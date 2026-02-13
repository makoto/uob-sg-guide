#!/usr/bin/env python3
"""
Generate walkability scores per subzone using BEH-NWI framework.

Components (Rundle et al. 2019):
  1. Population density — from existing subzone summary CSV
  2. Intersection density — walk network nodes with degree >= 3
  3. Transit access — network distance to nearest MRT exit (inverted)
  4. Destination accessibility — weighted avg network distance to 9 POI categories

Usage:
    /opt/anaconda3/envs/zensvi/bin/python3 generate_walkability.py
    /opt/anaconda3/envs/zensvi/bin/python3 generate_walkability.py --district bishan
"""

import argparse
import json
import math
import os
import statistics

import networkx as nx
import osmnx as ox
from shapely.geometry import shape, Point

# ---------------------------------------------------------------------------
# Args
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="Generate walkability scores for a district")
parser.add_argument("--district", default="queenstown",
                    help="District name (default: queenstown)")
args = parser.parse_args()

DISTRICT = args.district.lower().replace(" ", "-")

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = os.path.dirname(os.path.abspath(__file__))
GEO = os.path.join(BASE, "docs", "geo")
GLOBAL = os.path.join(GEO, "global")
GOV_SG = os.path.join(GEO, "gov-sg")

SUBZONES_PATH = os.path.join(GEO, f"{DISTRICT}-subzones.geojson")
BOUNDARY_PATH = os.path.join(GEO, f"{DISTRICT}-boundary.geojson")
SUMMARY_CSV = os.path.join(GEO, f"{DISTRICT}-subzone-summary.csv")
OUT_PATH = os.path.join(GLOBAL, f"{DISTRICT}-walkability.json")
GRID_PATH = os.path.join(GLOBAL, f"{DISTRICT}-walkability-grid.geojson")
WALK_NET_PATH = os.path.join(GLOBAL, f"{DISTRICT}-walk-network.geojson")

# POI layers with weights for destination accessibility
POI_LAYERS = {
    "hawker_centres":   {"file": f"{DISTRICT}-hawker-centres.geojson",   "weight": 3},
    "supermarkets":     {"file": f"{DISTRICT}-supermarkets.geojson",     "weight": 3},
    "mrt_exits":        {"file": f"{DISTRICT}-mrt-exits.geojson",        "weight": 2},
    "parks":            {"file": f"{DISTRICT}-parks.geojson",            "weight": 2},
    "community_clubs":  {"file": f"{DISTRICT}-community-clubs.geojson",  "weight": 2},
    "chas_clinics":     {"file": f"{DISTRICT}-chas-clinics.geojson",     "weight": 2},
    "preschools":       {"file": f"{DISTRICT}-preschools.geojson",       "weight": 1},
    "gyms":             {"file": f"{DISTRICT}-gyms.geojson",             "weight": 1},
    "park_facilities":  {"file": f"{DISTRICT}-park-facilities.geojson",  "weight": 1},
}


def load_geojson(path):
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def polygon_area_km2(geom):
    """Rough area in km² using projected Shoelace formula."""
    if geom.geom_type == "MultiPolygon":
        return sum(polygon_area_km2(p) for p in geom.geoms)
    coords = list(geom.exterior.coords)
    ref_lat = sum(c[1] for c in coords) / len(coords)
    m_lng = 111_320 * math.cos(math.radians(ref_lat))
    m_lat = 111_320
    projected = [(c[0] * m_lng, c[1] * m_lat) for c in coords]
    n = len(projected)
    area = 0
    for i in range(n):
        j = (i + 1) % n
        area += projected[i][0] * projected[j][1]
        area -= projected[j][0] * projected[i][1]
    return abs(area) / 2.0 / 1e6


def derive_bbox(boundary_path):
    """Derive bounding box [west, south, east, north] from boundary GeoJSON."""
    gj = load_geojson(boundary_path)
    all_coords = []
    for feat in gj["features"]:
        geom = feat["geometry"]
        if geom["type"] == "Polygon":
            for ring in geom["coordinates"]:
                all_coords.extend(ring)
        elif geom["type"] == "MultiPolygon":
            for poly in geom["coordinates"]:
                for ring in poly:
                    all_coords.extend(ring)
    lons = [c[0] for c in all_coords]
    lats = [c[1] for c in all_coords]
    return [min(lons), min(lats), max(lons), max(lats)]


def main():
    print(f"=== Walkability Scoring ({DISTRICT}) ===\n")

    # Derive bbox from boundary
    BBOX = derive_bbox(BOUNDARY_PATH)
    WEST, SOUTH, EAST, NORTH = BBOX
    print(f"  BBOX: [{WEST:.4f}, {SOUTH:.4f}, {EAST:.4f}, {NORTH:.4f}]")

    # 1. Download walk network
    print("Downloading walk network via OSMnx...")
    G = ox.graph_from_bbox(
        bbox=(NORTH, SOUTH, EAST, WEST),
        network_type="walk", simplify=True, truncate_by_edge=True,
    )
    print(f"  Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}")

    # Use largest connected component
    G_undirected = ox.convert.to_undirected(G)
    largest_cc = max(nx.connected_components(G_undirected), key=len)
    G = G.subgraph(largest_cc).copy()
    print(f"  Largest component: {G.number_of_nodes()} nodes")

    # Project to UTM for metre distances
    G_proj = ox.project_graph(G)

    # 2. Load subzones
    print("Loading subzones...")
    subzones_gj = load_geojson(SUBZONES_PATH)
    subzones = []
    for feat in subzones_gj["features"]:
        geom = shape(feat["geometry"])
        subzones.append({
            "name": feat["properties"]["SUBZONE_N"],
            "geom": geom,
            "area_km2": polygon_area_km2(geom),
        })
    print(f"  {len(subzones)} subzones")

    # 3. Map walk network nodes to subzones (using WGS84 graph for containment)
    print("Mapping nodes to subzones...")
    node_subzone = {}  # node_id -> subzone index
    for node_id, data in G.nodes(data=True):
        pt = Point(data["x"], data["y"])
        for i, sz in enumerate(subzones):
            if sz["geom"].contains(pt):
                node_subzone[node_id] = i
                break

    for i, sz in enumerate(subzones):
        sz["node_ids"] = [n for n, si in node_subzone.items() if si == i]
        print(f"  {sz['name']}: {len(sz['node_ids'])} nodes")

    # 4. Intersection density (nodes with degree >= 3)
    print("\nComputing intersection density...")
    for sz in subzones:
        intersections = sum(1 for n in sz["node_ids"] if G.degree(n) >= 3)
        area = sz["area_km2"]
        sz["intersection_density"] = round(intersections / area, 1) if area > 0 else None

    # 5. Load POIs and snap to nearest walk network node
    print("Loading POIs and snapping to walk network...")
    poi_nodes = {}  # category -> list of projected node ids
    for cat, info in POI_LAYERS.items():
        path = os.path.join(GOV_SG, info["file"])
        if not os.path.exists(path):
            print(f"  WARNING: {path} not found, skipping")
            poi_nodes[cat] = []
            continue
        gj = load_geojson(path)
        nodes = []
        for feat in gj["features"]:
            geom = feat["geometry"]
            if geom is None or geom["type"] != "Point":
                continue
            lon, lat = geom["coordinates"][:2]
            nearest = ox.nearest_nodes(G, lon, lat)
            if nearest in largest_cc:
                nodes.append(nearest)
        poi_nodes[cat] = list(set(nodes))
        print(f"  {cat}: {len(poi_nodes[cat])} POIs snapped")

    # 6. Multi-source Dijkstra for each POI category
    print("\nRunning shortest-path computations...")
    # distances[cat][node] = shortest distance to nearest POI in category
    distances = {}
    for cat, sources in poi_nodes.items():
        if not sources:
            distances[cat] = {}
            continue
        # Use projected graph for metre distances
        proj_sources = [n for n in sources if n in G_proj]
        if not proj_sources:
            distances[cat] = {}
            continue
        dist_dict = nx.multi_source_dijkstra_path_length(
            G_proj.to_undirected(), proj_sources, weight="length"
        )
        distances[cat] = dist_dict
        print(f"  {cat}: computed distances for {len(dist_dict)} nodes")

    # 7. Compute per-subzone metrics
    print("\nComputing per-subzone metrics...")

    # Transit access: median network distance to nearest MRT exit
    for sz in subzones:
        nodes = sz["node_ids"]
        if not nodes or not distances.get("mrt_exits"):
            sz["transit_access_score"] = None
            continue
        dists = [distances["mrt_exits"][n] for n in nodes if n in distances["mrt_exits"]]
        if dists:
            median_dist = statistics.median(dists)
            sz["transit_access_score"] = round(1000.0 / median_dist, 4) if median_dist > 0 else None
        else:
            sz["transit_access_score"] = None

    # Destination accessibility: weighted average of 1/median_dist across categories
    for sz in subzones:
        nodes = sz["node_ids"]
        if not nodes:
            sz["destination_accessibility"] = None
            continue
        weighted_sum = 0.0
        total_weight = 0
        for cat, info in POI_LAYERS.items():
            if not distances.get(cat):
                continue
            dists = [distances[cat][n] for n in nodes if n in distances[cat]]
            if not dists:
                continue
            median_dist = statistics.median(dists)
            if median_dist > 0:
                weighted_sum += info["weight"] * (1000.0 / median_dist)
                total_weight += info["weight"]
        sz["destination_accessibility"] = round(weighted_sum / total_weight, 4) if total_weight > 0 else None

    # 8. Read population density from existing CSV
    print("Reading population density from CSV...")
    import csv
    pop_density = {}
    if os.path.exists(SUMMARY_CSV):
        with open(SUMMARY_CSV, encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                name = row["subzone_name"]
                val = row.get("population_density", "")
                pop_density[name] = float(val) if val and val not in ("", "None") else 0
    for sz in subzones:
        sz["population_density"] = pop_density.get(sz["name"], 0)

    # 9. Z-score standardise and compute composite walkability index
    print("Computing z-score composite walkability index...")
    components = ["population_density", "intersection_density",
                   "transit_access_score", "destination_accessibility"]

    # Filter subzones with all components available
    valid_subzones = [sz for sz in subzones if all(
        sz.get(c) is not None for c in components
    )]

    if len(valid_subzones) >= 2:
        for comp in components:
            vals = [sz[comp] for sz in valid_subzones]
            mean = statistics.mean(vals)
            stdev = statistics.stdev(vals) if len(vals) > 1 else 1
            for sz in valid_subzones:
                sz[f"{comp}_z"] = (sz[comp] - mean) / stdev if stdev > 0 else 0

        for sz in valid_subzones:
            sz["walkability_index"] = round(
                sum(sz[f"{c}_z"] for c in components), 4
            )
    else:
        print("  WARNING: Not enough valid subzones for z-score")

    # Set walkability_index to None for subzones missing components
    for sz in subzones:
        if sz not in valid_subzones:
            sz["walkability_index"] = None

    # 10. Write output JSON
    print(f"\nWriting {OUT_PATH}...")
    output = []
    for sz in subzones:
        output.append({
            "subzone_name": sz["name"],
            "intersection_density": sz.get("intersection_density"),
            "transit_access_score": sz.get("transit_access_score"),
            "destination_accessibility": sz.get("destination_accessibility"),
            "walkability_index": sz.get("walkability_index"),
        })
    with open(OUT_PATH, "w", encoding="utf-8") as f:
        json.dump(output, f, indent=2, ensure_ascii=False)
    size_kb = os.path.getsize(OUT_PATH) / 1024
    print(f"  -> {os.path.basename(OUT_PATH)} ({size_kb:.1f} KB)")

    # 11. Generate 100m walkability grid (like RS grid)
    print("\nGenerating 100m walkability grid...")
    from scipy.spatial import KDTree
    import numpy as np
    from shapely.geometry import box as shapely_box

    boundary_gj = load_geojson(BOUNDARY_PATH)
    boundary_geom = shape(boundary_gj["features"][0]["geometry"])

    # Build KDTree of WGS84 node positions for fast nearest-node lookup
    all_node_ids = list(G.nodes())
    all_node_lnglat = np.array([[G.nodes[n]["x"], G.nodes[n]["y"]] for n in all_node_ids])
    node_tree = KDTree(all_node_lnglat)

    # Identify intersection nodes (degree >= 3) in projected coords
    intersection_ids = [n for n in G.nodes() if G.degree(n) >= 3]
    int_proj_coords = np.array(
        [[G_proj.nodes[n]["x"], G_proj.nodes[n]["y"]] for n in intersection_ids]
    ) if intersection_ids else np.empty((0, 2))
    int_tree = KDTree(int_proj_coords) if len(int_proj_coords) > 0 else None

    # Pre-compute per-node transit_access and dest_access
    node_transit = {}
    node_dest = {}
    for n in all_node_ids:
        d = distances.get("mrt_exits", {}).get(n)
        if d is not None and d > 0:
            node_transit[n] = round(1000.0 / d, 4)
        ws, tw = 0.0, 0
        for cat, info in POI_LAYERS.items():
            d = distances.get(cat, {}).get(n)
            if d is not None and d > 0:
                ws += info["weight"] * (1000.0 / d)
                tw += info["weight"]
        if tw > 0:
            node_dest[n] = round(ws / tw, 4)

    # Build subzone polygons for point-in-polygon lookup of population density
    sz_polys = []
    for sz in subzones:
        sz_polys.append({
            "geom": sz["geom"],
            "pop_density": sz.get("population_density", 0),
        })

    # Grid parameters: 100m cells — use center latitude of bbox
    center_lat = (SOUTH + NORTH) / 2
    M_LAT = 111_320
    M_LNG = 111_320 * math.cos(math.radians(center_lat))
    CELL_M = 100
    dlat = CELL_M / M_LAT
    dlng = CELL_M / M_LNG
    INT_RADIUS = 300  # metres for intersection density kernel

    grid_features = []
    cell_id = 0
    lat = SOUTH
    while lat < NORTH:
        lng = WEST
        while lng < EAST:
            cx, cy = lng + dlng / 2, lat + dlat / 2
            if not boundary_geom.contains(Point(cx, cy)):
                lng += dlng
                continue

            # Find nearest walk node
            _, idx = node_tree.query([cx, cy])
            nearest = all_node_ids[idx]

            # Transit access and destination access from nearest node
            transit = node_transit.get(nearest, 0)
            dest = node_dest.get(nearest, 0)

            # Intersection density: count degree>=3 nodes within 300m
            if int_tree is not None and nearest in G_proj.nodes:
                px, py = G_proj.nodes[nearest]["x"], G_proj.nodes[nearest]["y"]
                nearby = int_tree.query_ball_point([px, py], INT_RADIUS)
                area_km2 = math.pi * (INT_RADIUS / 1000) ** 2
                int_dens = round(len(nearby) / area_km2, 1)
            else:
                int_dens = 0

            # Population density from containing subzone
            cell_pop_density = 0
            cell_pt = Point(cx, cy)
            for sp in sz_polys:
                if sp["geom"].contains(cell_pt):
                    cell_pop_density = sp["pop_density"]
                    break

            coords = [
                [round(lng, 5), round(lat, 5)],
                [round(lng + dlng, 5), round(lat, 5)],
                [round(lng + dlng, 5), round(lat + dlat, 5)],
                [round(lng, 5), round(lat + dlat, 5)],
                [round(lng, 5), round(lat, 5)],
            ]

            grid_features.append({
                "type": "Feature",
                "properties": {
                    "cell_id": cell_id,
                    "pop_density": cell_pop_density,
                    "transit_access": transit,
                    "dest_access": dest,
                    "int_density": int_dens,
                },
                "geometry": {"type": "Polygon", "coordinates": [coords]},
            })
            cell_id += 1
            lng += dlng
        lat += dlat

    # Z-score composite walkability across grid cells (4 components, matching subzone method)
    grid_components = ["pop_density", "int_density", "transit_access", "dest_access"]
    for comp in grid_components:
        vals = [f["properties"][comp] for f in grid_features]
        nonzero = [v for v in vals if v > 0]
        if len(nonzero) >= 2:
            m = statistics.mean(vals)
            s = statistics.stdev(vals) if len(vals) > 1 else 1
            for f in grid_features:
                f["properties"][f"_{comp}_z"] = (f["properties"][comp] - m) / s if s > 0 else 0
        else:
            for f in grid_features:
                f["properties"][f"_{comp}_z"] = 0

    for f in grid_features:
        p = f["properties"]
        p["walkability"] = round(
            sum(p.get(f"_{c}_z", 0) for c in grid_components), 4
        )
        # Remove temp z-score fields
        for c in grid_components:
            del p[f"_{c}_z"]

    grid_fc = {"type": "FeatureCollection", "features": grid_features}
    with open(GRID_PATH, "w", encoding="utf-8") as f:
        json.dump(grid_fc, f, ensure_ascii=False)
    size_kb = os.path.getsize(GRID_PATH) / 1024
    print(f"  -> {os.path.basename(GRID_PATH)} ({size_kb:.0f} KB, {len(grid_features)} cells)")

    # 12. Optionally save walk network GeoJSON
    print("\nExporting walk network GeoJSON...")
    try:
        edges = ox.convert.graph_to_gdfs(G, nodes=False)
        keep_cols = ["geometry", "highway", "name", "length"]
        drop = [c for c in edges.columns if c not in keep_cols]
        edges = edges.drop(columns=drop, errors="ignore")
        if "length" in edges.columns:
            edges["length"] = edges["length"].round(1)

        import tempfile
        tmp_path = WALK_NET_PATH + ".tmp"
        edges.to_file(tmp_path, driver="GeoJSON")

        with open(tmp_path, encoding="utf-8") as f:
            gj = json.load(f)
        os.remove(tmp_path)

        def round_coords(coords):
            if isinstance(coords[0], (int, float)):
                return [round(c, 5) for c in coords]
            return [round_coords(c) for c in coords]

        for feat in gj["features"]:
            feat["geometry"]["coordinates"] = round_coords(feat["geometry"]["coordinates"])

        with open(WALK_NET_PATH, "w", encoding="utf-8") as f:
            json.dump(gj, f, ensure_ascii=False)

        size_kb = os.path.getsize(WALK_NET_PATH) / 1024
        print(f"  -> {os.path.basename(WALK_NET_PATH)} ({size_kb:.0f} KB)")
        if size_kb > 3000:
            print("  WARNING: Walk network > 3 MB, consider not committing")
            os.remove(WALK_NET_PATH)
            print("  Removed (too large)")
    except Exception as e:
        print(f"  Walk network export failed: {e}")

    # Print summary
    print(f"\n=== Walkability Scores ({DISTRICT}) ===")
    print(f"{'Subzone':<35} {'IntDens':>8} {'Transit':>8} {'DestAcc':>8} {'WalkIdx':>8}")
    print("-" * 75)
    for sz in sorted(subzones, key=lambda s: s.get("walkability_index") or -999, reverse=True):
        def fmt(v):
            return f"{v:8.2f}" if v is not None else "    null"
        print(f"{sz['name']:<35} {fmt(sz.get('intersection_density'))}"
              f" {fmt(sz.get('transit_access_score'))}"
              f" {fmt(sz.get('destination_accessibility'))}"
              f" {fmt(sz.get('walkability_index'))}")

    print("\nDone!")


if __name__ == "__main__":
    main()
