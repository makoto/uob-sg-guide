#!/usr/bin/env python3
"""
Fetch global datasets (OSMnx, Google Earth Engine, Mapillary/ZenSVI)
for a district and save to docs/geo/global/.

Usage:
    /opt/anaconda3/envs/zensvi/bin/python3 fetch_global_layers.py
    /opt/anaconda3/envs/zensvi/bin/python3 fetch_global_layers.py --district bishan
    /opt/anaconda3/envs/zensvi/bin/python3 fetch_global_layers.py --district bishan osmnx gee
"""

import argparse
import json
import math
import os
import statistics

from dotenv import load_dotenv

load_dotenv()

# ---------------------------------------------------------------------------
# Args
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="Fetch global layers for a district")
parser.add_argument("--district", default="queenstown",
                    help="District name (default: queenstown)")
parser.add_argument("tasks", nargs="*", default=["buildings", "osmnx", "gee", "mapillary"],
                    help="Tasks to run: buildings osmnx gee raster_pngs rs_grid mapillary")
args = parser.parse_args()

DISTRICT = args.district.lower().replace(" ", "-")

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = os.path.dirname(os.path.abspath(__file__))
GEO = os.path.join(BASE, "docs", "geo")
GLOBAL = os.path.join(GEO, "global")
BOUNDARY_PATH = os.path.join(GEO, f"{DISTRICT}-boundary.geojson")
SUBZONES_PATH = os.path.join(GEO, f"{DISTRICT}-subzones.geojson")

os.makedirs(GLOBAL, exist_ok=True)


def load_geojson(path):
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def write_geojson(path, fc):
    with open(path, "w", encoding="utf-8") as f:
        json.dump(fc, f, ensure_ascii=False)
    size_kb = os.path.getsize(path) / 1024
    print(f"  -> {os.path.basename(path)} ({size_kb:.0f} KB)")


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


# Derive bbox from boundary
BBOX = derive_bbox(BOUNDARY_PATH)
WEST, SOUTH, EAST, NORTH = BBOX
print(f"District: {DISTRICT}, BBOX: [{WEST:.4f}, {SOUTH:.4f}, {EAST:.4f}, {NORTH:.4f}]")


# ===================================================================
# 1. OSMnx Street Network
# ===================================================================
def fetch_street_network():
    """Download drive+walk+bike street network and compute betweenness centrality."""
    import osmnx as ox
    import networkx as nx

    out_path = os.path.join(GLOBAL, f"{DISTRICT}-street-network.geojson")
    print("\n=== OSMnx Street Network ===")

    # Download drive network (main streets for corridor analysis)
    print("  Downloading drive network...")
    G = ox.graph_from_bbox(
        north=NORTH, south=SOUTH, east=EAST, west=WEST,
        network_type="drive", simplify=True, truncate_by_edge=True,
    )
    print(f"  Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}")

    # Compute approximate edge betweenness centrality (sampled, fast)
    print("  Computing betweenness centrality (approximate, k=100)...")
    Gu = ox.convert.to_undirected(G)
    k = min(100, Gu.number_of_nodes())
    bc = nx.edge_betweenness_centrality(Gu, k=k, weight="length")

    # Map centrality back to directed edges
    edge_bc = {}
    for edge_key, val in bc.items():
        u, v = edge_key[0], edge_key[1]
        edge_bc[(u, v)] = val
        edge_bc[(v, u)] = val

    # Convert to GeoJSON via geopandas
    edges = ox.convert.graph_to_gdfs(G, nodes=False)

    # Add centrality column
    bc_vals = []
    for idx in edges.index:
        u, v, _ = idx
        bc_vals.append(edge_bc.get((u, v), 0.0))
    edges["betweenness"] = bc_vals

    # Keep useful columns, drop heavy ones
    keep_cols = ["geometry", "highway", "name", "length", "betweenness"]
    drop = [c for c in edges.columns if c not in keep_cols]
    edges = edges.drop(columns=drop, errors="ignore")

    # Round values
    edges["length"] = edges["length"].round(1)
    edges["betweenness"] = edges["betweenness"].round(6)

    # Write GeoJSON with reduced coordinate precision (5 decimal places ≈ 1m)
    import tempfile
    tmp_path = out_path + ".tmp"
    edges.to_file(tmp_path, driver="GeoJSON")

    # Reduce coordinate precision to shrink file
    with open(tmp_path, encoding="utf-8") as f:
        gj = json.load(f)
    os.remove(tmp_path)

    def round_coords(coords):
        if isinstance(coords[0], (int, float)):
            return [round(c, 5) for c in coords]
        return [round_coords(c) for c in coords]

    for feat in gj["features"]:
        feat["geometry"]["coordinates"] = round_coords(feat["geometry"]["coordinates"])

    write_geojson(out_path, gj)

    size_kb = os.path.getsize(out_path) / 1024
    if size_kb > 5000:
        print("  WARNING: File > 5 MB, may need further simplification")

    return out_path


# ===================================================================
# 2. Google Earth Engine — Remote Sensing Zonal Statistics
# ===================================================================
def fetch_gee_rasters():
    """
    Compute zonal statistics per subzone from multiple GEE raster datasets:
    - Landsat 9 LST
    - Sentinel-2 NDVI & NDBI
    - GHSL Built-up Height
    - ALOS AW3D30 DSM
    - SRTM DEM
    - Meta Canopy Height
    """
    import ee
    import geopandas as gpd

    print("\n=== Google Earth Engine Remote Sensing ===")
    ee.Initialize(project="uobdubai")
    print("  GEE initialized")

    out_path = os.path.join(GLOBAL, f"{DISTRICT}-remote-sensing.geojson")

    # Load subzones
    subzones_gj = load_geojson(SUBZONES_PATH)
    subzones_gdf = gpd.read_file(SUBZONES_PATH)

    # Convert subzones to ee.FeatureCollection
    ee_features = []
    for feat in subzones_gj["features"]:
        props = {"subzone_name": feat["properties"]["SUBZONE_N"]}
        geom = feat["geometry"]
        if geom["type"] == "Polygon":
            ee_geom = ee.Geometry.Polygon(geom["coordinates"])
        elif geom["type"] == "MultiPolygon":
            ee_geom = ee.Geometry.MultiPolygon(geom["coordinates"])
        else:
            continue
        ee_features.append(ee.Feature(ee_geom, props))
    ee_subzones = ee.FeatureCollection(ee_features)

    # Define AOI
    aoi = ee.Geometry.Rectangle([WEST, SOUTH, EAST, NORTH])

    # --- Landsat 9 LST ---
    print("  Computing LST from Landsat 9...")
    l9 = (
        ee.ImageCollection("LANDSAT/LC09/C02/T1_L2")
        .filterBounds(aoi)
        .filterDate("2023-01-01", "2025-01-01")
        .filter(ee.Filter.lt("CLOUD_COVER", 20))
    )
    # Thermal band ST_B10: scale 0.00341802, offset 149.0 → Kelvin
    lst_composite = l9.median()
    lst_kelvin = lst_composite.select("ST_B10").multiply(0.00341802).add(149.0)
    lst_celsius = lst_kelvin.subtract(273.15).rename("lst_c")

    # --- Sentinel-2 NDVI & NDBI ---
    print("  Computing NDVI & NDBI from Sentinel-2...")
    s2 = (
        ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
        .filterBounds(aoi)
        .filterDate("2023-01-01", "2025-01-01")
        .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", 20))
    )
    s2_median = s2.median()
    ndvi = s2_median.normalizedDifference(["B8", "B4"]).rename("ndvi")
    ndbi = s2_median.normalizedDifference(["B11", "B8"]).rename("ndbi")

    # --- GHSL Built-up Height ---
    print("  Computing GHSL built-up height...")
    ghsl = ee.Image("JRC/GHSL/P2023A/GHS_BUILT_H/2018").select("built_height")

    # --- ALOS AW3D30 DSM ---
    print("  Computing ALOS DSM...")
    dsm = ee.ImageCollection("JAXA/ALOS/AW3D30/V4_1").mosaic().select("DSM")

    # --- SRTM DEM ---
    print("  Computing SRTM DEM...")
    dem = ee.Image("USGS/SRTMGL1_003").select("elevation")

    # --- Hansen Global Forest Change (tree canopy cover) ---
    print("  Computing tree canopy cover (Hansen)...")
    hansen = ee.Image("UMD/hansen/global_forest_change_2024_v1_12")
    canopy_pct = hansen.select("treecover2000").rename("canopy_cover_pct")
    # Loss mask: where tree loss occurred
    canopy_loss = hansen.select("loss").rename("canopy_loss")

    # Stack all bands
    stack = (
        lst_celsius
        .addBands(ndvi)
        .addBands(ndbi)
        .addBands(ghsl.rename("ghsl_height"))
        .addBands(dsm.rename("dsm"))
        .addBands(dem.rename("dem"))
        .addBands(canopy_pct)
        .addBands(canopy_loss)
    )

    # Reduce regions — mean per subzone
    print("  Running zonal statistics (reduceRegions)...")
    reduced = stack.reduceRegions(
        collection=ee_subzones,
        reducer=ee.Reducer.mean(),
        scale=30,  # 30m resolution
    )

    # Fetch results
    results = reduced.getInfo()
    print(f"  Got {len(results['features'])} subzone results")

    # Build output: merge zonal stats into subzone GeoJSON
    # Create lookup by subzone name
    stats_by_name = {}
    for feat in results["features"]:
        name = feat["properties"]["subzone_name"]
        stats_by_name[name] = feat["properties"]

    # Build output features using original subzone geometries
    out_features = []
    for feat in subzones_gj["features"]:
        name = feat["properties"]["SUBZONE_N"]
        stats = stats_by_name.get(name, {})

        props = {
            "subzone_name": name,
            "subzone_code": feat["properties"]["SUBZONE_C"],
            "lst_mean_c": round(stats.get("lst_c", 0) or 0, 2),
            "ndvi_mean": round(stats.get("ndvi", 0) or 0, 4),
            "ndbi_mean": round(stats.get("ndbi", 0) or 0, 4),
            "ghsl_height_mean": round(stats.get("ghsl_height", 0) or 0, 2),
            "dsm_mean": round(stats.get("dsm", 0) or 0, 2),
            "dem_mean": round(stats.get("dem", 0) or 0, 2),
            "canopy_cover_pct": round(stats.get("canopy_cover_pct", 0) or 0, 2),
            "canopy_loss_pct": round((stats.get("canopy_loss", 0) or 0) * 100, 2),
        }

        out_features.append({
            "type": "Feature",
            "properties": props,
            "geometry": feat["geometry"],
        })

    fc = {"type": "FeatureCollection", "features": out_features}
    write_geojson(out_path, fc)
    return out_path


# ===================================================================
# 2b. Google Earth Engine — Raster PNG Exports
# ===================================================================
def generate_raster_pngs():
    """
    Export 4 raster overlay PNGs (NDVI, LST, NDBI, canopy) from GEE
    using getThumbURL(). Each is 1024px wide, georeferenced to the
    district bounding box for use as a SingleTileImageryProvider overlay.
    """
    import ee
    import urllib.request

    print("\n=== GEE Raster PNG Export ===")
    ee.Initialize(project="uobdubai")
    print("  GEE initialized")

    raster_dir = os.path.join(GLOBAL, "rasters")
    os.makedirs(raster_dir, exist_ok=True)

    aoi = ee.Geometry.Rectangle([WEST, SOUTH, EAST, NORTH])

    # --- Build raster layers (same as fetch_gee_rasters) ---
    print("  Building Sentinel-2 NDVI & NDBI composites...")
    s2 = (
        ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
        .filterBounds(aoi)
        .filterDate("2023-01-01", "2025-01-01")
        .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", 20))
    )
    s2_median = s2.median()
    ndvi = s2_median.normalizedDifference(["B8", "B4"]).rename("ndvi")
    ndbi = s2_median.normalizedDifference(["B11", "B8"]).rename("ndbi")

    print("  Building Landsat 9 LST composite...")
    l9 = (
        ee.ImageCollection("LANDSAT/LC09/C02/T1_L2")
        .filterBounds(aoi)
        .filterDate("2023-01-01", "2025-01-01")
        .filter(ee.Filter.lt("CLOUD_COVER", 20))
    )
    lst_composite = l9.median()
    lst_kelvin = lst_composite.select("ST_B10").multiply(0.00341802).add(149.0)
    lst_celsius = lst_kelvin.subtract(273.15).rename("lst_c")

    print("  Building Hansen canopy cover...")
    hansen = ee.Image("UMD/hansen/global_forest_change_2024_v1_12")
    canopy = hansen.select("treecover2000").rename("canopy")

    # --- Export each as PNG via getThumbURL ---
    layers = [
        ("ndvi", ndvi, {"min": -0.2, "max": 0.8, "palette": ["d73027", "fc8d59", "fee08b", "d9ef8b", "91cf60", "1a9850"]}),
        ("lst", lst_celsius, {"min": 25, "max": 40, "palette": ["313695", "74add1", "fed976", "fd8d3c", "e31a1c", "800026"]}),
        ("ndbi", ndbi, {"min": -0.3, "max": 0.3, "palette": ["1a9850", "91cf60", "d9ef8b", "fee08b", "fc8d59", "d73027"]}),
        ("canopy", canopy, {"min": 0, "max": 100, "palette": ["f7fcb1", "addd8e", "78c679", "31a354", "006837"]}),
    ]

    for name, image, vis in layers:
        out_path = os.path.join(raster_dir, f"{DISTRICT}-{name}.png")
        print(f"  Exporting {name}...")
        url = image.getThumbURL({
            "min": vis["min"],
            "max": vis["max"],
            "palette": vis["palette"],
            "region": aoi,
            "dimensions": 1024,
            "format": "png",
        })
        urllib.request.urlretrieve(url, out_path)
        size_kb = os.path.getsize(out_path) / 1024
        print(f"    -> {DISTRICT}-{name}.png ({size_kb:.0f} KB)")

    print("  All 4 raster PNGs exported")


# ===================================================================
# 2c. Google Earth Engine — 100m RS Grid
# ===================================================================
def generate_rs_grid():
    """
    Create a 100m fishnet grid over the district boundary, compute
    per-cell NDVI, NDBI, LST, and canopy cover via GEE reduceRegions,
    and output as GeoJSON.
    """
    import ee
    from shapely.geometry import shape, box, mapping
    from shapely.ops import unary_union

    print("\n=== GEE 100m RS Grid ===")
    ee.Initialize(project="uobdubai")
    print("  GEE initialized")

    out_path = os.path.join(GLOBAL, f"{DISTRICT}-rs-grid.geojson")

    # Load boundary
    boundary_gj = load_geojson(BOUNDARY_PATH)
    boundary_geom = unary_union([
        shape(f["geometry"]) for f in boundary_gj["features"]
    ])

    # Create 100m fishnet grid (~0.0009° at equator)
    cell_size = 0.0009  # degrees, ~100m
    minx, miny, maxx, maxy = boundary_geom.bounds
    print(f"  Boundary bounds: {minx:.4f}, {miny:.4f}, {maxx:.4f}, {maxy:.4f}")

    cells = []
    x = minx
    while x < maxx:
        y = miny
        while y < maxy:
            cell = box(x, y, x + cell_size, y + cell_size)
            clipped = cell.intersection(boundary_geom)
            if not clipped.is_empty and clipped.area > 0:
                cells.append(clipped)
            y += cell_size
        x += cell_size
    print(f"  Generated {len(cells)} grid cells")

    # Build raster stack
    aoi = ee.Geometry.Rectangle([WEST, SOUTH, EAST, NORTH])

    print("  Building raster composites...")
    s2 = (
        ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
        .filterBounds(aoi)
        .filterDate("2023-01-01", "2025-01-01")
        .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", 20))
    )
    s2_median = s2.median()
    ndvi = s2_median.normalizedDifference(["B8", "B4"]).rename("ndvi")
    ndbi = s2_median.normalizedDifference(["B11", "B8"]).rename("ndbi")

    l9 = (
        ee.ImageCollection("LANDSAT/LC09/C02/T1_L2")
        .filterBounds(aoi)
        .filterDate("2023-01-01", "2025-01-01")
        .filter(ee.Filter.lt("CLOUD_COVER", 20))
    )
    lst_composite = l9.median()
    lst_kelvin = lst_composite.select("ST_B10").multiply(0.00341802).add(149.0)
    lst_celsius = lst_kelvin.subtract(273.15).rename("lst_c")

    hansen = ee.Image("UMD/hansen/global_forest_change_2024_v1_12")
    canopy = hansen.select("treecover2000").rename("canopy_pct")

    stack = ndvi.addBands(ndbi).addBands(lst_celsius).addBands(canopy)

    # Process in batches (GEE has limits on FeatureCollection size)
    batch_size = 500
    all_results = []
    for i in range(0, len(cells), batch_size):
        batch = cells[i:i + batch_size]
        print(f"  Processing cells {i+1}-{i+len(batch)} of {len(cells)}...")

        ee_features = []
        for j, cell in enumerate(batch):
            geom = cell
            if geom.geom_type == "Polygon":
                ee_geom = ee.Geometry.Polygon(list(mapping(geom)["coordinates"]))
            elif geom.geom_type == "MultiPolygon":
                ee_geom = ee.Geometry.MultiPolygon(list(mapping(geom)["coordinates"]))
            else:
                continue
            ee_features.append(ee.Feature(ee_geom, {"cell_id": i + j}))

        ee_fc = ee.FeatureCollection(ee_features)
        reduced = stack.reduceRegions(
            collection=ee_fc,
            reducer=ee.Reducer.mean(),
            scale=10,  # 10m for Sentinel-2 resolution
        )
        results = reduced.getInfo()
        all_results.extend(results["features"])

    # Build output GeoJSON
    print(f"  Got {len(all_results)} cell results")
    out_features = []
    for k, (cell, result) in enumerate(zip(cells, all_results)):
        props = result["properties"]
        out_features.append({
            "type": "Feature",
            "properties": {
                "cell_id": k,
                "ndvi": round(props.get("ndvi") or 0, 4) if props.get("ndvi") is not None else None,
                "ndbi": round(props.get("ndbi") or 0, 4) if props.get("ndbi") is not None else None,
                "lst_c": round(props.get("lst_c") or 0, 2) if props.get("lst_c") is not None else None,
                "canopy_pct": round(props.get("canopy_pct") or 0, 2) if props.get("canopy_pct") is not None else None,
            },
            "geometry": mapping(cell),
        })

    # Round coordinates to 5 decimal places
    def round_coords(coords):
        if isinstance(coords[0], (int, float)):
            return [round(c, 5) for c in coords]
        return [round_coords(c) for c in coords]

    for feat in out_features:
        feat["geometry"]["coordinates"] = round_coords(feat["geometry"]["coordinates"])

    fc = {"type": "FeatureCollection", "features": out_features}
    write_geojson(out_path, fc)
    print(f"  {len(out_features)} grid cells written")
    return out_path


# ===================================================================
# 3. OSM Buildings + HDB Enrichment
# ===================================================================
def fetch_buildings():
    """
    Download building footprints from OpenStreetMap via osmnx and enrich
    with HDB property data. Produces 3 output files:
      - {district}-buildings.geojson (all buildings)
      - {district}-buildings-osm-only.geojson (non-HDB)
      - {district}-buildings-hdb-enriched.geojson (HDB only)
    """
    import csv as csv_module
    import osmnx as ox
    import geopandas as gpd
    from shapely.geometry import mapping, shape

    print("\n=== OSM Buildings + HDB Enrichment ===")

    out_all = os.path.join(GLOBAL, f"{DISTRICT}-buildings.geojson")
    out_osm = os.path.join(GLOBAL, f"{DISTRICT}-buildings-osm-only.geojson")
    out_hdb = os.path.join(GLOBAL, f"{DISTRICT}-buildings-hdb-enriched.geojson")

    # Download building footprints from OSM
    print("  Downloading building footprints from OSM...")
    tags = {"building": True}
    buildings = ox.features_from_bbox(
        bbox=(NORTH, SOUTH, EAST, WEST), tags=tags
    )
    # Keep only polygons/multipolygons (not nodes)
    buildings = buildings[buildings.geometry.type.isin(["Polygon", "MultiPolygon"])]
    print(f"  {len(buildings)} building footprints downloaded")

    # Load boundary for clipping
    boundary_gj = load_geojson(BOUNDARY_PATH)
    boundary_geom = shape(boundary_gj["features"][0]["geometry"])

    # Clip to boundary
    buildings_clipped = buildings[buildings.geometry.intersects(boundary_geom)]
    print(f"  {len(buildings_clipped)} buildings within boundary")

    # Load HDB property data for enrichment
    hdb_csv = os.path.join(BASE, "data", "data-gov-sg", "hdb-property.csv")
    hdb_lookup = {}  # (block, street_upper) -> row
    if os.path.exists(hdb_csv):
        with open(hdb_csv, encoding="utf-8") as f:
            reader = csv_module.DictReader(f)
            for row in reader:
                key = (row["blk_no"].strip(), row["street"].strip().upper())
                hdb_lookup[key] = row
        print(f"  Loaded {len(hdb_lookup)} HDB property records")
    else:
        print(f"  WARNING: {hdb_csv} not found, skipping HDB enrichment")

    # Street abbreviation map for matching
    STREET_ABBREVS = {
        "AVE": "AVENUE", "CL": "CLOSE", "CRES": "CRESCENT", "CT": "COURT",
        "DR": "DRIVE", "HTS": "HEIGHTS", "LN": "LANE", "PL": "PLACE",
        "RD": "ROAD", "ST": "STREET", "TER": "TERRACE",
    }

    def normalise_street(name):
        parts = name.strip().upper().split()
        return " ".join(STREET_ABBREVS.get(p, p) for p in parts)

    # Compute height from building_levels or GHSL default
    DEFAULT_FLOOR_HEIGHT = 3.0  # metres per floor
    DEFAULT_HEIGHT = 10.0

    # Build feature list
    all_features = []
    hdb_features = []
    osm_features = []

    for idx, row in buildings_clipped.iterrows():
        geom = row.geometry
        osm_id = idx[1] if isinstance(idx, tuple) else idx

        # Extract OSM properties
        name = str(row.get("name", "") or "")
        building_type = str(row.get("building", "") or "")
        height_raw = row.get("height", None)
        levels_raw = row.get("building:levels", None)
        min_level = str(row.get("building:min_level", "") or "")
        roof_levels = str(row.get("roof:levels", "") or "")
        addr_num = str(row.get("addr:housenumber", "") or "")
        addr_street = str(row.get("addr:street", "") or "")
        addr_postcode = str(row.get("addr:postcode", "") or "")
        amenity = str(row.get("amenity", "") or "")

        # Compute height
        height_m = DEFAULT_HEIGHT
        height_source = "default"
        if height_raw and str(height_raw) not in ("", "nan", "None", "0"):
            try:
                height_m = float(str(height_raw).replace("m", "").strip())
                height_source = "osm_height"
            except ValueError:
                pass
        elif levels_raw and str(levels_raw) not in ("", "nan", "None", "0"):
            try:
                height_m = float(str(levels_raw)) * DEFAULT_FLOOR_HEIGHT
                height_source = "osm_levels"
            except ValueError:
                pass

        # Try HDB match
        hdb_match = False
        hdb_props = {}
        if addr_num and addr_street:
            osm_street_upper = addr_street.strip().upper()
            # Try direct match
            hdb_row = hdb_lookup.get((addr_num.strip(), osm_street_upper))
            if not hdb_row:
                # Try normalised street
                norm = normalise_street(addr_street)
                hdb_row = hdb_lookup.get((addr_num.strip(), norm))
            if hdb_row:
                hdb_match = True
                max_floor = int(hdb_row.get("max_floor_lvl", 0) or 0)
                if max_floor > 0 and height_source == "default":
                    height_m = max_floor * DEFAULT_FLOOR_HEIGHT
                    height_source = "hdb_floors"
                hdb_props = {
                    "hdb_year_completed": hdb_row.get("year_completed", ""),
                    "hdb_total_dwelling_units": hdb_row.get("total_dwelling_units", ""),
                    "hdb_max_floor_lvl": hdb_row.get("max_floor_lvl", ""),
                    "hdb_residential": hdb_row.get("residential", ""),
                    "hdb_commercial": hdb_row.get("commercial", ""),
                    "hdb_market_hawker": hdb_row.get("market_hawker", ""),
                }

        props = {
            "osm_id": int(osm_id) if not isinstance(osm_id, str) else osm_id,
            "building": building_type,
            "name": name,
            "height": str(height_raw) if height_raw and str(height_raw) not in ("nan", "None") else "0",
            "building_levels": str(levels_raw) if levels_raw and str(levels_raw) not in ("nan", "None") else "",
            "building_min_level": min_level,
            "roof_levels": roof_levels,
            "addr_housenumber": addr_num,
            "addr_street": addr_street,
            "addr_postcode": addr_postcode,
            "amenity": amenity,
            "residential": "yes" if building_type in ("residential", "apartments", "house") else "",
            "data_source": "osm+hdb" if hdb_match else "osm",
            "hdb_match": hdb_match,
            "height_m": round(height_m, 1),
            "height_source": height_source,
            **hdb_props,
        }

        # Fill missing HDB fields for non-HDB buildings
        if not hdb_match:
            for k in ["hdb_year_completed", "hdb_total_dwelling_units", "hdb_max_floor_lvl",
                       "hdb_residential", "hdb_commercial", "hdb_market_hawker"]:
                props.setdefault(k, "")

        feature = {
            "type": "Feature",
            "properties": props,
            "geometry": json.loads(gpd.GeoSeries([geom]).to_json())["features"][0]["geometry"],
        }

        all_features.append(feature)
        if hdb_match:
            hdb_features.append(feature)
        else:
            osm_features.append(feature)

    print(f"  Total: {len(all_features)} buildings ({len(hdb_features)} HDB, {len(osm_features)} OSM-only)")

    # Write outputs
    def write_fc(path, features):
        fc = {"type": "FeatureCollection", "features": features}
        with open(path, "w", encoding="utf-8") as f:
            json.dump(fc, f, ensure_ascii=False)
        size_kb = os.path.getsize(path) / 1024
        print(f"  -> {os.path.basename(path)}: {len(features)} features ({size_kb:.0f} KB)")

    write_fc(out_all, all_features)
    write_fc(out_osm, osm_features)
    write_fc(out_hdb, hdb_features)

    return out_all


# ===================================================================
# 4. Mapillary / ZenSVI — Street-Level Imagery Metadata
# ===================================================================
def fetch_mapillary_metadata():
    """
    Download Mapillary image metadata for district bbox using ZenSVI.
    Saves sampled point locations with image IDs (GVI computed in follow-up).
    """
    import csv as csv_module
    import random

    from zensvi.download import MLYDownloader
    import geopandas as gpd

    print("\n=== Mapillary / ZenSVI ===")
    api_key = os.environ.get("YOUR_OWN_MLY_API_KEY", "")
    if not api_key:
        print("  WARNING: No Mapillary API key found in .env, skipping")
        return None

    out_path = os.path.join(GLOBAL, f"{DISTRICT}-mapillary-gvi.geojson")
    tmp_dir = os.path.join(BASE, f"tmp_mapillary_{DISTRICT}")
    os.makedirs(tmp_dir, exist_ok=True)

    meta_csv = os.path.join(tmp_dir, "mly_pids.csv")

    # Download metadata if not already cached
    if not os.path.exists(meta_csv):
        boundary_gdf = gpd.read_file(BOUNDARY_PATH)
        tmp_shp = os.path.join(tmp_dir, "boundary.shp")
        boundary_gdf.to_file(tmp_shp)

        print(f"  Downloading Mapillary metadata for {DISTRICT}...")
        downloader = MLYDownloader(mly_api_key=api_key)
        downloader.download_svi(
            dir_output=tmp_dir,
            input_shp_file=tmp_shp,
            metadata_only=True,
            start_date="2020-01-01",
            end_date="2025-12-31",
        )
    else:
        print("  Using cached metadata CSV")

    if not os.path.exists(meta_csv):
        print("  WARNING: No metadata CSV found after download")
        return None

    # Read all rows
    rows = []
    with open(meta_csv, encoding="utf-8") as f:
        reader = csv_module.DictReader(f)
        for row in reader:
            rows.append(row)
    print(f"  {len(rows)} total image records")

    # Downsample to ~5000 points for reasonable GeoJSON size
    random.seed(42)
    sample_size = min(5000, len(rows))
    sampled = random.sample(rows, sample_size)
    print(f"  Sampled {sample_size} points")

    features = []
    for row in sampled:
        lon = float(row["lon"])
        lat = float(row["lat"])
        props = {
            "image_id": row["id"],
            "is_pano": row.get("is_pano") == "True",
        }
        features.append({
            "type": "Feature",
            "properties": props,
            "geometry": {"type": "Point", "coordinates": [round(lon, 5), round(lat, 5)]},
        })

    fc = {"type": "FeatureCollection", "features": features}
    write_geojson(out_path, fc)
    return out_path


# ===================================================================
# Main
# ===================================================================
if __name__ == "__main__":
    tasks = args.tasks

    if "buildings" in tasks:
        fetch_buildings()

    if "osmnx" in tasks:
        fetch_street_network()

    if "gee" in tasks:
        fetch_gee_rasters()

    if "raster_pngs" in tasks:
        generate_raster_pngs()

    if "rs_grid" in tasks:
        generate_rs_grid()

    if "mapillary" in tasks:
        fetch_mapillary_metadata()

    print("\nDone!")
