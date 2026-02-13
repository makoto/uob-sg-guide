#!/usr/bin/env python3
"""
Fetch global datasets (OSMnx, Google Earth Engine, Mapillary/ZenSVI)
for Queenstown and save to docs/geo/global/.

Run with the zensvi conda env:
    /opt/anaconda3/envs/zensvi/bin/python3 fetch_global_layers.py
"""

import json
import math
import os
import statistics

from dotenv import load_dotenv

load_dotenv()

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = os.path.dirname(os.path.abspath(__file__))
GEO = os.path.join(BASE, "docs", "geo")
GLOBAL = os.path.join(GEO, "global")
BOUNDARY_PATH = os.path.join(GEO, "queenstown-boundary.geojson")
SUBZONES_PATH = os.path.join(GEO, "queenstown-subzones.geojson")

os.makedirs(GLOBAL, exist_ok=True)

# Queenstown bbox: [west, south, east, north]
BBOX = [103.7502, 1.2550, 103.8166, 1.3188]
WEST, SOUTH, EAST, NORTH = BBOX


def load_geojson(path):
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def write_geojson(path, fc):
    with open(path, "w", encoding="utf-8") as f:
        json.dump(fc, f, ensure_ascii=False)
    size_kb = os.path.getsize(path) / 1024
    print(f"  -> {os.path.basename(path)} ({size_kb:.0f} KB)")


# ===================================================================
# 1. OSMnx Street Network
# ===================================================================
def fetch_street_network():
    """Download drive+walk+bike street network and compute betweenness centrality."""
    import osmnx as ox
    import networkx as nx

    out_path = os.path.join(GLOBAL, "queenstown-street-network.geojson")
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

    out_path = os.path.join(GLOBAL, "queenstown-remote-sensing.geojson")

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
# 3. Mapillary / ZenSVI — Street-Level Imagery Metadata
# ===================================================================
def fetch_mapillary_metadata():
    """
    Download Mapillary image metadata for Queenstown bbox using ZenSVI.
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

    out_path = os.path.join(GLOBAL, "queenstown-mapillary-gvi.geojson")
    tmp_dir = os.path.join(BASE, "tmp_mapillary")
    os.makedirs(tmp_dir, exist_ok=True)

    meta_csv = os.path.join(tmp_dir, "mly_pids.csv")

    # Download metadata if not already cached
    if not os.path.exists(meta_csv):
        boundary_gdf = gpd.read_file(BOUNDARY_PATH)
        tmp_shp = os.path.join(tmp_dir, "boundary.shp")
        boundary_gdf.to_file(tmp_shp)

        print("  Downloading Mapillary metadata for Queenstown...")
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
    import sys

    tasks = sys.argv[1:] if len(sys.argv) > 1 else ["osmnx", "gee", "mapillary"]

    if "osmnx" in tasks:
        fetch_street_network()

    if "gee" in tasks:
        fetch_gee_rasters()

    if "mapillary" in tasks:
        fetch_mapillary_metadata()

    print("\nDone!")
