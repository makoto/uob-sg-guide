#!/usr/bin/env python3
"""
Generate per-subzone summary for a district.

Spatially joins multiple datasets (points, lines, CSV) to subzone polygons
and produces:
  - docs/geo/{district}-subzone-summary.geojson
  - docs/geo/{district}-subzone-summary.csv

Usage:
    python3 generate_subzone_summary.py                    # default: queenstown
    python3 generate_subzone_summary.py --district bishan
"""

import argparse
import csv
import json
import math
import os
import statistics

from shapely.geometry import shape, Point, LineString, MultiLineString
from shapely.ops import unary_union

# ---------------------------------------------------------------------------
# Parse args
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="Generate per-subzone summary for a district")
parser.add_argument("--district", default="queenstown",
                    help="District name (default: queenstown)")
args = parser.parse_args()

DISTRICT = args.district.lower().replace(" ", "-")

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = os.path.dirname(os.path.abspath(__file__))
GEO = os.path.join(BASE, "docs", "geo")
DATA = os.path.join(BASE, "data", "data-gov-sg")

SUBZONES_PATH = os.path.join(GEO, f"{DISTRICT}-subzones.geojson")
BUILDINGS_PATH = os.path.join(GEO, "global", f"{DISTRICT}-buildings.geojson")
REMOTE_SENSING_PATH = os.path.join(GEO, "global", f"{DISTRICT}-remote-sensing.geojson")
MAPILLARY_GVI_PATH = os.path.join(GEO, "global", f"{DISTRICT}-mapillary-gvi.geojson")
WALKABILITY_PATH = os.path.join(GEO, "global", f"{DISTRICT}-walkability.json")

POINT_LAYERS = {
    "hawker_centres": os.path.join(DATA, "hawker-centres.geojson"),
    "parks": os.path.join(DATA, "parks.geojson"),
    "mrt_exits": os.path.join(DATA, "mrt-station-exits.geojson"),
    "supermarkets": os.path.join(DATA, "supermarkets.geojson"),
    "gyms": os.path.join(DATA, "gyms.geojson"),
    "community_clubs": os.path.join(DATA, "community-clubs.geojson"),
    "sport_facilities": os.path.join(DATA, "sport-facilities.geojson"),
}

LINE_LAYERS = {
    "cycling_path": os.path.join(DATA, "cycling-path-network.geojson"),
    "park_connector": os.path.join(DATA, "park-connector-loop.geojson"),
    "nparks_track": os.path.join(DATA, "nparks-tracks.geojson"),
}

POPULATION_CSV = os.path.join(DATA, "population-census-2020.csv")
RESALE_CSV = os.path.join(DATA, "resale-flat-prices.csv")

OUT_GEOJSON = os.path.join(GEO, f"{DISTRICT}-subzone-summary.geojson")
OUT_CSV = os.path.join(GEO, f"{DISTRICT}-subzone-summary.csv")

# Resale town name mapping: district name → HDB resale CSV town name
# Most districts match directly; exceptions listed here
RESALE_TOWN_MAP = {
    "queenstown": "QUEENSTOWN",
    "bishan": "BISHAN",
    "outram": "CENTRAL AREA",
    "tampines": "TAMPINES",
    "newton": None,  # No HDB resale data for Newton
}

# Census name exceptions — subzone names that don't match simple title-casing
CENSUS_NAME_EXCEPTIONS = {
    "NATIONAL UNIVERSITY OF S'PORE": "National University Of S'pore",
    "PEARL'S HILL": "Pearl's Hill",
    "PEOPLE'S PARK": "People's Park",
    "MONK'S HILL": "Monk's Hill",
}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Approximate metres-per-degree at Singapore latitude (~1.3 N)
M_PER_DEG_LAT = 111_320
M_PER_DEG_LNG = 111_320 * math.cos(math.radians(1.3))


def load_geojson(path):
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def polygon_area_km2(geom):
    """Rough area in km² using the Shoelace formula in projected coords."""
    if geom.geom_type == "MultiPolygon":
        return sum(polygon_area_km2(p) for p in geom.geoms)
    coords = list(geom.exterior.coords)
    # Convert degrees to metres
    ref_lat = sum(c[1] for c in coords) / len(coords)
    m_lng = 111_320 * math.cos(math.radians(ref_lat))
    m_lat = 111_320
    projected = [(c[0] * m_lng, c[1] * m_lat) for c in coords]
    # Shoelace
    n = len(projected)
    area = 0
    for i in range(n):
        j = (i + 1) % n
        area += projected[i][0] * projected[j][1]
        area -= projected[j][0] * projected[i][1]
    return abs(area) / 2.0 / 1e6  # m² → km²


def line_length_m(line_geom):
    """Approximate length in metres for a LineString in WGS84."""
    coords = list(line_geom.coords)
    total = 0.0
    for i in range(len(coords) - 1):
        dx = (coords[i + 1][0] - coords[i][0]) * M_PER_DEG_LNG
        dy = (coords[i + 1][1] - coords[i][1]) * M_PER_DEG_LAT
        total += math.sqrt(dx * dx + dy * dy)
    return total


def subzone_to_census_name(subzone_name):
    """Convert UPPER CASE subzone name to census-style title case."""
    if subzone_name in CENSUS_NAME_EXCEPTIONS:
        return CENSUS_NAME_EXCEPTIONS[subzone_name]
    return subzone_name.title()


# ---------------------------------------------------------------------------
# 1. Load subzones
# ---------------------------------------------------------------------------
print(f"Loading {DISTRICT} subzones...")
subzones_gj = load_geojson(SUBZONES_PATH)
subzones = []
for feat in subzones_gj["features"]:
    geom = shape(feat["geometry"])
    props = feat["properties"]
    subzones.append({
        "name": props["SUBZONE_N"],
        "code": props["SUBZONE_C"],
        "geom": geom,
        "area_km2": polygon_area_km2(geom),
        "original_feature": feat,
    })
print(f"  {len(subzones)} subzones loaded")

# Build combined boundary for point filtering
boundary_geom = unary_union([sz["geom"] for sz in subzones])

# ---------------------------------------------------------------------------
# 2. Point-in-polygon counts
# ---------------------------------------------------------------------------
print("Counting points per subzone...")
# Initialise counts
for sz in subzones:
    for layer_name in POINT_LAYERS:
        sz[layer_name] = 0
    sz["mrt_stations"] = set()  # unique station names

for layer_name, path in POINT_LAYERS.items():
    if not os.path.exists(path):
        print(f"  WARNING: {path} not found, skipping")
        continue
    gj = load_geojson(path)
    count = 0
    for feat in gj["features"]:
        geom = feat["geometry"]
        if geom is None or geom["type"] != "Point":
            continue
        pt = Point(geom["coordinates"][:2])
        for sz in subzones:
            if sz["geom"].contains(pt):
                sz[layer_name] += 1
                if layer_name == "mrt_exits":
                    station = feat["properties"].get("STATION_NA", "")
                    if station:
                        sz["mrt_stations"].add(station)
                count += 1
                break
    print(f"  {layer_name}: {count} points in {DISTRICT}")

# Convert station sets to counts
for sz in subzones:
    sz["mrt_station_count"] = len(sz["mrt_stations"])
    del sz["mrt_stations"]

# ---------------------------------------------------------------------------
# 3. Line-in-polygon length (clip + sum)
# ---------------------------------------------------------------------------
print("Computing line lengths per subzone...")
for sz in subzones:
    sz["cycling_path_km"] = 0.0
    sz["park_connector_km"] = 0.0
    sz["nparks_track_km"] = 0.0

for layer_name, path in LINE_LAYERS.items():
    if not os.path.exists(path):
        print(f"  WARNING: {path} not found, skipping")
        continue
    gj = load_geojson(path)
    key = f"{layer_name}_km"
    total = 0.0
    for feat in gj["features"]:
        geom_raw = feat["geometry"]
        if geom_raw is None:
            continue
        line = shape(geom_raw)
        if not (line.geom_type == "LineString" or line.geom_type == "MultiLineString"):
            continue
        for sz in subzones:
            try:
                clipped = sz["geom"].intersection(line)
            except Exception:
                continue
            if clipped.is_empty:
                continue
            length_m = 0
            if clipped.geom_type == "LineString":
                length_m = line_length_m(clipped)
            elif clipped.geom_type == "MultiLineString":
                for part in clipped.geoms:
                    length_m += line_length_m(part)
            elif clipped.geom_type == "GeometryCollection":
                for part in clipped.geoms:
                    if part.geom_type == "LineString":
                        length_m += line_length_m(part)
                    elif part.geom_type == "MultiLineString":
                        for sub in part.geoms:
                            length_m += line_length_m(sub)
            sz[key] += length_m / 1000.0
            total += length_m / 1000.0
    print(f"  {layer_name}: {total:.2f} km total in {DISTRICT}")

# ---------------------------------------------------------------------------
# 3b. Polygon-area aggregation (clip + sum area)
# ---------------------------------------------------------------------------
print("Computing polygon areas per subzone...")

POLYGON_AREA_LAYERS = {
    "green_space": os.path.join(DATA, "nparks-parks-nature-reserves.geojson"),
    "abc_waters": os.path.join(DATA, "abc-waters.geojson"),
}

for sz in subzones:
    sz["green_space_area_km2"] = 0.0
    sz["abc_waters_area_km2"] = 0.0

for layer_name, path in POLYGON_AREA_LAYERS.items():
    if not os.path.exists(path):
        print(f"  WARNING: {path} not found, skipping")
        continue
    gj = load_geojson(path)
    key = f"{layer_name}_area_km2"
    total = 0.0
    for feat in gj["features"]:
        geom_raw = feat["geometry"]
        if geom_raw is None:
            continue
        geom = shape(geom_raw)
        if not geom.intersects(boundary_geom):
            continue
        for sz in subzones:
            try:
                clipped = sz["geom"].intersection(geom)
            except Exception:
                continue
            if clipped.is_empty:
                continue
            if clipped.geom_type in ("Polygon", "MultiPolygon"):
                area = polygon_area_km2(clipped)
                sz[key] += area
                total += area
    print(f"  {layer_name}: {total:.4f} km² total in {DISTRICT}")

# ---------------------------------------------------------------------------
# 4. Population (Census 2020)
# ---------------------------------------------------------------------------
print("Joining population data...")

for sz in subzones:
    sz["population_total"] = 0
    sz["population_male"] = 0
    sz["population_female"] = 0
    sz["population_elderly_65plus"] = 0

if os.path.exists(POPULATION_CSV):
    # Build name mapping: census row name → subzone name
    # Census uses title case e.g. "Commonwealth", subzone uses upper "COMMONWEALTH"
    census_rows = {}
    with open(POPULATION_CSV, encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            census_rows[row["Number"].strip()] = row

    def safe_int(val):
        val = val.strip().replace(",", "")
        if val in ("-", "", "na"):
            return 0
        return int(val)

    for sz in subzones:
        census_name = subzone_to_census_name(sz["name"])
        if census_name and census_name in census_rows:
            row = census_rows[census_name]
            sz["population_total"] = safe_int(row["Total_Total"])
            sz["population_male"] = safe_int(row["Males_Total"])
            sz["population_female"] = safe_int(row["Females_Total"])

            # Sum 65+ age brackets
            elderly = 0
            for age_col in ["65_69", "70_74", "75_79", "80_84", "85_89", "90andOver"]:
                key = f"Total_{age_col}"
                if key in row:
                    elderly += safe_int(row[key])
            sz["population_elderly_65plus"] = elderly
            print(f"  {sz['name']}: pop={sz['population_total']}, elderly={sz['population_elderly_65plus']}")
        else:
            print(f"  WARNING: census data not found for {sz['name']} (tried: {census_name})")
else:
    print(f"  WARNING: {POPULATION_CSV} not found, skipping population data")

# ---------------------------------------------------------------------------
# 5. Resale flat prices — join via buildings geometry
# ---------------------------------------------------------------------------
print("Joining resale flat prices...")

for sz in subzones:
    sz["resale_median_price"] = None
    sz["resale_transaction_count"] = 0

resale_town = RESALE_TOWN_MAP.get(DISTRICT, DISTRICT.upper())

if resale_town and os.path.exists(RESALE_CSV) and os.path.exists(BUILDINGS_PATH):
    # Street name abbreviation map (resale CSV → building GeoJSON)
    STREET_ABBREVS = {
        "AVE": "AVENUE", "CL": "CLOSE", "CRES": "CRESCENT", "CT": "COURT",
        "DR": "DRIVE", "HTS": "HEIGHTS", "LN": "LANE", "PL": "PLACE",
        "RD": "ROAD", "ST": "STREET", "TER": "TERRACE",
        "C'WEALTH": "COMMONWEALTH", "QUEEN'S": "QUEEN'S",
    }

    def normalise_street(name):
        """Expand common HDB resale CSV abbreviations to match OSM names."""
        parts = name.strip().upper().split()
        return " ".join(STREET_ABBREVS.get(p, p) for p in parts)

    # Load resale transactions for district, 2023-2025
    resale_by_block_street = {}
    with open(RESALE_CSV, encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row["town"] != resale_town:
                continue
            year = int(row["month"].split("-")[0])
            if year < 2023:
                continue
            key = (row["block"].strip(), normalise_street(row["street_name"]))
            resale_by_block_street.setdefault(key, []).append(float(row["resale_price"]))

    print(f"  {sum(len(v) for v in resale_by_block_street.values())} {DISTRICT} resale transactions (2023-2025)")

    # Load buildings to map block+street → subzone
    print("  Mapping resale prices to subzones via building geometry...")
    buildings_gj = load_geojson(BUILDINGS_PATH)

    block_to_subzone = {}
    for feat in buildings_gj["features"]:
        props = feat["properties"]
        block = (props.get("addr_housenumber") or "").strip()
        street = (props.get("addr_street") or "").strip().upper()
        if not block or not street:
            continue
        key = (block, street)
        if key in block_to_subzone:
            continue
        pt = shape(feat["geometry"]).centroid
        for sz in subzones:
            if sz["geom"].contains(pt):
                block_to_subzone[key] = sz["name"]
                break

    # Aggregate prices per subzone
    resale_prices_by_subzone = {}
    matched_txns = 0
    for (block, street), prices in resale_by_block_street.items():
        sz_name = block_to_subzone.get((block, street))
        if sz_name:
            resale_prices_by_subzone.setdefault(sz_name, []).extend(prices)
            matched_txns += len(prices)

    print(f"  {matched_txns} transactions matched to subzones")

    all_prices = [p for prices in resale_by_block_street.values() for p in prices]
    district_median = statistics.median(all_prices) if all_prices else 0

    for sz in subzones:
        prices = resale_prices_by_subzone.get(sz["name"], [])
        if prices:
            sz["resale_median_price"] = round(statistics.median(prices))
            sz["resale_transaction_count"] = len(prices)

    print(f"  District-wide median: ${district_median:,.0f}")
elif not resale_town:
    print(f"  No HDB resale town mapping for {DISTRICT}, skipping resale data")
elif not os.path.exists(BUILDINGS_PATH):
    print(f"  WARNING: {BUILDINGS_PATH} not found, skipping resale data")
else:
    print(f"  WARNING: {RESALE_CSV} not found, skipping resale data")

# ---------------------------------------------------------------------------
# 6. Building stats
# ---------------------------------------------------------------------------
print("Computing building stats per subzone...")

for sz in subzones:
    sz["building_count"] = 0
    sz["hdb_count"] = 0
    sz["building_heights"] = []

if os.path.exists(BUILDINGS_PATH):
    # Load buildings (tracked in git, always available)
    buildings_gj = load_geojson(BUILDINGS_PATH)

    for feat in buildings_gj["features"]:
        props = feat["properties"]
        centroid = shape(feat["geometry"]).centroid
        for sz in subzones:
            if sz["geom"].contains(centroid):
                sz["building_count"] += 1
                if props.get("hdb_match"):
                    sz["hdb_count"] += 1
                h = props.get("height_m")
                if h is not None and h > 0:
                    sz["building_heights"].append(h)
                break
else:
    print(f"  WARNING: {BUILDINGS_PATH} not found, skipping building stats")

for sz in subzones:
    heights = sz["building_heights"]
    sz["mean_height_m"] = round(statistics.mean(heights), 1) if heights else None
    sz["max_height_m"] = round(max(heights), 1) if heights else None
    del sz["building_heights"]

# HDB age and dwelling stats per subzone
for sz in subzones:
    sz["hdb_years"] = []
    sz["total_dwelling_units"] = 0

if os.path.exists(BUILDINGS_PATH):
    for feat in buildings_gj["features"]:
        props = feat["properties"]
        if not props.get("hdb_match"):
            continue
        centroid = shape(feat["geometry"]).centroid
        for sz in subzones:
            if sz["geom"].contains(centroid):
                year = props.get("hdb_year_completed")
                if year is not None:
                    try:
                        sz["hdb_years"].append(int(year))
                    except (ValueError, TypeError):
                        pass
                units = props.get("hdb_total_dwelling_units")
                if units is not None:
                    try:
                        sz["total_dwelling_units"] += int(units)
                    except (ValueError, TypeError):
                        pass
                break

for sz in subzones:
    years = sz["hdb_years"]
    sz["avg_hdb_year"] = round(statistics.mean(years)) if years else None
    del sz["hdb_years"]

# ---------------------------------------------------------------------------
# 7. Derived density metrics
# ---------------------------------------------------------------------------
print("Computing density metrics...")
for sz in subzones:
    area = sz["area_km2"]
    sz["population_density"] = round(sz["population_total"] / area) if area > 0 else 0

    amenity_total = (
        sz["hawker_centres"]
        + sz["parks"]
        + sz["supermarkets"]
        + sz["gyms"]
        + sz["community_clubs"]
    )
    sz["amenity_count"] = amenity_total
    sz["amenity_density"] = round(amenity_total / area, 1) if area > 0 else 0
    sz["dwelling_density"] = round(sz["total_dwelling_units"] / area) if area > 0 else 0
    sz["green_space_pct"] = round((sz["green_space_area_km2"] / area) * 100, 1) if area > 0 else 0

# ---------------------------------------------------------------------------
# 8. Remote sensing metrics (from fetch_global_layers.py GEE output)
# ---------------------------------------------------------------------------
print("Joining remote sensing data...")

RS_FIELDS = [
    "lst_mean_c", "ndvi_mean", "ndbi_mean", "ghsl_height_mean",
    "dsm_mean", "dem_mean", "canopy_cover_pct", "canopy_loss_pct",
]

for sz in subzones:
    for f in RS_FIELDS:
        sz[f] = None

if os.path.exists(REMOTE_SENSING_PATH):
    rs_gj = load_geojson(REMOTE_SENSING_PATH)
    rs_by_name = {}
    for feat in rs_gj["features"]:
        name = feat["properties"].get("subzone_name")
        if name:
            rs_by_name[name] = feat["properties"]
    for sz in subzones:
        rs = rs_by_name.get(sz["name"], {})
        for f in RS_FIELDS:
            val = rs.get(f)
            if val is not None:
                sz[f] = val
    print(f"  Merged remote sensing for {len(rs_by_name)} subzones")
else:
    print(f"  WARNING: {REMOTE_SENSING_PATH} not found, skipping")

# ---------------------------------------------------------------------------
# 8b. Mapillary GVI — mean Green View Index per subzone
# ---------------------------------------------------------------------------
print("Joining Mapillary GVI data...")

for sz in subzones:
    sz["gvi_mean"] = None

if os.path.exists(MAPILLARY_GVI_PATH):
    gvi_gj = load_geojson(MAPILLARY_GVI_PATH)
    gvi_by_subzone = {}
    for feat in gvi_gj["features"]:
        gvi = feat["properties"].get("gvi")
        if gvi is None:
            continue
        geom = feat["geometry"]
        if geom is None or geom["type"] != "Point":
            continue
        pt = Point(geom["coordinates"][:2])
        for sz in subzones:
            if sz["geom"].contains(pt):
                gvi_by_subzone.setdefault(sz["name"], []).append(gvi)
                break
    for sz in subzones:
        vals = gvi_by_subzone.get(sz["name"], [])
        if vals:
            sz["gvi_mean"] = round(statistics.mean(vals), 4)
    matched = sum(1 for sz in subzones if sz["gvi_mean"] is not None)
    print(f"  GVI computed for {matched} subzones")
else:
    print(f"  WARNING: {MAPILLARY_GVI_PATH} not found, skipping GVI")

# ---------------------------------------------------------------------------
# 8c. Walkability scoring (from generate_walkability.py output)
# ---------------------------------------------------------------------------
print("Joining walkability data...")

WALK_FIELDS = [
    "intersection_density", "transit_access_score",
    "destination_accessibility", "walkability_index",
]

for sz in subzones:
    for f in WALK_FIELDS:
        sz[f] = None

if os.path.exists(WALKABILITY_PATH):
    with open(WALKABILITY_PATH, encoding="utf-8") as f:
        walk_data = json.load(f)
    walk_by_name = {entry["subzone_name"]: entry for entry in walk_data}
    for sz in subzones:
        entry = walk_by_name.get(sz["name"], {})
        for f in WALK_FIELDS:
            val = entry.get(f)
            if val is not None:
                sz[f] = val
    matched = sum(1 for sz in subzones if sz["walkability_index"] is not None)
    print(f"  Merged walkability for {matched} subzones")
else:
    print(f"  WARNING: {WALKABILITY_PATH} not found, skipping walkability")

# ---------------------------------------------------------------------------
# 9. Assemble output
# ---------------------------------------------------------------------------
print("Writing output files...")

# Define output field order
FIELDS = [
    "subzone_name", "subzone_code", "area_km2",
    "population_total", "population_male", "population_female",
    "population_elderly_65plus", "population_density",
    "hawker_centres", "parks", "supermarkets", "gyms", "community_clubs",
    "amenity_count", "amenity_density",
    "mrt_exits", "mrt_station_count",
    "cycling_path_km", "park_connector_km", "nparks_track_km",
    "sport_facilities",
    "green_space_area_km2", "green_space_pct", "abc_waters_area_km2",
    "building_count", "hdb_count", "mean_height_m", "max_height_m",
    "avg_hdb_year", "total_dwelling_units", "dwelling_density",
    "resale_median_price", "resale_transaction_count",
    "lst_mean_c", "ndvi_mean", "ndbi_mean", "ghsl_height_mean",
    "dsm_mean", "dem_mean", "canopy_cover_pct", "canopy_loss_pct",
    "gvi_mean",
    "intersection_density", "transit_access_score",
    "destination_accessibility", "walkability_index",
]

features_out = []
csv_rows = []

for sz in subzones:
    props = {
        "subzone_name": sz["name"],
        "subzone_code": sz["code"],
        "area_km2": round(sz["area_km2"], 3),
        "population_total": sz["population_total"],
        "population_male": sz["population_male"],
        "population_female": sz["population_female"],
        "population_elderly_65plus": sz["population_elderly_65plus"],
        "population_density": sz["population_density"],
        "hawker_centres": sz["hawker_centres"],
        "parks": sz["parks"],
        "supermarkets": sz["supermarkets"],
        "gyms": sz["gyms"],
        "community_clubs": sz["community_clubs"],
        "amenity_count": sz["amenity_count"],
        "amenity_density": sz["amenity_density"],
        "mrt_exits": sz["mrt_exits"],
        "mrt_station_count": sz["mrt_station_count"],
        "cycling_path_km": round(sz["cycling_path_km"], 2),
        "park_connector_km": round(sz["park_connector_km"], 2),
        "nparks_track_km": round(sz["nparks_track_km"], 2),
        "sport_facilities": sz["sport_facilities"],
        "green_space_area_km2": round(sz["green_space_area_km2"], 4),
        "green_space_pct": sz["green_space_pct"],
        "abc_waters_area_km2": round(sz["abc_waters_area_km2"], 4),
        "building_count": sz["building_count"],
        "hdb_count": sz["hdb_count"],
        "mean_height_m": sz["mean_height_m"],
        "max_height_m": sz["max_height_m"],
        "avg_hdb_year": sz["avg_hdb_year"],
        "total_dwelling_units": sz["total_dwelling_units"],
        "dwelling_density": sz["dwelling_density"],
        "resale_median_price": sz["resale_median_price"],
        "resale_transaction_count": sz["resale_transaction_count"],
        "lst_mean_c": sz["lst_mean_c"],
        "ndvi_mean": sz["ndvi_mean"],
        "ndbi_mean": sz["ndbi_mean"],
        "ghsl_height_mean": sz["ghsl_height_mean"],
        "dsm_mean": sz["dsm_mean"],
        "dem_mean": sz["dem_mean"],
        "canopy_cover_pct": sz["canopy_cover_pct"],
        "canopy_loss_pct": sz["canopy_loss_pct"],
        "gvi_mean": sz["gvi_mean"],
        "intersection_density": sz["intersection_density"],
        "transit_access_score": sz["transit_access_score"],
        "destination_accessibility": sz["destination_accessibility"],
        "walkability_index": sz["walkability_index"],
    }

    # GeoJSON feature: reuse original geometry
    features_out.append({
        "type": "Feature",
        "properties": props,
        "geometry": sz["original_feature"]["geometry"],
    })

    csv_rows.append(props)

# Sort by subzone code
features_out.sort(key=lambda f: f["properties"]["subzone_code"])
csv_rows.sort(key=lambda r: r["subzone_code"])

# Write GeoJSON
geojson_out = {
    "type": "FeatureCollection",
    "features": features_out,
}
with open(OUT_GEOJSON, "w", encoding="utf-8") as f:
    json.dump(geojson_out, f, indent=2, ensure_ascii=False)
print(f"  -> {OUT_GEOJSON}")

# Write CSV
with open(OUT_CSV, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=FIELDS)
    writer.writeheader()
    writer.writerows(csv_rows)
print(f"  -> {OUT_CSV}")

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
print(f"\n=== Summary ({DISTRICT}) ===")
total_pop = sum(sz["population_total"] for sz in subzones)
total_buildings = sum(sz["building_count"] for sz in subzones)
total_hdb = sum(sz["hdb_count"] for sz in subzones)
total_amenities = sum(sz["amenity_count"] for sz in subzones)
print(f"  Subzones: {len(subzones)}")
print(f"  Total population: {total_pop:,}")
print(f"  Total buildings: {total_buildings:,}")
print(f"  Total HDB blocks: {total_hdb}")
print(f"  Total amenity points: {total_amenities}")
print("Done!")
