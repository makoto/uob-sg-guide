#!/usr/bin/env python3
"""
Extract boundary and subzones for a planning area from URA master plan data.

Usage:
    python3 extract_district.py BISHAN
    python3 extract_district.py QUEENSTOWN

Input:  data/ura/planning-area-boundary.geojson
        data/ura/subzone-boundary.geojson
Output: docs/geo/{district}-boundary.geojson
        docs/geo/{district}-subzones.geojson
"""

import json
import os
import sys

BASE = os.path.dirname(os.path.abspath(__file__))
URA = os.path.join(BASE, "data", "ura")
GEO = os.path.join(BASE, "docs", "geo")

PLANNING_AREA_PATH = os.path.join(URA, "planning-area-boundary.geojson")
SUBZONE_PATH = os.path.join(URA, "subzone-boundary.geojson")


def load_geojson(path):
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def write_geojson(path, features):
    fc = {"type": "FeatureCollection", "features": features}
    with open(path, "w", encoding="utf-8") as f:
        json.dump(fc, f, ensure_ascii=False)
    size_kb = os.path.getsize(path) / 1024
    print(f"  -> {os.path.basename(path)}: {len(features)} features ({size_kb:.1f} KB)")


def extract_district(district_name):
    district_upper = district_name.upper()
    district_lower = district_name.lower().replace(" ", "-")

    print(f"Extracting {district_upper}...")

    # 1. Extract planning area boundary
    pa_gj = load_geojson(PLANNING_AREA_PATH)
    boundary_features = [
        f for f in pa_gj["features"]
        if f["properties"]["PLN_AREA_N"] == district_upper
    ]
    if not boundary_features:
        print(f"  ERROR: Planning area '{district_upper}' not found")
        print("  Available:", sorted(set(
            f["properties"]["PLN_AREA_N"] for f in pa_gj["features"]
        )))
        sys.exit(1)

    boundary_path = os.path.join(GEO, f"{district_lower}-boundary.geojson")
    write_geojson(boundary_path, boundary_features)

    # 2. Extract subzones
    sz_gj = load_geojson(SUBZONE_PATH)
    subzone_features = [
        f for f in sz_gj["features"]
        if f["properties"]["PLN_AREA_N"] == district_upper
    ]

    subzones_path = os.path.join(GEO, f"{district_lower}-subzones.geojson")
    write_geojson(subzones_path, subzone_features)

    # 3. Print bounding box
    all_coords = []
    for feat in boundary_features:
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
    bbox = [min(lons), min(lats), max(lons), max(lats)]
    print(f"  Bounding box: [{bbox[0]:.4f}, {bbox[1]:.4f}, {bbox[2]:.4f}, {bbox[3]:.4f}]")
    print(f"  Subzone names: {[f['properties']['SUBZONE_N'] for f in subzone_features]}")
    print("Done!")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 extract_district.py DISTRICT_NAME")
        print("Example: python3 extract_district.py BISHAN")
        sys.exit(1)
    extract_district(sys.argv[1])
