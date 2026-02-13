#!/usr/bin/env python3
"""
Filter island-wide GeoJSON datasets to a district boundary.

Reads source files from data/data-gov-sg/, clips/filters to the district
boundary polygon, and writes compact GeoJSON files to docs/geo/gov-sg/.

Usage:
    python3 filter_district_layers.py                  # default: queenstown
    python3 filter_district_layers.py --district bishan
"""

import argparse
import json
import os

from shapely.geometry import shape, mapping, Point

BASE = os.path.dirname(os.path.abspath(__file__))
GEO = os.path.join(BASE, "docs", "geo")
GOV_SG = os.path.join(GEO, "gov-sg")
DATA = os.path.join(BASE, "data", "data-gov-sg")

# (source_file, output_suffix) — district prefix is prepended at runtime
POINT_LAYERS = [
    ("hawker-centres.geojson", "hawker-centres.geojson"),
    ("parks.geojson", "parks.geojson"),
    ("mrt-station-exits.geojson", "mrt-exits.geojson"),
    ("supermarkets.geojson", "supermarkets.geojson"),
    ("gyms.geojson", "gyms.geojson"),
    ("community-clubs.geojson", "community-clubs.geojson"),
    ("pre-schools.geojson", "preschools.geojson"),
    ("chas-clinics.geojson", "chas-clinics.geojson"),
    ("park-facilities.geojson", "park-facilities.geojson"),
]

LINE_LAYERS = [
    ("cycling-path-network.geojson", "cycling-paths.geojson"),
    ("park-connector-loop.geojson", "park-connectors.geojson"),
]

POLYGON_LAYERS = [
    ("ura-height-control.geojson", "ura-height-control.geojson"),
]


def load_geojson(path):
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def write_geojson(path, features):
    fc = {"type": "FeatureCollection", "features": features}
    with open(path, "w", encoding="utf-8") as f:
        json.dump(fc, f, ensure_ascii=False)
    size_kb = os.path.getsize(path) / 1024
    print(f"  -> {os.path.basename(path)}: {len(features)} features ({size_kb:.0f} KB)")


def filter_points(source_path, boundary_geom):
    gj = load_geojson(source_path)
    results = []
    for feat in gj["features"]:
        geom = feat.get("geometry")
        if geom is None or geom["type"] != "Point":
            continue
        pt = Point(geom["coordinates"][:2])
        if boundary_geom.contains(pt):
            results.append(feat)
    return results


def filter_and_clip_lines(source_path, boundary_geom):
    gj = load_geojson(source_path)
    results = []
    for feat in gj["features"]:
        geom_raw = feat.get("geometry")
        if geom_raw is None:
            continue
        line = shape(geom_raw)
        if line.geom_type not in ("LineString", "MultiLineString"):
            continue
        if not line.intersects(boundary_geom):
            continue
        try:
            clipped = boundary_geom.intersection(line)
        except Exception:
            continue
        if clipped.is_empty:
            continue
        # Normalize geometry collections to just line parts
        if clipped.geom_type == "GeometryCollection":
            from shapely.geometry import MultiLineString as MLS
            line_parts = [
                g for g in clipped.geoms
                if g.geom_type in ("LineString", "MultiLineString")
            ]
            if not line_parts:
                continue
            clipped = MLS(line_parts) if len(line_parts) > 1 else line_parts[0]
        results.append({
            "type": "Feature",
            "properties": feat["properties"],
            "geometry": mapping(clipped),
        })
    return results


def filter_polygons(source_path, boundary_geom):
    gj = load_geojson(source_path)
    results = []
    for feat in gj["features"]:
        geom_raw = feat.get("geometry")
        if geom_raw is None:
            continue
        geom = shape(geom_raw)
        if not geom.intersects(boundary_geom):
            continue
        try:
            clipped = boundary_geom.intersection(geom)
        except Exception:
            continue
        if clipped.is_empty:
            continue
        results.append({
            "type": "Feature",
            "properties": feat["properties"],
            "geometry": mapping(clipped),
        })
    return results


def main():
    parser = argparse.ArgumentParser(description="Filter island-wide layers to a district boundary")
    parser.add_argument("--district", default="queenstown",
                        help="District name (default: queenstown)")
    args = parser.parse_args()

    district = args.district.lower().replace(" ", "-")
    boundary_path = os.path.join(GEO, f"{district}-boundary.geojson")

    print(f"Loading {district} boundary...")
    boundary_gj = load_geojson(boundary_path)
    boundary_geom = shape(boundary_gj["features"][0]["geometry"])
    print(f"  Boundary loaded ({boundary_geom.geom_type})")

    os.makedirs(GOV_SG, exist_ok=True)

    total_size = 0

    print("\nFiltering point layers...")
    for source_name, output_suffix in POINT_LAYERS:
        source_path = os.path.join(DATA, source_name)
        if not os.path.exists(source_path):
            print(f"  WARNING: {source_name} not found, skipping")
            continue
        features = filter_points(source_path, boundary_geom)
        out_path = os.path.join(GOV_SG, f"{district}-{output_suffix}")
        write_geojson(out_path, features)
        total_size += os.path.getsize(out_path)

    print("\nFiltering & clipping line layers...")
    for source_name, output_suffix in LINE_LAYERS:
        source_path = os.path.join(DATA, source_name)
        if not os.path.exists(source_path):
            print(f"  WARNING: {source_name} not found, skipping")
            continue
        features = filter_and_clip_lines(source_path, boundary_geom)
        out_path = os.path.join(GOV_SG, f"{district}-{output_suffix}")
        write_geojson(out_path, features)
        total_size += os.path.getsize(out_path)

    print("\nFiltering polygon layers...")
    for source_name, output_suffix in POLYGON_LAYERS:
        source_path = os.path.join(DATA, source_name)
        if not os.path.exists(source_path):
            print(f"  WARNING: {source_name} not found, skipping")
            continue
        features = filter_polygons(source_path, boundary_geom)
        out_path = os.path.join(GOV_SG, f"{district}-{output_suffix}")
        write_geojson(out_path, features)
        total_size += os.path.getsize(out_path)

    print(f"\nTotal output: {total_size / 1024:.0f} KB")
    print("Done!")


if __name__ == "__main__":
    main()
