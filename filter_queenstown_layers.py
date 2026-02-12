#!/usr/bin/env python3
"""
Filter island-wide GeoJSON datasets to Queenstown boundary.

Reads source files from data/data-gov-sg/, clips/filters to the Queenstown
boundary polygon, and writes compact GeoJSON files to docs/geo/.
"""

import json
import os

from shapely.geometry import shape, mapping, Point

BASE = os.path.dirname(os.path.abspath(__file__))
GEO = os.path.join(BASE, "docs", "geo")
DATA = os.path.join(BASE, "data", "data-gov-sg")
BOUNDARY_PATH = os.path.join(GEO, "queenstown-boundary.geojson")

POINT_LAYERS = [
    ("hawker-centres.geojson", "queenstown-hawker-centres.geojson"),
    ("parks.geojson", "queenstown-parks.geojson"),
    ("mrt-station-exits.geojson", "queenstown-mrt-exits.geojson"),
    ("supermarkets.geojson", "queenstown-supermarkets.geojson"),
    ("gyms.geojson", "queenstown-gyms.geojson"),
    ("community-clubs.geojson", "queenstown-community-clubs.geojson"),
    ("pre-schools.geojson", "queenstown-preschools.geojson"),
    ("chas-clinics.geojson", "queenstown-chas-clinics.geojson"),
    ("park-facilities.geojson", "queenstown-park-facilities.geojson"),
]

LINE_LAYERS = [
    ("cycling-path-network.geojson", "queenstown-cycling-paths.geojson"),
    ("park-connector-loop.geojson", "queenstown-park-connectors.geojson"),
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


def main():
    print("Loading Queenstown boundary...")
    boundary_gj = load_geojson(BOUNDARY_PATH)
    boundary_geom = shape(boundary_gj["features"][0]["geometry"])
    print(f"  Boundary loaded ({boundary_geom.geom_type})")

    total_size = 0

    print("\nFiltering point layers...")
    for source_name, output_name in POINT_LAYERS:
        source_path = os.path.join(DATA, source_name)
        if not os.path.exists(source_path):
            print(f"  WARNING: {source_name} not found, skipping")
            continue
        features = filter_points(source_path, boundary_geom)
        out_path = os.path.join(GEO, output_name)
        write_geojson(out_path, features)
        total_size += os.path.getsize(out_path)

    print("\nFiltering & clipping line layers...")
    for source_name, output_name in LINE_LAYERS:
        source_path = os.path.join(DATA, source_name)
        if not os.path.exists(source_path):
            print(f"  WARNING: {source_name} not found, skipping")
            continue
        features = filter_and_clip_lines(source_path, boundary_geom)
        out_path = os.path.join(GEO, output_name)
        write_geojson(out_path, features)
        total_size += os.path.getsize(out_path)

    print(f"\nTotal output: {total_size / 1024:.0f} KB")
    print("Done!")


if __name__ == "__main__":
    main()
