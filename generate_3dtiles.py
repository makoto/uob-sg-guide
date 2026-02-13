"""
Generate 3D Tiles (B3DM) from buildings GeoJSON for a district.

Produces:
  docs/3dtiles/{district}/buildings.b3dm  — Batched 3D Model with extruded buildings
  docs/3dtiles/{district}/tileset.json    — 3D Tiles manifest

Usage:
    python3 generate_3dtiles.py                    # default: queenstown
    python3 generate_3dtiles.py --district bishan

Coordinate pipeline:
  WGS84 (lon,lat,height) → ECEF (via pyproj, accurate)
  → subtract ECEF center → local ECEF offsets
  → identity rotation tileset transform adds center back
  No glTF node transform — vertices are raw ECEF offsets.
"""

import argparse
import json
import math
import numpy as np
import mapbox_earcut
from pathlib import Path
from pyproj import Transformer
from py3dtiles.tileset.content import B3dm
from py3dtiles.tileset.content.batch_table import BatchTable
from py3dtiles.tileset.content.b3dm_feature_table import B3dmFeatureTable

# ── Args ───────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description="Generate 3D Tiles from buildings GeoJSON")
parser.add_argument("--district", default="queenstown",
                    help="District name (default: queenstown)")
args = parser.parse_args()

DISTRICT = args.district.lower().replace(" ", "-")

# ── Config ──────────────────────────────────────────────────────────────────
GEOJSON_PATH = Path(f"docs/geo/global/{DISTRICT}-buildings.geojson")
OUT_DIR = Path(f"docs/3dtiles/{DISTRICT}")

# WGS84 → ECEF transformer
wgs84_to_ecef = Transformer.from_crs("EPSG:4326", "EPSG:4978", always_xy=True)

# ── Helpers ─────────────────────────────────────────────────────────────────

def lonlat_to_ecef(lon, lat, height=0.0):
    """Convert WGS84 (lon,lat,height) to ECEF (x,y,z) using pyproj."""
    x, y, z = wgs84_to_ecef.transform(lon, lat, height)
    return [x, y, z]


def triangulate_polygon(exterior, holes=None):
    """Triangulate a 2D polygon using mapbox_earcut.
    exterior: list of [x,y] pairs (closed or open)
    holes: list of lists of [x,y] pairs
    Returns: triangle index array reshaped to (N,3)
    """
    if len(exterior) > 1 and exterior[0] == exterior[-1]:
        exterior = exterior[:-1]

    coords = list(exterior)
    ring_lengths = [len(exterior)]

    if holes:
        for hole in holes:
            if len(hole) > 1 and hole[0] == hole[-1]:
                hole = hole[:-1]
            ring_lengths.append(len(hole))
            coords.extend(hole)

    flat = np.array(coords, dtype=np.float64)
    rings = np.array(ring_lengths, dtype=np.uint32)
    indices = mapbox_earcut.triangulate_float64(flat, rings)
    return indices.reshape(-1, 3)


# ── Main ────────────────────────────────────────────────────────────────────

def main():
    if not GEOJSON_PATH.exists():
        print(f"WARNING: {GEOJSON_PATH} not found, skipping 3D tiles for {DISTRICT}")
        return

    print(f"Reading GeoJSON for {DISTRICT}...")
    with open(GEOJSON_PATH) as f:
        geojson = json.load(f)

    features = geojson["features"]
    print(f"  {len(features)} features")

    # ── Pass 1: compute ECEF center from all building centroids ─────────
    print("Computing ECEF center...")
    cx_sum, cy_sum, cz_sum = 0.0, 0.0, 0.0
    count = 0
    for feature in features:
        geom = feature["geometry"]
        if geom["type"] == "Polygon":
            ring = geom["coordinates"][0]
        elif geom["type"] == "MultiPolygon":
            ring = geom["coordinates"][0][0]
        else:
            continue
        # Centroid of first ring
        lons = [p[0] for p in ring]
        lats = [p[1] for p in ring]
        clon = sum(lons) / len(lons)
        clat = sum(lats) / len(lats)
        ex, ey, ez = lonlat_to_ecef(clon, clat, 0)
        cx_sum += ex
        cy_sum += ey
        cz_sum += ez
        count += 1

    ecef_center = np.array([cx_sum / count, cy_sum / count, cz_sum / count])
    print(f"  ECEF center: ({ecef_center[0]:.3f}, {ecef_center[1]:.3f}, {ecef_center[2]:.3f})")

    # ── Pass 2: build geometry ──────────────────────────────────────────
    print("Building geometry...")
    all_positions = []
    all_normals = []
    all_triangles = []
    batch_ids = []

    bt_props = {
        "name": [],
        "building": [],
        "height_m": [],
        "height_source": [],
        "data_source": [],
        "building_levels": [],
        "hdb_year_completed": [],
    }

    vertex_offset = 0
    batch_id = 0
    skipped = 0

    for feature in features:
        props = feature["properties"]
        geom = feature["geometry"]
        height = float(props.get("height_m", 10) or 10)

        if geom["type"] == "Polygon":
            polygons = [geom["coordinates"]]
        elif geom["type"] == "MultiPolygon":
            polygons = geom["coordinates"]
        else:
            skipped += 1
            continue

        feature_verts = []
        feature_normals = []
        feature_tris = []

        for poly_coords in polygons:
            exterior = poly_coords[0]
            holes = poly_coords[1:] if len(poly_coords) > 1 else None

            # Convert exterior ring to ECEF at base (h=0) and top (h=height)
            ext_base_ecef = [lonlat_to_ecef(p[0], p[1], 0) for p in exterior]
            ext_top_ecef = [lonlat_to_ecef(p[0], p[1], height) for p in exterior]

            # Remove closing vertex
            if len(ext_base_ecef) > 1 and ext_base_ecef[0] == ext_base_ecef[-1]:
                ext_base_ecef = ext_base_ecef[:-1]
                ext_top_ecef = ext_top_ecef[:-1]

            if len(ext_base_ecef) < 3:
                continue

            # Subtract ECEF center → local offsets
            ext_base_local = [(np.array(p) - ecef_center).tolist() for p in ext_base_ecef]
            ext_top_local = [(np.array(p) - ecef_center).tolist() for p in ext_top_ecef]

            # Prepare holes
            holes_base_local = None
            holes_top_local = None
            if holes:
                holes_base_local = []
                holes_top_local = []
                for hole in holes:
                    hb = [lonlat_to_ecef(p[0], p[1], 0) for p in hole]
                    ht = [lonlat_to_ecef(p[0], p[1], height) for p in hole]
                    if len(hb) > 1 and hb[0] == hb[-1]:
                        hb = hb[:-1]
                        ht = ht[:-1]
                    if len(hb) >= 3:
                        holes_base_local.append([(np.array(p) - ecef_center).tolist() for p in hb])
                        holes_top_local.append([(np.array(p) - ecef_center).tolist() for p in ht])

            # Triangulate footprint in 2D (lon, lat) — topology only
            ext_2d = [[p[0], p[1]] for p in exterior]
            if len(ext_2d) > 1 and ext_2d[0] == ext_2d[-1]:
                ext_2d = ext_2d[:-1]
            holes_2d = None
            if holes:
                holes_2d = []
                for hole in holes:
                    h2d = [[p[0], p[1]] for p in hole]
                    if len(h2d) > 1 and h2d[0] == h2d[-1]:
                        h2d = h2d[:-1]
                    if len(h2d) >= 3:
                        holes_2d.append(h2d)

            try:
                foot_tris = triangulate_polygon(ext_2d, holes_2d)
            except Exception:
                continue

            if len(foot_tris) == 0:
                continue

            # Build all ring vertices (for indexing consistency with earcut)
            all_base = list(ext_base_local)
            all_top = list(ext_top_local)
            if holes_base_local:
                for hb, ht in zip(holes_base_local, holes_top_local):
                    all_base.extend(hb)
                    all_top.extend(ht)

            n_ring = len(all_base)
            local_offset = len(feature_verts)

            # ─ Top face ─
            for p in all_top:
                feature_verts.append(p)
            # Compute up normal at center (radially outward from Earth)
            up = ecef_center / np.linalg.norm(ecef_center)
            up_list = up.tolist()
            for _ in range(n_ring):
                feature_normals.append(up_list)

            for tri in foot_tris:
                feature_tris.append([
                    local_offset + tri[0],
                    local_offset + tri[1],
                    local_offset + tri[2],
                ])

            # ─ Side walls ─
            rings_for_walls = [ext_base_local]
            rings_top_for_walls = [ext_top_local]
            if holes_base_local:
                rings_for_walls.extend(holes_base_local)
                rings_top_for_walls.extend(holes_top_local)

            for ring_b, ring_t in zip(rings_for_walls, rings_top_for_walls):
                n = len(ring_b)
                for i in range(n):
                    j = (i + 1) % n

                    base_off = len(feature_verts)

                    b0 = ring_b[i]
                    b1 = ring_b[j]
                    t0 = ring_t[i]
                    t1 = ring_t[j]

                    feature_verts.extend([b0, b1, t1, t0])

                    # Wall normal (cross product of two edges)
                    v0, v1, v2 = np.array(b0), np.array(b1), np.array(t1)
                    e1 = v1 - v0
                    e2 = v2 - v0
                    wn = np.cross(e1, e2)
                    wn_len = np.linalg.norm(wn)
                    if wn_len > 1e-10:
                        wn = wn / wn_len
                    wall_normal = wn.tolist()
                    feature_normals.extend([wall_normal] * 4)

                    feature_tris.append([base_off, base_off + 1, base_off + 2])
                    feature_tris.append([base_off, base_off + 2, base_off + 3])

        if len(feature_verts) == 0:
            skipped += 1
            continue

        n_verts = len(feature_verts)
        all_positions.extend(feature_verts)
        all_normals.extend(feature_normals)
        batch_ids.extend([batch_id] * n_verts)
        for tri in feature_tris:
            all_triangles.append([
                tri[0] + vertex_offset,
                tri[1] + vertex_offset,
                tri[2] + vertex_offset,
            ])
        vertex_offset += n_verts

        for key in bt_props:
            val = props.get(key, "") or ""
            if key == "height_m":
                bt_props[key].append(float(val) if val else 0.0)
            else:
                bt_props[key].append(str(val))
        batch_id += 1

    print(f"  {batch_id} buildings, {skipped} skipped")
    print(f"  {len(all_positions)} vertices, {len(all_triangles)} triangles")

    # ── Create B3DM ─────────────────────────────────────────────────────
    positions = np.array(all_positions, dtype=np.float32)
    normals = np.array(all_normals, dtype=np.float32)
    triangles = np.array(all_triangles, dtype=np.uint32).flatten()
    batchids = np.array(batch_ids, dtype=np.uint32)

    bt = BatchTable()
    for key, values in bt_props.items():
        bt.add_property_as_json(key, values)

    # Feature table with BATCH_LENGTH (required for CesiumJS to read batch table)
    ft = B3dmFeatureTable()
    ft.set_batch_length(batch_id)

    print(f"Creating B3DM (BATCH_LENGTH={batch_id}, no node transform)...")
    b3dm = B3dm.from_numpy_arrays(
        positions,
        triangles,
        batch_table=bt,
        feature_table=ft,
        normal=normals,
        batchids=batchids,
        transform=None,  # No glTF node transform
    )

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    b3dm_path = OUT_DIR / "buildings.b3dm"
    b3dm.save_as(b3dm_path)
    print(f"  Saved {b3dm_path} ({b3dm_path.stat().st_size / 1024 / 1024:.1f} MB)")

    # ── Create tileset.json ─────────────────────────────────────────────

    # Bounding region
    all_lons, all_lats = [], []
    for f in features:
        g = f["geometry"]
        if g["type"] == "Polygon":
            ring = g["coordinates"][0]
        elif g["type"] == "MultiPolygon":
            ring = g["coordinates"][0][0]
        else:
            continue
        for p in ring:
            all_lons.append(p[0])
            all_lats.append(p[1])

    west = math.radians(min(all_lons))
    south = math.radians(min(all_lats))
    east = math.radians(max(all_lons))
    north = math.radians(max(all_lats))

    # Identity rotation + ECEF center translation (column-major)
    transform = [
        1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        ecef_center[0], ecef_center[1], ecef_center[2], 1,
    ]

    tileset = {
        "asset": {
            "version": "1.0",
            "generator": f"{DISTRICT}-liveability-study",
            "gltfUpAxis": "Z",
        },
        "geometricError": 500,
        "root": {
            "boundingVolume": {
                "region": [west, south, east, north, 0, 200]
            },
            "geometricError": 100,
            "refine": "ADD",
            "content": {
                "uri": "buildings.b3dm"
            },
            "transform": transform,
        }
    }

    tileset_path = OUT_DIR / "tileset.json"
    with open(tileset_path, "w") as f:
        json.dump(tileset, f, indent=2)
    print(f"  Saved {tileset_path}")

    print("Done!")


if __name__ == "__main__":
    main()
