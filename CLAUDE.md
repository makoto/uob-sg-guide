# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Rules

- **Never commit unless explicitly asked.** Do not run `git add` or `git commit` unless the user specifically requests it.

## Project Overview

This is an academic research project investigating **perceived liveability in Singapore's Queenstown district**. The goal is to connect subjective perceptions of liveability with objective spatial indicators and translate findings into planning insights.

## Study Focus

- **Study area:** Queenstown (assigned to this group; other groups cover Bishan, Outram, Tampines, Newton)
- **Spatial boundary:** Must be defined from the URA Master Plan

## Research Framework

Two dimensions to investigate:
1. **Subjective:** How people perceive and feel about liveability in Queenstown
2. **Objective (pick ≥2):** Social intensity (POIs), accessibility (street/transit/green-blue networks), active mobility (walking/cycling), microclimate (weather), public facilities (amenities/services)

## Key Data Sources

- Open data: https://data.gov.sg/
- URA Map platform: https://eservice.ura.gov.sg/maps/
- OneMap (including 3D): https://www.onemap.gov.sg/
- Master Plan: https://www.uradraftmasterplan.gov.sg/regional-plans/
- Transport data: https://datamall.lta.gov.sg/content/datamall/en.html
- Street view imagery: https://www.mapillary.com/app

## Deliverables

1. A **geospatial map** showing distribution of perceived liveability and selected urban factors
2. A **report with analysis** (~1000 words)

## Repository Structure

### Data Pipeline

Island-wide GeoJSON datasets live in `data/data-gov-sg/` (gitignored, ~11.5 MB). Five Python scripts process them into committed outputs under `docs/`:

| Script | Input | Output | Run in CI |
|---|---|---|---|
| `generate_subzone_summary.py` | Subzones + point/line/CSV + remote sensing + walkability | `docs/geo/queenstown-subzone-summary.geojson` + `.csv` | Yes |
| `filter_queenstown_layers.py` | Boundary + 11 island-wide GeoJSONs | 11 `docs/geo/gov-sg/queenstown-*.geojson` files (~373 KB) | Yes |
| `generate_3dtiles.py` | Buildings GeoJSON | `docs/3dtiles/tileset.json` + B3DM | Yes |
| `fetch_global_layers.py` | OSMnx, GEE, Mapillary APIs | 3 `docs/geo/global/queenstown-*.geojson` files | Manual (requires conda env + API keys) |
| `generate_walkability.py` | OSMnx walk network + subzones | `docs/geo/global/queenstown-walkability.json` (~3 KB) + `queenstown-walkability-grid.geojson` (~704 KB, 2,178 cells) | Manual (requires conda env) |

The first three run in `.github/workflows/deploy.yml` on push to `main`. `fetch_global_layers.py` and `generate_walkability.py` are run manually with the `zensvi` conda environment (osmnx, networkx, shapely). `fetch_global_layers.py` additionally requires GEE authentication + a Mapillary API key.

### Processed Layers in `docs/geo/`

Files are organised into subdirectories by data provenance:

- **`docs/geo/`** (root) — boundary, subzones, subzone-summary GeoJSON + CSV
- **`docs/geo/gov-sg/`** — 12 layers from Singapore government open data (data.gov.sg / URA)
- **`docs/geo/global/`** — 9 layers from globally available sources (OSM buildings, street network, remote sensing, Mapillary, walkability grid)
- **`docs/geo/global/rasters/`** — 4 satellite raster PNGs from GEE (NDVI, LST, NDBI, canopy; 1024px each)

`gov-sg/` point layers (filtered by boundary containment):
- `queenstown-hawker-centres.geojson` (9), `queenstown-parks.geojson` (9), `queenstown-mrt-exits.geojson` (26), `queenstown-supermarkets.geojson` (19), `queenstown-gyms.geojson` (4), `queenstown-community-clubs.geojson` (4), `queenstown-preschools.geojson` (81), `queenstown-chas-clinics.geojson` (43), `queenstown-park-facilities.geojson` (370)

`gov-sg/` line layers (clipped to boundary):
- `queenstown-cycling-paths.geojson` (110), `queenstown-park-connectors.geojson` (36)

`gov-sg/` polygon layers:
- `queenstown-ura-height-control.geojson` (7)

`global/` layers:
- `queenstown-buildings.geojson` (8,671), `queenstown-buildings-osm-only.geojson` (8,363), `queenstown-buildings-hdb-enriched.geojson` (308)
- `queenstown-street-network.geojson` (3,788 edges with betweenness centrality from OSMnx drive network)
- `queenstown-remote-sensing.geojson` (15 subzone polygons with LST, NDVI, NDBI, GHSL height, DSM, DEM, canopy cover from GEE)
- `queenstown-rs-grid.geojson` (2,402 clickable 100m grid cells with per-cell NDVI, NDBI, LST, canopy % from GEE)
- `rasters/queenstown-ndvi.png`, `queenstown-lst.png`, `queenstown-ndbi.png`, `queenstown-canopy.png` (1024px georeferenced satellite imagery from GEE)
- `queenstown-mapillary-gvi.geojson` (5,000 sampled Mapillary image point locations)
- `queenstown-walkability.json` (15 subzones with BEH-NWI walkability scores: intersection_density, transit_access_score, destination_accessibility, walkability_index)
- `queenstown-walkability-grid.geojson` (2,178 clickable 100m grid cells with per-cell pop_density, transit_access, dest_access, int_density, walkability)

### Tech Catalogue (`docs/tech-catalogue.html`)

Documents all technologies, libraries, and frameworks used in the project: Python geospatial stack (Shapely, GeoPandas, OSMnx, NetworkX, GEE, ZenSVI, py3dtiles, pyproj, mapbox_earcut, NumPy), frontend (CesiumJS, CARTO basemaps), APIs (GEE, Mapillary, Overpass, Cesium ion), data formats (GeoJSON, CSV, 3D Tiles, PNG), CI/CD (GitHub Actions, GitHub Pages), script reference, and environment setup.

### 3D Viewer (`docs/3dtiles/viewer.html`)

CesiumJS-based viewer with:
- 3D building tileset (height-colored)
- Boundary + subzone overlays (on by default)
- **Layer panel** (top-right): 26 toggleable overlays grouped by 12 categories (Boundaries, Food & Daily Needs, Transit, Green & Recreation, Active Mobility, Community, Housing, Planning, Street Network, Remote Sensing, Street-Level)
- **Provenance filter bar**: two toggle buttons (Gov.sg / Global) that filter layer checkboxes by data source. Each layer row shows a source badge (SG or globe icon).
- Layers are lazy-loaded on first checkbox toggle via `GeoJsonDataSource`
- **HDB blocks** (Housing group): 308 dots coloured by construction era (red pre-1980, orange 1980-1999, green 2000+). Click shows year completed, dwelling units, and use flags.
- **Height control zones** (Planning group): semi-transparent polygons with storey-limit labels from `HT_CTL_TXT`.
- **Street network** (Street Network group): 3,788 road edges coloured by betweenness centrality (OSMnx drive network) with dynamic legend. Plus 4 walkability grid layers using 100m cells (like the remote sensing grid): Walkability grid (green), Intersection density grid (blue), Transit access grid (purple), Destination access grid (orange). Each uses a normalized colour ramp with dynamic legend; clicking a grid cell shows all 4 walkability metrics plus population density.
- **Remote sensing** (Remote Sensing group): 4 raster image overlays (NDVI, LST, NDBI, canopy) via `SingleTileImageryProvider`, plus a clickable 100m grid (2,402 cells) coloured by NDVI with 6-step green ramp and dynamic legend. The old subzone-level remote sensing GeoJSON remains for choropleth use.
- **Mapillary images** (Street-Level group): 5,000 sampled street-level image point locations. Click shows a link to view the image on Mapillary.
- Basemap toggle (dark/light CARTO)
- **Choropleth heatmap**: dropdown (26 options incl. Off) to colour subzones by one of 25 metrics (YlOrRd 5-step ramp, alpha 0.55). Original 15 metrics: population density, elderly share, amenity density, MRT stations, cycling paths, park connectors, resale flat price, avg building height, HDB blocks, total buildings, max building height, resale transactions, avg HDB year built, dwelling units, dwelling density. Plus 6 remote sensing metrics: Vegetation (NDVI), Built-up (NDBI), Land Surface Temp, Tree canopy cover, GHSL building height, Surface elevation. Plus 4 walkability metrics: Intersection density, Transit access, Destination access, Walkability (BEH-NWI). Data sourced from `queenstown-subzone-summary.geojson`, lazy-loaded and cached. Legend updates with formatted min/max per metric.
