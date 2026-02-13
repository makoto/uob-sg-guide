# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Rules

- **Never commit unless explicitly asked.** Do not run `git add` or `git commit` unless the user specifically requests it.

## Project Overview

This is an academic research project investigating **perceived liveability in Singapore's Queenstown district**, with multi-district support for comparison. The goal is to connect subjective perceptions of liveability with objective spatial indicators and translate findings into planning insights.

## Study Focus

- **Study area:** Queenstown (primary), with multi-district support for Bishan, Outram, Tampines, Newton
- **Spatial boundary:** Defined from the URA Master Plan via `extract_district.py`
- **Supported districts:** Queenstown (15 subzones), Bishan (3), Outram (4), Tampines (5), Newton (6)

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

Island-wide GeoJSON datasets live in `data/data-gov-sg/` (gitignored, ~11.5 MB). Six Python scripts process them into committed outputs under `docs/`. All scripts that accept `--district` default to `queenstown`.

| Script | Input | Output | Run in CI |
|---|---|---|---|
| `extract_district.py` | URA master plan subzones GeoJSON | `docs/geo/{district}-boundary.geojson` + `{district}-subzones.geojson` | Yes |
| `generate_subzone_summary.py` | Subzones + point/line/CSV + remote sensing + walkability | `docs/geo/{district}-subzone-summary.geojson` + `.csv` | Yes |
| `filter_district_layers.py` | Boundary + 11 island-wide GeoJSONs | 11 `docs/geo/gov-sg/{district}-*.geojson` files | Yes |
| `generate_3dtiles.py` | Buildings GeoJSON | `docs/3dtiles/{district}/tileset.json` + B3DM | Yes |
| `fetch_global_layers.py` | OSMnx, GEE, Mapillary APIs | `docs/geo/global/{district}-*.geojson` files | Manual (requires conda env + API keys). Default tasks: `buildings osmnx gee mapillary`. Additional tasks: `raster_pngs` (GEE raster PNG overlays), `rs_grid` (100m RS grid cells) |
| `generate_walkability.py` | OSMnx walk network + subzones | `docs/geo/global/{district}-walkability.json` + `{district}-walkability-grid.geojson` | Manual (requires conda env) |

The CI workflow (`.github/workflows/deploy.yml`) loops over all 5 districts for the first four scripts. `fetch_global_layers.py` and `generate_walkability.py` are run manually with the `zensvi` conda environment (osmnx, networkx, shapely). `fetch_global_layers.py` additionally requires GEE authentication + a Mapillary API key. `fetch_global_layers.py` has 4 default tasks: `buildings` (OSM footprints + HDB enrichment), `osmnx` (street network), `gee` (remote sensing), `mapillary` (street-level imagery). Two additional tasks: `raster_pngs` (GEE raster PNG overlays) and `rs_grid` (100m RS grid cells). Buildings, street network, remote sensing (subzone-level + RS grid + raster PNGs), and walkability data now exist for all 5 districts. Mapillary exists for Queenstown, Bishan, and Outram (not yet Tampines/Newton).

### Processed Layers in `docs/geo/`

Files are organised into subdirectories by data provenance. All filenames use `{district}-` prefix (e.g., `queenstown-`, `bishan-`).

- **`docs/geo/`** (root) — per-district boundary, subzones, subzone-summary GeoJSON + CSV (5 districts)
- **`docs/geo/gov-sg/`** — 11 layers per district from Singapore government open data (55 files total for 5 districts)
- **`docs/geo/global/`** — OSM buildings, street network, remote sensing (subzone-level + 100m RS grid), and walkability (scores + grid) for all 5 districts. Mapillary for Queenstown, Bishan, Outram.
- **`docs/geo/global/rasters/`** — 4 satellite raster PNGs per district from GEE (NDVI, LST, NDBI, canopy; 1024px each) for all 5 districts
- **`docs/3dtiles/{district}/`** — per-district 3D tileset directories for all 5 districts (queenstown: 15 MB, bishan: 4.8 MB, outram: 1.4 MB, tampines: 9.8 MB, newton: 1.3 MB)

`gov-sg/` point layers per district (filtered by boundary containment):
- `{district}-hawker-centres.geojson`, `{district}-parks.geojson`, `{district}-mrt-exits.geojson`, `{district}-supermarkets.geojson`, `{district}-gyms.geojson`, `{district}-community-clubs.geojson`, `{district}-preschools.geojson`, `{district}-chas-clinics.geojson`, `{district}-park-facilities.geojson`

`gov-sg/` line layers (clipped to boundary):
- `{district}-cycling-paths.geojson`, `{district}-park-connectors.geojson`

`gov-sg/` polygon layers:
- `{district}-ura-height-control.geojson`

`global/` building layers (all 5 districts):
- `{district}-buildings.geojson` — queenstown (8,671), bishan (3,766), outram (1,589), tampines (2,910), newton (759)
- `{district}-buildings-osm-only.geojson` — OSM-only buildings (non-HDB)
- `{district}-buildings-hdb-enriched.geojson` — HDB-enriched buildings

`global/` street network (all 5 districts):
- `{district}-street-network.geojson` — edges with betweenness centrality from OSMnx drive network: queenstown (3,788), bishan (1,363), outram (711), tampines (4,144), newton (838)

`global/` remote sensing — subzone-level (all 5 districts):
- `{district}-remote-sensing.geojson` — subzone polygons with LST, NDVI, NDBI, GHSL height, DSM, DEM, canopy cover from GEE: queenstown (15), bishan (3), outram (4), tampines (5), newton (6)

`global/` remote sensing — grid + rasters (all 5 districts):
- `{district}-rs-grid.geojson` — clickable 100m grid cells with per-cell NDVI, NDBI, LST, canopy % from GEE: queenstown (2,394), bishan (840), outram (175), tampines (2,230), newton (265)
- `rasters/{district}-ndvi.png`, `{district}-lst.png`, `{district}-ndbi.png`, `{district}-canopy.png` (1024px georeferenced satellite imagery from GEE, all 5 districts)

`global/` Mapillary (Queenstown, Bishan, Outram):
- `{district}-mapillary-gvi.geojson` — sampled street-level image point locations: queenstown (5,000), bishan (5,000), outram (5,000). Not yet available for Tampines or Newton.

`global/` walkability (all 5 districts):
- `{district}-walkability.json` — subzone-level BEH-NWI walkability scores (intersection_density, transit_access_score, destination_accessibility, walkability_index): queenstown (15), bishan (3), outram (4), tampines (5), newton (6)
- `{district}-walkability-grid.geojson` — clickable 100m grid cells with per-cell pop_density, transit_access, dest_access, int_density, walkability: queenstown (2,178), bishan (769), outram (135), tampines (2,085), newton (210)

### Tech Catalogue (`docs/tech-catalogue.html`)

Documents all technologies, libraries, and frameworks used in the project: Python geospatial stack (Shapely, GeoPandas, OSMnx, NetworkX, GEE, ZenSVI, py3dtiles, pyproj, mapbox_earcut, NumPy), frontend (CesiumJS, CARTO basemaps), APIs (GEE, Mapillary, Overpass, Cesium ion), data formats (GeoJSON, CSV, 3D Tiles, PNG), CI/CD (GitHub Actions, GitHub Pages), script reference, and environment setup.

### 3D Viewer (`docs/3dtiles/viewer.html`)

CesiumJS-based viewer with multi-district support:
- **District selector** dropdown (Queenstown, Bishan, Outram, Tampines, Newton) with `switchDistrict()` function that reloads all layers, camera position, and raster bounds per district. `layerFile()` helper dynamically swaps district name in file paths. Missing data handled gracefully (console warnings, no errors).
- 3D building tileset per district (`docs/3dtiles/{district}/`)
- Boundary + subzone overlays (on by default)
- **Layer panel** (top-right): 26 toggleable overlays grouped by 12 categories (Boundaries, Food & Daily Needs, Transit, Green & Recreation, Active Mobility, Community, Housing, Planning, Street Network, Remote Sensing, Street-Level)
- **Provenance filter bar**: two toggle buttons (Gov.sg / Global) that filter layer checkboxes by data source. Each layer row shows a source badge (SG or globe icon).
- Layers are lazy-loaded on first checkbox toggle via `GeoJsonDataSource`
- **HDB blocks** (Housing group): 308 dots coloured by construction era (red pre-1980, orange 1980-1999, green 2000+). Click shows year completed, dwelling units, and use flags.
- **Height control zones** (Planning group): semi-transparent polygons with storey-limit labels from `HT_CTL_TXT`.
- **Street network** (Street Network group): road edges coloured by betweenness centrality (OSMnx drive network) with dynamic legend. Available for all 5 districts (queenstown: 3,788; bishan: 1,363; outram: 711; tampines: 4,144; newton: 838). Plus 4 walkability grid layers using 100m cells (like the remote sensing grid): Walkability grid (green), Intersection density grid (blue), Transit access grid (purple), Destination access grid (orange). Each uses a normalized colour ramp with dynamic legend; clicking a grid cell shows all 4 walkability metrics plus population density. Walkability grids available for all 5 districts.
- **Remote sensing** (Remote Sensing group): 4 raster image overlays (NDVI, LST, NDBI, canopy) via `SingleTileImageryProvider`, plus a clickable 100m RS grid coloured by NDVI with 6-step green ramp and dynamic legend. Both raster PNGs and RS grid are available for all 5 districts. Subzone-level remote sensing GeoJSON also available for all 5 districts for choropleth use.
- **Mapillary images** (Street-Level group): sampled street-level image point locations. Available for Queenstown, Bishan, Outram (5,000 points each). Not yet available for Tampines or Newton. Click shows a link to view the image on Mapillary.
- Basemap toggle (dark/light CARTO)
- **Choropleth heatmap**: dropdown (26 options incl. Off) to colour subzones by one of 25 metrics (YlOrRd 5-step ramp, alpha 0.55). Original 15 metrics: population density, elderly share, amenity density, MRT stations, cycling paths, park connectors, resale flat price, avg building height, HDB blocks, total buildings, max building height, resale transactions, avg HDB year built, dwelling units, dwelling density. Plus 6 remote sensing metrics: Vegetation (NDVI), Built-up (NDBI), Land Surface Temp, Tree canopy cover, GHSL building height, Surface elevation. Plus 4 walkability metrics: Intersection density, Transit access, Destination access, Walkability (BEH-NWI). Data sourced from `{district}-subzone-summary.geojson` (per-district), lazy-loaded and cached. Legend updates with formatted min/max per metric.
