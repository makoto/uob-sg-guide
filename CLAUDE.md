# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Rules

- **Never commit unless explicitly asked.** Do not run `git add` or `git commit` unless the user specifically requests it.

## Project Overview

This is an academic research project investigating **perceived liveability in Singapore's Queenstown district**, with multi-district support for comparison. The goal is to connect subjective perceptions of liveability with objective spatial indicators and translate findings into planning insights.

## Study Focus

- **Study area:** Queenstown (primary), with multi-district support for Bishan, Outram, Tampines, Newton
- **Spatial boundary:** Defined from the URA Master Plan via `src/extract_district.py`
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

## Deliverables

1. A **geospatial map** showing distribution of perceived liveability and selected urban factors
2. A **report with analysis** (~1000 words)

## Repository Structure

### Data Pipeline

Island-wide GeoJSON datasets live in `data/data-gov-sg/` (gitignored, ~11.5 MB). Six Python scripts in the `src/` directory process them into committed outputs under `docs/`. All scripts that accept `--district` default to `queenstown`.

| Script | Input | Output | Run in CI |
|---|---|---|---|
| `src/extract_district.py` | URA master plan subzones GeoJSON | `docs/geo/{district}-boundary.geojson` + `{district}-subzones.geojson` | Yes |
| `src/generate_subzone_summary.py` | Subzones + point/line/CSV + remote sensing + walkability | `docs/geo/{district}-subzone-summary.geojson` + `.csv` | Yes |
| `src/filter_district_layers.py` | Boundary + 16 island-wide GeoJSONs | 16 `docs/geo/gov-sg/{district}-*.geojson` files | Yes |
| `src/generate_3dtiles.py` | Buildings GeoJSON | `docs/3dtiles/{district}/tileset.json` + B3DM | Yes |
| `src/fetch_global_layers.py` | OSMnx, GEE, Mapillary APIs | `docs/geo/global/{district}-*.geojson` files | Manual (requires conda env + API keys). Default tasks: `buildings osmnx gee mapillary`. Additional tasks: `raster_pngs` (GEE raster PNG overlays), `rs_grid` (100m RS grid cells) |
| `src/generate_walkability.py` | OSMnx walk network + subzones + GEE SRTM elevation | `docs/geo/global/{district}-walkability.json` + `{district}-walkability-grid.geojson` | Manual (requires conda env + GEE auth) |

The CI workflow (`.github/workflows/deploy.yml`) loops over all 5 districts for the first four scripts (calling `python src/scriptname.py`). `src/fetch_global_layers.py` and `src/generate_walkability.py` are run manually with the `zensvi` conda environment (osmnx, networkx, shapely). `src/fetch_global_layers.py` additionally requires GEE authentication + a Mapillary API key. `src/fetch_global_layers.py` has 4 default tasks: `buildings` (OSM footprints + HDB enrichment), `osmnx` (street network), `gee` (remote sensing), `mapillary` (street-level imagery). Two additional tasks: `raster_pngs` (GEE raster PNG overlays) and `rs_grid` (100m RS grid cells). Buildings, street network, remote sensing (subzone-level + RS grid + raster PNGs), and walkability data now exist for all 5 districts. Mapillary data exists for Queenstown, Bishan, and Outram (not yet Tampines/Newton). Note: the Mapillary layer has been removed from the viewer but the data files and script task remain.

### Private Directory

The `private/` directory (gitignored) contains literature review files and site visit reports that are not published to GitHub Pages.

### Processed Layers in `docs/geo/`

Files are organised into subdirectories by data provenance. All filenames use `{district}-` prefix (e.g., `queenstown-`, `bishan-`).

- **`docs/geo/`** (root) — per-district boundary, subzones, subzone-summary GeoJSON + CSV (5 districts)
- **`docs/geo/gov-sg/`** — 16 layers per district from Singapore government open data (80 files total for 5 districts)
- **`docs/geo/global/`** — OSM buildings, street network, remote sensing (subzone-level + 100m RS grid), and walkability (scores + grid) for all 5 districts. Mapillary data for Queenstown, Bishan, Outram (not shown in viewer).
- **`docs/geo/global/rasters/`** — 4 satellite raster PNGs per district from GEE (NDVI, LST, NDBI, canopy; 1024px each) for all 5 districts
- **`docs/3dtiles/{district}/`** — per-district 3D tileset directories for all 5 districts (queenstown: 15 MB, bishan: 4.8 MB, outram: 1.4 MB, tampines: 9.8 MB, newton: 1.3 MB)

`gov-sg/` point layers per district (filtered by boundary containment):
- `{district}-hawker-centres.geojson`, `{district}-parks.geojson`, `{district}-mrt-exits.geojson`, `{district}-supermarkets.geojson`, `{district}-gyms.geojson`, `{district}-community-clubs.geojson`, `{district}-preschools.geojson`, `{district}-chas-clinics.geojson`, `{district}-park-facilities.geojson`, `{district}-sport-facilities.geojson`

`gov-sg/` line layers (clipped to boundary):
- `{district}-cycling-paths.geojson`, `{district}-park-connectors.geojson`, `{district}-nparks-tracks.geojson`

`gov-sg/` polygon layers:
- `{district}-ura-height-control.geojson`, `{district}-nparks-parks-nature-reserves.geojson`, `{district}-abc-waters.geojson`

`global/` building layers (all 5 districts):
- `{district}-buildings.geojson` — queenstown (8,671), bishan (3,766), outram (1,589), tampines (2,910), newton (759)
- `{district}-buildings-osm-only.geojson` — OSM-only buildings (non-HDB)
- `{district}-buildings-hdb-enriched.geojson` — HDB-enriched buildings

`global/` street network (all 5 districts):
- `{district}-street-network.geojson` — edges with betweenness centrality from OSMnx drive network: queenstown (3,788), bishan (1,363), outram (711), tampines (4,144), newton (838)

`global/` remote sensing — subzone-level (all 5 districts):
- `{district}-remote-sensing.geojson` — subzone polygons with LST, NDVI, NDBI, GHSL height, DSM, DEM, canopy cover from GEE: queenstown (15), bishan (3), outram (4), tampines (5), newton (6)

`global/` remote sensing — grid + rasters (all 5 districts):
- `{district}-rs-grid.geojson` — clickable 100m grid cells with per-cell NDVI, NDBI, LST, canopy %, DEM elevation, slope from GEE: queenstown (2,394), bishan (840), outram (175), tampines (2,230), newton (265). 6 properties: ndvi, ndbi, lst_c, canopy_pct, dem_m (SRTM 30m), slope_deg (terrain slope).
- `rasters/{district}-ndvi.png`, `{district}-lst.png`, `{district}-ndbi.png`, `{district}-canopy.png` (1024px georeferenced satellite imagery from GEE, all 5 districts)

`global/` Mapillary (Queenstown, Bishan, Outram — data files exist but layer removed from viewer):
- `{district}-mapillary-gvi.geojson` — sampled street-level image point locations: queenstown (5,000), bishan (5,000), outram (5,000). Not yet available for Tampines or Newton.

`global/` walkability (all 5 districts):
- `{district}-walkability.json` — subzone-level BEH-NWI walkability scores (intersection_density, transit_access_score, destination_accessibility, walkability_index, transit_access_slope, dest_access_slope, walkability_slope): queenstown (15), bishan (3), outram (4), tampines (5), newton (6). Slope-penalised variants use Tobler's hiking function with SRTM elevation.
- `{district}-walkability-grid.geojson` — clickable 100m grid cells with per-cell pop_density, transit_access, dest_access, int_density, walkability, transit_slope, dest_slope, slope_walkability: queenstown (2,178), bishan (769), outram (135), tampines (2,085), newton (210). Walk network edges include elevation_start, gradient, tobler_time.

### Tech Catalogue (`docs/tech-catalogue.html`)

Documents all technologies, libraries, and frameworks used in the project: Python geospatial stack (Shapely, GeoPandas, OSMnx, NetworkX, GEE, ZenSVI, py3dtiles, pyproj, mapbox_earcut, NumPy), frontend (CesiumJS, deck.gl, MapLibre GL JS, CARTO basemaps, OpenTopoMap, Bing Maps satellite via Cesium Ion), APIs (GEE, Mapillary, Overpass, Cesium ion), data formats (GeoJSON, CSV, 3D Tiles, PNG), CI/CD (GitHub Actions, GitHub Pages), script reference, and environment setup.

### Glossary & Data Dictionary (`docs/glossary.html`)

Definitions of acronyms and technical terms used across the project, organised in 5 sections: Remote Sensing & Terrain (DEM, DSM, SRTM, NDVI, NDBI, LST, canopy, GHSL, GEE, slope, RS grid, NIR/SWIR), Walkability & Network Analysis (BEH-NWI, intersection density, transit access, destination accessibility, Tobler's hiking function, slope-adjusted walkability, betweenness centrality, Dijkstra, walk network, UNA, gravity accessibility), Geospatial & Data Formats (GeoJSON, 3D Tiles/B3DM, CesiumJS, OSM/OSMnx, ECEF, WGS84, EPSG codes, POI, GVI), Singapore-Specific Terms (HDB, URA, MRT, CHAS, NParks, ABC Waters, LTA, subzone, planning area, hawker centre), and Statistics & Visualisation (z-score, choropleth, colour ramp, CNT/RATE/IDX/SAT badges, feature count).

### User Manual (`docs/user-manual.html`)

End-user guide for the deck.gl map viewer. Covers map navigation (pan, zoom, rotate, tilt), the district selector (5 districts), basemap styles (3 CARTO options), the layer panel structure (10 collapsible groups, provenance filter, feature counts, lazy loading), a full reference table of all 33 toggleable layers across 10 groups, click-popup content by layer type, the 39 choropleth metrics across 6 sub-categories (with type badge definitions: CNT, RATE, IDX, SAT), the 5 RS grid colour modes (NDVI, DEM, Slope, LST, NDBI), the 5 UNA gravity accessibility layers (Queenstown only), legend descriptions, and a tips/troubleshooting table. Includes 9 annotated screenshots in `docs/images/` (viewer overview, info panel, layer panel, expanded group, choropleth panel, choropleth map, RS grid modes, click popup, district switch). Styled to match `glossary.html` and `tech-catalogue.html`.

### 3D Viewer (`docs/3dtiles/viewer.html`)

CesiumJS-based viewer with multi-district support:
- **District selector** dropdown (Queenstown, Bishan, Outram, Tampines, Newton) with `switchDistrict()` function that reloads all layers, camera position, and raster bounds per district. `layerFile()` helper dynamically swaps district name in file paths. Missing data handled gracefully (console warnings, no errors).
- 3D building tileset per district (`docs/3dtiles/{district}/`), toggleable via "3D Buildings" checkbox in the District group
- Boundary + subzone overlays (on by default)
- **Layer panel** (top-right, max-width 280px): 27 toggleable overlays grouped by 10 collapsible categories (District, Food & Daily Needs, Transit, Green & Recreation, Active Mobility, Community, Housing, Planning, Street Network, Remote Sensing). The District group (formerly Boundaries) contains 3 items: District boundary, Subzones, and 3D Buildings (all on by default). Each group header is clickable with a chevron indicator (▶/▼) and smooth CSS max-height animation. District group is expanded by default; all others are collapsed.
- **Feature counts**: Each layer row shows its feature count in parentheses, e.g. "Hawker centres (9) SG". Group headers show total feature count for the group. Zero-count layers are visually dimmed. Counts are stored in `layerCounts` per district within the `DISTRICTS` config object.
- **Provenance filter bar**: two toggle buttons (Gov.sg / Global) that filter layer checkboxes by data source. Each layer row shows a source badge (SG or globe icon).
- Layers are lazy-loaded on first checkbox toggle via `GeoJsonDataSource`
- **Data catalogue link**: Footer link at the bottom of the layer panel ("Data Catalogue & Downloads") linking to `../data-catalogue.html` (opens in a new tab).
- **HDB blocks** (Housing group): 308 dots coloured by construction era (red pre-1980, orange 1980-1999, green 2000+). Click shows year completed, dwelling units, and use flags.
- **NParks parks & nature reserves** (Green & Recreation group): park/reserve boundary polygons from NParks. Available for all 5 districts (queenstown: 12; bishan: 20; outram: 5; tampines: 13; newton: 1).
- **ABC Waters** (Green & Recreation group): Active Beautiful Clean Waters programme area polygons. Available for all 5 districts (queenstown: 5; bishan: 3; outram: 0; tampines: 3; newton: 0).
- **NParks tracks** (Green & Recreation group): hiking/walking track lines inside parks, clipped to boundary. Available for all 5 districts (queenstown: 1,559; bishan: 481; outram: 192; tampines: 449; newton: 0).
- **Sport facilities** (Community group): sport facility point locations. Available for all 5 districts (queenstown: 1; bishan: 2; outram: 1; tampines: 1; newton: 0).
- **Height control zones** (Planning group): semi-transparent polygons with storey-limit labels from `HT_CTL_TXT`.
- **Street network** (Street Network group): road edges coloured by betweenness centrality (OSMnx drive network) with dynamic legend. Available for all 5 districts (queenstown: 3,788; bishan: 1,363; outram: 711; tampines: 4,144; newton: 838). Plus an Intersection density grid layer using 100m cells: coloured by intersection density computed via spatial kernel method (counts walk-network nodes with degree>=3 within a 300m radius, divided by circle area). Uses a blue colour ramp with dynamic legend; clicking a grid cell shows the intersection density value and method description. Available for all 5 districts.
- **Remote sensing** (Remote Sensing group): 4 raster image overlays (NDVI, LST, NDBI, canopy) via `SingleTileImageryProvider`, plus a clickable 100m RS grid with a **colour mode selector** offering 5 toggleable modes (NDVI, DEM, Slope, LST, NDBI). Small toggle buttons appear below the RS grid checkbox in the layer panel. Clicking a mode re-colours all grid cells using a 6-step dynamic colour ramp with min/max values computed from the data. The legend updates to show a gradient bar with the current field name and value range. Default mode is NDVI. RS grid cells include 6 properties: NDVI, NDBI, LST, canopy %, DEM elevation (SRTM 30m), and slope. Both raster PNGs and RS grid are available for all 5 districts. Subzone-level remote sensing GeoJSON also available for all 5 districts for choropleth use.
- Basemap toggle cycling through 4 options: Dark (CARTO dark_all, default), Light (CARTO light_all), Satellite (Bing Maps via Cesium Ion asset 2, lazy-loaded), Terrain (OpenTopoMap)
- **Choropleth heatmap**: collapsible "Choropleth" accordion section containing 5 sub-categories: Demographics & Housing, Amenities & Infrastructure, Buildings, Remote Sensing, and Walkability. Each metric row shows a coloured type badge (CNT = count, RATE = rate/density, IDX = index/score, SAT = satellite-derived). Click a metric to activate it; an "Off" button at the top deactivates the choropleth. Internally, `CHOROPLETH_CATEGORIES` array (metrics derived via `flatMap`). 33 metrics total across the 5 sub-categories (YlOrRd 5-step ramp, alpha 0.55). Original 15 metrics: population density, elderly share, amenity density, MRT stations, cycling paths, park connectors, resale flat price, avg building height, HDB blocks, total buildings, max building height, resale transactions, avg HDB year built, dwelling units, dwelling density. Plus 5 green/blue infrastructure metrics in Amenities & Infrastructure: Green space coverage (%), Green space area (km2), ABC Waters area (km2), NParks tracks (km), Sport facilities. Plus 6 remote sensing metrics: Vegetation (NDVI), Built-up (NDBI), Land Surface Temp, Tree canopy cover, GHSL building height, Surface elevation. Plus 7 walkability metrics: Intersection density, Transit access, Destination access, Walkability (BEH-NWI), Walkability (slope-adjusted), Transit access (slope), Destination access (slope). The 3 slope-adjusted metrics use Tobler's hiking function with SRTM elevation to penalise steep gradients. Data sourced from `{district}-subzone-summary.geojson` (per-district), lazy-loaded and cached. Legend updates with formatted min/max per metric.
- **Disabled Cesium widgets**: `navigationHelpButton`, `sceneModePicker`, and `fullscreenButton` are disabled to reduce UI clutter.

### deck.gl Viewer (`docs/viewer-deckgl.html`)

Lightweight alternative viewer using deck.gl + MapLibre GL JS. Located at `docs/viewer-deckgl.html` (not inside `3dtiles/`). The existing CesiumJS viewer is not modified; both viewers coexist.

- **Stack**: deck.gl v9 CDN bundle (~700 KB) + MapLibre GL JS v4 for basemaps. No build tools, no Cesium Ion token required.
- **Buildings**: Extruded natively from GeoJSON via `GeoJsonLayer` with `extruded: true` (no 3D Tiles/B3DM needed). Height-based colour coding matching the CesiumJS viewer.
- **District selector** dropdown (Queenstown, Bishan, Outram, Tampines, Newton) with `switchDistrict()` function and `FlyToInterpolator` camera animation.
- **33 toggleable layers** grouped by 10 collapsible categories (same as CesiumJS viewer): District, Food & Daily Needs, Transit, Green & Recreation, Active Mobility, Community, Housing, Planning, Street Network, Remote Sensing. Additionally includes a Walkability grid layer and 5 UNA gravity accessibility layers (when data available).
- **Layer panel** (top-right): same UI structure as CesiumJS viewer -- collapsible groups with chevron indicators, feature counts per layer, zero-count dimming, provenance filter bar (Gov.sg / Global), source badges (SG / globe icon), RS grid colour mode selector (NDVI, DEM, Slope, LST, NDBI), and Data Catalogue footer link.
- **39 choropleth metrics** in 6 sub-categories: Demographics & Housing (8), Amenities & Infrastructure (9), Buildings (3), Remote Sensing (6), Walkability (7), UNA Gravity (6). Same YlOrRd 5-step ramp with type badges (CNT/RATE/IDX/SAT). The 6 UNA Gravity metrics are additional compared to the CesiumJS viewer's 33.
- **3 CARTO basemaps**: Dark (dark-matter, default), Light (positron), Voyager. No satellite option (no free token-less provider available).
- **Data paths**: Relative to `docs/` -- e.g., `geo/queenstown-boundary.geojson`, `geo/global/queenstown-buildings.geojson`. The `layerFile()` helper swaps district names in paths.
- **Raster overlays**: 4 satellite PNGs (NDVI, LST, NDBI, canopy) rendered via `BitmapLayer` with georeferenced bounds per district.
- **Click interaction**: `onMapClick` handler with `buildPopupHtml()` for contextual popups (buildings, HDB, RS grid, walkability grid, street network, height control, choropleth, gravity layers).
- **Lazy loading**: GeoJSON layers fetched on first checkbox toggle with URL-keyed cache (`dataCache`). Default layers (boundary, subzones, buildings) loaded on init.
