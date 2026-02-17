# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Rules

- **Never commit unless explicitly asked.** Do not run `git add` or `git commit` unless the user specifically requests it.

## Project Overview

This is an academic research project investigating **perceived liveability in Singapore's Queenstown district**, with multi-district support for comparison. The goal is to connect subjective perceptions of liveability with objective spatial indicators and translate findings into planning insights.

## Study Focus

- **Study area:** Queenstown (primary), with multi-district support for Marymount (Bishan), Outram, Bukit Merah, Newton
- **Spatial boundary:** Defined from the URA Master Plan via `src/extract_district.py`
- **Supported districts:** Queenstown (15 subzones), Marymount/Bishan (3), Outram (4), Bukit Merah (17), Newton (6)

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

Island-wide GeoJSON datasets live in `data/data-gov-sg/` (gitignored, ~11.5 MB). Five Python scripts in the `src/` directory process them into committed outputs under `docs/`. All scripts that accept `--district` default to `queenstown`.

| Script | Input | Output | Run in CI |
|---|---|---|---|
| `src/extract_district.py` | URA master plan subzones GeoJSON | `docs/geo/{district}-boundary.geojson` + `{district}-subzones.geojson` | Yes |
| `src/generate_subzone_summary.py` | Subzones + point/line/CSV + remote sensing + walkability | `docs/geo/{district}-subzone-summary.geojson` + `.csv` | Yes |
| `src/filter_district_layers.py` | Boundary + 16 island-wide GeoJSONs | 16 `docs/geo/gov-sg/{district}-*.geojson` files | Yes |
| `src/fetch_global_layers.py` | OSMnx, GEE, Mapillary APIs | `docs/geo/global/{district}-*.geojson` files | Manual (requires conda env + API keys). Default tasks: `buildings osmnx gee mapillary`. Additional tasks: `raster_pngs` (GEE raster PNG overlays), `rs_grid` (100m RS grid cells) |
| `src/generate_walkability.py` | OSMnx walk network + subzones + GEE SRTM elevation | `docs/geo/global/{district}-walkability.json` + `{district}-walkability-grid.geojson` | Manual (requires conda env + GEE auth) |

The CI workflow (`.github/workflows/deploy.yml`) loops over all 5 districts for the first three scripts (calling `python src/scriptname.py`). `src/fetch_global_layers.py` and `src/generate_walkability.py` are run manually with the `zensvi` conda environment (osmnx, networkx, shapely). `src/fetch_global_layers.py` additionally requires GEE authentication + a Mapillary API key. `src/fetch_global_layers.py` has 4 default tasks: `buildings` (OSM footprints + HDB enrichment), `osmnx` (street network), `gee` (remote sensing), `mapillary` (street-level imagery). Two additional tasks: `raster_pngs` (GEE raster PNG overlays) and `rs_grid` (100m RS grid cells). Buildings, street network, remote sensing (subzone-level + RS grid + raster PNGs), and walkability data now exist for all 5 districts. Mapillary data exists for Queenstown, Bishan, and Outram (not yet Bukit Merah/Newton). Note: the Mapillary layer has been removed from the viewer but the data files and script task remain.

### Private Directory

The `private/` directory (gitignored) contains literature review files and site visit reports that are not published to GitHub Pages.

### Processed Layers in `docs/geo/`

Files are organised into subdirectories by data provenance. All filenames use `{district}-` prefix (e.g., `queenstown-`, `bishan-`).

- **`docs/geo/`** (root) — per-district boundary, subzones, subzone-summary GeoJSON + CSV (5 districts)
- **`docs/geo/gov-sg/`** — 16 layers per district from Singapore government open data (80 files total for 5 districts)
- **`docs/geo/global/`** — OSM buildings, street network, remote sensing (subzone-level + 100m RS grid), and walkability (scores + grid) for all 5 districts. Mapillary data for Queenstown, Bishan, Outram (not shown in viewer).
- **`docs/geo/global/rasters/`** — 4 satellite raster PNGs per district from GEE (NDVI, LST, NDBI, canopy; 1024px each) for all 5 districts
`gov-sg/` point layers per district (filtered by boundary containment):
- `{district}-hawker-centres.geojson`, `{district}-parks.geojson`, `{district}-mrt-exits.geojson`, `{district}-supermarkets.geojson`, `{district}-gyms.geojson`, `{district}-community-clubs.geojson`, `{district}-preschools.geojson`, `{district}-chas-clinics.geojson`, `{district}-park-facilities.geojson`, `{district}-sport-facilities.geojson`

`gov-sg/` line layers (clipped to boundary):
- `{district}-cycling-paths.geojson`, `{district}-park-connectors.geojson`, `{district}-nparks-tracks.geojson`

`gov-sg/` polygon layers:
- `{district}-ura-height-control.geojson`, `{district}-nparks-parks-nature-reserves.geojson`, `{district}-abc-waters.geojson`

`global/` building layers (all 5 districts):
- `{district}-buildings.geojson` — queenstown (8,671), bishan (3,766), outram (1,589), bukit-merah (1,874), newton (759)
- `{district}-buildings-osm-only.geojson` — OSM-only buildings (non-HDB)
- `{district}-buildings-hdb-enriched.geojson` — HDB-enriched buildings

`global/` street network (all 5 districts):
- `{district}-street-network.geojson` — edges with betweenness centrality from OSMnx drive network: queenstown (3,788), bishan (1,363), outram (711), bukit-merah (2,447), newton (838)

`global/` remote sensing — subzone-level (all 5 districts):
- `{district}-remote-sensing.geojson` — subzone polygons with LST, NDVI, NDBI, GHSL height, DSM, DEM, canopy cover from GEE: queenstown (15), bishan (3), outram (4), bukit-merah (17), newton (6)

`global/` remote sensing — grid + rasters (all 5 districts):
- `{district}-rs-grid.geojson` — clickable 100m grid cells with per-cell NDVI, NDBI, LST, canopy %, DEM elevation, slope from GEE: queenstown (2,394), bishan (840), outram (175), bukit-merah (1,616), newton (265). 6 properties: ndvi, ndbi, lst_c, canopy_pct, dem_m (SRTM 30m), slope_deg (terrain slope).
- `rasters/{district}-ndvi.png`, `{district}-lst.png`, `{district}-ndbi.png`, `{district}-canopy.png` (1024px georeferenced satellite imagery from GEE, all 5 districts)

`global/` Mapillary (Queenstown, Bishan, Outram — data files exist but layer removed from viewer):
- `{district}-mapillary-gvi.geojson` — sampled street-level image point locations: queenstown (5,000), bishan (5,000), outram (5,000). Not yet available for Bukit Merah or Newton.

`global/` walkability (all 5 districts):
- `{district}-walkability.json` — subzone-level BEH-NWI walkability scores (intersection_density, transit_access_score, destination_accessibility, walkability_index, transit_access_slope, dest_access_slope, walkability_slope): queenstown (15), bishan (3), outram (4), bukit-merah (17), newton (6). Slope-penalised variants use Tobler's hiking function with SRTM elevation.
- `{district}-walkability-grid.geojson` — clickable 100m grid cells with per-cell pop_density, transit_access, dest_access, int_density, walkability, transit_slope, dest_slope, slope_walkability: queenstown (2,178), bishan (769), outram (135), bukit-merah (1,460), newton (210). Walk network edges include elevation_start, gradient, tobler_time.

### Tech Catalogue (`docs/tech-catalogue.html`)

Documents all technologies, libraries, and frameworks used in the project: Python geospatial stack (Shapely, GeoPandas, OSMnx, NetworkX, GEE, ZenSVI, NumPy), frontend (deck.gl, MapLibre GL JS, CARTO basemaps), APIs (GEE, Mapillary, Overpass), data formats (GeoJSON, CSV, PNG), CI/CD (GitHub Actions, GitHub Pages), script reference, and environment setup.

### Glossary & Data Dictionary (`docs/glossary.html`)

Definitions of acronyms and technical terms used across the project, organised in 5 sections: Remote Sensing & Terrain (DEM, DSM, SRTM, NDVI, NDBI, LST, canopy, GHSL, GEE, slope, RS grid, NIR/SWIR), Walkability & Network Analysis (BEH-NWI, intersection density, transit access, destination accessibility, Tobler's hiking function, slope-adjusted walkability, betweenness centrality, Dijkstra, walk network, UNA, gravity accessibility), Geospatial & Data Formats (GeoJSON, OSM/OSMnx, WGS84, EPSG codes, POI, GVI), Singapore-Specific Terms (HDB, URA, MRT, CHAS, NParks, ABC Waters, LTA, subzone, planning area, hawker centre), and Statistics & Visualisation (z-score, choropleth, colour ramp, CNT/RATE/IDX/SAT badges, feature count).

### User Manual (`docs/user-manual.html`)

End-user guide for the deck.gl map viewer. Covers map navigation (pan, zoom, rotate, tilt), the district selector (5 districts), basemap styles (3 CARTO options), the layer panel structure (10 collapsible groups, provenance filter, subzone filter, feature counts, lazy loading), a full reference table of all 29 toggleable layers across 10 groups (34 with `?experimental=true`), the subzone filter (per-subzone checkboxes, All/None toggles, X/Y selected badge, spatial filtering of layers by containment), click-popup content by layer type, the 33 choropleth metrics across 5 sub-categories (39 metrics across 6 sub-categories with `?experimental=true`; type badge definitions: CNT, RATE, IDX, SAT), the 5 RS grid colour modes (NDVI, DEM, Slope, LST, NDBI), the 5 UNA gravity accessibility layers (Queenstown only, hidden behind `?experimental=true`), shareable URL state (Copy link button, 7 URL parameters: `d`, `layers`, `sz`, `choro`, `rs`, `base`, `experimental`), legend descriptions, and a tips/troubleshooting table. Includes 9 annotated screenshots in `docs/images/` (viewer overview, info panel, layer panel, expanded group, choropleth panel, choropleth map, RS grid modes, click popup, district switch). Styled to match `glossary.html` and `tech-catalogue.html`.

### deck.gl Viewer (`docs/viewer-deckgl.html`)

Primary map viewer using deck.gl + MapLibre GL JS. Located at `docs/viewer-deckgl.html`.

- **Stack**: deck.gl v9 CDN bundle (~700 KB) + MapLibre GL JS v4 for basemaps. No build tools, no API token required.
- **Buildings**: Extruded natively from GeoJSON via `GeoJsonLayer` with `extruded: true`. Height-based colour coding.
- **District selector** dropdown (Queenstown, Marymount (Bishan), Outram, Bukit Merah, Newton) with `switchDistrict()` function and `FlyToInterpolator` camera animation.
- **29 toggleable layers** by default, grouped by 10 collapsible categories: District, Food & Daily Needs, Transit, Green & Recreation, Active Mobility, Community, Housing, Planning, Street Network, Remote Sensing. The District group contains 4 items: District boundary, Subzones, Buildings (extruded), and Other districts (gray outlines + name labels for the 4 non-selected districts, on by default, re-fetched on district switch; click shows district name with prompt to use district selector). Additionally includes a Walkability grid layer. The 5 UNA gravity accessibility layers are hidden by default and available via `?experimental=true` URL parameter (34 layers total when enabled).
- **Layer panel** (top-right): collapsible groups with chevron indicators, feature counts per layer, zero-count dimming, provenance filter bar (Gov.sg / Global), source badges (SG / globe icon), RS grid colour mode selector (NDVI, DEM, Slope, LST, NDBI), and Data Catalogue footer link.
- **Subzone filter**: collapsible "Subzone Filter" section in the layer panel (after provenance bar). Users can check/uncheck individual subzones to spatially filter all point, line, polygon, grid, building, HDB, and gravity layers by subzone containment (point-in-polygon). "All" / "None" quick-toggle buttons and an "X/Y selected" count badge. Unselected subzones are visually dimmed on the map (dark overlay + faded labels). Boundary, subzones outline, other districts, and raster overlays are not filtered. Choropleth only colours selected subzones. Filter resets to "all" on district switch.
- **33 choropleth metrics** by default in 5 sub-categories: Demographics & Housing (8), Amenities & Infrastructure (9), Buildings (3), Remote Sensing (6), Walkability (7). Same YlOrRd 5-step ramp with type badges (CNT/RATE/IDX/SAT). With `?experimental=true`, the UNA Gravity sub-category (6 metrics) is also shown, totalling 39 metrics across 6 sub-categories.
- **3 CARTO basemaps**: Dark (dark-matter, default), Light (positron), Voyager. No satellite option (no free token-less provider available).
- **Data paths**: Relative to `docs/` -- e.g., `geo/queenstown-boundary.geojson`, `geo/global/queenstown-buildings.geojson`. The `layerFile()` helper swaps district names in paths.
- **Raster overlays**: 4 satellite PNGs (NDVI, LST, NDBI, canopy) rendered via `BitmapLayer` with georeferenced bounds per district.
- **Click interaction**: `onMapClick` handler with `buildPopupHtml()` for contextual popups (buildings, HDB, RS grid, walkability grid, street network, height control, choropleth, gravity layers, other districts).
- **Lazy loading**: GeoJSON layers fetched on first checkbox toggle with URL-keyed cache (`dataCache`). Default layers (boundary, subzones, buildings, other districts) loaded on init.
- **Shareable URL state**: The viewer serializes its state into URL query parameters via `history.replaceState` on every state change, enabling shareable links. A "Copy link" button in the info panel (below basemap toggle) copies the current URL to the clipboard. URL parameters (only non-default values are included):
  - `d` — district key (e.g. `d=bishan`; omitted when Queenstown since it is the default)
  - `layers` — comma-separated visible layer IDs (e.g. `layers=boundary,subzones,buildings,hawker`; omitted when default set)
  - `sz` — comma-separated selected subzone names for the subzone filter (e.g. `sz=QUEENSWAY,TANGLIN%20HALT`; omitted when all selected)
  - `choro` — active choropleth metric key (e.g. `choro=population_density`; omitted when off)
  - `rs` — RS grid colour mode (e.g. `rs=lst`; omitted when default `ndvi`)
  - `base` — basemap index 0-2 (e.g. `base=1` for Light; omitted when 0/Dark)
  - `experimental` — `true` to enable UNA gravity layers (existing parameter)
  Key functions: `serializeState()` builds the query string, `restoreStateFromUrl()` parses URL params on page load and sets state before panel/layer init, `updateUrl()` called on every state change (layer toggle, subzone filter, choropleth, basemap, RS grid mode, district switch). `loadDefaultLayers()` loads all URL-specified visible layers (not just hardcoded defaults). `populateSubzoneFilter()` respects pre-set `selectedSubzones` from URL. Choropleth activation is deferred until after default layers finish loading if `choro` param was set.
