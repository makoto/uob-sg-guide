# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

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

Island-wide GeoJSON datasets live in `data/data-gov-sg/` (gitignored, ~11.5 MB). Three Python scripts process them into committed outputs under `docs/`:

| Script | Input | Output | Run in CI |
|---|---|---|---|
| `generate_subzone_summary.py` | Subzones + point/line/CSV datasets | `docs/geo/queenstown-subzone-summary.geojson` + `.csv` | Yes |
| `filter_queenstown_layers.py` | Boundary + 11 island-wide GeoJSONs | 11 `docs/geo/queenstown-*.geojson` files (~373 KB) | Yes |
| `generate_3dtiles.py` | Buildings GeoJSON | `docs/3dtiles/tileset.json` + B3DM | Yes |

All three run in `.github/workflows/deploy.yml` on push to `main`.

### Processed Layers in `docs/geo/`

Point layers (filtered by boundary containment):
- `queenstown-hawker-centres.geojson` (9), `queenstown-parks.geojson` (9), `queenstown-mrt-exits.geojson` (26), `queenstown-supermarkets.geojson` (19), `queenstown-gyms.geojson` (4), `queenstown-community-clubs.geojson` (4), `queenstown-preschools.geojson` (81), `queenstown-chas-clinics.geojson` (43), `queenstown-park-facilities.geojson` (370)

Line layers (clipped to boundary):
- `queenstown-cycling-paths.geojson` (110), `queenstown-park-connectors.geojson` (36)

### 3D Viewer (`docs/3dtiles/viewer.html`)

CesiumJS-based viewer with:
- 3D building tileset (height-colored)
- Boundary + subzone overlays (on by default)
- **Layer panel** (top-right): 13 toggleable overlays grouped by category (Food & Daily Needs, Transit, Green & Recreation, Active Mobility, Community)
- Layers are lazy-loaded on first checkbox toggle via `GeoJsonDataSource`
- Basemap toggle (dark/light CARTO)
