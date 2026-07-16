#!/usr/bin/env bash

# mk_coordinate_list.sh

# Example for Asama:
# uv run python ../../scripts/kml_to_csv.py \
#   --kml_in "asama.kml" \
#   --EPSG 6676 \
#   --md_out "asama.md" \
#   --fill_elevation_from_gsi \
#   --verbose

# Example for GeoJSON from GSI Maps (chiriin chizu):
# uv run python ../../../scripts/geojson_to_csv.py \
#   --geojson_in "tokachi02.geojson" \
#   --EPSG 6680 \
#   --fill_elevation_from_gsi \
#   --csv_out "tokachi02_table_epsg6680_gsi.csv" \
#   --verbose

EPSG_CODE=$epsg

# GSI API only
uv run python ../../../scripts/$kml_or_geojson_script \
  --$input_flag "$detector_input_basename" \
  --EPSG $$EPSG_CODE \
  --fill_elevation_from_gsi \
  --csv_out "${site_name}_table_epsg$${EPSG_CODE}_gsi.csv" \
  --verbose

# DEM only
uv run python ../../../scripts/kml_to_csv.py \
  --kml_in "$detector_input_basename" \
  --EPSG $$EPSG_CODE \
  --dem "$dem_relative_path" \
  --csv_out "${site_name}_table_epsg$${EPSG_CODE}_g2zbin.csv" \
  --verbose
