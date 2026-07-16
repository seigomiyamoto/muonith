#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Convert GSI Maps GeoJSON (FeatureCollection of Points) to analysis table (CSV + Markdown).

Usage:
  # GSI API elevation
  python3 geojson_to_csv.py \
    --geojson_in "tokachi.geojson" \
    --EPSG 6680 \
    --fill_elevation_from_gsi \
    --verbose

  # DEM only
  python3 geojson_to_csv.py \
    --geojson_in "tokachi.geojson" \
    --EPSG 6680 \
    --dem_only \
    --dem_tif "../dem/tokachi01_3x3km.tif" \
    --verbose
"""

import json
import sys
import argparse
import numpy as np

from kml_to_csv import (
  fill_missing_elevation_with_gsi,
  fill_elevation_from_geotiff,
  compute_table,
  df_to_markdown_table,
)

try:
  from IPython.display import display
except ImportError:
  def display(x):
    print(x)


def parse_geojson_points(path):
  """Extract Point features from a GeoJSON FeatureCollection file."""
  with open(path, "r", encoding="utf-8") as f:
    data = json.load(f)

  features = data.get("features", [])
  if data.get("type") == "Feature":
    features = [data]

  points = []
  for feat in features:
    geom = feat.get("geometry")
    if geom is None or geom.get("type") != "Point":
      continue

    coords = geom.get("coordinates", [])
    if len(coords) < 2:
      continue

    try:
      lon = float(coords[0])
      lat = float(coords[1])
      alt = float(coords[2]) if len(coords) >= 3 else np.nan
    except (ValueError, TypeError):
      continue

    props = feat.get("properties", {})
    name = props.get("name") or props.get("title") or ""
    if isinstance(name, str):
      name = name.strip()

    points.append({"name": name, "lon": lon, "lat": lat, "elev": alt})

  return points


def main():
  parser = argparse.ArgumentParser(
    description="Convert GSI Maps GeoJSON points to analysis table (CSV + Markdown)."
  )
  parser.add_argument("--geojson_in", required=True, help="Input GeoJSON file path")
  parser.add_argument("--EPSG", required=True, type=int, help="Output EPSG code (e.g., 6680)")
  parser.add_argument("--csv_out", required=False, help="Output CSV file path (*.csv)")
  parser.add_argument("--md_out", required=False, help="Output Markdown file path (*.md)")
  parser.add_argument("--fill_elevation_from_gsi", action="store_true",
                      help="Fetch ALL elevations from GSI API (ignores GeoJSON elevations)")
  parser.add_argument("--throttle_sec", type=float, default=0.5,
                      help="Throttle seconds per API call (default=0.5)")
  parser.add_argument("--timeout", type=float, default=5.0,
                      help="HTTP timeout for API calls (seconds)")
  parser.add_argument("--verbose", action="store_true", help="Verbose log for elevation filling")
  parser.add_argument("--dem_tif", type=str, default=None,
                      help="Path to GeoTIFF DEM file for fallback elevation lookup")
  parser.add_argument("--dem_only", action="store_true",
                      help="Use DEM only for elevation (ignore GeoJSON elevations, skip GSI API)")
  args = parser.parse_args()

  points = parse_geojson_points(args.geojson_in)
  if not points:
    print("[ERROR] No Point features found in the GeoJSON file.")
    sys.exit(1)

  if args.verbose:
    print(f"[INFO] Parsed {len(points)} point(s) from {args.geojson_in}")

  # --dem_only mode: use DEM only, skip GSI API
  if args.dem_only:
    if not args.dem_tif:
      print("[ERROR] --dem_only requires --dem_tif")
      sys.exit(1)
    for p in points:
      p["elev"] = np.nan
    if args.verbose:
      print("[DEBUG] dem_only: fetching all elevations from DEM")
    still_failed = fill_elevation_from_geotiff(points, args.dem_tif, verbose=args.verbose)
    if still_failed > 0:
      print(f"[WARN] {still_failed} points could not get elevation from DEM, trying GSI API fallback...")
      failures = fill_missing_elevation_with_gsi(
        points,
        throttle_sec=args.throttle_sec,
        timeout=args.timeout,
        verbose=args.verbose
      )
      if failures:
        print(f"[WARN] GSI API fallback also failed for {len(failures)} points")
      else:
        print(f"[INFO] All {still_failed} missing elevations fetched from GSI API")
    else:
      print(f"[INFO] All {len(points)} elevations fetched from DEM")
  else:
    if args.fill_elevation_from_gsi:
      for p in points:
        p["elev"] = np.nan
      if args.verbose:
        print("[DEBUG] fill_elevation_from_gsi: cleared all elevations, fetching from GSI API")
      failures = fill_missing_elevation_with_gsi(
        points,
        throttle_sec=args.throttle_sec,
        timeout=args.timeout,
        verbose=args.verbose
      )
      if failures:
        print(f"[WARN] GSI API failed for {len(failures)} points:")
        for i, fp in enumerate(failures, 1):
          print(f"  {i:02d}: name='{fp['name']}' lon={fp['lon']} lat={fp['lat']}")

        # GeoTIFF fallback
        if args.dem_tif:
          print(f"[INFO] Trying fallback from GeoTIFF ({args.dem_tif})...")
          still_failed = fill_elevation_from_geotiff(points, args.dem_tif, verbose=args.verbose)
          if still_failed > 0:
            print(f"[WARN] {still_failed} points still missing after DEM fallback")
          else:
            print(f"[INFO] All {len(failures)} elevations fetched from DEM")

  df = compute_table(points, args.EPSG)

  base = args.geojson_in.rsplit(".geojson", 1)[0]
  out_csv = args.csv_out if args.csv_out else f"{base}_table_epsg{args.EPSG}.csv"
  out_md = args.md_out if args.md_out else f"{base}_table_epsg{args.EPSG}.md"

  df.to_csv(out_csv, index=False, encoding="utf-8-sig")
  md = df_to_markdown_table(df)
  with open(out_md, "w", encoding="utf-8") as f:
    f.write(md)

  print(f"[OK] Output:")
  print(f"  CSV: {out_csv} ({len(df)} rows)")
  print(f"  MD : {out_md}")
  display(df)


if __name__ == "__main__":
  main()
