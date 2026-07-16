#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
# usage :
# Treat 0.0 as missing and fill it in (for the Asama case)
python3 kml_to_csv.py \
  --kml_in "asama.kml" \
  --EPSG 6676 \
  --md_out "asama.md" \
  --fill_elevation_from_gsi \
  --treat_zero_elev_as_missing \
  --throttle_sec 0.5 \
  --verbose

# Ignore KML elevations and fetch all points via the API
python3 kml_to_csv.py \
  --kml_in "asama.kml" \
  --EPSG 6676 \
  --fill_elevation_from_gsi \
  --force_gsi_for_all \
  --throttle_sec 0.5 \
  --verbose

# Fall back to DEM when the GSI API fails (auto-detects .tif / .g2zbin)
python3 kml_to_csv.py \
  --kml_in "iwojima.kml" \
  --EPSG 6670 \
  --fill_elevation_from_gsi \
  --force_gsi_for_all \
  --dem "path/to/dem.tif" \
  --verbose

# Get elevation from DEM (GSI API is used as fallback)
python3 kml_to_csv.py \
  --kml_in "tarumae.kml" \
  --EPSG 6680 \
  --dem "path/to/dem.g2zbin" \
  --verbose
"""

import xml.etree.ElementTree as ET
try:
  import pandas as pd
except ImportError:
  sys.exit(
    "ERROR: pandas is required.\n"
    "  pip install pandas  (or activate venv: source .venv/bin/activate)"
  )
import numpy as np
import math
import sys
import os
try:
  from pyproj import Transformer
except ImportError:
  sys.exit(
    "ERROR: pyproj is required for kml_to_csv.py\n"
    "  pip install pyproj  (or activate venv: source .venv/bin/activate)"
  )
import argparse
import re
import time
from typing import Optional, Tuple
try:
  from IPython.display import display
except ImportError:
  def display(x):
    print(x)
try:
  import requests
except ImportError:
  sys.exit(
    "ERROR: requests is required for kml_to_csv.py\n"
    "  pip install requests  (or activate venv: source .venv/bin/activate)"
  )

from g2zbin_io import axis_centers, read_g2zbin

try:
  import rasterio
  HAS_RASTERIO = True
except ImportError:
  HAS_RASTERIO = False
  import warnings
  warnings.warn(
    "rasterio module not found. DEM elevation lookup will be unavailable.\n"
    "  Hint: activate the venv — source .venv/bin/activate",
    stacklevel=2,
  )

# ===========================================
# Static variables (hardcoded items)
# ===========================================

# KML namespace
KML_NS = {
  "kml": "http://www.opengis.net/kml/2.2",
  "gx": "http://www.google.com/kml/ext/2.2"
}

# Priority list of keywords considered to indicate the reference point.
# A placemark whose name contains any of these is taken as the reference point.
REFERENCE_KEYWORDS = [
  "vent", "crater", "summit", "center", "main", "reference", "target"
]

# GSI (Geospatial Information Authority of Japan) elevation API (JSON)
GSI_ELEVATION_API_URL = "https://cyberjapandata2.gsi.go.jp/general/dem/scripts/getelevation.php"
GSI_ELEVATION_OUTTYPE = "JSON"

# ===========================================
# Function definitions
# ===========================================

def parse_kml_points(path):
  """Extract Point placemarks from a KML file."""
  tree = ET.parse(path)
  root = tree.getroot()
  placemarks = root.findall(".//kml:Placemark", KML_NS)

  points = []
  for pm in placemarks:
    name_el = pm.find("kml:name", KML_NS)
    name = name_el.text.strip() if (name_el is not None and name_el.text) else ""

    point_el = pm.find(".//kml:Point", KML_NS)
    if point_el is None:
      continue
    coords_el = point_el.find("kml:coordinates", KML_NS)
    if coords_el is None or not coords_el.text or not coords_el.text.strip():
      continue

    raw = coords_el.text.strip().split()
    first = raw[0].split(",")
    try:
      lon = float(first[0])
      lat = float(first[1])
      alt = float(first[2]) if len(first) >= 3 else np.nan
      points.append({"name": name, "lon": lon, "lat": lat, "elev": alt})
    except Exception:
      continue
  return points


def get_elevation_gsi(lon: float, lat: float, timeout: float = 5.0,
                      session: Optional[requests.Session] = None) -> Tuple[Optional[float], Optional[str]]:
  """
  Get elevation from the GSI (Geospatial Information Authority of Japan) elevation API.
  Returns: (elevation[m] or None, hsrc or None)
  Spec: https://maps.gsi.go.jp/development/elevation_s.html
  """
  s = session or requests.Session()
  params = {"lon": lon, "lat": lat, "outtype": GSI_ELEVATION_OUTTYPE}
  try:
    r = s.get(GSI_ELEVATION_API_URL, params=params, timeout=timeout)
    r.raise_for_status()
    data = r.json()
  except Exception:
    return None, None

  elev = data.get("elevation")
  hsrc = data.get("hsrc")

  # Spec: "-----" is returned on error
  if isinstance(elev, str) and set(elev) == {"-"}:
    return None, None

  try:
    elev_f = float(elev)
  except Exception:
    return None, None

  return elev_f, hsrc


def fill_missing_elevation_with_gsi(points, throttle_sec: float = 0.5,
                                    timeout: float = 5.0, verbose: bool = False):
  """
  Fill in values via the GSI elevation API for points whose elev is NaN.
  Returns: a simple list of points that failed to fetch (empty list if none failed)
  throttle_sec: interval between API calls (seconds)
  """
  s = requests.Session()
  failures = []
  for p in points:
    if np.isnan(p.get("elev", np.nan)):
      elev, hsrc = get_elevation_gsi(p["lon"], p["lat"], timeout=timeout, session=s)
      if elev is not None:
        p["elev"] = elev
        p["hsrc"] = hsrc
        if verbose:
          print(f'filled: name="{p.get("name","")}" elev={elev} hsrc={hsrc}')
      else:
        p["elev"] = -9999.0
        p["hsrc"] = None
        if verbose:
          print(f'failed : name="{p.get("name","")}" elev=-9999 (API failure)')
        failures.append({"name": p.get("name",""), "lon": p["lon"], "lat": p["lat"]})
      time.sleep(throttle_sec)
  return failures


def get_elevation_from_geotiff(lon: float, lat: float, dem_dataset, transformer_to_dem) -> Optional[float]:
  """
  Get elevation from a GeoTIFF.
  lon, lat: WGS84 latitude/longitude
  dem_dataset: rasterio dataset object
  transformer_to_dem: pyproj Transformer (EPSG:4326 -> DEM CRS)
  Returns: elevation[m] or None (out of range or NoData)
  """
  # Convert lon/lat to the DEM's coordinate system
  x_dem, y_dem = transformer_to_dem.transform(lon, lat)

  # Compute pixel coordinates
  inv_transform = ~dem_dataset.transform
  col, row = inv_transform * (x_dem, y_dem)
  col_int, row_int = int(round(col)), int(round(row))

  # Range check
  if not (0 <= row_int < dem_dataset.height and 0 <= col_int < dem_dataset.width):
    return None

  # Get the elevation value
  elev = dem_dataset.read(1)[row_int, col_int]

  # NoData check
  if dem_dataset.nodata is not None and elev == dem_dataset.nodata:
    return None

  return float(elev)


def fill_elevation_from_geotiff(points, dem_path: str, verbose: bool = False):
  """
  Fill in elevation via GeoTIFF for points whose elev is -9999.
  Returns: number of points that failed to fetch
  """
  if not HAS_RASTERIO:
    if verbose:
      print("[WARN] rasterio is not installed, skipping GeoTIFF elevation fill")
    return len([p for p in points if p.get("elev") == -9999.0])

  failed_count = 0
  with rasterio.open(dem_path) as src:
    dem_crs = src.crs
    if verbose:
      print(f"[DEBUG] DEM CRS: {dem_crs}, size: {src.width}x{src.height}")

    # Create a transformer for EPSG:4326 -> DEM CRS
    transformer_to_dem = Transformer.from_crs(4326, dem_crs, always_xy=True)

    for p in points:
      e = p.get("elev")
      if e == -9999.0 or (isinstance(e, float) and np.isnan(e)):
        elev = get_elevation_from_geotiff(p["lon"], p["lat"], src, transformer_to_dem)
        if elev is not None:
          p["elev"] = elev
          p["hsrc"] = f"GeoTIFF"
          if verbose:
            print(f'filled from DEM: name="{p.get("name","")}" elev={elev:.2f}')
        else:
          failed_count += 1
          if verbose:
            print(f'failed from DEM: name="{p.get("name","")}" (out of bounds or NoData)')

  return failed_count


def get_elevation_from_g2zbin(lon, lat, grid2d_info, Z, xi, yi, transformer_to_dem):
  """
  Get elevation from a g2zbin grid.
  lon, lat: WGS84 coordinates
  grid2d_info: dict from read_g2zbin
  Z: 2D numpy array (ny, nx)
  xi, yi: 1D coordinate arrays
  transformer_to_dem: pyproj Transformer (EPSG:4326 -> DEM CRS)
  Returns: elevation[m] or None (out of bounds)
  """
  x_dem, y_dem = transformer_to_dem.transform(lon, lat)

  ax = grid2d_info["x_axis"]
  ay = grid2d_info["y_axis"]

  # Containing-bin index (v2 min is the lower bin edge, so floor is nearest)
  ix = int((x_dem - ax["min"]) // ax["interval"])
  iy = int((y_dem - ay["min"]) // ay["interval"])

  if not (0 <= ix < ax["nbin"] and 0 <= iy < ay["nbin"]):
    return None

  return float(Z[iy, ix])


def fill_elevation_from_g2zbin(points, dem_path, epsg_dem, verbose=False):
  """
  Fill missing elevations (NaN or -9999) from a g2zbin DEM file.
  Returns: number of points that could not get elevation.
  """
  grid2d_info, Z = read_g2zbin(dem_path)
  ax = grid2d_info["x_axis"]
  ay = grid2d_info["y_axis"]
  # Bin centers (v2 min/max are bin edges)
  xi = axis_centers(ax)
  yi = axis_centers(ay)

  transformer_to_dem = Transformer.from_crs(4326, epsg_dem, always_xy=True)

  if verbose:
    print(f"[DEBUG] g2zbin DEM: {dem_path}")
    print(f"[DEBUG]   x: min={ax['min']:.1f} max={ax['max']:.1f} nbin={ax['nbin']} interval={ax['interval']:.1f}")
    print(f"[DEBUG]   y: min={ay['min']:.1f} max={ay['max']:.1f} nbin={ay['nbin']} interval={ay['interval']:.1f}")

  failed_count = 0
  for p in points:
    e = p.get("elev")
    if e == -9999.0 or (isinstance(e, float) and np.isnan(e)):
      elev = get_elevation_from_g2zbin(
        p["lon"], p["lat"], grid2d_info, Z, xi, yi, transformer_to_dem
      )
      if elev is not None:
        p["elev"] = elev
        p["hsrc"] = "g2zbin"
        if verbose:
          print(f'filled from g2zbin: name="{p.get("name","")}" elev={elev:.2f}')
      else:
        failed_count += 1
        if verbose:
          print(f'failed from g2zbin: name="{p.get("name","")}" (out of bounds)')

  return failed_count


def compute_table(points, epsg_out):
  """Compute projected coordinates and relative metrics."""
  if len(points) == 0:
    cols = [
      "det_name", "name", "lon", "lat",
      f"X (m, EPSG:{epsg_out})", f"Y (m, EPSG:{epsg_out})",
      "elevation(m)", "azimuth(deg)", "horizontal_distance(m)",
      "elevation_diff(m)", "elevation_angle(deg)", "tangent_elevation_angle"
    ]
    return pd.DataFrame(columns=cols)

  transformer = Transformer.from_crs(4326, epsg_out, always_xy=True)
  Xs, Ys = transformer.transform(
    [p["lon"] for p in points],
    [p["lat"] for p in points]
  )
  for i, p in enumerate(points):
    p["X"], p["Y"] = float(Xs[i]), float(Ys[i])

  # Determine the reference point. Azimuth, horizontal distance and elevation
  # difference of every detector are measured against it, so a wrong reference
  # point silently corrupts every row: fail instead of falling back.
  ref_idx = next(
    (i for i, p in enumerate(points)
     if any(kw in (p["name"] or "") for kw in REFERENCE_KEYWORDS)),
    None
  )
  if ref_idx is None:
    raise SystemExit(
      "No reference point found. Name one placemark after any of the following "
      "keywords (substring match, case-sensitive): "
      + ", ".join(REFERENCE_KEYWORDS)
      + ".\n  Placemark names found: "
      + ", ".join(repr(p["name"]) for p in points)
    )
  ref = points[ref_idx]

  # Generate table rows
  rows = []
  rows.append({
    "det_name": "-----",
    "name": ref["name"] or "reference_point",
    "lon": round(ref["lon"], 8),
    "lat": round(ref["lat"], 8),
    f"X (m, EPSG:{epsg_out})": round(ref["X"], 1),
    f"Y (m, EPSG:{epsg_out})": round(ref["Y"], 1),
    "elevation(m)": None if np.isnan(ref.get("elev", np.nan)) else round(ref["elev"], 1),
    "azimuth(deg)": 0,
    "horizontal_distance(m)": 0,
    "elevation_diff(m)": 0,
    "elevation_angle(deg)": 0,
    "tangent_elevation_angle": 0
  })

  det_counter = 0
  for i, p in enumerate(points):
    if i == ref_idx:
      continue

    name = p["name"] or f"candidate{det_counter:02d}"
    m = re.search(r"(\d+)", name)
    det_name = f"det_{int(m.group(1)):02d}" if m else f"det_{det_counter:02d}"
    det_counter += 1

    dx = ref["X"] - p["X"]
    dy = ref["Y"] - p["Y"]
    horiz = math.hypot(dx, dy)
    az_deg = math.degrees(math.atan2(dx, dy))

    if not np.isnan(ref.get("elev", np.nan)) and not np.isnan(p.get("elev", np.nan)):
      dz = ref["elev"] - p["elev"]
      tan_elev = dz / horiz if horiz != 0 else np.nan
      elev_deg = math.degrees(math.atan(tan_elev)) if horiz != 0 else np.nan
      dz_out, elev_deg_out, tan_out = (
        round(dz, 1),
        round(elev_deg, 6),
        round(tan_elev, 6)
      )
    else:
      dz_out = elev_deg_out = tan_out = ""

    rows.append({
      "det_name": det_name,
      "name": name,
      "lon": round(p["lon"], 8),
      "lat": round(p["lat"], 8),
      f"X (m, EPSG:{epsg_out})": round(p["X"], 1),
      f"Y (m, EPSG:{epsg_out})": round(p["Y"], 1),
      "elevation(m)": None if np.isnan(p.get("elev", np.nan)) else round(p["elev"], 1),
      "azimuth(deg)": round(az_deg, 3),
      "horizontal_distance(m)": round(horiz, 1),
      "elevation_diff(m)": dz_out,
      "elevation_angle(deg)": elev_deg_out,
      "tangent_elevation_angle": tan_out
    })

  cols = [
    "det_name","name","lon","lat",
    f"X (m, EPSG:{epsg_out})",f"Y (m, EPSG:{epsg_out})",
    "elevation(m)","azimuth(deg)","horizontal_distance(m)",
    "elevation_diff(m)","elevation_angle(deg)","tangent_elevation_angle"
  ]
  return pd.DataFrame(rows, columns=cols)


def df_to_markdown_table(df):
  """
  Markdown table generator that does not depend on Pandas or tabulate.
  - Header row
  - Separator row (left-aligned)
  - Data rows
  """
  def _fmt(v):
    if v is None:
      return ""
    if isinstance(v, float):
      s = f"{v:.8f}".rstrip("0").rstrip(".")
      return s
    return str(v)

  cols = list(df.columns)
  header = "| " + " | ".join(cols) + " |"
  sep = "| " + " | ".join(["---"] * len(cols)) + " |"

  lines = [header, sep]
  for _, row in df.iterrows():
    fields = [_fmt(row[c]) for c in cols]
    lines.append("| " + " | ".join(fields) + " |")
  return "\n".join(lines) + "\n"


def main():
  parser = argparse.ArgumentParser(description="Convert KML points to analysis table (CSV + Markdown, optional GSI elevation fill).")
  parser.add_argument("--kml_in", required=True, help="Input KML file path")
  parser.add_argument("--EPSG", required=True, type=int, help="Output EPSG code (e.g., 6676)")
  parser.add_argument("--csv_out", required=False, help="Output CSV file path (*.csv)")
  parser.add_argument("--md_out", required=False, help="Output Markdown file path (*.md)")
  parser.add_argument("--fill_elevation_from_gsi", action="store_true", help="Fetch ALL elevations from GSI API (ignores KML elevations)")
  parser.add_argument("--throttle_sec", type=float, default=0.5, help="Throttle seconds per API call (default=0.5)")
  parser.add_argument("--timeout", type=float, default=5.0, help="HTTP timeout for API calls (seconds)")
  parser.add_argument("--verbose", action="store_true", help="Verbose log for elevation filling")
  parser.add_argument("--dem", type=str, default=None,
                      help="Path to DEM file (.tif/.tiff or .g2zbin). Ignores KML elevations and uses DEM only (auto-detect format by extension)")
  args = parser.parse_args()

  # Resolve DEM format from extension
  args._dem_is_g2zbin = False
  if args.dem:
    ext = os.path.splitext(args.dem)[1].lower()
    if ext in (".tif", ".tiff"):
      pass
    elif ext == ".g2zbin":
      args._dem_is_g2zbin = True
    else:
      print(f"[ERROR] Unsupported DEM format: {ext} (expected .tif, .tiff, or .g2zbin)")
      sys.exit(1)

  points = parse_kml_points(args.kml_in)

  # Helper: fill from DEM (GeoTIFF or g2zbin)
  def _fill_from_dem(pts, verbose):
    if args._dem_is_g2zbin:
      return fill_elevation_from_g2zbin(pts, args.dem, args.EPSG, verbose=verbose)
    else:
      return fill_elevation_from_geotiff(pts, args.dem, verbose=verbose)

  # --dem mode: use DEM for all elevations, GSI API as fallback
  if args.dem:
    for p in points:
      p["elev"] = np.nan
    if args.verbose:
      print(f"[DEBUG] Fetching all elevations from DEM: {args.dem}")
    still_failed = _fill_from_dem(points, verbose=args.verbose)
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
  elif args.fill_elevation_from_gsi:
    # Clear all elevations to force GSI API fetch
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

  df = compute_table(points, args.EPSG)

  base = args.kml_in.rsplit(".kml", 1)[0]
  out_csv = args.csv_out if args.csv_out else f"{base}_table_epsg{args.EPSG}.csv"
  out_md = args.md_out if args.md_out else f"{base}_table_epsg{args.EPSG}.md"

  df.to_csv(out_csv, index=False, encoding="utf-8-sig")
  md = df_to_markdown_table(df)
  with open(out_md, "w", encoding="utf-8") as f:
    f.write(md)

  print("✅ Output complete:")
  print(f"  CSV: {out_csv} ({len(df)} rows)")
  print(f"  MD : {out_md}")
  display(df)


if __name__ == "__main__":
  main()
