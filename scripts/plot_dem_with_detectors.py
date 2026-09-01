#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# scripts/plot_dem_with_detectors.py
"""Plot a DEM (g2zbin or GeoTIFF) together with detector positions from a CSV file or JSON5 config."""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Sequence, Tuple

import matplotlib.patheffects as path_effects
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource
import numpy as np

from plot_dem import load_g2zbin_as_grid


# ---------------------------------------------------------------------------
# JSON5 helpers (restored from pre-aaa39a33 for --config support)
# ---------------------------------------------------------------------------

def strip_json_comments(src: str) -> str:
  """Remove // and /* */ style comments while keeping string literals intact."""
  result: List[str] = []
  i = 0
  in_string = False
  escape = False
  quote_char = ""
  length = len(src)

  while i < length:
    ch = src[i]
    if in_string:
      result.append(ch)
      if escape:
        escape = False
      elif ch == "\\":
        escape = True
      elif ch == quote_char:
        in_string = False
      i += 1
      continue

    if ch in ("'", '"'):
      in_string = True
      quote_char = ch
      result.append(ch)
      i += 1
      continue

    if ch == "/" and i + 1 < length:
      nxt = src[i + 1]
      if nxt == "/":
        i += 2
        while i < length and src[i] not in ("\n", "\r"):
          i += 1
        continue
      if nxt == "*":
        i += 2
        while i + 1 < length and not (src[i] == "*" and src[i + 1] == "/"):
          i += 1
        i += 2
        continue

    result.append(ch)
    i += 1

  return "".join(result)


def load_json5(path: Path) -> Dict[str, Any]:
  """Load a JSON5-ish file, falling back to comment stripping if json5 is unavailable."""
  text = path.read_text(encoding="utf-8")
  try:
    import json5 as _json5  # type: ignore
  except ImportError:
    _json5 = None

  if _json5:
    return _json5.loads(text)

  if str(path).endswith(".json5"):
    raise ImportError(
      f"Cannot parse JSON5 file without json5 module: {path}\n"
      "  Install: pip install json5\n"
      "  Or activate venv: source .venv/bin/activate"
    )
  return json.loads(strip_json_comments(text))


def resolve_path(base_dir: Path, candidate: str) -> Path:
  path = Path(candidate)
  if not path.is_absolute():
    path = (base_dir / path).resolve()
  return path


def load_detectors_from_config(det_files: Sequence[str], config_dir: Path) -> List[Dict[str, Any]]:
  """Load detector positions from JSON5 det_files listed in config."""
  detectors: List[Dict[str, Any]] = []
  for det_rel in det_files:
    det_path = resolve_path(config_dir, det_rel)
    det_cfg = load_json5(det_path)
    det_params = det_cfg.get("DETECTOR_PARAMETERS")
    if not det_params:
      raise KeyError(f"DETECTOR_PARAMETERS not found in: {det_path}")
    detectors.append({
        "name": det_params.get("name", det_path.stem),
        "x": det_params["x"],
        "y": det_params["y"],
    })
  return detectors


# ---------------------------------------------------------------------------


def read_geotiff(path: str):
  """Read a GeoTIFF DEM and return (xi, yi, Z) ready for plotting.

  Returns:
    (xi, yi, Z) where xi/yi are 1-D coordinate arrays and Z is 2-D elevation.
  """
  import rasterio

  with rasterio.open(path) as src:
    Z = src.read(1).astype(np.float64)
    nodata = src.nodata
    if nodata is not None:
      Z[Z == nodata] = np.nan
    transform = src.transform
    nrows, ncols = Z.shape
    # Pixel-center coordinates
    xi = np.array([transform.c + (col + 0.5) * transform.a for col in range(ncols)])
    yi = np.array([transform.f + (row + 0.5) * transform.e for row in range(nrows)])
    # Ensure yi is ascending (imshow origin="lower" expects it)
    if yi[0] > yi[-1]:
      yi = yi[::-1]
      Z = Z[::-1, :]

  return xi, yi, Z


def load_detectors_from_csv(csv_path: Path) -> Tuple[List[Dict[str, Any]], Tuple[float, float] | None]:
  """
  Load detector positions and reference point from a CSV file.

  Returns:
    (detectors, target_xy)
    - detectors: list of dicts with 'name', 'x', 'y'
    - target_xy: (x, y) of reference point (det_name == "-----"), or None
  """
  detectors: List[Dict[str, Any]] = []
  target_xy: Tuple[float, float] | None = None

  with open(csv_path, "r", encoding="utf-8-sig") as f:
    reader = csv.DictReader(f)

    # Find X/Y column names (they contain "X (m" and "Y (m")
    x_col = None
    y_col = None
    for col in reader.fieldnames or []:
      if "X (m" in col:
        x_col = col
      elif "Y (m" in col:
        y_col = col

    if not x_col or not y_col:
      raise KeyError(f"CSV must have columns containing 'X (m' and 'Y (m'. Found: {reader.fieldnames}")

    for row in reader:
      det_name = row.get("det_name", "").strip()
      name = row.get("name", "").strip()

      try:
        x = float(row[x_col])
        y = float(row[y_col])
      except (ValueError, KeyError):
        continue

      if det_name == "center" or re.fullmatch(r"-+", det_name):
        # Reference point (summit, center, or dashes like "-----", "--", etc.)
        target_xy = (x, y)
      else:
        detectors.append({
          "name": det_name or name,
          "x": x,
          "y": y,
        })

  return detectors, target_xy


def crop_to_window(xi, yi, Z, x_lo: float, x_hi: float, y_lo: float, y_hi: float):
  """Crop the DEM grid to the given window (meters).

  Args:
    xi: 1-D x coordinate array (meters).
    yi: 1-D y coordinate array (meters).
    Z: 2-D elevation array shaped (len(yi), len(xi)).
    x_lo, x_hi, y_lo, y_hi: Window edges, in meters.

  Returns:
    (xi, yi, Z) restricted to the window.
  """
  ix = np.where((xi >= x_lo) & (xi <= x_hi))[0]
  iy = np.where((yi >= y_lo) & (yi <= y_hi))[0]
  if ix.size == 0 or iy.size == 0:
    raise ValueError(
        f"Crop window is empty: x {x_lo:.1f}..{x_hi:.1f}, y {y_lo:.1f}..{y_hi:.1f} "
        f"does not overlap the DEM (x: {xi.min():.1f}..{xi.max():.1f}, "
        f"y: {yi.min():.1f}..{yi.max():.1f}).")
  return xi[ix], yi[iy], Z[iy[0]:iy[-1] + 1, ix[0]:ix[-1] + 1]


def crop_around(xi, yi, Z, center_xy: Tuple[float, float], half_range: float):
  """Crop the DEM grid to center_xy +/- half_range (meters)."""
  cx, cy = center_xy
  return crop_to_window(xi, yi, Z,
                        cx - half_range, cx + half_range,
                        cy - half_range, cy + half_range)


def add_range_circles(ax, center_xy: Tuple[float, float], radii: Sequence[float],
                      scale: float = 1.0, color: str = "black"):
  """Draw circles of the given radii (meters) centered on center_xy."""
  cx = center_xy[0] / scale
  cy = center_xy[1] / scale
  # A pale outline keeps the line readable where it crosses dark terrain.
  halo = [path_effects.withStroke(linewidth=2.5, foreground="white", alpha=0.8)]
  for radius in radii:
    r = radius / scale
    circle = plt.Circle((cx, cy), r, fill=False, color=color,
                        linewidth=1.0, linestyle="--", alpha=0.9, zorder=7,
                        path_effects=halo)
    ax.add_patch(circle)
    label = f"{radius / 1000.0:g} km"
    ax.text(cx, cy + r, label, color=color, fontsize=7,
            ha="center", va="bottom", zorder=8,
            bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.7))


def plot_background(ax, xi, yi, Z, cmap: str, zmin: float | None, zmax: float | None,
                    contour: bool, cint: float, xyunit: str = "m",
                    shade: bool = False, shade_alpha: float = 0.2,
                    shade_azdeg: float = 315.0, shade_altdeg: float = 45.0,
                    shade_vert_exag: float = 2.0):
  scale = 1000.0 if xyunit == "km" else 1.0
  xi_scaled = xi / scale
  yi_scaled = yi / scale

  extent = [float(xi_scaled.min()), float(xi_scaled.max()),
            float(yi_scaled.min()), float(yi_scaled.max())]
  im = ax.imshow(Z, extent=extent, origin="lower", aspect="equal",
                 vmin=zmin, vmax=zmax, cmap=cmap, interpolation="nearest")
  if shade:
    dx = float(np.nanmedian(np.diff(np.sort(np.unique(xi)))))
    dy = float(np.nanmedian(np.diff(np.sort(np.unique(yi)))))
    ls = LightSource(azdeg=shade_azdeg, altdeg=shade_altdeg)
    # LightSource.hillshade treats row 0 as the north edge, but this grid is drawn
    # with origin="lower", so row 0 is the south edge. Shading it as-is mirrors the
    # light north-south: with the default azdeg=315 it arrives from the lower left
    # and ridges read as valleys. Flip before shading and flip the result back so
    # shade_azdeg keeps its plain meaning (315 = light from the upper left).
    hillshade = ls.hillshade(Z[::-1, :], vert_exag=shade_vert_exag, dx=dx, dy=dy)[::-1, :]
    ax.imshow(hillshade, extent=extent, origin="lower", aspect="equal",
              cmap="gray", alpha=shade_alpha, interpolation="nearest")
  ax.set_xlabel(f"Easting ({xyunit})")
  ax.set_ylabel(f"Northing ({xyunit})")
  cbar = plt.colorbar(im, ax=ax, label="Elevation (meters above sea level)", extend="both")

  if contour:
    zmin_eff = zmin if zmin is not None else float(np.nanmin(Z))
    zmax_eff = zmax if zmax is not None else float(np.nanmax(Z))
    start = np.floor(zmin_eff / cint) * cint
    stop = np.ceil(zmax_eff / cint) * cint + 0.5 * cint
    levels = np.arange(start, stop, cint)
    cs = ax.contour(xi_scaled, yi_scaled, Z, levels=levels, colors="black",
                    linewidths=0.5, alpha=0.6)
    ax.clabel(cs, fmt="%.0f", fontsize=6, inline=True)

  return extent, cbar, scale


def name_label_position(det_x: float, det_y: float,
                        tgt_x: float | None, tgt_y: float | None,
                        offset: float, fallback_dx: float, fallback_dy: float):
  """Place a detector's name just past it, on the target-to-detector line.

  Args:
    det_x, det_y: Detector position in plot units.
    tgt_x, tgt_y: Target position in plot units, or None when there is no target.
    offset: How far past the detector to put the label, in plot units.
    fallback_dx, fallback_dy: Offset used when there is no target to aim away from.

  Returns:
    (x, y, ha, va) for the text call. With a target, the name box is centered on
    the point, so its middle sits on the target-to-detector line.
  """
  if tgt_x is None or tgt_y is None:
    return det_x + fallback_dx, det_y + fallback_dy, "left", "bottom"

  vx = det_x - tgt_x
  vy = det_y - tgt_y
  length = float(np.hypot(vx, vy))
  if length == 0.0:
    return det_x + fallback_dx, det_y + fallback_dy, "left", "bottom"

  ux = vx / length
  uy = vy / length
  # Centered on both axes so the middle of the name box lands on the line itself,
  # not one of its corners.
  return det_x + ux * offset, det_y + uy * offset, "center", "center"


def add_detectors(ax, detectors: List[Dict[str, Any]],
                  target_xy: Tuple[float, float] | None, extent, show_distance: bool,
                  scale: float = 1.0, label_offset: float | None = None):
  xspan = extent[1] - extent[0]
  yspan = extent[3] - extent[2]
  dx = max(xspan * 0.01, 2.0 / scale)
  dy = max(yspan * 0.01, 2.0 / scale)
  # label_offset arrives in meters; default to 3% of the plotted width.
  offset = (label_offset / scale) if label_offset is not None else xspan * 0.03
  # Pale outlines, same as the range circles, so the black marks stay visible
  # where they cross dark terrain.
  line_halo = [path_effects.withStroke(linewidth=2.2, foreground="white", alpha=0.8)]
  marker_halo = [path_effects.withStroke(linewidth=2.5, foreground="white", alpha=0.9)]

  for det in detectors:
    det_x = det["x"] / scale
    det_y = det["y"] / scale
    if target_xy:
      tgt_x = target_xy[0] / scale
      tgt_y = target_xy[1] / scale
      # Opaque on purpose: a see-through line lets the white outline underneath
      # show through and the line reads as white instead of black.
      ax.plot([det_x, tgt_x], [det_y, tgt_y],
              color="black", linewidth=1.0, zorder=4,
              path_effects=line_halo)
    ax.scatter(det_x, det_y, color="black", s=35, zorder=5,
               path_effects=marker_halo)
    label_x, label_y, ha, va = name_label_position(
        det_x, det_y,
        target_xy[0] / scale if target_xy else None,
        target_xy[1] / scale if target_xy else None,
        offset, dx, dy)
    ax.text(label_x, label_y, det["name"],
            color="black", fontsize=8, ha=ha, va=va,
            zorder=6, bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.7))
    if show_distance and target_xy:
      tgt_x = target_xy[0] / scale
      tgt_y = target_xy[1] / scale
      midx = (det_x + tgt_x) * 0.5
      midy = (det_y + tgt_y) * 0.5
      # Distance is always in meters
      dist = float(np.hypot(det["x"] - target_xy[0], det["y"] - target_xy[1]))
      ax.text(midx, midy, f"{dist:.0f}", color="black", fontsize=7,
              ha="center", va="center", zorder=6,
              bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.6))

  if target_xy:
    tgt_x = target_xy[0] / scale
    tgt_y = target_xy[1] / scale
    ax.scatter([tgt_x], [tgt_y], marker="*", s=160,
               color="black", zorder=6, path_effects=marker_halo)


def main():
  parser = argparse.ArgumentParser(
      description="Plot a DEM together with detector positions (from CSV or JSON5 config).")
  parser.add_argument("--config", default=None,
                      help="Path to JSON5 config file. Extracts DEM path, detectors, "
                           "and target coords. Overridable by --dem, --csv, --xcnt, --ycnt.")
  parser.add_argument("--csv", default=None,
                      help="Path to CSV file with detector coordinates (from kml_to_csv.py)")
  parser.add_argument("--dem", default=None,
                      help="Path to DEM file (.g2zbin or .tif/.tiff GeoTIFF)")
  parser.add_argument("--output", default=None,
                      help="Optional PNG path. Defaults to fig_<csv_name>.png")
  parser.add_argument("--cmap", default="terrain",
                      help="Matplotlib colormap name (default: terrain)")
  parser.add_argument("--zmin", type=float, default=None, help="Optional color scale min")
  parser.add_argument("--zmax", type=float, default=None, help="Optional color scale max")
  parser.add_argument("--contour", action="store_true", help="Overlay contour lines")
  parser.add_argument("--cint", type=float, default=10.0, help="Contour interval (default: 10)")
  parser.add_argument("--shade", action="store_true",
                      help="Overlay grayscale hillshade on the elevation map")
  parser.add_argument("--shade-alpha", type=float, default=0.2,
                      help="Hillshade overlay opacity (default: 0.2)")
  parser.add_argument("--shade-azdeg", type=float, default=315.0,
                      help="Hillshade light azimuth in degrees (default: 315)")
  parser.add_argument("--shade-altdeg", type=float, default=45.0,
                      help="Hillshade light altitude in degrees (default: 45)")
  parser.add_argument("--shade-vert-exag", type=float, default=2.0,
                      help="Hillshade vertical exaggeration (default: 2.0)")
  parser.add_argument("--dpi", type=int, default=300, help="Output PNG DPI (default: 300)")
  parser.add_argument("--show_dist", "--show_distance", dest="show_dist", action="store_true",
                      help="Show detector-target distance labels (always in meters)")
  parser.add_argument("--xcnt", type=float, default=None,
                      help="Override target X coordinate from CSV (in meters)")
  parser.add_argument("--ycnt", type=float, default=None,
                      help="Override target Y coordinate from CSV (in meters)")
  parser.add_argument("--title", default=None,
                      help="Plot title (default: DEM filename with detectors)")
  parser.add_argument("--xyunit", choices=["m", "km"], default="m",
                      help="X/Y axis unit: 'm' (meters) or 'km' (kilometers). Default: m")
  parser.add_argument("--epsg", type=int, default=None,
                      help="EPSG code to display in the plot title (e.g. 6677)")
  parser.add_argument("--half-range", "--half_range", dest="half_range",
                      type=float, default=None,
                      help="Crop the DEM to the target point +/- this many meters")
  parser.add_argument("--xmin", type=float, default=None,
                      help="Left edge of the drawn area, in meters. Cuts further than "
                           "--half-range; the plot center moves accordingly.")
  parser.add_argument("--xmax", type=float, default=None,
                      help="Right edge of the drawn area, in meters")
  parser.add_argument("--ymin", type=float, default=None,
                      help="Bottom edge of the drawn area, in meters")
  parser.add_argument("--ymax", type=float, default=None,
                      help="Top edge of the drawn area, in meters")
  parser.add_argument("--circles", default=None,
                      help="Comma-separated circle radii in meters, drawn around the "
                           "target point (e.g. 500,1000,1500)")
  parser.add_argument("--circle-color", "--circle_color", dest="circle_color",
                      default="black",
                      help="Color of the range circles and their labels (default: black)")
  parser.add_argument("--label-offset", "--label_offset", dest="label_offset",
                      type=float, default=None,
                      help="How far past each detector to put its name, in meters, "
                           "measured along the target-to-detector line "
                           "(default: 3%% of the plotted width)")

  args = parser.parse_args()

  # --- Resolve inputs: --config OR (--csv AND --dem) ---
  if args.config:
    config_path = Path(args.config)
    if not config_path.exists():
      raise FileNotFoundError(f"Config file not found: {config_path}")
    config_dir = config_path.parent

    prm = load_json5(config_path)

    # DEM: from config, overridable by --dem
    if args.dem:
      dem_path = Path(args.dem)
    else:
      grid2d = prm.get("GRID2D_PILLAR_PARAMETERS")
      if not grid2d:
        raise KeyError("GRID2D_PILLAR_PARAMETERS not found in config.")
      path_dem = grid2d.get("path_dem")
      if not path_dem:
        raise KeyError("GRID2D_PILLAR_PARAMETERS.path_dem not found in config.")
      dem_path = resolve_path(config_dir, path_dem)

    # Detectors: from --csv if given, otherwise from config det_files
    if args.csv:
      csv_path = Path(args.csv)
      if not csv_path.exists():
        raise FileNotFoundError(f"CSV file not found: {csv_path}")
      detectors, target_xy_from_source = load_detectors_from_csv(csv_path)
    else:
      detectors_section = prm.get("DETECTOR_PARAMETER_LISTS")
      if not detectors_section or not detectors_section.get("det_files"):
        parser.error("No detector info: --csv not given and DETECTOR_PARAMETER_LISTS.det_files "
                     "not found in config. Provide --csv or fix the config.")
      det_files = detectors_section["det_files"]
      detectors = load_detectors_from_config(det_files, config_dir)
      target_xy_from_source = None

    # Target: --xcnt/--ycnt > CSV/source > config DEPTH_RESOLUTION_SWEEP
    if args.xcnt is not None and args.ycnt is not None:
      target_xy = (args.xcnt, args.ycnt)
    elif target_xy_from_source is not None:
      target_xy = target_xy_from_source
    else:
      sweep = prm.get("DEPTH_RESOLUTION_SWEEP", {})
      sweep_common = sweep.get("common", {})
      if "x_cnt_obj" in sweep_common and "y_cnt_obj" in sweep_common:
        target_xy = (float(sweep_common["x_cnt_obj"]), float(sweep_common["y_cnt_obj"]))
      else:
        target_xy = None

  elif args.csv and args.dem:
    csv_path = Path(args.csv)
    dem_path = Path(args.dem)

    if not csv_path.exists():
      raise FileNotFoundError(f"CSV file not found: {csv_path}")

    detectors, target_xy_from_csv = load_detectors_from_csv(csv_path)

    if args.xcnt is not None and args.ycnt is not None:
      target_xy = (args.xcnt, args.ycnt)
    else:
      target_xy = target_xy_from_csv

  else:
    parser.error("Either --config or both --csv and --dem are required.")

  if not dem_path.exists():
    raise FileNotFoundError(f"DEM file not found: {dem_path}")

  if not detectors:
    print("[WARN] No detectors found")

  if args.show_dist and target_xy is None:
    print("[WARN] --show_distance requires a target point, but none was found. "
          "Use --xcnt/--ycnt to specify, or set DEPTH_RESOLUTION_SWEEP in config, "
          "or mark a row as 'center'/'-----' in CSV.")

  # Load DEM
  suffix = dem_path.suffix.lower()
  if suffix in (".tif", ".tiff"):
    xi, yi, Z = read_geotiff(str(dem_path))
  elif suffix == ".g2zbin":
    xi, yi, Z = load_g2zbin_as_grid(str(dem_path))
    # Legacy g2zbin files mark NoData as -9999; mask it like the GeoTIFF branch
    Z = np.where(Z <= -9998.0, np.nan, Z)
  else:
    raise ValueError(
        f"Unsupported DEM format: {suffix}. Only .g2zbin and .tif/.tiff are supported.")

  # Check target and detectors against DEM extent
  x_lo, x_hi = float(xi.min()), float(xi.max())
  y_lo, y_hi = float(yi.min()), float(yi.max())
  has_extent_error = False

  if target_xy is not None:
    tx, ty = target_xy
    if not (x_lo <= tx <= x_hi) or not (y_lo <= ty <= y_hi):
      print(f"[ERROR] Target point ({tx:.1f}, {ty:.1f}) is outside the DEM extent "
            f"(x: {x_lo:.1f}..{x_hi:.1f}, y: {y_lo:.1f}..{y_hi:.1f}). "
            "Please check --xcnt/--ycnt or the config values.")
      has_extent_error = True

  for det in detectors:
    dx, dy = det["x"], det["y"]
    if not (x_lo <= dx <= x_hi) or not (y_lo <= dy <= y_hi):
      print(f"[ERROR] Detector '{det['name']}' ({dx:.1f}, {dy:.1f}) is outside the DEM extent "
            f"(x: {x_lo:.1f}..{x_hi:.1f}, y: {y_lo:.1f}..{y_hi:.1f}). "
            "Please check the detector coordinates.")
      has_extent_error = True

  if has_extent_error:
    sys.exit(1)

  circle_radii: List[float] = []
  if args.circles:
    circle_radii = [float(tok) for tok in args.circles.split(",") if tok.strip()]

  if (args.half_range is not None or circle_radii) and target_xy is None:
    parser.error("--half-range and --circles need a target point. "
                 "Give --xcnt/--ycnt, or a config/CSV that provides one.")

  # Window: start from the whole DEM, narrow it by --half-range around the target,
  # then let any explicit edge cut it further.
  win = [x_lo, x_hi, y_lo, y_hi]
  if args.half_range is not None:
    win = [target_xy[0] - args.half_range, target_xy[0] + args.half_range,
           target_xy[1] - args.half_range, target_xy[1] + args.half_range]
  for i, edge in enumerate((args.xmin, args.xmax, args.ymin, args.ymax)):
    if edge is not None:
      win[i] = edge
  if win != [x_lo, x_hi, y_lo, y_hi]:
    xi, yi, Z = crop_to_window(xi, yi, Z, *win)

  # Plot
  fig, ax = plt.subplots(figsize=(8, 6))
  extent, _, scale = plot_background(ax, xi, yi, Z, args.cmap, args.zmin, args.zmax,
                                     args.contour, args.cint, args.xyunit,
                                     args.shade, args.shade_alpha,
                                     args.shade_azdeg, args.shade_altdeg,
                                     args.shade_vert_exag)
  add_detectors(ax, detectors, target_xy, extent, args.show_dist, scale, args.label_offset)
  if circle_radii:
    add_range_circles(ax, target_xy, circle_radii, scale, args.circle_color)
  # Keep the view on the DEM window even when detectors or circles reach outside it
  ax.set_xlim(extent[0], extent[1])
  ax.set_ylim(extent[2], extent[3])

  title = args.title or f"{dem_path.name} with detectors"
  if args.epsg is not None:
    title += f" (EPSG:{args.epsg})"
  ax.set_title(title)
  ax.grid(color="grey", linestyle=":", linewidth=1.0, alpha=0.6)
  plt.tight_layout()

  # Output path
  if args.output:
    output_path = Path(args.output)
  elif args.csv:
    output_path = Path(args.csv).with_name(f"fig_{Path(args.csv).stem}.png")
  else:
    output_path = dem_path.with_name(f"{dem_path.stem}_detectors.png")

  fig.savefig(output_path, dpi=args.dpi)
  print(f"Saved plot to {output_path}")


if __name__ == "__main__":
  main()
