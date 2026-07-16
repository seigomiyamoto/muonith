#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot detection limit heatmaps for muography.

Creates separate figures for each delta_density showing log10(p-value)
as a function of anomaly size and depth.

Usage:
  python plot_detection_limit_heatmap.py --main_config prm_reso.json5 --plot_config heatmap_config.json5 --outdir figs
"""

import argparse
import sys
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
try:
  import pandas as pd
except ImportError:
  sys.exit(
    "ERROR: pandas is required.\n"
    "  pip install pandas  (or activate venv: source .venv/bin/activate)"
  )
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
try:
  from scipy.interpolate import griddata
except ImportError:
  sys.exit(
    "ERROR: scipy is required for plot_detection_limit_heatmap.py\n"
    "  pip install scipy  (or activate venv: source .venv/bin/activate)"
  )

# Default configuration
DEFAULT_CONFIG: dict[str, Any] = {
  # Colormap settings
  "colormap": "RdYlBu_r",
  "vmin": -3.0,  # log10(p) = -3 => p = 0.001
  "vmax": 0.0,   # log10(p) = 0  => p = 1
  "ngrad": None,  # Number of gradient levels (null = continuous)

  # Colorbar extend: "neither", "min", "max", "both"
  "extend": "neither",
  # Invert colorbar axis (swap top/bottom)
  "colorbar_invert": False,
  # Colors for out-of-range values (null = use colormap extremes)
  "over_color": None,
  "under_color": None,

  # NaN background color (areas with no data)
  "nan_color": "lightgray",

  # Number of grid cells to expand NaN mask for contour suppression (0 = disabled)
  "contour_nan_margin": 0,

  # Contour lines (p-value thresholds)
  "contours": [
    {"p_value": 0.3, "color": "white", "linestyle": "-", "linewidth": 1.5, "label": True},
    {"p_value": 0.1, "color": "black", "linestyle": "--", "linewidth": 1.0, "label": True},
  ],

  # Legend location: "auto" (determined by depth_axis_mode), or explicit position
  #   e.g., "upper right", "lower right", "upper left", "lower left", "best", etc.
  "legend_loc": "auto",
  "legend_fontsize": 8,

  # Grid settings
  "grid": {
    "enabled": True,
    "color": "gray",
    "alpha": 0.5,
    "linestyle": ":",
    "linewidth": 0.5,
  },

  # Draw diagonal line (depth = size) using grid style
  "show_diagonal": False,

  # Axis ranges (null = auto)
  "size_max": None,
  "depth_max": None,

  # Surface elevation [m] for secondary Y-axis.
  # "auto" = read from DEM at (x_cnt_obj, y_cnt_obj), null = disabled, float = manual
  "surface_elevation": None,
  # Label for the secondary elevation axis
  "elevation_label": "Elevation a.s.l. [m]",

  # Depth axis mode: "bottom" uses depth_anom_btm, "top" uses depth_anom_top,
  #   "detector" uses depth_anom_top with size_bars auto-aligned at max depth_anom_btm per detector
  "depth_axis_mode": "bottom",

  # Size bars: T-capped vertical bars showing anomaly extent
  "size_bars": {
    "enabled": False,
    "sizes": [],
    "align": "top",
    "align_depth": 0,
    "color": "magenta",
    "linewidth": 2.0,
    "cap_width": 10,
  },

  # Resolution (number of bins, null = auto from data step)
  "n_size_bins": None,
  "n_depth_bins": None,

  # Figure settings
  "dpi": 150,
  "figsize_per_panel": [5, 4],  # [width, height] per panel

  # Font sizes
  "title_fontsize": 11,
  "label_fontsize": 10,  # x, y, colorbar labels
  "tick_fontsize": 9,    # x, y, colorbar ticks

  # Parallel processing settings (for PNG generation)
  "parallel": {
    "enabled": False,
    "max_workers": 4,
  },

  # Animation GIF settings
  "animation": {
    "tf_exec": False,
    "sec_per_frame": 0.5,
    "loop": 0,
    "tool": "convert",
    "output_subdir": "gif",  # Subdirectory for GIF output (relative to PNG dir)
  },
}


def load_config(config_path: Path | None) -> dict[str, Any]:
  """Load configuration from JSON5 file, merged with defaults."""
  config = DEFAULT_CONFIG.copy()

  if config_path is not None and config_path.exists():
    import json5
    with open(config_path, "r") as f:
      user_config = json5.load(f)

    # Deep merge for nested dicts
    for key, value in user_config.items():
      if key in config and isinstance(config[key], dict) and isinstance(value, dict):
        config[key] = {**config[key], **value}
      else:
        config[key] = value

  return config


def get_depth_column(config: dict[str, Any]) -> str:
  """Return CSV column name based on depth_axis_mode config."""
  mode = config.get("depth_axis_mode", "bottom")
  if mode in ("top", "detector"):
    return "depth_anom_top"
  return "depth_anom_btm"


def get_depth_axis_mode(config: dict[str, Any]) -> str:
  """Return depth_axis_mode from config."""
  return config.get("depth_axis_mode", "bottom")


def get_surface_elevation_from_dem(main_config_path: Path) -> float | None:
  """Read surface elevation at target center from DEM binary file.

  Uses x_cnt_obj/y_cnt_obj from DEPTH_RESOLUTION_SWEEP.common and
  path_dem from GRID2D_PILLAR_PARAMETERS to find the nearest DEM point.

  Returns:
    Surface elevation [m], or None if data is unavailable.
  """
  import json5
  with open(main_config_path, "r") as f:
    cfg = json5.load(f)

  sweep = cfg.get("DEPTH_RESOLUTION_SWEEP", {}).get("common", {})
  x_cnt = sweep.get("x_cnt_obj")
  y_cnt = sweep.get("y_cnt_obj")
  if x_cnt is None or y_cnt is None:
    return None

  dem_rel = cfg.get("GRID2D_PILLAR_PARAMETERS", {}).get("path_dem")
  if dem_rel is None:
    return None

  dem_path = main_config_path.parent / dem_rel
  if not dem_path.exists():
    return None

  # DEM must be in g2zbin format (magic bytes "G2ZBIN\0\0")
  with open(dem_path, "rb") as f:
    head = f.read(8)

  if head[:6] != b"G2ZBIN":
    return None

  from g2zbin_io import axis_centers, read_g2zbin
  info, Z = read_g2zbin(dem_path)
  # Bin centers (v2 min/max are bin edges)
  xi = axis_centers(info["x_axis"])
  yi = axis_centers(info["y_axis"])
  # Find nearest grid point
  xg, yg = np.meshgrid(xi, yi)
  dist_sq = (xg - x_cnt)**2 + (yg - y_cnt)**2
  idx = np.unravel_index(np.argmin(dist_sq), dist_sq.shape)
  return float(Z[idx])


def get_detector_name(main_config_path: Path, detid: str) -> str | None:
  """Get detector name from main config and detector parameter file.

  Args:
    main_config_path: Path to the main JSON5 config (e.g., prm_reso.json5)
    detid: Detector ID string (e.g., "00", "01")

  Returns:
    Detector name from DETECTOR_PARAMETERS.name, or None if not found
  """
  if main_config_path is None or not main_config_path.exists():
    return None

  import json5
  with open(main_config_path, "r") as f:
    main_config = json5.load(f)

  det_files = main_config.get("DETECTOR_PARAMETER_LISTS", {}).get("det_files", [])
  det_idx = int(detid)
  if det_idx >= len(det_files):
    return None

  det_file_path = main_config_path.parent / det_files[det_idx]
  if not det_file_path.exists():
    return None

  with open(det_file_path, "r") as f:
    det_config = json5.load(f)

  return det_config.get("DETECTOR_PARAMETERS", {}).get("name")


def load_and_filter_data(csv_path: Path, signal_amp: float | None = None) -> pd.DataFrame:
  """Load CSV and optionally filter by signal_amp."""
  df = pd.read_csv(csv_path)

  required_cols = ["obsize", "ddens", "depth_anom_top", "depth_anom_btm", "p_val", "signal_amp"]
  missing = [c for c in required_cols if c not in df.columns]
  if missing:
    raise ValueError(f"Missing required columns: {missing}")

  if signal_amp is not None:
    df = df[np.isclose(df["signal_amp"], signal_amp)].copy()
    if df.empty:
      raise ValueError(f"No data found for signal_amp={signal_amp}")

  return df


def create_heatmap_data(df: pd.DataFrame, ddens_val: float,
                        size_bins: np.ndarray, depth_bins: np.ndarray,
                        depth_col: str = "depth_anom_btm") -> np.ndarray:
  """Create 2D grid of log10(p-value) for a specific ddens value."""
  sub = df[df["ddens"] == ddens_val].copy()

  if sub.empty:
    return np.full((len(depth_bins) - 1, len(size_bins) - 1), np.nan)

  sizes = sub["obsize"].values
  depths = sub[depth_col].values
  pvals = sub["p_val"].values

  # Compute log10(p-value), handling p=0 or p=1 edge cases
  with np.errstate(divide='ignore', invalid='ignore'):
    log_pvals = np.log10(np.clip(pvals, 1e-10, 1.0))

  # Create grid using interpolation
  size_centers = 0.5 * (size_bins[:-1] + size_bins[1:])
  depth_centers = 0.5 * (depth_bins[:-1] + depth_bins[1:])
  grid_size, grid_depth = np.meshgrid(size_centers, depth_centers)

  if len(sizes) >= 4:
    # Use linear interpolation for smooth heatmap
    grid_logp = griddata(
      (sizes, depths), log_pvals,
      (grid_size, grid_depth),
      method='linear'
    )
  else:
    # Fall back to nearest neighbor for sparse data
    grid_logp = griddata(
      (sizes, depths), log_pvals,
      (grid_size, grid_depth),
      method='nearest'
    )

  return grid_logp


def infer_step(values: np.ndarray) -> float:
  """Infer step size from unique sorted values."""
  unique = np.sort(np.unique(values))
  if len(unique) < 2:
    return 1.0
  diffs = np.diff(unique)
  # Use median to be robust against outliers
  return float(np.median(diffs[diffs > 0]))


def plot_single_heatmap(df: pd.DataFrame, ddens_val: float, output_path: Path,
                        config: dict[str, Any]) -> None:
  """Create a single heatmap figure for one ddens value."""
  depth_col = get_depth_column(config)

  # Determine axis ranges
  size_max = config["size_max"]
  depth_max = config["depth_max"]
  if size_max is None:
    size_max = df["obsize"].max() * 1.05
  if depth_max is None:
    depth_max = df[depth_col].max() * 1.05

  # Create bins for heatmap (auto-detect from data if null)
  n_size_bins = config.get("n_size_bins")
  n_depth_bins = config.get("n_depth_bins")

  if n_size_bins is None:
    size_step = infer_step(df["obsize"].values)
    n_size_bins = max(1, int(np.ceil(size_max / size_step)))
  if n_depth_bins is None:
    depth_step = infer_step(df[depth_col].values)
    n_depth_bins = max(1, int(np.ceil(depth_max / depth_step)))

  size_bins = np.linspace(0, size_max, n_size_bins + 1)
  depth_bins = np.linspace(0, depth_max, n_depth_bins + 1)

  # Color normalization
  vmin = config["vmin"]
  vmax = config["vmax"]
  ngrad = config.get("ngrad")

  if ngrad is not None and ngrad > 0:
    # Discrete colormap with ngrad levels
    cmap = plt.get_cmap(config["colormap"], ngrad).copy()
    norm = Normalize(vmin=vmin, vmax=vmax)
  else:
    # Continuous colormap
    cmap = plt.get_cmap(config["colormap"]).copy()
    norm = Normalize(vmin=vmin, vmax=vmax)
  cmap.set_bad(color=config["nan_color"])
  if config.get("over_color"):
    cmap.set_over(color=config["over_color"])
  if config.get("under_color"):
    cmap.set_under(color=config["under_color"])

  # Create figure
  figsize_per_panel = config["figsize_per_panel"]
  fig, ax = plt.subplots(figsize=(figsize_per_panel[0], figsize_per_panel[1]))

  # Set background color for NaN regions
  ax.set_facecolor(config["nan_color"])

  # Create heatmap data
  grid_logp = create_heatmap_data(df, ddens_val, size_bins, depth_bins, depth_col)

  # Plot heatmap
  size_centers = 0.5 * (size_bins[:-1] + size_bins[1:])
  depth_centers = 0.5 * (depth_bins[:-1] + depth_bins[1:])

  im = ax.pcolormesh(size_bins, depth_bins, grid_logp,
                     cmap=cmap, norm=norm, shading='flat')

  # Expand NaN mask to suppress noisy contour lines near data boundary
  contour_nan_margin = config.get("contour_nan_margin", 0)
  if contour_nan_margin > 0:
    from scipy.ndimage import binary_dilation
    nan_mask = np.isnan(grid_logp)
    expanded_nan = binary_dilation(nan_mask, iterations=contour_nan_margin)
    grid_logp_contour = grid_logp.copy()
    grid_logp_contour[expanded_nan] = np.nan
  else:
    grid_logp_contour = grid_logp

  # Add contour lines
  contours = config["contours"]
  legend_handles = []
  legend_labels = []
  if not np.all(np.isnan(grid_logp_contour)):
    for cont in contours:
      p_val = cont["p_value"]
      log_p = np.log10(p_val)
      color = cont.get("color", "white")
      linestyle = cont.get("linestyle", "-")
      linewidth = cont.get("linewidth", 1.5)
      try:
        cs = ax.contour(
          size_centers, depth_centers, grid_logp_contour,
          levels=[log_p],
          colors=[color],
          linestyles=[linestyle],
          linewidths=[linewidth]
        )
        if cont.get("label", True) and len(cs.allsegs[0]) > 0:
          # Add to legend instead of inline label
          handle = Line2D([0], [0], color=color, linestyle=linestyle,
                          linewidth=linewidth)
          legend_handles.append(handle)
          legend_labels.append(f'p={p_val}')
      except (ValueError, IndexError):
        pass  # No contour if threshold not in range

  # Grid
  grid_cfg = config["grid"]
  if grid_cfg["enabled"]:
    ax.grid(
      True,
      color=grid_cfg["color"],
      alpha=grid_cfg["alpha"],
      linestyle=grid_cfg["linestyle"],
      linewidth=grid_cfg["linewidth"]
    )

  # Diagonal line (depth = size), only meaningful in "bottom" mode
  mode = get_depth_axis_mode(config)
  if config.get("show_diagonal", False) and mode == "bottom":
    diag_max = min(size_max, depth_max)
    ax.plot([0, diag_max], [0, diag_max],
            color=grid_cfg["color"], alpha=grid_cfg["alpha"],
            linestyle=grid_cfg["linestyle"], linewidth=grid_cfg["linewidth"])

  # Size bars (T-capped vertical bars showing anomaly extent)
  sb_cfg = config.get("size_bars", {})
  if sb_cfg.get("enabled", False):
    mode = get_depth_axis_mode(config)
    if mode == "detector":
      # Auto-determine: align bottom at max depth_anom_btm in data
      align = "bottom"
      align_depth = float(df["depth_anom_btm"].max())
    else:
      align = sb_cfg.get("align", "top")
      align_depth = sb_cfg.get("align_depth", 0)
    sb_color = sb_cfg.get("color", "magenta")
    sb_lw = sb_cfg.get("linewidth", 2.0)
    # T-cap half-width: prefer a fraction of the x-axis range (size_max) so the
    # cap keeps a consistent on-screen size at any zoom; fall back to an absolute
    # width [m] for backward compatibility.
    cap_frac = sb_cfg.get("cap_width_frac")
    if cap_frac is not None:
      half_cap = size_max * cap_frac / 2
    else:
      half_cap = sb_cfg.get("cap_width", 10) / 2
    sizes = sb_cfg.get("sizes", [])
    if sizes == "auto":
      # Follow the x-axis major ticks; keep only interior ticks (exclude 0 and the size_max edge)
      ax.set_xlim(0, size_max)
      sizes = [float(t) for t in ax.get_xticks() if 0 < t < size_max]
    for L in sizes:
      if L <= 0 or L > size_max:
        continue
      if align == "top":
        y_top, y_bottom = align_depth, align_depth + L
      else:  # bottom
        y_top, y_bottom = align_depth - L, align_depth
      # Keep the whole bar inside [0, depth_max]: if the bottom would fall below
      # the axis, shift the bar up by the overflow (preserving the bar length L).
      y_lo, y_hi = min(y_top, y_bottom), max(y_top, y_bottom)
      if y_hi > depth_max:
        y_lo -= (y_hi - depth_max)
        y_hi = depth_max
      y_lo = max(y_lo, 0.0)
      y_top, y_bottom = y_lo, y_hi
      ax.plot([L, L], [y_top, y_bottom], color=sb_color, linewidth=sb_lw)
      ax.plot([L - half_cap, L + half_cap], [y_top, y_top], color=sb_color, linewidth=sb_lw)
      ax.plot([L - half_cap, L + half_cap], [y_bottom, y_bottom], color=sb_color, linewidth=sb_lw)
    legend_handles.append(Line2D([0], [0], color=sb_color, linewidth=sb_lw))
    legend_labels.append("Anomaly\nsize indicator")

  # Combined legend (contour lines + size bars)
  if legend_handles:
    legend_loc = config.get("legend_loc", "auto")
    if legend_loc == "auto":
      # top mode (depth_anom_top): data in upper-left → legend lower-right
      # bottom mode (depth_anom_btm): data in lower-left → legend upper-right
      legend_loc = "lower right" if depth_col == "depth_anom_top" else "upper right"
    legend_fontsize = config.get("legend_fontsize", 8)
    loc_aliases = {
      "top-left": "upper left", "top-right": "upper right",
      "bottom-left": "lower left", "bottom-right": "lower right",
      "top": "upper center", "bottom": "lower center",
    }
    legend_loc = loc_aliases.get(legend_loc, legend_loc)
    leg = ax.legend(legend_handles, legend_labels, loc=legend_loc, fontsize=legend_fontsize)
    leg.get_frame().set_facecolor(config["nan_color"])
    leg.get_frame().set_edgecolor('none')

  # Font sizes
  title_fontsize = config.get("title_fontsize", 11)
  label_fontsize = config.get("label_fontsize", 10)
  tick_fontsize = config.get("tick_fontsize", 9)

  # Labels
  samp = config.get("signal_amp")
  detector_name = config.get("detector_name")
  detid = config.get("detid")
  # Use detector_name if available, otherwise fall back to detid
  if detector_name is not None:
    detid_prefix = f"detid={detector_name}, "
  elif detid is not None:
    detid_prefix = f"detid={detid}, "
  else:
    detid_prefix = ""
  if samp is not None:
    ax.set_title(f"{detid_prefix}$\\Delta\\rho$ = {int(ddens_val)} kg/m$^3$, signal_amp={samp:.1E}",
                 fontsize=title_fontsize)
  else:
    ax.set_title(f"{detid_prefix}$\\Delta\\rho$ = {int(ddens_val)} kg/m$^3$", fontsize=title_fontsize)
  ax.invert_yaxis()  # Depth increases downward
  ax.set_xlim(0, size_max)
  ax.set_ylim(depth_max, 0)
  depth_label = "Depth from surface (anomaly top) [m]" if depth_col == "depth_anom_top" else "Depth from surface (anomaly bottom) [m]"
  ax.set_ylabel(depth_label, fontsize=label_fontsize)
  ax.set_xlabel("Anomaly size L [m]", fontsize=label_fontsize)
  ax.tick_params(axis='both', labelsize=tick_fontsize)

  # Secondary Y-axis for elevation (when surface_elevation is configured)
  surface_elevation = config.get("surface_elevation")
  if surface_elevation is not None:
    elev_forward = lambda d, se=surface_elevation: se - d
    elev_inverse = lambda e, se=surface_elevation: se - e
    sec_ax = ax.secondary_yaxis('right', functions=(elev_forward, elev_inverse))
    elevation_label = config.get("elevation_label", "Elevation a.s.l. [m]")
    sec_ax.set_ylabel(elevation_label, fontsize=label_fontsize)
    sec_ax.tick_params(labelsize=tick_fontsize)

  # Add colorbar
  extend = config.get("extend", "neither")
  colorbar_pad = 0.14 if surface_elevation is not None else 0.04
  cbar = fig.colorbar(im, ax=ax, orientation='vertical', fraction=0.046, pad=colorbar_pad,
                      extend=extend)
  side_label = "two-sided" if config.get("both_side", False) else "one-sided"
  cbar.set_label(f"p-value ({side_label})", fontsize=label_fontsize)
  # Format ticks as 10^x
  from matplotlib.ticker import FuncFormatter
  def exp_formatter(x, pos):
    if x == int(x):
      return f"$10^{{{int(x)}}}$"
    else:
      return f"$10^{{{x:.1f}}}$"
  cbar.ax.yaxis.set_major_formatter(FuncFormatter(exp_formatter))
  cbar.ax.tick_params(labelsize=tick_fontsize)
  if config.get("colorbar_invert", False):
    cbar.ax.invert_yaxis()

  fig.tight_layout()
  fig.savefig(output_path, dpi=config["dpi"], bbox_inches='tight')
  plt.close(fig)
  print(f"Saved: {output_path}")


def create_animation_gif(png_files: list[Path], output_path: Path,
                         config: dict[str, Any]) -> None:
  """Create animated GIF from PNG files using ImageMagick convert or gifski.

  Args:
    png_files: List of PNG file paths in desired frame order
    output_path: Output GIF file path
    config: Configuration dict containing animation settings
  """
  import subprocess

  if not png_files:
    print("No PNG files provided for animation.")
    return

  # Check all PNG files exist
  existing_pngs = [p for p in png_files if p.exists()]
  if len(existing_pngs) != len(png_files):
    missing = [p for p in png_files if not p.exists()]
    print(f"Warning: {len(missing)} PNG files not found, skipping animation.")
    return

  anim_cfg = config.get("animation", {})
  tool = anim_cfg.get("tool", "convert")
  sec_per_frame = anim_cfg.get("sec_per_frame", 0.5)
  loop = anim_cfg.get("loop", 0)
  dpi = anim_cfg.get("dpi", 150)

  # Convert sec_per_frame to delay (centiseconds for ImageMagick)
  delay = max(1, int(sec_per_frame * 100))

  png_str_files = [str(p) for p in existing_pngs]

  if tool == "gifski":
    # gifski accepts fractional fps (e.g., 0.5 for 2 sec/frame)
    fps = 1.0 / sec_per_frame if sec_per_frame > 0 else 2
    cmd = ["gifski", "--fps", str(fps), "-o", str(output_path), "--"] + png_str_files
  else:  # default: ImageMagick convert
    cmd = ["convert", "-density", str(dpi), "-delay", str(delay),
           "-loop", str(loop)] + png_str_files + [str(output_path)]

  print(f"Creating animation: {output_path.name}")
  try:
    subprocess.run(cmd, check=True, capture_output=True)
    print(f"Saved animation: {output_path}")
  except subprocess.CalledProcessError as e:
    print(f"Error creating animation: {e}")
    if e.stderr:
      print(f"  stderr: {e.stderr.decode()}")
  except FileNotFoundError:
    print(f"Error: {tool} not found. Please install ImageMagick or gifski.")


def create_all_animation_gifs(gif_tasks: list[tuple[list[Path], Path]],
                              config: dict[str, Any]) -> None:
  """Create animation GIFs sequentially.

  Note: gifski uses internal multi-threading, so external parallelization
  is not needed and may cause resource contention.

  Args:
    gif_tasks: List of (png_files, output_path) tuples
    config: Configuration dict containing animation settings
  """
  if not gif_tasks:
    return

  # Create output directories if they don't exist
  output_dirs = set(out.parent for _, out in gif_tasks)
  for out_dir in output_dirs:
    out_dir.mkdir(parents=True, exist_ok=True)

  print(f"Creating {len(gif_tasks)} animations...")
  for pngs, out in gif_tasks:
    create_animation_gif(pngs, out, config)


def plot_heatmaps(df: pd.DataFrame, output_base: Path, config: dict[str, Any]) -> list[Path]:
  """Create separate heatmap figures for each ddens value.

  Args:
    df: DataFrame containing heatmap data
    output_base: Base path for output files
    config: Configuration dict (may contain "vec_delta_density" for ordering)

  Returns:
    List of generated PNG file paths in ddens order
  """
  # Use vec_delta_density order if provided, otherwise use data order (no sorting)
  vec_ddens = config.get("vec_delta_density")
  data_ddens = set(df["ddens"].unique())

  if vec_ddens is not None:
    # Filter to only ddens values that exist in data, keep original order
    ddens_values = [d for d in vec_ddens if d in data_ddens]
  else:
    # Use data order without sorting
    ddens_values = list(df["ddens"].unique())

  generated_pngs: list[Path] = []

  if len(ddens_values) == 0:
    print("No ddens values found in data.")
    return generated_pngs

  for ddens_val in ddens_values:
    # Generate output filename with ddens
    stem = output_base.stem
    suffix = output_base.suffix
    ddens_str = f"{ddens_val:+06.0f}".replace("+", "p").replace("-", "m")
    output_path = output_base.parent / f"{stem}_dd{ddens_str}{suffix}"

    plot_single_heatmap(df, ddens_val, output_path, config)
    generated_pngs.append(output_path)

  return generated_pngs


def load_main_config(main_config_path: Path) -> dict[str, Any]:
  """Load main configuration from JSON5 file."""
  import json5
  with open(main_config_path, "r") as f:
    return json5.load(f)


def process_png_task(task: dict[str, Any]) -> tuple[list[Path], Path | None]:
  """Process a single PNG generation task.

  Args:
    task: Dict containing csv_path, output_base, config, samp

  Returns:
    Tuple of (generated_pngs, gif_output_path or None)
  """
  csv_path = task["csv_path"]
  output_base = task["output_base"]
  config = task["config"]
  samp = task["samp"]

  try:
    df = load_and_filter_data(csv_path, samp)
    generated_pngs = plot_heatmaps(df, output_base, config)

    # Determine GIF output path if animation is enabled
    anim_cfg = config.get("animation", {})
    if anim_cfg.get("tf_exec", False) and generated_pngs:
      # Use output_subdir if specified (relative to PNG directory)
      output_subdir = anim_cfg.get("output_subdir", "")
      if output_subdir:
        gif_dir = output_base.parent / output_subdir
      else:
        gif_dir = output_base.parent
      gif_output = gif_dir / f"{output_base.stem}_anim.gif"
      return (generated_pngs, gif_output)
    return (generated_pngs, None)
  except Exception as e:
    detid = config.get("detid", "??")
    print(f"Error processing det{detid} samp{samp}: {e}")
    return ([], None)


def process_all_png_tasks(png_tasks: list[dict[str, Any]],
                          config: dict[str, Any]) -> list[tuple[list[Path], Path]]:
  """Process all PNG generation tasks, optionally in parallel.

  Note: Uses ProcessPoolExecutor (not ThreadPoolExecutor) because
  matplotlib is not thread-safe.

  Args:
    png_tasks: List of task dicts
    config: Configuration dict containing parallel settings

  Returns:
    List of (generated_pngs, gif_output_path) tuples for animation
  """
  from concurrent.futures import ProcessPoolExecutor, as_completed

  if not png_tasks:
    return []

  parallel_cfg = config.get("parallel", {})
  parallel = parallel_cfg.get("enabled", False)
  max_workers = parallel_cfg.get("max_workers", 4)

  gif_tasks: list[tuple[list[Path], Path]] = []

  if parallel and len(png_tasks) > 1:
    print(f"Generating PNGs in parallel (max_workers={max_workers}, tasks={len(png_tasks)})")
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
      futures = {executor.submit(process_png_task, task): task for task in png_tasks}
      for future in as_completed(futures):
        try:
          generated_pngs, gif_output = future.result()
          if gif_output is not None and generated_pngs:
            gif_tasks.append((generated_pngs, gif_output))
        except Exception as e:
          print(f"Error in PNG task: {e}")
  else:
    for task in png_tasks:
      generated_pngs, gif_output = process_png_task(task)
      if gif_output is not None and generated_pngs:
        gif_tasks.append((generated_pngs, gif_output))

  return gif_tasks


def main():
  parser = argparse.ArgumentParser(
    description="Plot detection limit heatmaps for muography"
  )
  parser.add_argument("--main_config", type=Path, required=True,
                      help="Main JSON5 config file (e.g., prm_reso.json5)")
  parser.add_argument("--plot_config", type=Path, default=None,
                      help="Plot configuration JSON5 file (e.g., heatmap_config.json5)")
  parser.add_argument("--outdir", type=Path, default=Path("figs"),
                      help="Output directory (default: figs)")
  parser.add_argument("--datadir", type=Path, default=Path("tmp"),
                      help="Data directory containing CSV files (default: tmp)")
  parser.add_argument("--size_max", type=float, default=None,
                      help="Maximum size for x-axis (overrides config)")
  parser.add_argument("--depth_max", type=float, default=None,
                      help="Maximum depth for y-axis (overrides config)")
  parser.add_argument("--surface_elevation", type=str, default=None,
                      help='Surface elevation [m] for secondary Y-axis: '
                           'float value, "auto" (from DEM), or omit to disable')
  args = parser.parse_args()

  if not args.main_config.exists():
    print(f"Error: Main config file not found: {args.main_config}")
    sys.exit(1)

  # Load main config to get detector list and signal amplifiers
  main_cfg = load_main_config(args.main_config)

  det_files = main_cfg.get("DETECTOR_PARAMETER_LISTS", {}).get("det_files", [])
  sweep_common = main_cfg.get("DEPTH_RESOLUTION_SWEEP", {}).get("common", {})
  output_prefix = sweep_common.get("output_ascii_prefix", "depth_res")
  signal_amplifiers = sweep_common.get("signal_noise_amplifiers", [[1.0, 0.0]])
  vec_delta_density = sweep_common.get("vec_delta_density", [])
  both_side = sweep_common.get("both_side", False)

  # Extract signal_amp values (first element of each pair)
  signal_amps = [float(sa[0]) for sa in signal_amplifiers]

  if not det_files:
    print("Error: No detector files found in main config")
    sys.exit(1)

  # Load plot config
  plot_config = load_config(args.plot_config)
  if args.size_max is not None:
    plot_config["size_max"] = args.size_max
  if args.depth_max is not None:
    plot_config["depth_max"] = args.depth_max

  # Resolve surface_elevation: "auto" → DEM lookup, float → use as-is, None → disabled
  se = plot_config.get("surface_elevation")
  if args.surface_elevation is not None:
    se = args.surface_elevation
  if se == "auto":
    se = get_surface_elevation_from_dem(args.main_config)
    if se is not None:
      print(f"Surface elevation from DEM: {se:.1f} m")
    else:
      print("Warning: Could not read surface elevation from DEM, disabling elevation axis.")
  elif isinstance(se, str):
    se = float(se)
  plot_config["surface_elevation"] = se

  args.outdir.mkdir(parents=True, exist_ok=True)

  # Collect PNG generation tasks
  png_tasks: list[dict[str, Any]] = []

  # Loop over all detectors and signal amplifiers to collect tasks
  for det_idx, det_file in enumerate(det_files):
    detid = f"{det_idx:02d}"
    detector_name = get_detector_name(args.main_config, detid)

    for samp in signal_amps:
      # Build CSV path: {datadir}/{output_prefix}_det{detid}_signal_signifi-tmp.csv
      csv_path = args.datadir / f"{output_prefix}_det{detid}_signal_signifi-tmp.csv"

      if not csv_path.exists():
        print(f"Warning: CSV not found, skipping: {csv_path}")
        continue

      # Build output base filename
      samp_str = f"{samp:.1E}"
      output_base = args.outdir / f"fig_heatmap_det{detid}_samp{samp_str}.png"

      # Prepare config for this combination
      config = plot_config.copy()
      config["signal_amp"] = samp
      config["detid"] = detid
      config["vec_delta_density"] = vec_delta_density
      config["both_side"] = both_side
      if detector_name is not None:
        config["detector_name"] = detector_name

      png_tasks.append({
        "csv_path": csv_path,
        "output_base": output_base,
        "config": config,
        "samp": samp,
      })

  # Process PNG tasks (sequential or parallel based on config)
  gif_tasks = process_all_png_tasks(png_tasks, plot_config)

  # Create animation GIFs
  if gif_tasks:
    create_all_animation_gifs(gif_tasks, plot_config)


if __name__ == "__main__":
  main()
