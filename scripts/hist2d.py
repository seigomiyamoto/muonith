#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
hist2d.py - Unified 2D histogram / grid plotter (numpy-only).

Merges hist2d_v08.py (grid plotter) and hist2d_v10.py (histogram plotter)
into a single script with two modes:

  --mode histogram  (default)  Bin scatter data with np.histogram2d.
  --mode grid                  Reshape pre-gridded (x,y,z) data for display.

Features carried over from v10:
  - numpy-only (no pandas dependency)
  - CSV output, mask overlay, colormap/color-over/color-under
  - Log scale, contour with interval, DPI/title/output control

Features carried over from v08:
  - Grid mode (pivot-like reshape without pandas)
  - Auto-z (data-driven color range)
  - Allow-mismatch (grid shape tolerance)
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import argparse
import copy
import sys
import os
from pathlib import Path


# ============================================================
# Argument parsing
# ============================================================
def build_parser():
  """Build argument parser with both triplet and individual bin options."""
  parser = argparse.ArgumentParser(
    description='Unified 2D histogram / grid plotter (numpy-only)',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  # Input file
  parser.add_argument('filename', help='Input file with 3 columns: x, y, z')

  # Mode
  parser.add_argument('--mode', choices=['histogram', 'grid'], default='histogram',
                      help='Plot mode: histogram (bin scatter data) or grid (reshape pre-gridded data)')

  # Triplet bin spec (nargs=3, type=str to avoid argparse negative-number issues)
  parser.add_argument('--binx', nargs=3, metavar=('NBINX', 'XMIN', 'XMAX'),
                      help='X binning as triplet: NBINX XMIN XMAX')
  parser.add_argument('--biny', nargs=3, metavar=('NBINY', 'YMIN', 'YMAX'),
                      help='Y binning as triplet: NBINY YMIN YMAX')
  parser.add_argument('--binz', nargs=3, metavar=('NGRAD', 'VMIN', 'VMAX'),
                      help='Z (color) binning as triplet: NGRAD VMIN VMAX')

  # Individual bin spec
  parser.add_argument('--nbinx', type=int, default=50, help='Number of X bins')
  parser.add_argument('--xmin', type=float, default=-1.0, help='X minimum')
  parser.add_argument('--xmax', type=float, default=1.0, help='X maximum')
  parser.add_argument('--nbiny', type=int, default=50, help='Number of Y bins')
  parser.add_argument('--ymin', type=float, default=-1.0, help='Y minimum')
  parser.add_argument('--ymax', type=float, default=1.0, help='Y maximum')
  parser.add_argument('--ngrad', type=int, default=256, help='Number of color gradations')
  parser.add_argument('--vmin', type=float, default=0.0, help='Color scale minimum')
  parser.add_argument('--vmax', type=float, default=1.0, help='Color scale maximum')

  # Output settings
  parser.add_argument('--output', default=None,
                      help='Output filename (default: auto-generated from input name)')
  parser.add_argument('--dpi', type=int, default=300, help='Output DPI')
  parser.add_argument('--title', type=str, default=None, help='Plot title (default: filename)')

  # Colormap settings
  parser.add_argument('--colormap', type=str, default='jet', help='Matplotlib colormap name')
  parser.add_argument('--color-over', type=str, default='black', help='Color for values > vmax')
  parser.add_argument('--color-under', type=str, default='white', help='Color for values < vmin')

  # Display options
  parser.add_argument('--log', action='store_true', help='Use logarithmic color scale')
  parser.add_argument('--auto-z', action='store_true',
                      help='Determine z range from data (ignores --vmin/--vmax)')
  parser.add_argument('--contour', action='store_true', help='Draw contour lines')
  parser.add_argument('--contour-interval', type=float, default=None,
                      help='Contour line interval (if not set, use ngrad divisions)')
  parser.add_argument('--aspect', choices=['equal', 'auto'], default='equal',
                      help="Axes aspect: 'equal' (1:1 data units, default) or 'auto' (fill figure box)")
  parser.add_argument('--figsize', nargs=2, type=float, default=None, metavar=('W', 'H'),
                      help='Figure size in inches; overrides the aspect-derived size')
  parser.add_argument('--xlabel', type=str, default=None,
                      help="X axis label (default: 'X')")
  parser.add_argument('--ylabel', type=str, default=None,
                      help="Y axis label (default: 'Y')")
  parser.add_argument('--zlabel', type=str, default=None,
                      help="Colorbar label (default: 'z')")

  # CSV output control (histogram mode only)
  parser.add_argument('--no-csv', action='store_true', help='Disable CSV output (histogram mode)')

  # Mask overlay (histogram mode only)
  parser.add_argument('--mask', type=str, default=None,
                      help='Mask file (x, y, 0/1 format). Cells with 0 are overlaid with white')
  parser.add_argument('--mask-alpha', type=float, default=0.3,
                      help='Mask overlay transparency (0=transparent, 1=opaque)')

  # Grid mode options
  parser.add_argument('--allow-mismatch', action='store_true',
                      help='[grid mode] Allow grid shape mismatch: use specified range with data bin count')

  return parser


def resolve_bins(args):
  """Override individual bin values with triplet values if provided."""
  if args.binx is not None:
    args.nbinx = int(args.binx[0])
    args.xmin = float(args.binx[1])
    args.xmax = float(args.binx[2])
  if args.biny is not None:
    args.nbiny = int(args.biny[0])
    args.ymin = float(args.biny[1])
    args.ymax = float(args.biny[2])
  if args.binz is not None:
    args.ngrad = int(args.binz[0])
    args.vmin = float(args.binz[1])
    args.vmax = float(args.binz[2])


# ============================================================
# Data loading
# ============================================================
def load_data(filename):
  """Load whitespace-separated x y z data. Returns shape (N, 3) or (0, 3)."""
  try:
    data = np.loadtxt(filename, comments='#')
    if data.ndim == 1:
      data = data.reshape(1, -1)
    if data.size == 0:
      return np.empty((0, 3))
    return data
  except Exception as e:
    print(f"Warning: Failed to load data ({e}), generating empty plot")
    return np.empty((0, 3))


# ============================================================
# Histogram mode (from v10)
# ============================================================
def compute_histogram(data, args):
  """Compute 2D histogram using np.histogram2d with weights.

  Returns:
    H:      shape (nbinx, nbiny) weighted histogram.
    xedges: shape (nbinx+1,) bin edges along x.
    yedges: shape (nbiny+1,) bin edges along y.
  """
  if data.size == 0:
    H = np.zeros((args.nbinx, args.nbiny))
    xedges = np.linspace(args.xmin, args.xmax, args.nbinx + 1)
    yedges = np.linspace(args.ymin, args.ymax, args.nbiny + 1)
    return H, xedges, yedges

  x, y, val = data[:, 0], data[:, 1], data[:, 2]
  H, xedges, yedges = np.histogram2d(
    x, y,
    bins=[args.nbinx, args.nbiny],
    range=[[args.xmin, args.xmax], [args.ymin, args.ymax]],
    weights=val
  )
  return H, xedges, yedges


def save_csv(H, xedges, yedges, args, output_csv):
  """Save histogram data to space-separated CSV format."""
  xstep = (args.xmax - args.xmin) / args.nbinx
  ystep = (args.ymax - args.ymin) / args.nbiny

  with open(output_csv, mode='w') as f:
    f.write(f"{args.nbinx:d} {args.xmin:E} {args.xmax:E} {xstep:E}\n")
    f.write(f"{args.nbiny:d} {args.ymin:E} {args.ymax:E} {ystep:E}\n")
    bin_id = 0
    for ix in range(args.nbinx):
      xlower = args.xmin + ix * xstep
      xupper = xlower + xstep
      for iy in range(args.nbiny):
        ylower = args.ymin + iy * ystep
        yupper = ylower + ystep
        value = H[ix, iy]
        f.write(f"{bin_id:d} {xlower:E} {xupper:E} {ylower:E} {yupper:E} {value:E} 1\n")
        bin_id += 1
  print(f"CSV file saved as {output_csv}")


# ============================================================
# Grid mode (from v08, pandas-free)
# ============================================================
def build_grid_from_scatter(data):
  """Convert scatter (x, y, z) triplets to a 2D grid using numpy.

  Duplicates at the same (x, y) are averaged (equivalent to
  pd.pivot_table with aggfunc=np.mean).

  Returns:
    x_centers: shape (nx,) sorted unique x values.
    y_centers: shape (ny,) sorted unique y values.
    Z:         shape (nx, ny), NaN where no data.
  """
  x_all, y_all, z_all = data[:, 0], data[:, 1], data[:, 2]
  x_centers = np.unique(x_all)
  y_centers = np.unique(y_all)
  nx, ny = len(x_centers), len(y_centers)

  accum = np.zeros((nx, ny), dtype=float)
  count = np.zeros((nx, ny), dtype=int)
  ix = np.searchsorted(x_centers, x_all)
  iy = np.searchsorted(y_centers, y_all)
  np.add.at(accum, (ix, iy), z_all)
  np.add.at(count, (ix, iy), 1)

  Z = np.full((nx, ny), np.nan, dtype=float)
  mask = count > 0
  Z[mask] = accum[mask] / count[mask]
  return x_centers, y_centers, Z


def infer_edges_from_centers(centers):
  """Derive bin edges from sorted center coordinates.

  Interior edges are midpoints.  Exterior edges are extrapolated using the
  mean spacing.  A single element yields [c - 0.5, c + 0.5].
  """
  if centers.size == 1:
    return np.array([centers[0] - 0.5, centers[0] + 0.5])
  d = np.diff(centers)
  dmean = float(np.mean(d)) if d.size > 0 else 1.0
  return np.concatenate(([centers[0] - 0.5 * dmean],
                         0.5 * (centers[:-1] + centers[1:]),
                         [centers[-1] + 0.5 * dmean]))


def build_mesh_edges(x_centers, y_centers, args):
  """Build bin-edge arrays for grid mode.

  If nbinx/nbiny > 0, use np.linspace from specified range.
  Otherwise, infer edges from the data centers.

  Returns:
    x_edges: 1D edge array.
    y_edges: 1D edge array.
    mismatch: True if shape mismatch was detected (and --allow-mismatch used).
  """
  # X edges
  if args.nbinx > 0:
    x_edges = np.linspace(args.xmin, args.xmax, args.nbinx + 1)
  else:
    x_edges = infer_edges_from_centers(x_centers)

  # Y edges
  if args.nbiny > 0:
    y_edges = np.linspace(args.ymin, args.ymax, args.nbiny + 1)
  else:
    y_edges = infer_edges_from_centers(y_centers)

  # Check shape match
  nx_spec = len(x_edges) - 1
  ny_spec = len(y_edges) - 1
  nx_data = len(x_centers)
  ny_data = len(y_centers)
  mismatch = False

  if nx_spec != nx_data or ny_spec != ny_data:
    msg = (f"grid shape ({ny_spec}, {nx_spec}) does not match "
           f"data shape ({ny_data}, {nx_data})")
    if args.allow_mismatch:
      print(f"WARNING: {msg}. Applying bin-spec range with data grid (--allow-mismatch).")
      mismatch = True
    else:
      raise ValueError(msg + ". Rerun with --allow-mismatch to override.")

  return x_edges, y_edges, mismatch


def allow_mismatch_remap(x_centers, y_centers, Z, x_edges, y_edges):
  """Remap a data grid onto a target bin specification using np.digitize.

  Data cells outside the specified range are discarded (NaN).

  Returns:
    Z_new of shape (len(x_edges)-1, len(y_edges)-1).
  """
  nx_tgt = len(x_edges) - 1
  ny_tgt = len(y_edges) - 1
  Z_new = np.full((nx_tgt, ny_tgt), np.nan, dtype=float)

  x_bins = np.digitize(x_centers, x_edges) - 1
  y_bins = np.digitize(y_centers, y_edges) - 1

  for ix_src, xb in enumerate(x_bins):
    if xb < 0 or xb >= nx_tgt:
      continue
    for iy_src, yb in enumerate(y_bins):
      if yb < 0 or yb >= ny_tgt:
        continue
      Z_new[xb, yb] = Z[ix_src, iy_src]

  return Z_new


# ============================================================
# Visualization (unified)
# ============================================================
EPSILON = 1.0e-6


def setup_colormap(args):
  """Create a colormap copy with bad/under/over colors set."""
  cmap = copy.copy(plt.colormaps[args.colormap])
  cmap.set_bad('silver')
  cmap.set_under(args.color_under)
  cmap.set_over(args.color_over)
  return cmap


def compute_norm(data_values, args):
  """Compute matplotlib normalization.

  - --auto-z: vmin/vmax determined from data (ignoring NaN).
  - --log:    LogNorm (vmin clamped to EPSILON).
  - else:     BoundaryNorm with ngrad+1 boundaries.

  Returns:
    norm: a matplotlib Normalize subclass.
  """
  if args.auto_z:
    valid = data_values[np.isfinite(data_values)]
    if valid.size == 0:
      vmin_use, vmax_use = EPSILON, 1.0
    else:
      vmin_use = float(np.min(valid))
      vmax_use = float(np.max(valid))
      if vmin_use == 0.0:
        vmin_use = EPSILON
      vmax_use = max(vmax_use, vmin_use * (1.0 + EPSILON))
  else:
    vmin_use = args.vmin
    vmax_use = args.vmax

  if args.log:
    vmin_plot = max(vmin_use, EPSILON)
    vmax_plot = max(vmax_use, vmin_plot * (1.0 + EPSILON))
    print(f"[LOG] vmin={vmin_plot:E} vmax={vmax_plot:E}")
    return mpl.colors.LogNorm(vmin=vmin_plot, vmax=vmax_plot)
  else:
    print(f"vmin={vmin_use:E} vmax={vmax_use:E}")
    bounds = np.linspace(vmin_use, vmax_use, args.ngrad + 1)
    return mpl.colors.BoundaryNorm(bounds, 256, extend='both')


def compute_contour_levels(args):
  """Compute contour level array from args.

  Returns None if --contour is not set or levels cannot be determined.
  """
  if not args.contour:
    return None

  vmin_eff = args.vmin
  vmax_eff = args.vmax
  if args.auto_z:
    # Contour with auto-z: cannot determine levels a priori
    print("WARNING: --contour with --auto-z is not supported. Skipping contour.")
    return None

  if args.log:
    vmin_eff = max(vmin_eff, EPSILON)
    vmax_eff = max(vmax_eff, vmin_eff * (1.0 + EPSILON))
    if args.contour_interval is not None:
      log_min = np.log10(vmin_eff)
      log_max = np.log10(vmax_eff)
      n_levels = max(2, int((log_max - log_min) /
                             np.log10(1 + args.contour_interval / max(vmin_eff, EPSILON))))
      return np.logspace(log_min, log_max, n_levels)
    else:
      return np.logspace(np.log10(vmin_eff), np.log10(vmax_eff), args.ngrad + 1)
  else:
    if args.contour_interval is not None:
      levels = np.arange(vmin_eff, vmax_eff + args.contour_interval * 0.5,
                         args.contour_interval)
      if levels.size == 0:
        print("WARNING: Failed to build contour levels. Skipping contour.")
        return None
      return levels
    else:
      return np.linspace(vmin_eff, vmax_eff, args.ngrad + 1)


def determine_output_path(args):
  """Determine output file path.

  --output given:  figs/{output}
  --output None:   figs/fig_{stem}{_log}.png  (auto-generated)
  """
  fig_dir = Path("figs")
  fig_dir.mkdir(parents=True, exist_ok=True)
  if args.output is not None:
    return fig_dir / args.output
  stem = Path(args.filename).stem
  suffix = "_log" if args.log else ""
  return fig_dir / f"fig_{stem}{suffix}.png"


def render_plot(X, Y, Z_plot, x_edges, y_edges, norm, cmap, contour_levels, args,
                mask_mesh=None):
  """Render pcolormesh plot and save to file.

  Args:
    X, Y:            2D meshgrid arrays for pcolormesh (bin centers).
    Z_plot:          masked array, shape (ny, nx).
    x_edges, y_edges: 1D edge arrays for axis limits.
    norm:            matplotlib Normalize object.
    cmap:            matplotlib Colormap.
    contour_levels:  np.ndarray or None.
    args:            parsed arguments.
    mask_mesh:       (X, Y) meshgrid tuple for mask overlay, or None.

  Returns:
    output file path string.
  """
  xmin, xmax = float(x_edges.min()), float(x_edges.max())
  ymin, ymax = float(y_edges.min()), float(y_edges.max())

  # Figure size (equal aspect, 10-inch wide)
  hori_inch = 10
  aspect_ratio = (xmax - xmin) / (ymax - ymin) if (ymax - ymin) != 0 else 1.0
  vert_inch = hori_inch / aspect_ratio if aspect_ratio != 0 else hori_inch

  # Font sizes
  title_fontsize = 18
  label_fontsize = 16
  tick_fontsize = 14

  fig, ax = plt.subplots()
  ax.set_xlim(xmin, xmax)
  ax.set_ylim(ymin, ymax)
  ax.set_aspect(args.aspect)
  if args.figsize is not None:
    fig.set_size_inches(args.figsize[0], args.figsize[1])
  elif args.aspect == "equal":
    fig.set_size_inches(hori_inch, vert_inch)
  else:
    fig.set_size_inches(hori_inch, hori_inch * 0.75)

  if args.mode == 'grid':
    ax.set_facecolor('silver')

  # Main plot
  mesh = ax.pcolormesh(X, Y, Z_plot, cmap=cmap, norm=norm, shading='auto')

  # Contour lines
  if contour_levels is not None:
    try:
      ax.contour(X, Y, Z_plot, levels=contour_levels, colors='black', linewidths=0.6)
    except ValueError as e:
      print(f"WARNING: Failed to draw contour lines ({e}).")

  # Colorbar
  cbar = plt.colorbar(mesh, ax=ax, extend='both')
  cbar.set_label(args.zlabel if args.zlabel else 'z', fontsize=label_fontsize)
  cbar.ax.tick_params(labelsize=tick_fontsize)

  # Mask overlay (histogram mode only)
  if args.mask is not None and mask_mesh is not None:
    try:
      mask_data = np.loadtxt(args.mask)
      if mask_data.size != 0:
        X_m, Y_m = mask_mesh
        x_mask, y_mask, mask_val = mask_data[:, 0], mask_data[:, 1], mask_data[:, 2]
        H_mask, _, _ = np.histogram2d(
          x_mask, y_mask,
          bins=[args.nbinx, args.nbiny],
          range=[[args.xmin, args.xmax], [args.ymin, args.ymax]],
          weights=mask_val
        )
        overlay = np.where(H_mask == 0, 0, np.nan)
        white_cmap = mpl.colors.ListedColormap(['white'])
        ax.pcolormesh(X_m, Y_m, np.ma.masked_invalid(overlay).T,
                      cmap=white_cmap, shading='auto', alpha=args.mask_alpha)
        print(f"Mask overlay applied from {args.mask}")
    except Exception as e:
      print(f"Warning: Failed to load mask file ({e})")

  # Grid lines
  x_ticks = ax.get_xticks()
  y_ticks = ax.get_yticks()
  ax.vlines(x_ticks, ymin=ymin, ymax=ymax, colors='grey', linestyles='dotted', linewidth=1.0)
  ax.hlines(y_ticks, xmin=xmin, xmax=xmax, colors='grey', linestyles='dotted', linewidth=1.0)

  # Labels
  ax.set_xlabel(args.xlabel if args.xlabel else 'X', fontsize=label_fontsize)
  ax.set_ylabel(args.ylabel if args.ylabel else 'Y', fontsize=label_fontsize)
  title = args.title if args.title else os.path.basename(args.filename)
  ax.set_title(title, fontsize=title_fontsize)
  ax.tick_params(axis='both', which='major', labelsize=tick_fontsize)

  # Save
  output_path = determine_output_path(args)
  plt.savefig(output_path, dpi=args.dpi)
  plt.close(fig)
  return str(output_path)


# ============================================================
# Main
# ============================================================
def main():
  parser = build_parser()
  args = parser.parse_args()

  # Resolve triplet -> individual
  resolve_bins(args)

  # Default title
  if args.title is None:
    args.title = os.path.basename(args.filename)

  # Load data
  data = load_data(args.filename)

  if args.mode == 'histogram':
    # ---- Histogram mode (v10 base) ----
    H, xedges, yedges = compute_histogram(data, args)

    # CSV output (before masking)
    if not args.no_csv:
      output_base = args.output if args.output else f"fig_{Path(args.filename).stem}.png"
      output_csv = output_base.replace('.png', '.csv')
      save_csv(H, xedges, yedges, args, output_csv)

    # Mask zero values for display
    H_masked = np.ma.masked_where(H == 0, H)

    # Compute norm
    norm = compute_norm(H.ravel(), args)

    # Meshgrid (bin centers)
    X, Y = np.meshgrid(
      0.5 * (xedges[:-1] + xedges[1:]),
      0.5 * (yedges[:-1] + yedges[1:])
    )

    cmap = setup_colormap(args)
    contour_levels = compute_contour_levels(args)

    output_path = render_plot(X, Y, H_masked.T, xedges, yedges,
                              norm, cmap, contour_levels, args,
                              mask_mesh=(X, Y))

  elif args.mode == 'grid':
    # ---- Grid mode (v08 base, pandas-free) ----
    if data.size == 0:
      print("Error: No data loaded for grid mode.")
      return 1

    x_centers, y_centers, Z = build_grid_from_scatter(data)

    # v08 compat: treat zero z-values as NaN
    Z = np.where(Z == 0.0, np.nan, Z)

    # Build mesh edges
    x_edges, y_edges, mismatch = build_mesh_edges(x_centers, y_centers, args)

    if mismatch:
      Z = allow_mismatch_remap(x_centers, y_centers, Z, x_edges, y_edges)
      x_plot = 0.5 * (x_edges[:-1] + x_edges[1:])
      y_plot = 0.5 * (y_edges[:-1] + y_edges[1:])
    else:
      x_plot = x_centers
      y_plot = y_centers

    # Compute norm (use finite values for auto-z)
    Z_flat = Z.ravel()
    norm = compute_norm(Z_flat, args)

    # Mask NaN for display
    Z_masked = np.ma.masked_invalid(Z)

    # Meshgrid
    X, Y = np.meshgrid(x_plot, y_plot)

    cmap = setup_colormap(args)
    contour_levels = compute_contour_levels(args)

    output_path = render_plot(X, Y, Z_masked.T, x_edges, y_edges,
                              norm, cmap, contour_levels, args)

  print(f"SUCCESS!! The plot filename is {output_path}")
  return 0


if __name__ == '__main__':
  sys.exit(main())
