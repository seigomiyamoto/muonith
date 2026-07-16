#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Plot g2zbin files as 2D heatmaps with optional contours.

Reads .g2zbin files and visualizes them with
x on the horizontal axis, y on the vertical axis, and z as color.
Auto-saves to PNG.  Supports --zmin/--zmax color range,
--contour/--cint contour lines, and --rotate_ccw90 rotation.
"""

import os
import numpy as np
import matplotlib.pyplot as plt

from g2zbin_io import axis_centers, read_g2zbin


def load_g2zbin_as_grid(path):
  """Read a .g2zbin file and return (xi, yi, Z).

  Args:
    path: Path to the .g2zbin file.

  Returns:
    (xi, yi, Z) where xi/yi are 1D bin-center coordinate arrays
    (v2 min/max are bin edges) and Z is (ny, nx).
  """
  info, Z = read_g2zbin(path)
  xi = axis_centers(info["x_axis"])
  yi = axis_centers(info["y_axis"])
  return xi, yi, Z


def rotate_ccw90_grid(xi, yi, Z):
  """Rotate grid data 90 degrees counter-clockwise"""
  Zr = np.rot90(Z)
  xir = yi.copy()
  yir = xi.copy()
  return xir, yir, Zr


def plot_dem(input_path, title=None, zmin=None, zmax=None,
                contour=False, cint=None, rotate_ccw90=False, cmap="terrain"):
  ext = os.path.splitext(input_path)[1].lower()
  if ext != ".g2zbin":
    raise ValueError(f"Unsupported input format: {ext}. Only .g2zbin is supported.")
  xi, yi, Z = load_g2zbin_as_grid(input_path)

  # Rotation option
  if rotate_ccw90:
    xi, yi, Z = rotate_ccw90_grid(xi, yi, Z)

  fig, ax = plt.subplots(figsize=(8, 6))

  extent = [xi.min(), xi.max(), yi.min(), yi.max()]
  im = ax.imshow(Z, extent=extent, origin='lower', aspect='equal',
                 vmin=zmin, vmax=zmax, cmap=cmap, interpolation='nearest')
  ax.grid(color="grey", linestyle=":", linewidth=1.0, alpha=0.6)
  plt.colorbar(im, ax=ax, label="Z value")

  ax.set_xlabel("X (meters)" if not rotate_ccw90 else "Y (meters)")
  ax.set_ylabel("Y (meters)" if not rotate_ccw90 else "X (meters)")
  if title:
    ax.set_title(title)

  # Contour lines
  if contour:
    levels = None
    if cint:
      zmin_eff = zmin if zmin is not None else float(np.nanmin(Z))
      zmax_eff = zmax if zmax is not None else float(np.nanmax(Z))
      levels = np.arange(np.floor(zmin_eff / cint) * cint,
                         np.ceil(zmax_eff / cint) * cint + 0.5 * cint, cint)
    cs = ax.contour(xi, yi, Z, levels=levels, colors='black',
                    linewidths=0.5, alpha=0.6)
    ax.clabel(cs, fmt="%.0f", fontsize=6, inline=True)
    print(f"Contours drawn (pitch={cint if cint else 'auto'})")

  plt.tight_layout()

  png_path = os.path.splitext(input_path)[0] + ".png"
  plt.savefig(png_path, dpi=300)
  print(f"Figure saved: {png_path}")


if __name__ == "__main__":
  import argparse

  p = argparse.ArgumentParser(description="Visualize g2zbin as a 2D heatmap with optional contours")
  p.add_argument("--input", required=True, help="Input .g2zbin file path")
  p.add_argument("--title", default=None, help="Figure title (optional)")
  p.add_argument("--zmin", type=float, default=None, help="Minimum value of the color scale (optional)")
  p.add_argument("--zmax", type=float, default=None, help="Maximum value of the color scale (optional)")
  p.add_argument("--contour", action="store_true", help="Draw contour lines")
  p.add_argument("--cint", type=float, default=10, help="Contour interval (in Z units)")
  p.add_argument("--rotate_ccw90", action="store_true",
                 help="Display rotated 90 degrees counter-clockwise (no rotation by default)")
  p.add_argument("--cmap", default="terrain", help="matplotlib colormap name (default: terrain)")

  args = p.parse_args()

  plot_dem(
    input_path=args.input,
    title=args.title,
    zmin=args.zmin,
    zmax=args.zmax,
    contour=args.contour,
    cint=args.cint,
    rotate_ccw90=args.rotate_ccw90,
    cmap=args.cmap
  )
