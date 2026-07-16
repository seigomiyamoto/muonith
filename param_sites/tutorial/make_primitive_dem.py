#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Generate a simple primitive-solid DEM and write it as a .g2zbin grid.

Supported shapes (single-valued height field z(x, y), vertical axis):
  - cylinder        : flat top at base+height within radius, else base
  - square_prism    : flat top at base+height within half-width (Chebyshev), else base
  - cone            : base + height * (1 - r / radius) for r <= radius, else base
  - square_pyramid  : base + height * (1 - c / half_width) for c <= half_width, else base
  - gaussian        : base + height * exp(-0.5 * ((dx/sigma_x)^2 + (dy/sigma_y)^2))

The g2zbin grid is uniform, so it maps directly onto the format
(Grid2d header + flat z-value array). See scripts/g2zbin_io.py for the layout.

These DEMs are meant to be (re)generated locally by the user for tutorials;
the .g2zbin binaries are not committed to the repository.
"""
import argparse
import sys
from pathlib import Path

import numpy as np

# Resolve the shared g2zbin writer via repo-relative paths (no personal absolute paths).
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]                          # repo root
sys.path.insert(0, str(ROOT / "scripts"))
from g2zbin_io import write_g2zbin


SHAPES = ("cylinder", "square_prism", "cone", "square_pyramid", "gaussian")

# Which per-shape size flag each shape actually uses (for the --width override notice).
SHAPE_SIZE_FLAG = {
  "cylinder": "radius", "cone": "radius",
  "square_prism": "half_width", "square_pyramid": "half_width",
  "gaussian": "sigma_x",
}


def _resolve_size(name: str, explicit, width, factor: float, default: float,
                  notify: bool) -> float:
  """Resolve a horizontal size by precedence: explicit per-shape value >
  --width-derived (width * factor) > built-in default.

  notify: print a one-line notice only when this flag is the one the chosen
  shape uses AND both the explicit flag and --width were given.
  """
  if explicit is not None:
    if width is not None and notify:
      print(f"note: --{name}={explicit:g} overrides --width for this shape")
    return explicit
  if width is not None:
    return width * factor
  return default


def compute_height_field(
    shape: str, xx: np.ndarray, yy: np.ndarray,
    xcnt: float, ycnt: float, base: float, height: float,
    radius: float, half_width: float, sigma_x: float, sigma_y: float,
) -> np.ndarray:
  """Return the surface elevation z(x, y) for the requested primitive.

  xx, yy are meshgrid coordinate arrays (row-major, shape (nbiny, nbinx)).
  """
  dx = xx - xcnt
  dy = yy - ycnt
  z = np.full(xx.shape, base, dtype=np.float64)   # base plane everywhere

  if shape == "cylinder":
    r = np.sqrt(dx * dx + dy * dy)
    z[r <= radius] = base + height
  elif shape == "square_prism":
    cheb = np.maximum(np.abs(dx), np.abs(dy))     # Chebyshev distance
    z[cheb <= half_width] = base + height
  elif shape == "cone":
    r = np.sqrt(dx * dx + dy * dy)
    inside = r <= radius
    z[inside] = base + height * (1.0 - r[inside] / radius)
  elif shape == "square_pyramid":
    cheb = np.maximum(np.abs(dx), np.abs(dy))
    inside = cheb <= half_width
    z[inside] = base + height * (1.0 - cheb[inside] / half_width)
  elif shape == "gaussian":
    z = base + height * np.exp(-0.5 * ((dx / sigma_x) ** 2 + (dy / sigma_y) ** 2))
  else:
    raise ValueError(f"unknown shape: {shape}")

  return z


def main() -> None:
  parser = argparse.ArgumentParser(
    description="Primitive-solid DEM -> g2zbin grid"
  )
  parser.add_argument("--shape", choices=SHAPES, required=True)

  # Plane extent (m)
  parser.add_argument("--x_min", type=float, default=-500.0)
  parser.add_argument("--x_max", type=float, default=500.0)
  parser.add_argument("--y_min", type=float, default=-500.0)
  parser.add_argument("--y_max", type=float, default=500.0)

  # Grid spacing (m)
  parser.add_argument("--x_interval", type=float, default=5.0,
                      help="grid spacing in x direction")
  parser.add_argument("--y_interval", type=float, default=5.0,
                      help="grid spacing in y direction")

  # Shape placement and size (m)
  parser.add_argument("--xcnt", type=float, default=0.0, help="shape center x")
  parser.add_argument("--ycnt", type=float, default=0.0, help="shape center y")
  parser.add_argument("--base", type=float, default=0.0, help="base plane elevation")
  parser.add_argument("--height", type=float, default=300.0, help="peak height above base")
  # Unified horizontal-size knob: cylinder/cone diameter = W, square side = W,
  # gaussian 4*sigma (+/-2 sigma) = W. Explicit per-shape flags below override it.
  parser.add_argument("--width", type=float, default=None,
                      help="unified footprint width W (m): radius=half_width=W/2, sigma=W/4")
  parser.add_argument("--radius", type=float, default=None,
                      help="radius for cylinder / cone (overrides --width; default 200)")
  parser.add_argument("--half_width", type=float, default=None,
                      help="half-width for square_prism / square_pyramid (overrides --width; default 200)")
  parser.add_argument("--sigma_x", type=float, default=None,
                      help="gaussian sigma in x (overrides --width; default 150)")
  parser.add_argument("--sigma_y", type=float, default=None,
                      help="gaussian sigma in y (if omitted, same as sigma_x)")

  # Output
  parser.add_argument("--outbin", type=str, default=None,
                      help="output .g2zbin path (default: HERE/eg_<shape>-<dx>m.g2zbin)")
  parser.add_argument("--float64", action="store_true",
                      help="store z as float64 instead of float32 (default: float32)")

  args = parser.parse_args()

  dx = args.x_interval
  dy = args.y_interval
  if dx <= 0 or dy <= 0:
    raise ValueError("x_interval and y_interval must be positive")
  if args.width is not None and args.width <= 0:
    raise ValueError("width must be positive")

  # Resolve horizontal sizes by precedence: explicit per-shape flag > --width > default.
  # The override notice is printed only for the flag the chosen shape uses.
  used = SHAPE_SIZE_FLAG[args.shape]
  radius = _resolve_size("radius", args.radius, args.width, 0.5, 200.0, used == "radius")
  half_width = _resolve_size("half_width", args.half_width, args.width, 0.5, 200.0,
                             used == "half_width")
  sigma_x = _resolve_size("sigma_x", args.sigma_x, args.width, 0.25, 150.0, used == "sigma_x")
  sigma_y = args.sigma_y if args.sigma_y is not None else sigma_x
  if radius <= 0 or half_width <= 0:
    raise ValueError("radius and half_width must be positive")
  if sigma_x <= 0 or sigma_y <= 0:
    raise ValueError("sigma_x and sigma_y must be positive")

  # Auto-compute grid counts (endpoints included)
  nx = int(round((args.x_max - args.x_min) / dx)) + 1
  ny = int(round((args.y_max - args.y_min) / dy)) + 1

  # Build the uniform grid (row-major: Z[iy, ix])
  x = args.x_min + np.arange(nx, dtype=np.float64) * dx
  y = args.y_min + np.arange(ny, dtype=np.float64) * dy
  xx, yy = np.meshgrid(x, y)

  Z = compute_height_field(
    args.shape, xx, yy,
    args.xcnt, args.ycnt, args.base, args.height,
    radius, half_width, sigma_x, sigma_y,
  )

  # v2 semantics: grid points are bin centers, so min/max must be the outer
  # bin edges (half an interval outside the outermost points).
  grid2d_info = {
    "name": args.shape,
    "x_axis": {
      "name": "x", "nbin": nx,
      "min": float(x[0]) - 0.5 * dx, "max": float(x[-1]) + 0.5 * dx, "interval": dx,
    },
    "y_axis": {
      "name": "y", "nbin": ny,
      "min": float(y[0]) - 0.5 * dy, "max": float(y[-1]) + 0.5 * dy, "interval": dy,
    },
  }

  outbin = args.outbin
  if outbin is None:
    outbin = str(HERE / f"eg_{args.shape}-{int(round(dx))}m.g2zbin")

  use_f32 = not args.float64
  write_g2zbin(outbin, grid2d_info, Z, use_float32=use_f32)

  print(f"Wrote {nx * ny} cells to {outbin} ({nx}x{ny} grid, "
        f"z range [{Z.min():.3f}, {Z.max():.3f}] m)")


if __name__ == "__main__":
  main()
