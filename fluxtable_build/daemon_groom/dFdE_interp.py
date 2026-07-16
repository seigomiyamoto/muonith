#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Convert  (log10P, cos(theta), dF/dP) table  ->  (cos(theta), log10KE, log10(dF/dE))
- KE = sqrt(P^2 + m^2) - m    where  m = 0.105658 GeV
- dF/dE = (dF/dP) * (E / P)
- Target grid:
    cos(theta): 0.000 to 1.000 in 0.001 steps
    log10(KE): -1.0 to 7.0 in 0.01 steps
- 2D spline interpolation:
    * costTheta (x) direction: cubic (kx=3)
    * log10KE (y) direction: linear (ky=1)
- Extrapolation allowed in y
- Output is sorted by x then y
- Also provides cross-section plotting functions

usage: python3 make_muon_flux_ke_table.py <input file> <output file> -j <n_threads>
"""

import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RectBivariateSpline
from multiprocessing import cpu_count
from pathlib import Path
from types import SimpleNamespace
import json

MUON_MASS_GEV = 0.105658  # GeV/c²

# Source-flux area-unit correction (cm^-2 -> m^-2).
# The input daemonflux.txt header is labeled "(m^2 sr s GeV/c)^-1", but the
# actual values are per cm^2: the vertical ~1 GeV/c value 2.50e-3 matches the
# textbook PDG sea-level vertical muon momentum spectrum in cm^-2 s^-1 sr^-1 (GeV/c)^-1.
# muonith uses m^-2 (consistent with the Honda tables), so multiply by 1e4.
CM2_TO_M2 = 1.0e4

# ----------------------------------------------------------------------
def parse_args():
  parser = argparse.ArgumentParser(
    description="2D spline (cubic x, linear y) with extrapolation and cross-section plotting")
  parser.add_argument("--infile",  help="input flux table")
  parser.add_argument("--outfile", help="output file name")
  parser.add_argument("--n_threads", type=int, default=cpu_count(),
    help="number of parallel processes (unused)")
  parser.add_argument("--json", type=str, default=None,
    help="JSON file specifying costhz[] and logKE[]")
  parser.add_argument("--smoothing", type=float, default=5.0e-5,
    help="smoothing factor (scaled by number of points), default=5.0e-5")

  return parser.parse_args()
# ----------------------------------------------------------------------
def momentum_to_ke(P_GeV):
  """P [GeV/c] → kinetic energy [GeV]"""
  E_total = np.sqrt(P_GeV * P_GeV + MUON_MASS_GEV * MUON_MASS_GEV)
  return E_total - MUON_MASS_GEV, E_total

# ----------------------------------------------------------------------
def dFdP_to_dFdE(P_GeV, dFdP, E_total):
  """convert dF/dP to dF/dE"""
  return dFdP * (E_total / P_GeV)

# ----------------------------------------------------------------------
def load_flux_table(path):
  """read table, return (costTheta_array, log10P_array, log10KE_array, dFdE_grid)"""
  try:
    data = np.loadtxt(path, comments="#")
  except Exception as e:
    sys.exit(f"[ERROR] failed to read {path}: {e}")

  log10P, costTheta, dFdP = data.T
  dFdP = dFdP * CM2_TO_M2  # correct mislabeled source units: cm^-2 -> m^-2
  P_GeV = 10.0 ** log10P
  KE, E_total = momentum_to_ke(P_GeV)
  dFdE = dFdP_to_dFdE(P_GeV, dFdP, E_total)

  uniq_log10P     = np.unique(log10P)
  uniq_costTheta  = np.unique(costTheta)
  nx, ny          = len(uniq_log10P), len(uniq_costTheta)

  if nx * ny != dFdE.size:
    sys.exit("[ERROR] table does not form a regular grid")

  dFdE_grid = dFdE.reshape((ny, nx))  # shape = (costTheta, log10P)
  KE_grid, _ = momentum_to_ke(10.0 ** uniq_log10P)
  log10KE = np.log10(KE_grid.clip(min=1e-10))

  return uniq_costTheta, uniq_log10P, log10KE, dFdE_grid

# ----------------------------------------------------------------------
def make_target_axes():
  costTheta_tgt = np.arange(0.000, 1.000 + 1e-9, 0.001)
  log10KE_tgt   = np.arange(-1.0, 7.0   + 1e-9, 0.01)
  return costTheta_tgt, log10KE_tgt

# ----------------------------------------------------------------------
def build_interpolator(costTheta_src, log10KE_src, dFdE_grid, smoothing_scale):
  """bicubic x, linear y spline in (costTheta, log10KE) with extrapolation"""
  try:
    Z         = np.log10(np.clip(dFdE_grid, 1e-30, None))
    npts      = Z.size
    kx, ky    = 3, 1
    smoothing = smoothing_scale * npts

    print(f"[INFO] Building RectBivariateSpline:")
    print(f"       - input grid shape = {Z.shape}")
    print(f"       - smoothing s = {smoothing:.3e}")
    print(f"       - kx = {kx} (cubic in x), ky = {ky} (linear in y)")

    spline = RectBivariateSpline(costTheta_src, log10KE_src, Z, kx=kx, ky=ky, s=smoothing)

    return SimpleNamespace(
      interp=spline,
      kx=kx,
      ky=ky,
      s=smoothing
    )

  except Exception as e:
    sys.exit(f"[ERROR] failed to build spline: {e}")

# ----------------------------------------------------------------------
def write_output(outfile, spline, costTheta_tgt, log10KE_tgt):
  """evaluate spline on full grid and write sorted by x then y"""
  Z_out = spline.interp(costTheta_tgt, log10KE_tgt)  # shape = (len(x), len(y))

  # ======== ASCII output ========
  lines = []
  for i, ct in enumerate(costTheta_tgt):
    ct_fmt = 0.0 if abs(ct) < 1e-10 else ct
    for j, ly in enumerate(log10KE_tgt):
      ly_fmt = 0.0 if abs(ly) < 1e-10 else ly
      val    = Z_out[i, j]
      lines.append(f"{ct_fmt: .6E} {ly_fmt: .6E} {val: .6E}")

  try:
    Path(outfile).write_text("\n".join(lines) + "\n")
  except Exception as e:
    sys.exit(f"[ERROR] cannot write {outfile}: {e}")

# ----------------------------------------------------------------------
def plot_costhz_slice(spline, costTheta_src, log10KE_src, dFdE_grid, costhz_in):
  """plot z vs log10KE at fixed costhz_in"""
  mask = np.isclose(costTheta_src, costhz_in, atol=1e-6)
  if not mask.any():
    nearest = costTheta_src[np.argmin(np.abs(costTheta_src - costhz_in))]
    sys.exit(f"[ERROR] costhz_in={costhz_in} not found in original data. "
             f"Closest available: {nearest}")

  idx = np.where(mask)[0][0]
  orig   = np.log10(np.clip(dFdE_grid[idx, :], 1e-30, None))
  interp = spline.interp(costhz_in, log10KE_src, grid=False)

  fig_ = plt.figure()
  ax   = fig_.add_subplot(1, 1, 1)
  ax.plot(log10KE_src, interp, '-', label='interp',
          linewidth=0.5, color='black', alpha=1.0)
  ax.plot(log10KE_src, orig, '+', label='original',
          markersize=6, color='blue', alpha=0.3)
  ax.set_xlabel('log10KE')
  ax.set_ylabel('log10 dFdE')

  title = f"costTheta = {costhz_in:.3f}  |  kx={spline.kx}, ky={spline.ky}, s={spline.s:.2e}"
  ax.set_title(title)
  ax.grid(True, which='both', linestyle=':', color='silver', linewidth=0.8)
  ax.legend()

  Path("figs").mkdir(parents=True, exist_ok=True)
  fig_.savefig(f"figs/fig_logdFdE_costhz{costhz_in:.3f}.png", dpi=150)
# ----------------------------------------------------------------------
def plot_logP_slice(spline, costTheta_src, log10P_src, dFdE_grid, logp_in):
  """plot z vs cos(theta) at fixed log10P = logp_in (i.e. original axis)"""
  mask = np.isclose(log10P_src, logp_in, atol=1e-6)
  if not mask.any():
    nearest = log10P_src[np.argmin(np.abs(log10P_src - logp_in))]
    sys.exit(f"[ERROR] log10P_in={logp_in} not found in original data. Closest available: {nearest}")

  j = np.where(mask)[0][0]
  orig   = np.log10(np.clip(dFdE_grid[:, j], 1e-30, None))  # shape = (costTheta)
  log10KE_val = np.log10(np.sqrt((10**logp_in)**2 + MUON_MASS_GEV**2) - MUON_MASS_GEV)
  interp = spline.interp(costTheta_src, log10KE_val, grid=False)

  fig_ = plt.figure()
  ax   = fig_.add_subplot(1, 1, 1)
  ax.plot(costTheta_src, interp, '-', label=f'interp (log10KE≈{log10KE_val:.3f})',
          linewidth=0.5, color='black', alpha=1.0)
  ax.plot(costTheta_src, orig, '+', label=f'original (log10P={logp_in})',
          markersize=6, color='blue', alpha=0.3)
  ax.set_xlabel('cos(theta)')
  ax.set_ylabel('log10 dFdE')

  title = f"log10P = {logp_in:.2f}  →  log10KE ≈ {log10KE_val:.3f}  |  kx={spline.kx}, ky={spline.ky}, s={spline.s:.2e}"
  ax.set_title(title)
  ax.grid(True, which='both', linestyle=':', color='silver', linewidth=0.8)
  ax.legend()

  Path("figs").mkdir(parents=True, exist_ok=True)
  fig_.savefig(f"figs/fig_logdFdE_logP{logp_in:.2f}.png", dpi=150)
# ----------------------------------------------------------------------
def main():
  opts = parse_args()

  if not opts.infile or not opts.outfile or not opts.json:
    sys.exit("[ERROR] --infile, --outfile, and --json are all required.")

  # read the input table
  costTheta_src, log10P_src, log10KE_src, dFdE_grid = load_flux_table(opts.infile)
  costTheta_tgt, log10KE_tgt = make_target_axes()
  spline = build_interpolator(
    costTheta_src, log10KE_src, dFdE_grid, opts.smoothing)

  # write the interpolated table
  write_output(opts.outfile, spline, costTheta_tgt, log10KE_tgt)

  # read the JSON file
  try:
    with open(opts.json, 'r') as f:
      config = json.load(f)
  except Exception as e:
    sys.exit(f"[ERROR] Failed to read JSON file: {e}")

  # plot the costhz slices
  for cz in config.get("costhz", []):
    if not (0.0 <= cz <= 1.0):
      print(f"[WARN] Skipping invalid costhz = {cz}")
      continue
    plot_costhz_slice(spline, costTheta_src, log10KE_src, dFdE_grid, cz)

  # plot the logP slices
  for lp in config.get("logP", []):
    if not (log10P_src.min() <= lp <= log10P_src.max()):
      print(f"[WARN] Skipping log10P = {lp} (out of range)")
      continue
    plot_logP_slice(spline, costTheta_src, log10P_src, dFdE_grid, lp)

# ----------------------------------------------------------------------
if __name__ == "__main__":
  try:
    main()
  except KeyboardInterrupt:
    sys.exit("[INFO] interrupted by user")
