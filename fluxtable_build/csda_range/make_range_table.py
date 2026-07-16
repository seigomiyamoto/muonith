#!/usr/bin/env python3
"""Build a muon CSDA range table from a PDG muE_<material>.txt file.

Input : PDG muon energy-loss table (https://pdg.lbl.gov/2024/AtomicNuclearProperties/MUE/),
        column 1 = kinetic energy T [MeV], column 9 = CSDA range [g/cm^2].
        Non-numeric header/comment lines are skipped.
Method: node-exact cubic spline in log10(KE [GeV]) vs log10(range [g/cm^2]),
        sampled on a uniform log10(KE) grid, then converted to kg/m^2 (+1 in log10).
Output: two-column "%E %E" ASCII, identical in format to
        fluxtable/daemon_groom/logRangeKgm2_logGeV_hokan_70001.txt
        (column 1 = log10(KE [GeV]), column 2 = log10(CSDA range [kg/m^2])).

The default grid (-2 .. 5, step 1e-4) gives 70001 rows. The lower bound of 10 MeV
follows the PDG file header warning that results below 10 MeV are not dependable.
"""

import argparse
import sys

import numpy as np
from scipy.interpolate import CubicSpline


def read_pdg_muE(path):
  """Return (T_MeV, csda_gcm2) arrays parsed from a PDG muE table."""
  rows = []
  with open(path) as f:
    for line in f:
      tokens = line.split()
      if len(tokens) < 9:
        continue
      try:
        rows.append((float(tokens[0]), float(tokens[8])))
      except ValueError:
        continue
  if not rows:
    raise RuntimeError(f"no numeric (T, CSDA) rows found in {path}")
  arr = np.array(rows)
  # Keep the first occurrence of each T (marker lines can duplicate energies).
  _, idx = np.unique(arr[:, 0], return_index=True)
  arr = arr[np.sort(idx)]
  return arr[:, 0], arr[:, 1]


def build_table(t_mev, csda_gcm2, xmin, xmax, step):
  """Return (x, z): log10(KE GeV) grid and log10(range kg/m^2) values."""
  lx = np.log10(t_mev / 1000.0)
  lz = np.log10(csda_gcm2)
  if lx[0] > xmin or lx[-1] < xmax:
    raise RuntimeError(
      f"input covers log10KE {lx[0]:.3f}..{lx[-1]:.3f}, "
      f"cannot sample {xmin}..{xmax} without extrapolation")
  n = int(round((xmax - xmin) / step)) + 1
  x = xmin + step * np.arange(n)
  spline = CubicSpline(lx, lz)
  z = spline(x) + 1.0  # log10(g/cm^2) -> log10(kg/m^2)
  return x, z


def compare(x, z, ref_path):
  """Print max deviations against a reference table on the same grid."""
  ref = np.loadtxt(ref_path)
  if len(ref) != len(x) or abs(ref[0, 0] - x[0]) > 1e-9:
    print(f"compare: grid mismatch with {ref_path} "
          f"({len(ref)} rows vs {len(x)}); skipping", file=sys.stderr)
    return
  d = np.abs(z - ref[:, 1])
  above_1gev = x >= 0.0
  print(f"compare vs {ref_path}")
  print(f"  max|dz| overall      : {d.max():.3e} dex (at log10KE={x[np.argmax(d)]:+.4f})")
  print(f"  max|dz| log10KE >= 0 : {d[above_1gev].max():.3e} dex")
  print(f"  |dz| at both ends    : {d[0]:.3e} / {d[-1]:.3e} dex")


def main():
  ap = argparse.ArgumentParser(description=__doc__,
                               formatter_class=argparse.RawDescriptionHelpFormatter)
  ap.add_argument("input", help="PDG muE_<material>.txt")
  ap.add_argument("output", help="output range table (2-column ASCII)")
  ap.add_argument("--xmin", type=float, default=-2.0, help="grid start in log10(KE GeV)")
  ap.add_argument("--xmax", type=float, default=5.0, help="grid end in log10(KE GeV)")
  ap.add_argument("--step", type=float, default=1.0e-4, help="grid step in log10(KE GeV)")
  ap.add_argument("--compare", metavar="REF",
                  help="reference table to report max deviations against")
  args = ap.parse_args()

  t_mev, csda_gcm2 = read_pdg_muE(args.input)
  print(f"read {len(t_mev)} points from {args.input} "
        f"(T = {t_mev[0]:.3g} .. {t_mev[-1]:.3g} MeV)")
  x, z = build_table(t_mev, csda_gcm2, args.xmin, args.xmax, args.step)

  with open(args.output, "w") as f:
    for xi, zi in zip(x, z):
      f.write(f"{xi:E} {zi:E}\n")
  print(f"wrote {len(x)} rows to {args.output}")

  if args.compare:
    compare(x, z, args.compare)


if __name__ == "__main__":
  main()
