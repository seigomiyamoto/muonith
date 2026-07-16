#!/usr/bin/env python3
"""Generate an eff_table for a detector from a JSON5 efficiency model.

The central efficiency of each angular bin is built as a PRODUCT of two
per-axis fit functions (tx, ty), times a base plateau efficiency:

    eff_cnt = base_eff * g_tx(tx) * g_ty(ty) * (1 + noise)

The optional (1 + noise) is a reproducible random
factor from the 'randomize' block, the product of three independent components:
a per-tx-bin draw (acts on g_tx), a per-ty-bin draw (acts on g_ty), and a
per-(tx,ty)-bin draw (amp), each (1 + Normal(0, value)) seeded for repeatability.
It also applies to the k-of-n layer model, scaling eff_low/cnt/upp together so
the band moves with the scattered center. The absolute uncertainty has its OWN dedicated
per-axis functions, independent of the central efficiency:

    sigma_eff = sigma.base * 0.5*(h_tx(tx) + h_ty(ty))
    half      = min(sigma_eff, eff_cnt, 1 - eff_cnt)   (symmetric, kept in [0,1])
    eff_low   = eff_cnt - half     so 0.5*(eff_upp + eff_low) == eff_cnt and
    eff_upp   = eff_cnt + half     0 <= eff_low <= eff_upp <= 1

This lets the uncertainty depend on x/y/r however the detector demands (e.g.
grow at large angle where the acceptance tail is poorly known) instead of
being tied to the central value. The band feeds the existing C_N diagonal
common mode (DetectorPanelArray::get_vecxf_var_eff_all, tf_independent=false),
which reads sigma_eff = 0.5*(eff_upp - eff_low) per bin.

The angular grid (the 4 columns xlow xup ylow yup) is NOT re-derived here:
it is read from an existing eff_table (default: the sample00 template) so it
always matches the consumer geometry. Only the 3 efficiency columns are filled.

Detector-type differences are expressed purely in the JSON5 config; no code
change is needed to add a new detector. See make_eff_table_default.json5.

Usage:
  python3 make_eff_table.py --config <model.json5> --out <eff_table.txt>
                            [--grid <existing_eff_table>] [--dry-run]
"""

import argparse
import math
import sys
from pathlib import Path

import numpy as np

try:
  import json5
except ImportError:  # pragma: no cover
  sys.stderr.write(
    "json5 module not found. Install with: pip install json5\n")
  raise

DEFAULT_CONFIG_FILE = Path(__file__).parent / "make_eff_table_default.json5"
DEFAULT_GRID_FILE = Path(__file__).parent / "templates" / "eff_table_sample00.tmp"


def load_config(path):
  """Load a JSON5 efficiency-model config."""
  with open(path, "r") as f:
    return json5.loads(f.read())


def read_grid(path):
  """Read the angular grid (first 4 columns) from an existing eff_table.

  Returns (xlow, xup, ylow, yup) as float arrays. Bin centers are derived
  by the caller. Rows are kept in file order so the output row order matches
  the input grid exactly.
  """
  data = np.loadtxt(path, usecols=(0, 1, 2, 3), dtype=float)
  if data.ndim != 2 or data.shape[1] != 4:
    raise ValueError(
      f"grid file must have >=4 numeric columns: {path}")
  return data[:, 0], data[:, 1], data[:, 2], data[:, 3]


# --- per-axis fit functions -------------------------------------------------
# Each takes the coordinate array x and the params dict, returns a factor in
# [0, inf) (clipped to [0,1] only after the full product is formed).

def _fn_flat(x, params):
  """Constant 1 (no dependence on this axis)."""
  return np.ones_like(x, dtype=float)


def _fn_poly(x, params):
  """Polynomial sum_k coef[k] * x^k evaluated at x."""
  coef = params.get("coef", [1.0])
  out = np.zeros_like(x, dtype=float)
  for k, c in enumerate(coef):
    out += c * np.power(x, k)
  return out


def _fn_gauss(x, params):
  """Gaussian exp(-((x-mu)^2)/(2 sigma^2))."""
  mu = float(params.get("mu", 0.0))
  sigma = float(params.get("sigma", 1.0))
  if sigma == 0.0:
    raise ValueError("gauss fn requires non-zero sigma")
  return np.exp(-((x - mu) ** 2) / (2.0 * sigma * sigma))


def _fn_cos_power(x, params):
  """cos(theta)^n with theta = atan(|x|) (x interpreted as tan(theta)).

  Optional cutoff_deg zeroes the factor beyond that zenith angle.
  """
  n = float(params.get("n", 1.0))
  theta = np.arctan(np.abs(x))
  out = np.power(np.cos(theta), n)
  cutoff_deg = params.get("cutoff_deg", None)
  if cutoff_deg is not None:
    out = np.where(theta > math.radians(float(cutoff_deg)), 0.0, out)
  return out


def _fn_sigmoid(x, params):
  """Plateau-then-drop sigmoid 1/(1+exp((|x|-x0)/w)).

  ~1 for |x| < x0, drops past x0 over width w (w>0). Models an HV/acceptance
  plateau with a knee at x0.
  """
  x0 = float(params.get("x0", 1.0))
  w = float(params.get("w", 0.1))
  if w == 0.0:
    raise ValueError("sigmoid fn requires non-zero w")
  return 1.0 / (1.0 + np.exp((np.abs(x) - x0) / w))


def _fn_step(x, params):
  """Two-level step: inner for |x|<x0, outer for |x|>x0, smooth over width w.

  f(x) = outer + (inner - outer) / (1 + exp((|x| - x0)/w)).
  Models a central plateau (inner) that transitions to a floor/shelf (outer)
  past |x|=x0. Use inner>outer for an efficiency plateau-with-floor, or
  inner<outer for an uncertainty that grows outside the plateau.
  """
  x0 = float(params.get("x0", 1.0))
  w = float(params.get("w", 0.1))
  inner = float(params.get("inner", 1.0))
  outer = float(params.get("outer", 0.0))
  if w == 0.0:
    raise ValueError("step fn requires non-zero w")
  return outer + (inner - outer) / (1.0 + np.exp((np.abs(x) - x0) / w))


_FIT_FNS = {
  "flat": _fn_flat,
  "poly": _fn_poly,
  "gauss": _fn_gauss,
  "cos_power": _fn_cos_power,
  "sigmoid": _fn_sigmoid,
  "step": _fn_step,
}


def apply_axis(x, spec):
  """Evaluate one axis factor from a {fn, params} spec."""
  if spec is None:
    return np.ones_like(x, dtype=float)
  fn = spec.get("fn", "flat")
  if fn not in _FIT_FNS:
    raise ValueError(
      f"unknown fit fn '{fn}'. choose from {sorted(_FIT_FNS)}")
  return _FIT_FNS[fn](x, spec.get("params", {}))


def _kofn_tail(p, n, k):
  """P(X >= k) for X ~ Binomial(n, p), vectorized over the array p.

  This is the probability that at least k of n layers fire given a per-layer
  efficiency p. Monotone increasing in p, so a band on p maps to a band here.
  """
  p = np.clip(np.asarray(p, dtype=float), 0.0, 1.0)
  q = 1.0 - p
  out = np.zeros_like(p)
  for j in range(int(k), int(n) + 1):
    out += math.comb(int(n), j) * np.power(p, j) * np.power(q, int(n) - j)
  return out


def _build_randomize_factor(rz, txc, tyc):
  """Build a reproducible multiplicative factor on the central efficiency.

  Three independent random factors, each (1 + Normal(0, value)), multiply
  together so the model can scatter along tx, along ty, and per (tx,ty) bin:

    tx  : one draw per tx bin  -> shared down each tx column (acts on g_tx)
    ty  : one draw per ty bin  -> shared across each ty row  (acts on g_ty)
    amp : one draw per (tx,ty) bin -> independent per-cell scatter

  Each component is disabled when its value is <= 0. Distinct seed offsets keep
  the three draws independent yet reproducible. Returns an all-ones array (same
  length as txc) when nothing is enabled.
  """
  rz = rz or {}
  seed = int(rz.get("seed", 0))
  n = txc.size
  factor = np.ones(n, dtype=float)

  a_tx = float(rz.get("tx", 0.0))
  if a_tx > 0.0:
    _, inv = np.unique(np.round(txc, 6), return_inverse=True)
    draws = 1.0 + np.random.default_rng(seed).normal(0.0, a_tx, size=inv.max() + 1)
    factor = factor * draws[inv]

  a_ty = float(rz.get("ty", 0.0))
  if a_ty > 0.0:
    _, inv = np.unique(np.round(tyc, 6), return_inverse=True)
    draws = 1.0 + np.random.default_rng(seed + 1).normal(0.0, a_ty, size=inv.max() + 1)
    factor = factor * draws[inv]

  a_xy = float(rz.get("amp", 0.0))
  if a_xy > 0.0:
    factor = factor * (1.0 + np.random.default_rng(seed + 2).normal(0.0, a_xy, size=n))

  return factor


def build_eff_layers(layers, r, factor=None):
  """Build (eff_low, eff_cnt, eff_upp) for a k-of-n layer-coincidence model.

  The track-finding efficiency is P(>= k of n layers fire) given a per-layer
  efficiency p(r). The band comes from the per-layer uncertainty per_layer_sigma
  propagated through the binomial tail: eff_low/upp = tail(p -/+ sigma_p).

  If a randomize factor array is given, it multiplies eff_low/cnt/upp together,
  so the central value scatters bin-to-bin (incl. tx/ty structure) while the
  band moves with it and the ordering eff_low <= eff_cnt <= eff_upp is kept.
  """
  n = int(layers.get("n", 10))
  k = int(layers.get("k", 8))
  sig_p = float(layers.get("per_layer_sigma", 0.0))
  p = np.clip(apply_axis(r, layers.get("per_layer")), 0.0, 1.0)
  eff_cnt = _kofn_tail(p, n, k)
  eff_low = _kofn_tail(np.clip(p - sig_p, 0.0, 1.0), n, k)
  eff_upp = _kofn_tail(np.clip(p + sig_p, 0.0, 1.0), n, k)

  if factor is not None:
    eff_low = np.clip(eff_low * factor, 0.0, 1.0)
    eff_cnt = np.clip(eff_cnt * factor, 0.0, 1.0)
    eff_upp = np.clip(eff_upp * factor, 0.0, 1.0)
  return eff_low, eff_cnt, eff_upp


def build_sigma(cfg, txc, tyc):
  """Build the absolute per-bin uncertainty sigma_eff from its own functions.

  sigma_eff = sigma.base * 0.5*(h_tx(tx) + h_ty(ty)), the MEAN of dedicated tx /
  ty fit functions independent of the central efficiency (no radial axis). The
  mean (not product) lets sigma vanish at tx=ty=0 while still rising along each
  axis: with h(t)=c*t^2, sigma is 0 at the origin yet > 0 on either axis. With
  flat factors (h=1) it reduces to sigma.base, a uniform band.
  Returns a zero array when no 'sigma' block is present.
  """
  sig_cfg = cfg.get("sigma", None)
  if not sig_cfg:
    return np.zeros_like(txc, dtype=float)
  base = float(sig_cfg.get("base", 0.0))
  h_tx = apply_axis(txc, sig_cfg.get("tx"))
  h_ty = apply_axis(tyc, sig_cfg.get("ty"))
  sigma_eff = base * 0.5 * (h_tx + h_ty)  # MEAN over tx, ty (flat -> base)
  return np.clip(sigma_eff, 0.0, None)


def build_eff(cfg, xlow, xup, ylow, yup):
  """Build (eff_low, eff_cnt, eff_upp) arrays from the model config."""
  txc = 0.5 * (xlow + xup)  # bin center, tx axis
  tyc = 0.5 * (ylow + yup)  # bin center, ty axis
  r = np.sqrt(txc * txc + tyc * tyc)  # radial coordinate (tan units)

  # reproducible random factor: per-tx, per-ty, and per-(tx,ty) bin (see
  # _build_randomize_factor). All-ones when the randomize block is absent/zero.
  factor = _build_randomize_factor(cfg.get("randomize"), txc, tyc)

  # k-of-n layer-coincidence model (radial only): track-finding efficiency
  layers = cfg.get("layers")
  if layers:
    return build_eff_layers(layers, r, factor)

  base = float(cfg.get("base_eff", 1.0))
  g_tx = apply_axis(txc, cfg.get("tx"))
  g_ty = apply_axis(tyc, cfg.get("ty"))

  eff_cnt = base * g_tx * g_ty  # PRODUCT over tx, ty (no radial axis)

  # reproducible random factor (per-tx / per-ty / per-(tx,ty) bin)
  eff_cnt = eff_cnt * factor

  eff_cnt = np.clip(eff_cnt, 0.0, 1.0)

  # band = absolute uncertainty from its own dedicated x/y/r functions.
  # Kept SYMMETRIC about eff_cnt (NOT clipped to [0,1]) so that, per bin,
  #   0.5*(eff_upp + eff_low) == eff_cnt   and   0.5*(eff_upp - eff_low) == sigma_eff
  # both hold exactly. Near the bounds eff_upp may exceed 1 (or eff_low drop
  # below 0); those endpoints encode the band, not a literal efficiency -- the
  # central value eff_cnt itself stays in [0,1].
  sigma_eff = build_sigma(cfg, txc, tyc)
  # symmetric half-width, shrunk where the full sigma would leave [0,1], so that
  #   0.5*(eff_upp + eff_low) == eff_cnt   (center always mid-band)   AND
  #   0 <= eff_low <= eff_upp <= 1         (endpoints stay valid efficiencies).
  # The cost: where eff_cnt approaches 1 (or 0) the band narrows toward 0, i.e.
  # the reported uncertainty shrinks to the available head-room at the bound.
  half = np.minimum(sigma_eff, np.minimum(eff_cnt, 1.0 - eff_cnt))
  eff_low = eff_cnt - half
  eff_upp = eff_cnt + half
  return eff_low, eff_cnt, eff_upp


def write_table(path, xlow, xup, ylow, yup, eff_low, eff_cnt, eff_upp):
  """Write the 7-column eff_table, matching the existing fixed-width style."""
  cols = np.column_stack([xlow, xup, ylow, yup, eff_low, eff_cnt, eff_upp])
  with open(path, "w") as f:
    for row in cols:
      f.write(" ".join(f"{v:7.4f}" for v in row) + "\n")


def main():
  parser = argparse.ArgumentParser(
    description="Generate an eff_table from a JSON5 efficiency model.")
  parser.add_argument("--config", default=str(DEFAULT_CONFIG_FILE),
                      help="JSON5 efficiency-model config")
  parser.add_argument("--grid", default=str(DEFAULT_GRID_FILE),
                      help="existing eff_table to read the angular grid from")
  parser.add_argument("--out", required=True, help="output eff_table path")
  parser.add_argument("--dry-run", action="store_true",
                      help="compute and print stats without writing")
  args = parser.parse_args()

  cfg = load_config(args.config)
  xlow, xup, ylow, yup = read_grid(args.grid)
  eff_low, eff_cnt, eff_upp = build_eff(cfg, xlow, xup, ylow, yup)

  n = eff_cnt.size
  print(f"rows={n}  eff_cnt min/mean/max="
        f"{eff_cnt.min():.4f}/{eff_cnt.mean():.4f}/{eff_cnt.max():.4f}  "
        f"band half-width mean={0.5 * (eff_upp - eff_low).mean():.4f}")
  if not (np.all(eff_low <= eff_cnt) and np.all(eff_cnt <= eff_upp)):
    raise ValueError("ordering violated: need eff_low <= eff_cnt <= eff_upp")

  if args.dry_run:
    print("dry-run: not writing")
    return
  write_table(args.out, xlow, xup, ylow, yup, eff_low, eff_cnt, eff_upp)
  print(f"wrote {args.out}")


if __name__ == "__main__":
  main()
