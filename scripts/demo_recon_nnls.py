#!/usr/bin/env python3
"""Standalone density-reconstruction demo driven only by mat/manifest.json.

This script is the first external consumer of the mat/ export written by
exemdl::pipeline::save (see mat/manifest.json). It never touches the C++
code: everything it knows about the binary layout comes from the manifest.

Pipeline:
  1. Load manifest.json and sanity-check every .bin it describes
     (header rows/cols, file size, row/col table counts, NaN/Inf scan).
  2. Build a synthetic true density x_true (uniform background + one
     low-density block placed via the col table's ix,iy,iz grid indices).
  3. Synthesize observations b_syn = A @ x_true from the center matrix.
  4. Reconstruct x_hat by non-negative regularized least squares:
       minimize ||A x - b_syn||^2 + lam * ||x - x_bg||^2  s.t. x >= 0
     solved with scipy.optimize.lsq_linear on an augmented LinearOperator.
  5. Report residual / relative error and save cross-section figures
     comparing x_true and x_hat.

Usage:
  python3 scripts/demo_recon_nnls.py <mat_dir> [--out DIR] [--lams 0.01,0.1,1]

The loader functions (load_manifest / load_matrix) are import-safe and are
meant to be reused by future prototype reconstruction scripts.
"""

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

HEADER_BYTES = 16  # two uint64: rows, cols
BG_DENSITY = 2000.0  # synthetic background density [kg/m^3]
BLOB_DENSITY = 1000.0  # synthetic low-density block [kg/m^3]


def load_manifest(mat_dir):
  """Load manifest.json found in a mat/ directory and return it as a dict."""
  path = Path(mat_dir) / "manifest.json"
  with open(path, "r") as f:
    return json.load(f)


def load_matrix(path):
  """Open one matrix .bin as a read-only float32 column-major memmap.

  Layout (per manifest 'files' section): uint64 rows, uint64 cols header,
  then rows*cols float32 values in column-major order.
  Returns (memmap, rows, cols).
  """
  path = Path(path)
  header = np.fromfile(path, dtype=np.uint64, count=2)
  rows, cols = int(header[0]), int(header[1])
  mm = np.memmap(path, dtype=np.float32, mode="r", offset=HEADER_BYTES,
                 shape=(rows, cols), order="F")
  return mm, rows, cols


def check_bin_against_manifest(path, meta):
  """Validate one .bin against its manifest entry. Returns list of errors."""
  errors = []
  mm, rows, cols = load_matrix(path)
  if rows != int(meta["rows"]) or cols != int(meta["cols"]):
    errors.append("%s: header (%d,%d) != manifest (%s,%s)"
                  % (path.name, rows, cols, meta["rows"], meta["cols"]))
  expected_size = HEADER_BYTES + rows * cols * 4
  actual_size = path.stat().st_size
  if actual_size != expected_size:
    errors.append("%s: file size %d != expected %d"
                  % (path.name, actual_size, expected_size))
  n_nonfinite = 0
  vmin = np.inf
  vmax = -np.inf
  col_step = max(1, (1 << 24) // max(rows, 1))  # ~16M floats per chunk
  for j0 in range(0, cols, col_step):
    blk = np.asarray(mm[:, j0:j0 + col_step])
    finite = np.isfinite(blk)
    n_nonfinite += int((~finite).sum())
    if finite.any():
      vmin = min(vmin, float(blk[finite].min()))
      vmax = max(vmax, float(blk[finite].max()))
  if n_nonfinite > 0:
    errors.append("%s: %d NaN/Inf values found" % (path.name, n_nonfinite))
  print("  %s: rows=%d cols=%d size_ok nonfinite=%d min=%.6g max=%.6g"
        % (path.name, rows, cols, n_nonfinite, vmin, vmax))
  return errors


def sanity_check(manifest, mat_dir):
  """Run all format checks. Returns list of error strings (empty = PASS)."""
  errors = []
  mat_dir = Path(mat_dir)
  print("[1] Sanity checks against manifest.json")
  for name, meta in sorted(manifest["files"].items()):
    path = mat_dir / name
    if not path.is_file():
      errors.append("%s: listed in manifest but missing on disk" % name)
      continue
    errors.extend(check_bin_against_manifest(path, meta))
  n_rows_tbl = len(manifest["rows"]["entries"])
  n_cols_tbl = len(manifest["cols"]["entries"])
  # The row/col tables describe the grouped dNdD matrices. The per-detector COO
  # and element tables carry their own shapes, so they cannot be checked here.
  dndd_name = next((n for n in sorted(manifest["files"])
                    if n.startswith("mat_dNdD_grouped")), None)
  if dndd_name is None:
    errors.append("manifest lists no mat_dNdD_grouped_*.bin "
                  "to check the row/col tables against")
  else:
    dndd_meta = manifest["files"][dndd_name]
    if n_rows_tbl != int(dndd_meta["rows"]):
      errors.append("rows.entries count %d != %s rows %s"
                    % (n_rows_tbl, dndd_name, dndd_meta["rows"]))
    if n_cols_tbl != int(dndd_meta["cols"]):
      errors.append("cols.entries count %d != %s cols %s"
                    % (n_cols_tbl, dndd_name, dndd_meta["cols"]))
  uqiv_span = int(manifest["cols"]["uqiv_max"]) - int(manifest["cols"]["uqiv_min"]) + 1
  if uqiv_span != n_cols_tbl:
    errors.append("uqiv range %d != cols.entries count %d" % (uqiv_span, n_cols_tbl))
  print("  row table: %d entries, col table: %d entries (uqiv %s..%s)"
        % (n_rows_tbl, n_cols_tbl,
           manifest["cols"]["uqiv_min"], manifest["cols"]["uqiv_max"]))
  print("  detectors (angle_unit):",
        ", ".join("det%s=%s" % (d["detid"], d["angle_unit"])
                  for d in manifest["detectors"]))
  row0 = manifest["rows"]["entries"][0]
  print("  row[0] sample [detid,igroup,uqig,xmin,xmax,ymin,ymax]:", row0)
  return errors


def build_x_true(manifest):
  """Synthetic truth: uniform background + one low-density block.

  The block is placed with the col table's [uqiv, ix, iy, iz] entries so the
  demo exercises the grid mapping documented in the manifest.
  Returns (x_true, ixyz, blob_mask).
  """
  entries = np.asarray(manifest["cols"]["entries"], dtype=np.int64)
  ixyz = entries[:, 1:4]  # columns: ix, iy, iz
  lo = ixyz.min(axis=0)
  hi = ixyz.max(axis=0)
  center = (lo + hi) // 2
  span = np.maximum((hi - lo + 1) // 6, 1)  # block half-width per axis
  blob_mask = np.all(np.abs(ixyz - center) <= span, axis=1)
  x_true = np.full(len(entries), BG_DENSITY, dtype=np.float64)
  x_true[blob_mask] = BLOB_DENSITY
  print("[2] Synthetic truth: grid ix %d..%d, iy %d..%d, iz %d..%d; "
        "block center %s halfwidth %s -> %d of %d voxels at %.0f (bg %.0f)"
        % (lo[0], hi[0], lo[1], hi[1], lo[2], hi[2],
           center.tolist(), span.tolist(), int(blob_mask.sum()),
           len(entries), BLOB_DENSITY, BG_DENSITY))
  return x_true, ixyz, blob_mask


def solve_nnls_regularized(A64, b, x_bg, lam):
  """Solve min ||A x - b||^2 + lam ||x - x_bg||^2 s.t. x >= 0.

  Implemented by stacking sqrt(lam)*I under A as a LinearOperator so no
  augmented dense matrix is ever materialized.
  """
  from scipy.optimize import lsq_linear
  from scipy.sparse.linalg import LinearOperator

  m, n = A64.shape
  sq = np.sqrt(lam)

  def matvec(x):
    return np.concatenate([A64 @ x, sq * x])

  def rmatvec(y):
    return A64.T @ y[:m] + sq * y[m:]

  op = LinearOperator((m + n, n), matvec=matvec, rmatvec=rmatvec,
                      dtype=np.float64)
  b_aug = np.concatenate([b, sq * x_bg])
  t0 = time.time()
  res = lsq_linear(op, b_aug, bounds=(0.0, np.inf), method="trf",
                   lsmr_tol="auto", max_iter=30, verbose=0)
  print("    lam=%.3g: lsq_linear %s in %.1fs (%d iterations)"
        % (lam, "converged" if res.success else "stopped",
           time.time() - t0, res.nit))
  return res.x


def _slice_figure(ixyz, x_true, x_hat, iz):
  """Build one side-by-side (truth | reconstruction) figure for slice iz."""
  import matplotlib.pyplot as plt

  lo = ixyz.min(axis=0)
  shape = ixyz.max(axis=0) - lo + 1
  img_t = np.full((shape[1], shape[0]), np.nan)
  img_r = np.full((shape[1], shape[0]), np.nan)
  sel = ixyz[:, 2] == iz
  xs = ixyz[sel, 0] - lo[0]
  ys = ixyz[sel, 1] - lo[1]
  img_t[ys, xs] = x_true[sel]
  img_r[ys, xs] = x_hat[sel]
  fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
  for ax, img, title in ((axes[0], img_t, "x_true (original)"),
                         (axes[1], img_r, "x_hat (reconstruction)")):
    im = ax.imshow(img, origin="lower", vmin=0.0, vmax=1.25 * BG_DENSITY)
    ax.set_title("%s, iz=%d" % (title, iz))
    ax.set_xlabel("ix")
    ax.set_ylabel("iy")
    fig.colorbar(im, ax=ax, label="density [kg/m^3]")
  return fig


def save_figures(out_dir, ixyz, x_true, x_hat, blob_mask):
  """Save a PNG at the block slice and a GIF scanning all z slices."""
  import matplotlib
  matplotlib.use("Agg")
  import matplotlib.pyplot as plt
  from PIL import Image

  out_dir = Path(out_dir)
  out_dir.mkdir(parents=True, exist_ok=True)
  paths = []
  iz_blob = int(np.round(ixyz[blob_mask][:, 2].mean()))
  fig = _slice_figure(ixyz, x_true, x_hat, iz_blob)
  png_path = out_dir / ("demo_recon_nnls_iz%03d.png" % iz_blob)
  fig.savefig(png_path, dpi=130)
  plt.close(fig)
  paths.append(png_path)
  print("  figure saved:", png_path)

  frames = []
  for iz in range(int(ixyz[:, 2].min()), int(ixyz[:, 2].max()) + 1):
    fig = _slice_figure(ixyz, x_true, x_hat, iz)
    fig.canvas.draw()
    buf = np.asarray(fig.canvas.buffer_rgba())[:, :, :3]
    frames.append(Image.fromarray(buf.copy()))
    plt.close(fig)
  gif_path = out_dir / "demo_recon_nnls_zscan.gif"
  frames[0].save(gif_path, save_all=True, append_images=frames[1:],
                 duration=500, loop=0)
  paths.append(gif_path)
  print("  figure saved:", gif_path)
  return paths


def main(argv=None):
  parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
  parser.add_argument("mat_dir", help="mat/ directory containing manifest.json")
  parser.add_argument("--out", default=None,
                      help="output directory (default: <mat_dir>/../demo_recon_nnls)")
  parser.add_argument("--lams", default="0.01,0.1,1.0",
                      help="comma-separated regularization weights to try")
  args = parser.parse_args(argv)

  mat_dir = Path(args.mat_dir)
  out_dir = Path(args.out) if args.out else mat_dir.parent / "demo_recon_nnls"
  manifest = load_manifest(mat_dir)

  errors = sanity_check(manifest, mat_dir)
  if errors:
    print("SANITY CHECK FAILED:")
    for e in errors:
      print("  -", e)
    return 1
  print("  sanity checks: PASS")

  x_true, ixyz, blob_mask = build_x_true(manifest)

  mm, rows, cols = load_matrix(mat_dir / "mat_dNdD_grouped_center.bin")
  A64 = np.asarray(mm, dtype=np.float64)
  b_syn = A64 @ x_true
  print("[3] b_syn = A @ x_true: %d observations, |b| range %.4g..%.4g"
        % (rows, float(np.abs(b_syn).min()), float(np.abs(b_syn).max())))

  x_bg = np.full(cols, BG_DENSITY, dtype=np.float64)
  col_scale = float(np.mean(np.sum(A64 * A64, axis=0)))  # mean squared col norm
  print("[4] Reconstruction (lam scaled by mean column norm^2 = %.4g)" % col_scale)
  best = None
  for lam_rel in (float(s) for s in args.lams.split(",")):
    lam = lam_rel * col_scale
    x_hat = solve_nnls_regularized(A64, b_syn, x_bg, lam)
    resid = float(np.linalg.norm(A64 @ x_hat - b_syn) / np.linalg.norm(b_syn))
    relerr = float(np.linalg.norm(x_hat - x_true) / np.linalg.norm(x_true))
    blob_mean = float(x_hat[blob_mask].mean())
    bg_mean = float(x_hat[~blob_mask].mean())
    print("    lam_rel=%.3g: residual=%.3e rel_err=%.3f "
          "blob_mean=%.0f bg_mean=%.0f" % (lam_rel, resid, relerr,
                                           blob_mean, bg_mean))
    if best is None or relerr < best[1]:
      best = (lam_rel, relerr, x_hat)
  lam_rel, relerr, x_hat = best
  print("[5] Best lam_rel=%.3g (rel_err=%.3f); saving figures" % (lam_rel, relerr))
  save_figures(out_dir, ixyz, x_true, x_hat, blob_mask)
  print("DEMO PASS: manifest-driven load, synthetic reconstruction completed.")
  return 0


if __name__ == "__main__":
  sys.exit(main())
