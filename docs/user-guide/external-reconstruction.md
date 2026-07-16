# External Reconstruction

The inversion writes its observation matrix and voxel/observation bookkeeping to a
self-describing `mat/` directory. Everything an external tool needs to run its own
density reconstruction — the matrix layout, the meaning of every row and column, the
units — is recorded in `mat/manifest.json`. No access to the C++ code is required.

This page documents the `mat/` export and a worked demo,
[`demo_recon_nnls.py`](../reference/scripts.md#demo_recon_nnlspy), that reconstructs
density from `mat/` alone.

## The `mat/` export

The `mat/` directory is written under a run's output directory (the same directory that
holds the binary checkpoints — see [Output](output.md#binary-checkpoints)). It contains
the observation matrix in binary form plus a single `manifest.json` that describes it.

### `manifest.json` structure

`manifest.json` has five top-level keys:

| Key | Content |
|-----|---------|
| `format_version` | Manifest schema version (currently `1`). |
| `files` | One entry per binary file, keyed by file name (see [Binary layout](#binary-layout)). |
| `detectors` | Per-detector metadata: `detid`, `angle_unit`, `pos_m` (position, meters), `dir`, `n_element`. |
| `rows` | Row contract of the matrix (see [Rows and columns](#rows-and-columns)). |
| `cols` | Column contract of the matrix (see [Rows and columns](#rows-and-columns)). |

The observation matrix is exported as three prior brackets:

| File | Meaning |
|------|---------|
| `mat_dNdD_grouped_lower.bin`  | Observation matrix $\mathbf{A} = \partial N / \partial \rho$ at the lower prior bracket. |
| `mat_dNdD_grouped_center.bin` | Observation matrix at the central prior. This is the matrix used by the demo. |
| `mat_dNdD_grouped_upper.bin`  | Observation matrix at the upper prior bracket. |

### Binary layout

Every matrix `.bin` file uses the same layout, also stated verbatim in each `files`
entry (`dtype`, `order`, `header`, `rows`, `cols`):

1. **Header** — two `uint64` values in native endianness: the row count, then the
   column count.
2. **Data** — `rows × cols` `float32` values in **column-major** order (Eigen's default).

To read one matrix in Python: read the 16-byte header as two `uint64`, then read the
remaining `rows × cols` `float32` values column-major. The demo's `load_matrix` helper
does exactly this and is safe to import and reuse.

### Rows and columns

The matrix follows a fixed contract, so an external solver can map every row to a
physical observation and every column to a voxel:

- **Row $i$ = one observation bin** (a single detector's angle group). `rows.columns`
  lists the fields `[detid, igroup, uqig, xmin, xmax, ymin, ymax]`, and `rows.entries`
  gives one entry per row in matrix order.
- **Column $j$ = one voxel**. `cols.columns` lists `[uqiv, ix, iy, iz]`, and the voxel
  index satisfies `uqiv = uqiv_min + j` (`cols.uqiv_min` / `cols.uqiv_max` bound the
  range). `cols.entries` gives one entry per column in matrix order.

## What the matrix represents

The exported matrix $\mathbf{A}$ is the observation matrix $\partial N / \partial \rho$
of the standard MAP (maximum a posteriori) reconstruction. The full estimator, including
the covariance terms, is described in [Inversion](../concepts/inversion.md#map-estimation):

$$\vec{\rho}' = \vec{\rho}_0 + \left(\mathbf{A}^T \mathbf{C}_N^{-1} \mathbf{A} + \mathbf{C}_\rho^{-1}\right)^{-1} \mathbf{A}^T \mathbf{C}_N^{-1} \left(\vec{N}_{\text{obs}} - \vec{N}_0\right)$$

Here $\mathbf{A}$ maps a voxel density vector to expected muon counts, $\mathbf{C}_N$ is
the observation (count) covariance, and $\mathbf{C}_\rho$ is the prior density
covariance. An external tool is free to ignore this estimator entirely and apply any
other method to the same $\mathbf{A}$ and observations.

## Demo: `demo_recon_nnls.py`

[`demo_recon_nnls.py`](../reference/scripts.md#demo_recon_nnlspy) is a standalone,
NagaInv-free reconstruction driven only by `mat/manifest.json`. It demonstrates that the
export is enough to run an alternative solver — here, non-negative regularized least
squares via SciPy.

```bash
python3 scripts/demo_recon_nnls.py <mat_dir> [--out DIR] [--lams 0.01,0.1,1.0]
```

The script:

1. Loads `mat/manifest.json` and sanity-checks every `.bin` it lists (header shape
   against the manifest, file size, and a scan for `NaN`/`Inf`).
2. Builds a synthetic ground truth from the column table: a uniform background density
   with one low-density block placed by voxel index.
3. Synthesizes observed counts $\vec{b} = \mathbf{A}\,\vec{x}_{\text{true}}$ using the
   central matrix.
4. Solves non-negative regularized least squares
   $\min_{\vec{x} \ge 0} \lVert \mathbf{A}\vec{x} - \vec{b} \rVert^2 + \lambda \lVert \vec{x} - \vec{x}_{\text{bg}} \rVert^2$
   for each $\lambda$ in `--lams`, and keeps the result with the smallest relative error.

It writes a PNG of the reconstructed slice through the block and a GIF scanning all $z$
slices to `--out` (default: `<mat_dir>/../demo_recon_nnls`). It requires NumPy, SciPy,
and Matplotlib.

Concrete example (using the `mat/` directory produced by a run):

```bash
python3 scripts/demo_recon_nnls.py path/to/run_output/mat
```
