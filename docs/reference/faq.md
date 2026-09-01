# FAQ

Frequently asked questions and troubleshooting tips.

## Installation and Building

### Q: CMake cannot find Eigen3

Eigen3 is downloaded automatically via FetchContent. If you see this error,
ensure you are using the current `CMakeLists.txt` (which includes
FetchContent declarations for Eigen3).

See [Installation](../getting-started/installation.md) for full dependency setup.

### Q: OpenBLAS linking errors (Linux)

Ensure OpenBLAS is installed:

```bash
sudo apt install -y libopenblas-openmp-dev liblapacke-dev
```

!!! note "macOS uses Accelerate by default"
    On macOS, CMake links against Apple's **Accelerate framework** automatically.
    OpenBLAS is not required. If you intentionally want to use Homebrew OpenBLAS,
    set `CMAKE_PREFIX_PATH`:

    ```bash
    export CMAKE_PREFIX_PATH="$(brew --prefix openblas):${CMAKE_PREFIX_PATH}"
    ```

### Q: Compiler version too old

MUONITH requires a C++20 compiler: GCC 13+ or Apple Clang 17+.

=== "Ubuntu 22.04"

    ```bash
    sudo add-apt-repository ppa:ubuntu-toolchain-r/test
    sudo apt update
    sudo apt install g++-13 gcc-13
    ```

=== "Ubuntu 24.04"

    G++ 13 is included by default.

=== "macOS"

    Apple Clang 17+ ships with the Xcode Command Line Tools and is what
    `build.sh` uses:

    ```bash
    xcode-select --install
    ```

    Or use the [Nix development environment](../getting-started/nix-setup.md) which provides
    the compiler automatically.

### Q: Nix shell does not activate / `nix develop` fails

Ensure Nix flakes are enabled. With the Determinate Systems installer, flakes
are enabled by default. Otherwise, add to `~/.config/nix/nix.conf`:

```
experimental-features = nix-command flakes
```

See [Nix Setup](../getting-started/nix-setup.md) for full instructions.

### Q: What is `IN_NIX_SHELL`?

An environment variable set to `1` by the Nix devShell. `build.sh` checks it
to decide the compiler selection strategy: in a Nix shell, it uses the
Nix-provided `cc`/`c++` wrappers instead of the platform default.

## Running

### Q: How do I set the number of threads?

Set the OpenMP thread count via environment variable before running:

```bash
export OMP_NUM_THREADS=16
```

MUONITH internally sets Eigen threads to `OMP_NUM_THREADS / 2` to balance parallelism between ray tracing and linear algebra.

### Q: Can I re-run just the inversion without recomputing path lengths?

Yes. Save the Trace Path Lengths (Module 4) results by setting `tf_save_arrdet_g3vox: true` in `PATH_LENGTH_PARAMETERS`. On subsequent runs, set `tf_load_arrdet_g3vox: true` so Trace Path Lengths (Module 4) reuses the saved path lengths, then let Modules 5–8 run.

Parameter sweeps on inversion-side knobs (`sigma_rho`, `corr_length`, etc.) can reuse a saved `checkpoint_trace_path_lengths` via `--resume <checkpoint_dir>`, which skips Modules 1–4 entirely.

### Q: How do I stop the pipeline early?

Use `--end-stage <N>` (or set `end_stage` in JSON5) to stop after Module N. Valid values are 3, 4, 5, 6, 7, 8 (default 8). For example, `--end-stage 4` stops right after the ray-tracing step, which is convenient for building a `checkpoint_trace_path_lengths` to resume from later. See [CLI](cli.md#muonithexe) for the full option list.

### Q: How long does a typical run take?

Runtime depends heavily on grid size (number of voxels), number of detectors, and angular bin resolution. Trace Path Lengths (Module 4) (ray tracing / path length computation) is typically the bottleneck. Saving and loading binary checkpoints (`tf_save_arrdet_g3vox`, `tf_load_arrdet_g3vox`) avoids recomputation.

Note that the dense covariance matrix `mat_cov_dens` requires `num_voxels² × 4` bytes (float). For example, 50,000 voxels consume approximately 9.3 GB of memory for this matrix alone.

## Parameters

### Q: What values should I use for `sigma_rho` and `corr_length`?

These are the most important regularization parameters:

- **`sigma_rho`** (density prior standard deviation): Controls the trade-off between data fit and prior. Typical range: 100–1000 kg/m³.
    - Too small → reconstruction stays close to the prior
    - Too large → reconstruction overfits the data (noisy)

- **`corr_length`** (spatial correlation length): Controls spatial smoothness. Should be comparable to or larger than the voxel size. Typical range: 30–200 m.
    - Too small → noisy, checkerboard-like artifacts
    - Too large → over-smoothed, anomalies are blurred

Use `NAGAINV_PARAM_SWEEP` to systematically explore combinations.

### Q: What is `n_hit_det_min`?

The minimum number of distinct detector panels whose rays must pass through a voxel for it to be included in the inversion. Setting `n_hit_det_min = 3` means voxels seen by fewer than 3 detectors are excluded — they have insufficient angular coverage for reliable reconstruction.

### Q: What is voxel merging?

Merging combines adjacent voxels to reduce the grid resolution. For example, `merge_factor = 4` combines 4×4×4 = 64 voxels into one. This reduces the number of unknowns and the size of the covariance matrix, making the inversion faster and more stable.

The `merge_center` parameters control the alignment of the merged grid.

### Q: What does `tf_build_g3vox = false` do?

It skips the 3D voxel construction entirely. The pipeline stops after Build Geometry (Module 3), producing only 2D projected density maps from the DEM. Useful for quick single-view analysis without 3D reconstruction.

## Input Data

### Q: What coordinate system should my DEM use?

MUONITH uses a right-handed, z-up coordinate system in meters. Any projected coordinate system (e.g., UTM, local plane coordinates) works. Geographic coordinates (lat/lon) must be converted first.

If your coordinates have large absolute values (common with UTM), enable `tf_shift_x` and `tf_shift_y` to shift the origin.

### Q: What DEM resolution do I need?

The DEM resolution determines the horizontal resolution of the terrain model. Typical choices:

| Application | DEM resolution |
|-------------|---------------|
| Volcanic edifice (~2 km) | 5–10 m |
| Archaeological structure (~100 m) | 0.5–2 m |
| Underground cavity (~50 m) | 1–5 m |

The DEM should cover the full extent of terrain that muon rays pass through, which is often much larger than the target volume itself.

### Q: What is the `.g2zbin` format?

A binary format for 2D grids. The file begins with an 8-byte magic (`G2ZBIN\x00\x00`),
followed by a version byte, grid metadata (name, x-axis, y-axis), and z-values in
row-major order (float32 or float64). It is the single binary DEM format read by
MUONITH. Generate DEM files from GSI data or existing GeoTIFF files with the
external [muonith-gsi-dem](../auxiliary-tools/gsi-dem.md) utility.

## Output

### Q: How do I visualize the cross-section output?

Cross-section files contain (x, y, value) data. Example with Python:

```python
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("cross_section_z300.txt")
x, y, rho = data[:, 0], data[:, 1], data[:, 2]

nx = len(np.unique(x))
ny = len(np.unique(y))

plt.pcolormesh(
    x.reshape(ny, nx),
    y.reshape(ny, nx),
    rho.reshape(ny, nx),
    cmap="RdBu_r"
)
plt.colorbar(label="Density [kg/m³]")
plt.xlabel("x [m]")
plt.ylabel("y [m]")
plt.title("z = 300 m")
plt.axis("equal")
plt.show()
```

### Q: What does a negative `delta_prior` mean?

A negative $\vec{\rho}' - \vec{\rho}_0$ means the reconstructed density is lower than the prior. This could indicate a cavity, fracture zone, or low-density material at that location.

## Troubleshooting

### Q: "Matrix inversion failed" error

This typically occurs when the covariance matrix is ill-conditioned. Possible solutions:

1. Increase `sigma_rho` (stronger regularization)
2. Increase `corr_length` (smoother covariance)
3. Reduce the number of voxels (increase merge factors)
4. Use double precision: the code falls back to `mp_reconst_density_double()` for better numerical stability

### Q: Memory usage is too high

The dense covariance matrix scales as $O(n_v^2)$. For 50,000 voxels, this requires ~9.3 GB. Solutions:

1. Increase merge factors to reduce voxel count
2. Increase `n_hit_det_min` to exclude poorly-constrained voxels
3. Reduce the voxel grid extent (tighter `zmin`/`zmax`)

### Q: Path length computation is very slow

Trace Path Lengths (Module 4) (path length computation) is the most compute-intensive step. To speed it up:

1. Increase `OMP_NUM_THREADS`
2. Save results with `tf_save_arrdet_g3vox: true` and load on subsequent runs
3. Use coarser voxel resolution (larger merge factors)
4. Reduce the number of angular bins per detector
