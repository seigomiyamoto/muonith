# CLI Reference

This page lists the executable programs included in MUONITH and their command-line interfaces.

## Building

### Build Scripts

| Script | Description |
|--------|-------------|
| `bash bdebug.sh` | Debug build (`-O0`, logging enabled) |
| `bash brelease.sh` | Release build (`-O3`, NODEBUG defined) |
| `bash ccdebug.sh` | Clean debug build (clears ccache, removes `build-debug/`, rebuilds with DBGPRINT) |
| `bash ccrelease.sh` | Clean release build (clears ccache, removes `build-release/`, rebuilds) |
| `bash pdebug.sh` | Python-oriented debug build (DBGPRINT enabled) |

All scripts call `build.sh` internally, which handles:

- **Compiler selection**: Nix shell → Homebrew GCC (14→13) → Apple Clang → system GCC
- **Parallel builds**: Job count = stricter of (a) core headroom (`logical_cpus - max(2, ncpu/8)`) and (b) a memory cap of ~2 GB per compile job. Override with the `BUILD_JOBS` environment variable
- **ccache**: Enabled by default for faster incremental builds
- **OpenMP**: On non-Nix macOS, detects Homebrew `libomp` and passes `-DOpenMP_ROOT` to CMake

Compiled binaries are placed in `build-release/exec/` or `build-debug/exec/`.

### Environment Variables

| Variable | Description |
|----------|-------------|
| `OMP_NUM_THREADS` | Number of OpenMP threads for runtime parallelization |
| `BUILD_JOBS` | Override the build parallelism (number of compile jobs). Bypasses the core-headroom and memory-cap heuristics in `build.sh` |
| `IN_NIX_SHELL` | Set to `1` by the Nix devShell; changes compiler selection in `build.sh` |
| `CC` / `CXX` | Override compiler selection (e.g., `CC=gcc-14 CXX=g++-14 bash build.sh`) |
| `ENABLE_NODEBUG` | `ON`/`OFF` — define `NODEBUG` macro (default: `OFF`) |
| `ENABLE_DBGPRINT` | `ON`/`OFF` — define `DBGPRINT` macro (default: `OFF`) |

### Environment Setup

Before running, set up the thread count:

```bash
export OMP_NUM_THREADS=16    # Adjust to your CPU core count
```

## Main Executables

### `muonith.exe`

The primary executable for 3D density reconstruction and parameter sweep.

```
muonith.exe --json <config.json5> [--seed <seed>] [--resume <checkpoint_dir>] [--end-stage <N>]
```

| Option | Long form | Type | Required | Description |
|--------|-----------|------|----------|-------------|
| `-j` | `--json` | string | Yes | Path to JSON5 configuration file |
| `-s` | `--seed` | uint | No | Random seed (overrides JSON5 value) |
| `-r` | `--resume` | path | No | Resume a previously interrupted parameter sweep from the given checkpoint directory |
| `-E` | `--end-stage` | int | No | Stop the pipeline after completing Module N. Valid values: 3, 4, 5, 6, 7, 8. Default: 8 (all modules). Takes precedence over `end_stage` in JSON5. Incompatible with sweep mode (`tf_exec=true`) when N < 7. |

**Execution modes** (determined by JSON5 configuration):

- **Single run**: When `NAGAINV_PARAM_SWEEP` section is absent or `tf_exec=false`, runs the full pipeline once (modules 1-8).
- **Parameter sweep**: When `NAGAINV_PARAM_SWEEP.tf_exec=true`, runs modules 1-6 once, then sweeps module 7-8 over all parameter combinations.

**Examples:**

```bash
# Run full pipeline (modules 1-8, single run)
muonith.exe -j prm_muonith.json5

# Run with a specific seed
muonith.exe -j prm_muonith.json5 -s 12345

# Resume an interrupted parameter sweep (pass the checkpoint directory)
muonith.exe -j prm_muonith.json5 --resume path/to/checkpoint_dir
```

**Pipeline modules** (executed in order):

| Module | Name | Description |
|--------|------|-------------|
| 1 | Initialization | Load JSON5, setup logging, set random seed |
| 2 | Load Parameters | Parse all parameter sections |
| 3 | Build Geometry | Load DEM, create Grid2dPillar and Grid3dVoxel |
| 4 | Trace Path Lengths | Build detector array, compute path lengths |
| 5 | Compute Prior | Compute prior density-length and expected counts |
| 6 | Build Observation Matrix | Build observation matrix (dN/dD) |
| 7 | Invert Density | MAP density reconstruction |
| 8 | Analyze Errors and Output | Error analysis and final output |

See [Pipeline](../user-guide/pipeline.md) for detailed descriptions of each module.

The JSON5 file may include a `NAGAINV_PARAM_SWEEP` section specifying the parameter values to sweep. See [Parameter Files](../user-guide/parameter-files.md#nagainv_param_sweep-optional) for the configuration format.

### `depth_reso.exe`

Depth vs. resolution analysis tool. Evaluates what size and magnitude of density anomaly is detectable at a given depth for a specified detector configuration.

```
depth_reso.exe -j <config.json5> [-s <seed>]
```

| Option | Long form | Type | Required | Description |
|--------|-----------|------|----------|-------------|
| `-j` | `--json` | string | Yes | JSON5 configuration file |
| `-s` | `--seed` | uint | No | Random seed |

## Utility Executables

### `make_peneflux_dFdR.exe`

Generates muon flux table files (penetrating flux and dF/dR) from source data. These tables are required by the main pipeline.

```
make_peneflux_dFdR.exe <config.json> <SECTION_NAME>
```

See [Making Penetrating Muon Flux Tables](../user-guide/flux-tables/index.md) for the full walkthrough (inputs, sample config, outputs, and verification).

### `rotate_vector.exe`

Rotates one or more vectors by a yaw/pitch/roll triple and prints the resulting orientation (the rotated vector, its azimuth, and its elevation). It is a hands-on tool for checking the rotation conventions and for seeing Euler-angle degeneracies such as gimbal lock.

```
rotate_vector.exe <yaw_deg> <pitch_deg> <roll_deg> <LOCAL|GLOBAL> [vx vy vz ...]
```

| Position | Argument | Description |
|----------|----------|-------------|
| 1 | `yaw_deg` | Rotation about the Z axis, in degrees |
| 2 | `pitch_deg` | Rotation about the Y axis, in degrees |
| 3 | `roll_deg` | Rotation about the X axis, in degrees |
| 4 | `LOCAL` \| `GLOBAL` | Frame convention (see below) |
| 5+ | `vx vy vz ...` | Zero or more vectors to rotate, given as consecutive triples |

- **Angle order**: ZYX Euler angles. The three numbers are read as yaw, then pitch, then roll (yaw → Z, pitch → Y, roll → X). Positive angles are counter-clockwise.
- **Coordinate system**: right-handed, z-up. Input angles are in degrees.
- **Frame**: `GLOBAL` keeps the axes fixed in space (`R = Rx * Ry * Rz`); `LOCAL` rotates the axes with the object (`R = Rz * Ry * Rx`).
- **Default vectors**: when no vector is given, the three body axes `x = (1,0,0)`, `y = (0,1,0)`, `z = (0,0,1)` are rotated, so a full orientation (and any gimbal lock) is visible in one run.

**Examples:**

```bash
# Gimbal lock: pitch = 90 deg folds yaw and roll onto the same axis.
# With no vector given, the three body axes are rotated.
rotate_vector.exe 0 90 45 LOCAL

# Rotate a single vector (1,0,0) by yaw = 90 deg about Z (GLOBAL frame).
rotate_vector.exe 90 0 0 GLOBAL 1 0 0

# Rotate two vectors at once by roll = 90 deg about X.
rotate_vector.exe 0 0 90 GLOBAL 1 0 0 0 1 0
```

Invalid input (an unknown frame name, a non-numeric angle, or a vector list whose length is not a multiple of three) is reported and the program exits with status 1.

## CLI Precedence

When the same setting is specified in multiple places, the following priority applies:

```
CLI arguments  >  JSON5 file  >  Default values
```

For example, `-s 12345` on the command line overrides `"seed": 42` in the JSON5 file.

## Typical Workflow

```bash
# 1. Set up environment
export OMP_NUM_THREADS=16

# 2. Navigate to working directory
cd work/my_site/

# 3. Run full reconstruction (single run)
../../build-release/exec/muonith.exe -j prm_muonith.json5

# 4. Run parameter sweep (set tf_exec=true in NAGAINV_PARAM_SWEEP)
../../build-release/exec/muonith.exe -j prm_muonith.json5
```
