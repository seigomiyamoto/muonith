# Input Data

This page describes the input data files required by MUONITH and their formats.

## Required Inputs

| Input | Format | Description |
|-------|--------|-------------|
| DEM file | `.g2zbin` (binary) | Digital Elevation Model — terrain surface |
| Flux tables | `.g2zbin` (binary) | Muon flux and flux derivative tables |
| Detector parameter files | `.json5` | One per detector panel |
| Main configuration | `.json5` | Pipeline parameters |

## DEM (Digital Elevation Model)

The DEM defines the terrain surface that muons pass through.

### Format

Binary file (`.g2zbin`) — a compact, little-endian uniform-grid format. A fixed header (magic `"G2ZBIN\0\0"`, version, `z` precision, storage order) is followed by a `Grid2d` header that defines the x and y axes (name, bin count, min, max, interval), then a flat row-major array of `z` (elevation) values (float32 by default, or float64). Each point's x/y is reconstructed from the axes, so only `z` is stored. Conceptually it represents a regular grid of 3D points:

| Quantity | Stored? | Type | Unit | Description |
|----------|---------|------|------|-------------|
| x | reconstructed from x-axis | float/double | m | Easting coordinate |
| y | reconstructed from y-axis | float/double | m | Northing coordinate |
| z | stored (row-major array) | float32/float64 | m | Elevation (m a.s.l.) |

The grid structure (axis ranges and interval) comes directly from the `.g2zbin` header. `tolerance_ratio` validates that the grid spacing stays uniform within the given relative tolerance.

### Configuration

```json5
"GRID2D_PILLAR_PARAMETERS": {
  "path_dem": "dem/my_terrain.g2zbin",
  "initial_uniform_density": 2000,     // kg/m³
  "zmin": 0.0,                          // Base elevation [m]
  "tf_shift_x": true,                   // Shift to local coordinates
  "tf_shift_y": true,
  "tolerance_ratio": 0.01               // Grid spacing tolerance
}
```

### Coordinate Shift

![tf_shift_x/y behavior](../assets/images/Grid2d_tf_shift_xy.dio.png)

When `tf_shift_x` and `tf_shift_y` are `true`, the grid origin is shifted so that the minimum x and y values become the grid origin. This is useful when DEM coordinates are in a geographic system with large absolute values (e.g., UTM coordinates).

### Preparing a DEM

DEM data can be obtained from:

- National geographic surveys (GSI in Japan, USGS in the US)
- Satellite-based DEMs (SRTM, ALOS)
- Airborne LiDAR surveys

The data must be converted to the `.g2zbin` binary format. See [muonith-gsi-dem](../auxiliary-tools/gsi-dem.md) for a tool that converts GSI (Japan) DEM ZIP archives to `.g2zbin` format. The grid spacing determines the horizontal resolution of the terrain model.

## Muon Flux Tables

MUONITH requires pre-computed muon flux tables for two purposes:

1. **Forward modeling**: Predicting expected muon counts for a given density distribution
2. **Sensitivity matrix**: Computing how muon counts change with density perturbations

### Required Tables

| Table | Content | Axes |
|-------|---------|------|
| Penetrating muon flux | $\log_{10} F$ [1/(m² sr s)] | cos(θ_zenith) × density-length [kg/m²] |
| Flux derivative dF/dR | $dF/dR$ | density-length × cos(θ_zenith) |

### Configuration

Two flux table sections are required — one for the prior model and one for the observed (real) data. They often point to the same files:

```json5
"FLUX_RANGE_DATA_TABLE_PRIOR": {
  "pathin_log_peneflux": "fluxtable/costhz-kgm2-log10peneflux-allinone.g2zbin",
  "pathin_dFdR_R_costhz": "fluxtable/g2_dFdR_R_costhz.g2zbin",
  "tf_xcnt_peneflux": true,
  "tf_ycnt_peneflux": false,
  "tf_xcnt_dFdR": true,
  "tf_ycnt_dFdR": false
},

"FLUX_RANGE_DATA_TABLE_REAL": {
  // Same format as PRIOR — can point to the same or different files
}
```

The `tf_xcnt_*` and `tf_ycnt_*` flags control whether the lookup uses bin centers or bin edges for interpolation.

### Flux Models

MUONITH ships with flux tables derived from the **daemonflux** model (combined with the
Groom range-energy table), located in `fluxtable/daemon_groom/`. This is the default
model for generated station configs. Tables based on the Honda model are kept in
`fluxtable/honda_groom/` for comparison, but have no confirmed redistribution permission.

The bundled `.g2zbin` files are not raw daemonflux samples; they are pre-computed
tables already interpolated from the original data with a 2D cubic spline.

To regenerate these tables from the source data, or to introduce another flux model or
absorber material, see [Making Penetrating Muon Flux Tables](flux-tables/index.md).

## Detector Parameter Files

Each detector panel is configured in its own JSON5 file. See [Parameter Files](parameter-files.md#detector-parameter-files) for the complete format.

### Key Parameters

| Parameter | Unit | Description |
|-----------|------|-------------|
| `x`, `y`, `z` | m | Detector position |
| `yaw_deg` | degrees | Viewing direction (0° = East, 90° = North) |
| `nbinx`, `nbiny` | — | Angular binning (path-length computation elements — not the instrument's angular resolution, see [Detector Model](../concepts/detector-model.md#detectorelement)) |
| `txmin`–`txmax`, `tymin`–`tymax` | (unit depends on `angle_unit`) | Angular range |
| `length_hori`, `length_vert`, `length_dept` | m | Physical detector dimensions |
| `n_unit` | — | Number of detector copies |
| `days` | days | Exposure time |

### Detector Efficiency

When detector efficiency is modeled (`tf_eff_cn_diag`), each detector element carries, per angular bin `(tx, ty)`, a central efficiency `eff_cnt` and a band `[eff_low, eff_upp]`. Two input routes exist per detector: an analytic efficiency model (`path_eff_model`, JSON5, preferred) and a per-bin text table (`path_eff_table`, for measured or otherwise non-analytic efficiency). See [Detector Efficiency](detector-efficiency/index.md) for the file formats, the evaluation formulas, and the precedence rule.

### Detector Placement

Detectors should be placed where they have a clear line of sight to the target volume. Key considerations:

- **Elevation**: Detectors below the target observe more rock, but with lower muon flux
- **Distance**: Closer detectors have finer angular resolution but narrower coverage
- **Multi-view**: For 3D reconstruction, place detectors at different azimuthal angles around the target

## Observed Muon Count Data

For real-data analysis, observed muon counts are provided through the detector parameter framework. In synthetic tests, muon counts are generated internally by the forward model using the "real" flux table and the input density distribution.

## Directory Structure Example

A typical working directory:

```
work/my_site/
├── prm_rec.json5              # Main configuration
├── dem/
│   └── my_terrain.g2zbin      # DEM data
├── detparams/
│   ├── det_00.json5           # Detector 0
│   ├── det_01.json5           # Detector 1
│   └── ...
├── fluxtable/ → ../fluxtable  # Symlink to shared flux tables
└── logs/                      # Log output (auto-created)
```
