# Utility Scripts

The `scripts/` directory contains Python utility scripts for data conversion, visualization, parameter management, and inspection. All scripts require **Python 3** and accept `--help` for full option details. Dependencies (numpy, matplotlib, pandas, scipy, pyproj, rasterio, json5, etc.) are installed with `uv sync` from `pyproject.toml` + `uv.lock` at the repository root (see [Installation](../getting-started/installation.md#python-packages-optional)).

## Overview

| Script | Purpose |
|--------|---------|
| `g2zbin_io.py` | Library module for reading/writing g2zbin files |
| `kml_to_csv.py` | Extract KML placemarks to CSV/Markdown with elevation lookup |
| `geojson_to_csv.py` | Convert GeoJSON points to CSV/Markdown analysis table |
| `mk_detjson5_from_template.py` | Generate per-detector JSON5 files from a CSV and template |
| `auto_plot.py` | Auto-detect output files and run appropriate plotting scripts |
| `plot_Grid3dVoxel.py` | Render 3D voxel cross-section plots with GIF/PDF output |
| `plot_detection_limit_heatmap.py` | Plot detection limit heatmaps (p-value vs. size and depth) |
| `hist2d.py` | Unified 2D histogram / grid plotter |
| `plot_dem.py` | Visualize g2zbin DEM as a 2D heatmap |
| `plot_dem_with_detectors.py` | Plot DEM with detector positions overlay |
| `plot_g2bg.py` | Plot Grid2d background fields per detector |
| `plot_det_arrdet.py` | Plot txty/g2bg figures straight from `det/arrdet_*.bin` binary checkpoints |
| `make_primitive_dem.py` | Generate a primitive-solid DEM (cylinder/square_prism/cone/square_pyramid/gaussian) as g2zbin (in `param_sites/tutorial/`) |
| `init_work_site.py` | Scaffold a new per-site working directory from KML/CSV inputs and JSON5 templates |

A small number of shell helpers also live at the repository root:

| Script | Purpose |
|--------|---------|
| `setup_station.sh` | Station setup wrapper. Auto-checks `/usr/bin/time`, `python3`, `taskset`; sets up `.venv/` via `uv sync` from `pyproject.toml` + `uv.lock` (falls back to venv+pip when `uv` is absent); auto-copies `param_sites/template.json5` for new stations; respects `data_source` for shared data dirs; relaxes KML requirement when `skip_detparams: true`. Pass `--skip-env-check` to bypass Phase 0 / 0-5 environment checks. |

---

## Data Conversion

### g2zbin_io.py

!!! info "Library Module"
    This is a Python library module, not a standalone script. Import it in other scripts with `from g2zbin_io import read_g2zbin, write_g2zbin`.

Provides `read_g2zbin()` and `write_g2zbin()` functions for the g2zbin binary format. The format stores Grid2d metadata (axis names, bin counts, min/max, interval) in a compact header followed by a flat z-value array. Used internally by `plot_dem.py`, `kml_to_csv.py`, and other grid-processing scripts.

---

## Geospatial & Detector Setup

### kml_to_csv.py

Extract point placemarks from a KML file and produce a CSV and Markdown analysis table with projected coordinates, distances, and azimuths.

```bash
python3 scripts/kml_to_csv.py --kml_in input.kml --EPSG 6676 [options]
```

Supports elevation lookup from the GSI (Geospatial Information Authority of Japan) API, GeoTIFF DEM, or g2zbin DEM. Automatically identifies a reference point (vent/summit) using configurable keywords and computes relative geometry for each detector candidate.

| Option | Description |
|--------|-------------|
| `--kml_in` | Input KML file (required) |
| `--EPSG` | Output EPSG code (required) |
| `--fill_elevation_from_gsi` | Fetch elevations from the GSI API |
| `--dem` | Path to DEM file (`.tif` or `.g2zbin`); GSI API used as fallback |
| `--throttle_sec` | API call interval in seconds (default: 0.5) |

```bash
python3 scripts/kml_to_csv.py --kml_in asama.kml --EPSG 6676 --dem dem.g2zbin --verbose
```

### geojson_to_csv.py

Convert GeoJSON point features to the same CSV/Markdown analysis table format as `kml_to_csv.py`.

```bash
python3 scripts/geojson_to_csv.py --geojson_in input.geojson --EPSG 6680 [options]
```

Accepts GeoJSON FeatureCollection files exported from GSI Maps or other GIS tools. Shares the elevation-fill logic with `kml_to_csv.py`, supporting `--fill_elevation_from_gsi`, `--dem_tif`, and `--dem_only` modes.

```bash
python3 scripts/geojson_to_csv.py --geojson_in tokachi.geojson --EPSG 6680 --dem_only --dem_tif dem.tif
```

### mk_detjson5_from_template.py

Generate per-detector JSON5 parameter files from a template and a CSV table.

```bash
python3 scripts/mk_detjson5_from_template.py --template template-det.json5 --csv table.csv
```

Reads a CSV produced by `kml_to_csv.py` and fills in x, y, z, and yaw_deg values for each detector row while preserving comments in the JSON5 template. Column detection is automatic (handles various EPSG label formats) but can be overridden with `--xcol`/`--ycol`. Rows with missing numeric values are safely skipped.

| Option | Description |
|--------|-------------|
| `--template` | Path to template JSON5 file (required) |
| `--csv` | Input CSV file (required) |
| `--z_pos_offset` | Add an elevation offset to all z values |
| `--epsg` | Prefer columns matching this EPSG code |

```bash
python3 scripts/mk_detjson5_from_template.py --template template-det.json5 --csv asama_table_epsg6676.csv --z_pos_offset 10.5
```

### init_work_site.py

Scaffold a new per-site working directory from a station JSON5 (under
`param_sites/`) that points at a KML/GeoJSON for detector placement and a
DEM for terrain. Produces a site directory containing `prm_muonith.json5`,
per-depth/sweep subdirectories (e.g. `depth001/`, `swp001/`), and ready-to-use
`auto_plot.json5` configurations.

```bash
python3 scripts/init_work_site.py param_sites/<station>.json5 --run
```

Typical use is to bootstrap a new station (e.g. `tarumae-dome`, `tarumae-mid01`)
without hand-editing the full JSON5 tree. Reconstruction volume fields are
written in **meters** (`x_aabb_meters`, `cylinder_radius_x_meters`, …)
matching the current `GRID3D_VOXEL_PARAMETERS.reconst_voxels` schema.

**Custom CRS stations (`skip_detparams`)**: For stations whose DEM and
detector positions use a non-EPSG custom orthogonal coordinate system,
set `"skip_detparams": true` in the station JSON5. The KML/CSV/det
conversion pipeline is then bypassed, and the existing
`work/<station>/detparams/<label>/runcard_det_*.json5` files are reused
as-is. Required keys are relaxed (`epsg`, `detector_input`, `detector_spec`
become optional). The `detector_spec.label` field can be supplied as a
hint to pin which existing label directory to use when multiple are present.

---

## Visualization

### auto_plot.py

Automatically detect output files and orchestrate the appropriate plotting scripts.

```bash
python3 scripts/auto_plot.py [--config config.json5] [--dry-run] [--exclude PATTERN]
```

Scans the `tmp/` directory (and sweep subdirectories) for `.tmp` files, then dispatches hist2d, g3vox cross-section, and g2bg plotting tasks. hist2d tasks run in parallel; g3vox, g2bg, and det tasks run sequentially. Configuration is read from `auto_plot_default.json5` and can be overridden with a user config file.

A fourth task family, **det**, is disabled by default: when the config sets `det_exec: true`, `auto_plot.py` also scans `det/arrdet_*.bin` and `checkpoint_*/det/arrdet_*.bin` — both directly in the working directory and under the tmp directory (the sweep driver writes its checkpoint bundles there, e.g. `tmp/checkpoint_m4/det/`) — and dispatches [`plot_det_arrdet.py`](#plot_det_arrdetpy) for each binary. This produces the standard txty/g2bg figures directly from binary checkpoints, without needing the `tmp/*.tmp` text exports. Keep `det_exec` off when the text exports are present — the same figures would be produced twice. The swp001 template (`scripts/templates/swp001_auto_plot.json5.tpl`) sets `det_exec: true`, because the swp001 `prm_muonith.json5` it is generated with disables the text exports (`tf_out_txty_ascii` / `tf_out_g2bg_ascii` false). The `det` section's `out_dir` / `n_jobs` keys control the output directory and the drawing parallelism (see [`plot_det_arrdet.py`](#plot_det_arrdetpy)); the swp001 template sets `out_dir: "figs"`, `n_jobs: -2`, and excludes `arrdet_g2pil_*`.

| Option | Description |
|--------|-------------|
| `--config`, `-c` | Custom JSON5/JSON config to override defaults |
| `--dry-run`, `-n` | Print commands without executing |
| `--skip-hist2d` | Skip hist2d plots |
| `--skip-g3vox` | Skip g3vox plots |
| `--skip-g2bg` | Skip g2bg plots |
| `--skip-det` | Skip det plots (binary-direct `plot_det_arrdet.py` tasks) |
| `--exclude`, `-e` | Exclude files matching glob pattern (repeatable) |
| `--delete-png` | Delete intermediate PNG files after GIF/PDF creation |
| `--tmp-dir` | Directory containing `.tmp` files (default: `tmp`) |
| `--error-log` | File for error logs (default: `error_auto_plot-tmp.log`; empty string disables) |

```bash
python3 scripts/auto_plot.py --dry-run --exclude "arrdet_g3vox_prior_*"
```

### plot_Grid3dVoxel.py

Render 3D voxel data as cross-section heatmap slices with animated GIF and PDF output.

```bash
python3 scripts/plot_Grid3dVoxel.py --config=config.json [--use-data-bins]
```

Reads `g3vox_*_zcross_all.tmp` (text) or `g3vox_*_zcross_all.tmpbin` (binary) files and produces per-slice PNG images along the specified cross-section axis (x, y, or z). Supports detector position markers, multiple cross-sections, and configurable color scales. Uses multiprocessing for parallel slice rendering.

| Option | Description |
|--------|-------------|
| `--config` | Path to JSON configuration file |
| `--use-header-bins` | Use bin parameters from the input file header |
| `--use-data-bins` | Compute bin edges from actual data coordinates |
| `--delete-png` | Remove intermediate PNG files after GIF/PDF creation |

```bash
python3 scripts/plot_Grid3dVoxel.py --config=plot_g3vox.json --use-data-bins
```

### hist2d.py

Unified 2D histogram and grid plotter with two modes: histogram (bin scatter data) and grid (reshape pre-gridded data).

```bash
python3 scripts/hist2d.py input.tmp [--mode histogram|grid] [options]
```

Reads whitespace-separated (x, y, z) data and produces a color-mapped PNG. Supports log scale, contour lines, mask overlays, and CSV export. Typically invoked by `auto_plot.py` rather than directly.

| Option | Description |
|--------|-------------|
| `--mode` | `histogram` (default) or `grid` |
| `--nbinx`, `--xmin`, `--xmax` | X-axis binning |
| `--ngrad`, `--vmin`, `--vmax` | Color scale range |
| `--log` | Use logarithmic color scale |
| `--colormap` | Matplotlib colormap name (default: jet) |
| `--contour` | Draw contour lines |
| `--no-csv` | Disable CSV output |

```bash
python3 scripts/hist2d.py tmp/hist_data.tmp --nbinx 320 --xmin -1.6 --xmax 1.6 --log
```

### plot_detection_limit_heatmap.py

Create detection limit heatmaps showing log10(p-value) as a function of anomaly size and depth for each delta-density value.

```bash
python3 scripts/plot_detection_limit_heatmap.py --main_config prm_reso.json5 [--plot_config heatmap.json5] [--outdir figs]
```

Reads significance CSV files produced by the depth-resolution sweep and generates separate heatmap figures per delta-density. Supports configurable contour lines, a secondary elevation axis, size-bar indicators, and optional animated GIF output.

```bash
python3 scripts/plot_detection_limit_heatmap.py --main_config prm_reso.json5 --plot_config heatmap_config.json5 --surface_elevation auto
```

### plot_dem.py

Visualize g2zbin DEM files as 2D heatmaps with optional contour lines.

```bash
python3 scripts/plot_dem.py --input dem.g2zbin [--contour] [--cint 10]
```

Supports color range control, contour line intervals, optional 90-degree rotation, and configurable colormaps. Saves output as PNG alongside the input file.

| Option | Description |
|--------|-------------|
| `--input` | Input `.g2zbin` file (required) |
| `--zmin`, `--zmax` | Color scale range |
| `--contour` | Draw contour lines |
| `--cint` | Contour interval (default: 10) |
| `--cmap` | Matplotlib colormap (default: terrain) |

```bash
python3 scripts/plot_dem.py --input dem.g2zbin --contour --cint 50 --cmap terrain
```

### plot_dem_with_detectors.py

Plot a DEM together with detector positions loaded from a CSV file or JSON5 configuration.

```bash
python3 scripts/plot_dem_with_detectors.py --config prm.json5 [--csv det_table.csv]
```

Overlays detector markers, a reference-point star, and optional distance labels on top of a DEM heatmap with contours. Accepts DEM in g2zbin or GeoTIFF format. Coordinates can be displayed in meters or kilometers.

| Option | Description |
|--------|-------------|
| `--config` | JSON5 config (extracts DEM path, detectors, target coordinates) |
| `--csv` | CSV with detector coordinates (from `kml_to_csv.py`) |
| `--dem` | Path to DEM file (overrides config) |
| `--show_dist` | Show detector-to-target distance labels |
| `--xyunit` | Axis unit: `m` or `km` (default: m) |

```bash
python3 scripts/plot_dem_with_detectors.py --config prm.json5 --contour --cint 50 --xyunit km
```

### plot_g2bg.py

Plot Grid2d background fields (signal, noise, density estimates, efficiency) per detector.

```bash
python3 scripts/plot_g2bg.py --jsonfile=config.json
```

Reads `g2bg_*_det*.tmp` files and generates PNG images for each configured field and detector. Supports reference-data differencing, virtual fields (e.g., S/N ratio differences), log scale, and parallel processing. Produces per-detector and per-field PDF/GIF aggregations.

```bash
python3 scripts/plot_g2bg.py --jsonfile=plot_g2bg.json
```

### plot_det_arrdet.py

Plot per-element quantities straight from `det/arrdet_*.bin` binary checkpoints (binary-direct, no intermediate text file).

```bash
python3 scripts/plot_det_arrdet.py <arrdet.bin | checkpoint_dir | det_dir> \
    [--field PL,dens,signal] [--det 0,1,2] [--out DIR] [--jobs N] \
    [--plot-config auto_plot.json5] [--check-mat <mat_dir>] [--list-fields]
```

Reads a serialized `DetectorPanelArray` checkpoint (`det/arrdet_g3vox_input.bin` or `det/arrdet_g2pil_naive.bin`, ~160 MB each) field by field, then reproduces the standard txty (hist2d) and g2bg figures with the rendering parameters from the same `auto_plot.json5` the normal flow uses, so the output naming and style match the established figures. Useful when only binary checkpoints are available — for example a `--end-stage 4` bundle, or runs where the `tmp/*.tmp` text exports were disabled or predate the `mat/` element-table export.

Also invoked automatically by [`auto_plot.py`](#auto_plotpy) as the **det** task family when the config sets `det_exec: true` (default: false). Figures are written to the `det.out_dir` of the auto_plot config (relative paths are anchored at the config's folder; the swp001 template sets `"figs"`), else to `plot_det_arrdet/` next to the binary; `--out` overrides both. Drawing runs in parallel: `--jobs` (else `det.n_jobs`, default -2) sets the number of workers — positive = that many, 1 = one by one, 0 = all CPU cores, negative = leave that many cores free.

Limitations: checkpoints older than `PIPELINE_VERSION` 7 are rejected with a clear error (v7 added the group efficiencies to the binary, so the `eff_*` figures are drawn on this path too) — regenerate them with the current code. For very large checkpoints the binary read can be slow — run the script manually if the `auto_plot.py` per-task time limit (600 s for det tasks) is exceeded.

| Option | Description |
|--------|-------------|
| `--field` | Comma-separated txty fields to plot (default: `PL,dens,signal`; see `--list-fields`) |
| `--det` | Comma-separated detector IDs to keep (default: all) |
| `--out` | Output directory (default: `det.out_dir` from the auto_plot config, else `<source_dir>/plot_det_arrdet`) |
| `--jobs` | Figure-drawing workers (default: `det.n_jobs` from the auto_plot config, else -2 = leave 2 CPU cores free) |
| `--plot-config` | `auto_plot.json5` for rendering parameters (default: nearest one above the binary, else `auto_plot_default.json5`) |
| `--check-mat` | Validate the binary read against the `mat/` element table row for row |
| `--g2bg` / `--no-g2bg` | Force g2bg group figures on/off (default follows the config `g2bg_exec` switch) |
| `--no-hist2d` | Never draw the per-element txty (hist2d) figures |
| `--list-fields` | Print the plottable field names and exit |

```bash
python3 scripts/plot_det_arrdet.py checkpoint_trace_path_lengths --check-mat run_output/mat
```

---

## Inspection & Comparison

### demo_recon_nnls.py

Standalone density reconstruction driven only by the exported `mat/manifest.json`, with no dependency on the C++ code.

```bash
python3 scripts/demo_recon_nnls.py <mat_dir> [--out DIR] [--lams 0.01,0.1,1.0]
```

Loads `mat/manifest.json`, sanity-checks every `.bin` it lists, builds a synthetic ground truth (uniform background with one low-density block), synthesizes observed counts from the central observation matrix, and solves non-negative regularized least squares over each `--lams` weight, keeping the best by relative error. Writes a slice PNG and a z-scan GIF to `--out` (default `<mat_dir>/../demo_recon_nnls`). Requires NumPy, SciPy, and Matplotlib. See [External Reconstruction](../user-guide/external-reconstruction.md) for the `mat/` layout.

```bash
python3 scripts/demo_recon_nnls.py path/to/run_output/mat
```

---

## Test Data Generation

### make_primitive_dem.py

Generate a simple primitive-solid DEM (single-valued height field) and write it as a g2zbin binary file. Lives in `param_sites/tutorial/`. Supersedes the former `mk_gaussian_g2zbin.py` (the Gaussian surface is now the `gaussian` shape). The generated `.g2zbin` binaries are not committed to the repository; users regenerate them locally for tutorials.

```bash
python3 param_sites/tutorial/make_primitive_dem.py --shape cone \
    --x_min -500 --x_max 500 --y_min -500 --y_max 500 \
    --x_interval 5 --y_interval 5 --radius 200 --height 300 \
    --outbin param_sites/tutorial/eg_cone-5m.g2zbin
```

Supported shapes (`--shape`): `cylinder`, `square_prism`, `cone`, `square_pyramid`, `gaussian`. Each computes a uniform grid and writes the height field in g2zbin format.

| Option | Description |
|--------|-------------|
| `--shape` | One of cylinder / square_prism / cone / square_pyramid / gaussian (required) |
| `--x_min`, `--x_max` | X-axis range |
| `--y_min`, `--y_max` | Y-axis range |
| `--x_interval`, `--y_interval` | Grid spacing |
| `--xcnt`, `--ycnt` | Shape center (default 0, 0) |
| `--base` | Base plane elevation (default 0) |
| `--height` | Peak height above base (default 300) |
| `--radius` | Radius for cylinder / cone (default 200) |
| `--half_width` | Half-width for square_prism / square_pyramid (default 200) |
| `--sigma_x`, `--sigma_y` | Gaussian widths (sigma_y defaults to sigma_x) |
| `--outbin` | Output g2zbin file path (default: `eg_<shape>-<dx>m.g2zbin`) |
| `--float64` | Store z as float64 instead of float32 (default: float32) |

For a full worked example that generates a primitive DEM, adds an internal density
anomaly, and runs the inversion end to end, see
[Synthetic DEM Examples](../user-guide/synthetic-dem-examples.md).
