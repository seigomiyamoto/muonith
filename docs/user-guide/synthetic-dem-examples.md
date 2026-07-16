# Synthetic DEM Examples

This page walks through two ready-to-run demo stations built on **primitive-solid
DEMs** (Digital Elevation Models). Each station is a synthetic mountain of a simple
shape that carries a single rectangular density anomaly inside it. Running the
station produces a synthetic observed muon signal from that "true" density and then
inverts it, so you can see the whole forward-then-inverse pipeline end to end without
any real survey data or licensing constraints.

Two demo stations ship with the repository:

| Station | Shape | Config |
|---------|-------|--------|
| `eg-cylinder` | flat-topped cylinder | `param_sites/eg-cylinder.json5` |
| `eg-square-pyramid` | square pyramid | `param_sites/eg-square-pyramid.json5` |

Both reuse the tutorial detector layout, so no detector data has to be generated.

!!! note "The DEM binary is never committed"
    The `.g2zbin` DEM is **generated at run time** by `run_prg.sh` (via the
    `dem_generator` block in the station config) and written into
    `work/<station>/dem/`. It is not stored in the repository. Every run
    regenerates it, so the example is fully reproducible from the config alone.

## Background: what a primitive DEM is

A DEM here is a single-valued height field `z(x, y)` — one surface elevation per
horizontal grid point — stored in the repository's `g2zbin` binary format. The
generator `param_sites/tutorial/make_primitive_dem.py` can produce five shapes:

| `--shape` | Height field |
|-----------|--------------|
| `cylinder` | flat top at `base + height` within a radius, else `base` |
| `square_prism` | flat top within a half-width (Chebyshev square), else `base` |
| `cone` | `base + height * (1 - r / radius)` inside the radius |
| `square_pyramid` | `base + height * (1 - c / half_width)` inside the half-width |
| `gaussian` | `base + height * exp(-0.5 * ((dx/sigma_x)^2 + (dy/sigma_y)^2))` |

A single **unified width** knob keeps shapes comparable: `--width W` sets the
cylinder/cone diameter to `W`, the square side to `W`, and the Gaussian to
`4*sigma` (`±2 sigma`) across `W`. Concretely `radius = half_width = W/2` and
`sigma = W/4`. Explicit per-shape flags (`--radius`, `--half_width`, `--sigma_x`,
`--sigma_y`) override `--width` when both are given. See the
[`make_primitive_dem.py` reference](../reference/scripts.md#make_primitive_dempy)
for the full option table.

The demo stations do not call the generator by hand — the `dem_generator` block in
each station config records the shape and size, and `run_prg.sh` invokes the
generator for you before every run.

## Example 1: `eg-cylinder`

A flat-topped cylinder with a unified footprint width `W = 500 m` (radius 250 m),
`base = 40 m`, `height = 360 m` (summit ≈ 400 m), centered on the tutorial detector
ring.

### Run it

```bash
# 1. Scaffold the work directory (detectors come from the tutorial station).
bash setup_station.sh eg-cylinder

# 2. Regenerate the DEM and run the inversion.
cd work/eg-cylinder/swp001 && bash run_prg.sh true 42
```

`run_prg.sh` takes `<run_exe: true|false> <seed>`: `true` runs the executable and
`42` is the random seed. How far the pipeline runs is taken from the station config
(`prm_muonith.json5`), so you do not have to specify it on the command line.

### What happens

First the DEM is regenerated (the binary is not in the repo):

```text
Wrote 90601 cells to ../dem/eg_cylinder-5m.g2zbin (301x301 grid, z range [40.000, 400.000] m)
```

Then the internal density anomaly is loaded and applied to the 3D voxel grid:

```text
read cuboid parameters...num of cuboid = 1
adding cuboid density structure: cuboid_demo
```

The [processing pipeline](pipeline.md) then runs to completion:

```text
Sweep completed: 1 configurations executed
```

### Outputs

Plots land in `work/eg-cylinder/swp001/figs/` (density cross-sections, per-detector
path-length and signal maps, and the reconstructed density). Cross-section binaries
and CSVs are written under the run directory. See [Output](output.md) for the full
list and naming scheme.

## Example 2: `eg-square-pyramid`

Same footprint width (`W = 500 m`, so a 500 m base side), `base` and `height` as
`eg-cylinder`, but the surface slopes linearly to the summit instead of being flat
on top — this changes how the terrain thickness varies across each detector's field
of view.

```bash
bash setup_station.sh eg-square-pyramid
cd work/eg-square-pyramid/swp001 && bash run_prg.sh true 42
```

The run proceeds exactly as above; only the terrain shape (and therefore the DEM
`eg_square_pyramid-5m.g2zbin` and the forward signal) differs.

## The internal density anomaly

Each demo station places one rectangular density anomaly inside the mountain via the
`cuboid_params` block. This is what the inversion tries to recover. In both stations
it is a single dense box near the center:

```json5
"cuboid_params": [
  {
    "tf_exec": true,
    "name": "cuboid_demo",
    "delta_density": 600,        // +600 kg/m^3 relative to the surrounding rock
    "xcnt": 130.0, "ycnt": 150.0, "zcnt": 200.0,   // center (m)
    "xlen": 200.0, "ylen": 200.0, "zlen": 120.0,   // side lengths (m)
    "theta_x_deg": 0.0, "theta_y_deg": 0.0, "theta_z_deg": 0.0,  // no rotation
    "rotation_type": "LOCAL",
  },
]
```

`delta_density` is added to whatever background density the voxel already has, only
inside the box. The box can be rotated with `theta_x_deg` / `theta_y_deg` /
`theta_z_deg` (yaw/pitch/roll); `rotation_type` selects `LOCAL` or `GLOBAL` axes,
matching the existing `CylinderParameters` convention. On the C++ side this is read
by `Grid3dVoxel::CuboidParameters` and applied by `add_density_structure`, so the
forward model computes the synthetic observed signal from a density field that
contains the box.

## Trying other shapes

`make_primitive_dem.py` also supports `cone`, `square_prism`, and `gaussian`. There
is no pre-built station for these, but you can add one by copying
`param_sites/eg-cylinder.json5` to `param_sites/eg-<shape>.json5` and changing
`dem_generator.shape` (and the size fields if desired). `setup_station.sh` and
`run_prg.sh` then work the same way — no other changes are needed.
