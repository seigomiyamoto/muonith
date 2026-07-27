# Detector Model

MUONITH models muon detectors as a hierarchy of geometric elements that define the angular binning and spatial coverage of each observation.

## Class Hierarchy

```
DetectorPanelArray
  └── DetectorPanel          ← one physical detector
        └── DetectorElement  ← one angular bin (path-length computation element)
```

## DetectorElement

![DetectorElement concept — 3D geometry, angular space, and muon path quantities](../assets/images/chart_DetectorElement_concept.drawio.png)

A DetectorElement represents the **minimum angular element of the path-length computation** — a single directional bin on which path lengths, fluxes, and efficiencies are evaluated.

!!! warning "Bin size is a computational choice, not the instrument's angular resolution"
    The angular grid (`nbinx`, `nbiny`, `txmin`–`txmax`, `tymin`–`tymax` in the
    [detector parameter file](../user-guide/parameter-files.md#detector-parameter-files))
    only sets how finely the computation samples the terrain. What a bin can *actually*
    resolve is limited by three independent factors:

    - **Count statistics** — finer bins collect fewer muons each, so per-bin Poisson noise
      grows. The [automatic bin grouping](bin-grouping.md) merges bins where
      the signal is sparse to recover statistical significance.
    - **Detector geometry** — the physical channel size and the separation of the detection
      planes set the hardware angular response; bins finer than this response do not add
      information.
    - **Angular resolution effects** — multiple scattering in the rock and in the detector
      smears the arrival directions.

    Choose the binning fine enough to sample the terrain, and judge what is actually
    resolvable from the statistics and the detectability analysis
    ([Depth vs. Spatial Resolution](depth-resolution.md)) — not from the bin count.

| Property | Type | Description |
|----------|------|-------------|
| `unique_index` | `Uqid` (int) | Globally unique identifier |
| `detid` | `Detid` (int) | Parent detector panel ID |
| `id_in_this_detector` | `Inthis` (int) | Position within the parent panel |
| `ray3d` | `Ray3d` | Center position and look direction |
| `txmin`, `txmax` | `float` | Angular bounds (horizontal) |
| `tymin`, `tymax` | `float` | Angular bounds (vertical) |
| `effective_area_m2` | `float` | Detector collection area [m²] |
| `solid_angle` | `float` | Angular coverage [sr] |
| `exposure_time_sec` | `float` | Measurement duration [s] |

### Computed Quantities

After ray tracing, each element stores:

| Quantity | Unit | Description |
|----------|------|-------------|
| PL | m | Path length through rock |
| DL | kg/m² | Density length (integrated density along path) |
| Penetrating muon flux | 1/(m² sr s) | Expected muon intensity at the given DL |
| Signal / Noise | counts | Expected muon counts and statistical uncertainty |
| `proj_density` | kg/m³ | Projected (line-averaged) density; set to -1.0 if invalid |

!!! note "Efficiency uncertainty in the projected-density error band"
    The 2D projected-density error band (its upper/lower bounds per angular bin) widens when the detector efficiency uncertainty is enabled via `tf_eff_cn_diag`. The same efficiency-uncertainty variance added to the observation covariance diagonal $\mathbf{C}_N$ (see [Inversion](inversion.md#observation-covariance-mathbfc_n)) also feeds the grouped projected-density bounds. With the flag off, the band reflects count statistics only.

### Workflow

A typical DetectorElement lifecycle:

1. **Construct** — Set position, angular bounds, and exposure parameters
2. **Ray trace** — Compute PL and DL through Grid2dPillar (terrain)
3. **Flux lookup** — Determine expected muon flux from flux table using DL and zenith angle
4. **Signal estimation** — Calculate expected muon counts: flux × area × solid angle × time
5. **Inversion** — Use in the observation equation for density reconstruction

## DetectorPanel

![DetectorPanel concept](../assets/images/cls_DetectorPanel.dio.png)

![DetectorPanel — 3D panel geometry and DetectorElement grid](../assets/images/chart_DetectorPanel_concept.drawio.png)

A DetectorPanel represents a **complete physical muon detector** consisting of a 2D grid of angular elements.

| Property | Description |
|----------|-------------|
| `name` | Detector identifier |
| `ray3d` | Detector center position and orientation |
| `v3_det_length` | Physical dimensions [horizontal, vertical, depth] in meters |
| `n_unit` | Number of detector unit copies |
| `days` | Exposure time in days |
| `angle_unit` | Angular unit for this panel (Tangent / Degree / Radian) |

Internally, a DetectorPanel contains a **Grid2dBinGroup** — a 2D grid that discretizes the detector's angular field of view into bins. Each bin corresponds to one DetectorElement.

```
DetectorPanel
├── Grid2dBinGroup (angular bin grid)
│     ├── x-axis: horizontal angle bins
│     └── y-axis: vertical angle bins
└── vec_vec_DetectorElement[iy][ix]
      └── one element per angular bin
```

Where the signal is sparse, neighboring bins are adaptively merged into
groups, and a bidirectional group ↔ bin index map keeps the two layers
consistent. The grouping algorithm (`auto_divide_by_zsum`), a worked example,
and the index map `OneToManyUOBimap` have their own page:
[Bin Grouping](bin-grouping.md).

### Detector Orientation

The detector's viewing direction is defined by its Ray3d:

- **Position**: The physical location of the detector [m]
- **Direction**: The central look direction (unit vector)
- **Yaw**: 0° = North (+y), 90° = East (+x), measured clockwise when viewed from above

### From (tx, ty) to a ray direction

Each element's bin coordinates `(tx, ty)` are converted to a world-frame ray
direction in three steps: the bin **center** is taken in `angle_unit` space,
turned into a detector-local unit vector (slopes in Tangent mode, spherical
angles in Degree/Radian mode), and rotated by the panel's yaw/pitch/roll.
The resulting direction fixes the zenith angle for the flux lookup and thus
the expected count of the element.

The full derivation — including the Tangent-vs-Degree grid distortion and the
comparison with [muonith-path-view](../auxiliary-tools/path-view.md) — has its
own page: [From (tx, ty) to a Ray Direction](detector-angles.md).

## DetectorPanelArray

![DetectorPanelArray concept](../assets/images/cls_DetectorPanelArray.dio.png)

A container for **multiple detector panels**, providing unified access to all DetectorElement instances across all panels.

In a multi-view observation (used for 3D density reconstruction), each panel observes the target from a different position and angle. The DetectorPanelArray maintains:

- Global detector element indexing (Uqid)
- Coordinate transformations between detector-local and global frames
- Aggregated path length and observation data

### Shell path lengths (`ShellPL`)

In addition to path lengths through the voxel grid, Trace Path Lengths (Module 4) computes **density-independent** geometric path lengths through the terrain shells that surround the grid. These are bundled in a `ShellPL` struct:

| Field | Type | Description |
|-------|------|-------------|
| `upper` | `Eigen::VectorXf` | Path length through terrain above the voxel grid |
| `lower` | `Eigen::VectorXf` | Path length through terrain below the voxel grid |
| `lateral` | `Eigen::VectorXf` | Path length through terrain to the sides of the voxel grid |

Each vector is laid out **flat** across all detector elements of a `DetectorPanelArray`; per-panel slices are obtained with `segment(offset, n_ele)`. Because the values are density-independent, they are computed once in Trace Path Lengths (Module 4) and reused across every prior-density combination in Compute Prior (Module 5) (`density_quad = [prior, shell_upper, shell_lower, shell_lateral]`).

## Configuration

Detector parameters are specified in a JSON parameter file. Key settings include:

| Parameter | Description |
|-----------|-------------|
| Detector position | (x, y, z) coordinates in meters |
| Look direction | Yaw and elevation angles |
| Angular range | Min/max angles for horizontal and vertical axes |
| Angular bin size | Resolution of angular discretization |
| Physical size | Detector dimensions |
| Exposure time | Observation duration |
| Number of units | Detector copies (affects effective area) |

See [Parameter Files](../user-guide/parameter-files.md) for the complete parameter file format.
