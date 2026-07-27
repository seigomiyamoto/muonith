# Bin Grouping

A DetectorPanel discretizes its angular field of view with a **Grid2dBinGroup** —
a 2D grid of angular bins, one [DetectorElement](detector-model.md#detectorelement)
per bin. Where the muon signal is sparse, neighboring bins are merged into groups
to recover statistical significance. This page explains the grouping algorithm and
the index map that binds grouped and individual bins.

The panel structure holding this grid is described in
[Detector Model — DetectorPanel](detector-model.md#detectorpanel).

![Grid2dBinGroup — grouped 2D grid bins for signal statistics](../assets/images/chart_Grid2dBinGroup_concept.drawio.png)

## Automatic Bin Grouping

When `auto_divide_by_zsum` is enabled, the angular bins are recursively subdivided based on signal distribution:

![auto_divide_by_zsum flow](../assets/images/auto_divide_by_signal_sum.dio.png)

The example below shows the effect on forward-modeled signal data (synthetic 11-detector example, detector 01). The left panel is the per-element signal field in angular (tangent) coordinates with a logarithmic color scale; the right panel is the same signal after grouping, where each cell is one Grid2dBinGroup group — coarse where the signal is sparse, fine where the signal is dense.

![Omuro DEM with the 11 synthetic detectors placed around the dome](../assets/images/omuro_dem_detectors.png)

*Geographic setting of the synthetic example: the Omuro DEM (5 m resolution, elevation in meters above sea level) with the 11 detectors placed around the dome. Numbers along the lines indicate each detector's horizontal distance to the summit. The two panels below show the per-element and grouped signal for detector 01 (`det_2018_B`), located on the western foot of the dome.*

| Per-element signal (before grouping) | Grouped bins (after grouping) |
|:-:|:-:|
| ![Per-element signal field in angular coordinates](../assets/images/bingroup_signal_per_element_det01.png) | ![Adaptively grouped Grid2dBinGroup bins](../assets/images/bingroup_grouped_det01.png) |

## OneToManyUOBimap — Bin Index Mapping

![OneToManyUOBimap — bidirectional one-to-many map for bin group indexing](../assets/images/chart_OneToManyUOBimap_concept.drawio.png)

Grid2dBinGroup uses **OneToManyUOBimap** to maintain bidirectional mappings between grouped bins and individual bins. This data structure provides:

| Direction | Lookup | Complexity |
|-----------|--------|------------|
| One → Many | `get_vec_Many(One)` — retrieve all individual bins belonging to a group | O(k) where k = group size |
| Many → One | `getOne(Many)` — find which group an individual bin belongs to | O(1) |

**Key constraints:**

- Each `Many` value maps to exactly **one** `One` key (uniqueness enforced on insert)
- Each `One` key can map to **multiple** `Many` values
- `insertOrOverwrite()` allows reassignment of a `Many` value to a different `One` key
