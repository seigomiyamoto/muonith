# Muography (Muon Radiography)

## What are cosmic-ray muons?

Cosmic rays from outer space constantly bombard the Earth's atmosphere.
When high-energy protons and nuclei collide with atmospheric molecules,
they produce showers of secondary particles. Among these, **muons** are
particularly useful for imaging: they are highly penetrating, naturally
abundant at ground level (~10,000 per m² per minute on a horizontal
surface), and travel in approximately straight lines through hundreds
of meters of rock.

## Principle of muography

Muography exploits the fact that the number of muons surviving passage
through matter depends on the **density length** (density × path length)
along their trajectory. A denser region absorbs more muons, just as bones
absorb more X-rays than soft tissue.

By placing a muon detector on one side of a geological structure and
counting the muons arriving from different directions, one obtains a
**projected density image** -- a 2D map of the average density along each
line of sight.

$$
X = \int \rho(l)\, dl
$$

$$
N_\text{muon}(S, \Delta \Omega, T) = T \times\int\limits_{\Delta\Omega} F(X,d\Omega)  S(d\Omega) d\Omega
$$

where:

- $X$ is the density length (integrated density along the muon path)
- $F(X,d\Omega)$ is the penetrating muon flux as a function of density length
- $S(d\Omega)$ is the effective detector area
- $\Delta\Omega$ is the solid angle of the angular bin
- $T$ is the exposure time

## From 2D to 3D: multi-directional muography

A single viewing direction provides only a projected (2D) density image.
To reconstruct **3D density distributions**, observations from multiple
directions are required -- analogous to how CT scans combine many X-ray
projections.

Recent advances in portable nuclear-emulsion detectors have made it
feasible to deploy many detectors around a target structure, enabling
multi-directional muography.

However, the number of candidate detector positions, viewing angles, and
exposure-time configurations grows rapidly. **MUONITH** is designed to
efficiently evaluate and compare these configurations before committing to
fieldwork, letting researchers weigh the trade-offs up front rather than
prescribing a single "best" setup.

## Key physical quantities

| Quantity | Symbol | Unit | Description |
|---|---|---|---|
| Path length | $l$ | m | Geometric length of a muon's path through a voxel |
| Density length | $\rho \cdot l$ | kg/m² | Integrated density along a path |
| Penetrating muon flux | $F$ | 1/(m² sr s) | Number of muons per unit area, solid angle, and time |
| Expected muon count | $N_\text{muon}$ | counts | Predicted number of detected muons |

## Applications

Muography has been applied to:

- **Volcanoes** -- Imaging magma conduits, hydrothermal systems, and lava
  domes (e.g., [Showa-shinzan](https://en.wikipedia.org/wiki/Sh%C5%8Dwa-shinzan), Izu-Omuroyama, La Soufrière de Guadeloupe)
- **Pyramids** -- Discovering hidden chambers (ScanPyramids project)
- **Underground cavities** -- Detecting tunnels and voids
- **Nuclear reactors** -- Imaging fuel debris (Fukushima Daiichi)
- **Mining and civil engineering** -- Mapping density variations in the
  subsurface

## Further reading

- Tanaka, H.K.M. (2019). "Muography (muon radiography)." *Nature Reviews Methods Primers*.
- Nagahara, S. et al. (2022). "Three-dimensional density structure..." *Earth, Planets and Space*. [DOI:10.1186/s40623-022-01609-0](https://doi.org/10.1186/s40623-022-01609-0)
- Rosas-Carbajal, M. et al. (2017). "Three-dimensional density structure of La Soufrière de Guadeloupe..." *Journal of Geophysical Research*.
