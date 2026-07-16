# MUONITH

**MUOgraphic Numerical Inversion & Tomography Harness**

A fast C++ simulation tool for multi-directional muography observation design
and 3D density reconstruction.

---

## What is muography?

Muography uses naturally occurring cosmic-ray muons to image the internal
bulk-density structure of large objects such as volcanoes ([Tanaka et al., 2007](https://doi.org/10.1029/2007GL031389); [Nishiyama et al., 2014](https://doi.org/10.5194/gi-3-29-2014); [Nishiyama et al., 2022](https://doi.org/10.1002/9781119722748.ch3)), pyramids ([Morishima et al., 2017](https://doi.org/10.1038/nature24647)),
underground cavities, and scoria cones such as [Mount Ōmuro](https://en.wikipedia.org/wiki/Mount_%C5%8Cmuro_(Shizuoka)) ([Miyamoto et al., 2022](https://doi.org/10.5194/gi-11-127-2022); [Nagahara et al., 2022](https://doi.org/10.1186/s40623-022-01609-0)). Just as medical X-rays reveal bones inside a body,
muons penetrate geological structures and their attenuation reveals
density variations within.

While conventional muography with 2--3 viewing directions provides
projected density images, portable nuclear-emulsion detectors now enable
observations from many directions, making detailed 3D density
reconstruction feasible. MUONITH is designed to rapidly evaluate and
compare such multi-directional observation configurations. It is fast
enough to run on a commodity laptop, so researchers can evaluate
observation designs at their own desk — without a dedicated compute
cluster.

## Key features

### Single-view analysis

For each detector viewing direction, MUONITH computes:

- Path-length distributions through a 3D voxel model
- Expected muon counts based on flux models
- Projected average density and its statistical uncertainty

### Spatial resolution vs. depth

Exclusion-plot analysis answers the question:
*"What size and magnitude of density anomaly is detectable at a given depth,
for a given detector area and exposure time?"*

This enables systematic comparison of candidate detector sites and
observation parameters.

### 3D density reconstruction

Given multi-directional muon-count data, MUONITH performs linear inversion
(based on [Nagahara et al., 2022](https://doi.org/10.1186/s40623-022-01609-0))
to reconstruct 3D bulk-density distributions. The tool also exports the
system matrix **A**, prior density vector, and observed counts, allowing
users to apply their own reconstruction methods.

## Why MUONITH?

Existing muography simulation tools fall into two categories:

| Approach | Strengths | Limitations |
|---|---|---|
| **Script-based tools**  | Flexible, user-friendly | Too slow for large parameter scans |
| **Monte Carlo tools**  | Detailed physics and noise modeling | Computationally heavy |

MUONITH fills the gap as a **fast approximate engine** for broad exploration
of observation designs:

- **~100x faster (/thread)** and **~1/4 peak memory** compared to a previous
  GNU Octave implementation for the same 3D inversion problem
  ([Nagahara et al., 2022](https://doi.org/10.1186/s40623-022-01609-0))
- Fast enough to explore **many observation-design variations** on a single
  machine — down to a commodity laptop — rather than a dedicated cluster
- Best candidates can then be validated with full Monte Carlo codes

### Technical highlights

- C++20 with a simple ray-tracing algorithm
- OpenMP parallelization for independent rays and detector elements
- Platform-native BLAS/LAPACK (Apple Accelerate on macOS, OpenBLAS on Linux)
- Sparse matrix storage to keep memory usage low with large 3D grids
- Reproducible builds via Nix development environment

### Platform support

| Platform | Status |
|---|---|
| Linux (Ubuntu 22.04 / 24.04, including WSL2) | Fully supported |
| Linux (Docker) | Fully supported |
| macOS (Apple Clang + Accelerate) | Fully supported |
| macOS (Nix environment) | Fully supported |

## Getting started

- [Installation](getting-started/installation.md) -- Build dependencies and
  compilation instructions
- [Nix Setup](getting-started/nix-setup.md) -- Optional: isolated reproducible
  development environment
- [Quick Start](getting-started/quickstart.md) -- Run your first simulation
