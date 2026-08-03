# MUONITH

**MUOgraphic Numerical Inversion & Tomography Harness**

A fast C++ simulation tool for multi-directional muography observation design
and 3D density reconstruction.

**Manual: <https://seigomiyamoto.github.io/muonith/>** -- the same pages are in
`docs/` and can be read offline.

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
optimize such multi-directional observation configurations.

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
- Scan **thousands of parameter configurations** within half a day
  on a single desktop machine
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

- [Installation](docs/getting-started/installation.md) -- Build dependencies and
  compilation instructions
- [Nix Setup](docs/getting-started/nix-setup.md) -- Optional: isolated reproducible
  development environment
- [Quick Start](docs/getting-started/quickstart.md) -- Run your first simulation

Once the build succeeds, run the bundled synthetic `tutorial` station in two
steps: generate the working directory, then run the reconstruction and plots:

```bash
bash setup_station.sh tutorial
cd work/tutorial/swp001
bash run_prg.sh true 42
```

## Related tools

MUONITH is developed alongside two standalone companion tools:

- [**muonith-gsi-dem**](https://github.com/seigomiyamoto/muonith-gsi-dem) --
  converts Geospatial Information Authority of Japan (GSI) DEM ZIP archives into
  the `.g2zbin` grids MUONITH reads.
- [**muonith-path-view**](https://github.com/seigomiyamoto/muonith-path-view) --
  visualizes detector viewsheds (line-of-sight paths) over a DEM.

Each tool also runs on its own; the bundled `tutorial` example needs neither of
them.

## Citation

If you use MUONITH in academic work, please cite the software. Archived releases
are deposited on Zenodo:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21694789.svg)](https://doi.org/10.5281/zenodo.21694789)

The DOI above always resolves to the latest release. To cite the exact version
you used, take the version-specific DOI from that Zenodo record; for v1.0.0 it is
[10.5281/zenodo.21694790](https://doi.org/10.5281/zenodo.21694790). Machine-readable
metadata is in [CITATION.cff](CITATION.cff).

## License

MUONITH source code and project-authored documentation are licensed under the
[BSD 3-Clause License](LICENSE).

Third-party libraries, flux tables, terrain data, detector example data, and
figures may be governed by separate licenses or attribution requirements. See
[NOTICE](NOTICE) and [LICENSES/](LICENSES/) for details.
