# Third-Party License Inventory

This directory records third-party license information for MUONITH. The
top-level `LICENSE` applies to MUONITH-authored source code and documentation;
it does not relicense third-party libraries, flux tables, terrain data, or
figures.

| Component | Use | License / terms to verify before redistribution |
| --- | --- | --- |
| Eigen 3.4.0 | Header-only linear algebra | Mainly MPL-2.0; Eigen also contains some files under other licenses. Consider building with `EIGEN_MPL2_ONLY` for public releases. |
| nlohmann/json 3.11.3 | JSON parser | MIT |
| fmt 11.2.0 | Formatting library | MIT |
| spdlog 1.15.3 | Logging library | MIT |
| OpenMP runtime | Parallel execution | Runtime-dependent; GCC libgomp uses GPL with GCC Runtime Library Exception, LLVM OpenMP uses Apache-2.0 WITH LLVM-exception. |
| BLAS/LAPACK/OpenBLAS/Accelerate | Linear algebra backend | Backend-dependent; OpenBLAS and reference LAPACK are BSD-style, Accelerate is supplied by macOS. |
| daemonflux-derived flux data | Flux tables or derived tables | daemonflux BSD-3-Clause; keep copyright notice, license terms, and disclaimer. |
| Groom/PDG range data | Muon range table inputs or derivatives | Attribute Groom et al. and PDG Review of Particle Physics; record interpolation or re-gridding. |
| MathJax 3.2.2 | Equation rendering in the built documentation site; vendored under `docs/assets/mathjax/` | Apache-2.0; the es5 build embeds the Speech Rule Engine, also Apache-2.0. |

## Full license texts in this directory

| File | Component |
| --- | --- |
| `Eigen-MPL-2.0.txt` | Eigen (MUONITH builds with `EIGEN_MPL2_ONLY`, so only the MPL-2.0 part is used) |
| `nlohmann_json-MIT.txt` | nlohmann/json |
| `fmt-MIT.txt` | fmt |
| `spdlog-MIT.txt` | spdlog |
| `OpenBLAS-BSD-3-Clause.txt` | OpenBLAS (Linux BLAS/LAPACK backend) |
| `daemonflux-BSD-3-Clause.txt` | daemonflux (source of the shipped flux tables) |
| `MathJax-Apache-2.0.txt` | MathJax (vendored in the documentation site, together with its embedded Speech Rule Engine) |

The OpenMP runtime and macOS Accelerate are supplied by the toolchain or the
operating system and are not redistributed here.

Before publishing a source archive, binary package, container image, or public
data bundle, copy the exact license texts for the redistributed third-party
materials into this directory or a generated notice bundle.
