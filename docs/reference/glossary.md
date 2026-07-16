# Glossary

This page is a one-stop index of the abbreviations and symbols used across the
MUONITH documentation. Each entry gives a short meaning and points to the page
that defines it in full. The linear-inversion matrix symbols are collected in
the [Appendix notation table](appendix.md#notation); this page links to them
rather than repeating the table.

## Abbreviations and quantities

| Term | Meaning | Unit | Defined in |
|---|---|---|---|
| PL | Path length: geometric length of a ray through matter | m | [Coordinate System — Units](../concepts/coordinate-system.md#units) |
| DL / density length | Density integrated along the ray path | kg/m² | [Coordinate System — Units](../concepts/coordinate-system.md#units) |
| tx, ty | Angular bin coordinates of a detector element (horizontal, vertical) | tangent / degree / radian (per `angle_unit`) | [Detector Model — From (tx, ty) to a ray direction](../concepts/detector-model.md#from-tx-ty-to-a-ray-direction) |
| uqid (`Uqid`) | Globally unique detector-element index | — | [Coordinate System — Key Types](../concepts/coordinate-system.md#key-types) |
| uqiv (`Uqiv`) | Globally unique (linearized) voxel index | — | [Coordinate System — Key Types](../concepts/coordinate-system.md#key-types) |
| shell / `ShellPL` | Density-independent path length through the terrain shells around the voxel grid (upper / lower / lateral) | m | [Detector Model — Shell path lengths](../concepts/detector-model.md#shell-path-lengths-shellpl) |
| reconst voxel | Voxel solved as an unknown in the inversion | — | [Grid System](../concepts/grid-system.md) |
| non-rec voxel | Voxel not solved as an unknown, but still contributing to ray attenuation through the prior | — | [Grid System](../concepts/grid-system.md) |

## Inversion symbols

The linear-inversion symbols — path-length matrix $\mathbf{L}$, observation
matrix $\mathbf{A}$, observation covariance $\mathbf{C}_N$, prior covariance
$\mathbf{C}_\rho$, posterior covariance $\mathbf{C}_{\rho'}$, prior and
posterior density $\boldsymbol{\rho}_0$ and $\boldsymbol{\rho}'$, density
uncertainty $\sigma_\rho$, correlation length $l_c$, density length $X$, and
muon counts $\mathbf{N}_{\text{obs}}$ and $\mathbf{N}_0$ — are listed with their
dimensions and input / intermediate / output role in the
[Appendix — Notation](appendix.md#notation) table. In code and on some pages the
covariances appear with a suffix, `C_N` and `C_rho`; these denote the same
$\mathbf{C}_N$ and $\mathbf{C}_\rho$.
