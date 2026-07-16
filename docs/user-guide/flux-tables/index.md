# Making Penetrating Muon Flux Tables

MUONITH ships pre-computed muon flux tables (see [Input Data](../input-data.md#muon-flux-tables)
for how the pipeline consumes them). This section describes how those tables are made,
how to reproduce the shipped tables, and how to introduce a different flux model or a
different absorber material.

!!! note "Standard-rock assumption"
    The shipped tables and this entire walkthrough assume **standard rock**
    (⟨Z/A⟩ = 0.5, density 2.65 g/cm³, mean excitation energy I = 136.4 eV) as the
    absorber composition. Density itself is *not* assumed: the pipeline works in
    density-length R = density × path length [kg/m²], so local density enters through
    the voxel densities. Only the *composition* per unit mass is fixed. For media whose
    composition differs substantially from standard rock (water, ice, iron, ...),
    regenerate the range table for that material (see [Range Table](range-table.md))
    and consult Lechmann et al. (2018), *The effect of rock composition on muon
    tomography measurements*, Solid Earth 9, 1517,
    [doi:10.5194/se-9-1517-2018](https://doi.org/10.5194/se-9-1517-2018).

## What is produced

Two tables are produced, both as `.g2zbin` binary grids (plus `.txt` dumps):

| Table | Content | Axes | Scale |
|-------|---------|------|-------|
| Penetrating flux F | $\log_{10} F$ [1/(m² sr s)] | cos(θ_zenith) × density-length R [kg/m²] | log10 |
| Flux derivative dF/dR | $dF/dR$ [1/(m² sr s)/(kg/m²)] | density-length R × cos(θ_zenith) | linear (all values ≤ 0) |

## Build pipeline

Two independent input branches are prepared and then joined by a single builder:

![Flux table build pipeline](../../assets/images/flux_table_pipeline.drawio.png)

- **Range branch (CSDA)** — the muon range-energy table derived from the PDG
  (Groom-Mokhov-Striganov) CSDA range data. Model-independent: every flux model
  reuses the same range table. See [Range Table (CSDA)](range-table.md).
- **Flux-model branch** — the sea-level differential flux dF/dE on a uniform grid.
  The shipped tables use daemonflux; you can bring your own model.
  See [Flux Model (dF/dE)](flux-model.md).
- **Build & install** — `make_peneflux_dFdR.exe` integrates the two inputs into the
  F and dF/dR tables and they are installed under `fluxtable/<flux_groom>/`.
  See [Build and Install](build-and-install.md).

## When you need this

- You want to **reproduce** the shipped `fluxtable/daemon_groom/` tables from source data.
- You want a **different binning** (range axis, cosθ axis) of the same model.
- You want to introduce **another flux model** (your own dF/dE table).
- You want a **different absorber material** (another PDG `muE_<material>.txt`).
