# Detector Efficiency

Detector efficiency is the probability (0 to 1) that a muon crossing the detector is actually counted. A real detector never counts every muon — scintillator response, electronics dead time, and readout thresholds all lose some events, and the loss usually depends on the arrival direction. If this loss is not modeled, the reconstruction interprets "fewer counts" as "more rock", biasing the density result.

MUONITH therefore lets each detector carry an efficiency value per angular bin. Every bin `(tx, ty)` — the same angular grid defined by `nbinx`/`txmin`/`txmax` and `nbiny`/`tymin`/`tymax` in the [detector parameter file](../parameter-files.md#detector-parameter-files) — holds three numbers. (These bins are the computational elements of the path-length calculation, not the instrument's angular resolution; see the note in [Detector Model](../../concepts/detector-model.md#detectorelement).)

| Value | Meaning |
|-------|---------|
| `eff_cnt` | Central (best-estimate) efficiency of the bin |
| `eff_low`, `eff_upp` | Lower / upper bound of the efficiency uncertainty band |

The half-width of the band feeds the observation covariance and the 2D projected-density error band when efficiency modeling is enabled (`tf_eff_cn_diag`; see [Parameter Reference](../../reference/parameter-reference.md#nagainv_parameters-array) and the [Appendix](../../reference/appendix.md#efficiency-uncertainty-in-mathbfc_n) for the exact formulas).

## Two Ways to Provide the Efficiency

The detector parameter file selects the input route with two keys:

| Key | Input | Status |
|-----|-------|--------|
| `path_eff_model` | Analytic efficiency model (one small JSON5 file) | **Preferred** |
| `path_eff_table` | Text table with one row per angular bin | For measured or otherwise non-analytic efficiency |

Rules:

- If `path_eff_model` is set to anything other than `"none"`, the model is used and `path_eff_table` is ignored.
- If `path_eff_model` is `"none"` (or absent), the text table at `path_eff_table` is read.
- If neither is provided, the run stops with a file-open error — one of the two routes is required.

## Which Route to Use

- **[Efficiency Model (JSON5)](efficiency-model.md)** — a ~20-line file of coefficients, evaluated at the center of every bin of the detector's own grid. It cannot go out of alignment with that grid, and one file serves detectors with different binnings. Use it whenever the efficiency can be written as a formula.
- **[Efficiency Table (text)](efficiency-table.md)** — one row per angular bin. The model is a separable product of two axis functions, so an arbitrary per-bin map — a measured one, for instance — generally cannot be expressed as a model. The table is the only input that accepts one.

## Quick Start (Efficiency Model)

1. Write a model file, e.g. `data/<site>/detparams/eff_model/my_model.json5`. The smallest useful model is a constant efficiency with a constant uncertainty:

    ```json5
    {
      "base_eff": 0.95,
      "sigma": { "base": 0.05 }
    }
    ```

    This gives every bin `eff_cnt` = 0.95 and a band of ±0.05.

2. Point the detector parameter file at it:

    ```json5
    "path_eff_model": "../detparams/eff_model/my_model.json5"
    ```

3. Run as usual. `setup_station.sh` links `data/<site>/detparams/eff_model/` into `work/<site>/detparams/eff_model/`, so the relative path above resolves from the run directory (see [File Path Resolution](../parameter-files.md#file-path-resolution)). The log confirms the route with `EfficiencyModel::load: path=...` and `Set efficiency for N elements from the efficiency model`.
