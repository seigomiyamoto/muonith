# Efficiency Table (text)

The table route is selected by the `path_eff_table` key of the detector parameter file, and only when `path_eff_model` is `"none"` or absent (see [Overview](index.md#two-ways-to-provide-the-efficiency)).

Use it when the efficiency cannot be written as a formula. The [efficiency model](efficiency-model.md) is a separable product `base_eff * g_tx(tx) * g_ty(ty)`, so an arbitrary per-bin efficiency map — a measured one, for instance — generally cannot be expressed as a model. The table is the only input that accepts one.

## Format

A whitespace-separated text file with one row per angular bin and seven columns:

```
tx_min tx_max ty_min ty_max eff_low eff_cnt eff_upp
```

Every bin of the detector grid gets its own row. The rows must match that grid exactly: a row whose bin center falls outside the detector's angular range is dropped with an `out of range` warning (the first five are logged), and a table generated for one binning cannot be reused for another. The model, by contrast, is evaluated on the detector's own grid and cannot go out of alignment.

The last three columns are the band of the bin. `eff_cnt` scales the expected muon count; the half-width implied by `eff_low` and `eff_upp` propagates into the observation covariance $\mathbf{C}_N$ when `tf_eff_cn_diag` is enabled. See [Efficiency Uncertainty in $\mathbf{C}_N$](../../reference/appendix.md#efficiency-uncertainty-in-mathbfc_n) for that derivation. Writing `eff_low` = `eff_cnt` = `eff_upp` gives a zero-width band, which adds nothing to $\mathbf{C}_N$.

## Using a User-Defined Table

1. Read the angular grid off the detector parameter file you intend to use: `nbinx`, `txmin`, `txmax`, `nbiny`, `tymin`, `tymax`. Your table must have exactly `nbinx * nbiny` rows covering that rectangle with uniform bins.

2. Write the table to `data/_shared/eff_table/<your_table>.txt`. One row per bin, seven numeric columns, whitespace-separated. Blank lines are skipped; **every other line must parse as seven numbers**, so do not add `#` comment lines or a header row — the reader stops with a format error on the first line it cannot parse.

    ```
    -1.60 -1.59 0.00 0.01 0.93 0.95 0.97
    -1.59 -1.58 0.00 0.01 0.93 0.95 0.97
    ...
    ```

    Column order within a row is fixed; row order is free, because each row is placed by its own bin center.

3. Point the detector parameter file at the table **and disable the model**, which otherwise takes precedence:

    ```json5
    , "path_eff_table" : "../detparams/eff_table/<your_table>.txt"
    , "path_eff_model" : "none"
    ```

    Omitting the second line is the most common mistake: a runcard that still names a model file reads the model and never opens your table.

4. Run `setup_station.sh <site>`. It links `data/_shared/eff_table/` into `work/<site>/detparams/eff_table/`, so the relative path above resolves from the run directory, and it then checks your table's grid against every detector parameter file. A mismatch in `nbinx`, `nbiny`, `txmin`, `txmax`, `tymin`, or `tymax` stops the setup with `Detector angle grid and eff_table grid differ`, before any physics runs.

5. Run as usual. The log confirms the route with `Set efficiency for N elements from <path>`, where `N` equals `nbinx * nbiny` when the grids agree.

`scripts/make_eff_table.py` writes tables in this format from a model file; use it as a worked example of the layout, or as a starting point for a script that emits your measured values.
