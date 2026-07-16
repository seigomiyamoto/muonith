# Range Table (CSDA)

The range branch of the [build pipeline](index.md) provides the muon range-energy
relation: how much matter (in kg/m²) a muon of a given kinetic energy can penetrate.
The builder inverts this relation to obtain the cutoff energy E_cut(R) during the
flux integration.

## The shipped table

`fluxtable/daemon_groom/logRangeKgm2_logGeV_hokan_70001.txt`:

| Property | Value |
|----------|-------|
| Column 1 | log10(kinetic energy [GeV]), −2.0000 to 5.0000, step 10⁻⁴ (70001 rows) |
| Column 2 | log10(CSDA range [kg/m²]) |
| Material | standard rock (⟨Z/A⟩ = 0.5, ρ = 2.65 g/cm³, I = 136.4 eV) |
| Source | PDG muon table `muE_standard_rock.txt` (Groom-Mokhov-Striganov) |

CSDA (Continuous Slowing Down Approximation) range is the deterministic path length
obtained by integrating the mean energy loss; at TeV energies the actual penetration
depth fluctuates because radiative losses are stochastic, and the CSDA value is used
here as the deterministic representative.

The lower bound of 10 MeV (log10 KE = −2) follows the PDG file itself, whose header
warns that results below 10 MeV are not dependable.

## Source data

The Particle Data Group publishes muon energy-loss tables for about 350 materials in
a common ASCII format:

```
https://pdg.lbl.gov/2024/AtomicNuclearProperties/MUE/muE_<material>.txt
```

Column 1 is the kinetic energy T [MeV] and column 9 is the CSDA range [g/cm²]
(non-numeric header/comment lines interleave the data and must be skipped).
The files carry no explicit license; cite Groom, Mokhov, Striganov,
*Atomic Data and Nuclear Data Tables* **78**, 183 (2001) and the PDG page when
redistributing derived tables. The raw PDG file itself is therefore **not** committed
to this repository; fetch it with the script below.

## Regenerating the table

Both scripts live in `fluxtable_build/csda_range/`:

```bash
cd fluxtable_build/csda_range
bash fetch_pdg_muE.sh                 # downloads muE_standard_rock.txt (md5-checked)
python3 make_range_table.py muE_standard_rock.txt out_logRangeKgm2_logGeV_70001.txt
```

`make_range_table.py` interpolates log10(CSDA range) versus log10(KE) with a
node-exact cubic spline, samples log10 KE = −2..5 in 10⁻⁴ steps (70001 rows), and
converts g/cm² to kg/m² (+1 in log10). The output format is identical to the shipped
table.

To compare a regenerated table against the shipped one:

```bash
python3 make_range_table.py muE_standard_rock.txt out.txt \
  --compare ../../fluxtable/daemon_groom/logRangeKgm2_logGeV_hokan_70001.txt
```

Expected agreement with the shipped table (which was produced from an earlier PDG
release with the same method): |Δ| < 10⁻³ dex above 1 GeV, growing to about
6×10⁻³ dex at the 10 MeV end.

## Other materials

Pass any other PDG file (e.g. `muE_water_liquid.txt`, `muE_iron_Fe.txt`) to
`make_range_table.py` to build a range table for that material, then follow
[Build and Install](build-and-install.md) with the new table. For the physics impact
of composition differences on muography, see Lechmann et al. (2018),
*The effect of rock composition on muon tomography measurements*, Solid Earth 9, 1517,
[doi:10.5194/se-9-1517-2018](https://doi.org/10.5194/se-9-1517-2018).
