# Appendix

This appendix collects implementation-level formulas that are useful when
checking MUONITH outputs against the code. It is intentionally narrower than the
concept pages: use the concept pages for the workflow, and this appendix for the
exact equations and source-code correspondence.

## Efficiency Uncertainty in $\mathbf{C}_N$

MUONITH represents the observation covariance $\mathbf{C}_N$ as a diagonal
matrix. With efficiency uncertainty disabled, the diagonal element for grouped
bin $g$ is the Poisson variance of the grouped muon count:

$$
(\mathbf{C}_N)_{gg} = N_g .
$$

When `tf_eff_cn_diag` is enabled, detector-efficiency uncertainty is added to
that diagonal entry:

$$
(\mathbf{C}_N)_{gg} = N_g + \sigma_{\varepsilon,g}^{2}.
$$

Let $k \in g$ be the angular elements merged into grouped bin $g$. For angular
element $k$, the efficiency table gives lower, central, and upper efficiencies
$\varepsilon_{\mathrm{low},k}$, $\varepsilon_k$, and
$\varepsilon_{\mathrm{upp},k}$. MUONITH uses the band half-width

$$
\delta\varepsilon_k =
\frac{1}{2}\left(\varepsilon_{\mathrm{upp},k}
- \varepsilon_{\mathrm{low},k}\right)
$$

as the efficiency uncertainty. With $b_k$ denoting the efficiency-free base
count of element $k$, the central expected count of that element is
$N_k=\varepsilon_k b_k$. In the grouped-relative notation used in the design
notes,

$$
s_k = \frac{\delta\varepsilon_k}{\varepsilon_k}, \qquad
N_g = \sum_{k\in g} N_k, \qquad
s_g = \frac{\sum_{k\in g}s_k N_k}{N_g}.
$$

Because $s_k N_k = \delta\varepsilon_k b_k$, the default grouped efficiency
variance can be written equivalently as

$$
\sigma_{\varepsilon,g}^{2}
= (s_g N_g)^2
= \left(\sum_{k\in g}\delta\varepsilon_k b_k\right)^2 .
$$

This is the correlated-in-group form: all angular elements in the grouped bin
move in the same direction. If `tf_eff_cn_diag_independent` is enabled,
MUONITH instead uses the independent-element form:

$$
\sigma_{\varepsilon,g}^{2}
= \sum_{k\in g}\left(\delta\varepsilon_k b_k\right)^2 .
$$

The same variance is also used when computing the grouped 2D projected-density
error band:

$$
\sigma_g = \sqrt{N_g + \sigma_{\varepsilon,g}^{2}} .
$$

## Why $\mathbf{C}_N$ Stays Diagonal

The diagonal treatment is a modeling choice, not a statement that correlated
efficiency errors cannot exist. The efficiency table gives each angular
element's lower, central, and upper efficiencies, so it defines the width of the
element-level uncertainty. It does not define correlations between angular bins,
detectors, or calibration-source groups.

MUONITH therefore treats the uncertainty consumed by `tf_eff_cn_diag` as closed
within each grouped angular bin. Correlated components inside the grouped bin
are folded into the diagonal variance term above. Cross-bin or cross-detector
common calibration errors would require off-diagonal covariance terms, and are
outside the current fast evaluation model. This keeps the implementation aligned
with MUONITH's main use: comparing observation layouts and relative density
contrast before data taking. Absolute-density studies that rely on
detector-wide common calibration should provide a separate correlation model or
evaluate that systematic outside the current diagonal $\mathbf{C}_N$ path.

## Source-Code Mapping

| Quantity | Code location | Meaning |
|---|---|---|
| $\delta\varepsilon_k$ | `src/cls_DetectorPanel.cpp`, `calc_var_eff_group` | Half-width of `eff_low` and `eff_upp` |
| $b_k$ | `src/cls_DetectorPanel.cpp`, `calc_var_eff_group` | Efficiency-free count from `DetectorElement::calc_signal()` |
| $\sigma_{\varepsilon,g}^{2}$ | `src/cls_DetectorPanel.cpp`, `calc_var_eff_group` | Grouped efficiency variance |
| $(\mathbf{C}_N)_{gg}$ Poisson term | `src/cls_NagaInv.cpp`, `mp_build_mat_cov_muon` | Diagonal Poisson variance |
| Efficiency variance addition | `src/cls_NagaInv.cpp`, `add_var_eff_to_cov_muon` | Adds `vecxf_var_eff(g)` to the diagonal |
| 2D projected-density band | `src/cls_DetectorPanel.cpp`, `calc_set_proj_dens_grouped` | Uses `sqrt(N_g + var_eff_g)` |
| Main gate | `include/cls_NagaInvParameters.hpp`, `get_tf_eff_cn_diag` | Enables the efficiency variance term |
| Independent-mode gate | `include/cls_NagaInvParameters.hpp`, `get_tf_eff_cn_diag_independent` | Switches from square-of-sum to sum-of-squares |

## Scope

This term is an additional uncertainty on observed counts. It does not replace
Poisson statistics, and it does not introduce a dense covariance matrix.
Correlated errors beyond a grouped bin are outside the current diagonal
implementation path described above.

## 3D Density Reconstruction

This section collects the linearized observation-matrix derivation, the MAP
posterior equations, the notation list, and the dense-covariance memory limit
behind MUONITH's 3D density reconstruction. The workflow and worked example are
in [3D Inversion](../concepts/inversion.md); this appendix is the equation and
source-code reference. The formulation follows Nagahara et al. (2022).

### Observation Matrix

Let $g$ index the grouped angular bins (the observation rows, the same grouped
bins used in the efficiency section above) and $j$ index the voxels. The density
length along grouped bin $g$ is linear in the voxel density vector
$\boldsymbol{\rho}$ (dimension $n_v$):

$$
X_g = \sum_j L_{gj}\, \rho_j \qquad [\mathrm{kg/m^2}],
$$

where $L_{gj}$ is the geometric path length [m] of bin $g$ through voxel $j$,
computed by the DDA ray tracer. The expected muon count is a nonlinear function
of the density length, $N_g = N_g(X_g)$, obtained by folding the CSDA flux table
with effective area, solid angle, and exposure time. Linearizing around the
prior density length $X_{0,g} = \sum_j L_{gj}\,\rho_{0,j}$ gives the observation
(muon-count sensitivity) matrix $\mathbf{A}$:

$$
A_{gj} = \frac{\partial N_g}{\partial \rho_j}
\approx L_{gj} \left.\frac{dN_g}{dX_g}\right|_{X_{0,g}},
$$

which is Eq. (18) of Nagahara et al. (2022). The derivative $dN_g/dX_g$ is the
differential sensitivity of the penetrating-muon flux, evaluated from finite
differences of the CSDA table. Because $\mathbf{A}$ is linear in $L$, it is
linear under the adaptive bin merging (see the efficiency section above): the
observation matrix of a merged bin is the plain sum of the rows of its sub-bins,
so merging does not degrade numerical accuracy.

### MAP Posterior

With a Gaussian prior [Tarantola, 2005], following the joint muon–gravity
formulation of Nishiyama et al. (2014b), the posterior covariance and the
posterior density are

$$
\mathbf{C}_{\rho'} = \left(\mathbf{A}^\top \mathbf{C}_N^{-1} \mathbf{A}
+ \mathbf{C}_\rho^{-1}\right)^{-1},
$$

$$
\delta\boldsymbol{\rho} = \mathbf{C}_{\rho'}\, \mathbf{A}^\top \mathbf{C}_N^{-1}
\left(\mathbf{N}_{\text{obs}} - \mathbf{N}_0\right),
$$

$$
\boldsymbol{\rho}' = \boldsymbol{\rho}_0 + \delta\boldsymbol{\rho}.
$$

The prior density covariance uses an exponential (isotropic by default)
correlation,

$$
C_{\rho,jj'} = \sigma_\rho^2 \exp\!\left(-\frac{d_{jj'}}{l_c}\right),
$$

with density uncertainty $\sigma_\rho$ [kg/m³] and correlation length $l_c$ [m].
MUONITH also implements a separable anisotropic form with independent horizontal
and vertical correlation lengths.

### Notation

| Symbol | Dimensions | Meaning | Role |
|---|---|---|---|
| $\mathbf{L}$ | $n_d \times n_v$ | Path-length matrix; $L_{gj}$ is the geometric path length [m] of bin $g$ through voxel $j$ (DDA ray tracing) | Input |
| $\boldsymbol{\rho}_0$ | $n_v \times 1$ | Prior density [kg/m³] | Input |
| $\sigma_\rho$ | scalar | Density uncertainty [kg/m³] | Input |
| $l_c$ | scalar | Spatial correlation length [m] (horizontal/vertical can be set independently) | Input |
| $\mathbf{N}_{\text{obs}}$ | $n_d \times 1$ | Observed muon counts | Input |
| $\mathbf{C}_N$ | $n_d \times n_d$ | Muon-count error covariance (diagonal; Poisson term, optionally plus efficiency variance) | Input |
| $\left.dN_g/dX_g\right\rvert_{X_{0,g}}$ | $n_d \times 1$ | Differential muon-flux sensitivity at the prior density length, precomputed from the flux table | Input |
| $X_g$ | scalar | Density length of bin $g$ [kg/m²] | Intermediate |
| $\mathbf{A}$ | $n_d \times n_v$ | Observation matrix; $A_{gj} = \partial N_g/\partial\rho_j \approx L_{gj}\,(dN_g/dX_g)\rvert_{X_{0,g}}$ | Intermediate |
| $\mathbf{C}_\rho$ | $n_v \times n_v$ | Prior density covariance (smoothness constraint) | Intermediate |
| $\mathbf{N}_0$ | $n_d \times 1$ | Expected muon counts at the prior model | Intermediate |
| $\boldsymbol{\rho}'$ | $n_v \times 1$ | Reconstructed density (posterior mean) | Output |
| $\mathbf{C}_{\rho'}$ | $n_v \times n_v$ | Posterior density covariance | Output |

Here $n_d$ is the number of detector elements (grouped observation rows) and
$n_v$ is the number of voxels.

### Numerical Constraints

The posterior covariance inverse is solved by LU decomposition
(`LAPACKE_sgetrf` + `LAPACKE_sgetri`). Single and double precision are both
available; single precision is the default, with double precision recommended
for ill-conditioned problems. The dense covariance matrix scales as $O(n_v^2)$
in memory, so `MAX_VOXEL_FOR_DENSE_COV = 50000` (about 9.3 GB for a float
matrix) is the safe limit. Beyond it, coarsen the resolution or restrict the
reconstruction sub-volume.

### Source-Code Mapping

| Quantity | Code location | Meaning |
|---|---|---|
| $\mathbf{L}$ | `vec_spmat_PL` (Trace Path Lengths, Module 4) | Density-independent path-length matrices |
| $dN_g/dX_g$, $\mathbf{N}_0$ | Compute Prior (Module 5) | Prior density length and expected counts |
| $\mathbf{A}$ | Build Observation Matrix (Module 6), `calc_dNdD` namespace | Observation matrix assembly |
| $\boldsymbol{\rho}'$, $\mathbf{C}_{\rho'}$ | Invert Density (Module 7) | Density reconstruction (LU via LAPACK) |
