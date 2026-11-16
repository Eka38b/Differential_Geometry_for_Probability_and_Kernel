# Survey 01 methodology

## Objective and admissible data

Study fixed-data fibers of the MathOverflow problem in the specialization $S=0$. The input matrix $C$ is orthogonal and fixes $e_{1}$. Both diagonal lists and both strictly positive superdiagonal lists are frozen. Only strictly lower-triangular entries are unknown.

The retained 69 inputs have dimensions between 2 and 20. They are selected compatible families, not a uniform or probabilistic sample of all possible fixed data. The inputs were originally generated and searched on 7 November 2026. The present export filters the historical records to the public problem's assumptions and recalculates the displayed totals.

## Witness construction

| Family | Retained inputs | Construction and purpose |
| --- | ---: | --- |
| Regular controls | 9 | Skew-tridiagonal normal matrices with positive superdiagonals and $C=I$; affine uniqueness controls |
| Cyclic controls | 13 | A cyclic permutation matrix, optional scalar shift, and a signed reversal fixing the first coordinate; singular Jacobians with exact isolation from cut budgets |
| Four-dimensional controls | 2 | The explicit $F(t)$ family in the problem note, with $a=3,b=5$, optionally shifted; exactly two solutions |
| Central spectral inputs | 12 | Real normal spectral blocks paired by symmetry, followed by Arnoldi reduction; avoid relying only on permutation overlap matrices |
| Mixed spectral inputs | 12 | Varied real/complex normal spectral blocks; separate Arnoldi reductions construct the dual Hessenberg form |
| Coupled zero-pivot inputs | 21 | Signed reversal on coordinates 2 through $n$, zero diagonal, and a jointly constructed normal witness with no zero superdiagonal |

For coupled inputs, $n=4,6,\ldots,16$ and three target values $a_{2}=1.05,1.15,1.30$ are used with $a_{1}=1$. The construction first solves the linear relation $D^{T}+CDC^{T}=0$ in an eligible Hessenberg basis. It then fits the quadratic normality equations together with the two target-superdiagonal equations. Simultaneous diagonal sign normalization makes the witness superdiagonals positive. Its resulting scalar lists and ordered $C$ are then frozen before the fixed-data root search.

The auxiliary witness fit is input generation; it is not counted as a root-search call with the scalar data already fixed. Spectral starting weights, seeds, target values, and construction diagnostics that were recorded are retained in the JSON metadata. Data files preserve the actual matrices rather than requiring rounded values to be retyped.

## Linear elimination

Eliminate $E=-C^{T}D^{T}C$. Write $D$ as its fixed upper part plus a linear combination of all strictly lower matrix units. Matching the prescribed upper part of $E$ gives a linear system. An SVD computes a particular compatible matrix and a Frobenius-orthonormal nullspace basis:

$$
D(z)=D_{0}+\sum_{\ell=1}^{r}z_{\ell}X_{\ell}.
$$

No lower coordinate is removed merely because it appears zero in a witness. The historical linear relative threshold was $10^{-11}$, audited at $10^{-9}$ and $10^{-13}$. The remaining quadratic system is normality of $D$; the dual relation makes normality of $E$ equivalent, and both are independently validated.

## Historical nonlinear search

The historical campaign used trust-region least squares with analytic Jacobians. A feasible witness was polished first. The search then used independent starts, with the budget capped at 64 for the larger coupled inputs; smaller active cases used up to 128 starts. The retained records contain 2,176 multistart calls in total.

Starts used the ellipsoid implied by the weighted cut identity when available:

$$
\sum_{i>j}(i-j)D_{ij}^{2}=\sum_{k=1}^{n-1}a_{k}^{2}.
$$

Substitution of the affine parameterization gives a quadratic ellipsoid in $z$. A Cholesky factor maps random unit directions to this ellipsoid; some starts are scaled radially. If the radius is numerically negligible, Gaussian parameter starts are used instead.

The historical residual is $D^{T}D-DD^{T}$, packed on and above the diagonal with off-diagonal weights $\sqrt{2}$ and normalized by matrix scale. The shared replay solver uses its negative, which has the same zeros. Its analytic derivative is obtained by differentiating each matrix product; this avoids finite-difference derivatives inside the solver.

## Validation, rank, and controls

Historical root acceptance required relative normality and Hessenberg residuals below $10^{-10}$ and prescribed-entry discrepancy below $10^{-8}$. Root clustering used $2\times10^{-5}$ times matrix scale. Jacobian rank used relative threshold $10^{-7}$ and absolute floor $10^{-10}$, with neighboring thresholds checked by the original audit.

Analytically known control roots were validated separately after numerical recovery was recorded. The data distinguish supplied witnesses, multistart roots, cut-budget controls, and analytic controls. Root occurrences are counted separately for different fixed inputs; they do not count individual least-squares convergences.

The retained cyclic examples include singular isolated roots. Their first-column bottom entry exhausts every cut budget, so all other strictly lower entries vanish. This exact argument establishes isolation independently of a singular Jacobian. No predictor/corrector continuation attempts were recorded in this retained subset: the saved continuation count is zero. It would be incorrect to describe it as an independent successful search against continuous motion.

## Included implementation and replay

`Code/fiber.py` implements the full fixed-data equations and analytic Jacobian. `Code/replay.py`, called by this directory's `run.py`, reads each frozen input, recomputes its full affine reduction, alternates trust-region and Levenberg–Marquardt starts, and writes new roots and diagnostics to the requested output directory.

The modern replay has explicit new initialization and clustering settings: Gaussian starts at matrix scale, default 24 calls for positive affine dimension, full-equation acceptance $2\times10^{-10}$, and clustering $2\times10^{-6}$. These are documented changes, not a claim to reproduce the historical trajectory or exact recovery count. Use the archived JSON for the original budgets and origins.

## Saved data and interpretation

[The data guide](../../Data/01_zero_source/README.md) lists the NPZ fields. The summary gives dimension, affine dimension, saved roots, calls, and validation error. The JSON preserves more detailed historical ranks, origins, and construction metadata.

The largest affine dimension is 171, largely from deliberately singular controls. All recovered coupled roots were locally full rank; 11 saved root occurrences in other controls have rank-deficient normality Jacobians. Neither fact proves the unrestricted finiteness statement. Generic inputs that already have zero affine dimension chiefly test the linear machinery, not the ability of normality to remove large fibers.
