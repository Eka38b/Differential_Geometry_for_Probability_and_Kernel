# Survey 03 methodology

## Objective

Test genuinely nonzero sources in the public MathOverflow problem while preserving the geometric shared-tangent condition $Se_{1}=0$. Both positive superdiagonal lists and both diagonal lists are fixed before solving. This is a deliberately constructed compatible family, not a sample of all possible sources or all possible normal pairs.

## Compatible witness construction

For each dimension $n=4,5,6,8,10$, choose positive coefficients independently and uniformly from $[0.5,1.5]$. Define the skew-tridiagonal witness

$$
D=E=J(a),\qquad J_{i,i+1}=a_{i},\quad J_{i+1,i}=-a_{i},\quad J_{ii}=0.
$$

Choose $C=I_{2}\oplus Q$ and set

$$
S=-J(a)+CJ(a)C^{T}.
$$

Then $D,E$ are normal, all superdiagonals are positive, and the dual affine relation holds. Since $Ce_{1}=e_{1}$ and $Ce_{2}=e_{2}$, while $J(a)e_{1}=-a_{1}e_{2}$, one has $Se_{1}=0$ exactly at the symbolic construction level. Floating-point evaluation preserves it to numerical precision. Zero sources are explicitly rejected.

Four overlap families are used, with two seeds in each dimension:

| Family | Choice of $Q$ | Purpose |
| --- | --- | --- |
| Reflection | $I-2vv^{T}$ for a normalized random $v$ | Structured compatible source |
| Zero-pivot reflection | Same, with the first component of $v$ chosen to equal the norm of the rest before normalization | Force the first diagonal entry of $Q$ to zero |
| Swap | Exchange the first and last coordinates of the complementary block | Insert exact zero pivots and retain some affine freedom |
| Orthogonal rotation family | QR-normalized random Gaussian matrix | Avoid only symmetric or permutation overlaps |

There are 40 retained inputs, ten per overlap family. The seed is a deterministic function of dimension, family, and replicate, recorded per case. The generator does not independently vary diagonal lists: these witnesses have zero diagonals. It also supplies equal superdiagonal lists. These limitations restrict the conclusions to the chosen family.

## Freeze data and eliminate linear constraints

After constructing $C,S,D,E$, freeze the complete input. Eliminate

$$
E(D)=C^{T}(S-D^{T})C.
$$

Perturb every strictly lower entry of $D$. Matching the full prescribed upper part of $E(D)$ gives linear constraints. Append the linear commutation equation

$$
[S,D+D^{T}]=0.
$$

Under the affine relation, this equation together with unstarred normality is equivalent to normality of both matrices. The SVD nullspace uses threshold $10^{-10}\max(1,\sigma_{1})$; ranks are also recorded at relative scales $10^{-8}$ and $10^{-12}$.

The resulting affine dimensions are zero for 34 inputs. The six positive-dimensional cases are the swap inputs in dimensions 6, 8, and 10, with dimensions 1, 3, and 5 respectively. Reconstructing the affine space from the full equations prevents the witness's skew symmetry from being imposed on unknown solutions.

## Nonlinear algorithms

The residual is the upper triangle of $DD^{T}-D^{T}D$, with off-diagonal weights $\sqrt{2}$, divided by $\max(1,\|D_{0}\|^{2})$. The derivative in direction $X$ is

$$
XD^{T}+DX^{T}-X^{T}D-D^{T}X.
$$

The supplied witness is accepted separately and labeled as such. For each positive-dimensional input, run 24 independent Gaussian starts in affine coordinates. Scale the start standard deviation by $\max(1,\|D_{0}\|)/\sqrt{r}$. Alternate trust-region reflective and Levenberg–Marquardt least squares, with analytic Jacobians, at most 600 function evaluations per call, and termination tolerances $3\times10^{-13}$.

No nonlinear starts are spent on the zero-dimensional affine cases. The full run therefore contains 144 local solver calls. The root count includes the initial witness, and repeated convergence to it is not a new solution.

## Acceptance and isolation diagnostics

Validate every candidate against both normality equations, the original affine relation, all prescribed upper entries, both Hessenberg zero patterns, orthogonality, $Ce_{1}=e_{1}$, skewness of $S$, and $Se_{1}=0$. Each scaled residual must be at most $2\times10^{-10}$. Positive superdiagonals are required with a $10^{-8}$ numerical lower threshold.

Cluster matrices within $2\times10^{-6}$ times matrix scale. At saved roots, compute the normality Jacobian restricted to the affine coordinates. Numerical rank uses relative threshold $10^{-7}$ and absolute floor $10^{-10}$. Save its rank, nullity, and smallest-to-largest singular-value ratio.

This run found no singular saved root in the positive-dimensional affine cases. Thus no predictor/corrector continuation probes were initiated or claimed. No interval proof of isolation is supplied for these newly generated source inputs.

## Outputs and scope of the conclusion

The script saves `C`, `S`, diagonal and superdiagonal lists, the witness-centered affine basis, roots and dual roots, and linear singular values in each NPZ. JSON records the seeds, family, source norms, leading-minor diagnostics, solver calls, origins, rank sweeps, and full-equation errors. The CSV is a compact index.

One root per input was stored, namely the supplied witness up to numerical recovery. No additional root was found in the six searched affine fibers. The data test nonzero shared-tangent sources and some zero-pivot overlaps, but do not rule out undiscovered solutions, continuous components elsewhere, or different behavior in nonskew witness families. They also do not claim realization in a fixed Fisher statistical model.
