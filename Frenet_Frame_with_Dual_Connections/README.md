# Dual Frenet frames: algebra, numerical experiments, and information geometry

Research archive of Yunju Lee (Eka38b). Compiled 16 November 2026.

This project asks whether two orthonormal Frenet frames associated with conjugate affine connections can be reconstructed from their diagonal coefficients, curvature lists, and relative orientation. Its algebraic core is a finiteness problem for pairs of real normal Hessenberg matrices with fixed upper entries. The public problem is [MathOverflow question 516790](https://mathoverflow.net/q/516790), cited below.

The archive records three stages: **A)** the geometric motivation and algebraic formulation; **B)** the numerical approaches already performed, with their code, saved inputs, results, and certificates; **C)** the interpretation of normality through information-geometric alpha connections. The small information-geometry checks dated 16 November are new consistency checks, separate from the earlier matrix surveys.

**Status.** General finiteness remains an open claim in this project. Linear uniqueness under a principal-minor condition is proved. The proposed strengthening that all solutions differ only by signs is false, with certified counterexamples. Neither finite root recovery nor local isolation of recovered roots proves global finiteness.

## Archive contents

| Path | Contents |
| --- | --- |
| [numerics/m_survey/](numerics/m_survey/) | First M.1/M.2 survey: code, fixed inputs, roots, audit, report, and logs |
| [numerics/m_survey2/](numerics/m_survey2/) | Alternative solvers and families; magnitude counterexamples and exact certificates |
| [numerics/mprime_probe/](numerics/mprime_probe/) | Noncyclic M′ continuation probe, saved sources, inputs, and diagnostics |
| [info_geometry/alpha_checks.py](info_geometry/alpha_checks.py) | New checks using a categorical statistical model and a nonmetric normal example |
| [info_geometry/alpha_checks_results.json](info_geometry/alpha_checks_results.json) | Recorded matrices, residuals, and step-size comparisons |
| [notes/Matrix_Problems_M_Mprime.pdf](notes/Matrix_Problems_M_Mprime.pdf) | Earlier settings, named claims/propositions, and proofs |
| [notes/Dual_Frenet_Frames_Discussion_Note.pdf](notes/Dual_Frenet_Frames_Discussion_Note.pdf) | Geometric reconstruction note and the shared-tangent obstruction |
| [references.bib](references.bib) | Supplied MathOverflow citation |

The numerical directories retain their original READMEs and recorded runs. Their historical M′ formulation allows more sources than the geometric formulation below; this distinction matters when interpreting their results. Relative links resolve after extracting the complete archive. The checks made during compilation are recorded in [verification_20261116.json](verification_20261116.json).

## A. From dual Frenet frames to an algebraic problem

### A.1 Geometric starting point

Let $(M,g,\nabla,\nabla^*)$ be a positive-definite dualistic structure:

$$
Xg(Y,Z)=g(\nabla_XY,Z)+g(Y,\nabla_X^*Z).
$$

Along a smooth unit-speed curve $\gamma(s)$, write $T=\dot\gamma$. Construct two full, individually orthonormal frames $N=(N_1,\ldots,N_n)^T$ and $N^*=(N_1^*,\ldots,N_n^*)^T$, sharing $N_1=N_1^*=T$. For the unstarred frame use the projection recursion

$$
W_{i+1}=\nabla_TN_i-\sum_{j=1}^i g(\nabla_TN_i,N_j)N_j,
\qquad a_i=\|W_{i+1}\|_g>0,\qquad N_{i+1}=W_{i+1}/a_i,
$$

and use the analogous recursion for $\nabla^*$. Work on intervals where all required residuals are nonzero. Full frames are essential: cross-Gram matrices between different proper osculating subspaces need not be orthogonal.

Define

$$
D_{ij}=g(\nabla_TN_i,N_j),\quad
E_{ij}=g(\nabla_T^*N_i^*,N_j^*),\quad
C_{ij}=g(N_i,N_j^*).
$$

Here $E$ denotes the earlier $D^*$; the star labels the second connection/frame and never means transpose, adjoint, or complex conjugation. Superscript $T$ means transpose. With these row-index coefficient conventions, $\nabla_TN=DN$.

Projection gives **lower Hessenberg** matrices, explicitly $D_{ij}=E_{ij}=0$ for $j>i+1$. Their superdiagonals are the positive curvature lists $a,a^*$. The diagonal lists $d,d^*$ and the entries below the diagonal are generally nonzero.

For a metric connection, orthonormality gives $D+D^T=0$, so the curvature list determines the skew-tridiagonal Frenet matrix. For general dual connections,

$$
D_{ij}+D_{ji}=-(\nabla_Tg)(N_i,N_j).
$$

The missing symmetric part is therefore the main reconstruction issue.

### A.2 The relative-frame relation

Differentiating $C$ with duality gives

$$
\dot C=DC+CE^T,\qquad
D^T+CEC^T=S,\qquad S:=C\dot C^T.
$$

Since $C$ is orthogonal, $S^T=-S$. The common tangent imposes

$$
Ce_1=e_1,\qquad Se_1=0.
$$

In particular, $C=1\oplus Q$ and $S=0\oplus S_\perp$. The second condition follows by differentiating the common-tangent identity; it must be retained when transferring the algebra back to actual moving frames.

The proposed curve data are $d(s),d^*(s),a(s),a^*(s),C(s)$, with $S(s)=C(s)\dot C(s)^T$. The algebraic problem treats their values at a fixed time. Coefficient reconstruction is distinct from ambient realization or uniqueness of curves up to ambient isometry. Even pointwise finiteness alone does not select a unique smooth branch along an interval.

### A.3 Fixed-data matrix formulation

Fix $n\ge2$, an **ordered** real orthogonal matrix $C$ with $Ce_1=e_1$, diagonal lists $d,d^*\in\mathbb R^n$, and superdiagonal lists $a,a^*$. Fix $S^T=-S$ for M′. The unknowns satisfy

$$
\begin{aligned}
D_{ii}=d_i,&\quad D_{i,i+1}=a_i,&\quad D_{ij}=0&& (j>i+1),\\
E_{ii}=d_i^*,&\quad E_{i,i+1}=a_i^*,&\quad E_{ij}=0&& (j>i+1).
\end{aligned}
$$

Only the strictly lower-triangular entries vary. The geometric/public problem uses $a_i,a_i^*>0$. Some historical block and singular controls allow zero superdiagonals; those controls are explicitly broader than the positive-curvature subclass. Reordering columns of a fixed $C$ changes the problem.

| Label | Affine relation | Normality | Additional assumptions | Status in this archive |
| --- | --- | --- | --- | --- |
| Claim M.1 | $D^T+CEC^T=0$ | Both matrices normal | (W) | General finiteness open |
| Claim M.2 | $D^T+CEC^T=0$ | Both matrices normal | No (W) | General finiteness open |
| Proposition M.3 | $D^T+CEC^T=S$ | Not required | All proper leading principal minors of $C$ nonzero | Compatible affine completion unique; criterion is also necessary for affine uniqueness |
| Proposition M′.4 | $D^T+CEC^T=S$ | Both matrices normal | $e_1$ cyclic for $S$ | At most one pair; outside the shared-tangent subclass |
| Claim M′.5 | $D^T+CEC^T=S$ | Both matrices normal | Neither (W) nor cyclicity assumed | General finiteness open |

Normality means $[D,D^T]=[E,E^T]=0$, where $[X,Y]=XY-YX$. In M, either normality condition implies the other. The original condition $DL=-LE^T$, $L=CE$, is equivalent to normality under the M dual relation.

The optional weighted relation is

$$
ACE=A^*C^TD,\qquad
A=\operatorname{diag}(1,a_1,\ldots,a_{n-1}),\quad
A^*=\operatorname{diag}(1,a_1^*,\ldots,a_{n-1}^*).\tag{W}
$$

Its geometric origin has not been established. It requires $d_1=d_1^*=0$ and $a_1=a_1^*$, so it is a substantive extra restriction. Without (W), necessary compatibility conditions still include $d_1^*=-d_1$ and $\sum_i(d_i+d_i^*)=0$.

### A.4 Reductions and established facts

**Eliminate the dual matrix:**

$$
E=C^T(S-D^T)C.
$$

All prescribed-entry conditions are linear in the unknown lower entries of $D$. For M′, imposing both normality conditions is equivalent to imposing normality of $D$ together with

$$
[S,D+D^T]=0.
$$

This commutation equation is linear and can be included before nonlinear solving. The earlier proposed condition

$$
-[S,D+D^T]=[D,D^T]+C[E,E^T]C^T
$$

is an identity under the affine relation, not an additional hypothesis.

**Principal-minor uniqueness.** For compatible fixed upper entries, the affine completion is unique iff $\det C[1:k,1:k]\ne0$ for every $1\le k<n$. Differences obey $X^T+CYC^T=0$ with strictly lower-triangular $X,Y$. Nonzero principal pivots give a triangular factorization forcing $X=Y=0$; a zero pivot produces a nonzero homogeneous direction. This criterion applies before normality or (W); after adding them, zero pivots need not imply nonuniqueness.

**Compactness.** Normality and the Hessenberg pattern imply the cut identities

$$
\sum_{i\le k<j}D_{ji}^2=a_k^2,\qquad 1\le k<n.
$$

Thus $|D_{ji}|\le\min_{i\le k<j}|a_k|$ for $i<j$, and the fixed-data real solution set is compact and algebraic. Compactness does not exclude positive-dimensional components. Full column rank of the reduced normality Jacobian at **every** real solution would imply finiteness; checking only recovered roots does not establish this premise. Singular Jacobians can occur at isolated roots.

**Cyclicity and its limitation.** Let

$$
V(S)=[e_1,Se_1,\ldots,S^{n-1}e_1],\quad
p(\lambda)=\det(\lambda I-S),\quad
q(\lambda)=\det(\lambda I-S[2:n,2:n]).
$$

For real skew-symmetric $S$, $\operatorname{Res}(p,q)=\det V(S)^2$. Proposition M′.4 applies when this is positive. Actual shared-tangent frames satisfy $Se_1=0$, so $e_1$ is never cyclic for $n>1$; here $p(\lambda)=\lambda q(\lambda)$. That proposition does not solve the geometric problem.

Fixing both diagonal lists is essential. If they are free and (W) is absent, $D\mapsto D+tI$, $E\mapsto E-tI$ already produces a continuous family. These examples do not contradict the fixed-diagonal claims.

## B. Numerical approaches already performed

### B.1 Common computational strategy

The searches eliminate $E$, solve the linear constraints with an SVD nullspace, and parameterize

$$
D(z)=D_0+\sum_{\ell=1}^r z_\ell X_\ell.
$$

They then solve quadratic normality equations using analytic Jacobians. For example, the derivative of $DD^T-D^TD$ in direction $X$ is

$$
XD^T+DX^T-X^TD-D^TX.
$$

Recovered roots are checked against the full original equations, prescribed entries, orthogonality, and (W) where requested. Root clustering, Jacobian singular values, and predictor/corrector continuation supply separate diagnostics. Witness construction is followed by freezing all scalar data before independent root searches. Failed witness generation and failed root searches are not treated as proofs of infeasibility.

### B.2 First M.1/M.2 survey — 7 November 2026

| Quantity | M.1, with (W) | M.2, without (W) |
| --- | ---: | ---: |
| Problem-condition combinations | 59 | 81 |
| Combinations with positive affine dimension $r$ | 37 | 47 |
| Largest matrix size | 20 | 20 |
| Largest affine dimension | 81 | 171 |
| Multistart attempts | 2,816 | 3,456 |

The 140 combinations reuse some fixed inputs under both conditions; they are not 140 independent random matrices. Families included regular affine controls, spectral/Arnoldi constructions, 21 irreducibly coupled inputs with nonzero superdiagonals, independent four-dimensional blocks, and singular cyclic controls.

Every recovered coupled root had a full-column-rank normality Jacobian: 66 root occurrences counted separately by input and condition, with minimum singular-value ratio about $4.62\times10^{-3}$. No continuous fixed-diagonal family was detected.

Exact block controls yielded $2^k$ solutions for $k=1,\ldots,5$, so finite multiplicity can grow with dimension. Some multistart runs missed known roots. Singular cyclic controls were isolated by the cut identities despite Jacobian nullity. A separate free-diagonal positive control successfully continued a known two-dimensional family in all four requested probes; it was excluded from the fixed-diagonal totals.

Acceptance used relative normality error below $10^{-10}$; clustering used $2\times10^{-5}$ times matrix scale. The master seed was `20261107`. The audit reported maximum normality error $7.34\times10^{-11}$ and dual-relation error $8.80\times10^{-16}$.

Entry points: [survey.py](numerics/m_survey/survey.py), [audit.py](numerics/m_survey/audit.py), [make_report.py](numerics/m_survey/make_report.py). Results: [summary.csv](numerics/m_survey/results/summary.csv), [results.json](numerics/m_survey/results/results.json), and `results/inputs/*.npz`. The local [README](numerics/m_survey/README.md) gives thresholds and full reproduction commands.

### B.3 Survey 2: alternative formulations and magnitude counterexamples — 7 November 2026

The second survey tested whether the earlier sign-branch pattern was intrinsic. It used trust-region and Levenberg–Marquardt searches, direct univariate polynomial candidates, complex total-degree homotopy for projected systems of affine dimension at most four, and simultaneous synthesis of two normal matrices sharing the same fixed data.

| Scope | Recorded result |
| --- | --- |
| Fixed inputs | 124, dimensions 4 through 12 |
| Local multistart calls | 10,752 |
| Implemented complex homotopy paths | 576 |
| Validated recorded root occurrences | 139 |
| Certified different-magnitude examples | Dimensions 6 and 8 |

**The sign-only strengthening is false.** In each certified example two solutions share $C,d,d^*,a,a^*$ but one corresponding squared entry differs by exactly $1/5$. Both have positive superdiagonals. The six-dimensional example does not satisfy (W); the eight-dimensional example does, so magnitude uniqueness fails under both M.2 and M.1.

For the eight-dimensional input, independent local algorithms recovered four roots in two magnitude classes. The certificates prove existence and local isolation of two roots in different classes, not an exhaustive count of four. Rational contraction bounds on radius-$10^{-35}$ boxes establish exact existence; separate nonsingular-minor bounds certify isolation in full affine dimensions 10 and 21. Decimal matrices are approximations to algebraic data defined by the polynomial systems and rational boxes.

Entry points: [survey2.py](numerics/m_survey2/survey2.py), [paired_probe.py](numerics/m_survey2/paired_probe.py), [certify_counterexamples.py](numerics/m_survey2/certify_counterexamples.py), and [run_all.py](numerics/m_survey2/run_all.py). See [final_results/summary.csv](numerics/m_survey2/final_results/summary.csv), [certificates/verification_summary.json](numerics/m_survey2/certificates/verification_summary.json), and the local [README](numerics/m_survey2/README.md). Failed targeted synthesis attempts are preserved and are not feasible counterexamples.

### B.4 M′ noncyclic continuation probe

This experiment retained fixed $C,S$, both diagonal lists, both superdiagonal lists, the affine relation, and both normality conditions, with and without (W). It used nonzero sources for which $e_1$ was noncyclic.

| Quantity | Recorded result |
| --- | ---: |
| Distinct fixed inputs | 67 |
| Problem-condition combinations | 102 |
| Dimensions | 4, 6, 8, 10, 12, 14 |
| Combinations with positive affine dimension | 28 |
| Combinations with singular tested Jacobian | 16 |
| Fixed-diagonal continuation probes | 152 |
| Successful fixed-diagonal moving probes | 0 |
| Successful probes for separate free-diagonal control | 36 of 36 |

Families comprised reflections, spectral sources with repeated frequencies or low rank, embedded independent blocks, and embedded cyclic controls. Continuation tested signed displacements of $0.01$ and $0.03$ times matrix scale, with full-equation acceptance below $2\times10^{-10}$. The singular controls were deliberately isolated examples, not evidence of surviving continuous freedom.

**Scope correction for geometric use:** inspection of all 103 saved input files, including the positive control, found none satisfying $Se_1\approx0$ at tolerance $10^{-10}$. Thus these are tests of the broader noncyclic algebraic M′ problem, not direct numerical tests of the nonzero-source shared-tangent subclass. Noncyclicity is weaker than $Se_1=0$.

Code: [probe.py](numerics/mprime_probe/probe.py). Recorded [summary.csv](numerics/mprime_probe/results/summary.csv), [results.json](numerics/mprime_probe/results/results.json), and the [README](numerics/mprime_probe/README.md) include inputs, roots, affine bases, and individual continuation outcomes.

### B.5 Reproduction and checks made for this archive

Each directory has its own dependency file or installation instructions. Recorded runs are included; the long searches need not be repeated to inspect them. From the indicated experiment directory:

```bash
# numerics/m_survey
python -m pip install -r requirements.txt
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python survey.py --starts 128 --out results_reproduced
python audit.py --results results_reproduced
```

```bash
# numerics/m_survey2: verify the existing exact certificates
python -m pip install -r requirements.txt
python certify_counterexamples.py --out certificates --verify-only

# Optional full rerun; the existing results are in final_results/
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python run_all.py --out results_reproduced
```

```bash
# numerics/mprime_probe; writes to its results/ directory
python -m pip install numpy scipy
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python probe.py
```

On 16 November 2026, the first survey's saved-matrix audit and Survey 2's audit and exact certificate verification were rerun successfully. The long multistart and continuation campaigns were not rerun. Historical recovery counts above refer to the included recorded campaigns.

## C. What alpha connections contribute

### C.1 Fix the sign convention and the cubic tensor

Let $g$ be the Fisher metric and $\Theta$ the Amari–Chentsov tensor. Define $\Theta^\sharp$ by $g(\Theta^\sharp(X,Y),Z)=\Theta(X,Y,Z)$, and use

$$
\nabla^{(\alpha)}=\nabla^{\mathrm{LC}}-\frac{\alpha}{2}\Theta^\sharp,
\qquad
\nabla^{(-\alpha)}=\nabla^{\mathrm{LC}}+\frac{\alpha}{2}\Theta^\sharp.
$$

This standard convention is stated, for example, by [Tarama and Françoise (2022)](https://doi.org/10.3390/psf2022005021), Section 2. Alpha connections are torsion free, conjugate in opposite-alpha pairs, and satisfy $(\nabla_X^{(\alpha)}g)(Y,Z)=\alpha\Theta(X,Y,Z)$. Conjugacy and torsion-freeness alone allow an arbitrary symmetric cubic tensor; fixing the statistical model fixes the tensor relevant here.

### C.2 Normality becomes a commutator condition

In the unstarred Frenet frame define

$$
\Omega_{ij}=g(\nabla_T^{\mathrm{LC}}N_i,N_j),\qquad
H_{ij}=\Theta(T,N_i,N_j).
$$

Then $\Omega^T=-\Omega$, $H^T=H$, and

$$
D=\Omega-\frac{\alpha}{2}H,\qquad
[D,D^T]=-\alpha[\Omega,H].
$$

For $\alpha\ne0$, normality is therefore exactly $[\Omega,H]=0$. Combining this with the dual affine relation gives the useful simultaneous-normality characterization

$$
[D,D^T]=[E,E^T]=0
\quad\Longleftrightarrow\quad
[\Omega,H]=0\ \text{and}\ [S,H]=0.
$$

These equivalences are derived here by separating the symmetric and skew-symmetric parts: with $B=-\alpha H/2$, the matrices in the common frame are $D=\Omega+B$ and $CEC^T=\Omega+S-B$.

This gives a geometric meaning and computational residual for the extra normality hypothesis. The Fisher cubic operator must commute with both the Levi–Civita frame rotation and the relative-frame source. In an eigenbasis of $H$, skew rotations can occur only inside repeated-eigenvalue blocks; if $H$ has simple spectrum, simultaneous normality forces $\Omega=S=0$ at that point. At $\alpha=0$ both coefficient matrices are skew-symmetric and automatically normal.

Alpha connections **do not automatically imply normality in their Frenet frames**. A Levi–Civita-parallel orthonormal frame makes the coefficient matrix symmetric, but generally is not the recursively constructed Frenet frame and need not have first vector $T$. Under a moving orthogonal change $\widetilde N=UN$,

$$
\widetilde D=\dot U U^T+UDU^T.
$$

The derivative term means arbitrary moving-frame changes need not preserve normality. Constant orthogonal changes do preserve it.

If the ambient cubic tensor is supplied, $D+D^T=-\alpha H$ directly supplies the missing lower entries. That is useful for numerical verification in a chosen model, but supplies more information than the invariant-data finiteness problem.

### C.3 New numerical information-geometry checks — 16 November 2026

The earlier material contained geometric derivations and examples rather than a numerical alpha-Frenet campaign. The present archive adds a small reproducible check on the interior of the three-outcome categorical simplex, where

$$
g_p(X,Y)=\sum_k\frac{X_kY_k}{p_k},\qquad
\Theta_p(X,Y,Z)=\sum_k\frac{X_kY_kZ_k}{p_k^2},\qquad \sum_kX_k=0.
$$

The probability curve is

$$
p(t)=\bigl(0.25+0.10t+0.03t^2,\ 0.35-0.04t+0.02t^2,\ 0.40-0.06t-0.05t^2\bigr).
$$

The code uses $T=p'(t)/\|p'(t)\|_g$ and differentiation along arc length. It constructs both Frenet frames independently for $\alpha=\pm1$, and the Levi–Civita frame for $\alpha=0$, at $t=-0.2,0,0.2$. Analytic curve derivatives and five-point frame differences with steps $10^{-4}$ and $5\times10^{-5}$ give 18 evaluated cases.

| Check | Recorded result |
| --- | --- |
| Largest identity/orthonormality residual | $1.34\times10^{-10}$ |
| Largest change in $D$ on halving the difference step | $1.81\times10^{-11}$ |
| Largest $\alpha=0$ normality residual | $1.34\times10^{-10}$ |
| Normality residual norms for $\alpha=\pm1$ | Approximately 4.98–9.57 |

The checks cover $D=\Omega-\alpha H/2$, the commutator identity, the overlap differential equation, the dual affine relation, and $Se_1=0$. The nonzero-alpha cases satisfy the dual geometric identities but fail normality, illustrating why the latter must be imposed separately. This is a two-dimensional tangent-space example: $S=0$ locally, so it does not test the additional nonzero-source condition $[S,H]=0$.

A positive example from the discussion note is also checked. On Euclidean $\mathbb R^2$, take $\nabla=\partial+K$, $\nabla^*=\partial-K$, with

$$
K(X,Y)=\sqrt2\,a\bigl(X_xY_x\partial_x-X_yY_y\partial_y\bigr).
$$

Along $\gamma(s)=s(1,1)/\sqrt2$, the recursive frames give

$$
D=E=\begin{pmatrix}0&a\\a&0\end{pmatrix},\quad
C=\operatorname{diag}(1,-1),\quad S=0.
$$

These connections are nonmetric, torsion free, dual, and flat, and both coefficient matrices are normal. The script uses $a=1.2$ and returns zero algebraic normality and dual-relation residuals. This is an abstract statistical-structure example; no particular likelihood family is asserted for it.

From the archive root, reproduce these focused checks with:

```bash
python -m pip install -r info_geometry/requirements.txt
python info_geometry/alpha_checks.py
```

### C.4 Next numerical question

A direct geometric search should generate or fit compatible data with $C=1\oplus Q$, $S=0\oplus S_\perp$, $S=C\dot C^T$, positive Frenet curvatures, and both commutators zero in the actual Frenet frames. In a fixed statistical model, geometric compatibility is stronger than the pointwise matrix equations. Tangent dimension at least three is needed for nonzero shared-tangent sources. This nonzero-source alpha-Frenet search has not been performed in the archived studies.

## Citation

Use the supplied citation for the public algebraic problem:

```bibtex
@MISC {516790,
    TITLE = {Dual Frenet moving frames: a finiteness problem for normal Hessenberg pairs},
    AUTHOR = {Eka38b (https://mathoverflow.net/users/602357/eka38b)},
    HOWPUBLISHED = {MathOverflow},
    NOTE = {URL:https://mathoverflow.net/q/516790 (version: 2026-11-10)},
    EPRINT = {https://mathoverflow.net/q/516790},
    URL = {https://mathoverflow.net/q/516790}
}
```

The included reports and code document exploratory work assisted by ChatGPT. Named open claims, proved reductions, exact certificates, historical numerical evidence, and the newly added consistency checks have different scopes and should be cited or described accordingly.
