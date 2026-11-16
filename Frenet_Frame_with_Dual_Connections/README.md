# Dual Frenet moving frames: the MathOverflow finiteness problem

This archive concerns only [MathOverflow question 516790](https://mathoverflow.net/q/516790): whether fixed compatible data determine finitely many pairs of real normal Hessenberg matrices arising from dual Frenet frames. It records the algebraic formulation, numerical methods and data, and the interpretation of the normality assumption through information-geometric alpha connections.

Updated 16 November 2026. General finiteness is unresolved in this archive. Numerical recovery counts are lower bounds, except where a separate exact argument is supplied.

## The problem

Fix real matrices $C,S$ and real lists $d,d^{*},a,a^{*}$ satisfying

$$
C^{T}C=I,\qquad Ce_{1}=e_{1},\qquad S^{T}=-S,\qquad Se_{1}=0,
\qquad a_{i}>0,\quad a_{i}^{*}>0.
$$

Find real matrices $D,E$ with fixed upper entries

$$
D_{ii}=d_{i},\quad E_{ii}=d_{i}^{*},\quad
D_{i,i+1}=a_{i},\quad E_{i,i+1}=a_{i}^{*},\quad
D_{ij}=E_{ij}=0\quad(j>i+1),
$$

such that

$$
D^{T}+CEC^{T}=S,\qquad
DD^{T}=D^{T}D,\qquad EE^{T}=E^{T}E.
$$

Only strictly lower-triangular entries are unknown. The superscript star labels the second list; transpose is written $T$. The matrices are lower Hessenberg in this explicit convention. Assuming a pair exists, **must the solution set be finite?** Can a nonconstant continuous family occur with all prescribed superdiagonals positive?

See [Problem.md](docs/Problem.md) for the dual-frame derivation, compatibility conditions, affine reduction, principal-minor uniqueness criterion, cut bounds, and the distinction between finiteness and uniqueness. The zero-source case $S=0$ is included as a specialization of this same question.

## Directory guide

| Directory | Purpose |
| --- | --- |
| [docs/](docs/) | Algebraic formulation, information geometry, provenance, and Markdown conventions |
| [Code/](Code/) | Shared fixed-data solver, alternative zero-source methods, paired synthesis, certificate verification, and audits |
| [Surveys/](Surveys/) | Four experiments, each with a README, detailed methodology, results, and an entry point |
| [Data/](Data/) | Frozen inputs, saved matrices, CSV tables, JSON diagnostics, exact certificates, and numerical geometry checks |
| [references.bib](references.bib) | The supplied MathOverflow citation |

Begin with [Data/README.md](Data/README.md) and [Data/index.csv](Data/index.csv) to find data suitable for inspection or reuse. Relative links resolve after extracting the complete archive.

## Numerical studies

| Study | Scope | Archived result | Detailed document |
| --- | --- | --- | --- |
| [01: zero source](Surveys/01_zero_source/README.md) | 69 compatible fixed inputs; $S=0$; positive superdiagonals | 83 saved root occurrences; 2,176 multistart calls in the retained historical subset | [Methodology](Surveys/01_zero_source/Methodology.md) |
| [02: magnitude multiplicity](Surveys/02_magnitude/README.md) | 122 compatible fixed inputs; alternative algorithms and paired synthesis | 127 saved root occurrences; certified pairs in dimensions 6 and 8 with different entry magnitudes | [Methodology](Surveys/02_magnitude/Methodology.md) |
| [03: nonzero shared-tangent sources](Surveys/03_shared_tangent_source/README.md) | 40 newly generated inputs satisfying $S\ne0$ and $Se_{1}=0$ | Six inputs have positive affine dimension; 144 local solver calls found no additional root beyond the witnesses | [Methodology](Surveys/03_shared_tangent_source/Methodology.md) |
| [04: alpha connections](Surveys/04_alpha_connections/README.md) | Independent Frenet construction in a categorical statistical model | 18 consistency checks; illustrates why normality is an extra restriction | [Methodology](Surveys/04_alpha_connections/Methodology.md) |

The first two datasets are curated from the recorded 7 November 2026 campaigns. Their counts include only the fixed-data, positive-superdiagonal cases retained here. The third and fourth studies were performed on 16 November 2026. Geometric diagnostic matrices that fail normality are clearly separated from feasible solution datasets.

The certified examples disprove entrywise-magnitude uniqueness. They do not disprove finiteness. Full-rank Jacobians support isolation of the specified roots; they do not exclude other undiscovered roots or components.

## Reproduction

Install Python dependencies from the archive root:

```bash
python -m pip install -r requirements.txt
```

Audit saved solution datasets and verify the exact certificates:

```bash
python Code/audit_data.py
python Code/certify_counterexamples.py --out Data/02_magnitude/certificates --verify-only
```

Replay a small selection of fixed inputs without overwriting archived results:

```bash
python Surveys/01_zero_source/run.py --out Reproduced/01 --starts 24 --limit 5
python Surveys/02_magnitude/run.py --out Reproduced/02 --starts 24 --limit 5 --alternative
```

Generate the shared-tangent source survey and repeat the geometry checks:

```bash
python Surveys/03_shared_tangent_source/run.py --out Reproduced/03 --starts 24
python Surveys/04_alpha_connections/run.py --out Reproduced/04/results.json
```

The replay scripts search the saved inputs using the included solver implementations. They do not claim bit-for-bit reproduction of every historical generation attempt or search trajectory. Each methodology document states its construction rules, equations, algorithms, initialization, thresholds, validation, saved fields, and limitations.

## Information geometry

With the convention

$$
\nabla^{(\alpha)}=\nabla^{\mathrm{LC}}-\frac{\alpha}{2}\Theta^{\sharp},
$$

the Frenet coefficient matrix has the decomposition

$$
D=\Omega-\frac{\alpha}{2}H,\qquad
H_{ij}=\Theta(T,N_{i},N_{j}),\qquad \Omega^{T}=-\Omega.
$$

For $\alpha\ne0$, simultaneous normality is equivalent to

$$
[\Omega,H]=0,\qquad [S,H]=0.
$$

These conditions interpret the extra hypothesis in the MathOverflow question; they are not consequences of alpha-connection duality alone. See [Information_Geometry.md](docs/Information_Geometry.md) for the derivation, frame dependence, and numerical implications.

## Citation

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
