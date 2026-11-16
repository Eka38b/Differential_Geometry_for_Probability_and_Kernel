# Survey 02 methodology

## Objective and retained scope

Test entrywise-magnitude uniqueness inside the MathOverflow solution fiber, without adding constraints beyond the public question. All data are fixed during a root search, every superdiagonal is positive, and $S=0$. Different entry magnitudes would disprove this strengthening of uniqueness, but would remain compatible with finiteness.

The retained 7 November 2026 campaign contains 122 fixed inputs in dimensions 4, 5, 6, 7, 8, 10, and 12. Its families comprise one irreducible four-dimensional sign control, 42 generic spectral inputs, 42 balanced-second-moment inputs, 35 balanced-second-and-third-moment inputs, and two targeted paired inputs. Tables and counts here concern this retained subset only.

## Spectral and moment constructions

Construct a real normal matrix from one-dimensional real spectral blocks and two-dimensional blocks representing complex conjugate eigenvalues. Positive weights specify a starting vector. Arnoldi reduction produces the Hessenberg coefficient matrix. A spectral reflection, transformed into the Arnoldi basis and combined with alternating signs, gives the overlap matrix and the dual relation. Positive normalization fixes the superdiagonal convention before data are frozen.

Generic spectra typically leave no affine freedom. To create less trivial fibers, balance centered holomorphic spectral moments. The second-moment construction rescales imaginary parts so that the weighted centered square moment vanishes. The second-and-third-moment construction solves linear equalities in squared imaginary parts, with positive lower bounds, using linear programming. Failed positive constructions are not evidence of infeasible fibers; they are failed witness-generation attempts.

The saved metadata include spectral real parts, imaginary parts, weights, seeds, and the chosen moment conditions. The retained inputs typically have affine dimensions zero, one, or three; the targeted inputs have dimensions 10 and 21.

## Affine reduction and local algorithms

Eliminate $E=-C^{T}D^{T}C$ and impose all upper entries of the eliminated matrix. An SVD with relative threshold $10^{-10}$ gives a full basis for the admissible strictly lower perturbations. No sign pattern or fixed absolute values are imposed on these parameters.

The normality residual is quadratic. Pack its symmetric upper triangle with off-diagonal weights $\sqrt{2}$ and normalize by the squared Frobenius scale. Alternate trust-region reflective and Levenberg–Marquardt least squares, with analytic Jacobians. Starts are Gaussian in the affine coordinates with random scale multipliers $0.2$, $1$, or $2$ times the matrix scale divided by the square root of affine dimension.

The supplied normal witness is stored separately as an origin. The retained campaign has 10,496 local multistart calls. It uses up to 600 function evaluations per local call and termination tolerances near $3\times10^{-13}$. Full-equation acceptance is $2\times10^{-10}$; root clustering is $2\times10^{-6}$ times matrix scale.

## Univariate candidates

When affine dimension is one, evaluate the full residual at parameter values zero, one, and minus one to recover its constant, linear, and quadratic coefficient vectors. Select a nontrivial component polynomial, solve it, retain real candidates, and check every original residual equation. A root of the selected component alone is not sufficient.

## Projected complex homotopy

For affine dimensions two through four, choose random square projections of the full normality residual. Connect the projected polynomial system to the start system with component equations $z_{i}^{2}=1$, using a random complex phase. Track every start sign vector with a predictor from the differentiated homotopy equation and Newton correction at adaptive steps.

The implementation begins with step $0.035$, limits it to $0.08$, shrinks it after failed correction, and rejects excessive norm growth or insufficient progress. A projected real endpoint is accepted only after full original-equation validation and local polishing. The retained run tracked 576 paths, all reaching endpoints. This is not a rigorous complete homotopy solver: square projection, singularities, numerical tracking, and endpoint filtering prevent an exhaustive-count conclusion.

## Paired synthesis and independent fixed-data searches

Paired synthesis uses a different parameterization: construct two Hessenberg matrices simultaneously with common upper entries and the signed-reversal overlap. Decompose them as $A+H$ and $A-H$, with $H$ strictly lower triangular and no first-column freedom. Impose normality of both and a squared-entry difference of $1/5$.

The synthesis uses analytic Jacobians, bounded trust-region least squares, positive bounds on variable superdiagonals, and a finite attempt budget. Once a pair is found, freeze its shared scalar data and perform independent fixed-data local searches. The examples at dimensions 6 and 8 gave two and four recovered roots, respectively. Failed synthesis attempts do not supply counterexamples or infeasibility proofs.

`Code/paired_synthesis.py` retains this construction. The witness-generation diagnostic records, including unsuccessful attempts, are stored separately from accepted input files.

## Exact rational certification

Each certified example is defined by a square system of 13 rational quadratic equations, including the exact squared-entry separation. The unknowns include common scalar data and the two lower-entry parameter lists. A rational center and rational approximate inverse define

$$
T(z)=z-YF(z).
$$

On a closed infinity-norm box of radius $\rho=10^{-35}$, compute exact rational bounds

$$
\beta=\|YF(z_{0})\|_{\infty},\quad
\delta=\|I-YJ(z_{0})\|_{\infty},\quad
L=\delta+\rho\eta.
$$

Here $\eta$ bounds the Jacobian variation after multiplication by $Y$. Since the original equations are quadratic, their Jacobian is affine and constant Hessian bounds suffice. The checks $L<1$ and $\beta+L\rho<\rho$ prove that $T$ is a contraction mapping the box into itself. This proves existence of exact real algebraic data and two solutions.

Symbolic checks connect the reduced equations to every matrix normality equation, the dual relation, all prescribed entries, and the shared tangent. The boxes also preserve positive superdiagonals. Separate rational inverse bounds certify full-rank minors of the fixed-data normality Jacobians, in affine dimensions 10 and 21, proving local isolation of the specified roots. They do not count all roots in those fibers.

## Code, outputs, and limits

`Code/alternative_search.py` retains the spectral constructions, two local methods, univariate candidates, and projected homotopy. The current replay reads the frozen inputs; it does not rerun the complete historical generation campaign. `Code/certify_counterexamples.py` verifies the exact stored certificates. `Data/02_magnitude/` contains summary CSV, detailed JSON, NPZ matrices, certificate boxes, and synthesis diagnostics.

The result is an exact counterexample to magnitude uniqueness within the public problem. It establishes neither a continuous family nor the maximal number of isolated solutions. Seeded numerical searches can recover different sets on different platforms.
