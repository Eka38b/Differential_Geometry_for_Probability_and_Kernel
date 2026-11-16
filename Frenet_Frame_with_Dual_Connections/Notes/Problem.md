# The algebraic problem posed at MathOverflow

This note concerns the fixed-data finiteness question in [MathOverflow 516790](https://mathoverflow.net/q/516790). Both diagonal lists and both positive superdiagonal lists are fixed throughout.

## From a pair of Frenet frames

Let $(M,g,\nabla,\nabla^{*})$ be a positive-definite dualistic structure:

$$
Xg(Y,Z)=g(\nabla_{X}Y,Z)+g(Y,\nabla_{X}^{*}Z).
$$

Along a smooth unit-speed curve with tangent $T$, construct two full, individually orthonormal frames. Their common first vector is $N_{1}=N_{1}^{*}=T$. For each connection, subtract the projections onto the previously constructed frame vectors and normalize the nonzero residual. For the unstarred connection,

$$
W_{i+1}=\nabla_{T}N_{i}-\sum_{j=1}^{i}g(\nabla_{T}N_{i},N_{j})N_{j},\qquad
a_{i}=\|W_{i+1}\|_{g}>0,\qquad N_{i+1}=W_{i+1}/a_{i}.
$$

Work on an interval where every required residual is nonzero. Full frames are assumed: the cross-Gram matrix between different proper osculating subspaces need not be orthogonal.

Define row-index coefficient matrices and the overlap matrix by

$$
D_{ij}=g(\nabla_{T}N_{i},N_{j}),\quad
E_{ij}=g(\nabla_{T}^{\*}N_{i}^{\*},N_{j}^{\*}),\quad
C_{ij}=g(N_{i},N_{j}^{\*}).
$$

The projection recursion gives the prescribed lower-Hessenberg pattern. Orthonormality does not generally make $D$ skew-symmetric:

$$
D_{ij}+D_{ji}=-(\nabla_{T}g)(N_{i},N_{j}).
$$

The missing symmetric part is therefore part of the reconstruction problem. Differentiating the overlap matrix with duality gives

$$
\dot C=DC+CE^{T},\qquad
D^{T}+CEC^{T}=S,\qquad S=C\dot C^{T}.
$$

Orthogonality gives $S^{T}=-S$. The common tangent gives $Ce_{1}=e_{1}$, and differentiating this identity gives $Se_{1}=0$. Hence the shared-tangent restriction is stronger than merely assuming a noncyclic source.

## Fixed data and the question

At one time, prescribe $C,S,d,d^{\*},a,a^{\*}$. The matrices $C,S$ satisfy the conditions above, and every entry of $a,a^{*}$ is positive. The matrices $D,E$ have those diagonals and superdiagonals and vanish for $j>i+1$. Only strictly lower entries vary. Impose the dual affine relation and normality of both matrices.

Assuming compatibility, is the fiber of pairs $(D,E)$ finite? A positive answer would still permit multiple isolated pairs. Pointwise finiteness would also require further work before it implied smooth branch selection or a fundamental theorem for ambient curves.

Necessary compatibility conditions include

$$
d_{1}^{\*}=-d_{1},\qquad \sum_{i=1}^{n}(d_{i}+d_{i}^{\*})=0.
$$

They do not guarantee existence. The curve-data interpretation supplies $S$ from the derivative of $C$; the algebraic question freezes its value.

## Affine reduction and normality

Eliminate the second matrix:

$$
E=C^{T}(S-D^{T})C.
$$

Its prescribed upper entries give linear equations in the strictly lower entries of $D$. Furthermore,

$$
C[E,E^{T}]C^{T}=-[S,D+D^{T}]-[D,D^{T}],
\qquad [A,B]=AB-BA.
$$

Consequently both normality conditions are equivalent to

$$
[D,D^{T}]=0,\qquad [S,D+D^{T}]=0.
$$

The second equation is linear and belongs in the affine reduction. An SVD nullspace then gives

$$
D(z)=D_{0}+\sum_{\ell=1}^{r}z_{\ell}X_{\ell},
$$

leaving quadratic normality equations. For zero source, only the fixed-entry equations are needed in this affine step.

## A linear uniqueness criterion

Before imposing normality, a compatible affine completion is unique exactly when every proper leading principal minor of $C$ is nonzero. This condition is about the ordered matrix $C$.

Indeed, differences of two completions satisfy

$$
X^{T}+CYC^{T}=0
$$

with $X,Y$ strictly lower triangular. Orthogonality and nonzero principal pivots give a factorization $C=UL$, with $U$ invertible upper triangular and $L$ invertible lower triangular. Thus

$$
U^{-1}X^{T}U=-LYL^{-1}.
$$

The two sides are respectively strictly upper and strictly lower triangular, so both vanish. Conversely, if the leading block $C[1:k,1:k]$ is singular, choose nonzero vectors $r,\ell$ supported on the first $k$ coordinates with that block annihilating $r$ and its transpose annihilating $\ell$. Then $Cr$ and $C^{T}\ell$ are supported on the last coordinates, and

$$
X=-(Cr)\ell^{T},\qquad Y=(C^{T}\ell)r^{T}
$$

give a nonzero strictly lower homogeneous direction. Normality may still remove this affine freedom. Thus a zero pivot identifies a potentially difficult case, not a counterexample to finiteness.

## Cut identities and compactness

Sum the first $k$ diagonal normality equations. Terms internal to the first $k$ coordinates cancel; the only upper entry crossing the cut is $a_{k}$. Therefore

$$
\sum_{i \le k < j}D_{ji}^{2}=a_{k}^{2},\qquad 1 \le k < n.
$$

Hence

$$
|D_{ji}|\le\min_{i\le k < j}a_{k}\qquad(i < j).
$$

The fixed-data solution fiber is closed, bounded, and real algebraic. It is compact. This does not exclude a curve or a higher-dimensional component. Full column rank of the reduced normality Jacobian at every real solution would imply that every point is isolated, hence that this compact fiber is finite. Numerical rank checks at recovered roots do not verify the universal premise.

A singular Jacobian does not establish continuous freedom. For example, a cyclic permutation matrix with all superdiagonals equal to one and bottom-left entry one exhausts every cut budget; all other lower entries must vanish even if the first-order normality equations are singular.

## Exact multiplicity with positive superdiagonals

A four-dimensional zero-source example uses

```math
C=\left( \begin{array}{rrrr}
1&0&0&0\cr
0&0&0&-1\cr
0&0&1&0\cr
0&-1&0&0
\end{array} \right),\qquad
F(t)=\left(\begin{array}{rrrr}
0&a&0&0\cr
0&0&b&0\cr
0&t&0&b\cr
a&0&t&0
\end{array}\right).
```

For $D=E=F(t)$, the affine relation holds and normality reduces to $t^{2}=b^{2}-a^{2}$. At $a=3,b=5$ there are exactly two pairs, $t=\pm4$; at $a=b=3$ the unique root is singular. All superdiagonals remain positive. Survey 2 supplies certified higher-dimensional examples whose solutions also have different entry magnitudes, so sign branches are not a general description of multiplicity.

## Shared-tangent consequences

Since $Se_{1}=0$, the Krylov space generated by $e_{1}$ under $S$ is just the line through $e_{1}$. No theorem requiring that vector to be cyclic can resolve this question when $n>1$. The nonzero-source survey in this archive enforces $Se_{1}=0$ explicitly.

The present archive preserves only fixed diagonal lists. Changing them permits translations $D\mapsto D+tI$, $E\mapsto E-tI$, which concern a different question and are not used as finiteness counterexamples here.
