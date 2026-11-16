# Alpha connections and the normality condition

This note explains the extra normality hypothesis in the MathOverflow problem. It uses the same row-index convention for Frenet coefficients and the same shared tangent. It does not assume that alpha connections automatically produce normal coefficient matrices.

## Statistical structure and signs

Let $g$ be the Fisher metric and $\Theta$ the symmetric Amari–Chentsov cubic tensor. Define $\Theta^{\sharp}$ by

$$
g(\Theta^{\sharp}(X,Y),Z)=\Theta(X,Y,Z).
$$

Use the convention

$$
\nabla^{(\alpha)}=\nabla^{\mathrm{LC}}-\frac{\alpha}{2}\Theta^{\sharp},\qquad
\nabla^{(-\alpha)}=\nabla^{\mathrm{LC}}+\frac{\alpha}{2}\Theta^{\sharp}.
$$

Then the two connections are torsion free, dual, and satisfy

$$
(\nabla_{X}^{(\alpha)}g)(Y,Z)=\alpha\Theta(X,Y,Z).
$$

Conjugacy and torsion-freeness alone leave the symmetric cubic tensor undetermined. A chosen statistical model fixes the tensor, making the nonmetric part concrete.

## Symmetric and skew parts in the Frenet frame

In the unstarred Frenet frame set

$$
\Omega_{ij}=g(\nabla_{T}^{\mathrm{LC}}N_{i},N_{j}),\qquad
H_{ij}=\Theta(T,N_{i},N_{j}),\qquad B=-\frac{\alpha}{2}H.
$$

The Levi–Civita part is skew-symmetric and the cubic part symmetric:

$$
\Omega^{T}=-\Omega,\qquad H^{T}=H,\qquad D=\Omega+B.
$$

Direct multiplication gives

$$
[D,D^{T}]=2[\Omega,B]=-\alpha[\Omega,H].
$$

For $\alpha\ne0$, unstarred normality is therefore equivalent to $[\Omega,H]=0$. The dual affine relation yields, in the same frame,

$$
CEC^{T}=S-D^{T}=\Omega+S-B.
$$

Normality is preserved by this orthogonal conjugation. Applying the same symmetric/skew calculation to the right-hand side proves

$$
[D,D^{T}]=[E,E^{T}]=0
\quad\Longleftrightarrow\quad
[\Omega,H]=0\ \text{and}\ [S,H]=0
\qquad(\alpha\ne0).
$$

This is a useful interpretation and numerical residual for the MathOverflow hypothesis. It constrains how the Fisher cubic operator can coexist with frame rotation and the relative-frame source. In an eigenbasis of $H$, a commuting skew-symmetric matrix can rotate only within repeated-eigenvalue blocks. If $H$ has simple spectrum, simultaneous normality forces both $\Omega$ and $S$ to vanish at that point.

At $\alpha=0$, the coefficient matrices in either orthonormal frame are skew-symmetric and automatically normal. The nonzero-alpha equivalence should not be divided by alpha at zero.

## Frame dependence

If a moving orthonormal frame is changed by $\widetilde N=UN$, then

$$
\widetilde D=\dot U U^{T}+UDU^{T}.
$$

The derivative term means that arbitrary moving-frame changes need not preserve normality. Constant orthogonal changes do preserve it. A Levi–Civita-parallel frame makes the alpha coefficient matrix symmetric, but generally is not the recursively constructed Frenet frame and need not have first vector $T$. It cannot replace that frame when checking the public question.

Supplying the ambient cubic tensor also supplies the missing symmetric part:

$$
D+D^{T}=-\alpha H,\qquad D_{ij}=-D_{ji}-\alpha H_{ij}\quad(i>j).
$$

This is convenient for verification in a specified statistical model. It gives more input than the fixed invariant-data finiteness problem.

## What has been checked numerically

[Survey 04](../Surveys/04_alpha_connections/README.md) constructs the two Frenet frames independently on the interior of the three-outcome categorical simplex. It verifies the alpha decomposition, the commutator identity, and the dual overlap equations. The nonzero-alpha examples fail normality while satisfying the other geometric identities. They are diagnostic examples, not feasible normal-pair solutions.

The same script checks a positive nonmetric example on Euclidean two-space. With $\nabla=\partial+K$, $\nabla^{*}=\partial-K$, and

$$
K(X,Y)=\sqrt{2}\,a\bigl(X_{x}Y_{x}\partial_{x}-X_{y}Y_{y}\partial_{y}\bigr),
$$

the unit-speed line with tangent $(1,1)/\sqrt{2}$ gives

$$
D=E=\begin{pmatrix}0&a\\a&0\end{pmatrix},\qquad
C=\begin{pmatrix}1&0\\0&-1\end{pmatrix},\qquad S=0.
$$

The cubic form associated with $K$ is symmetric, and the constant coordinate connection matrices commute, so the connections are dual, torsion free, and flat. Both Frenet matrices are symmetric and normal with positive superdiagonals. This is an abstract statistical structure; no particular likelihood family is asserted.

The categorical example has tangent dimension two, so its source is zero on a regular frame interval. A nonzero shared-tangent source requires tangent dimension at least three. Survey 03 tests such sources algebraically, but does not assert their realization in a fixed Fisher statistical model. A full nonzero-source numerical alpha-Frenet search remains a separate task.
