# Survey 04 methodology

## Purpose and relation to the public problem

Verify the alpha-connection decomposition of the Frenet matrix and check whether normality follows automatically in a concrete statistical model. The relevant matrices are calculated from independently constructed Frenet frames, not chosen by making an arbitrary orthonormal frame parallel. This is a diagnostic calculation; normality failures are part of its output.

## Model and curve

Use the interior of the three-outcome categorical simplex. A point is $p=(p_{1},p_{2},p_{3})$, with positive entries summing to one. Tangent vectors have component sum zero. The Fisher metric and cubic tensor are

$$
g_{p}(X,Y)=\sum_{k=1}^{3}\frac{X_{k}Y_{k}}{p_{k}},\qquad
\Theta_{p}(X,Y,Z)=\sum_{k=1}^{3}\frac{X_{k}Y_{k}Z_{k}}{p_{k}^{2}}.
$$

Take

$$
p(t)=\bigl(0.25+0.10t+0.03t^{2},\ 0.35-0.04t+0.02t^{2},\ 0.40-0.06t-0.05t^{2}\bigr).
$$

Evaluate at $t=-0.2,0,0.2$. The probabilities remain strictly positive on the evaluation neighborhoods. Compute the first and second curve derivatives analytically. Set speed $v=\|p'(t)\|_{g}$ and $T=p'(t)/v$; differentiation along arc length is $v^{-1}d/dt$.

## Alpha covariant differentiation

Use the convention $\nabla^{(\alpha)}=\nabla^{\mathrm{LC}}-\alpha\Theta^{\sharp}/2$. In probability components, the connection correction to the ordinary directional derivative is

$$
\frac{1+\alpha}{2}\left(g_{p}(X,Y)p-\left(\frac{X_{k}Y_{k}}{p_{k}}\right)_{k=1}^{3}\right).
$$

At alpha minus one this is the mixture connection; at alpha plus one it is the exponential connection. Use this explicit correction independently of the cubic decomposition later checked.

## Frenet construction and numerical derivatives

For each of $\alpha=-1,0,1$, calculate $\nabla_{T}^{(\alpha)}T$ using the analytic derivative of the normalized tangent. Subtract its component along $T$. Normalize the residual with its positive Fisher norm to form $N_{2}$. Reject a residual below $10^{-9}$ as a degenerate frame.

The resulting two-dimensional frame has first vector $T$ and positive curvature. For the dual pair, construct the opposite-alpha frame independently. Differentiating these frame functions requires the derivative of the second vector; approximate the complete frame derivative with the five-point formula

$$
f'(t)\approx\frac{f(t-2h)-8f(t-h)+8f(t+h)-f(t+2h)}{12h}.
$$

Use $h=10^{-4}$ and $h=5\times10^{-5}$ and divide by speed to obtain arc-length derivatives. There are three times, three alpha choices, and two step sizes, for 18 evaluated cases. This is a finite-difference consistency check, not exact arithmetic or a rigorous error bound.

## Independently compared quantities

Compute $D$ and $E$ from covariant frame derivatives. Compute $\Omega$ by applying the Levi–Civita connection to the same unstarred frame, and compute $H$ directly by contraction of the categorical cubic tensor. Compute $C$ from the cross-Gram matrix and $S=C\dot C^{T}$ from the overlap derivative.

Check orthonormality, skewness of $\Omega$, and

$$
D=\Omega-\frac{\alpha}{2}H,\qquad
[D,D^{T}]=-\alpha[\Omega,H],\qquad
\dot C=DC+CE^{T},\qquad D^{T}+CEC^{T}=S.
$$

Also check skewness of $S$ and $Se_{1}=0$. The script rejects an identity error exceeding $10^{-7}$ and rejects a change in $D$ above $10^{-7}$ on halving the difference step. The reported raw residual norms are much smaller than these conservative thresholds.

Record the normality norm separately; it is not required to vanish in this diagnostic experiment. Failure of normality is precisely what tests whether the extra hypothesis is automatic.

## Positive nonmetric example

The script also evaluates the explicit dual flat Euclidean example in the information-geometry note with $a=1.2$. It has $D=E$ symmetric with off-diagonal entries $a$, overlap diagonal entries 1 and minus 1, and zero source. The direct algebraic checks give zero normality and dual-relation residuals.

This example shows that nonmetric dual structures can satisfy the hypotheses. It is an abstract statistical-structure example, not an asserted likelihood realization.

## Saved data and limitations

The JSON stores all 18 sets of $D,E,C,S,\Omega,H$, curvature, alpha, time, difference step, residuals, and step-change diagnostics. The explicit positive example and environment versions are also stored. No random seeds or statistical fitting are involved.

The categorical tangent space has dimension two. On a regular interval the two common-tangent orthonormal frames differ only by a discrete choice of normal, so the relative source is zero. The calculation verifies neither the nonzero-source commutator condition in a nontrivial source example nor global solution finiteness. It demonstrates the need to check normality in the actual alpha-Frenet frames.
