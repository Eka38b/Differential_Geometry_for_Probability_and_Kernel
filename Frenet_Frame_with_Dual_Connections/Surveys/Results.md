# Survey 04 results

Performed 16 November 2026 using Python 3.12.14 and NumPy 2.3.5.

| Quantity | Result |
| --- | --- |
| Evaluated time/alpha/step combinations | 18 |
| Largest identity or orthonormality error | About $1.34\times10^{-10}$ |
| Largest change in $D$ after halving the step | About $1.81\times10^{-11}$ |
| Largest alpha-zero normality residual | About $1.34\times10^{-10}$ |
| Normality norms at alpha plus/minus one | Approximately 4.98 to 9.57 |
| Positive nonmetric example | Zero algebraic normality and dual-relation residuals |

The categorical nonzero-alpha frames obey the dual geometric identities but are not normal. They are diagnostic failures of an automatic-normality assertion, not feasible solutions of the normal-pair problem. The positive example demonstrates that the added restriction can hold. The full matrices and separate residuals are in [Data/04_alpha_connections/results.json](../../Data/04_alpha_connections/results.json).
