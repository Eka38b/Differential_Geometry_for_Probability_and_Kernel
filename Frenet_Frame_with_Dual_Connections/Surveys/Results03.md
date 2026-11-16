# Survey 03 results

Performed 16 November 2026 using Python 3.12.14, NumPy 2.3.5, and SciPy 1.17.0.

| Quantity | Result |
| --- | ---: |
| Fixed inputs | 40 |
| Inputs with nonzero source and shared tangent | 40 |
| Inputs with a numerically zero leading principal minor | 20 |
| Inputs with positive affine dimension | 6 |
| Largest affine dimension | 5 |
| Local solver calls | 144 |
| Saved root occurrences, including witnesses | 40 |
| Additional roots found | 0 |
| Maximum full-equation residual | About $1.46\times10^{-15}$ |

| Nontrivial affine family | Replicates | Affine dimension | Saved roots per input |
| --- | ---: | ---: | ---: |
| Swap overlap, dimension 6 | 2 | 1 | 1 |
| Swap overlap, dimension 8 | 2 | 3 | 1 |
| Swap overlap, dimension 10 | 2 | 5 | 1 |

The six restricted normality Jacobians had full column rank. Their minimum smallest-to-largest singular-value ratio was about 0.554. All other inputs were already unique at the affine level, subject to compatibility. These outcomes support local isolation for this tested family, not universal finiteness or exhaustive numerical classification.
