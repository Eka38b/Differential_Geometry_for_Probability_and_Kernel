# Survey 01: the zero-source specialization

This survey studies the MathOverflow problem with $S=0$. All retained inputs fix both diagonal lists and both positive superdiagonal lists, and impose both normality equations. The dataset contains 69 compatible inputs from the recorded 7 November 2026 campaign.

Read [Methodology.md](Methodology.md) for construction, affine reduction, algorithms, thresholds, and interpretation. [Results.md](Results.md) records the retained results. The reusable inputs and matrices are in [Data/01_zero_source/](../../Data/01_zero_source/README.md).

From the archive root, replay a small subset into a new directory:

```bash
python Surveys/01_zero_source/run.py --out Reproduced/01 --starts 24 --limit 5
```

Omit the limit to search all retained inputs. The replay solver is adapted to the current public formulation; historical counts refer to the saved campaign, not to this new replay. It initializes from supplied compatible witnesses and does not claim exhaustive root recovery.
