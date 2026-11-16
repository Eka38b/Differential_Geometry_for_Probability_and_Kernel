# Numerical implementation guide

All implementations in this directory concern the fixed-data problem in MathOverflow question 516790 or its certified zero-source specialization.

| File | Function |
| --- | --- |
| [fiber.py](fiber.py) | Full affine reduction, source commutation, analytic Jacobian, two local algorithms, and validation |
| [replay.py](replay.py) | Search frozen inputs and write new results separately |
| [alternative_search.py](alternative_search.py) | Historical zero-source spectral construction, local algorithms, univariate candidates, and projected complex homotopy |
| [paired_synthesis.py](paired_synthesis.py) | Simultaneously construct two normal matrices with shared upper entries and a magnitude separation |
| [certify_counterexamples.py](certify_counterexamples.py) | Exact rational existence and local-isolation verification |
| [audit_data.py](audit_data.py) | Independent full-equation audit of saved solution datasets and an index generator |
| [check_markdown.py](check_markdown.py) | Check local links, brace structure, blocked math macro, and braced subscripts |

The [survey methodology documents](../Surveys/README.md) explain which code is used for historical records, modern fixed-input replay, source generation, and geometric checks. The adapted replay is a new numerical trajectory and does not promise to recover the same roots in the same order as the original campaign.
