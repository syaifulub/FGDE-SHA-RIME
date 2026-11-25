# Results summary

## Frozen validation

Ten previously unused CEC2006 problems were evaluated after the V2 configuration and projection settings were frozen: G12, G14, G15, G16, G17, G18, G19, G21, G22, and G23.

| Method | Feasible runs | Feasibility rate | Mean feasibility-first rank |
|---|---:|---:|---:|
| FGDE-SHA-RIME V2 | 171/250 | 68.4% | 1.95 |
| DE-Deb | 178/250 | 71.2% | 2.60 |
| RIME | 111/250 | 44.4% | 3.10 |
| FGDE-SHA-RIME pre-projection | 100/250 | 40.0% | 3.35 |
| CA-SHA-RIME-FR | 100/250 | 40.0% | 4.00 |

The Friedman test across the ten problem-level rank profiles gave p = 0.03817. After Holm correction, pairwise differences between V2 and DE-Deb (0.3320), RIME (0.2578), and the pre-projection configuration (0.09375) were not significant; the comparison with CA-SHA-RIME-FR remained significant (0.02344).

The results support improved feasibility relative to the pre-projection framework but do not establish general superiority over DE-Deb or RIME.
