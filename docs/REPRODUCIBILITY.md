# Reproducibility notes

## Experimental chronology

The repository separates development evidence from frozen validation.

1. **Base-development suite:** G1, G3, G4, G6, G8, and G9.
2. **Diagnostic multi-equality development:** G5 and G13.
3. **Frozen validation:** G12, G14, G15, G16, G17, G18, G19, G21, G22, and G23.

G5 and G13 were used to motivate the multi-equality projection and therefore are not treated as unseen validation problems.

## Frozen settings

- Population size: 20
- Maximum function evaluations: 5,000
- Seeds: 101, 202, ..., 2525
- Equality tolerance: 1e-4
- Projection trigger: violation stagnation >= 4 iterations
- Projection frequency: every third iteration while triggered
- Finite-difference step: 1e-3 of the variable range
- Maximum sensitivity-selected dimensions: 6
- Projection correction bound: infinity norm <= 0.5
- Projection evaluations replace ordinary trials, so the total budget remains 5,000 evaluations.

## Environment

The manuscript specifies Python 3.12 with NumPy 2.3.5, pandas 2.2.3, SciPy 1.17.0, pymoo 0.6.1.6, and Matplotlib 3.10.8.

The archived frozen-validation package used CEC2006 definitions transcribed from the pymoo 0.6.1.6 source definitions. Before creating a permanent DOI release, rerun the frozen scripts in the pinned environment and confirm that the archived headline results are reproduced.

## Expected headline results

- Final FGDE-SHA-RIME V2: 171/250 feasible runs (68.4%)
- Pre-projection FGDE-SHA-RIME: 100/250 (40.0%)
- DE-Deb: 178/250 (71.2%)
- RIME: 111/250 (44.4%)
- CA-SHA-RIME-FR: 100/250 (40.0%)
- Mean feasibility-first ranks: 1.95, 3.35, 2.60, 3.10, and 4.00 respectively
- Friedman p = 0.03817
- Holm-adjusted p versus DE-Deb = 0.3320
- Holm-adjusted p versus RIME = 0.2578
- Holm-adjusted p versus pre-projection FGDE-SHA-RIME = 0.09375
- Holm-adjusted p versus CA-SHA-RIME-FR = 0.02344

These results do not support a claim of general superiority beyond the evaluated CEC2006 benchmark family.
