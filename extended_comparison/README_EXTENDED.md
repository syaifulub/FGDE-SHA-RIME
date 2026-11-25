# Extended comparison (revision, Reviewer comment on baselines and budget)

Copy `extended_comparison/` and `data/extended_comparison/` into the repository root.

## Contents
- `baselines.py`: re-implementations of eps-DEag (Takahama & Sakai, 2010), ECHT-DE
  (Mallipeddi & Suganthan, 2010), and L-SHADE (Tanabe & Fukunaga, 2014) with Deb's rules.
  All methods share Eq. (2), delta_h = 1e-4, the 1e-12 feasibility threshold, the 25 seeds,
  and exact evaluation accounting (finite-difference probes included).
- `development_ablations.py`: the frozen implementation with the two development-only
  ablations SHA-RIME and CA-SHA-RIME enabled. Results for all other methods are identical to
  `FGDE_SHA_RIME_V2_Frozen.py`. Its SHA-RIME and CA-SHA-RIME runs reproduce the feasibility
  totals reported for the original development benchmark (146/150 and 137/150).
- `run_extended_comparison.py`: runs the comparison at a chosen budget (resumable).
- `analyze_extended_comparison.py`: recomputes `data/extended_comparison/statistical_summary.json`.

## Data
- `raw_runs_5000_external_baselines.csv`: 750 runs (3 baselines x 10 problems x 25 seeds).
  The 5,000-evaluation results of the five original methods are `data/frozen_validation/raw_runs.csv`.
- `raw_runs_50000_all_methods.csv`: 2,000 runs (8 methods x 10 problems x 25 seeds).
- `convergence_*.csv`: best feasible objective and minimum violation at checkpoints.
- `development_rerun_six_methods.csv`: rerun of DE-Deb, RIME, CA-SHA-RIME-FR,
  FGDE-SHA-RIME, SHA-RIME, and CA-SHA-RIME on G1, G3, G4, G6, G8, and G9.
- `statistical_summary.json`: feasibility, eight-method feasibility-first ranks, Friedman,
  Holm-adjusted Wilcoxon, and run-level McNemar results for both budgets.

## Environment
Executed with Python 3.13.18, NumPy 2.5.3, SciPy 1.18.1, pandas 3.0.5.
The final configuration's least-squares projection is sensitive to floating-point
differences between numerical-library builds; individual runs can differ across
environments without changing the statistical conclusions.

## Commands
    python extended_comparison/run_extended_comparison.py --nfe 50000 --output outputs/extended_50000
    python extended_comparison/run_extended_comparison.py --nfe 5000 --baselines-only --output outputs/extended_5000
    python extended_comparison/analyze_extended_comparison.py

## Remaining CEC2006 problems (G2, G7, G10, G11, G24)
- `cec2006_remaining_defs.py` (repository root): definitions transcribed from the CEC2006
  report and checked at the published best known solutions; reported optima are used as references.
- `extended_comparison/run_remaining_problems.py`: runs all eight methods at 5,000 and 50,000 evaluations.
- `data/remaining_problems/`: run-level records, convergence records, and `summary.json` (Table 12).
G20 is excluded because no feasible solution is known.

## Sensitivity analysis (Table 13)
- `extended_comparison/sensitivity/fgde_parameterized.py`: frozen implementation with 18 hand-set
  constants exposed in `PRM`; with `DEFAULT` it reproduces `FGDE_SHA_RIME_V2_Frozen.py` exactly.
- `run_sensitivity.py`: one-at-a-time variations (36 configurations + default) on the frozen suite.
- `analyze_sensitivity.py` and `data/sensitivity/`: run records (9,250 runs) and paired McNemar/Holm summary.

## Statistical analyses (Tables 14-16)
- `extended_comparison/statistics/run_statistical_analyses.py`: per-problem run-level McNemar and
  Wilcoxon tests with risk differences and rank-biserial correlations, Friedman post hoc tests with
  a control method (Holm), and the projection analysis restricted to equality-constrained problems.
- `data/statistics/`: outputs at 5,000 and 50,000 evaluations.

## Credit-assignment correction, violation normalization, runtime (Table 17)
- `extended_comparison/corrections/fgde_corrections.py`: options `fix_credit` and `norm_viol`
  (defaults reproduce the frozen implementation exactly).
- `run_corrections.py`, `measure_runtime.py`; outputs in `data/corrections/`.
