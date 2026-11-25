# Data dictionary

## `data/frozen_validation/raw_runs.csv`

One row per method-problem-seed run.

- `Seed`: random seed.
- `Problem`: CEC2006 problem identifier.
- `Method`: evaluated method.
- `Objective`: final objective value.
- `Violation`: final aggregate constraint violation.
- `Feasible`: 1 if final violation <= 1e-12, otherwise 0.
- `Absolute Error`: absolute error from the known reference optimum for feasible runs; missing for infeasible runs.
- `First Feasible NFE`: first function-evaluation count at which feasibility was reached; 5001 indicates no feasible candidate within the 5,000-evaluation budget.

## `data/frozen_validation/convergence.csv`

Convergence records at predefined evaluation checkpoints.

## `data/frozen_validation/summary.csv`

Problem-method descriptive summaries, including feasibility rate, median error, median violation, and median first-feasible evaluation count.

## `data/frozen_validation/ranks_lexicographic.csv`

Problem-level feasibility-first ranks used for the frozen validation suite.

## `data/frozen_validation/problem_rank_statistics.csv`

Friedman and pairwise Wilcoxon-Holm results based on problem-level rankings.

## `data/diagnostic/v2_diagnostic_fd25.csv`

Diagnostic G5/G13 runs used during development of the multi-equality projection. These are development data, not unseen validation data.

## `data/development/dev_v1_v2_check.csv`

V1-V2 comparison on the original six development problems.
