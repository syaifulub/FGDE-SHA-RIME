# FGDE-SHA-RIME

**FGDE-SHA-RIME: A Feasibility-Guided Hybrid Method for Continuous Constrained Optimization**

This repository contains the reference implementation, benchmark definitions, diagnostic-development records, frozen-validation results, convergence records, and statistical outputs associated with the FGDE-SHA-RIME MethodsX manuscript.

## Method overview

FGDE-SHA-RIME combines Differential Evolution and success-history adaptive RIME with:

- feasibility-guided comparison;
- feasible-solution archiving and periodic reinjection;
- dynamic epsilon relaxation;
- stagnation-triggered feasibility restoration;
- adaptive DE-RIME operator selection; and
- a budget-neutral finite-difference projection for multiple equality constraints.

The projection evaluations replace an equal number of ordinary population trials, so the reported budget remains exactly **5,000 function evaluations per run**.

## Experimental design

The study uses a three-stage design:

1. **Base-framework development and ablation:** G1, G3, G4, G6, G8, and G9.
2. **Diagnostic multi-equality development:** G5 and G13.
3. **Frozen validation after all V2 settings were fixed:** G12, G14, G15, G16, G17, G18, G19, G21, G22, and G23.

The frozen-validation suite was not used to modify the final V2 configuration.

## Headline frozen-validation results

| Method | Feasible runs | Feasibility rate | Mean feasibility-first rank |
|---|---:|---:|---:|
| **FGDE-SHA-RIME V2** | **171/250** | **68.4%** | **1.95** |
| DE-Deb | 178/250 | 71.2% | 2.60 |
| RIME | 111/250 | 44.4% | 3.10 |
| FGDE-SHA-RIME pre-projection | 100/250 | 40.0% | 3.35 |
| CA-SHA-RIME-FR | 100/250 | 40.0% | 4.00 |

The problem-level Friedman test gave **p = 0.03817**. After Holm correction, comparisons of V2 with DE-Deb and RIME were not significant. The repository therefore does **not** claim general superiority beyond the evaluated CEC2006 benchmark family.

## Repository structure

```text
FGDE-SHA-RIME/
├── FGDE_SHA_RIME_V2_Frozen.py
├── cec2006_unseen_defs.py
├── cec2006_development_defs.py
├── requirements.txt
├── environment.yml
├── CITATION.cff
├── data/
│   ├── development/
│   │   └── dev_v1_v2_check.csv
│   ├── diagnostic/
│   │   └── v2_diagnostic_fd25.csv
│   └── frozen_validation/
│       ├── raw_runs.csv
│       ├── convergence.csv
│       ├── summary.csv
│       ├── ranks_lexicographic.csv
│       ├── problem_rank_statistics.csv
│       └── ...
├── scripts/
│   ├── run_frozen_validation.py
│   └── verify_archived_results.py
└── docs/
    ├── REPRODUCIBILITY.md
    ├── DATA_DICTIONARY.md
    ├── RESULTS_SUMMARY.md
    └── MANUSCRIPT_V2_FROZEN_UPDATE.md
```

## Installation

Python 3.12 is recommended.

```bash
python -m venv .venv
```

Activate the environment, then install the pinned dependencies:

```bash
pip install -r requirements.txt
```

Alternatively:

```bash
conda env create -f environment.yml
conda activate fgde-sha-rime
```

## Quick verification of archived results

```bash
python scripts/verify_archived_results.py
```

This recomputes the overall feasibility totals and mean feasibility-first ranks from the archived CSV files.

## Smoke test

To run one frozen-validation case:

```bash
python scripts/run_frozen_validation.py --smoke
```

## Full frozen-validation run

```bash
python scripts/run_frozen_validation.py --output outputs/frozen_validation
```

The full run evaluates 10 problems × 5 methods × 25 seeds = **1,250 runs**.

## Reproducibility notes

See [`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md) before creating a permanent release. In particular, the archived frozen-validation package used CEC2006 definitions transcribed from the pymoo 0.6.1.6 source definitions. A final pinned-environment rerun is recommended before DOI archiving.

## Data interpretation

- Objective errors are defined only for feasible runs.
- Final numerical feasibility uses aggregate violation <= 1e-12.
- Equality residuals use tolerance 1e-4 when forming aggregate constraint violation.
- G5 and G13 are diagnostic development problems, not unseen validation cases.
- Problem-level rankings use a feasibility-first lexicographic rule for the frozen validation suite.

## Authors

- Syaiful Anam — Universitas Brawijaya
- Ahmad Zarkasi — Universitas Mulawarman
- Hilmi Aziz Bukhori — Universitas Brawijaya
- Andrea Tri Rian Dani — Universitas Mulawarman

## Citation

The manuscript citation should be added after publication. A preliminary `CITATION.cff` is included and should be updated with the final article DOI and GitHub repository URL.

## License

No software license has been assigned in this package. Choose and add a license before making the repository public. For reusable academic software, MIT or BSD-3-Clause are common options.
