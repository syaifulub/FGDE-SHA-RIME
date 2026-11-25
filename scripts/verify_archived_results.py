from __future__ import annotations
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
raw = pd.read_csv(ROOT / "data" / "frozen_validation" / "raw_runs.csv")
ranks = pd.read_csv(ROOT / "data" / "frozen_validation" / "ranks_lexicographic.csv")
stats = pd.read_csv(ROOT / "data" / "frozen_validation" / "problem_rank_statistics.csv")

print("Feasibility totals")
for method, grp in raw.groupby("Method"):
    n = len(grp)
    feasible = int(grp["Feasible"].sum())
    print(f"{method:24s} {feasible}/{n} ({100*feasible/n:.1f}%)")

print("\nMean feasibility-first ranks")
print(ranks.groupby("Method")["Lexico Rank"].mean().sort_values())

print("\nProblem-level statistics")
print(stats.to_string(index=False))
