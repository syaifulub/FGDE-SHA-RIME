from __future__ import annotations
import argparse
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import cec2006_unseen_defs as cec
from FGDE_SHA_RIME_V2_Frozen import SEEDS, run

METHODS = [
    "DE-Deb",
    "RIME",
    "CA-SHA-RIME-FR",
    "FGDE-SHA-RIME",
    "FGDE-SHA-RIME-V2",
]


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the frozen ten-problem validation suite.")
    parser.add_argument("--output", default="outputs/frozen_validation", help="Output directory")
    parser.add_argument("--smoke", action="store_true", help="Run only g12, FGDE-SHA-RIME-V2, seed 101")
    args = parser.parse_args()

    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)

    rows = []
    conv = []
    problems = ["g12"] if args.smoke else list(cec.REG.keys())
    methods = ["FGDE-SHA-RIME-V2"] if args.smoke else METHODS
    seeds = [101] if args.smoke else SEEDS

    for problem in problems:
        for method in methods:
            for seed in seeds:
                result, history = run(method, problem, seed)
                rows.append(result)
                conv.extend(history)
                print(problem, method, seed, result[5], result[4])

    with (out / "raw_runs.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Seed", "Problem", "Method", "Objective", "Violation", "Feasible", "Absolute Error", "First Feasible NFE"])
        w.writerows(rows)

    with (out / "convergence.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Seed", "Problem", "Method", "NFE", "Best Feasible Objective", "Minimum Violation"])
        w.writerows(conv)

    print(f"Wrote {len(rows)} runs to {out}")


if __name__ == "__main__":
    main()
