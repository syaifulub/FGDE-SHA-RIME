"""Extended comparison: 8 methods x 10 frozen-validation problems x 25 seeds at a given budget.

Usage:
  python extended_comparison/run_extended_comparison.py --nfe 50000 --output outputs/extended_50000
  python extended_comparison/run_extended_comparison.py --nfe 5000 --baselines-only --output outputs/extended_5000
Each (problem, method) pair is written to its own CSV, so interrupted runs can be resumed.
"""
from __future__ import annotations
import argparse, csv, os, sys, time, pathlib
HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent)); sys.path.insert(0, str(HERE))
import FGDE_SHA_RIME_V2_Frozen as M
import baselines as B

ORIGINAL = ["FGDE-SHA-RIME-V2", "FGDE-SHA-RIME", "CA-SHA-RIME-FR", "RIME", "DE-Deb"]
EXTERNAL = ["eps-DEag", "ECHT-DE", "L-SHADE-Deb"]
PROBLEMS = ["g12", "g14", "g15", "g16", "g17", "g18", "g19", "g21", "g22", "g23"]
HEAD = ["Seed", "Problem", "Method", "Objective", "Violation", "Feasible", "Absolute Error", "First Feasible NFE"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nfe", type=int, default=50000)
    ap.add_argument("--output", default="outputs/extended")
    ap.add_argument("--baselines-only", action="store_true")
    a = ap.parse_args()
    out = pathlib.Path(a.output); out.mkdir(parents=True, exist_ok=True)
    M.NFE = a.nfe; M.ITERS = (a.nfe - M.POP) // M.POP
    cps = sorted({40, 220, 520, 1020, 1520, 2020, min(5000, a.nfe), a.nfe})
    methods = ([] if a.baselines_only else ORIGINAL) + EXTERNAL
    for p in PROBLEMS:
        for m in methods:
            fn = out / f"{a.nfe}__{p}__{m}.csv"
            if fn.exists():
                continue
            t = time.time(); rows = []; conv = []
            for s in M.SEEDS:
                if m in EXTERNAL:
                    r, h, used = B.run(m, p, s, a.nfe, cps); assert used == a.nfe
                else:
                    r, h = M.run(m, p, s)
                rows.append(r); conv.extend(h)
            with open(fn, "w", newline="") as f:
                w = csv.writer(f); w.writerow(HEAD); w.writerows(rows)
            with open(out / f"conv_{a.nfe}__{p}__{m}.csv", "w", newline="") as f:
                w = csv.writer(f); w.writerow(["Seed", "Problem", "Method", "NFE", "Best Feasible Objective", "Minimum Violation"]); w.writerows(conv)
            print(a.nfe, p, m, round(time.time() - t, 1), flush=True)


if __name__ == "__main__":
    main()
