"""Standalone RIME-M and SHA-F benchmark on constrained CEC2006 problems.

The implementation is separated from the larger FGDE-SHA-RIME engine so that
the two internal comparators can be inspected and reproduced independently.
It preserves the reference protocol: NumPy default_rng, sequential in-place
updates, Deb feasibility rules, population 20, 5,000 evaluations, and equality
tolerance 1e-4.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys; sys.path.insert(0,"/home/claude/repro/code")
from cec2006_development_defs import get_problem


PROBLEMS = ["g1", "g3", "g4", "g6", "g8", "g9"]
METHODS = ["RIME-M", "SHA-F"]
SEEDS = [101 + 101 * i for i in range(25)]
POPULATION_SIZE = 20
MAX_EVALUATIONS = 5_000
EQUALITY_TOLERANCE = 1e-4
EXACT_FEASIBILITY_TOLERANCE = 1e-12


def evaluate(problem, x: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return objective, aggregate violation, and signed equality residuals."""
    result = problem.evaluate(
        np.atleast_2d(x),
        return_values_of=["F", "G", "H"],
        return_as_dictionary=True,
    )
    objective = np.asarray(result["F"], dtype=float).reshape(-1)
    violation = np.zeros(len(objective), dtype=float)

    inequalities = result.get("G")
    if inequalities is not None and np.size(inequalities):
        violation += np.maximum(np.asarray(inequalities, dtype=float), 0.0).sum(axis=1)

    equalities = result.get("H")
    if equalities is None or not np.size(equalities):
        signed_equalities = np.empty((len(objective), 0), dtype=float)
    else:
        signed_equalities = np.asarray(equalities, dtype=float)
        violation += np.maximum(
            np.abs(signed_equalities) - EQUALITY_TOLERANCE, 0.0
        ).sum(axis=1)

    return objective, violation, signed_equalities


def deb_better(
    trial_objective: float,
    trial_violation: float,
    parent_objective: float,
    parent_violation: float,
) -> bool:
    """Deb feasibility comparison for a minimization problem."""
    trial_feasible = trial_violation <= EXACT_FEASIBILITY_TOLERANCE
    parent_feasible = parent_violation <= EXACT_FEASIBILITY_TOLERANCE
    if trial_feasible and not parent_feasible:
        return True
    if parent_feasible and not trial_feasible:
        return False
    if trial_feasible and parent_feasible:
        return trial_objective < parent_objective
    return trial_violation < parent_violation


def positive_cauchy(rng: np.random.Generator, location: float, scale: float = 0.1) -> float:
    """Draw a positive Cauchy value and truncate it at one."""
    value = location + scale * rng.standard_cauchy()
    while value <= 0.0:
        value = location + scale * rng.standard_cauchy()
    return min(value, 1.0)


def lehmer_memory_update(memory: float, values: list[float], gains: list[float]) -> float:
    """Apply the weighted Lehmer-mean success-history update."""
    if not values:
        return memory
    successful = np.asarray(values, dtype=float)
    weights = np.asarray(gains, dtype=float)
    weights /= weights.sum()
    lehmer_mean = np.sum(weights * successful**2) / np.sum(weights * successful)
    return 0.9 * memory + 0.1 * float(lehmer_mean)


def _run(method: str, problem_name: str, seed: int) -> tuple[dict, list[dict]]:
    if method not in METHODS:
        raise ValueError(f"Unknown method: {method}")

    problem = get_problem(problem_name)
    rng = np.random.default_rng(seed)
    dimensions = problem.n_var
    lower, upper = np.asarray(problem.xl), np.asarray(problem.xu)
    span = upper - lower

    population = lower + rng.random((POPULATION_SIZE, dimensions)) * span
    objectives, violations, equalities = evaluate(problem, population)
    evaluations = POPULATION_SIZE
    iterations = (MAX_EVALUATIONS - POPULATION_SIZE) // POPULATION_SIZE

    objective_memory = 0.5
    violation_memory = 0.5
    best_feasible_so_far = np.inf
    convergence: list[dict] = []

    def record_checkpoint() -> None:
        feasible = violations <= EXACT_FEASIBILITY_TOLERANCE
        nonlocal best_feasible_so_far
        if np.any(feasible):
            best_feasible_so_far = min(
                best_feasible_so_far, float(np.min(objectives[feasible]))
            )
        convergence.append(
            {
                "Seed": seed,
                "Problem": problem_name,
                "Method": method,
                "NFE": evaluations,
                "Best_Feasible_Objective": (
                    best_feasible_so_far if np.isfinite(best_feasible_so_far) else np.nan
                ),
                "Minimum_Violation": float(np.min(violations)),
            }
        )

    record_checkpoint()

    for iteration in range(iterations):
        feasible_indices = np.flatnonzero(
            violations <= EXACT_FEASIBILITY_TOLERANCE
        )
        if len(feasible_indices):
            best_index = int(feasible_indices[np.argmin(objectives[feasible_indices])])
        else:
            best_index = int(np.argmin(violations))
        leader = population[best_index].copy()

        successful_objective_scales: list[float] = []
        objective_gains: list[float] = []
        successful_violation_scales: list[float] = []
        violation_gains: list[float] = []

        # Sequential, in-place updates are intentional and reproducibility-critical.
        for i in range(POPULATION_SIZE):
            donor_pool = [j for j in range(POPULATION_SIZE) if j != i]
            r1, r2, _ = rng.choice(donor_pool, 3, replace=False)

            progress = (iteration + 1) / iterations
            soft_factor = (1.0 - progress) * math.cos(math.pi * progress)
            trial = population[i].copy()
            soft_mask = rng.random(dimensions) < 0.5
            trial[soft_mask] = (
                leader[soft_mask]
                + soft_factor
                * (2.0 * rng.random(soft_mask.sum()) - 1.0)
                * span[soft_mask]
            )
            hard_mask = rng.random(dimensions) < (0.3 * progress)
            trial[hard_mask] = leader[hard_mask]

            channel = (
                "violation"
                if violations[i] > EXACT_FEASIBILITY_TOLERANCE
                else "objective"
            )
            if method == "RIME-M":
                scale = 0.5
            else:  # SHA-F
                memory = violation_memory if channel == "violation" else objective_memory
                scale = positive_cauchy(rng, memory)

            trial += scale * (population[r1] - population[r2])
            trial = np.clip(trial, lower, upper)
            trial_f, trial_v, trial_h = evaluate(problem, trial)
            evaluations += 1
            fy, vy, hy = float(trial_f[0]), float(trial_v[0]), trial_h[0]

            if deb_better(fy, vy, objectives[i], violations[i]):
                old_f, old_v = float(objectives[i]), float(violations[i])
                population[i] = trial
                objectives[i] = fy
                violations[i] = vy
                equalities[i] = hy

                if method == "SHA-F":
                    if channel == "violation" and vy < old_v:
                        successful_violation_scales.append(scale)
                        violation_gains.append(max(old_v - vy, 1e-15))
                    elif channel == "objective" and fy < old_f:
                        successful_objective_scales.append(scale)
                        objective_gains.append(max(old_f - fy, 1e-15))

        if method == "SHA-F":
            objective_memory = lehmer_memory_update(
                objective_memory, successful_objective_scales, objective_gains
            )
            violation_memory = lehmer_memory_update(
                violation_memory, successful_violation_scales, violation_gains
            )
        record_checkpoint()

    assert evaluations == MAX_EVALUATIONS
    feasible_indices = np.flatnonzero(violations <= EXACT_FEASIBILITY_TOLERANCE)
    if len(feasible_indices):
        result_index = int(feasible_indices[np.argmin(objectives[feasible_indices])])
        feasible = 1
    else:
        result_index = int(np.argmin(violations))
        feasible = 0

    reference_optimum = float(np.asarray(problem.pareto_front()).reshape(-1)[0])
    result = {
        "Seed": seed,
        "Problem": problem_name,
        "Method": method,
        "Objective": float(objectives[result_index]),
        "Violation": float(violations[result_index]),
        "Feasible": feasible,
        "Absolute_Error": (
            abs(float(objectives[result_index]) - reference_optimum)
            if feasible
            else np.nan
        ),
        "NFE": evaluations,
    }
    return result, convergence


def run_rime_m(problem_name: str, seed: int) -> tuple[dict, list[dict]]:
    """Run the RIME-M comparator."""
    return _run("RIME-M", problem_name, seed)


def run_sha_f(problem_name: str, seed: int) -> tuple[dict, list[dict]]:
    """Run the SHA-F comparator."""
    return _run("SHA-F", problem_name, seed)


def run_experiment(
    problems: list[str] | None = None,
    seeds: list[int] | None = None,
    output_directory: str | Path = "outputs/rime_m_sha_f",
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Run paired experiments and save raw, convergence, and summary outputs."""
    selected_problems = problems or PROBLEMS
    selected_seeds = seeds or SEEDS
    output = Path(output_directory)
    output.mkdir(parents=True, exist_ok=True)

    rows, histories = [], []
    runners = {"RIME-M": run_rime_m, "SHA-F": run_sha_f}
    for problem_name in selected_problems:
        for seed in selected_seeds:
            for method, runner in runners.items():
                result, trace = runner(problem_name, seed)
                rows.append(result)
                histories.extend(trace)
        print(f"Completed {problem_name}", flush=True)

    raw = pd.DataFrame(rows)
    convergence = pd.DataFrame(histories)
    summary = (
        raw.groupby(["Problem", "Method"], as_index=False)
        .agg(
            Runs=("Seed", "size"),
            Feasible_Rate=("Feasible", "mean"),
            Median_Violation=("Violation", "median"),
            Median_Absolute_Error=("Absolute_Error", "median"),
        )
    )

    raw.to_csv(output / "raw_runs.csv", index=False)
    convergence.to_csv(output / "convergence.csv", index=False)
    summary.to_csv(output / "summary.csv", index=False)
    metadata = {
        "problems": selected_problems,
        "methods": METHODS,
        "seeds": selected_seeds,
        "population_size": POPULATION_SIZE,
        "maximum_evaluations": MAX_EVALUATIONS,
        "equality_tolerance": EQUALITY_TOLERANCE,
        "exact_feasibility_tolerance": EXACT_FEASIBILITY_TOLERANCE,
        "sequential_updates": True,
    }
    (output / "metadata.json").write_text(json.dumps(metadata, indent=2))

    plot_convergence(convergence, output / "convergence.png")
    return raw, convergence, summary


def plot_convergence(convergence: pd.DataFrame, destination: str | Path) -> None:
    """Plot median minimum violation over paired runs for every problem."""
    problems = list(convergence["Problem"].unique())
    fig, axes = plt.subplots(2, 3, figsize=(14, 8), constrained_layout=True)
    for axis, problem_name in zip(axes.flat, problems):
        subset = convergence[convergence["Problem"] == problem_name]
        for method in METHODS:
            curve = (
                subset[subset["Method"] == method]
                .groupby("NFE")["Minimum_Violation"]
                .median()
            )
            axis.semilogy(curve.index, np.maximum(curve.values, 1e-16), label=method)
        axis.set_title(problem_name.upper())
        axis.set_xlabel("Function evaluations")
        axis.set_ylabel("Median minimum violation")
        axis.grid(alpha=0.25)
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2)
    fig.savefig(destination, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--quick",
        action="store_true",
        help="Run two problems and two seeds for a fast installation check.",
    )
    parser.add_argument("--output", default="outputs/rime_m_sha_f")
    args = parser.parse_args()
    problems = PROBLEMS[:2] if args.quick else PROBLEMS
    seeds = SEEDS[:2] if args.quick else SEEDS
    _, _, summary = run_experiment(problems, seeds, args.output)
    print("\nSUMMARY\n", summary.to_string(index=False))


if __name__ == "__main__":
    main()
