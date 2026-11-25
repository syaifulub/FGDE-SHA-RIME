# RIME-M and SHA-F (development comparators)

`rime_m_sha_f.py` is the authors' original RIME-M / SHA-F implementation (Colab notebook, cell 3),
with the pymoo problem loader replaced by `cec2006_development_defs.get_problem`.

Run: `python rime_m_sha_f.py` (6 development problems x 25 seeds, 5,000 FEs).
Reproduces Table 2 (feasible runs: RIME-M 140/150, SHA-F 137/150; identical penalized mean ranks)
and the Table 5 Wilcoxon results for RIME-M and SHA-F.
Output: `data/development_comparators/raw_runs.csv`.
