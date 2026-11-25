# Projection ablation and minimal variant (Table 18)

Post hoc ablation of the multi-equality projection, run after the configuration was frozen.

Variants (`run_abl.py`): P (FGDE-SHA-RIME-P), B (FGDE-SHA-RIME-B), random selection of probed
variables, single-equality projection, no equality interpolation, Minimal-P (DE-only search,
fixed F = 0.5 and mutation probability 0.5, with epsilon relaxation, archive, restoration and
projection), DE-Deb, and DE-Deb+P (DE-Deb with the same projection).

Problems: G5, G13 (`cec2006_diagnostic_defs.py`, implemented from the CEC2006 report), G3, G11,
G14, G15, G17, G21, G22, G23; P, B, Minimal-P and DE-Deb also on G1, G4, G6, G8, G9, G12, G16,
G18, G19. 25 seeds, 5,000 evaluations.

Run: `python run_abl.py <worker> <n_workers> out` then `python an_abl.py`.
Records: `data/projection_ablation/raw_runs.csv`; summary: `data/projection_ablation/abl_results.json`.
