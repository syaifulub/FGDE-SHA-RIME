# FGDE-SHA-RIME — Frozen V2 Validation Update

## Recommended positioning

The final method should be presented as a **two-stage development followed by frozen independent validation**. The original six-problem study (G1, G3, G4, G6, G8, G9) and its 13-configuration ablation remain a **base-framework development analysis**. G5 and G13 were subsequently used **only as diagnostic development problems** to address the weakness of the original single-equality interpolation on multiple simultaneous equalities. After the added projection mechanism was fixed, the implementation and all parameters were frozen. The frozen method was then evaluated on ten previously unused CEC2006 problems: G12, G14, G15, G16, G17, G18, G19, G21, G22, and G23.

This framing avoids claiming that G5/G13 are independent evidence and keeps the unseen ten-problem suite genuinely confirmatory with respect to the frozen V2 implementation.

---

## New method subsection — Budget-neutral finite-difference multi-equality projection

### Budget-neutral multi-equality projection

The base restoration mechanism repairs one equality residual at a time by interpolation between population members located on opposite sides of an equality surface. This mechanism was effective on the single-equality development problem G3, but diagnostic experiments on G5 and G13 showed that one-at-a-time interpolation was insufficient for problems containing several simultaneous equality constraints. A budget-neutral finite-difference projection was therefore added to the final frozen implementation.

Let \(x^V\) denote the population member with the smallest aggregate constraint violation and let \(H(x)=[h_1(x),\ldots,h_p(x)]^\top\) collect its signed equality residuals. For each equality, a robust population scale is defined as

\[
s_k=\max\left\{\operatorname{median}_{i=1,\ldots,N}|h_k(x_i)|,\;10\delta_h\right\}.
\]

Population positions are expressed in normalized decision coordinates relative to \(x^V\),

\[
z_{i,d}=\frac{x_{i,d}-x^V_d}{x_d^U-x_d^L},
\]

and scaled residual differences are

\[
r_{i,k}=\frac{h_k(x_i)-h_k(x^V)}{s_k}.
\]

A least-squares population sensitivity model is fitted as

\[
B=\arg\min_B\|ZB-R\|_F^2.
\]

The sensitivity of decision variable \(d\) is measured by \(\|B_{d,:}\|_2\). At most six variables having the largest sensitivities are retained for the projection. For a selected dimension \(d_j\), a signed forward perturbation of magnitude \(\eta=10^{-3}\) of the variable range is evaluated, reversing direction when required by the upper bound. With \(\Delta z_j=\pm\eta\), the scaled finite-difference Jacobian is

\[
J_{k,j}=\frac{h_k(x^V+\Delta z_j(x_{d_j}^U-x_{d_j}^L)e_{d_j})/s_k-h_k(x^V)/s_k}{\Delta z_j}.
\]

The normalized projection step is obtained from the least-squares system

\[
\Delta z=\arg\min_{u}\|Ju+\widetilde H(x^V)\|_2,
\qquad
\widetilde H_k(x^V)=\frac{h_k(x^V)}{s_k}.
\]

To prevent unstable extrapolation, the step is rescaled when necessary so that

\[
\|\Delta z\|_\infty\le 0.5.
\]

The projected point is clipped to the decision bounds and evaluated once. Finite-difference probes are also valid evaluated candidates; among the probes and the projected point, the best candidate according to the same feasibility-preserving epsilon comparison is retained when it improves the least-violating incumbent.

The projection is activated only when equality constraints are present, the global violation-stagnation counter has reached four iterations, and the current iteration is a multiple of three. The mechanism is **budget neutral**. If \(c\) evaluations are used for finite-difference probes and the projected candidate, only \(N-c\) ordinary population trials are evaluated in that iteration. Consequently, every complete iteration still consumes exactly \(N=20\) evaluations and the total budget remains exactly 5,000 evaluations.

---

## New validation subsection — Diagnostic development and frozen unseen validation

### Diagnostic equality-constrained development

G5 and G13 were used as diagnostic development problems because both contain three equality constraints and were not part of the original six-problem base-development suite. No claim of independent validation is made from these two problems. Under the original FGDE-SHA-RIME configuration, neither G5 nor G13 produced a feasible final solution in the preliminary holdout experiment. After adding the budget-neutral multi-equality projection and fixing its trigger, finite-difference step, dimensional cap, and evaluation accounting, the final development check used the same 25 seeds, population size 20, 5,000-evaluation budget, and equality tolerance \(10^{-4}\). The frozen V2 configuration reached feasibility in **25/25 runs on G5 and 25/25 runs on G13**. Median first-feasible counts were 2,380 and 3,040 evaluations, respectively. Median absolute objective errors were 2.989795 on G5 and 0.654323 on G13.

The projection did not reduce feasibility on the original six development problems. V1 and V2 both remained feasible in all 150 runs. Their median errors were identical on G1, G4, G6, G8, and G9. On G3, the V2 median error increased from 0.470895 to 0.612684 while feasibility remained 25/25, illustrating that stronger multi-equality restoration does not necessarily improve objective refinement on every equality-constrained problem.

### Frozen independent validation suite

After the projection mechanism and all parameters were fixed, the implementation was evaluated without further retuning on ten previously unused CEC2006 problems: G12, G14, G15, G16, G17, G18, G19, G21, G22, and G23. These problems cover inequality-only formulations, multiple equality constraints, mixed inequality-equality constraints, dimensions from 3 to 22, and up to 19 simultaneous equalities. Each method used 25 paired seeds, population size 20, exactly 5,000 candidate evaluations, and equality tolerance \(10^{-4}\).

Five configurations were compared: DE-Deb, the RIME-derived comparator, CA-SHA-RIME-FR, the pre-projection FGDE-SHA-RIME base configuration, and the final frozen FGDE-SHA-RIME configuration with multi-equality projection. Feasibility was evaluated separately from objective error. Problem-level ranking used a feasibility-first lexicographic rule: higher feasibility rate was preferred first, followed by lower median absolute error among feasible runs, and then lower median violation when objective error was undefined. This ranking avoids the method-specific infeasibility penalty used in the original development benchmark.

### Table X. Frozen independent validation of the final FGDE-SHA-RIME configuration

| Problem | Feasible runs | Feasibility rate | Median absolute error | Median violation | Median first-feasible NFE |
|---|---:|---:|---:|---:|---:|
| G12 | 25/25 | 100% | 0.000000 | 0 | 20 |
| G14 | 23/25 | 92% | 4.455861 | 0 | 3,700 |
| G15 | 19/25 | 76% | 0.502749 | 0 | 3,280 |
| G16 | 16/25 | 64% | 0.176122 | 0 | 4,840 |
| G17 | 21/25 | 84% | 112.367868 | 0 | 3,880 |
| G18 | 25/25 | 100% | 0.224013 | 0 | 4,280 |
| G19 | 25/25 | 100% | 4.887018 | 0 | 20 |
| G21 | 7/25 | 28% | 149.707526 | 0.434029 | >5,000 for median run |
| G22 | 0/25 | 0% | undefined | 1.13587e7 | >5,000 |
| G23 | 10/25 | 40% | 395.368001 | 0.763200 | >5,000 for median run |

Across the ten unseen problems, the final configuration produced **171 feasible outcomes in 250 runs (68.4%)**, compared with **100/250 (40.0%)** for the pre-projection FGDE-SHA-RIME configuration. DE-Deb produced 178/250 feasible outcomes (71.2%), RIME 111/250 (44.4%), and CA-SHA-RIME-FR 100/250 (40.0%). Thus, the projection substantially improved equality-constrained feasibility relative to the base hybrid, but did not exceed DE-Deb in overall feasibility.

The final frozen configuration obtained the smallest mean feasibility-first lexicographic rank (**1.95**), followed by DE-Deb (2.60), RIME (3.10), the pre-projection FGDE-SHA-RIME configuration (3.35), and CA-SHA-RIME-FR (4.00). The Friedman test across the ten problem-level ranks was significant, \(\chi_F^2=10.1376\), \(p=0.03817\). Pairwise Wilcoxon comparisons against the frozen method did not remain significant after Holm correction for DE-Deb (adjusted \(p=0.3320\)), RIME (adjusted \(p=0.2578\)), or the pre-projection FGDE-SHA-RIME configuration (adjusted \(p=0.09375\)). The comparison with CA-SHA-RIME-FR remained significant (adjusted \(p=0.02344\)). These results support improved robustness relative to the earlier internal configurations but do not establish general superiority over DE-Deb or RIME.

As a descriptive aggregate across the 250 paired seed-problem runs, the final projection configuration was uniquely feasible in 71 cases, whereas the pre-projection configuration was never uniquely feasible. An exact sign/binomial calculation over these discordant cases gives \(p=8.47\times10^{-22}\); this aggregate is reported descriptively because the primary confirmatory inference treats the benchmark problem, rather than the individual run, as the statistical unit. The largest feasibility changes occurred on G14 (23/25 versus 4/25), G15 (19/25 versus 5/25), G17 (21/25 versus 0/25), G21 (7/25 versus 0/25), and G23 (10/25 versus 0/25). No method reached final feasibility on G22 under the 5,000-evaluation budget, demonstrating that highly constrained high-dimensional equality systems remain difficult.

---

## Recommended revised abstract

Continuous constrained optimization requires balancing objective improvement and constraint satisfaction under a limited evaluation budget. FGDE-SHA-RIME combines Differential Evolution and success-history adaptive RIME with feasibility-guided comparison, feasible archiving, epsilon relaxation, and stagnation-triggered restoration. A controlled base-framework ablation was first conducted on six CEC2006 development problems. Diagnostic experiments on the multi-equality problems G5 and G13 subsequently motivated a budget-neutral finite-difference projection that simultaneously reduces several equality residuals without increasing the 5,000-evaluation budget. After this mechanism and all parameters were frozen, independent validation was conducted on ten previously unused CEC2006 problems using 25 paired runs per method. The final configuration reached feasibility in 171/250 unseen runs, compared with 100/250 for the pre-projection configuration, and obtained the smallest mean feasibility-first rank (1.95). The overall Friedman test was significant (p=0.0382), but pairwise differences from DE-Deb and RIME were not significant after Holm correction. The results support improved multi-equality feasibility and broader robustness while showing persistent difficulty on highly constrained problems and do not establish general superiority.

---

## Recommended revised scope/limitations paragraph

The validation was deliberately separated into development and confirmatory stages. The original six CEC2006 problems and G5/G13 were used during method development and therefore do not provide independent evidence of generalization. The final projection mechanism and all control settings were frozen before evaluation on the ten-problem unseen suite. Independent results improved substantially relative to the pre-projection configuration, particularly on several multi-equality problems, but performance remained problem dependent. DE-Deb achieved a slightly higher overall feasibility rate on the unseen suite, pairwise problem-rank differences between the final method and DE-Deb or RIME were not significant after Holm correction, and none of the evaluated methods reached feasibility on G22 within 5,000 evaluations. The method therefore should be interpreted as a reproducible feasibility-guided framework with improved multi-equality handling rather than as a generally superior constrained optimizer. Evaluation on newer constrained benchmark suites, larger budgets, higher dimensions, mixed-variable problems, and real-world constrained applications remains necessary.

---

## Manuscript consistency changes required

1. Keep the original six-problem study explicitly labeled **base-development validation**.
2. Rename the original 13-configuration study **base-framework controlled ablation** rather than implying that it analyzes the final projection-enhanced implementation.
3. Add the budget-neutral projection subsection to *Method details* and update Algorithm 1 to include the projection trigger before ordinary candidate trials.
4. Add G5/G13 as **diagnostic development**, not independent validation.
5. Add the ten-problem frozen unseen validation as the main evidence for generalization.
6. Use feasibility-first lexicographic ranking for the unseen suite; retain Eq. (49) only for historical compatibility with the original development benchmark and keep the robustness check showing that its original mean-rank conclusion was unchanged.
7. Do not claim superiority over DE-Deb or RIME. The strongest supported claim is improved feasibility/generalization relative to the pre-projection FGDE-SHA-RIME configuration and CA-SHA-RIME-FR, with problem-dependent performance.
