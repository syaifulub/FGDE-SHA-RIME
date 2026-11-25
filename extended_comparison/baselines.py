"""Re-implemented external baselines for the FGDE-SHA-RIME revision.

All methods share the evaluation routine, aggregate violation (Eq. 2, delta_h = 1e-4),
numerical feasibility threshold (v <= 1e-12), seeds, and exact function-evaluation
accounting used by FGDE_SHA_RIME_V2_Frozen.py. Every evaluated point, including
finite-difference probes, counts against the budget.

Methods
-------
eps-DEag      : epsilon-constrained DE with an archive and gradient-based mutation,
                following Takahama & Sakai (CEC 2010).
ECHT-DE       : DE with an ensemble of four constraint-handling techniques
                (feasibility rules, self-adaptive penalty, epsilon constraint,
                stochastic ranking), following Mallipeddi & Suganthan (2010).
L-SHADE-Deb   : L-SHADE (Tanabe & Fukunaga, 2014) with Deb's feasibility rules.
"""
from __future__ import annotations
import math, sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import numpy as np
import FGDE_SHA_RIME_V2_Frozen as M

FEAS = 1e-12


class Budget(Exception):
    pass


class Ctx:
    """Evaluation wrapper with exact budget accounting and convergence logging."""

    def __init__(self, problem, nfe, checkpoints, seed, pname, method):
        self.p, self.nfe, self.used = problem, nfe, 0
        self.cps = sorted(checkpoints); self.ci = 0
        self.best_feas = np.inf; self.min_v = np.inf; self.first_feas = nfe + 1
        self.hist = []; self.meta = (seed, pname, method)

    def remaining(self):
        return self.nfe - self.used

    def eval(self, X):
        X = np.atleast_2d(X)
        k = min(len(X), self.remaining())
        if k <= 0:
            raise Budget
        d = self.p.evaluate(X[:k])
        f = np.asarray(d['F'], float).reshape(-1)
        G = np.asarray(d['G'], float) if 'G' in d else np.zeros((k, 0))
        H = np.asarray(d['H'], float) if 'H' in d else np.zeros((k, 0))
        v = np.maximum(G, 0).sum(axis=1) + np.maximum(np.abs(H) - M.EQ_TOL, 0).sum(axis=1)
        self.last_G = G
        for j in range(k):
            self.used += 1
            if v[j] <= FEAS:
                if self.first_feas > self.nfe:
                    self.first_feas = self.used
                self.best_feas = min(self.best_feas, float(f[j]))
            self.min_v = min(self.min_v, float(v[j]))
            while self.ci < len(self.cps) and self.used >= self.cps[self.ci]:
                self.hist.append([*self.meta, self.cps[self.ci], self.best_feas, self.min_v]); self.ci += 1
        return f, v, H, k


def deb_key(f, v):
    """Sort key implementing Deb's rules (smaller is better)."""
    feas = v <= FEAS
    return np.where(feas, 0, 1), np.where(feas, f, v)


def deb_better(f1, v1, f2, v2):
    return M.deb_better(f1, v1, f2, v2)


def eps_better(f1, v1, f2, v2, eps):
    """Epsilon-level comparison (Takahama & Sakai): <=eps treated as feasible."""
    a1, a2 = (v1 <= eps), (v2 <= eps)
    if a1 and a2:
        return f1 < f2 if f1 != f2 else v1 < v2
    if v1 == v2:
        return f1 < f2
    return v1 < v2


def final_answer(ctx, X, f, v):
    feas = np.where(v <= FEAS)[0]
    idx = feas[np.argmin(f[feas])] if len(feas) else int(np.argmin(v))
    opt = float(np.asarray(ctx.p.pareto_front()).reshape(-1)[0])
    ok = v[idx] <= FEAS
    s, pname, method = ctx.meta
    return [s, pname, method, float(f[idx]), float(v[idx]), int(ok),
            abs(float(f[idx]) - opt) if ok else np.nan, ctx.first_feas]


def reflect(y, x, xl, xu):
    """Bound handling: midpoint between parent and violated bound."""
    lo, hi = y < xl, y > xu
    y = y.copy(); y[lo] = (xl[lo] + x[lo]) / 2; y[hi] = (xu[hi] + x[hi]) / 2
    return y


# ---------------------------------------------------------------------------
# eps-DEag
# ---------------------------------------------------------------------------

def run_epsdeag(ctx, rng):
    p = ctx.p; n = p.n_var; xl, xu = p.xl, p.xu; span = xu - xl
    N = 4 * n; Msz = 100 * n
    F, CR, Pg, Rg = 0.5, 0.9, 0.2, 3
    cpmin = 3.0
    A = xl + rng.random((Msz, n)) * span
    fA, vA, HA, k = ctx.eval(A)
    A, fA, vA, HA = A[:k], fA[:k], vA[:k], HA[:k]
    theta = max(1, int(0.2 * len(A)))
    eps0 = float(np.sort(vA)[theta - 1])
    # Control period Tc = 20% of the post-initialization budget (Tc = 1000 of ~5000
    # generations in the original CEC2010 setting), measured in function evaluations
    # because gradient-based mutation makes the cost per generation variable.
    used0 = ctx.used; Tc_fe = 0.2 * max(1, ctx.nfe - used0)
    cp = cpmin if eps0 <= 0 else max(cpmin, (-5 - math.log10(eps0)) / math.log10(0.05))
    # population = top N of the archive under the eps0-level comparison
    key = [(0 if vv <= eps0 else 1, ff if vv <= eps0 else vv) for ff, vv in zip(fA, vA)]
    order = sorted(range(len(A)), key=lambda i: key[i])
    A, fA, vA, HA = A[order], fA[order], vA[order], HA[order]
    Np = min(N, len(A))  # first Np archive slots form the population

    def eps_level():
        t = (ctx.used - used0) / Tc_fe
        return eps0 * (1 - t) ** cp if t < 1 else 0.0

    def child_of(i):
        r = rng.choice([j for j in range(len(A)) if j != i], 3, replace=False)
        v_ = A[r[0]] + F * (A[r[1]] - A[r[2]])
        y = A[i].copy(); j = rng.integers(n); L = 0
        while True:  # exponential crossover
            y[j] = v_[j]; j = (j + 1) % n; L += 1
            if L >= n or rng.random() >= CR:
                break
        return reflect(y, A[i], xl, xu)

    def grad_mutation(y, fy, vy, gy, hy):
        for _ in range(Rg):
            if vy <= FEAS:
                break
            idx_g = np.where(gy > 0)[0]
            dC = np.concatenate([gy[idx_g], hy])
            if dC.size == 0:
                break
            h0 = 1e-6 * np.maximum(span, 1e-12)
            Xp = np.repeat(y[None, :], n, axis=0)
            Xp[np.arange(n), np.arange(n)] = np.where(y + h0 <= xu, y + h0, y - h0)
            fp, vp, Hp, kk = ctx.eval(Xp)
            Gp = ctx.last_G
            if kk < n:
                raise Budget
            Cp = np.concatenate([Gp[:, idx_g], Hp], axis=1)
            step = (Xp - y)[np.arange(n), np.arange(n)]
            J = ((Cp - dC[None, :]) / step[:, None]).T  # (m x n)
            dx = -np.linalg.pinv(J) @ dC
            if not np.all(np.isfinite(dx)):
                break
            y = np.clip(y + dx, xl, xu)
            fyy, vyy, hyy, kk = ctx.eval(y)
            fy, vy, gy, hy = float(fyy[0]), float(vyy[0]), ctx.last_G[0], hyy[0]
        return y, fy, vy, hy

    cp_updated = False
    try:
        while True:
            if not cp_updated and (ctx.used - used0) >= 0.95 * Tc_fe:
                cp = 0.3 * cp + 0.7 * cpmin; cp_updated = True
            eps = eps_level()
            for i in range(Np):
                for attempt in range(2):  # parent generates another child if the first fails
                    y = child_of(i)
                    fy, vy, hy, _ = ctx.eval(y); gy = ctx.last_G[0]; fy, vy, hy = float(fy[0]), float(vy[0]), hy[0]
                    if vy > FEAS and rng.random() < Pg:
                        y, fy, vy, hy = grad_mutation(y, fy, vy, gy, hy)
                    if eps_better(fy, vy, fA[i], vA[i], eps) or (fy == fA[i] and vy == vA[i]):
                        A[i], fA[i], vA[i], HA[i] = y, fy, vy, hy
                        break
    except Budget:
        pass
    return final_answer(ctx, A[:Np], fA[:Np], vA[:Np])


# ---------------------------------------------------------------------------
# ECHT-DE
# ---------------------------------------------------------------------------

def run_echtde(ctx, rng):
    p = ctx.p; n = p.n_var; xl, xu = p.xl, p.xu; span = xu - xl
    Ns, F, CR = 50, 0.5, 0.9
    pops = []
    for _ in range(4):
        X = xl + rng.random((Ns, n)) * span
        f, v, H, k = ctx.eval(X)
        pops.append([X[:k], f[:k], v[:k]])
    gmax = max(1, ctx.remaining() // (4 * Ns))
    Tc = max(1, int(0.2 * gmax)); cp = 5.0
    v_init = pops[2][2]
    eps0 = float(np.sort(v_init)[max(0, int(0.2 * len(v_init)) - 1)])

    def select_sf(X, f, v, m):
        a, b = deb_key(f, v)
        o = np.lexsort((b, a))[:m]
        return o

    def select_eps(X, f, v, m, eps):
        rel = v <= eps
        o = np.lexsort((np.where(rel, f, v), np.where(rel, 0, 1)))[:m]
        return o

    def select_sr(X, f, v, m, pf):
        idx = list(range(len(f)))
        fl = f.tolist(); vl = v.tolist()
        for _ in range(len(idx)):
            swapped = False
            u = rng.random(len(idx) - 1).tolist()
            for j in range(len(idx) - 1):
                a, b = idx[j], idx[j + 1]
                both_feas = vl[a] <= FEAS and vl[b] <= FEAS
                if both_feas or u[j] < pf:
                    swap = fl[a] > fl[b]
                else:
                    swap = vl[a] > vl[b]
                if swap:
                    idx[j], idx[j + 1] = b, a; swapped = True
            if not swapped:
                break
        return np.array(idx[:m])

    def select_sp(X, f, v, m):
        feas = v <= FEAS; rf = feas.mean()
        fmin, fmax = f.min(), f.max()
        fn = (f - fmin) / (fmax - fmin) if fmax > fmin else np.zeros_like(f)
        vmax = v.max(); vn = v / vmax if vmax > 0 else np.zeros_like(v)
        d = vn if rf == 0 else np.sqrt(fn ** 2 + vn ** 2)
        Xp = np.zeros_like(v) if rf == 0 else vn
        Yp = np.where(feas, 0.0, fn)
        Fp = d + (1 - rf) * Xp + rf * Yp
        return np.argsort(Fp, kind='stable')[:m]

    t = 0
    try:
        while True:
            pf = 0.475 - (0.475 - 0.025) * min(1.0, t / max(1, gmax - 1))
            eps = eps0 * (1 - t / Tc) ** cp if t < Tc else 0.0
            offs = []
            for (X, f, v) in pops:
                m = len(X); Y = np.empty_like(X)
                for i in range(m):
                    r = rng.choice([j for j in range(m) if j != i], 3, replace=False)
                    donor = X[r[0]] + F * (X[r[1]] - X[r[2]])
                    mask = rng.random(n) < CR; mask[rng.integers(n)] = True
                    Y[i] = reflect(np.where(mask, donor, X[i]), X[i], xl, xu)
                fy, vy, Hy, k = ctx.eval(Y)
                offs.append([Y[:k], fy[:k], vy[:k]])
                if k < m:
                    raise Budget
            OY = np.vstack([o[0] for o in offs]); Of = np.concatenate([o[1] for o in offs]); Ov = np.concatenate([o[2] for o in offs])
            new = []
            for c, (X, f, v) in enumerate(pops):
                CX = np.vstack([X, OY]); Cf = np.concatenate([f, Of]); Cv = np.concatenate([v, Ov])
                m = len(X)
                if c == 0: o = select_sf(CX, Cf, Cv, m)
                elif c == 1: o = select_sp(CX, Cf, Cv, m)
                elif c == 2: o = select_eps(CX, Cf, Cv, m, eps)
                else: o = select_sr(CX, Cf, Cv, m, pf)
                new.append([CX[o], Cf[o], Cv[o]])
            pops = new; t += 1
    except Budget:
        pass
    AX = np.vstack([q[0] for q in pops]); Af = np.concatenate([q[1] for q in pops]); Av = np.concatenate([q[2] for q in pops])
    return final_answer(ctx, AX, Af, Av)


# ---------------------------------------------------------------------------
# L-SHADE with Deb's feasibility rules
# ---------------------------------------------------------------------------

def run_lshade_deb(ctx, rng):
    p = ctx.p; n = p.n_var; xl, xu = p.xl, p.xu; span = xu - xl
    Ninit = 18 * n; Nmin = 4; H = 6; pbest_rate = 0.11; arc_rate = 2.6
    X = xl + rng.random((Ninit, n)) * span
    f, v, Hh, k = ctx.eval(X); X, f, v = X[:k], f[:k], v[:k]
    MF = np.full(H, 0.5); MCR = np.full(H, 0.5); kpos = 0
    Arc = np.empty((0, n))
    try:
        while True:
            N = len(X)
            a, b = deb_key(f, v); order = np.lexsort((b, a))
            r = rng.integers(H, size=N)
            CRs = np.where(np.isnan(MCR[r]), 0.0, np.clip(MCR[r] + 0.1 * rng.standard_normal(N), 0, 1))
            Fs = np.empty(N)
            for i in range(N):
                Fi = MF[r[i]] + 0.1 * rng.standard_cauchy()
                while Fi <= 0:
                    Fi = MF[r[i]] + 0.1 * rng.standard_cauchy()
                Fs[i] = min(Fi, 1.0)
            npb = max(2, int(round(pbest_rate * N)))
            U = np.empty_like(X)
            pool = np.vstack([X, Arc]) if len(Arc) else X
            for i in range(N):
                pb = X[order[rng.integers(npb)]]
                r1 = rng.integers(N)
                while r1 == i:
                    r1 = rng.integers(N)
                r2 = rng.integers(len(pool))
                while r2 == i or r2 == r1:
                    r2 = rng.integers(len(pool))
                donor = X[i] + Fs[i] * (pb - X[i]) + Fs[i] * (X[r1] - pool[r2])
                mask = rng.random(n) < CRs[i]; mask[rng.integers(n)] = True
                U[i] = reflect(np.where(mask, donor, X[i]), X[i], xl, xu)
            fu, vu, Hu, k = ctx.eval(U)
            SF, SCR, W = [], [], []
            for i in range(k):
                if deb_better(fu[i], vu[i], f[i], v[i]) or (fu[i] == f[i] and vu[i] == v[i]):
                    if deb_better(fu[i], vu[i], f[i], v[i]):
                        both = vu[i] <= FEAS and v[i] <= FEAS
                        gain = (f[i] - fu[i]) if both else max(v[i] - vu[i], 0.0)
                        if gain <= 0:
                            gain = abs(f[i] - fu[i]) + 1e-15
                        SF.append(Fs[i]); SCR.append(CRs[i]); W.append(gain)
                        Arc = np.vstack([Arc, X[i]]) if len(Arc) else X[i][None, :]
                    X[i], f[i], v[i] = U[i], fu[i], vu[i]
            if k < N:
                raise Budget
            if SF:
                w = np.asarray(W) / np.sum(W); SF = np.asarray(SF); SCR = np.asarray(SCR)
                MF[kpos] = np.sum(w * SF ** 2) / np.sum(w * SF)
                if np.isnan(MCR[kpos]) or SCR.max() == 0:
                    MCR[kpos] = np.nan
                else:
                    MCR[kpos] = np.sum(w * SCR ** 2) / np.sum(w * SCR)
                kpos = (kpos + 1) % H
            # linear population size reduction
            Nnew = int(round(((Nmin - Ninit) / ctx.nfe) * ctx.used + Ninit))
            Nnew = max(Nmin, Nnew)
            if Nnew < len(X):
                a, b = deb_key(f, v); keep = np.lexsort((b, a))[:Nnew]
                X, f, v = X[keep], f[keep], v[keep]
            amax = int(round(arc_rate * len(X)))
            if len(Arc) > amax:
                Arc = Arc[rng.choice(len(Arc), amax, replace=False)]
    except Budget:
        pass
    return final_answer(ctx, X, f, v)


RUNNERS = {'eps-DEag': run_epsdeag, 'ECHT-DE': run_echtde, 'L-SHADE-Deb': run_lshade_deb}


def run(method, pname, seed, nfe, checkpoints):
    p = M.cec.get_problem(pname)
    ctx = Ctx(p, nfe, checkpoints, seed, pname, method)
    rng = np.random.default_rng(seed)
    res = RUNNERS[method](ctx, rng)
    assert ctx.used <= nfe
    return res, ctx.hist, ctx.used
