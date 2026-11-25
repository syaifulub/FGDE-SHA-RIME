import sys, glob, json
import numpy as np, pandas as pd
from scipy.stats import friedmanchisquare, wilcoxon, binomtest, rankdata

FINAL = 'FGDE-SHA-RIME-V2'
ORDER = ['FGDE-SHA-RIME-V2', 'FGDE-SHA-RIME', 'CA-SHA-RIME-FR', 'RIME', 'DE-Deb', 'eps-DEag', 'ECHT-DE', 'L-SHADE-Deb']


def lex_ranks(df):
    out = []
    for p, g in df.groupby('Problem'):
        s = g.groupby('Method').agg(fr=('Feasible', 'mean'), me=('Absolute Error', 'median'), mv=('Violation', 'median')).reset_index()
        s['me'] = s['me'].fillna(np.inf)
        key = list(zip(-s.fr, s.me, s.mv))
        order = sorted(set(key)); r0 = {k: i for i, k in enumerate(order)}
        s['rank'] = rankdata([r0[k] for k in key]); s['Problem'] = p
        out.append(s)
    return pd.concat(out)


def holm(ps):
    items = sorted(ps.items(), key=lambda kv: kv[1]); k = len(items); adj = {}; run = 0
    for i, (c, p) in enumerate(items):
        run = max(run, min(1, (k - i) * p)); adj[c] = run
    return adj


def analyze(df, label):
    res = {'label': label}
    meths = [m for m in ORDER if m in df.Method.unique()]
    feas = df.groupby('Method').Feasible.agg(['sum', 'count'])
    res['feasible'] = {m: [int(feas.loc[m, 'sum']), int(feas.loc[m, 'count'])] for m in meths}
    r = lex_ranks(df); piv = r.pivot(index='Problem', columns='Method', values='rank')[meths]
    res['mean_rank'] = piv.mean().round(3).to_dict()
    fr = friedmanchisquare(*[piv[c] for c in meths]); res['friedman'] = [float(fr.statistic), float(fr.pvalue)]
    raw = {}
    for c in meths:
        if c == FINAL: continue
        d = piv[FINAL] - piv[c]
        raw[c] = float(wilcoxon(piv[FINAL], piv[c]).pvalue) if (d != 0).any() else 1.0
    res['wilcoxon_raw'] = raw; res['wilcoxon_holm'] = holm(raw)
    pv = df.pivot_table(index=['Problem', 'Seed'], columns='Method', values='Feasible')
    mc = {}
    for c in meths:
        if c == FINAL: continue
        a = int(((pv[FINAL] == 1) & (pv[c] == 0)).sum()); b = int(((pv[FINAL] == 0) & (pv[c] == 1)).sum())
        mc[c] = [a, b, float(binomtest(a, a + b, .5).pvalue) if a + b else 1.0]
    res['mcnemar_pooled'] = mc
    # per-problem McNemar, Holm over all problem x comparator tests
    per = {}
    for p in sorted(df.Problem.unique()):
        for c in meths:
            if c == FINAL: continue
            sub = pv.loc[p]
            a = int(((sub[FINAL] == 1) & (sub[c] == 0)).sum()); b = int(((sub[FINAL] == 0) & (sub[c] == 1)).sum())
            per[(p, c)] = [a, b, float(binomtest(a, a + b, .5).pvalue) if a + b else 1.0]
    adj = holm({k: v[2] for k, v in per.items()})
    res['mcnemar_problem_sig'] = {f'{k[0]}|{k[1]}': per[k][:2] + [round(adj[k], 6)] for k in per if adj[k] < 0.05}
    tab_f = df.pivot_table(index='Problem', columns='Method', values='Feasible', aggfunc='sum')[meths]
    tab_e = df.groupby(['Problem', 'Method'])['Absolute Error'].median().unstack()[meths]
    res['per_problem_feasible'] = tab_f.astype(int).to_dict()
    res['per_problem_median_error'] = tab_e.to_dict()
    res['ranks'] = piv.to_dict()
    return res


if __name__ == '__main__':
    import pathlib
    root = pathlib.Path(__file__).resolve().parents[1] / 'data'
    archived = pd.read_csv(root / 'frozen_validation' / 'raw_runs.csv')
    b5 = pd.read_csv(root / 'extended_comparison' / 'raw_runs_5000_external_baselines.csv')
    d50 = pd.read_csv(root / 'extended_comparison' / 'raw_runs_50000_all_methods.csv')
    out = {'5000': analyze(pd.concat([archived, b5]), '5000'), '50000': analyze(d50, '50000')}
    json.dump(out, open(root / 'extended_comparison' / 'statistical_summary.json', 'w'), indent=1, default=str)
    for k, r in out.items():
        print('=== budget', k)
        print('feasible', r['feasible']); print('mean rank', r['mean_rank']); print('friedman', r['friedman'])
        print('wilcoxon holm', {c: round(p, 4) for c, p in r['wilcoxon_holm'].items()})
        print('mcnemar pooled', r['mcnemar_pooled'])
        print('per-problem sig', r['mcnemar_problem_sig'])
        print(pd.DataFrame(r['per_problem_feasible']))
        print(pd.DataFrame(r['per_problem_median_error']).round(4))
