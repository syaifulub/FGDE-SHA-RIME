import pandas as pd, numpy as np, glob, json
from scipy.stats import binomtest, wilcoxon
d=pd.concat([pd.read_csv(f) for f in glob.glob('out/*.csv')])
def holm(ps):
    ps=np.asarray(ps,float); o=np.argsort(ps); m=len(ps); adj=np.empty(m); run=0
    for r,i in enumerate(o): run=max(run,min(1,(m-r)*ps[i])); adj[i]=run
    return adj
def mcn(a,b):
    x=int(((a==1)&(b==0)).sum()); y=int(((a==0)&(b==1)).sum())
    return x,y,(binomtest(x,x+y,.5).pvalue if x+y else 1.0)
EQ=['g5','g13','g3','g11','g14','g15','g17','g21','g22','g23']
V=['P','B','P-random-dims','P-single-equality','P-no-interpolation','Minimal-P','DE-Deb','DE-Deb+P']
piv=d.pivot_table(index=['Problem','Seed'],columns='Variant',values='Feasible')
err=d.pivot_table(index=['Problem','Seed'],columns='Variant',values='Absolute Error',dropna=False)
tab={}
for p in EQ+['g12','g16','g18','g19','g1','g4','g6','g8','g9']:
    row={}
    for v in V:
        if v in piv.columns and not piv.loc[p][v].isna().all():
            fe=int(piv.loc[p][v].sum()); me=err.loc[p][v].median() if p in err.index.get_level_values(0) else np.nan
            row[v]=(fe, None if np.isnan(me) else float(me))
    tab[p]=row
# McNemar vs P for projection variants on eq problems (incl. B), and DE-Deb+P vs DE-Deb, Minimal-P vs P
tests=[]
for p in EQ:
    a=piv.loc[p]
    for v in ['B','P-random-dims','P-single-equality','P-no-interpolation','Minimal-P']:
        x,y,pv=mcn(a['P'],a[v]); tests.append(('P vs '+v,p,x,y,pv))
    x,y,pv=mcn(a['DE-Deb+P'],a['DE-Deb']); tests.append(('DE-Deb+P vs DE-Deb',p,x,y,pv))
adj=holm([t[4] for t in tests])
res=[dict(cmp=t[0],problem=t[1],win=t[2],loss=t[3],p=t[4],pholm=float(a)) for t,a in zip(tests,adj)]
# pooled
pooled={}
for v in V:
    sub=d[(d.Variant==v)&(d.Problem.isin(EQ))]; pooled[v]=int(sub.Feasible.sum())
for grp,probs in [('frozen_eq',['g14','g15','g17','g21','g22','g23']),('diag',['g5','g13']),('single',['g3','g11'])]:
    for v in V:
        pooled[f'{grp}:{v}']=int(d[(d.Variant==v)&(d.Problem.isin(probs))].Feasible.sum())
# pooled McNemar over eq problems
pm={}
for v in ['B','P-random-dims','P-single-equality','P-no-interpolation','Minimal-P']:
    a=piv.loc[EQ]; pm['P vs '+v]=mcn(a['P'],a[v])
a=piv.loc[EQ]; pm['DE-Deb+P vs DE-Deb']=mcn(a['DE-Deb+P'],a['DE-Deb'])
# minimal across all 19 problems
ALL=EQ+['g12','g16','g18','g19','g1','g4','g6','g8','g9']
a=piv.loc[ALL]; pm['P vs Minimal-P (all)']=mcn(a['P'],a['Minimal-P']); pm['P vs DE-Deb (all)']=mcn(a['P'],a['DE-Deb'])
tot={v:int(d[(d.Variant==v)&(d.Problem.isin(ALL))].Feasible.sum()) for v in ['P','B','Minimal-P','DE-Deb']}
# error comparisons P vs Minimal-P on jointly feasible, per problem
errt=[]
for p in ALL:
    if p not in err.index.get_level_values(0): continue
    e=err.loc[p][['P','Minimal-P']].dropna()
    dd=e['P']-e['Minimal-P']; dd=dd[np.abs(dd)>1e-8*(1+np.abs(e['P']))]
    if len(dd)>=6:
        w=wilcoxon(dd); wp=(dd<0).sum(); errt.append(dict(problem=p,n=int(len(dd)),p=float(w.pvalue),P_better=int(wp),median_P=float(e['P'].median()),median_min=float(e['Minimal-P'].median())))
ea=holm([e['p'] for e in errt])
for e,x in zip(errt,ea): e['pholm']=float(x)
json.dump(dict(tab=tab,tests=res,pooled=pooled,pooled_mcn={k:list(v) for k,v in pm.items()},tot=tot,err_min=errt),open('abl_results.json','w'),indent=1)
print(json.dumps(pooled,indent=0)); print(pm); print(tot)
for r in res:
    if r['pholm']<0.05: print(r)
for e in errt: print(e)
