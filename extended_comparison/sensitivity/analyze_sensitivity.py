import glob, json, pandas as pd, numpy as np
from scipy.stats import binomtest
fs=sorted(glob.glob('outputs/sensitivity/*.csv'))
d=pd.concat([pd.read_csv(f) for f in fs]); d['Value']=d['Value'].astype(str)
base=d[d.Param=='default'].set_index(['Problem','Seed'])
rows=[]; ps={}
for (k,v),g in d[d.Param!='default'].groupby(['Param','Value']):
    g=g.set_index(['Problem','Seed']).loc[base.index]
    a=int(((g.Feasible==1)&(base.Feasible==0)).sum()); b=int(((g.Feasible==0)&(base.Feasible==1)).sum())
    p=binomtest(a,a+b,.5).pvalue if a+b else 1.0
    # problem-wise median error comparison on problems where both have >=1 feasible
    me_g=g.groupby(level=0)['Absolute Error'].median(); me_b=base.groupby(level=0)['Absolute Error'].median()
    both=me_g.notna()&me_b.notna(); better=int((me_g[both]<me_b[both]-1e-12).sum()); worse=int((me_g[both]>me_b[both]+1e-12).sum())
    rows.append(dict(Param=k,Value=v,Feasible=int(g.Feasible.sum()),gain=a,loss=b,p=p,err_better=better,err_worse=worse,
                     per_problem=g.groupby(level=0).Feasible.sum().to_dict()))
    ps[(k,v)]=p
# Holm
items=sorted(ps.items(),key=lambda kv:kv[1]); K=len(items); adj={}; run=0
for i,(c,p) in enumerate(items): run=max(run,min(1,(K-i)*p)); adj[c]=run
for r in rows: r['p_holm']=adj[(r['Param'],r['Value'])]
res=dict(default=int(base.Feasible.sum()), default_per_problem=base.groupby(level=0).Feasible.sum().to_dict(), rows=rows, n_configs=len(rows))
json.dump(res,open('data/sensitivity/sensitivity_summary.json','w'),indent=1,default=str)
df=pd.DataFrame(rows).drop(columns='per_problem'); pd.set_option('display.width',200)
print('default',res['default'], res['default_per_problem']); print(df.sort_values(['Param','Value']).to_string())
