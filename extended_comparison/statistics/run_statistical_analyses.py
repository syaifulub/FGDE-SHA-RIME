import json, glob, numpy as np, pandas as pd
from scipy.stats import binomtest, wilcoxon, norm, rankdata
import sys, pathlib; H=pathlib.Path(__file__).resolve(); sys.path.insert(0,str(H.parents[1])); import analyze_extended_comparison as A
FINAL='FGDE-SHA-RIME-V2'
arch=pd.read_csv(str(H.parents[2]/'data/frozen_validation/raw_runs.csv'))
b5=pd.read_csv(str(H.parents[2]/'data/extended_comparison/raw_runs_5000_external_baselines.csv'))
d5=pd.concat([arch,b5]); d50=pd.read_csv(str(H.parents[2]/'data/extended_comparison/raw_runs_50000_all_methods.csv'))
def holm(d):
    items=sorted(d.items(),key=lambda kv:kv[1]); k=len(items); out={}; run=0
    for i,(c,p) in enumerate(items): run=max(run,min(1,(k-i)*p)); out[c]=run
    return out
def per_problem(d, comps, probs):
    fe={}; er={}; info={}
    for p in probs:
        for c in comps:
            a=d[(d.Problem==p)&(d.Method==FINAL)].set_index('Seed').sort_index(); b=d[(d.Problem==p)&(d.Method==c)].set_index('Seed').sort_index()
            g=int(((a.Feasible==1)&(b.Feasible==0)).sum()); l=int(((a.Feasible==0)&(b.Feasible==1)).sum())
            fe[(p,c)]=binomtest(g,g+l,.5).pvalue if g+l else 1.0
            rd=(a.Feasible.sum()-b.Feasible.sum())/25
            jf=(a.Feasible==1)&(b.Feasible==1); ea=a['Absolute Error'][jf]; eb=b['Absolute Error'][jf]
            diff=(eb-ea).values; tol=1e-8*(1+np.abs(ea.values))  # positive = final better
            diff=np.where(np.abs(diff)<=tol,0,diff); nz=diff[diff!=0]
            r=None; pe=None
            if len(nz)>=1:
                rk=rankdata(np.abs(nz)); wp=rk[nz>0].sum(); wm=rk[nz<0].sum(); r=(wp-wm)/(wp+wm)
                if len(nz)>=6: pe=wilcoxon(nz).pvalue
            if pe is not None: er[(p,c)]=pe
            info[(p,c)]=dict(gain=g,loss=l,RD=rd,n_joint=int(jf.sum()),n_nonzero=int(len(nz)),r=r)
    hf=holm(fe); he=holm(er)
    for k in info: info[k]['p_feas_holm']=hf[k]; info[k]['p_err_holm']=he.get(k)
    return info
def posthoc(mean_ranks, N, control):
    k=len(mean_ranks); se=np.sqrt(k*(k+1)/(6*N)); ps={}; zs={}
    for m,r in mean_ranks.items():
        if m==control: continue
        z=(r-mean_ranks[control])/se; zs[m]=z; ps[m]=2*norm.sf(abs(z))
    return {m:dict(z=zs[m],p=ps[m],p_holm=h) for m,h in holm(ps).items()}, se
probs=['g12','g14','g15','g16','g17','g18','g19','g21','g22','g23']
comps=['FGDE-SHA-RIME','CA-SHA-RIME-FR','RIME','DE-Deb','eps-DEag','ECHT-DE','L-SHADE-Deb']
out={}
out['pp5']={f'{p}|{c}':v for (p,c),v in per_problem(d5,comps,probs).items()}
out['pp50']={f'{p}|{c}':v for (p,c),v in per_problem(d50,comps,probs).items()}
# Friedman post hoc
dev_ranks={'FGDE-SHA-RIME':1.750,'DE-Deb':2.583,'RIME-M':3.917,'SHA-F':4.250,'CA-SHA-RIME-FR':4.667,'SHA-RIME':5.750,'CA-SHA-RIME':6.000,'RIME':7.083}
out['ph_dev'],se=posthoc(dev_ranks,6,'FGDE-SHA-RIME')
r5=A.lex_ranks(arch).pivot(index='Problem',columns='Method',values='rank').mean().to_dict()
out['ph_frozen5'],_=posthoc(r5,10,FINAL)
for b,d in [('ph_ext5',d5),('ph_ext50',d50)]:
    mr=A.lex_ranks(d).pivot(index='Problem',columns='Method',values='rank').mean().to_dict(); out[b],_=posthoc(mr,10,FINAL)
# projection restricted to equality problems
eq=['g14','g15','g17','g21','g22','g23']
for b,d in [('5000',d5),('50000',d50)]:
    info=per_problem(d,['FGDE-SHA-RIME'],eq)
    # holm within the six feasibility tests only
    fp={}
    for p in eq:
        a=d[(d.Problem==p)&(d.Method==FINAL)].set_index('Seed').sort_index(); c=d[(d.Problem==p)&(d.Method=='FGDE-SHA-RIME')].set_index('Seed').sort_index()
        g=int(((a.Feasible==1)&(c.Feasible==0)).sum()); l=int(((a.Feasible==0)&(c.Feasible==1)).sum()); fp[p]=binomtest(g,g+l,.5).pvalue if g+l else 1.0
        info[(p,'FGDE-SHA-RIME')].update(fa=int(a.Feasible.sum()),fb=int(c.Feasible.sum()),p_raw=fp[p])
    hh=holm(fp)
    for p in eq: info[(p,'FGDE-SHA-RIME')]['p_holm6']=hh[p]
    sub=d[d.Problem.isin(eq)&d.Method.isin([FINAL,'FGDE-SHA-RIME'])]
    rk=A.lex_ranks(sub).pivot(index='Problem',columns='Method',values='rank')
    w=wilcoxon(rk[FINAL],rk['FGDE-SHA-RIME']) if (rk[FINAL]!=rk['FGDE-SHA-RIME']).any() else None
    pv=sub.pivot_table(index=['Problem','Seed'],columns='Method',values='Feasible'); g=int(((pv[FINAL]==1)&(pv['FGDE-SHA-RIME']==0)).sum()); l=int(((pv[FINAL]==0)&(pv['FGDE-SHA-RIME']==1)).sum())
    out['proj'+b]=dict(per={p:info[(p,'FGDE-SHA-RIME')] for p in eq}, problem_wilcoxon=None if w is None else float(w.pvalue), better=int((rk[FINAL]<rk['FGDE-SHA-RIME']).sum()), pooled=[g,l,binomtest(g,g+l,.5).pvalue])
json.dump(out,open(str(H.parents[2]/'data/statistics/statistical_analyses.json'),'w'),indent=1,default=str)
for k in ['ph_dev','ph_frozen5','ph_ext5','ph_ext50']:
    print(k,{m:(round(v['z'],2),round(v['p_holm'],4)) for m,v in out[k].items()})
for b in ['5000','50000']:
    print('proj',b,{p:(v['fa'],v['fb'],round(v['p_holm6'],6),v['RD']) for p,v in out['proj'+b]['per'].items()}, out['proj'+b]['problem_wilcoxon'], out['proj'+b]['better'], out['proj'+b]['pooled'])
pp=pd.DataFrame(out['pp5']).T; pp['sigF']=pp.p_feas_holm.astype(float)<.05; pp['sigE']=pp.p_err_holm.astype(float)<.05
print(pp[['RD','r','n_joint','p_feas_holm','p_err_holm']].to_string())
print('5k sig feas', pp.sigF.sum(), 'sig err', pp.sigE.sum())
pp50=pd.DataFrame(out['pp50']).T; print('50k sig feas',(pp50.p_feas_holm.astype(float)<.05).sum(),'sig err',(pd.to_numeric(pp50.p_err_holm)<.05).sum())
