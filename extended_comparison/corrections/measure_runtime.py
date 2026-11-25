import sys, time, json, numpy as np
import pathlib
H=pathlib.Path(__file__).resolve(); sys.path.insert(0,str(H.parents[2])); sys.path.insert(0,str(H.parent)); sys.path.insert(0,str(H.parents[1]))
import fgde_corrections as F, baselines as B
F.PRM=dict(F.DEFAULT); F.NFE=5000
probs=['g12','g14','g15','g16','g17','g18','g19','g21','g22','g23']; seeds=[101,202,303,404,505]
stat={'calls':0,'evals':0,'time':0.0}
orig=F.fd_projection_step
def wrapped(*a,**k):
    t=time.perf_counter(); r=orig(*a,**k); stat['time']+=time.perf_counter()-t; stat['calls']+=1; stat['evals']+=r[1]; return r
F.fd_projection_step=wrapped
out={}
for m in ['FGDE-SHA-RIME-V2','FGDE-SHA-RIME','CA-SHA-RIME-FR','RIME','DE-Deb','eps-DEag','ECHT-DE','L-SHADE-Deb']:
    per={}
    for p in probs:
        ts=[]; s0=dict(stat)
        for s in seeds:
            t=time.perf_counter()
            if m in B.RUNNERS: B.run(m,p,s,5000,[5000])
            else: F.run(m,p,s)
            ts.append(time.perf_counter()-t)
        per[p]=dict(mean_s=float(np.mean(ts)), proj_calls=(stat['calls']-s0['calls'])/len(seeds), proj_evals=(stat['evals']-s0['evals'])/len(seeds), proj_time=(stat['time']-s0['time'])/len(seeds))
    out[m]=per; print(m, round(np.mean([v['mean_s'] for v in per.values()]),3), flush=True)
json.dump(out,open('runtime.json','w'),indent=1)
