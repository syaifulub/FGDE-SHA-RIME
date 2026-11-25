import sys, csv, os, time, json
import pathlib
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parent))
import fgde_parameterized as F
# usage: python run_sensitivity.py <worker> <n_workers> <output_dir>
worker,nw,outdir=int(sys.argv[1]),int(sys.argv[2]),sys.argv[3]; os.makedirs(outdir,exist_ok=True)
LEVELS=dict(eps_q=(.5,.95),eps_frac=(.5,1.0),eps_cp=(1,4),reinj=(10,40),proj_trig=(2,8),proj_freq=(1,6),fd_step=(1e-4,1e-2),max_dims=(3,10),dz_bound=(.25,1.0),cauchy=(.05,.2),pbest=(.1,.4),cr=(.5,1.0),soft_p=(.3,.7),hard_max=(.15,.6),restore_trig=(4,16),sigma0=(.05,.2),lr=(.05,.2),pop=(10,40))
jobs=[('default',None)]+[(k,v) for k,vs in LEVELS.items() for v in vs]
probs=['g12','g14','g15','g16','g17','g18','g19','g21','g22','g23']
for j,(k,v) in enumerate(jobs):
    if j%nw!=worker: continue
    fn=f'{outdir}/{k}__{v}.csv'
    if os.path.exists(fn): continue
    F.PRM=dict(F.DEFAULT)
    if k!='default': F.PRM[k]=v
    t=time.time(); rows=[]
    for p in probs:
        for s in F.SEEDS:
            r,_=F.run('FGDE-SHA-RIME-V2',p,s); rows.append([k,v]+r)
    with open(fn+'.tmp','w',newline='') as f:
        w=csv.writer(f); w.writerow(['Param','Value',"Seed","Problem","Method","Objective","Violation","Feasible","Absolute Error","First Feasible NFE"]); w.writerows(rows)
    os.replace(fn+'.tmp',fn); print(k,v,round(time.time()-t,1),flush=True)
