import sys, csv, os, time
import pathlib
H=pathlib.Path(__file__).resolve(); sys.path.insert(0,str(H.parents[2])); sys.path.insert(0,str(H.parent))
import fgde_corrections as F, cec2006_development_defs as DV, cec2006_remaining_defs as RM
F.cec.REG.update(DV.REG); F.cec.REG.update(RM.REG)
# usage: python run_corrections.py <worker> <n_workers> <output_dir>
worker,nw,outdir=int(sys.argv[1]),int(sys.argv[2]),sys.argv[3]; os.makedirs(outdir,exist_ok=True)
FROZ=['g12','g14','g15','g16','g17','g18','g19','g21','g22','g23']; REM=['g2','g7','g10','g11','g24']; DEV=['g1','g3','g4','g6','g8','g9']
METH=['FGDE-SHA-RIME-V2','FGDE-SHA-RIME','CA-SHA-RIME-FR']
jobs=[]
for p in FROZ+REM+DEV:
    for m in METH: jobs.append(('fix',5000,p,m))
for p in FROZ: jobs.append(('fix',50000,p,'FGDE-SHA-RIME-V2'))
for p in FROZ:
    for m in ['FGDE-SHA-RIME-V2','DE-Deb','RIME']: jobs.append(('norm',5000,p,m))
for p in FROZ:
    for m in ['FGDE-SHA-RIME-V2','DE-Deb']: jobs.append(('norm',50000,p,m))
for k,(var,nfe,p,m) in enumerate(jobs):
    if k%nw!=worker: continue
    fn=f'{outdir}/{var}__{nfe}__{p}__{m}.csv'
    if os.path.exists(fn): continue
    F.NFE=nfe; F.PRM=dict(F.DEFAULT, fix_credit=(var=='fix'), norm_viol=(var=='norm'))
    t=time.time(); rows=[]
    for s in F.SEEDS:
        r,_=F.run(m,p,s); rows.append(r)
    with open(fn+'.tmp','w',newline='') as f:
        w=csv.writer(f); w.writerow(["Seed","Problem","Method","Objective","Violation","Feasible","Absolute Error","First Feasible NFE"]); w.writerows(rows)
    os.replace(fn+'.tmp',fn); print(var,nfe,p,m,round(time.time()-t,1),flush=True)
