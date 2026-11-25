import sys, csv, os, time
import pathlib
HERE=pathlib.Path(__file__).resolve().parent; sys.path.insert(0,str(HERE.parent)); sys.path.insert(0,str(HERE))
import FGDE_SHA_RIME_V2_Frozen as M, cec2006_remaining_defs as R, baselines as B
M.cec.REG.update(R.REG)
# usage: python run_remaining_problems.py <worker> <n_workers> <output_dir>
worker,nw,outdir=int(sys.argv[1]),int(sys.argv[2]),sys.argv[3]; os.makedirs(outdir,exist_ok=True)
probs=['g2','g7','g10','g11','g24']
ORIG=["FGDE-SHA-RIME-V2","FGDE-SHA-RIME","CA-SHA-RIME-FR","RIME","DE-Deb"]; EXT=['eps-DEag','ECHT-DE','L-SHADE-Deb']
jobs=[(n,p,m) for n in [5000,50000] for p in probs for m in ORIG+EXT]
for k,(nfe,p,m) in enumerate(jobs):
    if k%nw!=worker: continue
    fn=f'{outdir}/{nfe}__{p}__{m}.csv'
    if os.path.exists(fn): continue
    M.NFE=nfe; M.ITERS=(nfe-M.POP)//M.POP
    cps=sorted({40,220,520,1020,1520,2020,min(5000,nfe),nfe})
    t=time.time(); rows=[]; conv=[]
    for s in M.SEEDS:
        if m in EXT:
            r,h,used=B.run(m,p,s,nfe,cps); assert used==nfe
        else:
            r,h=M.run(m,p,s)
        rows.append(r); conv.extend(h)
    with open(fn+'.tmp','w',newline='') as f:
        w=csv.writer(f); w.writerow(["Seed","Problem","Method","Objective","Violation","Feasible","Absolute Error","First Feasible NFE"]); w.writerows(rows)
    with open(f'{outdir}/conv_{nfe}__{p}__{m}.csv','w',newline='') as f:
        w=csv.writer(f); w.writerow(["Seed","Problem","Method","NFE","Best Feasible Objective","Minimum Violation"]); w.writerows(conv)
    os.replace(fn+'.tmp',fn); print(nfe,p,m,round(time.time()-t,1),flush=True)
