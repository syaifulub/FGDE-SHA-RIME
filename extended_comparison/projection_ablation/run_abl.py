import sys, csv, os, time
import os; ROOT=os.path.abspath(os.path.join(os.path.dirname(__file__),'..','..')); sys.path.insert(0,ROOT); sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import fgde_abl as A, cec2006_development_defs as DV, cec2006_remaining_defs as RM, cec2006_diagnostic_defs as DG
for R in (DV.REG,RM.REG,DG.REG): A.cec.REG.update(R)
worker,nw,outdir=int(sys.argv[1]),int(sys.argv[2]),sys.argv[3]
VAR={'P':('FGDE-SHA-RIME-V2',{}),'B':('FGDE-SHA-RIME',{}),
     'P-random-dims':('FGDE-SHA-RIME-V2',{'proj_dims':'random'}),
     'P-single-equality':('FGDE-SHA-RIME-V2',{'proj_single':True}),
     'P-no-interpolation':('FGDE-SHA-RIME-V2',{'no_interp':True}),
     'Minimal-P':('FGDE-SHA-RIME-V2',{'de_only':True,'fixed_params':True}),
     'DE-Deb':('DE-Deb',{}),'DE-Deb+P':('DE-Deb',{'proj_on_de':True})}
EQ=['g5','g13','g3','g11','g14','g15','g17','g21','g22','g23']
NONEQ=['g12','g16','g18','g19','g1','g4','g6','g8','g9']
jobs=[(v,p) for p in EQ for v in VAR]+[(v,p) for p in NONEQ for v in ['P','B','Minimal-P','DE-Deb']]
for k,(v,p) in enumerate(jobs):
    if k%nw!=worker: continue
    fn=f'{outdir}/{v}__{p}.csv'
    if os.path.exists(fn): continue
    m,ov=VAR[v]; A.PRM=dict(A.DEFAULT,**ov); t=time.time(); rows=[]
    for s in A.SEEDS:
        r,_=A.run(m,p,s); rows.append([v]+r)
    with open(fn+'.tmp','w',newline='') as f:
        w=csv.writer(f); w.writerow(["Variant","Seed","Problem","Method","Objective","Violation","Feasible","Absolute Error","First Feasible NFE"]); w.writerows(rows)
    os.replace(fn+'.tmp',fn); print(v,p,round(time.time()-t,1),flush=True)
