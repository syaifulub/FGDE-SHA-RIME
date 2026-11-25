"""Frozen implementation with exposed constants plus two revision options:
  fix_credit=True : discard the sampled scale factor and attempt count when restoration replaces a trial.
  norm_viol=True  : order infeasible candidates by constraint violations normalized by their initial maxima
                    (feasibility still decided on the raw sum).
With DEFAULT settings, results are identical to FGDE_SHA_RIME_V2_Frozen.py."""
from __future__ import annotations
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))
import math, json, time
import numpy as np
import cec2006_unseen_defs as cec
DEFAULT=dict(eps_q=.8,eps_frac=.8,eps_cp=2,reinj=20,proj_trig=4,proj_freq=3,fd_step=1e-3,max_dims=6,dz_bound=.5,cauchy=.1,pbest=.2,cr=.9,soft_p=.5,hard_max=.3,restore_trig=8,sigma0=.10,lr=.1,pop=20,fix_credit=False,norm_viol=False)
PRM=dict(DEFAULT)

POP,NFE,EQ_TOL=20,5000,1e-4
ITERS=(NFE-POP)//POP
SEEDS=[101+101*i for i in range(25)]


WEIGHTS=None
def violations(problem,X):
    d=problem.evaluate(np.atleast_2d(X),return_values_of=['F','G','H'],return_as_dictionary=True)
    f=np.asarray(d.get('F'),float).reshape(-1); g=d.get('G'); h=d.get('H')
    parts=[]
    if g is not None and np.size(g): parts.append(np.maximum(np.asarray(g,float),0))
    H=np.asarray(h,float) if h is not None and np.size(h) else np.empty((len(f),0))
    if H.size: parts.append(np.maximum(np.abs(H)-EQ_TOL,0))
    V=np.hstack(parts) if parts else np.zeros((len(f),0))
    return f,V,H
def evaluate(problem,X):
    f,V,H=violations(problem,X)
    raw=V.sum(axis=1)
    if WEIGHTS is None: return f,raw,H
    # normalized violation orders infeasible candidates; feasibility is still decided on the raw sum
    cv=np.where(raw<=1e-12,raw,np.maximum(V@WEIGHTS,2e-12))
    return f,cv,H

def deb_better(f1,v1,f2,v2):
    return (v1<=1e-12<v2) or (v1<=1e-12 and v2<=1e-12 and f1<f2) or (v1>1e-12 and v2>1e-12 and v1<v2)

def eps_better(f1,v1,f2,v2,eps):
    exact1,exact2=v1<=1e-12,v2<=1e-12
    if exact1 and not exact2:return True
    if exact2 and not exact1:return False
    if exact1 and exact2:return f1<f2
    relaxed1,relaxed2=v1<=eps,v2<=eps
    return (relaxed1 and not relaxed2) or (relaxed1 and relaxed2 and f1<f2) or (not relaxed1 and not relaxed2 and v1<v2)

def bounded(X,xl,xu): return np.clip(X,xl,xu)

def fd_projection_step(p,X,f,cv,H,idx,span,xl,xu,step_rel=1e-3,max_dims=6,dzb=.5):
    q=H.shape[1]; n=X.shape[1]
    if q==0:return None,0,[]
    scale=np.median(np.abs(H),axis=0); scale=np.maximum(scale,10*EQ_TOL)
    Z=(X-X[idx])/np.where(span>0,span,1.0); R=(H-H[idx])/scale
    try:
        B=np.linalg.lstsq(Z,R,rcond=1e-8)[0]
        sens=np.linalg.norm(B,axis=1)
        md=min(n,max_dims,max(1,min(q+1,max_dims)))
        dims=np.argsort(sens)[::-1][:md]
    except Exception:
        dims=np.arange(min(n,max_dims))
    h0=H[idx]/scale; J=np.zeros((q,len(dims))); probes=[]
    for jj,d in enumerate(dims):
        sign=1.0
        if X[idx,d]+step_rel*span[d]>xu[d]:sign=-1.0
        xp=X[idx].copy(); xp[d]=np.clip(xp[d]+sign*step_rel*span[d],xl[d],xu[d])
        actual=(xp[d]-X[idx,d])/span[d] if span[d]>0 else 0.0
        if abs(actual)<1e-14:continue
        fp,vp,hp=evaluate(p,xp); fp=float(fp[0]); vp=float(vp[0]); hp=hp[0]
        probes.append((xp,fp,vp,hp)); J[:,jj]=(hp/scale-h0)/actual
    if not probes:return None,0,[]
    try:dz=np.linalg.lstsq(J,-h0,rcond=1e-8)[0]
    except Exception:return None,len(probes),probes
    if not np.all(np.isfinite(dz)):return None,len(probes),probes
    mx=float(np.max(np.abs(dz))) if len(dz) else 0.0
    if mx>dzb:dz*=dzb/mx
    y=X[idx].copy()
    for val,d in zip(dz,dims):y[d]+=val*span[d]
    y=np.clip(y,xl,xu); fy,vy,hy=evaluate(p,y)
    return (y,float(fy[0]),float(vy[0]),hy[0]),len(probes)+1,probes

def run(method,pname,seed):
    P=PRM; POP=P['pop']; ITERS=(NFE-POP)//POP
    p=cec.get_problem(pname); rng=np.random.default_rng(seed); n=p.n_var; xl,xu=p.xl,p.xu; span=xu-xl
    global WEIGHTS
    WEIGHTS=None
    X=xl+rng.random((POP,n))*span; f,cv,H=evaluate(p,X)
    if P['norm_viol']:
        _,V0,_=violations(p,X); mx=V0.max(axis=0)
        WEIGHTS=np.where(mx>0,1.0/np.where(mx>0,mx,1.0),1.0)
        f,cv,H=evaluate(p,X)
    eps0=float(np.quantile(cv,P['eps_q'])); mf_obj=mf_cv=.5; pm_obj=pm_cv=.5; p_de=.5
    hist=[]; best_feas=np.inf; archive_x=None; archive_f=np.inf; archive_h=None
    first_feas=NFE+1
    feas0=np.where(cv<=1e-12)[0]
    if len(feas0):
        ai=feas0[np.argmin(f[feas0])]; archive_x=X[ai].copy(); archive_f=float(f[ai]); archive_h=H[ai].copy(); first_feas=POP
    best_cv_seen=float(cv.min()); cv_stag=0
    is_fgde=method in ['FGDE-SHA-RIME','FGDE-SHA-RIME-V2']
    uses_ca=method in ['CA-SHA-RIME-FR','FGDE-SHA-RIME','FGDE-SHA-RIME-V2']
    for t in range(ITERS):
        feasible=np.where(cv<=1e-12)[0]; feas_ratio=len(feasible)/POP
        b=int(feasible[np.argmin(f[feasible])]) if len(feasible) else int(np.argmin(cv)); best=X[b].copy()
        eps=eps0*max(0,1-t/(P['eps_frac']*ITERS))**P['eps_cp'] if uses_ca else 0.0
        if uses_ca and archive_x is not None:
            best=archive_x.copy()
            if t>0 and t%P['reinj']==0:
                infeas=np.where(cv>1e-12)[0]
                if len(infeas):
                    wi=infeas[np.argmax(cv[infeas])]; X[wi],f[wi],cv[wi]=archive_x.copy(),archive_f,0.0
                    if H.shape[1]:H[wi]=archive_h.copy()
        proj_cost=0; process_indices=list(range(POP))
        if method=='FGDE-SHA-RIME-V2' and H.shape[1]>0 and cv_stag>=P['proj_trig'] and t%P['proj_freq']==0:
            idxp=int(np.argmin(cv)); final,proj_cost,probes=fd_projection_step(p,X,f,cv,H,idxp,span,xl,xu,step_rel=P['fd_step'],max_dims=P['max_dims'],dzb=P['dz_bound'])
            candidates=list(probes)
            if final is not None:candidates.append(final)
            bestcand=None
            for cand in candidates:
                if eps_better(cand[1],cand[2],f[idxp],cv[idxp],eps):
                    if bestcand is None or eps_better(cand[1],cand[2],bestcand[1],bestcand[2],eps):bestcand=cand
            if bestcand is not None:X[idxp],f[idxp],cv[idxp],H[idxp]=bestcand
            proj_cost=min(proj_cost,POP-1); process_indices=list(rng.permutation(POP)[:POP-proj_cost])
            feasible=np.where(cv<=1e-12)[0]; feas_ratio=len(feasible)/POP
            if len(feasible):
                b=int(feasible[np.argmin(f[feasible])]); best=archive_x.copy() if archive_x is not None else X[b].copy()
        suc_obj=[];gain_obj=[];suc_cv=[];gain_cv=[];used_obj=won_obj=used_cv=won_cv=0;used_de=won_de=used_rime=won_rime=0
        for i in process_indices:
            ids=rng.choice([j for j in range(POP) if j!=i],3,replace=False); op=None
            if method=='DE-Deb':
                donor=X[ids[0]]+.5*(X[ids[1]]-X[ids[2]]); mask=rng.random(n)<.9; mask[rng.integers(n)]=True
                y=np.where(mask,donor,X[i]);Fused=None;channel=None
            elif is_fgde and cv[i]<=1e-12:
                channel='obj';Fused=mf_obj+P['cauchy']*rng.standard_cauchy()
                while Fused<=0:Fused=mf_obj+P['cauchy']*rng.standard_cauchy()
                Fused=min(Fused,1.);p_de_eff=float(np.clip(.5*p_de+.5*(.1+.8*feas_ratio),.1,.9))
                if rng.random()<p_de_eff:
                    op='de';used_de+=1;used_obj+=1
                    feas_order=feasible[np.argsort(f[feasible])]
                    top=feas_order[:max(1,int(np.ceil(P['pbest']*len(feas_order))))]
                    pbest=X[int(rng.choice(top))]
                    donor=X[i]+Fused*(pbest-X[i])+Fused*(X[ids[0]]-X[ids[1]])
                    mask=rng.random(n)<P['cr'];mask[rng.integers(n)]=True;y=np.where(mask,donor,X[i])
                else:
                    op='rime';used_rime+=1;Fused=None
                    soft=(1-(t+1)/ITERS)*math.cos(math.pi*(t+1)/ITERS);y=X[i].copy();mask=rng.random(n)<P['soft_p']
                    y[mask]=best[mask]+soft*(2*rng.random(mask.sum())-1)*span[mask];hard=rng.random(n)<(P['hard_max']*(t+1)/ITERS);y[hard]=best[hard]
            else:
                soft=(1-(t+1)/ITERS)*math.cos(math.pi*(t+1)/ITERS);y=X[i].copy();mask=rng.random(n)<P['soft_p']
                y[mask]=best[mask]+soft*(2*rng.random(mask.sum())-1)*span[mask];hard=rng.random(n)<(P['hard_max']*(t+1)/ITERS);y[hard]=best[hard]
                channel='cv' if cv[i]>eps else 'obj';Fused=None
                mutate=False
                if method=='RIME':mutate=False
                elif method=='CA-SHA-RIME-FR' or is_fgde:mutate=rng.random()<(pm_cv if channel=='cv' else pm_obj)
                if mutate:
                    mem=mf_cv if channel=='cv' else mf_obj;Fused=mem+P['cauchy']*rng.standard_cauchy()
                    while Fused<=0:Fused=mem+P['cauchy']*rng.standard_cauchy()
                    Fused=min(Fused,1.);y=y+Fused*(X[ids[0]]-X[ids[1]])
                    if channel=='cv':used_cv+=1
                    else:used_obj+=1
            if uses_ca and cv[i]>1e-12 and cv_stag>=P['restore_trig']:
                if P['fix_credit'] and Fused is not None:
                    if channel=='cv':used_cv-=1
                    else:used_obj-=1
                    Fused=None
                repaired=False
                if H.shape[1]:
                    q=int(np.argmax(np.abs(H[i])));pos=np.where(H[:,q]>0)[0];neg=np.where(H[:,q]<0)[0]
                    if len(pos) and len(neg):
                        aidx=pos[np.argmin(cv[pos])];bidx=neg[np.argmin(cv[neg])];ha,hb=H[aidx,q],H[bidx,q]
                        if abs(ha-hb)>1e-14:
                            tau=ha/(ha-hb);y=X[aidx]+tau*(X[bidx]-X[aidx]);repaired=True
                if not repaired:
                    leader=X[int(np.argmin(cv))];sigma=P['sigma0']*(1-t/max(1,ITERS-1))+.005
                    y=X[i]+rng.uniform(.4,1.0)*(leader-X[i])+sigma*rng.normal(size=n)*span
            y=bounded(y,xl,xu);fy_arr,vy_arr,hy_arr=evaluate(p,y);fy,vy,hy=fy_arr[0],vy_arr[0],hy_arr[0]
            accept=eps_better(fy,vy,f[i],cv[i],eps) if uses_ca else deb_better(fy,vy,f[i],cv[i])
            if accept:
                oldf,oldv=f[i],cv[i];X[i],f[i],cv[i],H[i]=y,fy,vy,hy
                if op=='de' and oldv<=1e-12 and vy<=1e-12 and fy<oldf:won_de+=1
                if op=='rime' and oldv<=1e-12 and vy<=1e-12 and fy<oldf:won_rime+=1
                if Fused is not None and (method=='CA-SHA-RIME-FR' or is_fgde):
                    if channel=='cv' and oldv>vy:suc_cv.append(Fused);gain_cv.append(max(oldv-vy,1e-15));won_cv+=1
                    elif channel=='obj' and fy<oldf:suc_obj.append(Fused);gain_obj.append(max(oldf-fy,1e-15));won_obj+=1
        def update_mem(mem,vals,gains):
            if not vals:return mem
            v=np.asarray(vals);w=np.asarray(gains);w/=w.sum();return (1-P['lr'])*mem+P['lr']*float(np.sum(w*v*v)/np.sum(w*v))
        if method=='CA-SHA-RIME-FR' or is_fgde:
            mf_obj=update_mem(mf_obj,suc_obj,gain_obj);mf_cv=update_mem(mf_cv,suc_cv,gain_cv)
            pm_obj=float(np.clip((1-P['lr'])*pm_obj+P['lr']*(won_obj/max(used_obj,1)),.1,.9));pm_cv=float(np.clip((1-P['lr'])*pm_cv+P['lr']*(won_cv/max(used_cv,1)),.1,.9))
        if is_fgde and used_de+used_rime:
            sr_de=won_de/max(used_de,1);sr_rime=won_rime/max(used_rime,1);target=(sr_de+.01)/(sr_de+sr_rime+.02);p_de=float(np.clip((1-P['lr'])*p_de+P['lr']*target,.1,.9))
        feas=np.where(cv<=1e-12)[0]
        if len(feas):
            if first_feas>NFE:first_feas=(t+1)*POP+POP
            fi=feas[np.argmin(f[feas])];best_feas=min(best_feas,float(f[fi]))
            if f[fi]<archive_f:archive_x=X[fi].copy();archive_f=float(f[fi]);archive_h=H[fi].copy()
        cur=float(cv.min())
        if cur<best_cv_seen-1e-14:best_cv_seen=cur;cv_stag=0
        else:cv_stag+=1
        if (t+1) in [1,10,25,50,75,100,ITERS]:hist.append([seed,pname,method,(t+1)*POP+POP,best_feas,cur])
    WEIGHTS_SAVED=WEIGHTS; WEIGHTS=None
    if WEIGHTS_SAVED is not None:
        f,cv,H=evaluate(p,X)
    feas=np.where(cv<=1e-12)[0]
    if len(feas):idx=feas[np.argmin(f[feas])]
    else:idx=int(np.argmin(cv))
    optimum=float(np.asarray(p.pareto_front()).reshape(-1)[0])
    return [seed,pname,method,float(f[idx]),float(cv[idx]),int(cv[idx]<=1e-12),abs(float(f[idx])-optimum) if cv[idx]<=1e-12 else np.nan,first_feas],hist

if __name__=='__main__':
    for p in cec.REG:
        for m in ['DE-Deb','RIME','CA-SHA-RIME-FR','FGDE-SHA-RIME','FGDE-SHA-RIME-V2']:
            t=time.time();r,h=run(m,p,101);print(p,m,round(time.time()-t,3),r)
