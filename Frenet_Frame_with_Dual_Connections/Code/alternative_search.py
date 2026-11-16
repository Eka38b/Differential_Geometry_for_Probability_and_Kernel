"""Zero-source algorithms used in the 7 November 2026 magnitude survey.

Adapted here for replay of the curated MathOverflow-compatible inputs.
"""
from __future__ import annotations
import argparse
import csv
import itertools
import json
import platform
import time
from dataclasses import dataclass, field
from pathlib import Path
import numpy as np
import scipy
from scipy.linalg import block_diag, null_space
from scipy.optimize import least_squares, linprog, brentq

SEED = 202611072
RCOND = 1e-10
ROOT_TOL = 2e-10
CLUSTER_TOL = 2e-6
MAG_TOL = 1e-7

def norm(x): return float(np.linalg.norm(x))
def comm(D): return D.T @ D - D @ D.T
def dual(D, C): return -C.T @ D.T @ C

@dataclass
class Case:
    name: str
    family: str
    C: np.ndarray
    D: np.ndarray
    metadata: dict = field(default_factory=dict)

class Fiber:
    """Witness-centered nullspace, with every strictly lower entry eligible."""
    def __init__(self, case):
        self.case=case; self.C=case.C; self.D0=case.D.copy(); self.E0=dual(self.D0,self.C)
        self.n=n=len(self.D0); self.positions=list(zip(*np.tril_indices(n,-1)))
        B=np.zeros((len(self.positions),n,n))
        for k,(i,j) in enumerate(self.positions): B[k,i,j]=1
        Y=-self.C.T[None] @ B.transpose(0,2,1) @ self.C[None]
        self.mask=np.triu_indices(n)
        M=Y[:,self.mask[0],self.mask[1]].T
        _,self.linear_sv,vh=np.linalg.svd(M,full_matrices=False)
        self.linear_rank=int(np.sum(self.linear_sv>RCOND*max(1.,self.linear_sv[0])))
        self.X=np.einsum('pk,pij->kij',vh[self.linear_rank:].T,B)
        self.r=len(self.X)
        self.tri=np.triu_indices(n); self.weights=np.where(self.tri[0]==self.tri[1],1.,np.sqrt(2.))
        self.scale=max(1.,norm(self.D0)**2)
    def matrix(self,z): return self.D0+np.einsum('k,kij->ij',z,self.X)
    def fun(self,z):
        F=comm(self.matrix(z)); return F[self.tri]*self.weights/self.scale
    def jac(self,z):
        D=self.matrix(z); X=self.X
        J=(X.transpose(0,2,1)@D+D.T@X-X@D.T-D@X.transpose(0,2,1))
        return J[:,self.tri[0],self.tri[1]].T*self.weights[:,None]/self.scale
    def validate(self,D):
        E=dual(D,self.C); a=self.D0; b=self.E0
        return {'normal':norm(comm(D))/max(1.,norm(D)**2),
                'dual_normal':norm(comm(E))/max(1.,norm(E)**2),
                'fixed_D':norm(np.triu(D-a))/max(1.,norm(a)),
                'fixed_E':norm(np.triu(E-b))/max(1.,norm(b)),
                'hessenberg_D':norm(np.triu(D,2))/max(1.,norm(D)),
                'hessenberg_E':norm(np.triu(E,2))/max(1.,norm(E)),
                'original':norm(D.T+self.C@E@self.C.T)/max(1.,norm(D)),
                'orthogonality':norm(self.C@self.C.T-np.eye(self.n))}
    def coordinates(self,D): return np.einsum('kij,ij->k',self.X,D-self.D0)

def arnoldi_normal(A,v,Q0,name,family,metadata):
    n=len(A); v=v/norm(v); basis=[]; w=v
    for j in range(n):
        if norm(w)<1e-10: raise ValueError('Krylov breakdown')
        u=w/norm(w); basis.append(u); w=A.T@u
        for _ in range(2):
            for q in basis: w-=q*(q@w)
    U=np.array(basis).T; D=U.T@A@U; Q=U.T@Q0@U
    S=np.diag((-1.)**np.arange(n)); C=Q@S
    D[np.triu_indices(n,2)]=0
    case=Case(name,family,C,D,metadata)
    if np.min(np.diag(D,1))<1e-8: raise ValueError('unreduced normalization failed')
    return case

def spectral(n,seed,moments=0):
    """Mixed real/complex eigenvalues, positive weights, nonconstant diagonals.

    moments=2 balances the second centered holomorphic moment.
    moments=3 balances the second AND third, by a positive linear program in
    squared imaginary parts. This avoids forcing a zero diagonal or a fixed
    signed-reversal C.
    """
    rng=np.random.default_rng(seed)
    nc=max(1,n//2-(seed%2 if n>=6 else 0)); nr=n-2*nc; m=nc+nr
    re=rng.normal(loc=rng.uniform(-.8,.8),scale=1.,size=m)
    wt=rng.uniform(.2,1.5,size=m); wt/=wt.sum()
    mu=wt@re; centered=re-mu
    b=rng.uniform(.2,1.6,size=nc)
    V=wt@(centered**2)
    if moments==2: b*=np.sqrt(V/np.sum(wt[:nc]*b*b))
    if moments==3:
        if nc<2: raise ValueError('too few complex planes')
        eq=np.array([wt[:nc],3*wt[:nc]*centered[:nc]])
        rhs=np.array([V,wt@(centered**3)])
        lower=.02*V/max(.1,wt[:nc].sum())
        sols=[]
        for _ in range(4):
            lp=linprog(rng.normal(size=nc),A_eq=eq,b_eq=rhs,bounds=[(lower,None)]*nc,method='highs')
            if lp.success: sols.append(lp.x)
        if not sols: raise ValueError('positive moment balancing infeasible')
        b=np.sqrt(np.mean(sols,axis=0))
    blocks=[]; signs=[]; v=[]
    for j in range(nc):
        blocks.append(np.array([[re[j],-b[j]],[b[j],re[j]]]))
        signs += [1.,-1.]; v += [np.sqrt(wt[j]),0.]
    for j in range(nc,m):
        blocks.append(np.array([[re[j]]]))
        signs += [1.]; v += [np.sqrt(wt[j])]
    A=block_diag(*blocks); Q0=np.diag(signs)
    family={0:'generic_spectral',2:'balanced_moment2',3:'balanced_moments23'}[moments]
    return arnoldi_normal(A,np.array(v),Q0,f'{family}_n{n}_s{seed}',family,
        {'seed':seed,'complex_planes':nc,'real_eigenvalues':nr,'balanced_moments':moments,
         'spectral_real_parts':re.tolist(),'spectral_imaginary_parts':b.tolist(),'weights':wt.tolist()})

def magnitude_gap(A,B): return norm(np.abs(A)-np.abs(B))/max(1.,norm(A),norm(B))

def complex_homotopy(fib,rng,projection_count=2):
    """Total-degree paths for random square projections of full quadratics.

    All endpoints must pass every original equation. Path failures and
    projection-extraneous roots are recorded; no completeness claim is made.
    """
    r=fib.r; candidates=[]; stats={'paths':0,'paths_completed':0,'real_endpoints':0,'extraneous_real_endpoints':0}
    if not 2<=r<=4: return candidates,stats
    m=len(fib.fun(np.zeros(r)))
    for projection in range(projection_count):
        Q,_=np.linalg.qr(rng.normal(size=(m,r))); P=Q.T
        gamma=np.exp(1j*rng.uniform(.25,2.8))
        def G(z): return P@fib.fun(z)
        def JG(z): return P@fib.jac(z)
        for signs in itertools.product([-1.,1.],repeat=r):
            stats['paths']+=1; z=np.array(signs,dtype=complex); t=0.; step=.035; ok=True; nsteps=0
            while t<1.-1e-13 and nsteps<1800:
                nsteps+=1; dt=min(step,1.-t)
                J=(1-t)*gamma*np.diag(2*z)+t*JG(z)
                try: velocity=np.linalg.solve(J,-(G(z)-gamma*(z*z-1)))
                except np.linalg.LinAlgError: ok=False;break
                zz=z+dt*velocity; tt=t+dt; converged=False
                for it in range(14):
                    H=(1-tt)*gamma*(zz*zz-1)+tt*G(zz)
                    if norm(H)<2e-10: converged=True;break
                    JJ=(1-tt)*gamma*np.diag(2*zz)+tt*JG(zz)
                    try: correction=np.linalg.solve(JJ,H)
                    except np.linalg.LinAlgError: break
                    zz-=correction
                    if not np.all(np.isfinite(zz)) or norm(zz)>1e7: break
                if converged and norm(zz)<1e6:
                    z=zz;t=tt;step=min(.08,step*(1.35 if it<5 else 1.))
                else:
                    step*=.4
                    if step<2e-7: ok=False;break
            if not ok or t<1.-1e-12: continue
            stats['paths_completed']+=1
            if norm(z.imag)<1e-6*max(1.,norm(z.real)):
                stats['real_endpoints']+=1
                if norm(fib.fun(z.real))<1e-7: candidates.append((z.real,'complex_homotopy'))
                else: stats['extraneous_real_endpoints']+=1
    return candidates,stats

def univariate(fib):
    z=np.zeros(1); f0=fib.fun(z); fp=fib.fun(np.ones(1)); fm=fib.fun(-np.ones(1))
    f1=(fp-fm)/2; f2=(fp+fm)/2-f0
    co=np.array([f2,f1,f0]).T; row=co[np.argmax(np.linalg.norm(co,axis=1))]
    while len(row)>1 and abs(row[0])<1e-13: row=row[1:]
    roots=np.roots(row)
    return [(np.array([float(q.real)]),'univariate_polynomial') for q in roots if abs(q.imag)<1e-8]

def solve_case(case,out,starts=96,homotopy=True):
    begin=time.time(); fib=Fiber(case); rng=np.random.default_rng(case.metadata.get('seed',SEED)+31)
    roots=[]; origins=[]; calls={'trf':0,'lm':0}; homstats={}; recoveries={}
    def accept(z,origin):
        D=fib.matrix(z); errors=fib.validate(D)
        if max(errors.values())>ROOT_TOL: return False
        match=next((i for i,X in enumerate(roots) if norm(D-X)<=CLUSTER_TOL*max(1.,norm(D),norm(X))),None)
        if match is None:
            match=len(roots);roots.append(D);origins.append(origin)
        recoveries.setdefault(origin,set()).add(match)
        return True
    accept(np.zeros(fib.r),'supplied_witness')
    if fib.r==1:
        for z,origin in univariate(fib): accept(z,origin)
    if fib.r:
        radius=max(1.,norm(case.D)); scale=radius/np.sqrt(fib.r)
        for j in range(starts):
            z0=rng.normal(scale=scale*rng.choice([.2,1.,2.]),size=fib.r)
            method='trf' if j%2==0 else 'lm'; calls[method]+=1
            res=least_squares(fib.fun,z0,jac=fib.jac,method=method,max_nfev=600,
                              ftol=3e-13,xtol=3e-13,gtol=3e-13)
            accept(res.x,method)
        if homotopy:
            candidates,homstats=complex_homotopy(fib,rng)
            for z,origin in candidates:
                res=least_squares(fib.fun,z,jac=fib.jac,method='lm',max_nfev=300,
                                  ftol=2e-14,xtol=2e-14,gtol=2e-14)
                accept(res.x,origin)
    gaps=[]
    for i,j in itertools.combinations(range(len(roots)),2):
        D,E=roots[i],roots[j]
        gaps.append({'roots':[i,j],'D':magnitude_gap(D,E),'Dstar':magnitude_gap(dual(D,case.C),dual(E,case.C))})
    maximum=max([max(x['D'],x['Dstar']) for x in gaps],default=0.)
    jac_ranks=[]; margins=[]
    for D in roots:
        sv=np.linalg.svd(fib.jac(fib.coordinates(D)),compute_uv=False) if fib.r else np.zeros(0)
        jac_ranks.append(int(np.sum(sv>max(1e-10,1e-7*sv[0]))) if len(sv) else 0)
        margins.append(float(sv[-1]/sv[0]) if len(sv) and sv[0] else None)
    row={'case_id':case.name,'family':case.family,'n':fib.n,'linear_dof':fib.r,
         'roots_found':len(roots),'origins':origins,'solver_calls':calls,'homotopy':homstats,
         'solver_recoveries':{k:sorted(v) for k,v in recoveries.items()},
         'magnitude_pairs':gaps,'maximum_magnitude_gap':maximum,'non_sign_candidate':maximum>MAG_TOL,
         'jacobian_ranks':jac_ranks,'jacobian_singular_ratios':margins,
         'nonconstant_diagonal':float(np.ptp(np.diag(case.D)))>1e-7,
         'minimum_superdiagonal':float(np.min(np.abs(np.diag(case.D,1)))),
         'maximum_validation_error':max(max(fib.validate(D).values()) for D in roots),
         'metadata':case.metadata,'seconds':time.time()-begin}
    np.savez_compressed(out/'inputs'/f'{case.name}.npz',C=case.C,D0=case.D,E0=fib.E0,
                        X=fib.X,roots=np.array(roots),linear_singular_values=fib.linear_sv)
    return row

