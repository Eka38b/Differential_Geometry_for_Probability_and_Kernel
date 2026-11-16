#!/usr/bin/env python3
"""Targeted synthesis of two normal solutions with unequal magnitudes.

Shared scalar data are constructed jointly, then frozen. This is a search for
a counterexample, not a random sample or an exhaustive existence test.
"""
import argparse,io,json,time,zipfile
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from alternative_search import Case,Fiber,comm,norm,magnitude_gap

def construct(n,seed,attempts=8,target=.2):
    rng=np.random.default_rng(seed)
    C=np.zeros((n,n));C[0,0]=1
    for j in range(1,n):C[n-j,j]=(-1.)**j
    ap=[(i,i) for i in range(1,n)]+[(i,i+1) for i in range(1,n-1)]+list(zip(*np.tril_indices(n,-1)))
    ap=[(i,j) for i,j in ap if j>0]
    hp=[(i,j) for i,j in zip(*np.tril_indices(n,-1)) if j>0]
    BA=np.zeros((len(ap),n,n));BH=np.zeros((len(hp),n,n))
    for k,(i,j) in enumerate(ap):BA[k,i,j]=1
    for k,(i,j) in enumerate(hp):BH[k,i,j]=1
    base=np.zeros((n,n));base[0,1]=base[-1,0]=1
    p=len(ap); basisP=np.r_[BA,BH];basisM=np.r_[BA,-BH]
    tri=np.triu_indices(n); weights=np.where(tri[0]==tri[1],1.,np.sqrt(2.))
    idxA=ap.index((2,1));idxH=p+hp.index((2,1))
    def matrices(z):
        A=base+np.einsum('k,kij->ij',z[:p],BA);H=np.einsum('k,kij->ij',z[p:],BH)
        return A+H,A-H
    def fun(z):
        P,M=matrices(z)
        return np.r_[comm(P)[tri]*weights,comm(M)[tri]*weights,4*z[idxA]*z[idxH]-target]
    def jac(z):
        P,M=matrices(z);rows=[]
        for D,B in [(P,basisP),(M,basisM)]:
            J=B.transpose(0,2,1)@D+D.T@B-B@D.T-D@B.transpose(0,2,1)
            rows.append(J[:,tri[0],tri[1]].T*weights[:,None])
        t=np.zeros(len(z));t[idxA]=4*z[idxH];t[idxH]=4*z[idxA]
        return np.vstack((*rows,t))
    lo=np.full(p+len(hp),-6.);hi=np.full_like(lo,6.)
    for k,(i,j) in enumerate(ap):
        if j==i+1:lo[k]=.15;hi[k]=5.
    best=None
    for attempt in range(attempts):
        z=rng.normal(scale=.4,size=len(lo))
        for k,(i,j) in enumerate(ap):
            if j==i+1:z[k]=rng.uniform(1.05,1.5)
        result=least_squares(fun,z,jac=jac,bounds=(lo,hi),max_nfev=1800,
                xtol=3e-13,ftol=3e-13,gtol=3e-13)
        error=norm(fun(result.x));P,M=matrices(result.x)
        record={'n':n,'seed':seed,'attempt':attempt,'error':error,'magnitude_gap':magnitude_gap(P,M),'nfev':result.nfev}
        if best is None or error<best[0]:best=(error,P,M,record)
        print(json.dumps(record),flush=True)
        if error<2e-9:break
    return C,best

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',default='paired_results');ap.add_argument('--attempts',type=int,default=8)
    ap.add_argument('--dimensions',nargs='+',type=int,default=[4,6,8,10]);args=ap.parse_args()
    out=Path(args.out);out.mkdir(exist_ok=True,parents=True);records=[]
    for n in args.dimensions:
        C,best=construct(n,202611072+n,args.attempts);error,P,M,record=best
        record['feasible_pair']=error<2e-9
        if record['feasible_pair']:
            case=Case(f'paired_n{n}','targeted_pair',C,P)
            f=Fiber(case);record['second_validation']=f.validate(M);record['linear_dof']=f.r
        np.savez_compressed(out/f'paired_n{n}.npz',C=C,Dplus=P,Dminus=M)
        records.append(record);(out/'summary.json').write_text(json.dumps(records,indent=2))
if __name__=='__main__':main()
