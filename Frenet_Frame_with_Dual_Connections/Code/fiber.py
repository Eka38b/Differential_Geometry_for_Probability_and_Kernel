"""Fixed-data solver for the problem in MathOverflow question 516790.

All inputs have real orthogonal C fixing e1, skew S with S e1=0,
fixed diagonals and positive superdiagonals, and two normal Hessenberg matrices.
"""
from pathlib import Path
import json
import numpy as np
from scipy.optimize import least_squares


def norm(a):
    return float(np.linalg.norm(a))


def normal(d):
    return d @ d.T - d.T @ d


def validate(c, s, d, e, d0, e0):
    n=len(c);eye=np.eye(n);u=eye[:,0]
    return dict(orthogonality=norm(c@c.T-eye),common_tangent=norm(c@u-u),
                skew_source=norm(s+s.T)/max(1.,norm(s)),source_tangent=norm(s@u)/max(1.,norm(s)),
                dual_relation=norm(d.T+c@e@c.T-s)/max(1.,norm(d),norm(e),norm(s)),
                normality_d=norm(normal(d))/max(1.,norm(d)**2),
                normality_e=norm(normal(e))/max(1.,norm(e)**2),
                fixed_d=norm(np.triu(d-d0))/max(1.,norm(d0)),
                fixed_e=norm(np.triu(e-e0))/max(1.,norm(e0)),
                hessenberg_d=norm(np.triu(d,2))/max(1.,norm(d)),
                hessenberg_e=norm(np.triu(e,2))/max(1.,norm(e)))


class Fiber:
    def __init__(self,c,s,d0,e0=None,rcond=1e-10):
        self.C=c;self.S=s;self.D0=d0;self.n=n=len(c)
        self.E0=self.dual(d0) if e0 is None else e0
        if min(np.min(np.diag(d0,1)),np.min(np.diag(self.E0,1)))<=1e-8:
            raise ValueError('Strictly positive superdiagonals required')
        if max(validate(c,s,d0,self.dual(d0),d0,self.E0).values())>2e-9:
            raise ValueError('The supplied witness fails the public problem constraints')
        positions=list(zip(*np.tril_indices(n,-1)))
        b=np.zeros((len(positions),n,n))
        for k,(i,j) in enumerate(positions):b[k,i,j]=1.
        y=-c.T[None]@b.transpose(0,2,1)@c[None]
        tri=np.triu_indices(n)
        prescribed=y[:,tri[0],tri[1]].T
        sc=max(1.,norm(s))
        source=np.array([(s@(x+x.T)-(x+x.T)@s).ravel()/sc for x in b]).T
        linear=np.vstack([prescribed,source])
        _,sv,vh=np.linalg.svd(linear,full_matrices=False)
        self.linear_sv=sv
        rank=int(np.sum(sv>rcond*max(1.,sv[0]))) if len(sv) else 0
        self.basis=np.einsum('pk,pij->kij',vh[rank:].T,b)
        self.r=len(self.basis)
        self.tri=tri;self.weights=np.where(tri[0]==tri[1],1.,np.sqrt(2.))
        self.scale=max(1.,norm(d0)**2)
        self.rank_sweep={str(t):int(np.sum(sv>t*max(1.,sv[0]))) for t in [1e-8,1e-10,1e-12]}

    def dual(self,d):
        return self.C.T@(self.S-d.T)@self.C

    def matrix(self,z):
        return self.D0+np.einsum('k,kij->ij',z,self.basis)

    def fun(self,z):
        return normal(self.matrix(z))[self.tri]*self.weights/self.scale

    def jac(self,z):
        d=self.matrix(z);x=self.basis;xt=x.transpose(0,2,1)
        j=x@d.T+d@xt-xt@d-d.T@x
        return j[:,self.tri[0],self.tri[1]].T*self.weights[:,None]/self.scale

    def errors(self,d):
        return validate(self.C,self.S,d,self.dual(d),self.D0,self.E0)


def solve(fiber,starts=24,seed=20261116):
    rng=np.random.default_rng(seed);roots=[];origins=[];calls=dict(trf=0,lm=0)
    radius=max(1.,norm(fiber.D0))
    def accept(z,origin):
        d=fiber.matrix(z)
        if max(fiber.errors(d).values())>2e-10:return
        if any(norm(d-x)<2e-6*max(1.,norm(d),norm(x)) for x in roots):return
        roots.append(d);origins.append(origin)
    accept(np.zeros(fiber.r),'supplied_witness')
    if fiber.r:
        for k in range(starts):
            method='trf' if k%2==0 else 'lm';calls[method]+=1
            z=rng.normal(size=fiber.r)*radius/np.sqrt(fiber.r)
            result=least_squares(fiber.fun,z,jac=fiber.jac,method=method,max_nfev=600,
                                 ftol=3e-13,xtol=3e-13,gtol=3e-13)
            accept(result.x,method)
    ranks=[];nullities=[];smallest=[]
    for d in roots:
        if not fiber.r:ranks.append(0);nullities.append(0);smallest.append(None);continue
        z=np.einsum('kij,ij->k',fiber.basis,d-fiber.D0)
        sv=np.linalg.svd(fiber.jac(z),compute_uv=False)
        rank=int(np.sum(sv>max(1e-10,1e-7*sv[0])))
        ranks.append(rank);nullities.append(fiber.r-rank)
        smallest.append(float(sv[-1]/sv[0]) if sv[0] else 0.)
    return roots,dict(affine_dimension=fiber.r,saved_root_count=len(roots),
                      root_origins=origins,solver_calls=calls,jacobian_ranks=ranks,
                      jacobian_nullities=nullities,jacobian_singular_ratios=smallest,
                      linear_rank_sweep=fiber.rank_sweep,
                      maximum_validation_error=max(max(fiber.errors(d).values())for d in roots))


def save_case(path,fiber,roots):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(path,C=fiber.C,S=fiber.S,D0=fiber.D0,E0=fiber.E0,
                        d=np.diag(fiber.D0),dstar=np.diag(fiber.E0),
                        a=np.diag(fiber.D0,1),astar=np.diag(fiber.E0,1),
                        basis=fiber.basis,roots=np.array(roots),
                        dual_roots=np.array([fiber.dual(d)for d in roots]),linear_singular_values=fiber.linear_sv)
