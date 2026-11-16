#!/usr/bin/env python3
"""Exact rational contraction certificates for algebraic counterexamples.

The examples are defined by rational quadratic systems, not by treating
rounded matrix entries as exact. Rational boxes certify existence. Further
interval bounds certify full-rank fixed-data normality Jacobians at both roots.
"""
import argparse,json
from pathlib import Path
import mpmath as mp
import numpy as np
import sympy as s
from scipy.linalg import qr

def model(n):
    if n==6:
        b,c,p,q,r=s.symbols('b c p q r');xs=s.symbols('x1 y1 z1 u1 x2 y2 z2 u2');vs=[b,c,p,q,r,*xs]
        def eq(x,y,z,u):return [1-b*b+p*p+2*x*x+y*y,b*p-b*q+q*x+u*x+y*z,
            -b*z+c*x+r*y+x*z,-b*u+c*y+q*x,b*b-c*c+u*u+z*z,c*q-c*r-q*z+r*z+u*z]
        def matrix(x,y,z,u):return s.Matrix([[0,1,0,0,0,0],[0,p,b,0,0,0],[0,x,q,c,0,0],
            [0,y,z,r,c,0],[0,x,u,z,q,b],[1,p,x,y,x,p]])
        P=matrix(*xs[:4]);M=matrix(*xs[4:]);F=eq(*xs[:4])+eq(*xs[4:])+[xs[0]**2-xs[4]**2-s.Rational(1,5)]
    elif n==8:
        b,c,d=s.symbols('b c d');xs=s.symbols('x1 y1 z1 u1 v1 x2 y2 z2 u2 v2');vs=[b,c,d,*xs]
        def eq(x,y,z,u,v):return [1-b*b+u*u+2*x*x,-b*y+c*x+u*z+v*x,-b*v+d*u+x*y,
            b*b-c*c+v*v+y*y,-c*z+d*y+v*z,c*c-d*d+v*v-y*y+z*z]
        def matrix(x,y,z,u,v):return s.Matrix([[0,1,0,0,0,0,0,0],[0,0,b,0,0,0,0,0],[0,x,0,c,0,0,0,0],
            [0,0,y,0,d,0,0,0],[0,u,0,z,0,d,0,0],[0,0,v,0,z,0,c,0],[0,x,0,v,0,y,0,b],[1,0,x,0,u,0,x,0]])
        P=matrix(*xs[:5]);M=matrix(*xs[5:]);F=eq(*xs[:5])+eq(*xs[5:])+[xs[0]**2-xs[5]**2-s.Rational(1,5)]
    else:raise ValueError(n)
    C=s.zeros(n);C[0,0]=1
    for j in range(1,n):C[n-j,j]=(-1)**j
    return vs,s.Matrix(F),P,M,C

def exact_structure_checks(n,vs,F,P,M,C):
    assert C*C.T==s.eye(n)
    assert C[:,0]==s.eye(n)[:,0]
    assert all(P[i,j]==0 and M[i,j]==0 for i in range(n) for j in range(i+2,n))
    E1=-C.T*P.T*C;E2=-C.T*M.T*C
    assert all(s.expand(E1[i,i+1]-P[i,i+1])==0 and s.expand(E2[i,i+1]-M[i,i+1])==0 for i in range(n-1))
    assert all(s.expand(P[i,j]-M[i,j])==0 and s.expand(E1[i,j]-E2[i,j])==0 for i in range(n) for j in range(i,n))
    assert all(E1[i,j]==0 and E2[i,j]==0 for i in range(n) for j in range(i+2,n))
    allowed={s.expand(f) for f in list(F[:-1])+[-f for f in F[:-1]]+[s.Integer(0)]}
    for D in [P,M]:
        assert all(s.expand(v) in allowed for v in D.T*D-D*D.T)
    return {'fixed_data_identical_symbolically':True,'dual_hessenberg_symbolically':True,
            'normality_reduces_to_certified_equations':True,
            'squared_entry_difference':'1/5'}

def rational(x):return s.Rational(str(mp.nstr(x,65)))
def maxrow(A):return max(sum(abs(A[i,j]) for j in range(A.cols)) for i in range(A.rows))
def textnum(x):return str(s.N(x,12))

def inverse_approx(J):
    MM=mp.matrix([[mp.mpf(str(s.N(J[i,j],110))) for j in range(J.cols)] for i in range(J.rows)])
    inv=MM**-1
    return s.Matrix([[rational(inv[i,j]) for j in range(J.cols)] for i in range(J.rows)])

def derivative_bound(J,vs):
    # Every J entry is affine, hence these derivatives are exact constants.
    return s.Matrix([sum(abs(s.diff(J[i,j],v)) for j in range(J.cols) for v in vs) for i in range(J.rows)])

def contraction(vs,F,center,Y,radius):
    rep=dict(zip(vs,center));J=F.jacobian(vs);Jc=J.subs(rep);Fc=F.subs(rep)
    beta=max(abs(v) for v in Y*Fc)
    delta=maxrow(s.eye(len(vs))-Y*Jc)
    eta=max(abs(v) for v in abs(Y)*derivative_bound(J,vs))
    lipschitz=delta+radius*eta
    passed=bool(lipschitz<1 and beta+lipschitz*radius<radius)
    return {'passed':passed,'beta':textnum(beta),'delta':textnum(delta),'eta':textnum(eta),
            'lipschitz_bound':textnum(lipschitz),'radius':str(radius),
            'self_map_margin':textnum(radius-beta-lipschitz*radius)}

def fiber_jacobian(D):
    n=D.rows;positions=[(k,l) for k in range(n) for l in range(1,k)]
    rows=[(i,j) for i in range(n) for j in range(i,n)]
    J=s.zeros(len(rows),len(positions))
    for row,(i,j) in enumerate(rows):
        for col,(k,l) in enumerate(positions):
            J[row,col]=(D[k,j] if i==l else 0)+(D[k,i] if j==l else 0)-(D[j,l] if i==k else 0)-(D[i,l] if j==k else 0)
    return J,positions

def isolation_certificate(vs,D,center,radius,stored=None):
    J,positions=fiber_jacobian(D);rep=dict(zip(vs,center));Jc=J.subs(rep)
    if stored is None:
        _,_,piv=qr(np.array(Jc,float).T,mode='economic',pivoting=True)
        selected=list(map(int,piv[:J.cols]));block=J[selected,:];Y=inverse_approx(block.subs(rep))
        stored={'rows':selected,'inverse':[[str(v) for v in row] for row in Y.tolist()]}
    else:
        selected=stored['rows'];block=J[selected,:];Y=s.Matrix([[s.Rational(v) for v in row] for row in stored['inverse']])
    delta=maxrow(s.eye(J.cols)-Y*block.subs(rep))
    eta=max(abs(v) for v in abs(Y)*derivative_bound(block,vs));bound=delta+radius*eta
    result={'passed':bool(bound<1),'rank':J.cols,'invertibility_bound':textnum(bound),'rows':selected}
    return result,stored

def initial(n,source):
    raw=np.load(source);P=raw['Dplus'];M=raw['Dminus']
    if n==6:return [P[1,2],P[2,3],P[1,1],P[2,2],P[3,3],P[2,1],P[3,1],P[3,2],P[4,2],M[2,1],M[3,1],M[3,2],M[4,2]]
    return [P[1,2],P[2,3],P[3,4],P[2,1],P[3,2],P[4,3],P[4,1],P[5,2],M[2,1],M[3,2],M[4,3],M[4,1],M[5,2]]

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',default='certificates');ap.add_argument('--source',default='paired_pilot')
    ap.add_argument('--verify-only',action='store_true');args=ap.parse_args();out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
    mp.mp.dps=110;summary=[]
    for n in [6,8]:
        vs,F,P,M,C=model(n);structure=exact_structure_checks(n,vs,F,P,M,C);path=out/f'certificate_n{n}.json'
        if args.verify_only:
            cert=json.loads(path.read_text());center=list(map(s.Rational,cert['center']));Y=s.Matrix([[s.Rational(v) for v in row] for row in cert['inverse']]);radius=s.Rational(cert['radius'])
        else:
            fun=s.lambdify(vs,list(F),'mpmath');Jfun=s.lambdify(vs,F.jacobian(vs),'mpmath')
            root=mp.findroot(fun,tuple(map(float,initial(n,Path(args.source)/f'paired_n{n}.npz'))),J=Jfun,tol=mp.mpf('1e-100'),maxsteps=80)
            center=[rational(v) for v in root];radius=s.Rational(1,10**35)
            Y=inverse_approx(F.jacobian(vs).subs(dict(zip(vs,center))))
            cert={'n':n,'variables':list(map(str,vs)),'equations':list(map(str,F)),
                  'center':list(map(str,center)),'radius':str(radius),'inverse':[[str(v) for v in row] for row in Y.tolist()]}
        existence=contraction(vs,F,center,Y,radius);assert existence['passed']
        center_rep=dict(zip(vs,center))
        super_lower_bounds=[s.simplify(P[i,i+1].subs(center_rep)-radius*sum(abs(s.diff(P[i,i+1],v)) for v in vs)) for i in range(n-1)]
        assert min(super_lower_bounds)>0
        iso=[];saved=[]
        for i,D in enumerate([P,M]):
            result,stored=isolation_certificate(vs,D,center,radius,cert.get('isolation',[None,None])[i])
            assert result['passed'];iso.append(result);saved.append(stored)
        cert['isolation']=saved;cert['existence_result']=existence;cert['structure']=structure;cert['isolation_results']=iso
        rep=dict(zip(vs,center));Pnum=np.array(P.subs(rep),float);Mnum=np.array(M.subs(rep),float)
        np.savez_compressed(out/f'counterexample_n{n}.npz',C=np.array(C,float),Dplus=Pnum,Dminus=Mnum,
            Eplus=np.array((-C.T*P.T*C).subs(rep),float),Eminus=np.array((-C.T*M.T*C).subs(rep),float))
        if not args.verify_only:path.write_text(json.dumps(cert,indent=2))
        row={'n':n,'existence':existence,'structure':structure,'local_isolation':iso,
             'minimum_superdiagonal':float(min(Pnum[i,i+1] for i in range(n-1)))}
        summary.append(row);print(json.dumps(row),flush=True)
    (out/'verification_summary.json').write_text(json.dumps(summary,indent=2))

if __name__=='__main__':main()
