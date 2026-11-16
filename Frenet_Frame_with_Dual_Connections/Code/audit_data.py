"""Audit every stored solution against the MathOverflow constraints."""
import argparse,csv,json
from pathlib import Path
import numpy as np
from fiber import Fiber,validate,norm


def main():
    root=Path(__file__).resolve().parents[1]
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,default=root/'Data/validation.json')
    args=ap.parse_args();maximum={};counts=[];index=[];derivative_errors=[]
    for folder in ['01_zero_source','02_magnitude','03_shared_tangent_source']:
        base=root/'Data'/folder;record=json.loads((base/'results.json').read_text())
        roots_checked=0
        for k,row in enumerate(record['rows']):
            path=base/row['data_file'];z=np.load(path,allow_pickle=False)
            c,s,d0,e0=z['C'],z['S'],z['D0'],z['E0']
            assert min(np.min(z['a']),np.min(z['astar']))>1e-8
            assert np.allclose(np.diag(d0),z['d'],rtol=0,atol=1e-12)
            assert np.allclose(np.diag(e0),z['dstar'],rtol=0,atol=1e-12)
            assert np.allclose(np.diag(d0,1),z['a'],rtol=0,atol=1e-12)
            assert np.allclose(np.diag(e0,1),z['astar'],rtol=0,atol=1e-12)
            for d,e in zip(z['roots'],z['dual_roots']):
                errors=validate(c,s,d,e,d0,e0)
                for name,value in errors.items():maximum[name]=max(maximum.get(name,0.),value)
                if max(errors.values())>2e-9:raise AssertionError((path,errors))
                roots_checked+=1
            if folder=='03_shared_tangent_source' and row['affine_dimension']>0:
                f=Fiber(c,s,d0,e0);rng=np.random.default_rng(100+k)
                point=rng.normal(size=f.r)*.01;direction=rng.normal(size=f.r);direction/=norm(direction)
                h=1e-5
                approx=(f.fun(point+h*direction)-f.fun(point-h*direction))/(2*h)
                discrepancy=norm(approx-f.jac(point)@direction)
                derivative_errors.append(discrepancy)
                if discrepancy>1e-7:raise AssertionError(discrepancy)
            index.append(dict(survey=folder,case_id=row['case_id'],family=row['family'],n=len(c),
                              affine_dimension=row['affine_dimension'],saved_root_count=len(z['roots']),
                              nonzero_source=bool(norm(s)>1e-8),file=str(path.relative_to(root/'Data'))))
        counts.append(dict(survey=folder,inputs=len(record['rows']),saved_pairs=roots_checked))
    with (root/'Data/index.csv').open('w')as out:
        writer=csv.DictWriter(out,fieldnames=list(index[0]));writer.writeheader();writer.writerows(index)
    result=dict(date='2026-11-16',passed=True,acceptance_tolerance=2e-9,surveys=counts,
                inputs=sum(x['inputs']for x in counts),saved_pairs=sum(x['saved_pairs']for x in counts),
                maximum_errors=maximum,analytic_jacobian_direction_checks=len(derivative_errors),
                maximum_jacobian_discrepancy=max(derivative_errors,default=0.))
    args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
