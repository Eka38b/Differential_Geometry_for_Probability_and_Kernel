"""Replay frozen MathOverflow-compatible inputs without changing archived data."""
import argparse,csv,json,platform
from pathlib import Path
import numpy as np
import scipy
from fiber import Fiber,solve,save_case


def main():
    root=Path(__file__).resolve().parents[1]
    ap=argparse.ArgumentParser()
    ap.add_argument('--survey',choices=['01_zero_source','02_magnitude'],required=True)
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--starts',type=int,default=24)
    ap.add_argument('--limit',type=int)
    ap.add_argument('--alternative',action='store_true')
    args=ap.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    source=root/'Data'/args.survey
    records=json.loads((source/'results.json').read_text())['rows']
    if args.limit:records=records[:args.limit]
    rows=[]
    for i,record in enumerate(records):
        z=np.load(source/record['data_file'])
        if args.alternative:
            if np.linalg.norm(z['S'])>1e-12:raise ValueError('Alternative algorithms require zero source')
            from alternative_search import Case,solve_case
            (args.out/'inputs').mkdir(exist_ok=True)
            row=solve_case(Case(record['case_id'],record['family'],z['C'],z['D0'],{'seed':20261116+i}),args.out,args.starts)
            row.update(affine_dimension=row['linear_dof'],saved_root_count=row['roots_found'])
        else:
            f=Fiber(z['C'],z['S'],z['D0'],z['E0'])
            roots,row=solve(f,args.starts,20261116+i)
            row.update(case_id=record['case_id'],family=record['family'],n=len(z['C']),
                       data_file=f"inputs/{record['case_id']}.npz")
            save_case(args.out/row['data_file'],f,roots)
        rows.append(row)
        print(json.dumps({k:row.get(k)for k in ['case_id','affine_dimension','saved_root_count','maximum_validation_error']}),flush=True)
    (args.out/'results.json').write_text(json.dumps(dict(date='2026-11-16',replay=True,
        historical_run_reproduced_bit_for_bit=False,versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),rows=rows),indent=2)+'\n')


if __name__=='__main__':main()
