"""Bounded parallel batch dispatcher, one dataset owned by each worker."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import argparse,json
from concurrent.futures import ProcessPoolExecutor,as_completed
from run_dataset import run
HERE=Path(__file__).resolve().parent

def one(task):
 campaign,cell,seed_index,out=task
 return run(campaign,cell,seed_index,out)

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--campaign',choices=['power','robustness'],required=True);ap.add_argument('--batch-index',type=int,required=True);ap.add_argument('--workers',type=int,default=16);ap.add_argument('--output',required=True);a=ap.parse_args()
 mapping=json.loads((HERE/'batch_mapping.json').read_text())[a.campaign];batch=mapping[a.batch_index]
 with ProcessPoolExecutor(max_workers=a.workers) as pool:
  futures={pool.submit(one,(a.campaign,batch['cell_index'],i,a.output)):i for i in batch['seed_indices']}
  failures=[]
  for f in as_completed(futures):
   try:
    result=f.result()
    if not result['complete']:failures.append({'seed_index':futures[f],'unresolved':True})
   except Exception as e:failures.append({'seed_index':futures[f],'error':str(e)})
 if failures:raise RuntimeError(json.dumps(failures))
if __name__=='__main__':main()
