"""Export just the previously executed bridge reference identities (development only)."""
import hashlib,json
from pathlib import Path
root=Path('/N/slate/tg11/ebmcal_v2')
source=root/'runs/dev_calibration_a'
out=root/'analysis/bridge_reference';out.mkdir(parents=True,exist_ok=True)
seeds={31000200,31000201}
files=[]
for p in sorted(source.glob('cell_*/chunk_*/fit_records.jsonl')):
 rows=[]
 with p.open() as f:
  for s in f:
   try:r=json.loads(s)
   except json.JSONDecodeError:
    if s.endswith("\n"):raise
    break # Only a torn final append from an active writer may be left out.
   if r.get('seed') in seeds and (r.get('scheme') is None or r.get('perm_id') in {0,1,2}):rows.append(r)
 if rows:
  dest=out/p.relative_to(source);dest.parent.mkdir(parents=True,exist_ok=True)
  text=''.join(json.dumps(r,allow_nan=False)+'\n' for r in rows);dest.write_text(text)
  for name in ['manifest.json','progress.json','rows.jsonl']:
   original=p.with_name(name)
   if original.exists():dest.with_name(name).write_bytes(original.read_bytes())
  files.append({'path':str(dest.relative_to(out)),'records':len(rows),'sha256':hashlib.sha256(text.encode()).hexdigest()})
(out/'export_manifest.json').write_text(json.dumps({'run_id':'dev_calibration_a','selected_seeds':sorted(seeds),'permutation_ids':[0,1,2],'include_observed':True,'source':str(source),'files':files},indent=2)+'\n')
print(json.dumps({'files':len(files),'records':sum(x['records'] for x in files)}))
