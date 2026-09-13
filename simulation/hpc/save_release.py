"""Create a named, immutable local v2 checkpoint including a checksummed manifest."""
import argparse,hashlib,io,json,tarfile
from datetime import datetime,timezone
from pathlib import Path
ap=argparse.ArgumentParser();ap.add_argument('--label',required=True);ap.add_argument('--weekly-used',type=int,required=True)
ap.add_argument('--include-complete-mechanism-logs',action='store_true',help='Include the completed four-cell, 16-chunk-per-cell mechanism raw logs; other permutation logs remain excluded.')
args=ap.parse_args()
root=Path(__file__).resolve().parents[2]
if not args.label.replace('_','').replace('-','').isalnum():raise SystemExit('Use a simple unique release label')
out=root/'releases'/args.label
if out.exists():raise SystemExit('Release exists; choose a new label')
out.mkdir(parents=True)
folders=['scripts/v2','configs/v2','hpc/v2','report/v2','output/pdf','docs/review_20260908','docs/review_v2','runs/v2/snapshots']
folders += [str(p.relative_to(root)) for p in (root/'results').glob('v2_*') if p.is_dir()]
folders += ['runs/v2/hpc_control','runs/v2/campaign_versions','runs/v2/toolkit_versions']
files=[]
for folder in folders:
 for p in (root/folder).rglob('*'):
  if p.is_file() and '__pycache__' not in p.parts and p.suffix not in {'.aux','.log','.out','.pyc'}:files.append(p)
files += [p for p in (root/'docs').glob('*v2*.md')]
files += [root/'docs/模拟阶段审阅_20260908.md']
files += [root/'README.md', root/'results/README.md']
files += [root/'docs/HPC指南.md']
files += [root/'lit/前期文献调研.md', root/'docs/关键修正_20260906.md']
files += list((root/'runs/v2').glob('*.json'))
files += list((root/'runs/v2').glob('usage_observations.jsonl'))
for p in (root/'runs/v2/hpc_results').glob('*/cell_*/chunk_*/*'):
 if p.name in {'rows.jsonl','manifest.json','progress.json','replay_audit.jsonl'}:files.append(p)
mechanism_logs=[]
if args.include_complete_mechanism_logs:
 mechanism_root=root/'runs/v2/hpc_results/confirm_mechanism_b'
 mechanism_logs=sorted(mechanism_root.glob('cell_*/chunk_*/fit_records.jsonl'))
 expected={mechanism_root/f'cell_{cell}'/f'chunk_{chunk}'/'fit_records.jsonl' for cell in range(4) for chunk in range(16)}
 if set(mechanism_logs)!=expected or any(json.loads((p.parent/'progress.json').read_text()).get('done') is not True for p in mechanism_logs):
  raise SystemExit('Complete mechanism logs require every frozen cell/chunk and done progress; no partial raw inclusion')
 files += mechanism_logs
files=sorted(set(p for p in files if p.is_file()))
manifest={'label':args.label,'created_utc':datetime.now(timezone.utc).isoformat(),'weekly_used_percent':args.weekly_used,'status':'working_version_with_explicit_pending_experiments','files':{},'capture_policy':'Each entry is read once, hashed and archived from the same bytes. Cross-file capture is an interval, not an atomic run snapshot.','included_complete_mechanism_fit_logs':len(mechanism_logs),'excluded_large_evidence':'Permutation fit logs remain in the main workspace/HPC output roots. This checkpoint includes dataset summaries, frozen source/config, provenance, and active-job state; completed B0 mechanism raw logs are included only when explicitly requested and counted separately.'}
archive=out/(args.label+'.tar.gz')
with tarfile.open(archive,'w:gz') as tar:
 for p in files:
  payload=p.read_bytes(); name=str(p.relative_to(root))
  info=tar.gettarinfo(str(p),arcname=name);info.size=len(payload)
  manifest['files'][name]={'bytes':len(payload),'sha256':hashlib.sha256(payload).hexdigest()}
  tar.addfile(info,io.BytesIO(payload))
 manifest['capture_finished_utc']=datetime.now(timezone.utc).isoformat()
 payload=(json.dumps(manifest,indent=2,ensure_ascii=False)+'\n').encode()
 (out/'manifest.json').write_bytes(payload)
 info=tarfile.TarInfo('RELEASE_MANIFEST.json');info.size=len(payload)
 tar.addfile(info,io.BytesIO(payload))
receipt={'archive':str(archive.relative_to(root)),'sha256':hashlib.sha256(archive.read_bytes()).hexdigest(),'bytes':archive.stat().st_size,'files':len(files),'created_utc':manifest['created_utc']}
(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt))
