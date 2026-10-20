#!/usr/bin/env python3
"""Verify and restore the exact released analysis inputs/results; run no analyses."""
import argparse,concurrent.futures,gzip,hashlib,json,os,pathlib,shutil,tempfile
ROOT=pathlib.Path(__file__).resolve().parents[1]
def digest(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(2**20),b''):h.update(b)
 return h.hexdigest()
def safe(p):
 out=ROOT/p
 if pathlib.Path(p).is_absolute()or '..'in pathlib.Path(p).parts:raise ValueError('Unsafe manifest path')
 return out

def main():
 p=argparse.ArgumentParser(__doc__);p.add_argument('--verify-only',action='store_true');p.add_argument('--workers',type=int,default=4);a=p.parse_args();rows=json.loads((ROOT/'provenance/data_manifest.json').read_text())['files']
 assets={q['path']:q for row in rows for q in row['assets']}
 def check(q):
  f=safe(q['path'])
  if f.stat().st_size!=q['bytes']or digest(f)!=q['sha256']:raise ValueError('Data asset checksum mismatch: '+q['path'])
 with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers)as pool:list(pool.map(check,assets.values()))
 if a.verify_only:print(json.dumps({'status':'PASS','assets':len(assets),'logical_files':len(rows),'analyses_run':False}));return
 def restore(row):
  dst=safe(row['path']);dst.parent.mkdir(parents=True,exist_ok=True)
  if dst.exists():
   if digest(dst)!=row['sha256']:raise ValueError('Existing file differs; preserve it and choose a clean checkout: '+row['path'])
   return
  with tempfile.TemporaryDirectory(prefix='.restore-',dir=dst.parent)as td:
   archive=pathlib.Path(td)/'packed';output=pathlib.Path(td)/'restored'
   with archive.open('wb')as out:
    for q in row['assets']:
     with safe(q['path']).open('rb')as f:shutil.copyfileobj(f,out,2**20)
   if row['compression']=='gzip':
    with gzip.open(archive,'rb')as f,output.open('wb')as out:shutil.copyfileobj(f,out,2**20)
   else:archive.rename(output)
   if output.stat().st_size!=row['bytes']or digest(output)!=row['sha256']:raise ValueError('Restored bytes mismatch: '+row['path'])
   os.replace(output,dst)
 with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers)as pool:list(pool.map(restore,rows))
 # Add the legacy location needed by the executed numerical sources.
 run=ROOT/'runs/revision_run_01';run.mkdir(parents=True,exist_ok=True)
 ammo=run/'ammo';target=ROOT/'ammo_analysis_package_20261014/ammo'
 if not ammo.exists():ammo.symlink_to(os.path.relpath(target,ammo.parent),target_is_directory=True)
 # Final model summaries used by the S6/S13 drawing scripts.
 gam=ROOT/'analyses/driver_models/results';compat=run/'outputs/oct16_revision/all_driver_models_postprocess';compat.mkdir(parents=True,exist_ok=True)
 for f in gam.iterdir():
  if not f.is_file():continue
  out=compat/f.name.removesuffix('.gz')
  if out.exists():continue
  if f.suffix=='.gz':
   with gzip.open(f,'rb')as src,out.open('wb')as dst:shutil.copyfileobj(src,dst,2**20)
  else:out.symlink_to(os.path.relpath(f,out.parent))
 cfg={'output_root':str(ROOT/'data/mapping_work/outputs'),'existing_pca_scores':str(ROOT/'revision_oct17_mapping/inputs/PCA_scores_and_group_assignments.csv'),'existing_pca_variance':str(ROOT/'revision_oct17_mapping/inputs/PCA_variance.csv')}
 out=ROOT/'data/mapping_plot_config.json';out.write_text(json.dumps(cfg,indent=2))
 print(json.dumps({'status':'PASS','assets':len(assets),'restored_files':len(rows),'root':str(ROOT),'analyses_run':False}))
if __name__=='__main__':main()
