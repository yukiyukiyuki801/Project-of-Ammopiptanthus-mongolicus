#!/usr/bin/env python3
"""Fetch explicitly selected pinned upstream software, not data or analysis jobs."""
import argparse,base64,hashlib,json,pathlib,subprocess,urllib.request
ROOT=pathlib.Path(__file__).resolve().parents[1]
def main():
 p=argparse.ArgumentParser(__doc__);p.add_argument('--pyrona',action='store_true');p.add_argument('--ecospat-functions',action='store_true');p.add_argument('--hierfstat-functions',action='store_true');a=p.parse_args()
 if not(a.pyrona or a.ecospat_functions or a.hierfstat_functions):p.error('Select a named software source')
 if a.pyrona:
  rec=next(x for x in json.loads((ROOT/'provenance/upstream_software.json').read_text())['repositories']if x['path'].endswith('/pyRona'));dst=ROOT/rec['path']
  if not dst.exists():dst.parent.mkdir(parents=True,exist_ok=True);subprocess.run(['git','clone',rec['url'],str(dst)],check=True);subprocess.run(['git','-C',str(dst),'checkout','--detach',rec['commit']],check=True)
  got=subprocess.check_output(['git','-C',str(dst),'rev-parse','HEAD'],text=True).strip()
  if got!=rec['commit']:raise ValueError('Existing pyRona checkout differs; not overwritten')
 if a.ecospat_functions:
  for rec in json.loads((ROOT/'provenance/ecospat_source_blobs.json').read_text()):
   dst=ROOT/rec['destination'];dst.parent.mkdir(parents=True,exist_ok=True)
   if dst.exists():data=dst.read_bytes()
   else:
    req=urllib.request.Request(rec['blob_api_url'],headers={'Accept':'application/vnd.github+json','User-Agent':'Ammopiptanthus-analysis-reproduction'});raw=json.load(urllib.request.urlopen(req,timeout=60));data=base64.b64decode(raw['content'])
   if hashlib.sha256(data).hexdigest()!=rec['sha256']:raise ValueError('Upstream source SHA256 mismatch')
   dst.write_bytes(data)
 if a.hierfstat_functions:
  out=ROOT/'runs/revision_run_01/outputs/fst_subset_comparison/primary_sources'
  manifest=json.loads((out/'SOURCE_MANIFEST.json').read_text())
  for rec in manifest['files']:
   dst=out/rec['file']
   data=dst.read_bytes() if dst.exists() else urllib.request.urlopen(rec['url'],timeout=60).read()
   if hashlib.sha256(data).hexdigest()!=rec['sha256']:raise ValueError('hierfstat source SHA256 mismatch')
   dst.write_bytes(data)
 print('Selected software fetched and pinned sources verified; no analyses run.')
if __name__=='__main__':main()
