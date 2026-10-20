#!/usr/bin/env python3
"""Materialize portable reproduction paths; never run an analysis or submit jobs."""
import argparse,csv,hashlib,json,pathlib,shutil

def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def main():
 p=argparse.ArgumentParser(__doc__)
 p.add_argument('--repository-root',type=pathlib.Path,default=pathlib.Path(__file__).resolve().parent)
 p.add_argument('--work-root',required=True,type=pathlib.Path)
 p.add_argument('--raw-fastq-root',required=True,type=pathlib.Path)
 p.add_argument('--nanus-reference',required=True,type=pathlib.Path)
 p.add_argument('--mongolicus-reference',required=True,type=pathlib.Path)
 p.add_argument('--skip-reference-bytecheck',action='store_true',help='Configure absent downloads only; subsequent workflow still verifies references')
 a=p.parse_args();root=a.repository_root.resolve();base=root/'revision_oct17_mapping';work=a.work_root.resolve();raw=a.raw_fastq_root.resolve()
 tokens={'${REPOSITORY_ROOT}':str(root),'${WORK_ROOT}':str(work),'${RAW_FASTQ_ROOT}':str(raw)}
 def fill(x):
  if isinstance(x,dict):return {k:fill(v) for k,v in x.items()}
  if isinstance(x,list):return [fill(v) for v in x]
  if isinstance(x,str):
   for t,v in tokens.items():x=x.replace(t,v)
  return x
 def write(path,obj):path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(obj,indent=2)+'\n')
 for path in sorted(base.rglob('*.template.json')):write(path.with_name(path.name.replace('.template.json','.json')),fill(json.loads(path.read_text())))
 cfgpath=base/'parent/config.frozen.json';cfg=json.loads(cfgpath.read_text());auditpath=base/'references/REFERENCE_AUDIT.json';audit=json.loads(auditpath.read_text())
 for arm,source in [('nanus',a.nanus_reference.resolve()),('mongolicus',a.mongolicus_reference.resolve())]:
  expected=cfg['references'][arm]['compressed_sha256']
  if not a.skip_reference_bytecheck:
   if not source.is_file():raise FileNotFoundError(source)
   if sha(source)!=expected:raise ValueError(f'{arm}: exact reference SHA256 mismatch; do not substitute another assembly/compression')
  cfg['references'][arm]['source_gz']=str(source)
  matches=[x for x in audit['references'] if x['sha256']==expected];assert len(matches)==1;matches[0]['path']=str(source)
 write(auditpath,audit);write(cfgpath,cfg)
 text=(base/'workflow/samples.template.tsv').read_text()
 for t,v in tokens.items():text=text.replace(t,v)
 (base/'workflow/samples.tsv').write_text(text)
 # Lock portable copies to their new content/paths. Never reuse old output checkpoints.
 locked=[cfgpath,auditpath,base/'workflow/samples.tsv',base/'workflow/run_workflow.py',base/'workflow/pca_comparison.R',base/'parallel_call/call_shards.py',base/'parallel_call/config.json']+sorted((base/'parallel_call').glob('*.bed'))
 guard={'scope':'Fresh portable reproduction of exact recorded settings; paths and executable lookup adapted','sha256':{str(f.resolve()):sha(f) for f in locked}}
 write(base/'parallel_call/FROZEN_SOURCES.json',guard)
 missing_tools=[t for t in cfg['tools'] if not shutil.which(cfg['tools'][t])]
 print(json.dumps({'configuration_written':str(cfgpath),'reference_checks':'SKIPPED' if a.skip_reference_bytecheck else 'PASS','tools_not_on_PATH':missing_tools,'analysis_run':False,'jobs_submitted':False},indent=2))
if __name__=='__main__':main()
