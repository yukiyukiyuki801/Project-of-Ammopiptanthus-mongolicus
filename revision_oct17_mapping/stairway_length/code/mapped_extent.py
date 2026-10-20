#!/usr/bin/env python3
"""Observed mapped SNP-region extent, not historical Stairway L or callable length."""
import argparse,bisect,csv,hashlib,heapq,json,mmap,os,re,shutil,subprocess,sys,tempfile,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3];BASE=ROOT/'revision_oct17_mapping/stairway_length'
CIGAR=re.compile(r'(\d+)([MIDNSHP=X])')
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def dump(p,d):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_name(p.name+'.partial');tmp.write_text(json.dumps(d,indent=2)+'\n');tmp.replace(p)
def blocks(start,cigar):
 """start is zero based. Insertions/clips consume no reference; D/N split."""
 parts=CIGAR.findall(cigar);assert parts and ''.join(a+b for a,b in parts)==cigar,('Invalid CIGAR',cigar)
 pos=start;left=None
 for text,op in parts:
  n=int(text);assert n>0
  if op in 'M=X':
   if left is None:left=pos
   pos+=n
  elif op in 'DN':
   if left is not None:yield left,pos;left=None
   pos+=n
  elif op in 'ISHP':pass
 if left is not None:yield left,pos

def merge_intervals(rows):
 current=None;previous=None
 for chrom,start,end in rows:
  assert 0<=start<end
  if previous is not None:assert (chrom,start,end)>=previous,'Unsorted intervals'
  previous=(chrom,start,end)
  if current and chrom==current[0] and start<=current[2]:
   assert start>=current[1],'Unsorted intervals';current=(chrom,current[1],max(end,current[2]))
  else:
   if current:
    assert (chrom,start)>=current[:2],'Unsorted intervals';yield current
   current=(chrom,start,end)
 if current:yield current

def bed_rows(path):
 with open(path) as f:
  for line in f:
   x=line.rstrip().split('\t');yield x[0],int(x[1]),int(x[2])
class Fasta:
 def __init__(self,p):
  self.p=Path(p);self.handle=self.p.open('rb');self.data=mmap.mmap(self.handle.fileno(),0,access=mmap.ACCESS_READ);self.index={}
  for l in Path(str(p)+'.fai').read_text().splitlines():
   c,n,o,b,w=l.split('\t')[:5];self.index[c]=tuple(map(int,(n,o,b,w)))
 def sequence(self,c,s,e):
  n,o,b,w=self.index[c];assert 0<=s<=e<=n
  if s==e:return b''
  a=o+(s//b)*w+s%b;z=o+((e-1)//b)*w+(e-1)%b+1
  x=self.data[a:z].replace(b'\n',b'').replace(b'\r',b'').upper();assert len(x)==e-s;return x
 def acgt(self,c,s,e):
  count=0
  for a in range(s,e,1024*1024):
   x=self.sequence(c,a,min(e,a+1024*1024));count+=sum(x.count(k) for k in (b'A',b'C',b'G',b'T'))
  return count
 def close(self):self.data.close();self.handle.close()
def select_regions(rows,targets):
 """Keep entire covered component if it contains a target. Never bridge gaps."""
 for c,s,e in rows:
  positions=targets.get(c,[]);i=bisect.bisect_left(positions,s)
  if i<len(positions) and positions[i]<e:yield c,s,e

def config(path):
 cfg=json.loads(Path(path).read_text());assert cfg['measurement']=='observed_mapped_positive_coverage_extent';assert cfg['reference']=='nanus';assert cfg['min_mapq']==20 and cfg['exclude_flags']==0x904
 out=Path(cfg['output_root']);assert out.is_absolute() and '/scratch/' in str(out),'Heavy products must use scratch'
 workflow=json.loads(Path(cfg['workflow_config']).read_text());assert workflow['execution_frozen']
 return cfg,workflow,out

def fingerprint(cfgpath,cfg,w,stage):
 return {'stage':stage,'source_sha256':sha(__file__),'config_sha256':sha(cfgpath),'workflow_config_sha256':sha(cfg['workflow_config']),'manifest_sha256':sha(w['manifest']),'registry_sha256':sha(cfg['registry'])}
def valid_saved(marker,fp):
 if not marker.exists():return False
 d=json.loads(marker.read_text());assert d['validated'] and d['fingerprint']==fp,'Changed checkpoint inputs; use a new output root'
 for p,h in d['outputs'].items():assert Path(p).is_file() and sha(p)==h,('Changed product',p)
 return True

def sample_stage(cfgpath,index):
 cfg,w,out=config(cfgpath);rows=list(csv.DictReader(open(w['manifest']),delimiter='\t'));assert len(rows)==217 and len({r['sample_id'] for r in rows})==217
 assert 0<=index<len(rows);sample=rows[index]['sample_id'];fp=fingerprint(cfgpath,cfg,w,'sample_'+sample)
 maproot=Path(w['output_root']);checkpoint=maproot/'checkpoints'/f'map_nanus_{sample}.json';d=json.loads(checkpoint.read_text());bam=maproot/'mapping/nanus'/sample/'aligned.bam'
 assert d['validated'] and d['stage']==f'map_nanus_{sample}'
 assert d['fingerprint']['config']==sha(cfg['workflow_config']) and d['fingerprint']['manifest']==sha(w['manifest'])
 assert d['fingerprint']['script']==sha(ROOT/'revision_oct17_mapping/workflow/run_workflow.py'),'Mapping source changed'
 expected=d['outputs'].get(str(bam.resolve()));assert expected,'BAM absent from validated mapping checkpoint'
 fp.update(mapping_checkpoint_sha256=sha(checkpoint),bam_sha256=expected)
 assert sha(bam)==expected,'Input BAM checksum mismatch'
 dest=out/'samples';dest.mkdir(parents=True,exist_ok=True);marker=dest/(sample+'.json')
 if valid_saved(marker,fp):print('VALIDATED '+sample);return
 bed=dest/(sample+'.bed');assert not bed.exists(),'Uncheckpointed BED exists; inspect before restart'
 lock=dest/(sample+'.lock');fd=os.open(lock,os.O_CREAT|os.O_EXCL|os.O_WRONLY);os.close(fd)
 tmp=Path(tempfile.mkdtemp(prefix=sample+'.partial.',dir=out));procs=[]
 try:
  tool=w['tools']['samtools'];subprocess.run([tool,'quickcheck','-v',str(bam)],check=True)
  sortedfile=tmp/'blocks.sorted.tsv';env={**os.environ,'LC_ALL':'C'}
  with sortedfile.open('w') as sortedout,(tmp/'samtools.stderr').open('w') as err:
   sorter=subprocess.Popen(['/usr/bin/sort','-T',str(tmp),'-S','512M','--parallel=2','-k1,1','-k2,2n','-k3,3n'],stdin=subprocess.PIPE,stdout=sortedout,text=True,env=env);procs.append(sorter)
   sam=subprocess.Popen([tool,'view','-q','20','-F','2308',str(bam)],stdout=subprocess.PIPE,stderr=err,text=True);procs.append(sam)
   reads=0;parts=0
   for line in sam.stdout:
    fields=line.split('\t',6);assert len(fields)>=6;chrom=fields[2];start=int(fields[3])-1
    for s,e in blocks(start,fields[5]):sorter.stdin.write(f'{chrom}\t{s}\t{e}\n');parts+=1
    reads+=1
   sam.stdout.close();assert sam.wait()==0,'samtools view failed';sorter.stdin.close();assert sorter.wait()==0,'sort failed'
  merged=tmp/'merged.bed';intervals=0;bases=0
  with merged.open('w') as f:
   for c,s,e in merge_intervals(bed_rows(sortedfile)):f.write(f'{c}\t{s}\t{e}\n');intervals+=1;bases+=e-s
  assert reads>0 and intervals>0,'Empty mapping extent';merged.replace(bed)
  dump(marker,{'validated':True,'fingerprint':fp,'outputs':{str(bed):sha(bed)},'sample_id':sample,'eligible_primary_alignments':reads,'aligned_blocks_before_union':parts,'union_intervals':intervals,'observed_reference_bases_including_ambiguous':bases,'completed_unix':time.time(),'interpretation':'Per-sample positive mapped coverage; no genotype or base-quality callability threshold.'})
  shutil.rmtree(tmp)
 except BaseException as exc:
  for p in procs:
   if p.poll() is None:p.terminate()
  dump(dest/(sample+'.error.json'),{'error':repr(exc),'partial_directory':str(tmp),'time':time.time()});raise
 finally:lock.unlink(missing_ok=True)

def aggregate(cfgpath):
 cfg,w,out=config(cfgpath);registry=json.loads(Path(cfg['registry']).read_text());assert registry['status']=='VERIFIED_MEMBERSHIP'
 fp=fingerprint(cfgpath,cfg,w,'aggregate');markers={};beds={};rows=list(csv.DictReader(open(w['manifest']),delimiter='\t'));ids={r['sample_id'] for r in rows};assert len(ids)==217
 for s in ids:
  p=out/'samples'/(s+'.json');d=json.loads(p.read_text());assert d['validated'] and d['sample_id']==s
  expected=fingerprint(cfgpath,cfg,w,'sample_'+s)
  assert all(d['fingerprint'][k]==v for k,v in expected.items()),('Sample source differs',s)
  for f,h in d['outputs'].items():assert sha(f)==h;beds[s]=Path(f)
  markers[s]=sha(p)
 fp['sample_checkpoint_hashes']=markers
 ref=Path(w['output_root'])/'references/nanus/genome.fa';refmark=Path(w['output_root'])/'checkpoints/reference_nanus.json';rd=json.loads(refmark.read_text());assert rd['validated'] and rd['fingerprint']['config']==sha(cfg['workflow_config'])
 assert rd['fingerprint']['script']==sha(ROOT/'revision_oct17_mapping/workflow/run_workflow.py'),'Reference source changed'
 for p in (ref,Path(str(ref)+'.fai')):assert sha(p)==rd['outputs'][str(p.resolve())]
 fp['reference_checkpoint_sha256']=sha(refmark);fp['all_target_reference_bases_validated']=True;dest=out/'aggregate';marker=dest/'MEASUREMENT.json'
 if valid_saved(marker,fp):print('VALIDATED aggregate');return
 assert not dest.exists(),'Uncheckpointed aggregate exists; inspect before restart';lock=out/'aggregate.lock';fd=os.open(lock,os.O_CREAT|os.O_EXCL|os.O_WRONLY);os.close(fd)
 tmp=Path(tempfile.mkdtemp(prefix='aggregate.partial.',dir=out));fa=Fasta(ref)
 try:
  groupbeds={};outputs=[]
  for group,samples in registry['groups'].items():
   assert len(samples)==len(set(samples)) and set(samples)<=ids
   p=tmp/(group+'.union.bed')
   with p.open('w') as f:
    streams=[bed_rows(beds[s]) for s in samples]
    for c,a,b in merge_intervals(heapq.merge(*streams)):f.write(f'{c}\t{a}\t{b}\n')
   groupbeds[group]=p;outputs.append(p)
  results=[]
  for t in registry['targets']:
   assert t['samples']==registry['groups'][t['membership_group']]
   targetpath=Path(cfg['registry']).parent/t['target_bed'];allelepath=Path(cfg['registry']).parent/t['alleles_tsv']
   assert sha(targetpath)==t['target_bed_sha256'] and sha(allelepath)==t['alleles_sha256']
   targets={};refsites={}
   for c,pos,refbase,alt in (l.split('\t') for l in allelepath.read_text().splitlines()):
    pos=int(pos);assert fa.sequence(c,pos-1,pos).decode()==refbase,('Target REF mismatch',c,pos);assert len(refbase)==1 and len(alt)==1
    refsites[(c,pos-1)]=refbase
   for c,s,e in bed_rows(targetpath):
    assert e==s+1 and (c,s) in refsites;targets.setdefault(c,[]).append(s)
   for c in targets:targets[c]=sorted(set(targets[c]))
   assert sum(map(len,targets.values()))==t['positive_SFS_contributor_sites']
   name=targetpath.stem;path=tmp/(name+'.mapped_regions.bed');missing=tmp/(name+'.uncovered_targets.bed');covered=set();intervals=0;total=0;acgt=0
   with path.open('w') as f:
    for c,a,b in select_regions(bed_rows(groupbeds[t['membership_group']]),targets):
     f.write(f'{c}\t{a}\t{b}\n');intervals+=1;total+=b-a;acgt+=fa.acgt(c,a,b)
     positions=targets[c];i=bisect.bisect_left(positions,a);j=bisect.bisect_left(positions,b)
     covered.update((c,p) for p in positions[i:j])
   assert intervals>0 and covered,'No historical target SNP observed; inspect inputs'
   missing.write_text(''.join(f'{c}\t{p}\t{p+1}\n' for c,poss in sorted(targets.items()) for p in poss if (c,p) not in covered));outputs.extend([path,missing])
   results.append({'target_bed':t['target_bed'],'membership_group':t['membership_group'],'samples':len(t['samples']),'projection_chromosomes':t['projection_chromosomes'],'target_snps':len(refsites),'covered_target_snps':len(covered),'uncovered_target_snps':len(refsites)-len(covered),'mapped_regions':intervals,'observed_reference_extent_bp':total,'observed_ACGT_extent_bp':acgt,'ambiguous_reference_bp_excluded':total-acgt,'historical_L_replaced':False,'callable_denominator_verified':False})
  fp['all_target_reference_bases_validated']=True
  (tmp/'RESULTS.json').write_text(json.dumps(results,indent=2)+'\n');outputs.append(tmp/'RESULTS.json')
  with (tmp/'RESULTS.csv').open('w') as f:
   writer=csv.DictWriter(f,fieldnames=list(results[0]));writer.writeheader();writer.writerows(results)
  outputs.append(tmp/'RESULTS.csv');finalhash={str(dest/p.name):sha(p) for p in outputs};tmp.replace(dest)
  dump(marker,{'validated':True,'fingerprint':fp,'outputs':finalhash,'results':results,'interpretation':'New observed mapped extent of regions containing historical SFS-positive SNPs. Not historic input L, not invariant-inclusive genotype-callable sequence. No demographic rerun or recalibration.','completed_unix':time.time()})
 except BaseException as exc:dump(out/'aggregate.error.json',{'error':repr(exc),'partial_directory':str(tmp)});raise
 finally:fa.close();lock.unlink(missing_ok=True)
def main():
 p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--stage',choices=['sample','aggregate'],required=True);p.add_argument('--index',type=int);a=p.parse_args()
 if a.stage=='sample':sample_stage(a.config,a.index if a.index is not None else int(os.environ['SLURM_ARRAY_TASK_ID']))
 else:aggregate(a.config)
if __name__=='__main__':main()
