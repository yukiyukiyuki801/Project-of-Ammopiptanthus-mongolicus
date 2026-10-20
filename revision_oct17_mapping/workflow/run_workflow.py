#!/usr/bin/env python3
"""Parent-submitted Slurm worker; no submission and no automatic stage chaining."""
import argparse,csv,gzip,hashlib,json,os,re,shutil,subprocess,sys,tempfile,time
from pathlib import Path
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
def digest(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''): h.update(b)
 return h.hexdigest()
def dump(p,x):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
 t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(x,indent=2)+'\n');t.replace(p)
def command(args,stdout=None):
 print(json.dumps([str(x) for x in args]),flush=True)
 subprocess.run([str(x) for x in args],check=True,stdout=stdout)
def count(tool,args): return int(subprocess.check_output([tool,'view','-c',*map(str,args)],text=True).strip())
def parse_manifest(path,old,require_reads=True):
 rows=list(csv.DictReader(open(path),delimiter='\t'))
 assert len(rows)==217,'Expected exactly 217 samples, combine lanes first; do not drop samples'
 assert {'sample_id','population','class','layout','r1','r2'} <= set(rows[0])
 expected={r['sample_id']:r for r in csv.DictReader(open(old))}
 assert len(expected)==217 and set(expected)=={r['sample_id'] for r in rows},'Sample identity differs from existing PCA'
 assert len({r['sample_id'] for r in rows})==217
 paths=[]
 for r in rows:
  assert re.fullmatch(r'[A-Za-z0-9_.-]+',r['sample_id']), 'Unsafe sample ID'
  assert r['layout']=='PAIRED', 'Only verified paired-end reads implemented'
  assert r['population']==expected[r['sample_id']]['population']
  assert r['class']==expected[r['sample_id']]['group']
  assert r['r1']!=r['r2']
  for key in ('r1','r2'):
   p=Path(r[key]);assert p.is_absolute(), 'FASTQ paths must be absolute'
   if require_reads: assert p.is_file() and p.stat().st_size>0,str(p)
   paths.append(str(p))
 assert len(paths)==len(set(paths)),'FASTQ path reused across samples or mates'
 # Freeze original PCA sample order, independently of manifest row order.
 order={s:i for i,s in enumerate(expected)}
 return sorted(rows,key=lambda r:order[r['sample_id']])
def verify_site_annotations(bcftools,path,nsamples):
 """Check refreshed biallelic INFO against retained GT, not caller metadata."""
 p=subprocess.Popen([bcftools,'query','-f','%REF\t%ALT\t%INFO/AC\t%INFO/AN\t%INFO/AF\t%INFO/NS[\t%GT:%DP]\n',str(path)],stdout=subprocess.PIPE,text=True)
 n=0;max_error=0
 try:
  for line in p.stdout:
   fields=line.rstrip('\n').split('\t');ref,alt,ac,an,af,ns=fields[:6]
   assert len(ref)==1 and len(alt)==1 and ref!=alt,'Non-SNP retained'
   assert len(fields[6:])==nsamples,'Wrong number of genotypes'
   gt_ac=0;gt_an=0;gt_ns=0
   for sample in fields[6:]:
    gt,dp=sample.split(':');alleles=re.split(r'[/|]',gt)
    assert len(alleles)==2,'Non-diploid GT'
    called=[a for a in alleles if a!='.'];assert all(a in ('0','1') for a in called)
    if called:
     assert len(called)==2 and dp!='.' and int(dp)>=4,'Unmasked low/missing-depth GT'
     gt_ac+=sum(int(a) for a in called);gt_an+=len(called);gt_ns+=1
   assert int(ac)==gt_ac and int(an)==gt_an and int(ns)==gt_ns,'Stale INFO allele/sample counts'
   assert gt_an>0 and gt_ns/nsamples>=0.8
   err=abs(float(af)-gt_ac/gt_an);assert err<1e-5,'Stale INFO AF';max_error=max(max_error,err);n+=1
  assert p.wait()==0 and n>0,'Annotation verification failed or no sites retained'
 except BaseException:
  p.kill();p.wait();raise
 return {'validated':True,'sites_checked':n,'pure_biallelic_SNPs':True,'INFO_recomputed_from_masked_genotypes':True,'maximum_AF_rounding_error':max_error,'DP_mask_and_completeness_checked':True}
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--config',required=True);ap.add_argument('--stage',required=True,choices=['preflight','reference','clean','map','call','filter','pca','summarize']);ap.add_argument('--sample');ap.add_argument('--reference',choices=['nanus','mongolicus']);a=ap.parse_args()
 cfgpath=Path(a.config).resolve();cfg=json.loads(cfgpath.read_text());out=Path(cfg['output_root']);out.mkdir(parents=True,exist_ok=True)
 rows=[] if a.stage=='reference' else parse_manifest(cfg['manifest'],cfg['existing_pca_scores']);tools=cfg['tools'];threads=int(cfg['threads']);assert threads>0
 refaudit=json.loads(Path(cfg['reference_audit']).read_text()); refs=cfg['references']
 for arm in refs:
  match=[x for x in refaudit['references'] if x['sha256']==refs[arm]['compressed_sha256'] and Path(x['path']).resolve()==Path(refs[arm]['source_gz']).resolve()]
  assert len(match)==1,'Configured reference must match audited exact source'
 selected=[r for r in rows if r['sample_id']==a.sample] if a.sample else []
 if a.stage in ('clean','map'): assert len(selected)==1,'Supply valid --sample'
 if a.stage in ('map','reference'): assert a.reference,'Supply --reference'
 needed={'preflight':list(tools),'reference':['bwa','samtools'],'clean':['fastp'],'map':['bwa','samtools'],'call':['bcftools'],'filter':['bcftools','vcftools','plink'],'pca':['Rscript'],'summarize':[]}[a.stage]
 for name in needed: assert shutil.which(tools[name]),f'Missing executable: {name} {tools[name]}'
 if a.stage=='preflight':
  versions={}
  for name in needed:
   p=subprocess.run([tools[name],'--version'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
   versions[name]=p.stdout[:5000]
  dump(out/'preflight.json',{'validated_inputs':True,'analysis_executed':False,'sample_count':len(rows),'class_counts':{g:sum(r['class']==g for r in rows) for g in ('core','peripheral')},'versions_raw':versions,'configuration_sha256':digest(cfgpath)})
  return
 assert cfg.get('execution_frozen') is True,'Parent must freeze explicit cleaning/caller settings and input manifest before execution'
 assert cfg['caller']=='bcftools_mpileup_call','Only explicit bcftools sensitivity implemented; not the published two-caller intersection'
 scratch=Path(cfg['scratch_root']);scratch.mkdir(parents=True,exist_ok=True)
 key=a.stage+('_'+a.reference if a.reference else '')+('_'+a.sample if a.sample else '')
 marker=out/'checkpoints'/f'{key}.json'
 # Script/config/manifest changes invalidate checkpoint; referenced products also checksum-verified.
 fingerprint={'script':digest(__file__),'config':digest(cfgpath),'manifest':digest(cfg['manifest']) if a.stage!='reference' else 'NOT_REQUIRED_FOR_INDEX','stage':key,'pca_script':digest(HERE/'pca_comparison.R'),'existing_pca_scores':digest(cfg['existing_pca_scores']),'existing_pca_variance':digest(cfg['existing_pca_variance']),'reference_audit':digest(cfg['reference_audit'])}
 if marker.exists():
  old=json.loads(marker.read_text())
  assert old['fingerprint']==fingerprint,'Checkpoint inputs changed; use a new output root instead of overwrite'
  assert all(Path(p).is_file() and digest(p)==h for p,h in old['outputs'].items()),'Checkpoint product missing/changed'
  print('VALIDATED CHECKPOINT: '+key);return
 outputs=[]
 def must(stage):
  p=out/'checkpoints'/f'{stage}.json';assert p.exists(),f'Missing dependency {p}'
  d=json.loads(p.read_text());assert d['validated'] and d['fingerprint']['config']==fingerprint['config'] and (stage.startswith('reference_') or d['fingerprint']['manifest']==fingerprint['manifest'])
  assert all(d['fingerprint'].get(k)==fingerprint[k] for k in ('script','pca_script','existing_pca_scores','existing_pca_variance','reference_audit')),'Dependency source fingerprint changed'
  assert all(Path(q).is_file() and digest(q)==h for q,h in d['outputs'].items()),'Dependency product changed'
 def bc(*args):command([tools['bcftools'],*args])
 with tempfile.TemporaryDirectory(prefix=key+'_',dir=scratch) as tmp:
  tmp=Path(tmp)
  if a.stage=='reference':
   rr=refs[a.reference];src=Path(rr['source_gz']);assert digest(src)==rr['compressed_sha256']
   dest=out/'references'/a.reference;dest.mkdir(parents=True,exist_ok=True);fa=dest/'genome.fa'
   assert not fa.exists(),'Uncheckpointed reference exists: inspect before restarting'
   with gzip.open(src,'rb') as f,open(fa,'wb') as g:shutil.copyfileobj(f,g)
   command([tools['samtools'],'faidx',fa]);command([tools['bwa'],'index',fa])
   outputs=[fa,Path(str(fa)+'.fai'),*[Path(str(fa)+e) for e in ('.amb','.ann','.bwt','.pac','.sa')]]
  elif a.stage=='clean':
   r=selected[0];dest=out/'clean'/a.sample;dest.mkdir(parents=True,exist_ok=True)
   assert not (dest/'R1.fastq.gz').exists(),'Uncheckpointed cleaning exists; inspect before restart'
   extra=cfg['fastp_arguments'];assert isinstance(extra,list) and all(isinstance(x,str) for x in extra)
   command([tools['fastp'],'--in1',r['r1'],'--in2',r['r2'],'--out1',dest/'R1.fastq.gz','--out2',dest/'R2.fastq.gz','--thread',min(16,threads),'--json',dest/'fastp.json','--html',dest/'fastp.html',*extra])
   report=json.loads((dest/'fastp.json').read_text());assert report['summary']['after_filtering']['total_reads']>0
   assert report['summary']['after_filtering']['total_reads']%2==0
   info={'source_fastq_sha256':{k:digest(r[k]) for k in ('r1','r2')},'clean_pairs':report['summary']['after_filtering']['total_reads']//2}
   dump(dest/'source_identity.json',info);outputs=list(dest.iterdir())
  elif a.stage=='map':
   must('clean_'+a.sample);must('reference_'+a.reference)
   fa=out/'references'/a.reference/'genome.fa';dest=out/'mapping'/a.reference/a.sample;dest.mkdir(parents=True,exist_ok=True);bam=dest/'aligned.bam'
   assert not bam.exists(),'Uncheckpointed BAM exists; inspect before restart'
   rg='@RG\\tID:'+a.sample+'\\tSM:'+a.sample+'\\tPL:ILLUMINA'
   cmd=[tools['bwa'],'mem','-t',str(threads),'-R',rg,str(fa),str(out/'clean'/a.sample/'R1.fastq.gz'),str(out/'clean'/a.sample/'R2.fastq.gz')]
   print(json.dumps(cmd),flush=True)
   with open(dest/'bwa.stderr.txt','w') as log:
    p=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=log)
    try:
     q=subprocess.run([tools['samtools'],'sort','-@',str(max(1,threads//2)),'-m',cfg['sort_memory_per_thread'],'-T',str(tmp/'sort'),'-o',str(bam),'-'],stdin=p.stdout)
     p.stdout.close();rc=p.wait();assert q.returncode==0 and rc==0,'BWA/samtools pipe failed'
    except BaseException:p.kill();p.wait();raise
   command([tools['samtools'],'quickcheck','-v',bam]);command([tools['samtools'],'index','-@',threads,bam])
   total=count(tools['samtools'],['-F',2304,bam]);mapped=count(tools['samtools'],['-F',2308,bam]);proper=count(tools['samtools'],['-f',2,'-F',2308,bam]);mq=count(tools['samtools'],['-q',20,'-F',2308,bam])
   clean=json.loads((out/'clean'/a.sample/'source_identity.json').read_text());assert total==clean['clean_pairs']*2,'Primary alignment denominator differs from cleaned read count'
   # Stream only covered bases, never materialize a genome-sized depth table.
   proc=subprocess.Popen([tools['samtools'],'depth','-q','20','-Q','20','-g','0x600','-G','0x904',str(bam)],stdout=subprocess.PIPE,text=True)
   hist={};nbase=0;sumdepth=0
   for line in proc.stdout:
    depth=int(line.rsplit('\t',1)[1])
    if depth==0:continue
    hist[depth]=hist.get(depth,0)+1;nbase+=1;sumdepth+=depth
   assert proc.wait()==0 and nbase>0
   def quantile(frac):
    threshold=(nbase-1)*frac;cum=0
    for d,n in sorted(hist.items()):
     cum+=n
     if cum>threshold:return d
   genome=sum(int(line.split('\t')[1]) for line in open(str(fa)+'.fai'));covered4=sum(n for d,n in hist.items() if d>=4)
   metrics={'sample_id':a.sample,'population':selected[0]['population'],'class':selected[0]['class'],'reference':a.reference,'clean_pairs':clean['clean_pairs'],'primary_records':total,'mapped_primary':mapped,'mapped_fraction':mapped/total,'properly_paired_primary_fraction':proper/total,'mapq20_primary_fraction':mq/total,'covered_bases_q20Q20':nbase,'covered_bases_ge4':covered4,'covered_base_mean_depth':sumdepth/nbase,'covered_base_median_depth_lower':quantile(.5),'assembly_bases':genome,'assembly_breadth_ge1':nbase/genome,'assembly_breadth_ge4':covered4/genome,'whole_assembly_mean_depth':sumdepth/genome,'depth_denominator_warning':'Observed covered bases differ between assemblies; not homologous-locus callability or genotype missingness'}
   dump(dest/'metrics.json',metrics);dump(dest/'depth_histogram.json',hist);outputs=[bam,Path(str(bam)+'.bai'),dest/'metrics.json',dest/'depth_histogram.json',dest/'bwa.stderr.txt']
  elif a.stage=='call':
   must('reference_mongolicus')
   for r in rows:must('map_mongolicus_'+r['sample_id'])
   dest=out/'variants';dest.mkdir(parents=True,exist_ok=True);bams=dest/'bam.list';bams.write_text(''.join(str(out/'mapping'/'mongolicus'/r['sample_id']/'aligned.bam')+'\n' for r in rows));fa=out/'references'/'mongolicus'/'genome.fa'
   raw=dest/'raw.vcf.gz';assert not raw.exists(),'Uncheckpointed calls exist'
   # Explicit per-sample depth cap, MAPQ/base-quality filters, and diploid multiallelic caller.
   pile=tmp/'pileup.bcf';bc('mpileup','--threads',threads,'-f',fa,'-b',bams,'-q','20','-Q','20','-d',cfg['pileup_max_depth_per_sample'],'-a','FORMAT/DP,FORMAT/AD','-Ob','-o',pile)
   bc('call','--threads',threads,'-m','-v','--ploidy','2','-Oz','-o',raw,pile);bc('index','-t',raw)
   outputs=[raw,Path(str(raw)+'.tbi'),bams]
  elif a.stage=='filter':
   must('call');dest=out/'variants';raw=dest/'raw.vcf.gz';masked=dest/'dp_masked.vcf.gz'
   masked_stale=tmp/'masked_stale.bcf'
   bc('+setGT',raw,'-Ob','-o',masked_stale,'--','-t','q','-n','.','-i','FMT/DP<4 | FMT/DP="."')
   bc('+fill-tags',masked_stale,'-Oz','-o',masked,'--','-t','AC,AN,AF,NS');bc('index','-t',masked)
   prefix=dest/'site_filtered';command([tools['vcftools'],'--gzvcf',masked,'--min-alleles','2','--max-alleles','2','--remove-indels','--max-missing','0.8','--mac','3','--maf','0.01','--minQ','100','--recode','--recode-INFO-all','--out',prefix])
   sites=dest/'site_filtered.vcf.gz';bc('view','-v','snps','-m2','-M2','-Oz','-o',sites,str(prefix)+'.recode.vcf');bc('index','-t',sites)
   verify=verify_site_annotations(tools['bcftools'],sites,len(rows));dump(dest/'site_annotation_validation.json',verify)
   ids=dest/'sites_ids.vcf.gz';bc('annotate','--set-id','%CHROM:%POS:%REF:%FIRST_ALT','-Oz','-o',ids,sites);bc('index','-t',ids)
   prune=dest/'ld';command([tools['plink'],'--vcf',ids,'--double-id','--allow-extra-chr','--set-missing-var-ids','@:#','--indep-pairwise','50','10','0.8','--out',prune])
   retained=dest/'pca.vcf.gz';bc('view','-i',f'ID=@{prune}.prune.in','-Oz','-o',retained,ids);bc('index','-t',retained)
   # Explicit source counts and filtered genotype missingness, not mapping missingness.
   counts={}
   for name,v in [('raw',raw),('site_filtered',sites),('ld_pruned',retained)]:
    result=subprocess.run([tools['bcftools'],'index','-n',str(v)],capture_output=True,text=True,check=True);counts[name]=int(result.stdout)
   assert counts['ld_pruned']>=2,'Insufficient retained loci for PCA'
   with open(dest/'dosage.tsv','w') as f:bcquery=[tools['bcftools'],'query','-f','%CHROM\t%POS\t%ID[\t%GT]\n',retained];command(bcquery,stdout=f)
   with open(dest/'sample_order.txt','w') as f:command([tools['bcftools'],'query','-l',retained],stdout=f)
   dump(dest/'filter_counts.json',counts);outputs=[masked,Path(str(masked)+'.tbi'),sites,Path(str(sites)+'.tbi'),ids,Path(str(ids)+'.tbi'),retained,Path(str(retained)+'.tbi'),Path(str(prune)+'.prune.in'),dest/'dosage.tsv',dest/'sample_order.txt',dest/'filter_counts.json',Path(str(prefix)+'.log'),dest/'site_annotation_validation.json']
  elif a.stage=='pca':
   must('filter');dest=out/'pca';dest.mkdir(parents=True,exist_ok=True)
   command([tools['Rscript'],HERE/'pca_comparison.R',cfgpath])
   outputs=list(dest.iterdir());assert (dest/'validation.json').exists()
  elif a.stage=='summarize':
   metrics=[]
   for arm in refs:
    for r in rows:
     must('map_'+arm+'_'+r['sample_id']);metrics.append(json.loads((out/'mapping'/arm/r['sample_id']/'metrics.json').read_text()))
   dest=out/'summary';dest.mkdir(parents=True,exist_ok=True);p=dest/'paired_mapping_depth.csv'
   with open(p,'w') as f:
    w=csv.DictWriter(f,fieldnames=list(metrics[0]));w.writeheader();w.writerows(metrics)
   outputs=[p]
 assert outputs and all(p.is_file() and p.stat().st_size>0 for p in outputs),'Empty/missing stage output'
 dump(marker,{'validated':True,'stage':key,'fingerprint':fingerprint,'outputs':{str(p.resolve()):digest(p) for p in outputs},'completed_unix':time.time(),'scientific_interpretation_validated':False})
if __name__=='__main__':main()
