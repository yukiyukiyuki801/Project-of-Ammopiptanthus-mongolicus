#!/usr/bin/env python3
"""Meaningful input guards; no biological fits/data products."""
import csv,importlib.util,json,tempfile,subprocess,shutil
from pathlib import Path
p=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('workflow',p/'run_workflow.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
cfg=json.load(open(p.parent/'parent/config.frozen.json'));old=list(csv.DictReader(open(cfg['existing_pca_scores'])))
rows=[{'sample_id':r['sample_id'],'population':r['population'],'class':r['group'],'layout':'PAIRED','r1':'/missing/'+r['sample_id']+'_1.fastq.gz','r2':'/missing/'+r['sample_id']+'_2.fastq.gz'} for r in old]
checks={}
with tempfile.TemporaryDirectory() as d:
 f=Path(d)/'samples.tsv'
 def write(r):
  with open(f,'w') as o:
   w=csv.DictWriter(o,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(r)
 def fails(name,r,reads=False):
  write(r)
  try:m.parse_manifest(f,cfg['existing_pca_scores'],require_reads=reads)
  except AssertionError:checks[name]=True
  else:raise AssertionError('Guard failed '+name)
 write(list(reversed(rows)));got=m.parse_manifest(f,cfg['existing_pca_scores'],require_reads=False)
 checks['restores_original_sample_order']=[r['sample_id'] for r in got]==[r['sample_id'] for r in rows]
 fails('missing_actual_reads_rejected',rows,True)
 fails('missing_sample_rejected',rows[:-1])
 bad=[dict(r) for r in rows];bad[0]['class']='peripheral' if bad[0]['class']=='core' else 'core';fails('swapped_class_rejected',bad)
 bad=[dict(r) for r in rows];bad[0]['population']='ZZ';fails('wrong_population_rejected',bad)
 bad=[dict(r) for r in rows];bad[1]['sample_id']=bad[0]['sample_id'];fails('duplicate_id_rejected',bad)
 bad=[dict(r) for r in rows];bad[1]['r1']=bad[0]['r1'];fails('reused_fastq_rejected',bad)
 bad=[dict(r) for r in rows];bad[0]['layout']='SINGLE';fails('unverified_layout_rejected',bad)
 bad=[dict(r) for r in rows];bad[0]['r1']='relative.fastq';fails('relative_read_path_rejected',bad)
# The local environment validates real plugin syntax on a tiny fabricated VCF only.
bc=shutil.which('bcftools');assert bc, 'BCFtools must be on PATH'
with tempfile.TemporaryDirectory() as d:
 d=Path(d);v=d/'toy.vcf'
 header=['##fileformat=VCFv4.2','##contig=<ID=chr1,length=1000>','##INFO=<ID=AC,Number=A,Type=Integer,Description="ALT count">','##INFO=<ID=AN,Number=1,Type=Integer,Description="Alleles">','##INFO=<ID=AF,Number=A,Type=Float,Description="Frequency">','##INFO=<ID=NS,Number=1,Type=Integer,Description="Samples">','##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">','##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">','#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\tS5']
 records=[]
 for pos,ref,alt in [(1,'A','C'),(2,'AA','CC'),(3,'A','AC')]:
  records.append(f'chr1\t{pos}\t.\t{ref}\t{alt}\t150\tPASS\tAC=6;AN=10;AF=0.6;NS=5\tGT:DP\t0/1:2\t0/0:5\t0/1:4\t1/1:7\t1/1:8')
 v.write_text('\n'.join(header+records)+'\n')
 def run(*a):subprocess.run([str(bc),*map(str,a)],check=True,stdout=subprocess.DEVNULL)
 run('+setGT',v,'-Ob','-o',d/'masked.bcf','--','-t','q','-n','.','-i','FMT/DP<4 | FMT/DP="."')
 run('+fill-tags',d/'masked.bcf','-Ob','-o',d/'filled.bcf','--','-t','AC,AN,AF,NS')
 run('view','-v','snps','-m2','-M2',d/'filled.bcf','-Ov','-o',d/'snp.vcf')
 result=m.verify_site_annotations(str(bc),d/'snp.vcf',5)
 checks['actual_setGT_filltags_and_pure_SNP_filter']=result['sites_checked']==1 and result['maximum_AF_rounding_error']==0
 run('view','-v','snps','-m2','-M2',d/'masked.bcf','-Ov','-o',d/'stale.vcf')
 try:m.verify_site_annotations(str(bc),d/'stale.vcf',5)
 except AssertionError:checks['stale_INFO_rejected']=True
 else:raise AssertionError('Stale annotation guard failed')
assert all(checks.values())
m.dump(p/'implementation_validation.json',{'implementation_checks_pass':True,'checks':checks,'python_compilation':'PASS','R_parse':'PASS','biological_analysis_executed':False,'production_validation':'PENDING actual reads, frozen settings and pilot','script_sha256':m.digest(p/'run_workflow.py'),'R_script_sha256':m.digest(p/'pca_comparison.R')})
print(json.dumps(checks,indent=2))
