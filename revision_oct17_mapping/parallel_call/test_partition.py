#!/usr/bin/env python3
"""Tiny synthetic calling equivalence; no study data, jobs or source mutation."""
import importlib.util,pathlib,tempfile,subprocess,json,shutil
base=pathlib.Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('shards',base/'call_shards.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
bc=shutil.which('bcftools');sam=shutil.which('samtools');assert bc and sam, 'BCFtools and SAMtools must be on PATH'
def run(*args):return subprocess.check_output(list(map(str,args)),stderr=subprocess.PIPE,text=True)
with tempfile.TemporaryDirectory(prefix='parallel_test_') as tmp:
 p=pathlib.Path(tmp);ref=p/'ref.fa';seq='ACGT'*20
 ref.write_text(''.join('>'+c+'\n'+seq+'\n' for c in ('chr1','scaffold1','empty')));run(sam,'faidx',ref)
 bams=[]
 for sample in ('a','b','c'):
  rows=['@HD\tVN:1.6\tSO:coordinate']+['@SQ\tSN:'+c+'\tLN:80' for c in ('chr1','scaffold1','empty')]+['@RG\tID:'+sample+'\tSM:'+sample]
  for chrom in ('chr1','scaffold1'):
   for i in range(16):
    bases=list(seq)
    if i<(12 if sample=='a' else 8 if sample=='b' else 0):bases[29]='A' if bases[29]!='A' else 'T'
    rows.append('\t'.join([sample+chrom+str(i),'0',chrom,'1','60','80M','*','0','0',''.join(bases),'I'*80,'RG:Z:'+sample]))
  source=p/(sample+'.sam');source.write_text('\n'.join(rows)+'\n');bam=p/(sample+'.bam');run(sam,'sort','-o',bam,source);run(sam,'index',bam);bams.append(bam)
 listing=p/'bams.list';listing.write_text(''.join(str(x)+'\n' for x in bams))
 cfg={'shards':[{'index':0,'contigs':[['chr1',80]],'bed':str(p/'0.bed')},{'index':1,'contigs':[['scaffold1',80],['empty',80]],'bed':str(p/'1.bed')}],'total_contigs':3,'assembly_bases':240}
 for sh in cfg['shards']:pathlib.Path(sh['bed']).write_text(''.join(f'{c}\t0\t{n}\n' for c,n in sh['contigs']))
 fai=mod.validate_layout(cfg,str(ref)+'.fai')
 def call(dest,bed=None):
  mp=[bc,'mpileup','-f',ref,'-b',listing,'-q','20','-Q','20','-d','10000','-a','FORMAT/DP,FORMAT/AD','--seed','0','-Ou']
  if bed:mp+=['-R',bed]
  with open(p/'errors.log','ab') as log:
   first=subprocess.Popen(list(map(str,mp)),stdout=subprocess.PIPE,stderr=log)
   second=subprocess.Popen(list(map(str,[bc,'call','-m','-v','--ploidy','2','-Oz','-o',dest])),stdin=first.stdout,stderr=log)
   first.stdout.close();r2=second.wait();r1=first.wait();assert r1==r2==0
  run(bc,'index','-t',dest)
 full=p/'full.vcf.gz';call(full)
 part=[]
 for s in cfg['shards']:
  f=p/(str(s['index'])+'.vcf.gz');call(f,s['bed']);mod.validate_vcf(str(bc),f,['a','b','c'],set(c for c,n in s['contigs']),fai);part.append(f)
 joined=p/'joined.vcf.gz';run(bc,'concat','-Oz','-o',joined,*part);run(bc,'index','-t',joined)
 v=mod.validate_vcf(str(bc),joined,['a','b','c'],set(c for c,n in fai),fai)
 fmt='%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%INFO[\t%GT:%DP:%AD:%PL]\n'
 assert run(bc,'query','-f',fmt,full)==run(bc,'query','-f',fmt,joined)
 assert v['variant_records']>=2
 # An entirely empty contig must remain a valid shard, not silently dropped.
 emptybed=p/'empty.bed';emptybed.write_text('empty\t0\t80\n');empty=p/'empty.vcf.gz';call(empty,emptybed)
 assert mod.validate_vcf(str(bc),empty,['a','b','c'],{'empty'},fai)['variant_records']==0
 # Intentional overlap/omission and sample mismatch must fail closed.
 bad=json.loads(json.dumps(cfg));bad['shards'][1]['contigs'][0][0]='chr1'
 try:mod.validate_layout(bad,str(ref)+'.fai')
 except AssertionError:pass
 else:raise AssertionError('Duplicate/missingcontig accepted')
 try:mod.validate_vcf(str(bc),full,['b','a','c'],set(c for c,n in fai),fai)
 except AssertionError:pass
 else:raise AssertionError('Sampleorder mismatch accepted')
 result={'status':'PASS','synthetic_only':True,'checks':['Whole-contig plan exactly covers reference in order','Joint sample order and fullcontigheaders verified','Wholegenome versus fullcontigparallel calls: identicalCHROM/POS/REF/ALT/QUAL/INFO/GT/DP/AD/PL on tiny synthetic alignment','Orderedconcat count/sorting/duplicates verified','Zero-variantcontig shard accepted','Overlapping or missingcontig plan refused','Different sample order refused'],'study_data_used':False,'test_reference_contigs':3,'synthetic_variant_records':v['variant_records']}
 (base/'TEST_VALIDATION.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
