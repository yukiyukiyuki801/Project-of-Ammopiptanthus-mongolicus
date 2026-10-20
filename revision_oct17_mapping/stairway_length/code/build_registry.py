#!/usr/bin/env python3
import ast,csv,json,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]; BASE=ROOT/'revision_oct17_mapping/stairway_length'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 audit=json.loads((BASE/'SAVED_INPUT_AUDIT.json').read_text()); meta=list(csv.DictReader((ROOT/'revision_oct17_mapping/inputs/PCA_scores_and_group_assignments.csv').open()))
 groups={'all':[r['sample_id'] for r in meta],'population_core':[r['sample_id'] for r in meta if r['group']=='core'],'population_peripheral':[r['sample_id'] for r in meta if r['group']=='peripheral'],'adaptive_peripheral':[r['sample_id'] for r in meta if r['group']=='peripheral' and r['archived_group']!='peripheral'],'genetic_nonadaptive':[r['sample_id'] for r in meta if r['archived_group'] in ('core','peripheral')]}
 for g in ('core','peripheral','adaptive_1','adaptive_2'):groups['genetic_'+g]=[r['sample_id'] for r in meta if r['archived_group']==g]
 vcf=ROOT/'ammo_analysis_package_20261014/cluster_recovery/personal2/chong/Ammo/vcftools/Amon.LDfilter_08.vcf';counts={g:{} for g in groups}
 for line in vcf.open():
  if line.startswith('#CHROM'):
   ids=line.rstrip().split('\t')[9:];idx={g:[ids.index(s) for s in ss] for g,ss in groups.items()}
  elif not line.startswith('#'):
   f=line.rstrip().split('\t');key=(f[0],int(f[1]));gts=[x.split(':')[0].replace('|','/').split('/') for x in f[9:]]
   for g,ix in idx.items():
    al=[a for i in ix for a in gts[i] if a!='.'];counts[g][key]=(al.count('0'),al.count('1'))
 targets=[]
 for t in audit['target_sets']:
  p=ROOT/t['source_datadict'];d={}
  for l in p.read_text().splitlines():
   name,s=l.split('{',1);c,pos=name.rsplit('-',1);x=ast.literal_eval('{'+s);d[(c,int(pos))]=(tuple(next(iter(x['calls'].values()))),x['segregating'])
  shared=set(d)&set(counts['all']);assert len(shared)>10000,'Insufficient independent shared-coordinate membership evidence'
  matches=[g for g in groups if all(counts[g][k]==d[k][0] for k in shared)]
  assert len(matches)==1,(p,matches,len(shared));g=matches[0]
  refout=BASE/(t['target_bed'].removesuffix('.bed')+'.alleles.tsv')
  sites=[(f[0],int(f[2])) for f in (l.split('\t') for l in (BASE/t['target_bed']).read_text().splitlines())]
  refout.write_text(''.join(f'{c}\t{pos}\t{d[(c,pos)][1][0]}\t{d[(c,pos)][1][1]}\n' for c,pos in sites))
  targets.append({**t,'membership_group':g,'samples':groups[g],'membership_evidence_shared_sites':len(shared),'membership_mismatches':0,'alleles_tsv':refout.name,'alleles_sha256':sha(refout)})
 result={'status':'VERIFIED_MEMBERSHIP','membership_validation_vcf':str(vcf),'membership_validation_vcf_sha256':sha(vcf),'groups':groups,'targets':targets,'note':'Population and genetic groups intentionally distinct. Candidate memberships are verified by identical archived REF/ALT counts at every coordinate shared with inherited 38124 VCF; extra archived coordinates do not independently validate additional historical genotypes.'}
 (BASE/'MEASUREMENT_REGISTRY.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'targets':len(targets),'groups':{g:len(s) for g,s in groups.items()}}))
if __name__=='__main__':main()
