#!/usr/bin/env python3
"""Validate keyed inputs, preserve called genotypes and stage numerical matrices."""
import csv
import hashlib
import json
import os
from collections import Counter
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
B=ROOT/'ammo_analysis_package_20261014'
RUN=Path(os.environ['AMMO_REVISION_RUN'])
OUT=RUN/'outputs/inputs';OUT.mkdir(parents=True,exist_ok=True)
def read(path):return list(csv.DictReader(path.open()))
def write(name,rows):
 with (OUT/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def digest(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  while chunk:=f.read(1048576):h.update(chunk)
 return h.hexdigest()
def parse_vcf(path,expected):
 keys=[];genotypes=[];quality=[]
 with path.open() as f:
  for line in f:
   if line.startswith('##'):continue
   fields=line.rstrip().split('\t')
   if line.startswith('#CHROM'):
    sample_ids=fields[9:];assert set(sample_ids)==set(expected)
    order=[sample_ids.index(s) for s in expected];continue
   assert len(fields[3])==len(fields[4])==1 and fields[3] in 'ACGT' and fields[4] in 'ACGT'
   keys.append((fields[0],int(fields[1]),fields[3],fields[4]))
   fmt=fields[8].split(':');gtindex=fmt.index('GT');dosages=[]
   for j in order:
    gt=fields[9+j].split(':')[gtindex].replace('|','/').split('/')
    assert len(gt)==2
    if '.' in gt:dosages.append(9)
    else:
     assert all(x in ('0','1') for x in gt)
     dosages.append(sum(map(int,gt)))
   genotypes.append(dosages)
 assert len(set(keys))==len(keys)
 return keys,np.asarray(genotypes,dtype=np.int8).T,sample_ids

manifest=read(B/'ammo/rona/validated_sample_manifest.csv')
ids=[r['sample_id'] for r in manifest]
assert len(ids)==len(set(ids))==217
assert [int(r['vcf_index']) for r in manifest]==list(range(1,218))
classes=read(ROOT/'recovery_audit/published_2022_population_classification.csv')
classmap={r['population']:r for r in classes}
coords=read(B/'ammo/gf/ammo_xy_19pop.csv')
popmap={r['code']:r for r in coords};pops=sorted(popmap)
counts=Counter(r['population'] for r in manifest)
assert len(pops)==19 and set(counts)==set(pops)==set(classmap)
for pop in pops:
 assert popmap[pop]['group1']==classmap[pop]['published_2022_group']
 assert counts[pop]==int(classmap[pop]['sample_count'])
population_rows=[dict(population=p,n=counts[p],longitude=float(popmap[p]['long']),latitude=float(popmap[p]['lat']),altitude=float(popmap[p]['Alt']),group=popmap[p]['group1'],genetic_cluster='UNVERIFIED') for p in pops]
write('population_metadata.csv',population_rows)
sample_rows=[dict(sample_id=r['sample_id'],friendly_id=r['renamed_sample_id'],population=r['population'],group=popmap[r['population']]['group1'],vcf_index=int(r['vcf_index']),longitude=float(r['longitude']),latitude=float(r['latitude']),altitude=float(r['field_elevation'])) for r in manifest]
write('sample_metadata.csv',sample_rows)

adaptive_path=B/'ammo/ammo_outliers_3overlap.recode.vcf'
full_path=B/'cluster_recovery/personal2/chong/Ammo/vcftools/Amon.LDfilter_08.vcf'
adaptive,geno,adapt_order=parse_vcf(adaptive_path,ids)
full,fullgeno,full_order=parse_vcf(full_path,ids)
assert adapt_order==ids and full_order==ids
assert len(adaptive)==1637 and len(full)==38124
fullindex={k:i for i,k in enumerate(full)}
assert set(adaptive)<=set(full)
adapt_in_full=[fullindex[k] for k in adaptive]
assert np.array_equal(geno,fullgeno[:,adapt_in_full])
archived=np.loadtxt(B/'ammo/rona/ammo_outliers_3overlap.recode.lfmm_imputed.lfmm',dtype=np.int8)
assert archived.shape==geno.shape
assert np.isin(archived,[0,1,2]).all()
called=geno!=9
changed=int(np.sum(archived[called]!=geno[called]))
assert changed==0,'Archived imputation changed called genotypes or has a wrong row order'
np.savetxt(OUT/'adaptive_called.tsv',geno,fmt='%d',delimiter='\t')
np.savetxt(OUT/'adaptive_imputed.tsv',archived,fmt='%d',delimiter='\t')
np.save(OUT/'full_genotypes.npy',fullgeno)
np.save(OUT/'nonadaptive_indices.npy',np.asarray([i for i,k in enumerate(full) if k not in set(adaptive)],dtype=np.int32))
write('adaptive_loci.csv',[dict(locus_index=i+1,chrom=k[0],pos=k[1],ref=k[2],alt=k[3],full_locus_index=fullindex[k]+1) for i,k in enumerate(adaptive)])
write('full_loci.csv',[dict(locus_index=i+1,chrom=k[0],pos=k[1],ref=k[2],alt=k[3],is_adaptive=k in set(adaptive)) for i,k in enumerate(full)])

# pyRona groups genotypes by contiguous blocks, not by independent ID joins.
grouped=sorted(range(217),key=lambda i:(manifest[i]['population'],i))
write('pyrona_grouped_sample_order.csv',[dict(grouped_index=j+1,vcf_index=i+1,sample_id=ids[i],population=manifest[i]['population']) for j,i in enumerate(grouped)])
np.savetxt(OUT/'adaptive_imputed_population_grouped.tsv',archived[grouped],fmt='%d',delimiter='\t')
blocks=sum(i==0 or manifest[i]['population']!=manifest[i-1]['population'] for i in range(217))
assert blocks>19

env=np.loadtxt(B/'ammo/rona/ammo_current_corrected.txt')
expected=np.array([[float(r[k]) for k in ['bio4','bio5','bio6','bio15','bio17']] for r in manifest])
assert env.shape==(217,5) and np.allclose(env,expected,rtol=0,atol=1e-5)
np.savetxt(OUT/'current_environment_manifest.tsv',expected,delimiter='\t',fmt='%.8f')
fstrows=read(B/'ammo/gf/ammo_outlier_pairwise_fst.csv')
assert len(fstrows)==19 and set(fstrows[0])==set(pops)|{'pop'}
fstmap={r['pop']:r for r in fstrows};assert set(fstmap)==set(pops)
fst=np.array([[float(fstmap[p][q]) for q in pops] for p in pops])
assert np.isfinite(fst).all() and np.allclose(fst,fst.T,atol=1e-12)
assert np.allclose(np.diag(fst),0) and fst.max()<=1
write('adaptive_fst.csv',[dict(population=p,**{q:fst[i,j] for j,q in enumerate(pops)}) for i,p in enumerate(pops)])
sampleqc=[]
for i,r in enumerate(manifest):sampleqc.append(dict(sample_id=ids[i],population=r['population'],adaptive_missing_fraction=float(np.mean(geno[i]==9)),full_missing_fraction=float(np.mean(fullgeno[i]==9))))
write('genotype_missingness_by_sample.csv',sampleqc)
lineage=[dict(dataset='adaptive',source=str(adaptive_path.relative_to(ROOT)),sha256=digest(adaptive_path),samples=217,loci=1637,status='verified_sample_locus_and_dosage_match'),dict(dataset='full',source=str(full_path.relative_to(ROOT)),sha256=digest(full_path),samples=217,loci=38124,status='verified_original_coordinate_full_set'),dict(dataset='archived_imputed_adaptive',source='ammo_analysis_package_20261014/ammo/rona/ammo_outliers_3overlap.recode.lfmm_imputed.lfmm',sha256=digest(B/'ammo/rona/ammo_outliers_3overlap.recode.lfmm_imputed.lfmm'),samples=217,loci=1637,status='all_called_genotypes_preserved; historical_imputation_parameters_unverified')]
write('dataset_lineage.csv',lineage)
summary=dict(validated=True,samples=217,populations=19,core_populations=9,peripheral_populations=10,adaptive_loci=1637,full_loci=38124,nonadaptive_loci=len(full)-len(adaptive),adaptive_missing_calls=int(np.sum(~called)),adaptive_missing_fraction=float(np.mean(~called)),imputed_called_genotype_changes=changed,raw_population_contiguous_blocks=blocks,pyrona_grouping_corrected=True,negative_pairwise_fst=int(np.sum(fst<0)),fst_min=float(fst.min()),fst_max=float(fst.max()),unresolved=['Original archived neutral locus crosswalk','Demographic sample-cluster assignments and SFS/callable-L provenance','Functional-variant metric definitions','Historical imputation parameters'],scope='Primary metadata/genotype prerequisites verified; not completion of all A01 lineage questions.')
(OUT/'input_validation.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2),flush=True)
