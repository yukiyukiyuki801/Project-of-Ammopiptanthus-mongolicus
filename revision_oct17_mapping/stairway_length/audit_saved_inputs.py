#!/usr/bin/env python3
"""Read-only recovery of archived Stairway inputs; no demographic fitting."""
import ast,csv,hashlib,json,math,re
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]; OUT=ROOT/'revision_oct17_mapping/stairway_length'
BASE=ROOT/'ammo_analysis_package_20261014/cluster_recovery/personal2/chong/Ammo/Demography'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def bp_parse(p):
 d={}
 for l in p.read_text().splitlines():
  l=l.split('#')[0].strip()
  if ':' in l:
   k,v=l.split(':',1);d[k.strip()]=v.strip()
 return d
sfs=[]
for p in (BASE/'easySFS').rglob('*.sfs'):
 ls=p.read_text().splitlines()
 try:
  n=int(ls[0].split()[0])-1; x=list(map(float,ls[1].split()));v=x[1:n//2+1]
  if len(x)!=n+1:continue
  sfs.append((p,n,v))
 except (ValueError,IndexError):continue
rows=[]; sources={}; targetsets={}
for p in sorted((BASE/'stariwayplot').rglob('*.blueprint')):
 d=bp_parse(p); n=int(d['nseq']); v=list(map(float,d['SFS'].split())); L=int(d['L'])
 matches=[f for f,k,x in sfs if k==n and len(x)==len(v) and max(abs(a-b) for a,b in zip(x,v))<1e-8]
 entry={'blueprint':str(p.relative_to(ROOT)),'blueprint_sha256':sha(p),'population':d['popid'],'projection_chromosomes':n,'historical_L_bp':L,'projected_segregating_sum':sum(v),'L_equals_round_SFS_times150':L==round(sum(v))*150,'matching_saved_sfs':[]}
 for f in matches:
  dd=f.parent.parent/'datadict.txt'; obj={'path':str(f.relative_to(ROOT)),'sha256':sha(f),'datadict_exists':dd.exists()}
  if dd.exists():
   ds=sha(dd);obj.update(datadict=str(dd.relative_to(ROOT)),datadict_sha256=ds)
   key=ds+'_'+str(n)
   if key not in targetsets:
    rec=[]
    for line in dd.read_text().splitlines():
     coord,content=line.split('{',1); chrom,pos=coord.rsplit('-',1); calls=ast.literal_eval('{'+content)['calls']; assert len(calls)==1
     group,(ref,alt)=next(iter(calls.items())); rec.append((chrom,int(pos),int(ref),int(alt)))
    contributing=[r for r in rec if r[2]>0 and r[3]>0 and r[2]+r[3]>=n]
    fn='targets_'+ds[:12]+'_proj'+str(n)+'.bed'
    (OUT/fn).write_text(''.join(f'{c}\t{pos-1}\t{pos}\t{ref}\t{alt}\n' for c,pos,ref,alt in sorted(contributing)))
    targetsets[key]={'source_datadict':str(dd.relative_to(ROOT)),'sha256':ds,'projection_chromosomes':n,'all_datadict_sites':len(rec),'positive_SFS_contributor_sites':len(contributing),'target_bed':fn,'target_bed_sha256':sha(OUT/fn),'meaning':'BED columns chrom/start0/end1/REFcount/ALTcount; both alleles present and called chromosomes >= projection. Positive expected polymorphism after projection; not one realized random SNP subset.'}
   obj['target_set']=targetsets[key]
  entry['matching_saved_sfs'].append(obj)
 rows.append(entry)
(OUT/'SAVED_INPUT_AUDIT.json').write_text(json.dumps({'blueprints':rows,'target_sets':list(targetsets.values()),'all_23_L_match_rounded_projected_SFS_times150':all(r['L_equals_round_SFS_times150'] for r in rows),'interpretation':'Arithmetic identity is established; historical derivation script and historical BAM length measurement are not recovered. New remapping lengths cannot be described as historical numeric inputs without additional evidence.'},indent=2)+'\n')
with (OUT/'HISTORICAL_L.csv').open('w') as o:
 w=csv.DictWriter(o,fieldnames=['blueprint','population','projection_chromosomes','historical_L_bp','projected_segregating_sum','L_equals_round_SFS_times150','matched_sfs_files']);w.writeheader()
 for r in rows:w.writerow({**{k:r[k] for k in w.fieldnames[:-1]},'matched_sfs_files':len(r['matching_saved_sfs'])})
print(json.dumps({'blueprints':len(rows),'all_rounded_identity':all(r['L_equals_round_SFS_times150'] for r in rows),'sfs_unmatched':[r['blueprint'] for r in rows if not r['matching_saved_sfs']],'unique_target_sets':len(targetsets)}))
