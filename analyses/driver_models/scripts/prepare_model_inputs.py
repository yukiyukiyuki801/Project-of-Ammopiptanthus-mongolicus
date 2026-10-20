#!/usr/bin/env python3
"""Select explicitly audited predictors and preserve the existing spatial block partition."""
import csv,sys,json,hashlib
from pathlib import Path
if len(sys.argv)!=4 or '--help' in sys.argv:
    print('Usage: python prepare_model_inputs.py audited_population_case.csv output.csv geographic_blocks.csv')
    sys.exit(0 if '--help' in sys.argv else 2)
src=Path(sys.argv[1]); dest=Path(sys.argv[2]);blocks_src=Path(sys.argv[3])
blocks={}
for r in csv.DictReader(blocks_src.open()):
    p,b=r['population'],r['heldout_block']
    if p in blocks and blocks[p]!=b:raise ValueError('Inconsistent saved block')
    blocks[p]=b
mapping={'climate_change':'climate_change','private_alleles':'Private_allele_archived_full','Ho':'Ho_full','He':'He_unbiased_full','Fis':'Fis_archived_full','Pi':'Pi_archived_full','Lon':'Lon','Lat':'Lat','Alt':'Alt','DEL':'DEL_ratio','LoF':'LOF_ratio'}
rows=[]
for r in csv.DictReader(src.open()):
    z={'case_id':'_'.join(r[k] for k in ['ssp','period','bio']),'population':r['population'],'response':r['rona'],'block':blocks[r['population']]}
    z.update({new:r[old] for new,old in mapping.items()});rows.append(z)
assert len(rows)==380 and len({r['case_id'] for r in rows})==20
dest.parent.mkdir(parents=True,exist_ok=True)
with dest.open('w',newline='')as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
manifest={'input':str(src),'input_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'output':str(dest),'output_sha256':hashlib.sha256(dest.read_bytes()).hexdigest(),'predictor_mapping':mapping,'block_source':str(blocks_src),'block_populations':blocks,'n_predictors':11,'n_subsets':2048,'definition_limits':['Private,Fis,Pi are original saved full-panel metrics; Ho/He are validated full-panel recomputations.','Pi and He are near-redundant but not numerically identical; both included as requested, coefficient/inclusion claims must discuss collinearity.','DEL and LoF are corrected reference-relative synonymous-normalized ratios; not measured fitness load.','Every response is fixed previously estimated RONA: downstream nested selection does not refit candidate detection/RONA.']}
Path(str(dest)+'.manifest.json').write_text(json.dumps(manifest,indent=2))
print(json.dumps({'rows':len(rows),'cases':20,'predictors':11,'subsets':2048,'output':str(dest)}))
