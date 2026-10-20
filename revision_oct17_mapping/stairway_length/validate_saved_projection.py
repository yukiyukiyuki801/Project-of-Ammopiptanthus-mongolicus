#!/usr/bin/env python3
"""Verify saved projected spectra against saved allele-count dictionaries."""
import ast,collections,json
from pathlib import Path
import numpy as np
from scipy.stats import hypergeom
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'revision_oct17_mapping/stairway_length'
audit=json.loads((OUT/'SAVED_INPUT_AUDIT.json').read_text());out=[]
for row in audit['target_sets']:
 count=collections.Counter()
 for l in (ROOT/row['source_datadict']).read_text().splitlines():
  _,d=l.split('{',1); calls=ast.literal_eval('{'+d)['calls']; count[tuple(next(iter(calls.values())))]+=1
 n=row['projection_chromosomes']; unfolded=np.zeros(n+1); contributing=0
 for (ref,alt),m in count.items():
  if ref+alt<n:continue
  unfolded += m*hypergeom.pmf(np.arange(n+1),ref+alt,alt,n)
  if min(ref,alt)>0:contributing+=m
 folded=unfolded[:n//2+1].copy()
 for k in range(n//2+1):
  if n-k!=k:folded[k]+=unfolded[n-k]
 matches=[]
 for bp in audit['blueprints']:
  for saved in bp['matching_saved_sfs']:
   if saved.get('target_set',{}).get('target_bed')==row['target_bed']:
    p=ROOT/saved['path']; vec=np.array(list(map(float,p.read_text().splitlines()[1].split())))[1:n//2+1]
    matches.append({'saved_sfs':saved['path'],'max_abs_bin_difference':float(np.max(np.abs(folded[1:]-vec)))})
 maxerr=max(x['max_abs_bin_difference'] for x in matches)
 assert maxerr<1e-7,(row['source_datadict'],maxerr)
 assert contributing==row['positive_SFS_contributor_sites']
 out.append({'target_bed':row['target_bed'],'projection_chromosomes':n,'contributing_coordinates':contributing,'projected_polymorphic_expected_sum':float(folded[1:].sum()),'max_abs_bin_difference':maxerr,'saved_matches':matches})
(OUT/'PROJECTION_RECOVERY_VALIDATION.json').write_text(json.dumps({'status':'PASS','method':'Independent hypergeometric projection of saved per-group REF/ALT counts; fold complementary counts; mask monomorphic bin. Saved results audit only, not demographic rerun.','targets':out},indent=2)+'\n')
print(json.dumps({'status':'PASS','target_sets':len(out),'max_abs_bin_difference':max(r['max_abs_bin_difference'] for r in out)}))
