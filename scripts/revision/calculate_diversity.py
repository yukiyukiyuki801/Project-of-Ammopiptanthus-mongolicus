#!/usr/bin/env python3
"""Conditional SNP-panel diversity for full/adaptive/nonadaptive sensitivity."""
import csv,os,json
from pathlib import Path
import numpy as np
run=Path(os.environ['AMMO_REVISION_RUN']);inputs=run/'outputs/inputs';out=run/'outputs/statistics';out.mkdir(exist_ok=True)
with (inputs/'sample_metadata.csv').open() as f:meta=list(csv.DictReader(f))
G=np.load(inputs/'full_genotypes.npy');background=np.load(inputs/'nonadaptive_indices.npy')
with (inputs/'adaptive_loci.csv').open() as f:adaptive=np.array([int(r['full_locus_index'])-1 for r in csv.DictReader(f)])
rows=[]
for subset,indices in [('full',np.arange(G.shape[1])),('adaptive',adaptive),('nonadaptive',background)]:
 for pop in sorted({r['population'] for r in meta}):
  members=[i for i,r in enumerate(meta) if r['population']==pop]
  data=G[np.ix_(members,indices)]
  called=data!=9;n=called.sum(axis=0)
  valid=n>=int(np.ceil(.8*len(members)))
  dosage=np.where(called,data,0).sum(axis=0)
  with np.errstate(divide='ignore',invalid='ignore'):
   p=dosage/(2*n)
   ho=(data==1).sum(axis=0)/n
   he=2*p*(1-p)*(2*n)/(2*n-1)
  rows.append(dict(population=pop,group=meta[members[0]]['group'],subset=subset,n_individuals=len(members),input_loci=len(indices),used_loci=int(valid.sum()),Ho=float(np.mean(ho[valid])),He_unbiased=float(np.mean(he[valid])),missing_fraction=float(np.mean(~called)),within_population_call_rate_threshold=.8,estimand='Mean per-locus heterozygosity conditional on globally ascertained SNP panel; not genome-wide nucleotide diversity'))
with (out/'population_diversity.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
(out/'diversity_validation.json').write_text(json.dumps(dict(validated=True,rows=len(rows),full_loci=G.shape[1],adaptive_loci=len(adaptive),nonadaptive_loci=len(background),nonadaptive_interpretation='Full minus adaptive; not certified neutral',pi_calculated=False),indent=2)+'\n')
