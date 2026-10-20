#!/usr/bin/env python3
"""Numerical check against unmodified pyRona; no graphics are produced."""
import csv,json,os
from pathlib import Path
import numpy as np
from pyRona import file_parser as fp
from pyRona.pyRona import RonaClass, calculate_rona
run=Path(os.environ['AMMO_REVISION_RUN']);out=run/'outputs/rona'
with (out/'rona_population_long.csv').open() as f:rows=list(csv.DictReader(f))
target=[r for r in rows if r['bio']=='bio4' and r['ssp']=='ssp126' and r['period']=='2041-2070' and r['gcm']=='ENSEMBLE']
pops=[r['population'] for r in target]; expected=np.array([float(r['rona']) for r in target])
pvals=np.loadtxt(run/'outputs/lfmm/associations_K4.csv',delimiter=',')
associations=[(str(i+1),'bio4') for i,p in enumerate(pvals[:,0]) if p<1]
names,alleles=fp.lfmm_to_pop_allele_freqs(str(run/'outputs/inputs/adaptive_imputed_population_grouped.tsv'),str(out/'pyrona_grouped_current.tsv'),associations,popnames=True)
assert names==pops
present=fp.parse_lfmm_envfile(str(out/'pyrona_grouped_current.tsv'))[0]
with (out/'population_climate_change.csv').open() as f:climate=list(csv.DictReader(f))
future=np.array([float(r['future']) for r in climate if r['bio']=='bio4' and r['ssp']=='ssp126' and r['period']=='2041-2070' and r['gcm']=='ENSEMBLE'])
RonaClass.POP_NAMES=np.array(pops)
rona=RonaClass('bio4');excluded=0
for locus,frequency in alleles.items():
 if np.ptp(frequency)==0 or not np.isfinite(frequency).all():excluded+=1;continue
 weight=np.corrcoef(present,frequency)[0,1]**2
 if not np.isfinite(weight) or weight<=0:excluded+=1;continue
 calculate_rona(locus,rona,present.copy(),future.copy(),frequency.copy(),None,False,'absdiff')
rona.basic_stats(True)
actual=np.asarray(rona.avg_ronas)
assert actual.shape==expected.shape
assert np.allclose(actual,expected,rtol=1e-9,atol=1e-11),(actual,expected)
report=dict(validated=True,populations=len(pops),loci_used=len(rona.pop_ronas),excluded_loci=excluded,max_abs_difference=float(np.max(abs(actual-expected))),reference='Installed unmodified pyRona functions, absdiff with weighted means; population-contiguous paired genotype/environment rows',graphics_generated=False)
(out/'pyrona_numeric_validation.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report),flush=True)
