#!/usr/bin/env python3
"""Convert completed targeted SIFT outputs and exact selected consequences to R inputs."""
import csv,json,pathlib,re,collections,datetime,hashlib
import numpy as np
ROOT=pathlib.Path.cwd(); BASE=ROOT/'runs/revision_run_01/outputs/functional_reconstruction'; OUT=BASE/'postprocess'
def read(p):return list(csv.DictReader(p.open()))
def write(name,rows):
 with (OUT/name).open('w') as h:
  w=csv.DictWriter(h,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def main():
 execution=json.loads((BASE/'sift_execution.json').read_text());assert execution['validated'], 'Engine completion must be audited by parent first'
 OUT.mkdir(exist_ok=True)
 variants=read(BASE/'consequences_by_variant.csv');assert len(variants)==38124
 samples=read(ROOT/'runs/revision_run_01/outputs/inputs/sample_metadata.csv');samples=sorted(samples,key=lambda x:int(x['vcf_index']));assert [int(x['vcf_index']) for x in samples]==list(range(1,218))
 D=np.load(BASE/'variant_dosages.npy');assert D.shape==(38124,217) and set(np.unique(D))<={-1,0,1,2}
 predictions={};audit=[];unknown=[]
 expected={p.stem:set(p.read_text().split()) for p in (BASE/'substitutions').glob('*.subst')}
 for tid,substitutions in sorted(expected.items()):
  path=BASE/'sift_predictions'/f'{tid}.SIFTprediction';found=set();warnings=[];parsed=[]
  if path.exists():
   for line in path.read_text().splitlines():
    if not line.strip():continue
    if line.startswith('WARNING!'):warnings.append(line);continue
    fields=line.split()
    if len(fields)!=6 or not re.fullmatch('[A-Z][0-9]+[A-Z]',fields[0]) or fields[1] not in ('DELETERIOUS','TOLERATED'):unknown.append((str(path),line));continue
    sub,label,score,median,npos,ntotal=fields;score=float(score);median=float(median);npos=int(npos);ntotal=int(ntotal)
    assert sub in substitutions and sub not in found and 0<=score<=1 and 0<=median<=4.33 and 0<=npos<=ntotal
    # Emitted labels reflect unrounded score; printed 0.05 alone cannot reproduce the decision.
    assert (label=='DELETERIOUS' and score<=.05) or (label=='TOLERATED' and score>=.05)
    found.add(sub);parsed.append(dict(transcript=tid,substitution=sub,prediction=label,printed_score=score,median_sequence_information=median,sequences_at_position=npos,total_sequences=ntotal,printed_threshold_boundary=score==.05))
   assert found==substitutions, f'Incomplete substitution output: {tid}'
   for row in parsed:
    pos=re.search(r'\d+',row['substitution']).group();w=[x for x in warnings if re.search(r'(?<!\d)'+pos+r'(?!\d)',x)]
    row['reference_or_position_warning']=' | '.join(w)
    row['strict_quality']=row['median_sequence_information']<=3.25 and row['sequences_at_position']>=2 and not w
    predictions[(tid,row['substitution'])]=row
  audit.append(dict(transcript=tid,requested_substitutions=len(substitutions),prediction_file=path.exists(),parsed_substitutions=len(found),warnings=' | '.join(warnings),missing_status='' if path.exists() else 'No output; engine skips queries lacking selected alignments; not tolerated'))
 assert not unknown,unknown[:5]
 assert predictions,'No scored substitutions'
 write('SIFT_predictions_parsed.csv',list(predictions.values()));write('SIFT_protein_coverage.csv',audit)
 chosen=[];sitecounts=[];populations=sorted(set(s['population'] for s in samples))
 for i,v in enumerate(variants):
  assert int(v['index'])==i
  effect=v['consequence'];label=effect;score={}
  if effect=='missense':
   key=(v['transcript'],v['aa_ref']+v['aa_position']+v['aa_alt']);score=predictions.get(key,{})
   label='SIFT_'+score['prediction'].lower() if score else 'missense_unscored'
  classes=[]
  if effect=='synonymous':classes=['synonymous']
  elif effect in ('stop_gained','splice_donor','splice_acceptor'):classes=['strict_LOF']
  elif effect=='missense':
   classes=['missense_all',label]
   if score and score['strict_quality']:classes.append(label+'_strict_quality')
  if classes:
   for pop in populations:
    cols=[j for j,s in enumerate(samples) if s['population']==pop];calls=D[i,cols];called=calls[calls>=0]
    sitecounts.append(dict(index=i,population=pop,group=samples[cols[0]]['group'],class_labels=';'.join(classes),population_individuals=len(cols),called_individuals=len(called),called_chromosomes=2*len(called),ALT_copies=int(called.sum()),homozygous_ALT_sites=int((called==2).sum())))
  chosen.append({**v,'SIFT_class':label,'SIFT_prediction':score.get('prediction',''),'SIFT_printed_score':score.get('printed_score',''),'SIFT_median_information':score.get('median_sequence_information',''),'SIFT_sequences_at_position':score.get('sequences_at_position',''),'SIFT_strict_quality':score.get('strict_quality','')})
 write('selected_variant_SIFT_consequences.csv',chosen);write('population_variant_counts.csv',sitecounts)
 # Per-individual dosage ratios retain fixed ALT sites; observed called genotypes only.
 rows=[]
 for j,s in enumerate(samples):
  r={k:s[k] for k in ('sample_id','population','group')}
  for cls in ('synonymous','strict_LOF','SIFT_deleterious','SIFT_tolerated','SIFT_deleterious_strict_quality','SIFT_tolerated_strict_quality'):
   ids=[]
   for i,v in enumerate(chosen):
    hit=(cls=='synonymous' and v['consequence']=='synonymous') or (cls=='strict_LOF' and v['consequence'] in ('stop_gained','splice_donor','splice_acceptor')) or (v['SIFT_class']==cls) or (cls.endswith('_strict_quality') and v['SIFT_class']==cls.removesuffix('_strict_quality') and v['SIFT_strict_quality'] is True)
    if hit:ids.append(i)
   calls=D[ids,j];calls=calls[calls>=0];r[cls+'_called_sites']=len(calls);r[cls+'_ALT_dosage']=int(calls.sum());r[cls+'_homozygous_ALT']=int((calls==2).sum())
  rows.append(r)
 write('individual_functional_counts.csv',rows)
 marker=dict(validated=True,completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),reviewer_basis='R3.M1 and explicit user request for population segregating-site and ALT-dosage ratio alternatives',predicted_proteins=sum(x['prediction_file'] for x in audit),requested_proteins=len(audit),parsed_substitutions=len(predictions),selected_missense_variants=sum(v['consequence']=='missense' for v in chosen),selected_missense_scored=sum(v['SIFT_prediction']!='' for v in chosen),selected_missense_strict_quality=sum(v['SIFT_strict_quality'] is True for v in chosen),rules='Most-severe canonical chosen variant row ONLY; SIFT emitted label uses raw score<0.05; equality tolerated. Printed score has two decimals, so use emitted label. All valid predictions primary; strict-quality sensitivity median<=3.25,>=2 sequences at position,no warning. Missing is unscored.',limits='REF-relative ALT, no ancestral polarization or fitness load. SwissProt reconstructed database differs from published UniRef100. No invariant/callable genome denominator.',script_sha256=hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest())
 (OUT/'score_parsing_validation.json').write_text(json.dumps(marker,indent=2)+'\n')
 print(json.dumps(marker,indent=2))
if __name__=='__main__':main()
