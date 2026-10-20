#!/usr/bin/env python3
"""Variant-panel consequence reconstruction; no ancestral or fitness inference."""
import bisect,csv,gzip,hashlib,json,pathlib,re,collections,datetime
import numpy as np
ROOT=pathlib.Path.cwd(); RUN=ROOT/'runs/revision_run_01'; OUT=RUN/'outputs/functional_reconstruction'; OUT.mkdir(parents=True,exist_ok=True)
BASE=ROOT/'ammo_analysis_package_20261014/cluster_recovery/personal2/sun/sift4g/ammo'
VCF=ROOT/'ammo_analysis_package_20261014/cluster_recovery/personal2/chong/Ammo/vcftools/Amon.LDfilter_08.vcf'
GTF=BASE/'gene-annotation-src/ammo.gtf.gz'; REF=BASE/'chr-src/ammo.fasta.gz'; PROT=BASE/'all_prot.fasta'
COMP=str.maketrans('ACGTNacgtn','TGCANtgcan')
def rc(s):return s.translate(COMP)[::-1]
# Standard nuclear genetic code, codon order T,C,A,G at each position.
CODE=dict(zip((a+b+c for a in 'TCAG' for b in 'TCAG' for c in 'TCAG'),'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'))
def translate(s):return ''.join(CODE.get(s[i:i+3],'X') for i in range(0,len(s)-2,3))
def fasta(path):
 op=gzip.open if str(path).endswith('.gz') else open
 with op(path,'rt') as h:
  name=None;seq=[]
  for line in h:
   if line.startswith('>'):
    if name is not None:yield name,''.join(seq).upper()
    name=line[1:].split()[0];seq=[]
   else:seq.append(line.strip())
  if name is not None:yield name,''.join(seq).upper()
def write(name,rows,fields=None):
 with open(OUT/name,'w') as h:
  w=csv.DictWriter(h,fieldnames=fields or list(rows[0]));w.writeheader();w.writerows(rows)
variants=[];bychrom=collections.defaultdict(list);doses=[]
with open(VCF) as h:
 for line in h:
  if line.startswith('##'):continue
  f=line.rstrip().split('\t')
  if line.startswith('#CHROM'):samples=f[9:];continue
  assert len(f[3])==len(f[4])==1 and f[3] in 'ACGT' and f[4] in 'ACGT', 'Expected biallelic SNP panel'
  i=len(variants);variants.append(dict(index=i,chrom=f[0],position=int(f[1]),ref=f[3],alt=f[4]));bychrom[f[0]].append(i)
  gtidx=f[8].split(':').index('GT');row=[]
  for s in f[9:]:
   gt=re.split('[/|]',s.split(':')[gtidx]);row.append(sum(map(int,gt)) if len(gt)==2 and all(g in ['0','1'] for g in gt) else -1)
  doses.append(row)
assert len(variants)==38124 and len(samples)==217
D=np.array(doses,dtype=np.int8);del doses
pos={c:[variants[i]['position'] for i in ids] for c,ids in bychrom.items()}
for c in pos:assert pos[c]==sorted(pos[c])
def hits(c,a,b):
 ids=bychrom.get(c,[]);p=pos.get(c,[])
 return ids[bisect.bisect_left(p,a):bisect.bisect_right(p,b)]
tx={}
with gzip.open(GTF,'rt') as h:
 for line in h:
  if line.startswith('#'):continue
  f=line.rstrip().split('\t')
  if len(f)!=9 or f[2] not in ['CDS','exon']:continue
  attrs=dict(re.findall(r'(\w+) "([^"]+)"',f[8]));tid=attrs.get('transcript_id');gene=attrs.get('gene_id')
  assert tid and gene
  t=tx.setdefault(tid,dict(transcript=tid,gene=gene,chrom=f[0],strand=f[6],CDS=[],exon=[]))
  t[f[2]].append((int(f[3]),int(f[4]),int(f[7]) if f[7]!='.' else None))
proteins=dict(fasta(PROT));target_genes=set()
for t in tx.values():
 for a,b,_ in t['CDS']:
  if hits(t['chrom'],a-2,b+2):target_genes.add(t['gene']);break
selected={k:v for k,v in tx.items() if v['gene'] in target_genes and v['CDS']}
contigtx=collections.defaultdict(list)
for tid,t in selected.items():contigtx[t['chrom']].append(tid)
annotations=collections.defaultdict(list);audit=[];ref_checked=0
for chrom,seq in fasta(REF):
 for vi in bychrom.get(chrom,[]):
  v=variants[vi];assert seq[v['position']-1]==v['ref'],f'Reference mismatch {v}'
  ref_checked+=1
 for tid in contigtx.get(chrom,[]):
  t=selected[tid];ordered=sorted(t['CDS'],reverse=t['strand']=='-');firstphase=ordered[0][2];phase_ok=firstphase in [0,1,2];length=-firstphase if phase_ok else 0
  for a,b,phase in ordered:
   if (a,b,phase)!=ordered[0]:phase_ok &= phase==(3-length%3)%3
   length+=b-a+1
  cds=''.join(seq[a-1:b] if t['strand']=='+' else rc(seq[a-1:b]) for a,b,_ in ordered)
  cds=cds[firstphase:] if phase_ok else cds
  protein=translate(cds);expected=proteins.get(tid,'').rstrip('*');translation_ok=bool(expected) and protein.rstrip('*')==expected and '*' not in protein.rstrip('*') and 'X' not in protein
  valid=bool(phase_ok and translation_ok and firstphase==0 and len(cds)%3==0)
  t.update(valid=valid,protein=expected,cds_length=len(cds))
  audit.append(dict(transcript=tid,gene=t['gene'],strand=t['strand'],cds_bases=len(cds),first_phase=firstphase,phase_consistent=bool(phase_ok),protein_exact_match=translation_ok,valid_primary=valid))
  if not valid:continue
  offset=-firstphase;codonvariants=collections.defaultdict(list);local=[]
  for a,b,_ in ordered:
   for vi in hits(chrom,a,b):
    v=variants[vi];idx=offset+(v['position']-a if t['strand']=='+' else b-v['position']);codonidx=idx//3;codon=cds[codonidx*3:codonidx*3+3];alt=v['alt'] if t['strand']=='+' else rc(v['alt']);ref=v['ref'] if t['strand']=='+' else rc(v['ref'])
    assert codon[idx%3]==ref
    altcodon=codon[:idx%3]+alt+codon[idx%3+1:];aa=CODE[codon];alt_aa=CODE[altcodon]
    effect='synonymous' if aa==alt_aa else 'stop_gained' if alt_aa=='*' else 'stop_lost' if aa=='*' else 'start_lost' if codonidx==0 and aa=='M' else 'missense'
    local.append(dict(index=vi,transcript=tid,gene=t['gene'],consequence=effect,aa_ref=aa,aa_alt=alt_aa,aa_position=codonidx+1,codon_ref=codon,codon_alt=altcodon,multiple_variant_codon=False))
    codonvariants[codonidx].append(vi)
   offset+=b-a+1
  for row in local:row['multiple_variant_codon']=len(set(codonvariants[row['aa_position']-1]))>1;annotations[row['index']].append(row)
  # Essential intronic2bp only where both exons overlap annotated CDS.
  exons=sorted(t['exon']);cds_intervals=t['CDS']
  for left,right in zip(exons,exons[1:]):
   if left[1]>=right[0]-4:continue
   if not all(any(a<=ex[1] and b>=ex[0] for a,b,_ in cds_intervals) for ex in [left,right]):continue
   donor=(left[1]+1,left[1]+2);acceptor=(right[0]-2,right[0]-1)
   motifs=(seq[donor[0]-1:donor[1]],seq[acceptor[0]-1:acceptor[1]])
   canonical=motifs==('GT','AG') if t['strand']=='+' else motifs==('CT','AC')
   if not canonical:continue
   for bounds,kind in [(donor,'splice_donor' if t['strand']=='+' else 'splice_acceptor'),(acceptor,'splice_acceptor' if t['strand']=='+' else 'splice_donor')]:
    for vi in hits(chrom,*bounds):annotations[vi].append(dict(index=vi,transcript=tid,gene=t['gene'],consequence=kind,aa_ref='',aa_alt='',aa_position='',codon_ref='',codon_alt='',multiple_variant_codon=False))
assert ref_checked==38124
bygene=collections.defaultdict(list)
for t in selected.values():bygene[t['gene']].append(t)
canonical={g:sorted([t for t in ts if t['valid']],key=lambda t:(-t['cds_length'],t['transcript']))[0]['transcript'] for g,ts in bygene.items() if any(t['valid'] for t in ts)}
severity={c:i for i,c in enumerate(['stop_gained','splice_donor','splice_acceptor','stop_lost','start_lost','missense','synonymous'])}
chosen=[];allrows=[];queries=collections.defaultdict(set);site_classes=[]
for vi,v in enumerate(variants):
 rr=[r for r in annotations.get(vi,[]) if canonical.get(r['gene'])==r['transcript']]
 for r in annotations.get(vi,[]):allrows.append({**v,**r,'canonical_transcript':canonical.get(r['gene'])==r['transcript']})
 if rr:
  r=min(rr,key=lambda r:severity[r['consequence']]);effect=r['consequence']
  ambiguous=any(x['multiple_variant_codon'] for x in rr)
  effect='multiple_variant_codon_unresolved' if ambiguous else effect
  chosen.append({**v,**r,'consequence':effect,'overlapping_canonical_gene_count':len(set(x['gene'] for x in rr))})
  if not ambiguous:
   for x in rr:
    if x['consequence']=='missense':queries[x['transcript']].add(x['aa_ref']+str(x['aa_position'])+x['aa_alt'])
 else:effect='not_assigned_valid_canonical_coding_or_splice';chosen.append({**v,'transcript':'','gene':'','consequence':effect,'aa_ref':'','aa_alt':'','aa_position':'','codon_ref':'','codon_alt':'','multiple_variant_codon':False,'overlapping_canonical_gene_count':0})
 site_classes.append(effect)
write('transcript_validation.csv',audit);write('consequences_by_transcript.csv',allrows);write('consequences_by_variant.csv',chosen)
subdir=OUT/'substitutions';subdir.mkdir(exist_ok=True)
with open(OUT/'target_proteins.fasta','w') as h:
 for tid in sorted(queries):
  h.write('>'+tid+'\n'+selected[tid]['protein']+'\n');(subdir/(tid+'.subst')).write_text('\n'.join(sorted(queries[tid]))+'\n')
metadata={r['sample_id']:r for r in csv.DictReader(open(RUN/'outputs/inputs/sample_metadata.csv'))};burdens=[];C=np.array(site_classes)
for j,sample in enumerate(samples):
 row={'sample_id':sample,'population':metadata[sample]['population'],'group':metadata[sample]['group']}
 for label,classes in {'stop_gained':['stop_gained'],'essential_splice':['splice_donor','splice_acceptor'],'stop_lost':['stop_lost'],'start_lost':['start_lost'],'missense':['missense'],'synonymous':['synonymous'],'strict_predicted_LOF':['stop_gained','splice_donor','splice_acceptor']}.items():
  mask=np.isin(C,classes);calls=D[mask,j];calls=calls[calls>=0]
  row[label+'_called_sites']=int(len(calls));row[label+'_ALT_dosage']=int(calls.sum());row[label+'_homozygous_ALT_sites']=int((calls==2).sum())
 denom=row['synonymous_called_sites']
 row['strict_predicted_LOF_ALT_dosage_per_called_synonymous_allele_opportunity']=row['strict_predicted_LOF_ALT_dosage']/(2*denom) if denom else ''
 row['polarization']='Reference-relative ALT; ancestral state unverified';burdens.append(row)
write('reference_relative_burden_by_sample.csv',burdens)
np.save(OUT/'variant_dosages.npy',D)
settings={'validated':True,'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'reviewer_basis':'R3.M1 exact original request for LOF consequences, derived-allele evidence and LOF/DEL calculation definitions','vcf_snps':len(variants),'samples':len(samples),'all_reference_bases_match':ref_checked,'target_genes':len(target_genes),'valid_canonical_genes':len(canonical),'sift_target_proteins':len(queries),'sift_substitutions':sum(map(len,queries.values())),'consequence_counts':dict(collections.Counter(site_classes)),'rules':'Longest valid protein-matching phase-consistent transcript per gene, lexical ID tie-break; most severe consequence across overlapping chosen genes; multi-SNP codons excluded from primary consequence classes because haplotype unresolved; canonical GT-AG coding-adjacent intron2bp only; no indels present','limits':'Variant-only panel; no invariant/callable denominator, no BAM/DP/GQ mapping audit, no ancestral polarization, no measured fitness load. Stoploss/startloss separate from strict predicted LOF. SIFT not yet run; missense is NOT deleterious classification. SwissProt latest available is new reconstruction, not historical database reproduction.','sources':{str(p.relative_to(ROOT)):{'bytes':p.stat().st_size} for p in [VCF,GTF,REF,PROT]},'source_sha256':hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()}
(OUT/'consequence_validation.json').write_text(json.dumps(settings,indent=2)+'\n');print(json.dumps(settings,indent=2))
