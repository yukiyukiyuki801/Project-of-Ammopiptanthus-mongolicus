#!/usr/bin/env python3
"""Parent-queued whole-contig joint calling; no scheduler submission."""
import argparse, collections, csv, datetime, hashlib, json, os, pathlib, re, subprocess, sys, traceback

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[1]

def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(8*1024*1024), b''):
            h.update(b)
    return h.hexdigest()

def dump(path, data):
    path = pathlib.Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(data, indent=2)+'\n'); tmp.replace(path)

def command(args):
    print(json.dumps(list(map(str, args))), flush=True)
    subprocess.run(list(map(str, args)), check=True)

def output(args):
    return subprocess.check_output(list(map(str, args)), text=True)

def assert_products(marker):
    assert marker['validated'] is True
    for path, digest in marker['outputs'].items():
        assert pathlib.Path(path).is_file() and sha(path)==digest, 'Changed product: '+path

def original_fingerprint(cfg, original):
    return dict(script=sha(original['original_workflow_source']), config=sha(original['original_config']),
                manifest=sha(cfg['manifest']), stage='call',
                pca_script=sha(ROOT/'revision_oct17_mapping/workflow/pca_comparison.R'),
                existing_pca_scores=sha(cfg['existing_pca_scores']),
                existing_pca_variance=sha(cfg['existing_pca_variance']),
                reference_audit=sha(cfg['reference_audit']))

def validate_layout(config, fai):
    expected = [(v[0], int(v[1])) for l in pathlib.Path(fai).read_text().splitlines() if (v:=l.split('\t'))]
    actual = [(c,n) for s in config['shards'] for c,n in s['contigs']]
    assert actual == expected and len(actual)==config['total_contigs']
    assert len(set(c for c,n in actual))==len(actual)
    assert sum(n for c,n in actual)==config['assembly_bases']
    assert [s['index'] for s in config['shards']]==list(range(len(config['shards'])))
    for s in config['shards']:
        assert pathlib.Path(s['bed']).suffix=='.bed'
        assert pathlib.Path(s['bed']).read_text()==''.join(f'{c}\t0\t{n}\n' for c,n in s['contigs'])
    return expected

def validate_vcf(bc, path, samples, allowed, fai):
    assert output([bc,'query','-l',path]).splitlines()==samples, 'Sample order mismatch'
    # Empty variant shards are legitimate; still require a complete valid header.
    h=output([bc,'view','-h',path])
    header_contigs=[]
    for line in h.splitlines():
        if line.startswith('##contig=<ID='):
            match=re.match(r'##contig=<ID=([^,>]+),length=([0-9]+)',line)
            assert match, 'Unrecognized contig header'
            header_contigs.append((match[1],int(match[2])))
    assert header_contigs==fai, 'Reference header order/length mismatch'
    ranks={c:i for i,(c,n) in enumerate(fai)}; lengths=dict(fai)
    p=subprocess.Popen([str(bc),'query','-f','%CHROM\t%POS\t%REF\t%ALT\n',str(path)],stdout=subprocess.PIPE,text=True)
    prev=(-1,-1); seen=set(); counts=collections.Counter()
    try:
        for line in p.stdout:
            c,pos,ref,alt=line.rstrip('\n').split('\t'); pos=int(pos)
            assert c in allowed and 1<=pos<=lengths[c], 'Record outside assigned contigs'
            at=(ranks[c],pos); assert at>=prev, 'Unsorted records'
            if at!=prev: seen.clear()
            key=(ref,alt); assert key not in seen, 'Duplicate variant'; seen.add(key)
            prev=at; counts[c]+=1
        assert p.wait()==0, 'VCF decode failed'
    except BaseException:
        p.kill(); p.wait(); raise
    assert int(output([bc,'index','-n',path]).strip())==sum(counts.values())
    return dict(validated=True, variant_records=sum(counts.values()), contig_records=dict(counts), sample_count=len(samples), sample_order_verified=True)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--stage',choices=['setup','shard','assemble'],required=True)
    ap.add_argument('--index',type=int); ap.add_argument('--check-only',action='store_true'); a=ap.parse_args()
    config=json.loads((HERE/'config.json').read_text())
    frozen=json.loads((HERE/'FROZEN_SOURCES.json').read_text())
    for path,digest in frozen['sha256'].items(): assert sha(path)==digest, 'Frozen source changed: '+path
    cfg=json.loads(pathlib.Path(config['original_config']).read_text()); out=pathlib.Path(cfg['output_root'])
    dest=pathlib.Path(config['output_root']); dest.mkdir(parents=True,exist_ok=True)
    ref=out/'references/mongolicus/genome.fa'; fai=validate_layout(config,str(ref)+'.fai')
    bc=cfg['tools']['bcftools']; rows=list(csv.DictReader(open(cfg['manifest']),delimiter='\t'))
    samples=[r['sample_id'] for r in rows]; assert len(samples)==217 and len(set(samples))==217
    fingerprint=original_fingerprint(cfg,config)
    if a.check_only: print('Layout and frozen-source guards PASS; no calling performed'); return
    key=a.stage+(f'_{a.index:02d}' if a.stage=='shard' else '')
    marker=dest/(key+'.json')
    if marker.exists():
        m=json.loads(marker.read_text()); assert m['source_guard']==frozen; assert_products(m)
        print('Validated parallel checkpoint: '+key); return
    setup_path=dest/'setup.json'; bams=dest/'bam.list'; outputs=[]; extra={}
    if a.stage=='setup':
        for keyname in ['reference_mongolicus']+['map_mongolicus_'+s for s in samples]:
            m=json.loads((out/'checkpoints'/f'{keyname}.json').read_text())
            for field,value in fingerprint.items():
                if field=='stage' or (field=='manifest' and keyname.startswith('reference_')): continue
                assert m['fingerprint'][field]==value, 'Input checkpoint fingerprint changed'
            assert_products(m)
        bam_paths=[out/'mapping/mongolicus'/s/'aligned.bam' for s in samples]
        for path in bam_paths: assert pathlib.Path(str(path)+'.bai').is_file()
        bams.write_text(''.join(str(x)+'\n' for x in bam_paths))
        outputs=[bams]; extra={'validated_mapping_samples':217,'all_reference_contigs_partitioned_once':True,'total_contigs':len(fai)}
    else:
        setup=json.loads(setup_path.read_text()); assert setup['source_guard']==frozen; assert_products(setup)
        if a.stage=='shard':
            assert a.index is not None and 0<=a.index<len(config['shards'])
            shard=config['shards'][a.index]; sd=dest/f'shard_{a.index:02d}'; sd.mkdir(exist_ok=True)
            path=sd/'calls.vcf.gz'; partial=sd/'calls.partial.vcf.gz'
            assert not path.exists() and not partial.exists(), 'Uncheckpointed shard output; inspect before restart'
            cmd1=[bc,'mpileup','--threads',str(config['compression_threads']),'-f',str(ref),'-b',str(bams),'-q','20','-Q','20','-d',str(cfg['pileup_max_depth_per_sample']),'-a','FORMAT/DP,FORMAT/AD','--seed',str(config['seed']),'-R',shard['bed'],'-Ou']
            cmd2=[bc,'call','--threads',str(config['compression_threads']),'-m','-v','--ploidy','2','-Oz','-o',str(partial)]
            dump(sd/'commands.json',{'mpileup':cmd1,'call':cmd2,'joint_samples':217,'whole_contig_regions':shard['contigs'],'parameters_match_original':True,'compression_and_pipe_only_changed':True})
            print(json.dumps({'mpileup':cmd1,'call':cmd2}),flush=True)
            with open(sd/'mpileup.stderr.txt','w') as e1, open(sd/'call.stderr.txt','w') as e2:
                first=subprocess.Popen(cmd1,stdout=subprocess.PIPE,stderr=e1)
                try:
                    second=subprocess.Popen(cmd2,stdin=first.stdout,stderr=e2)
                    first.stdout.close(); rc2=second.wait(); rc1=first.wait()
                    assert rc1==0 and rc2==0, f'mpileup/call failed: {rc1}/{rc2}'
                except BaseException:
                    first.kill(); first.wait(); raise
            command([bc,'index','-t',partial])
            v=validate_vcf(bc,partial,samples,set(c for c,n in shard['contigs']),fai)
            partial.replace(path); pathlib.Path(str(partial)+'.tbi').replace(str(path)+'.tbi')
            dump(sd/'validation.json',v)
            outputs=[path,pathlib.Path(str(path)+'.tbi'),sd/'validation.json',sd/'commands.json',sd/'mpileup.stderr.txt',sd/'call.stderr.txt']; extra=v
        else:
            if config.get('portable_fresh_execution'):
                assert config.get('original_call_job') is None
                retirement={'original_job': None, 'confirmed_terminal': True, 'reason': 'Fresh reproduction; no monolithic caller was submitted'}
            else:
                retirement=json.loads(pathlib.Path(config['old_call_retirement_record']).read_text())
                assert retirement['original_job']==config['original_call_job'] and retirement['confirmed_terminal'] is True
            assert retirement['use_parallel_replacement'] is True
            assert not (out/'checkpoints/call.json').exists(), 'Existing completed original call; do not replace'
            paths=[]; total=0
            for shard in config['shards']:
                m=json.loads((dest/f"shard_{shard['index']:02d}.json").read_text())
                assert m['source_guard']==frozen; assert_products(m)
                paths.append(dest/f"shard_{shard['index']:02d}/calls.vcf.gz"); total+=m['details']['variant_records']
            listing=dest/'ordered_shard_vcfs.list'; listing.write_text(''.join(str(x)+'\n' for x in paths))
            merged=dest/'assembled.vcf.gz'; assert not merged.exists(), 'Uncheckpointed assembly exists'
            command([bc,'concat','--threads','4','-f',listing,'-Oz','-o',merged]); command([bc,'index','-t',merged])
            v=validate_vcf(bc,merged,samples,set(c for c,n in fai),fai); assert v['variant_records']==total and total>0
            # Read-only REF check; normalized variants are discarded, never substituted.
            command([bc,'norm','-c','e','-f',ref,'-Ou','-o','/dev/null',merged])
            raw=out/'variants/raw.vcf.gz'; assert not raw.exists(), 'Canonical output exists; do not overwrite'
            merged.replace(raw); pathlib.Path(str(merged)+'.tbi').replace(str(raw)+'.tbi')
            canonical_bams=out/'variants/bam.list'
            if canonical_bams.exists(): assert canonical_bams.read_bytes()==bams.read_bytes()
            else: canonical_bams.write_bytes(bams.read_bytes())
            provenance=dest/'ASSEMBLY_VALIDATION.json'
            dump(provenance,dict(v,shards=len(paths),all_contigs_partitioned_once=True,reference_contigs=len(fai),same_joint_sample_order=True,original_call_job_retired=retirement,implementation='Full-contig parallel bcftools mpileup/call, ordered concat; original worker/config not modified',execution_source_sha256=sha(__file__),source_guard=frozen,bitwise_identity_to_monolithic_run_not_claimed=True))
            products=[raw,pathlib.Path(str(raw)+'.tbi'),canonical_bams,provenance,HERE/'config.json',HERE/'FROZEN_SOURCES.json']
            # Compatible input fingerprint for the unmodified downstream filter. Truthful
            # adapter identity is explicit; this does not claim original worker execution.
            dump(out/'checkpoints/call.json',dict(validated=True,stage='call',fingerprint=fingerprint,outputs={str(p.resolve()):sha(p) for p in products},completed_unix=__import__('time').time(),scientific_interpretation_validated=False,execution_adapter=dict(path=str(HERE/'call_shards.py'),sha256=sha(__file__),implementation='parallel whole-contig joint calling',original_single_call_worker_executed=False),validation=str(provenance)))
            outputs=products+[listing]; extra=v
    dump(marker,dict(validated=True,stage=key,source_guard=frozen,outputs={str(p.resolve()):sha(p) for p in outputs},details=extra,completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))

if __name__=='__main__':
    try: main()
    except BaseException:
        print(traceback.format_exc(),file=sys.stderr); raise
