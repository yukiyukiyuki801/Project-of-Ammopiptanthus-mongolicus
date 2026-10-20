#!/usr/bin/env bash
set -euo pipefail
OUT_DIR="$PWD/runs/revision_run_01/outputs/functional_reconstruction"
mkdir -p "$OUT_DIR/sift_predictions"
sha256sum public_dependencies/uniprot_sprot.fasta "$OUT_DIR/target_proteins.fasta" > "$OUT_DIR/sift_input_sha256.txt"
public_dependencies/software/sift4g/bin/sift4g -q "$OUT_DIR/target_proteins.fasta" -d "$PWD/public_dependencies/uniprot_sprot.fasta" --subst "$OUT_DIR/substitutions" --out "$OUT_DIR/sift_predictions" -t "${AMMO_THREADS:-8}"
"${AMMO_PYTHON:-python3}" - <<'PY'
import json,pathlib,datetime
p=pathlib.Path('runs/revision_run_01/outputs/functional_reconstruction')
expected=sum(1 for l in (p/'target_proteins.fasta').open() if l.startswith('>'))
files=list((p/'sift_predictions').glob('*.SIFTprediction'))
assert len(files)==expected,(len(files),expected)
assert all(f.stat().st_size>0 for f in files)
(p/'sift_execution.json').write_text(json.dumps({'validated':True,'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'expected_proteins':expected,'prediction_files':len(files),'status':'Execution complete; score parsing and scientific validation still required','database_provenance':'public_dependencies/DOWNLOAD_PROVENANCE.txt; downloaded2026-10-13; new reconstruction, not historical database'},indent=2)+'\n')
PY
