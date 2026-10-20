# Reference-mapping and folded-SFS sequence-length reproduction

This module contains the executed two-reference mapping workflow, the whole-contig parallel caller that actually produced the new-reference variants, read retrieval, descriptive reporting, and bounded mapped-region measurement code. No analysis is run by installation or configuration. Use fresh output directories; old machine-specific checkpoint hashes cannot be reused after relocation.

The public reads are PRJNA1545357 / SRP746694. `revision_oct17_mapping/read_access/download/runs.txt` lists 222 runs. The 217-individual crosswalk and five excluded runs are separate TSVs. `original_fastq_checksums.tsv` records the submitted original FASTQ names, sizes and MD5 values. The downloader requests fresh NCBI locator URLs and downloads full-quality paired FASTQs. It does not use quality-reduced SRA Lite. No expiring signed URLs, FASTQs or BAMs are bundled.

## References and requirements

`references/PUBLIC_REFERENCE_SOURCES.json` gives the exact compressed-file SHA256/MD5 and the permanent publisher download URL for A. mongolicus AMO v1.0. The publisher assembly has 498 sequences/843,074,137 bases and differs from NCBI GCA_039655645.1, which cannot be silently substituted. The exact bundled A. nanus compressed file is SHA256 `5d99d8b199bc23c54400199e62c44c9efb62abddbdb5c21cea73444cc8a5f7f2` (1,099 sequences/823,736,225 bases). A verified public URL for that exact file was not recovered; supply the matching source file. This is an explicit remaining upstream-access dependency, not a missing analysis result.

Required executables on PATH: Python3, gzip, curl8.x, bwa0.7.17, samtools, fastp0.23.2, bcftools, vcftools, plink1.9 and Rscript. See the release environment/tool-version metadata for exact recorded versions. R PCA/reporting scripts require their declared R packages. `curl` can be overridden with the CURL environment variable. Portable copies change only paths/default executable lookup and the obsolete monolithic-job retirement guard; original executed Python sources are included under `provenance/executed_sources/` with source hashes.

## Configure a fresh working copy

Run from the repository root containing `revision_oct17_mapping/`:

```bash
python3 prepare_mapping.py --work-root /large_volume/ammo_reproduction --raw-fastq-root /large_volume/ammo_reproduction/raw_fastq --nanus-reference /references/ammo.fasta.gz --mongolicus-reference /references/AMO_genome_v1.0.fasta.gz
```

This writes local absolute paths from portable templates and recomputes source guards. It does not download, map, call or submit jobs. The local generated `config.frozen.json`, `config.json`, `FROZEN_SOURCES.json` and `samples.tsv` should remain untracked because paths are machine-specific. `--skip-reference-bytecheck` permits configuration before downloads; this does not waive subsequent worker checksum checks.

Retrieve any project run with:

```bash
python3 revision_oct17_mapping/read_access/download/download_fastq.py --run SRR41233426 --root /large_volume/ammo_reproduction/raw_fastq
```

Use `--index 0..221` for an externally managed download array; the script itself never submits jobs. Each paired FASTQ must pass size, MD5 and gzip CRC. Names and run directories must match the generated sample TSV.

## Execute individual validated stages

The workflow is deliberately stage-specific. Run preflight, independently index each reference, then clean each sample once and map the same cleaned pair against each reference. Do not use the old monolithic `--stage call` path for reproducing the published revision.

```bash
python3 revision_oct17_mapping/workflow/run_workflow.py --config revision_oct17_mapping/parent/config.frozen.json --stage preflight
python3 revision_oct17_mapping/workflow/run_workflow.py --config revision_oct17_mapping/parent/config.frozen.json --stage reference --reference nanus
python3 revision_oct17_mapping/workflow/run_workflow.py --config revision_oct17_mapping/parent/config.frozen.json --stage reference --reference mongolicus
python3 revision_oct17_mapping/workflow/run_workflow.py --config revision_oct17_mapping/parent/config.frozen.json --stage clean --sample aa
python3 revision_oct17_mapping/workflow/run_workflow.py --config revision_oct17_mapping/parent/config.frozen.json --stage map --sample aa --reference nanus
python3 revision_oct17_mapping/workflow/run_workflow.py --config revision_oct17_mapping/parent/config.frozen.json --stage map --sample aa --reference mongolicus
```

Repeat clean/map for all 217 sample IDs using the resource manager of your choice. Only after every mongolicus BAM checkpoint passes:

```bash
python3 revision_oct17_mapping/parallel_call/call_shards.py --stage setup
python3 revision_oct17_mapping/parallel_call/call_shards.py --stage shard --index 0
# Run all ten indices 0..9 independently, then assemble once all pass.
python3 revision_oct17_mapping/parallel_call/call_shards.py --stage assemble
python3 revision_oct17_mapping/workflow/run_workflow.py --config revision_oct17_mapping/parent/config.frozen.json --stage filter
python3 revision_oct17_mapping/workflow/run_workflow.py --config revision_oct17_mapping/parent/config.frozen.json --stage pca
python3 revision_oct17_mapping/workflow/run_workflow.py --config revision_oct17_mapping/parent/config.frozen.json --stage summarize
```

Nine chromosome shards plus a tenth shard of remaining whole scaffolds cover all 498 reference sequences exactly once. Calling remains joint across the same217 diploid samples, using bcftools mpileup MAPQ20/BQ20/depthcap10000, FORMAT/DP+AD, seed0, then multiallelicvariant calling. Genotypes DP<4/missingDP are masked; INFO allelecounts are refreshed before biallelic-SNP completeness/MAC/MAF/QUAL/LD filtering. Exact filtering flags are in `workflow/run_workflow.py`; the final count is checkpointed rather than forced to equal the historical nanus panel. Nanus PCA is reused from the existing38,124-SNP data; no new nanus SNP calling is part of this sensitivity.

Run reporting only after these stages validate:

```bash
Rscript revision_oct17_mapping/reporting/summarize_mapping.R revision_oct17_mapping/parent/config.frozen.json revision_oct17_mapping/inputs/sample_metadata.csv revision_oct17_mapping/read_access/crosswalk/crosswalk217.tsv results/mapping
python3 revision_oct17_mapping/reporting/export_workbook.py results/mapping
```

Mapping/PCA similarity supports broad population-structure sensitivity only. It does not validate cross-reference functional-load or offset equivalence, ancestry reconstruction or genotype callability at unrepresented invariant sites.

## Folded-SFS mapped length and retained histories

Stairway histories were reused from saved folded-SFS results, not refitted in this revision. The inventory expands exact final six summary/blueprint paths and saved700-deme EEMS inputs. No EEMS MCMC should be launched. Mapped-length code measures the union of positive MAPQ20 aligned M/=/X blocks containing the specified SFS SNPs, with overlapping reference positions counted once. It does not define invariant-site callability or estimate new ancestral states. The explicit historical length arithmetic/projection audit and later author confirmation are retained separately. `stairway_length/code/config.template.json` controls per-sample measurement and aggregation; inspect source registry dependencies before executing. Historical easySFS data dictionaries/blueprints and reference-matched target BEDs are required to regenerate input membership from first principles; compact targets and verified membership registry are included.


## Exact nanus reference supplied in this repository

After `python3 tools/prepare_data.py`, the exact compressed reference is restored at `ammo_analysis_package_20261014/cluster_recovery/personal2/sun/sift4g/ammo/chr-src/ammo.fasta.gz` (SHA256 `5d99d8b199bc23c54400199e62c44c9efb62abddbdb5c21cea73444cc8a5f7f2`). It is stored in bounded binary parts with complete SHA256 verification, rather than assigning an unverified assembly accession. Use that restored path for `--nanus-reference`. The absence of a separate verified public original URL does not prevent access to these exact released bytes.
