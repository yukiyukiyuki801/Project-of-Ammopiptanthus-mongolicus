# Ammopiptanthus mongolicus genomic vulnerability analyses

Code and analysis data for **Distinct patterns of genomic vulnerability to climate change in core and peripheral populations of Ammopiptanthus mongolicus**.

Maintainer and code contributor: **Pei-Wei Sun** ([yukiyukiyuki801](https://github.com/yukiyukiyuki801)). This version replaces the earlier May2026 pipeline with the code, numerical inputs, results and figure sources used by the revised manuscript. The earlier version remains available in Git history.

## Start here

Clone the repository, install the dependencies listed under `environment/`, and restore the exact input/result files:

```bash
git clone https://github.com/yukiyukiyuki801/Project-of-Ammopiptanthus-mongolicus.git
cd Project-of-Ammopiptanthus-mongolicus
python3 tools/prepare_data.py --verify-only
python3 tools/prepare_data.py
source tools/activate.sh
```

Restoration checks SHA256 for every stored asset and every reconstructed file. It performs no statistical fitting, figure rendering, downloads, or job submission. Existing files with different contents are preserved and cause an error. The data store includes compressed text and bounded binary parts; `provenance/data_manifest.json` maps each asset to its exact analysis path. Data restoration requires additional disk space beyond the Git checkout.

Run commands from the repository root. Existing numerical results are included, so plots and summaries can be inspected without refitting every model. Execute only the analysis needed. No command automatically submits the whole workflow.

## Analyses

| Analysis | Code and entry points | Inputs and retained results after data restoration |
|---|---|---|
| Samples, SNPs and population classes | `scripts/revision/audit_inputs.py` | Exact 38,124-SNP and 1,637-candidate VCFs, 217-sample crosswalk, 19-population metadata; `runs/revision_run_01/outputs/inputs` |
| Climate extraction | `scripts/revision/extract_climate.R` | Current/ensemble CHELSA crops and five individual GCMs; `outputs/climate` |
| Genetic diversity | `scripts/revision/calculate_diversity.py`, `population_statistics.R` | Called genotypes and original supplied population summaries; `outputs/statistics` |
| Functional consequences and population ratios | `reconstruct_functional_consequences.py`, `run_functional_sift.sh`, `summarize_functional_scores.py`, `summarize_functional_populations.R` | Exact annotation/protein/reference inputs, SIFT substitutions/predictions, consequence tables, frequencies and population comparisons; `outputs/functional_reconstruction` |
| Pairwise differentiation | `scripts/revision/compare_fst_subsets.R` | Candidate/full/noncandidate SNP definitions and numerical FST comparison; `outputs/fst_subset_comparison` |
| Genetic PCA | `scripts/revision/calculate_genetic_pca.R` | Existing nanus-panel genotype matrix, scores and variance; `outputs/figures/genetic_pca` |
| Demographic histories | Latest `revision_oct18/figure_updates/render_S2.R`; mapped-length code in `revision_oct17_mapping/stairway_length` | Saved folded-SFS blueprints, summary trajectories and mapped-region inputs/results; no demographic refitting during this revision |
| RONA and LFMM | `fit_lfmm.R`, `calculate_rona.R`, `verify_rona_reference.py`, `imputation_sensitivity.R`, `population_statistics.R` | Called/imputed candidate genotypes, current/future environments, associations, RONA and sensitivity tables |
| RONA driver models | [analyses/driver_models](analyses/driver_models/README.md) | Exact 380-row input, final all-subset linear/GAM code, correlation screens, nested validation, aggregate results and 40 selected model objects |
| Climatic niche | `reconstruct_niche.R` | Shared background, population climates, density grids and all null-test results; `outputs/niche_reconstruction` |
| Species distribution | `reconstruct_sdm.R` | Occurrences, climate, elevation, MaxEnt settings, fitted models, suitability and habitat transitions; `outputs/sdm_reconstruction` |
| GDM and matching searches | `fit_gdm.R`, `landscape_search.R`, `exact_gdm_grid_search.cpp`, `summarize_landscape_populations.R`, `project_required_grid_offsets.R` | Candidate FST, climate, fitted models, full matching domain, population/grid offsets and endpoint summaries |
| Classifications and distance limits | `classification_sensitivity.R`, `oct16_distance_classification.R` | Published population classes, occurrence-based alternatives, and unrestricted/10/50/100/200km matching summaries |
| EEMS | `scripts/revision/extract_existing_eems.R` and saved-result figure sources | Retained original posterior summaries and plotting data; no new EEMS MCMC |
| Reference mapping and new-reference calling/PCA | [README_MAPPING.md](README_MAPPING.md) | 222 public runs, verified 217-individual subset, exact two-reference workflow, ten whole-contig calling shards, filtering, QC, PCA and mapped-length measurements |
| Figures and supplementary tables | [docs/figures.md](docs/figures.md), `tools/render_saved_figure.R` | All6 main and20 supplementary figures in PDF/PNG, exact source data and renderers, final28-table workbook |

Unless a full path is shown, numerical sources in this table are under `scripts/revision/`, and `outputs/...` means `runs/revision_run_01/outputs/...`. The final GAM module supersedes the earlier fixed-model diagnostics; obsolete models are not presented as the current analysis.

## Software and execution

R package versions used are recorded in `environment/R_versions_used.json` and retained `sessionInfo.txt` files. `environment/packages-R.txt` and `requirements-python.txt` identify dependencies; these lists do not automatically install them. The final GAM module also provides its own environment specification. Mapping additionally requires fastp, BWA, SAMtools, BCFtools, VCFtools and PLINK; the exact settings and version provenance are in the mapping module. MaxEnt requires Java/dismo. SIFT4G requires its compiled executable and the recorded Swiss-Prot reference protein database.

Fetch selected pinned upstream sources only when needed:

```bash
python3 tools/fetch_software.py --pyrona --ecospat-functions --hierfstat-functions
```

These fetches retain upstream software authorship. The ecospat functions are verified against the exact Git blob and SHA256 used; the three hierfstat oracle functions are pinned by commit and SHA256; pyRona is checked out at the recorded commit. Other tools and databases are obtained from their own providers using `provenance/upstream_software.json`, the retained download checksums and mapping documentation.

[docs/commands.md](docs/commands.md) gives the numerical execution order. The final driver-model module and mapping module have their own command documentation. Most plotting code operates on saved results. Optional local figure-QA helpers are not bundled; absent checks are explicitly `NOT_RUN`, rather than reported as successful.

Validation for this release covers file checksums, restoration, code syntax, identities, table/figure provenance and targeted configuration/input checks. Preparing this release did not refit statistical models. Complete fresh-environment reproduction of every analysis has not been performed.

## Data sources

The original RAD reads are in NCBI SRA **PRJNA1545357 / SRP746694**. The mapping module records all222 runs and the verified217-sample crosswalk; raw FASTQs and BAMs are obtained or regenerated externally. The previously published filtered SNP deposit is [Mendeley Data, 10.17632/9g8fd4t2nj.2](https://doi.org/10.17632/9g8fd4t2nj.2). The exact variant-only SNP inputs used here are included separately and are not asserted to be identical to that filtered release. No invariant sites are inferred from these variant VCFs.

Climate inputs come from CHELSA2.1, and elevation from WorldClim2.1. The exact cropped analysis inputs are supplied with provenance. The exact A.nanus reference bytes associated with the original coordinates are included in bounded parts because an independent public URL for those exact bytes was not established; restoration verifies the recorded SHA256. The A.mongolicus AMO v1.0 reference is retrieved from its publisher using the URL/checksum in the mapping module. Do not substitute the similarly named NCBI assembly without checking sequence identity. Original source/data attribution remains applicable; the MIT license covers this repository's code and does not replace third-party data or software terms.

The repository excludes manuscript drafts, reviewer correspondence, credentials, installed environments, raw sequencing downloads and BAMs. `provenance/` records source-to-release mappings, data checksums and software sources. See `CITATION.cff` for the repository software citation.
