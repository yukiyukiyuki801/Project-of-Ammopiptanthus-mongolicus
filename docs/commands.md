# Numerical commands and dependencies

Run from the repository root after restoring data and sourcing `tools/activate.sh`. Commands below are individual entry points; they are not an instruction to rerun every completed stage. Use fresh output locations or preserve supplied products before numerical refitting.

| Analysis | Command | Upstream inputs |
|---|---|---|
| Input validation | `python3 scripts/revision/audit_inputs.py` | Exact VCFs, sample/environment manifests, population classes and supplied imputed matrix |
| Climate | `Rscript scripts/revision/extract_climate.R` | Validated samples, current/ensemble crops and individual-GCM rasters |
| Diversity | `python3 scripts/revision/calculate_diversity.py` | Full genotype matrix and locus/sample keys |
| LFMM | `Rscript scripts/revision/fit_lfmm.R` | Candidate genotypes and climates |
| RONA | `Rscript scripts/revision/calculate_rona.R` | LFMM output, candidate genotypes and climates; pinned pyRona for the numerical oracle |
| Imputation sensitivity | `Rscript scripts/revision/imputation_sensitivity.R` | Called/imputed candidate genotypes and climates |
| GDM fits | `Rscript scripts/revision/fit_gdm.R` | Candidate FST and current/future population climates |
| Full landscape matching | `Rscript scripts/revision/landscape_search.R` | GDM fits, regional climate and habitat geometry |
| Matching summaries | `Rscript scripts/revision/summarize_landscape_populations.R` | Completed matching tables |
| Population tests | `Rscript scripts/revision/population_statistics.R` | RONA, GDM and diversity outputs |
| Classification sensitivity | `Rscript scripts/revision/classification_sensitivity.R` | Populations, occurrences, matching outputs and habitat polygons |
| Distance limits | `Rscript scripts/revision/oct16_distance_classification.R` | Existing full-domain offsets |
| Functional consequences | `python3 scripts/revision/reconstruct_functional_consequences.py` | Full variant VCF, exact reference, GTF and protein FASTA |
| SIFT | `AMMO_THREADS=8 bash scripts/revision/run_functional_sift.sh` | Target proteins, substitutions and recorded Swiss-Prot database; compiled SIFT4G executable |
| Functional score parsing | `python3 scripts/revision/summarize_functional_scores.py` | Completed SIFTprediction files and consequence/dosage tables |
| Population functional ratios | `Rscript scripts/revision/summarize_functional_populations.R` | Parsed variant consequences and sample/population metadata |
| FST panel comparison | `Rscript scripts/revision/compare_fst_subsets.R` | Full/candidate/noncandidate genotypes and supplied differentiation matrix |
| Niche tests | `Rscript scripts/revision/reconstruct_niche.R` | Current climate, populations, and verified ecospat source functions |
| Species distribution | `Rscript scripts/revision/reconstruct_sdm.R` | Climate, occurrences, elevation, MaxEnt Java/dismo |

Driver-model commands, correlation screening and nested validation are documented in `analyses/driver_models/README.md`. Two-reference mapping, ten-shard joint calling, filtering, PCA and folded-SFS length commands are in `README_MAPPING.md`. Those upstream workflows require suitable CPU, RAM and temporary disk; no scheduler jobs are submitted by restoration or configuration.

The supplied results, fitted models and plot inputs allow inspection and redrawing without repeating expensive fits. Historical demographic summaries and EEMS are retained inputs, not inference reruns. For figure reconstruction use the current sources and commands in `docs/figures.md`, rather than a superseded renderer.
