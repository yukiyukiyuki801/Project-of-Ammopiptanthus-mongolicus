# Population-level RONA driver models

This module releases the completed October 16 revision analysis of 19 Ammopiptanthus mongolicus populations, 20 fixed RONA responses and 11 potential drivers. It contains the final correlation-screened production fitter, its postprocessor, numerical inputs, complete aggregate result tables and compact selected model objects. No models were fitted when preparing this release.

## Inputs and interpretation

`data/model_inputs_380.csv` is the exact production input: 20 cases × 19 populations. Cases combine SSP126/SSP585, 2041–2070/2071–2100 and bio4/bio5/bio6/bio15/bio17. `data/population_case_saved_inputs_380.csv` supplies the audited precursor rows; `data/geographic_blocks.csv` preserves the previously fixed three-block partition.

The response is previously estimated RONA. Potential drivers, in input-column order, are climate change, private alleles, Ho, He, Fis, Pi, longitude, latitude, altitude, DEL/SYN and LoF/SYN. Ho and unbiased He were recomputed from the validated full SNP panel. Private alleles, Fis and Pi use the original saved full-panel statistics; they are not a wholly recomputed diversity panel. Pi is variant-panel diversity, and private counts are unrarefied. DEL/SYN and LoF/SYN are reference-relative functional indices, not measured fitness load or independently polarized deleterious alleles.

The correlation screen is response-independent and repeated within each training fold. Retention priority is `climate_change, Ho, He, DEL, LoF, private_alleles, Fis, Alt, Lat, Lon, Pi`. A predictor is rejected when its absolute Spearman correlation with an already retained predictor exceeds 0.8. This gives original variables priority over added Pi and records the actual competing variable and signed correlation. It does not establish causal independence. The full-data screen retains ten drivers in 17 cases and nine in three cases: Pi is excluded in all cases and latitude additionally in three. Screened-out variables receive NA inclusion weights, not evidence of biological unimportance.

All 2^11 = 2,048 additive subsets remain in the traceability universe for each model family; only subsets eligible after screening can compete. No variable, including climate, is mandatory. The null model is included. Families are Gaussian linear regression and Gaussian penalized GAMs; each included GAM predictor receives a thin-plate smooth with k=3, ML estimation and `select=TRUE`. Ineligible fits remain explicitly flagged. No interactions are added.

AICc uses the effective parameter count from `logLik`, including Gaussian scale. GAM AICc is an approximate small-sample criterion. Full-data ranks and within-family Akaike weights describe conditional support. Parameter intervals and smooth-test p values condition on selection and omit model-search uncertainty.

Primary validation repeats screening, training-only scaling, subset search and training-AICc selection inside each outer leave-one-population-out or fixed geographic-block fold. Joint selection compares the two families' training-AICc winners; it does not choose a family using held-out error. The mean-only baseline is recalculated from training data. Upstream candidate discovery and RONA estimation remain fixed: this is conditional downstream validation, not nested end-to-end validation.

## Current numerical results

The full-data joint search selects a GAM in 16 of 20 cases and a linear model in four. Joint nested selection beats the training-mean baseline in 15/20 LOPO and 14/20 geographic-block comparisons; GAM-only selection has the same counts, and linear-only selection improves on the baseline in 12/20 and 13/20, respectively. Individual-case magnitudes and exceptions are in `results/nested_selection_validation_performance.csv`. High in-sample fit should not be substituted for these held-out results.

`results/top1_combined_full_models.csv` is the full-data joint winner table. `results/nested_selection_predictor_stability.csv` reports fold selection conditional on eligibility; it is not causal importance. `results/top1_GAM_conditional_curve_source.csv` and `results/top1_GAM_partial_residual_source.csv` contain numerical sources for response-curve plots. Saved top-one models and standardizations are in `selected_models/<case>/`. Read an RDS model with `readRDS()` in a compatible R/mgcv environment; the model's predictors must receive its saved training means and SDs before prediction.

## Reproduction commands

Run from this module directory with Python 3 and R. Production used R 4.4.3, mgcv 1.9-4 and geosphere 1.6-8; `results/sessionInfo.txt` preserves the full environment, including its historical BLAS path. Python input preparation needs only the standard library.

Regenerate the modelling input from the supplied precursor and block lookup:

```sh
python scripts/prepare_model_inputs.py data/population_case_saved_inputs_380.csv rebuilt_inputs.csv data/geographic_blocks.csv
```

The resulting CSV must be byte-identical to `data/model_inputs_380.csv`; the original production SHA256 is `df13e415e90048d72e72b8b24f2e1c22216eb1859d2e2c4206c6800efc56f79b`.

The following commands explicitly request an expensive new fit. They are provided for reproducibility and were not run for this release. Each case is independent and can be scheduled separately. Use a new writable output directory; case indices use R's sorted case IDs.

```sh
Rscript --vanilla scripts/fit_all_subsets_production.R data/model_inputs_380.csv rerun_cases 1 both
# Repeat for case indices 2 through 20, or schedule one index per task.
Rscript --vanilla scripts/postprocess_driver_models.R rerun_cases rerun_summary 20
```

The fitter preserves completion/resume checks only when the saved input rows agree. Do not set `AMMO_SUBSET_BENCHMARK=1` for production. Postprocessing requires all 20 complete case outputs and computes diagnostics without refitting models; it will not finalize a partial set. It recalculates permutation Moran tests (1,999 permutations) and GAM k-checks (400 replicates) with fixed seeds.

For compact distribution, the three large aggregate CSVs are losslessly compressed as `*.csv.gz`. Read them directly with Python's `gzip` or R's `read.csv(gzfile(path))`. Decompressed bytes match the completed source tables. The complete production per-case top-ten cache is intentionally not distributed. `selected_models` retains one winner per family per case, not every file required by the postprocessor; to rerun the full postprocessor, first generate complete `rerun_cases` outputs with the fitter. Already computed postprocess results are provided directly.

## Provenance and validation limits

`SOURCE_MANIFEST.json` records original relative paths, original SHA256 hashes, release hashes and every adaptation. The two R scripts are byte-identical to the completed production sources; input preparation changes only the block-file path to a required CLI argument and adds usage/help. `provenance/original_model_inputs_manifest.json` records the historical predictor mapping and block assignment; its historical paths are provenance, not runnable dependencies.

This release was checked for Python syntax/CLI usage, R syntax, input byte-for-byte regeneration, table dimensions, completion markers and file hashes. No scientific fit was launched and no new statistical results were generated. A fully fresh environment installation and a complete end-to-end rerun were not performed. `RELEASE_VALIDATION.json` records the checks actually completed.
