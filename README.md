# Ammopiptanthus genomic vulnerability under climate change

This repository contains scripts used to reproduce the main climate-change genomic offset analyses for:

**Peripheral populations buffer genomic vulnerability of core populations to climate change in *Ammopiptanthus mongolicus***

## Contents

| Script | Purpose |
|---|---|
| `00_config.R` | Shared paths, scenario names, climate variables, and package checks. |
| `01_prepare_lfmm2_for_pyrona.R` | Converts the adaptive SNP VCF to LFMM format, imputes missing genotypes, runs LFMM2, and writes `associations.csv` for pyRona. |
| `01_run_pyrona.sh` | Runs pyRona for SSP126/SSP585 and 2041-2070/2071-2100. |
| `02_rona_maps.R` | Plots population-level RONA values over climate-change rasters. |
| `03_sdm_ensemble.R` | Fits ensemble species distribution models and predicts current/future habitat suitability. |
| `04_gdm_offsets.R` | Fits GDM and estimates local, forward, and reverse genomic offsets. |
| `05_core_periphery_offset_tests.R` | Compares local, forward, and reverse offsets between core and peripheral regions. |
| `06_gene_flow_source_summary.R` | Classifies GDM best-matching source regions as core, peripheral, or other. |

## Input data

The scripts expect an `ammo` data directory containing the derived SNP, climate, range, RONA, and GDM input files. By default:

```bash
/home/liao/DataVol1/personal2/sun/ammo
```

To use another location:

```bash
export AMMO_DIR=/path/to/ammo
```

Required files include:

| File or folder | Description |
|---|---|
| `ammo/ammo_outliers_3overlap.recode.vcf` | Adaptive SNP VCF. |
| `ammo/rona/ammo_current.txt` | Current environmental covariates for pyRona/LFMM2. |
| `ammo/rona/ammo_rona_ssp126_2040-2070.txt` | Future covariates for SSP126, 2041-2070. |
| `ammo/rona/ammo_rona_ssp126_2070-2100.txt` | Future covariates for SSP126, 2071-2100. |
| `ammo/rona/ammo_rona_ssp585_2040-2070.txt` | Future covariates for SSP585, 2041-2070. |
| `ammo/rona/ammo_rona_ssp585_2070-2100.txt` | Future covariates for SSP585, 2071-2100. |
| `ammo/rona/var_name.txt` | Covariate names for pyRona. |
| `ammo/gf/ammo_outlier_pairwise_fst.csv` | Pairwise adaptive genetic differentiation matrix for GDM. |
| `ammo/gf/ammo_xy_19pop.csv` | Population coordinates and core/peripheral labels. |
| `ammo/gf/ammo_binary.shp` | Species range mask. |
| `ammo/gf/ammo_core.shp` | Core-region mask. |
| `ammo/gf/ammo_peri.shp` | Peripheral-region mask. |
| `ammo/gf/chelsa2_crop_selected/` | Current and future CHELSA climate rasters. |
| `ammo/rona/*_rona_xy.csv` | Population-level RONA estimates for mapping. |
| `ammo/rona/*_diff.asc` | Climate-change raster differences for RONA maps. |

For species distribution modelling, provide the thinned occurrence shapefile:

```bash
export AMMO_OCCURRENCE_SHP=/path/to/ammo_thin_10km.shp
```

## Run order

Run from the repository root:

```bash
Rscript 01_prepare_lfmm2_for_pyrona.R
PYRONA=/path/to/pyRona.py bash 01_run_pyrona.sh
Rscript 02_rona_maps.R
Rscript 03_sdm_ensemble.R
Rscript 04_gdm_offsets.R
Rscript 05_core_periphery_offset_tests.R
Rscript 06_gene_flow_source_summary.R
```

`04_gdm_offsets.R` can be computationally intensive. Set the number of cores with:

```bash
AMMO_CORES=8 Rscript 04_gdm_offsets.R
```

Outputs are written to `outputs/` by default. To change this:

```bash
export AMMO_ANALYSIS_OUT=/path/to/output
```

## Main outputs

| Script | Outputs |
|---|---|
| `01_prepare_lfmm2_for_pyrona.R` | Imputed LFMM genotype file and `associations.csv`. |
| `01_run_pyrona.sh` | RONA result text files and pyRona PDFs. |
| `02_rona_maps.R` | RONA map PDFs. |
| `03_sdm_ensemble.R` | Current and future habitat suitability rasters. |
| `04_gdm_offsets.R` | Forward and reverse GDM offset tables. |
| `05_core_periphery_offset_tests.R` | Core/peripheral offset summaries and Wilcoxon tests. |
| `06_gene_flow_source_summary.R` | Source-region summary table and figure. |

## Software

R packages:

```text
LEA, nFactors, raster, sf, gdm, geosphere, foreach, doParallel, ggplot2, dplyr, sdm
```

External software:

```text
pyRona
```

## Data availability

RAD-seq data and derived SNP datasets are available at Mendeley Data: DOI `10.17632/9g8fd4t2nj.2`. Climate rasters were obtained from CHELSA.

