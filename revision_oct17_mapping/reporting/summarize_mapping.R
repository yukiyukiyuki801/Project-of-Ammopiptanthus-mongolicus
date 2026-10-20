#!/usr/bin/env Rscript
# Descriptive reporting only; no tests, additional models, or PCA refitting.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==4)
config_path <- normalizePath(args[1]); metadata_path <- normalizePath(args[2])
crosswalk_path <- normalizePath(args[3]); destination <- args[4]
library(jsonlite)
cfg <- fromJSON(config_path); root <- cfg$output_root
needed <- c(file.path(root,'summary/paired_mapping_depth.csv'),
 file.path(root,'variants/filter_counts.json'),file.path(root,'pca/mongolicus_PCA_scores.csv'),
 file.path(root,'pca/mongolicus_PCA_variance.csv'),file.path(root,'pca/new_arm_genotype_missingness.csv'),
 file.path(root,'pca/validation.json'),file.path(root,'preflight.json'),
 cfg$existing_pca_scores,cfg$existing_pca_variance,config_path,metadata_path,crosswalk_path,
 cfg$reference_audit)
stopifnot(all(file.exists(needed)),cfg$execution_frozen)
# Validate existing stage success markers; do not generate incomplete result summaries.
for (s in c('filter','pca','summarize')) {
 marker <- fromJSON(file.path(root,'checkpoints',paste0(s,'.json')))
 stopifnot(isTRUE(marker$validated))
}
qc <- read.csv(needed[1],stringsAsFactors=FALSE,check.names=FALSE)
meta <- read.csv(metadata_path,stringsAsFactors=FALSE)
cross <- read.delim(crosswalk_path,stringsAsFactors=FALSE)
old <- read.csv(cfg$existing_pca_scores,stringsAsFactors=FALSE)
new <- read.csv(needed[3],stringsAsFactors=FALSE)
oldvar <- read.csv(cfg$existing_pca_variance); newvar <- read.csv(needed[4])
missing <- read.csv(needed[5]); counts <- fromJSON(needed[2]); pca_validation <- fromJSON(needed[6])
stopifnot(nrow(qc)==434,nrow(meta)==217,nrow(cross)==217,nrow(old)==217,nrow(new)==217,
 !anyDuplicated(meta$sample_id),!anyDuplicated(cross$sample_id),!anyDuplicated(old$sample_id),
 !anyDuplicated(new$sample_id),!anyDuplicated(missing$sample_id),nrow(missing)==217,
 setequal(meta$sample_id,old$sample_id),setequal(meta$sample_id,new$sample_id),
 setequal(meta$sample_id,missing$sample_id),setequal(meta$sample_id,cross$sample_id),
 identical(sort(unique(qc$reference)),sort(c('nanus','mongolicus'))),
 !anyDuplicated(paste(qc$sample_id,qc$reference)),
 all(c('sample_id','population','class','reference') %in% names(qc)))
stopifnot(identical(as.integer(table(meta$group)[c('core','peripheral')]),c(104L,113L)))
cm <- meta[match(cross$sample_id,meta$sample_id),]
stopifnot(all(cross$group==cm$group),all(cross$population==cm$population),
 all(cross$friendly_id==cm$friendly_id),!anyDuplicated(cross$run_accession),
 isTRUE(pca_validation$validated))
for (arm in c('nanus','mongolicus')) {
 q <- qc[qc$reference==arm,]; m <- meta[match(q$sample_id,meta$sample_id),]
 stopifnot(nrow(q)==217,setequal(q$sample_id,meta$sample_id),
  all(q$population==m$population),all(q$class==m$group))
}
for (scores in list(old,new)) {
 m <- meta[match(scores$sample_id,meta$sample_id),]
 stopifnot(all(scores$group==m$group),all(is.finite(scores$PC1)),all(is.finite(scores$PC2)))
}
for (v in list(oldvar,newvar)) stopifnot(all(is.finite(v$variance_fraction)),
 all(v$variance_fraction>=0),sum(v$variance_fraction)<=1+1e-7,nrow(v)>=2)
stopifnot(counts$raw>=counts$site_filtered,counts$site_filtered>=counts$ld_pruned,
 counts$ld_pruned>=pca_validation$retained_loci,pca_validation$retained_loci>=2,
 all(is.finite(missing$missing_fraction)),all(missing$missing_fraction>=0 & missing$missing_fraction<=1))
metrics <- c('clean_pairs','primary_records','mapped_primary','mapped_fraction',
 'properly_paired_primary_fraction','mapq20_primary_fraction','covered_bases_q20Q20',
 'covered_bases_ge4','covered_base_mean_depth','covered_base_median_depth_lower',
 'assembly_bases','assembly_breadth_ge1','assembly_breadth_ge4','whole_assembly_mean_depth')
stopifnot(all(metrics %in% names(qc)))
for (nm in metrics) stopifnot(is.numeric(qc[[nm]]),all(is.finite(qc[[nm]])),all(qc[[nm]]>=0))
for (nm in grep('fraction|breadth',metrics,value=TRUE)) stopifnot(all(qc[[nm]]<=1))
stopifnot(all(qc$primary_records==2*qc$clean_pairs),
 all(abs(qc$mapped_fraction-qc$mapped_primary/qc$primary_records)<1e-12),
 all(qc$covered_bases_ge4<=qc$covered_bases_q20Q20),
 all(abs(qc$covered_base_mean_depth*qc$covered_bases_q20Q20-
  qc$whole_assembly_mean_depth*qc$assembly_bases)<1e-6*pmax(1,qc$covered_base_mean_depth*qc$covered_bases_q20Q20)))
qn <- qc[qc$reference=='nanus',]; qm <- qc[qc$reference=='mongolicus',]
qm <- qm[match(qn$sample_id,qm$sample_id),]
stopifnot(identical(qn$clean_pairs,qm$clean_pairs))
idx <- match(qc$sample_id,meta$sample_id); j <- match(qc$sample_id,cross$sample_id)
qc$friendly_id <- meta$friendly_id[idx]; qc$run_accession <- cross$run_accession[j]
qc$biosample <- cross$biosample[j]
qc <- qc[order(match(qc$sample_id,meta$sample_id),match(qc$reference,c('nanus','mongolicus'))),]
qc <- qc[,c('sample_id','friendly_id','population','class','run_accession','biosample',
 'reference',metrics,'depth_denominator_warning')]
summarize_numeric <- function(d,arm,group,level) do.call(rbind,lapply(metrics,function(nm){
 x <- d[[nm]];data.frame(reference=arm,stratum=group,stratum_level=level,metric=nm,
 n_individuals=length(x),mean=mean(x),sd=sd(x),median=median(x),min=min(x),max=max(x))
}))
descriptive <- do.call(rbind,lapply(c('nanus','mongolicus'),function(arm){
 d <- qc[qc$reference==arm,]
 rbind(summarize_numeric(d,arm,'overall','all_individuals'),
  summarize_numeric(d[d$class=='core',],arm,'core','population_class'),
  summarize_numeric(d[d$class=='peripheral',],arm,'peripheral','population_class'))
}))
# One paired descriptive change per individual, preserving the matched design.
# No p values or inference about independent reference assemblies are produced.
delta_metrics <- c('mapped_fraction','properly_paired_primary_fraction','mapq20_primary_fraction',
 'covered_base_mean_depth','whole_assembly_mean_depth','assembly_breadth_ge1','assembly_breadth_ge4')
delta <- data.frame(sample_id=qn$sample_id,class=qn$class)
for (nm in delta_metrics) delta[[paste0(nm,'_mongolicus_minus_nanus')]] <- qm[[nm]]-qn[[nm]]
delta_summary <- do.call(rbind,lapply(c('overall','core','peripheral'),function(g){
 d <- if(g=='overall')delta else delta[delta$class==g,]
 do.call(rbind,lapply(names(d)[-(1:2)],function(nm){x<-d[[nm]];data.frame(stratum=g,
  metric=nm,n_individuals=length(x),mean=mean(x),sd=sd(x),median=median(x),min=min(x),max=max(x))}))
}))
variance <- rbind(data.frame(reference='nanus',oldvar),data.frame(reference='mongolicus',newvar))
definitions <- data.frame(metric=c(metrics,'new_arm_genotype_missing_fraction'),
 definition=c('Number of paired reads retained by shared fastp cleaning.',
 'Primary alignment records excluding secondary and supplementary records; unmapped records included.',
 'Mapped primary alignment records.','Mapped primary records / all primary records.',
 'Properly paired mapped primary records / all primary records.',
 'Mapped primary records with MAPQ >=20 / all primary records.',
 'Number of assembly bases observed at depth >=1 after MAPQ >=20 and base quality >=20 filters.',
 'Number of assembly bases observed at depth >=4 under the same filters.',
 'Mean retained depth over bases with nonzero retained coverage; observed bases can differ between assemblies.',
 'Lower median retained depth over bases with nonzero retained coverage.',
 'Total reference assembly length in bp.','Observed depth >=1 bases / complete assembly length.',
 'Observed depth >=4 bases / complete assembly length.',
 'Sum of retained per-base depths / complete assembly length; uncovered bases contribute zero.',
 'Fraction of missing genotypes per individual in the new LD-pruned PCA panel, before mean imputation.'),
 units=c(rep('records_or_pairs',3),rep('fraction',3),rep('bp',2),rep('read_depth',2),'bp',
 rep('fraction',2),'read_depth','fraction'))
preflight <- fromJSON(file.path(root,'preflight.json'))
bwa_raw <- suppressWarnings(system2(cfg$tools$bwa,stdout=TRUE,stderr=TRUE))
versions <- preflight$versions_raw; versions$bwa_actual_invocation <- paste(bwa_raw,collapse='\n')
notes <- c('Descriptive statistics summarize 217 matched individuals (104 core and 113 peripheral) in each reference arm. SD is the sample standard deviation among individuals, not among independent populations.',
 'Class means weight individuals equally; no population-level independent-replicate test is implied.',
 'Depth and breadth denominators refer to each complete assembly or its observed covered bases; they are not matched homologous RAD regions and are not genotype missingness.',
 'The two reference arms use the same cleaned reads. Existing A. nanus PCA scores were reused; new A. mongolicus variants were called with bcftools, not the historical two-caller intersection.',
 'No new inferential tests, PCA alignment, Procrustes analysis, functional-load or genomic-offset reruns were performed by this reporting script.',
 'PCA agreement must be assessed separately from the plotted point distributions; these summaries do not automatically establish robustness.',
 'Fractions are stored unrounded (0-1); multiply by 100 for percentages. Paired differences of fractions become percentage points when multiplied by 100.')
sha256 <- function(p) {
 z <- system2('sha256sum',shQuote(normalizePath(p)),stdout=TRUE)
 stopifnot(length(z)==1);sub(' .*','',z)
}
source_files <- unique(c(needed,file.path(root,'checkpoints',paste0(c('filter','pca','summarize'),'.json'))))
provenance <- data.frame(path=source_files,sha256=vapply(source_files,sha256,''))
report <- list(validated=TRUE,completed_utc=format(Sys.time(),tz='UTC',usetz=TRUE),
 biological_result_status='DESCRIPTIVE_RESULTS_ONLY; PCA visual and manuscript interpretation pending',
 n_individuals=217,n_populations=length(unique(meta$population)),n_core=104,n_peripheral=113,
 mapping_records=nrow(qc),variant_counts=counts,new_pca_retained_loci=pca_validation$retained_loci,
 mapping_summaries=descriptive,paired_change_summaries=delta_summary,pca_variance=variance,
 new_arm_missingness=list(mean=mean(missing$missing_fraction),sd=sd(missing$missing_fraction),
 median=median(missing$missing_fraction),min=min(missing$missing_fraction),max=max(missing$missing_fraction)),
 inference_performed=FALSE,notes=notes,settings=cfg,tool_versions=versions,source_manifest=provenance)
stopifnot(!dir.exists(destination) || length(list.files(destination,all.files=FALSE))==0)
dir.create(destination,recursive=TRUE,showWarnings=FALSE)
tables <- list(sample_QC=qc,class_overall_summary=descriptive,paired_changes=delta,
 paired_change_summary=delta_summary,variant_counts=data.frame(stage=names(counts),sites=unlist(counts)),
 PCA_variance=variance,new_genotype_missing=missing,metric_definitions=definitions,
 reporting_notes=data.frame(note=notes),provenance=provenance)
for(nm in names(tables)) write.csv(tables[[nm]],file.path(destination,paste0(nm,'.csv')),row.names=FALSE,na='')
write_json(report,file.path(destination,'NUMERICAL_SUMMARY.json'),pretty=TRUE,auto_unbox=TRUE,digits=NA,na='null')
write_json(versions,file.path(destination,'TOOL_VERSIONS.json'),pretty=TRUE,auto_unbox=TRUE)
writeLines(capture.output(sessionInfo()),file.path(destination,'R_sessionInfo.txt'))
cat('Reporting CSV and JSON complete; package workbook with export_workbook.py. No interpretation or final figure QA claimed.\n')
