# Numerical PCA of the exact retained SNP matrix; no figure rendering or clustering.
source("scripts/revision/common.R")
directory <- file.path(output_dir,"figures","genetic_pca")
dir.create(directory,recursive=TRUE,showWarnings=FALSE)
qa_python <- Sys.getenv("AMMO_PYTHON", unset="python3")
# Python is used solely to convert an already validated numpy matrix to text.
conversion <- paste("import numpy as np; from pathlib import Path;",
  "p=Path('runs/revision_run_01/outputs'); a=np.load(p/'inputs/full_genotypes.npy');",
  "assert a.shape==(217,38124); np.savetxt(p/'figures/genetic_pca/full_genotypes_conversion.tsv',a,fmt='%d',delimiter='\\t')")
stopifnot(system2(qa_python,c("-c",shQuote(conversion)))==0)
g <- as.matrix(read.table(file.path(directory,"full_genotypes_conversion.tsv"),header=FALSE))
stopifnot(identical(dim(g),c(217L,38124L)),all(g %in% c(0,1,2,9)))
g[g==9] <- NA_real_
dosage_mean <- colMeans(g,na.rm=TRUE);p <- dosage_mean/2
scale <- sqrt(2*p*(1-p));keep <- is.finite(scale) & scale>0
missing <- which(is.na(g),arr.ind=TRUE)
g[missing] <- dosage_mean[missing[,2]]
scaled <- sweep(sweep(g[,keep,drop=FALSE],2,dosage_mean[keep],"-"),2,scale[keep],"/")
stopifnot(all(is.finite(scaled)))
fit <- prcomp(scaled,center=FALSE,scale.=FALSE,rank.=2)
variance <- fit$sdev^2/sum(fit$sdev^2)
crosswalk <- read.csv(file.path(output_dir,"figures","demography","archived_individual_group_crosswalk.csv"))
crosswalk <- crosswalk[match(sample_meta$sample_id,crosswalk$sample_id),]
stopifnot(identical(crosswalk$sample_id,sample_meta$sample_id),nrow(fit$x)==217)
d <- data.frame(sample_meta,archived_group=crosswalk$archived_genetic_group,PC1=fit$x[,1],PC2=fit$x[,2])
stopifnot(all(table(d$group)[c("core","peripheral")]==c(104,113)),
  all(table(d$archived_group)[c("core","peripheral","adaptive_1","adaptive_2")]==c(101,69,23,24)))
write_csv(d,file.path(directory,"PCA_scores_and_group_assignments.csv"))
write_csv(data.frame(component=seq_along(variance),variance_fraction=variance),file.path(directory,"PCA_variance.csv"))
write_csv(data.frame(locus_index=seq_along(keep),mean_dosage=dosage_mean,scale=scale,retained=keep),file.path(directory,"PCA_locus_standardization.csv"))
saveRDS(fit,file.path(directory,"PCA_reconstructed.rds"))

save_session(directory)
marker(directory,list(individuals=217,input_loci=38124,retained_loci=sum(keep),excluded_constant_or_invalid=sum(!keep),missing_dosages_imputed=nrow(missing),source="Exact numerical block of plot_genetic_pca_reconstruction.R; graphics and optional QA removed"))
