source(file.path("scripts", "revision", "common.R"))
library(LEA)
directory <- file.path(output_dir,"lfmm")
dir.create(directory,recursive=TRUE,showWarnings=FALSE)
Y <- as.matrix(read.table(file.path(output_dir,"inputs","adaptive_imputed.tsv")))
X <- as.matrix(read.table(file.path(output_dir,"inputs","current_environment_manifest.tsv")))
stopifnot(identical(dim(Y),c(217L,1637L)),identical(dim(X),c(217L,5L)),all(is.finite(X)),all(Y %in% 0:2))
writeLines(capture.output(list(fit_arguments=formals(LEA::lfmm2),test_arguments=formals(LEA::lfmm2.test),test_body=body(LEA::lfmm2.test))),file.path(directory,"LEA_implementation_audit.txt"))
manifest <- read.csv(file.path(output_dir,"inputs","sample_metadata.csv"))
stopifnot(identical(manifest$sample_id,sample_meta$sample_id))
summaries <- list()
for(k in c(4L,2L,6L)) {
  model_file <- file.path(directory,paste0("model_K",k,".rds"))
  assoc_file <- file.path(directory,paste0("associations_K",k,".csv"))
  checkpoint <- file.path(directory,paste0("K",k,"_complete.json"))
  if(file.exists(checkpoint) && file.exists(model_file) && file.exists(assoc_file)) {
    message("Validated checkpoint exists for K=",k)
    p <- as.matrix(read.csv(assoc_file,header=FALSE))
  } else {
    set.seed(20261013+k)
    message("Fitting LFMM2 K=",k)
    mod <- LEA::lfmm2(input=Y,env=X,K=k)
    pv <- LEA::lfmm2.test(object=mod,input=Y,env=X,linear=TRUE)
    p <- pv$pvalues
    original_dim <- dim(p)
    if(identical(dim(p),c(5L,1637L))) p <- t(p)
    stopifnot(identical(dim(p),c(1637L,5L)),all(is.finite(p)),all(p>=0 & p<=1))
    saveRDS(list(model=mod,test=pv,K=k,sample_ids=sample_meta$sample_id,bio_names=bio_names,
      seed=20261013+k,imputation="archived validated imputation, all called genotypes preserved; held fixed for K sensitivity"),model_file)
    write.table(p,assoc_file,sep=",",row.names=FALSE,col.names=FALSE)
    jsonlite::write_json(list(validated=TRUE,K=k,original_pvalue_dimensions=original_dim,output_dimensions=dim(p),
      test_call="LEA::lfmm2.test(linear=TRUE); all other arguments at recorded version defaults"),checkpoint,pretty=TRUE,auto_unbox=TRUE)
  }
  stopifnot(identical(dim(p),c(1637L,5L)),all(is.finite(p)))
  q <- apply(p,2,p.adjust,method="BH")
  write.table(q,file.path(directory,paste0("qvalues_K",k,".csv")),sep=",",row.names=FALSE,col.names=FALSE)
  summaries[[as.character(k)]] <- data.frame(K=k,bio=bio_names,p_lt_1=colSums(p<1),p_lt_005=colSums(p<.05),BH_q_lt_005=colSums(q<.05),median_p=apply(p,2,median))
}
write_csv(do.call(rbind,summaries),file.path(directory,"association_summary.csv"))
save_session(directory)
marker(directory,list(primary_K=4,sensitivity_K=c(2,6),samples=217,loci=1637,covariates=5,
  discovery_scope="Conditional analysis of inherited candidate SNPs, not a new unbiased genome scan",
  K_sensitivity_scope="LFMM latent-factor K changes; archived genotype imputation held fixed and validated against called genotypes"))
