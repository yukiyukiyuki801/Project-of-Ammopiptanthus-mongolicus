source(file.path("scripts","revision","common.R"))
library(LEA)
directory <- file.path(output_dir,"imputation_sensitivity")
dir.create(directory,recursive=TRUE,showWarnings=FALSE)
raw <- as.matrix(read.table(file.path(output_dir,"inputs","adaptive_called.tsv")))
archived <- as.matrix(read.table(file.path(output_dir,"inputs","adaptive_imputed.tsv")))
file <- file.path(directory,"adaptive_called.lfmm")
write.table(raw,file,quote=FALSE,row.names=FALSE,col.names=FALSE)
set.seed(20261013)
project <- LEA::snmf(file,K=4,entropy=TRUE,repetitions=10,project="new",seed=20261013)
entropy <- LEA::cross.entropy(project,K=4)
best <- which.min(entropy)
set.seed(20261013)
LEA::impute(project,file,method="mode",K=4,run=best)
fresh <- as.matrix(LEA::read.lfmm(paste0(file,"_imputed.lfmm")))
stopifnot(identical(dim(fresh),dim(raw)),all(fresh %in% 0:2),all(fresh[raw!=9]==raw[raw!=9]))
write_csv(data.frame(repetition=seq_along(entropy),cross_entropy=entropy,selected=seq_along(entropy)==best),file.path(directory,"snmf_imputation_repetitions.csv"))
write_csv(data.frame(sample_id=sample_meta$sample_id,population=sample_meta$population,
  missing_calls=rowSums(raw==9),changed_imputations=rowSums((raw==9) & fresh!=archived)),file.path(directory,"imputation_difference_by_sample.csv"))
current <- as.matrix(read.table(file.path(output_dir,"inputs","current_environment_manifest.tsv")))
future <- read.csv(file.path(output_dir,"climate","future_sample_climate.csv"))
pops <- pop_meta$population
group_mean <- function(x) do.call(rbind,lapply(pops,function(p) colMeans(x[sample_meta$population==p,,drop=FALSE],na.rm=TRUE)))
C <- group_mean(current)
called <- raw;called[called==9] <- NA_real_
datasets <- list(archived_imputation=archived,fresh_seeded_imputation=fresh,observed_calls_no_imputation=called)
results <- list();ri <- 0L
for(source in names(datasets)) {
  F <- group_mean(datasets[[source]])/2
  amplitude <- apply(F,2,function(x) if(all(is.finite(x))) max(x)-min(x) else NA_real_)
  for(ssp in c("ssp126","ssp585")) for(period in c("2041-2070","2071-2100")) for(gcm in unique(future$gcm)) {
    d <- future[future$ssp==ssp & future$period==period & future$gcm==gcm,]
    d <- d[match(sample_meta$sample_id,d$sample_id),]
    xf <- group_mean(as.matrix(d[,bio_names]))
    for(b in seq_along(bio_names)) {
      x <- C[,b];xc <- x-mean(x);fc <- sweep(F,2,colMeans(F),"-")
      slope <- as.numeric(crossprod(xc,fc))/sum(xc^2)
      intercept <- colMeans(F)-slope*mean(x)
      now <- outer(x,slope)+matrix(intercept,19,1637,byrow=TRUE)
      projected <- outer(xf[,b],slope)+matrix(intercept,19,1637,byrow=TRUE)
      r2 <- as.numeric(crossprod(xc,fc))^2/(sum(xc^2)*colSums(fc^2))
      relative <- sweep(abs(abs(F-now)-abs(F-projected)),2,amplitude,"/")
      valid <- is.finite(amplitude) & amplitude>0 & is.finite(r2) & r2>0 & apply(is.finite(relative),2,all)
      value <- as.numeric(relative[,valid,drop=FALSE] %*% r2[valid]/sum(r2[valid]))
      stopifnot(all(is.finite(value)))
      ri <- ri+1L
      results[[ri]] <- data.frame(population=pops,group=pop_meta$group,bio=bio_names[b],ssp=ssp,period=period,gcm=gcm,genotype_source=source,rona=value,used_loci=sum(valid),excluded_loci=sum(!valid),association_scope="All inherited 1637 candidates; p<1 primary rule; raw-call sensitivity has no missing-data LFMM test")
    }
  }
  message("RONA imputation sensitivity complete for ",source)
}
result <- do.call(rbind,results)
write_csv(result,file.path(directory,"rona_imputation_sensitivity.csv"))
old <- read.csv(file.path(output_dir,"rona","rona_population_long.csv"))
key <- function(d) paste(d$population,d$bio,d$ssp,d$period,d$gcm,sep="|")
baseline <- result[result$genotype_source=="archived_imputation",]
idx <- match(key(baseline),key(old));stopifnot(!anyNA(idx),max(abs(baseline$rona-old$rona[idx]))<1e-12)
summary <- list();si <- 0L
for(source in c("fresh_seeded_imputation","observed_calls_no_imputation")) for(bio in bio_names) for(ssp in unique(result$ssp)) for(period in unique(result$period)) {
  a <- result[result$genotype_source==source & result$bio==bio & result$ssp==ssp & result$period==period & result$gcm=="ENSEMBLE",]
  b <- baseline[baseline$bio==bio & baseline$ssp==ssp & baseline$period==period & baseline$gcm=="ENSEMBLE",]
  b <- b[match(a$population,b$population),]
  si <- si+1L
  summary[[si]] <- data.frame(genotype_source=source,bio=bio,ssp=ssp,period=period,max_abs_RONA_difference=max(abs(a$rona-b$rona)),population_rank_rho=cor(a$rona,b$rona,method="spearman"),
    difference_peripheral_minus_core=mean(a$rona[a$group=="peripheral"])-mean(a$rona[a$group=="core"]),baseline_difference=mean(b$rona[b$group=="peripheral"])-mean(b$rona[b$group=="core"]))
}
write_csv(do.call(rbind,summary),file.path(directory,"imputation_effect_stability.csv"))
save_session(directory)
marker(directory,list(fresh_snmf_K=4,repetitions=10,seed=20261013,method="mode",best_repetition=best,missing_calls=sum(raw==9),called_genotype_changes=sum(fresh[raw!=9]!=raw[raw!=9]),changed_missing_imputations=sum(fresh[raw==9]!=archived[raw==9]),RONA_rows=nrow(result),
  comparison="Archived validated imputation versus fresh seeded mode imputation versus observed-call population frequencies; no silent change to primary results"))
