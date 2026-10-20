source(file.path("scripts","revision","common.R"))
directory <- file.path(output_dir,"rona")
dir.create(directory,recursive=TRUE,showWarnings=FALSE)
G <- as.matrix(read.table(file.path(output_dir,"inputs","adaptive_imputed.tsv")))
current <- as.matrix(read.table(file.path(output_dir,"inputs","current_environment_manifest.tsv")))
future <- read.csv(file.path(output_dir,"climate","future_sample_climate.csv"))
pops <- pop_meta$population
group_mean <- function(x) do.call(rbind,lapply(pops,function(p) colMeans(x[sample_meta$population==p,,drop=FALSE])))
F <- group_mean(G)/2
C <- group_mean(current)
stopifnot(identical(dim(F),c(19L,1637L)),identical(dim(C),c(19L,5L)))
write_csv(data.frame(population=pops,F,check.names=FALSE),file.path(directory,"allele_frequency_by_population.csv"))
write_csv(data.frame(population=pops,C,check.names=FALSE),file.path(directory,"current_population_means.csv"))
order_group <- order(sample_meta$population,seq_len(nrow(sample_meta)))
write.table(data.frame(sample_meta$population[order_group],current[order_group,,drop=FALSE]),file.path(directory,"pyrona_grouped_current.tsv"),sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
amplitude <- apply(F,2,function(x) max(x)-min(x))
results <- list();counter <- 0L;marker_counts <- list();mcounter <- 0L
climate_long <- list();ccounter <- 0L
for(ssp in c("ssp126","ssp585")) for(period in c("2041-2070","2071-2100")) for(gcm in unique(future$gcm)) {
  subset <- future[future$ssp==ssp & future$period==period & future$gcm==gcm,]
  stopifnot(nrow(subset)==217,!anyDuplicated(subset$sample_id))
  subset <- subset[match(sample_meta$sample_id,subset$sample_id),]
  stopifnot(identical(subset$sample_id,sample_meta$sample_id))
  future_means <- group_mean(as.matrix(subset[,bio_names]))
  for(b in seq_along(bio_names)) {
    x <- C[,b];xf <- future_means[,b]
    xc <- x-mean(x);fc <- sweep(F,2,colMeans(F),"-")
    slope <- as.numeric(crossprod(xc,fc))/sum(xc^2)
    intercept <- colMeans(F)-slope*mean(x)
    pred_now <- outer(x,slope)+matrix(intercept,19,1637,byrow=TRUE)
    pred_future <- outer(xf,slope)+matrix(intercept,19,1637,byrow=TRUE)
    r2 <- as.numeric(crossprod(xc,fc))^2/(sum(xc^2)*colSums(fc^2))
    rel <- sweep(abs(abs(F-pred_now)-abs(F-pred_future)),2,amplitude,"/")
    valid <- amplitude>0 & is.finite(r2) & r2>0 & apply(is.finite(rel),2,all)
    ccounter <- ccounter+1L
    climate_long[[ccounter]] <- data.frame(population=pops,group=pop_meta$group,bio=bio_names[b],ssp=ssp,period=period,gcm=gcm,current=x,future=xf,signed_change=xf-x,absolute_change=abs(xf-x))
    for(k in c(4L,2L,6L)) {
      pvalues <- as.matrix(read.csv(file.path(output_dir,"lfmm",paste0("associations_K",k,".csv")),header=FALSE))[,b]
      qvalues <- p.adjust(pvalues,"BH")
      for(filter in c("inherited_p_lt_1","p_lt_0.05","BH_q_lt_0.05")) {
        selected <- valid & switch(filter,inherited_p_lt_1=pvalues<1,p_lt_0.05=pvalues<.05,BH_q_lt_0.05=qvalues<.05)
        n <- sum(selected)
        values <- if(n>0) as.numeric(rel[,selected,drop=FALSE] %*% r2[selected]/sum(r2[selected])) else rep(NA_real_,19)
        counter <- counter+1L
        results[[counter]] <- data.frame(population=pops,group=pop_meta$group,bio=bio_names[b],ssp=ssp,period=period,gcm=gcm,K=k,association_filter=filter,rona=values,selected_loci=n,excluded_constant_or_invalid=sum(!valid),
          rona_type="absdiff",weighting="population_regression_R_squared",imputation="archived_validated_fixed")
        if(ssp=="ssp126" && period=="2041-2070" && gcm=="ENSEMBLE") {
          mcounter <- mcounter+1L
          marker_counts[[mcounter]] <- data.frame(bio=bio_names[b],K=k,association_filter=filter,selected_loci=n,excluded_constant_or_invalid=sum(!valid))
        }
      }
    }
  }
  message("RONA projected ",ssp," ",period," ",gcm)
}
result <- do.call(rbind,results)
climate <- do.call(rbind,climate_long)
stopifnot(nrow(result)==19*5*24*3*3)
primary <- result[result$K==4 & result$association_filter=="inherited_p_lt_1",]
stopifnot(all(is.finite(primary$rona)))
write_csv(result,file.path(directory,"rona_all_sensitivities.csv"))
write_csv(primary,file.path(directory,"rona_population_long.csv"))
write_csv(climate,file.path(directory,"population_climate_change.csv"))
write_csv(do.call(rbind,marker_counts),file.path(directory,"association_marker_counts.csv"))
save_session(directory)
code <- system2(Sys.getenv("AMMO_PYTHON", unset="python3"),file.path(root,"scripts","revision","verify_rona_reference.py"))
stopifnot(code==0)
marker(directory,list(primary_rows=nrow(primary),all_sensitivity_rows=nrow(result),reference_numeric_check="passed pyRona default absdiff/R2-weighted calculation and grouped sample parsing",primary_K=4,primary_association_filter="p<1 as inherited command; conditional on preselected candidates",zero_association_handling="NA, never zero",constant_population_frequency_loci="excluded with counts; amplitude-normalized statistic undefined at amplitude zero",uncertainty_note="Locus dispersion is not independent population-level uncertainty"))
