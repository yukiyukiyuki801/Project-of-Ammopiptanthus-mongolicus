source(file.path("scripts","revision","common.R"))
directory <- file.path(output_dir,"landscape_statistics")
dir.create(directory,recursive=TRUE,showWarnings=FALSE)
writeLines(c("Reviewer scope: R2.M3 explicitly requests populations or spatial blocks rather than 1000 autocorrelated pixels, with effect sizes and confidence intervals.",
  "R1.40 requests explicit scenario comparisons. This stage only summarizes existing full-habitat population matching results. No GDM fitting or landscape search is repeated.",
  "Both saved scoring definitions are retained; geographic-penalty scoring uses inherited coordinate-degree distances and is not silently treated as a metric-distance model.",
  "Class tests are unpaired (9 core vs10 peripheral); time/scenario contrasts are paired by population. Bootstrap intervals assume exchangeable populations and do not remove spatial dependence."),file.path(directory,"reviewer_scope.txt"))
x <- read.csv(file.path(output_dir,"landscape","full_habitat_domain_population_matches.csv"))
x <- x[x$radius_km=="unrestricted",]
stopifnot(nrow(x)==19*24*2*2)
keys <- c("population","group","ssp","period","gcm","search_scoring")
forward <- x[x$direction=="forward",c(keys,"offset","local_offset")]
reverse <- x[x$direction=="reverse",c(keys,"offset","local_offset")]
names(forward)[names(forward)=="offset"] <- "forward"
names(reverse)[names(reverse)=="offset"] <- "reverse"
joined <- merge(forward,reverse,by=keys,suffixes=c(".f",".r"),sort=FALSE)
stopifnot(nrow(joined)==19*24*2,all(abs(joined$local_offset.f-joined$local_offset.r)<1e-10))
joined$local <- joined$local_offset.f
joined$local_offset.f <- joined$local_offset.r <- NULL
write_csv(joined,file.path(directory,"population_offsets_wide.csv"))
R <- 2000L
effects <- list();i <- 0L
for(score in unique(joined$search_scoring)) for(ssp in c("ssp126","ssp585")) for(period in c("2041-2070","2071-2100")) for(gcm in unique(joined$gcm)) {
  d <- joined[joined$search_scoring==score & joined$ssp==ssp & joined$period==period & joined$gcm==gcm,]
  stopifnot(nrow(d)==19,!anyDuplicated(d$population))
  for(metric in c("local","forward","reverse")) {
    a <- d[[metric]][d$group=="core"];b <- d[[metric]][d$group=="peripheral"]
    stopifnot(length(a)==9,length(b)==10,all(is.finite(c(a,b))))
    ci <- c(NA_real_,NA_real_);p <- rb <- NA_real_
    if(gcm=="ENSEMBLE") {
      boot <- replicate(R,mean(sample(b,length(b),replace=TRUE))-mean(sample(a,length(a),replace=TRUE)))
      ci <- unname(quantile(boot,c(.025,.975)))
      test <- wilcox.test(b,a,paired=FALSE,exact=FALSE,correct=TRUE)
      p <- test$p.value;rb <- 2*as.numeric(test$statistic)/(length(a)*length(b))-1
    }
    i <- i+1L
    effects[[i]] <- data.frame(search_scoring=score,ssp=ssp,period=period,gcm=gcm,metric=metric,
      n_core=9,n_peripheral=10,mean_core=mean(a),mean_peripheral=mean(b),
      difference_peripheral_minus_core=mean(b)-mean(a),CI_low=ci[1],CI_high=ci[2],
      rank_biserial=rb,p_unpaired_rank_sum=p)
  }
}
effects <- do.call(rbind,effects);effects$q_BH <- NA_real_
family <- interaction(effects$search_scoring,effects$ssp,effects$period,drop=TRUE)
for(f in unique(family)) {
  ind <- which(family==f & effects$gcm=="ENSEMBLE")
  stopifnot(length(ind)==3)
  effects$q_BH[ind] <- p.adjust(effects$p_unpaired_rank_sum[ind],method="BH")
}
write_csv(effects,file.path(directory,"population_class_effects.csv"))
contrasts <- list();i <- 0L
for(score in unique(joined$search_scoring)) {
  defs <- list(
    list(name="ssp126_late_minus_mid",a=c("ssp126","2041-2070"),b=c("ssp126","2071-2100")),
    list(name="ssp585_late_minus_mid",a=c("ssp585","2041-2070"),b=c("ssp585","2071-2100")),
    list(name="mid_ssp585_minus_ssp126",a=c("ssp126","2041-2070"),b=c("ssp585","2041-2070")),
    list(name="late_ssp585_minus_ssp126",a=c("ssp126","2071-2100"),b=c("ssp585","2071-2100")))
  for(def in defs) {
    a <- joined[joined$search_scoring==score & joined$gcm=="ENSEMBLE" & joined$ssp==def$a[1] & joined$period==def$a[2],]
    b <- joined[joined$search_scoring==score & joined$gcm=="ENSEMBLE" & joined$ssp==def$b[1] & joined$period==def$b[2],]
    a <- a[match(pop_meta$population,a$population),];b <- b[match(pop_meta$population,b$population),]
    stopifnot(identical(a$population,b$population),nrow(a)==19)
    for(metric in c("local","forward","reverse")) {
      delta <- b[[metric]]-a[[metric]]
      ci <- unname(quantile(replicate(R,mean(sample(delta,19,replace=TRUE))),c(.025,.975)))
      test <- wilcox.test(b[[metric]],a[[metric]],paired=TRUE,exact=FALSE,correct=TRUE)
      i <- i+1L
      contrasts[[i]] <- data.frame(search_scoring=score,contrast=def$name,metric=metric,n_paired_populations=19,
        mean_paired_difference=mean(delta),median_paired_difference=median(delta),CI_low=ci[1],CI_high=ci[2],
        p_paired_signed_rank=test$p.value)
    }
  }
}
contrasts <- do.call(rbind,contrasts);contrasts$q_BH <- NA_real_
family <- interaction(contrasts$search_scoring,contrasts$contrast,drop=TRUE)
for(f in unique(family)) {ind <- which(family==f);contrasts$q_BH[ind] <- p.adjust(contrasts$p_paired_signed_rank[ind],"BH")}
write_csv(contrasts,file.path(directory,"paired_scenario_contrasts.csv"))
save_session(directory)
marker(directory,list(review_comments=c("R2.M3","R1.40"),source="Saved full-habitat matching output; no refitting or new searches",population_rows=nrow(joined),class_effect_rows=nrow(effects),scenario_contrast_rows=nrow(contrasts),
  bootstrap_repetitions=R,seed=20261013,multiple_testing="BH across three offset metrics within scoring x scenario-period for class tests, or within scoring x paired contrast for scenario tests",spatial_limitation="Population-level tests avoid pixel pseudoreplication; intervals/tests still assume exchangeable populations and are not spatially adjusted."))
