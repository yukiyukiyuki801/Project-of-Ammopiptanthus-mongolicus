source(file.path("scripts","revision","common.R"))
directory <- file.path(output_dir,"statistics")
dir.create(directory,recursive=TRUE,showWarnings=FALSE)
stopifnot(system2(Sys.getenv("AMMO_PYTHON", unset="python3"),file.path(root,"scripts","revision","calculate_diversity.py"))==0)
rona <- read.csv(file.path(output_dir,"rona","rona_population_long.csv"))
climate <- read.csv(file.path(output_dir,"rona","population_climate_change.csv"))
gdm <- read.csv(file.path(output_dir,"gdm","population_finite_domain_offsets.csv"))
gdm <- gdm[gdm$model=="geo_TRUE" & gdm$search_scoring=="environment_only",]
diversity <- read.csv(file.path(directory,"population_diversity.csv"))
effects <- list();ei <- 0L;spatial <- list();si <- 0L
# Coordinate-only blocks are a sensitivity diagnostic. A block deletion that
# removes all members of a class is explicitly non-estimable, never filled in.
coords <- scale(pop_meta[,c("longitude","latitude")])
set.seed(20261013)
blocks <- kmeans(coords,centers=3,nstart=30)$cluster
write_csv(data.frame(population=pop_meta$population,geographic_block=blocks),file.path(directory,"geographic_block_definition.csv"))
compare <- function(data,value,analysis,variable,ssp,period,gcm,domain) {
  x <- data[data$group=="core",value];y <- data[data$group=="peripheral",value]
  stopifnot(length(x)==9,length(y)==10,all(is.finite(c(x,y))))
  delta <- mean(y)-mean(x)
  if(gcm=="ENSEMBLE") {
    boot <- replicate(2000,mean(sample(y,length(y),TRUE))-mean(sample(x,length(x),TRUE)))
    interval <- unname(quantile(boot,c(.025,.975)))
    test <- wilcox.test(y,x,paired=FALSE,exact=FALSE,correct=TRUE)
    p <- test$p.value
  } else {interval <- c(NA,NA);p <- NA_real_}
  ei <<- ei+1L
  effects[[ei]] <<- data.frame(analysis=analysis,variable=variable,ssp=ssp,period=period,gcm=gcm,domain=domain,n_core=9,n_peripheral=10,mean_core=mean(x),mean_peripheral=mean(y),difference_peripheral_minus_core=delta,CI_low=interval[1],CI_high=interval[2],
    median_difference=median(y)-median(x),rank_biserial=mean(sign(outer(y,x,"-"))),p_unpaired_rank_sum=p,bootstrap_assumption="Exchangeable populations within class; spatial sensitivity reported separately")
  if(gcm=="ENSEMBLE") for(block in sort(unique(blocks))) {
    omit <- pop_meta$population[blocks==block]
    kept <- data[!data$population %in% omit,]
    a <- kept[kept$group=="core",value];b <- kept[kept$group=="peripheral",value]
    si <<- si+1L
    spatial[[si]] <<- data.frame(analysis=analysis,variable=variable,ssp=ssp,period=period,omitted_block=block,n_core=length(a),n_peripheral=length(b),difference=if(length(a)>0 && length(b)>0) mean(b)-mean(a) else NA_real_,estimable=length(a)>0 && length(b)>0,interpretation="Coordinate-defined block-deletion sensitivity; not a spatially corrected significance test")
  }
}
for(ssp in unique(rona$ssp)) for(period in unique(rona$period)) for(gcm in unique(rona$gcm)) {
  for(bio in bio_names) {
    d <- rona[rona$ssp==ssp & rona$period==period & rona$gcm==gcm & rona$bio==bio,]
    compare(d,"rona","RONA",bio,ssp,period,gcm,"19_sampled_populations")
  }
  d <- gdm[gdm$ssp==ssp & gdm$period==period & gdm$gcm==gcm,]
  for(metric in c("local","forward","reverse")) compare(d,metric,"GDM",metric,ssp,period,gcm,"19_sampled_population_candidate_domain")
}
effects <- do.call(rbind,effects)
effects$q_BH <- NA_real_
for(analysis in unique(effects$analysis)) for(ssp in unique(effects$ssp)) for(period in unique(effects$period)) {
  idx <- effects$analysis==analysis & effects$ssp==ssp & effects$period==period & effects$gcm=="ENSEMBLE"
  effects$q_BH[idx] <- p.adjust(effects$p_unpaired_rank_sum[idx],"BH")
}
write_csv(effects,file.path(directory,"population_group_effects.csv"))
write_csv(do.call(rbind,spatial),file.path(directory,"spatial_block_deletion_sensitivity.csv"))
correlations <- list();ci <- 0L
for(ssp in unique(rona$ssp)) for(period in unique(rona$period)) for(bio in bio_names) {
  d <- rona[rona$ssp==ssp & rona$period==period & rona$bio==bio & rona$gcm=="ENSEMBLE",]
  c <- climate[climate$ssp==ssp & climate$period==period & climate$bio==bio & climate$gcm=="ENSEMBLE",]
  c <- c[match(d$population,c$population),]
  stopifnot(identical(d$population,c$population))
  for(type in c("signed_change","absolute_change")) {
    x <- c[[type]];y <- d$rona
    test <- cor.test(x,y,method="spearman",exact=FALSE)
    boot <- replicate(2000,{i <- sample(seq_along(x),length(x),TRUE);suppressWarnings(cor(x[i],y[i],method="spearman"))})
    interval <- quantile(boot,c(.025,.975),na.rm=TRUE)
    ci <- ci+1L
    correlations[[ci]] <- data.frame(ssp=ssp,period=period,bio=bio,change_type=type,n=19,rho=unname(test$estimate),CI_low=interval[1],CI_high=interval[2],p=test$p.value,min_change=min(x),max_change=max(x),constant_bootstrap_resamples=sum(!is.finite(boot)))
  }
}
correlations <- do.call(rbind,correlations);correlations$q_BH <- NA_real_
for(ssp in unique(correlations$ssp)) for(period in unique(correlations$period)) for(type in unique(correlations$change_type)) {
  idx <- correlations$ssp==ssp & correlations$period==period & correlations$change_type==type
  correlations$q_BH[idx] <- p.adjust(correlations$p[idx],"BH")
}
write_csv(correlations,file.path(directory,"rona_climate_correlations.csv"))

# Fixed low-dimensional formulas, never stepwise selection. Full/background
# marker definitions are kept separate; no unvalidated load term is included.
formulae <- list(intercept=y~1,climate=y~climate_change,diversity=y~Ho,climate_and_diversity=y~climate_change+Ho)
predictions <- list();coefficients <- list();models <- list();mi <- 0L;pi <- 0L;coi <- 0L;influence <- list();ii <- 0L
for(ssp in unique(rona$ssp)) for(period in unique(rona$period)) for(bio in bio_names) for(marker_set in c("full","nonadaptive")) {
  d <- rona[rona$ssp==ssp & rona$period==period & rona$bio==bio & rona$gcm=="ENSEMBLE",]
  c <- climate[climate$ssp==ssp & climate$period==period & climate$bio==bio & climate$gcm=="ENSEMBLE",]
  div <- diversity[diversity$subset==marker_set,]
  data <- data.frame(population=d$population,y=d$rona,climate_change=c$absolute_change[match(d$population,c$population)],Ho=div$Ho[match(d$population,div$population)])
  stopifnot(nrow(data)==19,all(is.finite(as.matrix(data[,-1]))))
  for(name in names(formulae)) {
    model <- lm(formulae[[name]],data)
    model_key <- paste(ssp,period,bio,marker_set,name,sep="_")
    models[[model_key]] <- model
    ct <- summary(model)$coefficients
    interval <- confint(model)
    coi <- coi+1L
    coefficients[[coi]] <- data.frame(model=model_key,term=rownames(ct),estimate=ct[,1],SE=ct[,2],CI_low=interval[,1],CI_high=interval[,2],p=ct[,4],residual_df=df.residual(model),adjusted_R2=summary(model)$adj.r.squared,predictor_correlation=cor(data$climate_change,data$Ho))
    for(j in seq_len(19)) {
      fit <- lm(formulae[[name]],data[-j,])
      pred <- as.numeric(predict(fit,newdata=data[j,,drop=FALSE]))
      pi <- pi+1L
      predictions[[pi]] <- data.frame(model=model_key,ssp=ssp,period=period,bio=bio,marker_set=marker_set,candidate=name,heldout_population=data$population[j],observed=data$y[j],predicted=pred,training_mean=mean(data$y[-j]))
      beta <- coef(fit)
      ii <- ii+1L
      influence[[ii]] <- data.frame(model=model_key,omitted_population=data$population[j],term=names(beta),estimate=unname(beta))
    }
  }
}
predictions <- do.call(rbind,predictions)
stopifnot(all(is.finite(predictions$predicted)))
performance <- do.call(rbind,lapply(split(predictions,predictions$model),function(d) data.frame(model=d$model[1],heldout_n=nrow(d),RMSE=sqrt(mean((d$predicted-d$observed)^2)),MAE=mean(abs(d$predicted-d$observed)),baseline_RMSE=sqrt(mean((d$training_mean-d$observed)^2)),prediction_bias=mean(d$predicted-d$observed),negative_prediction_count=sum(d$predicted<0))))
write_csv(predictions,file.path(directory,"heldout_model_predictions.csv"))
write_csv(performance,file.path(directory,"model_validation_performance.csv"))
write_csv(do.call(rbind,coefficients),file.path(directory,"model_coefficients.csv"))
write_csv(do.call(rbind,influence),file.path(directory,"leave_one_population_coefficient_stability.csv"))
saveRDS(models,file.path(directory,"linear_models.rds"))
save_session(directory)
marker(directory,list(group_effect_rows=nrow(effects),correlation_rows=nrow(correlations),models=length(models),heldout_prediction_rows=nrow(predictions),
  statistics_scope="First population-based numerical batch. Exchangeability-based intervals plus coordinate-block sensitivity; no claim of full spatial correction.",
  GDM_scope="Local offsets and sampled-population candidate-domain forward/reverse sensitivity; full-landscape forward/reverse main analysis still pending",
  model_scope="Prespecified linear candidates, full/nonadaptive heterozygosity sensitivity; no load terms, no stepwise search; residual spatial and GAM diagnostics next stage",
  causality="Associations and conditional prediction of estimated RONA do not establish causal effects or fitness validity"))
