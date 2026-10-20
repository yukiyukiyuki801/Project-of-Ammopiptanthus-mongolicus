source(file.path("scripts","revision","common.R"))
library(gdm)
directory <- file.path(output_dir,"gdm")
dir.create(directory,recursive=TRUE,showWarnings=FALSE)
current <- read.csv(file.path(output_dir,"climate","current_population_climate.csv"))
future <- read.csv(file.path(output_dir,"climate","future_population_climate.csv"))
current <- current[match(pop_meta$population,current$population),]
stopifnot(identical(current$population,pop_meta$population))
fstdata <- read.csv(file.path(output_dir,"inputs","adaptive_fst.csv"),check.names=FALSE)
fstdata <- fstdata[match(pop_meta$population,fstdata$population),]
fst <- as.matrix(fstdata[,pop_meta$population])
pair_frame <- function(first,second,use_actual_coordinates=TRUE) {
  n <- nrow(first)
  stopifnot(nrow(second)==n)
  result <- data.frame(distance=rep(0,n),weights=rep(1,n),
    s1.xCoord=first$longitude,s1.yCoord=first$latitude,
    s2.xCoord=if(use_actual_coordinates) second$longitude else first$longitude,
    s2.yCoord=if(use_actual_coordinates) second$latitude else first$latitude)
  cbind(result,setNames(first[,bio_names],paste0("s1.",bio_names)),setNames(second[,bio_names],paste0("s2.",bio_names)))
}
comb <- t(combn(seq_len(19),2))
site_pairs <- pair_frame(current[comb[,1],],current[comb[,2],])
site_pairs$distance <- fst[comb]
write_csv(data.frame(population_1=pop_meta$population[comb[,1]],population_2=pop_meta$population[comb[,2]],site_pairs),file.path(directory,"site_pairs.csv"))
writeLines(capture.output(list(gdm_transform=body(gdm::gdm.transform),gdm_arguments=formals(gdm::gdm))),file.path(directory,"gdm_implementation_audit.txt"))
models <- list(); diagnostics <- list(); curves <- list();ci <- 0L
for(geographic in c(TRUE,FALSE)) {
  tag <- if(geographic) "geo_TRUE" else "geo_FALSE"
  set.seed(20261013)
  model <- gdm::gdm(site_pairs,geo=geographic)
  stopifnot(inherits(model,"gdm"),all(is.finite(model$coefficients)))
  models[[tag]] <- model
  saveRDS(model,file.path(directory,paste0("model_",tag,".rds")))
  fitted <- as.numeric(predict(model,site_pairs))
  diagnostics[[tag]] <- data.frame(model=tag,metric=c("training_RMSE","explained_deviance"),value=c(sqrt(mean((fitted-site_pairs$distance)^2)),model$explained))
  start <- 1L
  for(j in seq_along(model$predictors)) {
    finish <- start+model$splines[j]-1L
    ci <- ci+1L
    curves[[ci]] <- data.frame(model=tag,predictor=model$predictors[j],fitted_spline_height=sum(model$coefficients[start:finish]),meaning="Fitted contribution, not independently validated predictive importance")
    start <- finish+1L
  }
}
write_csv(do.call(rbind,diagnostics),file.path(directory,"training_diagnostics.csv"))
write_csv(do.call(rbind,curves),file.path(directory,"fitted_predictor_contributions.csv"))
# Leave out every pair touching the test population, avoiding pairwise leakage.
heldout <- list();hi <- 0L
for(tag in names(models)) for(hold in seq_len(19)) {
  train <- !(comb[,1]==hold | comb[,2]==hold)
  set.seed(20261013+hold)
  fit <- tryCatch(gdm::gdm(site_pairs[train,],geo=tag=="geo_TRUE"),error=function(e) NULL)
  pred <- if(is.null(fit)) rep(NA_real_,sum(!train)) else tryCatch(as.numeric(predict(fit,site_pairs[!train,])),error=function(e) rep(NA_real_,sum(!train)))
  hi <- hi+1L
  heldout[[hi]] <- data.frame(model=tag,heldout_population=pop_meta$population[hold],population_1=pop_meta$population[comb[!train,1]],population_2=pop_meta$population[comb[!train,2]],observed=site_pairs$distance[!train],predicted=pred,training_mean=mean(site_pairs$distance[train]),fold_valid=all(is.finite(pred)))
  message("GDM validation ",tag," held out ",pop_meta$population[hold])
}
heldout <- do.call(rbind,heldout)
write_csv(heldout,file.path(directory,"heldout_population_predictions.csv"))
stopifnot(all(heldout$fold_valid))
distance_km <- geosphere::distm(as.matrix(pop_meta[,c("longitude","latitude")]))/1000
score <- function(model,first,second,penalty=FALSE) {
  p <- as.numeric(predict(model,pair_frame(first,second,penalty)))
  stopifnot(all(is.finite(p)),all(p<1))
  eta <- -log1p(-p)-model$intercept
  stopifnot(all(eta> -1e-8))
  1-exp(-pmax(0,eta))
}
results <- list();ri <- 0L
for(ssp in c("ssp126","ssp585")) for(period in c("2041-2070","2071-2100")) for(gcm in unique(future$gcm)) {
  fut <- future[future$ssp==ssp & future$period==period & future$gcm==gcm,]
  fut <- fut[match(pop_meta$population,fut$population),]
  stopifnot(identical(fut$population,pop_meta$population))
  for(tag in names(models)) for(penalty in c(FALSE,TRUE)) {
    model <- models[[tag]]
    for(i in seq_len(19)) {
      local <- score(model,current[i,],fut[i,],penalty)
      forward <- score(model,current[rep(i,19),],fut,penalty)
      reverse <- score(model,fut[rep(i,19),],current,penalty)
      pick <- function(values) {
        tied <- which(abs(values-min(values))<=1e-10)
        tied[order(distance_km[i,tied],pop_meta$population[tied])][1]
      }
      f <- pick(forward);r <- pick(reverse)
      stopifnot(forward[f]<=local+1e-8,reverse[r]<=local+1e-8)
      ri <- ri+1L
      results[[ri]] <- data.frame(population=pop_meta$population[i],group=pop_meta$group[i],ssp=ssp,period=period,gcm=gcm,model=tag,search_scoring=if(penalty) "geographic_penalty" else "environment_only",candidate_domain="19_sampled_population_coordinates",local=local,forward=forward[f],reverse=reverse[r],
        forward_future_destination=pop_meta$population[f],reverse_current_candidate_source=pop_meta$population[r],forward_distance_km=distance_km[i,f],reverse_distance_km=distance_km[i,r],forward_ties=sum(abs(forward-min(forward))<=1e-10),reverse_ties=sum(abs(reverse-min(reverse))<=1e-10))
    }
  }
  message("GDM finite-domain offsets ",ssp," ",period," ",gcm)
}
result <- do.call(rbind,results)
stopifnot(nrow(result)==19*24*2*2)
write_csv(result,file.path(directory,"population_finite_domain_offsets.csv"))
save_session(directory)
marker(directory,list(populations=19,site_pairs=171,heldout_prediction_rows=nrow(heldout),offset_rows=nrow(result),
  models=names(models),scope="Validated GDM fits, local population offsets, sampled-population analogue sensitivity; full landscape-domain forward/reverse search remains pending",geo_handling="geo_TRUE fit and geo_FALSE sensitivity; temporal matching with geographic distance zero or actual-distance penalty, both reported",intercept_handling="removed in link space; roundoff-only nonnegative truncation",ties="within 1e-10, then nearest distance and lexical population ID"))
