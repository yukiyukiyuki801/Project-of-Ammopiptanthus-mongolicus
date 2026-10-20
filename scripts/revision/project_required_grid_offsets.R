# Missing original-design Figure 3 / S9 / S10 numerics only; no model fitting.
source('scripts/revision/common.R')
library(gdm); library(terra); library(geosphere); library(Rcpp)
directory <- file.path(output_dir,'figure_revision_original_design','required_grid_projection')
dir.create(directory,recursive=TRUE,showWarnings=FALSE)
args <- commandArgs(TRUE); validate_only <- any(c('--validate','--validate-only','--benchmark') %in% args)
# Resolve the installed PROJ development library instead of author-machine paths.
stopifnot(nzchar(Sys.which("pkg-config")))
proj_cflags <- system2("pkg-config", c("--cflags","proj"), stdout=TRUE)
proj_libs <- system2("pkg-config", c("--libs","proj"), stdout=TRUE)
stopifnot(is.null(attr(proj_cflags,"status")),is.null(attr(proj_libs,"status")))
Sys.setenv(PKG_CPPFLAGS=paste(proj_cflags,collapse=" "),PKG_LIBS=paste(proj_libs,collapse=" "))
Rcpp::sourceCpp(file.path(root,'scripts/revision/exact_gdm_grid_search.cpp'),rebuild=TRUE)
model_path <- file.path(output_dir,'gdm','model_geo_TRUE.rds'); model <- readRDS(model_path)
stopifnot(model$geo,identical(as.integer(model$splines),rep(3L,6)),all(model$coefficients>=0))
envmodel <- model;envmodel$geo <- FALSE;envmodel$predictors <- model$predictors[-1];envmodel$coefficients <- model$coefficients[-(1:3)];envmodel$knots <- model$knots[-(1:3)];envmodel$splines <- model$splines[-1]
domain_path <- file.path(output_dir,'landscape','current_candidate_domain.csv'); domain <- read.csv(domain_path)
n_grid <- sum(domain$candidate_type=='habitat_grid_cell');stopifnot(nrow(domain)==n_grid+19)
transform_domain <- function(d) cbind(as.matrix(gdm.transform(envmodel,d[,bio_names])),as.matrix(d[,c('longitude','latitude')]))
cur <- transform_domain(domain)
base <- file.path(bundle,'ammo','gf','chelsa2_crop_selected')
pair_frame <- function(first,second,penalty){n<-nrow(second);cbind(data.frame(distance=rep(0,n),weights=rep(1,n),s1.xCoord=rep(first$longitude,n),s1.yCoord=rep(first$latitude,n),s2.xCoord=if(penalty)second$longitude else rep(first$longitude,n),s2.yCoord=if(penalty)second$latitude else rep(first$latitude,n)),setNames(first[rep(1,n),bio_names],paste0('s1.',bio_names)),setNames(second[,bio_names],paste0('s2.',bio_names)))}
score <- function(first,second,penalty) {raw<-as.numeric(predict(model,pair_frame(first,second,penalty))); -expm1(-pmax(0,-log1p(-raw)-model$intercept))}
checks <- list(); ci <- 0L; times <- list()
for(ssp in c('ssp126','ssp585')) for(period in c('2041-2070','2071-2100')) {
 tag <- paste(ssp,period,'ENSEMBLE',sep='_'); short <- if(period=='2041-2070')'2040' else '2070'
 future <- rast(file.path(base,paste0(ssp,'_',short),paste0('b',bio_numbers,'.asc')));names(future)<-bio_names
 fv <- terra::extract(future,domain[,c('longitude','latitude')],method='simple',ID=FALSE)
 stopifnot(all(is.finite(as.matrix(fv))))
 futdomain <- cbind(domain[,c('longitude','latitude')],fv); fut <- transform_domain(futdomain)
 saved <- read.csv(file.path(output_dir,'landscape',paste0(tag,'_matches.csv')))
 # Validate deterministic nontrivial small target grid and all saved focal populations.
 small <- unique(as.integer(round(seq(1,nrow(domain),length.out=521)))); qsmall <- unique(as.integer(round(seq(1,n_grid,length.out=17))))
 for(penalty in c(TRUE,FALSE)) for(direction in c('forward','reverse')) {
   targets <- if(direction=='forward')fut else cur; queries <- if(direction=='forward')cur else fut
   td <- if(direction=='forward')futdomain else domain; qd <- if(direction=='forward')domain else futdomain
   ans <- exact_gdm_search(targets[small,,drop=FALSE],queries[qsmall,,drop=FALSE],model$knots[1:3],model$coefficients[1:3],penalty)
   for(i in seq_along(qsmall)){
     truth<-score(qd[qsmall[i],,drop=FALSE],td[small,,drop=FALSE],penalty); minimum<-min(truth); tie<-which(abs(truth-minimum)<=1e-10)
     dist<-distGeo(as.matrix(qd[qsmall[i],c('longitude','latitude')]),as.matrix(td[small[tie],c('longitude','latitude')]))/1000
     w<-tie[order(dist,small[tie])][1];stopifnot(abs(ans[i,1]-minimum)<2e-12,small[ans[i,2]]==small[w],abs(ans[i,3]-min(dist))<1e-8,ans[i,4]==length(tie))
   }
   popidx<-n_grid+seq_len(19); elapsed<-system.time(popans<-exact_gdm_search(targets,queries[popidx,,drop=FALSE],model$knots[1:3],model$coefficients[1:3],penalty))[['elapsed']]
   sc<-if(penalty)'inherited_geographic_penalty_degrees' else 'environment_only'; old<-saved[saved$direction==direction&saved$search_scoring==sc&saved$radius_km=='unrestricted',];old<-old[match(pop_meta$population,old$population),]
   stopifnot(nrow(old)==19,!anyNA(old$offset),max(abs(popans[,1]-old$offset))<2e-12,all(popans[,2]==old$endpoint_index),max(abs(popans[,3]-old$distance_km))<1e-8,all(popans[,4]==old$tie_count))
   ci<-ci+1L;checks[[ci]]<-data.frame(scenario=tag,scoring=sc,direction=direction,small_queries=length(qsmall),small_targets=length(small),population_queries=19,score_max_error=max(abs(popans[,1]-old$offset)),endpoint_match=TRUE,tie_count_match=TRUE,elapsed_seconds=elapsed,mean_leaf_evaluations=mean(popans[,5]))
   if(!validate_only && penalty){
      csv_path <- file.path(directory,paste(tag,direction,sc,'grid.csv',sep='_')); checkpoint <- sub('grid.csv$','grid_complete.json',csv_path)
      if(file.exists(csv_path) && file.exists(checkpoint) && isTRUE(jsonlite::read_json(checkpoint)$validated)){message('Reusing completed projection ',csv_path);next}
      gc();elapsed<-system.time(ans<-exact_gdm_search(targets,queries[seq_len(n_grid),,drop=FALSE],model$knots[1:3],model$coefficients[1:3],penalty))[['elapsed']]
      local_eta<-rowSums(abs(cur[seq_len(n_grid),1:5,drop=FALSE]-fut[seq_len(n_grid),1:5,drop=FALSE]));local<- -expm1(-local_eta)
      result<-data.frame(candidate_index=seq_len(n_grid),longitude=domain$longitude[seq_len(n_grid)],latitude=domain$latitude[seq_len(n_grid)],ssp=ssp,period=period,direction=direction,search_scoring=sc,local_offset=local,offset=ans[,1],endpoint_index=as.integer(ans[,2]),distance_km=ans[,3],tie_count=as.integer(ans[,4]),endpoint_longitude=domain$longitude[ans[,2]],endpoint_latitude=domain$latitude[ans[,2]],leaf_evaluations=ans[,5])
      stopifnot(all(result$offset<=result$local_offset+1e-8));write_csv(result,csv_path)
      jsonlite::write_json(list(validated=TRUE,scenario=tag,direction=direction,scoring=sc,grid_cells=n_grid,candidates=nrow(domain),seconds=elapsed,no_model_refitting=TRUE),checkpoint,pretty=TRUE,auto_unbox=TRUE)
      times[[length(times)+1L]]<-data.frame(scenario=tag,direction=direction,scoring=sc,seconds=elapsed,queries=n_grid,evaluations=sum(ans[,5]))
   }
 }
}
write_csv(do.call(rbind,checks),file.path(directory,'algorithm_validation.csv'))
if(length(times))write_csv(do.call(rbind,times),file.path(directory,'runtime.csv'))
jsonlite::write_json(list(validated=TRUE,validation_only=validate_only,scenarios=4,production_searches=8,validation_searches=16,grid_cells=n_grid,candidates=nrow(domain),no_model_refitting=TRUE,score_tolerance=1e-10,numerical_pruning_safety=2e-14),file.path(directory,if(validate_only)'algorithm_validation.json' else 'validation.json'),pretty=TRUE,auto_unbox=TRUE)
save_session(directory)
