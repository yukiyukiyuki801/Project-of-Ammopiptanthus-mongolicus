#!/usr/bin/env Rscript
# Reconstruct missing SDM numerical sources under explicit user authorization.
# Source bundle and population class labels are never changed.
source('scripts/revision/common.R')
options(java.parameters='-Xmx6g', timeout=1200)
suppressPackageStartupMessages({library(terra); library(dismo); library(jsonlite)})
directory <- file.path(output_dir,'sdm_reconstruction')
dir.create(directory,recursive=TRUE,showWarnings=FALSE)
log_path <- file.path(directory,'progress.log')
log <- function(...) { s <- paste(format(Sys.time(),tz='UTC'),paste(...,collapse=' '));message(s);cat(s,'\n',file=log_path,append=TRUE) }
seed <- 20261014L
predictors <- c(bio_names,'elevation')
settings <- list(authorization='User explicitly requested reconstructing missing niche/SDM results with supplied data; missing numerical sources prevent redrawing S11 and Fig6.',
 scope='New documented reconstruction, not exact historical replication; source outputs and population classes unchanged.',
 occurrence_source='recovery_audit/published_2022_occurrences.csv; 37 supplied published occurrences',
 occurrence_thinning='Retain first record per occupied native current CHELSA cell; no outcome-based filtering.',
 background='10000 uniformly sampled valid cells across supplied current CHELSA crop (~100–110 E, 36–44 N), excluding occupied cells; inherited study extent is operational calibration domain, not demonstrated accessible area.',
 predictors=predictors,elevation='WorldClim v2.1 30 arcsecond country tiles CHN and MNG; static for future; bilinear alignment onto native CHELSA grid',
 elevation_source='https://www.worldclim.org/data/worldclim21.html',
 model='MaxEnt 3.4.4 via dismo; auto features, default regularization multiplier 1; logistic output; no tuning or optional model benchmark.',
 maxent_args=c('autofeature=true','betamultiplier=1','outputformat=logistic','randomseed=false','threads=1','removeduplicates=false','responsecurves=false','jackknife=true','writeplotdata=true'),
 seed=seed,cross_validation='Seeded balanced random ten-fold presence holdout; background partitioned into matching folds; finite sample AUC is presence versus background, not presence versus confirmed absence; random folds do not establish spatial transferability.',
 threshold='Maximum training sensitivity plus specificity (MTSS) on final all-occurrence model versus calibration background; threshold fixed for all futures.',
 projection='Five individual GCMs × two SSPs × two periods; native future climate bilinear-resampled onto current grid; constant elevation; arithmetic mean/SD/min/max over five suitability predictions. No averaging climates before nonlinear projection.',
 extrapolation='MaxEnt default clamping; per-pixel number of predictors outside observed training/background ranges saved separately; unsuitable/clamped surfaces are not realized colonization.',
 uncertainty='GCM spread is model spread, not independent biological replication; ten random folds have few held-out occurrences and do not prove robust transferability.')
write_json(settings,file.path(directory,'settings_prespecified.json'),pretty=TRUE,auto_unbox=TRUE)
log('Settings saved before model fitting.')
# Country tiles avoid a 323 MB global download. Downloads are durable and checked.
dem_dir <- file.path(directory,'elevation_source');dir.create(dem_dir,showWarnings=FALSE)
dem_paths <- file.path(dem_dir,paste0(c('CHN','MNG'),'_wc2.1_30s_elev.tif'))
for(i in seq_along(dem_paths)) if(!file.exists(dem_paths[i])) {
 url <- paste0('https://geodata.ucdavis.edu/climate/worldclim/2_1/tiles/iso/',basename(dem_paths[i]))
 tmp <- paste0(dem_paths[i],'.part');log('Downloading',url)
 status <- system2('curl',c('--fail','--location','--retry','3','--connect-timeout','30','--max-time','1200','--output',shQuote(tmp),shQuote(url)))
 stopifnot(status==0L,file.info(tmp)$size>1000000);stopifnot(file.rename(tmp,dem_paths[i]))
}
write_csv(data.frame(file=dem_paths,bytes=file.info(dem_paths)$size,md5=unname(tools::md5sum(dem_paths))),file.path(directory,'elevation_source_checksums.csv'))
current <- rast(file.path(bundle,'ammo/gf/chelsa2_crop_selected/current',paste0('b',bio_numbers,'.asc')));names(current)<-bio_names
if(!nzchar(crs(current)))crs(current)<-'EPSG:4326'
template <- current[[1]]
# Crop source tiles before aligning; no global grid is materialized.
dems <- lapply(dem_paths,function(p)crop(rast(p),ext(template)+.05))
elevation <- resample(do.call(merge,dems),template,method='bilinear');names(elevation)<-'elevation'
writeRaster(elevation,file.path(directory,'elevation_aligned.tif'),overwrite=TRUE)
current <- c(current,elevation);names(current)<-predictors
occ <- read.csv(file.path(root,'recovery_audit/published_2022_occurrences.csv'))
occ$cell <- cellFromXY(template,occ[,c('longitude','latitude')]);occ$retained <- !is.na(occ$cell)&!duplicated(occ$cell)
occ_values <- terra::extract(current,occ[,c('longitude','latitude')],ID=FALSE)
occ$retained <- occ$retained&complete.cases(occ_values)
write_csv(cbind(occ,occ_values),file.path(directory,'occurrence_audit.csv'))
presence <- as.data.frame(occ_values[occ$retained,,drop=FALSE]);stopifnot(nrow(presence)>=20)
all_values <- values(current,mat=TRUE)
valid <- which(apply(all_values,1,function(x)all(is.finite(x))))
set.seed(seed);bg_cells <- sample(setdiff(valid,occ$cell[occ$retained]),min(10000L,length(setdiff(valid,occ$cell[occ$retained]))))
background <- as.data.frame(all_values[bg_cells,,drop=FALSE]);names(background)<-predictors
write_csv(cbind(cell=bg_cells,as.data.frame(xyFromCell(template,bg_cells)),background),file.path(directory,'background_sample.csv'))
write_csv(as.data.frame(cor(rbind(presence,background))),file.path(directory,'predictor_correlations.csv'))
set.seed(seed);pfold <- sample(rep(seq_len(10),length.out=nrow(presence)));bfold <- sample(rep(seq_len(10),length.out=nrow(background)))
write_csv(data.frame(source_row=which(occ$retained),label=occ$published_label[occ$retained],fold=pfold),file.path(directory,'presence_fold_assignments.csv'))
jar <- system.file('java/maxent.jar',package='dismo');stopifnot(file.exists(jar),dismo::maxent(silent=TRUE))
write_csv(data.frame(file=jar,md5=unname(tools::md5sum(jar)),bytes=file.info(jar)$size),file.path(directory,'maxent_jar.csv'))
base_args <- settings$maxent_args
fit_model <- function(p,b,path,jackknife=FALSE) {
 args <- base_args;args[grepl('^jackknife=',args)]<-paste0('jackknife=',tolower(jackknife))
 dismo::maxent(rbind(p,b),p=c(rep(1L,nrow(p)),rep(0L,nrow(b))),path=path,args=args)
}
predict_model <- function(model,x) as.numeric(predict(model,x,args=c('outputformat=logistic','doclamp=true')))
metrics<-list();roc<-list();held<-list()
for(f in seq_len(10)) {
 path<-file.path(directory,paste0('fold_',sprintf('%02d',f)));rds<-paste0(path,'.rds')
 if(file.exists(rds)){model<-readRDS(rds);log('Reusing completed fold',f)} else {
  log('Fitting fold',f);model<-fit_model(presence[pfold!=f,,drop=FALSE],background[bfold!=f,,drop=FALSE],path);saveRDS(model,rds)
 }
 p<-predict_model(model,presence[pfold==f,,drop=FALSE]);b<-predict_model(model,background[bfold==f,,drop=FALSE]);ev<-dismo::evaluate(p=p,a=b)
 metrics[[f]]<-data.frame(fold=f,n_presence=length(p),n_background=length(b),AUC=ev@auc)
 roc[[f]]<-data.frame(fold=f,threshold=ev@t,FPR=1-ev@TNR,TPR=ev@TPR)
 held[[f]]<-rbind(data.frame(fold=f,type='presence',suitability=p),data.frame(fold=f,type='background',suitability=b))
 write_csv(do.call(rbind,metrics),file.path(directory,'cross_validation.csv'))
}
write_csv(do.call(rbind,roc),file.path(directory,'cross_validation_ROC.csv'));write_csv(do.call(rbind,held),file.path(directory,'heldout_predictions.csv'))
model_path<-file.path(directory,'final_model.rds')
if(file.exists(model_path)){final<-readRDS(model_path);log('Reusing final model')}else{log('Fitting final model and required jackknife');final<-fit_model(presence,background,file.path(directory,'final_maxent'),TRUE);saveRDS(final,model_path)}
train_p<-predict_model(final,presence);train_b<-predict_model(final,background);ev<-dismo::evaluate(p=train_p,a=train_b)
threshold<-as.numeric(dismo::threshold(ev,'spec_sens'));stopifnot(length(threshold)==1,is.finite(threshold))
write_csv(data.frame(method='maximum training sensitivity plus specificity',threshold=threshold,training_AUC=ev@auc),file.path(directory,'threshold.csv'))
write_csv(data.frame(threshold=ev@t,FPR=1-ev@TNR,TPR=ev@TPR),file.path(directory,'final_training_ROC.csv'))
write_csv(data.frame(metric=rownames(final@results),value=final@results[,1]),file.path(directory,'maxent_final_results_long.csv'))
# Marginal conditional responses retain the final model, averaging other predictors at presence means.
response <- list();means<-colMeans(presence)
for(v in predictors){z<-as.data.frame(matrix(rep(means,each=121),nrow=121));names(z)<-predictors;z[[v]]<-seq(min(background[[v]],presence[[v]]),max(background[[v]],presence[[v]]),length.out=121);response[[v]]<-data.frame(predictor=v,value=z[[v]],suitability=predict_model(final,z))}
write_csv(do.call(rbind,response),file.path(directory,'conditional_response_curves.csv'))
# Chunk predictions avoid a full dense pair matrix; 1.15 million cells × 6 predictors only.
project <- function(x,path) {
 if(file.exists(path)){log('Reusing completed projection',basename(path));return(rast(path))}
 vals<-values(x,mat=TRUE);out<-rep(NA_real_,nrow(vals));idx<-which(complete.cases(vals))
 for(start in seq.int(1,length(idx),by=50000L)){ii<-idx[start:min(start+49999L,length(idx))];out[ii]<-predict_model(final,as.data.frame(vals[ii,,drop=FALSE]))}
 stopifnot(all(is.finite(out[idx])),all(out[idx]>=0&out[idx]<=1))
 r<-setValues(template,out);names(r)<-'suitability';writeRaster(r,path,overwrite=TRUE,wopt=list(gdal='COMPRESS=DEFLATE'));r
}
current_prediction<-project(current,file.path(directory,'current_suitability.tif'))
writeRaster(current_prediction>=threshold,file.path(directory,'current_suitable_MTSS.tif'),overwrite=TRUE)
train_range<-apply(rbind(presence,background),2,range)
novelty <- function(x,path){v<-values(x,mat=TRUE);count<-rowSums(sweep(v,2,train_range[1,],'<')|sweep(v,2,train_range[2,],'>'));writeRaster(setValues(template,count),path,overwrite=TRUE,wopt=list(gdal='COMPRESS=DEFLATE'))}
area<-values(cellSize(template,unit='km'),mat=FALSE);cur<-values(current_prediction,mat=FALSE);cur_binary<-cur>=threshold
sums<-list();pop_predictions<-list();gcm_names<-c('GFDL-ESM4','IPSL-CM6A-LR','MPI-ESM1-2-HR','MRI-ESM2-0','UKESM1-0-LL')
for(ssp in c('ssp126','ssp585'))for(period in c('2041-2070','2071-2100')) {
 ensemble_paths<-character()
 for(gcm in gcm_names){log('Projecting',ssp,period,gcm);prefix<-file.path(directory,paste(ssp,period,gcm,sep='_'))
  future<-rast(file.path(root,'public_dependencies/climate_individual/CHELSA_v2.1_regional',period,gcm,ssp,paste0(bio_names,'.tif')));names(future)<-bio_names
  future<-c(resample(future,template,method='bilinear'),elevation);names(future)<-predictors
  pred<-project(future,paste0(prefix,'_suitability.tif'));ensemble_paths<-c(ensemble_paths,paste0(prefix,'_suitability.tif'))
  novelty(future,paste0(prefix,'_novel_predictor_count.tif'))
  fp<-values(pred,mat=FALSE);valid_pair<-is.finite(cur)&is.finite(fp);suitable<-fp>=threshold
  sums[[length(sums)+1]]<-data.frame(ssp=ssp,period=period,gcm=gcm,threshold=threshold,current_area_km2=sum(area[cur_binary&valid_pair]),future_area_km2=sum(area[suitable&valid_pair]),stable_area_km2=sum(area[cur_binary&suitable&valid_pair]),loss_area_km2=sum(area[cur_binary&!suitable&valid_pair]),gain_area_km2=sum(area[!cur_binary&suitable&valid_pair]),shared_valid_area_km2=sum(area[valid_pair]))
  pop_predictions[[length(pop_predictions)+1]]<-cbind(pop_meta,ssp=ssp,period=period,gcm=gcm,suitability=terra::extract(pred,pop_meta[,c('longitude','latitude')],ID=FALSE)[,1])
  write_csv(do.call(rbind,sums),file.path(directory,'habitat_area_per_GCM.csv'))
 }
 stack<-rast(ensemble_paths);ensemble_prefix<-file.path(directory,paste(ssp,period,'ENSEMBLE',sep='_'))
 for(stat in c('mean','sd','min','max')){r<-app(stack,fun=get(stat));writeRaster(r,paste0(ensemble_prefix,'_',stat,'.tif'),overwrite=TRUE,wopt=list(gdal='COMPRESS=DEFLATE'))}
 en<-rast(paste0(ensemble_prefix,'_mean.tif'));writeRaster(en>=threshold,paste0(ensemble_prefix,'_MTSS.tif'),overwrite=TRUE)
 writeRaster((current_prediction>=threshold)+2*(en>=threshold),paste0(ensemble_prefix,'_habitat_transition.tif'),overwrite=TRUE)
}
write_csv(do.call(rbind,pop_predictions),file.path(directory,'population_future_suitability.csv'))
write_csv(cbind(pop_meta,suitability=terra::extract(current_prediction,pop_meta[,c('longitude','latitude')],ID=FALSE)[,1]),file.path(directory,'population_current_suitability.csv'))
save_session(directory)
marker(directory,list(analysis='Explicitly authorized missing-source SDM reconstruction',original_occurrences=nrow(occ),retained_occurrences=nrow(presence),background_cells=nrow(background),folds=10,future_GCM_scenarios=20,predictors=predictors,threshold_MTSS=threshold,mean_heldout_AUC=mean(do.call(rbind,metrics)$AUC),scope='Reconstructed numerical source, not historical model recovery',population_class_labels_changed=FALSE,figure_rendering_pending=TRUE))
log('SDM reconstruction completed and validated.')
