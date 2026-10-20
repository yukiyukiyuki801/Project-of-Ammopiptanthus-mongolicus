#!/usr/bin/env Rscript
# Read saved selected models only: this script never fits or selects additional models.
suppressPackageStartupMessages(library(mgcv))
suppressPackageStartupMessages(library(geosphere))
args<-commandArgs(TRUE);if(length(args)<2)stop('Usage: postprocess_driver_models.R model_root output_dir [expected_cases=20]')
root<-normalizePath(args[1]);out<-args[2];expected_cases<-if(length(args)>2)as.integer(args[3]) else 20L
dir.create(out,recursive=TRUE,showWarnings=FALSE)
casedirs<-sort(list.dirs(root,recursive=FALSE,full.names=TRUE));casedirs<-casedirs[file.exists(file.path(casedirs,'validation.json'))]
if(length(casedirs)!=expected_cases)stop('Expected ',expected_cases,' validated cases; found ',length(casedirs),'. No partial finalization.')
read<-function(p)read.csv(p,check.names=FALSE,stringsAsFactors=FALSE)
aggregate_files<-c('model_fit_audit.csv','all_model_validation_performance.csv','nested_selected_heldout_predictions.csv','nested_selection_validation_performance.csv','correlation_screen_decisions.csv','correlation_screen_summary.csv','heldout_fit_eligibility_counts.csv')
for(fn in aggregate_files){z<-do.call(rbind,lapply(casedirs,function(cd)read(file.path(cd,fn))));write.csv(z,file.path(out,fn),row.names=FALSE)}
ranks<-do.call(rbind,lapply(casedirs,function(cd)read(file.path(cd,'combined_full_model_ranking.csv'))));write.csv(ranks,file.path(out,'all_full_model_rankings.csv'),row.names=FALSE)
inclusions<-do.call(rbind,lapply(casedirs,function(cd)do.call(rbind,lapply(c('linear','gam'),function(fm)read(file.path(cd,paste0(fm,'_inclusion_weights.csv')))))));write.csv(inclusions,file.path(out,'all_inclusion_weights.csv'),row.names=FALSE)
# Moran definition follows the existing revision diagnostic: binary symmetric union3NN.
moran_test<-function(e,W,seed){
  e<-e-mean(e);n<-length(e);fun<-function(z)if(sum(z*z)<1e-20)NA_real_ else n/sum(W)*as.numeric(crossprod(z,W%*%z))/sum(z*z)
  observed<-fun(e);set.seed(seed);null<-replicate(1999,fun(sample(e)));expected<- -1/(n-1)
  p<-if(is.finite(observed))(1+sum(abs(null-expected)>=abs(observed-expected),na.rm=TRUE))/2000 else NA_real_
  c(I=observed,p=p)
}
diags<-list();smooths<-list();concur<-list();curves<-list();points<-list();intervals<-list();tops<-list();dd<-0L;ss<-0L;cc<-0L;gg<-0L;pp<-0L;ii<-0L;tt<-0L
for(ci in seq_along(casedirs)){
 cd<-casedirs[ci];d<-read(file.path(cd,'case_input.csv'));case<-d$case_id[1]
 dist<-geosphere::distm(as.matrix(d[c('Lon','Lat')]));W<-matrix(0,nrow(d),nrow(d));for(j in seq_len(nrow(d))){v<-dist[j,];v[j]<-Inf;W[j,order(v)[1:3]]<-1};W<-pmax(W,t(W))
 for(fm in c('linear','gam')){
  rank<-read(file.path(cd,paste0(fm,'_full_model_ranking.csv')));rank<-rank[rank$status=='eligible'&is.finite(rank$aicc),];rank<-rank[order(rank$aicc,rank$model_id),];if(!nrow(rank))stop('No eligible full model')
  win<-rank[1,];id<-win$model_id;prefix<-paste0(fm,'_',id);fit<-readRDS(file.path(cd,paste0(prefix,'_conditional_fit.rds')));std<-read(file.path(cd,paste0(prefix,'_predictor_standardization.csv')))
  tt<-tt+1L;tops[[tt]]<-win
  e<-as.numeric(residuals(fit));h<-tryCatch(as.numeric(hatvalues(fit)),error=function(err)as.numeric(fit$hat));if(length(h)!=nrow(d))h<-rep(NA_real_,nrow(d))
  cooks<-tryCatch(as.numeric(cooks.distance(fit)),error=function(err)rep(NA_real_,nrow(d)));if(length(cooks)!=nrow(d))cooks<-rep(NA_real_,nrow(d))
  m<-moran_test(e,W,1061600+ci*2+match(fm,c('linear','gam')));sm<-summary(fit)
  dd<-dd+1L;diags[[dd]]<-data.frame(case_id=case,family=fm,model_id=id,n=nrow(d),n_predictors=win$n_predictors,total_EDF=if(inherits(fit,'gam'))sum(fit$edf) else fit$rank,residual_EDF=df.residual(fit),adjusted_R2=if(inherits(fit,'gam'))sm$r.sq else sm$adj.r.squared,deviance_explained=if(inherits(fit,'gam'))sm$dev.expl else 1-sum(e^2)/sum((d$response-mean(d$response))^2),max_leverage=if(any(is.finite(h)))max(h,na.rm=TRUE) else NA_real_,max_Cooks_D=if(any(is.finite(cooks)))max(cooks,na.rm=TRUE) else NA_real_,Moran_I=m['I'],Moran_permutation_p=m['p'],Moran_permutations=1999,spatial_weights='Symmetric union of3nearest geographic neighbours; binary weights',interpretation='Postselection conditional diagnostic; low power at19populations; not proof of independence')
  int<-read(file.path(cd,paste0(prefix,'_conditional_parameter_intervals.csv')));ii<-ii+1L;intervals[[ii]]<-cbind(case_id=case,family=fm,model_id=id,int)
  if(inherits(fit,'gam')){
   set.seed(2061600+ci);kc<-mgcv::k.check(fit,subsample=5000,n.rep=400)
   st<-as.data.frame(sm$s.table);st$smooth<-rownames(st);rownames(st)<-NULL
   kk<-as.data.frame(kc);kk$smooth<-rownames(kk);rownames(kk)<-NULL
   joined<-merge(st,kk,by='smooth',all=TRUE);ss<-ss+1L;smooths[[ss]]<-cbind(case_id=case,family=fm,model_id=id,joined)
   cv<-mgcv::concurvity(fit,full=TRUE);cj<-expand.grid(diagnostic=rownames(cv),term=colnames(cv),stringsAsFactors=FALSE);cj$value<-as.vector(cv);cc<-cc+1L;concur[[cc]]<-cbind(case_id=case,model_id=id,cj)
   vars<-std$predictor;standardized<-d
   for(v in vars){q<-std[std$predictor==v,];standardized[[v]]<-(d[[v]]-q$training_mean)/q$training_SD}
   for(v in vars){
    q<-std[std$predictor==v,];grid<-seq(min(d[[v]]),max(d[[v]]),length.out=100);nd<-as.data.frame(matrix(0,nrow=100,ncol=length(vars),dimnames=list(NULL,vars)));nd[[v]]<-(grid-q$training_mean)/q$training_SD
    pr<-predict(fit,newdata=nd,type='response',se.fit=TRUE,unconditional=TRUE);critical<-qt(.975,df=df.residual(fit));ref<-as.numeric(predict(fit,newdata=nd[1,,drop=FALSE]*0))
    gg<-gg+1L;curves[[gg]]<-data.frame(case_id=case,family=fm,model_id=id,predictor=v,x_original=grid,conditional_mean=as.numeric(pr$fit),SE=as.numeric(pr$se.fit),lower95=as.numeric(pr$fit-critical*pr$se.fit),upper95=as.numeric(pr$fit+critical*pr$se.fit),difference_from_predictor_mean=as.numeric(pr$fit)-ref,other_predictors='held at full-case means; standardized zero',interval='pointwise95% mean-response CI; includes available smoothing uncertainty; conditional on selected formula, no search correction')
    pd<-as.data.frame(matrix(0,nrow=nrow(d),ncol=length(vars),dimnames=list(NULL,vars)));pd[[v]]<-standardized[[v]];pfit<-as.numeric(predict(fit,newdata=pd));pp<-pp+1L;points[[pp]]<-data.frame(case_id=case,model_id=id,predictor=v,population=d$population,x_original=d[[v]],observed_RONA=d$response,conditional_partial_residual=pfit+e,conditional_mean_at_x=pfit,residual=e)
   }
  }
 }
}
write_list<-function(x,fn){if(length(x))write.csv(do.call(rbind,x),file.path(out,fn),row.names=FALSE) else write.csv(data.frame(),file.path(out,fn),row.names=FALSE)}
topall<-do.call(rbind,tops);topcombined<-do.call(rbind,lapply(split(topall,topall$case_id),function(z)z[order(z$aicc,match(z$family,c('linear','gam')),z$model_id)[1],,drop=FALSE]));write.csv(topcombined,file.path(out,'top1_combined_full_models.csv'),row.names=FALSE)
write_list(tops,'top1_each_family_full_models.csv');write_list(diags,'top1_conditional_diagnostics.csv');write_list(intervals,'top1_conditional_parameter_intervals.csv');write_list(smooths,'top1_GAM_smooth_EDF_kcheck.csv');write_list(concur,'top1_GAM_concurvity.csv');write_list(curves,'top1_GAM_conditional_curve_source.csv');write_list(points,'top1_GAM_partial_residual_source.csv')
# Selection stability uses folds as units, not duplicate geographic held-out populations.
sel<-read(file.path(out,'nested_selected_heldout_predictions.csv'));folds<-unique(sel[c('case_id','family','scope','fold','model_id','predictors','selected_family')]);write.csv(folds,file.path(out,'nested_selected_fold_formulas.csv'),row.names=FALSE)
screenall<-read(file.path(out,'correlation_screen_decisions.csv'))
stability<-list();j<-0L
for(key in split(folds,interaction(folds$case_id,folds$family,folds$scope,drop=TRUE))){
 vars<-names(d)[!names(d)%in%c('case_id','population','response','block')]
 for(v in vars){included<-vapply(strsplit(key$predictors,' + ',fixed=TRUE),function(z)v%in%z,logical(1));eligible<-screenall[screenall$case_id==key$case_id[1]&screenall$scope==key$scope[1]&screenall$predictor==v,];n_eligible<-sum(eligible$retained);j<-j+1L;stability[[j]]<-data.frame(case_id=key$case_id[1],family=key$family[1],scope=key$scope[1],predictor=v,n_selected_folds=sum(included),n_folds=nrow(key),n_screen_eligible_folds=n_eligible,selected_fraction_all_folds=if(n_eligible)mean(included) else NA_real_,selected_fraction_when_eligible=if(n_eligible)sum(included)/n_eligible else NA_real_,interpretation='Selection frequency under correlated candidate competition; NA if always excluded; not causal importance')}
}
write_list(stability,'nested_selection_predictor_stability.csv')
combinedfolds<-folds[folds$family=='combined',];familyfreq<-do.call(rbind,lapply(split(combinedfolds,interaction(combinedfolds$case_id,combinedfolds$scope,drop=TRUE)),function(z)data.frame(case_id=z$case_id[1],scope=z$scope[1],selected_family=c('linear','gam'),n_folds=c(sum(z$selected_family=='linear'),sum(z$selected_family=='gam')),total_folds=nrow(z),fraction=c(mean(z$selected_family=='linear'),mean(z$selected_family=='gam')))));write.csv(familyfreq,file.path(out,'nested_selected_family_stability.csv'),row.names=FALSE)
status<-read(file.path(out,'model_fit_audit.csv'));status_counts<-aggregate(list(n_models=status$model_id),status[c('case_id','family','status')],length);write.csv(status_counts,file.path(out,'full_model_eligibility_counts.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'));writeLines(sprintf('{"validated":true,"n_cases":%s,"no_new_model_fits":true,"top1_models":%s,"Moran_permutations":1999,"GAM_kcheck_replicates":400}',length(casedirs),length(diags)),file.path(out,'validation.json'))
cat('POSTPROCESS_COMPLETE',length(casedirs),'cases\n')
