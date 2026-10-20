#!/usr/bin/env Rscript
# All additive driver subsets: no predictor is mandatory. Run one response case per Slurm array task.
suppressPackageStartupMessages(library(mgcv))
args <- commandArgs(TRUE)
if(length(args)<3) stop('Usage: fit_all_subsets.R input.csv output_dir case_index [linear|gam|both]')
input <- normalizePath(args[1]); outroot <- args[2]; case_index <- as.integer(args[3]); family_arg <- if(length(args)>3)args[4] else 'both'
benchmark <- Sys.getenv('AMMO_SUBSET_BENCHMARK','0')=='1'
resume <- Sys.getenv('AMMO_SUBSET_RESUME','1')=='1'
families <- switch(family_arg, linear='linear', gam='gam', both=c('linear','gam'),stop('Unknown family'))
dir.create(outroot,recursive=TRUE,showWarnings=FALSE)
dat <- read.csv(input,check.names=FALSE,stringsAsFactors=FALSE)
meta <- c('case_id','population','response','block')
if(!all(meta %in% names(dat)))stop('Required input columns: case_id,population,response,block and numeric predictors')
predictors <- setdiff(names(dat),c(meta,'bio','ssp','period','marker_set'))
if(!length(predictors)||length(predictors)>15)stop('Unexpected predictor count')
if(!all(vapply(dat[predictors],is.numeric,logical(1))))stop('All drivers must be numeric')
cases <- sort(unique(dat$case_id)); if(case_index<1||case_index>length(cases))stop('Bad case index')
case_id <- cases[case_index]; d <- dat[dat$case_id==case_id,,drop=FALSE]; d <- d[order(d$population),]
if(nrow(d)!=19||anyDuplicated(d$population))stop('Exactly19 unique populations required')
if(any(!is.finite(as.matrix(d[c('response',predictors)]))))stop('Nonfinite data; no model-specific complete-case deletion allowed')
if(length(unique(d$block))!=3)stop('Require the three fixed previously validated geographic blocks')
out <- file.path(outroot,gsub('[^A-Za-z0-9_-]','_',case_id));dir.create(out,showWarnings=FALSE,recursive=TRUE)
if(resume && file.exists(file.path(out,'validation.json')) && !benchmark){saved<-read.csv(file.path(out,'case_input.csv'),check.names=FALSE);if(!isTRUE(all.equal(saved,d,check.attributes=FALSE)))stop('Completed case inputs changed; use a new output directory');cat('SKIP_VALIDATED',case_id,'\n');quit(status=0)}
write.csv(d,file.path(out,'case_input.csv'),row.names=FALSE)
# Prespecified biological/original-factor order: all original10 before optional Pi.
# This filter never sees response values and is recomputed independently in training folds.
screen_priority<-c('climate_change','Ho','He','DEL','LoF','private_alleles','Fis','Alt','Lat','Lon','Pi')
if(!setequal(screen_priority,predictors))stop('Unexpected predictor universe for fixed correlation-screen priority')
correlation_screen<-function(train,scope,fold){
  retained<-character();decisions<-list()
  for(v in screen_priority){
    rhos<-if(length(retained))vapply(retained,function(k)suppressWarnings(cor(train[[v]],train[[k]],method='spearman')),numeric(1)) else numeric()
    blocked<-names(rhos)[is.finite(rhos)&abs(rhos)>.8]
    keep<-!length(blocked)
    decisions[[length(decisions)+1L]]<-data.frame(case_id,scope,fold=as.character(fold),n_training=nrow(train),predictor=v,priority=match(v,screen_priority),retained=keep,blocking_predictors=paste(blocked,collapse=';'),blocking_Spearman_rho=paste(sprintf('%.17g',rhos[blocked]),collapse=';'),rule='absolute Spearman rho >0.8; fixed original-factor/biological priority; response not used')
    if(keep)retained<-c(retained,v)
  }
  list(retained=retained,decisions=do.call(rbind,decisions))
}
fullscreen<-correlation_screen(d,'full','all');screenlogs<-list(fullscreen$decisions)
subsets <- lapply(0:(2^length(predictors)-1),function(mask)predictors[as.logical(intToBits(mask)[seq_along(predictors)])])
ids <- sprintf('M%0*d',nchar(length(subsets)),seq_along(subsets)-1L)
universe <- data.frame(model_id=ids,n_predictors=lengths(subsets),predictors=vapply(subsets,paste,collapse=' + ',FUN.VALUE=''),stringsAsFactors=FALSE)
universe$full_screen_eligible<-vapply(subsets,function(v)all(v%in%fullscreen$retained),logical(1))
universe$full_screen_retained<-universe$full_screen_eligible
write.csv(universe,file.path(out,'all_subset_universe.csv'),row.names=FALSE)
if(benchmark){keep <- unique(c(1L,2L,which(lengths(subsets)==2)[1:3],which(lengths(subsets)==4)[1:2],which(lengths(subsets)==7)[1:2],length(subsets)));keep<-keep[is.finite(keep)];subsets<-subsets[keep];ids<-ids[keep];universe<-universe[keep,];out<-file.path(out,'benchmark');dir.create(out,showWarnings=FALSE)}
# ML makes differing fixed/smooth structures comparable. REML scores are not used to select different subsets.
fit_one <- function(train,test,vars,family,details=FALSE,allowed_predictors=predictors){
  warn <- character();fail <- function(why)list(status=why,aicc=NA_real_,aic=NA_real_,K=NA_real_,edf=NA_real_,rank=NA_integer_,condition=NA_real_,pred=rep(NA_real_,nrow(test)),warning=paste(unique(warn),collapse=' | '))
  if(any(!vars%in%allowed_predictors))return(fail('correlation_screen_excluded_predictor'))
  if(length(vars)){
    mu <- vapply(train[vars],mean,numeric(1));sdv <- vapply(train[vars],sd,numeric(1))
    if(any(!is.finite(sdv)|sdv<1e-12))return(fail('constant_predictor_in_training'))
    train[vars] <- as.data.frame(sweep(sweep(as.matrix(train[vars]),2,mu,'-'),2,sdv,'/'))
    test[vars] <- as.data.frame(sweep(sweep(as.matrix(test[vars]),2,mu,'-'),2,sdv,'/'))
  }
  X <- if(length(vars))model.matrix(reformulate(vars),train) else matrix(1,nrow(train),1); qrX <- qr(X,tol=1e-9)
  if(qrX$rank<ncol(X))return(fail('rank_deficient_linear_predictor_matrix'))
  cond <- kappa(X,exact=TRUE)
  if(family=='gam' && length(vars) && (1+2*length(vars))>nrow(train)-1)return(fail('smooth_basis_exceeds_training_information'))
  rhs <- if(!length(vars))'1' else if(family=='gam')paste(sprintf('s(%s,k=3,bs="tp")',vars),collapse=' + ') else paste(vars,collapse=' + ')
  form <- as.formula(paste('response ~',rhs))
  fit <- tryCatch(withCallingHandlers(if(family=='gam'&&length(vars))gam(form,data=train,method='ML',select=TRUE) else lm(form,data=train),warning=function(w){warn<<-c(warn,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)e)
  if(inherits(fit,'error'))return(fail(paste0('fit_error: ',conditionMessage(fit))))
  if(inherits(fit,'gam')&&!isTRUE(fit$converged))return(fail('gam_not_converged'))
  ll <- logLik(fit); K <- as.numeric(attr(ll,'df')); aic <- -2*as.numeric(ll)+2*K
  if(!is.finite(aic)||nrow(train)-K-1<=0)return(fail('AICc_undefined_at_effective_parameter_count'))
  pr <- tryCatch(as.numeric(predict(fit,newdata=test)),error=function(e)rep(NA_real_,nrow(test)))
  if(any(!is.finite(pr)))return(fail('nonfinite_prediction'))
  ans <- list(status='eligible',aicc=aic+2*K*(K+1)/(nrow(train)-K-1),aic=aic,K=K,edf=if(inherits(fit,'gam'))sum(fit$edf) else fit$rank,rank=qrX$rank,condition=cond,pred=pr,warning=paste(unique(warn),collapse=' | '))
  if(details){ans$fit <- fit;ans$standardization<-if(length(vars))data.frame(predictor=vars,training_mean=mu,training_SD=sdv) else data.frame(predictor=character(),training_mean=numeric(),training_SD=numeric())}
  ans
}
rows <- list();selectedrows<-list();diagrows<-list();perfs<-list();foldaudit<-list();counter<-0L;sc<-0L;dc<-0L;fc<-0L
for(family in families){
  fullresults <- vector('list',length(subsets));foldpred <- matrix(NA_real_,nrow(d),length(subsets))
  for(i in seq_along(subsets)){
    f <- fit_one(d,d,subsets[[i]],family,allowed_predictors=fullscreen$retained);fullresults[[i]]<-f
    counter<-counter+1L; rows[[counter]]<-data.frame(case_id,family,model_id=ids[i],scope='full',fold='all',n_train=nrow(d),status=f$status,aicc=f$aicc,aic=f$aic,K=f$K,edf=f$edf,rank=f$rank,condition=f$condition,near_collinear=is.finite(f$condition)&&f$condition>30,warning=f$warning)
  }
  for(cv in c('LOPO','geographic_block')){
    sse<-sae<-nvalid<-rep(0,length(subsets));baseline_sse<-0;expected<-0
    folds <- if(cv=='LOPO')seq_len(nrow(d)) else sort(unique(d$block))
    if(benchmark)folds<-folds[1]
    for(fold in folds){
      testids <- if(cv=='LOPO')fold else which(d$block==fold);train<-d[-testids,,drop=FALSE];test<-d[testids,,drop=FALSE];res<-vector('list',length(subsets))
      foldscreen<-correlation_screen(train,cv,fold);if(family==families[1])screenlogs[[length(screenlogs)+1L]]<-foldscreen$decisions
      for(i in seq_along(subsets)){
        f<-fit_one(train,test,subsets[[i]],family,allowed_predictors=foldscreen$retained);res[[i]]<-f
        if(all(is.finite(f$pred))){nvalid[i]<-nvalid[i]+nrow(test);sse[i]<-sse[i]+sum((f$pred-test$response)^2);sae[i]<-sae[i]+sum(abs(f$pred-test$response))}
        if(cv=='LOPO')foldpred[testids,i]<-f$pred
      }
      baseline_sse<-baseline_sse+sum((mean(train$response)-test$response)^2);expected<-expected+nrow(test)
      st<-table(vapply(res,`[[`,character(1),'status'));fc<-fc+1L;foldaudit[[fc]]<-data.frame(case_id,family,scope=cv,fold=as.character(fold),status=names(st),n_models=as.integer(st))
      scores<-vapply(res,function(x)if(x$status=='eligible')x$aicc else Inf,numeric(1));winner<-which.min(scores);if(!any(is.finite(scores)))stop('No eligible model including intercept-only in training fold')
      sc<-sc+1L;selectedrows[[sc]]<-data.frame(case_id,family,scope=cv,fold=as.character(fold),model_id=ids[winner],predictors=universe$predictors[winner],n_eligible=sum(is.finite(scores)),population=test$population,observed=test$response,predicted=res[[winner]]$pred,training_mean=mean(train$response),training_AICc=scores[winner],selection='lowest_training_AICc_then_model_id')
    }
    pp<-data.frame(case_id,family,scope=cv,model_id=ids,n=nvalid,n_expected=expected,RMSE=ifelse(nvalid==expected,sqrt(sse/expected),NA_real_),MAE=ifelse(nvalid==expected,sae/expected,NA_real_),baseline_RMSE=sqrt(baseline_sse/expected));perfs[[length(perfs)+1L]]<-pp
    # Durable checkpoint after each family and validation scope.
    write.csv(do.call(rbind,rows),file.path(out,'model_fit_audit.csv'),row.names=FALSE)
    write.csv(do.call(rbind,perfs),file.path(out,'all_model_validation_performance.csv'),row.names=FALSE)
    write.csv(do.call(rbind,foldaudit),file.path(out,'heldout_fit_eligibility_counts.csv'),row.names=FALSE)
    write.csv(do.call(rbind,selectedrows),file.path(out,'nested_selected_heldout_predictions.csv'),row.names=FALSE)
  }
  scores<-vapply(fullresults,function(x)if(x$status=='eligible')x$aicc else Inf,numeric(1));best<-which.min(scores)
  ord<-order(scores);top<-ord[seq_len(min(10,sum(is.finite(scores))))]
  weights<-exp(-.5*(scores-min(scores)));weights[!is.finite(weights)]<-0;weights<-weights/sum(weights)
  selection<-universe;selection$status<-vapply(fullresults,`[[`,character(1),'status');selection$case_id<-case_id;selection$family<-family;selection$aicc<-scores;selection$delta_aicc<-scores-min(scores);selection$akaike_weight<-weights
  selection$LOPO_RMSE<-apply(foldpred,2,function(p)if(all(is.finite(p)))sqrt(mean((p-d$response)^2)) else NA_real_)
  write.csv(selection,file.path(out,paste0(family,'_full_model_ranking.csv')),row.names=FALSE)
  importance<-data.frame(case_id,family,predictor=predictors,full_screen_retained=predictors%in%fullscreen$retained,summed_Akaike_weight=vapply(predictors,function(v)if(v%in%fullscreen$retained)sum(weights[vapply(subsets,function(s)v%in%s,logical(1))]) else NA_real_,numeric(1)),interpretation=ifelse(predictors%in%fullscreen$retained,'relative support within screened candidate set; not causal importance','screened out for redundancy; importance not estimable'))
  write.csv(importance,file.path(out,paste0(family,'_inclusion_weights.csv')),row.names=FALSE)
  for(i in top){
    f<-fit_one(d,d,subsets[[i]],family,details=TRUE,allowed_predictors=fullscreen$retained);fit<-f$fit
    write.csv(f$standardization,file.path(out,paste0(family,'_',ids[i],'_predictor_standardization.csv')),row.names=FALSE)
    beta<-coef(fit);vv<-vcov(fit);se<-sqrt(pmax(0,diag(vv)));rdf<-df.residual(fit);critical<-qt(.975,df=rdf)
    ci<-data.frame(parameter=names(beta),estimate=as.numeric(beta),SE=as.numeric(se),lower95=as.numeric(beta-critical*se),upper95=as.numeric(beta+critical*se),residual_df=rdf,interval='95% pointwise conditional on selected formula; no model-search correction',parameter_type=if(inherits(fit,'gam'))'GAM basis coefficient on standardized predictors; not a linear driver effect' else 'linear coefficient on standardized predictors')
    write.csv(ci,file.path(out,paste0(family,'_',ids[i],'_conditional_parameter_intervals.csv')),row.names=FALSE)
    saveRDS(fit,file.path(out,paste0(family,'_',ids[i],'_conditional_fit.rds')))
    hat<-tryCatch(hatvalues(fit),error=function(e)fit$hat)
    cook<-tryCatch(cooks.distance(fit),error=function(e)rep(NA_real_,nrow(d)))
    dc<-dc+1L;diagrows[[dc]]<-data.frame(case_id,family,model_id=ids[i],population=d$population,residual=as.numeric(residuals(fit)),fitted=as.numeric(fitted(fit)),leverage=hat,Cooks_D=cook)
    if(length(subsets[[i]])>1){X<-model.matrix(reformulate(subsets[[i]]),d);C<-cor(X[,-1,drop=FALSE]);vv<-tryCatch(diag(solve(C)),error=function(e)rep(NA_real_,ncol(C)));write.csv(data.frame(predictor=colnames(C),linear_VIF=vv),file.path(out,paste0(family,'_',ids[i],'_linear_predictor_VIF.csv')),row.names=FALSE)}
    if(inherits(fit,'gam')){write.csv(summary(fit)$s.table,file.path(out,paste0(family,'_',ids[i],'_conditional_smooth_tests.csv')));write.csv(concurvity(fit,full=TRUE),file.path(out,paste0(family,'_',ids[i],'_concurvity.csv')))} else write.csv(coef(summary(fit)),file.path(out,paste0(family,'_',ids[i],'_conditional_coefficients.csv')))
  }
}
write.csv(do.call(rbind,diagrows),file.path(out,'top_models_conditional_diagnostics.csv'),row.names=FALSE)
sel<-do.call(rbind,selectedrows)
sel$selected_family<-sel$family
# Joint family/subset selection uses training AICc only, never outer held-out errors.
# Within-family minimization followed by minimization of those two scores equals
# minimization over the full union, with linear preferred for exact-score ties.
if(all(c('linear','gam')%in%families)){
  combined<-do.call(rbind,lapply(split(sel,interaction(sel$scope,sel$fold,drop=TRUE)),function(z){
    fam<-unique(z[c('family','training_AICc','model_id')]);fam<-fam[order(fam$training_AICc,match(fam$family,c('linear','gam')),fam$model_id),]
    chosen<-z[z$family==fam$family[1],,drop=FALSE];chosen$selected_family<-chosen$family;chosen$family<-'combined';chosen$selection<-'lowest_training_AICc_across_all_family_subsets_then_linear_tie';chosen
  }))
  sel<-rbind(sel,combined)
}
write.csv(sel,file.path(out,'nested_selected_heldout_predictions.csv'),row.names=FALSE)
# Combined full-data ranking is descriptive. Family-specific Akaike weights remain
# unchanged; no joint weights are claimed for the duplicate null model.
fullrank<-do.call(rbind,lapply(families,function(fm)read.csv(file.path(out,paste0(fm,'_full_model_ranking.csv')),check.names=FALSE)))
fullrank<-fullrank[order(fullrank$aicc,match(fullrank$family,c('linear','gam')),fullrank$model_id),]
fullrank$duplicate_null<-fullrank$n_predictors==0 & fullrank$family=='gam'
write.csv(fullrank,file.path(out,'combined_full_model_ranking.csv'),row.names=FALSE)
metric<-function(z)data.frame(n=sum(is.finite(z$predicted)),n_expected=nrow(z),RMSE=if(all(is.finite(z$predicted)))sqrt(mean((z$predicted-z$observed)^2)) else NA_real_,MAE=if(all(is.finite(z$predicted)))mean(abs(z$predicted-z$observed)) else NA_real_,baseline_RMSE=sqrt(mean((z$training_mean-z$observed)^2)))
perfsel<-do.call(rbind,lapply(split(sel,interaction(sel$family,sel$scope,drop=TRUE)),function(z)cbind(z[1,c('case_id','family','scope')],metric(z))))
write.csv(perfsel,file.path(out,'nested_selection_validation_performance.csv'),row.names=FALSE)
screenall<-do.call(rbind,screenlogs)
screensummary<-do.call(rbind,lapply(split(screenall,interaction(screenall$scope,screenall$fold,drop=TRUE)),function(z)data.frame(case_id,scope=z$scope[1],fold=z$fold[1],n_training=z$n_training[1],n_retained=sum(z$retained),n_screen_eligible_subsets=2^sum(z$retained),retained_predictors=paste(z$predictor[z$retained],collapse=';'),excluded_predictors=paste(z$predictor[!z$retained],collapse=';'))))
write.csv(screensummary,file.path(out,'correlation_screen_summary.csv'),row.names=FALSE)
write.csv(screenall,file.path(out,'correlation_screen_decisions.csv'),row.names=FALSE)
write.csv(cor(d[predictors],method='spearman'),file.path(out,'predictor_Spearman_correlation.csv'))
write.csv(cor(d[predictors]),file.path(out,'predictor_correlation.csv'))
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
writeLines(sprintf('{"validated":true,"case_id":"%s","n_populations":19,"n_predictors":%s,"n_subsets":%s,"families":"%s","validation":"conditional_upstream_RONA_fixed_nested_downstream_AICc_selection"}',case_id,length(predictors),length(subsets),paste(families,collapse=',')),file.path(out,'validation.json'))
cat('COMPLETE',case_id,'\n')
