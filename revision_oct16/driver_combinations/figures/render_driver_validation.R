#!/usr/bin/env Rscript
# Replace S13 from saved downstream nested all-subset validation; never fit models.
# R-only quantitative grid, native190x230mm, complete20cases, no forced climate baseline.
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(jsonlite)})
args <- commandArgs(TRUE)
if(length(args)<2) stop('Usage: Rscript render_driver_validation.R postprocess_dir output_dir [--synthetic-layout]')
in_dir <- args[1]; figure_dir <- args[2]; mock <- '--synthetic-layout' %in% args
root <- normalizePath(getwd()); dir.create(figure_dir,recursive=TRUE,showWarnings=FALSE)
skill <- Sys.getenv("AMMO_FIGURE_QA_DIR", unset="")
source(file.path(skill,'scripts/panel_alignment.R'));source('scripts/revision/measure_panels.R')
python <- Sys.getenv("AMMO_PYTHON", unset="python3")
scenarios <- c('ssp126_2041-2070','ssp126_2071-2100','ssp585_2041-2070','ssp585_2071-2100')
bios <- c('bio4','bio5','bio6','bio15','bio17')
keys <- unlist(lapply(scenarios,function(z)paste(z,bios,sep='_')))
row_labels <- setNames(unlist(lapply(c('126/M','126/L','585/M','585/L'),function(z)paste(z,toupper(bios)))),keys)
families <- c('linear','gam','combined');display <- c(linear='Linear subsets',gam='GAM subsets',combined='Combined subsets')
shapes <- c(linear=1,gam=2,combined=16)
drivers <- c('climate_change','Ho','He','DEL','LoF','private_alleles','Fis','Alt','Lat','Lon')
driver_labels <- c(climate_change='Climate change',Ho='Ho',He='He',DEL='DEL / SYN',LoF='LoF / SYN',private_alleles='Private alleles',Fis='FIS',Alt='Elevation',Lat='Latitude',Lon='Longitude')
read <- function(fn)read.csv(file.path(in_dir,fn),check.names=FALSE,stringsAsFactors=FALSE)
if(mock){
 performance <- expand.grid(case_id=keys,family=families,scope=c('LOPO','geographic_block'),stringsAsFactors=FALSE)
 performance$RMSE <- exp(seq(log(.15),log(4),length.out=nrow(performance)));performance$baseline_RMSE<-1;performance$n<-19;performance$n_expected<-19
 stability <- expand.grid(case_id=keys,family=families,scope='LOPO',predictor=c(drivers,'Pi'),stringsAsFactors=FALSE)
 stability$n_folds<-19;stability$n_screen_eligible_folds<-ifelse(stability$predictor=='Pi',0,19);stability$n_selected_folds<-ifelse(stability$predictor=='Pi',0,seq_len(nrow(stability))%%20)
} else {
 marker <- jsonlite::read_json(file.path(in_dir,'validation.json'),simplifyVector=TRUE)
 stopifnot(isTRUE(marker$validated),marker$n_cases==20)
 performance <- read('nested_selection_validation_performance.csv');stability <- read('nested_selection_predictor_stability.csv')
}
stopifnot(setequal(unique(performance$case_id),keys),setequal(unique(performance$family),families))
# Scope names are fixed by the analytical source. Fail rather than plot an unknown protocol.
stopifnot(setequal(unique(performance$scope),c('LOPO','geographic_block')))
performance$RMSE_ratio <- performance$RMSE/performance$baseline_RMSE
stopifnot(nrow(performance)==120,all(is.finite(performance$RMSE_ratio)),all(performance$RMSE_ratio>0),all(performance$n==19),all(performance$n_expected==19))
stopifnot(!anyDuplicated(performance[c('case_id','family','scope')]))
performance$family <- factor(performance$family,levels=families)
performance$y <- match(performance$case_id,rev(keys))+(as.integer(performance$family)-2)*.18
lopo <- performance[performance$scope=='LOPO',];block<-performance[performance$scope=='geographic_block',]
s <- stability[stability$scope=='LOPO',]
stopifnot(all(s$n_selected_folds>=0),all(s$n_screen_eligible_folds>=s$n_selected_folds),!anyDuplicated(s[c('case_id','family','predictor')]))
# A wholly excluded predictor is not assigned a zero importance score.
excluded <- aggregate(cbind(n_selected_folds,n_screen_eligible_folds)~predictor,s,sum)
if(any(!excluded$predictor%in%drivers & excluded$n_screen_eligible_folds>0)) stop('Additional driver is eligible in LOPO folds; return source for explicit panel-c update before excluding it')
s <- s[s$predictor%in%drivers,]
stopifnot(nrow(s)==20*3*10,setequal(unique(s$case_id),keys))
pool <- aggregate(cbind(n_selected_folds,n_screen_eligible_folds,n_folds)~predictor+family,s,sum)
pool$selected_fraction <- with(pool,ifelse(n_screen_eligible_folds>0,n_selected_folds/n_screen_eligible_folds,NA_real_))
pool$family<-factor(pool$family,levels=families);pool$y<-match(pool$predictor,rev(drivers))+(as.integer(pool$family)-2)*.16
stopifnot(nrow(pool)==30,all(pool$n_folds==380),all(is.finite(pool$selected_fraction)),all(pool$selected_fraction>=0&pool$selected_fraction<=1))
ratios <- performance$RMSE_ratio
limits <- 10^c(floor(log10(min(c(ratios,1)))),ceiling(log10(max(c(ratios,1)))))
if(limits[1]==limits[2])limits<-limits*c(.1,10)
breaks <- 10^seq(log10(limits[1]),log10(limits[2])); if(length(breaks)>5)breaks<-breaks[unique(round(seq(1,length(breaks),length.out=5)))]
base <- theme_classic(base_size=9,base_family='sans')+theme(axis.text=element_text(size=8.1,colour='black'),axis.title=element_text(size=8.1),legend.text=element_text(size=8.1),legend.title=element_blank(),legend.key.width=unit(5,'mm'),legend.key.height=unit(5,'mm'),plot.tag=element_text(size=11,face='bold',hjust=1,vjust=0),plot.tag.position=c(0,1),plot.margin=margin(8,7,7,7),legend.position='bottom',panel.grid=element_blank(),plot.background=element_rect(fill='white',colour=NA))
shape_scale <- function()scale_shape_manual(values=shapes,breaks=families,labels=display,drop=FALSE)
validation <- function(d,tag,show_y,xlab)ggplot(d,aes(RMSE_ratio,y,shape=family))+
 geom_vline(xintercept=1,linewidth=.3,linetype=2,colour='grey55')+geom_point(size=1.6,stroke=.4)+shape_scale()+
 scale_x_log10(limits=limits,breaks=breaks,labels=function(x)vapply(x,function(z)format(z,scientific=FALSE,trim=TRUE,digits=2),character(1)),expand=expansion(mult=c(.05,.05)))+
 scale_y_continuous(breaks=seq_along(rev(keys)),labels=if(show_y)row_labels[rev(keys)] else NULL,limits=c(.3,20.7),expand=expansion(mult=0))+
 labs(x=xlab,y=NULL,tag=tag)+base+guides(shape=guide_legend(nrow=1,byrow=TRUE,override.aes=list(size=2)))
pa<-validation(lopo,'a',TRUE,'LOPO RMSE /\nmean baseline')
pb<-validation(block,'b',FALSE,'Geographic RMSE /\nmean baseline')+theme(axis.ticks.y=element_blank())+guides(shape='none')
pc<-ggplot(pool,aes(selected_fraction,y,shape=family))+geom_point(size=1.6,stroke=.4)+shape_scale()+
 geom_text(data=data.frame(y=seq_along(rev(drivers)),label=driver_labels[rev(drivers)]),aes(x=-.96,y=y,label=label),inherit.aes=FALSE,hjust=0,size=8.1/ggplot2::.pt)+
 scale_x_continuous(limits=c(-1,1),breaks=c(0,.5,1),expand=expansion(mult=c(.05,.05)))+
 scale_y_continuous(breaks=seq_along(rev(drivers)),labels=NULL,limits=c(.3,10.7),expand=expansion(mult=0))+
 labs(x='Selected / eligible\nLOPO folds',y=NULL,tag='c')+base+theme(axis.ticks.y=element_blank())+guides(shape='none')
p <- (pa|pb|pc)+plot_layout(guides='collect') & theme(legend.position='bottom')
caption <- paste0('Figure S13. Conditional prediction and selection stability of alternative RONA-driver models. ',
 '(a) Leave-one-population-out (LOPO) prediction, with correlation screening, predictor standardization and all-subset AICc selection repeated using each training fold. ',
 '(b) Prediction after deleting each of three coordinate-defined geographic blocks (3, 12 and 4 populations), with the same training-only procedures. ',
 'Open circles, open triangles and filled circles identify selection among linear subsets, GAM subsets and the combined family set, respectively. Climate is optional in every search. ',
 'In (a–b), each point pools held-out predictions for all 19 populations and gives root mean squared error (RMSE) divided by the corresponding training-mean baseline RMSE. The horizontal scale is logarithmic; the dashed line at 1 marks baseline performance. Small vertical offsets separate selection families. ',
 '(c) Number of LOPO training folds selecting each predictor divided by the number in which that predictor survived correlation screening, pooled across all 20 response cases. The denominator is at most 380 folds per predictor and selection family; exact denominators accompany the source data. A predictor excluded in every fold is not assigned zero selection frequency. ',
 'These descriptive frequencies reflect nominal inclusion in selected formulas under competition among candidate subsets, not causal importance. A selected formula may contain smooths penalized nearly to zero; frequencies do not measure effect magnitude. The 20 cases share populations and predictors, and training folds overlap; no independence-based confidence interval is shown. ',
 'Within each training set, redundant predictors (absolute Spearman correlation > 0.8) are screened using the prespecified original-factor and biological priority. Eligible subsets include an intercept-only model and all combinations surviving screening and estimability constraints. ',
 'GAM smooths use k = 3, maximum likelihood and selection penalization; subset ranking uses training AICc. DEL/SYN and LoF/SYN are reference-relative functional-burden indices, not measured fitness load. ',
 'Rows in (a–b) identify SSP126 or SSP585 and mid-century (M; 2041–2070) or late-century (L; 2071–2100). BIO4: temperature seasonality; BIO5: maximum temperature of the warmest month; BIO6: minimum temperature of the coldest month; BIO15: precipitation seasonality; BIO17: precipitation of the driest quarter. ',
 'Validation is nested for downstream screening and model selection only: inherited candidate discovery and upstream RONA estimation remain fixed. Prediction of RONA does not validate population fitness or establish biological causation. Table S1 provides model rankings, uncertainty conditional on selected formulas, diagnostics and selection details.')
stem<-file.path(figure_dir,if(mock)'SYNTHETIC_LAYOUT_ONLY_S13' else 'Figure_S13')
writeLines(caption,paste0(stem,'.caption.txt'))
write.csv(lopo,paste0(stem,'_LOPO.csv'),row.names=FALSE);write.csv(block,paste0(stem,'_geographic_block.csv'),row.names=FALSE);write.csv(pool,paste0(stem,'_selection_frequency.csv'),row.names=FALSE);write.csv(excluded,paste0(stem,'_all_predictor_eligibility.csv'),row.names=FALSE)
require_patchwork_panel_alignment(p,paste0(stem,'.alignment-layout.json'),paste0(stem,'.alignment.json'),width_in=190/25.4,height_in=230/25.4,tolerance_pt=1.5,gutter_tolerance_pt=1.5,audit_script=file.path(skill,'scripts/audit_panel_alignment.py'),python=python,strict=TRUE)
# Native pdf keeps editable standard-font text and exposes per-label glyph geometry;
# cairo_pdf on this HPC build merges distant labels into false overlapping trace boxes.
grDevices::pdf(paste0(stem,'.pdf'),width=190/25.4,height=230/25.4,family='sans',useDingbats=FALSE);print(p);dev.off()
svglite::svglite(paste0(stem,'.svg'),width=190/25.4,height=230/25.4);print(p);dev.off()
ggsave(paste0(stem,'.png'),p,width=190,height=230,units='mm',dpi=600,type='cairo',bg='white')
saveRDS(list(plot=p,LOPO=lopo,geographic_block=block,selection_frequency=pool,synthetic=mock),paste0(stem,'.rds'))
stopifnot(system2(python,c(shQuote(file.path(skill,'scripts/audit_pdf_text.py')),shQuote(paste0(stem,'.pdf')),'--min-pt','8','--json'),stdout=paste0(stem,'.text-audit.json'))==0)
stopifnot(system2(python,c(shQuote(file.path(skill,'scripts/audit_figure_collisions.py')),shQuote(paste0(stem,'.pdf')),'--json-out',shQuote(paste0(stem,'.collision-audit.json'))))==0)
writeLines(capture.output(sessionInfo()),file.path(figure_dir,'sessionInfo.txt'))
if(!mock){
 used<-c('validation.json','nested_selection_validation_performance.csv','nested_selection_predictor_stability.csv')
 hashes<-vapply(file.path(in_dir,used),function(fn)strsplit(system2('sha256sum',shQuote(fn),stdout=TRUE),' ',fixed=TRUE)[[1]][1],character(1))
 jsonlite::write_json(list(validated=TRUE,no_model_fitting=TRUE,figure_count=1L,number='S13',LOPO_rows=nrow(lopo),geographic_rows=nrow(block),selection_rows=nrow(pool),source_SHA256=as.list(setNames(hashes,used)),manual_visual_QA='Pending parent review'),file.path(figure_dir,'validation.json'),auto_unbox=TRUE,pretty=TRUE)
}
cat(if(mock)'SYNTHETIC_LAYOUT_ONLY_COMPLETE' else 'S13_RENDER_COMPLETE','\n')
