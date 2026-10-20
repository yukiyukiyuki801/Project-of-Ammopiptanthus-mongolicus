# Three reviewer-requested evidence displays; exclusively saved numerical results.
# Native pdf retains independent editable text spans; cairo_pdf batches distant labels here.
# Contract: S13 conditional model validation; S14 five-GCM spread; S15 geographic sensitivity.
# Quantitative grids, final 190-mm width, >=8-pt text, native editable PDF and 600-dpi PNG.
source('scripts/revision/common.R')
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(jsonlite)})
figure_dir <- file.path(output_dir,'figure_revision_original_design','supplement_S13_S15')
dir.create(figure_dir,recursive=TRUE,showWarnings=FALSE)
skill <- Sys.getenv("AMMO_FIGURE_QA_DIR", unset="")
source(file.path(skill,'scripts/panel_alignment.R'))
source('scripts/revision/measure_panels.R')
python <- Sys.getenv("AMMO_PYTHON", unset="python3")
read_saved <- function(path) read.csv(file.path(output_dir,path))
case_levels <- c('ssp126_2041-2070','ssp126_2071-2100','ssp585_2041-2070','ssp585_2071-2100')
case_short <- setNames(c('126/M','126/L','585/M','585/L'),case_levels)
variable_levels <- c('bio4','bio5','bio6','bio15','bio17','local','forward','reverse')
var_labels <- setNames(c('BIO4','BIO5','BIO6','BIO15','BIO17','Local','Forward','Reverse'),variable_levels)
condition_levels <- unlist(lapply(case_levels,function(z) paste(z,variable_levels,sep='|')))
condition_labels <- setNames(unlist(lapply(case_levels,function(z) paste(case_short[z],var_labels))),condition_levels)
condition <- function(d,variable='variable') {
  key <- paste(paste(d$ssp,d$period,sep='_'),d[[variable]],sep='|')
  d$condition <- factor(key,levels=rev(condition_levels[condition_levels %in% key]))
  d
}
base_theme <- theme_classic(base_size=9,base_family='sans')+theme(
 axis.text=element_text(size=8,colour='black'),axis.title=element_text(size=9),
 legend.text=element_text(size=8),legend.title=element_blank(),
 plot.tag=element_text(size=11,face='bold'),plot.tag.position=c(0,1),
 plot.margin=margin(9,9,6,9),legend.position='bottom',legend.key.width=unit(5,'mm'),
 panel.grid=element_blank(),plot.background=element_rect(fill='white',colour=NA),legend.background=element_rect(fill='white',colour=NA))
export <- function(p,id,height,data,caption) {
 stem <- file.path(figure_dir,paste0('Figure_S',id))
 require_patchwork_panel_alignment(p,paste0(stem,'.alignment-layout.json'),paste0(stem,'.alignment.json'),width_in=190/25.4,height_in=height/25.4,tolerance_pt=1.5,gutter_tolerance_pt=1.5,audit_script=file.path(skill,'scripts/audit_panel_alignment.py'),python=python,strict=TRUE)
 grDevices::pdf(paste0(stem,'.pdf'),width=190/25.4,height=height/25.4,useDingbats=FALSE);print(p);dev.off()
 svglite::svglite(paste0(stem,'.svg'),width=190/25.4,height=height/25.4);print(p);dev.off()
 ggsave(paste0(stem,'.png'),p,width=190,height=height,units='mm',dpi=600,type='cairo',bg='white')
 saveRDS(list(plot=p,data=data),paste0(stem,'.rds'))
 for(n in names(data)) write_csv(data[[n]],paste0(stem,'_',n,'.csv'))
 writeLines(caption,paste0(stem,'.caption.txt'))
 stopifnot(system2(python,c(shQuote(file.path(skill,'scripts/audit_pdf_text.py')),shQuote(paste0(stem,'.pdf')),'--min-pt','8','--json'),stdout=paste0(stem,'.text-audit.json'))==0)
 audit_status <- system2(python,c(shQuote(file.path(skill,'scripts/audit_figure_collisions.py')),shQuote(paste0(stem,'.pdf')),'--json-out',shQuote(paste0(stem,'.collision-audit.json'))))
 stopifnot(audit_status==0)
 invisible(stem)
}
# Historical S13 branch omitted: final S13 is the all-driver nested validation figure.
# S14: every individual GCM retained; no new tests or population bootstrap intervals.
r <- read_saved('statistics/population_group_effects.csv');r <- r[r$analysis=='RONA',]
g <- read_saved('landscape_statistics/population_class_effects.csv');g <- g[g$search_scoring=='inherited_geographic_penalty_degrees',];g$variable <- g$metric
a14 <- condition(r);b14 <- condition(g)
stopifnot(nrow(a14)==120,nrow(b14)==72)
effect_plot <- function(d,tag) {
 individual <- d[d$gcm!='ENSEMBLE',];ensemble <- d[d$gcm=='ENSEMBLE',]
 gcm_order <- sort(unique(individual$gcm));individual$offset <- match(individual$gcm,gcm_order)*.10-.30
 ggplot()+geom_vline(xintercept=0,linewidth=.35,linetype=2,colour='grey55')+
 geom_point(data=individual,aes(difference_peripheral_minus_core,as.numeric(condition)+offset),colour='grey55',size=1.65)+
 geom_point(data=ensemble,aes(difference_peripheral_minus_core,as.numeric(condition)),shape=18,size=2.5)+
 scale_y_continuous(breaks=seq_along(levels(d$condition)),labels=condition_labels[levels(d$condition)],expand=expansion(add=.7))+
 labs(x='Peripheral - core mean difference',y=NULL,tag=tag)+base_theme
}
p14 <- effect_plot(a14,'a')|effect_plot(b14,'b')
caption14 <- 'Figure S14. Variation among five climate-model projections in population-class comparisons. Peripheral-minus-core mean differences are shown for (a) RONA and (b) GDM local, forward and reverse offsets under the primary inherited geographic scoring. Each grey point is one of five GCMs (GFDL-ESM4, IPSL-CM6A-LR, MPI-ESM1-2-HR, MRI-ESM2-0 and UKESM1-0-LL); black diamonds are the corresponding ensemble estimates. Small vertical offsets separate the five GCM estimates within each row. Positive values indicate a higher mean for peripheral populations. Each class comparison uses nine core and ten peripheral populations. The dashed line marks zero difference. Grey-point spread represents differences among alternative climate projections, not a population bootstrap confidence interval or five independent biological replicates; population-level intervals and tests are reported separately in Tables S3 and S4. SSP126/SSP585 and M/L denote 2041–2070/2071–2100. BIO4: temperature seasonality; BIO5: maximum temperature of the warmest month; BIO6: minimum temperature of the coldest month; BIO15: precipitation seasonality; BIO17: precipitation of the driest quarter. GDM comparisons retain the inherited geographic penalty in longitude/latitude degrees. Suitable-area variation among GCMs is reported in Table S9.'
# S15: original class vs continuous coordinate measure; all 32 ensemble associations.
classes <- read_saved('classification/population_class_and_geographic_marginality.csv')
sensitivity <- read_saved('classification/continuous_geographic_marginality_sensitivity.csv')
sensitivity <- condition(sensitivity[sensitivity$gcm=='ENSEMBLE',])
stopifnot(nrow(classes)==19,nrow(sensitivity)==32,all(sensitivity$n==19))
classes$population <- factor(classes$population,levels=classes$population[order(classes$geographic_marginality_km)])
classes$group <- factor(classes$group,levels=c('core','peripheral'))
p15a <- ggplot(classes,aes(geographic_marginality_km,population,colour=group))+
 geom_vline(xintercept=0,linewidth=.3,linetype=2,colour='grey55')+
 geom_point(size=2.2)+scale_colour_manual(values=c(core='#F8766D',peripheral='#00BFC4'),labels=c('Core','Peripheral'))+
 labs(x='Geographic marginality (km)',y=NULL,tag='a')+base_theme
p15b <- ggplot(sensitivity,aes(rho,condition))+
 geom_vline(xintercept=0,linewidth=.3,linetype=2,colour='grey55')+
 geom_segment(aes(x=CI_low,xend=CI_high,yend=condition),linewidth=.4,colour='grey45')+
 geom_point(size=1.75)+scale_x_continuous(limits=c(-1,1),breaks=c(-1,0,1))+
 scale_y_discrete(labels=condition_labels)+labs(x='Spearman correlation',y=NULL,tag='b')+base_theme
p15 <- (p15a|p15b)+plot_layout(guides='collect') & theme(legend.position='bottom')
caption15 <- 'Figure S15. Sensitivity to continuous geographic marginality. (a) Geographic marginality of the 19 populations, coloured by their published core/peripheral class. Marginality is the negative signed interior distance (km) from a convex hull built from 37 published occurrence coordinates; more negative values are farther inside the sampled envelope and values near zero lie near its boundary. (b) Saved Spearman correlations between geographic marginality and RONA (BIO4, BIO5, BIO6, BIO15 and BIO17) or GDM offsets (local, forward and reverse) for the four ensemble scenario-period combinations. Points are correlations across 19 populations; horizontal intervals are the saved population-bootstrap 95% intervals. The dashed line marks zero correlation. SSP126/SSP585 and M/L denote 2041–2070/2071–2100. BIO4: temperature seasonality; BIO5: maximum temperature of the warmest month; BIO6: minimum temperature of the coldest month; BIO15: precipitation seasonality; BIO17: precipitation of the driest quarter. These are exploratory sensitivities: the occurrence hull is a sampled geographic envelope rather than a certified species boundary, and uneven sampling, shared occurrence provenance and spatial dependence remain. GDM associations here use the saved finite population-candidate domain with environment-only scoring from the geographic GDM fit, and should not be conflated with the full-landscape matching estimates in the main GDM figures. Full estimates, intervals and adjusted probabilities are given in Table S10.'
# Caption/source metadata are written before export so document assembly can proceed.
for(i in c(14,15)) writeLines(get(paste0('caption',i)),file.path(figure_dir,paste0('Figure_S',i,'.caption.txt')))
export(p14,14,200,list(RONA=a14,GDM=b14),caption14)
export(p15,15,235,list(populations=classes,associations=sensitivity),caption15)
save_session(figure_dir)
marker(figure_dir,list(figures=length(14:15),scope='Reviewer-requested GCM spread and continuous marginality; saved results only',no_model_fitting=TRUE,no_new_tests=TRUE,manual_visual_QA='Pending parent review'))
