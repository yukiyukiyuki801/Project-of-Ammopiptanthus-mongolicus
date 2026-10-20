# Original two-panel functional comparison design; saved data only.
source('scripts/revision/common.R')
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(jsonlite)})
directory <- file.path(output_dir,'figure_revision_original_design','supplement_S5')
dir.create(directory,recursive=TRUE,showWarnings=FALSE)
skill <- Sys.getenv("AMMO_FIGURE_QA_DIR", unset="")
source(file.path(skill,'scripts/panel_alignment.R'))
source('scripts/revision/measure_panels.R')
python <- Sys.getenv("AMMO_PYTHON", unset="python3")
d <- read.csv(file.path(output_dir,'functional_reconstruction/postprocess/population_ratios.csv'))
tests <- read.csv(file.path(output_dir,'functional_reconstruction/postprocess/population_class_comparisons.csv'))
classes <- c('SIFT_deleterious','strict_LOF')
stopifnot(all(classes %in% d$numerator_class))
d <- d[d$metric=='summed_ALT_frequency_ratio_common_mask' & d$numerator_class %in% classes,]
stopifnot(nrow(d)==38,all(table(d$numerator_class,d$group)==rep(c(9,10),each=2)))
write_csv(d,file.path(directory,'Figure_S5_source.csv'))
plots <- lapply(seq_along(classes),function(i) {
 x <- d[d$numerator_class==classes[i],];x$group <- factor(x$group,levels=c('core','peripheral'))
 s <- do.call(rbind,lapply(levels(x$group),function(g) data.frame(group=g,mean=mean(x$value[x$group==g]),SD=sd(x$value[x$group==g]))))
 s$group <- factor(s$group,levels=levels(x$group))
 test <- tests[tests$numerator_class==classes[i]&tests$metric=='summed_ALT_frequency_ratio_common_mask',]
 stopifnot(nrow(test)==1)
 ymax <- max(s$mean+s$SD)*1.22
 ggplot(s,aes(group,mean,fill=group))+geom_col(width=.7)+
   geom_errorbar(aes(ymin=pmax(0,mean-SD),ymax=mean+SD),width=.15,linewidth=.35)+
   scale_fill_manual(values=c(core='#F8766D',peripheral='#00BFC4'))+
   scale_y_continuous(limits=c(0,ymax),expand=expansion(mult=c(0,.02)))+
   annotate('segment',x=1,xend=2,y=ymax*.91,yend=ymax*.91,linewidth=.3)+
   annotate('text',x=1.5,y=ymax*.96,label=sprintf('q = %.3f',test$BH_q),size=8/2.845)+
   labs(x=NULL,y=if(i==1)'Deleterious / synonymous' else 'LoF / synonymous',tag=letters[i])+theme_classic(base_size=9)+
   theme(axis.text=element_text(size=8,colour='black'),axis.title=element_text(size=9),legend.position='none',plot.tag=element_text(size=12),plot.margin=margin(8,10,5,8))
})
p <- wrap_plots(plots,nrow=1)
stem <- file.path(directory,'Figure_S5')
require_patchwork_panel_alignment(p,paste0(stem,'.alignment-layout.json'),paste0(stem,'.alignment.json'),width_in=180/25.4,height_in=85/25.4,tolerance_pt=1.5,gutter_tolerance_pt=1.5,audit_script=file.path(skill,'scripts/audit_panel_alignment.py'),python=python,strict=TRUE)
grDevices::pdf(paste0(stem,'.pdf'),width=180/25.4,height=85/25.4,useDingbats=FALSE);print(p);dev.off()
ggsave(paste0(stem,'.png'),p,width=180,height=85,units='mm',dpi=600,type='cairo',bg='white');saveRDS(p,paste0(stem,'.rds'))
stopifnot(system2(python,c(shQuote(file.path(skill,'scripts/audit_pdf_text.py')),shQuote(paste0(stem,'.pdf')),'--min-pt','8','--json'),stdout=paste0(stem,'.text-audit.json'))==0)
stopifnot(system2(python,c(shQuote(file.path(skill,'scripts/audit_figure_collisions.py')),shQuote(paste0(stem,'.pdf')),'--json-out',shQuote(paste0(stem,'.collision-audit.json'))))==0)
writeLines(c('Figure S5. Reference-relative functional variant ratios. (a) SIFT-predicted deleterious missense/synonymous; (b) strict predicted LoF/synonymous. The ratios divide summed ALT-allele frequencies at each consequence class by the synonymous sum, using the common mask requiring at least80% called individuals at every population. Bars show means and error bars ±1SD across nine core and ten peripheral populations. Colors follow the original figure: salmon core, turquoise peripheral. The saved Welch-test q values were BH adjusted across21 class–metric comparisons. These are RAD-panel reference-relative indices, not independently polarized derived-allele burden or fitness load. Other estimators, historical annotation fractions, scoring coverage and mask sensitivity remain in Tables S28–S39.'),paste0(stem,'.caption.txt'))
marker(directory,list(scope='Original S5 two-bar-panel layout; corrected saved population ratios only',figures=1,PDF=paste0(stem,'.pdf'),PNG=paste0(stem,'.png'),no_new_tests=TRUE,manual_visual_QA='Pending'))
