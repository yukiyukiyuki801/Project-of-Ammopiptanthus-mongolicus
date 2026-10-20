#!/usr/bin/env Rscript
# Data-only rendering: no statistical refit or coordinate transformation.
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1L || length(args) > 2L) stop('Usage: Rscript render_S20.R workflow-config.json [--validate-only]')
validate_only <- length(args)==2L && identical(args[2], '--validate-only')
if (length(args)==2L && !validate_only) stop('Unknown option')
suppressPackageStartupMessages({library(jsonlite);library(ggplot2);library(patchwork)})
cfg <- fromJSON(args[1]); pca <- file.path(cfg$output_root,'pca')
required <- c(cfg$existing_pca_scores,cfg$existing_pca_variance,file.path(pca,'mongolicus_PCA_scores.csv'),file.path(pca,'mongolicus_PCA_variance.csv'),file.path(pca,'validation.json'))
if (!all(file.exists(required))) stop(paste('Actual saved PCA inputs are missing:',paste(required[!file.exists(required)],collapse=', ')))
check <- fromJSON(required[5]);stopifnot(isTRUE(check$validated), check$matched_individuals==217)
old <- read.csv(required[1],stringsAsFactors=FALSE);new <- read.csv(required[3],stringsAsFactors=FALSE)
oldvar <- read.csv(required[2]);newvar <- read.csv(required[4])
check_scores <- function(d) {
 stopifnot(all(c('sample_id','population','group','PC1','PC2') %in% names(d)),nrow(d)==217L,!anyDuplicated(d$sample_id),!anyNA(d$sample_id),all(is.finite(d$PC1)),all(is.finite(d$PC2)),length(unique(d$population))==19L,all(d$group %in% c('core','peripheral')),sum(d$group=='core')==104L,sum(d$group=='peripheral')==113L)
}
check_variance <- function(d) {
 stopifnot(all(c('component','variance_fraction') %in% names(d)),nrow(d)>=2L,identical(as.integer(d$component[1:2]),1:2),all(is.finite(d$variance_fraction)),all(d$variance_fraction>=0),all(d$variance_fraction<=1),abs(sum(d$variance_fraction)-1)<1e-6)
}
check_scores(old);check_scores(new);check_variance(oldvar);check_variance(newvar)
stopifnot(setequal(old$sample_id,new$sample_id));new <- new[match(old$sample_id,new$sample_id),]
stopifnot(identical(old$sample_id,new$sample_id),identical(old$population,new$population),identical(old$group,new$group))
new$PC1 <- -new$PC1; new$PC2 <- -new$PC2
if (validate_only) {cat('Real-input guards PASS; no figure rendered\n');quit(status=0)}
script_arg <- grep('^--file=',commandArgs(),value=TRUE)
base <- dirname(normalizePath(sub('^--file=','',script_arg[1])))
dest <- file.path(base,'S20');dir.create(dest,showWarnings=FALSE)
stem <- file.path(dest,'Figure_S20_reference_PCA')
if (file.exists(paste0(stem,'.pdf'))) stop('Existing rendered PDF requires explicit parent-approved new output or removal after job completion')
skill <- Sys.getenv("AMMO_FIGURE_QA_DIR", unset="")
python <- Sys.getenv("AMMO_PYTHON", unset="python3")
source(file.path(base,'panel_alignment_viewports.R'))
make_panel <- function(d,v,tag) {
 ggplot(d,aes(PC1,PC2,colour=group,shape=group)) + geom_point(size=1.7,alpha=.8) +
 scale_colour_manual(values=c(core='#F8766D',peripheral='#00BFC4'),breaks=c('core','peripheral'),labels=c('Core (104)','Peripheral (113)')) +
 scale_shape_manual(values=c(core=16,peripheral=17),breaks=c('core','peripheral'),labels=c('Core (104)','Peripheral (113)')) +
 labs(x=sprintf('PC1 (%.2f%%)',100*v[1]),y=sprintf('PC2 (%.2f%%)',100*v[2]),tag=tag,colour=NULL,shape=NULL) +
 theme_classic(base_size=11,base_family='Helvetica') +
 theme(axis.text=element_text(size=10,colour='black'),axis.title=element_text(size=11),legend.text=element_text(size=10),legend.position='bottom',plot.tag=element_text(size=13,face='bold'),plot.margin=margin(7,7,7,7),plot.background=element_rect(fill='white',colour=NA))
}
figure <- (make_panel(old,oldvar$variance_fraction,'a') | make_panel(new,newvar$variance_fraction,'b')) + plot_layout(guides='collect',widths=c(1,1)) & theme(legend.position='bottom')
w <- 183/25.4;h <- 94/25.4
require_patchwork_panel_alignment(figure,manifest_path=paste0(stem,'.alignment-layout.json'),report_path=paste0(stem,'.alignment.json'),overlay_svg=paste0(stem,'.alignment.svg'),width_in=w,height_in=h,panel_ids=c('a','b'),row_groups=list(c('a','b')),tolerance_pt=1.5,gutter_tolerance_pt=1.5,audit_script=file.path(skill,'scripts','audit_panel_alignment.py'),python=python,strict=TRUE)
ggsave(paste0(stem,'.pdf'),figure,width=183,height=94,units='mm',device=grDevices::cairo_pdf,family='Helvetica')
ggsave(paste0(stem,'.png'),figure,width=183,height=94,units='mm',dpi=400,type='cairo',bg='white')
saveRDS(figure,paste0(stem,'.rds'))
run_audit <- function(script,arguments,stdout='') {
 status <- system2(python,c(shQuote(if(script=='audit_figure_collisions_chars.py') file.path(base,script) else file.path(skill,'scripts',script)),arguments),stdout=stdout)
 if (!identical(status,0L)) stop(paste('Figure audit blocked delivery:',script,'exit',status))
}
run_audit('audit_pdf_text.py',c(shQuote(paste0(stem,'.pdf')),'--min-pt','8','--json'),stdout=paste0(stem,'.text-audit.json'))
run_audit('audit_figure_collisions_chars.py',c('--strict',shQuote(paste0(stem,'.pdf')),'--json-out',shQuote(paste0(stem,'.collision-audit.json')),'--overlay-pdf',shQuote(paste0(stem,'.collision-audit.pdf'))))
write.csv(old,file.path(dest,'panel_a_source_scores.csv'),row.names=FALSE)
write.csv(new,file.path(dest,'panel_b_source_scores.csv'),row.names=FALSE)
write.csv(oldvar,file.path(dest,'panel_a_source_variance.csv'),row.names=FALSE)
write.csv(newvar,file.path(dest,'panel_b_source_variance.csv'),row.names=FALSE)
write_json(list(rendered=TRUE,input_individuals_per_panel=217,excluded_individuals=0,pca_refitted=FALSE,coordinate_transforms='Panel b: PC1 and PC2 multiplied by minus one for visual orientation; no refit',width_mm=183,height_mm=94,manual_visual_QA='PENDING',interpretation='PENDING',inputs=data.frame(path=required,md5=unname(tools::md5sum(required)))),file.path(dest,'render_validation.json'),pretty=TRUE,auto_unbox=TRUE)
cat('Rendered PDF and PNG; manual panel QA and interpretation remain PENDING\n')
