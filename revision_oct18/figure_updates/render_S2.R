#!/usr/bin/env Rscript
# Original S1/S2 data-driven remakes; archived results only, no model fitting.
source('scripts/revision/common.R')
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(sf);library(terra);library(ggrepel);library(jsonlite)})
directory<-file.path(root,'revision_oct18/figure_updates/S2');dir.create(directory,recursive=TRUE,showWarnings=FALSE)
skill<-Sys.getenv("AMMO_FIGURE_QA_DIR", unset="");source(file.path(skill,'scripts/panel_alignment.R'));source('scripts/revision/measure_S1_S2_panels.R');python<-Sys.getenv("AMMO_PYTHON", unset="python3")
cols<-c(core='#F8766D',peripheral='#00BFC4');group_cols<-c(core='#0072B2',peripheral='#D55E00',adaptive_1='#CC79A7',adaptive_2='#009E73')
theme_set(theme_bw(base_size=9,base_family='sans')+theme(axis.text=element_text(size=8,color='black'),axis.title=element_text(size=9),plot.tag=element_text(size=12),legend.text=element_text(size=8),legend.title=element_blank(),legend.position='bottom',panel.grid.minor=element_blank(),plot.background=element_rect(fill='white',colour=NA),legend.background=element_blank(),legend.box.background=element_blank(),plot.margin=margin(5,5,5,5)))
# S2 three original panels: reconstructed genetic PCA + two preserved history sets.
pca_dir<-file.path(output_dir,'figures/genetic_pca');scores<-read.csv(file.path(pca_dir,'PCA_scores_and_group_assignments.csv'));var<-read.csv(file.path(pca_dir,'PCA_variance.csv'));stopifnot(nrow(scores)==217)
pca<-ggplot(scores,aes(PC1,PC2,color=group,shape=archived_group))+geom_point(size=1.2,alpha=.85)+scale_color_manual(values=cols,labels=c(core='Core pop.',peripheral='Peripheral pop.'))+scale_shape_manual(values=c(core=16,peripheral=17,adaptive_1=15,adaptive_2=18),labels=c(core='Core',peripheral='Peripheral',adaptive_1='Adaptive 1',adaptive_2='Adaptive 2'))+labs(x=sprintf('PC1 (%.1f%%)',100*var$variance_fraction[1]),y=sprintf('PC2 (%.1f%%)',100*var$variance_fraction[2]),tag='a')+guides(colour=guide_legend(ncol=1,order=1),shape=guide_legend(ncol=2,order=2,override.aes=list(colour='black')))+theme(legend.box='vertical',legend.spacing.y=unit(1,'mm'),legend.key.width=unit(3,'mm'),legend.margin=margin(0,0,0,0))
manifest<-read.csv(file.path(output_dir,'figures/demography/history_source_manifest.csv'));selected<-manifest[manifest$page==1,];stopifnot(nrow(selected)==6)
curves<-read.csv(gzfile(file.path(output_dir,'figures/demography/lossless_history_plot_data.csv.gz')))
curves<-merge(curves,selected[,c('source','group')],by='source',sort=FALSE);stopifnot(setequal(unique(curves$source),selected$source))
curves$panel<-ifelse(grepl('core_peripheral/neutral/',curves$source),'b','c')
stopifnot(all(curves$year>0),all(curves$Ne>=0),all(curves$high>=curves$low))
# Saved percentile envelopes are retained. All archived curve vertices preserved
# in the lossless display data, with no inference or smoothing.
history<-function(d,colors,labels,tag){ggplot(d,aes(year,Ne,color=group,fill=group))+geom_ribbon(aes(ymin=low,ymax=high),alpha=.18,colour=NA)+geom_line(linewidth=.5)+scale_color_manual(values=colors,labels=labels)+scale_fill_manual(values=colors,labels=labels)+scale_x_log10(breaks=c(.001,1,1e3,1e5),labels=c('0.001','1','1k','100k'))+scale_y_continuous(breaks=scales::breaks_pretty(n=3),labels=scales::label_number(scale_cut=scales::cut_short_scale()))+labs(x='Years before present',y=expression(italic(N)[e]),tag=tag)+guides(color=guide_legend(ncol=2),fill=guide_legend(ncol=2))+theme(axis.title.y=element_text(size=12))}
pool<-history(curves[curves$panel=='b',],c(Core=cols['core'][[1]],all_Peri=cols['peripheral'][[1]]),c(Core='Core',all_Peri='Peripheral'),'b')
group<-history(curves[curves$panel=='c',],group_cols,c(core='Core',peripheral='Periph.',adaptive_1='Adaptive 1',adaptive_2='Adaptive 2'),'c')
fig2<-wrap_plots(pca,pool,group,nrow=1)
write_csv(scores,file.path(directory,'S2_genetic_PCA_points.csv'));write_csv(curves,file.path(directory,'S2_saved_history_plot_data.csv'));write_csv(selected,file.path(directory,'S2_history_source_manifest.csv'))
export<-function(p,stem,w,h,ids=NULL,rows=NULL,columns=NULL){stem<-file.path(directory,stem);require_patchwork_panel_alignment(p,paste0(stem,'.alignment-layout.json'),paste0(stem,'.alignment.json'),width_in=w/25.4,height_in=h/25.4,panel_ids=ids,row_groups=rows,column_groups=columns,tolerance_pt=1.5,gutter_tolerance_pt=1.5,audit_script=file.path(skill,'scripts/audit_panel_alignment.py'),python=python,strict=TRUE)
 ggsave(paste0(stem,'.pdf'),p,width=w,height=h,units='mm',device=grDevices::pdf,useDingbats=FALSE,bg='white');ggsave(paste0(stem,'.png'),p,width=w,height=h,units='mm',dpi=300,type='cairo',bg='white');saveRDS(p,paste0(stem,'.rds'))
 stopifnot(system2(python,c(shQuote(file.path(skill,'scripts/audit_pdf_text.py')),shQuote(paste0(stem,'.pdf')),'--min-pt','8','--json'),stdout=paste0(stem,'.text-audit.json'))==0)
 status<-system2(python,c(shQuote(file.path(skill,'scripts/audit_figure_collisions.py')),shQuote(paste0(stem,'.pdf')),'--json-out',shQuote(paste0(stem,'.collision-audit.json'))));stopifnot(status==0)
}
export(fig2,'S2_original_design',190,100)

jsonlite::write_json(list(rendered=TRUE,model_refitted=FALSE,shape_legend=TRUE,individuals=217,manual_QA="PENDING"),file.path(directory,"validation.json"),pretty=TRUE,auto_unbox=TRUE)
