#!/usr/bin/env Rscript
# Original-design S3/S4 rendering from validated saved results. No inference rerun.
source('scripts/revision/common.R')
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(jsonlite)})
directory<-file.path(output_dir,'figure_revision_original_design/supplement_S1_S4');dir.create(directory,recursive=TRUE,showWarnings=FALSE)
input<-file.path(output_dir,'niche_reconstruction')
stopifnot(fromJSON(file.path(input,'validation.json'))$validated)
skill<-Sys.getenv("AMMO_FIGURE_QA_DIR", unset="")
source(file.path(skill,'scripts/panel_alignment.R'));source('scripts/revision/measure_S1_S4_panels.R')
python<-Sys.getenv("AMMO_PYTHON", unset="python3")
cols<-c(core='#F8766D',peripheral='#00BFC4')
theme_set(theme_classic(base_size=9,base_family='sans')+theme(axis.text=element_text(size=8,color='black'),axis.title=element_text(size=9),plot.tag=element_text(size=12,face='plain'),legend.text=element_text(size=8),legend.title=element_blank(),legend.position='bottom',plot.background=element_rect(fill='white',colour=NA),legend.background=element_blank(),legend.box.background=element_blank(),plot.margin=margin(6,6,6,6)))
obj<-readRDS(file.path(input,'niche_plot_objects.rds'))
core<-read.csv(file.path(input,'core_density.csv'));peri<-read.csv(file.path(input,'peripheral_density.csv'))
stopifnot(identical(core[,1:2],peri[,1:2]),nrow(core)==10000)
score<-read.csv(file.path(input,'population_pca_scores.csv'))
variance<-read.csv(file.path(input,'pca_variance.csv'))
summary<-read.csv(file.path(input,'test_summary.csv'))
clim<-read.csv(file.path(input,'population_climate.csv'))
stopifnot(nrow(score)==19,nrow(clim)==19,nrow(summary)==6,all(summary$n==999))
# Display overlap uses all saved corrected-density pixels. No kernel/test refitting.
overlap<-core[,1:2];overlap$core_density<-core$corrected_density;overlap$peripheral_density<-peri$corrected_density
overlap$class<-ifelse(overlap$core_density>0&overlap$peripheral_density>0,'Overlap',ifelse(overlap$core_density>0,'Core','Peripheral'))
overlap$alpha<-sqrt(pmax(overlap$core_density,overlap$peripheral_density));overlap<-overlap[overlap$alpha>0,]
limits_x<-range(c(overlap$PC1,score$PC1));limits_y<-range(c(overlap$PC2,score$PC2))
limits_x<-limits_x+c(-1,1)*diff(limits_x)*.08;limits_y<-limits_y+c(-1,1)*diff(limits_y)*.08;limits_y[1]<-limits_y[1]-.15*diff(limits_y)
p1<-ggplot(overlap,aes(PC1,PC2))+geom_raster(aes(fill=class,alpha=alpha),interpolate=FALSE,show.legend=TRUE)+scale_alpha_identity()+scale_fill_manual(values=c(Core=unname(cols['core']),Peripheral=unname(cols['peripheral']),Overlap='#7B7B7B'),breaks=c('Core','Peripheral','Overlap'),limits=c('Core','Peripheral','Overlap'),drop=FALSE)+guides(fill=guide_legend(override.aes=list(alpha=1)))+coord_cartesian(xlim=limits_x,ylim=limits_y,expand=FALSE)+labs(x=sprintf('PC1 (%.1f%%)',variance$proportion[1]*100),y=sprintf('PC2 (%.1f%%)',variance$proportion[2]*100),tag='a')
p1<-p1+annotate('text',x=limits_x[2]-.03*diff(limits_x),y=limits_y[1]+.025*diff(limits_y),label=sprintf('D = %.3f\nI = %.3f',obj$observed$D,obj$observed$I),hjust=1,vjust=0,size=8/ggplot2::.pt)
# The two-axis correlation circle is rotation * saved component SD, not a new PCA.
load<-data.frame(variable=toupper(rownames(obj$pca$rotation)),PC1=obj$pca$rotation[,1]*obj$pca$sdev[1],PC2=obj$pca$rotation[,2]*obj$pca$sdev[2])
stopifnot(nrow(load)==5,all(load$PC1^2+load$PC2^2<=1+1e-8))
circle<-data.frame(x=cos(seq(0,2*pi,length.out=361)),y=sin(seq(0,2*pi,length.out=361)));circle$xend<-c(circle$x[-1],circle$x[1]);circle$yend<-c(circle$y[-1],circle$y[1]);load$label_x<-load$PC1/sqrt(load$PC1^2+load$PC2^2)*1.22;load$label_y<-load$PC2/sqrt(load$PC1^2+load$PC2^2)*1.22
p2<-ggplot()+geom_segment(data=circle,aes(x,y,xend=xend,yend=yend),linewidth=.35)+annotate('segment',x=-1,xend=1,y=0,yend=0,linewidth=.2,color='grey60')+annotate('segment',x=0,xend=0,y=-1,yend=1,linewidth=.2,color='grey60')+geom_segment(data=load,aes(x=0,y=0,xend=PC1*.88,yend=PC2*.88),arrow=grid::arrow(length=unit(1.5,'mm')),linewidth=.35)+geom_text(data=load,aes(x=label_x,y=label_y,label=variable),size=8/ggplot2::.pt)+coord_fixed(xlim=c(-1.5,1.5),ylim=c(-1.5,1.5),expand=FALSE)+scale_x_continuous(breaks=c(-1,0,1))+scale_y_continuous(breaks=c(-1,0,1))+labs(x='PC1 correlation',y='PC2 correlation',tag='b')
# Descriptive 1D density curves recreate c–g geometry only. Each population
# contributes one coordinate value; native climate units in saved source data.
# Bandwidth=nrd0,512points,cut3 match base R density defaults, shared grid perBIO.
# These visual summaries are NOT the prevalence-corrected bivariate niche tests.
density_rows<-list();density_panels<-list()
units<-c(bio4='BIO4 (temperature SD, C)',bio5='BIO5 (C)',bio6='BIO6 (C)',bio15='BIO15 (CV)',bio17='BIO17 (mm)')
for(i in seq_along(bio_names)){
 bio<-bio_names[i];v<-clim[[bio]];if(bio=='bio4')v<-v/100
 bandwidth<-vapply(c('core','peripheral'),function(g)stats::bw.nrd0(v[clim$group==g]),numeric(1));bounds<-range(v)+c(-1,1)*3*max(bandwidth)
 d<-do.call(rbind,lapply(c('core','peripheral'),function(g){z<-stats::density(v[clim$group==g],bw=bandwidth[g],n=512,from=bounds[1],to=bounds[2]);data.frame(bio=bio,group=g,x=z$x,density=z$y,bandwidth=unname(bandwidth[g]))}))
 density_rows[[i]]<-d
 density_panels[[i]]<-ggplot(d,aes(x,density,color=group))+geom_line(linewidth=.55)+scale_color_manual(values=cols,labels=c(core='Core',peripheral='Peripheral'))+labs(x=units[[bio]],y='Population density',tag=letters[i+2])+theme(legend.position='none')
}
# Original topology: two upper panels, three middle panels, two lower panels.
# All climate-density panels share equal grid spans and physical sizes.
fig3<-wrap_plots(list(wrap_plots(p1,p2,ncol=2),wrap_plots(density_panels,ncol=3)),ncol=1,heights=c(1.15,2))
write_csv(overlap,file.path(directory,'S3_overlap_display.csv'));write_csv(load,file.path(directory,'S3_correlation_circle.csv'));write_csv(do.call(rbind,density_rows),file.path(directory,'S3_descriptive_density_curves.csv'));write_csv(clim,file.path(directory,'S3_population_climate_source.csv'))
# Original observed-vs-null histogram grammar, corrected six-test scope.
# D upper row, I lower row; tests across columns. No lettered-panel titles.
null_panels<-list();null_rows<-list();k<-0
names_test<-c(equivalency='Equivalency',similarity_core_fixed='Core fixed',similarity_peripheral_fixed='Peripheral fixed')
for(metric in c('D','I'))for(test in names(names_test)){
 k<-k+1;null<-read.csv(file.path(input,paste0(test,'_null.csv')));vals<-null[[metric]];stopifnot(all(is.finite(vals)),all(vals>=-1e-12),all(vals<=1+1e-12));vals[pmin(abs(vals),1)<1e-12]<-0
 s<-summary[summary$test==test&summary$metric==metric,];stopifnot(nrow(s)==1)
 d<-data.frame(test=test,metric=metric,overlap=vals);null_rows[[k]]<-d
 null_panels[[k]]<-ggplot(d,aes(overlap))+geom_histogram(binwidth=.025,boundary=0,fill='grey65',colour='white',linewidth=.2)+geom_vline(xintercept=s$observed,color='red',linewidth=.65)+scale_x_continuous(limits=c(0,1),breaks=c(0,.5,1),expand=expansion(mult=c(.01,.01)))+scale_y_continuous(breaks=scales::breaks_pretty(n=3),expand=expansion(mult=c(0,.08)))+labs(x=paste(names_test[[test]],if(metric=='D')"Schoener's D" else "Warren's I",sep='\n'),y='Randomizations',tag=letters[k])
}
fig4<-wrap_plots(null_panels,ncol=3)
write_csv(do.call(rbind,null_rows),file.path(directory,'S4_null_values.csv'));write_csv(summary,file.path(directory,'S4_test_summary.csv'))
export<-function(p,stem,w,h,ids=NULL,rows=NULL,columns=NULL){stem<-file.path(directory,stem)
 require_patchwork_panel_alignment(p,paste0(stem,'.alignment-layout.json'),paste0(stem,'.alignment.json'),width_in=w/25.4,height_in=h/25.4,panel_ids=ids,row_groups=rows,column_groups=columns,tolerance_pt=1.5,gutter_tolerance_pt=1.5,audit_script=file.path(skill,'scripts/audit_panel_alignment.py'),python=python,strict=TRUE)
 ggsave(paste0(stem,'.pdf'),p,width=w,height=h,units='mm',device=grDevices::pdf,useDingbats=FALSE,bg='white');ggsave(paste0(stem,'.png'),p,width=w,height=h,units='mm',dpi=300,type='cairo',bg='white');saveRDS(p,paste0(stem,'.rds'))
 stopifnot(system2(python,c(shQuote(file.path(skill,'scripts/audit_pdf_text.py')),shQuote(paste0(stem,'.pdf')),'--min-pt','8','--json'),stdout=paste0(stem,'.text-audit.json'))==0)
 status<-system2(python,c(shQuote(file.path(skill,'scripts/audit_figure_collisions.py')),shQuote(paste0(stem,'.pdf')),'--json-out',shQuote(paste0(stem,'.collision-audit.json'))));stopifnot(status==0)
}
export(fig3,'S3_original_design',190,230,letters[1:7],list(c('a','b'),c('c','d','e'),c('f','g')),list(c('c','f'),c('d','g')))
export(fig4,'S4_original_design',190,145)
writeLines(c('Figure S3. Reconstructed population climate niches, preserving the original seven-panel composition. (a) Saved environmental-prevalence-corrected density in the shared two-axis environmental PCA. Salmon: core; turquoise: peripheral; gray: cells occupied by both estimated niches; opacity increases with square-root relative density. The display window retains all nonzero density cells and all19 population PCA coordinates. Saved overlap: Schoener D=0.12527 and Warren I=0.35181. (b) Variable correlations with the same saved PCA axes (loadings multiplied by component SD), unit circle and arrows. (c–g) Descriptive one-dimensional Gaussian kernel densities of the19 equally weighted population climate values, nine core and ten peripheral. R density defaults (nrd0 bandwidth per class,512points, shared plotting bounds spanning both groups plus three maximum bandwidths) are used only for these display curves; they are not the prevalence-corrected two-axis test estimator. BIO4:temperature seasonality (monthly temperature SD, divided by100 for C); BIO5:maximum temperature of warmest month(C); BIO6:minimum temperature of coldest month(C); BIO15:precipitation seasonality(CV points); BIO17:precipitation of driest quarter(mm). New reconstruction, not recovery of historical exact niche objects.',
'Figure S4. Observed versus randomized niche overlap, preserving gray histograms and red observed lines. Columns: equivalency, core-fixed similarity, peripheral-fixed similarity. Top row: Schoener D; bottom row: Warren I. Gray bar heights count the999saved randomizations per test in overlap bins of width0.025. Red vertical lines show the observed overlap. Equivalency randomizes pooled population labels while preserving group sizes and uses a lower-tail probability. Directional similarity holds the named niche fixed and shifts the other within the common environmental background, using upper-tail probabilities. Six-panel arrangement is the necessary exception to the original four-panel arrangement because the two directional tests are distinct. Exact raw plus-one probabilities and six-test BH-adjusted values are preserved in S4_test_summary.csv and the corresponding supplementary table. No tests rerun.'),file.path(directory,'captions_draft.txt'))
write_json(list(validated=TRUE,scope='Render saved S3/S4 results with original design; no PCA, niche test or biological inference refitted',reviewer_basis=c('R1.63','R1.66','R3.m16','R3.m20','R3.m21','R3.m22'),original_R_template='ammo_analysis_package_20261014/documents/revision_materials/ammo.R lines2248–2319:ecospat density comparisons; original panel architecture from immutable EMF references',S3_panels=7,S4_panels=6,manual_QA='PENDING'),file.path(directory,'validation.json'),pretty=TRUE,auto_unbox=TRUE)
save_session(directory)
