#!/usr/bin/env Rscript
# Draw saved GAM predictions only; no model fitting or subset selection.
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(jsonlite)})
a<-commandArgs(TRUE);stopifnot(length(a)==2);src<-a[1];out<-a[2];dir.create(out,recursive=TRUE,showWarnings=FALSE)
skill<-Sys.getenv("AMMO_FIGURE_QA_DIR", unset="");python<-Sys.getenv("AMMO_PYTHON", unset="python3");source(file.path(skill,'scripts/panel_alignment.R'));source('scripts/revision/measure_panels.R')
v<-read_json(file.path(src,'validation.json'),simplifyVector=TRUE);stopifnot(isTRUE(v$validated),v$n_cases==20)
d<-read.csv(file.path(src,'top1_GAM_conditional_curve_source.csv'));d$bio<-sub('.*_','',d$case_id);d$scenario<-sub('_bio[0-9]+$','',d$case_id)
bios<-c('bio4','bio5','bio6','bio15','bio17');drivers<-c('climate_change','Ho','He','DEL','LoF','private_alleles','Fis','Alt','Lat','Lon')
labels<-c('Absolute climate\nchange','Observed\nheterozygosity','Expected\nheterozygosity','DEL / SYN\n','LoF / SYN\n','Private\nallele count','FIS\n','Elevation\n(m)','Latitude\n(degrees N)','Longitude\n(degrees E)');names(labels)<-drivers
scenarios<-c('ssp126_2041-2070','ssp126_2071-2100','ssp585_2041-2070','ssp585_2071-2100');scenario_labels<-c('SSP126 mid-century','SSP126 late-century','SSP585 mid-century','SSP585 late-century');colors<-setNames(c('#0072B2','#56B4E9','#D55E00','#CC79A7'),scenarios);lines<-setNames(c('solid','dashed','solid','dashed'),scenarios);d$scenario<-factor(d$scenario,levels=scenarios)
stopifnot(all(is.finite(as.matrix(d[c('x_original','conditional_mean','lower95','upper95')]))),all(d$lower95<=d$conditional_mean),all(d$upper95>=d$conditional_mean))
base<-theme_classic(base_size=8.5,base_family='sans')+theme(axis.text=element_text(size=8,colour='black'),axis.title=element_text(size=8.5),plot.title=element_text(size=8.5,hjust=.5,face='plain'),legend.text=element_text(size=8.5),legend.title=element_blank(),legend.key.width=unit(7,'mm'),plot.margin=margin(4,5,4,5),panel.grid=element_blank(),plot.background=element_rect(fill='white',colour=NA),legend.background=element_rect(fill='white',colour=NA),legend.box.background=element_blank())
panels<-list();empty<-character();summary<-list();k<-0
for(i in seq_along(bios)){
 z<-d[d$bio==bios[i],];yr<-range(c(z$lower95,z$upper95));ylim<-yr+c(-1,1)*max(diff(yr),.01)*.04
 for(j in seq_along(drivers)){
  k<-k+1;zz<-z[z$predictor==drivers[j],];id<-paste(bios[i],drivers[j],sep='_');if(!nrow(zz))empty<-c(empty,id)
  # Empty cells preserve the grid; they do not denote zero effect.
  if(!nrow(zz))zz<-data.frame(x_original=c(0,1),conditional_mean=NA_real_,lower95=NA_real_,upper95=NA_real_,scenario=factor(scenarios[1],levels=scenarios))
  q<-ggplot(zz,aes(x_original,conditional_mean,color=scenario,fill=scenario,linetype=scenario,group=scenario))+geom_ribbon(aes(ymin=lower95,ymax=upper95),alpha=.10,linewidth=0,color=NA,na.rm=TRUE)+geom_line(linewidth=.4,na.rm=TRUE)+scale_color_manual(values=colors,breaks=scenarios,labels=scenario_labels,drop=FALSE)+scale_fill_manual(values=colors,breaks=scenarios,labels=scenario_labels,drop=FALSE)+scale_linetype_manual(values=lines,breaks=scenarios,labels=scenario_labels,drop=FALSE)+scale_x_continuous(n.breaks=if(drivers[j]=='Fis')2 else 3,labels=function(x)format(signif(x,3),trim=TRUE,scientific=FALSE))+scale_y_continuous(limits=ylim,n.breaks=3)+labs(x=NULL,y=if(j==1)paste0(toupper(bios[i]),'\nConditional RONA')else NULL,title=if(i==1)labels[drivers[j]]else NULL)+base+theme(legend.position='bottom')+guides(color=guide_legend(nrow=1),fill='none',linetype=guide_legend(nrow=1))
  if(k>1)q<-q+guides(color='none',fill='none',linetype='none')
  if(!nrow(z[z$predictor==drivers[j],]))q<-q+theme(axis.line=element_line(colour='transparent'),axis.ticks=element_line(colour='transparent'),axis.text=element_text(colour='transparent'))
  panels[[k]]<-q
 }
 summary[[bios[i]]]<-list(n_panels=length(unique(z$predictor)),predictors=drivers[drivers%in%z$predictor])
}
p<-wrap_plots(panels,ncol=10,guides='collect')&theme(legend.position='bottom');stem<-file.path(out,'Figure_S6');w<-420;h<-260
write_patchwork_panel_layout(p,paste0(stem,'.alignment-layout.json'),w/25.4,h/25.4,panel_ids=unlist(lapply(bios,function(b)paste(b,drivers,sep='_'))))
stopifnot(system2(python,c(shQuote(file.path(skill,'scripts/audit_panel_alignment.py')),shQuote(paste0(stem,'.alignment-layout.json')),'--json-out',shQuote(paste0(stem,'.alignment.json'))))==0)
# Native PDF retains editable glyphs; cairo_pdf was avoided because its PDF geometry causes false whole-path collision alarms.
grDevices::pdf(paste0(stem,'.pdf'),width=w/25.4,height=h/25.4,family='sans',useDingbats=FALSE);print(p);dev.off()
ggsave(paste0(stem,'.png'),p,width=w,height=h,units='mm',dpi=400,type='cairo',bg='white');svglite::svglite(paste0(stem,'.svg'),width=w/25.4,height=h/25.4);print(p);dev.off();saveRDS(d,paste0(stem,'_data.rds'));write.csv(d,paste0(stem,'_source.csv'),row.names=FALSE)
cap<-paste('Figure S6. Conditional RONA predictions from selected GAMs. Rows correspond to BIO4 (temperature seasonality), BIO5 (maximum temperature of the warmest month), BIO6 (minimum temperature of the coldest month), BIO15 (precipitation seasonality), and BIO17 (precipitation of the driest quarter). Columns identify predictors; x axes use original predictor units. Curves hold all other predictors at their population means. Colours identify SSP126 mid-century (blue), SSP126 late-century (light blue), SSP585 mid-century (orange-red), and SSP585 late-century (purple); solid and dashed lines distinguish periods. Ribbons show pointwise 95% intervals for the conditional mean, including available smoothing-parameter uncertainty but excluding subset-selection uncertainty. Curves appear only for predictors in each scenario-specific selected formula; empty cells and absent curves do not denote zero effects. Nominally retained smooths can be penalized almost to zero. All panels in a BIO row share the response scale, including scenario baseline differences, so curves are not centered effect sizes. Models use 19 populations and are conditional associations rather than causal or independently validated fitness effects. DEL/SYN and LoF/SYN are reference-relative functional-burden indices. Table S1 and Figure S13 provide model details and nested predictive validation. The figure is supplied as one landscape sheet to retain readable axes for every selected predictor.')
writeLines(cap,paste0(stem,'.caption.txt'))
stopifnot(system2(python,c(shQuote(file.path(skill,'scripts/audit_pdf_text.py')),shQuote(paste0(stem,'.pdf')),'--min-pt','8','--json'),stdout=paste0(stem,'.text-audit.json'))==0)
stopifnot(system2(python,c(shQuote(file.path(skill,'scripts/audit_figure_collisions.py')),shQuote(paste0(stem,'.pdf')),'--json-out',shQuote(paste0(stem,'.collision-audit.json'))))==0)
write_json(list(validated=TRUE,no_model_fitting=TRUE,figure_number='S6',n_plates=1,width_mm=w,height_mm=h,source_rows=nrow(d),empty_grid_cells=empty,panels=summary,manual_visual_QA='Pending parent review'),file.path(out,'validation.json'),auto_unbox=TRUE,pretty=TRUE)
