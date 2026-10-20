#!/usr/bin/env Rscript
# Saved-result rendering only, restoring original Figure2 and Figure6 designs.
source('scripts/revision/common.R')
suppressPackageStartupMessages({library(terra);library(sf);library(ggplot2);library(patchwork);library(jsonlite)})
directory<-file.path(output_dir,'figure_revision_original_design/main_2_6');dir.create(directory,recursive=TRUE,showWarnings=FALSE)
skill<-Sys.getenv("AMMO_FIGURE_QA_DIR", unset="")
source(file.path(skill,'scripts/panel_alignment.R'));source('scripts/revision/measure_panels.R')
python<-Sys.getenv("AMMO_PYTHON", unset="python3")
stopifnot(fromJSON(file.path(output_dir,'rona/validation.json'))$validated,fromJSON(file.path(output_dir,'sdm_reconstruction/validation.json'))$validated)
sf_use_s2(FALSE)
read_shape<-function(name){x<-st_read(file.path(bundle,'ammo/gf',name),quiet=TRUE);if(is.na(st_crs(x)))st_crs(x)<-4326;st_transform(st_make_valid(st_transform(x,32648)),4326)}
habitat<-read_shape('ammo_binary.shp')
core<-read_shape('ammo_core.shp');peri<-read_shape('ammo_peri.shp')
country<-st_read(file.path(bundle,'ammo/gf/chelsa2_crop_selected/current/ammo_extent_country.shp'),quiet=TRUE)
if(is.na(st_crs(country)))st_crs(country)<-4326
country<-st_transform(st_make_valid(country),4326)
# Crop rendering geometries before plotting, so clipped off-map path segments
# cannot extend into axis labels or neighboring panels in the exported PDF.
map_bbox<-st_bbox(c(xmin=100,ymin=36,xmax=110,ymax=44),crs=st_crs(4326))
country<-suppressWarnings(st_crop(country,map_bbox))
habitat_draw<-suppressWarnings(st_crop(habitat,map_bbox))
core<-suppressWarnings(st_crop(core,map_bbox));peri<-suppressWarnings(st_crop(peri,map_bbox))
cases<-data.frame(ssp=c('ssp126','ssp126','ssp585','ssp585'),period=c('2041-2070','2071-2100','2041-2070','2071-2100'),short=c('2040','2070','2040','2070'))
rona<-read.csv(file.path(output_dir,'rona/rona_population_long.csv'));rona<-rona[rona$gcm=='ENSEMBLE',]
stopifnot(nrow(rona)==380,!anyDuplicated(rona[,c('population','bio','ssp','period')]))
avg<-aggregate(rona~population+bio,data=rona,FUN=mean)
avg<-merge(avg,pop_meta[,c('population','longitude','latitude','group')],by='population',sort=FALSE)
stopifnot(nrow(avg)==95,all(is.finite(avg$rona)),all(avg$rona>=0))
write_csv(avg,file.path(directory,'Figure2_RONA_four_case_mean.csv'))
# Confirm original RONA pooling, without reusing old numerical values.
old_directory<-file.path(bundle,'ammo/rona');oldavg<-read.csv(file.path(old_directory,'avg_avg_rona_xy.csv'),row.names=1)
oldcases<-lapply(seq_len(nrow(cases)),function(i)read.csv(file.path(old_directory,paste0(cases$ssp[i],'_',gsub('-','_',cases$period[i]),'_rona_xy.csv')),row.names=1))
old_error<-max(vapply(bio_numbers,function(b){v<-paste0('Bio',b);max(abs(oldavg[,v]-Reduce('+',lapply(oldcases,function(d)d[rownames(oldavg),v]))/4))},numeric(1)))
stopifnot(old_error<1e-12)
units<-c(bio4='degrees C SD',bio5='degrees C',bio6='degrees C',bio15='CV points',bio17='mm')
# Labels are legend identifiers/units, not panel titles.
legend_names<-c(bio4='BIO4 change (C SD)',bio5='BIO5 change (C)',bio6='BIO6 change (C)',bio15='BIO15 change (CV)',bio17='BIO17 change (mm)')
theme_set(theme_bw(base_size=9,base_family='sans')+theme(panel.grid=element_blank(),axis.text=element_text(size=8,color='black'),axis.text.x=element_text(margin=margin(t=5)),axis.text.y=element_text(margin=margin(r=5)),axis.title=element_text(size=8),plot.tag=element_text(size=12,face='plain'),legend.text=element_text(size=8),legend.title=element_text(size=8),plot.margin=margin(5,5,5,5)))
# Native pixels are color-mapped, not spatially interpolated or thinned.
raster_layer<-function(r,palette,lim){v<-values(r,mat=FALSE);e<-as.vector(ext(r));idx<-pmin(length(palette),pmax(1,round((v-lim[1])/diff(lim)*(length(palette)-1))+1));col<-palette[idx];col[!is.finite(v)]<-'white';im<-matrix(col,nrow=nrow(r),ncol=ncol(r),byrow=TRUE);annotation_raster(im,xmin=e[1],xmax=e[2],ymin=e[3],ymax=e[4],interpolate=FALSE)}
plots<-list();audit<-list()
for(i in seq_along(bio_names)){
 bio<-bio_names[i];b<-bio_numbers[i];base<-file.path(bundle,'ammo/gf/chelsa2_crop_selected');cur<-rast(file.path(base,'current',paste0('b',b,'.asc')))
 fut<-lapply(seq_len(4),function(j)rast(file.path(base,paste0(cases$ssp[j],'_',cases$short[j]),paste0('b',b,'.asc'))))
 stopifnot(all(vapply(fut,function(x)compareGeom(cur,x,crs=FALSE),logical(1))))
 delta<-Reduce('+',lapply(fut,function(x)x-cur))/4
 original<-rast(file.path(old_directory,paste0('Bio',b,'_avg_avg_diff.asc')))
 stopifnot(compareGeom(original,delta,crs=FALSE))
 av<-values(original,mat=FALSE);dv<-values(delta,mat=FALSE);ok<-is.finite(av)&is.finite(dv)
 err<-max(abs(av[ok]-dv[ok]));stopifnot(any(ok),err<1e-4)
 # BIO4 conventional SD*100 is converted to degreesC SD solely for display.
 if(bio=='bio4')delta<-delta/100
 masked<-mask(delta,vect(habitat));v<-values(masked,mat=FALSE);lim<-range(v,na.rm=TRUE);stopifnot(diff(lim)>0)
 palette<-grDevices::colorRampPalette(c('lightblue','yellow','tomato'))(256)
 pt<-avg[avg$bio==bio,];stopifnot(nrow(pt)==19,sum(pt$group=='core')==9,sum(pt$group=='peripheral')==10)
 legend_data<-data.frame(x=c(100,100),y=c(36,36),value=lim)
 p<-ggplot()+raster_layer(masked,palette,lim)+geom_sf(data=habitat_draw,fill=NA,colour='black',linewidth=.15)+geom_sf(data=country,fill=NA,colour='grey70',linewidth=.2)+geom_point(data=legend_data,aes(x,y,fill=value),shape=22,alpha=0,show.legend=TRUE)+geom_point(data=pt,aes(longitude,latitude,size=rona),shape=21,fill=NA,colour='red',stroke=.5)+scale_size_continuous(name='RONA',range=c(1,4),limits=c(0,.6),breaks=c(.1,.3,.5),guide='none')+scale_fill_gradientn(colours=palette,limits=lim,breaks=pretty(lim,n=2),name=legend_names[[bio]])+guides(fill=guide_colourbar(title.position='top',barwidth=unit(36,'mm'),barheight=unit(2,'mm')))+coord_sf(default_crs=st_crs(4326),xlim=c(100,110),ylim=c(36,44),expand=FALSE)+scale_x_continuous(breaks=c(101,105,109))+scale_y_continuous(breaks=c(36,40,44))+labs(x='Longitude',y='Latitude',tag=letters[i])+theme(legend.position='bottom',legend.margin=margin(0,0,0,0))
 plots[[i]]<-p
 audit[[i]]<-data.frame(bio=bio,original_average_raster_max_abs_error=err,comparison_cells=sum(ok),display_divisor=if(bio=='bio4')100 else 1,masked_cells=sum(is.finite(v)),markers=19,core=9,peripheral=10,RONA_min=min(pt$rona),RONA_max=max(pt$rona),climate_min=lim[1],climate_max=lim[2])
}
# Original vacant sixth grid slot hosts one common circle-size legend.
key<-data.frame(x=0,y=seq(0,1,length.out=3),rona=c(.1,.3,.5))
plots[[6]]<-ggplot(key,aes(x,y,size=rona))+geom_point(shape=21,colour='red',fill=NA,stroke=.5,alpha=0,show.legend=TRUE)+scale_size_continuous(name='RONA',range=c(1,4),limits=c(0,.6),breaks=c(.1,.3,.5))+scale_x_continuous(limits=c(-1,1))+scale_y_continuous(limits=c(-1,2))+theme_void()+theme(legend.position='right',legend.title=element_text(size=9),legend.text=element_text(size=8))+guides(size=guide_legend(override.aes=list(colour='red',fill=NA,alpha=1)))
# Hide duplicated marks in the spare slot; the scale key remains.
plots[[6]]<-plots[[6]]+geom_blank()+theme(plot.margin=margin(5,5,5,5))
fig2<-wrap_plots(plots,ncol=3)
write_csv(do.call(rbind,audit),file.path(directory,'Figure2_source_audit.csv'))
# Figure6 exact original TIFF legend colors sampled without editing the source.
pal<-read.csv(file.path(directory,'original_Figure6_legend_palette.csv'));stopifnot(nrow(pal)==256)
fig6plots<-lapply(seq_len(4),function(i){case<-cases[i,];r<-rast(file.path(output_dir,'sdm_reconstruction',paste0(case$ssp,'_',case$period,'_ENSEMBLE_mean.tif')));v<-values(r,mat=FALSE);stopifnot(min(v,na.rm=TRUE)>=0,max(v,na.rm=TRUE)<=1)
 d<-data.frame(x=c(100,100),y=c(36,36),value=c(0,1))
 ggplot()+raster_layer(r,pal$colour,c(0,1))+geom_sf(data=peri,fill='grey',colour=NA,alpha=.45)+geom_sf(data=core,fill='green',colour=NA,alpha=.30)+geom_point(data=d,aes(x,y,fill=value),shape=22,alpha=0,show.legend=TRUE)+scale_fill_gradientn(colours=pal$colour,limits=c(0,1),breaks=c(0,.25,.5,.75,1),name='Suitability')+coord_sf(default_crs=st_crs(4326),xlim=c(100,110),ylim=c(36,44),expand=FALSE)+labs(title=paste(toupper(case$ssp),gsub('-','–',case$period)),x=NULL,y=NULL)+theme_void(base_size=9)+theme(plot.title=element_text(size=10,hjust=.5),legend.position='right',legend.text=element_text(size=8),legend.title=element_text(size=9),plot.margin=margin(5,5,5,5))+guides(fill=guide_colourbar(barwidth=unit(3,'mm'),barheight=unit(42,'mm')))
})
fig6<-wrap_plots(fig6plots,ncol=2,guides='collect')&theme(legend.position='right')
export<-function(p,stem,w,h,alignment=TRUE){stem<-file.path(directory,stem)
 if(alignment)require_patchwork_panel_alignment(p,paste0(stem,'.alignment-layout.json'),paste0(stem,'.alignment.json'),width_in=w/25.4,height_in=h/25.4,tolerance_pt=1.5,gutter_tolerance_pt=1.5,audit_script=file.path(skill,'scripts/audit_panel_alignment.py'),python=python,strict=TRUE)
 ggsave(paste0(stem,'.pdf'),p,width=w,height=h,units='mm',device=grDevices::pdf,useDingbats=FALSE,bg='white')
 ggsave(paste0(stem,'.png'),p,width=w,height=h,units='mm',dpi=300,type='cairo',bg='white')
 saveRDS(p,paste0(stem,'.rds'))
 stopifnot(system2(python,c(shQuote(file.path(skill,'scripts/audit_pdf_text.py')),shQuote(paste0(stem,'.pdf')),'--min-pt','8','--json'),stdout=paste0(stem,'.text-audit.json'))==0)
 stopifnot(system2(python,c(shQuote(file.path(skill,'scripts/audit_figure_collisions.py')),shQuote(paste0(stem,'.pdf')),'--json-out',shQuote(paste0(stem,'.collision-audit.json'))))==0)
}
# The spare sixth slot is only a common legend, not a comparable map.
# Measure the actual final layout; explicit groups audit the five map panels.
require_patchwork_panel_alignment(fig2,file.path(directory,'Figure2.alignment-layout.json'),file.path(directory,'Figure2.alignment.json'),width_in=190/25.4,height_in=165/25.4,panel_ids=c('a','b','c','d','e','key'),row_groups=list(c('a','b','c'),c('d','e')),column_groups=list(c('a','d'),c('b','e')),exemptions=list(list(panels='key',checks='all',reason='Vacant original sixth grid slot contains only the common RONA legend; it is not a geographic data panel.')),tolerance_pt=1.5,gutter_tolerance_pt=1.5,audit_script=file.path(skill,'scripts/audit_panel_alignment.py'),python=python,strict=TRUE)
export(fig2,'Figure2_original_design',190,165,FALSE)
export(fig6,'Figure6_original_design',190,155)
writeLines(c('Figure2. Corrected population RONA and projected climatic change. Panels(a–e) showBIO4, temperature seasonality (monthly temperature SD); BIO5, maximum temperature of the warmest month; BIO6, minimum temperature of the coldest month; BIO15, precipitation seasonality (coefficient of variation); andBIO17, precipitation of the driest quarter. Hollow red circle size shows population RONA. Raster color shows signed future-minus-current climatic change. Both are arithmetic means acrossSSP126/SSP585 and2041–2070/2071–2100; this descriptive aggregation does not substitute for scenario-specific tests. BIO4change is divided by100 for degreesC SD; BIO5/BIO6 are degreesC, BIO15CV points, BIO17mm. All19populations retained at actual coordinates; no jitter. Climatic pixels masked to inherited habitat; no landscape RONA interpolation. Circle sizes use a nonzero display floor and should be interpreted against the key. Corrected numbers replace historical values; no analyses refitted for this figure.',
'Figure6. Future potential suitability from the completed reconstructed MaxEnt model. SSP126top andSSP585bottom;2041–2070left and2071–2100right. Colors preserve the original figure legend from low white/lightblue to high green/yellow/orange. Predictions are arithmetic means of five GCM-specific suitability projections. Semitransparent green/gray overlays are inherited current core/peripheral habitat polygons and do not redefine the published population classes. Projections use the same reconstructed model and static elevation. They describe potential suitability within the regional crop, not observed distributions, migration or extinction. No model fitting occurred for this rendering.'),file.path(directory,'captions_draft.txt'))
save_session(directory)
marker(directory,list(scope='Original Figure2/Figure6 design restored from corrected saved results; no model fitting',original_RONA_four_case_mean_max_error=old_error,Figure2_maps=5,Figure2_population_markers=95,Figure6_future_maps=4,manual_visual_QA='PENDING',PDFs=c('Figure2_original_design.pdf','Figure6_original_design.pdf')))
