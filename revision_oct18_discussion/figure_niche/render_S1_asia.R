#!/usr/bin/env Rscript
# Original S1/S2 data-driven remakes; archived results only, no model fitting.
source('scripts/revision/common.R')
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(sf);library(terra);library(ggrepel);library(jsonlite)})
directory<-file.path(root,'revision_oct18_discussion/figure_niche/S1');dir.create(directory,recursive=TRUE,showWarnings=FALSE)
skill<-Sys.getenv("AMMO_FIGURE_QA_DIR", unset="");source(file.path(skill,'scripts/panel_alignment.R'));source('scripts/revision/measure_S1_S2_panels.R');python<-Sys.getenv("AMMO_PYTHON", unset="python3")
cols<-c(core='#F8766D',peripheral='#00BFC4');group_cols<-c(core='#0072B2',peripheral='#D55E00',adaptive_1='#CC79A7',adaptive_2='#009E73')
theme_set(theme_bw(base_size=9,base_family='sans')+theme(axis.text=element_text(size=8,color='black'),axis.title=element_text(size=9),plot.tag=element_text(size=12),legend.text=element_text(size=8),legend.title=element_blank(),legend.position='bottom',panel.grid.minor=element_blank(),plot.background=element_rect(fill='white',colour=NA),legend.background=element_blank(),legend.box.background=element_blank(),plot.margin=margin(5,5,5,5)))
# S2 three original panels: reconstructed genetic PCA + two preserved history sets.
pca_dir<-file.path(output_dir,'figures/genetic_pca');scores<-read.csv(file.path(pca_dir,'PCA_scores_and_group_assignments.csv'));var<-read.csv(file.path(pca_dir,'PCA_variance.csv'));stopifnot(nrow(scores)==217)
pca<-ggplot(scores,aes(PC1,PC2,color=group,shape=archived_group))+geom_point(size=1.2,alpha=.85)+scale_color_manual(values=cols,labels=c(core='Core population',peripheral='Peripheral population'))+scale_shape_manual(values=c(core=16,peripheral=17,adaptive_1=15,adaptive_2=18),labels=c(core='Core group',peripheral='Peripheral group',adaptive_1='Adaptive 1',adaptive_2='Adaptive 2'))+labs(x=sprintf('PC1 (%.1f%%)',100*var$variance_fraction[1]),y=sprintf('PC2 (%.1f%%)',100*var$variance_fraction[2]),tag='a')+guides(colour=guide_legend(ncol=1,order=1),shape='none')
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
# S1 original two large maps with central regional overview and sampling zoom.
# Every component is drawn from coordinates/geometries/numerical suitability.
sf_use_s2(FALSE)
read_shape<-function(path){x<-st_read(path,quiet=TRUE);if(is.na(st_crs(x)))st_crs(x)<-4326;st_transform(st_make_valid(st_transform(x,32648)),4326)}
country<-read_shape(file.path(bundle,'ammo/gf/chelsa2_crop_selected/current/ammo_extent_country.shp'))
core<-read_shape(file.path(bundle,'ammo/gf/ammo_core.shp'));peri<-read_shape(file.path(bundle,'ammo/gf/ammo_peri.shp'))
context<-sf::st_read(file.path(root,'revision_oct18/figure_updates/map_sources/ne_110m_admin_0_countries.geojson'),quiet=TRUE)
# The original recovered overview geometry covers China. Dissolve county borders
# for the context overview; this is cartographic geometry processing only.
# Keep individual country boundaries and crop actual inset geometry before plotting.
context <- suppressWarnings(sf::st_crop(sf::st_make_valid(context),c(xmin=65,ymin=15,xmax=145,ymax=60)))
map_bbox<-st_bbox(c(xmin=100,ymin=36,xmax=110,ymax=44),crs=st_crs(4326));zoom_bbox<-st_bbox(c(xmin=105.8,ymin=38.5,xmax=107.5,ymax=40.8),crs=st_crs(4326));country<-suppressWarnings(st_crop(country,map_bbox));core<-suppressWarnings(st_crop(core,map_bbox));peri<-suppressWarnings(st_crop(peri,map_bbox));zoom_country<-suppressWarnings(st_crop(country,zoom_bbox))
# Separate native vector segments retain every boundary vertex and avoid PDF
# ring bounding boxes extending over unrelated text at other parts of the map.
edges<-function(x){q<-st_coordinates(st_cast(st_geometry(x),'MULTILINESTRING'));idcols<-setdiff(colnames(q),c('X','Y'));ids<-apply(q[,idcols,drop=FALSE],1,paste,collapse='|');k<-which(ids[-length(ids)]==ids[-1]);data.frame(x=q[k,'X'],y=q[k,'Y'],xend=q[k+1,'X'],yend=q[k+1,'Y'])}
country_edges<-edges(country);zoom_edges<-edges(zoom_country)
# Original in-map legend is an opaque map-key overlay; mask the underlying
# administrative strokes in that key box, avoiding invisible stroke collisions.
legend_mask<-function(d)d[!(pmax(d$x,d$xend)>=100.1&pmin(d$x,d$xend)<=102.8&pmax(d$y,d$yend)>=36.2&pmin(d$y,d$yend)<=39.4),]
country_a_edges<-legend_mask(country_edges)
edge_layer<-function(d,color='grey60',lw=.15)geom_segment(data=d,aes(x,y,xend=xend,yend=yend),inherit.aes=FALSE,colour=color,linewidth=lw)
r<-rast(file.path(output_dir,'sdm_reconstruction/current_suitability.tif'));binary<-rast(file.path(output_dir,'sdm_reconstruction/current_suitable_MTSS.tif'));stopifnot(compareGeom(r,binary))
# Show full unthresholded MaxEnt suitability, as requested; no refit.
palette<-grDevices::colorRampPalette(c('#2C7BB6','#ABDDA4','#FFFFBF','#FDAE61','#D7191C'))(256)
raster_layer<-function(r){v<-values(r,mat=FALSE);e<-as.vector(ext(r));idx<-pmax(1,pmin(256,round(v*255)+1));col<-palette[idx];col[!is.finite(v)]<-'white';annotation_raster(matrix(col,nrow=nrow(r),ncol=ncol(r),byrow=TRUE),xmin=e[1],xmax=e[2],ymin=e[3],ymax=e[4],interpolate=FALSE)}
pop<-pop_meta;stopifnot(nrow(pop)==19,sum(pop$group=='core')==9)
# Readability correction: avoid repeating the11central labels on both full maps.
# Their true points remain on both maps; the original central zoom carries their
# labels. Eight outer populations are labelled on both full maps.
outer_labels<-data.frame(population=c('AE','AF','AYE','AZS','BJ','BW','YL','ZS'),label_x=c(102.8,104.4,103.4,104.7,103.1,107.6,107.1,104.5),label_y=c(42.5,40.1,41.4,40.95,37.35,42.2,37.8,37.0))
zoom_labels<-data.frame(population=c('AZA','AG','AA','AD','AB','AZM','EHA','EEY','WWQ','AC','EEA'),label_x=c(rep(106.1,6),rep(107.25,5)),label_y=c(40.65,40.15,39.9,39.6,39.3,38.75,40.65,40.4,40.12,39.72,38.75))
# Clip cartographic strokes underneath opaque text labels, keeping every
# underlying source geometry intact in the input manifest.
mask_labels<-function(edges,lab,xhalf=.6,yhalf=.24){for(i in seq_len(nrow(lab))){hit<-pmax(edges$x,edges$xend)>=lab$label_x[i]-xhalf&pmin(edges$x,edges$xend)<=lab$label_x[i]+xhalf&pmax(edges$y,edges$yend)>=lab$label_y[i]-yhalf&pmin(edges$y,edges$yend)<=lab$label_y[i]+yhalf;edges<-edges[!hit,]};edges}
country_a_edges<-mask_labels(country_a_edges,outer_labels);country_b_edges<-mask_labels(country_edges,outer_labels);zoom_edges<-mask_labels(zoom_edges,zoom_labels,.25,.10)
zoom_connections<-merge(pop,zoom_labels,by='population');zoom_connections$xend<-ifelse(zoom_connections$label_x<106.5,zoom_connections$label_x+.23,zoom_connections$label_x-.23);zoom_connections$yend<-zoom_connections$label_y
zoom_connections<-do.call(rbind,lapply(seq_len(nrow(zoom_connections)),function(i){d<-zoom_connections[i,];f<-seq(0,1,length.out=51);data.frame(x=head(d$longitude+f*(d$xend-d$longitude),-1),y=head(d$latitude+f*(d$yend-d$latitude),-1),xend=tail(d$longitude+f*(d$xend-d$longitude),-1),yend=tail(d$latitude+f*(d$yend-d$latitude),-1))}));zoom_connections<-mask_labels(zoom_connections,zoom_labels,.25,.10)
labels<-function(d){geom_label(data=d,aes(label_x,label_y,label=population),inherit.aes=FALSE,size=8/ggplot2::.pt,fill='white',colour='black',linewidth=0,label.padding=unit(.15,'mm'))}
map_base<-function(){ggplot()+coord_sf(crs=st_crs(4326),default_crs=st_crs(4326),xlim=c(100,110),ylim=c(36,44),expand=FALSE)+scale_x_continuous(breaks=c(101,105,109))+scale_y_continuous(breaks=c(36,40,44))+labs(x='Longitude',y='Latitude')+theme(panel.grid=element_blank())}
point_layer<-geom_point(data=pop,aes(longitude,latitude,color=group),size=1.8,stroke=.4)
# Original suitability thresholding is represented by the reconstructed MTSS mask;
# the threshold and source correction are explicit in the caption.
mapa<-map_base()+raster_layer(r)+edge_layer(country_a_edges)+point_layer+labels(outer_labels)+scale_colour_manual(values=cols,labels=c(core='Core',peripheral='Peripheral'))+labs(tag='a')
# Add a compact numeric suitability key in the original vacant southwest space.
key<-data.frame(xmin=100.25,xmax=100.55,ymin=seq(36.3,38.3,length.out=5),ymax=seq(36.3,38.3,length.out=5)+.4,colour=palette[round(seq(0,1,length.out=5)*255)+1],label=seq(0,1,length.out=5))
mapa<-mapa+annotate('rect',xmin=100.1,xmax=102.8,ymin=36.2,ymax=39.4,fill='white',colour=NA)+geom_rect(data=key,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),inherit.aes=FALSE,fill=key$colour,colour=NA)+geom_text(data=key,aes(x=100.8,y=(ymin+ymax)/2,label=label),inherit.aes=FALSE,size=8/ggplot2::.pt,hjust=0)+annotate('text',x=100.25,y=39.1,label='Suitability',size=8/ggplot2::.pt,hjust=0)
mapb<-map_base()+geom_sf(data=peri,fill='grey85',colour=NA)+geom_sf(data=core,fill='#B8DAB5',colour=NA)+edge_layer(country_b_edges)+point_layer+labels(outer_labels)+scale_color_manual(values=cols,labels=c(core='Core',peripheral='Peripheral'))+labs(tag='b')
overview<-ggplot()+geom_sf(data=context,fill='grey95',color='grey50',linewidth=.2)+annotate('rect',xmin=100,xmax=110,ymin=36,ymax=44,fill=NA,colour='red',linewidth=.5)+coord_sf(xlim=c(65,145),ylim=c(15,60),expand=FALSE)+theme_void()+theme(plot.background=element_rect(fill='white',colour=NA),plot.margin=margin(3,3,3,3))
zoom_pop<-pop[pop$longitude>=105.8&pop$longitude<=107.5&pop$latitude>=38.5&pop$latitude<=40.8,]
zoom<-ggplot()+geom_sf(data=zoom_country,fill='grey95',colour=NA)+edge_layer(zoom_edges)+geom_point(data=zoom_pop,aes(longitude,latitude,color=group),size=1.7)+edge_layer(zoom_connections,'grey50',.2)+labels(zoom_labels)+scale_color_manual(values=cols,guide='none')+coord_sf(crs=st_crs(4326),default_crs=st_crs(4326),xlim=c(105.8,107.5),ylim=c(38.5,40.8),expand=FALSE)+theme_void()+theme(plot.background=element_rect(fill='white',colour=NA),panel.border=element_rect(fill=NA,colour='grey50',linewidth=.3),plot.margin=margin(3,3,3,3))
fig1<-wrap_plots(list(mapa,wrap_plots(overview,zoom,ncol=1,heights=c(.65,1.35)),mapb),nrow=1,widths=c(1.4,.7,1.4),guides='collect')&theme(legend.position='bottom')
write_csv(outer_labels,file.path(directory,'S1_outer_label_positions.csv'));write_csv(zoom_labels,file.path(directory,'S1_zoom_label_positions.csv'));write_csv(pop,file.path(directory,'S1_population_coordinates_classes.csv'));write_csv(zoom_pop,file.path(directory,'S1_zoom_population_subset.csv'));write_csv(data.frame(layer=c('suitability','threshold_mask','context','regional_boundaries','core_polygon','peripheral_polygon'),source=c('outputs/sdm_reconstruction/current_suitability.tif','outputs/sdm_reconstruction/current_suitable_MTSS.tif','cluster_recovery/personal2/chong/Ammo/ENM/soilgrids/china_mongolia_3.shp','ammo/gf/chelsa2_crop_selected/current/ammo_extent_country.shp','ammo/gf/ammo_core.shp','ammo/gf/ammo_peri.shp')),file.path(directory,'S1_cartographic_sources.csv'))
export<-function(p,stem,w,h,ids=NULL,rows=NULL,columns=NULL){stem<-file.path(directory,stem);require_patchwork_panel_alignment(p,paste0(stem,'.alignment-layout.json'),paste0(stem,'.alignment.json'),width_in=w/25.4,height_in=h/25.4,panel_ids=ids,row_groups=rows,column_groups=columns,tolerance_pt=1.5,gutter_tolerance_pt=1.5,audit_script=file.path(skill,'scripts/audit_panel_alignment.py'),python=python,strict=TRUE)
 ggsave(paste0(stem,'.pdf'),p,width=w,height=h,units='mm',device=grDevices::pdf,useDingbats=FALSE,bg='white');ggsave(paste0(stem,'.png'),p,width=w,height=h,units='mm',dpi=300,type='cairo',bg='white');saveRDS(p,paste0(stem,'.rds'))
 stopifnot(system2(python,c(shQuote(file.path(skill,'scripts/audit_pdf_text.py')),shQuote(paste0(stem,'.pdf')),'--min-pt','8','--json'),stdout=paste0(stem,'.text-audit.json'))==0)
 status<-system2(python,c(shQuote(file.path(skill,'scripts/audit_figure_collisions.py')),shQuote(paste0(stem,'.pdf')),'--json-out',shQuote(paste0(stem,'.collision-audit.json'))));stopifnot(status==0)
}

export(fig1,'S1_original_design',190,125,c('a','overview','zoom','b'),list(c('a','b')),NULL)

writeLines("Figure S1. Sampling locations and habitat context. (a) Full current MaxEnt suitability, without removing predictions below MTSS, and the 19 sampled populations. Salmon and turquoise points denote the published nine core and ten peripheral populations, respectively. Upper inset: country boundaries in eastern and northern Asia, with a red rectangle marking the study region (100–110°E, 36–44°N). Lower inset: the 11 central populations. The remaining eight populations are labeled on both full maps. (b) Core (light green) and peripheral (gray) habitat polygons. Population classification follows Yang et al. (2022), whereas habitat classes describe matching locations; these classifications are distinct. Coordinates are WGS84.",file.path(directory,"caption.txt"))
jsonlite::write_json(list(rendered=TRUE,models_refitted=FALSE,suitability_threshold_masked=FALSE,country_context="Natural Earth 1:110m countries cropped to65–145E,15–60N",population_n=19,manual_QA="PENDING"),file.path(directory,"validation.json"),pretty=TRUE,auto_unbox=TRUE)
