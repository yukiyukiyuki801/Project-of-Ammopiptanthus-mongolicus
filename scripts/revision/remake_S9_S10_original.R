#!/usr/bin/env Rscript
# Original supplemental designs from validated current-model grid values.
# No model fitting, analogue searches or source-image reuse.
source('scripts/revision/common.R')
# Use the PROJ/GDAL resources of the installed R spatial environment.
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(sf);library(geosphere)})
directory<-file.path(output_dir,'figure_revision_original_design/supplement_S9_S10')
dir.create(directory,recursive=TRUE,showWarnings=FALSE)
input_dir<-file.path(output_dir,'figure_revision_original_design/required_grid_projection')
stopifnot(isTRUE(jsonlite::read_json(file.path(input_dir,'validation.json'))$validated))
cases<-data.frame(ssp=c('ssp126','ssp126','ssp585','ssp585'),period=rep(c('2041-2070','2071-2100'),2))
paths<-unlist(lapply(seq_len(4),function(i)file.path(input_dir,paste0(cases$ssp[i],'_',cases$period[i],'_ENSEMBLE_',c('forward','reverse'),'_inherited_geographic_penalty_degrees_grid.csv'))))
jsonlite::write_json(list(input_files=paths,md5=as.list(tools::md5sum(paths)),no_model_refitting=TRUE,no_search=TRUE),file.path(directory,'inputmanifest.json'),pretty=TRUE,auto_unbox=TRUE)
grids<-lapply(seq_len(4),function(i){f<-read.csv(paths[(i-1)*2+1]);r<-read.csv(paths[(i-1)*2+2]);stopifnot(nrow(f)==157745,nrow(r)==157745,identical(f$candidate_index,r$candidate_index),max(abs(f$local_offset-r$local_offset))<1e-12);list(forward=f,reverse=r)})
boundary<-sf::st_read(file.path(bundle,'ammo/gf/chelsa2_crop_selected/current/ammo_extent_country.shp'),quiet=TRUE)
if(is.na(sf::st_crs(boundary)))sf::st_crs(boundary)<-4326
boundary<-suppressWarnings(sf::st_crop(boundary,c(xmin=100,xmax=110,ymin=36,ymax=44)))
skill<-Sys.getenv("AMMO_FIGURE_QA_DIR", unset="")
source(file.path(skill,'scripts/panel_alignment.R'))
source('scripts/revision/measure_panels.R')

# Measure actual rendered grid viewports. The upstream helper's direct
# convertWidth(null units) reports zero sizes outside the panel viewport.
# Keep the same manifest schema and auditor/tolerance; measure on the same
# native PDF device used for export instead of accepting unauditable geometry.
write_patchwork_panel_layout <- function(plot,manifest_path,width_in,height_in,
  panel_ids=NULL,row_groups=NULL,column_groups=NULL,exemptions=list()) {
  grob <- patchwork::patchworkGrob(plot)
  rows <- grob$layout[grepl("^panel(?:-[0-9]+|;.*)?$",grob$layout$name,perl=TRUE),,drop=FALSE]
  rows <- rows[order(rows$t,rows$l),,drop=FALSE]
  if(nrow(rows)<2) {
    print(grob$layout[,c("name","t","l","b","r")]);stop("Cannot identify two plot-area viewport cells")
  }
  if(is.null(panel_ids)) panel_ids <- letters[seq_len(nrow(rows))]
  stopifnot(length(panel_ids)==nrow(rows))
  # Only comparable geographic data panels are measured; headers, row labels,
  # shared keys and polar summaries serve distinct roles.
  keep <- grepl("^map_",panel_ids)
  rows <- rows[keep,,drop=FALSE];panel_ids <- panel_ids[keep]
  path <- tempfile(fileext=".pdf")
  grDevices::pdf(path,width=width_in,height=height_in,family="sans",useDingbats=FALSE)
  on.exit({grDevices::dev.off();unlink(path)},add=TRUE)
  grid::grid.newpage();grid::grid.draw(grob);grid::grid.force()
  panels <- lapply(seq_len(nrow(rows)),function(i) {
    row <- rows[i,]
    viewport <- paste0(row$name,".",row$t,"-",row$l,"-",row$b,"-",row$r)
    grid::upViewport(0)
    grid::seekViewport(viewport)
    if(grepl(";",row$name,fixed=TRUE)) {
      grid::downViewport("layout")
      grid::downViewport("panel.1-1-1-1")
    }
    bottom <- grid::deviceLoc(grid::unit(0,"npc"),grid::unit(0,"npc"),valueOnly=TRUE)
    top <- grid::deviceLoc(grid::unit(1,"npc"),grid::unit(1,"npc"),valueOnly=TRUE)
    bbox <- c(bottom$x,bottom$y,top$x,top$y)*72
    stopifnot(all(is.finite(bbox)),bbox[3]>bbox[1],bbox[4]>bbox[2])
    list(id=panel_ids[i],bbox_pt=unname(bbox),grid_id="patchwork-grid-1",row_start=as.integer(row$t-1),row_stop=as.integer(row$b),col_start=as.integer(row$l-1),col_stop=as.integer(row$r))
  })
  manifest <- list(schema_version=1L,backend="r-patchwork",figure=list(width_pt=width_in*72,height_pt=height_in*72),
    panels=panels,row_groups=.nature_alignment_groups(row_groups),column_groups=.nature_alignment_groups(column_groups),
    exemptions=exemptions,measurement="Actual panel viewport device coordinates on final-size native PDF; no tolerance exemptions")
  jsonlite::write_json(manifest,manifest_path,pretty=TRUE,auto_unbox=TRUE,null="null",digits=12)
  invisible(manifest)
}

python<-Sys.getenv("AMMO_PYTHON", unset="python3")
theme_set(theme_bw(base_size=8,base_family='sans')+theme(panel.grid=element_blank(),axis.text=element_text(size=8),axis.text.x=element_text(size=8,margin=margin(t=5)),axis.text.y=element_text(size=8,margin=margin(r=5)),axis.title=element_text(size=8),plot.title=element_blank(),plot.subtitle=element_text(size=8),plot.tag=element_text(size=12),plot.tag.position=c(.02,.96),plot.margin=margin(4,4,4,4),legend.position='none'))
# Native PDF core fonts retain editable searchable text; no whole-figure rasterization.
# The numeric raster is drawn by R; boundaries, labels, keys and axes remain vectors.
map_plot<-function(d,value,palette,limits,subtitle=' ',tag=NULL){
  colfun<-scales::col_numeric(palette,domain=limits,na.color=NA_character_)
  cols<-colfun(value);step<-1/120;x0<-min(d$longitude);y0<-min(d$latitude)
  ix<-round((d$longitude-x0)/step)+1L;iy<-round((d$latitude-y0)/step)+1L
  stopifnot(max(abs(d$longitude-(x0+(ix-1)*step)))<1e-7,max(abs(d$latitude-(y0+(iy-1)*step)))<1e-7)
  img<-matrix(NA_character_,nrow=max(iy),ncol=max(ix));img[cbind(max(iy)-iy+1L,ix)]<-cols
  ggplot()+annotation_raster(img,xmin=x0-step/2,xmax=x0+(max(ix)-.5)*step,ymin=y0-step/2,ymax=y0+(max(iy)-.5)*step,interpolate=FALSE)+geom_sf(data=boundary,fill=NA,colour='black',linewidth=.25)+coord_sf(xlim=c(100,110),ylim=c(36,44),expand=FALSE,datum=sf::st_crs(4326))+scale_x_continuous(breaks=c(102,105,108),labels=scales::label_number(accuracy=1))+scale_y_continuous(breaks=c(36,40,44),labels=scales::label_number(accuracy=1))+labs(x='Longitude (°E)',y='Latitude (°N)',subtitle=subtitle,tag=tag)+theme(panel.background=element_rect(fill='white',colour=NA))
}
label_plot<-function(label,angle=0){ggplot()+annotate('text',x=0,y=0,label=label,angle=angle,size=8/ggplot2::.pt)+coord_cartesian(xlim=c(-1,1),ylim=c(-1,1))+theme_void()+theme(plot.margin=margin(0,0,0,0))}
key_plot<-function(limits,palette,label,breaks=NULL){
  d<-data.frame(x=seq(limits[1],limits[2],length.out=256),y=0)
  ggplot(d,aes(x,y,fill=x))+geom_raster()+scale_fill_gradientn(colours=palette,limits=limits)+scale_x_continuous(breaks=if(is.null(breaks))scales::breaks_pretty(n=4)(limits)else breaks,expand=expansion(0))+scale_y_continuous(expand=expansion(0))+coord_cartesian(ylim=c(-.5,.5),clip='off')+labs(x=label,y=NULL)+theme_classic(base_size=8)+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line=element_blank(),axis.title.x=element_text(size=8,margin=margin(t=2)),axis.text.x=element_text(size=8),legend.position='none',plot.margin=margin(2,10,2,10))
}
export_plot<-function(plot,stem,width,height,ids,rows,columns){
  prefix<-file.path(directory,stem)
  require_patchwork_panel_alignment(plot,paste0(prefix,'.alignment-layout.json'),paste0(prefix,'.alignment.json'),width_in=width/25.4,height_in=height/25.4,panel_ids=ids,row_groups=rows,column_groups=columns,audit_script=file.path(skill,'scripts/audit_panel_alignment.py'),python=python,strict=TRUE)
  ggsave(paste0(prefix,'.pdf'),plot,width=width,height=height,units='mm',device=grDevices::pdf,useDingbats=FALSE)
  ggsave(paste0(prefix,'.png'),plot,width=width,height=height,units='mm',dpi=300,type='cairo',bg='white')
  svglite::svglite(paste0(prefix,'.svg'),width=width/25.4,height=height/25.4);print(plot);dev.off()
  saveRDS(plot,paste0(prefix,'.rds'))
  text_status<-system2(python,c(shQuote(file.path(skill,'scripts/audit_pdf_text.py')),shQuote(paste0(prefix,'.pdf')),'--min-pt','7','--json'),stdout=paste0(prefix,'.text-audit.json'))
  collision_status<-system2(python,c(shQuote(file.path(skill,'scripts/audit_figure_collisions.py')),shQuote(paste0(prefix,'.pdf')),'--json-out',shQuote(paste0(prefix,'.collision-audit.json')),'--overlay-pdf',shQuote(paste0(prefix,'.collision-audit.pdf'))))
  stopifnot(text_status==0,collision_status==0)
}
# S9 preserves the original exact3×4metric/scenario map matrix. Shared scale
# provides direct comparisons across all twelve maps; no per-map rescaling.
if(!'--s10-only'%in%commandArgs(TRUE)){
max_offset<-max(unlist(lapply(grids,function(g)c(g$forward$offset,g$reverse$offset,g$forward$local_offset))))
offset_palette<-c('darkmagenta','yellow')
plots9<-lapply(seq_len(4),function(i)label_plot(paste(toupper(cases$ssp[i]),cases$period[i],sep='\n')))
ids9<-paste0('header_',1:4);map_ids<-matrix('',3,4)
for(i in seq_len(3)){
  plots9[[length(plots9)+1]]<-label_plot(c('Forward','Reverse','Local')[i],angle=90);ids9<-c(ids9,paste0('row_',i))
  for(j in seq_len(4)){
    f<-grids[[j]]$forward;r<-grids[[j]]$reverse
    value<-switch(i,f$offset,r$offset,f$local_offset)
    plots9[[length(plots9)+1]]<-map_plot(f,value,offset_palette,c(0,max_offset),subtitle=NULL)
    map_ids[i,j]<-paste0('map_',i,'_',j);ids9<-c(ids9,map_ids[i,j])
  }
}
plots9[[length(plots9)+1]]<-key_plot(c(0,max_offset),offset_palette,'Genomic offset');ids9<-c(ids9,'shared_key')
fig9<-wrap_plots(plots9,design='#ABCD\nEFGHI\nJKLMN\nOPQRS\n#TTTT',widths=c(.16,1,1,1,1),heights=c(.24,1,1,1,.28))
export_plot(fig9,'Figure_S9',190,165,ids9,lapply(1:3,function(i)map_ids[i,]),lapply(1:4,function(j)map_ids[,j]))
writeLines('Figure S9. Forward, reverse and local genomic offsets (rows) under SSP126 and SSP585 during 2041–2070 and2071–2100 (columns). Each panel contains all157,745 current-habitat focal grid cells. Values are projections from the current saved GDM model, using the inherited geographic penalty in degrees and unrestricted fixed current-habitat candidates plus sampled coordinates. Forward comparisons select future destinations; reverse comparisons select present-day analogues; local comparisons retain location. All twelve panels share the same darkmagenta-to-yellow offset scale, with no within-panel rescaling. Numeric surfaces correspond to the RGB channels in Figure3. These are model-based matching offsets, not measured migration or fitness.',file.path(directory,'Figure_S9_caption.txt'))
}
# S10 reproduces the original forward-only bearing/rose/distance layout.
bearing_palette<-c('#440154FF','#21908CFF','#FDE725FF')
distance_max<-max(vapply(grids,function(g)max(g$forward$distance_km),numeric(1)))
# Original gradient2: darkmagenta low, orange at100km, yellow high, globally shared.
distance_palette<-scales::gradient_n_pal(c('darkmagenta','orange','yellow'),values=c(0,100/distance_max,1))(seq(0,1,length.out=256))
bearings<-list();roses<-list();audit<-list()
for(i in seq_len(4)){
  f<-grids[[i]]$forward
  b<-geosphere::bearing(as.matrix(f[,c('longitude','latitude')]),as.matrix(f[,c('endpoint_longitude','endpoint_latitude')]))
  zero<-f$distance_km<1e-9;b[zero]<-NA_real_
  d<-data.frame(candidate_index=f$candidate_index,longitude=f$longitude,latitude=f$latitude,endpoint_longitude=f$endpoint_longitude,endpoint_latitude=f$endpoint_latitude,bearing_degrees=b,distance_km=f$distance_km)
  con<-gzfile(file.path(directory,paste0('S10_',cases$ssp[i],'_',cases$period[i],'_plot_data.csv.gz')),'wt');write.csv(d,con,row.names=FALSE,na='');close(con)
  bins<-cut(b,breaks=seq(-180,180,by=10),include.lowest=TRUE,right=FALSE)
  counts<-as.numeric(table(bins));rose<-data.frame(angle=seq(-175,175,by=10),n=counts)
  stopifnot(sum(rose$n)==sum(!zero),all(is.finite(b[!zero])))
  bearings[[i]]<-b;roses[[i]]<-rose
  audit[[i]]<-data.frame(ssp=cases$ssp[i],period=cases$period[i],focal_cells=nrow(f),undefined_zero_distance=sum(zero),rose_cells=sum(rose$n))
}
rose_max<-max(unlist(lapply(roses,function(d)d$n)))
plots10<-list();ids10<-character();bearing_ids<-distance_ids<-character()
for(i in seq_len(4)){
  f<-grids[[i]]$forward;subtitle<-paste(toupper(cases$ssp[i]),cases$period[i])
  plots10[[length(plots10)+1]]<-map_plot(f,bearings[[i]],bearing_palette,c(-180,180),subtitle=subtitle,tag=if(i==1)'a'else NULL)
  bearing_ids<-c(bearing_ids,paste0('map_bearing_',i));ids10<-c(ids10,tail(bearing_ids,1))
  rose<-ggplot(roses[[i]],aes(angle,n,fill=angle))+geom_col(width=10,colour=NA)+scale_fill_gradientn(colours=bearing_palette,limits=c(-180,180))+scale_x_continuous(limits=c(-180,180),breaks=c(-180,-90,0,90),labels=NULL,expand=expansion(0))+scale_y_continuous(limits=c(0,rose_max*1.5),breaks=rose_max,labels=NULL,expand=expansion(0))+geom_text(data=data.frame(angle=c(-180,-90,0,90),n=rose_max*1.28,label=c('S','W','N','E')),aes(label=label),inherit.aes=TRUE,size=8/ggplot2::.pt)+coord_polar(start=pi,clip='off')+labs(x=NULL,y=NULL)+theme_minimal(base_size=8)+theme(axis.text=element_text(size=8),panel.grid.major=element_line(colour='#BDBDBD',linewidth=.25),panel.grid.major.x=element_blank(),panel.grid.minor=element_blank(),legend.position='none',plot.margin=margin(2,2,2,2))
  plots10[[length(plots10)+1]]<-rose;ids10<-c(ids10,paste0('rose_',i))
  plots10[[length(plots10)+1]]<-map_plot(f,f$distance_km,distance_palette,c(0,distance_max),subtitle=subtitle,tag=if(i==1)'b'else NULL)
  distance_ids<-c(distance_ids,paste0('map_distance_',i));ids10<-c(ids10,tail(distance_ids,1))
}
plots10[[length(plots10)+1]]<-key_plot(c(-180,180),bearing_palette,'Bearing (° clockwise from north)',breaks=c(-180,-90,0,90,180));ids10<-c(ids10,'bearing_key')
plots10[[length(plots10)+1]]<-key_plot(c(0,distance_max),distance_palette,'Matching distance (km)');ids10<-c(ids10,'distance_key')
fig10<-wrap_plots(plots10,design='ABC\nDEF\nGHI\nJKL\nMMN',widths=c(1,.42,1),heights=c(1,1,1,1,.24))
export_plot(fig10,'Figure_S10',190,250,ids10,lapply(1:4,function(i)c(bearing_ids[i],distance_ids[i])),list(bearing_ids,distance_ids))
write_csv(do.call(rbind,audit),file.path(directory,'S10_bearing_data_audit.csv'))
write_csv(do.call(rbind,lapply(1:4,function(i)cbind(cases[i,],roses[[i]]))),file.path(directory,'S10_rose_histogram_counts.csv'))
writeLines(paste0('Figure S10. Forward matching endpoint directions and distances across four scenarios, ordered SSP1262041–2070, SSP1262071–2100, SSP5852041–2070, SSP5852071–2100 (rows). (a)Bearing maps with adjacent rose histograms; (b)great-circle distances from present focal cells to their selected future endpoints. Bearings are calculated directly from the saved endpoint coordinates, in degrees clockwise from north (N=0°, E=90°, W=−90°, S=±180°). Roses use10°bins and identical radial count scales across scenarios; their outer radial ring corresponds to ',rose_max,' focal cells per bin. Each scenario contains157,745 focal current-habitat grid cells. Zero-distance matches, if present, have undefined bearing and are omitted only from bearings/roses; counts are recorded in the data audit. All distances and all focal cells remain in distance maps. Bearing and distance color scales are shared across scenarios. These model-selected future destinations do not estimate actual dispersal, migration or gene flow.'),file.path(directory,'Figure_S10_caption.txt'))
save_session(directory)
marker(directory,list(rendered=c('S9','S10'),rendered_this_job=if('--s10-only'%in%commandArgs(TRUE))'S10' else c('S9','S10'),source='T044 validated current-fit full-grid projections',no_model_refitting=TRUE,no_search=TRUE,manual_review_pending=TRUE,width_mm=190,height_mm=c(165,250),focal_cells_per_scenario=157745))
