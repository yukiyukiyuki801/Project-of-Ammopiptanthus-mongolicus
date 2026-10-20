#!/usr/bin/env Rscript
# Restore author-provided Figures 3 and 5; no fitting or analogue search.
source('scripts/revision/common.R')
suppressPackageStartupMessages({library(gdm);library(terra);library(sf);library(ggplot2);library(patchwork)})
directory <- file.path(root,'revision_oct18/figure_updates/main5')
dir.create(file.path(directory,'scope'),recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(directory,'contracts'),recursive=TRUE,showWarnings=FALSE)
base <- file.path(bundle,'ammo/gf')
cases <- data.frame(ssp=c('ssp126','ssp126','ssp585','ssp585'),short=c('2040','2070','2040','2070'),period=rep(c('2041-2070','2071-2100'),2))
args<-commandArgs(TRUE)
# T044 supplies validated current-model grid projections. Never substitute the
# incompatible legacy map archive. This script renders saved values only.
projection_dir<-file.path(output_dir,'figure_revision_original_design/required_grid_projection')
projection_validation<-jsonlite::read_json(file.path(projection_dir,'validation.json'))
stopifnot(isTRUE(projection_validation$validated))
compatible<-TRUE
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
  # Shared categorical key is not a comparable data panel. Measure every pie.
  keep <- panel_ids != "class_key"
  rows <- rows[keep,,drop=FALSE]; panel_ids <- panel_ids[keep]
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
theme_set(theme_bw(base_size=9,base_family='sans')+theme(panel.grid=element_blank(),axis.title=element_text(size=9),axis.text=element_text(size=8),plot.title=element_blank(),plot.subtitle=element_text(size=8),plot.tag=element_text(size=13,face='plain'),plot.tag.position=c(0,1),plot.margin=margin(6,6,6,6),legend.title=element_text(size=8),legend.text=element_text(size=8)))
audit_export<-function(plot,stem,width_mm,height_mm,ids,rows,cols,exemptions=list()) {
  prefix<-file.path(directory,stem)
  require_patchwork_panel_alignment(plot,paste0(prefix,'.alignment-layout.json'),paste0(prefix,'.alignment.json'),width_in=width_mm/25.4,height_in=height_mm/25.4,panel_ids=ids,row_groups=rows,column_groups=cols,exemptions=exemptions,audit_script=file.path(skill,'scripts/audit_panel_alignment.py'),python=python,strict=TRUE)
  ggsave(paste0(prefix,'.pdf'),plot,width=width_mm,height=height_mm,units='mm',device=grDevices::pdf,useDingbats=FALSE)
  svglite::svglite(paste0(prefix,'.svg'),width=width_mm/25.4,height=height_mm/25.4);print(plot);dev.off()
  ggsave(paste0(prefix,'.png'),plot,width=width_mm,height=height_mm,units='mm',dpi=300,type='cairo')
  saveRDS(plot,paste0(prefix,'.rds'))
  text_status<-system2(python,c(shQuote(file.path(skill,'scripts/audit_pdf_text.py')),shQuote(paste0(prefix,'.pdf')),'--min-pt','7','--json'),stdout=paste0(prefix,'.text-audit.json'))
  collision_status<-system2(python,c(shQuote(file.path(skill,'scripts/audit_figure_collisions.py')),shQuote(paste0(prefix,'.pdf')),'--json-out',shQuote(paste0(prefix,'.collision-audit.json')),'--overlay-pdf',shQuote(paste0(prefix,'.collision-audit.pdf'))))
  stopifnot(text_status==0,collision_status==0)
}
if(!'--figure5-only'%in%args) {
  grids<-lapply(seq_len(nrow(cases)),function(i){
    ca<-cases[i,];stem<-paste(ca$ssp,ca$period,'ENSEMBLE',sep='_')
    f<-read.csv(file.path(projection_dir,paste0(stem,'_forward_inherited_geographic_penalty_degrees_grid.csv')))
    r<-read.csv(file.path(projection_dir,paste0(stem,'_reverse_inherited_geographic_penalty_degrees_grid.csv')))
    stopifnot(nrow(f)==157745,nrow(r)==157745,identical(f$candidate_index,r$candidate_index),max(abs(f$local_offset-r$local_offset))<1e-12)
    list(forward=f,reverse=r)
  })
  all_channels<-do.call(rbind,lapply(grids,function(g)cbind(local=g$forward$local_offset,forward=g$forward$offset,reverse=g$reverse$offset)))
  channel_limits<-apply(all_channels,2,quantile,probs=c(.2,.7),na.rm=TRUE)
  write_csv(data.frame(metric=colnames(channel_limits),lower=channel_limits[1,],upper=channel_limits[2,]),file.path(directory,'Figure_3_RGB_channel_limits.csv'))
  boundary<-sf::st_read(file.path(base,'chelsa2_crop_selected/current/ammo_extent_country.shp'),quiet=TRUE)
  if(is.na(sf::st_crs(boundary)))sf::st_crs(boundary)<-4326
  boundary<-suppressWarnings(sf::st_crop(boundary,c(xmin=100,xmax=110,ymin=36,ymax=44)))
  # Spatial positions and all157745cells retained, no interpolation or smoothing.
  rgb_maps<-lapply(seq_len(nrow(cases)),function(i){
    ca<-cases[i,];f<-grids[[i]]$forward;r<-grids[[i]]$reverse
    z<-cbind(local=f$local_offset,forward=f$offset,reverse=r$offset)
    stretch<-sapply(seq_len(3),function(j){lim<-channel_limits[,j];stopifnot(diff(lim)>0);pmin(1,pmax(0,(z[,j]-lim[1])/diff(lim)))})
    dat<-data.frame(x=f$longitude,y=f$latitude,colour=grDevices::rgb(stretch[,1],stretch[,2],stretch[,3]))
    # Integer index annotation_raster avoids false irregular-grid assumptions.
    step<-1/120;x0<-min(dat$x);y0<-min(dat$y);ix<-round((dat$x-x0)/step)+1L;iy<-round((dat$y-y0)/step)+1L
    stopifnot(max(abs(dat$x-(x0+(ix-1)*step)))<1e-7,max(abs(dat$y-(y0+(iy-1)*step)))<1e-7)
    image<-matrix(NA_character_,nrow=max(iy),ncol=max(ix));image[cbind(max(iy)-iy+1L,ix)]<-dat$colour
    ggplot()+annotation_raster(image,xmin=x0-step/2,xmax=x0+(max(ix)-.5)*step,ymin=y0-step/2,ymax=y0+(max(iy)-.5)*step,interpolate=FALSE)+geom_sf(data=boundary,fill=NA,colour='black',linewidth=.35)+coord_sf(xlim=c(100,110),ylim=c(36,44),expand=FALSE,datum=NA)+scale_x_continuous(breaks=c(100,102.5,105,107.5,110))+scale_y_continuous(breaks=c(36,38,40,42,44))+labs(x='Longitude (°E)',y='Latitude (°N)',tag=letters[i])+theme(panel.background=element_rect(fill='#EBEBEB',colour=NA),panel.border=element_rect(fill=NA,colour='black',linewidth=.4))
  })
  vertices<-data.frame(x=c(-.866,0,.866,.866,0,-.866,-.866),y=c(.5,1,.5,-.5,-1,-.5,.5))
  dots<-data.frame(x=vertices$x[1:6],y=vertices$y[1:6],colour=c('#FFFF00','#00FF00','#00FFFF','#0000FF','#FF00FF','#FF0000'))
  key<-ggplot()+geom_path(data=vertices,aes(x,y),colour='#C9C9C9',linewidth=1.1,linetype='dotted')+geom_point(data=dots,aes(x,y,fill=colour),shape=21,size=5,stroke=.4)+geom_point(aes(x=0,y=0),size=5,colour='black')+geom_text(aes(x=-.9,y=-1.28,label='Local'),size=2.8)+geom_text(aes(x=0,y=1.43,label='Forward'),size=2.8)+geom_text(aes(x=.9,y=-1.28,label='Reverse'),size=2.8)+scale_fill_identity()+coord_equal(xlim=c(-1.55,1.55),ylim=c(-1.6,1.7),clip='off')+labs(tag='e')+theme_void(base_size=9)+theme(plot.tag=element_text(size=13),plot.tag.position=c(0,1),plot.margin=margin(6,6,6,6))
  design<-'AB#\nCDE'
  fig3<-wrap_plots(c(rgb_maps,list(key)),design=design,widths=c(1,1,.65))
  audit_export(fig3,'Figure_3_original_design',190,150,letters[1:5],list(c('a','b'),c('c','d')),list(c('a','c'),c('b','d')))
  writeLines('Figure 3. RGB composites of local, forward and reverse offsets across four scenarios: (a) SSP126 2041–2070, (b) SSP126 2071–2100, (c) SSP585 2041–2070, (d) SSP585 2071–2100. (e) Red=local, green=forward, blue=reverse; mixed colors combine these channels. Each metric is linearly stretched between its pooled 20th and 70th percentiles across all four scenarios, with identical metric-specific limits applied to every panel and clipping outside those limits. Colors are a display composite, not a summed offset. Projections use the current saved GDM model, inherited geographic penalty in degrees and unrestricted fixed current-habitat candidate grid plus sampled coordinates. All 157,745 focal grid cells are retained without interpolation. Forward endpoints are future destinations; reverse endpoints are present analogues. Neither estimates migration. Separate metric maps are provided in Figure S9.',file.path(directory,'Figure_3_caption_DRAFT.txt'))
}
# Figure5 preserves the original sixteen-pie arrangement while correcting the
# old endpoint-raster count collapse. Existing numerical results only.
if(!'--figure3-only'%in%args){
percent<-read.csv(file.path(root,'revision_oct18/figure_updates/endpoint_percentages.csv'))
percent<-percent[percent$gcm=='ENSEMBLE'&percent$radius_km=='unrestricted'&percent$search_scoring=='inherited_geographic_penalty_degrees',]
stopifnot(all(percent$denominator%in%c(9,10)))
row_specs<-data.frame(direction=c('forward','forward','reverse','reverse'),group=c('core','peripheral','core','peripheral'),panel=letters[1:4])
pies<-list();retained<-list();k<-0L
palette<-c(core='#F8766D',peripheral='#00BFC4')
for(i in seq_len(nrow(row_specs)))for(j in seq_len(nrow(cases))){
  rr<-row_specs[i,];ca<-cases[j,]
  d<-percent[percent$group==rr$group&percent$direction==rr$direction&percent$ssp==ca$ssp&percent$period==ca$period,]
  stopifnot(nrow(d)>0,abs(sum(d$percentage_of_focal_populations)-100)<1e-8,sum(d$n)==if(rr$group=='core')9 else 10)
  stopifnot(all(d$endpoint_habitat_class%in%names(palette)))
  d$endpoint_habitat_class<-factor(d$endpoint_habitat_class,levels=names(palette));d<-d[order(d$endpoint_habitat_class),]
  d$label<-sprintf('%.1f%%',d$percentage_of_focal_populations)
  k<-k+1L;retained[[k]]<-d
  # Cartesian wedge polygons retain exact proportional pie areas and place
  # narrow-slice labels outside, avoiding text crossing the white wedge edge.
  # Small labels have explicit leaders; no category or percentage is hidden.
  d$fraction<-d$n/sum(d$n);d$end<-cumsum(d$fraction)*2*pi+pi/2;d$start<-c(pi/2,head(d$end,-1));d$middle<-(d$start+d$end)/2
  polygons<-do.call(rbind,lapply(seq_len(nrow(d)),function(z){angle<-seq(d$start[z],d$end[z],length.out=max(20,ceiling(400*d$fraction[z])));data.frame(x=c(0,cos(angle),0),y=c(0,sin(angle),0),category=d$endpoint_habitat_class[z],wedge=z)}))
  if(nrow(d)==1){angle<-seq(0,2*pi,length.out=401);polygons<-data.frame(x=cos(angle),y=sin(angle),category=d$endpoint_habitat_class,wedge=1)}
  labels<-data.frame(x=.54*cos(d$middle),y=.54*sin(d$middle),label=d$label,small=d$fraction<.25)
  leaders<-data.frame(x=cos(d$middle[labels$small]),y=sin(d$middle[labels$small]),xend=1.00*sign(cos(d$middle[labels$small])),yend=1.08*sin(d$middle[labels$small]))
  labels$x[labels$small]<-1.06*sign(cos(d$middle[labels$small]));labels$y[labels$small]<-1.08*sin(d$middle[labels$small]);labels$hjust<-ifelse(labels$small,ifelse(labels$x>0,0,1),.5)
  if(nrow(d)==2 && all(d$fraction>.4 & d$fraction<.6)){labels$x<-.52*sign(cos(d$middle));labels$y<-.22*sign(sin(d$middle))}
  if(nrow(d)==1){labels$x<-0;labels$y<-0}
  pies[[k]]<-ggplot()+geom_polygon(data=polygons,aes(x,y,group=wedge,fill=category),colour='white',linewidth=.4,show.legend=FALSE)+geom_segment(data=leaders,aes(x,y,xend=xend,yend=yend),linewidth=.25,colour='#555555')+geom_text(data=labels,aes(x,y,label=label,hjust=hjust),size=2.65)+scale_fill_manual(values=palette,breaks=names(palette),labels=c(core='Core',peripheral='Peripheral'),limits=names(palette),drop=FALSE,name=NULL)+coord_equal(xlim=c(-1.65,1.65),ylim=c(-1.25,1.25),expand=FALSE,clip='off')+labs(subtitle=if(i==1)paste(if(ca$ssp=='ssp126')'SSP126'else'SSP585',ca$period,sep='\n')else' ',tag=if(j==1)rr$panel else NULL)+theme_void(base_size=8)+theme(plot.subtitle=element_text(size=8,hjust=.5),plot.tag=element_text(size=13),plot.tag.position=c(0,1),legend.position='bottom',legend.text=element_text(size=8),legend.key.size=unit(3,'mm'),legend.spacing.y=unit(0,'mm'),plot.margin=margin(5,2,5,2))
}
fig5<-wrap_plots(pies,ncol=4)+plot_annotation(caption='■ Core    ■ Peripheral    ■ Other',theme=theme(plot.caption=element_blank()))
# Exactly one shared class legend, rendered independently of observed wedges.
legend_data<-data.frame(x=c(.30,.62),y=0,category=names(palette),label=c('Core','Peripheral'))
key5<-ggplot(legend_data,aes(x,y))+geom_point(aes(colour=category),shape=15,size=3,show.legend=FALSE)+geom_text(aes(label=label),nudge_x=.025,hjust=0,size=2.8)+scale_colour_manual(values=palette)+coord_cartesian(xlim=c(0,1),ylim=c(-1,1),expand=FALSE)+theme_void()
fig5<-wrap_plots(c(pies,list(key5)),design='ABCD\nEFGH\nIJKL\nMNOP\nQQQQ',heights=c(1,1,1,1,.22))
ids<-paste0('pie',seq_len(16));row_groups<-lapply(0:3,function(i)ids[i*4+1:4]);col_groups<-lapply(1:4,function(j)ids[seq(j,16,by=4)])
audit_export(fig5,'Figure_5_original_design',190,170,c(ids,'class_key'),row_groups,col_groups)
write_csv(do.call(rbind,retained),file.path(directory,'Figure_5_source_percentages.csv'))
writeLines(c('Figure5. Habitat classes of selected matching endpoints under four future scenarios. Rows(a–b)show forward future-destination endpoints for core and peripheral focal populations,respectively; rows(c–d)show reverse present-day candidate-source endpoints for core and peripheral focal populations,respectively. Columns show SSP1262041–2070,SSP1262071–2100,SSP5852041–2070,andSSP5852071–2100. Colors preserve the original design: core coral,peripheral turquoise,other grey.','Each pie summarizes one selected endpoint per focal population,with n=9core or n=10peripheral populations. Results reuse validated fixed-habitat searches plus exact sampled coordinates,with the inherited geographic penalty in degrees and unrestricted radius. Percentages are population counts,not habitat area,actual migration or source-preference probabilities. Forward endpoints are future destinations,not present sources. These corrected denominators replace historical endpoint-raster fractions that collapsed many-to-one focal matches. Numerical rounding may affect displayed sums.'),file.path(directory,'Figure_5_caption_DRAFT.txt'))
}
save_session(directory)
jsonlite::write_json(list(figure5_rendered=!'--figure3-only'%in%args,figure3_rendered=!'--figure5-only'%in%args,figure3_current_model_compatible=compatible,scope='Saved-result rendering only. Manual visual review remains required before integration.'),file.path(directory,'render_status.json'),pretty=TRUE,auto_unbox=TRUE)

# These markers validate successful rendering and automated gates, not final
# publication adoption. Parent performs final-size visual review separately.
for(number in c(3L,5L)){
  enabled<-if(number==3L) !'--figure5-only'%in%args else !'--figure3-only'%in%args
  if(enabled){
    stem<-file.path(directory,paste0('Figure_',number,'_original_design'))
    expected<-paste0(stem,c('.pdf','.png','.svg','.rds','.alignment.json','.text-audit.json','.collision-audit.json'))
    stopifnot(all(file.exists(expected)),all(file.info(expected)$size>0))
    report<-jsonlite::read_json(paste0(stem,'.collision-audit.json'))
    stopifnot(isTRUE(report$auditable),report$summary$fail==0)
    jsonlite::write_json(list(validated=TRUE,figure=number,no_model_refitting=TRUE,render_only=TRUE,
      width_mm=190,height_mm=if(number==3L)150 else 170,
      output_files=expected,manual_review_pending=TRUE,
      scope=if(number==3L)'T044 saved current-model grid offsets; original RGB design; shared pooled metric-specific display scales' else 'Validated19-focal-population endpoint percentages; original16-pie design',
      collision_failures=report$summary$fail,collision_warnings=report$summary$warn),
      file.path(directory,if(number==3L)'validation_main3.json' else 'validation.json'),pretty=TRUE,auto_unbox=TRUE)
  }
}
