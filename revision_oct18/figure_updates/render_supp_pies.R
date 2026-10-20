#!/usr/bin/env Rscript
# Render authorized distance-limited endpoint compositions from saved population counts.
source("scripts/revision/common.R")
suppressPackageStartupMessages({library(ggplot2);library(patchwork)})
skill <- Sys.getenv("AMMO_FIGURE_QA_DIR", unset="")
source(file.path(skill,"scripts/panel_alignment.R"))
directory <- file.path(root,"revision_oct18/figure_updates/distance_pies")
dir.create(directory,recursive=TRUE,showWarnings=FALSE)
if("--first-case" %in% commandArgs(TRUE))marker <- function(directory,value){value$validated <- TRUE;jsonlite::write_json(value,file.path(directory,"validation_precheck.json"),pretty=TRUE,auto_unbox=TRUE,na="null")}
input_directory <- file.path(output_dir,"oct16_revision","distance_classification")
validation <- jsonlite::read_json(file.path(input_directory,"validation.json"))
stopifnot(isTRUE(validation$validated))
percent <- read.csv(file.path(root,"revision_oct18/figure_updates/endpoint_percentages.csv"))
percent <- percent[percent$search_scoring=="inherited_geographic_penalty_degrees",]
stopifnot(nrow(percent)>0,all(percent$gcm=="ENSEMBLE"),all(percent$denominator%in%c(9,10)))
stopifnot(all(percent$n[percent$endpoint_habitat_class=="overlap"]==0))
write_csv(percent,file.path(directory,"Figure_S16_source_percentages.csv"))
python <- Sys.getenv("AMMO_PYTHON", unset="python3")
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

cases <- data.frame(ssp=c("ssp126","ssp126","ssp585","ssp585"),period=rep(c("2041-2070","2071-2100"),2))
radii <- c("10","50","100","200","unrestricted")
row_specs <- data.frame(direction=c("forward","forward","reverse","reverse"),group=c("core","peripheral","core","peripheral"),panel=letters[1:4])
palette <- c(core="#F8766D",peripheral="#00BFC4")
plots <- list();reports <- list()
case_indices <- if("--first-case" %in% commandArgs(TRUE))1L else seq_len(nrow(cases))
for(case_i in case_indices) {
 ca <- cases[case_i,];pies <- list();k <- 0L
 for(i in seq_len(nrow(row_specs))) for(j in seq_along(radii)) {
  rr <- row_specs[i,]
  d <- percent[percent$ssp==ca$ssp & percent$period==ca$period & percent$group==rr$group & percent$direction==rr$direction & percent$radius_km==radii[j] & percent$n>0,]
  stopifnot(abs(sum(d$percentage_of_focal_populations)-100)<1e-8,sum(d$n)==if(rr$group=="core")9 else 10)
  d$endpoint_habitat_class <- factor(d$endpoint_habitat_class,levels=names(palette));d <- d[order(d$endpoint_habitat_class),]
  d$label <- sprintf("%.1f%%",d$percentage_of_focal_populations)
  d$fraction <- d$n/sum(d$n);d$end <- cumsum(d$fraction)*2*pi+pi/2;d$start <- c(pi/2,head(d$end,-1));d$middle <- (d$start+d$end)/2
  polygons <- do.call(rbind,lapply(seq_len(nrow(d)),function(z){angle <- seq(d$start[z],d$end[z],length.out=max(20,ceiling(400*d$fraction[z])));data.frame(x=c(0,cos(angle),0),y=c(0,sin(angle),0),category=d$endpoint_habitat_class[z],wedge=z)}))
  if(nrow(d)==1){angle <- seq(0,2*pi,length.out=401);polygons <- data.frame(x=cos(angle),y=sin(angle),category=d$endpoint_habitat_class,wedge=1)}
  labels <- data.frame(x=.54*cos(d$middle),y=.54*sin(d$middle),label=d$label,small=d$fraction<1)
  leaders <- data.frame(x=cos(d$middle[labels$small]),y=sin(d$middle[labels$small]),xend=1.02*sign(cos(d$middle[labels$small])),yend=1.12*sin(d$middle[labels$small]))
  labels$x[labels$small] <- 1.14*sign(cos(d$middle[labels$small]));labels$y[labels$small] <- 1.12*sin(d$middle[labels$small])
  labels$hjust <- ifelse(labels$small,ifelse(labels$x>0,0,1),.5)
  if(nrow(d)==1){labels$x <- 0;labels$y <- 0}
  k <- k+1L
  pies[[k]] <- ggplot()+geom_polygon(data=polygons,aes(x,y,group=wedge,fill=category),colour="white",linewidth=.4,show.legend=FALSE)+geom_segment(data=leaders,aes(x,y,xend=xend,yend=yend),linewidth=.25,colour="#555555")+geom_text(data=labels,aes(x,y,label=label,hjust=hjust),size=2.82)+scale_fill_manual(values=palette,limits=names(palette),drop=FALSE)+coord_equal(xlim=c(-2.5,2.5),ylim=c(-1.5,1.5),expand=FALSE,clip="off")+labs(subtitle=if(i==1)if(radii[j]=="unrestricted")"Unlimited" else paste(radii[j],"km") else " ",tag=if(j==1)rr$panel else NULL)+theme_void(base_size=8)+theme(plot.subtitle=element_text(size=8,hjust=.5),plot.tag=element_text(size=12),plot.tag.position=c(0,1),plot.margin=margin(7,3,7,3))
 }
 key_data <- data.frame(x=c(.30,.62),y=0,category=names(palette),label=c("Core","Peripheral"))
 key <- ggplot(key_data,aes(x,y))+geom_point(aes(colour=category),shape=15,size=3,show.legend=FALSE)+geom_text(aes(label=label),nudge_x=.025,hjust=0,size=2.82)+scale_colour_manual(values=palette)+coord_cartesian(xlim=c(0,1),ylim=c(-1,1),expand=FALSE)+theme_void()
 plot <- wrap_plots(c(pies,list(key)),design="ABCDE\nFGHIJ\nKLMNO\nPQRST\nUUUUU",heights=c(1,1,1,1,.18))
 stem <- file.path(directory,paste0("Figure_S",15L+case_i))
 ids <- paste0("pie",seq_len(20));row_groups <- lapply(0:3,function(i)ids[i*5+1:5]);col_groups <- lapply(1:5,function(j)ids[seq(j,20,by=5)])
 require_patchwork_panel_alignment(plot,paste0(stem,".alignment-layout.json"),paste0(stem,".alignment.json"),width_in=190/25.4,height_in=200/25.4,panel_ids=c(ids,"class_key"),row_groups=row_groups,column_groups=col_groups,audit_script=file.path(skill,"scripts/audit_panel_alignment.py"),python=python,strict=TRUE)
 ggsave(paste0(stem,".pdf"),plot,width=190,height=200,units="mm",device=grDevices::pdf,useDingbats=FALSE,bg="white")
 ggsave(paste0(stem,".png"),plot,width=190,height=200,units="mm",dpi=300,type="cairo",bg="white")
 svglite::svglite(paste0(stem,".svg"),width=190/25.4,height=200/25.4);print(plot);dev.off()
 saveRDS(plot,paste0(stem,".rds"))
 text_status <- system2(python,c(shQuote(file.path(skill,"scripts/audit_pdf_text.py")),shQuote(paste0(stem,".pdf")),"--min-pt","7","--json"),stdout=paste0(stem,".text-audit.json"))
 collision_status <- system2(python,c(shQuote(file.path(skill,"scripts/audit_figure_collisions.py")),shQuote(paste0(stem,".pdf")),"--json-out",shQuote(paste0(stem,".collision-audit.json")),"--overlay-pdf",shQuote(paste0(stem,".collision-audit.pdf"))))
 stopifnot(text_status==0,collision_status==0)
 collision <- jsonlite::read_json(paste0(stem,".collision-audit.json"));stopifnot(collision$summary$fail==0)
 plots[[case_i]] <- plot
 reports[[case_i]] <- list(ssp=ca$ssp,period=ca$period,stem=stem,collision_fail=collision$summary$fail,collision_warn=collision$summary$warn)
 writeLines(paste0("Figure S",15L+case_i,". Distance-limited matching endpoints under ",if(ca$ssp=="ssp126")"SSP1-2.6"else"SSP5-8.5"," in ",ca$period,". Rows (a–b) show forward matching locations for core and peripheral focal populations, respectively; rows (c–d) show reverse present-day candidate source locations for core and peripheral focal populations. Columns restrict geographic distance to 10, 50, 100, 200 km or no limit. Salmon denotes core and turquoise peripheral endpoint class. Exact sampling coordinates use the published population class; other habitat-grid endpoints use polygon class. Each pie counts one endpoint per focal population (nine core or ten peripheral). Matching includes geographic distance in the score. Local self matches are eligible; endpoint percentages do not estimate migration probabilities."),paste0(stem,".caption.txt"))
}

jsonlite::write_json(list(rendered=TRUE,figures=16:19,models_refitted=FALSE,source="endpoint_percentages.csv",manual_QA="PENDING"),file.path(directory,"validation.json"),pretty=TRUE,auto_unbox=TRUE)
