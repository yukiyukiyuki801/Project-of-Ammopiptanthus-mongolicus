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
