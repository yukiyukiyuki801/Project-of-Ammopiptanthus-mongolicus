# Actual leaf viewports, including nested original two-plus-three layouts.
write_patchwork_panel_layout<-function(plot,manifest_path,width_in,height_in,panel_ids=NULL,row_groups=NULL,column_groups=NULL,exemptions=list()){
 g<-patchwork::patchworkGrob(plot);path<-tempfile(fileext='.pdf');grDevices::pdf(path,width=width_in,height=height_in,useDingbats=FALSE);on.exit({grDevices::dev.off();unlink(path)},add=TRUE)
 grid::grid.newpage();grid::grid.draw(g);grid::grid.force()
 v<-grid::grid.ls(viewports=TRUE,grobs=FALSE,print=FALSE)
 idx<-which(grepl('^panel(-[0-9]+)?\\.[0-9]+-[0-9]+-[0-9]+-[0-9]+$',v$name))
 full<-paste(v$vpPath[idx],v$name[idx],sep='::');idx<-idx[!duplicated(full)];full<-unique(full)
 boxes<-lapply(full,function(x){parts<-strsplit(x,'::',fixed=TRUE)[[1]];parts<-parts[parts!='ROOT'&nzchar(parts)];grid::upViewport(0);grid::downViewport(do.call(grid::vpPath,as.list(parts)));b<-grid::deviceLoc(grid::unit(0,'npc'),grid::unit(0,'npc'),valueOnly=TRUE);t<-grid::deviceLoc(grid::unit(1,'npc'),grid::unit(1,'npc'),valueOnly=TRUE);c(b$x,b$y,t$x,t$y)*72})
 boxes<-do.call(rbind,boxes);boxes<-boxes[order(-round(boxes[,4],2),boxes[,1]),,drop=FALSE];if(nrow(boxes)==4){cx<-(boxes[,1]+boxes[,3])/2;left<-which.min(cx);right<-which.max(cx);middle<-setdiff(seq_len(4),c(left,right));middle<-middle[order(-boxes[middle,4])];boxes<-boxes[c(left,middle,right),,drop=FALSE]}
 if(is.null(panel_ids))panel_ids<-letters[seq_len(nrow(boxes))];stopifnot(nrow(boxes)==length(panel_ids))
 # Topology of seven-panel original S3 uses separate two- and three-column grids.
 panels<-lapply(seq_len(nrow(boxes)),function(i){if(nrow(boxes)==4){main<-i%in%c(1,4);return(list(id=panel_ids[i],bbox_pt=unname(boxes[i,]),grid_id=if(main)'main-maps' else paste0('inset-',i),row_start=0L,row_stop=1L,col_start=if(i==4)2L else 0L,col_stop=if(i==4)3L else 1L))};list(id=panel_ids[i],bbox_pt=unname(boxes[i,]),grid_id='three-panels',row_start=0L,row_stop=1L,col_start=as.integer(i-1),col_stop=as.integer(i))})
 jsonlite::write_json(list(schema_version=1L,backend='r-patchwork',figure=list(width_pt=width_in*72,height_pt=height_in*72),panels=panels,row_groups=.nature_alignment_groups(row_groups),column_groups=.nature_alignment_groups(column_groups),exemptions=exemptions,measurement='Actual final-size native PDF leaf panel viewports, with original nested two-column versus three-column topology'),manifest_path,pretty=TRUE,auto_unbox=TRUE,null='null',digits=12)
}
