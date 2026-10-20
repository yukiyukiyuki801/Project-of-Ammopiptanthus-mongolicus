#!/usr/bin/env Rscript
# Export saved posterior plotting coordinates without fitting any model.
args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(if(length(args)) args[[1]] else '.', mustWork=TRUE)
out <- if(length(args)>1) args[[2]] else file.path(root,'runs/revision_run_01/outputs/eems_existing_audit')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
input <- file.path(root,'ammo_analysis_package_20261014/ammo/ammo_all_deme700/ammo_all_deme700-rEEMSplots-rdist.RData')
e <- new.env(parent=baseenv()); load(input,envir=e)
stopifnot(all(c('xym.values','xyq.values','B.component','G.component','W.component') %in% ls(e)))
summary <- list()
for(n in c('xym.values','xyq.values','B.component','G.component','W.component')) {
  d <- as.data.frame(e[[n]])
  write.csv(d,file.path(out,paste0(n,'.csv')),row.names=FALSE)
  summary[[n]] <- list(rows=nrow(d),columns=names(d))
  if(n %in% c('xym.values','xyq.values')) {
    stopifnot(ncol(d)==3,all(is.finite(as.matrix(d))),!anyDuplicated(d[,1:2]))
    summary[[n]]$longitude_range <- range(d[[1]])
    summary[[n]]$latitude_range <- range(d[[2]])
    summary[[n]]$value_range <- range(d[[3]])
    summary[[n]]$regular_grid_complete <- nrow(d)==length(unique(d[[1]]))*length(unique(d[[2]]))
  }
}
stopifnot(identical(e[['xym.values']][,1:2],e[['xyq.values']][,1:2]))
jsonlite::write_json(list(validated=TRUE,source=input,source_md5=unname(tools::md5sum(input)),
 scope='Existing 700-deme plotting-data export only; no model fit and no convergence claim',
 scale='Official plotting source labels these log10 rates after spatial mean centering under defaults; exact original plotting call not recovered',
 objects=summary),file.path(out,'extraction_validation.json'),pretty=TRUE,auto_unbox=TRUE)
