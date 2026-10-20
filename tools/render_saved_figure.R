#!/usr/bin/env Rscript
# Run the paper's drawing code from saved data, with optional separately licensed QA.
# This wrapper does not fit models or change plot geometry, text, colours or data.
a <- commandArgs(trailingOnly=TRUE)
if (!length(a)) stop('Usage: Rscript tools/render_saved_figure.R FIGURE [renderer arguments]')
repo <- normalizePath(getwd())
manifest <- jsonlite::read_json(file.path(repo,'provenance/figure_inventory.json'), simplifyVector=FALSE)
idx <- which(vapply(manifest$figures,function(x)identical(as.character(x$figure),a[1]),logical(1)))
if(length(idx)!=1) stop('Unknown figure identifier')
f <- manifest$figures[[idx]]; src <- file.path(repo,f$script)
qa_dir <- Sys.getenv('AMMO_FIGURE_QA_DIR',unset='')
qa_files <- c('panel_alignment.R','audit_panel_alignment.py','audit_pdf_text.py','audit_figure_collisions.py')
qa_available <- nzchar(qa_dir) && all(file.exists(file.path(qa_dir,'scripts',qa_files)))
if(!qa_available) warning('Optional figure QA NOT_RUN: separately supplied nature-figure helpers are unavailable. Saved scientific results and figure drawing continue.')
Sys.setenv(AMMO_REVISION_RUN=file.path(repo,'runs/revision_run_01'))
# Runtime call guards apply only to optional QA. Parsing retains all drawing AST.
qa_pattern <- 'audit_panel_alignment|audit_pdf_text|audit_figure_collisions|audit_figure_collisions_chars|panel_alignment\\.R|\\.collision-audit|\\.text-audit|\\.alignment'
qa_names <- c('require_patchwork_panel_alignment','write_patchwork_panel_layout')
qa_vars <- character()
contains_qa <- function(e) grepl(qa_pattern,paste(deparse(e),collapse=' '))
scan_vars <- function(e) {
 if(!is.call(e))return(invisible(NULL))
 n<-paste(deparse(e[[1]]),collapse='')
 if(n %in% c('<-','=') && length(e)>=3 && is.symbol(e[[2]]) && contains_qa(e[[3]]) && !identical(as.character(e[[3]][[1]])[1],'function'))qa_vars <<- unique(c(qa_vars,as.character(e[[2]])))
 if(length(e)>1)for(i in seq.int(2,length(e))){if(identical(e[[i]],quote(expr=)))next;scan_vars(e[[i]])}
 invisible(NULL)
}
expr <- parse(src);for(e in expr)scan_vars(e)
guard <- function(e) {
 if(!is.call(e))return(e)
 n<-paste(deparse(e[[1]]),collapse='')
 # Entire QA assertions are conditional; no synthetic success status is returned.
 qa_assert<-n=='stopifnot' && (contains_qa(e) || any(vapply(qa_vars,function(v)v %in% all.names(e),logical(1))))
 qa_source<-n=='source' && contains_qa(e)
 qa_read<-n %in% c('jsonlite::read_json','read_json','fromJSON','jsonlite::fromJSON') && contains_qa(e)
 qa_write<-n %in% c('jsonlite::write_json','write_json') && any(vapply(qa_vars,function(v)v %in% all.names(e),logical(1)))
 if(n %in% qa_names || qa_assert || qa_source || qa_read || qa_write || (n %in% c('system2','system') && contains_qa(e)))return(substitute(if(qa_available) ORIGINAL else invisible(NULL),list(ORIGINAL=e)))
 if(length(e)>1)for(i in seq.int(2,length(e))){if(identical(e[[i]],quote(expr=)))next;e[i]<-list(guard(e[[i]]))}
 e
}
# commandArgs is used by original renderers; give them only intended arguments.
args <- a[-1]; render_env<-new.env(parent=globalenv())
render_env$commandArgs<-function(trailingOnly=FALSE,...) if(trailingOnly) args else c('Rscript',paste0('--file=',src),'--args',args)
render_env$qa_available<-qa_available
# Independent minimal normalization needed by authored workspace measurement code.
render_env$.nature_alignment_groups<-function(groups)if(is.null(groups))NULL else lapply(groups,as.character)
for(e in expr)eval(guard(e),envir=render_env)
status_dir<-file.path(repo,'runs/revision_run_01/outputs/figure_optional_QA');dir.create(status_dir,recursive=TRUE,showWarnings=FALSE)
jsonlite::write_json(list(figure=a[1],status=if(qa_available)'OPTIONAL_CHECKS_REQUESTED' else 'NOT_RUN',reason=if(qa_available)'See actual individual audit outputs; this wrapper does not certify PASS.' else 'External QA helpers absent; original plotting calls ran without optional alignment/text/collision checks.'),file.path(status_dir,paste0('Figure_',a[1],'.json')),pretty=TRUE,auto_unbox=TRUE)
