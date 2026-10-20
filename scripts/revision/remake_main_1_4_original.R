# Original-design remakes from validated saved results only. No model fitting.
source(file.path('scripts','revision','common.R'))
library(ggplot2)
library(patchwork)
base_out <- file.path(output_dir,'figure_revision_original_design','main_1_4')
dir.create(base_out,recursive=TRUE,showWarnings=FALSE)
skill <- Sys.getenv("AMMO_FIGURE_QA_DIR", unset="")
source(file.path(skill,'scripts','panel_alignment.R'))
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
  viewport_listing <- grid::grid.ls(viewports=TRUE,grobs=FALSE,print=FALSE)
  panels <- lapply(seq_len(nrow(rows)),function(i) {
    row <- rows[i,]
    viewport <- paste0(row$name,".",row$t,"-",row$l,"-",row$b,"-",row$r)
    grid::upViewport(0)
    grid::seekViewport(viewport)
    if(grepl(";",row$name,fixed=TRUE)) {
      grid::downViewport("layout")
      parent_path <- paste0("ROOT::layout::",viewport,"::layout")
      choices <- which(viewport_listing$vpPath==parent_path &
                       grepl("^panel(?:-[0-9]+)*[.]",viewport_listing$name,perl=TRUE))
      if(length(choices)!=1L) stop("Cannot identify unique actual facet viewport: ",viewport)
      grid::downViewport(viewport_listing$name[choices])
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

qa_python <- Sys.getenv("AMMO_PYTHON", unset="python3")
palette <- c(core='#F8766D',peripheral='#00BFC4')
periods <- c('2041-2070','2071-2100')
ssps <- c('ssp126','ssp585')
base_theme <- theme_bw(base_size=9,base_family='sans') +
  theme(axis.text=element_text(size=8,colour='black'),axis.title=element_text(size=9),
        axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=.5),
        strip.text=element_text(size=9),strip.background=element_rect(fill='grey85',colour=NA),
        panel.grid.minor=element_blank(),panel.grid.major=element_line(colour='grey90',linewidth=.3),
        legend.title=element_blank(),legend.text=element_text(size=9),
        legend.key.size=grid::unit(4,'mm'),plot.margin=margin(3,3,3,3),
        plot.tag=element_text(size=12,face='plain'),plot.tag.position=c(0,1))
theme_set(base_theme)
export_figure <- function(plot,name,width,height,caption) {
  stem <- file.path(base_out,name)
  require_patchwork_panel_alignment(plot,manifest_path=paste0(stem,'.alignment-layout.json'),
    report_path=paste0(stem,'.alignment.json'),overlay_svg=paste0(stem,'.alignment.svg'),
    width_in=width/25.4,height_in=height/25.4,
    row_groups=if(name=='Figure_4') unname(split(letters[1:12],rep(1:6,each=2))) else NULL,
    column_groups=if(name=='Figure_4') lapply(1:6,function(i) letters[c(i,i+6)]) else NULL,tolerance_pt=1.5,gutter_tolerance_pt=1.5,
    audit_script=file.path(skill,'scripts','audit_panel_alignment.py'),python=qa_python,strict=TRUE)
  grDevices::pdf(paste0(stem,'.pdf'),width=width/25.4,height=height/25.4,family='sans',useDingbats=FALSE)
  print(plot);grDevices::dev.off()
  dir.create(file.path(base_out,'qa'),showWarnings=FALSE)
  svglite::svglite(file.path(base_out,'qa',paste0(name,'.editable-text.svg')),width=width/25.4,height=height/25.4)
  print(plot);grDevices::dev.off()
  ggsave(paste0(stem,'.png'),plot,width=width,height=height,units='mm',dpi=600,type='cairo')
  writeLines(caption,paste0(stem,'.caption.txt'))
  saveRDS(plot,paste0(stem,'.rds'))
  stopifnot(system2(qa_python,c(shQuote(file.path(skill,'scripts','audit_pdf_text.py')),
    shQuote(paste0(stem,'.pdf')),'--min-pt','8','--json'),stdout=paste0(stem,'.text-audit.json'))==0)
  collision_status <- system2(qa_python,c(shQuote(file.path(skill,'scripts','audit_figure_collisions.py')),
    shQuote(paste0(stem,'.pdf')),'--json-out',shQuote(paste0(stem,'.collision-audit.json')),
    '--overlay-pdf',shQuote(paste0(stem,'.collision-audit.pdf'))))
  writeLines(as.character(collision_status),paste0(stem,'.collision-exit-code.txt'))
  stopifnot(collision_status==0)
}
# Figure 1: the original avg_avg aggregation is the arithmetic mean across
# the four scenario-period ensemble RONA values for each population and BIO.
r <- read.csv(file.path(output_dir,'rona','rona_population_long.csv'))
r <- r[r$gcm=='ENSEMBLE' & r$K==4 & r$association_filter=='inherited_p_lt_1',]
stopifnot(nrow(r)==19*5*4,!anyDuplicated(r[c('population','bio','ssp','period')]))
stopifnot(all(table(r$population,r$bio)==4),all(is.finite(r$rona)))
avg <- aggregate(rona~population+group+bio,r,mean)
stopifnot(nrow(avg)==95,all(table(avg$bio,avg$group)==rep(c(9,10),each=5)))
avg$group <- factor(avg$group,levels=names(palette))
write_csv(avg,file.path(base_out,'Figure_1_source_population_four_scenario_means.csv'))
fig1_panels <- lapply(bio_names,function(bio) {
  d <- avg[avg$bio==bio,];d$BIO <- factor(toupper(bio))
  ggplot(d,aes(group,rona,fill=group)) +
    geom_boxplot(width=.65,linewidth=.4,outlier.size=1.1) +
    facet_wrap(~BIO) + scale_fill_manual(values=palette,drop=FALSE) +
    scale_y_continuous(breaks=scales::breaks_pretty(3),expand=expansion(mult=c(.04,.08))) +
    labs(x=NULL,y=if(bio=='bio4') 'RONA' else NULL)
})
fig1 <- wrap_plots(fig1_panels,ncol=5,guides='collect') & theme(legend.position='right')
export_figure(fig1,'Figure_1',190,70,c(
 'Figure 1. RONA distributions in core and peripheral populations for BIO4, BIO5, BIO6, BIO15 and BIO17.',
 'For each population and climatic variable, RONA is the arithmetic mean of the four SSP-by-period ensemble projections (SSP126 and SSP585; 2041–2070 and 2071–2100). Each ensemble projection averages five GCMs. Population identities are corrected; nine core and ten peripheral populations are biological units. Boxes show medians and interquartile ranges; whiskers extend to the most extreme observations within 1.5 times the interquartile range, and observations outside the whiskers are plotted individually. These averaged distributions are descriptive; no significance from individual scenarios is transferred to this average.',
 'BIO4: temperature seasonality; BIO5: maximum temperature of the warmest month; BIO6: minimum temperature of the coldest month; BIO15: precipitation seasonality; BIO17: precipitation of the driest quarter.'))
# Figure 4: retain the original three metric blocks and mean +/- SD bars.
# Full-habitat matching outputs are summarized at 19 focal populations.
g <- read.csv(file.path(output_dir,'landscape_statistics','population_offsets_wide.csv'))
g <- g[g$gcm=='ENSEMBLE' & g$search_scoring=='inherited_geographic_penalty_degrees',]
stopifnot(nrow(g)==19*4,!anyDuplicated(g[c('population','ssp','period')]));g$group <- factor(g$group,levels=names(palette))
metrics <- c('forward','reverse','local')
rows <- do.call(rbind,lapply(metrics,function(metric) {
 d <- g[c('population','group','ssp','period')];d$metric <- metric;d$value <- g[[metric]];d
}))
summary <- aggregate(value~group+ssp+period+metric,rows,function(x)c(mean=mean(x),sd=sd(x),n=length(x)))
summary <- cbind(summary[1:4],as.data.frame(summary$value));names(summary)[5:7] <- c('mean','sd','n')
stopifnot(all(summary$n[summary$group=='core']==9),all(summary$n[summary$group=='peripheral']==10))
tests <- read.csv(file.path(output_dir,'landscape_statistics','population_class_effects.csv'))
tests <- tests[tests$gcm=='ENSEMBLE' & tests$search_scoring=='inherited_geographic_penalty_degrees',]
stopifnot(nrow(tests)==12)
write_csv(rows,file.path(base_out,'Figure_4_source_population_offsets.csv'))
write_csv(summary,file.path(base_out,'Figure_4_source_bar_mean_sd.csv'))
write_csv(tests,file.path(base_out,'Figure_4_saved_BH_tests.csv'))
# Twelve plot areas preserve the original three paired metric blocks. Comparable
# row groups are period pairs within each metric/SSP, and comparable columns
# are SSP pairs; larger inter-block gutters intentionally contain y-axis labels.
# The original six-column geometry is retained:
# within each of three blocks, early/late periods are left/right; SSP126 is
# the upper row and SSP585 the lower row. Necessary short facet identifiers
# remain; there are no figure titles or metric headings above lettered panels.
fig4_panels <- list()
for(ssp in ssps) for(metric in metrics) for(period in periods) {
 d <- summary[summary$ssp==ssp & summary$period==period & summary$metric==metric,]
 test <- tests[tests$ssp==ssp & tests$period==period & tests$metric==metric,]
 stopifnot(nrow(d)==2,nrow(test)==1)
 lim <- max(summary$mean[summary$metric==metric]+summary$sd[summary$metric==metric])*1.5
 d$scenario <- paste(toupper(ssp),if(period==periods[1]) '2041-2070' else '2071-2100',sep='\n')
 label <- sprintf('q=%.3f',test$q_BH)
 p <- ggplot(d,aes(group,mean,fill=group)) + geom_col(width=.65,linewidth=.25,colour='grey30') +
   geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=.16,linewidth=.35) +
   facet_wrap(~scenario) + scale_fill_manual(values=palette,drop=FALSE) +
   scale_y_continuous(limits=c(0,lim),breaks=scales::breaks_pretty(3),labels=function(x) sprintf('%.3f',x),expand=expansion(mult=c(0,0))) +
   annotate('segment',x=1,xend=2,y=lim*.82,yend=lim*.82,linewidth=.3) +
   annotate('text',x=1.5,y=lim*.91,label=label,size=8/ggplot2::.pt) +
   labs(x=NULL,y=if(period==periods[1]) 'GDM offset' else NULL,
        tag=if(ssp==ssps[1] && period==periods[1]) letters[match(metric,metrics)] else NULL) +
   theme(strip.text=element_text(size=8),plot.tag.position=c(.03,.98),
         panel.grid.major.x=element_blank(),plot.background=element_blank(),
         axis.text.y=if(period==periods[1]) element_text(size=8,colour='black') else element_blank(),
         axis.ticks.y=if(period==periods[1]) element_line() else element_blank())
 fig4_panels[[length(fig4_panels)+1L]] <- p
}
fig4 <- wrap_plots(fig4_panels,ncol=6,guides='collect') & theme(legend.position='right')
export_figure(fig4,'Figure_4',190,115,c(
 'Figure 4. Population-level comparisons of (a) forward, (b) reverse and (c) local GDM offsets.',
 'Each metric block contains four scenario-period combinations: SSP126 above SSP585, and 2041–2070 left of 2071–2100. Bars show means across nine core or ten peripheral focal populations; error bars show one standard deviation. Salmon denotes core and turquoise peripheral populations. Results use unrestricted full-habitat candidate matching with the inherited geographic penalty in longitude/latitude degrees. This score is not a dispersal-distance model. Environmental-only matching results are retained separately in supporting results.',
 'Brackets show saved population-level rank-sum comparisons with Benjamini–Hochberg-adjusted probabilities across the three offset metrics within each scoring definition, SSP and period. Numeric q values replace significance-star legends. Population identities and test values are corrected; map cells and GCMs are not biological replicates.'))
save_session(base_out)
marker(base_out,list(figures=c('Figure_1.pdf','Figure_4.pdf'),scope='Original-design rendering of validated saved numeric results; no refits or new tests',manual_visual_QA='Pending parent inspection'))
