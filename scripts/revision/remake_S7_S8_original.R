# Original-design S7/S8; render validated saved results only, no model fitting.
source(file.path('scripts','revision','common.R'))
base_out <- file.path(output_dir,'figure_revision_original_design','supplement_S7_S8')
dir.create(base_out,recursive=TRUE,showWarnings=FALSE)
skill <- Sys.getenv("AMMO_FIGURE_QA_DIR", unset="")
python <- Sys.getenv("AMMO_PYTHON", unset="python3")
fit <- readRDS(file.path(output_dir,'gdm','model_geo_TRUE.rds'))
curves <- read.csv(file.path(output_dir,'figures','gdm_functions','ispline_source_data.csv'))
curves <- curves[curves$model=='geo_TRUE',]
contrib <- read.csv(file.path(output_dir,'gdm','fitted_predictor_contributions.csv'))
contrib <- contrib[contrib$model=='geo_TRUE',]
stopifnot(nrow(contrib)==6,length(fit$observed)==171)
write_csv(curves,file.path(base_out,'S7_saved_isplines.csv'))
write_csv(contrib,file.path(base_out,'S7_saved_spline_heights.csv'))
write_csv(data.frame(observed=fit$observed,predicted=fit$predicted,ecological=fit$ecological),file.path(base_out,'S7_saved_fit_diagnostics.csv'))
d <- read.csv(file.path(output_dir,'figures','gdm_transforms','transformed_grid_and_RGB.csv'))
stopifnot(nrow(d)==157745,all(is.finite(as.matrix(d[c('longitude','latitude')]))))
# Exact numerical inputs copied, rather than recomputing the saved PCA/signs.
file.copy(file.path(output_dir,'figures','gdm_transforms','transformed_grid_and_RGB.csv'),file.path(base_out,'S8_saved_transformed_grid_and_RGB.csv'),overwrite=TRUE)
record <- function(id) {
 f<-par('fig');p<-par('plt'); w<-dev.size('in')*72
 bbox<-c(f[1]+p[1]*diff(f[1:2]),f[3]+p[3]*diff(f[3:4]),f[1]+p[2]*diff(f[1:2]),f[3]+p[4]*diff(f[3:4]))*rep(w,2)
 list(id=id,bbox_pt=unname(bbox))
}
render <- function(name,w,h,draw,rows,cols) {
 stem<-file.path(base_out,name)
 pdf(paste0(stem,'.pdf'),width=w/25.4,height=h/25.4,pointsize=8,family='sans',useDingbats=FALSE)
 geometry<-draw();dev.off()
 jsonlite::write_json(list(schema_version=1,backend='r-base',figure=list(width_pt=w/25.4*72,height_pt=h/25.4*72),panels=geometry,row_groups=rows,column_groups=cols,measurement='Native base-R par(fig)/par(plt) final-device rectangles'),paste0(stem,'.alignment-layout.json'),auto_unbox=TRUE,pretty=TRUE,digits=12)
 stopifnot(system2(python,c(shQuote(file.path(skill,'scripts','audit_panel_alignment.py')),shQuote(paste0(stem,'.alignment-layout.json')),'--json-out',shQuote(paste0(stem,'.alignment.json'))))==0)
 svglite::svglite(paste0(stem,'.svg'),width=w/25.4,height=h/25.4,pointsize=8);draw();dev.off()
 png(paste0(stem,'.png'),width=w,height=h,units='mm',res=600,type='cairo',pointsize=8);draw();dev.off()
 stopifnot(system2(python,c(shQuote(file.path(skill,'scripts','audit_pdf_text.py')),shQuote(paste0(stem,'.pdf')),'--min-pt','8','--json'),stdout=paste0(stem,'.text-audit.json'))==0)
 status<-system2(python,c(shQuote(file.path(skill,'scripts','audit_figure_collisions.py')),shQuote(paste0(stem,'.pdf')),'--json-out',shQuote(paste0(stem,'.collision-audit.json')),'--overlay-pdf',shQuote(paste0(stem,'.collision-audit.pdf'))))
 writeLines(as.character(status),paste0(stem,'.collision-exit-code.txt'))
 stopifnot(status==0)
}
drawS7 <- function() {
 layout(rbind(c(1,2,3,4),c(1,5,6,7),c(0,8,9,0)),widths=c(1.7,1,1,1))
 par(cex=1,cex.axis=1,cex.lab=1,cex.main=1,font.main=1,mgp=c(1.7,.45,0),tcl=-.25,mar=c(3.5,3.4,2.2,.7),xaxs='i',yaxs='i')
 o<-c('bio15','bio17','bio4','bio5','bio6','Geographic');q<-contrib[match(o,contrib$predictor),]
 barplot(q$fitted_spline_height,names.arg=c('BIO15','BIO17','BIO4','BIO5','BIO6','Geographic'),col='steelblue',border=NA,las=2,xlim=c(0,7.5),ylim=c(0,max(q$fitted_spline_height)*1.15),ylab='Fitted spline height',cex.names=1)
 box();mtext('a',side=3,adj=0,line=.45,cex=1.5);g<-list(record('height'))
 par(mar=c(3.2,3.2,2.2,.8))
 plot(fit$ecological,fit$observed,pch=16,cex=.35,col='blue',xlab='Predicted distance',ylab='Observed dissimilarity',axes=FALSE);axis(1);axis(2);box()
 xx<-seq(0,max(fit$ecological),length.out=200);lines(xx,1-exp(-xx),lwd=1)
 mtext('b',side=3,adj=0,line=.45,cex=1.5);g[[2]]<-record('diagnostic1')
 plot(fit$predicted,fit$observed,pch=16,cex=.35,col='blue',xlab='Predicted dissimilarity',ylab='Observed dissimilarity',axes=FALSE);axis(1);axis(2);box();abline(0,1)
 g[[3]]<-record('diagnostic2')
 for(pred in c('bio4','bio17','bio5','Geographic','bio15','bio6')) {
  z<-curves[curves$predictor==pred,];stopifnot(nrow(z)>0)
  plot(z$x,z$transformed,type='l',lwd=1,col='black',xlab=if(pred=='Geographic') 'Distance (degrees)' else toupper(pred),ylab='Fitted transformation',axes=FALSE)
  axis(1,at=pretty(z$x,3));axis(2,at=pretty(z$transformed,3));box();g[[length(g)+1]]<-record(pred)
 }
 g
}
# Scalar transform colors reproduce the original blue-cyan-green-yellow-red map sequence.
ramp <- grDevices::colorRampPalette(c('blue','cyan','green','yellow','red'))(100)
cell<-1/120;xi<-round((d$longitude-min(d$longitude))/cell)+1L;yi<-round((d$latitude-min(d$latitude))/cell)+1L
xr<-min(d$longitude)+seq(0,max(xi)-1)*cell;yr<-min(d$latitude)+seq(0,max(yi)-1)*cell
stopifnot(max(abs(xr[xi]-d$longitude))<1e-7,max(abs(yr[yi]-d$latitude))<1e-7)
make_matrix<-function(v){z<-matrix(NA_real_,length(xr),length(yr));z[cbind(xi,yi)]<-v;z}
mapvars<-c('Geographic_X','Geographic_Y','BIO15','BIO17','BIO4','BIO5','BIO6')
mapnames<-c('Longitude','Latitude','BIO15','BIO17','BIO4','BIO5','BIO6')
matrices<-lapply(mapvars,function(v)make_matrix(d[[v]]))
rgbmatrix<-matrix('#FFFFFF',length(xr),length(yr));rgbmatrix[cbind(xi,yi)]<-d$RGB
# RGB raster stored as a genuine geographic layer; axes and panel letters remain vector/text.
drawS8 <- function() {
 layout(rbind(c(1,2,3,8,8),c(4,5,6,8,8),c(7,0,0,8,8)),widths=c(1,1,1,1.35,1.35))
 par(cex=1,cex.axis=1,cex.lab=1,cex.main=1,font.main=1,mgp=c(1.4,.3,0),tcl=-.2,mar=c(2.3,2.2,2.8,2.6))
 g<-list()
 for(i in seq_along(mapvars)) {
  v<-d[[mapvars[i]]];lim<-range(v);z<-matrices[[i]]
  image(xr,yr,z,col=ramp,zlim=lim,useRaster=TRUE,xlab='',ylab='',axes=FALSE,xlim=c(99.5,110.5),ylim=c(36.8,43.5))
  axis(1,at=c(100,105,110));axis(2,at=c(37,40,43),las=1);box();mtext(mapnames[i],side=3,line=.3)
  if(i==1)mtext('a',side=3,adj=0,line=1.2,cex=1.5)
  # Original narrow adjacent colour scale, with numeric endpoints in fitted units.
  u<-par('usr');dx<-diff(u[1:2]);dy<-diff(u[3:4]);bottom<-u[3]+.15*dy;top<-u[3]+.85*dy
  par(xpd=NA);yy<-seq(bottom,top,length.out=101)
  rect(u[2]+.06*dx,yy[-101],u[2]+.11*dx,yy[-1],col=ramp,border=NA)
  text(u[2]+.16*dx,c(bottom,top),labels=format(signif(lim,2),trim=TRUE),adj=0,cex=1);par(xpd=FALSE)
  g[[i]]<-record(mapvars[i])
 }
 par(fig=c(.53,.995,.24,.965),new=TRUE,mar=c(2.3,2.2,2.8,.5))
 plot.new();plot.window(xlim=c(99.5,110.5),ylim=c(36.8,43.5),xaxs='i',yaxs='i')
 rasterImage(as.raster(t(rgbmatrix[,ncol(rgbmatrix):1,drop=FALSE])),min(xr)-cell/2,min(yr)-cell/2,max(xr)+cell/2,max(yr)+cell/2,interpolate=FALSE)
 axis(1,at=c(100,105,110));axis(2,at=c(37,40,43),las=1);box();mtext('b',side=3,adj=0,line=1.2,cex=1.5)
 g[[8]]<-record('RGB');g
}
render('Figure_S7',190,150,drawS7,rows=list(c('diagnostic1','diagnostic2','bio4'),c('bio17','bio5','Geographic'),c('bio15','bio6')),cols=list(c('diagnostic1','bio17','bio15'),c('diagnostic2','bio5','bio6'),c('bio4','Geographic')))
render('Figure_S8',190,110,drawS8,rows=list(mapvars[1:3],mapvars[4:6]),cols=list(mapvars[c(1,4,7)],mapvars[c(2,5)],mapvars[c(3,6)]))
writeLines(c('Figure S7. Saved GDM fit for 171 candidate-SNP population pairs. (a) Fitted spline heights; these are fitted model contributions, not permutation-based predictive importance or causal effects. (b) Original diagnostic and I-spline arrangement: observed versus ecological distance, observed versus predicted dissimilarity, BIO4; BIO17, BIO5, geographic distance; BIO15, BIO6. Lines show the saved fitted transformations without uncertainty intervals. Geographic distance uses inherited coordinate-degree units.','BIO4: temperature seasonality; BIO5: maximum temperature of the warmest month; BIO6: minimum temperature of the coldest month; BIO15: precipitation seasonality; BIO17: precipitation of the driest quarter.'),file.path(base_out,'Figure_S7.caption.txt'))
writeLines(c('Figure S8. Saved fitted GDM transformations across 157,745 current-habitat grid cells. (a) Longitude, latitude and five climatic transformations, in original order; each blue-to-red scale spans the individual observed fitted-transformation range. Axes show longitude and latitude in degrees. (b) Saved centered, unscaled PCA of the seven transformed columns: independently min-max-scaled PC1, PC2 and PC3 control red, green and blue. Similar colors indicate similarity along these display axes; they do not identify genetic clusters. This PCA composite is distinct from the local/forward/reverse offset RGB shown in Figure 3.','BIO4: temperature seasonality; BIO5: maximum temperature of the warmest month; BIO6: minimum temperature of the coldest month; BIO15: precipitation seasonality; BIO17: precipitation of the driest quarter. No GDM refit or PCA rerun was undertaken.'),file.path(base_out,'Figure_S8.caption.txt'))
saveRDS(list(curves=curves,contributions=contrib,diagnostics=data.frame(observed=fit$observed,predicted=fit$predicted,ecological=fit$ecological),mapdata=d),file.path(base_out,'source_plotting_data.rds'))
save_session(base_out)
