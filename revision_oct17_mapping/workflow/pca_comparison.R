#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1)
library(jsonlite);library(ggplot2);library(patchwork)
cfg<-fromJSON(args[1]);o<-cfg$output_root;dest<-file.path(o,'pca');dir.create(dest,showWarnings=FALSE)
old<-read.csv(cfg$existing_pca_scores,stringsAsFactors=FALSE)
oldvar<-read.csv(cfg$existing_pca_variance)
samples<-readLines(file.path(o,'variants','sample_order.txt'))
stopifnot(length(samples)==217,!anyDuplicated(samples),setequal(samples,old$sample_id))
x<-read.delim(file.path(o,'variants','dosage.tsv'),header=FALSE,stringsAsFactors=FALSE,check.names=FALSE)
stopifnot(ncol(x)==220,nrow(x)>=2)
gt<-as.matrix(x[,4:ncol(x)]);dimnames(gt)<-NULL
valid<-c('0/0','0|0','0/1','1/0','0|1','1|0','1/1','1|1','./.','.|.','.')
stopifnot(all(gt%in%valid))
g<-matrix(NA_real_,nrow=nrow(gt),ncol=ncol(gt));g[gt%in%c('0/0','0|0')]<-0;g[gt%in%c('0/1','1/0','0|1','1|0')]<-1;g[gt%in%c('1/1','1|1')]<-2
missingness<-data.frame(sample_id=samples,missing_fraction=colMeans(is.na(g)))
write.csv(missingness,file.path(dest,'new_arm_genotype_missingness.csv'),row.names=FALSE)
g<-t(g);colnames(g)<-x$V3;rownames(g)<-samples
mu<-colMeans(g,na.rm=TRUE);p<-mu/2;scale<-sqrt(2*p*(1-p));keep<-is.finite(scale)&scale>0
stopifnot(sum(keep)>=2);g<-g[,keep,drop=FALSE];mu<-mu[keep];scale<-scale[keep]
for(j in seq_len(ncol(g)))g[is.na(g[,j]),j]<-mu[j]
z<-sweep(sweep(g,2,mu,'-'),2,scale,'/')
fit<-prcomp(z,center=FALSE,scale.=FALSE,rank.=2);variance<-fit$sdev^2/sum(fit$sdev^2)
new<-old[match(samples,old$sample_id),];new$PC1<-fit$x[,1];new$PC2<-fit$x[,2]
stopifnot(all(is.finite(new$PC1)),all(table(new$group)[c('core','peripheral')]==c(104,113)))
write.csv(new,file.path(dest,'mongolicus_PCA_scores.csv'),row.names=FALSE)
write.csv(data.frame(component=seq_along(variance),variance_fraction=variance),file.path(dest,'mongolicus_PCA_variance.csv'),row.names=FALSE)
write.csv(data.frame(variant=colnames(g),mean_dosage=mu,scale=scale),file.path(dest,'locus_standardization.csv'),row.names=FALSE)
saveRDS(list(fit=fit,source_variants=x$V3[keep],sample_order=samples),file.path(dest,'mongolicus_PCA.rds'))
makeplot<-function(d,v,tag){
 ggplot(d,aes(PC1,PC2,colour=group,shape=group))+geom_point(size=1.7,alpha=.8)+
 scale_colour_manual(values=c(core='#F8766D',peripheral='#00BFC4'),breaks=c('core','peripheral'),labels=c('Core (104)','Peripheral (113)'))+
 scale_shape_manual(values=c(core=16,peripheral=17),breaks=c('core','peripheral'),labels=c('Core (104)','Peripheral (113)'))+
 labs(x=sprintf('PC1 (%.2f%%)',100*v[1]),y=sprintf('PC2 (%.2f%%)',100*v[2]),tag=tag,colour=NULL,shape=NULL)+
 theme_classic(base_size=11)+theme(axis.text=element_text(colour='black'),legend.position='bottom',plot.tag=element_text(face='bold',size=13))
}
fig<-(makeplot(old,oldvar$variance_fraction,'a')|makeplot(new,variance,'b'))+plot_layout(guides='collect')&theme(legend.position='bottom')
ggsave(file.path(dest,'Reference_PCA_comparison.pdf'),fig,width=7.2,height=3.7,device=pdf)
ggsave(file.path(dest,'Reference_PCA_comparison.png'),fig,width=7.2,height=3.7,dpi=400,bg='white')
writeLines(c('Reference-mapping sensitivity: population structure. (a) Existing reconstructed PCA from the inherited A. nanus-aligned 38,124-SNP panel. (b) Newly called and filtered A. mongolicus-aligned SNP panel. Colours and shapes denote published population-level core/peripheral assignments; identical 217 individuals are shown. Each PCA uses its own loci, axes and explained variance; axis signs are arbitrary. Both analyses use per-locus mean imputation, centring at 2p, and scaling by sqrt(2p(1-p)). Differing reference assemblies and reconstructed calling workflows mean this is a qualitative panel comparison, not an isolated estimate of reference bias. Similar structure would not establish robustness of functional-load or genomic-offset estimates.'),file.path(dest,'caption.txt'))
write_json(list(validated=TRUE,matched_individuals=217,retained_loci=ncol(g),manual_visual_QA='PENDING',analysis_interpretation='PENDING'),file.path(dest,'validation.json'),pretty=TRUE,auto_unbox=TRUE)
