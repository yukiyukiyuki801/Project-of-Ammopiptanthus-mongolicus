# New, explicitly documented PCA-env reconstruction; never an exact historical rerun.
source("scripts/revision/common.R")
suppressPackageStartupMessages({library(terra); library(ade4); library(jsonlite)})
d <- file.path(output_dir, "niche_reconstruction")
dir.create(d, recursive=TRUE, showWarnings=FALSE)
seed <- 20261014L
set.seed(seed)
source_dir <- file.path(root,"public_dependencies/software/ecospat_verified")
source(file.path(source_dir,"ecospat.nichedynamic.R"))
source(file.path(source_dir,"ecospat.nicheoverlap.R"))
replications <- 999L
grid_resolution <- 100L
background_n <- 1000L
settings <- list(analysis="New PCA-env reconstruction, not recovery of historical niche tests",
 authorization="User explicitly requested rerunning missing niche results with their data; original S4 numerical grids/nulls absent; reviewer R1.24 asks metric/null interpretation and R1.66 asks readable niche figure colour coding.",
 seed=seed, repetitions=replications, grid_resolution=grid_resolution,
 variables=bio_names, background_n=background_n,
 background="One shared regional accessible-background approximation: random equal-cell sample of complete five-variable current CHELSA raster cells across the supplied 100-110 E / 36-44 N regional extent. This is an explicit geographic assumption, not a recovered accessible-area estimate; no suitability threshold or group-specific climate clipping.",
 units="19 equally weighted verified population coordinates: nine published core, ten peripheral. The 217 individual records are not treated as independent ecological occurrences. Other recovered 2022 occurrences lack verified population-class identities and are not assigned classes here.",
 pca="Centred and scaled PCA fitted once on the shared background; two axes retained; population values projected using fitted centring/scaling/loadings.",
 kernel="Official ecospat source functions, adehabitat bivariate Gaussian href kernel, th.sp=0, th.env=0; R=100; environmental-prevalence correction enabled for D/I.",
 equivalency="Pooled population-occurrence labels randomly reassigned preserving 9/10 group sizes; lower-tail overlap alternative (observed overlap lower than random relabelling).",
 similarity="Directional rand.type=2 tests: first niche fixed, second niche randomly shifted relative to its background; higher-tail overlap alternative. Both directions have the same shared background; direction does not indicate dispersal or adaptation.",
 intersection=0, threads=1,
 limitations="Sparse spatially structured sample (19 populations), assumed shared rectangular regional background, and two-axis projection restrict inference. Nonsignificant equivalency is not evidence that niches are identical; significant similarity is not evidence of adaptation, gene flow or rescue.")
write_json(settings,file.path(d,"settings.json"),pretty=TRUE,auto_unbox=TRUE)
rasters <- file.path(bundle,"ammo/gf/chelsa2_crop_selected/current",paste0("b",bio_numbers,".asc"))
stopifnot(all(file.exists(rasters)),nrow(pop_meta)==19,
 sum(pop_meta$group=="core")==9,sum(pop_meta$group=="peripheral")==10,
 !anyDuplicated(pop_meta$population))
r <- rast(rasters); names(r) <- bio_names
v <- values(r,mat=TRUE)
valid <- which(rowSums(is.finite(v))==length(bio_names))
stopifnot(length(valid)>=background_n)
bg_cells <- sort(sample(valid,background_n,replace=FALSE))
bg <- as.data.frame(v[bg_cells,,drop=FALSE]); names(bg) <- bio_names
bg_xy <- xyFromCell(r,bg_cells)
write_csv(cbind(cell=bg_cells,longitude=bg_xy[,1],latitude=bg_xy[,2],bg),file.path(d,"background_climate.csv"))
pop_clim <- extract(r,as.matrix(pop_meta[,c("longitude","latitude")]))
stopifnot(all(is.finite(as.matrix(pop_clim))),ncol(pop_clim)==5)
names(pop_clim) <- bio_names
write_csv(cbind(pop_meta,pop_clim),file.path(d,"population_climate.csv"))
rm(v); gc()
pca <- prcomp(bg,center=TRUE,scale.=TRUE)
bg_scores <- as.data.frame(predict(pca,bg)[,1:2,drop=FALSE])
pop_scores <- as.data.frame(predict(pca,pop_clim)[,1:2,drop=FALSE])
# Global limits include projected occurrence points, even if outside sampled background extremes.
glob <- rbind(bg_scores,pop_scores)
core <- pop_scores[pop_meta$group=="core",,drop=FALSE]
peri <- pop_scores[pop_meta$group=="peripheral",,drop=FALSE]
write_csv(cbind(pop_meta,pop_scores),file.path(d,"population_pca_scores.csv"))
write_csv(data.frame(variable=rownames(pca$rotation),pca$rotation),file.path(d,"pca_loadings.csv"))
write_csv(data.frame(axis=seq_along(pca$sdev),proportion=pca$sdev^2/sum(pca$sdev^2)),file.path(d,"pca_variance.csv"))
z1 <- ecospat.grid.clim.dyn(glob=glob,glob1=bg_scores,sp=core,R=grid_resolution,th.sp=0,th.env=0)
z2 <- ecospat.grid.clim.dyn(glob=glob,glob1=bg_scores,sp=peri,R=grid_resolution,th.sp=0,th.env=0)
observed <- ecospat.niche.overlap(z1,z2,cor=TRUE)
stopifnot(all(is.finite(unlist(observed))),all(unlist(observed)>=0),all(unlist(observed)<=1))
saveRDS(list(pca=pca,background_scores=bg_scores,population_scores=pop_scores,
 core=terra::wrap(z1$z.cor),peripheral=terra::wrap(z2$z.cor),observed=observed,settings=settings),file.path(d,"niche_plot_objects.rds"))
for (nm in c("core","peripheral")) {
 z <- if(nm=="core") z1 else z2
 writeRaster(z$z.cor,file.path(d,paste0(nm,"_density.tif")),overwrite=TRUE)
 xy <- xyFromCell(z$z.cor,seq_len(ncell(z$z.cor)))
 write_csv(data.frame(PC1=xy[,1],PC2=xy[,2],
  corrected_density=as.vector(values(z$z.cor))),file.path(d,paste0(nm,"_density.csv")))
}
write_csv(data.frame(metric=names(observed),value=unlist(observed)),file.path(d,"observed_overlap.csv"))
test_defs <- list(equivalency=list(kind="equivalency",a=z1,b=z2,alternative="lower",seed=seed+1L),
 similarity_core_fixed=list(kind="similarity",a=z1,b=z2,alternative="higher",seed=seed+2L),
 similarity_peripheral_fixed=list(kind="similarity",a=z2,b=z1,alternative="higher",seed=seed+3L))
summaries <- list()
for(nm in names(test_defs)) {
 def <- test_defs[[nm]]; target <- file.path(d,paste0(nm,".rds"))
 if(file.exists(target)) {
  checkpoint <- readRDS(target)
  stopifnot(identical(checkpoint$settings,settings),identical(checkpoint$seed,def$seed))
  ans <- checkpoint$result
 } else {
  set.seed(def$seed)
  cat("Starting",nm,"with",replications,"randomizations\n"); flush.console()
  if(def$kind=="equivalency") ans <- ecospat.niche.equivalency.test(def$a,def$b,
   rep=replications,intersection=0,overlap.alternative=def$alternative,ncores=1)
  else ans <- ecospat.niche.similarity.test(def$a,def$b,
   rep=replications,intersection=0,rand.type=2,overlap.alternative=def$alternative,ncores=1)
  stopifnot(nrow(ans$sim)==replications,all(is.finite(ans$sim$D)),all(is.finite(ans$sim$I)))
  saveRDS(list(settings=settings,seed=def$seed,result=ans),target)
 }
 write_csv(ans$sim,file.path(d,paste0(nm,"_null.csv")))
 for(metric in c("D","I")) {
  obs <- unname(ans$obs[[metric]])
  null <- ans$sim[[metric]]
  # Explicit plus-one p-values; one prespecified tail for each test type.
  p <- if(def$alternative=="lower") (1+sum(null<=obs))/(replications+1) else (1+sum(null>=obs))/(replications+1)
  stopifnot(is.finite(p),abs(p-ans[[paste0("p.",metric)]])<1e-12)
  summaries[[paste(nm,metric)]] <- data.frame(test=nm,metric=metric,observed=obs,
   alternative=def$alternative,p=p,null_q025=unname(quantile(null,.025)),null_q975=unname(quantile(null,.975)),n=replications)
 }
}
summary <- do.call(rbind,summaries)
summary$q_BH_six_reported_tests <- p.adjust(summary$p,"BH")
write_csv(summary,file.path(d,"test_summary.csv"))
source_files <- c(rasters,file.path(output_dir,"inputs/population_metadata.csv"),
 file.path(source_dir,c("ecospat.nichedynamic.R","ecospat.nicheoverlap.R")))
write_csv(data.frame(path=source_files,md5=unname(tools::md5sum(source_files))),file.path(d,"input_manifest.csv"))
save_session(d)
marker(d,list(n_populations=19,n_core=9,n_peripheral=10,background_cells=background_n,
 available_complete_raster_cells=length(valid),repetitions=replications,seed=seed,
 observed=observed,scope="New documented population-level niche reconstruction; original historical test objects unavailable; no SDM fitting or EEMS inference."))
