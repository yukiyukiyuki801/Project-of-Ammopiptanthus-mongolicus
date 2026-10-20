source(file.path("scripts", "revision", "common.R"))
library(terra)
directory <- file.path(output_dir, "climate")
dir.create(directory, recursive=TRUE, showWarnings=FALSE)
rasters <- file.path(bundle,"ammo","gf","chelsa2_crop_selected")
read_stack <- function(folder, native=FALSE) {
  files <- file.path(folder, if(native) paste0(bio_names,".tif") else paste0("b",bio_numbers,".asc"))
  stopifnot(all(file.exists(files)))
  x <- terra::rast(files);names(x) <- bio_names
  x
}
extract_at <- function(x, meta, method="simple") {
  result <- terra::extract(x, meta[,c("longitude","latitude")], method=method, ID=FALSE)
  stopifnot(nrow(result)==nrow(meta), all(is.finite(as.matrix(result))))
  result
}
current <- read_stack(file.path(rasters,"current"))
current_sample <- extract_at(current,sample_meta)
current_pop <- extract_at(current,pop_meta)
manifest <- as.matrix(read.table(file.path(output_dir,"inputs","current_environment_manifest.tsv")))
stopifnot(ncol(manifest)==5)
checks <- data.frame(bio=bio_names, max_abs_difference=apply(abs(as.matrix(current_sample)-manifest),2,max),
                     median_raster=apply(current_sample,2,median),median_manifest=apply(manifest,2,median))
write_csv(checks,file.path(directory,"current_raster_vs_manifest.csv"))
write_csv(cbind(sample_meta,current_sample),file.path(directory,"current_sample_climate.csv"))
write_csv(cbind(pop_meta,current_pop),file.path(directory,"current_population_climate.csv"))
write_csv(cbind(sample_meta,setNames(as.data.frame(manifest),bio_names)),file.path(directory,"current_manifest_climate.csv"))
# Gross unit mismatches are an inference blocker. Small grid/extraction
# differences remain explicit; the verified manifest is used for LFMM/RONA.
stopifnot(all(abs(checks$median_raster-checks$median_manifest)<pmax(1,abs(checks$median_manifest)*0.1)))
future_sample <- list();future_pop <- list();archive_comparison <- list()
index <- 0L
gcms <- c("GFDL-ESM4","IPSL-CM6A-LR","MPI-ESM1-2-HR","MRI-ESM2-0","UKESM1-0-LL")
for (ssp in c("ssp126","ssp585")) for(short in c("2040","2070")) {
  period <- if(short=="2040") "2041-2070" else "2071-2100"
  ensemble <- read_stack(file.path(rasters,paste0(ssp,"_",short)))
  for(gcm in c("ENSEMBLE",gcms)) {
    native <- gcm!="ENSEMBLE"
    stack <- if(native) read_stack(file.path(root,"public_dependencies","climate_individual","CHELSA_v2.1_regional",period,gcm,ssp),TRUE) else ensemble
    es <- extract_at(stack,sample_meta); ep <- extract_at(stack,pop_meta)
    index <- index+1L
    future_sample[[index]] <- cbind(sample_meta,ssp=ssp,period=period,gcm=gcm,es)
    future_pop[[index]] <- cbind(pop_meta,ssp=ssp,period=period,gcm=gcm,ep)
    message("Extracted ",ssp," ",period," ",gcm)
    if(!native) {
      path <- file.path(bundle,"ammo","rona",paste0("ammo_rona_",ssp,"_",short,if(short=="2040") "-2070.txt" else "-2100.txt"))
      old <- read.table(path)
      stopifnot(nrow(old)==217,ncol(old)==6)
      # Neither row count nor old population labels establish individual ID.
      archive_comparison[[length(archive_comparison)+1L]] <- data.frame(ssp=ssp,period=period,bio=bio_names,
          positional_RMSE=sqrt(colMeans((as.matrix(old[,-1])-as.matrix(es))^2)),
          population_label_mismatches=sum(old[,1]!=sample_meta$population),
          replacement="reextracted_by_verified_individual_coordinates")
    }
  }
}
sample_result <- do.call(rbind,future_sample);pop_result <- do.call(rbind,future_pop)
stopifnot(nrow(sample_result)==217*24,nrow(pop_result)==19*24)
write_csv(sample_result,file.path(directory,"future_sample_climate.csv"))
write_csv(pop_result,file.path(directory,"future_population_climate.csv"))
write_csv(do.call(rbind,archive_comparison),file.path(directory,"archived_future_positional_comparison.csv"))
save_session(directory)
jsonlite::write_json(list(validated=TRUE,sample_future_rows=nrow(sample_result),population_future_rows=nrow(pop_result),
  extraction="terra simple cell lookup at explicit sample/population coordinates",current_raster_manifest_max_differences=checks,
  current_for_LFMM_RONA="verified individual manifest; raster comparison retained",native_GCM_alignment="native crop cells; archived ensemble on archived grid; sample-coordinate extraction from each",future_sample_identity="sample_id carried from verified metadata; no positional use of unlabeled archived future tables"),
  file.path(directory,"extraction_validation.json"),pretty=TRUE,auto_unbox=TRUE)
