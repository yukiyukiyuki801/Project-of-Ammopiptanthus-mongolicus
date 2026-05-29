args <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])))
source(file.path(script_dir, "00_config.R"))
need(c("sdm", "raster", "sf"))

occ_shp <- Sys.getenv("AMMO_OCCURRENCE_SHP", file.path(ammo_dir, "ammo_thin_10km.shp"))
if (!file.exists(occ_shp)) stop("Set AMMO_OCCURRENCE_SHP to the thinned occurrence shapefile.", call. = FALSE)

read_stack <- function(folder) {
  files <- list.files(folder, pattern = "\\.asc$", full.names = TRUE)
  files <- files[match(bio_vars, tools::file_path_sans_ext(basename(files)))]
  if (anyNA(files)) stop("Missing one or more SDM raster predictors in ", folder, call. = FALSE)
  raster::stack(files)
}

species <- raster::shapefile(occ_shp)
current <- read_stack(file.path(gf_dir, "chelsa2_crop_selected", "current"))
d <- sdmData(Species ~ ., train = species, predictors = current,
             bg = list(n = 500, method = "gRandom", remove = TRUE))

fit <- sdm(Species ~ ., data = d,
           methods = c("glm", "maxent", "svm", "gam", "mda", "mlp"),
           replication = "cv", cv.folds = 10)
saveRDS(fit, file.path(out_dir, "sdm_ensemble_fit.rds"))

predict_ensemble <- function(folder, label) {
  stk <- read_stack(folder)
  pred <- ensemble(fit, newdata = stk, setting = list(method = "weighted", stat = "TSS", opt = 2))
  raster::writeRaster(pred, file.path(out_dir, paste0("sdm_", label, ".asc")), format = "ascii", overwrite = TRUE)
}

predict_ensemble(file.path(gf_dir, "chelsa2_crop_selected", "current"), "current")
for (i in seq_len(nrow(scenarios))) {
  label <- paste(scenarios$scenario[i], scenarios$period[i], sep = "_")
  predict_ensemble(file.path(gf_dir, "chelsa2_crop_selected", label), label)
}
