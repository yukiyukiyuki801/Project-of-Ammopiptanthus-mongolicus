args <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])))
source(file.path(script_dir, "00_config.R"))
need(c("raster", "sf", "gdm", "geosphere", "foreach", "doParallel"))

remove_intercept <- function(mod, pred) {
  adjust <- -log(1 - pred) - mod$intercept
  1 - exp(-adjust)
}

stack_climate <- function(folder) raster::stack(list.files(folder, pattern = "\\.asc$", full.names = TRUE))

fst <- read.csv(file.path(gf_dir, "ammo_outlier_pairwise_fst.csv"))
fst[fst < 0] <- 0
pops <- read.csv(file.path(gf_dir, "ammo_xy_19pop.csv"))
stopifnot(all(fst$pop == pops$code))

range_shp <- raster::shapefile(file.path(gf_dir, "ammo_binary.shp"))
current <- stack_climate(file.path(gf_dir, "chelsa2_crop_selected", "current"))
pred <- data.frame(pop = pops$code, long = pops$long, lat = pops$lat,
                   raster::extract(current, pops[, c("long", "lat")]))
site_pair <- formatsitepair(fst, bioFormat = 3, siteColumn = "pop",
                            XColumn = "long", YColumn = "lat", predData = pred)
mod <- gdm(na.omit(site_pair), geo = TRUE)
saveRDS(mod, file.path(out_dir, "gdm_offset_fit.rds"))

offset_one <- function(future_folder, label, cores = as.integer(Sys.getenv("AMMO_CORES", "4"))) {
  future <- stack_climate(future_folder)
  current_cells <- na.omit(as.data.frame(raster::mask(current, range_shp), xy = TRUE))
  future_cells <- na.omit(as.data.frame(raster::mask(future, range_shp), xy = TRUE))
  current_pop <- split(data.frame(distance = 1, weights = 1, current_cells), seq_len(nrow(current_cells)))
  future_dat <- data.frame(distance = 1, weights = 1, future_cells)

  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  fwd <- foreach::foreach(i = seq_along(current_pop), .combine = rbind, .packages = c("gdm", "geosphere")) %dopar% {
    one <- current_pop[[i]]
    setup <- cbind(one, future_cells)
    colnames(setup) <- c("distance", "weights", "s1.xCoord", "s1.yCoord", paste0("s1.", bio_vars),
                         "s2.xCoord", "s2.yCoord", paste0("s2.", bio_vars))
    dat <- setup[, c("distance", "weights", "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord",
                     paste0("s1.", bio_vars), paste0("s2.", bio_vars))]
    pred_fst <- remove_intercept(mod, predict(mod, dat, time = FALSE))
    cd <- data.frame(dat[, c("s2.xCoord", "s2.yCoord")], predFst = pred_fst)
    coord <- one[, c("x", "y")]
    best <- cd[cd$predFst == min(cd$predFst), ]
    best$dists <- geosphere::distGeo(coord, best[, 1:2])
    best <- best[which.min(best$dists), ]
    local <- cd[cd$s2.xCoord == coord$x & cd$s2.yCoord == coord$y, "predFst"]
    c(x1 = coord$x, y1 = coord$y, local = local, forwardFst = best$predFst,
      predDist = best$dists, bearing = geosphere::bearing(coord, best[, 1:2]),
      x2 = best$s2.xCoord, y2 = best$s2.yCoord)
  }

  rev <- foreach::foreach(i = seq_len(nrow(future_dat)), .combine = rbind, .packages = c("gdm", "geosphere")) %dopar% {
    one <- future_dat[i, ]
    setup <- cbind(one, current_cells)
    colnames(setup) <- c("distance", "weights", "s1.xCoord", "s1.yCoord", paste0("s1.", bio_vars),
                         "s2.xCoord", "s2.yCoord", paste0("s2.", bio_vars))
    dat <- setup[, c("distance", "weights", "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord",
                     paste0("s1.", bio_vars), paste0("s2.", bio_vars))]
    pred_fst <- remove_intercept(mod, predict(mod, dat, time = FALSE))
    cd <- data.frame(dat[, c("s2.xCoord", "s2.yCoord")], predFst = pred_fst)
    coord <- one[, c("x", "y")]
    best <- cd[cd$predFst == min(cd$predFst), ]
    best$dists <- geosphere::distGeo(coord, best[, 1:2])
    best <- best[which.min(best$dists), ]
    local <- cd[cd$s2.xCoord == coord$x & cd$s2.yCoord == coord$y, "predFst"]
    c(x1 = coord$x, y1 = coord$y, local = local, reverseFst = best$predFst,
      predDist = best$dists, bearing = geosphere::bearing(coord, best[, 1:2]),
      x2 = best$s2.xCoord, y2 = best$s2.yCoord)
  }

  write.csv(fwd, file.path(out_dir, paste0("ammo_", label, "_fwdGDM.csv")), row.names = FALSE)
  write.csv(rev, file.path(out_dir, paste0("ammo_", label, "_revGDM.csv")), row.names = FALSE)
}

for (i in seq_len(nrow(scenarios))) {
  label <- paste(scenarios$scenario[i], scenarios$period[i], sep = "_")
  offset_one(file.path(gf_dir, "chelsa2_crop_selected", label), label)
}
