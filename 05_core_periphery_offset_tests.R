args <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])))
source(file.path(script_dir, "00_config.R"))
need(c("raster", "sf", "dplyr"))

sample_offsets <- function(label, region, n = 10000) {
  shp <- sf::st_read(file.path(gf_dir, paste0("ammo_", region, ".shp")), quiet = TRUE)
  shp_sp <- as(shp, "Spatial")
  fwd <- read.csv(file.path(gf_dir, paste0("ammo_", label, "_fwdGDM.csv")))
  rev <- read.csv(file.path(gf_dir, paste0("ammo_", label, "_revGDM.csv")))

  rasters <- list(
    forward = raster::rasterFromXYZ(fwd[, c("x1", "y1", "forwardFst")]),
    reverse = raster::rasterFromXYZ(rev[, c("x1", "y1", "reverseFst")]),
    local = raster::rasterFromXYZ(rev[, c("x1", "y1", "local")])
  )

  dplyr::bind_rows(lapply(names(rasters), function(metric) {
    r <- raster::mask(raster::crop(rasters[[metric]], shp_sp), shp_sp)
    vals <- na.omit(raster::values(r))
    data.frame(label = label, region = region, metric = metric,
               value = sample(vals, n, replace = length(vals) < n))
  }))
}

all_samples <- dplyr::bind_rows(lapply(seq_len(nrow(scenarios)), function(i) {
  label <- paste(scenarios$scenario[i], scenarios$period[i], sep = "_")
  dplyr::bind_rows(sample_offsets(label, "core"), sample_offsets(label, "peri"))
}))

write.csv(all_samples, file.path(out_dir, "core_periphery_offset_samples.csv"), row.names = FALSE)

summary_tab <- all_samples |>
  group_by(label, metric, region) |>
  summarise(mean = mean(value), sd = sd(value), median = median(value), .groups = "drop")
write.csv(summary_tab, file.path(out_dir, "core_periphery_offset_summary.csv"), row.names = FALSE)

tests <- all_samples |>
  group_by(label, metric) |>
  summarise(p_value = wilcox.test(value[region == "core"], value[region == "peri"])$p.value,
            core_mean = mean(value[region == "core"]),
            peri_mean = mean(value[region == "peri"]),
            .groups = "drop")
write.csv(tests, file.path(out_dir, "core_periphery_offset_wilcox_tests.csv"), row.names = FALSE)
