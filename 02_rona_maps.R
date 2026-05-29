args <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])))
source(file.path(script_dir, "00_config.R"))
need(c("raster", "sf", "ggplot2"))

range_shp <- sf::st_read(file.path(gf_dir, "ammo_binary.shp"), quiet = TRUE)
range_sp <- as(range_shp, "Spatial")

plot_rona <- function(variable, scenario, period) {
  period_long <- if (period == "2040") "2041_2070" else "2071_2100"
  points_path <- file.path(rona_dir, paste0(scenario, "_", period_long, "_rona_xy.csv"))
  diff_path <- file.path(rona_dir, paste0(variable, "_", scenario, "_", period_long, "_diff.asc"))
  if (!file.exists(points_path) || !file.exists(diff_path)) return(NULL)

  rona <- read.csv(points_path, row.names = 1)
  r <- raster::mask(raster::crop(raster::raster(diff_path), range_sp), range_sp)
  rdf <- as.data.frame(as(r, "SpatialPixelsDataFrame"))
  value_col <- names(rdf)[1]

  p <- ggplot() +
    geom_raster(data = rdf, aes(x = x, y = y, fill = .data[[value_col]])) +
    geom_sf(data = range_shp, fill = NA, color = "black", linewidth = 0.3) +
    geom_point(data = rona, aes(x = x, y = y, size = .data[[variable]]),
               shape = 21, fill = "red", color = "black", alpha = 0.8) +
    scale_fill_gradient2(low = "lightblue", mid = "yellow", high = "tomato",
                         midpoint = mean(range(rdf[[value_col]], na.rm = TRUE)), name = "Climate change") +
    scale_size_continuous(name = "RONA") +
    coord_sf(expand = FALSE) +
    theme_bw() +
    theme(panel.grid = element_blank())

  ggsave(file.path(out_dir, paste0("rona_", variable, "_", scenario, "_", period, ".pdf")), p, width = 7, height = 5)
}

for (v in rona_vars) {
  for (i in seq_len(nrow(scenarios))) plot_rona(v, scenarios$scenario[i], scenarios$period[i])
}
