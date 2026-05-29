args <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])))
source(file.path(script_dir, "00_config.R"))
need(c("sf", "dplyr", "ggplot2"))

core <- sf::st_read(file.path(gf_dir, "ammo_core.shp"), quiet = TRUE)
peri <- sf::st_read(file.path(gf_dir, "ammo_peri.shp"), quiet = TRUE)
sf::st_crs(core) <- 4326
sf::st_crs(peri) <- 4326

classify_points <- function(df, x = "x2", y = "y2") {
  pts <- sf::st_as_sf(df, coords = c(x, y), crs = 4326, remove = FALSE)
  in_core <- lengths(sf::st_intersects(pts, core)) > 0
  in_peri <- lengths(sf::st_intersects(pts, peri)) > 0
  ifelse(in_core, "core", ifelse(in_peri, "peripheral", "other"))
}

summarise_sources <- function(label, sink_region, direction) {
  path <- file.path(gf_dir, paste0("ammo_", label, "_", sink_region, "_", direction, "GDM.csv"))
  if (!file.exists(path)) return(NULL)
  df <- read.csv(path, row.names = 1)
  df$source_region <- classify_points(df, "x2", "y2")
  df |>
    count(source_region, name = "n") |>
    mutate(label = label, sink_region = sink_region, direction = direction,
           percent = 100 * n / sum(n))
}

source_summary <- dplyr::bind_rows(lapply(seq_len(nrow(scenarios)), function(i) {
  label <- paste(scenarios$scenario[i], scenarios$period[i], sep = "_")
  dplyr::bind_rows(
    summarise_sources(label, "core", "fwd"),
    summarise_sources(label, "peri", "fwd"),
    summarise_sources(label, "core", "rev"),
    summarise_sources(label, "peri", "rev")
  )
}))

write.csv(source_summary, file.path(out_dir, "gene_flow_source_region_summary.csv"), row.names = FALSE)

p <- ggplot(source_summary, aes(x = label, y = percent, fill = source_region)) +
  geom_col(width = 0.75) +
  facet_grid(direction ~ sink_region) +
  ylab("Predicted source cells (%)") +
  xlab(NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(out_dir, "gene_flow_source_region_summary.pdf"), p, width = 8, height = 5)
