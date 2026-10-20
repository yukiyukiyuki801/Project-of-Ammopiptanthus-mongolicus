run_dir <- Sys.getenv("AMMO_REVISION_RUN")
stopifnot(nzchar(run_dir), dir.exists(run_dir))
output_dir <- file.path(run_dir, "outputs")
root <- normalizePath(file.path(run_dir, "..", ".."))
bundle <- file.path(root, "ammo_analysis_package_20261014")
set.seed(20261013)
options(stringsAsFactors = FALSE, warn = 1)
sample_meta <- read.csv(file.path(output_dir, "inputs", "sample_metadata.csv"))
pop_meta <- read.csv(file.path(output_dir, "inputs", "population_metadata.csv"))
bio_numbers <- c(4,5,6,15,17)
bio_names <- paste0("bio", bio_numbers)
write_csv <- function(x, path) write.csv(x, path, row.names=FALSE, na="")
marker <- function(directory, value) {
  value$validated <- TRUE
  jsonlite::write_json(value, file.path(directory, "validation.json"), pretty=TRUE, auto_unbox=TRUE, na="null")
}
save_session <- function(directory) writeLines(capture.output(sessionInfo()), file.path(directory,"sessionInfo.txt"))
