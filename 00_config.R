`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || is.na(x)) y else x

if (!exists("script_dir", inherits = FALSE)) {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1] %||% NA_character_)
  script_dir <- if (!is.na(file_arg)) dirname(normalizePath(file_arg)) else getwd()
}

ammo_dir <- Sys.getenv("AMMO_DIR", "/home/liao/DataVol1/personal2/sun/ammo")
gf_dir <- file.path(ammo_dir, "gf")
rona_dir <- file.path(ammo_dir, "rona")
out_dir <- Sys.getenv("AMMO_ANALYSIS_OUT", file.path(script_dir, "outputs"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

scenarios <- expand.grid(
  scenario = c("ssp126", "ssp585"),
  period = c("2040", "2070"),
  stringsAsFactors = FALSE
)

bio_vars <- c("b4", "b5", "b6", "b15", "b17")
rona_vars <- c("Bio4", "Bio5", "Bio6", "Bio15", "Bio17")

need <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install required R packages: ", paste(missing, collapse = ", "), call. = FALSE)
  invisible(lapply(pkgs, library, character.only = TRUE))
}

read_required <- function(path, fun = read.csv, ...) {
  if (!file.exists(path)) stop("Missing required file: ", path, call. = FALSE)
  fun(path, ...)
}

message("Using ammo_dir: ", ammo_dir)
message("Writing outputs to: ", out_dir)
