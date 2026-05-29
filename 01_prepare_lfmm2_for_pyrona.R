args <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])))
source(file.path(script_dir, "00_config.R"))
need(c("LEA", "nFactors"))

original_vcf <- Sys.getenv("AMMO_ADAPTIVE_VCF", file.path(ammo_dir, "ammo_outliers_3overlap.recode.vcf"))
target_lfmm <- file.path(rona_dir, "ammo_outliers_3overlap.recode.lfmm")
current_env <- file.path(rona_dir, "ammo_current.txt")
association_table <- file.path(rona_dir, "associations.csv")
k_value <- Sys.getenv("LFMM_K", "4")

read_required(original_vcf, function(path, ...) path)
read_required(current_env, function(path, ...) path)

if (k_value == "estimate") {
  k_value <- "estimate"
} else {
  k_value <- as.integer(k_value)
}

env_file_compat <- function(env_file) {
  if (!grepl("\\.env$", env_file)) {
    env_vars <- read.table(env_file, sep = "\t", header = FALSE)[, -1, drop = FALSE]
    new_filename <- sub("\\.[^.]*$", ".env", env_file)
    write.table(env_vars, file = new_filename, sep = "\t",
                col.names = FALSE, row.names = FALSE, quote = FALSE)
    env_file <- new_filename
  }
  env_file
}

choose_k <- function(genetic_data, k) {
  pc <- prcomp(genetic_data)
  pdf(file.path(out_dir, "lfmm_pca_k_check.pdf"), width = 5, height = 4)
  plot(pc$sdev[seq_len(min(20, length(pc$sdev)))]^2,
       xlab = "PC", ylab = "Variance explained", pch = 19)
  if (identical(k, "estimate")) {
    k <- as.integer(as.character(nFactors::nCng(pc$sdev^2, cor = FALSE, details = FALSE)$nFactors))
  }
  points(k, pc$sdev[k]^2, type = "h", lwd = 3, col = "blue")
  dev.off()
  k
}

message("Converting VCF to LFMM: ", original_vcf)
LEA::vcf2lfmm(original_vcf, target_lfmm, force = TRUE)

genetic_data <- LEA::read.lfmm(target_lfmm)
k_value <- choose_k(genetic_data, k_value)
message("Using LFMM K = ", k_value)

message("Imputing missing genotypes with snmf")
project <- LEA::snmf(target_lfmm, K = k_value, entropy = TRUE,
                     repetitions = 10, project = "new")
best <- which.min(LEA::cross.entropy(project, K = k_value))
LEA::impute(project, target_lfmm, method = "mode", K = k_value, run = best)

imputed_lfmm <- paste0(target_lfmm, "_imputed.lfmm")
env_file <- env_file_compat(current_env)

message("Running LFMM2 association model")
imputed_data <- LEA::read.lfmm(imputed_lfmm)
mod <- LEA::lfmm2(input = imputed_data, env = env_file, K = k_value)
pv <- LEA::lfmm2.test(object = mod, input = imputed_data, env = env_file, linear = TRUE)
pvalues <- t(pv$pvalues)

write.table(pvalues, file = association_table, col.names = FALSE,
            row.names = FALSE, sep = ",", quote = FALSE)

message("Wrote imputed LFMM: ", imputed_lfmm)
message("Wrote association table: ", association_table)
