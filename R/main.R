# ==============================================================================
# ORCHESTRATOR - Bayesian spatio-temporal dengue model for São Paulo State
# ------------------------------------------------------------------------------
# Runs the core pipeline modules in order, in the global environment.
#
# Reproduces: the main INLA model, the DLNM exposure-lag-response surfaces,
# and the leakage-free rolling-origin predictive validation.
#
# The spatial structure is the Besag (ICAR) effect built from the municipal
# shapefile (SWITCH_GRAFO_VERSAO = 1 in R/00_config.R).
#
# USAGE
#   1. Place the processed cache file in  data/  (see README, "Data").
#   2. From the repository root, run:   Rscript R/main.R
#   All paths are relative to the repository root; nothing is hard-coded.
# ==============================================================================

rm(list = ls()); gc()

# Repository root = parent of this R/ directory. Works with Rscript and source().
if (requireNamespace("here", quietly = TRUE)) {
  DIR_BASE_CODIGO <- here::here()
} else {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
  DIR_BASE_CODIGO <- if (length(file_arg))
    normalizePath(file.path(dirname(file_arg), "..")) else normalizePath("..")
}
DIR_MOD_INT <- file.path(DIR_BASE_CODIGO, "R")

run_module <- function(file) {
  path <- file.path(DIR_MOD_INT, file)
  message("\n==================================================================")
  message("\u25B6  RUNNING MODULE: ", file)
  message("==================================================================")
  tryCatch(
    source(path, encoding = "UTF-8", local = FALSE),
    error = function(e) {
      stop("\n\u274C FAILED in module '", file, "':\n   ", conditionMessage(e),
           "\n   (Subsequent modules were NOT executed.)", call. = FALSE)
    }
  )
}

run_module("00_config.R")
run_module("R00b_figures.R")
run_module("01_data.R")
run_module("02_features_dlnm_lags.R")
run_module("03_model.R")
run_module("04_outputs.R")               # metrics table, DLNM curves, forest plots
run_module("External_1_rolling_origin.R")        # leakage-free rolling-origin validation
run_module("External_2_predictive_intervals.R")  # out-of-sample predictive intervals + MAE
run_module("05_diagnostics.R")           # external-module dispatcher (idle by default)

message("\n\u2705 Pipeline finished (SWITCH_GRAFO_VERSAO = ", SWITCH_GRAFO_VERSAO, ").")