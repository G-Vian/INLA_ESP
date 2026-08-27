# ==============================================================================
# 05 - EXTERNAL-MODULE DISPATCHER
# ------------------------------------------------------------------------------
# Runs the external analysis modules that ship with this repository. In the
# published pipeline this dispatcher launched many optional sub-analyses; here
# it is reduced to the two that reproduce the predictive validation reported in
# the paper:
#
#   External 1 - leakage-free rolling-origin cross-validation
#   External 2 - out-of-sample posterior predictive intervals (coverage, MAE)
#
# Both are sourced by R/main.R directly, so by default this module only reports
# what has run. If you prefer to launch them from here instead of from main.R,
# set the switches below and remove the corresponding lines from main.R.
#
# NOTE: the DLNM exposure-lag-response surfaces and the spatial relative-risk
# maps are produced inside module 03 / the model diagnostics. If you maintain a
# separate figure module for them, add it to R/ and source it here.
# ==============================================================================

DIR_EXTERNAL <- if (exists("DIR_MOD_INT")) DIR_MOD_INT else
  file.path(if (exists("DIR_BASE_CODIGO")) DIR_BASE_CODIGO else getwd(), "R")

run_external <- function(file) {
  path <- file.path(DIR_EXTERNAL, file)
  if (!file.exists(path)) {
    message("[05] external module not found, skipping: ", file)
    return(invisible(FALSE))
  }
  message("[05] running external module: ", file)
  source(path, encoding = "UTF-8", local = FALSE)
  invisible(TRUE)
}

# By default these are launched from main.R (so they are not re-run here).
# Flip to TRUE to have this dispatcher launch them instead.
RUN_EXTERNAL_FROM_DISPATCHER <- get0("RUN_EXTERNAL_FROM_DISPATCHER",
                                     ifnotfound = FALSE)

if (isTRUE(RUN_EXTERNAL_FROM_DISPATCHER)) {
  run_external("External_1_rolling_origin.R")
  run_external("External_2_predictive_intervals.R")
} else {
  message("[05] external modules are launched from main.R; dispatcher idle.")
}
