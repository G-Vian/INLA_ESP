# ==============================================================================
# 04 - RESULTS OUTPUTS: model-selection metrics, DLNM curves, forest plots
# ------------------------------------------------------------------------------
# Writes to disk the results reported in the paper that the fitting modules
# compute but do not otherwise save:
#   (1) the stepwise model-selection table (WAIC, DIC, log-score, pseudo-R2)
#   (2) the DLNM overall exposure-response curves for each climate variable
#   (3) forest plots and a table of the fixed-effect relative risks
#
# Consumes objects produced by module 03 (mod_final, tabela_res, f_vencedora)
# and module 02 (cb_objects_global). Everything is written to figures/.
# ==============================================================================

# ---- required packages -------------------------------------------------------
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
})

message("\n=== MODULE 04: results outputs ===")

DIR_OUT_04 <- get0("DIR_OUTPUT", ifnotfound = "figures")
dir.create(DIR_OUT_04, recursive = TRUE, showWarnings = FALSE)

`%||%` <- function(a, b) if (is.null(a)) b else a

# ==============================================================================
# 1. MODEL-SELECTION METRICS (stepwise): WAIC / DIC / log-score / pseudo-R2
# ==============================================================================
if (exists("tabela_res") && is.data.frame(tabela_res) && nrow(tabela_res) > 0) {
  write.csv(tabela_res, file.path(DIR_OUT_04, "TABLE_model_selection_metrics.csv"),
            row.names = FALSE)
  message("[04] TABLE_model_selection_metrics.csv written (",
          nrow(tabela_res), " model steps).")
  
  # tidy plain-text version
  sink(file.path(DIR_OUT_04, "TABLE_model_selection_metrics.txt"))
  cat("Stepwise model-selection metrics\n")
  cat("Lower WAIC / DIC / log-score = better; higher pseudo-R2 = better.\n\n")
  print(knitr::kable(tabela_res, format = "markdown", digits = 3))
  sink()
} else {
  message("[04] tabela_res not found — model-selection table skipped.")
}

# Also report the final model's headline metrics on their own.
if (exists("mod_final") && !is.null(mod_final)) {
  ls_final <- {
    cpo <- mod_final$cpo$cpo
    if (!is.null(cpo)) sum(-log(cpo[cpo > 0 & !is.na(cpo)])) else NA_real_
  }
  final_metrics <- data.frame(
    DIC      = tryCatch(mod_final$dic$dic,   error = function(e) NA_real_),
    WAIC     = tryCatch(mod_final$waic$waic, error = function(e) NA_real_),
    LogScore = ls_final,
    MargLik  = tryCatch(mod_final$mlik[1],   error = function(e) NA_real_)
  )
  write.csv(final_metrics, file.path(DIR_OUT_04, "TABLE_final_model_metrics.csv"),
            row.names = FALSE)
  message("[04] TABLE_final_model_metrics.csv written.")
}

# ==============================================================================
# 2. DLNM OVERALL EXPOSURE-RESPONSE CURVES (one per climate variable)
# ==============================================================================
if (exists("mod_final") && !is.null(mod_final) &&
    exists("cb_objects_global") && length(cb_objects_global) > 0 &&
    requireNamespace("dlnm", quietly = TRUE)) {
  
  fx <- mod_final$summary.fixed
  q_hi_dlnm <- get0("CONFIG_DLNM_PLOT_QUANTIL", ifnotfound = 0.99)
  var_labels <- c(temp_min = "Minimum temperature",
                  precip_med = "Mean precipitation",
                  temp_max = "Maximum temperature")
  
  dlnm_all <- list()
  for (v in names(cb_objects_global)) {
    idx <- grep(paste0("^basis_", v, "\\."), rownames(fx))
    if (length(idx) == 0) next
    betas <- fx$mean[idx]
    vc    <- diag(fx$sd[idx]^2)                 # diagonal (conservative)
    x     <- dados_final_modelo[[v]]
    if (is.null(x)) next
    hi   <- as.numeric(stats::quantile(x, q_hi_dlnm, na.rm = TRUE))
    grid <- seq(min(x, na.rm = TRUE), hi, length.out = 60)
    cen  <- stats::median(x, na.rm = TRUE)
    pr <- tryCatch(
      dlnm::crosspred(cb_objects_global[[v]], coef = betas, vcov = vc,
                      model.link = "log", at = grid, cen = cen),
      error = function(e) NULL)
    if (is.null(pr)) next
    
    dfc <- data.frame(
      Variable  = v,
      Exposure  = as.numeric(names(pr$allRRfit)),
      RR        = as.numeric(pr$allRRfit),
      RR_low    = as.numeric(pr$allRRlow),
      RR_high   = as.numeric(pr$allRRhigh))
    dlnm_all[[v]] <- dfc
    
    lab <- var_labels[v] %||% v
    pD <- ggplot(dfc, aes(Exposure, RR)) +
      geom_hline(yintercept = 1, linetype = "dashed", colour = "grey45") +
      geom_ribbon(aes(ymin = RR_low, ymax = RR_high), alpha = 0.18,
                  fill = "#2166AC") +
      geom_line(colour = "#2166AC", linewidth = 1) +
      labs(title = paste0("DLNM overall exposure-response: ", lab),
           subtitle = "Cumulative relative risk over the lag window (ref. = median)",
           x = lab, y = "Relative Risk (RR)") +
      theme_bw(base_size = 12) +
      theme(panel.grid.minor = element_blank(),
            plot.title = element_text(face = "bold"),
            axis.title = element_text(face = "bold"))
    ggsave(file.path(DIR_OUT_04, paste0("FIG_DLNM_", v, ".png")),
           pD, width = 8, height = 5, dpi = 600, bg = "white")
    
    # 2D exposure-lag-response heatmap (matRRfit), restricted to P01-P99 to
    # avoid extrapolation into the sparse tails of the monthly distribution.
    pr_hm <- tryCatch(
      dlnm::crosspred(cb_objects_global[[v]], coef = betas, vcov = vc,
                      model.link = "log", at = grid, bylag = 0.2, cen = cen),
      error = function(e) NULL)
    if (!is.null(pr_hm) && !is.null(pr_hm$matRRfit)) {
      rr_mat <- pr_hm$matRRfit
      hm <- data.frame(
        Exposure = as.numeric(rownames(rr_mat)[row(rr_mat)]),
        Lag      = as.numeric(sub("lag", "", colnames(rr_mat)[col(rr_mat)])),
        RR       = as.numeric(rr_mat))
      rng <- range(rr_mat, na.rm = TRUE)
      lo <- max(0.3, rng[1]); hi2 <- min(4, rng[2])
      wp <- (1 - lo) / (hi2 - lo)                     # white anchored at RR = 1
      pH <- ggplot(hm, aes(Lag, Exposure, fill = RR)) +
        geom_raster(interpolate = TRUE) +
        geom_hline(yintercept = cen, linetype = "dashed",
                   colour = "black", linewidth = 0.4) +
        scale_fill_gradientn(
          colours = c("#2166AC","#4393C3","#92C5DE","#D1E5F0","white",
                      "#FDDBC7","#F4A582","#D6604D","#B2182B"),
          values  = c(seq(0, wp, length.out = 4), wp,
                      seq(wp, 1, length.out = 5)[-1]),
          limits  = c(lo, hi2), oob = scales::squish,
          name = "RR") +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        labs(title = paste0("Exposure-lag-response surface: ", lab),
             subtitle = "Restricted to the 1st-99th percentiles (ref. = median)",
             x = "Lag (months)", y = lab) +
        theme_bw(base_size = 12) +
        theme(panel.grid = element_blank(),
              plot.title = element_text(face = "bold"),
              axis.title = element_text(face = "bold"),
              legend.key.height = unit(1.4, "cm"))
      ggsave(file.path(DIR_OUT_04, paste0("FIG_DLNM_heatmap_", v, ".png")),
             pH, width = 8, height = 6, dpi = 600, bg = "white")
    }
  }
  
  if (length(dlnm_all) > 0) {
    write.csv(bind_rows(dlnm_all),
              file.path(DIR_OUT_04, "TABLE_DLNM_curves.csv"), row.names = FALSE)
    message("[04] DLNM curves written for: ",
            paste(names(dlnm_all), collapse = ", "), ".")
  } else {
    message("[04] no DLNM curve could be built (check basis_ terms).")
  }
} else {
  message("[04] cb_objects_global / dlnm unavailable — DLNM curves skipped.")
}

# ==============================================================================
# 3. FOREST PLOT OF FIXED-EFFECT RELATIVE RISKS
# ==============================================================================
if (exists("mod_final") && !is.null(mod_final)) {
  fx <- mod_final$summary.fixed
  # keep interpretable fixed effects: drop the DLNM basis columns and intercept
  keep <- !grepl("^basis_", rownames(fx)) &
    !grepl("Intercept|intercept|\\(Intercept\\)", rownames(fx))
  rr <- data.frame(
    Term = rownames(fx)[keep],
    RR   = exp(fx$mean[keep]),
    LI   = exp(fx$`0.025quant`[keep]),
    LS   = exp(fx$`0.975quant`[keep]),
    stringsAsFactors = FALSE)
  
  if (nrow(rr) > 0) {
    rr$Credible <- (rr$LI > 1) | (rr$LS < 1)
    write.csv(rr, file.path(DIR_OUT_04, "TABLE_fixed_effects_RR.csv"),
              row.names = FALSE)
    
    rr_plot <- rr %>% arrange(RR) %>% mutate(Term = factor(Term, levels = Term))
    pF <- ggplot(rr_plot, aes(x = RR, y = Term, colour = Credible)) +
      geom_vline(xintercept = 1, linetype = "dashed", colour = "grey45") +
      geom_point(size = 2.6) +
      geom_errorbar(aes(xmin = LI, xmax = LS), orientation = "y",
                    width = 0.25, linewidth = 0.7) +
      scale_colour_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "grey55"),
                          labels = c(`TRUE` = "credible", `FALSE` = "not credible"),
                          name = NULL) +
      labs(title = "Fixed-effect relative risks (final model)",
           subtitle = "Posterior mean and 95% credible interval; dashed line at RR = 1",
           x = "Relative Risk (RR)", y = NULL) +
      theme_bw(base_size = 12) +
      theme(panel.grid.minor = element_blank(),
            plot.title = element_text(face = "bold"),
            axis.title = element_text(face = "bold"),
            legend.position = "bottom")
    ggsave(file.path(DIR_OUT_04, "FIG_forest_fixed_effects.png"),
           pF, width = 8, height = 5.5, dpi = 600, bg = "white")
    message("[04] forest plot + TABLE_fixed_effects_RR.csv written (",
            nrow(rr), " terms).")
  } else {
    message("[04] no interpretable fixed effect found — forest plot skipped.")
  }
}

message("=== MODULE 04 done. Outputs in ", DIR_OUT_04, " ===\n")