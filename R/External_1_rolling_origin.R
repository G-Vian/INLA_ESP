# ==============================================================================
# EXTERNAL 1 - ROLLING-ORIGIN CROSS-VALIDATION (LEAKAGE-FREE)
# ------------------------------------------------------------------------------
# For each target month t, refits the model on months <= t-1 and predicts month
# t. The autoregressive predictors for t are t-1 and t-2, both inside the
# training window, so no information outside the training set enters any
# prediction. Checkpoints one RDS per origin so an interrupted run resumes.
# Produces per-origin predictions with SD_PRED, consumed by External 2.
# Original working comments are kept in Portuguese below.
# ==============================================================================

# ==============================================================================
# MONTHLY ROLLING-ORIGIN CROSS-VALIDATION (LEAKAGE-FREE)
# ------------------------------------------------------------------------------
# OBJETIVO
#   Predictive validation with a monthly sliding origin. For each target month t:
#     - train EXCLUSIVELY on ID_TEMPO <= t-1
#     - predict month t
#     - advance the origin by one month and refit
#
#   The autoregressive predictors for month t are months t-1 and t-2:
#     when predicting month t, lag-1 is month t-1 and lag-2 is month t-2,
#     both of which belong ENTIRELY to the training set (ID_TEMPO <= t-1).
#     No observation outside the training window is used at any point,
#     so the design is leakage-free by construction.
#
#
#
# ROBUSTEZ
#   - One checkpoint per origin (1 RDS per month). If the job stops, resubmit:
#     completed origins are skipped automatically.
#   - Warm start: the hyperparameter mode of the previous origin is reused
#     (control.mode), substantially reducing the time per fit.
#   - Checkpoints store only small objects (predictions, RR, theta),
#     nunca o objeto INLA completo.
#
# PIPELINE DEPENDENCIES
#   dados_final_modelo  (with ID_TEMPO, DATA, ANO, MES, CASOS, POP_TOTAL,
#                        COD_GEO, DRS_ID)
#   f_vencedora         (final model formula)
#   ctrl_inla_params, ctrl_family_params, ctrl_fixed_params, CONFIG_FAMILIA
#   DIR_OUTPUT
#
# OUTPUTS  (in DIR_OUTPUT/teste_rolling_mensal/)
#   TABELA_ROLLING_MENSAL_METRICAS.txt / .csv
#   TABELA_ROLLING_MENSAL_POR_ANO.csv
#   TABELA_ROLLING_MENSAL_POR_DRS.csv
#   TABELA_ROLLING_MENSAL_POR_PORTE.csv
#   TABELA_ROLLING_MENSAL_RR_EVOLUCAO.csv
#   PREDICOES_ROLLING_MENSAL.csv
#   FIG_RollingOrigin_Painel_Publicacao.png   (main figure, 4 panels)
#   FIG_RollingOrigin_SerieEstado.png
#   FIG_RollingOrigin_RR_Estabilidade.png
# ==============================================================================

message("\n======================================================================")
message("  MÓDULO 24: ROLLING-ORIGIN MENSAL (LEAKAGE-FREE)")
message("======================================================================\n")

# Operador auxiliar (definido no topo porque e usado dentro do bloco abaixo)
if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a

# ------------------------------------------------------------------------------
# 0. GUARDA DE ENTRADA
# ------------------------------------------------------------------------------
.rolling_pode_rodar <- exists("dados_final_modelo") &&
  exists("f_vencedora") &&
  isTRUE(get0("SWITCH_TESTE_ROLLING_MENSAL", ifnotfound = FALSE))

if (!.rolling_pode_rodar) {
  
  message("Módulo 24 ignorado (SWITCH_TESTE_ROLLING_MENSAL desligado ou ",
          "objetos ausentes).")
  
} else {
  
  suppressPackageStartupMessages({
    library(dplyr); library(tidyr); library(ggplot2); library(tibble)
  })
  
  # ----------------------------------------------------------------------------
  # 1. CONFIGURATION (all with get0 defaults — nothing required in the config)
  # ----------------------------------------------------------------------------
  CFG <- list(
    # Number of initial months reserved for training before the first prediction.
    # 36 = three years.
    burnin        = get0("CONFIG_ROLLING_BURNIN_MESES",  ifnotfound = 36L),
    
    # Restrict origins to a range of years (NULL = whole series).
    # e.g. c(2020, 2024) to run only the five most recent years.
    anos_teste    = get0("CONFIG_ROLLING_ANOS_TESTE",    ifnotfound = NULL),
    
    # Warm start from the previous origin's mode.
    warmstart     = get0("CONFIG_ROLLING_WARMSTART",     ifnotfound = TRUE),
    
    # Diagonal added for numerical stability.
    diagonal      = get0("CONFIG_ROLLING_DIAGONAL",      ifnotfound = 0.1),
    
    # Quantile used to define an "epidemic month" for AUC / hit rate.
    q_surto       = get0("CONFIG_ROLLING_QUANTIL_SURTO", ifnotfound = 0.75),
    
    # If TRUE, overwrite the global df_cv/metrics_cv objects.
    # Default FALSE to avoid overwriting them silently.
    ponte_df_cv   = get0("CONFIG_ROLLING_PONTE_DF_CV",   ifnotfound = FALSE),
    
    # Reprocessar do zero, ignorando checkpoints existentes.
    forcar        = get0("CONFIG_ROLLING_FORCAR_RECALC", ifnotfound = FALSE)
  )
  
  DIR_ROLL <- file.path(DIR_OUTPUT, "teste_rolling_mensal")
  DIR_CKPT <- file.path(DIR_ROLL, "checkpoints")
  dir.create(DIR_CKPT, recursive = TRUE, showWarnings = FALSE)
  
  if (isTRUE(CFG$forcar)) {
    unlink(list.files(DIR_CKPT, full.names = TRUE, pattern = "\\.rds$"))
    message("[cfg] CONFIG_ROLLING_FORCAR_RECALC = TRUE -> checkpoints apagados.")
  }
  
  # ----------------------------------------------------------------------------
  # 2. ORIGIN DEFINITION
  # ----------------------------------------------------------------------------
  if (!"ID_TEMPO" %in% names(dados_final_modelo))
    stop("Módulo 24: coluna ID_TEMPO ausente em dados_final_modelo.")
  
  mapa_tempo <- dados_final_modelo %>%
    distinct(ID_TEMPO, DATA, ANO, MES) %>%
    arrange(ID_TEMPO)
  
  id_min <- min(mapa_tempo$ID_TEMPO)
  id_max <- max(mapa_tempo$ID_TEMPO)
  
  ids_alvo <- mapa_tempo$ID_TEMPO[mapa_tempo$ID_TEMPO >= (id_min + CFG$burnin)]
  
  if (!is.null(CFG$anos_teste)) {
    anos_ok  <- mapa_tempo$ANO[match(ids_alvo, mapa_tempo$ID_TEMPO)]
    ids_alvo <- ids_alvo[anos_ok >= min(CFG$anos_teste) &
                           anos_ok <= max(CFG$anos_teste)]
  }
  
  if (length(ids_alvo) == 0)
    stop("Módulo 24: nenhuma origem elegível. Revise burn-in / anos_teste.")
  
  message(sprintf("[plano] %d origens mensais | de %s ate %s | burn-in de %d meses",
                  length(ids_alvo),
                  format(mapa_tempo$DATA[match(min(ids_alvo), mapa_tempo$ID_TEMPO)], "%Y-%m"),
                  format(mapa_tempo$DATA[match(max(ids_alvo), mapa_tempo$ID_TEMPO)], "%Y-%m"),
                  CFG$burnin))
  message(sprintf("[plano] warm start: %s | checkpoints em: %s",
                  ifelse(CFG$warmstart, "ATIVO", "desativado"), DIR_CKPT))
  
  ckpt_path <- function(id) file.path(DIR_CKPT, sprintf("origem_%04d.rds", id))
  
  ja_feitos <- ids_alvo[file.exists(vapply(ids_alvo, ckpt_path, character(1)))]
  se_fazer  <- setdiff(ids_alvo, ja_feitos)
  
  if (length(ja_feitos) > 0)
    message(sprintf("[ckpt] %d origens ja concluidas serao puladas | %d restantes.",
                    length(ja_feitos), length(se_fazer)))
  
  # ----------------------------------------------------------------------------
  # 3. RECUPERA theta MAIS RECENTE PARA WARM START (sobrevive a queda do job)
  # ----------------------------------------------------------------------------
  theta_prev <- NULL
  if (isTRUE(CFG$warmstart) && length(ja_feitos) > 0) {
    ult <- max(ja_feitos)
    theta_prev <- tryCatch(readRDS(ckpt_path(ult))$theta, error = function(e) NULL)
    if (!is.null(theta_prev))
      message("[ckpt] warm start recuperado da origem ", ult, ".")
  }
  
  # ----------------------------------------------------------------------------
  # 4. MAIN LOOP — ONE ORIGIN PER ITERATION
  # ----------------------------------------------------------------------------
  t_inicio_geral <- Sys.time()
  
  for (k in seq_along(se_fazer)) {
    
    id_alvo <- se_fazer[k]
    linha_t <- mapa_tempo[mapa_tempo$ID_TEMPO == id_alvo, ]
    
    message(sprintf("\n   [%d/%d] origem -> prevendo %s (ID_TEMPO=%d) | treino: ID_TEMPO <= %d",
                    k, length(se_fazer),
                    format(linha_t$DATA, "%Y-%m"), id_alvo, id_alvo - 1L))
    
    t0 <- Sys.time()
    
    # --- Recorte temporal: nada posterior ao mes-alvo entra na base -----------
    dt <- dados_final_modelo %>% filter(ID_TEMPO <= id_alvo)
    idx_alvo <- which(dt$ID_TEMPO == id_alvo)
    
    # A resposta do mes-alvo vira NA: essas linhas NAO entram na verossimilhanca.
    # The lags of these rows point to t-1 and t-2, which ARE in the training set.
    dt$CASOS[idx_alvo] <- NA
    
    ctrl <- ctrl_inla_params
    ctrl$diagonal <- CFG$diagonal
    
    ctrl_mode <- NULL
    if (isTRUE(CFG$warmstart) && !is.null(theta_prev))
      ctrl_mode <- list(theta = theta_prev, restart = TRUE)
    
    m <- tryCatch({
      INLA::inla(
        f_vencedora,
        family            = CONFIG_FAMILIA,
        data              = dt,
        control.family    = ctrl_family_params,
        control.inla      = ctrl,
        control.fixed     = ctrl_fixed_params,
        control.mode      = ctrl_mode,
        control.predictor = list(compute = TRUE, link = 1),
        control.compute   = list(dic = FALSE, waic = FALSE, cpo = FALSE)
      )
    }, error = function(e) {
      message("      [ERRO] INLA falhou nesta origem: ", e$message)
      NULL
    })
    
    # If the warm start hurt, try once without it before giving up.
    if (is.null(m) && !is.null(ctrl_mode)) {
      message("      [retry] refazendo sem warm start...")
      m <- tryCatch({
        INLA::inla(
          f_vencedora, family = CONFIG_FAMILIA, data = dt,
          control.family = ctrl_family_params, control.inla = ctrl,
          control.fixed = ctrl_fixed_params,
          control.predictor = list(compute = TRUE, link = 1),
          control.compute = list(dic = FALSE, waic = FALSE, cpo = FALSE)
        )
      }, error = function(e) NULL)
    }
    
    if (is.null(m)) {
      saveRDS(list(id_tempo = id_alvo, ok = FALSE, theta = theta_prev),
              ckpt_path(id_alvo))
      message("      [ckpt] origem marcada como FALHA e registrada.")
      next
    }
    
    # --- Target-month predictions ---------------------------------------------
    obs_alvo <- dados_final_modelo$CASOS[dados_final_modelo$ID_TEMPO == id_alvo]
    
    pred_df <- data.frame(
      ID_TEMPO  = id_alvo,
      DATA      = dt$DATA[idx_alvo],
      ANO       = dt$ANO[idx_alvo],
      MES       = dt$MES[idx_alvo],
      COD_GEO   = dt$COD_GEO[idx_alvo],
      DRS_ID    = if ("DRS_ID" %in% names(dt)) dt$DRS_ID[idx_alvo] else NA_integer_,
      POP_TOTAL = if ("POP_TOTAL" %in% names(dt)) dt$POP_TOTAL[idx_alvo] else NA_real_,
      OBS       = obs_alvo,
      PRED      = m$summary.fitted.values$mean[idx_alvo],
      # Posterior sd of the fitted mean; required to build predictive intervals.
      SD_PRED   = m$summary.fitted.values$sd[idx_alvo],
      LI        = m$summary.fitted.values$`0.025quant`[idx_alvo],
      LS        = m$summary.fitted.values$`0.975quant`[idx_alvo],
      N_TREINO  = id_alvo - 1L,
      stringsAsFactors = FALSE
    )
    
    # --- Relative risks for this origin ---------------------------------------
    rr_df <- tryCatch({
      base_rr <- m$summary.fixed %>%
        tibble::rownames_to_column("Termo") %>%
        mutate(RR    = exp(mean),
               LI_rr = exp(`0.025quant`),
               LS_rr = exp(`0.975quant`)) %>%
        select(Termo, RR, LI_rr, LS_rr)
      
      for (nm in c("CASOS_LAG1_LOG", "CASOS_LAG2_LOG")) {
        if (!is.null(m$summary.random[[nm]])) {
          extra <- m$summary.random[[nm]] %>%
            mutate(Termo = nm, RR = exp(mean),
                   LI_rr = exp(`0.025quant`), LS_rr = exp(`0.975quant`)) %>%
            select(Termo, RR, LI_rr, LS_rr)
          base_rr <- bind_rows(base_rr, extra)
        }
      }
      base_rr %>% mutate(ID_TEMPO = id_alvo, DATA = linha_t$DATA,
                         ANO = linha_t$ANO)
    }, error = function(e) NULL)
    
    # --- Hyperparameters (for stability diagnostics) --------------------------
    hyper_df <- tryCatch({
      m$summary.hyperpar %>%
        tibble::rownames_to_column("Hiper") %>%
        select(Hiper, mean, sd) %>%
        mutate(ID_TEMPO = id_alvo)
    }, error = function(e) NULL)
    
    theta_prev <- tryCatch(m$mode$theta, error = function(e) theta_prev)
    
    # --- CHECKPOINT (small objects only) --------------------------------------
    saveRDS(
      list(id_tempo = id_alvo, data = linha_t$DATA, ok = TRUE,
           pred = pred_df, rr = rr_df, hyper = hyper_df,
           theta = theta_prev,
           segundos = as.numeric(difftime(Sys.time(), t0, units = "secs"))),
      ckpt_path(id_alvo)
    )
    
    dur <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2)
    resta <- length(se_fazer) - k
    message(sprintf("      OK em %.2f min | MAE do mes: %.2f | restam ~%d origens (~%.1f h)",
                    dur,
                    mean(abs(pred_df$OBS - pred_df$PRED), na.rm = TRUE),
                    resta, resta * dur / 60))
    
    rm(m, dt); gc(verbose = FALSE)
  }
  
  message(sprintf("\n[loop] concluido em %.2f h.",
                  as.numeric(difftime(Sys.time(), t_inicio_geral, units = "hours"))))
  
  # ----------------------------------------------------------------------------
  # 5. CONSOLIDATION FROM CHECKPOINTS
  # ----------------------------------------------------------------------------
  message("\n--- Consolidando checkpoints ---")
  
  arquivos <- sort(list.files(DIR_CKPT, pattern = "^origem_\\d+\\.rds$",
                              full.names = TRUE))
  ckpts <- lapply(arquivos, function(f) tryCatch(readRDS(f), error = function(e) NULL))
  ckpts <- Filter(function(x) !is.null(x) && isTRUE(x$ok), ckpts)
  
  if (length(ckpts) == 0) {
    
    message("Módulo 24: nenhum checkpoint valido. Nada a consolidar.")
    
  } else {
    
    df_roll   <- bind_rows(lapply(ckpts, `[[`, "pred"))
    df_rr     <- bind_rows(lapply(ckpts, `[[`, "rr"))
    df_hyper  <- bind_rows(lapply(ckpts, `[[`, "hyper"))
    tempos    <- vapply(ckpts, function(x) x$segundos %||% NA_real_, numeric(1))
    
    n_falhas <- length(arquivos) - length(ckpts)
    if (n_falhas > 0)
      message(sprintf("[aviso] %d origens falharam e foram excluidas das metricas.",
                      n_falhas))
    
    message(sprintf("[dados] %d origens | %d linhas municipio-mes preditas.",
                    dplyr::n_distinct(df_roll$ID_TEMPO), nrow(df_roll)))
    
    # --- Categoria de porte populacional (se ausente no pipeline) ------------
    if (!"CATEGORIA_POP" %in% names(df_roll) && !all(is.na(df_roll$POP_TOTAL))) {
      df_roll <- df_roll %>%
        mutate(CATEGORIA_POP = cut(
          POP_TOTAL,
          breaks = c(-Inf, 20000, 100000, 500000, Inf),
          labels = c("1. Small (<20k)", "2. Medium (20-100k)",
                     "3. Large (100-500k)", "4. Metropolis (>500k)")))
    }
    
    # --------------------------------------------------------------------------
    # 6. METRICS
    # --------------------------------------------------------------------------
    TH <- as.numeric(quantile(df_roll$OBS, CFG$q_surto, na.rm = TRUE))
    
    calc_metricas <- function(d) {
      d %>% summarise(
        N_OBS           = dplyr::n(),
        MAE             = mean(abs(OBS - PRED), na.rm = TRUE),
        RMSE            = sqrt(mean((OBS - PRED)^2, na.rm = TRUE)),
        # Legacy definition (denominator = population), kept for comparison:
        MAE_PCT         = (mean(abs(OBS - PRED), na.rm = TRUE) /
                             mean(POP_TOTAL, na.rm = TRUE)) * 100,
        # Relative error against the mean observed count.
        MAE_REL_OBS_PCT = (mean(abs(OBS - PRED), na.rm = TRUE) /
                             mean(OBS, na.rm = TRUE)) * 100,
        BIAS            = mean(PRED - OBS, na.rm = TRUE),
        COBERTURA_PCT   = mean(OBS >= floor(LI) & OBS <= ceiling(LS),
                               na.rm = TRUE) * 100,
        SPEARMAN        = suppressWarnings(
          cor(OBS, PRED, method = "spearman", use = "complete.obs")),
        AUC             = as.numeric(tryCatch(
          pROC::roc(ifelse(OBS > TH, 1, 0), PRED, quiet = TRUE)$auc,
          error = function(e) NA_real_)),
        HIT_RATE_PCT    = sum(OBS > TH & PRED > TH, na.rm = TRUE) /
          max(1, sum(OBS > TH, na.rm = TRUE)) * 100,
        FALSE_ALARM_PCT = sum(OBS <= TH & PRED > TH, na.rm = TRUE) /
          max(1, sum(OBS <= TH, na.rm = TRUE)) * 100,
        MAX_OBS         = max(OBS,  na.rm = TRUE),
        MAX_PRED        = max(PRED, na.rm = TRUE),
        .groups = "drop"
      )
    }
    
    met_global   <- calc_metricas(df_roll)
    met_ano      <- df_roll %>% group_by(ANO) %>% calc_metricas()
    met_mes      <- df_roll %>% group_by(MES) %>% calc_metricas()
    met_porte    <- if ("CATEGORIA_POP" %in% names(df_roll))
      df_roll %>% group_by(CATEGORIA_POP) %>% calc_metricas() else NULL
    met_drs      <- if (!all(is.na(df_roll$DRS_ID)))
      df_roll %>% group_by(DRS_ID) %>% calc_metricas() else NULL
    
    # Human-readable DRS name, if the pipeline provides it
    if (!is.null(met_drs) && exists("df_geo_drs")) {
      met_drs <- met_drs %>%
        left_join(df_geo_drs %>% distinct(DRS_ID, NOME_DRS), by = "DRS_ID")
    }
    
    # --------------------------------------------------------------------------
    # 7. EXPORT
    # --------------------------------------------------------------------------
    wcsv <- function(x, nome) if (!is.null(x))
      tryCatch(write.csv(x, file.path(DIR_ROLL, nome), row.names = FALSE),
               error = function(e) message("   [aviso] falha ao salvar ", nome))
    
    wcsv(met_global, "TABELA_ROLLING_MENSAL_METRICAS.csv")
    wcsv(met_ano,    "TABELA_ROLLING_MENSAL_POR_ANO.csv")
    wcsv(met_mes,    "TABELA_ROLLING_MENSAL_POR_MES.csv")
    wcsv(met_porte,  "TABELA_ROLLING_MENSAL_POR_PORTE.csv")
    wcsv(met_drs,    "TABELA_ROLLING_MENSAL_POR_DRS.csv")
    wcsv(df_rr,      "TABELA_ROLLING_MENSAL_RR_EVOLUCAO.csv")
    wcsv(df_hyper,   "TABELA_ROLLING_MENSAL_HIPERPARAMETROS.csv")
    wcsv(df_roll,    "PREDICOES_ROLLING_MENSAL.csv")
    
    sink(file.path(DIR_ROLL, "TABELA_ROLLING_MENSAL_METRICAS.txt"))
    cat("================================================================\n")
    cat("  MONTHLY ROLLING-ORIGIN CROSS-VALIDATION (LEAKAGE-FREE)\n")
    cat("================================================================\n\n")
    cat("Design: for each target month t, the model is refitted using only\n")
    cat("ID_TEMPO <= t-1 and predicts month t. The autoregressive predictors\n")
    cat("for month t are months t-1 and t-2, both inside the training set.\n")
    cat("No observation outside the training window is used at any point.\n\n")
    cat(sprintf("Origins evaluated : %d\n", dplyr::n_distinct(df_roll$ID_TEMPO)))
    cat(sprintf("Period            : %s to %s\n",
                format(min(df_roll$DATA), "%Y-%m"),
                format(max(df_roll$DATA), "%Y-%m")))
    cat(sprintf("Outbreak threshold: %.2f cases (quantile %.2f of observed)\n",
                TH, CFG$q_surto))
    cat(sprintf("Mean fit time     : %.2f min per origin\n\n",
                mean(tempos, na.rm = TRUE) / 60))
    cat("--- GLOBAL METRICS ---\n")
    cat(knitr::kable(met_global, format = "markdown", digits = 3), sep = "\n")
    cat("\n\n--- BY YEAR ---\n")
    cat(knitr::kable(met_ano, format = "markdown", digits = 3), sep = "\n")
    cat("\n\n--- BY CALENDAR MONTH ---\n")
    cat(knitr::kable(met_mes, format = "markdown", digits = 3), sep = "\n")
    if (!is.null(met_porte)) {
      cat("\n\n--- BY POPULATION SIZE ---\n")
      cat(knitr::kable(met_porte, format = "markdown", digits = 3), sep = "\n")
    }
    if (!is.null(met_drs)) {
      cat("\n\n--- BY REGIONAL HEALTH DEPARTMENT (DRS) ---\n")
      cat(knitr::kable(met_drs, format = "markdown", digits = 3), sep = "\n")
    }
    cat("\n\nNOTE ON RELATIVE ERROR:\n")
    cat("  MAE_PCT         = MAE / mean(POP_TOTAL) * 100   (legacy definition)\n")
    cat("  MAE_REL_OBS_PCT = MAE / mean(OBS)       * 100   (standard definition)\n")
    cat("  State the formula explicitly in the manuscript.\n")
    sink()
    
    message("   -> tabelas exportadas em ", DIR_ROLL)
    
    # --------------------------------------------------------------------------
    # 8. PUBLICATION FIGURE (4-panel)
    # --------------------------------------------------------------------------
    tema_pub <- theme_bw(base_size = 12) +
      theme(
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
        panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.7),
        axis.text        = element_text(colour = "black"),
        axis.title       = element_text(face = "bold"),
        plot.title       = element_text(face = "bold", size = rel(1.0)),
        plot.subtitle    = element_text(colour = "grey35", size = rel(0.82)),
        legend.position  = "bottom",
        legend.title     = element_blank(),
        legend.background = element_blank(),
        plot.margin      = margin(8, 10, 8, 8)
      )
    
    COR_OBS  <- "#2C3E50"
    COR_PRED <- "#B2182B"
    COR_AUX  <- "#2166AC"
    
    # --- Panel A: state-aggregated series -------------------------------------
    df_estado <- df_roll %>%
      group_by(DATA) %>%
      summarise(OBS = sum(OBS, na.rm = TRUE),
                PRED = sum(PRED, na.rm = TRUE), .groups = "drop")
    
    pA <- ggplot(df_estado, aes(x = DATA)) +
      geom_line(aes(y = OBS,  colour = "Observed"),  linewidth = 0.85) +
      geom_line(aes(y = PRED, colour = "Predicted"), linewidth = 0.85,
                linetype = "22") +
      scale_colour_manual(values = c(Observed = COR_OBS, Predicted = COR_PRED)) +
      scale_y_continuous(labels = scales::comma) +
      scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
      labs(title = "A. Observed versus predicted monthly cases, São Paulo State",
           subtitle = "Each point predicted from a model refitted on all data up to the preceding month",
           x = NULL, y = "Monthly cases") +
      tema_pub
    
    # --- Panel B: municipality-month scatter (log-log) ------------------------
    df_sc <- df_roll %>% filter(OBS > 0, PRED > 0)
    pB <- ggplot(df_sc, aes(x = OBS, y = PRED)) +
      geom_point(alpha = 0.06, size = 0.5, colour = COR_AUX) +
      geom_abline(slope = 1, intercept = 0, colour = COR_PRED,
                  linetype = "dashed", linewidth = 0.6) +
      scale_x_log10(labels = scales::comma) +
      scale_y_log10(labels = scales::comma) +
      annotate("text", x = min(df_sc$OBS), y = max(df_sc$PRED),
               hjust = 0, vjust = 1, size = 3.4,
               label = sprintf("Spearman = %.3f\nAUC = %.3f",
                               met_global$SPEARMAN, met_global$AUC)) +
      labs(title = "B. Municipality-month agreement",
           subtitle = "Log-log scale; dashed line is the 1:1 identity",
           x = "Observed cases", y = "Predicted cases") +
      tema_pub
    
    # --- Painel C: erro relativo por ano -------------------------------------
    pC <- ggplot(met_ano, aes(x = ANO, y = MAE_REL_OBS_PCT)) +
      geom_col(fill = COR_AUX, alpha = 0.85, width = 0.65) +
      geom_text(aes(label = sprintf("%.1f", MAE_REL_OBS_PCT)),
                vjust = -0.4, size = 2.9) +
      scale_x_continuous(breaks = met_ano$ANO) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
      labs(title = "C. Relative error by predicted year",
           subtitle = "No systematic upward trend indicates no dependence on specific historical periods",
           x = NULL, y = "Relative MAE (%)") +
      tema_pub +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    # --- Painel D: estabilidade dos RR ao longo das origens -------------------
    termos_foco <- c("CASOS_LAG1_LOG", "CASOS_LAG2_LOG")
    if (exists("vars_social_on")) termos_foco <- c(termos_foco, vars_social_on)
    termos_foco <- intersect(termos_foco, unique(df_rr$Termo))
    
    pD <- NULL
    if (length(termos_foco) > 0) {
      df_rr_plot <- df_rr %>% filter(Termo %in% termos_foco)
      pD <- ggplot(df_rr_plot, aes(x = DATA, y = RR, colour = Termo)) +
        geom_hline(yintercept = 1, linetype = "dashed", colour = "grey45") +
        geom_ribbon(aes(ymin = LI_rr, ymax = LS_rr, fill = Termo),
                    alpha = 0.13, colour = NA) +
        geom_line(linewidth = 0.8) +
        scale_colour_viridis_d(option = "turbo", end = 0.85) +
        scale_fill_viridis_d(option = "turbo", end = 0.85) +
        scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
        labs(title = "D. Stability of relative risks across rolling origins",
             subtitle = "Each point re-estimated from an expanding training window",
             x = NULL, y = "Relative Risk (RR)") +
        tema_pub +
        theme(legend.text = element_text(size = 7.5))
    }
    
    # --- Montagem -------------------------------------------------------------
    paineis <- Filter(Negate(is.null), list(pA, pB, pC, pD))
    salvou <- FALSE
    
    if (requireNamespace("patchwork", quietly = TRUE) && length(paineis) == 4) {
      fig <- patchwork::wrap_plots(pA, pB, pC, pD, ncol = 2) +
        patchwork::plot_annotation(
          title = "Monthly rolling-origin cross-validation (leakage-free design)",
          subtitle = sprintf("%d origins, %s to %s | model refitted at every origin using only prior data",
                             dplyr::n_distinct(df_roll$ID_TEMPO),
                             format(min(df_roll$DATA), "%b %Y"),
                             format(max(df_roll$DATA), "%b %Y")),
          theme = theme(
            plot.title    = element_text(face = "bold", size = 15),
            plot.subtitle = element_text(colour = "grey35", size = 10)))
      ggsave(file.path(DIR_ROLL, "FIG_RollingOrigin_Painel_Publicacao.png"),
             fig, width = 15.5, height = 11, dpi = 600, bg = "white")
      salvou <- TRUE
    }
    
    if (!salvou) {
      message("   [aviso] 'patchwork' indisponivel — salvando paineis separados.")
      ggsave(file.path(DIR_ROLL, "FIG_RollingOrigin_A_serie.png"),
             pA, width = 11, height = 5.2, dpi = 600, bg = "white")
      ggsave(file.path(DIR_ROLL, "FIG_RollingOrigin_B_dispersao.png"),
             pB, width = 7, height = 6, dpi = 600, bg = "white")
      ggsave(file.path(DIR_ROLL, "FIG_RollingOrigin_C_erro_ano.png"),
             pC, width = 8, height = 5, dpi = 600, bg = "white")
      if (!is.null(pD))
        ggsave(file.path(DIR_ROLL, "FIG_RollingOrigin_D_rr.png"),
               pD, width = 9, height = 5.5, dpi = 600, bg = "white")
    }
    
    # Stand-alone high-resolution figures
    ggsave(file.path(DIR_ROLL, "FIG_RollingOrigin_SerieEstado.png"),
           pA, width = 12, height = 5.5, dpi = 700, bg = "white")
    if (!is.null(pD))
      ggsave(file.path(DIR_ROLL, "FIG_RollingOrigin_RR_Estabilidade.png"),
             pD, width = 10, height = 6, dpi = 700, bg = "white")
    
    message("   -> figuras salvas em ", DIR_ROLL)
    
    # --------------------------------------------------------------------------
    # 9. OBJETOS GLOBAIS PARA USO POSTERIOR
    # --------------------------------------------------------------------------
    df_rolling_mensal      <<- df_roll
    metrics_rolling_mensal <<- met_global
    rr_rolling_mensal      <<- df_rr
    
    # OPTIONAL bridge to the diagnostics module (off by default:
    # overwriting df_cv silently is undesirable).
    if (isTRUE(CFG$ponte_df_cv)) {
      df_cv      <<- df_roll
      metrics_cv <<- met_global
      message("   [ponte] df_cv e metrics_cv sobrescritos pelo rolling mensal.")
    }
    
    message("\n   Módulo 24 concluído.")
    message(sprintf("   MAE            : %.3f", met_global$MAE))
    message(sprintf("   MAE rel (OBS)  : %.2f%%", met_global$MAE_REL_OBS_PCT))
    message(sprintf("   MAE_PCT legado : %.3f%%", met_global$MAE_PCT))
    message(sprintf("   AUC            : %.3f", met_global$AUC))
    message(sprintf("   Spearman       : %.3f", met_global$SPEARMAN))
    message(sprintf("   Hit rate       : %.2f%%", met_global$HIT_RATE_PCT))
    message(sprintf("   False alarm    : %.2f%%", met_global$FALSE_ALARM_PCT))
  }
}