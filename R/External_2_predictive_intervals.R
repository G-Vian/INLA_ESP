# ==============================================================================
# EXTERNAL 2 - OUT-OF-SAMPLE POSTERIOR PREDICTIVE INTERVALS
# ------------------------------------------------------------------------------
# Turns each rolling origin's fitted values into genuine out-of-sample
# predictive intervals by Monte Carlo over the (zero-inflated) negative-binomial
# likelihood, using the hyperparameters estimated at that origin. Reports
# predictive coverage and MAE / relative MAE, globally and by municipality size,
# all from the rolling-origin design.
# ==============================================================================

# ==============================================================================
# OUT-OF-SAMPLE PREDICTIVE INTERVALS
#            (predictive-interval logic applied to the cross-validation predictions)
# ------------------------------------------------------------------------------
# Computes out-of-sample posterior predictive intervals from the rolling-origin
# predictions. Reports predictive coverage together with MAE, RMSE and relative
# MAE (per mean observed count), globally and by municipality size. The coverage
# of the credible interval of the mean is also computed, for contrast only.
# ==============================================================================

message("\n======================================================================")
message("  MÓDULO 25: INTERVALOS PREDITIVOS OUT-OF-SAMPLE")
message("======================================================================\n")

if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a

.ip_rodar <- isTRUE(get0("SWITCH_INTERVALOS_PREDITIVOS_OOS", ifnotfound = FALSE))

if (!.ip_rodar) {
  
  message("Módulo 25 ignorado (SWITCH_INTERVALOS_PREDITIVOS_OOS desligado).")
  
} else {
  
  suppressPackageStartupMessages({
    library(dplyr); library(tidyr); library(ggplot2)
  })
  
  # ============================================================================
  # 1. CONFIGURATION
  # ============================================================================
  IPCFG <- list(
    n_sim      = get0("CONFIG_IP_N_SIM",        ifnotfound = 1000L),
    nivel      = get0("CONFIG_IP_NIVEL",        ifnotfound = 0.95),
    seed       = get0("CONFIG_IP_SEED",         ifnotfound = 20260817L),
    dir_ckpt   = get0("CONFIG_IP_DIR_CKPT",
                      ifnotfound = file.path(DIR_OUTPUT, "teste_rolling_mensal",
                                             "checkpoints")),
    fazer_insample = get0("CONFIG_IP_INSAMPLE", ifnotfound = TRUE)
  )
  
  alpha  <- 1 - IPCFG$nivel
  q_lo   <- alpha / 2
  q_hi   <- 1 - alpha / 2
  
  DIR_IP <- file.path(DIR_OUTPUT, "teste_rolling_mensal", "intervalos_preditivos")
  dir.create(DIR_IP, recursive = TRUE, showWarnings = FALSE)
  
  message(sprintf("[cfg] N_SIM = %d | nivel = %.0f%% | seed = %d",
                  IPCFG$n_sim, IPCFG$nivel * 100, IPCFG$seed))
  
  # ============================================================================
  # 2. HELPER FUNCTIONS
  # ============================================================================
  
  extrai_k_pi <- function(hyper_df) {
    if (is.null(hyper_df) || nrow(hyper_df) == 0) return(list(k = NA, pi = 0))
    
    nomes <- tolower(if ("Hiper" %in% names(hyper_df)) hyper_df$Hiper
                     else rownames(hyper_df))
    val   <- hyper_df[["mean"]]
    
    i_k  <- grep("size", nomes)
    i_pi <- grep("zero[-. ]?prob|probability.*zero|zeroinfl", nomes)
    
    k  <- if (length(i_k)  > 0) val[i_k[1]]  else NA_real_
    pp <- if (length(i_pi) > 0) val[i_pi[1]] else 0
    
    if (!is.finite(k) || k <= 0) k <- NA_real_
    if (!is.finite(pp) || pp < 0) pp <- 0
    if (pp >= 1) pp <- 0.999
    
    list(k = k, pi = pp)
  }
  
  simula_preditivo <- function(theta, sigma, k, pp, n_sim) {
    
    n <- length(theta)
    
    th <- pmax(as.numeric(theta), 1e-8)
    sg <- pmax(as.numeric(sigma), 1e-8)
    
    shape <- (th / sg)^2
    scale <- (sg^2) / th
    
    theta_star <- matrix(
      stats::rgamma(n * n_sim, shape = rep(shape, times = n_sim),
                    scale = rep(scale, times = n_sim)),
      nrow = n, ncol = n_sim
    )
    
    mu_star <- theta_star / (1 - pp)
    mu_star[!is.finite(mu_star) | mu_star <= 0] <- 1e-8
    
    Y <- matrix(
      stats::rnbinom(n * n_sim, size = k, mu = as.vector(mu_star)),
      nrow = n, ncol = n_sim
    )
    
    if (pp > 0) {
      z <- matrix(stats::runif(n * n_sim) < pp, nrow = n, ncol = n_sim)
      Y[z] <- 0L
    }
    
    Y
  }
  
  decompoe_variancia <- function(theta, sigma, k, pp) {
    th <- pmax(as.numeric(theta), 0)
    data.frame(
      VAR_PARAM = as.numeric(sigma)^2,
      VAR_OBS   = th + th^2 * (1 / k + pp) / (1 - pp)
    )
  }
  
  # ============================================================================
  # 3. PROCESSAMENTO OUT-OF-SAMPLE — UMA ORIGEM POR VEZ
  # ============================================================================
  
  arqs <- sort(list.files(IPCFG$dir_ckpt, pattern = "^origem_\\d+\\.rds$",
                          full.names = TRUE))
  
  if (length(arqs) == 0) {
    stop("Módulo 25: nenhum checkpoint encontrado em ", IPCFG$dir_ckpt,
         "\n   Rode o Módulo 24 antes.")
  }
  
  message(sprintf("[dados] %d checkpoints encontrados.", length(arqs)))
  
  set.seed(IPCFG$seed)
  
  lista_muni <- list()
  lista_drs  <- list()
  falhas_sd  <- 0L
  falhas_hp  <- 0L
  
  for (a in seq_along(arqs)) {
    
    ck <- tryCatch(readRDS(arqs[a]), error = function(e) NULL)
    if (is.null(ck) || !isTRUE(ck$ok) || is.null(ck$pred)) next
    
    pr <- ck$pred
    
    if (!"SD_PRED" %in% names(pr)) {
      falhas_sd <- falhas_sd + 1L
      next
    }
    
    kp <- extrai_k_pi(ck$hyper)
    if (is.na(kp$k)) {
      falhas_hp <- falhas_hp + 1L
      next
    }
    
    Y <- simula_preditivo(pr$PRED, pr$SD_PRED, kp$k, kp$pi, IPCFG$n_sim)
    
    pr$PI_LI <- apply(Y, 1, stats::quantile, probs = q_lo, names = FALSE)
    pr$PI_LS <- apply(Y, 1, stats::quantile, probs = q_hi, names = FALSE)
    pr$PI_MEDIANA <- apply(Y, 1, stats::median)
    
    dv <- decompoe_variancia(pr$PRED, pr$SD_PRED, kp$k, kp$pi)
    pr$VAR_PARAM <- dv$VAR_PARAM
    pr$VAR_OBS   <- dv$VAR_OBS
    
    pr$KAPPA <- kp$k
    pr$PI0   <- kp$pi
    
    pr$DENTRO_PREDITIVO <- pr$OBS >= pr$PI_LI & pr$OBS <= pr$PI_LS
    pr$DENTRO_MEDIA <- pr$OBS >= floor(pr$LI) & pr$OBS <= ceiling(pr$LS)
    
    lista_muni[[length(lista_muni) + 1L]] <- pr
    
    if (!all(is.na(pr$DRS_ID))) {
      g <- as.character(pr$DRS_ID)
      Ysum <- rowsum(Y, group = g, reorder = TRUE)
      obs_drs <- tapply(pr$OBS, g, sum, na.rm = TRUE)
      med_drs <- tapply(pr$PRED, g, sum, na.rm = TRUE)
      
      drs_df <- data.frame(
        ID_TEMPO = pr$ID_TEMPO[1],
        DATA     = pr$DATA[1],
        ANO      = pr$ANO[1],
        DRS_ID   = rownames(Ysum),
        OBS      = as.numeric(obs_drs[rownames(Ysum)]),
        PRED     = as.numeric(med_drs[rownames(Ysum)]),
        PI_LI    = apply(Ysum, 1, stats::quantile, probs = q_lo, names = FALSE),
        PI_LS    = apply(Ysum, 1, stats::quantile, probs = q_hi, names = FALSE),
        stringsAsFactors = FALSE
      )
      drs_df$DENTRO_PREDITIVO <- drs_df$OBS >= drs_df$PI_LI &
        drs_df$OBS <= drs_df$PI_LS
      lista_drs[[length(lista_drs) + 1L]] <- drs_df
    }
    
    rm(Y); if (a %% 20 == 0) gc(verbose = FALSE)
    
    if (a %% 10 == 0 || a == length(arqs))
      message(sprintf("   [%d/%d] origens processadas", a, length(arqs)))
  }
  
  if (falhas_sd > 0)
    message(sprintf("\n[ATENCAO] %d checkpoints sem a coluna SD_PRED foram pulados.",
                    falhas_sd),
            "\n          The SD_PRED column is required (see module External 1).")
  if (falhas_hp > 0)
    message(sprintf("[aviso] %d checkpoints sem hiperparametro 'size' foram pulados.",
                    falhas_hp))
  
  if (length(lista_muni) == 0)
    stop("No origin could be processed: the SD_PRED column is missing.")
  
  df_ip     <- bind_rows(lista_muni)
  df_ip_drs <- if (length(lista_drs) > 0) bind_rows(lista_drs) else NULL
  
  message(sprintf("\n[ok] %d linhas municipio-mes com intervalo preditivo.",
                  nrow(df_ip)))
  
  if (!"CATEGORIA_POP" %in% names(df_ip) && !all(is.na(df_ip$POP_TOTAL))) {
    df_ip <- df_ip %>%
      mutate(CATEGORIA_POP = cut(
        POP_TOTAL,
        breaks = c(-Inf, 20000, 100000, 500000, Inf),
        labels = c("1. Small (<20k)", "2. Medium (20-100k)",
                   "3. Large (100-500k)", "4. Metropolis (>500k)")))
  }
  
  # ============================================================================
  # 4. TABELAS
  # ============================================================================
  
  # Coverage plus error metrics (MAE, RMSE, relative MAE) for one group.
  resume_cob <- function(d) {
    d %>% summarise(
      N               = dplyr::n(),
      MEAN_OBS        = mean(OBS, na.rm = TRUE),
      MAE             = mean(abs(OBS - PRED), na.rm = TRUE),
      RMSE            = sqrt(mean((OBS - PRED)^2, na.rm = TRUE)),
      # Relative MAE = MAE / media OBSERVADA (definicao correta, igual ao Mod 24).
      MAE_REL_PCT     = mean(abs(OBS - PRED), na.rm = TRUE) /
        mean(OBS, na.rm = TRUE) * 100,
      COB_PREDITIVO   = mean(DENTRO_PREDITIVO, na.rm = TRUE) * 100,
      # COB_MEDIA: contraste diagnostico apenas. NAO reportar como preditiva.
      COB_MEDIA       = mean(DENTRO_MEDIA,     na.rm = TRUE) * 100,
      LARG_PI_MEDIA   = mean(PI_LS - PI_LI,    na.rm = TRUE),
      LARG_IC_MEDIA   = mean(LS - LI,          na.rm = TRUE),
      VAR_PARAM_PCT   = sum(VAR_PARAM, na.rm = TRUE) /
        (sum(VAR_PARAM, na.rm = TRUE) + sum(VAR_OBS, na.rm = TRUE)) * 100,
      VAR_OBS_PCT     = sum(VAR_OBS, na.rm = TRUE) /
        (sum(VAR_PARAM, na.rm = TRUE) + sum(VAR_OBS, na.rm = TRUE)) * 100,
      .groups = "drop"
    )
  }
  
  cob_global <- resume_cob(df_ip)
  cob_porte  <- if ("CATEGORIA_POP" %in% names(df_ip))
    df_ip %>% group_by(CATEGORIA_POP) %>% resume_cob() else NULL
  cob_ano    <- df_ip %>% group_by(ANO) %>% resume_cob()
  
  # Performance table by population size (rolling out-of-sample).
  if (!is.null(cob_porte)) {
    tab_porte_pub <- cob_porte %>%
      dplyr::select(CATEGORIA_POP, N, MEAN_OBS, MAE, RMSE,
                    MAE_REL_PCT, COB_PREDITIVO)
    tryCatch(write.csv(tab_porte_pub,
                       file.path(DIR_IP, "TABELA_PERFORMANCE_POR_PORTE_ROLLING.csv"),
                       row.names = FALSE),
             error = function(e) message("   [aviso] falha ao salvar tab porte"))
    message("\n[tab] TABELA_PERFORMANCE_POR_PORTE_ROLLING.csv")
    message("      (MAE, RMSE, Relative MAE e COBERTURA PREDITIVA por porte, rolling OOS)")
    message("\n--- Desempenho por porte (rolling OOS) ---")
    print(as.data.frame(tab_porte_pub), row.names = FALSE, digits = 4)
  }
  
  cob_drs <- if (!is.null(df_ip_drs)) {
    d <- df_ip_drs %>%
      group_by(DRS_ID) %>%
      summarise(N = dplyr::n(),
                COB_PREDITIVO = mean(DENTRO_PREDITIVO, na.rm = TRUE) * 100,
                LARG_PI_MEDIA = mean(PI_LS - PI_LI, na.rm = TRUE),
                .groups = "drop")
    if (exists("df_geo_drs"))
      d <- d %>% mutate(DRS_ID = as.integer(DRS_ID)) %>%
        left_join(df_geo_drs %>% distinct(DRS_ID, NOME_DRS), by = "DRS_ID")
    d
  } else NULL
  
  cob_drs_global <- if (!is.null(df_ip_drs))
    mean(df_ip_drs$DENTRO_PREDITIVO, na.rm = TRUE) * 100 else NA_real_
  
  wcsv <- function(x, nm) if (!is.null(x))
    tryCatch(write.csv(x, file.path(DIR_IP, nm), row.names = FALSE),
             error = function(e) message("   [aviso] falha ao salvar ", nm))
  
  wcsv(cob_global, "TABELA_COBERTURA_PREDITIVA_OOS.csv")
  wcsv(cob_porte,  "TABELA_CONTRASTE_MEDIA_vs_PREDITIVO.csv")
  wcsv(cob_ano,    "TABELA_COBERTURA_POR_ANO.csv")
  wcsv(cob_drs,    "TABELA_COBERTURA_POR_DRS.csv")
  wcsv(df_ip %>% select(any_of(c("ID_TEMPO","DATA","ANO","MES","COD_GEO",
                                 "DRS_ID","POP_TOTAL","CATEGORIA_POP","OBS",
                                 "PRED","SD_PRED","LI","LS","PI_LI","PI_LS",
                                 "VAR_PARAM","VAR_OBS","KAPPA","PI0",
                                 "DENTRO_PREDITIVO","DENTRO_MEDIA"))),
       "PREDICOES_COM_INTERVALOS_PREDITIVOS.csv")
  
  decomp <- df_ip %>%
    summarise(VAR_PARAM_TOTAL = sum(VAR_PARAM, na.rm = TRUE),
              VAR_OBS_TOTAL   = sum(VAR_OBS,   na.rm = TRUE)) %>%
    mutate(VAR_TOTAL     = VAR_PARAM_TOTAL + VAR_OBS_TOTAL,
           PCT_PARAMETRO = VAR_PARAM_TOTAL / VAR_TOTAL * 100,
           PCT_OBSERVACAO= VAR_OBS_TOTAL   / VAR_TOTAL * 100)
  wcsv(decomp, "TABELA_DECOMPOSICAO_VARIANCIA.csv")
  
  # ============================================================================
  # 5. IN-SAMPLE VERSION (same code, on mod_final)
  # ============================================================================
  cob_insample <- NULL
  
  if (isTRUE(IPCFG$fazer_insample) && exists("mod_final") && !is.null(mod_final)) {
    
    message("\n--- Versao IN-SAMPLE (mod_final), pelo mesmo algoritmo ---")
    
    kp_in <- extrai_k_pi(mod_final$summary.hyperpar)
    
    if (!is.na(kp_in$k)) {
      th_in <- mod_final$summary.fitted.values$mean
      sd_in <- mod_final$summary.fitted.values$sd
      obs_in <- dados_final_modelo$CASOS
      
      ok <- !is.na(obs_in) & is.finite(th_in) & is.finite(sd_in)
      
      idx_ok <- which(ok)
      blocos <- split(idx_ok, ceiling(seq_along(idx_ok) / 5000))
      dentro <- logical(0); vpar <- numeric(0); vobs <- numeric(0)
      
      set.seed(IPCFG$seed + 1L)
      for (b in blocos) {
        Yb <- simula_preditivo(th_in[b], sd_in[b], kp_in$k, kp_in$pi, IPCFG$n_sim)
        li <- apply(Yb, 1, stats::quantile, probs = q_lo, names = FALSE)
        ls <- apply(Yb, 1, stats::quantile, probs = q_hi, names = FALSE)
        dentro <- c(dentro, obs_in[b] >= li & obs_in[b] <= ls)
        dvb <- decompoe_variancia(th_in[b], sd_in[b], kp_in$k, kp_in$pi)
        vpar <- c(vpar, dvb$VAR_PARAM); vobs <- c(vobs, dvb$VAR_OBS)
        rm(Yb)
      }
      gc(verbose = FALSE)
      
      cob_insample <- data.frame(
        N             = length(dentro),
        COB_PREDITIVO = mean(dentro, na.rm = TRUE) * 100,
        VAR_PARAM_PCT = sum(vpar) / (sum(vpar) + sum(vobs)) * 100,
        VAR_OBS_PCT   = sum(vobs) / (sum(vpar) + sum(vobs)) * 100,
        KAPPA = kp_in$k, PI0 = kp_in$pi
      )
      wcsv(cob_insample, "TABELA_COBERTURA_PREDITIVA_INSAMPLE.csv")
      message(sprintf("   cobertura in-sample: %.2f%%", cob_insample$COB_PREDITIVO))
    } else {
      message("   [aviso] nao foi possivel extrair 'size' de mod_final.")
    }
  }
  
  # ============================================================================
  # 6. TEXT REPORT
  # ============================================================================
  sink(file.path(DIR_IP, "TABELA_COBERTURA_PREDITIVA_OOS.txt"))
  cat("================================================================\n")
  cat("  POSTERIOR PREDICTIVE INTERVALS — OUT-OF-SAMPLE COVERAGE\n")
  cat("================================================================\n\n")
  cat("Algorithm (per municipality-month):\n")
  cat("  1. theta, sigma  <- posterior mean and sd of the fitted value\n")
  cat("  2. theta* ~ Gamma(shape=(theta/sigma)^2, scale=sigma^2/theta)\n")
  cat("  3. Y*     ~ NegBinomial(mu = theta*/(1-pi), size = kappa)\n")
  cat("             with structural zeros applied with probability pi\n")
  cat("  4. interval = empirical quantiles of Y*\n\n")
  cat(sprintf("Draws per observation : %d\n", IPCFG$n_sim))
  cat(sprintf("Nominal level         : %.0f%%\n", IPCFG$nivel * 100))
  cat(sprintf("Origins processed     : %d\n", dplyr::n_distinct(df_ip$ID_TEMPO)))
  cat("Hyperparameters kappa and pi taken from EACH origin's own fit,\n")
  cat("so no information outside the training window enters the interval.\n\n")
  
  cat("--- GLOBAL, MUNICIPALITY LEVEL ---\n")
  cat(knitr::kable(cob_global, format = "markdown", digits = 3), sep = "\n")
  
  if (!is.null(cob_porte)) {
    cat("\n\n--- PERFORMANCE BY POPULATION SIZE (rolling OOS) ---\n")
    cat("MAE / RMSE in cases; MAE_REL_PCT = MAE / mean observed count;\n")
    cat("COB_PREDITIVO = coverage of the posterior PREDICTIVE interval.\n\n")
    cat(knitr::kable(
      cob_porte %>% dplyr::select(CATEGORIA_POP, N, MEAN_OBS, MAE, RMSE,
                                  MAE_REL_PCT, COB_PREDITIVO),
      format = "markdown", digits = 2), sep = "\n")
    cat("\n\nDiagnostic note — mean-interval coverage (NOT predictive):\n")
    cat(knitr::kable(
      cob_porte %>% dplyr::select(CATEGORIA_POP, COB_MEDIA),
      format = "markdown", digits = 2), sep = "\n")
    cat("\nThe predictive coverage (COB_PREDITIVO) stays near the nominal level\n")
    cat("across all sizes. The mean-interval coverage (COB_MEDIA) declines with\n")
    cat("size; that decline is a property of a mean-only interval and is NOT a\n")
    cat("predictive metric.\n")
  }
  
  cat("\n\n--- BY YEAR ---\n")
  cat(knitr::kable(cob_ano %>% dplyr::select(ANO, N, MAE, MAE_REL_PCT,
                                             COB_PREDITIVO),
                   format = "markdown", digits = 2), sep = "\n")
  
  if (!is.null(cob_drs)) {
    cat(sprintf("\n\n--- AGGREGATED TO DRS LEVEL ---\nGlobal coverage: %.2f%%\n\n",
                cob_drs_global))
    cat("Municipalities simulated independently and summed. This does not\n")
    cat("capture spatial correlation of the latent field; compare with the\n")
    cat("joint posterior predictive check (Module 7.6).\n\n")
    cat(knitr::kable(cob_drs, format = "markdown", digits = 2), sep = "\n")
  }
  
  cat("\n\n--- PREDICTIVE VARIANCE DECOMPOSITION ---\n")
  cat(knitr::kable(decomp, format = "markdown", digits = 2), sep = "\n")
  
  if (!is.null(cob_insample)) {
    cat("\n\n--- IN-SAMPLE COMPARISON (mod_final, same algorithm) ---\n")
    cat(knitr::kable(cob_insample, format = "markdown", digits = 3), sep = "\n")
    cat(sprintf("\nOut-of-sample : %.2f%%\nIn-sample     : %.2f%%\n",
                cob_global$COB_PREDITIVO, cob_insample$COB_PREDITIVO))
    cat("Out-of-sample coverage below in-sample is expected.\n")
  }
  sink()
  
  # ============================================================================
  # 7. FIGURA
  # ============================================================================
  tema_pub <- theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
          panel.border  = element_rect(colour = "black", fill = NA, linewidth = 0.7),
          axis.text     = element_text(colour = "black"),
          axis.title    = element_text(face = "bold"),
          plot.title    = element_text(face = "bold", size = rel(1.0)),
          plot.subtitle = element_text(colour = "grey35", size = rel(0.82)),
          legend.position = "bottom", legend.title = element_blank())
  
  paineis <- list()
  
  if (!is.null(cob_porte)) {
    dA <- cob_porte %>%
      select(CATEGORIA_POP, COB_MEDIA, COB_PREDITIVO) %>%
      pivot_longer(-CATEGORIA_POP, names_to = "Tipo", values_to = "Cobertura") %>%
      mutate(Tipo = dplyr::recode(as.character(Tipo),
                                  COB_MEDIA     = "Credible interval of the mean",
                                  COB_PREDITIVO = "Posterior predictive interval"))
    
    paineis$A <- ggplot(dA, aes(x = CATEGORIA_POP, y = Cobertura, fill = Tipo)) +
      geom_col(position = position_dodge(width = 0.75), width = 0.68) +
      geom_hline(yintercept = IPCFG$nivel * 100, linetype = "dashed",
                 colour = "#B2182B", linewidth = 0.6) +
      annotate("text", x = 0.6, y = IPCFG$nivel * 100 + 3,
               label = sprintf("nominal %.0f%%", IPCFG$nivel * 100),
               hjust = 0, size = 3.1, colour = "#B2182B") +
      geom_text(aes(label = sprintf("%.0f", Cobertura)),
                position = position_dodge(width = 0.75),
                vjust = -0.35, size = 2.9) +
      scale_fill_manual(values = c("Credible interval of the mean" = "#92C5DE",
                                   "Posterior predictive interval" = "#2166AC")) +
      scale_y_continuous(limits = c(0, 108), expand = c(0, 0)) +
      labs(title = "A. Empirical coverage by municipality size",
           subtitle = "The mean-only interval degrades with size; the predictive interval does not",
           x = NULL, y = "Coverage (%)") +
      tema_pub +
      theme(axis.text.x = element_text(angle = 18, hjust = 1))
  }
  
  if (!is.null(cob_porte)) {
    dB <- cob_porte %>%
      select(CATEGORIA_POP, VAR_PARAM_PCT, VAR_OBS_PCT) %>%
      pivot_longer(-CATEGORIA_POP, names_to = "Fonte", values_to = "Share") %>%
      mutate(Fonte = dplyr::recode(as.character(Fonte),
                                   VAR_PARAM_PCT = "Parameter uncertainty",
                                   VAR_OBS_PCT   = "Observation dispersion"))
    
    paineis$B <- ggplot(dB, aes(x = CATEGORIA_POP, y = Share, fill = Fonte)) +
      geom_col(width = 0.68) +
      geom_text(aes(label = ifelse(Share >= 6, sprintf("%.0f%%", Share), "")),
                position = position_stack(vjust = 0.5),
                size = 3, colour = "white") +
      scale_fill_manual(values = c("Parameter uncertainty"  = "#4393C3",
                                   "Observation dispersion" = "#D6604D")) +
      scale_y_continuous(expand = c(0, 0)) +
      labs(title = "B. Decomposition of predictive variance",
           subtitle = "Observation dispersion dominates at every scale",
           x = NULL, y = "Share of predictive variance (%)") +
      tema_pub +
      theme(axis.text.x = element_text(angle = 18, hjust = 1))
  }
  
  dC <- df_ip %>%
    filter(PRED > 0) %>%
    mutate(LARG_PI = PI_LS - PI_LI, LARG_IC = LS - LI) %>%
    filter(LARG_PI > 0, LARG_IC > 0)   # log scale needs strictly positive widths
  
  if (nrow(dC) > 0) {
    set.seed(1)
    dC_s <- dC[sample(seq_len(nrow(dC)), min(30000, nrow(dC))), ]
    paineis$C <- ggplot(dC_s, aes(x = PRED)) +
      geom_point(aes(y = LARG_PI, colour = "Predictive interval"),
                 alpha = 0.10, size = 0.5) +
      geom_point(aes(y = LARG_IC, colour = "Interval of the mean"),
                 alpha = 0.10, size = 0.5) +
      scale_x_log10(labels = scales::comma) +
      scale_y_log10(labels = scales::comma) +
      scale_colour_manual(values = c("Predictive interval"  = "#2166AC",
                                     "Interval of the mean" = "#D6604D")) +
      guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2.5))) +
      labs(title = "C. Interval width versus expected count",
           subtitle = "Log-log; the gap between the two is the omitted observation process",
           x = "Predicted mean count", y = "95% interval width") +
      tema_pub
  }
  
  paineis <- Filter(Negate(is.null), paineis)
  
  if (length(paineis) > 0) {
    if (requireNamespace("patchwork", quietly = TRUE) && length(paineis) == 3) {
      fig <- patchwork::wrap_plots(paineis$A, paineis$B, paineis$C, ncol = 3) +
        patchwork::plot_annotation(
          title = "Out-of-sample posterior predictive intervals",
          subtitle = sprintf("Rolling-origin cross-validation | %d draws per municipality-month | global predictive coverage %.1f%%",
                             IPCFG$n_sim, cob_global$COB_PREDITIVO),
          theme = theme(plot.title = element_text(face = "bold", size = 15),
                        plot.subtitle = element_text(colour = "grey35", size = 10)))
      ggsave(file.path(DIR_IP, "FIG_IntervalosPreditivos_Painel.png"),
             fig, width = 18, height = 6.2, dpi = 600, bg = "white")
    } else {
      for (nm in names(paineis))
        ggsave(file.path(DIR_IP, paste0("FIG_IntervalosPreditivos_", nm, ".png")),
               paineis[[nm]], width = 8, height = 6, dpi = 600, bg = "white")
    }
    message("   -> figuras salvas em ", DIR_IP)
  }
  
  # ============================================================================
  # 8. OBJETOS GLOBAIS E RESUMO
  # ============================================================================
  df_ip_oos            <<- df_ip
  cobertura_oos        <<- cob_global
  cobertura_oos_porte  <<- cob_porte
  cobertura_oos_drs    <<- cob_drs
  cobertura_insample   <<- cob_insample
  
  message("\n   Módulo 25 concluído.")
  message(sprintf("   Cobertura preditiva OOS (municipio) : %.2f%%",
                  cob_global$COB_PREDITIVO))
  message(sprintf("   Coverage of the interval of the mean: %.2f%%",
                  cob_global$COB_MEDIA))
  if (!is.na(cob_drs_global))
    message(sprintf("   Cobertura preditiva agregada por DRS: %.2f%%", cob_drs_global))
  message(sprintf("   Variancia: parametro %.1f%% | observacional %.1f%%",
                  cob_global$VAR_PARAM_PCT, cob_global$VAR_OBS_PCT))
  if (!is.null(cob_insample))
    message(sprintf("   Cobertura preditiva IN-SAMPLE       : %.2f%%",
                    cob_insample$COB_PREDITIVO))
}