# ==============================================================================
# 03 - MODEL BUILDING AND FINAL INLA FIT
# ------------------------------------------------------------------------------
# Defines the ZINB likelihood, priors and INLA controls; runs the stepwise
# model-building sequence (DIC/WAIC/log-score) and produces the final fit
# (mod_final, f_vencedora, dados_final_modelo) consumed downstream.
# Original working comments are kept in Portuguese below.
# ==============================================================================

# ---- required packages -------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyr)
  library(tibble)
  library(INLA)
  library(dplyr)
})


# 5. MODEL SELECTION AND METRIC EVOLUTION (STEPWISE)
# ==============================================================================
message("\n--- 5. Seleção de Modelos e Evolução de Métricas (Stepwise) ---")
CONFIG_FAMILIA <- "zeroinflatednbinomial1"

# [REGRA RESTAURADA] Conecta as Priors da Mesa de Som ao Motor do INLA
if (SWITCH_USAR_PC_PRIORS) {
  ctrl_family_params <- list(hyper = list(
    theta1 = list(prior = "pc.prec", param = c(PRIOR_ZINB_OVERDISP_U, PRIOR_ZINB_OVERDISP_ALPHA)), 
    theta2 = list(prior = "gaussian", param = c(PRIOR_ZINB_PROB0_MEAN, PRIOR_ZINB_PROB0_PREC))
  ))
} else {
  ctrl_family_params <- list() 
}

pc_prec_param <- list(prior = "pc.prec", param = c(PRIOR_PREC_BESAG_U, PRIOR_PREC_BESAG_ALPHA))
part_trend_str     <- if (SWITCH_TENDENCIA_TEMPORAL) paste0("f(ID_TEMPO, model='rw1', scale.model=TRUE, hyper=list(prec=list(prior='pc.prec', param=c(", PRIOR_PREC_RW1_TREND_U, ", ", PRIOR_PREC_RW1_TREND_ALPHA, "))))") else ""
part_ano_global    <- if (SWITCH_EFEITO_ANO) paste0("f(ANO_ID, model='iid', hyper=list(prec=list(prior='pc.prec', param=c(", PRIOR_PREC_IID_ANO_U, ", ", PRIOR_PREC_IID_ANO_ALPHA, "))))") else ""
part_interacao_ano <- if (SWITCH_INTERACAO_ANO) paste0("f(ID_AREA_ANO, model='iid', hyper=list(prec=list(prior='pc.prec', param=c(", PRIOR_PREC_IID_ANO_U, ", ", PRIOR_PREC_IID_ANO_ALPHA, "))))") else ""
part_olre_peq      <- if (SWITCH_OLRE_PEQUENAS) paste0("f(ID_OBS_PEQ, model='iid', hyper=list(prec=list(prior='pc.prec', param=c(", PRIOR_PREC_IID_OLRE_U, ", ", PRIOR_PREC_IID_OLRE_ALPHA, "))))") else ""
part_seas_global   <- paste0("f(MES_GLOBAL, model='rw2', cyclic=TRUE, scale.model=TRUE, hyper=list(prec=list(prior='pc.prec', param=c(", PRIOR_PREC_RW2_SEAS_U, ", ", PRIOR_PREC_RW2_SEAS_ALPHA, "))))")
part_seas_drs      <- paste0("f(MES_DRS, model='rw2', cyclic=TRUE, scale.model=TRUE, replicate=DRS_ID, hyper=list(prec=list(prior='pc.prec', param=c(", PRIOR_PREC_RW2_SEAS_U, ", ", PRIOR_PREC_RW2_SEAS_ALPHA, "))))")

part_ar_temp_comps <- c()
if (SWITCH_AUTOREGRESSAO_LAG1) part_ar_temp_comps <- c(part_ar_temp_comps, paste0("f(CASOS_LAG1_LOG, model='linear', mean.linear=0, prec.linear=", PRIOR_PREC_LAGS, ")"))
if (SWITCH_AUTOREGRESSAO_LAG2) part_ar_temp_comps <- c(part_ar_temp_comps, paste0("f(CASOS_LAG2_LOG, model='linear', mean.linear=0, prec.linear=", PRIOR_PREC_LAGS, ")"))
if (SWITCH_AR1_LATENTE_GLOBAL) part_ar_temp_comps <- c(part_ar_temp_comps, paste0("f(ID_TEMPO, model='ar1', hyper=list(rho=list(prior='pc.cor1', param=c(", PRIOR_AR1_RHO_U, ", ", PRIOR_AR1_RHO_ALPHA, ")), prec=list(prior='pc.prec', param=c(1, 0.01))))"))
if (SWITCH_AR1_LATENTE_CIDADE) part_ar_temp_comps <- c(part_ar_temp_comps, paste0("f(ID_TEMPO, model='ar1', replicate=ID_AREA, hyper=list(rho=list(prior='pc.cor1', param=c(", PRIOR_AR1_RHO_U, ", ", PRIOR_AR1_RHO_ALPHA, ")), prec=list(prior='pc.prec', param=c(1, 0.01))))"))
part_ar_temp <- if (length(part_ar_temp_comps) > 0) paste(part_ar_temp_comps, collapse = " + ") else ""

part_imun_str      <- if (SWITCH_IMUNIDADE_LOCAL) "f(inla.group(IMUNIDADE_LOCAL, n=20), model='rw1')" else ""
part_ar_esp        <- if (SWITCH_AUTOREGRESSAO_ESPACIAL) "PRESSAO_ESPACIAL_LAG1_LOG" else ""
# --- EFEITO ESPACIAL: besag (original) OU generic0 ponderado (Grafos #1/#2/#3) ---
# A escolha vem de SWITCH_GRAFO_VERSAO (modulo_00); os objetos C_ESPACIAL/
# A_ESPACIAL/e_ESPACIAL/RANKDEF_ESPACIAL are prepared in module 01.
if (!SWITCH_GRAFO) {
  part_besag <- ""
} else if (as.character(SWITCH_GRAFO_VERSAO) == "1") {
  # Besag (ICAR) with the internal graph 'g'
  part_besag <- paste0("f(ID_AREA, model='besag', graph=g, scale.model=TRUE, constr=TRUE, hyper=list(prec=pc_prec_param))")
} else {
  # Alternative graphs: generic0 with a weighted C matrix and sum-to-zero constraints.
  # Cmatrix, extraconstr and rankdef are passed as names of global objects
  # o INLA resolve em tempo de ajuste (C_ESPACIAL etc. existem no ambiente global).
  part_besag <- paste0(
    "f(ID_AREA, model='generic0', Cmatrix=C_ESPACIAL, constr=FALSE, ",
    "extraconstr=list(A=A_ESPACIAL, e=e_ESPACIAL), rankdef=RANKDEF_ESPACIAL, ",
    "hyper=list(prec=pc_prec_param))"
  )
}

# HYBRID FORMULA ASSEMBLY (joins the DLNM bases and the linear covariates)
vars_clima_on <- c()
if (length(vars_clima_dlnm) > 0) vars_clima_on <- c(vars_clima_on, grep("^basis_", names(dados_final_modelo), value = TRUE))
if (length(vars_clima_simples) > 0) vars_clima_on <- c(vars_clima_on, paste0("LAG_", vars_clima_simples))

if (SWITCH_DADOS_SOCIAIS) vars_social_on <- intersect(names(dados_final_modelo), names(SELECAO_SOCIAL)[unlist(SELECAO_SOCIAL)]) else vars_social_on <- character(0)

part_clim <- if (length(vars_clima_on) > 0) paste(vars_clima_on, collapse = " + ") else ""
part_soc  <- if (length(vars_social_on) > 0) paste(vars_social_on, collapse = " + ") else ""
part_offset <- if (SWITCH_OFFSET_POP) "offset(LOG_OFFSET)" else ""

montar_formula <- function(base, comps) {
  valid <- comps[comps != ""]
  offset_str <- if (part_offset != "") paste("+", part_offset) else ""
  f <- as.formula(paste(base, "~ 1", offset_str, if (length(valid) > 0) paste("+", paste(valid, collapse = " + ")) else ""))
  environment(f) <- globalenv()   # garante que o INLA ache C_ESPACIAL, g, etc.
  f
}

lista_formulas <- list()
comps_atuais <- c()
lista_formulas[["M00: Base (Intercepto)"]] = montar_formula("CASOS", comps_atuais)

if (SWITCH_TENDENCIA_TEMPORAL) { comps_atuais <- c(comps_atuais, part_trend_str); lista_formulas[["M01: +Tendência Temporal"]] = montar_formula("CASOS", comps_atuais) }
comps_ruido <- c()
if (SWITCH_EFEITO_ANO) comps_ruido <- c(comps_ruido, part_ano_global)
if (SWITCH_INTERACAO_ANO) comps_ruido <- c(comps_ruido, part_interacao_ano)
if (SWITCH_OLRE_PEQUENAS) comps_ruido <- c(comps_ruido, part_olre_peq)
if (length(comps_ruido) > 0) { comps_atuais <- c(comps_atuais, comps_ruido); lista_formulas[["M02: +Efeitos Aleatórios (Ano/OLRE)"]] = montar_formula("CASOS", comps_atuais) }
if (SWITCH_SAZONALIDADE) {
  comps_atuais <- c(comps_atuais, part_seas_global); lista_formulas[["M03: +Sazonalidade Global"]] = montar_formula("CASOS", comps_atuais)
  if (SWITCH_SAZONALIDADE_DRS) { comps_atuais <- c(comps_atuais, part_seas_drs); lista_formulas[["M04: +Sazonalidade DRS"]] = montar_formula("CASOS", comps_atuais) }
}
if (part_ar_temp != "") { comps_atuais <- c(comps_atuais, part_ar_temp); lista_formulas[["M05: +AR Temporal"]] = montar_formula("CASOS", comps_atuais) }
if (SWITCH_IMUNIDADE_LOCAL) { comps_atuais <- c(comps_atuais, part_imun_str); lista_formulas[["M06: +Imunidade Local"]] = montar_formula("CASOS", comps_atuais) }
if (SWITCH_AUTOREGRESSAO_ESPACIAL) { comps_atuais <- c(comps_atuais, part_ar_esp); lista_formulas[["M07: +AR Espacial (Cluster)"]] = montar_formula("CASOS", comps_atuais) }
if (length(vars_clima_on) > 0) { comps_atuais <- c(comps_atuais, part_clim); lista_formulas[["M08: +Clima (DLNM / Linear)"]] = montar_formula("CASOS", comps_atuais) }
if (length(vars_social_on) > 0) { comps_atuais <- c(comps_atuais, part_soc); lista_formulas[["M09: +Sociais"]] = montar_formula("CASOS", comps_atuais) }
if (part_besag != "") { comps_atuais <- c(comps_atuais, part_besag)
  .rotulo_esp <- if (as.character(SWITCH_GRAFO_VERSAO) == "1") "Besag original" else paste0("generic0 Grafo v", SWITCH_GRAFO_VERSAO)
  lista_formulas[[paste0("M10: +Espaço (", .rotulo_esp, ")")]] = montar_formula("CASOS", comps_atuais) }

tabela_res <- data.frame()
mod_final <- NULL

if (SWITCH_MODO_RAPIDO) {
  ctrl_inla_params <- list(strategy = "gaussian", int.strategy = "eb", huge = TRUE, diagonal = 0.001)
} else {
  ctrl_inla_params <- list(strategy = "simplified.laplace", int.strategy = "ccd", huge = TRUE, tolerance = 0.01, diagonal = 0.001)
}

prec_list <- list(default = PRIOR_PREC_FIXED)
if (SWITCH_AUTOREGRESSAO_LAG1) prec_list[["CASOS_LAG1_LOG"]] <- PRIOR_PREC_LAGS
if (SWITCH_AUTOREGRESSAO_LAG2) prec_list[["CASOS_LAG2_LOG"]] <- PRIOR_PREC_LAGS

ctrl_fixed_params <- list(mean.intercept = 0, prec.intercept = PRIOR_PREC_INTERCEPT, mean = list(default = 0), prec = prec_list)

message("   -> Calculando linha de base (Modelo Nulo) para o R2...")
mod_null <- tryCatch({
  inla(lista_formulas[["M00: Base (Intercepto)"]], family = CONFIG_FAMILIA, data = dados_final_modelo, control.family = ctrl_family_params, control.compute = list(mlik = T, cpo = F), control.inla = ctrl_inla_params, control.fixed = ctrl_fixed_params)
}, error = function(e) NULL)

if (is.null(mod_null)) {
  mod_null <- inla(lista_formulas[["M00: Base (Intercepto)"]], family = "poisson", data = dados_final_modelo, control.compute = list(mlik = T, cpo = F), control.inla = ctrl_inla_params, control.fixed = ctrl_fixed_params)
}
logL0 <- mod_null$mlik[1]

for (n in names(lista_formulas)) {
  message("🚀 Rodando: ", n)
  is_complex <- grepl("Espaço|AR Temporal", n)
  ctrl_local <- ctrl_inla_params
  if (is_complex) ctrl_local$diagonal <- 0.1 
  
  m <- tryCatch({
    inla(lista_formulas[[n]], family = CONFIG_FAMILIA, data = dados_final_modelo, control.family = ctrl_family_params, control.inla = ctrl_local, control.fixed = ctrl_fixed_params, control.compute = list(dic = T, waic = T, cpo = T, mlik = T))
  }, error = function(e) NULL)
  
  if (!is.null(m)) {
    r2 <- 1 - exp(-(2 / nrow(dados_final_modelo)) * (m$mlik[1] - logL0))
    ls_val <- NA
    if (!is.null(m$cpo$cpo)) { cpo_vals <- m$cpo$cpo; ls_val <- sum(-log(cpo_vals[cpo_vals > 0 & !is.na(cpo_vals)])) }
    tabela_res <- bind_rows(tabela_res, data.frame(Modelo = n, DIC = m$dic$dic, WAIC = m$waic$waic, LogScore = ls_val, R2_LR = r2))
    message("✅ Sucesso.")
  }
}

tabela_res <- tabela_res %>% mutate(Ganho_DIC = DIC - lag(DIC), Ganho_WAIC = WAIC - lag(WAIC), Ganho_LogScore = LogScore - lag(LogScore), Ganho_R2 = R2_LR - lag(R2_LR))

if (nrow(tabela_res) > 1) {
  df_plot_melt <- tabela_res %>% select(Modelo, DIC, WAIC, LogScore, R2_LR) %>% pivot_longer(cols = c(DIC, WAIC, LogScore, R2_LR), names_to = "Metrica", values_to = "Valor")
  # Panel labels
  df_plot_melt$Metrica <- factor(df_plot_melt$Metrica, levels = c("WAIC", "DIC", "LogScore", "R2_LR"), labels = c("WAIC (Lower is Better)", "DIC (Lower is Better)", "Log-Score (Lower is Better)", "Pseudo-R² (Higher is Better)"))
  df_plot_melt$Modelo <- factor(df_plot_melt$Modelo, levels = tabela_res$Modelo)
  
  g_fit <- ggplot(df_plot_melt, aes(x = Modelo, y = Valor, group = 1)) + geom_line(color = "#2c3e50", linewidth = 1) + geom_point(color = "#e74c3c", size = 3) + facet_wrap( ~ Metrica, scales = "free_y", ncol = 2) + theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), strip.text = element_text(face = "bold", size = 10)) + 
    # Title and axis labels
    labs(title = "Model Evolution: DIC, WAIC, Log-Score and R²", subtitle = "Impact of adding each structural component sequentially", x = "Model Step", y = "Metric Value")
  
  ggsave(file.path(DIR_OUTPUT, "PLOT_Model_Evolution.png"), g_fit, width = 12, height = 8)
}

message("🚀 Recalculando modelo final com métricas completas...")
f_vencedora <- lista_formulas[[length(lista_formulas)]]
ctrl_final <- ctrl_inla_params
if (grepl("Espaço|AR Temporal", names(lista_formulas)[length(lista_formulas)])) ctrl_final$diagonal <- 0.1
mod_final <- inla(f_vencedora, family = CONFIG_FAMILIA, data = dados_final_modelo, control.family = ctrl_family_params, control.inla = ctrl_final, control.fixed = ctrl_fixed_params, control.compute = list(dic = T, waic = T, cpo = T, mlik = T, config = TRUE))

# ==============================================================================
