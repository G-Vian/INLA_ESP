# ==============================================================================
# 02 - FEATURE ENGINEERING: DLNM CROSS-BASES AND AUTOREGRESSIVE LAGS
# ------------------------------------------------------------------------------
# Builds the distributed-lag non-linear cross-bases for minimum temperature and
# mean precipitation, and the log-lagged case predictors. Runs the dataset
# audit before the model is fitted.
# Original working comments are kept in Portuguese below.
# ==============================================================================

# ---- required packages -------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyr)
  library(dlnm)
  library(dplyr)
  library(RcppRoll)
  library(splines)
})


# 3. ENGENHARIA DE FEATURES
# ==============================================================================
message("\n--- 3. Engenharia de Features ---")
if (exists("nb")) {
  df_edges <- data.frame(ID_AREA = rep(1:length(nb), sapply(nb, length)), VIZINHO_ID = unlist(nb)) %>% filter(VIZINHO_ID > 0)
  self_edges <- data.frame(ID_AREA = 1:length(nb), VIZINHO_ID = 1:length(nb))
  df_edges <- bind_rows(df_edges, self_edges) %>% distinct()
  casos_pop_t <- dados_base %>% select(ID_TEMPO, ID_AREA, CASOS, POP_TOTAL)
  vizinhos_casos <- df_edges %>% left_join(casos_pop_t, by = c("VIZINHO_ID" = "ID_AREA"), relationship = "many-to-many") %>% group_by(ID_AREA, ID_TEMPO) %>% summarise(CASOS_CLUSTER = sum(CASOS, na.rm = TRUE), POP_CLUSTER = sum(POP_TOTAL, na.rm = TRUE), .groups = "drop") %>% mutate(TAXA_ESPACIAL = (CASOS_CLUSTER / POP_CLUSTER) * 100000)
  dados_base <- dados_base %>% left_join(vizinhos_casos, by = c("ID_AREA", "ID_TEMPO")) %>% mutate(TAXA_ESPACIAL = replace_na(TAXA_ESPACIAL, 0))
} else { dados_base$TAXA_ESPACIAL <- 0 }

# --- SPLIT: DLNM CLIMATE VS LINEAR TERMS ---
todas_vars_clima <- names(SELECAO_CLIMA)[unlist(SELECAO_CLIMA)]
vars_clima_dlnm <- character(0)
vars_clima_simples <- character(0)

if (SWITCH_DLNM_CLIMA) {
  vars_clima_dlnm <- todas_vars_clima
  # If a variable's DLNM is OFF, move it to the linear-term list (lag 0)
  if (!SWITCH_DLNM_NDVI && "INDICE_VEGETACAO" %in% vars_clima_dlnm) {
    vars_clima_dlnm <- setdiff(vars_clima_dlnm, "INDICE_VEGETACAO")
    vars_clima_simples <- "INDICE_VEGETACAO"
  }
} else {
  vars_clima_simples <- todas_vars_clima
}
keys_clima_ativos <- c(vars_clima_dlnm, vars_clima_simples)

# 1. Apply the cross-basis only to the DLNM variables
# Global list holding the crossbasis objects (needed later for crosspred)
cb_objects_global <- list()

if (length(vars_clima_dlnm) > 0) {
  dados_final_modelo <- dados_base
  
  for (v in vars_clima_dlnm) {
    serie_completa <- dados_base[[v]]
    if (is.null(serie_completa) || all(is.na(serie_completa))) next
    
    lag_atual <- if (v == "ONI") CONFIG_DLNM_MAXLAG_ONI else CONFIG_DLNM_MAXLAG
    
    # Knots computed ONCE on the global series — fixed across municipalities
    knots_var <- quantile(serie_completa, probs = c(0.10, 0.50, 0.90), na.rm = TRUE)
    
    # Per-municipality cross-basis (correct lags) with global knots (consistent basis)
    resultado_cb <- dados_base %>%
      arrange(COD_GEO, DATA) %>%
      group_by(COD_GEO) %>%
      do({
        sub <- .
        serie_mun <- sub[[v]]
        # Imputa NAs dos primeiros meses (Jan/Fev/Mar 2010 sem clima de lookback)
        # These months are dropped later by drop_na — imputation does not affect the model
        serie_mun[is.na(serie_mun)] <- mean(serie_completa, na.rm = TRUE)
        serie_mun[!is.finite(serie_mun)] <- mean(serie_completa, na.rm = TRUE)
        cb_mun <- dlnm::crossbasis(
          serie_mun,
          lag    = lag_atual,
          argvar = list(fun            = CONFIG_DLNM_FUN_VAR,
                        df             = CONFIG_DLNM_DF_VAR,
                        knots          = knots_var,
                        Boundary.knots = range(serie_completa, na.rm = TRUE)),
          arglag = list(fun = CONFIG_DLNM_FUN_LAG,
                        df  = CONFIG_DLNM_DF_LAG)
        )
        df_cb <- as.data.frame(cb_mun)
        names(df_cb) <- paste0("basis_", v, ".v", 1:ncol(df_cb))
        bind_cols(sub %>% select(COD_GEO, DATA), df_cb)
      }) %>%
      ungroup()
    
    dados_final_modelo <- dados_final_modelo %>%
      left_join(resultado_cb, by = c("COD_GEO", "DATA"))
    
    
    # Store a reference object using the FULL series (all municipalities).
    # Reason: the prediction grid uses the min/max of the full dataset.
    # If cb_objects_global held a single municipality, crosspred would fail
    # silently when the grid extrapolated beyond that municipality's range.
    # Impute NAs on the full series before building the reference object
    serie_completa_imputada <- serie_completa
    media_v <- mean(serie_completa, na.rm = TRUE)
    serie_completa_imputada[is.na(serie_completa_imputada)] <- media_v
    serie_completa_imputada[!is.finite(serie_completa_imputada)] <- media_v
    
    cb_objects_global[[v]] <- dlnm::crossbasis(
      serie_completa_imputada,
      lag    = lag_atual,
      argvar = list(fun            = CONFIG_DLNM_FUN_VAR,
                    df             = CONFIG_DLNM_DF_VAR,
                    knots          = knots_var,
                    Boundary.knots = range(serie_completa, na.rm = TRUE)),
      arglag = list(fun = CONFIG_DLNM_FUN_LAG,
                    df  = CONFIG_DLNM_DF_LAG)
    )
  }
} else {
  dados_final_modelo <- dados_base
}
# 2. Apply a linear lag to the simple (non-DLNM) variables
if (length(vars_clima_simples) > 0) {
  dados_final_modelo <- dados_final_modelo %>% group_by(COD_GEO) %>% mutate(across(all_of(vars_clima_simples), ~ dplyr::lag(., CONFIG_LAGS[[cur_column()]]), .names = "LAG_{.col}")) %>% ungroup()
}
dados_final_modelo <- dados_final_modelo %>%
  group_by(COD_GEO) %>% arrange(DATA) %>%
  mutate(
    TAXA_CASOS = (CASOS / POP_TOTAL) * 100000,
    CASOS_LAG1_LOG = log(dplyr::lag(TAXA_CASOS, 1) + 1),
    CASOS_LAG2_LOG = log(dplyr::lag(TAXA_CASOS, 2) + 1),
    CASOS_LAG3_LOG = log(dplyr::lag(TAXA_CASOS, 3) + 1),
    CASOS_LAG4_LOG = log(dplyr::lag(TAXA_CASOS, 4) + 1),
    PRESSAO_ESPACIAL_LAG1 = dplyr::lag(TAXA_ESPACIAL, 1),
    PRESSAO_ESPACIAL_LAG1_LOG = log(PRESSAO_ESPACIAL_LAG1 + 1),
    SOMA_12M = RcppRoll::roll_sum(TAXA_CASOS, n = 12, align = "right", fill = NA, na.rm = TRUE),
    SURTO_ANTERIOR = dplyr::lag(SOMA_12M, 12),
    IMUNIDADE_LOCAL = ifelse(is.na(as.numeric(scale(SURTO_ANTERIOR))), 0, as.numeric(scale(SURTO_ANTERIOR))),
    MES_EPI = ifelse(MES >= 9, MES - 8, MES + 4), MES_GLOBAL = MES_EPI, MES_DRS = MES_EPI
  ) %>% 
  ungroup() %>% 
  filter(ANO %in% CONFIG_ANOS) %>%
  # ==============================================================================
# START FILTER (drops Jan and Feb of the first year)
# Ensures the model only predicts from March 2010 onward,
# utilizando Janeiro e Fevereiro de 2010 como os Lags 1 e 2.
# ==============================================================================
filter(!(ANO == min(CONFIG_ANOS) & MES %in% 1:CONFIG_DLNM_MAXLAG))


# ==============================================================================
# GUARD AGAINST NAs IN THE COVARIATES
# CASOS_LAG1_LOG and CASOS_LAG2_LOG are added to the cleaning list.
# INLA does not accept empty predictors.
# ==============================================================================
colunas_para_limpar <- c("CASOS", "POP_TOTAL", "IMUNIDADE_LOCAL", "PRESSAO_ESPACIAL_LAG1_LOG", "CASOS_LAG1_LOG", "CASOS_LAG2_LOG")

if (length(vars_clima_dlnm) > 0) colunas_para_limpar <- c(colunas_para_limpar, grep("^basis_", names(dados_final_modelo), value = TRUE))
if (length(vars_clima_simples) > 0) colunas_para_limpar <- c(colunas_para_limpar, paste0("LAG_", vars_clima_simples))

dados_final_modelo <- dados_final_modelo %>% 
  drop_na(any_of(colunas_para_limpar)) %>% 
  filter(if_all(where(is.numeric), ~ !is.infinite(.)))

dados_final_modelo$ID_AREA_TEMPO <- 1:nrow(dados_final_modelo)
if (SWITCH_OLRE_PEQUENAS) dados_final_modelo$ID_OBS_PEQ <- ifelse(dados_final_modelo$POP_TOTAL < 20000, 1:nrow(dados_final_modelo), NA)

# Alert for "ghost" municipalities misaligned in the graph
municipios_nos_dados <- unique(dados_final_modelo$ID_AREA)
municipios_no_grafo  <- 1:g$n
municipios_faltando  <- setdiff(municipios_no_grafo, municipios_nos_dados)
if (length(municipios_faltando) > 0) {
  warning(paste("⚠️ ALERTA:", length(municipios_faltando), "municípios no grafo ficaram sem dados após o drop_na! O INLA fará a interpolação espacial (smoothing) para esses nós fantasmas."))
}

source("/home/g.vian/Pesquisa_Epidemic/Estado_SP_Dengue_INLA_2/INLA_CLIMA_SOCIAIS_MENSAL/modulos/modulo_testes_dataset.R", encoding = "UTF-8")

# ==============================================================================
