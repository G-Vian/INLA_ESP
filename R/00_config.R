# ==============================================================================
# 00 - CONFIGURATION
# ------------------------------------------------------------------------------
# Central configuration for the São Paulo dengue spatio-temporal model.
# Restricted to the settings used by the FINAL model reported in the paper.
#
# Some switches for model variants that are NOT used in the final specification
# are kept and set to FALSE: several modules read them with if(SWITCH_...) and
# would error if the variable were simply removed. They are grouped and marked
# below as "kept FALSE for compatibility".
#
# All input paths are RELATIVE to the repository root (DIR_BASE_CODIGO, set by
# R/main.R). Nothing is hard-coded to an absolute location.
# ==============================================================================

CONFIG_MODELO_TIPO <- "ORIGINAL_RESGATE2024_DRS_V75_AR1_ESTRUTURADO"

# ------------------------------------------------------------------------------
# 1. FINAL-MODEL STRUCTURE SWITCHES
# ------------------------------------------------------------------------------
SWITCH_AUTOREGRESSAO_LAG1 <- TRUE     # log-lag 1 on cases
SWITCH_AUTOREGRESSAO_LAG2 <- TRUE     # log-lag 2 on cases

SWITCH_SAZONALIDADE       <- TRUE     # global cyclic RW2 seasonality
SWITCH_SAZONALIDADE_DRS   <- TRUE     # region-specific (DRS) cyclic seasonality

SWITCH_GRAFO         <- TRUE          # spatial effect on
SWITCH_GRAFO_VERSAO  <- 1             # 1 = Besag ICAR, queen contiguity (paper)
SWITCH_ESCALAR_GRAFO <- TRUE

SWITCH_DADOS_SOCIAIS  <- TRUE
SWITCH_OFFSET_POP     <- TRUE         # log-population offset
SWITCH_USAR_PC_PRIORS <- TRUE
SWITCH_DLNM_CLIMA     <- TRUE

# Validation: rolling-origin is the main design (module 24); module 25 builds
# the out-of-sample predictive intervals. Leave-one-year-out is NOT used.
SWITCH_TESTE_ROLLING_MENSAL      <- TRUE
SWITCH_INTERVALOS_PREDITIVOS_OOS <- TRUE

# --- Model variants NOT in the final specification (kept FALSE for compatibility)
SWITCH_TENDENCIA_TEMPORAL     <- FALSE
SWITCH_EFEITO_ANO             <- FALSE
SWITCH_INTERACAO_ANO          <- FALSE
SWITCH_AR1_LATENTE_GLOBAL     <- FALSE
SWITCH_AR1_LATENTE_CIDADE     <- FALSE
SWITCH_AUTOREGRESSAO_ESPACIAL <- FALSE
SWITCH_IMUNIDADE_LOCAL        <- FALSE
SWITCH_OLRE_PEQUENAS          <- FALSE
SWITCH_MODO_RAPIDO            <- TRUE   # simplified integration for speed
SWITCH_DLNM_NDVI              <- FALSE  # vegetation index not used
SWITCH_NDVI_MEDIA_ANUAL       <- FALSE
SWITCH_NDVI_MEDIA_TUDO        <- FALSE

# ------------------------------------------------------------------------------
# 2. PRIORS
#    PC priors are specified as P(sd > U) = alpha for the precisions.
# ------------------------------------------------------------------------------
# Fixed effects and autoregressive log-lags
PRIOR_PREC_FIXED      <- 0.05
PRIOR_PREC_LAGS       <- 0.05
PRIOR_PREC_INTERCEPT  <- 0.1
PRIOR_AR              <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)))

# Zero-inflated negative binomial hyperpriors
PRIOR_ZINB_OVERDISP_U     <- 1
PRIOR_ZINB_OVERDISP_ALPHA <- 0.3
PRIOR_ZINB_PROB           <- list(prior = "gaussian", param = c(0, 0.4))
PRIOR_ZINB_PROB0_MEAN     <- 0      # gaussian prior mean for the zero-inflation prob
PRIOR_ZINB_PROB0_PREC     <- 0.4    # gaussian prior precision

# AR1 correlation (rho) — PC prior P(rho > U) = alpha
PRIOR_AR1_RHO_U           <- 0
PRIOR_AR1_RHO_ALPHA       <- 0.9

# RW1 temporal trend precision (used only if SWITCH_TENDENCIA_TEMPORAL = TRUE)
PRIOR_PREC_RW1_TREND_U     <- 1
PRIOR_PREC_RW1_TREND_ALPHA <- 0.01

# Structured spatial (Besag/ICAR) precision
PRIOR_PREC_BESAG_U     <- 1
PRIOR_PREC_BESAG_ALPHA <- 0.01
pc_prec_param <- list(prior = "pc.prec",
                      param = c(PRIOR_PREC_BESAG_U, PRIOR_PREC_BESAG_ALPHA))
PRIOR_PREC_RW <- list(prior = "pc.prec", param = c(1, 0.5))

# Random-walk / IID precisions read by module 03
PRIOR_PREC_RW2_SEAS_U     <- 1
PRIOR_PREC_RW2_SEAS_ALPHA <- 0.5
PRIOR_PREC_IID_ANO_U      <- 1
PRIOR_PREC_IID_ANO_ALPHA  <- 0.01
PRIOR_PREC_IID_OLRE_U     <- 1
PRIOR_PREC_IID_OLRE_ALPHA <- 0.01

# ------------------------------------------------------------------------------
# 3. DLNM CROSS-BASIS SETTINGS
# ------------------------------------------------------------------------------
CONFIG_DLNM_MAXLAG     <- 3
CONFIG_DLNM_MAXLAG_ONI <- 6           # unused (ONI off); kept for module compatibility
CONFIG_DLNM_FUN_VAR    <- "ns"
CONFIG_DLNM_FUN_LAG    <- "ns"
CONFIG_DLNM_DF_VAR     <- 3
CONFIG_DLNM_DF_LAG     <- 3

# ------------------------------------------------------------------------------
# 4. STUDY WINDOW AND VARIABLE SELECTION
# ------------------------------------------------------------------------------
CONFIG_ANOS    <- 2010:2024
CONFIG_ANOS_CV <- 2010:2024
CONFIG_ESTADO  <- "35"                # São Paulo State (IBGE code)

# Climate exposures entered as DLNM cross-bases
SELECAO_CLIMA <- list(temp_min = TRUE, precip_med = TRUE)

# Socioeconomic indicators available in the processed data. Set each to TRUE to
# include it as a fixed effect, or FALSE to leave it out. The five set to TRUE
# below reproduce the final model reported in the paper; the focal deprivation
# indicator is swapped per specification inside module 03. The full list of
# indicators that can be toggled is:
#   PROP_PRETOS_PARDOS  proportion of Black and mixed-race residents
#   PROP_SEM_LIXO       proportion of households without garbage collection
#   PROP_SEM_AGUA       proportion of households without piped water
#   PROP_SEM_ESGOTO     proportion of households without sewage
#   DENSIDADE_DEMO      population density
#   MEDIA_MORADORES     mean residents per household
#   RENDA_MEDIA         mean household income
#   PROP_POBREZA        proportion in poverty
#   PROP_BAIXA_EDUC     proportion with low education
#   INDICE_GINI         Gini index of income inequality
SELECAO_SOCIAL <- list(
  PROP_PRETOS_PARDOS = TRUE,
  PROP_SEM_LIXO      = FALSE,
  PROP_SEM_AGUA      = TRUE,
  PROP_SEM_ESGOTO    = FALSE,
  DENSIDADE_DEMO     = TRUE,
  MEDIA_MORADORES    = TRUE,
  RENDA_MEDIA        = FALSE,
  PROP_POBREZA       = FALSE,
  PROP_BAIXA_EDUC    = FALSE,
  INDICE_GINI        = TRUE
)

# Non-DLNM lag map read by module 02 (only used for variants that are off here)
CONFIG_LAGS <- list()

# ------------------------------------------------------------------------------
# 5. ROLLING-ORIGIN VALIDATION SETTINGS (module 24)
# ------------------------------------------------------------------------------
CONFIG_ROLLING_BURNIN_MESES <- 36L           # 3-year burn-in before first origin
CONFIG_ROLLING_ANOS_TESTE   <- c(2010, 2024) # target-year range for origins
CONFIG_ROLLING_WARMSTART    <- TRUE
CONFIG_ROLLING_QUANTIL_SURTO <- 0.75         # outbreak-month threshold (quantile)

# ------------------------------------------------------------------------------
# 6. PATHS (relative to the repository root)
# ------------------------------------------------------------------------------
if (!exists("DIR_BASE_CODIGO"))
  DIR_BASE_CODIGO <- normalizePath(getwd())

DIR_DATA   <- file.path(DIR_BASE_CODIGO, "data")
DIR_OUTPUT <- file.path(DIR_BASE_CODIGO, "figures")
DIR_CACHE  <- DIR_DATA

# Processed cache bundling the analysis-ready inputs (see README, "Data").
# Distributed separately because of its size; place it in data/.
FILE_CACHE  <- file.path(DIR_CACHE, "DADOS_BASE_MENSAL_V75.RData")

# Rolling-origin checkpoints
DIR_ROLLING_CKPT <- file.path(DIR_DATA, "rolling_checkpoints")

if (!dir.exists(DIR_OUTPUT))       dir.create(DIR_OUTPUT, recursive = TRUE)
if (!dir.exists(DIR_ROLLING_CKPT)) dir.create(DIR_ROLLING_CKPT, recursive = TRUE)

# ------------------------------------------------------------------------------
# 7. INPUT CHECK
# ------------------------------------------------------------------------------
if (!file.exists(FILE_CACHE)) {
  message("\n", strrep("!", 78))
  message("  PROCESSED CACHE NOT FOUND:")
  message("    ", FILE_CACHE)
  message("  Download it (see README, \"Data\") and place it in data/ before")
  message("  running the pipeline.")
  message(strrep("!", 78), "\n")
}

message("[config] Final-model configuration loaded ",
        "(SWITCH_GRAFO_VERSAO = ", SWITCH_GRAFO_VERSAO, ").")