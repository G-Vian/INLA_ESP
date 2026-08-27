# ==============================================================================
# 01 - DATA LOADING AND PANEL CONSTRUCTION
# ------------------------------------------------------------------------------
# Loads the processed cache (FILE_CACHE) and assembles the balanced monthly
# municipality panel used by the model. When the cache is present it is loaded
# directly; the raw-data reconstruction path is not exercised in this
# repository (the cache is provided — see README, "Data").
# Original working comments are kept in Portuguese below.
# ==============================================================================

# ---- required packages -------------------------------------------------------
suppressPackageStartupMessages({
  library(stringr)
  library(lubridate)
  library(tidyr)
  library(INLA)
  library(sf)
  library(spdep)
  library(Matrix)
  library(ISOweek)
  library(dplyr)
})


# ==============================================================================
# 2. PROCESSAMENTO E DADOS
# ==============================================================================
if (file.exists(FILE_CACHE)) {
  message("💾 CACHE ENCONTRADO! Avaliando integridade...")
  load(FILE_CACHE)
  if (exists("dados_base") && "LOG_OFFSET" %in% names(dados_base)) {
    message("✅ Cache aprovado! Carregando dados mensais...")
    if (exists("nb")) {
      graph_file <- file.path(tempdir(), "map.graph")
      nb2INLA(graph_file, nb)
      g <- inla.read.graph(graph_file)
    }
  } else {
    message("⚠️ Cache desatualizado. Reconstruindo...")
    suppressWarnings(rm(dados_base, mapa_sf, nb, g, df_clima_mensal, df_geo_drs, df_dengue_raw))
  }
}

if (!exists("dados_base")) {
  message("⚠️ INICIANDO PROCESSAMENTO COMPLETO...")
  carregar_rdata <- function(caminho) {
    env <- new.env(); load(caminho, envir = env); obj <- env[[ls(env)[1]]]
    if (inherits(obj, "data.table") || is.data.frame(obj)) { cols <- names(obj); for (j in cols) if (inherits(obj[[j]], "IDate")) obj[[j]] <- as.Date(obj[[j]]) }
    return(obj)
  }
  
  df_dengue_raw <- carregar_rdata(PATH_DENGUE) %>% as_tibble()
  df_clima_raw  <- carregar_rdata(PATH_CLIMA) %>% as_tibble() 
  df_social     <- read_csv(PATH_SOCIAL, show_col_types = FALSE)
  df_oni        <- read_csv(PATH_ONI, show_col_types = FALSE) %>% as_tibble()
  df_veg        <- read_csv(PATH_VEG, show_col_types = FALSE) %>% 
    mutate(COD_MUN_6 = as.character(COD_MUN_6)) %>% 
    as_tibble()
  
  message(">> Processando Tabela de Regiões de Saúde (geo_full)...")
  geo_full_raw <- carregar_rdata(PATH_GEO_FULL) %>% as_tibble()
  df_geo_drs <- geo_full_raw %>%
    select(ID_MUN, DRS) %>%
    rename(COD_MUN_FULL = ID_MUN, NOME_DRS = DRS) %>%
    mutate(COD_GEO_6 = str_trim(str_sub(as.character(COD_MUN_FULL), 1, 6)), NOME_DRS = str_trim(as.character(NOME_DRS)), DRS_ID = as.numeric(as.factor(NOME_DRS))) %>% distinct(COD_GEO_6, .keep_all = TRUE)
  
  df_dengue_mensal <- df_dengue_raw %>%
    filter(str_sub(COD_GEO, 1, 2) == CONFIG_ESTADO) %>%
    mutate(DATA_OBJ = ISOweek::ISOweek2date(sprintf("%04d-W%02d-1", ANO_EPI, SEMANA_EPI)), ANO = year(DATA_OBJ), MES = month(DATA_OBJ)) %>%
    group_by(COD_GEO, ANO, MES) %>% summarise(CASOS = sum(CASOS, na.rm = TRUE), .groups = "drop")
  
  df_clima_mensal <- df_clima_raw %>% filter(str_sub(COD_GEO, 1, 2) == CONFIG_ESTADO) %>% mutate(COD_GEO = as.character(COD_GEO), ANO = as.numeric(ANO), MES = as.numeric(MES))
  
  municipios <- unique(df_dengue_mensal$COD_GEO)
  ano_inicio_grid <- min(CONFIG_ANOS) 
  grid_tempo <- expand.grid(COD_GEO = municipios, ANO = ano_inicio_grid:max(CONFIG_ANOS), MES = 1:12) %>% as_tibble()
  
  dados_base <- grid_tempo %>%
    left_join(df_dengue_mensal, by = c("COD_GEO", "ANO", "MES")) %>% mutate(CASOS = replace_na(CASOS, 0)) %>% left_join(df_clima_mensal, by = c("COD_GEO", "ANO", "MES")) %>% left_join(df_oni, by = c("ANO", "MES"))
  
  df_social_pop <- df_social %>% select(COD_MUN_6, POP_TOTAL) %>% mutate(COD_MUN_6 = as.character(COD_MUN_6))
  # Filter to São Paulo only before computing the mean and SD used by scale()
  df_social_norm <- df_social %>% 
    filter(as.character(COD_MUN_6) %in% as.character(str_sub(municipios, 1, 6))) %>%
    select(-POP_TOTAL) %>% 
    mutate(COD_MUN_6 = as.character(COD_MUN_6)) %>% 
    mutate(across(where(is.numeric) & !starts_with("COD_"), ~ as.numeric(scale(.))))
  
  
  dados_base <- dados_base %>%
    mutate(COD_GEO_6 = str_trim(str_sub(as.character(COD_GEO), 1, 6))) %>%
    left_join(df_social_norm, by = c("COD_GEO_6" = "COD_MUN_6")) %>%
    left_join(df_social_pop, by = c("COD_GEO_6" = "COD_MUN_6")) %>%
    left_join(df_geo_drs, by = "COD_GEO_6") %>%
    left_join(df_veg, by = c("COD_GEO_6" = "COD_MUN_6", "ANO" = "Ano", "MES" = "Mes")) %>%
    group_by(ANO, MES) %>% 
    mutate(INDICE_VEGETACAO = replace_na(INDICE_VEGETACAO, mean(INDICE_VEGETACAO, na.rm = TRUE))) %>%
    ungroup() %>%
    mutate(INDICE_VEGETACAO = replace_na(INDICE_VEGETACAO, mean(INDICE_VEGETACAO, na.rm = TRUE))) %>%
    mutate(POP_TOTAL = replace_na(POP_TOTAL, 1), POP_TOTAL = ifelse(POP_TOTAL <= 0, 1, POP_TOTAL), LOG_OFFSET = log(POP_TOTAL / 100000), CATEGORIA_POP_ID = case_when(POP_TOTAL < 20000 ~ 1, POP_TOTAL < 100000 ~ 2, POP_TOTAL < 500000 ~ 3, TRUE ~ 4), DRS_ID = replace_na(DRS_ID, 1))
  
  mapa_sf <- read_sf(PATH_SHAPE) %>% mutate(CD_MUN = as.character(CD_MUN), code_join = str_sub(CD_MUN, 1, 6), ID_AREA = 1:n()) %>% st_make_valid()
  dados_base <- dados_base %>% left_join(mapa_sf %>% st_drop_geometry() %>% select(code_join, ID_AREA), by = c("COD_GEO_6" = "code_join")) %>% mutate(DATA = make_date(ANO, MES, 1), ID_TEMPO = as.numeric(as.factor(DATA)), MES_CICLICO = MES, ANO_ID = as.numeric(as.factor(ANO)), ID_AREA_ANO = as.numeric(as.factor(paste0(COD_GEO, "_", ANO)))) %>% arrange(COD_GEO, DATA)
  
  nb <- poly2nb(mapa_sf); n_comp <- n.comp.nb(nb)
  if (n_comp$nc > 1) { coords <- st_coordinates(st_centroid(mapa_sf)); knn <- knearneigh(coords, k = 1); nb <- union.nb(nb, knn2nb(knn)) }
  nb <- make.sym.nb(nb)
  graph_file <- file.path(tempdir(), "map.graph")
  nb2INLA(graph_file, nb); g <- inla.read.graph(graph_file)
  
  save(dados_base, mapa_sf, nb, g, df_clima_mensal, df_geo_drs, df_dengue_raw, file = FILE_CACHE)
}

# ==============================================================================
# 2.B  SELECTABLE SPATIAL STRUCTURE (SWITCH_GRAFO_VERSAO)
# ------------------------------------------------------------------------------
# Version 1 -> uses the Besag graph 'g' built above (default behaviour).
# Versions 2/3/4 -> load the matrix C = D - W from an external graph and its constraints,
# para uso com model="generic0" (ponderado) no modulo_03.
#
# CRITICAL SAFETY: row i of matrix C corresponds to municipality i IN THE ORDER
# in which the graph script read the shapefile (id_inla = 1:n). Here ID_AREA is also
# 1:n() in read_sf(PATH_SHAPE) order. If the two shapefiles are the same
# arquivo, a ordem bate. Conferimos isso explicitamente abaixo (stopifnot) e
# we log everything to a dedicated file, since a misalignment would leave the model
# SILENTLY wrong (swapped neighbourhoods).
# ==============================================================================
C_ESPACIAL <- NULL; A_ESPACIAL <- NULL; e_ESPACIAL <- NULL; RANKDEF_ESPACIAL <- NULL

# Dedicated graph log (append: grows on each run)
DIR_LOG_GRAFO <- if (exists("DIR_OUTPUT")) DIR_OUTPUT else getwd()
PATH_LOG_GRAFO <- file.path(DIR_LOG_GRAFO, "LOG_INLA_GRAFO.txt")
log_grafo <- function(...) {
  msg <- paste0(...)
  message(msg)                                   # mostra no console
  cat(msg, "\n", file = PATH_LOG_GRAFO, append = TRUE)  # e grava no arquivo
}

cat("\n", strrep("=", 78), "\n",
    "LOG DO GRAFO ESPACIAL  |  ", as.character(Sys.time()), "\n",
    strrep("=", 78), "\n", file = PATH_LOG_GRAFO, append = TRUE)

.vsel <- as.character(SWITCH_GRAFO_VERSAO)
log_grafo(">> SWITCH_GRAFO_VERSAO = ", .vsel,
          c("1"=" (ORIGINAL: besag binário, grafo interno)",
            "2"=" (GRAFO #1: rodoviário puro, generic0)",
            "3"=" (GRAFO #2: estradas + mobilidade, generic0)",
            "4"=" (GRAFO #3: estradas + mobilidade + atalhos, generic0)",
            "5"=" (GRAFO #4: Grafo #3 + camada aérea, generic0)")[.vsel])
log_grafo(">> SWITCH_GRAFO (liga/desliga efeito espacial) = ", SWITCH_GRAFO)

if (isTRUE(SWITCH_GRAFO) && .vsel != "1") {
  
  if (is.null(GRAFOS_ARQUIVOS[[.vsel]]))
    stop("SWITCH_GRAFO_VERSAO=", .vsel, " inválido (use 1, 2, 3, 4 ou 5).")
  
  arq <- GRAFOS_ARQUIVOS[[.vsel]]
  log_grafo("Carregando matriz do grafo a partir de:")
  log_grafo("   C  : ", arq$C)
  log_grafo("   A  : ", arq$A)
  log_grafo("   e  : ", arq$e)
  
  faltam <- unlist(arq)[!file.exists(unlist(arq))]
  if (length(faltam))
    stop("Arquivo(s) do grafo não encontrado(s):\n  ", paste(faltam, collapse = "\n  "),
         "\nConfira DIR_GRAFOS e os nomes em GRAFOS_ARQUIVOS (modulo_00).")
  
  C_ESPACIAL <- readRDS(arq$C)
  A_ESPACIAL <- readRDS(arq$A)
  e_ESPACIAL <- readRDS(arq$e)
  
  # --- dimensions ---
  nC <- nrow(C_ESPACIAL); nMapa <- nrow(mapa_sf)
  log_grafo("Dimensão da matriz C : ", nC, " x ", ncol(C_ESPACIAL))
  log_grafo("Nº de municípios (mapa_sf/ID_AREA): ", nMapa)
  if (nC != nMapa)
    stop("INCOMPATIBILIDADE DE TAMANHO: C tem ", nC, " linhas mas o mapa tem ",
         nMapa, " municípios. Os dois shapefiles precisam ter os mesmos municípios.")
  
  # --- AUTOMATIC ORDER CHECK (hard stopifnot) --------------------------------
  # Row i of matrix C corresponds to municipality i IN THE ORDER the graph script
  # read SP_Municipios_2024.shp (id_inla = row_number()). Here ID_AREA is also
  # 1:n() na mesma leitura. Se for o MESMO arquivo lido da mesma forma, as ordens
  # 1:n(). To PROVE the orders match (rather than trust mapa_sf, which may be cached),
  # we RE-READ the shapefile from disk NOW and compare order municipality by municipality.
  cd_mapa <- as.character(mapa_sf$CD_MUN)
  log_grafo("Verificação de ordem dos municípios (CD_MUN):")
  log_grafo("   1º município (linha 1)  : ", cd_mapa[1],  " - ", mapa_sf$NM_MUN[1])
  log_grafo("   último município (linha ", nMapa, "): ",
            cd_mapa[nMapa], " - ", mapa_sf$NM_MUN[nMapa])
  
  ordem_ok <- FALSE
  if (exists("PATH_SHAPE") && file.exists(PATH_SHAPE)) {
    # Re-read the shapefile from scratch, with the SAME logic as the pipeline:
    # row_number() sobre read_sf, sem reordenar -> id_inla 1:n.
    cd_ref <- tryCatch({
      ref <- sf::read_sf(PATH_SHAPE)
      as.character(ref$CD_MUN)
    }, error = function(e) { log_grafo("   [aviso] não consegui reler PATH_SHAPE: ", conditionMessage(e)); NULL })
    
    if (!is.null(cd_ref)) {
      cd_ref  <- trimws(as.character(cd_ref))
      cd_mapa_cmp <- trimws(as.character(cd_mapa))
      log_grafo("   Releitura de PATH_SHAPE para conferência: ", length(cd_ref), " municípios.")
      # Mesmo tamanho?
      if (length(cd_ref) != length(cd_mapa_cmp)) {
        stop("VERIFICAÇÃO DE ORDEM FALHOU: o shapefile relido tem ", length(cd_ref),
             " municípios, mas mapa_sf (em uso) tem ", length(cd_mapa_cmp),
             ".\nO mapa_sf provavelmente veio de um CACHE com shapefile diferente. ",
             "Apague o FILE_CACHE e rode de novo.")
      }
      # Identical order element by element?
      divergencias <- which(cd_ref != cd_mapa_cmp)
      if (length(divergencias) > 0) {
        i1 <- divergencias[1]
        stop("VERIFICAÇÃO DE ORDEM FALHOU: a ordem de mapa_sf (em uso) NÃO bate com ",
             "a ordem do shapefile em PATH_SHAPE.\n",
             "Primeira divergência na posição ", i1, ": mapa_sf tem CD_MUN=", cd_mapa_cmp[i1],
             " mas o shapefile tem CD_MUN=", cd_ref[i1], ".\n",
             "Total de posições divergentes: ", length(divergencias), " de ", length(cd_ref), ".\n",
             ">> A matriz C assume a ordem do shapefile; com mapa_sf desalinhado, o ",
             "efeito espacial ficaria trocado. Apague o FILE_CACHE (DADOS_BASE_*.RData) ",
             "e rode de novo para reconstruir mapa_sf a partir do shapefile correto.")
      }
      ordem_ok <- TRUE
      log_grafo("   ✅ ORDEM CONFERIDA: mapa_sf bate 100% com o shapefile (",
                length(cd_ref), "/", length(cd_ref), " municípios na mesma posição).")
      log_grafo("   => A linha i da matriz C corresponde ao município i dos dados. OK.")
    }
  } else {
    log_grafo("   [aviso] PATH_SHAPE não encontrado; conferência cruzada pulada.")
  }
  
  # Hard guard: if re-reading was possible and anything mismatched, it already aborted.
  # If re-reading was not possible, we record that order could not be proven.
  if (!ordem_ok)
    log_grafo("   [ATENÇÃO] Ordem NÃO foi verificada automaticamente (releitura indisponível). ",
              "Confira manualmente que PATH_SHAPE é o MESMO shapefile usado no script do grafo.")
  
  # Matrix sanity check: symmetric and row sums ~ 0 (Laplacian)
  somaLin <- Matrix::rowSums(C_ESPACIAL)
  simetrica <- Matrix::isSymmetric(C_ESPACIAL)
  log_grafo("C é simétrica? ", simetrica)
  log_grafo("Soma das linhas de C (deve ~0 p/ Laplaciano): min=",
            round(min(somaLin), 6), " | max=", round(max(somaLin), 6))
  if (!simetrica) log_grafo("AVISO: C não é exatamente simétrica - verifique o grafo.")
  if (max(abs(somaLin)) > 1e-6)
    log_grafo("AVISO: linhas de C não somam zero - confirme que é C = D - W.")
  
  # rankdef = number of connected components = number of rows of Aconstr
  RANKDEF_ESPACIAL <- nrow(A_ESPACIAL)
  log_grafo("Restrições soma-zero (Aconstr): ", RANKDEF_ESPACIAL, " linha(s)  => rankdef = ", RANKDEF_ESPACIAL)
  log_grafo("Comprimento de econstr: ", length(e_ESPACIAL),
            if (length(e_ESPACIAL) == RANKDEF_ESPACIAL) "  (OK, casa com Aconstr)" else "  (AVISO: não casa com Aconstr!)")
  
  # --- ESCALONAMENTO DA MATRIZ (equivale ao scale.model=TRUE do besag) ---------
  # inla.scale.model() standardises C so that the GEOMETRIC marginal variance of
  # the effect is ~1, RESPECTING the sum-to-zero constraints. Without it, the scale of
  # the precision hyperparameter depends on each graph's geometry/weights, and the
  # same pc.prec prior would mean different things across graphs -> comparison
  # de DIC/WAIC enviesada. Com isto, os 4 grafos ficam na MESMA base.
  if (isTRUE(SWITCH_ESCALAR_GRAFO)) {
    diag_antes <- mean(Matrix::diag(C_ESPACIAL))
    C_ESPACIAL <- INLA::inla.scale.model(
      C_ESPACIAL,
      constr = list(A = A_ESPACIAL, e = e_ESPACIAL)
    )
    diag_depois <- mean(Matrix::diag(C_ESPACIAL))
    log_grafo("Escalonamento inla.scale.model: APLICADO (scale.model equivalente).")
    log_grafo("   diagonal média de C: antes=", round(diag_antes, 4),
              " -> depois=", round(diag_depois, 4),
              "  (mudança de escala é esperada)")
  } else {
    log_grafo("Escalonamento inla.scale.model: DESLIGADO (SWITCH_ESCALAR_GRAFO=FALSE).")
    log_grafo("   ATENÇÃO: sem escalonar, a prior de precisão e o DIC NÃO são",
              " diretamente comparáveis entre grafos de escalas diferentes.")
  }
  
  # Type coercions for INLA (avoids Matrix class warnings)
  C_ESPACIAL <- as(C_ESPACIAL, "TsparseMatrix")
  log_grafo("Matriz do grafo pronta para generic0.")
  
  # --- OVERWRITE 'nb' with THIS graph's neighbourhood ---
  # The weight matrix is W = -C off the diagonal (since C = D - W). The binary
  # adjacency (W != 0) defines the neighbours. NOTE: the scaling above multiplies
  # C's values by constants but PRESERVES the non-zero pattern, so the
  # ADJACENCY (who neighbours whom) is identical - 'nb' is correct either way.
  #
  # 'nb' reflete o grafo selecionado.
  nb_interno_backup <- if (exists("nb")) nb else NULL   # guarda o poly2nb original
  Wadj <- -C_ESPACIAL
  Matrix::diag(Wadj) <- 0
  Aadj <- as(Wadj != 0, "lMatrix")                      # TRUE onde há aresta
  nb <- spdep::mat2listw(Aadj, style = "B", zero.policy = TRUE)$neighbours
  attr(nb, "region.id") <- as.character(seq_len(nrow(C_ESPACIAL)))
  .grau <- spdep::card(nb)
  log_grafo("'nb' redefinido a partir da matriz do grafo (para os diagnósticos):")
  log_grafo("   nº de nós     : ", length(nb))
  log_grafo("   nº de arestas : ", sum(.grau) / 2)
  log_grafo("   grau médio    : ", round(mean(.grau), 2),
            " | grau máximo: ", max(.grau))
  log_grafo("   componentes   : ", spdep::n.comp.nb(nb)$nc, "\n")
  
} else if (isTRUE(SWITCH_GRAFO) && .vsel == "1") {
  log_grafo("Usando estrutura ORIGINAL: efeito 'besag' com o grafo interno (g).")
  log_grafo("Grafo interno: nós=", length(nb), " | componentes=",
            spdep::n.comp.nb(nb)$nc, "\n")
} else {
  log_grafo("SWITCH_GRAFO = FALSE: nenhum efeito espacial será incluído.\n")
}

message("\n--- AUDITORIA DE DADOS ---")
message("Total de Casos na Base: ", sum(dados_base$CASOS[dados_base$ANO %in% CONFIG_ANOS], na.rm=TRUE))
if (SWITCH_OFFSET_POP) {
  message(">> OFFSET POPULACIONAL ATIVADO! População média (cidade/mês): ", round(mean(dados_base$POP_TOTAL, na.rm=TRUE), 0))
}
# ==============================================================================
# 2.5 OPTIONAL ANNUAL/GLOBAL TRANSFORMATION FOR A LINEAR VARIABLE
# ==============================================================================
if ("INDICE_VEGETACAO" %in% names(dados_base)) {
  if (SWITCH_NDVI_MEDIA_TUDO) {
    message("\n>> TRANSFORMAÇÃO NDVI: Usando a Média Histórica Total por Município (Valor Fixo).")
    dados_base <- dados_base %>%
      group_by(COD_GEO) %>%
      mutate(INDICE_VEGETACAO = mean(INDICE_VEGETACAO, na.rm = TRUE)) %>%
      ungroup()
    
    # Safety guard: a fixed value cannot have a lag or DLNM
    SWITCH_DLNM_NDVI <- FALSE
    CONFIG_LAGS$INDICE_VEGETACAO <- 0
    message("   ⚠️ Segurança: DLNM do NDVI foi desligado e Lag = 0 (Evita crash matemático por falta de variância temporal).")
    
  } else if (SWITCH_NDVI_MEDIA_ANUAL) {
    message("\n>> TRANSFORMAÇÃO NDVI: Usando a Média Anual por Município (Valor fixo por ano).")
    dados_base <- dados_base %>%
      group_by(COD_GEO, ANO) %>%
      mutate(INDICE_VEGETACAO = mean(INDICE_VEGETACAO, na.rm = TRUE)) %>%
      ungroup()
    
    # Safety guard: an annual value has no monthly variation for a short lag or DLNM
    SWITCH_DLNM_NDVI <- FALSE
    CONFIG_LAGS$INDICE_VEGETACAO <- 0
    message("   ⚠️ Segurança: DLNM do NDVI foi desligado e Lag = 0 (Evita crash matemático por falta de variância mensal).")
  }
}

# ==============================================================================