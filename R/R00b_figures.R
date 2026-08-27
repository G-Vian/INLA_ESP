# ==============================================================================
# R00b - PUBLICATION FIGURE THEME
# ------------------------------------------------------------------------------
# Sets the shared ggplot2 theme, palette and save defaults for publication
# figures.
# Original working comments are kept in Portuguese below.
# ==============================================================================

# ==============================================================================
if (!isTRUE(get0("CONFIG_FIG_PUBLICACAO", ifnotfound = TRUE))) {
  message("[figuras] CONFIG_FIG_PUBLICACAO = FALSE — tema antigo mantido.")
} else {

  suppressPackageStartupMessages(library(ggplot2))

  # ============================================================================
  # 1. PARAMETROS
  # ============================================================================
  # ESCALA controla o quanto o texto cresce em relacao a figura.
  #   1.00 = no scaling
  #   1.35 = recomendado (o que usamos na versao revisada)
  #   1.60 = para figuras que vao entrar em COLUNA UNICA
  FIG_ESCALA <- as.numeric(get0("CONFIG_FIG_ESCALA", ifnotfound = 1.35))

  # Largura final em polegadas. 7.0 = pagina inteira da Acta Tropica.
  # 3.4 = coluna unica. Desenhar aqui e o que faz a fonte ficar certa.
  FIG_LARGURA_PADRAO <- as.numeric(get0("CONFIG_FIG_LARGURA", ifnotfound = 7.5))
  FIG_ALTURA_PADRAO  <- as.numeric(get0("CONFIG_FIG_ALTURA",  ifnotfound = 5.0))
  FIG_DPI            <- as.numeric(get0("CONFIG_FIG_DPI",     ifnotfound = 600))

  BASE <- 11 * FIG_ESCALA   # ~14.9 pt

  message(sprintf(paste0("[figuras] modo publicacao ON | escala %.2f | ",
                         "fonte base %.1f pt | tela padrao %.1f x %.1f pol | %d dpi"),
                  FIG_ESCALA, BASE, FIG_LARGURA_PADRAO, FIG_ALTURA_PADRAO,
                  as.integer(FIG_DPI)))

  # ============================================================================
  # 2. PALETA SEGURA PARA DALTONISMO (Okabe-Ito)
  # ============================================================================
  # Substitui a paleta anterior mantendo os MESMOS nomes, para nao quebrar
  # nenhuma chamada a REV_CORES$... espalhada pelos modulos.
  REV_CORES <<- list(
    neutro   = "#000000",
    azul     = "#0072B2",
    azul_cl  = "#56B4E9",
    vermelho = "#D55E00",
    verm_cl  = "#E69F00",
    laranja  = "#E69F00",
    verde    = "#009E73",
    cinza    = "#666666",
    roxo     = "#CC79A7"
  )

  REV_PALETA <<- unname(unlist(REV_CORES[c("azul", "vermelho", "verde",
                                           "roxo", "azul_cl", "verm_cl",
                                           "cinza", "neutro")]))

  # ============================================================================
  # 3. TEMA
  # ============================================================================
  # Tamanhos ABSOLUTOS, nao rel(): rel() se perde quando um painel e montado
  # com patchwork e o base_size do conjunto muda.
  rev_tema <<- function(base_size = BASE) {
    ggplot2::theme_bw(base_size = base_size) +
      ggplot2::theme(
        panel.grid.minor  = ggplot2::element_blank(),
        panel.grid.major  = ggplot2::element_line(colour = "grey90",
                                                  linewidth = 0.35),
        panel.border      = ggplot2::element_rect(colour = "black", fill = NA,
                                                  linewidth = 0.9),
        axis.text         = ggplot2::element_text(colour = "black",
                                                  size = base_size * 0.90),
        axis.title        = ggplot2::element_text(face = "bold",
                                                  size = base_size * 1.00),
        axis.ticks        = ggplot2::element_line(colour = "black",
                                                  linewidth = 0.55),
        axis.ticks.length = ggplot2::unit(3.5, "pt"),
        strip.text        = ggplot2::element_text(face = "bold",
                                                  size = base_size * 0.95),
        strip.background  = ggplot2::element_rect(fill = "grey94",
                                                  colour = "black",
                                                  linewidth = 0.6),
        plot.title        = ggplot2::element_text(face = "bold",
                                                  size = base_size * 1.15),
        plot.subtitle     = ggplot2::element_text(colour = "grey30",
                                                  size = base_size * 0.88),
        plot.caption      = ggplot2::element_text(colour = "grey30",
                                                  size = base_size * 0.78),
        legend.position   = "bottom",
        legend.text       = ggplot2::element_text(size = base_size * 0.90),
        legend.title      = ggplot2::element_text(size = base_size * 0.95,
                                                  face = "bold"),
        legend.key.size   = ggplot2::unit(1.05, "lines"),
        legend.background = ggplot2::element_blank(),
        plot.background   = ggplot2::element_rect(fill = "white", colour = NA),
        plot.margin       = ggplot2::margin(9, 11, 9, 9)
      )
  }

  # Vale tambem para os modulos PRINCIPAIS, que nao chamam rev_tema().
  ggplot2::theme_set(rev_tema())

  # ============================================================================
  # 4. GEOMS MAIS GROSSOS
  # ============================================================================
  # Linha de 0.5 e ponto de 1.5 desaparecem quando a figura e reduzida.
  suppressWarnings({
    try(ggplot2::update_geom_defaults("line",      list(linewidth = 0.95)), silent = TRUE)
    try(ggplot2::update_geom_defaults("path",      list(linewidth = 0.95)), silent = TRUE)
    try(ggplot2::update_geom_defaults("point",     list(size = 2.4)),       silent = TRUE)
    try(ggplot2::update_geom_defaults("errorbar",  list(linewidth = 0.75)), silent = TRUE)
    try(ggplot2::update_geom_defaults("errorbarh", list(linewidth = 0.75)), silent = TRUE)
    try(ggplot2::update_geom_defaults("hline",     list(linewidth = 0.75)), silent = TRUE)
    try(ggplot2::update_geom_defaults("vline",     list(linewidth = 0.75)), silent = TRUE)
    try(ggplot2::update_geom_defaults("text",      list(size = BASE / 3.2)), silent = TRUE)
    try(ggplot2::update_geom_defaults("label",     list(size = BASE / 3.2)), silent = TRUE)
  })

  # ============================================================================
  # 5. SALVAMENTO
  # ============================================================================
  # A tela ENCOLHE (7.5x5 no lugar de 10x6). E contraintuitivo, mas e o que
  # aumenta a fonte efetiva na pagina impressa. Chamadas antigas que passam
  # largura/altura explicitos sao reescaladas na mesma proporcao, para nao
  # ter que editar 35 modulos.
  .fig_ajusta <- function(largura, altura) {
    # Aspect ratio for two-panel layouts.
    # algo diferente disso.
    esc <- FIG_LARGURA_PADRAO / 10
    c(largura * esc, altura * esc)
  }

  rev_salva_fig <<- function(plot, arquivo, subdir, largura = 10, altura = 6,
                             dpi = FIG_DPI) {
    d   <- .fig_ajusta(largura, altura)
    cam <- file.path(rev_dir_saida(subdir), arquivo)
    tryCatch({
      ggplot2::ggsave(cam, plot, width = d[1], height = d[2],
                      dpi = dpi, bg = "white", limitsize = FALSE)
      rev_info(sprintf("[fig ] %s  (%.1f x %.1f pol, %d dpi)",
                       arquivo, d[1], d[2], as.integer(dpi)))
    }, error = function(e)
      rev_aviso("falha ao salvar ", arquivo, ": ", e$message))
    invisible(cam)
  }

  rev_painel <<- function(lista_plots, arquivo, subdir, ncol = 2,
                          titulo = NULL, subtitulo = NULL,
                          largura = 15, altura = 11, dpi = FIG_DPI) {
    lista_plots <- Filter(Negate(is.null), lista_plots)
    if (length(lista_plots) == 0) return(invisible(NULL))

    if (requireNamespace("patchwork", quietly = TRUE) && length(lista_plots) > 1) {
      fig <- patchwork::wrap_plots(lista_plots, ncol = ncol) +
        patchwork::plot_annotation(
          title = titulo, subtitle = subtitulo,
          tag_levels = NULL,
          theme = ggplot2::theme(
            plot.title    = ggplot2::element_text(face = "bold",
                                                  size = BASE * 1.25),
            plot.subtitle = ggplot2::element_text(colour = "grey30",
                                                  size = BASE * 0.90)))
      rev_salva_fig(fig, arquivo, subdir, largura, altura, dpi)
    } else {
      for (i in seq_along(lista_plots))
        rev_salva_fig(lista_plots[[i]],
                      sub("\\.png$", paste0("_p", i, ".png"), arquivo),
                      subdir, largura / ncol,
                      altura / max(1, ceiling(length(lista_plots) / ncol)),
                      dpi)
    }
  }

  # ============================================================================
  # 6. ESCALAS DE COR PADRAO
  # ============================================================================
  # Only override the default for DISCRETE scales. Continuous scales
  # (mapas, heatmaps) ficam como estao: viridis e derivados ja sao seguros.
  options(
    ggplot2.discrete.colour = function(...)
      ggplot2::scale_colour_manual(..., values = REV_PALETA),
    ggplot2.discrete.fill = function(...)
      ggplot2::scale_fill_manual(..., values = REV_PALETA)
  )

  message("[figuras] tema, geoms, paleta e salvamento redefinidos.")
}
