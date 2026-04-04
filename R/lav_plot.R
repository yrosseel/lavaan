lav_plot <- function(
  model = NULL,
  infile = NULL,
  varlv = FALSE,
  placenodes = NULL,
  edgelabelsbelow = NULL,
  group_covar_indicators = FALSE,
  common_opts = list(sloped_labels = TRUE,
                     mlovcolors = c("lightgreen", "lightblue"),
                     italic = TRUE,
                     lightness = 1,
                     auto_subscript = TRUE),
  rplot = list(outfile = "",
               addgrid = TRUE),
  tikz = list(outfile = "",
              cex = 1.3,
              standalone = FALSE),
  svg = list(outfile = "",
             stroke_width = 2L,
             font_size = 20L,
             idx_font_size = 15L,
             dy = 5L,
             font_family = "Latin Modern Math, arial, Aerial, sans",
             standalone = FALSE)
) {
  tmp <- lav_model_plotinfo(model,
                            infile = infile,
                            varlv = varlv)
  tmp <- lav_plotinfo_positions(tmp,
                            placenodes = placenodes,
                            edgelabelsbelow = edgelabelsbelow,
                            group_covar_indicators = group_covar_indicators)
  mc <- match.call()
  create_rplot <- !is.null(mc$rplot) || (is.null(mc$tikz) && is.null(mc$svg))
  if (create_rplot)
    do.call(lav_plotinfo_rgraph, c(list(plotinfo = tmp),
                                   common_opts,
                                   rplot))
  if (!is.null(mc$tikz))
    do.call(lav_plotinfo_tikzcode, c(list(plotinfo = tmp),
                             common_opts,
                             tikz))
  if (!is.null(mc$svg))
    do.call(lav_plotinfo_svgcode, c(list(plotinfo = tmp),
                             common_opts,
                             svg))
}
