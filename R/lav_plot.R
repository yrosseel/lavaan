lavplot <- lav_plot <- function(
  model = NULL,
  infile = NULL,
  varlv = FALSE,
  placenodes = NULL,
  edgelabelsbelow = NULL,
  group.covar.indicators = FALSE,
  common.opts = list(sloped.labels = TRUE,
                     mlovcolors = c("lightgreen", "lightblue"),
                     italic = TRUE,
                     lightness = 1,
                     auto.subscript = TRUE),
  rplot = list(outfile = "",
               addgrid = TRUE),
  tikz = list(outfile = "",
              cex = 1.3,
              standalone = FALSE),
  svg = list(outfile = "",
             stroke.width = 2L,
             font.size = 20L,
             idx.font.size = 15L,
             dy = 5L,
             font.family = "Latin Modern Math, arial, Aerial, sans",
             standalone = FALSE)
) {
  tmp <- lav_get_model_info(model,
                            infile = infile,
                            varlv = varlv)
  tmp <- lav_position_nodes(tmp,
                            placenodes = placenodes,
                            edgelabelsbelow = edgelabelsbelow,
                            group.covar.indicators = group.covar.indicators)
  mc <- match.call()
  create_rplot <- !is.null(mc$rplot) || (is.null(mc$tikz) && is.null(mc$svg))
  if (create_rplot)
    do.call(lav_make_rplot, c(list(nodes.edges = tmp),
                                   common.opts,
                                   rplot))
  if (!is.null(mc$tikz))
    do.call(lav_make_tikz, c(list(nodes.edges = tmp),
                             common.opts,
                             tikz))
  if (!is.null(mc$svg))
    do.call(lav_make_svg, c(list(nodes.edges = tmp),
                             common.opts,
                             svg))
}
