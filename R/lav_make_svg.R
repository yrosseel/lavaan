lav_make_svg <- function(nodes.edges,
                         outfile = "",
                         sloped.labels = TRUE,
                         standalone = FALSE,
                         stroke.width = 2L,
                         font.size = 20L,
                         idx.font.size = 15L,
                         dy = 5L,
                         mlovcolors = c("lightgreen", "lightblue"),
                         lightness = 1,
                         font.family = "Latin Modern Math, arial, Aerial, sans",
                         italic = TRUE,
                         auto.subscript = TRUE
                         ) {
  textattr <- paste0('fill="black" font-size="', font.size,
              '" font-family="', font.family, '" ',
              ifelse(italic, 'font-style="italic"',''))
  tmpcol <- col2rgb(mlovcolors)
  wovcol <- paste(as.hexmode(tmpcol[, 1L]), collapse = "")
  bovcol <- paste(as.hexmode(tmpcol[, 2L]), collapse = "")
  node_elements_svg <- function(nodetiepe, noderadius, waar, stroke.width) {
    # define form, color and anchors for a node
    localradius <- noderadius
    if (nodetiepe == "varlv") localradius <- noderadius * .8
    ovxy <- localradius * sqrt(0.5)
    cvxy <- localradius * c(0.5, sqrt(0.75))
    constxy <- cvxy
    drawit <- switch(nodetiepe,
                     lv = ,
                     varlv = paste0('<circle cx="', waar[1], '" cy="', waar[2],
                                    '" r="', localradius,
                                    '" stroke-width="', stroke.width,
                                    '" stroke="black" fill="white"/>'),
                     ov = paste0('<rect width="', 2 * ovxy, '" height="',
                                 2 * ovxy, '" x="', waar[1] - ovxy, '" y="',
                                 waar[2] - ovxy,
                                 '" stroke-width="', stroke.width,
                                 '" stroke="black" fill="white" />'),
                     wov =  paste0('<rect width="', 2 * ovxy, '" height="',
                                   2 * ovxy, '" x="', waar[1] - ovxy, '" y="',
                                   waar[2] - ovxy, '" rx="', ovxy / 3, '" ry="',
                                   ovxy / 3, '" stroke-width="', stroke.width,
                                   '" stroke="black" fill="#', wovcol,
                                   '" />'),
                     bov =  paste0('<rect width="', 2 * ovxy, '" height="',
                                   2 * ovxy, '" x="', waar[1] - ovxy, '" y="',
                                   waar[2] - ovxy, '" rx="', ovxy / 3, '" ry="',
                                   ovxy / 3, '" stroke-width="', stroke.width,
                                   '" stroke="black" fill="#', bovcol,
                                   '" />'),
                     cv = paste0('<polygon points="',
                                 waar[1] - cvxy[1], ',', waar[2] - cvxy[2], ' ',
                                 waar[1] + cvxy[1], ',', waar[2] - cvxy[2], ' ',
                                 waar[1] + localradius, ',', waar[2], ' ',
                                 waar[1] + cvxy[1], ',', waar[2] + cvxy[2], ' ',
                                 waar[1] - cvxy[1], ',', waar[2] + cvxy[2], ' ',
                                 waar[1] - localradius, ',', waar[2],
                                 '" stroke-width="', stroke.width,
                                 '" stroke="black" fill="none" />'),
                     const = paste0('<polygon points="',
                                    waar[1], ',', waar[2] - localradius, ' ',
                                    waar[1] + constxy[2], ',', waar[2] + constxy[1], ' ',
                                    waar[1] - constxy[2], ',', waar[2] + constxy[1],
                                    '" stroke-width="', stroke.width,
                                    '" stroke="black" fill="none" />')
    )
    n <- c(waar[1], switch(nodetiepe,
                           lv = , varlv = , const = waar[2] - localradius,
                           ov = , wov = , bov = waar[2] - ovxy,
                           cv = waar[2] - cvxy[2]))
    s <- c(waar[1], switch(nodetiepe,
                           lv = , varlv = waar[2] + localradius,
                           ov = , wov = , bov = waar[2] + ovxy,
                           cv = waar[2]  + cvxy[2],
                           const = waar[2] + constxy[1]))
    e <- switch(nodetiepe,
                lv = , varlv = , cv = waar + c(localradius, 0),
                ov = , wov = , bov = waar + c(ovxy, 0),
                const = waar + c(constxy[2], constxy[1]))
    w <- switch(nodetiepe,
                lv = , varlv = , cv = waar + c(-localradius, 0),
                ov = , wov = , bov = waar + c(-ovxy, 0),
                const = waar + c(-constxy[2], constxy[1]))
    ne <- switch(nodetiepe,
                 lv = , varlv = , ov = , wov = ,
                 bov = waar + ovxy * c(1, -1),
                 cv = waar + c(cvxy[1], -cvxy[2]),
                 const = e)
    nw <- switch(nodetiepe,
                 lv = , varlv = , ov = , wov = ,
                 bov = waar + ovxy * c(-1, -1),
                 cv = waar + c(-cvxy[1], -cvxy[2]),
                 const = w)
    se <- switch(nodetiepe,
                 lv = , varlv = , ov = , wov = ,
                 bov = waar + ovxy * c(1, 1),
                 cv = waar + cvxy,
                 const = e)
    sw <- switch(nodetiepe,
                 lv = , varlv = , ov = , wov = ,
                 bov = waar + ovxy * c(-1, 1),
                 cv = waar + c(-cvxy[1L], cvxy[2L]),
                 const = w)
    list(drawit = drawit, n = n, ne = ne, e = e,
         se = se, s = s, sw = sw, w = w, nw = nw)
  }
  get_file_extension <- function(path) {
    if (path == "") return("")
    delen <- strsplit(path, ".", fixed = TRUE)[[1]]
    if (length(delen) > 1L) return(tolower(delen[length(delen)]))
    return("")
  }
  if (is.character(outfile) && outfile != "") {
    stopifnot(standalone || get_file_extension(outfile) == "svg",
              !standalone || get_file_extension(outfile) %in% c("htm", "html"))
  }
  mlrij <- nodes.edges$mlrij
  if (is.null(mlrij))
    lav_msg_stop(gettext(
      "nodes.edges hasn't been processed by lav_position_nodes!"))
  if (outfile == "") outfile <- stdout()
  if (is.character(outfile)) {
    zz <- file(outfile, open = "w")
    closezz <- TRUE
  } else {
    zz <- outfile
    closezz <- FALSE
  }
  nodes <- nodes.edges$nodes
  edges <- nodes.edges$edges
  nodedist <- 100
  noderadius <- 0.3
  rijen <- max(nodes$rij)
  kolommen <- max(nodes$kolom)
  nodes$rij <- nodes$rij + 1
  nodes$kolom <- nodes$kolom + 1
  if (lightness != 1) {
    mlrij <- lightness * mlrij
    nodes$kolom <- lightness * nodes$kolom
    nodes$rij <- lightness * nodes$rij
    edges$controlpt.kol <- lightness * edges$controlpt.kol
    edges$controlpt.rij <- lightness * edges$controlpt.rij
  }
  if (standalone) {
    writeLines(c(
      '<!DOCTYPE html>',
      '<html>',
      '<body>',
      '<h2>SVG diagram created by lavplot R package</h2>'),
      zz)
  }
  writeLines(c(
    paste0('<svg width="', lightness * (kolommen + 3) * nodedist, '" height="',
           lightness * (rijen + 3) * nodedist,
           '" version="1.1" xmlns="http://www.w3.org/2000/svg"',
           ' xmlns:xlink="http://www.w3.org/1999/xlink">'),
    '<rect width="100%" height="100%" fill="white" />',
    '<defs>',
    '  <marker id="arr" markerWidth="6" markerHeight="6"',
    '          refX="6" refY="2.5" orient="auto">',
    '    <path d="M 0 0 L 6 2.5 L 0 5 L 2 2.5 z" fill="black" />',
    '  </marker>',
    '  <marker id="sarr" markerWidth="6" markerHeight="6"',
    '          refX="0" refY="2.5" orient="auto">',
    '    <path d="M 0 2.5 L 6 0 L 4 2.5 L 6 5 z" fill="black" />',
    '  </marker>',
    '</defs>'),
     zz)
  plot_edge <- function(van, naar, label = "", dubbel = FALSE,
                        control = NA_real_, below = FALSE,
                        id = 0) {
    labele <- lav_format_label(label,
                               idx.font.size = idx.font.size,
                               dy = dy,
                               auto.subscript = auto.subscript)$svg
    dirvec <- naar - van
    theta <- atan2(naar[2] - van[2], naar[1] - van[1])
    if (is.na(control[1L])) { # line
      if (van[1L] <= naar[1L]) {
        writeLines(paste0('<path id="L', id, '" d="M ', van[1L],
                          ' ', van[2L], ' L ', naar[1L], " ", naar[2L],
                        '" stroke-width="', stroke.width, '" stroke="black" ',
                        ifelse(dubbel,'marker-start="url(#sarr)" ', ''),
                        'marker-end="url(#arr)" />'), zz)
      } else {
        writeLines(paste0('<path d="M ', van[1L],
                          ' ', van[2L], ' L ', naar[1L], " ", naar[2L],
                          '" stroke-width="', stroke.width, '" stroke="black" ',
                          ifelse(dubbel,'marker-start="url(#sarr)" ', ''),
                          'marker-end="url(#arr)" />'), zz)
        writeLines(paste0('<path id="L', id, '" d="M ', naar[1L],
                          ' ', naar[2L], ' L ', van[1L], " ", van[2L],
                          '" stroke-width="0" stroke="none" fill="none" />'),
                   zz)
      }
      midden <- (van + naar) * 0.5
    } else {  # path Q (quadratic Bézier)
      if (van[1L] <= naar[1L]) {
        writeLines(paste0('<path id="L', id, '" d="M ', van[1L], ' ',
                          van[2L], ' Q ', control[1L], ' ', control[2L],
                          ' ', naar[1L], " ", naar[2L],
                          '" stroke-width="', stroke.width, '" stroke="black" fill="none" ',
                          ifelse(dubbel,'marker-start="url(#sarr)" ', ''),
                          'marker-end="url(#arr)" />'), zz)
      } else {
        writeLines(paste0('<path d="M ', van[1L], ' ',
                          van[2L], ' Q ', control[1L], ' ', control[2L],
                          ' ', naar[1L], " ", naar[2L],
                          '" stroke-width="', stroke.width, '" stroke="black" fill="none" ',
                          ifelse(dubbel,'marker-start="url(#sarr)" ', ''),
                          'marker-end="url(#arr)" />'), zz)
        writeLines(paste0('<path id="L', id, '" d="M ', naar[1L], ' ',
                          naar[2L], ' Q ', control[1L], ' ', control[2L],
                          ' ', van[1L], " ", van[2L],
                          '" stroke-width="0" stroke="none" fill="none" />'),
                   zz)
      }
      midden <- 0.25 * (van + naar) + 0.5 * control
    }
    if (label != "") {
      if (sloped.labels) {
        writeLines(
          c(paste0('<text ', textattr, ' text-anchor="middle">'),
            paste0('<textPath xlink:href="#L', id, '" startOffset="50%">',
                   labele, '</textPath>'),
            '</text>'), zz)
      } else {
        if (below) {
        if (theta >= 0 && theta < pi / 2) {
          extra <- 'dy="30"'
        } else if (theta >= pi / 2) {
          extra <- 'dy="30" text-anchor="end"'
        } else if (theta < -pi/2) {
          extra <- 'dy="30"'
        } else {
          extra <- 'dy="0" text-anchor="end"'
        }
      } else {
        if (theta >= 0 && theta < pi / 2) {
          extra <- 'text-anchor="end"'
        } else if (theta >= pi / 2) {
          extra <- ' '
        } else if (theta < -pi/2) {
          extra <- 'text-anchor="end"'
        } else {
          extra <- ' '
        }
        writeLines(paste0('<text x="', midden[1L], '" y="', midden[2L],
                          '" ', textattr, ' ', extra, '>', labele, '</text>'),
                   zz)
        }
      }
    }
  }
  plot_var <- function(waar, noderadius, label = "", side = "n") {
    labele <- lav_format_label(label,
                               idx.font.size = idx.font.size,
                               dy = dy,
                               auto.subscript = auto.subscript)$svg
    thetarange <- c(pi / 6, 11 * pi / 6)
    if (side == "s") thetarange <- thetarange + 3 * pi / 2
    if (side == "e") thetarange <- thetarange + pi
    if (side == "n") thetarange <- thetarange + pi / 2
    localradius <- noderadius * 0.8
    middelpt <- switch(side,
                       n = c(0, -localradius),
                       w = c(-localradius, 0),
                       s = c(0, localradius),
                       e = c(localradius, 0))
    middelpt <- middelpt + waar
    # cirkelsegment
    straal <- localradius
    xs <- middelpt[1] + cos(thetarange) * straal
    ys <- middelpt[2] + sin(thetarange) * straal
    writeLines(paste0(
      '<path d="M ', xs[1L], ' ', ys[1L], ' A ', straal, ' ', straal ,
      ' 0 1,1 ', xs[2L], ' ', ys[2L] , '" stroke-width="', stroke.width,
      '" stroke="black" fill="none" ',
      'marker-start="url(#sarr)" marker-end="url(#arr)" />'
    ), zz)
    # label
    if (label != "") {
      writeLines(paste0('<text x="', middelpt[1L], '" y="', middelpt[2L],
                        '" text-anchor="middle" ',textattr, '>', labele,
                        '</text>'), zz)
    }
  }
  plot_node <- function(waar, tiepe, label = "") {
    labele <- lav_format_label(label,
                               idx.font.size = idx.font.size,
                               dy = dy,
                               auto.subscript = auto.subscript)$svg
    elems <- node_elements_svg(tiepe, nodedist * noderadius, waar, stroke.width)
    writeLines(c(
      elems$drawit,
      paste0('<text x="', waar[1], '" y="', waar[2], '" ',
              textattr, ' dominant-baseline="central" text-anchor="middle">',
             labele, '</text>')
    ), zz)
  }

  if (mlrij > 0L) {
    mlrij <- mlrij + lightness
    writeLines(paste0('<path d="M 1 ', mlrij * nodedist, ' L ',
                       (max(nodes$kolom) + lightness) * nodedist,
                      ' ', mlrij * nodedist, '" stroke="black"/>'),
               zz)
  }
  yrange <- nodedist * range(nodes$rij)
  xrange <- nodedist * range(nodes$kolom)
  midxy <- c(mean(xrange), mean(yrange))
  for (j in seq.int(nrow(edges))) {
    if (edges$naar[j] != edges$van[j]) {
      van <- which(nodes$id == edges$van[j])
      naar <- which(nodes$id == edges$naar[j])
      adrvan <- c(nodedist * nodes$kolom[van], nodedist * nodes$rij[van])
      elems <- node_elements_svg(nodes$tiepe[van], nodedist * noderadius,
                                 adrvan, stroke.width)
      adrvan <- elems[[edges$vananker[j]]]
      adrnaar <- c(nodedist * nodes$kolom[naar], nodedist * nodes$rij[naar])
      elems <- node_elements_svg(nodes$tiepe[naar], nodedist * noderadius,
                                 adrnaar, stroke.width)
      adrnaar <- elems[[edges$naaranker[j]]]
      if (is.na(edges$controlpt.rij[j])) {
        plot_edge(adrvan, adrnaar, edges$label[j],
                  dubbel = (edges$tiepe[j] == "~~"),
                  below = edges$labelbelow[j], id = j)
      } else {
        controlpt <- nodedist * c(edges$controlpt.kol[j] + 1,
                                  edges$controlpt.rij[j] + 1)
        plot_edge(adrvan, adrnaar, edges$label[j],
                  dubbel = (edges$tiepe[j] == "~~"),
                  below = edges$labelbelow[j],
                  control = controlpt,
                  id = j
        )
      }
    } else {
      van <- which(nodes$id == edges$van[j])
      adrvan <- c(nodedist * nodes$kolom[van], nodedist * nodes$rij[van])
      elems <- node_elements_svg(nodes$tiepe[van], nodedist * noderadius, adrvan, stroke.width)
      adrvan <- elems[[edges$vananker[j]]]
      plot_var(adrvan, noderadius * nodedist, edges$label[j], edges$vananker[j])
    }
  }
  for (j in seq.int(nrow(nodes))) {
    plot_node(nodedist * c(nodes$kolom[j], nodes$rij[j]),
              nodes$tiepe[j],
              nodes$naam[j])
  }
  writeLines("</svg>", zz)
  if (standalone) writeLines(c("</body>", "</html>"), zz)
  if (closezz) close(zz)
  return(invisible(NULL))
}
