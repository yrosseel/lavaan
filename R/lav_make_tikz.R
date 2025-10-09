beziersq2beziersc <- function(x) {
  # x contains the beziers points P1, P, P2 (P = control) for quadratic beziers
  # x is a matrix with 2 rows and 3 columns
  # returns rtval which contains P1, C1, C2, P2 so that the cubic beziers
  # with control points C1 and C2 is as 'high' as the quadratic one
  # rtval is a matrix with 2 rows and 4 columns
  matrix(c(x[ , 1L], x[ , 1L] / 3 + 2 * x[ ,2L] / 3,
         x[ , 3L] / 3 + 2 * x[ ,2L] / 3, x [ ,3L]), nrow = 2)
}
getcoord <- function(nodeid, anker, nodes, maxrij) {
  nodenr <- which(nodes$id == nodeid)
  middelpunt <- c(nodes$kolom[nodenr], maxrij - nodes$rij[nodenr])
  delta <- switch(anker, n = c(0, 0.3), ne = c(0.3, 0.3), e = c(0.3, 0),
                  se = c(0.3, -0.3), s = c(0, -0.3), sw = c(-0.3, -0.3),
                  w = c(-0.3, 0), nw = c(-0.3, 0.3))
  middelpunt + delta
}
lav_make_tikz <- function(nodes.edges,
                          outfile = "",
                          cex = 1.3,
                          sloped.labels = TRUE,
                          standalone = FALSE,
                          mlovcolors = c("lightgreen", "lightblue"),
                          lightness = 1,
                          italic = TRUE,
                          auto.subscript = TRUE
                          ) {
  tmpcol <- col2rgb(mlovcolors)
  wovcol <- paste(round(tmpcol[, 1L]/255, 2), collapse = ",")
  bovcol <- paste(round(tmpcol[, 2L]/255, 2), collapse = ",")
  nodenaam <- function(nm, blk) {
    if (blk > 0L) return(gsub("_", "", paste0("B", blk, nm)))
    return(gsub("_", "", nm))
    }
  mlrij <- nodes.edges$mlrij
  if (is.null(mlrij)) 
    lav_msg_stop(gettext(
      "nodes.edges hasn't been processed by lav_position_nodes!"))
  nodes <- nodes.edges$nodes
  edges <- nodes.edges$edges
  if (lightness != 1) {
    mlrij <- lightness * mlrij
    nodes$kolom <- lightness * nodes$kolom
    nodes$rij <- lightness * nodes$rij
    edges$controlpt.kol <- lightness * edges$controlpt.kol
    edges$controlpt.rij <- lightness * edges$controlpt.rij
  }
  if (is.character(outfile)) {
    zz <- file(outfile, open = "w")
    closezz <- TRUE
  } else {
    zz <- outfile
    closezz <- FALSE
  }
  if (standalone) writeLines(c(
    "\\documentclass{article}",
    "\\usepackage{amsmath, amssymb}",
    "\\usepackage{amsfonts}",
    "\\usepackage[utf8]{inputenc}",
    "\\usepackage[english]{babel}",
    "\\usepackage{color}",
    "\\usepackage{tikz}"), zz)
  commstyle <- paste0("draw, minimum size=", round(6 * cex), "mm")
  writeLines (c(
    "\\usetikzlibrary {shapes.geometric}",
    paste0("\\definecolor{wovcol}{rgb}{", wovcol, "}"),
    paste0("\\definecolor{bovcol}{rgb}{", bovcol, "}"),
    "\\tikzset{",
    ">=stealth,",
    paste0("x={(", cex, "cm,0cm)}, y={(0cm,", cex, "cm)},"),
    paste0("lv/.style={circle, ", commstyle, ", thick},"),
    paste0("varlv/.style={circle, draw, minimum size=", round(4 * cex), "mm, semithick},"),
    paste0("cv/.style={regular polygon, regular polygon sides=6, ", commstyle, ", thick},"),
    paste0("ov/.style={rectangle, ", commstyle,", thick},"),
    paste0("wov/.style={rectangle, rounded corners, fill=wovcol, ", commstyle, ", thick},"),
    paste0("bov/.style={rectangle, rounded corners, fill=bovcol, ", commstyle, ", thick},"),
    paste0("const/.style={regular polygon, regular polygon sides=3, ", commstyle, ", thick}"),
    "}"), zz)
  if (standalone) writeLines("\\begin{document}", zz)
  writeLines("\\begin{tikzpicture}", zz)
  maxrij <- max(nodes$rij)
  maxcol <- max(nodes$kolom)
  if (mlrij > 0L) {
    writeLines(paste("\\draw (0, ", maxrij - mlrij, ") -- (", maxcol, ",", maxrij - mlrij,
    ");", sep = ""), zz)
  }
  for (j in seq.int(nrow(nodes))) {
    xpos <- nodes$kolom[j]
    ypos <- maxrij - nodes$rij[j]
    writeLines(paste(
      "\\node[", nodes$tiepe[j], "] (", nodenaam(nodes$naam[j], nodes$blok[j]),
      ") at (", xpos, ",", ypos, ") {",
      lav_format_label(nodes$naam[j], italic = italic,
                    auto.subscript = auto.subscript)$tikz, "};", sep = ""), zz)
  }
  varlv <-any(nodes$tiepe == "varlv")
  for (j in seq.int(nrow(edges))) {
    van <- which(nodes$id == edges$van[j])
    vannaam <- nodenaam(nodes$naam[van], nodes$blok[van])
    naar <- which(nodes$id == edges$naar[j])
    naarnaam <- nodenaam(nodes$naam[naar], nodes$blok[naar])
    nodelabel <- lav_format_label(edges$label[j], italic = italic,
                                  auto.subscript = auto.subscript)$tikz
    if (van == naar) { # self
      if (nodes$kolom[van] == 1L) {
        writeLines(paste("\\path[<->] (", vannaam,
                         ") edge [in=160, out=-160, looseness=8] node[right] {",
                         nodelabel, "} (",
                         vannaam, ");",
                         sep = ""), zz)
      } else if (nodes$rij[van] == maxrij) {
        writeLines(paste("\\path[<->] (", vannaam,
                         ") edge [in=-110, out=-70, looseness=8] node[above] {",
                         nodelabel, "} (",
                         vannaam, ");",
                         sep = ""), zz)
      } else if (nodes$kolom[van] == maxcol) {
        writeLines(paste("\\path[<->] (", vannaam,
                         ") edge [in=20, out=-20, looseness=8] node[left] {",
                         nodelabel, "} (",
                         vannaam, ");",
                         sep = ""), zz)
      } else {
        writeLines(paste("\\path[<->] (", vannaam,
                         ") edge [in=110, out=70, looseness=8] node[below] {",
                         nodelabel, "} (",
                         vannaam, ");",
                         sep = ""), zz)
      }
    } else {
      anchorv <- switch(edges$vananker[j],
                        n = ".north", e = ".east", s = ".south", w = ".west")
      anchorn <- switch(edges$naaranker[j],
                        n = ".north", e = ".east", s = ".south", w = ".west")
      if (is.na(edges$controlpt.kol[j])) {
        pathtype <- " -- "
      } else {
        vanadr <- getcoord(edges$van[j], edges$vananker[j], nodes, maxrij)
        naaradr <- getcoord(edges$naar[j], edges$naaranker[j], nodes, maxrij)
        controlq <- c(edges$controlpt.kol[j], maxrij - edges$controlpt.rij[j])
        beziersc <- beziersq2beziersc(
          matrix(c(vanadr, controlq, naaradr), nrow = 2L)
        )
        pathtype <- paste0(" .. controls (", beziersc[1L, 2L] , ",",
                           beziersc[2L, 2L], ") and (", beziersc[1L, 3L] , ",",
                           beziersc[2L, 3L], ") .. ")
      }
      thelabel <- lav_format_label(edges$label[j], italic = italic,
                                   auto.subscript = auto.subscript)$tikz
      if (thelabel != "") {
        thelabel <- paste0("node[pos=0.5,",
                           ifelse(edges$labelbelow[j], "below", "above"),
                           ifelse(sloped.labels, ",sloped", ""),
                           "] {", thelabel, "} ")
      }
      pijl <- ifelse(edges$tiepe[j] %in% c("~~", "~~~"), "<->", "->")
      writeLines(paste0("\\draw[", pijl, "] (", vannaam, anchorv, ")",
                pathtype, "(", naarnaam, anchorn, ") ",
                thelabel, ";", sep = ""), zz)
    }
  }
  writeLines("\\end{tikzpicture}", zz)
  if(standalone) writeLines("\\end{document}", zz)
  if (closezz) close(zz)
  return(invisible(NULL))
}
