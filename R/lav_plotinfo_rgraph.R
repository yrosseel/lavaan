lav_points_beziers <- function(x, y = NULL, col = par("col"), lwd = par("lwd")) {
  if (is.null(y)) {
    if (dim(x)[1L] == 2) {
      Px <- x[1L, ]
      Py <- x[2L, ]
    } else {
      Px <- x[, 1L]
      Py <- x[, 2L]
    }
  }
  else {
    Px <- x
    Py <- y
  }
  stopifnot(length(Px) == length(Py))
  t <- seq(0, 1, length.out = 50)
  if (length(Px) == 3L) {
    PuntenX <- (1 - t) ^ 2 * Px[1] + 2 * (1-t) * t * Px[2] + t ^ 2 * Px[3]
    PuntenY <- (1 - t) ^ 2 * Py[1] + 2 * (1-t) * t * Py[2] + t ^ 2 * Py[3]
  } else {
    PuntenX <- (1 - t) ^ 3 * Px[1] + 3 * t * (1 - t) ^ 2 * Px[2] +
      3 * t ^ 2 * (1 - t) * Px[3] + t ^ 3 * Px[4]
    PuntenY <- (1 - t) ^ 3 * Py[1] + 3 * t * (1 - t) ^ 2 * Py[2] +
      3 * t ^ 2 * (1 - t) * Py[3] + t ^ 3 * Py[4]
  }
  lines(PuntenX, PuntenY, col = col, lwd = lwd)
}

lav_plotinfo_rgraph <- function(plotinfo,
                           sloped.labels = TRUE,
                           outfile = "",
                           addgrid = TRUE,
                           mlovcolors = c("lightgreen", "lightblue"),
                           lightness = 1,
                           italic = TRUE,
                           auto.subscript = TRUE
                           ) {
  font <- ifelse(italic, 3L, 1L)
  node_elements <- function(nodetiepe, noderadius) {
    # define form, color and anchors for a node type
    thetas <- switch(nodetiepe,
                     lv = ,
                     varlv = seq(0, 2 * pi, length.out = 50L),
                     ov = seq(pi / 4, 2 * pi, by = pi / 2),
                     wov = ,
                     bov =c(
                       seq(pi / 4 - pi / 10, pi / 4 + pi / 10, by = pi / 60),
                       seq(3 * pi / 4 - pi / 10, 3 * pi / 4 + pi / 10, by = pi / 60),
                       seq(5 * pi / 4 - pi / 10, 5 * pi / 4 + pi / 10, by = pi / 60),
                       seq(7 * pi / 4 - pi / 10, 7 * pi / 4 + pi / 10, by = pi / 60)
                     ),
                     cv = seq(0, 2 * pi, by = pi / 3),
                     const = seq(pi / 2, 2 * pi, by = 2 * pi / 3 )
    )
    localradius <- noderadius
    if (nodetiepe == "varlv") localradius <- noderadius * .8
    drawx <- localradius * cos(thetas)
    drawy <- localradius * sin(thetas)
    wovbovflat <- max(drawx)
    boxcol <- switch(nodetiepe,
                     lv = NA_integer_,
                     varlv = NA_integer_,
                     ov = NA_integer_,
                     wov = mlovcolors[1L],
                     bov = mlovcolors[2L],
                     cv = NA_integer_,
                     const = NA_integer_)
    n <- c(0, switch(nodetiepe,
                     lv = , varlv = , const = localradius,
                     ov = localradius * sqrt(2) / 2,
                     wov = , bov = wovbovflat,
                     cv = localradius * sqrt(3) / 2))
    s <- c(0, switch(nodetiepe,
                     lv = , varlv = -localradius,
                     ov = -localradius * sqrt(2) / 2,
                     wov = , bov = -wovbovflat,
                     cv = -localradius * sqrt(3) / 2,
                     const = -localradius * 0.5))
    e <- switch(nodetiepe,
                lv = , varlv = , cv = c(localradius, 0),
                ov = c(localradius * sqrt(2) / 2, 0),
                wov = , bov = c(wovbovflat, 0),
                const = c(localradius * sqrt(3) / 2, -localradius * 0.5))
    w <- c(-e[1L], e[2L])
    ne <- switch(nodetiepe,
                 lv = , varlv = , ov = , wov = ,
                 bov = localradius * sqrt(0.5) * c(1, 1),
                 cv = localradius * c(0.5, sqrt(3) / 2),
                 const = e)
    nw <- c(-ne[1L], ne[2L])
    se <- switch(nodetiepe,
                 lv = , varlv = , ov = , wov = ,
                 bov = localradius * sqrt(0.5) * c(1, -1),
                 cv = localradius * c(-0.5, sqrt(3) / 2),
                 const = e)
    sw <- c(-se[1L], se[2L])
    list(drawx = drawx, drawy = drawy, boxcol = boxcol, n = n, ne = ne, e = e,
         se = se, s = s, sw = sw, w = w, nw = nw)
  }

  vecrotate <- function(vec, angle) {
    c(cos(angle)*vec[1]+sin(angle)*vec[2],
      -sin(angle)*vec[1]+cos(angle)*vec[2])
  }
  plot_arrow <- function(tip, dirvec) {
    unitvec <- dirvec / sqrt(sum(dirvec ^ 2))
    arrowangle <- pi * 25 / 180
    arrowinset <- 0.4
    args <- rbind(tip,
                  tip + vecrotate(-unitvec * arrowlength, arrowangle),
                  tip - unitvec * arrowlength * (1 - arrowinset),
                  tip + vecrotate(-unitvec * arrowlength, -arrowangle))
    polygon(args, col = "black", border = NA)
  }
  plot_edge <- function(van, naar, label = "", dubbel = FALSE,
                        control = NA_real_, below = FALSE, txtcex = 0.9) {
    labele <- lav_label_code(label, auto.subscript = auto.subscript)$r
    dirvec <- naar - van
    theta <- atan2(naar[2] - van[2], naar[1] - van[1])
    srt <- ifelse(sloped.labels, 180 * theta / pi, 0)
    if (srt > 90) srt <- srt - 180
    if (srt < -90) srt <- srt + 180
    if (is.na(control[1L])) {
      args <- rbind(van, naar)
      lines(args, lwd = 2)
      plot_arrow(naar, dirvec)
      if (dubbel) plot_arrow(van, -dirvec)
      midden <- (van + naar) * 0.5
    } else {
      # gebogen lijn (quadratic lav_points_beziers)
      lav_points_beziers(rbind(van, control, naar), lwd = 2)
      midden <- (van + naar) / 4 + control / 2
      plot_arrow(naar, naar - control)
      if (dubbel) plot_arrow(van, van - control)
    }
    if (label != "") {
      if (below) {
        if (theta >= 0 && theta < 90) {
          text(midden[1L], midden[2L], labele, adj = c(0, 1),
               srt = srt, cex = txtcex, font = font)
        } else if (theta >= 90) {
          text(midden[1L], midden[2L], labele, adj = c(1, 1),
               srt = srt, cex = txtcex, font = font)
        } else if (theta < -90) {
          text(midden[1L], midden[2L], labele, adj = c(0, 1),
               srt = srt, cex = txtcex, font = font)
        } else {
          text(midden[1L], midden[2L], labele, adj = c(1, 1),
               srt = srt, cex = txtcex, font = font)
        }
      } else {
        if (theta >= 0 && theta < 90) {
          text(midden[1L], midden[2L], labele, adj = c(1, 0),
               srt = srt, cex = txtcex, font = font)
        } else if (theta >= 90) {
          text(midden[1L], midden[2L], labele, adj = c(0, 0),
               srt = srt, cex = txtcex, font = font)
        } else if (theta < -90) {
          text(midden[1L], midden[2L], labele, adj = c(1, 0),
               srt = srt, cex = txtcex, font = font)
        } else {
          text(midden[1L], midden[2L], labele, adj = c(0, 0),
               srt = srt, cex = txtcex, font = font)
        }
      }
    }
  }
  plot_var <- function(waar, noderadius, label = "", side = "n", txtcex = 0.9) {
    labele <- lav_label_code(label, auto.subscript = auto.subscript)$r
    thetarange <- c(pi / 6, 11 * pi / 6)
    if (side == "s") thetarange <- thetarange + pi / 2
    if (side == "e") thetarange <- thetarange + pi
    if (side == "n") thetarange <- thetarange + 3 * pi / 2
    localradius <- noderadius * 0.8
    middelpt <- switch(side,
                       n = c(0, localradius),
                       w = c(-localradius, 0),
                       s = c(0, -localradius),
                       e = c(localradius, 0))
    middelpt <- middelpt + waar
    # cirkelsegment
    thetas <- seq(thetarange[1L], thetarange[2L], length.out = 40)
    straal <- localradius
    xs <- middelpt[1] + cos(thetas) * straal
    ys <- middelpt[2] + sin(thetas) * straal
    lines(xs, ys)
    # pijlen
    plot_arrow(c(xs[1], ys[1]), c(sin(thetarange[1]), -cos(thetarange[1])))
    plot_arrow(c(xs[40], ys[40]), c(-sin(thetarange[2]), cos(thetarange[2])))
    # label
    if (label != "")
      text(middelpt[1L], middelpt[2L], labele, adj = 0.5, cex = txtcex * 0.8,
           font = font)
  }
  plot_node <- function(waar, tiepe, label = "", txtcex = 0.9) {
    labele <- lav_label_code(label, auto.subscript = auto.subscript)$r
    elems <- node_elements(tiepe, noderadius)
    x <- waar[1] + elems$drawx
    y <- waar[2] + elems$drawy
    polygon(x, y, col = elems$boxcol, lwd = 1)
    text(waar[1L], waar[2L], labele, adj = 0.5, cex = txtcex, font = font)
  }
  mlrij <- plotinfo$mlrij
  if (is.null(mlrij))
    lav_msg_stop(gettext(
      "plotinfo hasn't been processed by lav_plotinfo_positions!"))
  nodes <- plotinfo$nodes
  edges <- plotinfo$edges
  noderadius <- 0.3
  arrowlength <- noderadius / 3
  rijen <- max(nodes$rij)
  kolommen <- max(nodes$kolom)
  if (outfile != "") png(outfile, 960, 960, "px")
  opar <- par(mar = c(1L,1L, 1L,1L) + 0.1)
  plot.default(x = c(0, lightness * kolommen + 1),
               c(0, lightness * rijen + 1), type="n",
               xlab = "", ylab = "", axes = FALSE, asp = 1)
  if (addgrid) {
    abline(v = seq.int(kolommen) * lightness,
           h = seq.int(rijen) * lightness,
           lwd = 1,
           lty = "dotted", col = "lightgray")
    text(seq.int(kolommen) * lightness, 0.3,
         labels=seq.int(kolommen), adj = 1, cex = 0.7)
    text(0.3, seq.int(rijen) * lightness,
         labels=seq.int(rijen, 1), adj = 1, cex = 0.7)
  }
  if (mlrij > 0L) abline(h = rijen - mlrij + 1, lwd = 2)
  yrange <- rijen - range(nodes$rij) + 1
  xrange <- range(nodes$kolom)
  midxy <- c(mean(xrange), mean(yrange))
  for (j in seq.int(nrow(edges))) {
    if (edges$naar[j] != edges$van[j]) {
      van <- which(nodes$id == edges$van[j])
      naar <- which(nodes$id == edges$naar[j])
      adrvan <- lightness *
        c(nodes$kolom[van], rijen - nodes$rij[van] + 1)
      elems <- node_elements(nodes$tiepe[van], noderadius)
      adrvan <- adrvan + elems[[edges$vananker[j]]]
      adrnaar <- lightness *
        c(nodes$kolom[naar], rijen - nodes$rij[naar] + 1)
      elems <- node_elements(nodes$tiepe[naar], noderadius)
      adrnaar <- adrnaar + elems[[edges$naaranker[j]]]
      if (is.na(edges$controlpt.rij[j])) {
        plot_edge(adrvan, adrnaar, edges$label[j],
                  dubbel = (edges$tiepe[j] == "~~"),
                  below = edges$labelbelow[j])
      } else {
        controlpt <- lightness * c(edges$controlpt.kol[j],
                       rijen - edges$controlpt.rij[j] + 1)
        plot_edge(adrvan, adrnaar, edges$label[j],
                  dubbel = (edges$tiepe[j] == "~~"),
                  below = edges$labelbelow[j],
                  control = controlpt
                  )
      }
    } else {
      van <- which(nodes$id == edges$van[j])
      adrvan <- lightness * c(nodes$kolom[van], rijen - nodes$rij[van] + 1)
      elems <- node_elements(nodes$tiepe[van], noderadius)
      adrvan <- adrvan + elems[[edges$vananker[j]]]
      plot_var(adrvan, noderadius, edges$label[j], edges$vananker[j])
    }
  }
  for (j in seq.int(nrow(nodes))) {
    plot_node(lightness *
                c(nodes$kolom[j], rijen - nodes$rij[j] + 1),
              nodes$tiepe[j],
              nodes$naam[j])
  }
  par(opar)
  if (outfile != "") dev.off()
  return(invisible(NULL))
}
