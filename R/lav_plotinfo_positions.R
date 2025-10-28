lav_points_normalform <- function(P1, P2) {
  # compute normal form of a line
  # the constant in the normal form will always be <= 0
  xy1 <- matrix(c(P1, 1, P2, 1), byrow = TRUE, ncol = 3)
  a <- det(xy1[, 2:3])
  b <- -det(xy1[, c(1L, 3L)])
  c <- det(xy1[, 1:2])
  fac <- 1 / sqrt(a * a + b * b)
  if (c > 0) fac <- -fac
  fac * c(a, b, c)
}
lav_pointslope_normalform <- function(P, slope) {
  # compute normal form of a line defined by intercept a and slope b
  # a vertical line is defined by x = a and b = Inf or -Inf
  if (is.infinite(slope)) return(c(1, 0, -P[1]))
  coeffs <- c(-slope, 1, slope * P[1] - P[2]) / sqrt(1 + slope ^ 2)
  if (coeffs[3] > 0) coeffs <- -coeffs
  coeffs
}
lav_lines_intersection <- function(Line1, Line2) {
  # get the intersection point of two lines for which coefficients of
  # equations are given
  m <- matrix(c(Line1[1:2], Line2[1:2]), byrow = TRUE, ncol = 2L)
  b <- matrix(c(-Line1[3], -Line2[3]), ncol = 1L)
  if (abs(det(m)) < 1e-12) return(NA_real_)
  as.vector(solve(m) %*% b)
}
lav_edge_bezierscontrolpoint <- function(van, naar, maxrij, maxcol) {
  middelpunt <- c(maxrij + 1, maxcol + 1) / 2
  delta <- sqrt(sum((van - naar) ^ 2)) /
    sqrt(sum(c(maxrij - 1, maxcol - 1) ^ 2))
  lijn <- lav_points_normalform(van, naar)
  middenlijn <- (van + naar) / 2
  lijnmidden <- sum(lijn * c(middelpunt, 1))
  if (lijnmidden > 0) {
    lijn[3] <- lijn[3] + 0.5 + delta
  } else {
    lijn[3] <- lijn[3] - 0.5 - delta
  }
  orthoslope <- lijn[2]/lijn[1]
  loodlijn <- lav_pointslope_normalform(middenlijn, orthoslope)
  lav_lines_intersection(lijn, loodlijn)
}
lav_edge_bezierscp_corner <- function(van, naar, wvannaar, maxrij, maxcol) {
  dif <- (abs(van[1L] - naar[1L]) + abs(van[2L] - naar[2L]) - 1) /
    (maxrij + maxcol - 2)
  p <- switch(wvannaar,
              ne = , en = c(1, maxcol) + c(-dif, dif),
              nw = , wn = c(1, 1) + c(-dif, -dif),
              se = , es = c(maxrij, maxcol) + c(dif, dif),
              sw = , ws = c(maxrij, 1) + c(dif, -dif)
  )
  2 * (p - 0.25 * (van + naar))
}
lav_edgenode_distanceok <- function(edgevan, edgenaar, node) {
  # coordinates edgevan, edgenaar and node are integer,
  # so if node lies outside the rectangle defined by
  # edgevan and edgenaar the distance is at least 1, which is OK
  if ((edgevan[1L] - node[1L]) * (edgenaar[1L]-node[1L]) > 0.001) return (TRUE)
  if ((edgevan[2L] - node[2L]) * (edgenaar[2L]-node[2L]) > 0.001) return (TRUE)
  van <- edgevan - node
  naar <- edgenaar - node
  edgevec <- edgevan - edgenaar
  if (isTRUE(all.equal(edgevec, c(0, 0))) || isTRUE(all.equal(van, c(0,0))))
             return(FALSE)
  distancesquared <- sum(van * van) * (1 - sum(edgevec * (-van)) ** 2 /
                     (sum(edgevec * edgevec) * sum(van * van)))
  return(distancesquared > 0.25)
}
lav_nodeindex_position <- function(j, nodes, edges) {
  nodestocheck <- which(!is.na(nodes$rij))    # index in nodes
  positionednodes <- nodes$id[nodestocheck]   # id of the nodes to check
  check_position <- function(testpos) {
    if (any(testpos < 0)) return(FALSE)
    isok <- TRUE
    for (i in seq_len(nrow(edges))) {
      if (edges$van[i] != edges$naar[i] &&
          !is.na(nodes$rij[edges$van[i]]) &&
          !is.na(nodes$rij[edges$naar[i]])) {
        if (!lav_edgenode_distanceok(
          c(nodes$rij[edges$van[i]], nodes$kolom[edges$van[i]]),
          c(nodes$rij[edges$naar[i]], nodes$kolom[edges$naar[i]]),
          testpos
        )) {
          isok <- FALSE
          break
        }
      }
    }
    if (!isok) return(FALSE)
    for (i in seq_len(nrow(edges))) {
      if (edges$van[i] != edges$naar[i] &&
          !is.na(nodes$rij[edges$van[i]]) &&
          edges$naar[i] == nodes$naam[j]) {
        isok <- any(FALSE == vapply(seq_along(nodestocheck$rij), function(k) {
          lav_edgenode_distanceok(
            c(nodes$rij[edges$van[i]], nodes$kolom[edges$van[i]]),
            testpos,
            c(nodes$rij[k], nodes$kolom[k]))
        }, TRUE)
        )
        if (!isok) break
      }
    }
    if (!isok) return(FALSE)
    for (i in seq_len(nrow(edges))) {
      if (edges$van[i] != edges$naar[i] &&
          !is.na(nodes$rij[edges$naar[i]]) &&
          edges$van[i] == nodes$naam[j]) {
        isok <- any(FALSE == vapply(seq_along(nodestocheck$rij), function(k) {
          lav_edgenode_distanceok(
            c(nodes$rij[edges$naar[i]], nodes$kolom[edges$naar[i]]),
            testpos,
            c(nodes$rij[k], nodes$kolom[k]))
        }, TRUE)
        )
        if (!isok) break
      }
    }
    return(isok)
  }
  if (length(positionednodes) == 0) {
    availcols <- availrows <- c(1L, 3L)
  } else {
    availrows <- range(nodes$rij, na.rm = TRUE)
    availcols <- range(nodes$kolom, na.rm = TRUE)
    availrows[1L] <- availrows[1] - 1L
    availrows[2L] <- availrows[2L] + 1L
    availcols[1L] <- availcols[1] - 1L
    availcols[2L] <- availcols[2L] + 1L
  }
  connectednodes <- intersect(positionednodes,
              union(edges$naar[edges$van == j], edges$van[edges$naar == j]))
  if (length(connectednodes) == 0L) {
    bestposition <- as.integer(c(mean(availrows), mean(availcols)))
  } else if (length(connectednodes) == 1L) {
    bestposition <- as.integer(c(nodes$rij[nodes$id == connectednodes],
                               nodes$kolom[nodes$id == connectednodes]))
  } else {
    bestposition <- as.integer(round(c(
      mean(nodes$rij[nodes$id %in% connectednodes]),
      mean(nodes$kolom[nodes$id %in% connectednodes])
    )))
  }
  found <- FALSE
  for (coordist in 0:5) {
    if (coordist == 0L) {
      dir1s <- 0L
    } else {
      dir1s <- c(-coordist, coordist)
    }
    dir2s <- seq.int(-coordist, coordist)
    dir2s <- dir2s[order(abs(dir2s))]
    for (dir1 in dir1s) {
      for (dir2 in dir2s) {
        testposition <- bestposition + c(dir1, dir2)
        if (check_position(testposition)) return(testposition)
        if (dir1 != dir2) {
          testposition <- bestposition + c(dir2, dir1)
          if (check_position(testposition)) return(testposition)
        }
      }
    }
  }
  return(bestposition) # is very unlikely
}
lav_plotinfo_anchors <- function(nodes, edges, maxrij, maxkol) {
  if (all(!is.na(edges$vananker))) return(edges)
  adaptedges <- which(is.na(edges$vananker))
  breaks <- c(-pi - 0.01, -7 * pi / 8, -5 * pi / 8, -3 * pi / 8, -pi / 8,
              pi / 8, 3 * pi / 8, 5 * pi / 8, 7 * pi / 8, pi + 0.01)
  winds <- c("w", "sw", "s", "se", "e", "ne", "n", "nw", "w")
  for (i in adaptedges) {
    nodevan <- which(nodes$id == edges$van[i])
    nodenaar <- which(nodes$id == edges$naar[i])
    rrij <- range(nodes$rij)
    rkol <- range(nodes$kol)
    if (edges$tiepe[i] == "~~~") {
      if (nodes$kolom[nodevan] == rkol[1L]) {
        edges$vananker[i] <- "w"
        edges$naaranker[i] <- "w"
      } else if (nodes$kolom[nodevan] == rkol[2L]) {
        edges$vananker[i] <- "e"
        edges$naaranker[i] <- "e"
      } else if (nodes$rij[nodevan] == rrij[1L]) {
        edges$vananker[i] <- "n"
        edges$naaranker[i] <- "n"
      } else if (nodes$rij[nodevan] == rrij[2L]) {
        edges$vananker[i] <- "s"
        edges$naaranker[i] <- "s"
      } else {
        edges$vananker[i] <- "n"
        edges$naaranker[i] <- "n"
      }
      next
    }
    hoek <- atan2(nodes$rij[nodevan] - nodes$rij[nodenaar],
                  nodes$kolom[nodenaar] - nodes$kolom[nodevan])
    wind <- cut(hoek, breaks, winds)
    edges$vananker[i] <- as.character(wind)
    if (hoek > 0) hoek <- hoek - pi else hoek <- hoek + pi
    wind <- cut(hoek, breaks, winds)
    edges$naaranker[i] <- as.character(wind)
    if (edges$tiepe[i] == "~~") {
      bc <- lav_edge_bezierscontrolpoint(c(nodes$rij[nodevan], nodes$kolom[nodevan]),
                           c(nodes$rij[nodenaar], nodes$kolom[nodenaar]),
                           maxrij, maxkol)
      edges$controlpt.rij[i] <- bc[1L]
      edges$controlpt.kol[i] <- bc[2L]
    }
  }
  edges
}
lav_plotinfo_positions <- function(plotinfo,
                               placenodes = NULL,
                               edgelabelsbelow = NULL,
                               group.covar.indicators = FALSE) {
  #### lav_plotinfo_positions MAIN ####
  nodes <- plotinfo$nodes
  nodes$rij <- NA_integer_
  nodes$kolom <- NA_integer_
  edges <- plotinfo$edges
  edges$vananker <- NA_character_
  edges$naaranker <- NA_character_
  edges$controlpt.kol <- NA_real_
  edges$controlpt.rij <- NA_real_
  edges$labelbelow <- FALSE
  if (length(nodes$rij) == 1L) { # Only 1 node !
    nodes$rij[1L] <- 1L
    nodes$kolom[1L] <- 1L
    return(list(nodes = nodes, edges = edges, mlrij = 0L))
  }
  if (any(nodes$blok > 0L)) { # Multilevel, only level:1 and level:2 accepted
    nodes1 <- nodes[nodes$blok >= 2L, ]
    edges1 <- edges[edges$van %in% nodes1$id, ]
    nodes1$blok <- 0L
    nodes2 <- nodes[nodes$blok == 1L, ]
    edges2 <- edges[edges$van %in% nodes2$id, ]
    nodes2$blok <- 0L
    result1 <- lav_plotinfo_positions(list(nodes = nodes1, edges = edges1))
    result2 <- lav_plotinfo_positions(list(nodes = nodes2, edges = edges2))
    rijen1 <- max(result1$nodes$rij)
    result2$nodes$rij <- result2$nodes$rij + rijen1 + 1L
    result2$edges$controlpt.rij <- result2$edges$controlpt.rij + rijen1 + 1L
    result1$nodes$blok <- 2L
    result2$nodes$blok <- 1L
    nodes <- rbind(result1$nodes, result2$nodes)
    edges <- rbind(result1$edges, result2$edges)
    return(list(nodes = nodes, edges = edges, mlrij = rijen1 + 1L))
  }
  structnodes <- setdiff(nodes$id, union(edges$naar[edges$tiepe == "=~"],
                                         edges$van[edges$tiepe == "<~"]))
  structnodeindicators <- lapply(structnodes, function(sn) {
    union(edges$naar[edges$tiepe == "=~" & edges$van == sn],
          edges$van[edges$tiepe == "<~" & edges$naar == sn])
  })

  # create partition of structnodes via column sidegroup
  nodes$sidegroup <- rep(0L, length(nodes$naam))
  for (j in seq_along(structnodes)) {
    k <- which(nodes$id == structnodes[j])
    if (j == 1L) {
      nodes$sidegroup[k] <- j
      next
    }
    for (jj in seq.int(1L, j - 1L)) {
      if (any(structnodeindicators[[j]] %in% structnodeindicators[[jj]])) {
        nodes$sidegroup[k] <- nodes$sidegroup[structnodes[jj]]
        next
      }
      if (group.covar.indicators) {
        if (any(outer(structnodeindicators[[j]], structnodeindicators[[jj]],
            function(ind1, ind2) {
              sapply(seq_along(ind1), function(i) {
              length(edges$id[edges$van == ind1[i] & edges$naar == ind2[i]]) +
              length(edges$id[edges$naar == ind1[i] & edges$van == ind2[i]]) >
              0L
              })
        }))) {
          nodes$sidegroup[k] <- nodes$sidegroup[structnodes[jj]]
        }
      }
      if (nodes$sidegroup[k] > 0L) break
    }
    if (nodes$sidegroup[k] == 0L) nodes$sidegroup[k] <- j
  }

  # modify voorkeur to be the same within each partition group
  for (groep in unique(nodes$sidegroup)) {
    if (groep > 0L) {
      if (any(nodes$voorkeur[nodes$sidegroup == groep] == "l")) {
        nodes$voorkeur[nodes$sidegroup == groep] = "l"
      }
      if (any(nodes$voorkeur[nodes$sidegroup == groep] == "r")) {
        nodes$voorkeur[nodes$sidegroup == groep] = "r"
      }
    }
  }
  # number of groups per side (l, m or r)
  sides <- c("l", "m", "r")
  nbgroup <- rep(0L, 3)
  nbindic <- rep(0L, 3)
  for (k in seq_along(sides)) {
    nbgroup[k] <- length(unique(nodes$sidegroup[nodes$voorkeur == sides[k]]))
    welke <- nodes$id[nodes$voorkeur == sides[k]]
    nbindic[k] <-
      length(unique(edges$naar[edges$tiepe=="=~" & edges$van %in% welke])) +
      length(unique(edges$van[edges$tiepe=="<~" & edges$naar %in% welke]))
  }
  # move last 'm' groups to bottom if possible
  bottomoccupied <- FALSE
  if (nbgroup[2L] > 1L && nbindic[2L] > 8L) {
    mgroepen <- nodes$sidegroup[nodes$voorkeur == "m"]
    mgroep <- mgroepen[length(mgroepen)]
    nodes$voorkeur[nodes$sidegroup == mgroep] <- "b"
    bottomoccupied <- TRUE
  } else {
    # move last 'l' group to bottom if possible
    if (nbgroup[1L] > 1L && nbindic[1L] > 8L) {
      lgroepen <- nodes$sidegroup[nodes$voorkeur == "l"]
      lgroep <- lgroepen[length(lgroepen)]
      nodes$voorkeur[nodes$sidegroup == lgroep] <- "b"
      bottomoccupied <- TRUE
    }
    # move first 'r' group to top or last r group to bottom, if possible
    if (nbgroup[3L] > 1L && nbindic[3L] > 8L) {
      rgroepen <- nodes$sidegroup[nodes$voorkeur == "r"]
      if (nbgroup[2L] == 0L) {
        rgroep <- lgroepen[1L]
        nodes$voorkeur[nodes$sidegroup == rgroep] <- "m"
      } else {
        if (!bottomoccupied) {
          rgroep <- lgroepen[length(rgroepen)]
          nodes$voorkeur[nodes$sidegroup == rgroep] <- "b"
          bottomoccupied <- TRUE
        }
      }
    }
  }
  nodes$sidegroup <- NULL
  nodes$indicatorside <- NA_character_
  varlvs <- which(nodes$tiepe == "varlv")
  varlv <- length(varlvs) > 0L
  d_lv <- ifelse(varlv, 3L, 2L)  # distance from border
  d_ind <- ifelse(varlv, 2L, 1L) # depends on presence varlv's
  # structural part + measurement
  lvcvs <- nodes$tiepe %in% c("lv", "cv")
  if (any(lvcvs)) {
    # count number of indicators on each side
    nbindic <- c(l=0L, r=0L, m=0L, b=0L)
    for (side in names(nbindic)) {
      indicatorsassigned <- integer(0L)
      strucs <- which(nodes$voorkeur == side & lvcvs)
      for (k in strucs) {
        lvs <- edges$tiepe == "=~" & edges$van == nodes$id[k]
        cvs <- edges$tiepe == "<~" & edges$naar ==  nodes$id[k]
        indics <- unique(c(edges$naar[lvs], edges$van[cvs]))
        indicstoplace <- setdiff(indics, indicatorsassigned)
        indicatorsassigned <- union(indicatorsassigned, indics)
        nbindic[side] <- nbindic[side] + length(indicstoplace)
      }
    }
    # convert nbindic to offsets for the different sides
    if (nbindic["l"] == 0L || nbindic["r"] == 0L ||
        nbindic["l"] == nbindic["r"]) {
      nbindic["l"] <- 0L
      nbindic["r"] <- 0L
    } else {
      if (nbindic["l"] > nbindic["r"]) {
        nbindic["r"] <- as.integer((nbindic["l"] - nbindic["r"]) / 2)
        nbindic["l"] <- 0L
      } else {
        nbindic["l"] <- as.integer((nbindic["r"] - nbindic["l"]) / 2)
        nbindic["r"] <- 0L
      }
    }
    if (nbindic["m"] == 0L || nbindic["b"] == 0L ||
        nbindic["m"] == nbindic["b"]) {
      nbindic["m"] <- 0L
      nbindic["b"] <- 0L
    } else {
      if (nbindic["m"] > nbindic["b"]) {
        nbindic["b"] <- as.integer((nbindic["m"] - nbindic["b"]) / 2)
        nbindic["m"] <- 0L
      } else {
        nbindic["m"] <- as.integer((nbindic["b"] - nbindic["m"]) / 2)
        nbindic["r"] <- 0L
      }
    }
    for (side in names(nbindic)) {
      curi <- d_lv + nbindic[side]
      indicatorsassigned <- integer(0L)
      strucs <- which(nodes$voorkeur == side & lvcvs)
      for (k in strucs) {
        lvs <- edges$tiepe == "=~" & edges$van == nodes$id[k]
        cvs <- edges$tiepe == "<~" & edges$naar ==  nodes$id[k]
        edges$vananker[lvs] <- switch(side, l="w", r="e", m="n", b="s")
        edges$naaranker[lvs] <- switch(side, l="e", r="w", m="s", b="n")
        edges$vananker[cvs] <- switch(side, l="e", r="w", m="s", b="n")
        edges$naaranker[cvs] <- switch(side, l="w", r="e", m="n", b="s")
        indics <- unique(c(edges$naar[lvs], edges$van[cvs]))
        indicstoplace <- setdiff(indics, indicatorsassigned)
        indicatorsassigned <- union(indicatorsassigned, indics)
        if (side == "l" || side == "r") {
          nodes$rij[k] <- curi + as.integer(length(indicstoplace)/2.01)
          nodes$kolom[k] <- ifelse(side=="l", d_lv, 101L - d_lv)
        } else {
          nodes$kolom[k] <- curi + as.integer(length(indicstoplace)/2.01)
          nodes$rij[k] <- ifelse(side=="m", d_lv, 101L - d_lv)
        }
        for (kk in seq_along(indicstoplace)) {
          kkk <- which(nodes$id == indicstoplace[kk])
          nodes$indicatorside[kkk] <- side
          if (side == "l" || side == "r") {
            nodes$rij[kkk] <- curi
            nodes$kolom[kkk] <- ifelse(side=="l", d_ind, 101L - d_ind)
          } else {
            nodes$kolom[kkk] <- curi
            nodes$rij[kkk] <- ifelse(side=="m", d_ind, 101L - d_ind)
          }

          if (varlv) { # variances of indicators?
            varlvedges <- which(edges$naar == indicstoplace[kk] &
                                  edges$van %in% varlvs)
            if (length(varlvedges) > 0L) {
              edges$vananker[varlvedges[1L]] <- switch(side, l="e", r="w", m = "s", b = "n")
              edges$naaranker[varlvedges[1L]] <- switch(side, l="w", r="e", m = "n", b = "s")
              lvvarid <- edges$van[varlvedges[1L]]
              lvvarindex <- match(lvvarid, nodes$id)
              nodes$kolom[lvvarindex] <- switch(side, l = 1L, r = 100L, m =, b = curi)
              nodes$rij[lvvarindex] <- switch(side, m = 1L, b = 100L, l =, r = curi)
              nodes$indicatorside[lvvarindex] <- side
            }
          }
          curi <- curi + 1L
        }
      }
    }
  } else {  #### only observed variables ####
    nbnodes <- sapply(c("l", "r", "m", "b"),
                      function(n) sum(nodes$voorkeur == n))
    # offsets for the different sides
    nodesoffset <- nbnodes
    if (nodesoffset["l"] == 0L || nodesoffset["r"] == 0L ||
        nodesoffset["l"] == nodesoffset["r"]) {
      nodesoffset["l"] <- 0L
      nodesoffset["r"] <- 0L
    } else {
      if (nodesoffset["l"] > nodesoffset["r"]) {
        nodesoffset["r"] <- as.integer((nodesoffset["l"] - nodesoffset["r"]) / 2)
        nodesoffset["l"] <- 0L
      } else {
        nodesoffset["l"] <- as.integer((nodesoffset["r"] - nodesoffset["l"]) / 2)
        nodesoffset["r"] <- 0L
      }
    }
    if (nodesoffset["m"] == 0L || nodesoffset["b"] == 0L ||
        nodesoffset["m"] == nodesoffset["b"]) {
      nodesoffset["m"] <- 0L
      nodesoffset["b"] <- 0L
    } else {
      if (nodesoffset["m"] > nodesoffset["b"]) {
        nodesoffset["b"] <- as.integer((nodesoffset["m"] - nodesoffset["b"]) / 2)
        nodesoffset["m"] <- 0L
      } else {
        nodesoffset["m"] <- as.integer((nodesoffset["b"] - nodesoffset["m"]) / 2)
        nodesoffset["r"] <- 0L
      }
    }
    if (nbnodes["l"] > 0L) {
      strucs <- which(nodes$voorkeur == "l")
      nodes$rij[strucs] <- 1L + nodesoffset["l"] + seq.int(length(strucs))
      nodes$kolom[strucs] <- 1L
    }
    if (nbnodes["r"] > 0L) {
      strucs <- which(nodes$voorkeur == "r")
      nodes$rij[strucs] <- 1L + nodesoffset["r"] + seq.int(length(strucs))
      nodes$kolom[strucs] <- 100L
    }
    if (nbnodes["m"] > 0L) {
      strucs <- which(nodes$voorkeur == "m")
      nodes$kolom[strucs] <- 1L + nodesoffset["m"] + seq.int(length(strucs))
      nodes$rij[strucs] <- 1L
    }
    if (nbnodes["b"] > 0L) {
      strucs <- which(nodes$voorkeur == "b")
      nodes$kolom[strucs] <- 1L + nodesoffset["b"] + seq.int(length(strucs))
      nodes$rij[strucs] <- 100L
    }
  }
  #### remove empty rows / cols ####
  if (all(nodes$voorkeur == "")) { # no regressions defined
    # do nothing, should be covered by lav_nodeindex_position subroutine
  } else {
    # remove the holes in rows and columns, always include row/col in the middle
    rijen <- sort(unique(c(d_lv + 1L, nodes$rij)))
    nodes$rij <- match(nodes$rij, rijen)
    kolommen <- sort(unique(c(d_lv + 1L, nodes$kolom)))
    nodes$kolom <- match(nodes$kolom, kolommen)
  }
  #### adapt rij, kolom to be 1: ... ####
  minrij <- min(nodes$rij, na.rm = TRUE)
  if (minrij != 1L) nodes$rij <- nodes$rij - minrij + 1L
  maxrij <- max(nodes$rij, na.rm = TRUE)
  minkol <- min(nodes$kolom, na.rm = TRUE)
  if (minkol != 1L) nodes$kolom <- nodes$kolom - minkol + 1L
  maxkol <- max(nodes$kolom, na.rm = TRUE)
  #### compress neighboring rows (columns) if possible ####
  rijen <- max(nodes$rij, na.rm = TRUE)
  if (rijen > 2) {
    for (k in c(rijen, 2)) {
      nodesk <- which(nodes$rij == k)
      nodesk1 <- which(nodes$rij == k - 1L)
      if (length(intersect(nodes$kolom[nodesk],
                          nodes$kolom[nodesk1])) == 0L) {
        nodesadapt <- which(nodes$rij >= k)
        nodes$rij[nodesadapt] <- nodes$rij[nodesadapt] - 1L
      }
    }
  }
  kolommen <- max(nodes$kolom, na.rm = TRUE)
  if (kolommen > 2) {
    for (k in c(kolommen, 2)) {
      nodesk <- which(nodes$kolom == k)
      nodesk1 <- which(nodes$kolom == k - 1L)
      if (length(intersect(nodes$rij[nodesk],
                           nodes$rij[nodesk1])) == 0L) {
        nodesadapt <- which(nodes$kolom >= k)
        nodes$kolom[nodesadapt] <- nodes$kolom[nodesadapt] - 1L
      }
    }
  }
  #### handle nodes which are not yet placed ####
  while (any(is.na(nodes$rij))) {
    j <- which(is.na(nodes$rij))[1L]
    x <- lav_nodeindex_position(j, nodes, edges)
    nodes$rij[j] <- x[1L]
    nodes$kolom[j] <- x[2L]
  }
  maxrij <- max(nodes$rij)
  maxkol <- max(nodes$kolom)
  #### place nodes demanded by user ? ####
  if (!is.null(placenodes)) {
    for (nn in names(placenodes)) {
      w <- which(nodes$naam == nn)
      if (length(w) == 0) {
        lav_msg_warn(gettextf("placenodes: node name %s not found!", nn))
      }
      nodes$rij[w] <- placenodes[[nn]][1L]
      nodes$kolom[w] <- placenodes[[nn]][2L]
    }
  }
  indicatorids <- union(union(edges$van[edges$tiepe == "<~"],
                              edges$naar[edges$tiepe == "=~"]),
                        nodes$id[nodes$tiepe == "varlv"])
  #### adapt anchors for covariances in first or last rows/columns ####
  adaptableedges <- which(edges$tiepe == "~~" & edges$van != edges$naar)
  for (i in adaptableedges) {
    if (all(edges$van[i] != indicatorids)) next;
    nodevan <- which(nodes$id == edges$van[i])
    nodenaar <- which(nodes$id == edges$naar[i])
    if (nodes$rij[nodevan] == nodes$rij[nodenaar]) {
      if (nodes$indicatorside[nodevan] == "m") {
        edges$vananker[i] <- "n"
        edges$naaranker[i] <- "n"
        edges$controlpt.kol[i] <- (nodes$kolom[nodevan] + nodes$kolom[nodenaar]) / 2
        edges$controlpt.rij[i] <- nodes$rij[nodevan] - 0.3 -
          abs(nodes$kolom[nodevan] - nodes$kolom[nodenaar]) / 2
      } else {
        edges$vananker[i] <- "s"
        edges$naaranker[i] <- "s"
        edges$controlpt.kol[i] <- (nodes$kolom[nodevan] + nodes$kolom[nodenaar]) / 2
        edges$controlpt.rij[i] <- nodes$rij[nodevan] + 0.3 +
          abs(nodes$kolom[nodevan] - nodes$kolom[nodenaar]) / 2
      }
    } else if (nodes$kolom[nodevan] == nodes$kolom[nodenaar]) {
      if (nodes$indicatorside[nodevan] == "l") {
        edges$vananker[i] <- "w"
        edges$naaranker[i] <- "w"
        edges$controlpt.rij[i] <- (nodes$rij[nodevan] + nodes$rij[nodenaar]) / 2
        edges$controlpt.kol[i] <- nodes$kolom[nodevan] - 0.3 -
          abs(nodes$rij[nodevan] - nodes$rij[nodenaar]) / 2
      } else {
        edges$vananker[i] <- "e"
        edges$naaranker[i] <- "e"
        edges$controlpt.rij[i] <- (nodes$rij[nodevan] + nodes$rij[nodenaar]) / 2
        edges$controlpt.kol[i] <- nodes$kolom[nodevan] + 0.3 +
          abs(nodes$rij[nodevan] - nodes$rij[nodenaar]) / 2
      }
    } else {
      if (nodes$indicatorside[nodevan] == "l") {
        wvan <- "w"
      } else if (nodes$indicatorside[nodevan] == "r") {
        wvan <- "e"
      } else if (nodes$indicatorside[nodevan] == "m") {
        wvan <- "n"
      } else {
        wvan <- "s"
      }
      if (nodes$indicatorside[nodenaar] == "l") {
        wnaar <- "w"
      } else if (nodes$indicatorside[nodenaar] == "r") {
        wnaar <- "e"
      } else if (nodes$indicatorside[nodenaar] == "m") {
        wnaar <- "n"
      } else {
        wnaar <- "s"
      }
      wvannaar <- paste0(wvan, wnaar)
      edges$vananker[i] <-
        switch(wvannaar,
               nn = , ne = , nw = "n",
               ns = ifelse(nodes$kolom[nodevan] < nodes$kolom[nodenaar], "e", "w"),
               ss = , se = , sw = "s",
               sn = ifelse(nodes$kolom[nodevan] < nodes$kolom[nodenaar], "e", "w"),
               ww = , wn = , ws = "w",
               we = ifelse(nodes$rij[nodevan] < nodes$rij[nodenaar], "s", "n"),
               ee = , en = , es = "e",
               ew = ifelse(nodes$rij[nodevan] < nodes$rij[nodenaar], "n", "s")
        )
      edges$naaranker[i] <-
        switch(wvannaar,
               nn = , en = , wn = "n",
               sn = "s",
               ss = , es = , ws = "s",
               ns = "n",
               ww = , nw = , sw = "w",
               ew = "e",
               ee = , ne = , se = "e",
               we = "w"
        )
      edges$controlpt.kol[i] <-
        switch(wvannaar,
               nn = , ss = (nodes$kolom[nodevan] + nodes$kolom[nodenaar]) / 2,
               ee = nodes$kolom[nodevan] + 0.3 +
                 abs(nodes$rij[nodevan] - nodes$rij[nodenaar]) / 2,
               ww = nodes$kolom[nodevan] - 0.3 -
                 abs(nodes$rij[nodevan] - nodes$rij[nodenaar]) / 2,
               ns = , sn = nodes$kolom[nodenaar],
               we = , ew = nodes$kolom[nodevan],
               lav_edge_bezierscp_corner(
                 c(nodes$rij[nodevan], nodes$kolom[nodevan]),
                 c(nodes$rij[nodenaar], nodes$kolom[nodenaar]),
                 wvannaar, maxrij, maxkol)[2L]
        )
      edges$controlpt.rij[i] <-
        switch(wvannaar,
               nn =  nodes$rij[nodevan] - 0.3 -
                 abs(nodes$kolom[nodevan] - nodes$kolom[nodenaar]) / 2,
               ss =  nodes$rij[nodevan] + 0.3 +
                 abs(nodes$kolom[nodevan] - nodes$kolom[nodenaar]) / 2,
               ee = , ww = (nodes$rij[nodevan] + nodes$rij[nodenaar]) / 2,
               ew = , we = nodes$rij[nodenaar],
               ns = , sn = nodes$rij[nodevan],
               lav_edge_bezierscp_corner(
                 c(nodes$rij[nodevan], nodes$kolom[nodevan]),
                 c(nodes$rij[nodenaar], nodes$kolom[nodenaar]),
                 wvannaar, maxrij, maxkol)[1L]
               )
    }
  }
  #### fill anchors structural edges ####
  edges <- lav_plotinfo_anchors(nodes, edges, maxrij, maxkol)
  #### labelsbelow demanded by user ? ####
  if (!is.null(edgelabelsbelow)) {
    for (i in seq_along(edgelabelsbelow)) {
      n1 <- which(nodes$naam == edgelabelsbelow[[i]][1L])
      if (length(n1) == 0) {
        lav_msg_warn(gettextf("edgelabelsbelow: node name %s not found!",
                              edgelabelsbelow[[i]][1L]))
      }
      n2 <- which(nodes$naam == edgelabelsbelow[[i]][2L])
      if (length(n2) == 0) {
        lav_msg_warn(gettextf("edgelabelsbelow: node name %s not found!",
                              edgelabelsbelow[[i]][2L]))
      }
      ed <- which(edges$van == nodes$id[n1] & edges$naar == nodes$id[n2])
      if (length(ed) == 0L) {
        ed <- which(edges$naar == nodes$id[n1] & edges$van == nodes$id[n2])
      }
      if (length(ed) == 0L) {
        lav_msg_warn(gettextf("edgelabelsbelow: edge %s -- %s not found!",
                              nodes$naam[n1], nodes$naam[n2]))
      }
      edges$labelbelow[ed] <- TRUE
    }
  }
  #### RETURN ####
  nodes$indicatorside <- NULL
  return(list(nodes = nodes, edges = edges, mlrij = 0L))
}
