# This function takes as input a plotinfo structure with two data.frames:
#   nodes: id, naam, tiepe, voorkeur, blok
#   edges: id, label, van, naar, tiepe
# , parameters to control the placement of nodes and edge labels, a
# switch to indicate whether indicators for which a covariance is explicit in
# the model should be grouped at the same side of the graph and a debug
# switch.
# It returns a plotinfo strucure with modified nodes and edges data.frames and
# an integer mlrij giving the position at which a line should be drawn for
# multilevel models.
# The data.frames have new columns defined as follows:
# nodes:
#   rij: index of the row where the node will be placed.
#   kolom: ndex of the column where the node will be placed.
# edges
#   vananker: character, anchor point for starting node.
#   naaranker: character, anchor point for destination node.
#   controlpt.kol: real, column position of control point if the
#     edge has to be drawn as a quadratic Beziers curve.
#   controlpt.row: real, row position of control point if the
#     edge has to be drawn as a quadratic Beziers curve.
#   labelbelow: logical, TRUE if label has to be positioned under
#     the line, initalized FALSE.
#
lav_plotinfo_positions <- function(
  plotinfo,
  placenodes = NULL,
  edgelabelsbelow = NULL,
  group.covar.indicators = FALSE,
  debug = FALSE
) {
  # add new columns to nodes and edges data.frames
  nodes <- plotinfo$nodes
  nodes$rij <- NA_integer_
  nodes$kolom <- NA_integer_
  edges <- plotinfo$edges
  edges$vananker <- NA_character_
  edges$naaranker <- NA_character_
  edges$controlpt.kol <- NA_real_
  edges$controlpt.rij <- NA_real_
  edges$labelbelow <- FALSE
  # if only one node, place it at (1, 1) and return
  if (length(nodes$rij) == 1L) {
    # Only 1 node, e.g. model = 'x ~~ x'
    nodes$rij[1L] <- 1L
    nodes$kolom[1L] <- 1L
    return(list(
      nodes = nodes,
      edges = edges,
      mlrij = 0L
    ))
  }
  # if there are multiple levels, split the plotinfo in two and call function
  # lav_plotinfo_positions_one separately for the two levels and combine the
  # results
  if (any(nodes$blok > 0L)) {
    # Multilevel, only level:1 and level:2 accepted
    nodes1 <- nodes[nodes$blok >= 2L, ]
    edges1 <- edges[edges$van %in% nodes1$id, ]
    edges1$van <- match(edges1$van, nodes1$id)
    edges1$naar <- match(edges1$naar, nodes1$id)
    nodes1$id <- seq_along(nodes1$tiepe)
    nodes1$blok <- 0L
    nodes2 <- nodes[nodes$blok == 1L, ]
    edges2 <- edges[edges$van %in% nodes2$id, ]
    edges2$van <- match(edges2$van, nodes2$id)
    edges2$naar <- match(edges2$naar, nodes2$id)
    nodes2$id <- seq_along(nodes2$tiepe)
    nodes2$blok <- 0L
    result1 <- lav_plotinfo_positions_one(
      list(nodes = nodes1, edges = edges1),
      placenodes,
      edgelabelsbelow,
      group.covar.indicators,
      debug
    )
    result2 <- lav_plotinfo_positions_one(
      list(nodes = nodes2, edges = edges2),
      placenodes,
      edgelabelsbelow,
      group.covar.indicators,
      debug
    )
    rijen1 <- max(result1$nodes$rij)
    result2$nodes$rij <- result2$nodes$rij + rijen1 + 1L
    result2$edges$controlpt.rij <- result2$edges$controlpt.rij + rijen1 + 1L
    result1$nodes$blok <- 2L
    result2$nodes$blok <- 1L
    result2$nodes$id <- result2$nodes$id + length(result1$nodes)
    result2$edges$van <- result2$edges$van + length(result1$nodes)
    result2$edges$naar <- result2$edges$naar + length(result1$nodes)
    nodes <- rbind(result1$nodes, result2$nodes)
    edges <- rbind(result1$edges, result2$edges)
    return(list(
      nodes = nodes,
      edges = edges,
      mlrij = rijen1 + 1L
    ))
  }
  # if there is only one level, call function lav_plotinfo_positions_one
  plotinfo <- list(nodes = nodes, edges = edges)
  lav_plotinfo_positions_one(
    plotinfo,
    placenodes,
    edgelabelsbelow,
    group.covar.indicators,
    debug
  )
}

# This function computes the positions for the nodes and the anchors and
# control points for the edges in the diagram for a single level.
lav_plotinfo_positions_one <- function(
  plotinfo,
  placenodes,
  edgelabelsbelow,
  group.covar.indicators,
  debug
) {
  nods <- plotinfo$nodes
  edgs <- plotinfo$edges
  # assign groupnumbers to the nodes, i.e. partition the nodes in groups
  # which belong together, e.g. indicators and the corresponding latent
  # variable
  nods <- lav_plotinfo_nodes_groups(plotinfo, group.covar.indicators)
  plotinfo <- list(nodes = nods, edges = edgs)
  # compute the information for the groups
  groups <- lav_plotinfo_groups(plotinfo)
  # order the groups in a matrix via topological sorting
  groups <- lav_groups_order(groups, plotinfo)
  # debug <- TRUE
  if (debug) {
    cat("debug start\nnodes\n")
    nods1 <- within(nods, {
      rm(blok, rij, kolom)
    })
    print(nods1)
    cat("edges\n")
    edgs1 <- within(edgs, {
      rm(vananker, naaranker, controlpt.kol, controlpt.rij)
    })
    print(edgs1)
    rm(edgs1, nods1)
    cat("matrix with groups after ordening\n")
    print(lav_groups_matrix(groups))
    cat("debug end\n")
  }
  for (g in seq_along(groups)) {
    group <- groups[[g]]
    group <- lav_group_order(group, plotinfo)
    if (debug && group$nb.nodes > 1L) {
      plot(group$offsets.lin, -group$offsets.out,
        pch = 16, xlim = c(0, 1 + max(group$offsets.lin)),
        axes = FALSE, main = paste("Group", g))
      arrows(0.5, 0, 0.5, -max(group$offsets.out), col = "blue")
      text(0.6, -0.6, "out", col = "blue")
      text(group$offsets.lin  + 0.2, -group$offsets.out,
           nods$naam[group$nodes.id])
    }
    groups[[g]] <- group
  }
  width <- sapply(groups, \(g) g$width.height[1L])
  height <- sapply(groups, \(g) g$width.height[2L])
  columns <- sapply(groups, \(g) g$matrixrowcol[2L])
  rows <- sapply(groups, \(g) g$matrixrowcol[1L])
  colnums <- sort(unique(columns))
  rownums <- sort(unique(rows))
  colwidths <- sapply(colnums, \(i) max(width[columns == i]))
  rowheights <- sapply(rownums, \(i) max(height[rows == i]))
  for (g in groups) {
    thisrow <- g$matrixrowcol[1L]
    thiscol <- g$matrixrowcol[2L]
    gtop <- 1L
    if (thisrow > 1L) gtop <- gtop + cumsum(rowheights)[thisrow - 1L]
    gleft <- 1L
    if (thiscol > 1L) gleft <- gleft + cumsum(colwidths)[thiscol - 1L]
    for (n in seq_along(g$nodes.id)) {
      k <- which(nods$id == g$nodes.id[n])
      if (g$loc == "l" || g$loc == "?") {
        nods$rij[k] <- gtop + g$offsets.lin[n]
        nods$kolom[k] <- gleft + colwidths[g$matrixrowcol[2L]] -
                         g$offsets.out[n] - 1L
      } else if (g$loc == "r") {
        nods$rij[k] <- gtop + g$offsets.lin[n]
        nods$kolom[k] <- gleft + g$offsets.out[n]
      } else if (g$loc == "t") {
        nods$rij[k] <- gtop + rowheights[1] - g$offsets.out[n] - 1L
        nods$kolom[k] <- gleft + g$offsets.lin[n]
      } else if (g$loc == "b") {
        nods$rij[k] <- gtop + g$offsets.out[n]
        nods$kolom[k] <- gleft + g$offsets.lin[n]
      }
    }
  }
  # compress graph by joining adjacent rows without common column-elements and
  # adjacent columns without common row-elements
  maxrow <- max(nods$rij)
  maxcol <- max(nods$kolom)
  for (j in seq(maxrow - 1L, 1L, -1L)) {
    col1 <- nods$kolom[nods$rij == j]
    col2 <- nods$kolom[nods$rij == j + 1L]
    if (all(is.na(match(col1, col2)))) {
      dimins <- which(nods$rij > j)
      nods$rij[dimins] <- nods$rij[dimins] - 1L
    }
  }
  for (j in seq(maxcol - 1L, 1L, -1L)) {
    col1 <- nods$rij[nods$kolom == j]
    col2 <- nods$rij[nods$kolom == j + 1L]
    if (all(is.na(match(col1, col2)))) {
      dimins <- which(nods$kolom > j)
      nods$kolom[dimins] <- nods$kolom[dimins] - 1L
    }
  }
  # check for debug
  if (debug) {
      plot(nods$kolom, -nods$rij,
        pch = 16, xlim = c(0, 1 + max(nods$kolom)), axes = FALSE,
        main = "Absolute positions nodes")
      text(nods$kolom + 0.3, -nods$rij, nods$naam)
  }
  #### place nodes demanded by user ? ####
  if (!is.null(placenodes)) {
    for (nn in names(placenodes)) {
      w <- which(nods$naam == nn)
      if (length(w) == 0) {
        lav_msg_warn(gettextf("placenodes: node name %s not found!", nn))
      }
      nods$rij[w] <- placenodes[[nn]][1L]
      nods$kolom[w] <- placenodes[[nn]][2L]
    }
  }
  #### place anchors ####
  groupmaxcol <- max(sapply(groups, \(g) g$matrixrowcol[2L]))
  groupmaxrow <- max(sapply(groups, \(g) g$matrixrowcol[1L]))
  for (j in seq_along(edgs$id)) {
    thisedge <- edgs[j, ]
    van <- nods[which(nods$id == thisedge$van), ]
    naar <- nods[which(nods$id == thisedge$naar), ]
    vangroup <- groups[[van$group]]
    if (thisedge$tiepe == "=~") { # define latent variable
      if (vangroup$matrixrowcol[[2L]] == 1L) {
        thisedge$vananker <- "w"
        thisedge$naaranker <- "e"
      } else if (vangroup$matrixrowcol[[2L]] == groupmaxcol) {
        thisedge$vananker <- "e"
        thisedge$naaranker <- "w"
      } else if (vangroup$matrixrowcol[[1L]] == 1L) {
        thisedge$vananker <- "n"
        thisedge$naaranker <- "s"
      } else if (vangroup$matrixrowcol[[1L]] == groupmaxrow) {
        thisedge$vananker <- "s"
        thisedge$naaranker <- "n"
      }
    } else if (thisedge$tiepe == "<~") { # define composite variable
      if (vangroup$matrixrowcol[[2L]] == 1L) {
        thisedge$vananker <- "e"
        thisedge$naaranker <- "w"
      } else if (vangroup$matrixrowcol[[2L]] == groupmaxcol) {
        thisedge$vananker <- "w"
        thisedge$naaranker <- "e"
      } else if (vangroup$matrixrowcol[[1L]] == 1L) {
        thisedge$vananker <- "s"
        thisedge$naaranker <- "n"
      } else if (vangroup$matrixrowcol[[1L]] == groupmaxrow) {
        thisedge$vananker <- "n"
        thisedge$naaranker <- "s"
      }
    } else if (thisedge$tiepe == "~." || thisedge$tiepe == "~") {
      # regression/varlv
      if (van$rij == naar$rij) {
        thisedge$vananker <- if (van$kolom < naar$kolom) "e" else "w"
        thisedge$naaranker <- if (van$kolom < naar$kolom) "w" else "e"
      } else if (van$kolom == naar$kolom) {
        thisedge$vananker <- if (van$rij < naar$rij) "s" else "n"
        thisedge$naaranker <- if (van$rij < naar$rij) "n" else "s"
      } else {
        thisedge <- lav_plotinfo_anchors(thisedge, van, naar,
                                         max(nods$kolom), max(nods$rij))
      }
    } else if (thisedge$tiepe == "~~~") { # (remaining) variance
        thisedge <- lav_plotinfo_anchors(thisedge, van, naar,
                                         max(nods$kolom), max(nods$rij))
    } else { # covariance
        thisedge <- lav_plotinfo_anchors(thisedge, van, naar,
                                         max(nods$kolom), max(nods$rij))
    }
    edgs[j, ] <- thisedge
  }
  #### labelsbelow demanded by user ? ####
  if (!is.null(edgelabelsbelow)) {
    for (i in seq_along(edgelabelsbelow)) {
      n1 <- which(nods$naam == edgelabelsbelow[[i]][1L])
      if (length(n1) == 0) {
        lav_msg_warn(gettextf(
          "edgelabelsbelow: node name %s not found!",
          edgelabelsbelow[[i]][1L]
        ))
      }
      n2 <- which(nods$naam == edgelabelsbelow[[i]][2L])
      if (length(n2) == 0) {
        lav_msg_warn(gettextf(
          "edgelabelsbelow: node name %s not found!",
          edgelabelsbelow[[i]][2L]
        ))
      }
      ed <- which(
        edgs$van == nods$id[n1] &
          edgs$naar == nods$id[n2]
      )
      if (length(ed) == 0L) {
        ed <- which(edgs$naar == nods$id[n1] & edgs$van == nods$id[n2])
      }
      if (length(ed) == 0L) {
        lav_msg_warn(
          gettextf(
            "edgelabelsbelow: edge %s -- %s not found!",
            nods$naam[n1],
            nods$naam[n2]
          )
        )
      }
      edgs$labelbelow[ed] <- TRUE
    }
  }
  #### RETURN ####
  list(nodes = nods, edges = edgs, mlrij = 0L)
}

lav_plotinfo_nodes_groups <- function(plotinfo, group.covar.indicators) {
  nodes <- plotinfo$nodes
  edges <- plotinfo$edges
  nodes$group <- seq.int(nrow(nodes)) # each node in its own group
  merge_groups <- function(group, i1, i2) { # merge group of i1 into group of i2
    group[group == group[i1]] <- group[i2]
    group
  }
  indicator.ids <- lav_plotinfo_indicator_ids(nodes, edges)
  # process edges to form groups
  for (j in seq.int(nrow(edges))) {
    van.id <- which(nodes$id == edges$van[j])
    naar.id <- which(nodes$id == edges$naar[j])
    merge.them <- FALSE
    if (edges$tiepe[j] == "=~") {
      if (nodes$tiepe[naar.id] != "lv" && nodes$tiepe[naar.id] != "cv") {
        merge.them <- TRUE
      }
    } else if (edges$tiepe[j] == "<~") {
      if (nodes$tiepe[van.id] != "lv" && nodes$tiepe[van.id] != "cv") {
        merge.them <- TRUE
      }
    } else if (edges$tiepe[j] == "~.") {
      merge.them <- TRUE
    } else if (edges$tiepe[j] == "~~") {
      if (group.covar.indicators &&
        any(edges$van[j] == indicator.ids) &&
        any(edges$naar[j] == indicator.ids)
      ) {
        merge.them <- TRUE
      }
    } else if (edges$tiepe[j] == "~") {
      if (nodes$tiepe[van.id] == "varlv") {
        merge.them <- TRUE
      }
    }
    if (merge.them) {
      nodes$group <- merge_groups(nodes$group, naar.id, van.id)
    }
  }
  nodes$group <- match(nodes$group, unique(nodes$group)) # renumber 1:n
  nodes
}
lav_plotinfo_indicator_ids <- function(nodes, edges) {
  # compute ids of observed indicators
  indicator.ids <- unique(c(
    edges$naar[edges$tiepe == "=~"],
    edges$van[edges$tiepe == "<~"]
  ))
  for (i in seq_along(indicator.ids)) {
    indicator.id <- which(nodes$id == indicator.ids[i])
    if (
      nodes$tiepe[indicator.id] == "lv" ||
        nodes$tiepe[indicator.id] == "cv"
    ) {
      indicator.ids[i] <- 0L
    }
  }
  indicator.ids <- indicator.ids[indicator.ids != 0L]
  indicator.ids
}
lav_plotinfo_groups <- function(plotinfo) {
  nodes <- plotinfo$nodes
  edges <- plotinfo$edges
  nb.of.groups <- max(nodes$group)
  indicator.ids <- lav_plotinfo_indicator_ids(nodes, edges)
  groups <- lapply(seq_len(nb.of.groups), function(g) {
    group <- list(
      id = g,
      nodes.id = integer(),
      edges.id = integer(),
      offsets.out = integer(),
      offsets.lin = integer(),
      nb.nodes = integer(),
      nb.indicators = 0L,
      measurement = FALSE,
      loc = "?",
      matrixrowcol = c(0L, 0L),
      width.height = c(1L, 1L),
      indic = ""
    )
    group$nodes.id <- which(nodes$group == g)
    group$nodes.id <- nodes$id[group$nodes.id]
    group$edges.id <- which(
      edges$van %in% group$nodes.id & edges$naar %in% group$nodes.id
    )
    group$edges.id <- edges$id[group$edges.id]
    group$nb.nodes <- length(group$nodes.id)
    group$nb.indicators <-
      as.integer(sum(nodes$group[nodes$id %in% indicator.ids] == g))
    group$measurement <- (group$nb.indicators > 0L)
    group
  })
  groups
}
lav_groups_matrix <- function(groups) {
  maxrow <- 1L
  maxcol <- 1L
  for (group in groups) {
    rowcol <- group$matrixrowcol
    if (rowcol[1L] > maxrow) {
      maxrow <- rowcol[1L]
    }
    if (rowcol[2L] > maxcol) {
      maxcol <- rowcol[2L]
    }
  }
  m <- matrix(0L, nrow = maxrow, ncol = maxcol)
  for (group in groups) {
    rowcol <- group$matrixrowcol
    m[rowcol[1L], rowcol[2L]] <- group$id
  }
  m
}

 # internal ordering of a group
 lav_group_order <- function(group.object, plotinfo) {
  edges <- plotinfo$edges
  group <- group.object
  group$offsets.out <- rep(0L, group$nb.nodes)
  group$offsets.lin <- rep(0L, group$nb.nodes)
  if (group$nb.nodes > 1L) {
    defined <- integer()
    definedby <- integer()
    bordernodes <- integer()
    for (j in group$edges.id) {
      thisedge <- edges[which(edges$id == j), ]
      if (thisedge$tiepe == "~.") {
        definedby <- c(definedby, thisedge$van)
        defined <- c(defined, thisedge$naar)
      } else {
        if (thisedge$tiepe == "=~") {
          definedby <- c(definedby, thisedge$naar)
          defined <- c(defined, thisedge$van)
          bordernodes <- c(bordernodes, thisedge$van)
        } else if (thisedge$tiepe == "<~") {
          definedby <- c(definedby, thisedge$van)
          defined <- c(defined, thisedge$naar)
          bordernodes <- c(bordernodes, thisedge$naar)
        } else { # include nodes not yet present indefined or definedby
          if (!(thisedge$van %in% definedby || thisedge$van %in% defined)) {
            definedby <- c(definedby, thisedge$van)
            defined <- c(defined, 999L)
          }
          if (!(thisedge$naar %in% definedby || thisedge$naar %in% defined)) {
            definedby <- c(definedby, thisedge$naar)
            defined <- c(defined, 999L)
          }
        }
      }
    }
    internaldf <- lav_graph_topological_matrix(defined, definedby,
                                               bordernodes = bordernodes)
    for (j in seq_along(group$nodes.id)) {
      group$offsets.lin[j] <- internaldf$rows[internaldf$nodes ==
                                              group$nodes.id[j]] - 1L
      group$offsets.out[j] <- max(internaldf$cols) -
          internaldf$cols[internaldf$nodes == group$nodes.id[j]]
    }
  }
  if (group$loc == "l" || group$loc == "r" || group$loc == "?") {
    group$width.height <- c(
      1L + max(group$offsets.out),
      1L + max(group$offsets.lin)
    )
  } else {
    group$width.height <- c(
      1L + max(group$offsets.lin),
      1L + max(group$offsets.out)
    )
  }
  group
 }

lav_groups_order <- function(groups, plotinfo) {
  nodes <- plotinfo$nodes
  edges <- plotinfo$edges
  dependencies <- data.frame(defined = integer(0), definedby = integer(0))
  add_dependency <- function(depends, def, defby) {
    if (any(depends$defined == def & depends$definedby == defby)) {
      return(depends)
    }
    rbind(depends, data.frame(defined = def, definedby = defby))
  }
  for (j in seq_along(edges$id)) {
    van.id <- which(nodes$id == edges$van[j])
    naar.id <- which(nodes$id == edges$naar[j])
    if (nodes$group[van.id] == nodes$group[naar.id]) {
      next
    }
    if (edges$tiepe[j] == "=~") {
      if (nodes$tiepe[naar.id] == "lv" || nodes$tiepe[naar.id] == "cv") {
        dependencies <- add_dependency(
          dependencies,
          nodes$group[van.id],
          nodes$group[naar.id]
        )
      }
    } else if (edges$tiepe[j] == "<~") {
      if (nodes$tiepe[van.id] == "lv" || nodes$tiepe[van.id] == "cv") {
        dependencies <- add_dependency(
          dependencies,
          nodes$group[naar.id],
          nodes$group[van.id]
        )
      }
    } else if (edges$tiepe[j] == "~") {
      dependencies <- add_dependency(
        dependencies,
        nodes$group[naar.id],
        nodes$group[van.id]
      )
    }
  }
  bordernodes <- integer(0)
  for (g in groups) {
    if (g$measurement) bordernodes <- c(bordernodes, g$id)
  }
  if (nrow(dependencies) == 0L) {
    for (g in groups) dependencies <- add_dependency(dependencies, 999, g$id)
  }
  groupmatrixdf <- lav_graph_topological_matrix(
    dependencies$defined,
    dependencies$definedby,
    bordernodes = bordernodes,
    warn = TRUE
  )
  for (g in seq_along(groups)) {
    group <- groups[[g]]
    gmcol <- which(groupmatrixdf$nodes == group$id)
    group$matrixrowcol <- as.integer(c(
      groupmatrixdf$rows[gmcol],
      groupmatrixdf$cols[gmcol]
    ))
    group$indic <- groupmatrixdf$indic[gmcol]
    groups[[g]] <- group
  }
  rm(groupmatrixdf)
  # Set loc = "l" for all measurement groups in first column.
  # Set loc = "r" for all measurement groups in last column.
  # Set loc = "t" for measurement groups in first row and
  #                               not first or last column.
  # Set loc = "b" for measurement groups in another row and
  #                               not first or last column.
  group.matrix <- lav_groups_matrix(groups)
  gmcols <- ncol(group.matrix)
  for (g in seq_along(groups)) {
    group <- groups[[g]]
    if (group$measurement) {
      if (group$matrixrowcol[2L] == 1L) {
        group$loc <- "l"
      } else if (group$matrixrowcol[2L] == gmcols) {
        group$loc <- "r"
      } else {
        if (group$matrixrowcol[1L] == 1L) {
          group$loc <- "t"
        } else {
          group$loc <- "b"
        }
      }
      groups[[g]] <- group
    }
  }
  groups
}

lav_points_normalform <- function(p1, p2) {
  # compute normal form of a line
  # the constant in the normal form will always be <= 0
  xy1 <- matrix(
    c(
      p1, 1, p2, 1
    ),
    byrow = TRUE,
    ncol = 3
  )
  a <- det(xy1[, 2:3])
  b <- -det(xy1[, c(
    1L, 3L
  )])
  c <- det(xy1[, 1:2])
  fac <- 1 / sqrt(a * a + b * b)
  if (c > 0) {
    fac <- -fac
  }
  fac *
    c(
      a, b, c
    )
}
lav_pointslope_normalform <- function(p, slope) {
  # compute normal form of a line defined by intercept a and slope b
  # a vertical line is defined by x = a and b = Inf or -Inf
  if (is.infinite(slope)) {
    return(c(
      1, 0, -p[1]
    ))
  }
  coeffs <- c(
    -slope, 1, slope * p[1] - p[2]
  ) /
    sqrt(1 + slope^2)
  if (coeffs[3] > 0) {
    coeffs <- -coeffs
  }
  coeffs
}
lav_lines_intersection <- function(line1, line2) {
  # get the intersection point of two lines for which coefficients of
  # equations are given
  m <- matrix(
    c(
      line1[1:2], line2[1:2]
    ),
    byrow = TRUE,
    ncol = 2L
  )
  b <- matrix(
    c(
      -line1[3], -line2[3]
    ),
    ncol = 1L
  )
  if (abs(det(m)) < 1e-12) {
    return(NA_real_)
  }
  as.vector(solve(m) %*% b)
}
lav_edge_bezierscontrolpoint <- function(van, naar, maxrij, maxcol) {
  middelpunt <- c(
    maxrij + 1, maxcol + 1
  ) /
    2
  delta <- sqrt(sum((van - naar)^2)) /
    sqrt(sum(
      c(
        maxrij - 1, maxcol - 1
      )^2
    ))
  lijn <- lav_points_normalform(van, naar)
  middenlijn <- (van + naar) / 2
  lijnmidden <- sum(
    lijn *
      c(
        middelpunt, 1
      )
  )
  if (lijnmidden > 0) {
    lijn[3] <- lijn[3] + 0.5 + delta
  } else {
    lijn[3] <- lijn[3] - 0.5 - delta
  }
  orthoslope <- lijn[2] / lijn[1]
  loodlijn <- lav_pointslope_normalform(middenlijn, orthoslope)
  lav_lines_intersection(lijn, loodlijn)
}
lav_edge_bezierscp_corner <- function(van, naar, wvannaar, maxrij, maxcol) {
  dif <- (abs(van[1L] - naar[1L]) + abs(van[2L] - naar[2L]) - 1) /
    (maxrij + maxcol - 2)
  p <- switch(wvannaar,
    ne = ,
    en = c(
      1, maxcol
    ) +
      c(
        -dif, dif
      ),
    nw = ,
    wn = c(
      1, 1
    ) +
      c(
        -dif, -dif
      ),
    se = ,
    es = c(
      maxrij, maxcol
    ) +
      c(
        dif, dif
      ),
    sw = ,
    ws = c(
      maxrij, 1
    ) +
      c(
        dif, -dif
      )
  )
  2 * (p - 0.25 * (van + naar))
}

lav_plotinfo_anchors <- function(edge, nodevan, nodenaar, maxkol, maxrij) {
  if (edge$tiepe == "~~~") {
    if (nodevan$kolom == 1L) {
      edge$vananker <- "w"
      edge$naaranker <- "w"
    } else if (nodevan$kolom == maxkol) {
      edge$vananker <- "e"
      edge$naaranker <- "e"
    } else if (nodevan$rij == 1L) {
      edge$vananker <- "n"
      edge$naaranker <- "n"
    } else if (nodevan$rij == maxrij) {
      edge$vananker <- "s"
      edge$naaranker <- "s"
    } else {
      edge$vananker <- "n"
      edge$naaranker <- "n"
    }
  } else if (edge$tiepe == "~~") {
    thrucorner <- FALSE
    if (edge$tiepe == "~~") {
      if (nodevan$kolom == nodenaar$kolom) {
        if (nodevan$kolom < maxkol / 2) {
          edge$vananker <- edge$naaranker <- "w"
        } else {
          edge$vananker <- edge$naaranker <- "e"
        }
      } else if (nodevan$rij == nodenaar$rij) {
        if (nodevan$rij < maxrij / 2) {
          edge$vananker <- edge$naaranker <- "n"
        } else {
          edge$vananker <- edge$naaranker <- "s"
        }
      } else {
        if (nodevan$rij < nodenaar$rij) {
          if (nodevan$kolom < nodenaar$kolom) {
            if (nodevan$kolom == 1L && nodenaar$rij == maxrij) {
              edge$vananker <- "w"
              edge$naaranker <- "s"
              thrucorner <- TRUE
            } else {
              edge$vananker <- "s"
              edge$naaranker <- "w"
            }
          } else {
            if (nodevan$kolom == maxkol && nodenaar$rij == maxrij) {
              edge$vananker <- "e"
              edge$naaranker <- "s"
              thrucorner = TRUE
            } else {
              edge$vananker <- "s"
              edge$naaranker <- "e"
            }
            }
        } else {
          if (nodevan$kolom < nodenaar$kolom) {
            if (nodevan$kolom == 1L && nodenaar$rij == 1L) {
              edge$vananker <- "w"
              edge$naaranker <- "n"
              thrucorner <- TRUE
            } else {
              edge$vananker <- "n"
              edge$naaranker <- "w"
            }
          } else {
            if (nodevan$kolom == maxkol && nodenaar$rij == 1L) {
              edge$vananker <- "e"
              edge$naaranker <- "n"
              thrucorner <- TRUE
            } else {
              edge$vananker <- "n"
              edge$naaranker <- "e"
            }
          }
        }
      }
      if (thrucorner) {
        bc <- lav_edge_bezierscp_corner(
                  c(
          nodevan$rij, nodevan$kolom
        ),
        c(
          nodenaar$rij, nodenaar$kolom
        ),
        paste0(edge$vananker, edge$naaranker),
        maxrij,
        maxkol
        )
      } else {
      bc <- lav_edge_bezierscontrolpoint(
        c(
          nodevan$rij, nodevan$kolom
        ),
        c(
          nodenaar$rij, nodenaar$kolom
        ),
        maxrij,
        maxkol
      )
      }
      edge$controlpt.rij <- bc[1L]
      edge$controlpt.kol <- bc[2L]
    }
  } else {
    breaks <- c(
    -pi - 0.01,
    -7 * pi / 8,
    -5 * pi / 8,
    -3 * pi / 8,
    -pi / 8,
    pi / 8,
    3 * pi / 8,
    5 * pi / 8,
    7 * pi / 8,
    pi + 0.01
    )
    winds <- c("w", "sw", "s", "se", "e", "ne", "n", "nw", "w")
    hoek <- atan2(nodevan$rij - nodenaar$rij, nodenaar$kolom - nodevan$kolom)
    wind <- cut(hoek, breaks, winds)
    edge$vananker <- as.character(wind)
    if (hoek > 0) {
      hoek <- hoek - pi
    } else {
      hoek <- hoek + pi
    }
    wind <- cut(hoek, breaks, winds)
    edge$naaranker <- as.character(wind)
  }
  edge
}
