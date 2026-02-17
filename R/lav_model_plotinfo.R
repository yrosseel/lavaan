# extract info from model
lav_model_plotinfo <- function(model = NULL, infile = NULL, varlv = FALSE) {
  if (is.null(model) == is.null(infile)) {
    lav_msg_stop(gettext("either model or infile must be specified"))
  }
  edge_label <- function(label, fix) {
    if (label == "" && fix == "") {
      return("")
    }
    if (label == "") {
      return(fix)
    }
    if (fix == "") {
      return(label)
    }
    paste(label, fix, sep = "=")
  }
  if (!is.null(infile)) {
    if (!file.exists(infile)) {
      lav_msg_stop(gettextf("Infile %s not found!", infile))
    }
    model <- readLines(infile)
  }
  if (inherits(model, "lavaan")) {
    model <- as.data.frame(model@ParTable)
    model <- model[model$user == 1, ]
    model$fixed <- round(model$est, 3L)
    if (max(model$block) == 1L) model$block <- rep(0L, length(model$block))
    model$est <- NULL
  }
  if (is.list(model) && !is.null(model$op) && !is.null(model$lhs) &&
    !is.null(model$rhs) && !is.null(model$label) &&
    !is.null(model$fixed)) {
    tbl <- as.data.frame(model)
  } else if (is.character(model)) {
    tbl <- lavParseModelString(
      model.syntax = model,
      as.data.frame. = TRUE
    )
  } else {
    lav_msg_stop(gettext(
      "model, or content of infile, not in interpretable form!"
    ))
  }
  if (is.null(tbl$block)) {
    tbl$block <- 0L
  } else if (all(tbl$block == tbl$block[1L])) {
    tbl$block <- 0L
  } else {
    blockmin <- min(tbl$block)
    tbl$block[tbl$block == blockmin] <- 999L
    tbl$block[tbl$block != 999L] <- 2L
    tbl$block[tbl$block == 999L] <- 1L
  }
  create_node_id <- function(blok, naam) {
    if (blok == 0L) {
      naam
    } else {
      paste(blok, naam, sep = ".")
    }
  }
  maxedges <- nrow(tbl)
  maxnodes <- 2 * maxedges
  nodes <- data.frame(
    id = character(maxnodes),
    naam = character(maxnodes),
    tiepe = character(maxnodes), # ov, lv, varlv, cv, wov, bov, const
    # cv: composite; wov = within; bov = between; const = intercept
    blok = integer(maxnodes)
  )
  edges <- data.frame(
    id = integer(maxedges),
    label = character(maxedges),
    van = character(maxedges),
    naar = character(maxedges),
    tiepe = character(maxedges) # lavaan operator, behalve ~~ self -> "~~~"
  )
  curnode <- 0L
  curedge <- 0L
  for (i in seq.int(nrow(tbl))) {
    if (tbl$op[i] == "=~") {
      #### =~ : is manifested by ####
      # lhs node
      jl <- match(create_node_id(tbl$block[i], tbl$lhs[i]),
                  nodes$id, nomatch = 0L)
      if (jl == 0L) {
        curnode <- curnode + 1L
        jl <- curnode
        nodes$naam[curnode] <- tbl$lhs[i]
        nodes$tiepe[curnode] <- "lv"
        nodes$blok[curnode] <- tbl$block[i]
        nodes$id[curnode] <- create_node_id(nodes$blok[curnode],
                                            nodes$naam[curnode])
      } else {
        nodes$tiepe[jl] <- "lv"
      }
      # rhs node
      jr <- match(create_node_id(tbl$block[i], tbl$rhs[i]),
                  nodes$id, nomatch = 0L)
      nodetype <- "ov"
      if (length(unique(tbl$block[tbl$rhs == tbl$rhs[i]])) > 1L) {
        nodetype <- switch(tbl$block[i],
          "wov",
          "bov"
        )
      }
      if (jr == 0L) {
        curnode <- curnode + 1L
        jr <- curnode
        nodes$naam[curnode] <- tbl$rhs[i]
        nodes$tiepe[curnode] <- nodetype
        nodes$blok[curnode] <- tbl$block[i]
        nodes$id[curnode] <- create_node_id(nodes$blok[curnode],
                                            nodes$naam[curnode])
      } else {
        if (nodes$tiepe[jr] == "") nodes$tiepe[jr] <- nodetype
      }
      # edge
      curedge <- curedge + 1L
      edges$id[curedge] <- curedge
      edges$label[curedge] <- edge_label(tbl$label[i], tbl$fixed[i])
      edges$van[curedge] <- nodes$id[jl]
      edges$naar[curedge] <- nodes$id[jr]
      edges$tiepe[curedge] <- tbl$op[i]
    } else if (tbl$op[i] == "<~") {
      #### <~ : is a result of ####
      # lhs node
      jl <- match(create_node_id(tbl$block[i], tbl$lhs[i]),
                  nodes$id, nomatch = 0L)
      if (jl == 0L) {
        curnode <- curnode + 1L
        jl <- curnode
        nodes$naam[curnode] <- tbl$lhs[i]
        nodes$tiepe[curnode] <- "cv"
        nodes$blok[curnode] <- tbl$block[i]
        nodes$id[curnode] <- create_node_id(nodes$blok[curnode],
                                            nodes$naam[curnode])
      } else {
        nodes$tiepe[jl] <- "cv"
      }
      # rhs node
      jr <- match(create_node_id(tbl$block[i], tbl$rhs[i]),
                  nodes$id, nomatch = 0L)
      nodetype <- "ov"
      if (length(unique(tbl$block[tbl$rhs == tbl$rhs[i]])) > 1L) {
        nodetype <- switch(tbl$block[i],
          "wov",
          "bov"
        )
      }
      if (jr == 0L) {
        curnode <- curnode + 1L
        jr <- curnode
        nodes$naam[curnode] <- tbl$rhs[i]
        nodes$tiepe[curnode] <- nodetype
        nodes$blok[curnode] <- tbl$block[i]
        nodes$id[curnode] <- create_node_id(nodes$blok[curnode],
                                            nodes$naam[curnode])
      } else {
        if (nodes$tiepe[jr] == "") nodes$tiepe[jr] <- nodetype
      }
      # edge
      curedge <- curedge + 1L
      edges$id[curedge] <- curedge
      edges$label[curedge] <- edge_label(tbl$label[i], tbl$fixed[i])
      edges$van[curedge] <- nodes$id[jr]
      edges$naar[curedge] <- nodes$id[jl]
      edges$tiepe[curedge] <- tbl$op[i]
    } else if (tbl$op[i] == "~") {
      #### ~ : is regressed on ####
      # lhs node
      jl <- match(create_node_id(tbl$block[i], tbl$lhs[i]),
                                 nodes$id, nomatch = 0L)
      nodetype <- "ov"
      if (length(unique(tbl$block[tbl$rhs == tbl$rhs[i]])) > 1L) {
        nodetype <- switch(tbl$block[i],
          "wov",
          "bov"
        )
      }
      if (jl == 0L) {
        curnode <- curnode + 1L
        jl <- curnode
        nodes$naam[curnode] <- tbl$lhs[i]
        nodes$tiepe[curnode] <- nodetype
        nodes$blok[curnode] <- tbl$block[i]
        nodes$id[curnode] <- create_node_id(nodes$blok[curnode],
                                            nodes$naam[curnode])
      } else {
        if (nodes$tiepe[jl] == "") nodes$tiepe[jl] <- nodetype
      }
      # rhs node
      jr <- match(create_node_id(tbl$block[i], tbl$rhs[i]),
                  nodes$id, nomatch = 0L)
      if (jr == 0L) {
        curnode <- curnode + 1L
        jr <- curnode
        nodes$naam[curnode] <- tbl$rhs[i]
        nodes$tiepe[curnode] <- nodetype
        nodes$blok[curnode] <- tbl$block[i]
        nodes$id[curnode] <- create_node_id(nodes$blok[curnode],
                                            nodes$naam[curnode])
      } else {
        if (nodes$tiepe[jr] == "") nodes$tiepe[jr] <- nodetype
      }
      # edge
      curedge <- curedge + 1L
      edges$id[curedge] <- curedge
      edges$label[curedge] <- edge_label(tbl$label[i], tbl$fixed[i])
      edges$van[curedge] <- nodes$id[jr]
      edges$naar[curedge] <- nodes$id[jl]
      edges$tiepe[curedge] <- tbl$op[i]
    } else if (tbl$op[i] == "~1") {
      #### ~1 : intercept ####
      # lhs node
      jl <- match(create_node_id(tbl$block[i], tbl$lhs[i]),
                  nodes$id, nomatch = 0L)
      nodetype <- "ov"
      if (length(unique(tbl$block[tbl$rhs == tbl$lhs[i] |
        tbl$lhs == tbl$lhs[i]])) > 1L) {
        nodetype <- switch(tbl$block[i],
          "wov",
          "bov"
        )
      }
      if (jl == 0L) {
        curnode <- curnode + 1L
        jl <- curnode
        nodes$naam[curnode] <- tbl$lhs[i]
        nodes$tiepe[curnode] <- nodetype
        nodes$blok[curnode] <- tbl$block[i]
        nodes$id[curnode] <- create_node_id(nodes$blok[curnode],
                                            nodes$naam[curnode])
      } else {
        if (nodes$tiepe[jl] == "") nodes$tiepe[jl] <- nodetype
      }
      # rhs node
      jr <- 0L
      curnode <- curnode + 1L
      jr <- curnode
      nodes$naam[curnode] <- paste0("1van", tbl$lhs[i])
      nodes$tiepe[curnode] <- "const"
      nodes$blok[curnode] <- tbl$block[i]
      nodes$id[curnode] <- create_node_id(nodes$blok[curnode],
                                          nodes$naam[curnode])
      # edge
      curedge <- curedge + 1L
      edges$id[curedge] <- curedge
      edges$label[curedge] <- edge_label(tbl$label[i], tbl$fixed[i])
      edges$van[curedge] <- nodes$id[jr]
      edges$naar[curedge] <- nodes$id[jl]
      edges$tiepe[curedge] <- "~"
    } else if (tbl$op[i] == "~~") {
      #### ~~ : is correlated with ####
      # lhs node
      jl <- match(create_node_id(tbl$block[i], tbl$lhs[i]),
                  nodes$id, nomatch = 0L)
      nodetype <- "ov"
      if (length(unique(tbl$block[tbl$rhs == tbl$lhs[i] |
        tbl$lhs == tbl$lhs[i]])) > 1L) {
        nodetype <- switch(tbl$block[i],
          "wov",
          "bov"
        )
      }
      if (jl == 0L) {
        curnode <- curnode + 1L
        jl <- curnode
        nodes$naam[curnode] <- tbl$lhs[i]
        nodes$tiepe[curnode] <- nodetype
        nodes$blok[curnode] <- tbl$block[i]
        nodes$id[curnode] <- create_node_id(nodes$blok[curnode],
                                            nodes$naam[curnode])
      } else {
        if (nodes$tiepe[jl] == "") nodes$tiepe[jl] <- nodetype
      }
      # rhs node
      jr <- match(create_node_id(tbl$block[i], tbl$rhs[i]),
                  nodes$id, nomatch = 0L)
      nodetype <- "ov"
      if (length(unique(tbl$block[tbl$rhs == tbl$lhs[i] |
        tbl$lhs == tbl$lhs[i]])) > 1L) {
        nodetype <- switch(tbl$block[i],
          "wov",
          "bov"
        )
      }
      if (jr == 0L) {
        curnode <- curnode + 1L
        jr <- curnode
        nodes$naam[curnode] <- tbl$rhs[i]
        nodes$tiepe[curnode] <- nodetype
        nodes$blok[curnode] <- tbl$block[i]
        nodes$id[curnode] <- create_node_id(nodes$blok[curnode],
                                            nodes$naam[curnode])
      } else {
        if (nodes$tiepe[jr] == "") nodes$tiepe[jr] <- nodetype
      }
      # edge
      curedge <- curedge + 1L
      edges$id[curedge] <- curedge
      edges$label[curedge] <- edge_label(tbl$label[i], tbl$fixed[i])
      if (varlv && jl == jr) { # prepare for handling varlv
        edges$label[curedge] <- paste(tbl$label[i], tbl$fixed[i], sep = "=")
      }
      edges$van[curedge] <- nodes$id[jr]
      edges$naar[curedge] <- nodes$id[jl]
      edges$tiepe[curedge] <- ifelse(jl == jr, "~~~", tbl$op[i])
    }
  }
  # aanpassingen voor varlv
  if (varlv && any(edges$tiepe == "~~~")) {
    #       wijzigen varianties
    welke <- which(edges$tiepe == "~~~")
    varlvnodes <- rep(0L, length(welke))
    lvnodes <- rep(0L, length(welke))
    for (ji in seq_along(welke)) {
      j <- welke[ji]
      curnode <- curnode + 1L
      nodes$naam[curnode] <- gsub("=.*$", "", edges$label[j])
      if (nodes$naam[curnode] == "")
        nodes$naam[curnode] <- paste0("var", edges$van[j])
      nodes$tiepe[curnode] <- "varlv"
      nodes$blok[curnode] <- nodes$blok[which(nodes$id == edges$van[j])[[1L]]]
      nodes$id[curnode] <- create_node_id(nodes$blok[curnode],
                                          nodes$naam[curnode])
      varlvnodes[ji] <- curnode
      lvnodes[ji] <- edges$van[j]
      edges$van[j] <- nodes$id[curnode]
      edges$tiepe[j] <- "~."
      edges$label[j] <- ifelse(grepl("=", edges$label[j]),
        gsub("^.*=", "", edges$label[j]), ""
      )
    }
    #       wijzigen covarianties
    for (j1 in seq_along(welke)) {
      for (j2 in seq_along(welke)) {
        if (j1 != j2 && any(edges$tiepe == "~~" &
          edges$van == lvnodes[j1] & edges$naar == lvnodes[j2])) {
          edg <- which(edges$tiepe == "~~" &
            edges$van == lvnodes[j1] & edges$naar == lvnodes[j2])[[1L]]
          edges$van[edg] <- varlvnodes[j1]
          edges$naar[edg] <- varlvnodes[j2]
        }
      }
    }
  }
  nodes <- nodes[seq.int(curnode), ]
  edges <- edges[seq.int(curedge), ]
  edges$label <- trimws(edges$label)
  if (any(grepl("[<>&]", c(edges$label, nodes$naam)))) {
    lav_msg_warn(gettext(
      "some labels contain '<', '>' or '&', which can result in errors!"
    ))
  }
  list(nodes = nodes, edges = edges)
}
