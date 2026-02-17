# miscellaneous functions related to graphs
# often using an adjancency matrix as input
#
# collecting functions that were in lav_utils.R
# YR 15 Oct 2025

# find index of 'ancestors' (predictors) for all nodes in a DAG
# given an adjacency matrix B (rows are y's, columns are x's)
#
# this (speedy!) version is written by Luc De Wilde
lav_graph_get_ancestors <- function(B = NULL) {
  B <- abs(B)
  nr <- nrow(B)
  out_env <- new.env(parent = emptyenv())

  # container to hold ancestor indices per node
  out.idx <- vector("list", length = nr)

  get_ancestors <- function(nr, callers) {
    if (any(callers == nr)) {
      lav_msg_warn(gettextf("Cycle detected for element nr %d !", nr))
      return(integer(0))
    }
    x <- get0(as.character(nr), envir = out_env, ifnotfound = NULL)
    if (!is.null(x)) return(x)
    retval <- integer(0L)
    x.direct.idx <- which(B[nr, ] != 0)
    for (j in seq_along(x.direct.idx)) {
      thisone <- x.direct.idx[j]
      retval <- c(retval, thisone)
      sub <- get_ancestors(thisone, c(callers, nr))
      if (all(sub != nr)) retval <- c(retval, sub)
    }
    retval <- sort.int(unique(retval))
    assign(as.character(nr), retval, envir = out_env)
    retval
  }

  # run over each node
  for (i in seq_len(nr)) {
    out.idx[[i]] <- get_ancestors(i, integer(0))
  } # all nodes
  out.idx
}


# cluster rows/cols that are linked/connected
# return list of connected nodes
# we assume A is the square/symmetric (binary) adjacency matrix of an
# undirected graph
#
# this version written by Luc De Wilde
lav_graph_get_connected_nodes <- function(A) {

  # make sure we have square symmetric matrix
  A <- as.matrix(A)

  # A must be square
  stopifnot(nrow(A) == ncol(A))

  # A must be symmetric
  stopifnot(isSymmetric(A))

  # set diagonal to zero (just in case)
  diag(A) <- 0L

  # make it logical
  A <- (A != 0)

  # number of cols/rows
  M <- ncol(A)

  # catch diagonal A
  if (all(lavaan::lav_matrix_vech(A, diagonal = FALSE) == 0L)) {
    return(seq_len(M))
  }

  visited <- rep(FALSE, M)         # track visited nodes
  membership <- integer(M)         # component id for each node
  component.id <- 0L               # current component id

  put_node_in_component <- function(node, componentid) {
    visited[node] <<- TRUE
    membership[node] <<- componentid
    toadd <- which(A[, node])
    for (n in toadd) {
      if (!visited[n]) put_node_in_component(n, componentid)
    }
  }

  for (node in seq_len(M)) {
    if (!visited[node]) {
      component.id <- component.id + 1L
      put_node_in_component(node, component.id)
    }
  }

  membership
}

# This routine tries to order the expressions such that all expressions that
# use a certain variable come after the expression that defines this variable.
# input adj.mat is an adjacency matrix
# when there is a 1 in element (r, c) this means that the variable
# defined in expression r is used in expression c.

# contributed by ldw (adapted from function defined by Kss2k, github issue #445)
lav_graph_order_adj_mat <- function(adj.mat, warn = TRUE) {
  adjmat <- adj.mat
  n <- nrow(adjmat)
  k <- 0
  testen <- 1:n
  while (TRUE) {
    found <- FALSE
    tests <- testen # which expressions are tested in the next loop
    for (i in tests) {
      if (all(adjmat[1:n, i] == 0)) { # var defined in expression i doesn't use
        # vars not already defined, select as next
        k <- k + 1            # increment order counter
        adjmat[i, ] <- 0      # remove usage indicator of defined variable
        adjmat[i, i] <- k     # keep order of expressions in diagonal element
        found <- TRUE
        testen <- testen[testen != i] # don't test this i in next loop
      }
    }
    if (k == n) return(order(diag(adjmat)))  # all done
    if (!found) {  # no definable var found
        if (warn) {
          lav_msg_warn(gettext("unable to sort `:=` parameters;",
                               "system of defined parameters contains a cycle"))
        }
        return(seq_len(n))    # cycle detected; return original order
    }
  }
}

  # This function performs the same task as lav_graph_order_adj_mat but
  # takes as input two vectors (character or integer) of the same length
  # defining the dependencies:
  # inputs are vectors defined and definedby of the same length, where
  # for each index j for these vectors
  # node definedby[j] is used in the definition of node defined[j].
  # This routine tries to order all nodes such that all node definitions that
  # use a certain node come after the definition of that node.

  # contributed by ldw as alternative for lav_graph_order_adj_mat
  # this function is faster, requires less memory and handles also
  # character vectors
lav_graph_topological_sort <- function(defined, definedby, warn = TRUE) {
  nodes <- unique(c(definedby, defined))
  n <- length(nodes)
  rv <- vector(mode(definedby), n)
  k <- 0L
  testnodes <- rep(TRUE, n)
  testedges <- rep(TRUE, length(defined))
  while (TRUE) {
    found <- FALSE
    tests <- which(testnodes)
    for (i in tests) {
      tocheck <- nodes[i]
      if (all(tocheck != defined[testedges])) {
        k <- k + 1L
        testnodes[i] <- FALSE
        rv[k] <- tocheck
        found <- TRUE
        testedges[definedby == tocheck] <- FALSE
      }
    }
    if (k == n) return(rv)  # all done
    if (!found) {  # no definable node found
      if (warn) {
        lav_msg_warn(gettext("unable to sort;",
                             "dependencies contain a cycle"))
      }
      return(nodes)    # cycle detected; return original order
    }
  }
}

# Topological grouping and placing of nodes in a matrix.
# Nodes which need to be placed at a border are given in argument bordernodes.
# This routine does a topological sort and returns a data.frame with the nodes,
# their position in a  matrix (rows, cols) and an indication of root
# (no dependencies, indic == "r") or
# leave (no other nodes depend on this one, indic == "l"), as follows :
# the first column contains all nodes without dependencies (*1) (*2)
# the second column contains nodes with only dependency
#                                       in the first column (*1)(*3)
# the third column contains nodes with only dependencies in the first and
#     second, with at least one in the second column and so on (*3)
# (*1) nodes with successors but no successors in the next column are promoted
#     to the column just before the least column of the successors, except
#     the nodes mentioned in bordernodes and in the first column.
# (*2) the rows in the first column are chosen so that they are in the
#     neighborhood of the mean of the rows from the nodes depending on them.(*4)
# (*3) the rows in the second to last column are chosen so that the sum of the
#     row-distances to connected nodes is minimized
# (*4) if there are nodes forced to the top border in columns 2 through
#     maxcol-1, the items in column1s 1 and maxcol cannot occupy the first row;
#     analogue for the bottom
lav_graph_topological_matrix <- function(
  defined,
  definedby,
  bordernodes = character(0),
  warn = TRUE) {
  rv <- lav_graph_topological_sort(defined, definedby, warn)
  n <- length(rv)
  rvrow <- integer(n)
  rvcol <- integer(n)
  rvindic <- character(n)
  # position nodes in matrix, column by column, and mark root and leave nodes
  # rows are assigned but not yet adapted
  colmax <- 0L
  bordersincol <- 0L
  topborderfixed <- FALSE
  bottomborderfixed <- FALSE
  for (i in seq.int(n)) {
    predecessors <- definedby[rv[i] == defined]
    followers <- defined[rv[i] == definedby]
    if (length(predecessors) == 0L) {
      rvcol[i] <- 1L
      colmax[1L] <- colmax[1L] + 1L
      rvrow[i] <- colmax[1L]
      rvindic[i] <- "r" # root
    } else {
      if (length(followers) == 0L) {
        rvindic[i] <- "l"
      }
      predecessor.ind <- match(predecessors, rv)
      rvcol[i] <- max(rvcol[predecessor.ind]) + 1L
      if (length(colmax) < rvcol[i]) {
        colmax[rvcol[i]] <- 1L
        bordersincol <- 0L
      } else {
        colmax[rvcol[i]] <- colmax[rvcol[i]] + 1L
      }
      rvrow[i] <- colmax[rvcol[i]]
      if (rv[i] %in% bordernodes && rvindic[i] != "l") {
        bordersincol <- bordersincol + 1L
        topborderfixed <- TRUE
        if (bordersincol == 2L) bottomborderfixed <- TRUE
        if (bordersincol > 2L) {
          rvcol[i]  <- rvcol[i] + 1L
          colmax[rvcol[i]] <- 1L
          rvrow[i] <- colmax[rvcol[i]]
        }
      }
    }
  }
  # increment columns ? (*)
  incremented <- TRUE
  while (incremented) {
    incremented <- FALSE
    for (i in seq_along(rv)) {
      if (rvindic[i] != "l" && !(rv[i] %in% bordernodes && rvindic[i] == "r")) {
        followers <- defined[rv[i] == definedby]
        followers.ind <- match(followers, rv)
        mincol <- min(rvcol[followers.ind])
        if (rvcol[i] < mincol - 1L) {
          incremented <- TRUE
          curcol <- rvcol[i]
          rvcol[i] <- mincol - 1L
          # adapt rows in curcol
          rvrow[rvcol == curcol & rvrow > rvrow[i]] <-
            rvrow[rvcol == curcol & rvrow > rvrow[i]] - 1L
          colmax[curcol] <- colmax[curcol] - 1L
          # adapt rows in newcol
          colmax[mincol - 1L] <- colmax[mincol - 1L] + 1L
          rvrow[i] <- colmax[mincol - 1L]
        }
      }
    }
  }
  # place bordernodes in columns 2 through maxcol-1 in the first or last row
  for (col in seq.int(2L, length(colmax) - 1L)) {
    bordernodes.incol <- which(rvcol == col & rv %in% bordernodes)
    if (length(bordernodes.incol) > 0L) {
      totop <- bordernodes.incol[1L]
      if (rvrow[totop] != 1L) {
        swappie <- which(rvcol == col & rvrow == 1L)
        rvrow[swappie] <- rvrow[totop]
        rvrow[totop] <- 1L
      }
    }
    if (length(bordernodes.incol) > 1L) {
      tobottom <- bordernodes.incol[2L]
      if (rvrow[tobottom] != colmax[col]) {
        swappie <- which(rvcol == col & rvrow == colmax[col])
        rvrow[swappie] <- rvrow[tobottom]
        rvrow[tobottom] <- max(colmax)
      }
    }
  }
  # help function to order doubles as integers in a specified range
  order_doubles_interval <- function(inorder, outrange) {
    stopifnot(length(inorder) < outrange[2] - outrange[1] + 2)
    if (length(inorder) == 1) {
      if (inorder >= outrange[1L] && inorder <= outrange[2L]) {
        return(as.integer(inorder))
      }
      return(as.integer(sum(outrange) / 2))
    }
    in.order <- order(inorder)
    in.order[in.order] <- seq_along(inorder)
    if (length(inorder) == outrange[2] - outrange[1] + 1) {
      return(as.integer(outrange[1]) + in.order - 1L)
    }
    as.integer(outrange[1] + (in.order - 1L) *
      (outrange[2] - outrange[1] + 0.99) / (length(in.order) - 1L))
  }
  # arrange nodes in first column (***)
  rowmax <- max(colmax)
  addrows <- 0L
  if (topborderfixed) addrows <- 1L
  if (bottomborderfixed) addrows <- 2L
  rowmax <- max(colmax, colmax[1L] + addrows, colmax[length(colmax)] + addrows)
  nodescol1 <- which(rvcol == 1L)
  optimalrows <- sapply(nodescol1, function(ci) {
    nextnodes.ind <- which(rv[ci] == definedby)
    if (length(nextnodes.ind) == 0L) {
      colmax[1L] / 2
    } else {
      nextnodes <- defined[nextnodes.ind]
      mean(rvrow[match(nextnodes, rv)])
    }
  })
  interval1l <- c(
    if (topborderfixed) 2L else 1L,
    if (bottomborderfixed) rowmax - 1L else rowmax
  )
  rvrow[nodescol1] <- order_doubles_interval(optimalrows, interval1l)

  # adapt rows in columns to match as close as possible the rows of the
  # connected nodes in prior column(s), except for
  # borders for all but the last column: (**)
  if (length(colmax) > 1L) {
    for (c in seq.int(2L, length(colmax))) {
      if (c == length(colmax)) {
        cnodes.ind <- which(rvcol == c)
        interval <- interval1l
      } else {
        cnodes.ind <- which(rvcol == c & !(rv %in% bordernodes))
        if (length(cnodes.ind) == 0L) next
        interval <- range(rvrow[cnodes.ind])
        interval[2L] <- max(interval[2L], colmax - 1L)
      }
      optimalrows <- sapply(cnodes.ind, function(ci) {
        prevnodes.ind <- which(rv[ci] == defined)
        if (length(prevnodes.ind) == 0L) {
          rvrow[ci]
        } else {
          prevnodes <- definedby[prevnodes.ind]
          mean(rvrow[match(prevnodes, rv)])
        }
      })
      rvrow[cnodes.ind] <- order_doubles_interval(optimalrows, interval)
    }
  }
  # order on column, then row, and return
  neworder <- order(rvcol * 1000L + rvrow)
  stopifnot(length(neworder) == length(unique(neworder)))
  data.frame(
    nodes = rv[neworder],
    rows = rvrow[neworder],
    cols = rvcol[neworder],
    indic = rvindic[neworder]
  )
}
