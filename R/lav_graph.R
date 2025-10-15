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
  OUTENV <- new.env(parent = emptyenv())

  # container to hold ancestor indices per node
  out.idx <- vector("list", length = nr)

  get_ancestors <- function(nr, callers) {
    if (any(callers == nr)) {
      lav_msg_warn(gettextf("Cycle detected for element nr %d !", nr))
      return(integer(0));
    }
    x = get0(as.character(nr), envir = OUTENV, ifnotfound = NULL)
    if (!is.null(x)) return(x);
    retval <- integer(0L);
    x.direct.idx <- which(B[nr,] != 0)
    for (j in seq_along(x.direct.idx)) {
      thisone <- x.direct.idx[j]
      retval <- c(retval, thisone)
      sub <- get_ancestors(thisone, c(callers, nr))
      if (all(sub != nr)) retval = c(retval, sub)
    }
    retval <- sort.int(unique(retval))
    assign(as.character(nr), retval, envir = OUTENV)
    return(retval)
  }

  # run over each node
  for (i in seq_len(nr)) {
    out.idx[[i]] <- get_ancestors(i, integer(0));
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
  A <- adj.mat
  n <- nrow(A)
  k <- 0
  testen <- 1:n
  while (TRUE) {
    found <- FALSE
    tests <- testen # which expressions are tested in the next loop
    for (i in tests) {
      if (all(A[1:n, i] == 0)) { # var defined in expression i doesn't use
        # vars not already defined, select as next
        k <- k + 1            # increment order counter
        A[i, ] <- 0           # remove usage indicator of defined variable
        A[i, i] <- k          # keep order of expressions in diagonal element
        found <- TRUE
        testen <- testen[testen != i] # don't test this i in next loop
      }
    }
    if (k == n) return (order(diag(A)))  # all done
    if (!found) {  # no definable var found or all done
        if (warn) {
          lav_msg_warn(gettext("unable to sort `:=` parameters;",
                               "system of defined parameters contains a cycle"))
        }
        return(seq_len(n))    # cycle detected; return original order
    }
  }
}

