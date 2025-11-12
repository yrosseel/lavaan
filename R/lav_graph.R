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
        adjmat[i, ] <- 0           # remove usage indicator of defined variable
        adjmat[i, i] <- k          # keep order of expressions in diagonal element
        found <- TRUE
        testen <- testen[testen != i] # don't test this i in next loop
      }
    }
    if (k == n) return (order(diag(adjmat)))  # all done
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
  # this function is faster, requires less memory and handles also character vectors
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
# This routine returns a data.frame with the nodes and their position in a
# matrix (rows, cols), as follows :
# the first column contains all nodes without dependencies
# the second column contains nodes with only dependency in the first column (*)
# the third column contains nodes with only dependencies in the first and
#     second, with at least one in the second column (*)
# and so on 
# (*) all nodes which have no successors are moved to the last column
# the rows in the second to last column are chosen so that they are in the 
# neighborhood of the mean of the rows from the nodes they depend on.
lav_graph_topological_matrix <- function(defined, definedby, warn = TRUE) {
  rv <- lav_graph_topological_sort(defined, definedby, warn)
  n <- length(rv)
  rvrow <- integer(n)
  rvcol <-  integer(n)
  colmax <- 0L
  lastcolmax <- 0L
  for (i in seq.int(n)) {
    if (i == 1L) {
      rvrow[i] <- rvcol[i] <- 1L 
      colmax <- 1L
    } else {
      predecessors <- definedby[rv[i] == defined]  
      followers <- defined[rv[i] == definedby]
      if (length(predecessors) == 0L) {
        rvcol[i] <- 1L
        colmax[1L] <- colmax[1L] + 1L
        rvrow[i] <- colmax[1L]
      } else {
        if (length(followers) == 0) { # put in last column, temporarily -1
          rvcol[i] <- -1L
          lastcolmax <- lastcolmax + 1L
          rvrow[i] <- lastcolmax
        } else {
          predecessor.ind <- match(predecessors, rv)
          rvcol[i] <- max(rvcol[predecessor.ind]) + 1L
          if (length(colmax) < rvcol[i]) {
            colmax[rvcol[i]] <- 1L
          } else {
            colmax[rvcol[i]] <- colmax[rvcol[i]] + 1L
          }
        }
        rvrow[i] <- colmax[rvcol[i]]
      }     
    }
  }
  # last column
  rvcol[rvcol == -1L] <- length(colmax) + 1L
  colmax[length(colmax) + 1L] <- lastcolmax
  # adapt rows in columns to match as close as possible the rows of the
  # connected nodes in the previous column(s) :
  if (length(colmax) > 1L) {
    rowmax <- max(colmax)
    for (c in seq.int(2L, length(colmax))) {
      cnodes.ind <- which(rvcol == c)
      optimalrows <- sapply(cnodes.ind, function (ci) {
        prevnodes.ind <- which(rv[ci] == defined)
        prevnodes <- definedby[prevnodes.ind]
        mean(rvrow[match(prevnodes, rv)])
      })
      or_order <- order(optimalrows)
      optimalrows <- optimalrows[or_order]
      for (j in seq_along(optimalrows)) {
        if (j == 1L) {
          optimalrows[j] <- as.integer(optimalrows[j]) 
        } else {
          if (as.integer(optimalrows[j]) <= optimalrows[j - 1L]) {
            optimalrows[j] <- optimalrows[j - 1L] + 1L
          } else {
            optimalrows[j] <- as.integer(optimalrows[j]) 
          }
        }
      }
      optimalrows[or_order] <- optimalrows
      rvrow[cnodes.ind] <- optimalrows
    }
  }
  # order on column, then row, and return  
  neworder <- order(rvcol * 1000L + rvrow)
  data.frame(nodes = rv[neworder], 
             rows = rvrow[neworder], 
             cols = rvcol[neworder])
}