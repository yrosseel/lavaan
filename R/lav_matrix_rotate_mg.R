# rotation algorithms for multiple groups
#
# YR 22 Dec 2025 -- initial version
#
# - rstarts only affects the 'initial' rotation

# main function to rotate a list of matrices 'Alist'
lav_matrix_rotate_mg <- function(Alist = NULL, # original matrices
                                 orthogonal = FALSE, # default is oblique
                                 method = "geomin", # default rot method
                                 method.args = list(
                                   geomin.epsilon = 0.01,
                                   orthomax.gamma = 1,
                                   cf.gamma = 0,
                                   oblimin.gamma = 0,
                                   promax.kappa = 4,
                                   target = matrix(0, 0, 0),
                                   target.mask = matrix(0, 0, 0)
                                 ),
                                 init.rotList = NULL, # initial rotation matrix
                                 init.ROT.check = TRUE, # check if init ROT is ok
                                 rstarts = 30L, # number of random starts
                                 # row.weights = "default", # row weighting
                                 # std.ov = FALSE, # rescale ov
                                 # ov.var = NULL, # ov variances
                                 # algorithm = "gpa", # rotation algorithm
                                 mg.algorithm = "pairwise",
                                 mg.agreement.method = "procrustes",
                                 mg.agreement.weight = 0.5,
                                 reflect = TRUE, # refect sign
                                 order.lv.by = "index", # how to order the lv's
                                 # gpa.tol = 0.00001, # stopping tol gpa
                                 tol = 1e-07, # stopping tol others
                                 keep.rep = FALSE, # store replications
                                 max.iter = 10000L) { # max iterations

  # check Alist
  if (!is.list(Alist)) {
    lav_msg_stop(gettext("Alist does not seem to be a list"))
  }
  ngroups <- length(Alist)

  # check A
  A <- Alist[[1]]
  if (!inherits(A, "matrix")) {
    lav_msg_stop(gettext("A does not seem to be a matrix"))
  }

  P <- nrow(A)
  M <- ncol(A)
  if (M < 2L) { # single dimension
    res <- list(
      lambdaList = Alist, rotList = NULL,
      orthogonal = orthogonal, method = "none",
      method.args = list(), row.weights = "none",
      algorithm = "none", iter = 0L, converged = TRUE,
      method.value = 0
    )
    return(res)
  }

  # method
  method <- tolower(method)

  # no promax (for now)
  if (method == "promax") {
    lav_msg_stop(gettext("mg-promax + agreement not supported"))
  }

  # check init.ROT per group
  if (!is.null(init.rotList) && init.ROT.check) {
    for(g in seq_len(ngroups)) {
      init.ROT <- init.rotList[[g]]
      if (!inherits(init.ROT, "matrix")) {
        lav_msg_stop(gettext("init.ROT does not seem to a matrix"))
      }
      if (nrow(init.ROT) != M) {
        lav_msg_stop(gettextf(
          "nrow(init.ROT) = %1$s does not equal ncol(A) = %2$s",
          nrow(init.ROT), M
        ))
      }
      if (nrow(init.ROT) != ncol(init.ROT)) {
        lav_msg_stop(gettextf(
          "nrow(init.ROT) = %1$s does not equal ncol(init.ROT) = %2$s",
          nrow(init.ROT), ncol(init.ROT)
        ))
      }
      # rotation matrix?
      if (!lav_matrix_rotate_check(init.ROT, orthogonal = orthogonal)) {
        lav_msg_stop(gettext("init.ROT does not look like a rotation matrix"))
      }
    } # group
  } # check init.rotList

  # determine method function name
  if (method %in% c(
    "cf-quartimax", "cf-varimax", "cf-equamax",
    "cf-parsimax", "cf-facparsim"
  )) {
    method.fname <- "lav_matrix_rotate_cf"
    method.args$cf.gamma <- switch(method,
      "cf-quartimax" = 0,
      "cf-varimax"   = 1 / P,
      "cf-equamax"   = M / (2 * P),
      "cf-parsimax"  = (M - 1) / (P + M - 2),
      "cf-facparsim" = 1
    )
  } else if (method %in% c("bi-quartimin", "biquartimin")) {
    method.fname <- "lav_matrix_rotate_biquartimin"
  } else if (method %in% c("bi-geomin", "bigeomin")) {
    method.fname <- "lav_matrix_rotate_bigeomin"
  } else if (method == "target.strict") {
    method.fname <- "lav_matrix_rotate_target"
  } else {
    method.fname <- paste("lav_matrix_rotate_", method, sep = "")
  }

  # check if rotation method exists
  check <- try(get(method.fname), silent = TRUE)
  if (inherits(check, "try-error")) {
    lav_msg_stop(gettext("unknown rotation method:"), method.fname)
  }

  # if target, check target matrix
  if (method == "target.strict" || method == "pst") {
    target <- method.args$target
    if (is.list(target)) {
      method.args$target <- target <- target[[1L]] # always take the first group
    }
    # check dimension of target/A
    if (nrow(target) != nrow(A)) {
      lav_msg_stop(gettext("nrow(target) != nrow(A)"))
    }
    if (ncol(target) != ncol(A)) {
      lav_msg_stop(gettext("ncol(target) != ncol(A)"))
    }
  }
  if (method == "pst") {
    target.mask <- method.args$target.mask
    if (is.list(target.mask)) {
      method.args$target.mask <- target.mask <- target.mask[[1L]]
    }
    # check dimension of target.mask/A
    if (nrow(target.mask) != nrow(A)) {
      lav_msg_stop(gettext("nrow(target.mask) != nrow(A)"))
    }
    if (ncol(target.mask) != ncol(A)) {
      lav_msg_stop(gettext("col(target.mask) != ncol(A)"))
    }
  }
  # we keep this here, so lav_matrix_rotate() can be used independently
  if (method == "target.strict" && anyNA(target)) {
    method <- "pst"
    method.fname <- "lav_matrix_rotate_pst"
    target.mask <- matrix(1, nrow = nrow(target), ncol = ncol(target))
    target.mask[is.na(target)] <- 0
    method.args$target.mask <- target.mask
  }

  # set orthogonal option
  if (missing(orthogonal)) {
    # the default is oblique, except for varimax, entropy and a few others
    if (method %in% c(
      "varimax", "entropy", "mccammon",
      "tandem1", "tandem2"
    )) {
      orthogonal <- TRUE
    } else {
      orthogonal <- FALSE
    }
  } else {
    if (!orthogonal && method %in% c(
      "varimax", "entropy", "mccammon",
      "tandem1", "tandem2"
    )) {
      lav_msg_warn(gettextf(
        "rotation method %s may not work with oblique rotation.",
        dQuote(method)
      ))
    }
  }

  # 0. initialize: rotate each group separately + reorder
  Alist.orig <- Alist
  orderList <- vector("list", length = ngroups)
  for (g in seq_len(ngroups)) {
    if (lav_verbose()) {
      cat("  group = ", g, "\n")
    }
    res0 <- lav_matrix_rotate(
      A = Alist[[g]],
      orthogonal = orthogonal,
      method = method,
      method.args = method.args,
      init.ROT = NULL,
      init.ROT.check = FALSE,
      rstarts = rstarts,
      row.weights = "none",
      std.ov = FALSE,
      ov.var = NULL,
      algorithm = "pairwise",
      reflect = reflect,
      order.lv.by = order.lv.by,
      tol = tol,
      max.iter = max.iter,
      group = g
    )
    if (g > 1L) {
      # reorder factors to align with the first group
      order.idx <- lav_efa_find_best_order(
        lambda_ref = Alist[[1]],
        lambda_target = res0$LAMBDA
      )
      Alist[[g]] <- res0$LAMBDA[, order.idx, drop = FALSE]
      orderList[[g]] <- order.idx
    } else {
      Alist[[g]] <- res0$LAMBDA
      orderList[[g]] <- res0$order.idx
    }
  } # group

  #  # set row.weights
  #  row.weights <- tolower(row.weights)
  #  if (row.weights == "default") {
  #    # the default is "none", except for varimax
  #    if (method %in% c("varimax", "promax")) {
  #      row.weights <- "kaiser"
  #    } else {
  #      row.weights <- "none"
  #    }
  #  }

  #   # check algorithm
  #   algorithm <- tolower(algorithm)
  #   if (algorithm %in% c("gpa", "pairwise", "none")) {
  #     # nothing to do
  #   } else {
  #     lav_msg_stop(gettext("algorithm must be gpa or pairwise"))
  #   }

  # 1. compute row weigths

  #   # 1.a cov -> cor?
  #   if (std.ov) {
  #     A <- A * 1 / sqrt(ov.var)
  #   }
  #
  #   if (row.weights == "none") {
  #     weights <- rep(1.0, P)
  #   } else if (row.weights == "kaiser") {
  #     weights <- lav_matrix_rotate_kaiser_weights(A)
  #   } else if (row.weights == "cureton-mulaik") {
  #     weights <- lav_matrix_rotate_cm_weights(A)
  #   } else {
  #     lav_msg_stop(gettext("row.weights can be none, kaiser or cureton-mulaik"))
  #   }
  #   A <- A * weights


  # 2. rotate

  # multiple random starts?
#  if (rstarts > 0L) {
#    REP <- lapply(seq_len(rstarts), function(rep) {
#      # random start (always orthogonal)
#      # init.ROT <- lav_matrix_rotate_gen(M = M, orthogonal = TRUE)
#      # init.ROT <- lav_matrix_rotate_gen(M = M, orthogonal = orthogonal)
#
#      if (lav_verbose()) {
#        cat("\n")
#        cat("rstart = ", sprintf("%4d", rep), " start:\n")
#      }
#
#      # mg-pairwise + agreement
#      rotList <- lav_matrix_rotate_pairwise_mg(
#        Alist = Alist,
#        orthogonal = orthogonal,
#        random.ROT = TRUE,
#        method.fname = method.fname,
#        method.args = method.args,
#        tol = tol,
#        max.iter = max.iter
#      )
#      info <- attr(rotList, "info")
#      attr(rotList, "info") <- NULL
#      res <- list(value = info$method.value, rotList = rotList)
#
#      if (lav_verbose()) {
#        cat(
#          "rstart = ", sprintf("%4d", rep),
#          " end; current crit = ", sprintf("%17.15f", res$value), "\n"
#        )
#      }
#      res
#    })
#    best.idx <- which.min(sapply(REP, "[[", "value"))
#    rotList <- REP[[best.idx]]$rotList
#    if (keep.rep) {
#      info <- list(method.value = REP[[best.idx]]$value, REP = REP)
#    } else {
#      info <- list(method.value = REP[[best.idx]]$value)
#    }
#  } else { # just a single rstart
    # initial rotation matrix
    # if (is.null(init.ROT)) {
    #   init.ROT <- diag(M)
    # }

    if (lav_verbose()) {
      cat("  mg-pairwise + agreement rotation:\n")
    }
    rotList <- lav_matrix_rotate_pairwise_mg(
      Alist = Alist,
      orthogonal = orthogonal,
      random.ROT = FALSE,
      method.fname = method.fname,
      method.args = method.args,
      mg.agreement.method = mg.agreement.method,
      mg.agreement.weight = mg.agreement.weight,
      tol = tol,
      max.iter = max.iter
    )
    info <- attr(rotList, "info")
    attr(rotList, "info") <- NULL
  #}
  lambdaList <- vector("list", length = ngroups)
  for (g in seq_len(ngroups)) {
    A <- Alist[[g]]
    ROT <- rotList[[g]]

    if (orthogonal) {
      # LAMBDA <- A %*% solve(t(ROT))
      # note: when ROT is orthogonal, solve(t(ROT)) == ROT
      LAMBDA <- A %*% ROT
    } else {
      # LAMBDA <- A %*% solve(t(ROT))
      LAMBDA <- t(solve(ROT, t(A)))
    }

    #  # 3. undo row weighting
    #  LAMBDA <- LAMBDA / weights

    #  # 3.b undo cov -> cor
    #  if (std.ov) {
    #    LAMBDA <- LAMBDA * sqrt(ov.var)
    #  }

    # put back in List
    lambdaList[[g]] <- LAMBDA
  } # group

  # re-compute final rotation matrix, using Alist.orig
  # because in step 0, we already rotated (per group)
  rotList <- vector("list", length = ngroups)
  if (orthogonal) {
    for (g in seq_len(ngroups)) {
      rotList[[g]] <- solve(
        crossprod(Alist.orig[[g]]),
        crossprod(Alist.orig[[g]], lambdaList[[g]])
      )
    }
  } else {
    for (g in seq_len(ngroups)) {
      # to be compatible with GPa
      ROTt.inv <- solve(
        crossprod(Alist.orig[[g]]),
        crossprod(Alist.orig[[g]], lambdaList[[g]])
      )
      rotList[[g]] <- solve(t(ROTt.inv))
    }
  }

  # 6. return results as a list
  res <- list(
    lambdaList = lambdaList, rotList = rotList,
    orderList = orderList,
    orthogonal = orthogonal, method = method,
    method.args = method.args, row.weights = "none"
  )

  # add method info
  res <- c(res, info)

  res
}


# pairwise rotation algorithm with direct line search
# (see lav_matrix_rotate_pairwise()
#
# but adapted to handle multiple groups + an agreement criterion
#
lav_matrix_rotate_pairwise_mg <- function(Alist = NULL, # original matrices
                                          orthogonal = FALSE,
                                          random.ROT = FALSE,
                                          method.fname = NULL, # crit function
                                          method.args = list(), # method args
                                          tol = 1e-8,
                                          mg.agreement.weight = 0.50,
                                          mg.agreement.method = "procrustes",
                                          scale = TRUE,
                                          max.iter = 1000L) {
  # first group is 'A'
  A <- Alist[[1]]
  ngroups <- length(Alist)

  # number of columns
  M <- ncol(A)

  # initial LAMBDA + PSI
  lambdaList <- Alist
  psiList <- rep(list(diag(M)), ngroups)

  # random initial rotation?
  if (random.ROT) {
    for (g in seq_len(ngroups)) {
      init.ROT <- lav_matrix_rotate_gen(M = M, orthogonal = TRUE)
      if (orthogonal) {
        lambdaList[[g]] <- Alist[[g]] %*% init.ROT
      } else {
        lambdaList[[g]] <- t(solve(init.ROT, t(Alist[[g]])))
        psiList[[g]] <- crossprod(init.ROT)
      }
    }
  }

  # using the current LAMBDA, evaluate the user-specified
  # rotation criteron; return Q (the criterion) only
  Q.current <- lav_matrix_rotate_mg_agreement(
    lambdaList = lambdaList, method.fname = method.fname,
    method.args = method.args,
    mg.agreement.method = mg.agreement.method,
    w = mg.agreement.weight, scale = scale
  )

  # if verbose, print
  if (lav_verbose()) {
    Q.mg <- attr(Q.current, "Q.mg")
    A.mg <- attr(Q.current, "A.mg")
    cat(
      "  iter = ", sprintf("%4d", 0),
      " Q = ", sprintf("%13.11f", Q.current),
      " Q.rot = ", sprintf("%13.11f", Q.mg),
      " Q.agr = ", sprintf("%13.11f", A.mg), "\n"
    )
  }

  # plane combinations
  if (orthogonal) {
    PLANE <- utils::combn(M, 2)
  } else {
    tmp <- utils::combn(M, 2)
    PLANE <- cbind(tmp, tmp[c(2, 1), , drop = FALSE])
  }

  # define objective function -- orthogonal
  objf_orth <- function(g = 1L, theta = 0, A = NULL, col1 = 0L, col2 = 0L) {
    # construct ROT
    ROT <- diag(M)
    ROT[col1, col1] <- base::cos(theta)
    ROT[col1, col2] <- base::sin(theta)
    ROT[col2, col1] <- -1 * base::sin(theta)
    ROT[col2, col2] <- base::cos(theta)

    # rotate
    lambdaList[[g]] <- A %*% ROT

    # evaluate criterion
    Q <- lav_matrix_rotate_mg_agreement(
      lambdaList = lambdaList, method.fname = method.fname,
      method.args = method.args, mg.agreement.method = mg.agreement.method,
      w = mg.agreement.weight, scale = scale
    )

    Q
  }

  # define objective function -- oblique
  objf_obliq <- function(g = 1, delta = 0, A = NULL, col1 = 0L, col2 = 0L,
                         psi12 = 0) {
    # construct ROT
    ROT <- diag(M)

    # gamma
    gamma2 <- 1 + (2 * delta * psi12) + (delta * delta)

    ROT[col1, col1] <- sqrt(abs(gamma2))
    ROT[col1, col2] <- -1 * delta
    ROT[col2, col1] <- 0
    ROT[col2, col2] <- 1

    # rotate
    lambdaList[[g]] <- A %*% ROT

    # evaluate criterion
    Q <- lav_matrix_rotate_mg_agreement(
      lambdaList = lambdaList, method.fname = method.fname,
      method.args = method.args, mg.agreement.method = mg.agreement.method,
      w = mg.agreement.weight, scale = scale
    )

    Q
  }

  # start iterations
  converged <- FALSE
  Q.old <- Q.current
  for (iter in seq_len(max.iter)) {
    # for each group
    for (g in seq_len(ngroups)) {
      # rotate - one cycle
      for (pl in seq_len(ncol(PLANE))) {
        # choose plane
        col1 <- PLANE[1, pl]
        col2 <- PLANE[2, pl]

        # optimize
        if (orthogonal) {
          out <- optimize(
            f = objf_orth, interval = c(-pi / 4, +pi / 4),
            g = g,
            A = lambdaList[[g]], col1 = col1, col2 = col2,
            maximum = FALSE, tol = .Machine$double.eps^0.25
          )
          # best rotation - for this plane
          theta <- out$minimum

          # construct ROT
          ROT <- diag(M)
          ROT[col1, col1] <- base::cos(theta)
          ROT[col1, col2] <- base::sin(theta)
          ROT[col2, col1] <- -1 * base::sin(theta)
          ROT[col2, col2] <- base::cos(theta)
        } else {
          psi12 <- psiList[[g]][col1, col2]
          out <- optimize(
            f = objf_obliq, interval = c(-1, +1),
            g = g,
            A = lambdaList[[g]], col1 = col1, col2 = col2,
            psi12 = psi12,
            maximum = FALSE, tol = .Machine$double.eps^0.25
          )

          # best rotation - for this plane
          delta <- out$minimum

          # construct ROT
          ROT <- diag(M)

          # gamma
          gamma2 <- 1 + (2 * delta * psi12) + (delta * delta)
          gamma <- sqrt(abs(gamma2))

          ROT[col1, col1] <- gamma
          ROT[col1, col2] <- -1 * delta
          ROT[col2, col1] <- 0
          ROT[col2, col2] <- 1
        }

        # rotate
        lambdaList[[g]] <- lambdaList[[g]] %*% ROT

        if (!orthogonal) {
          PSI <- psiList[[g]]
          # rotate PSI
          PSI[col1, ] <- (1 / gamma) * PSI[col1, ] + (delta / gamma) * PSI[col2, ]
          PSI[, col1] <- PSI[col1, ]
          PSI[col1, col1] <- 1
          psiList[[g]] <- PSI
        }
      } # all planes
    } # ngroups

    # check for convergence
    Q.current <- lav_matrix_rotate_mg_agreement(
      lambdaList = lambdaList, method.fname = method.fname,
      method.args = method.args, mg.agreement.method = mg.agreement.method,
      w = mg.agreement.weight, scale = scale
    )

    # absolute change in Q
    diff <- abs(Q.old - Q.current)

    # if verbose, print
    if (lav_verbose()) {
      Q.mg <- attr(Q.current, "Q.mg")
      A.mg <- attr(Q.current, "A.mg")
      cat(
        "  iter = ", sprintf("%4d", iter),
        " Q = ", sprintf("%13.11f", Q.current),
        " Q.rot = ", sprintf("%13.11f", Q.mg),
        " Q.agr = ", sprintf("%13.11f", A.mg),
        " change = ", sprintf("%13.11f", diff), "\n"
      )
    }

    if (!is.finite(diff)) {
      # something went wrong ... bail out
      lav_msg_warn(gettextf("mg-pairwise rotation algorithm failed!"))
      break
    } else if (diff < tol) {
      converged <- TRUE
      break
    } else {
      Q.old <- Q.current
    }
  } # iter

  # warn if no convergence
  if (!converged) {
    lav_msg_warn(gettextf(
      "mg-pairwise rotation algorithm did not converge after %s iterations",
      max.iter
    ))
  }

  # compute final rotation matrix
  rotList <- vector("list", length = ngroups)
  if (orthogonal) {
    for (g in seq_len(ngroups)) {
      rotList[[g]] <- solve(
        crossprod(Alist[[g]]),
        crossprod(Alist[[g]], lambdaList[[g]])
      )
    }
  } else {
    for (g in seq_len(ngroups)) {
      # to be compatible with GPa
      ROTt.inv <- solve(
        crossprod(Alist[[g]]),
        crossprod(Alist[[g]], lambdaList[[g]])
      )
      rotList[[g]] <- solve(t(ROTt.inv))
    }
  }

  # algorithm information
  info <- list(
    algorithm = "mg-pairwise",
    iter = iter,
    converged = converged,
    method.value = Q.current
  )

  attr(rotList, "info") <- info

  rotList
}
