# rotation algorithms
#
# YR  3 April 2019 -- gradient projection algorithm
# YR 21 April 2019 -- pairwise rotation algorithm
# YR 11 May   2020 -- order.idx is done in rotation matrix
#                     (suggested by Florian Scharf)
# YR 02 June  2024 -- add group argument, so target and target_mask can
#                     be a list

# main function to rotate a single matrix 'A'
lav_matrix_rotate <- function(a = NULL, # original matrix
                              orthogonal = FALSE, # default is oblique
                              method = "geomin", # default rot method
                              method_args = list(
                                geomin_epsilon = 0.01,
                                orthomax_gamma = 1,
                                cf_gamma = 0,
                                oblimin_gamma = 0,
                                promax_kappa = 4,
                                target = matrix(0, 0, 0),
                                target_mask = matrix(0, 0, 0)
                              ),
                              init_rot = NULL, # initial rotation matrix
                              init_rot_check = TRUE, # check if init ROT is ok
                              rstarts = 100L, # number of random starts
                              row_weights = "default", # row weighting
                              std_ov = FALSE, # rescale ov
                              ov_var = NULL, # ov variances
                              algorithm = "gpa", # rotation algorithm
                              reflect = TRUE, # reflect sign
                              order_lv_by = "index", # how to order the lv's
                              gpa_tol = 0.00001, # stopping tol gpa
                              tol = 1e-07, # stopping tol others
                              keep_rep = FALSE, # store replications
                              max_iter = 10000L, # max gpa iterations
                group = 1L) { # group number

  # check A
  if (!inherits(a, "matrix")) {
    lav_msg_stop(gettext("A does not seem to be a matrix"))
  }

  p <- nrow(a)
  m <- ncol(a)
  if (m < 2L) { # single dimension
    res <- list(
      mm_lambda = a, PHI = matrix(1, 1, 1), ROT = matrix(1, 1, 1),
      orthogonal = orthogonal, method = "none",
      method.args = list(), row_weights = "none",
      algorithm = "none", iter = 0L, converged = TRUE,
      method.value = 0
    )
    return(res)
  }

  # method
  method <- tolower(method)

  # if promax, skip everything, then call promax() later
  if (method == "promax") {
    # orig.algorithm <- algorithm
    # orig.rstarts <- rstarts

    algorithm <- "none"
    rstarts <- 0L
    init_rot <- NULL
    rot <- diag(m)
  }

  # check init.ROT
  if (!is.null(init_rot) && init_rot_check) {
    if (!inherits(init_rot, "matrix")) {
      lav_msg_stop(gettext("init.ROT does not seem to be a matrix"))
    }
    if (nrow(init_rot) != m) {
      lav_msg_stop(gettextf(
        "nrow(init.ROT) = %1$s does not equal ncol(A) = %2$s",
        nrow(init_rot), m))
    }
    if (nrow(init_rot) != ncol(init_rot)) {
      lav_msg_stop(gettextf(
        "nrow(init.ROT) = %1$s does not equal ncol(init.ROT) = %2$s",
        nrow(init_rot), ncol(init_rot)))
    }
    # rotation matrix?
    if (!lav_matrix_rotate_check(init_rot, orthogonal = orthogonal)) {
      lav_msg_stop(gettext("init.ROT does not look like a rotation matrix"))
    }
  }

  # determine method function name
  if (method %in% c(
    "cf-quartimax", "cf-varimax", "cf-equamax",
    "cf-parsimax", "cf-facparsim"
  )) {
    method_fname <- "lav_matrix_rotate_cf"
    method_args$cf_gamma <- switch(method,
      "cf-quartimax" = 0,
      "cf-varimax"   = 1 / p,
      "cf-equamax"   = m / (2 * p),
      "cf-parsimax"  = (m - 1) / (p + m - 2),
      "cf-facparsim" = 1
    )
  } else if (method %in% c("bi-quartimin", "biquartimin")) {
    method_fname <- "lav_matrix_rotate_biquartimin"
  } else if (method %in% c("bi-geomin", "bigeomin")) {
    method_fname <- "lav_matrix_rotate_bigeomin"
  } else if (method == "target.strict") {
    method_fname <- "lav_matrix_rotate_target"
  } else {
    method_fname <- paste("lav_matrix_rotate_", method, sep = "")
  }

  # check if rotation method exists
  check <- try(get(method_fname), silent = TRUE)
  if (inherits(check, "try-error")) {
    lav_msg_stop(gettext("unknown rotation method:"), method_fname)
  }

  # if target, check target matrix
  if (method == "target.strict" || method == "pst") {
    target <- method_args$target
  if (is.list(target)) {
    method_args$target <- target <- target[[group]]
  }
    # check dimension of target/A
    if (nrow(target) != nrow(a)) {
      lav_msg_stop(gettext("nrow(target) != nrow(A)"))
    }
    if (ncol(target) != ncol(a)) {
      lav_msg_stop(gettext("ncol(target) != ncol(A)"))
    }
  }
  if (method == "pst") {
    target_mask <- method_args$target_mask
  if (is.list(target_mask)) {
    method_args$target_mask <- target_mask <- target_mask[[group]]
  }
    # check dimension of target_mask/A
    if (nrow(target_mask) != nrow(a)) {
      lav_msg_stop(gettext("nrow(target_mask) != nrow(A)"))
    }
    if (ncol(target_mask) != ncol(a)) {
      lav_msg_stop(gettext("col(target_mask) != ncol(A)"))
    }
  }
  # we keep this here, so lav_matrix_rotate() can be used independently
  if (method == "target.strict" && anyNA(target)) {
    method <- "pst"
    method_fname <- "lav_matrix_rotate_pst"
    target_mask <- matrix(1, nrow = nrow(target), ncol = ncol(target))
    target_mask[is.na(target)] <- 0
    method_args$target_mask <- target_mask
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

  # set row_weights
  row_weights <- tolower(row_weights)
  if (row_weights == "default") {
    # the default is "none", except for varimax
    if (method %in% c("varimax", "promax")) {
      row_weights <- "kaiser"
    } else {
      row_weights <- "none"
    }
  }

  # check algorithm
  algorithm <- tolower(algorithm)
  if (algorithm %in% c("gpa", "pairwise", "none")) {
    # nothing to do
  } else {
    lav_msg_stop(gettext("algorithm must be gpa or pairwise"))
  }



  # 1. compute row weigths

  # 1.a cov -> cor?
  if (std_ov) {
    a <- a * 1 / sqrt(ov_var)
  }

  if (row_weights == "none") {
    weights <- rep(1.0, p)
  } else if (row_weights == "kaiser") {
    weights <- lav_matrix_rotate_kaiser_weights(a)
  } else if (row_weights == "cureton-mulaik") {
    weights <- lav_matrix_rotate_cm_weights(a)
  } else {
    lav_msg_stop(gettext("row_weights can be none, kaiser or cureton-mulaik"))
  }
  a <- a * weights


  # 2. rotate

  # multiple random starts?
  if (rstarts > 0L) {
    rep_1 <- sapply(seq_len(rstarts), function(rep) {
      # random start (always orthogonal)
      init_rot <- lav_matrix_rotate_gen(m = m, orthogonal = TRUE)
      # init.ROT <- lav_matrix_rotate_gen(M = M, orthogonal = orthogonal)

      if (lav_verbose()) {
        cat("\n")
        cat("rstart = ", sprintf("%4d", rep), " start:\n")
      }


      # choose rotation algorithm
      if (algorithm == "gpa") {
        rot <- lav_matrix_rotate_gpa(
          a = a, orthogonal = orthogonal,
          init_rot = init_rot,
          method_fname = method_fname,
          method_args = method_args,
          gpa_tol = gpa_tol,
          max_iter = max_iter
        )
        info <- attr(rot, "info")
        attr(rot, "info") <- NULL
        res <- c(info$method.value, lav_matrix_vec(rot))
      } else if (algorithm == "pairwise") {
        rot <- lav_matrix_rotate_pairwise(
          a = a,
          orthogonal = orthogonal,
          init_rot = init_rot,
          method_fname = method_fname,
          method_args = method_args,
          tol = tol,
          max_iter = max_iter
        )
        info <- attr(rot, "info")
        attr(rot, "info") <- NULL
        res <- c(info$method.value, lav_matrix_vec(rot))
      }

      if (lav_verbose()) {
        cat(
          "rstart = ", sprintf("%4d", rep),
          " end; current crit = ", sprintf("%17.15f", res[1]), "\n"
        )
      }
      res
    })
    best_idx <- which.min(rep_1[1, ])
    rot <- matrix(rep_1[-1, best_idx], nrow = m, ncol = m)
    if (keep_rep) {
      info <- list(method.value = rep_1[1, best_idx], REP = rep_1)
    } else {
      info <- list(method.value = rep_1[1, best_idx])
    }
  } else if (algorithm != "none") {
    # initial rotation matrix
    if (is.null(init_rot)) {
      init_rot <- diag(m)
    }

    # Gradient Projection Algorithm
    if (algorithm == "gpa") {
      rot <- lav_matrix_rotate_gpa(
        a = a, orthogonal = orthogonal,
        init_rot = init_rot,
        method_fname = method_fname,
        method_args = method_args,
        gpa_tol = gpa_tol,
        max_iter = max_iter
      )
    } else if (algorithm == "pairwise") {
      rot <- lav_matrix_rotate_pairwise(
        a = a,
        orthogonal = orthogonal,
        init_rot = init_rot,
        method_fname = method_fname,
        method_args = method_args,
        tol = tol,
        max_iter = max_iter
      )
    }
    info <- attr(rot, "info")
    attr(rot, "info") <- NULL
  }

  # final rotation
  if (orthogonal) {
    # mm_lambda <- A %*% solve(t(ROT))
    # note: when ROT is orthogonal, solve(t(ROT)) == ROT
    mm_lambda <- a %*% rot
    phi <- diag(ncol(mm_lambda)) # correlation matrix == I
  } else {
    # mm_lambda <- A %*% solve(t(ROT))
    mm_lambda <- t(solve(rot, t(a)))
    phi <- crossprod(rot) # correlation matrix
  }

  # 3. undo row weighting
  mm_lambda <- mm_lambda / weights

  # here, after re-weighted, we run promax if needed
  if (method == "promax") {
    lambda_orig <- mm_lambda

    # first, run 'classic' varimax using varimax() from the stats package
    # we split varimax from promax, so we can control the normalize flag
    normalize_flag <- row_weights == "kaiser"
    xx <- stats::varimax(x = mm_lambda, normalize = normalize_flag)

    # promax
    kappa <- method_args$promax_kappa
    out <- lav_matrix_rotate_promax(
      x = xx$loadings, m = kappa,
      varimax_rot = xx$rotmat
    )
    mm_lambda <- out$loadings
    phi <- solve(crossprod(out$rotmat))

    # compute 'ROT' to be compatible with GPa
    rott_inv <- solve(
      crossprod(lambda_orig),
      crossprod(lambda_orig, mm_lambda)
    )
    rot <- solve(t(rott_inv))

    info <- list(
      algorithm = "promax", iter = 0L, converged = TRUE,
      method.value = as.numeric(NA)
    )
  }

  # 3.b undo cov -> cor
  if (std_ov) {
    mm_lambda <- mm_lambda * sqrt(ov_var)
  }

  # 4.a reflect so that column sum is always positive
  if (reflect) {
    sum_1 <- colSums(mm_lambda)
    neg_idx <- which(sum_1 < 0)
    if (length(neg_idx) > 0L) {
      mm_lambda[, neg_idx] <- -1 * mm_lambda[, neg_idx, drop = FALSE]
      rot[, neg_idx] <- -1 * rot[, neg_idx, drop = FALSE]
      if (!orthogonal) {
        # recompute PHI
        phi <- crossprod(rot)
      }
    }
  }

  # 4.b reorder the columns
  if (order_lv_by == "sumofsquares") {
    l2 <- mm_lambda * mm_lambda
    order_idx <- base::order(colSums(l2), decreasing = TRUE)
  } else if (order_lv_by == "index") {
    # reorder using Asparouhov & Muthen 2009 criterion (see Appendix D)
    max_loading <- apply(abs(mm_lambda), 2, max)
    # 1: per factor, number of the loadings that are at least 0.8 of the
    #    highest loading of the factor
    # 2: mean of the index numbers
    average_index <- sapply(seq_len(ncol(mm_lambda)), function(i) {
      mean(which(abs(mm_lambda[, i]) >= 0.8 * max_loading[i]))
    })
    # order of the factors
    order_idx <- base::order(average_index)
  } else if (order_lv_by == "none") {
    order_idx <- seq_len(ncol(mm_lambda))
  } else {
    lav_msg_stop(gettext("order must be index, sumofsquares or none"))
  }

  # do the same in PHI
  mm_lambda <- mm_lambda[, order_idx, drop = FALSE]
  phi <- phi[order_idx, order_idx, drop = FALSE]

  # new in 0.6-6, also do this in ROT, so we won't have to do this
  # again upstream
  rot <- rot[, order_idx, drop = FALSE]

  # 6. return results as a list
  res <- list(
    mm_lambda = mm_lambda, PHI = phi, ROT = rot, order.idx = order_idx,
    orthogonal = orthogonal, method = method,
    method_args = method_args, row_weights = row_weights
  )

  # add method info
  res <- c(res, info)

  res
}


# Gradient Projection Algorithm (Jennrich 2001, 2002)
#
# - this is a translation of the SAS PROC IML code presented in the Appendix
#   of Bernaards & Jennrich (2005)
# - as the orthogonal and oblique algorithm are so similar, they are
#   combined in a single function
# - the default is oblique rotation
#
lav_matrix_rotate_gpa <- function(a = NULL, # original matrix
                                  orthogonal = FALSE, # default is oblique
                                  init_rot = NULL, # initial rotation
                                  method_fname = NULL, # criterion function
                                  method_args = list(), # optional method args
                                  gpa_tol = 0.00001,
                                  max_iter = 10000L) {
  # number of columns
  m <- ncol(a)

  # transpose of A (not needed for orthogonal)
  at <- t(a)

  # check init.ROT
  if (is.null(init_rot)) {
    rot <- diag(m)
  } else {
    rot <- init_rot
  }

  # set initial value of alpha to 1
  alpha <- 1

  # initial rotation
  if (orthogonal) {
    mm_lambda <- a %*% rot
  } else {
    mm_lambda <- t(solve(rot, at))
  }

  # using the current mm_lambda, evaluate the user-specified
  # rotation criterion; return Q (the criterion) and its gradient Gq
  q_1 <- do.call(
    method_fname,
    c(list(mm_lambda = mm_lambda), method_args, list(grad = TRUE))
  )
  gq <- attr(q_1, "grad")
  attr(q_1, "grad") <- NULL
  q_current <- q_1

  # compute gradient GRAD of f() at ROT from the gradient Gq of Q at mm_lambda
  # in a manner appropriate for orthogonal or oblique rotation
  if (orthogonal) {
    grad <- crossprod(a, gq)
  } else {
    grad <- -1 * solve(t(init_rot), crossprod(gq, mm_lambda))
  }

  # start iterations
  converged <- FALSE
  for (iter in seq_len(max_iter + 1L)) {
    # compute projection Gp of GRAD onto the linear manifold tangent at
    # ROT to the manifold of orthogonal or normal (for oblique) matrices
    #
    # this projection is zero if and only if ROT is a stationary point of
    # f() restricted to the orthogonal/normal matrices
    if (orthogonal) {
      mm <- crossprod(rot, grad)
      symm <- (mm + t(mm)) / 2
      gp <- grad - (rot %*% symm)
    } else {
      gp <- grad - t(t(rot) * colSums(rot * grad))
    }

    # check Frobenius norm of Gp
    frob <- sqrt(sum(gp * gp))

    # if verbose, print
    if (lav_verbose()) {
      cat(
        "  iter = ", sprintf("%4d", iter - 1),
        " Q = ", sprintf("%9.7f", q_current),
        " frob.log10 = ", sprintf("%10.7f", log10(frob)),
        " alpha = ", sprintf("%9.7f", alpha), "\n"
      )
    }

    if (frob < gpa_tol) {
      converged <- TRUE
      break
    }

    # update
    alpha <- 2 * alpha
    for (i in seq_len(1000)) { # make option?

      # step in the negative projected gradient direction
      # (note, the original algorithm in Jennrich 2001 used G, not Gp)
      x <- rot - alpha * gp

      if (orthogonal) {
        # use SVD to compute the projection ROTt of X onto the manifold
        # of orthogonal matrices
        svd_out <- svd(x)
        u <- svd_out$u
        v_1 <- svd_out$v
        rott <- u %*% t(v_1)
      } else {
        # compute the projection ROTt of X onto the manifold
        # of normal matrices
        v <- 1 / sqrt(apply(x^2, 2, sum))
        rott <- x %*% diag(v)
      }

      # rotate again
      if (orthogonal) {
        mm_lambda <- a %*% rott
      } else {
        mm_lambda <- t(solve(rott, at))
      }

      # evaluate criterion
      q_new <- do.call(method_fname, c(
        list(mm_lambda = mm_lambda),
        method_args, list(grad = TRUE)
      ))
      gq <- attr(q_new, "grad")
      attr(q_new, "grad") <- NULL

      # check stopping criterion
      if (q_new < q_current - 0.5 * frob * frob * alpha) {
        break
      } else {
        alpha <- alpha / 2
      }

      if (i == 1000) {
        lav_msg_warn(gettext("half-stepping failed in GPA"))
      }
    }

    # update
    rot <- rott
    q_current <- q_new

    if (orthogonal) {
      grad <- crossprod(a, gq)
    } else {
      grad <- -1 * solve(t(rot), crossprod(gq, mm_lambda))
    }
  } # iter

  # warn if no convergence
  if (!converged) {
    lav_msg_warn(gettextf(
      "GP rotation algorithm did not converge after %s iterations",
      max_iter
    ))
  }

  # algorithm information
  info <- list(
    algorithm = "gpa",
    iter = iter - 1L,
    converged = converged,
    method.value = q_current
  )

  attr(rot, "info") <- info

  rot
}


# pairwise rotation algorithm with direct line search
#
# based on Kaiser's (1959) algorithm and Jennrich and Sampson (1966) algorithm
# but to make it generic, a line search is used; inspired by Browne 2001
#
# - orthogonal: rotate one pair of columns (=plane) at a time
# - oblique: rotate 1 factor in one pair of columns (=plane) at a time
#            note: in the oblique case, (1,2) is not the same as (2,1)
# - BUT use optimize() to find the optimal angle (for each plane)
#   (see Browne, 2001, page 130)
# - repeat until the changes in the f() criterion are below tol
#

lav_matrix_rotate_pairwise <- function(a = NULL, # original matrix
                                       orthogonal = FALSE,
                                       init_rot = NULL,
                                       method_fname = NULL, # crit function
                                       method_args = list(), # method args
                                       tol = 1e-8,
                                       max_iter = 1000L) {
  # number of columns
  m <- ncol(a)

  # initial mm_lambda + PHI
  if (is.null(init_rot)) {
    mm_lambda <- a
    if (!orthogonal) {
      phi <- diag(m)
    }
  } else {
    if (orthogonal) {
      mm_lambda <- a %*% init_rot
    } else {
      mm_lambda <- t(solve(init_rot, t(a)))
      phi <- crossprod(init_rot)
    }
  }

  # using the current mm_lambda, evaluate the user-specified
  # rotation criterion; return Q (the criterion) only
  q_current <- do.call(method_fname, c(
    list(mm_lambda = mm_lambda),
    method_args, list(grad = FALSE)
  ))

  # if verbose, print
  if (lav_verbose()) {
    cat(
      "  iter = ", sprintf("%4d", 0),
      " Q = ", sprintf("%13.11f", q_current), "\n"
    )
  }

  # plane combinations
  if (orthogonal) {
    plane <- utils::combn(m, 2)
  } else {
    tmp <- utils::combn(m, 2)
    plane <- cbind(tmp, tmp[c(2, 1), , drop = FALSE])
  }

  # define objective function -- orthogonal
  objf_orth <- function(theta = 0, A = NULL, col1 = 0L, col2 = 0L) {
    # construct ROT
    rot <- diag(m)
    rot[col1, col1] <- base::cos(theta)
    rot[col1, col2] <- base::sin(theta)
    rot[col2, col1] <- -1 * base::sin(theta)
    rot[col2, col2] <- base::cos(theta)

    # rotate
    mm_lambda <- a %*% rot

    # evaluate criterion
    q_1 <- do.call(method_fname, c(
      list(mm_lambda = mm_lambda),
      method_args, list(grad = FALSE)
    ))

    q_1
  }

  # define objective function -- oblique
  objf_obliq <- function(delta = 0, a = NULL, col1 = 0L, col2 = 0L,
                         phi12 = 0) {
    # construct ROT
    rot <- diag(m)

    # gamma
    gamma2 <- 1 + (2 * delta * phi12) + (delta * delta)

    rot[col1, col1] <- sqrt(abs(gamma2))
    rot[col1, col2] <- -1 * delta
    rot[col2, col1] <- 0
    rot[col2, col2] <- 1

    # rotate
    mm_lambda <- a %*% rot

    # evaluate criterion
    q_1 <- do.call(method_fname, c(
      list(mm_lambda = mm_lambda),
      method_args, list(grad = FALSE)
    ))
    q_1
  }

  # start iterations
  converged <- FALSE
  q_old <- q_current
  for (iter in seq_len(max_iter)) {
    # rotate - one cycle
    for (pl in seq_len(ncol(plane))) {
      # choose plane
      col1 <- plane[1, pl]
      col2 <- plane[2, pl]

      # optimize
      if (orthogonal) {
        out <- optimize(
          f = objf_orth, interval = c(-pi / 4, +pi / 4),
          A = mm_lambda, col1 = col1, col2 = col2,
          maximum = FALSE, tol = .Machine$double.eps^0.25
        )
        # best rotation - for this plane
        theta <- out$minimum

        # construct ROT
        rot <- diag(m)
        rot[col1, col1] <- base::cos(theta)
        rot[col1, col2] <- base::sin(theta)
        rot[col2, col1] <- -1 * base::sin(theta)
        rot[col2, col2] <- base::cos(theta)
      } else {
        phi12 <- phi[col1, col2]
        out <- optimize(
          f = objf_obliq, interval = c(-1, +1),
          A = mm_lambda, col1 = col1, col2 = col2,
          phi12 = phi12,
          maximum = FALSE, tol = .Machine$double.eps^0.25
        )

        # best rotation - for this plane
        delta <- out$minimum

        # construct ROT
        rot <- diag(m)

        # gamma
        gamma2 <- 1 + (2 * delta * phi12) + (delta * delta)
        gamma <- sqrt(abs(gamma2))

        rot[col1, col1] <- gamma
        rot[col1, col2] <- -1 * delta
        rot[col2, col1] <- 0
        rot[col2, col2] <- 1
      }

      # rotate
      mm_lambda <- mm_lambda %*% rot

      if (!orthogonal) {
        # rotate PHI
        phi[col1, ] <- (1 / gamma) * phi[col1, ] + (delta / gamma) * phi[col2, ]
        phi[, col1] <- phi[col1, ]
        phi[col1, col1] <- 1
      }
    } # all planes

    # check for convergence
    q_current <- do.call(method_fname, c(
      list(mm_lambda = mm_lambda),
      method_args, list(grad = FALSE)
    ))

    # absolute change in Q
    diff <- abs(q_old - q_current)

    # if verbose, print
    if (lav_verbose()) {
      cat(
        "  iter = ", sprintf("%4d", iter),
        " Q = ", sprintf("%13.11f", q_current),
        " change = ", sprintf("%13.11f", diff), "\n"
      )
    }

    if (diff < tol) {
      converged <- TRUE
      break
    } else {
      q_old <- q_current
    }
  } # iter

  # warn if no convergence
  if (!converged) {
    lav_msg_warn(gettextf(
      "pairwise rotation algorithm did not converge after %s iterations",
      max_iter
    ))
  }

  # compute final rotation matrix
  if (orthogonal) {
    rot <- solve(crossprod(a), crossprod(a, mm_lambda))
  } else {
    # to be compatible with GPa
    rott_inv <- solve(crossprod(a), crossprod(a, mm_lambda))
    rot <- solve(t(rott_inv))
  }

  # algorithm information
  info <- list(
    algorithm = "pairwise",
    iter = iter,
    converged = converged,
    method.value = q_current
  )

  attr(rot, "info") <- info

  rot
}
