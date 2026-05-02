# rotation algorithms for multiple groups
#
# YR 22 Dec 2025 -- initial version
#
# - rstarts only affects the 'initial' rotation

# main function to rotate a list of matrices 'Alist'
lav_matrix_rotate_mg <- function(a_list = NULL, # original matrices
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
                                 init_rot_list = NULL, # initial rotation matrix
                                 init_rot_check = TRUE, #check if init ROT is ok
                                 rstarts = 30L, # number of random starts
                                 # row_weights = "default", # row weighting
                                 # std.ov = FALSE, # rescale ov
                                 # ov.var = NULL, # ov variances
                                 # algorithm = "gpa", # rotation algorithm
                                 mg_algorithm = "pairwise",
                                 mg_agreement_method = "procrustes",
                                 mg_agreement_weight = 0.5,
                                 reflect = TRUE, # refect sign
                                 order_lv_by = "index", # how to order the lv's
                                 # gpa.tol = 0.00001, # stopping tol gpa
                                 tol = 1e-07, # stopping tol others
                                 keep_rep = FALSE, # store replications
                                 max_iter = 10000L) { # max iterations

  # check Alist
  if (!is.list(a_list)) {
    lav_msg_stop(gettext("a_list does not seem to be a list"))
  }
  ngroups <- length(a_list)

  # check A
  a <- a_list[[1]]
  if (!inherits(a, "matrix")) {
    lav_msg_stop(gettext("a does not seem to be a matrix"))
  }

  p <- nrow(a)
  m <- ncol(a)
  if (m < 2L) { # single dimension
    res <- list(
      lambdaList = a_list, rotList = NULL,
      orthogonal = orthogonal, method = "none",
      method.args = list(), row_weights = "none",
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
  if (!is.null(init_rot_list) && init_rot_check) {
    for (g in seq_len(ngroups)) {
      init_rot <- init_rot_list[[g]]
      if (!inherits(init_rot, "matrix")) {
        lav_msg_stop(gettext("init.ROT does not seem to be a matrix"))
      }
      if (nrow(init_rot) != m) {
        lav_msg_stop(gettextf(
          "nrow(init.ROT) = %1$s does not equal ncol(A) = %2$s",
          nrow(init_rot), m
        ))
      }
      if (nrow(init_rot) != ncol(init_rot)) {
        lav_msg_stop(gettextf(
          "nrow(init.ROT) = %1$s does not equal ncol(init.ROT) = %2$s",
          nrow(init_rot), ncol(init_rot)
        ))
      }
      # rotation matrix?
      if (!lav_matrix_rotate_check(init_rot, orthogonal = orthogonal)) {
        lav_msg_stop(gettext("init.ROT does not look like a rotation matrix"))
      }
    } # group
  } # check init.rotList

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
      method_args$target <- target <- target[[1L]] # always take the first group
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
      method_args$target_mask <- target_mask <- target_mask[[1L]]
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

  # 0. initialize: rotate each group separately + reorder
  alist_orig <- a_list
  order_list <- vector("list", length = ngroups)
  for (g in seq_len(ngroups)) {
    if (lav_verbose()) {
      cat("  group = ", g, "\n")
    }
    res0 <- lav_matrix_rotate(
      a = a_list[[g]],
      orthogonal = orthogonal,
      method = method,
      method_args = method_args,
      init_rot = NULL,
      init_rot_check = FALSE,
      rstarts = rstarts,
      row_weights = "none",
      std_ov = FALSE,
      ov_var = NULL,
      algorithm = "pairwise",
      reflect = reflect,
      order_lv_by = order_lv_by,
      tol = tol,
      max_iter = max_iter,
      group = g
    )
    if (g > 1L) {
      # reorder factors to align with the first group
      order_idx <- lav_efa_find_best_order(
        lambda_ref = a_list[[1]],
        lambda_target = res0$LAMBDA
      )
      a_list[[g]] <- res0$LAMBDA[, order_idx, drop = FALSE]
      order_list[[g]] <- order_idx
    } else {
      a_list[[g]] <- res0$LAMBDA
      order_list[[g]] <- res0$order.idx
    }
  } # group

  #  # set row_weights
  #  row_weights <- tolower(row_weights)
  #  if (row_weights == "default") {
  #    # the default is "none", except for varimax
  #    if (method %in% c("varimax", "promax")) {
  #      row_weights <- "kaiser"
  #    } else {
  #      row_weights <- "none"
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
  #   if (row_weights == "none") {
  #     weights <- rep(1.0, P)
  #   } else if (row_weights == "kaiser") {
  #     weights <- lav_matrix_rotate_kaiser_weights(A)
  #   } else if (row_weights == "cureton-mulaik") {
  #     weights <- lav_matrix_rotate_cm_weights(A)
  #   } else {
  #     lav_msg_stop(gettext("row_weights can be none,
  #                    kaiser or cureton-mulaik"))
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
    rot_list <- lav_matrix_rotate_pairwise_mg(
      a_list = a_list,
      orthogonal = orthogonal,
      random_rot = FALSE,
      method_fname = method_fname,
      method_args = method_args,
      mg_agreement_method = mg_agreement_method,
      mg_agreement_weight = mg_agreement_weight,
      tol = tol,
      max_iter = max_iter
    )
    info <- attr(rot_list, "info")
    attr(rot_list, "info") <- NULL
  #}
  lambda_list <- vector("list", length = ngroups)
  for (g in seq_len(ngroups)) {
    a <- a_list[[g]]
    rot <- rot_list[[g]]

    if (orthogonal) {
      # LAMBDA <- A %*% solve(t(ROT))
      # note: when ROT is orthogonal, solve(t(ROT)) == ROT
      mm_lambda <- a %*% rot
    } else {
      # LAMBDA <- A %*% solve(t(ROT))
      mm_lambda <- t(solve(rot, t(a)))
    }

    #  # 3. undo row weighting
    #  LAMBDA <- LAMBDA / weights

    #  # 3.b undo cov -> cor
    #  if (std.ov) {
    #    LAMBDA <- LAMBDA * sqrt(ov.var)
    #  }

    # put back in List
    lambda_list[[g]] <- mm_lambda
  } # group

  # re-compute final rotation matrix, using Alist.orig
  # because in step 0, we already rotated (per group)
  rot_list <- vector("list", length = ngroups)
  if (orthogonal) {
    for (g in seq_len(ngroups)) {
      rot_list[[g]] <- solve(
        crossprod(alist_orig[[g]]),
        crossprod(alist_orig[[g]], lambda_list[[g]])
      )
    }
  } else {
    for (g in seq_len(ngroups)) {
      # to be compatible with GPa
      rott_inv <- solve(
        crossprod(alist_orig[[g]]),
        crossprod(alist_orig[[g]], lambda_list[[g]])
      )
      rot_list[[g]] <- solve(t(rott_inv))
    }
  }

  # 6. return results as a list
  res <- list(
    lambdaList = lambda_list, rotList = rot_list,
    orderList = order_list,
    orthogonal = orthogonal, method = method,
    method_args = method_args, row_weights = "none"
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
lav_matrix_rotate_pairwise_mg <- function(a_list = NULL, # original matrices
                                          orthogonal = FALSE,
                                          random_rot = FALSE,
                                          method_fname = NULL, # crit function
                                          method_args = list(), # method args
                                          tol = 1e-8,
                                          mg_agreement_weight = 0.50,
                                          mg_agreement_method = "procrustes",
                                          scale = TRUE,
                                          max_iter = 1000L) {
  # first group is 'A'
  a <- a_list[[1]]
  ngroups <- length(a_list)

  # number of columns
  m <- ncol(a)

  # initial LAMBDA + PSI
  lambda_list <- a_list
  psi_list <- rep(list(diag(m)), ngroups)

  # random initial rotation?
  if (random_rot) {
    for (g in seq_len(ngroups)) {
      init_rot <- lav_matrix_rotate_gen(m = m, orthogonal = TRUE)
      if (orthogonal) {
        lambda_list[[g]] <- a_list[[g]] %*% init_rot
      } else {
        lambda_list[[g]] <- t(solve(init_rot, t(a_list[[g]])))
        psi_list[[g]] <- crossprod(init_rot)
      }
    }
  }

  # using the current LAMBDA, evaluate the user-specified
  # rotation criteron; return Q (the criterion) only
  q_current <- lav_matrix_rotate_mg_agreement(
    lambda_list = lambda_list, method_fname = method_fname,
    method_args = method_args,
    mg_agreement_method = mg_agreement_method,
    w = mg_agreement_weight, scale = scale
  )

  # if verbose, print
  if (lav_verbose()) {
    q_mg <- attr(q_current, "Q.mg")
    a_mg <- attr(q_current, "A.mg")
    cat(
      "  iter = ", sprintf("%4d", 0),
      " Q = ", sprintf("%13.11f", q_current),
      " Q.rot = ", sprintf("%13.11f", q_mg),
      " Q.agr = ", sprintf("%13.11f", a_mg), "\n"
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
  objf_orth <- function(g = 1L, theta = 0, a = NULL, col1 = 0L, col2 = 0L) {
    # construct ROT
    rot <- diag(m)
    rot[col1, col1] <- base::cos(theta)
    rot[col1, col2] <- base::sin(theta)
    rot[col2, col1] <- -1 * base::sin(theta)
    rot[col2, col2] <- base::cos(theta)

    # rotate
    lambda_list[[g]] <- a %*% rot

    # evaluate criterion
    q_1 <- lav_matrix_rotate_mg_agreement(
      lambda_list = lambda_list, method_fname = method_fname,
      method_args = method_args, mg_agreement_method = mg_agreement_method,
      w = mg_agreement_weight, scale = scale
    )

    q_1
  }

  # define objective function -- oblique
  objf_obliq <- function(g = 1, delta = 0, A = NULL, col1 = 0L, col2 = 0L,
                         psi12 = 0) {
    # construct ROT
    rot <- diag(m)

    # gamma
    gamma2 <- 1 + (2 * delta * psi12) + (delta * delta)

    rot[col1, col1] <- sqrt(abs(gamma2))
    rot[col1, col2] <- -1 * delta
    rot[col2, col1] <- 0
    rot[col2, col2] <- 1

    # rotate
    lambda_list[[g]] <- a %*% rot

    # evaluate criterion
    q_1 <- lav_matrix_rotate_mg_agreement(
      lambda_list = lambda_list, method_fname = method_fname,
      method_args = method_args, mg_agreement_method = mg_agreement_method,
      w = mg_agreement_weight, scale = scale
    )

    q_1
  }

  # start iterations
  converged <- FALSE
  q_old <- q_current
  for (iter in seq_len(max_iter)) {
    # for each group
    for (g in seq_len(ngroups)) {
      # rotate - one cycle
      for (pl in seq_len(ncol(plane))) {
        # choose plane
        col1 <- plane[1, pl]
        col2 <- plane[2, pl]

        # optimize
        if (orthogonal) {
          out <- optimize(
            f = objf_orth, interval = c(-pi / 4, +pi / 4),
            g = g,
            a = lambda_list[[g]], col1 = col1, col2 = col2,
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
          psi12 <- psi_list[[g]][col1, col2]
          out <- optimize(
            f = objf_obliq, interval = c(-1, +1),
            g = g,
            A = lambda_list[[g]], col1 = col1, col2 = col2,
            psi12 = psi12,
            maximum = FALSE, tol = .Machine$double.eps^0.25
          )

          # best rotation - for this plane
          delta <- out$minimum

          # construct ROT
          rot <- diag(m)

          # gamma
          gamma2 <- 1 + (2 * delta * psi12) + (delta * delta)
          gamma <- sqrt(abs(gamma2))

          rot[col1, col1] <- gamma
          rot[col1, col2] <- -1 * delta
          rot[col2, col1] <- 0
          rot[col2, col2] <- 1
        }

        # rotate
        lambda_list[[g]] <- lambda_list[[g]] %*% rot

        if (!orthogonal) {
          mm_psi <- psi_list[[g]]
          # rotate PSI
          mm_psi[col1, ] <- (1 / gamma) * mm_psi[col1, ] +
            (delta / gamma) * mm_psi[col2, ]
          mm_psi[, col1] <- mm_psi[col1, ]
          mm_psi[col1, col1] <- 1
          psi_list[[g]] <- mm_psi
        }
      } # all planes
    } # ngroups

    # check for convergence
    q_current <- lav_matrix_rotate_mg_agreement(
      lambda_list = lambda_list, method_fname = method_fname,
      method_args = method_args, mg_agreement_method = mg_agreement_method,
      w = mg_agreement_weight, scale = scale
    )

    # absolute change in Q
    diff <- abs(q_old - q_current)

    # if verbose, print
    if (lav_verbose()) {
      q_mg <- attr(q_current, "Q.mg")
      a_mg <- attr(q_current, "A.mg")
      cat(
        "  iter = ", sprintf("%4d", iter),
        " Q = ", sprintf("%13.11f", q_current),
        " Q.rot = ", sprintf("%13.11f", q_mg),
        " Q.agr = ", sprintf("%13.11f", a_mg),
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
      q_old <- q_current
    }
  } # iter

  # warn if no convergence
  if (!converged) {
    lav_msg_warn(gettextf(
      "mg-pairwise rotation algorithm did not converge after %s iterations",
      max_iter
    ))
  }

  # compute final rotation matrix
  rot_list <- vector("list", length = ngroups)
  if (orthogonal) {
    for (g in seq_len(ngroups)) {
      rot_list[[g]] <- solve(
        crossprod(a_list[[g]]),
        crossprod(a_list[[g]], lambda_list[[g]])
      )
    }
  } else {
    for (g in seq_len(ngroups)) {
      # to be compatible with GPa
      rott_inv <- solve(
        crossprod(a_list[[g]]),
        crossprod(a_list[[g]], lambda_list[[g]])
      )
      rot_list[[g]] <- solve(t(rott_inv))
    }
  }

  # algorithm information
  info <- list(
    algorithm = "mg-pairwise",
    iter = iter,
    converged = converged,
    method.value = q_current
  )

  attr(rot_list, "info") <- info

  rot_list
}
