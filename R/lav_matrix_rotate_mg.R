# rotation algorithms for multiple groups
#
# YR 22 Dec 2025 -- initial version
#
# - rstarts only affects the 'initial' rotation

# main function to rotate a list of matrices 'Alist'
lav_mat_rotate_mg <- function(a_list = NULL, # original matrices
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
                                 reflect = TRUE, # reflect sign
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
      if (!lav_mat_rotate_check(init_rot, orthogonal = orthogonal)) {
        lav_msg_stop(gettext("init.ROT does not look like a rotation matrix"))
      }
    } # group
  } # check init.rotList

  # determine method function name + validate target/target_mask
  # (multigroup always uses the first group's target/target_mask)
  tmp_resolve <- lav_mat_rotate_resolve_method(
    method = method, method_args = method_args, p = p, m = m, group = 1L
  )
  method <- tmp_resolve$method
  method_fname <- tmp_resolve$method_fname
  method_args <- tmp_resolve$method_args

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
    res0 <- lav_mat_rotate(
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
        lambda_target = res0$mm_lambda
      )
      a_list[[g]] <- res0$mm_lambda[, order_idx, drop = FALSE]
      order_list[[g]] <- order_idx
    } else {
      a_list[[g]] <- res0$mm_lambda
      order_list[[g]] <- res0$order.idx
    }
  } # group

  # rotate: mg-pairwise + agreement
  # (note: rstarts only affects the per-group initial rotation above)
  if (lav_verbose()) {
    cat("  mg-pairwise + agreement rotation:\n")
  }
  rot_list <- lav_mat_rotate_pairwise_mg(
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
# (see lav_mat_rotate_pairwise()
#
# but adapted to handle multiple groups + an agreement criterion
#
lav_mat_rotate_pairwise_mg <- function(a_list = NULL, # original matrices
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
      init_rot <- lav_mat_rotate_gen(m = m, orthogonal = TRUE)
      if (orthogonal) {
        lambda_list[[g]] <- a_list[[g]] %*% init_rot
      } else {
        lambda_list[[g]] <- t(solve(init_rot, t(a_list[[g]])))
        psi_list[[g]] <- crossprod(init_rot)
      }
    }
  }

  # using the current LAMBDA, evaluate the user-specified
  # rotation criterion; return Q (the criterion) only
  q_current <- lav_mat_rotate_mg_agreement(
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
    rot <- lav_mat_givens_orth(m, col1, col2, theta)

    # rotate
    lambda_list[[g]] <- a %*% rot

    # evaluate criterion
    q_1 <- lav_mat_rotate_mg_agreement(
      lambda_list = lambda_list, method_fname = method_fname,
      method_args = method_args, mg_agreement_method = mg_agreement_method,
      w = mg_agreement_weight, scale = scale
    )

    q_1
  }

  # define objective function -- oblique
  # ('a' is the current loading matrix for group g, passed in by optimize())
  objf_obliq <- function(g = 1, delta = 0, a = NULL, col1 = 0L, col2 = 0L,
                         psi12 = 0) {
    # construct ROT
    rot <- lav_mat_givens_obliq(m, col1, col2, delta, psi12)

    # rotate
    lambda_list[[g]] <- a %*% rot

    # evaluate criterion
    q_1 <- lav_mat_rotate_mg_agreement(
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
          rot <- lav_mat_givens_orth(m, col1, col2, theta)
        } else {
          psi12 <- psi_list[[g]][col1, col2]
          out <- optimize(
            f = objf_obliq, interval = c(-1, +1),
            g = g,
            a = lambda_list[[g]], col1 = col1, col2 = col2,
            psi12 = psi12,
            maximum = FALSE, tol = .Machine$double.eps^0.25
          )

          # best rotation - for this plane
          delta <- out$minimum

          # construct ROT (rot[col1, col1] == gamma, the PSI scaling factor)
          rot <- lav_mat_givens_obliq(m, col1, col2, delta, psi12)
          gamma <- rot[col1, col1]
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
    q_current <- lav_mat_rotate_mg_agreement(
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
