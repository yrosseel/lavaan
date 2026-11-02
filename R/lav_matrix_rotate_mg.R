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
                                   simplimax_k = 0L,
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
                                 mg_algorithm = "gpa", # agreement algorithm
                                 mg_agreement_method = "procrustes",
                                 mg_agreement_epsilon = 0.01, # alignment
                                 mg_agreement_weight = 0.5,
                                 mg_agreement_group_weights = NULL,
                                 mg_agreement_scales = "equal", # or "free"
                                 mg_agreement_rotate = TRUE, # joint rotation
                                 mg_agreement_intercepts = "none", # or
                                                   # "fixed" / "free"
                                 nu_list = NULL, # intercepts (alignment)
                                 reflect = TRUE, # reflect sign
                                 order_lv_by = "index", # how to order the lv's
                                 pairwise_accelerate = TRUE, # sweep extrap.
                                 gpa_tol = 0.00001, # stopping tol (gradient)
                                 tol = 1e-07, # stopping tol others
                                 keep_rep = FALSE, # store replications
                                 max_iter = 10000L) { # max iterations

  # scale split across groups: "equal" (equal column sums of squares, the
  # scales are determined by the rotation) or "free" (the scales of the
  # factors are free parameters of the agreement objective, subject to
  # mean variance one across groups, as in multigroup factor rotation)
  mg_agreement_scales <- match.arg(mg_agreement_scales, c("equal", "free"))
  free_scales <- mg_agreement_scales == "free"
  # intercept alignment: the factor means are extra variables of the
  # agreement objective ("fixed": zero in the first group; "free": free
  # in every group), and the intercepts of the indicators (nu_list) are
  # aligned along with the loadings
  mg_agreement_intercepts <- match.arg(mg_agreement_intercepts,
                                       c("none", "fixed", "free"))
  if (mg_agreement_intercepts != "none" && is.null(nu_list)) {
    lav_msg_stop(gettext(
      "the alignment of the intercepts (mg_agreement_intercepts) requires
       the intercepts of the indicators (nu_list)."))
  }
  if (mg_agreement_intercepts == "none") {
    nu_list <- NULL
  }

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
  # bifactor rotation: the first column is the general factor and stays in
  # place when the groups are aligned to the first group
  bifactor_flag <- method %in% c(
    "bi-geomin", "bigeomin", "bi-quartimin", "biquartimin"
  )
  # a free scale split needs a criterion that is bounded below when the
  # loadings of a factor are scaled up (the agreement term grows with the
  # square of the loadings, so a criterion that decreases with the fourth
  # power, like orthomax/varimax, sends the scale of a factor in one group
  # to zero)
  if (free_scales) {
    scale_ok <- method %in% c(
      "quartimin", "geomin", "target", "target.strict", "pst", "simplimax",
      "bi-quartimin", "biquartimin", "bi-geomin", "bigeomin"
    ) || startsWith(method, "cf") ||
      (method == "oblimin" && isTRUE(method_args$oblimin_gamma == 0))
    if (!scale_ok) {
      lav_msg_stop(gettextf(
        "mg_agreement_scales = \"free\" requires a rotation criterion that is
         bounded below when the loadings of a factor are scaled up (oblimin
         with gamma = 0, quartimin, geomin, the Crawford-Ferguson family,
         target, pst, simplimax, bi-quartimin, bi-geomin); not %s.",
        dQuote(method)))
    }
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
  # factor correlations of the per-group solutions: the starting point of
  # the agreement sweep below (oblique rotation; identity if orthogonal)
  phi_list <- vector("list", length = ngroups)
  # per-group information about the initial rotation (for lavInspect())
  group_info <- vector("list", length = ngroups)
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
      pairwise_accelerate = pairwise_accelerate,
      reflect = reflect,
      order_lv_by = order_lv_by,
      gpa_tol = gpa_tol,
      tol = tol,
      keep_rep = keep_rep,
      max_iter = max_iter,
      group = g
    )
    group_info[[g]] <- list(
      algorithm = res0$algorithm, converged = res0$converged,
      iter = res0$iter, method.value = res0$method.value,
      REP = res0$REP
    )
    if (g > 1L) {
      # reorder factors to align with the first group
      order_idx <- lav_efa_find_best_order(
        lambda_ref = a_list[[1]],
        lambda_target = res0$mm_lambda,
        general_first = bifactor_flag
      )
      a_list[[g]] <- res0$mm_lambda[, order_idx, drop = FALSE]
      order_list[[g]] <- order_idx
      phi_list[[g]] <- res0$PHI[order_idx, order_idx, drop = FALSE]
    } else {
      a_list[[g]] <- res0$mm_lambda
      order_list[[g]] <- res0$order.idx
      phi_list[[g]] <- res0$PHI
    }
  } # group

  # rotate: all groups together + agreement
  # (note: rstarts only affects the per-group initial rotation above)
  # the gradient projection algorithm needs a smooth criterion; for
  # simplimax we use the (derivative-free) plane sweeps instead
  if (mg_algorithm == "gpa" && !lav_mat_rotate_is_smooth(method_fname)) {
    mg_algorithm <- "pairwise"
  }
  # mg_agreement_rotate = FALSE: the per-group (aligned) patterns are kept
  # and only the scale split is optimized for agreement (the alignment of
  # the factor variances); with equal scales nothing moves here
  if (!mg_agreement_rotate) {
    mg_algorithm <- "scales"
  }
  if (lav_verbose()) {
    cat("  mg-", mg_algorithm, " + agreement rotation:\n", sep = "")
  }
  if (mg_algorithm == "scales") {
    rot_list <- lav_mat_rotate_scales_mg(
      a_list = a_list,
      phi_list = phi_list,
      orthogonal = orthogonal,
      method_fname = method_fname,
      method_args = method_args,
      mg_agreement_method = mg_agreement_method,
      mg_agreement_epsilon = mg_agreement_epsilon,
      mg_agreement_weight = mg_agreement_weight,
      group_weights = mg_agreement_group_weights,
      scales = mg_agreement_scales,
      intercepts = mg_agreement_intercepts,
      nu_list = nu_list,
      rstarts = rstarts
    )
  } else if (mg_algorithm == "gpa") {
    rot_list <- lav_mat_rotate_gpa_mg(
      a_list = a_list,
      phi_list = phi_list,
      orthogonal = orthogonal,
      method_fname = method_fname,
      method_args = method_args,
      mg_agreement_method = mg_agreement_method,
      mg_agreement_epsilon = mg_agreement_epsilon,
      mg_agreement_weight = mg_agreement_weight,
      group_weights = mg_agreement_group_weights,
      scales = mg_agreement_scales,
      intercepts = mg_agreement_intercepts,
      nu_list = nu_list,
      gpa_tol = gpa_tol,
      max_iter = max_iter
    )
  } else {
    rot_list <- lav_mat_rotate_pairwise_mg(
      a_list = a_list,
      phi_list = phi_list,
      orthogonal = orthogonal,
      random_rot = FALSE,
      method_fname = method_fname,
      method_args = method_args,
      mg_agreement_method = mg_agreement_method,
      mg_agreement_epsilon = mg_agreement_epsilon,
      mg_agreement_weight = mg_agreement_weight,
      group_weights = mg_agreement_group_weights,
      scales = mg_agreement_scales,
      intercepts = mg_agreement_intercepts,
      nu_list = nu_list,
      accelerate = pairwise_accelerate,
      gpa_tol = gpa_tol,
      tol = tol,
      max_iter = max_iter
    )
  }
  info <- attr(rot_list, "info")
  attr(rot_list, "info") <- NULL
  # the aligned factor means (per group), if any
  mean_list <- NULL
  if (!is.null(info$means)) {
    mean_list <- lapply(seq_len(ngroups), function(g) info$means[, g])
    info$means <- NULL
  }
  lambda_list <- vector("list", length = ngroups)
  for (g in seq_len(ngroups)) {
    a <- a_list[[g]]
    rot <- rot_list[[g]]

    if (orthogonal && !free_scales) {
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
  # (with free scales the rotation matrices include the diagonal scaling,
  # so the general form lambda = a solve(t(rot)) is used for orthogonal
  # rotation too)
  rot_list <- vector("list", length = ngroups)
  if (orthogonal && !free_scales) {
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
    orderList = order_list, group_info = group_info,
    orthogonal = orthogonal, method = method,
    method_args = method_args, row_weights = "none",
    mg_agreement_scales = mg_agreement_scales,
    mg_agreement_intercepts = mg_agreement_intercepts,
    meanList = mean_list
  )

  # add method info
  res <- c(res, info)

  res
}


# pairwise (plane) rotation algorithm for multiple groups + an agreement
# criterion: the plane sweeps (over all groups) are run by the engine
# shared with the single-group algorithm (lav_mat_rotate_pairwise_core(),
# in lav_matrix_rotate.R), driven by the gradient of the agreement
# objective (lav_mat_rotate_mg_agreement(grad = TRUE))
#
# phi_list: the factor correlation matrices that go with a_list (oblique
# rotation). The starting matrices are the per-group rotated solutions, so
# their factors are correlated; starting the oblique plane sweeps from an
# identity matrix instead (as before 0.7-3) mis-normalizes the factors
# (the rotated solution then had factor variances different from one)
#
lav_mat_rotate_pairwise_mg <- function(a_list = NULL, # original matrices
                                       phi_list = NULL, # factor correlations
                                       orthogonal = FALSE,
                                       random_rot = FALSE,
                                       method_fname = NULL, # crit function
                                       method_args = list(), # method args
                                       accelerate = TRUE,
                                       gpa_tol = 1e-5,
                                       tol = 1e-8,
                                       mg_agreement_weight = 0.50,
                                       mg_agreement_method = "procrustes",
                                       mg_agreement_epsilon = 0.01,
                                       group_weights = NULL,
                                       scale = TRUE,
                                       scales = "equal",
                                       intercepts = "none",
                                       nu_list = NULL,
                                       max_iter = 1000L) {
  ngroups <- length(a_list)

  # number of columns
  m <- ncol(a_list[[1L]])

  # free scale split: the plane sweeps (given the scales) alternate with
  # a minimization over the scales (given the rotations); the objective
  # is then evaluated on lambda_g diag(1/c_g) without the equal-SS
  # rescaling. The same alternation handles the factor means of the
  # intercept alignment
  free_scales <- identical(scales, "free")
  if (free_scales) {
    scale <- FALSE
  }
  align_means <- intercepts != "none" && !is.null(nu_list)
  a_mat <- NULL
  if (align_means) {
    a_mat <- matrix(0, m, ngroups)
  }

  # initial LAMBDA + PSI
  lambda_list <- a_list
  if (is.null(phi_list) || orthogonal) {
    psi_list <- rep(list(diag(m)), ngroups)
  } else {
    psi_list <- phi_list
  }

  # random initial rotation?
  if (random_rot) {
    for (g in seq_len(ngroups)) {
      init_rot <- lav_mat_rotate_gen(m = m, orthogonal = TRUE)
      if (orthogonal) {
        lambda_list[[g]] <- a_list[[g]] %*% init_rot
      } else {
        # lambda = a %*% solve(t(init_rot)): the factors transform as
        # init_rot^-1 f, so PSI becomes init_rot^-1 PSI init_rot^-T
        # (= t(init_rot) PSI init_rot for an orthogonal init_rot)
        lambda_list[[g]] <- t(solve(init_rot, t(a_list[[g]])))
        psi_list[[g]] <- t(init_rot) %*% psi_list[[g]] %*% init_rot
      }
    }
  }

  # objective: the agreement criterion (+ per-group gradients); with free
  # scales c_mat (m x ngroups, read from this environment) the loadings
  # are scaled first and the gradient is chained back
  c_mat <- NULL
  objf <- function(lambda_list, grad = FALSE) {
    if (free_scales) {
      lambda_list <- lapply(seq_len(ngroups), function(g) {
        t(t(lambda_list[[g]]) / c_mat[, g])
      })
    }
    mean_list <- NULL
    if (align_means) {
      mean_list <- lapply(seq_len(ngroups), function(g) a_mat[, g])
    }
    out <- lav_mat_rotate_mg_agreement(
      lambda_list = lambda_list, method_fname = method_fname,
      method_args = method_args, mg_agreement_method = mg_agreement_method,
      mg_agreement_epsilon = mg_agreement_epsilon,
      w = mg_agreement_weight, scale = scale,
      group_weights = group_weights,
      nu_list = nu_list, mean_list = mean_list, grad = grad
    )
    if (free_scales && grad && !is.null(attr(out, "grad"))) {
      gr <- attr(out, "grad")
      attr(out, "grad") <- lapply(seq_len(ngroups), function(g) {
        t(t(gr[[g]]) / c_mat[, g])
      })
    }
    out
  }

  if (!free_scales && !align_means) {
    res <- lav_mat_rotate_pairwise_core(
      lambda_list = lambda_list, phi_list = psi_list,
      orthogonal = orthogonal, objf = objf,
      use_grad = lav_mat_rotate_is_smooth(method_fname),
      accelerate = accelerate,
      gpa_tol = gpa_tol, tol = tol, max_iter = max_iter
    )
    lambda_list <- res$lambda_list
    iter <- res$iter
    converged <- res$converged
    q_current <- res$q
  } else {
    c_mat <- lav_mat_rotate_mg_scales_start(lambda_list)
    iter <- 0L
    converged <- FALSE
    q_current <- as.numeric(objf(lambda_list, grad = FALSE))
    for (outer in seq_len(max_iter)) {
      q_old <- q_current
      res <- lav_mat_rotate_pairwise_core(
        lambda_list = lambda_list, phi_list = psi_list,
        orthogonal = orthogonal, objf = objf,
        use_grad = lav_mat_rotate_is_smooth(method_fname),
        accelerate = accelerate,
        gpa_tol = gpa_tol, tol = tol, max_iter = max_iter
      )
      lambda_list <- res$lambda_list
      psi_list <- res$phi_list
      iter <- iter + res$iter
      c_mat <- lav_mat_rotate_mg_scales_update(
        lambda_list = lambda_list, c_mat = c_mat,
        method_fname = method_fname, method_args = method_args,
        mg_agreement_method = mg_agreement_method,
        mg_agreement_epsilon = mg_agreement_epsilon,
        w = mg_agreement_weight, group_weights = group_weights,
        free_scales = free_scales, nu_list = nu_list, a_mat = a_mat,
        mean_fixed_first = (intercepts == "fixed")
      )
      a_mat <- attr(c_mat, "a_mat")
      attr(c_mat, "a_mat") <- NULL
      q_current <- as.numeric(objf(lambda_list, grad = FALSE))
      if (lav_verbose()) {
        cat(
          "  outer = ", sprintf("%3d", outer),
          " Q = ", sprintf("%13.11f", q_current),
          " sweeps = ", sprintf("%4d", iter), "\n"
        )
      }
      if (abs(q_old - q_current) < tol) {
        converged <- res$converged
        break
      }
    }
    # the rotated (and scaled) loadings
    if (free_scales) {
      lambda_list <- lapply(seq_len(ngroups), function(g) {
        t(t(lambda_list[[g]]) / c_mat[, g])
      })
    }
  }

  # warn if no convergence
  if (!converged) {
    lav_msg_warn(gettextf(
      "mg-pairwise rotation algorithm did not converge after %s iterations",
      iter
    ))
  }

  # compute final rotation matrix
  rot_list <- vector("list", length = ngroups)
  if (orthogonal && !free_scales) {
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
  if (free_scales) {
    info$scales <- c_mat
  }
  if (align_means) {
    info$means <- a_mat
  }

  attr(rot_list, "info") <- info

  rot_list
}


# starting values for the free scale split: the equal-column-sums-of-squares
# rule, c_gj = sqrt(s_gj / mean_h s_hj), which has mean variance one; the
# loadings lambda_g diag(1/c_g) then have equal column sums of squares
# across groups (lav_mat_rotate_mg_agreement(scale = TRUE))
lav_mat_rotate_mg_scales_start <- function(lambda_list) {
  ss <- vapply(lambda_list, function(x) colSums(x * x),
               numeric(ncol(lambda_list[[1L]])))
  ss <- matrix(ss, nrow = ncol(lambda_list[[1L]]))
  sqrt(ss / rowMeans(ss))
}

# minimize the agreement objective over the scale split, given the rotated
# (unscaled) loadings: c_gj^2 = ngroups * exp(z_gj) / sum_h exp(z_hj)
# (positive, mean variance one by construction), with z_Gj = 0 so that the
# parametrization has no redundant direction; analytic gradient
lav_mat_rotate_mg_scales_update <- function(lambda_list, c_mat,
                                            method_fname, method_args,
                                            mg_agreement_method,
                                            mg_agreement_epsilon = 0.01,
                                            w, group_weights,
                                            free_scales = TRUE,
                                            nu_list = NULL, a_mat = NULL,
                                            mean_fixed_first = TRUE) {
  ngroups <- length(lambda_list)
  m <- ncol(lambda_list[[1L]])
  align_means <- !is.null(nu_list) && !is.null(a_mat)
  if (align_means) {
    a_free <- if (mean_fixed_first) seq_len(ngroups)[-1L] else seq_len(ngroups)
  } else {
    a_free <- integer(0L)
  }
  n_z <- if (free_scales) m * (ngroups - 1L) else 0L
  n_a <- m * length(a_free)

  unpack <- function(x) {
    c_new <- c_mat
    if (free_scales) {
      z <- cbind(matrix(x[seq_len(n_z)], m, ngroups - 1L), 0)
      p <- exp(z - apply(z, 1, max))
      p <- p / rowSums(p)
      c_new <- sqrt(ngroups * p)
    }
    a_new <- a_mat
    if (align_means) {
      a_new[, a_free] <- matrix(x[n_z + seq_len(n_a)], m, length(a_free))
    }
    list(c = c_new, a = a_new)
  }
  evalx <- function(x, grad = FALSE) {
    u <- unpack(x)
    ll <- lapply(seq_len(ngroups), function(g) {
      t(t(lambda_list[[g]]) / u$c[, g])
    })
    mean_list <- NULL
    if (align_means) {
      mean_list <- lapply(seq_len(ngroups), function(g) u$a[, g])
    }
    out <- lav_mat_rotate_mg_agreement(
      lambda_list = ll, method_fname = method_fname,
      method_args = method_args, mg_agreement_method = mg_agreement_method,
      mg_agreement_epsilon = mg_agreement_epsilon,
      w = w, scale = FALSE, group_weights = group_weights,
      nu_list = nu_list, mean_list = mean_list, grad = grad
    )
    attr(out, "ll") <- ll
    attr(out, "c") <- u$c
    out
  }
  objx <- function(x) {
    as.numeric(evalx(x, grad = FALSE))
  }
  gradx <- function(x) {
    out <- evalx(x, grad = TRUE)
    gx <- numeric(0L)
    if (free_scales) {
      ll <- attr(out, "ll")
      c_new <- attr(out, "c")
      p <- c_new * c_new / ngroups
      gr <- attr(out, "grad")
      # ctilde_gj = (lambda_g' grad_g)_jj; dF/dz_hj = -ctilde_hj/2 +
      # p_hj sum_g ctilde_gj / 2
      ct <- vapply(seq_len(ngroups), function(g) colSums(ll[[g]] * gr[[g]]),
                   numeric(m))
      ct <- matrix(ct, nrow = m)
      gz <- -ct / 2 + p * rowSums(ct) / 2
      gx <- c(gx, as.vector(gz[, -ngroups, drop = FALSE]))
    }
    if (align_means) {
      gm <- matrix(unlist(attr(out, "grad_mean")), m, ngroups)
      gx <- c(gx, as.vector(gm[, a_free, drop = FALSE]))
    }
    gx
  }

  x0 <- numeric(0L)
  if (free_scales) {
    z0 <- log(c_mat * c_mat / ngroups)
    x0 <- c(x0, as.vector(z0[, -ngroups, drop = FALSE] - z0[, ngroups]))
  }
  if (align_means) {
    x0 <- c(x0, as.vector(a_mat[, a_free, drop = FALSE]))
  }
  if (length(x0) == 0L) {
    attr(c_mat, "a_mat") <- a_mat
    attr(c_mat, "converged") <- TRUE
    attr(c_mat, "value") <- objx(x0)
    return(c_mat)
  }
  res <- try(nlminb(x0, objx, gradient = gradx,
                    control = list(iter.max = 1000L, rel.tol = 1e-10)),
             silent = TRUE)
  if (inherits(res, "try-error") || !is.finite(res$objective) ||
      res$objective > objx(x0)) {
    attr(c_mat, "a_mat") <- a_mat
    attr(c_mat, "converged") <- FALSE
    attr(c_mat, "value") <- objx(x0)
    return(c_mat)
  }
  u <- unpack(res$par)
  out <- u$c
  attr(out, "a_mat") <- u$a
  # converged: nlminb says so, or the gradient is small relative to the
  # objective (nlminb reports "singular convergence" at the optimum when
  # the relative tolerance is below what the loss allows)
  attr(out, "converged") <- (res$convergence == 0L ||
    sqrt(sum(gradx(res$par)^2)) < 1e-5 * max(1, abs(res$objective)))
  attr(out, "value") <- res$objective
  out
}


# scales-only agreement (mg_agreement_rotate = FALSE): the per-group
# rotated and aligned patterns a_list are kept, and only the split of each
# factor's scale over the groups is optimized (if scales = "free"); this is
# the alignment of the factor variances given the patterns. Returns the
# rotation matrices relative to a_list, as the other mg algorithms
lav_mat_rotate_scales_mg <- function(a_list = NULL,
                                     phi_list = NULL,
                                     orthogonal = FALSE,
                                     method_fname = NULL,
                                     method_args = list(),
                                     mg_agreement_weight = 0.50,
                                     mg_agreement_method = "procrustes",
                                     mg_agreement_epsilon = 0.01,
                                     group_weights = NULL,
                                     scales = "equal",
                                     intercepts = "none",
                                     nu_list = NULL,
                                     rstarts = 0L) {
  ngroups <- length(a_list)
  m <- ncol(a_list[[1L]])
  free_scales <- identical(scales, "free")
  align_means <- intercepts != "none" && !is.null(nu_list)
  converged <- TRUE
  lambda_list <- a_list
  a_mat <- NULL
  if (free_scales || align_means) {
    # from the equal-SS start (and zero means), plus rstarts random starts
    # (the alignment loss has local optima); keep the best
    c_start <- lav_mat_rotate_mg_scales_start(a_list)
    a_start <- if (align_means) matrix(0, m, ngroups) else NULL
    c_mat <- NULL
    for (i in seq_len(rstarts + 1L)) {
      if (i > 1L) {
        if (free_scales) {
          z <- matrix(rnorm(m * ngroups, sd = 0.5), m, ngroups)
          p <- exp(z) / rowSums(exp(z))
          c_start <- sqrt(ngroups * p)
        }
        if (align_means) {
          a_start <- matrix(rnorm(m * ngroups, sd = 0.5), m, ngroups)
          if (intercepts == "fixed") {
            a_start[, 1L] <- 0
          }
        }
      }
      c_i <- lav_mat_rotate_mg_scales_update(
        lambda_list = a_list, c_mat = c_start,
        method_fname = method_fname, method_args = method_args,
        mg_agreement_method = mg_agreement_method,
        mg_agreement_epsilon = mg_agreement_epsilon,
        w = mg_agreement_weight, group_weights = group_weights,
        free_scales = free_scales, nu_list = nu_list, a_mat = a_start,
        mean_fixed_first = (intercepts == "fixed")
      )
      if (is.null(c_mat) || isTRUE(attr(c_i, "value") < attr(c_mat, "value"))) {
        c_mat <- c_i
      }
    }
    converged <- isTRUE(attr(c_mat, "converged"))
    q_current <- attr(c_mat, "value")
    a_mat <- attr(c_mat, "a_mat")
    attr(c_mat, "converged") <- NULL
    attr(c_mat, "value") <- NULL
    attr(c_mat, "a_mat") <- NULL
    if (free_scales) {
      lambda_list <- lapply(seq_len(ngroups), function(g) {
        t(t(a_list[[g]]) / c_mat[, g])
      })
    }
  } else {
    q_current <- as.numeric(lav_mat_rotate_mg_agreement(
      lambda_list = lambda_list, method_fname = method_fname,
      method_args = method_args, mg_agreement_method = mg_agreement_method,
      mg_agreement_epsilon = mg_agreement_epsilon,
      w = mg_agreement_weight, scale = TRUE,
      group_weights = group_weights, grad = FALSE
    ))
  }
  if (!converged) {
    lav_msg_warn(gettext(
      "the minimization over the scale split (and the factor means) of the
       agreement rotation did not converge"))
  }

  rot_list <- vector("list", length = ngroups)
  for (g in seq_len(ngroups)) {
    rott_inv <- solve(
      crossprod(a_list[[g]]),
      crossprod(a_list[[g]], lambda_list[[g]])
    )
    if (orthogonal && !free_scales) {
      rot_list[[g]] <- rott_inv
    } else {
      rot_list[[g]] <- solve(t(rott_inv))
    }
  }
  info <- list(
    algorithm = "mg-scales",
    iter = 0L,
    converged = converged,
    method.value = q_current
  )
  if (free_scales) {
    info$scales <- c_mat
  }
  if (align_means) {
    info$means <- a_mat
  }
  attr(rot_list, "info") <- info

  rot_list
}


# gradient projection algorithm for the multiple-group agreement objective
# (lav_mat_rotate_mg_agreement()): the rotation matrices of all groups are
# updated jointly, with the Barzilai-Borwein step size and the non-monotone
# line search of lav_mat_rotate_gpa(). a_list/phi_list are the per-group
# (pre-rotated) loading matrices and their factor correlations, as for
# lav_mat_rotate_pairwise_mg(); returns the list of rotation matrices
# relative to a_list (attribute "info" as for the other algorithms)
lav_mat_rotate_gpa_mg <- function(a_list = NULL, # original matrices
                                  phi_list = NULL, # factor correlations
                                  orthogonal = FALSE,
                                  method_fname = NULL, # crit function
                                  method_args = list(), # method args
                                  gpa_tol = 1e-5,
                                  gpa_fwindow = 10L,
                                  mg_agreement_weight = 0.50,
                                  mg_agreement_method = "procrustes",
                                  mg_agreement_epsilon = 0.01,
                                  group_weights = NULL,
                                  scale = TRUE,
                                  scales = "equal",
                                  intercepts = "none",
                                  nu_list = NULL,
                                  max_iter = 10000L) {
  ngroups <- length(a_list)
  m <- ncol(a_list[[1L]])
  if (is.null(phi_list) || orthogonal) {
    phi_list <- rep(list(diag(m)), ngroups)
  }

  # intercept alignment: the factor means a_mat (m x ngroups) are extra
  # (unconstrained) variables; "fixed" keeps the first group at zero
  align_means <- intercepts != "none" && !is.null(nu_list)
  a_mat <- NULL
  a_mask <- NULL
  if (align_means) {
    a_mat <- matrix(0, m, ngroups)
    a_mask <- matrix(1, m, ngroups)
    if (intercepts == "fixed") {
      a_mask[, 1L] <- 0
    }
  }

  # free scale split: the scales c_gj (m x ngroups) are extra variables,
  # one vector per factor on the sphere sum_g c_gj^2 = ngroups (mean
  # variance one); the loadings are lambda_g diag(1/c_g), the objective
  # is evaluated without the equal-SS rescaling. The projected gradient,
  # the BB step and the line search run on the stacked (rot, c) vector
  free_scales <- identical(scales, "free")
  if (free_scales) {
    scale <- FALSE
    c_mat <- lav_mat_rotate_mg_scales_start(a_list)
  } else {
    c_mat <- NULL
  }

  # GPA parametrization: lambda_g = astar_g %*% solve(t(rot_g)), with
  # rot_g' rot_g = phi_g; we start at rot_g = chol(phi_g), so that
  # lambda_g = a_g at the start
  rot_list <- lapply(phi_list, chol)
  astar_list <- lapply(seq_len(ngroups), function(g) {
    a_list[[g]] %*% t(rot_list[[g]])
  })

  objf <- function(lambda_list, a_mat = NULL, grad = FALSE) {
    mean_list <- NULL
    if (align_means) {
      mean_list <- lapply(seq_len(ngroups), function(g) a_mat[, g])
    }
    lav_mat_rotate_mg_agreement(
      lambda_list = lambda_list, method_fname = method_fname,
      method_args = method_args, mg_agreement_method = mg_agreement_method,
      mg_agreement_epsilon = mg_agreement_epsilon,
      w = mg_agreement_weight, scale = scale,
      group_weights = group_weights,
      nu_list = nu_list, mean_list = mean_list, grad = grad
    )
  }
  # gradient with respect to the factor means (masked)
  grad_means <- function(q_val) {
    gm <- matrix(unlist(attr(q_val, "grad_mean")), m, ngroups)
    gm * a_mask
  }
  rotate_all <- function(rot_list, c_mat = NULL) {
    lapply(seq_len(ngroups), function(g) {
      if (orthogonal) {
        lam <- astar_list[[g]] %*% rot_list[[g]]
      } else {
        lam <- t(solve(rot_list[[g]], t(astar_list[[g]])))
      }
      if (free_scales) {
        lam <- t(t(lam) / c_mat[, g])
      }
      lam
    })
  }
  # gradient of the objective with respect to rot_g, from the gradient gq
  # with respect to lambda_g (with free scales: lambda_list and gq_list
  # are for the scaled loadings; chain back through diag(1/c_g))
  grad_all <- function(rot_list, lambda_list, gq_list, c_mat = NULL) {
    lapply(seq_len(ngroups), function(g) {
      gq <- gq_list[[g]]
      lam <- lambda_list[[g]]
      if (free_scales) {
        gq <- t(t(gq) / c_mat[, g])
        lam <- t(t(lam) * c_mat[, g])
      }
      if (orthogonal) {
        crossprod(astar_list[[g]], gq)
      } else {
        -1 * solve(t(rot_list[[g]]), crossprod(gq, lam))
      }
    })
  }
  # gradient with respect to the scales: dF/dc_gj = -(lambda_g' gq_g)_jj
  # / c_gj (scaled quantities), projected on the sphere per factor
  grad_scales <- function(lambda_list, gq_list, c_mat) {
    gc <- vapply(seq_len(ngroups), function(g) {
      -1 * colSums(lambda_list[[g]] * gq_list[[g]]) / c_mat[, g]
    }, numeric(m))
    gc <- matrix(gc, nrow = m)
    gc - c_mat * (rowSums(c_mat * gc) / ngroups)
  }
  # projection onto the tangent space of the orthogonal/normal matrices
  project_all <- function(rot_list, grad_list) {
    lapply(seq_len(ngroups), function(g) {
      rot <- rot_list[[g]]
      grad <- grad_list[[g]]
      if (orthogonal) {
        mm <- crossprod(rot, grad)
        grad - rot %*% ((mm + t(mm)) / 2)
      } else {
        grad - t(t(rot) * colSums(rot * grad))
      }
    })
  }

  lambda_list <- rotate_all(rot_list, c_mat)
  q_1 <- objf(lambda_list, a_mat, grad = TRUE)
  gq_list <- attr(q_1, "grad")
  if (is.null(gq_list)) {
    lav_msg_stop(gettext(
      "mg_agreement_algorithm = \"gpa\" requires a rotation criterion with
       a gradient."))
  }
  q_current <- as.numeric(q_1)
  grad_list <- grad_all(rot_list, lambda_list, gq_list, c_mat)
  ga <- NULL
  if (align_means) {
    ga <- grad_means(q_1)
  }

  alpha <- 1
  f_hist <- numeric(max_iter + 1L)
  rot_prev <- NULL
  gp_prev <- NULL
  c_prev <- NULL
  gpc_prev <- NULL
  gpc <- NULL
  a_prev <- NULL
  ga_prev <- NULL
  converged <- FALSE
  degenerate <- FALSE
  for (iter in seq_len(max_iter + 1L)) {
    gp_list <- project_all(rot_list, grad_list)
    if (free_scales) {
      gpc <- grad_scales(lambda_list, gq_list, c_mat)
    }
    frob <- sqrt(sum(unlist(gp_list)^2) + sum(gpc^2) + sum(ga^2))
    f_hist[iter] <- q_current
    if (lav_verbose()) {
      cat(
        "  iter = ", sprintf("%4d", iter - 1L),
        " Q = ", sprintf("%13.11f", q_current),
        " frob.log10 = ", sprintf("%10.7f", log10(frob)),
        " alpha = ", sprintf("%9.7f", alpha), "\n"
      )
    }
    if (frob < gpa_tol) {
      converged <- TRUE
      break
    }

    # Barzilai-Borwein step size (over all groups jointly)
    if (!is.null(rot_prev)) {
      d_rot <- c(unlist(rot_list) - unlist(rot_prev), c_mat - c_prev,
                 a_mat - a_prev)
      d_gp <- c(unlist(gp_list) - unlist(gp_prev), gpc - gpc_prev,
                ga - ga_prev)
      if (sum(d_gp * d_gp) > 0) {
        alpha_new <- sum(d_rot * d_rot) / abs(sum(d_rot * d_gp))
        if (is.finite(alpha_new)) {
          alpha <- alpha_new
        }
        alpha <- max(1e-10, min(alpha, 20))
      }
    } else {
      alpha <- 2 * alpha
    }

    # non-monotone line search
    target_f <- max(f_hist[max(1L, iter - gpa_fwindow + 1L):iter])
    accepted <- FALSE
    for (i in seq_len(11L)) {
      rott_list <- vector("list", ngroups)
      for (g in seq_len(ngroups)) {
        x <- rot_list[[g]] - alpha * gp_list[[g]]
        if (orthogonal) {
          svd_out <- svd(x)
          rott_list[[g]] <- svd_out$u %*% t(svd_out$v)
        } else {
          v <- 1 / sqrt(colSums(x * x))
          rott_list[[g]] <- x %*% diag(v, nrow = m)
        }
      }
      c_t <- NULL
      if (free_scales) {
        x_c <- c_mat - alpha * gpc
        if (any(x_c <= 0)) {
          # a scale must stay positive: shorter step
          alpha <- alpha / 2
          next
        }
        c_t <- x_c * sqrt(ngroups / rowSums(x_c * x_c))
      }
      a_t <- NULL
      if (align_means) {
        a_t <- a_mat - alpha * ga
      }
      lambda_t <- tryCatch(rotate_all(rott_list, c_t),
                           error = function(e) NULL)
      if (is.null(lambda_t)) {
        alpha <- alpha / 2
        next
      }
      q_new <- objf(lambda_t, a_t, grad = TRUE)
      if (is.finite(q_new) &&
          as.numeric(q_new) < target_f - 0.5 * frob * frob * alpha) {
        accepted <- TRUE
        break
      }
      alpha <- alpha / 2
    }
    if (!accepted) {
      # no acceptable step: stop here (not converged)
      break
    }

    rot_prev <- rot_list
    gp_prev <- gp_list
    c_prev <- c_mat
    gpc_prev <- gpc
    a_prev <- a_mat
    ga_prev <- ga
    rot_list <- rott_list
    c_mat <- c_t
    a_mat <- a_t
    lambda_list <- lambda_t
    q_current <- as.numeric(q_new)
    gq_list <- attr(q_new, "grad")
    if (align_means) {
      ga <- grad_means(q_new)
    }

    # oblique: collapsed factors?
    if (!orthogonal) {
      for (g in seq_len(ngroups)) {
        phi_g <- crossprod(rot_list[[g]])
        if (max(abs(phi_g[lower.tri(phi_g)])) > 1 - 1e-6) {
          degenerate <- TRUE
        }
      }
      if (degenerate) {
        break
      }
    }
    grad_list <- tryCatch(grad_all(rot_list, lambda_list, gq_list, c_mat),
      error = function(e) NULL)
    if (is.null(grad_list)) {
      break
    }
  } # iter

  if (degenerate) {
    converged <- FALSE
    q_current <- Inf
    lav_msg_warn(gettext(
      "the multiple-group agreement rotation is degenerate: two or more
       factors have collapsed."))
  } else if (!converged) {
    lav_msg_warn(gettextf(
      "mg-gpa rotation algorithm did not converge after %s iterations",
      iter - 1L
    ))
  }

  # rotation matrices relative to a_list (as lav_mat_rotate_pairwise_mg());
  # with free scales lambda_list includes the scaling, and the general
  # form applies to orthogonal rotation too
  rot_out <- vector("list", length = ngroups)
  for (g in seq_len(ngroups)) {
    rott_inv <- solve(
      crossprod(a_list[[g]]),
      crossprod(a_list[[g]], lambda_list[[g]])
    )
    if (orthogonal && !free_scales) {
      rot_out[[g]] <- rott_inv
    } else {
      rot_out[[g]] <- solve(t(rott_inv))
    }
  }

  info <- list(
    algorithm = "mg-gpa",
    iter = iter - 1L,
    converged = converged,
    degenerate = degenerate,
    method.value = q_current
  )
  if (free_scales) {
    info$scales <- c_mat
  }
  if (align_means) {
    info$means <- a_mat
  }
  attr(rot_out, "info") <- info

  rot_out
}
