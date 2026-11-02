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
                                 mg_algorithm = "pairwise",
                                 mg_agreement_method = "procrustes",
                                 mg_agreement_weight = 0.5,
                                 mg_agreement_group_weights = NULL,
                                 reflect = TRUE, # reflect sign
                                 order_lv_by = "index", # how to order the lv's
                                 pairwise_accelerate = TRUE, # sweep extrap.
                                 gpa_tol = 0.00001, # stopping tol (gradient)
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
  # bifactor rotation: the first column is the general factor and stays in
  # place when the groups are aligned to the first group
  bifactor_flag <- method %in% c(
    "bi-geomin", "bigeomin", "bi-quartimin", "biquartimin"
  )

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
      max_iter = max_iter,
      group = g
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
  if (lav_verbose()) {
    cat("  mg-", mg_algorithm, " + agreement rotation:\n", sep = "")
  }
  if (mg_algorithm == "gpa") {
    rot_list <- lav_mat_rotate_gpa_mg(
      a_list = a_list,
      phi_list = phi_list,
      orthogonal = orthogonal,
      method_fname = method_fname,
      method_args = method_args,
      mg_agreement_method = mg_agreement_method,
      mg_agreement_weight = mg_agreement_weight,
      group_weights = mg_agreement_group_weights,
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
      mg_agreement_weight = mg_agreement_weight,
      group_weights = mg_agreement_group_weights,
      accelerate = pairwise_accelerate,
      gpa_tol = gpa_tol,
      tol = tol,
      max_iter = max_iter
    )
  }
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
                                       group_weights = NULL,
                                       scale = TRUE,
                                       max_iter = 1000L) {
  ngroups <- length(a_list)

  # number of columns
  m <- ncol(a_list[[1L]])

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

  # objective: the agreement criterion (+ per-group gradients)
  objf <- function(lambda_list, grad = FALSE) {
    lav_mat_rotate_mg_agreement(
      lambda_list = lambda_list, method_fname = method_fname,
      method_args = method_args, mg_agreement_method = mg_agreement_method,
      w = mg_agreement_weight, scale = scale,
      group_weights = group_weights, grad = grad
    )
  }

  res <- lav_mat_rotate_pairwise_core(
    lambda_list = lambda_list, phi_list = psi_list,
    orthogonal = orthogonal, objf = objf,
    use_grad = lav_mat_rotate_is_smooth(method_fname),
    accelerate = accelerate,
    gpa_tol = gpa_tol, tol = tol, max_iter = max_iter
  )
  lambda_list <- res$lambda_list

  # warn if no convergence
  if (!res$converged) {
    lav_msg_warn(gettextf(
      "mg-pairwise rotation algorithm did not converge after %s iterations",
      res$iter
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
    iter = res$iter,
    converged = res$converged,
    method.value = res$q
  )

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
                                  group_weights = NULL,
                                  scale = TRUE,
                                  max_iter = 10000L) {
  ngroups <- length(a_list)
  m <- ncol(a_list[[1L]])
  if (is.null(phi_list) || orthogonal) {
    phi_list <- rep(list(diag(m)), ngroups)
  }

  # GPA parametrization: lambda_g = astar_g %*% solve(t(rot_g)), with
  # rot_g' rot_g = phi_g; we start at rot_g = chol(phi_g), so that
  # lambda_g = a_g at the start
  rot_list <- lapply(phi_list, chol)
  astar_list <- lapply(seq_len(ngroups), function(g) {
    a_list[[g]] %*% t(rot_list[[g]])
  })

  objf <- function(lambda_list, grad = FALSE) {
    lav_mat_rotate_mg_agreement(
      lambda_list = lambda_list, method_fname = method_fname,
      method_args = method_args, mg_agreement_method = mg_agreement_method,
      w = mg_agreement_weight, scale = scale,
      group_weights = group_weights, grad = grad
    )
  }
  rotate_all <- function(rot_list) {
    lapply(seq_len(ngroups), function(g) {
      if (orthogonal) {
        astar_list[[g]] %*% rot_list[[g]]
      } else {
        t(solve(rot_list[[g]], t(astar_list[[g]])))
      }
    })
  }
  # gradient of the objective with respect to rot_g, from the gradient gq
  # with respect to lambda_g
  grad_all <- function(rot_list, lambda_list, gq_list) {
    lapply(seq_len(ngroups), function(g) {
      if (orthogonal) {
        crossprod(astar_list[[g]], gq_list[[g]])
      } else {
        -1 * solve(t(rot_list[[g]]), crossprod(gq_list[[g]], lambda_list[[g]]))
      }
    })
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

  lambda_list <- rotate_all(rot_list)
  q_1 <- objf(lambda_list, grad = TRUE)
  gq_list <- attr(q_1, "grad")
  if (is.null(gq_list)) {
    lav_msg_stop(gettext(
      "mg_agreement_algorithm = \"gpa\" requires a rotation criterion with
       a gradient."))
  }
  q_current <- as.numeric(q_1)
  grad_list <- grad_all(rot_list, lambda_list, gq_list)

  alpha <- 1
  f_hist <- numeric(max_iter + 1L)
  rot_prev <- NULL
  gp_prev <- NULL
  converged <- FALSE
  degenerate <- FALSE
  for (iter in seq_len(max_iter + 1L)) {
    gp_list <- project_all(rot_list, grad_list)
    frob <- sqrt(sum(unlist(gp_list)^2))
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
      d_rot <- unlist(rot_list) - unlist(rot_prev)
      d_gp <- unlist(gp_list) - unlist(gp_prev)
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
      lambda_t <- tryCatch(rotate_all(rott_list), error = function(e) NULL)
      if (is.null(lambda_t)) {
        alpha <- alpha / 2
        next
      }
      q_new <- objf(lambda_t, grad = TRUE)
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
    rot_list <- rott_list
    lambda_list <- lambda_t
    q_current <- as.numeric(q_new)
    gq_list <- attr(q_new, "grad")

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
    grad_list <- tryCatch(grad_all(rot_list, lambda_list, gq_list),
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

  # rotation matrices relative to a_list (as lav_mat_rotate_pairwise_mg())
  rot_out <- vector("list", length = ngroups)
  for (g in seq_len(ngroups)) {
    rott_inv <- solve(
      crossprod(a_list[[g]]),
      crossprod(a_list[[g]], lambda_list[[g]])
    )
    if (orthogonal) {
      rot_out[[g]] <- rott_inv
    } else {
      rot_out[[g]] <- solve(t(rott_inv))
    }
  }

  attr(rot_out, "info") <- list(
    algorithm = "mg-gpa",
    iter = iter - 1L,
    converged = converged,
    degenerate = degenerate,
    method.value = q_current
  )

  rot_out
}
