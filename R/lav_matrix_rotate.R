# rotation algorithms
#
# YR  3 April 2019 -- gradient projection algorithm
# YR 21 April 2019 -- pairwise rotation algorithm
# YR 11 May   2020 -- order.idx is done in rotation matrix
#                     (suggested by Florian Scharf)
# YR 02 June  2024 -- add group argument, so target and target_mask can
#                     be a list

# main function to rotate a single matrix 'A'
lav_mat_rotate <- function(a = NULL, # original matrix
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
                              init_rot = NULL, # initial rotation matrix
                              init_rot_check = TRUE, # check if init ROT is ok
                              rstarts = 100L, # number of random starts
                              row_weights = "default", # row weighting
                              std_ov = FALSE, # rescale ov
                              ov_var = NULL, # ov variances
                              algorithm = "gpa", # rotation algorithm
                              pairwise_accelerate = TRUE, # sweep extrapolation
                              reflect = TRUE, # reflect sign
                              order_lv_by = "index", # how to order the lv's
                              gpa_tol = 0.00001, # stopping tol gpa
                              gpa_algorithm = "bb", # step size/projection gpa
                              gpa_fwindow = 0L, # line search window gpa
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
    if (!lav_mat_rotate_check(init_rot, orthogonal = orthogonal)) {
      lav_msg_stop(gettext("init.ROT does not look like a rotation matrix"))
    }
  }

  # determine method function name + validate target/target_mask
  tmp_resolve <- lav_mat_rotate_resolve_method(
    method = method, method_args = method_args, p = p, m = m, group = group
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

  # check gpa_algorithm here as well: inside the random-starts loop, a
  # user error would be caught by the failed-start safety net (and turn
  # into a misleading "all random starts failed" message)
  if (algorithm == "gpa") {
    gpa_algorithm <- tolower(gpa_algorithm)
    if (!gpa_algorithm %in% c("legacy", "bb", "cayley")) {
      lav_msg_stop(gettext("gpa_algorithm must be legacy, bb or cayley"))
    }
    if (gpa_algorithm == "cayley" && !orthogonal) {
      lav_msg_stop(gettext(
        "gpa_algorithm = \"cayley\" is only available for orthogonal rotation."
      ))
    }
  }



  # 1. compute row weigths

  # 1.a cov -> cor?
  if (std_ov) {
    a <- a * 1 / sqrt(ov_var)
  }

  if (row_weights == "none") {
    weights <- rep(1.0, p)
  } else if (row_weights == "kaiser") {
    weights <- lav_mat_rotate_kaiser_weights(a)
  } else if (row_weights == "cureton-mulaik") {
    weights <- lav_mat_rotate_cm_weights(a)
  } else {
    lav_msg_stop(gettext("row_weights can be none, kaiser or cureton-mulaik"))
  }
  a <- a * weights


  # 2. rotate

  # multiple random starts?
  if (rstarts > 0L) {
    rep_1 <- sapply(seq_len(rstarts), function(rep) {
      # random start (always orthogonal)
      init_rot <- lav_mat_rotate_gen(m = m, orthogonal = TRUE)
      # init.ROT <- lav_mat_rotate_gen(M = M, orthogonal = orthogonal)

      if (lav_verbose()) {
        cat("\n")
        cat("rstart = ", sprintf("%4d", rep), " start:\n")
      }


      # choose rotation algorithm; if a start fails (eg a singular
      # transformation matrix for a factor-collapsing criterion), discard
      # it (criterion = +Inf) instead of aborting the whole rotation
      res <- tryCatch({
        if (algorithm == "gpa") {
          rot <- lav_mat_rotate_gpa(
            a = a, orthogonal = orthogonal,
            init_rot = init_rot,
            method_fname = method_fname,
            method_args = method_args,
            gpa_tol = gpa_tol,
            gpa_algorithm = gpa_algorithm,
            gpa_fwindow = gpa_fwindow,
            max_iter = max_iter
          )
        } else if (algorithm == "pairwise") {
          rot <- lav_mat_rotate_pairwise(
            a = a,
            orthogonal = orthogonal,
            init_rot = init_rot,
            method_fname = method_fname,
            method_args = method_args,
            accelerate = pairwise_accelerate,
            gpa_tol = gpa_tol,
            tol = tol,
            max_iter = max_iter
          )
        }
        info <- attr(rot, "info")
        attr(rot, "info") <- NULL
        c(info$method.value, info$iter, as.numeric(isTRUE(info$converged)),
          lav_mat_vec(rot))
      }, error = function(e) c(Inf, NA_real_, 0, rep(NA_real_, m * m)))

      if (lav_verbose()) {
        cat(
          "rstart = ", sprintf("%4d", rep),
          " end; current crit = ", sprintf("%17.15f", res[1]), "\n"
        )
      }
      res
    })
    best_idx <- which.min(rep_1[1, ])
    if (!is.finite(rep_1[1, best_idx])) {
      lav_msg_stop(gettext(
        "rotation failed or degenerated (collapsed factors) for all random
        starts. The rotation criterion may be unbounded for these data;
        consider another rotation method, or fewer factors."))
    }
    n_bad <- sum(!is.finite(rep_1[1, ]))
    if (n_bad > 0L) {
      lav_msg_warn(gettextf(
        "rotation failed or degenerated (collapsed factors) for %1$s out of
        %2$s random starts; these starts were discarded.",
        n_bad, rstarts
      ))
    }
    rot <- matrix(rep_1[-(1:3), best_idx], nrow = m, ncol = m)
    # algorithm information of the best start; REP has one column per
    # start: the criterion value, the number of iterations, the
    # convergence flag (0/1) and the vectorized rotation matrix
    info <- list(
      algorithm = algorithm,
      iter = as.integer(rep_1[2, best_idx]),
      converged = isTRUE(rep_1[3, best_idx] == 1),
      method.value = rep_1[1, best_idx]
    )
    if (keep_rep) {
      info$REP <- rep_1
    }
  } else if (algorithm != "none") {
    # initial rotation matrix
    if (is.null(init_rot)) {
      init_rot <- diag(m)
    }

    # Gradient Projection Algorithm
    if (algorithm == "gpa") {
      rot <- lav_mat_rotate_gpa(
        a = a, orthogonal = orthogonal,
        init_rot = init_rot,
        method_fname = method_fname,
        method_args = method_args,
        gpa_tol = gpa_tol,
        gpa_algorithm = gpa_algorithm,
        gpa_fwindow = gpa_fwindow,
        max_iter = max_iter
      )
    } else if (algorithm == "pairwise") {
      rot <- lav_mat_rotate_pairwise(
        a = a,
        orthogonal = orthogonal,
        init_rot = init_rot,
        method_fname = method_fname,
        method_args = method_args,
        accelerate = pairwise_accelerate,
        gpa_tol = gpa_tol,
        tol = tol,
        max_iter = max_iter
      )
    }
    info <- attr(rot, "info")
    attr(rot, "info") <- NULL
    if (isTRUE(info$degenerate)) {
      lav_msg_warn(gettext(
        "the rotated solution is degenerate: two or more factors have
        collapsed (correlation of (nearly) one). The rotation criterion may
        be unbounded for these data; consider another rotation method, or
        fewer factors."))
    }
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
    out <- lav_mat_rotate_promax(
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

  # 4.b reorder the columns; for a bifactor rotation, the first column is
  #     the general factor and stays in place: the ordering rule is applied
  #     to the specific factors only
  bifactor_flag <- method %in% c(
    "bi-geomin", "bigeomin", "bi-quartimin", "biquartimin"
  )
  order_idx <- lav_efa_order_idx(mm_lambda, order_lv_by,
    general_first = bifactor_flag
  )

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
# - new in 0.7-2 (following GPArotation 2026-6): the step size alpha can be
#   set using the Barzilai-Borwein (1988) method (gpa_algorithm = "bb", now
#   the default), combined with a non-monotone line search (Grippo, Lampariello
#   & Lucidi, 1986) over a window of gpa_fwindow previous criterion values;
#   for orthogonal rotation only, the projection onto the orthogonal manifold
#   can be done using the Cayley transform (gpa_algorithm = "cayley") instead
#   of the SVD; gpa_algorithm = "legacy" gives the original (pre 0.7-2)
#   behavior: alpha doubling + monotone line search + SVD projection
#
lav_mat_rotate_gpa <- function(a = NULL, # original matrix
                                  orthogonal = FALSE, # default is oblique
                                  init_rot = NULL, # initial rotation
                                  method_fname = NULL, # criterion function
                                  method_args = list(), # optional method args
                                  gpa_tol = 0.00001,
                                  gpa_algorithm = "bb",
                                  gpa_fwindow = 0L, # 0 = auto
                                  max_iter = 10000L) {
  # number of columns
  m <- ncol(a)

  # check gpa_algorithm
  gpa_algorithm <- tolower(gpa_algorithm)
  if (!gpa_algorithm %in% c("legacy", "bb", "cayley")) {
    lav_msg_stop(gettext("gpa_algorithm must be legacy, bb or cayley"))
  }
  if (gpa_algorithm == "cayley" && !orthogonal) {
    lav_msg_stop(gettext(
      "gpa_algorithm = \"cayley\" is only available for orthogonal rotation."))
  }

  # width of the (non-monotone) line search window; 0 = auto
  gpa_fwindow <- as.integer(gpa_fwindow)
  if (gpa_fwindow < 1L) {
    gpa_fwindow <- if (gpa_algorithm == "legacy") 1L else 10L
  }

  # maximum number of step-halvings in the line search
  if (gpa_algorithm == "legacy") {
    ls_max_iter <- 1000L
  } else {
    ls_max_iter <- 11L
  }

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
  degenerate <- FALSE
  # history of criterion values (used by the non-monotone line search)
  f_hist <- numeric(max_iter + 1L)
  rot_prev <- NULL
  gp_prev <- NULL
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

    # store current criterion value
    f_hist[iter] <- q_current

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

    # determine (initial) step size alpha
    if (gpa_algorithm != "legacy" && !is.null(rot_prev)) {
      # Barzilai-Borwein step size; alpha is set directly, and we enter
      # the line search without doubling
      d_rot <- rot - rot_prev
      d_gp <- gp - gp_prev
      if (sum(d_gp * d_gp) > 0) {
        alpha_new <- sum(d_rot * d_rot) / abs(sum(d_rot * d_gp))
        if (is.finite(alpha_new)) {
          alpha <- alpha_new
        }
        alpha <- max(1e-10, min(alpha, 20))
      }
    } else {
      # original behavior: double the (accepted) step size
      alpha <- 2 * alpha
    }

    # reference value for the (non-monotone) line search: the worst
    # criterion value over the last gpa_fwindow iterations
    # (if gpa_fwindow == 1, this is just q_current: monotone line search)
    target_f <- max(f_hist[max(1L, iter - gpa_fwindow + 1L):iter])

    ls_ok <- FALSE # did any trial step produce a usable ROTt?
    for (i in seq_len(ls_max_iter)) {

      if (orthogonal && gpa_algorithm == "cayley") {
        # use the Cayley transform to move over the manifold of
        # orthogonal matrices along a descent curve
        w_skew <- (alpha / 2) * (tcrossprod(gp, rot) - tcrossprod(rot, gp))
        rott <- solve(diag(m) + w_skew, rot - w_skew %*% rot)
      } else {
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
      }

      # rotate again
      if (orthogonal) {
        mm_lambda <- a %*% rott
      } else {
        # for factor-collapsing criteria (eg oblimin with gamma > 0), a
        # trial ROTt can be (numerically) singular; treat this as a
        # failed step: halve alpha and try again (as alpha decreases,
        # ROTt approaches the current -- nonsingular -- ROT)
        mm_lambda <- tryCatch(t(solve(rott, at)), error = function(e) NULL)
        if (is.null(mm_lambda)) {
          alpha <- alpha / 2
          next
        }
      }
      rott_ok <- rott
      ls_ok <- TRUE

      # evaluate criterion
      q_new <- do.call(method_fname, c(
        list(mm_lambda = mm_lambda),
        method_args, list(grad = TRUE)
      ))
      gq <- attr(q_new, "grad")
      attr(q_new, "grad") <- NULL

      # check stopping criterion
      if (q_new < target_f - 0.5 * frob * frob * alpha) {
        break
      } else {
        alpha <- alpha / 2
      }

      if (i == ls_max_iter && gpa_algorithm == "legacy") {
        lav_msg_warn(gettext("half-stepping failed in GPA"))
      }
    }

    # no trial step was usable (all singular): give up (not converged)
    if (!ls_ok) {
      break
    }

    # update
    rot_prev <- rot
    gp_prev <- gp
    rot <- rott_ok
    q_current <- q_new

    # oblique: stop early if two factors have (numerically) collapsed;
    # some criteria (eg oblimin with gamma > 0) are unbounded below for
    # oblique rotation, and the algorithm then chases -Inf towards a
    # singular transformation matrix
    if (!orthogonal) {
      phi_1 <- crossprod(rot)
      if (max(abs(phi_1[lower.tri(phi_1)])) > 1 - 1e-6) {
        degenerate <- TRUE
        break
      }
    }

    if (orthogonal) {
      grad <- crossprod(a, gq)
    } else {
      # the (transposed) solve may still fail if the accepted ROT sits
      # right at the singularity threshold; stop here (not converged)
      grad <- tryCatch(
        -1 * solve(t(rot), crossprod(gq, mm_lambda)),
        error = function(e) NULL
      )
      if (is.null(grad)) {
        break
      }
    }
  } # iter

  # a non-converged oblique run that ends with (nearly) perfectly
  # correlated factors, or with a criterion value that has diverged far
  # below its starting value, is a factor-collapse trajectory (for an
  # unbounded criterion, the projected gradient never becomes small):
  # flag it
  if (!degenerate && !orthogonal && !converged) {
    phi_1 <- crossprod(rot)
    if (max(abs(phi_1[lower.tri(phi_1)])) > 0.999 ||
        q_current < -1000 * (1 + abs(f_hist[1]))) {
      degenerate <- TRUE
    }
  }

  # degenerate (collapsed) solution: not converged, and make sure this
  # solution is never selected when multiple random starts are used
  if (degenerate) {
    converged <- FALSE
    q_current <- Inf
  } else if (!converged) {
    # warn if no convergence
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
    degenerate = degenerate,
    method.value = q_current
  )

  attr(rot, "info") <- info

  rot
}


# pairwise (plane) rotation algorithm
#
# based on Kaiser's (1959) algorithm and Jennrich and Sampson (1966)
# algorithm, made generic (any criterion) as in Browne (2001)
#
# - orthogonal: rotate one pair of columns (= plane) at a time
# - oblique: rotate 1 factor in one pair of columns (= plane) at a time
#            note: in the oblique case, (1,2) is not the same as (2,1)
# - new in 0.7-3: the plane minimizations use the gradient of the
#   criterion (see lav_mat_rotate_pairwise_core()); the derivative-free
#   line search (optimize()) is only used for a non-smooth criterion
#   (simplimax) or a criterion without a gradient
#
lav_mat_rotate_pairwise <- function(a = NULL, # original matrix
                                    orthogonal = FALSE,
                                    init_rot = NULL,
                                    method_fname = NULL, # crit function
                                    method_args = list(), # method args
                                    accelerate = TRUE,
                                    gpa_tol = 1e-5,
                                    tol = 1e-8,
                                    max_iter = 1000L) {
  # number of columns
  m <- ncol(a)

  # initial mm_lambda + PHI
  if (is.null(init_rot)) {
    mm_lambda <- a
    phi <- diag(m)
  } else if (orthogonal) {
    mm_lambda <- a %*% init_rot
    phi <- diag(m)
  } else {
    mm_lambda <- t(solve(init_rot, t(a)))
    phi <- crossprod(init_rot)
  }

  # objective: the criterion (and its gradient) for a 1-element list
  objf <- function(lambda_list, grad = FALSE) {
    q <- do.call(method_fname, c(
      list(mm_lambda = lambda_list[[1L]]),
      method_args, list(grad = grad)
    ))
    if (grad) {
      gq <- attr(q, "grad")
      attr(q, "grad") <- if (is.null(gq)) NULL else list(gq)
    }
    q
  }

  res <- lav_mat_rotate_pairwise_core(
    lambda_list = list(mm_lambda), phi_list = list(phi),
    orthogonal = orthogonal, objf = objf,
    use_grad = lav_mat_rotate_is_smooth(method_fname),
    plane_quartic = lav_mat_rotate_is_plane_quartic(method_fname),
    accelerate = accelerate,
    gpa_tol = gpa_tol, tol = tol, max_iter = max_iter
  )
  mm_lambda <- res$lambda_list[[1L]]

  # warn if no convergence
  if (!res$converged) {
    lav_msg_warn(gettextf(
      "pairwise rotation algorithm did not converge after %s iterations",
      res$iter
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
    iter = res$iter,
    converged = res$converged,
    method.value = res$q
  )

  attr(rot, "info") <- info

  rot
}


# plane-sweep engine shared by lav_mat_rotate_pairwise() (single group) and
# lav_mat_rotate_pairwise_mg() (multiple groups + agreement criterion)
#
# lambda_list: the current loading matrices (one per group); phi_list: the
# corresponding factor correlation matrices (oblique rotation only);
# objf(lambda_list, grad): the objective, with (if grad = TRUE) the
# attribute "grad" holding the list of gradients with respect to each
# loading matrix, or NULL when no gradient is available
#
# Gradient-driven version (the default for a smooth criterion):
# - the directional derivative g0 of the objective along a plane rotation
#   is an inner product of the current gradient with the plane direction,
#   so it costs nothing; planes with |g0| below skip_tol are skipped
# - for an orthogonal rotation with a criterion whose plane objective is
#   a + b cos(4t) + c sin(4t) (plane_quartic = TRUE; the orthomax,
#   Crawford-Ferguson and oblimin families, see
#   lav_mat_rotate_is_plane_quartic()), the exact plane minimizer follows
#   from the values at t = 0 and t = +/- pi/8 (three evaluations per plane)
# - otherwise the plane objective is minimized along the plane by a
#   safeguarded secant (Illinois) iteration on the directional derivative:
#   a first (Newton) trial step from the curvature remembered for this
#   plane, bracket expansion until the derivative changes sign, then
#   secant steps; every evaluation yields the objective AND its gradient,
#   so a plane typically costs 3 to 6 evaluations instead of the 15 to 25
#   of a derivative-free line search; the best evaluated point is taken,
#   and only if it improves the objective
# - convergence: the largest |g0| over all planes (and groups) below
#   gpa_tol; the plane directions span the tangent space of the orthogonal
#   (or normal) matrices, so this is the same projected-gradient criterion
#   as used by the GPA algorithm
# - sweep-level acceleration (accelerate = TRUE): plane sweeps converge
#   linearly, so after each sweep the cumulative transformation matrix
#   M_k (lambda_k = lambda_0 M_k, per group) is extrapolated along the
#   last sweep displacement, M_k + s (M_k - M_{k-1}), with the Aitken
#   factor s = r / (1 - r) estimated from the contraction ratio r of the
#   last two displacements; the extrapolated point is projected back onto
#   the manifold (orthogonal: polar projection; oblique: unit-variance
#   factors) and accepted only if it lowers the objective (one extra
#   evaluation per sweep)
#
# Derivative-free fallback (use_grad = FALSE, or no gradient available):
# Brent's method (optimize()) per plane, as in the original algorithm;
# convergence when the change in the objective over a sweep is below tol.
# This is the path for a non-smooth criterion (simplimax), where the
# derivative-based iteration can stall at a kink
lav_mat_rotate_pairwise_core <- function(lambda_list = NULL,
                                         phi_list = NULL,
                                         orthogonal = FALSE,
                                         objf = NULL,
                                         use_grad = TRUE,
                                         plane_quartic = FALSE,
                                         accelerate = TRUE,
                                         gpa_tol = 1e-5,
                                         tol = 1e-8,
                                         max_iter = 1000L,
                                         plane_max_eval = 12L) {
  ngroups <- length(lambda_list)
  m <- ncol(lambda_list[[1L]])

  plane_quartic <- plane_quartic && orthogonal

  # starting point + cumulative transformation per group (for the
  # sweep-level extrapolation): lambda_g = lambda_0_g %*% m_g, and
  # phi_g = m_g^-1 phi_0_g m_g^-T
  lambda_0 <- lambda_list
  phi_0 <- phi_list
  m_list <- rep(list(diag(m)), ngroups)
  m_prev <- NULL # m_list at the end of the previous sweep
  d_prev_norm <- NA_real_ # norm of the previous sweep displacement
  n_extrap <- 0L

  # planes whose directional derivative is below skip_tol are left alone;
  # convergence is declared at gpa_tol (so that the final sweep still
  # polishes the planes between skip_tol and gpa_tol)
  skip_tol <- gpa_tol / 100

  # plane combinations
  if (orthogonal) {
    plane <- utils::combn(m, 2L)
  } else {
    tmp <- utils::combn(m, 2L)
    plane <- cbind(tmp, tmp[c(2L, 1L), , drop = FALSE])
  }
  nplanes <- ncol(plane)

  # largest plane step (angle, or shear)
  t_max <- if (orthogonal) pi / 4 else 1

  # apply a plane rotation to the two affected columns only; returns the
  # updated matrix and (oblique) the scaling factor gamma of column c1
  plane_apply <- function(mm_lambda, c1, c2, t, phi12 = 0) {
    l1 <- mm_lambda[, c1]
    l2 <- mm_lambda[, c2]
    if (orthogonal) {
      cs <- cos(t)
      sn <- sin(t)
      mm_lambda[, c1] <- cs * l1 - sn * l2
      mm_lambda[, c2] <- sn * l1 + cs * l2
      gam <- 1
    } else {
      gam <- sqrt(abs(1 + 2 * t * phi12 + t * t))
      mm_lambda[, c1] <- gam * l1
      mm_lambda[, c2] <- l2 - t * l1
    }
    list(mm_lambda = mm_lambda, gam = gam)
  }

  # directional derivative of the objective along the plane at step t,
  # given the gradient gq at that point and the columns l1/l2 BEFORE the
  # step
  plane_deriv <- function(gq, l1, l2, c1, c2, t, phi12 = 0, gam = 1) {
    if (orthogonal) {
      cs <- cos(t)
      sn <- sin(t)
      sum(gq[, c1] * (-sn * l1 - cs * l2)) + sum(gq[, c2] * (cs * l1 - sn * l2))
    } else {
      sum(gq[, c1] * l1) * (phi12 + t) / gam - sum(gq[, c2] * l1)
    }
  }

  # oblique: update PHI of group g after a plane step
  phi_update <- function(phi, c1, c2, t, gam) {
    phi[c1, ] <- (1 / gam) * phi[c1, ] + (t / gam) * phi[c2, ]
    phi[, c1] <- phi[c1, ]
    phi[c1, c1] <- 1
    phi
  }

  # initial objective (+ gradients)
  q_0 <- objf(lambda_list, grad = use_grad)
  grad_list <- attr(q_0, "grad")
  use_grad <- use_grad && !is.null(grad_list)
  q_current <- as.numeric(q_0)
  neval <- 1L

  if (lav_verbose()) {
    cat(
      "  iter = ", sprintf("%4d", 0L),
      " Q = ", sprintf("%13.11f", q_current), "\n"
    )
  }

  # per-plane (and per-group) curvature estimates (NA = unknown)
  curv_pl <- matrix(as.numeric(NA), nrow = nplanes, ncol = ngroups)

  converged <- FALSE
  for (iter in seq_len(max_iter)) {
    q_old <- q_current
    g_max <- 0
    n_accepted <- 0L

    for (g in seq_len(ngroups)) {
      for (pl in seq_len(nplanes)) {
        c1 <- plane[1L, pl]
        c2 <- plane[2L, pl]
        mm_lambda <- lambda_list[[g]]
        l1 <- mm_lambda[, c1]
        l2 <- mm_lambda[, c2]
        phi12 <- if (orthogonal) 0 else phi_list[[g]][c1, c2]

        if (!use_grad) {
          # derivative-free: Brent's method on the plane objective
          objf_plane <- function(t) {
            lambda_list[[g]] <- plane_apply(mm_lambda, c1, c2, t,
              phi12 = phi12)$mm_lambda
            as.numeric(objf(lambda_list, grad = FALSE))
          }
          out <- optimize(
            f = objf_plane, interval = c(-t_max, t_max),
            maximum = FALSE, tol = .Machine$double.eps^0.25
          )
          t <- out$minimum
          step <- plane_apply(mm_lambda, c1, c2, t, phi12 = phi12)
          lambda_list[[g]] <- step$mm_lambda
          m_list[[g]] <- plane_apply(m_list[[g]], c1, c2, t,
            phi12 = phi12)$mm_lambda
          if (!orthogonal) {
            phi_list[[g]] <- phi_update(phi_list[[g]], c1, c2, t, step$gam)
          }
          next
        }

        # directional derivative at the current point (free)
        g0 <- plane_deriv(grad_list[[g]], l1, l2, c1, c2, t = 0,
          phi12 = phi12, gam = 1)
        g_max <- max(g_max, abs(g0))
        if (abs(g0) < skip_tol) {
          next
        }

        # evaluate the objective (+ gradient) and the directional
        # derivative at plane step t
        eval_t <- function(t) {
          step <- plane_apply(mm_lambda, c1, c2, t, phi12 = phi12)
          lambda_list[[g]] <- step$mm_lambda
          q_t <- objf(lambda_list, grad = TRUE)
          gq <- attr(q_t, "grad")
          list(
            t = t, q = as.numeric(q_t), grad = gq, gam = step$gam,
            mm_lambda = step$mm_lambda,
            g = if (is.null(gq)) NA_real_ else {
              plane_deriv(gq[[g]], l1, l2, c1, c2, t = t,
                phi12 = phi12, gam = step$gam)
            }
          )
        }

        best <- NULL
        if (plane_quartic) {
          # the plane objective is a + b cos(4t) + c sin(4t): the value and
          # the derivative at t = 0 (both known: a + b and 4c) and the value
          # at t = pi/8 (a + c) give a, b, c; the minimizer is
          # 4t = atan2(c, b) + pi
          lambda_list[[g]] <- plane_apply(mm_lambda, c1, c2, pi / 8)$mm_lambda
          f_plus <- as.numeric(objf(lambda_list, grad = FALSE))
          neval <- neval + 1L
          c_q <- g0 / 4
          a_q <- f_plus - c_q
          b_q <- q_current - a_q
          t_q <- (atan2(c_q, b_q) + pi) / 4
          t_q <- ((t_q + pi / 4) %% (pi / 2)) - pi / 4
          ps <- eval_t(t_q)
          neval <- neval + 1L
          if (is.finite(ps$q) && ps$q < q_current) {
            best <- ps
          }
          # curvature at the minimizer: 16 sqrt(b^2 + c^2)
          curv_pl[pl, g] <- 16 * sqrt(b_q * b_q + c_q * c_q)
        } else {
        # 1. first trial step: Newton step from the remembered curvature,
        #    or a fixed fraction of the maximal step in the descent direction
        curv <- curv_pl[pl, g]
        if (is.finite(curv) && curv > 0) {
          t1 <- -g0 / curv
        } else {
          t1 <- -sign(g0) * 0.1 * t_max
        }
        t1 <- max(-t_max, min(t_max, t1))
        if (abs(t1) < 1e-6) {
          t1 <- -sign(g0) * 1e-6
        }
        pa <- list(t = 0, q = q_current, g = g0)
        pb <- eval_t(t1)
        neval <- neval + 1L
        if (is.finite(pb$q) && pb$q < q_current) {
          best <- pb
        }
        n_ev <- 1L

        # 2. bracket: expand (doubling) until the derivative changes sign
        while (is.finite(pb$g) && sign(pb$g) == sign(g0) &&
               abs(pb$t) < t_max && n_ev < plane_max_eval) {
          pa <- pb
          t_new <- max(-t_max, min(t_max, 2 * pb$t))
          pb <- eval_t(t_new)
          neval <- neval + 1L
          n_ev <- n_ev + 1L
          if (is.finite(pb$q) && (is.null(best) || pb$q < best$q)) {
            best <- pb
          }
        }

        # 3. secant (Illinois) iterations on the directional derivative
        #    within the bracket [pa, pb]
        if (is.finite(pb$g) && sign(pb$g) != sign(pa$g)) {
          side <- 0L
          while (n_ev < plane_max_eval) {
            if (abs(pb$g) < skip_tol || abs(pa$g) < skip_tol ||
                abs(pb$t - pa$t) < 1e-10) {
              break
            }
            t_s <- (pa$t * pb$g - pb$t * pa$g) / (pb$g - pa$g)
            if (!is.finite(t_s)) {
              break
            }
            ps <- eval_t(t_s)
            neval <- neval + 1L
            n_ev <- n_ev + 1L
            if (!is.finite(ps$q) || !is.finite(ps$g)) {
              break
            }
            if (is.null(best) || ps$q < best$q) {
              best <- ps
            }
            if (sign(ps$g) == sign(pb$g)) {
              pb <- ps
              if (side == -1L) {
                pa$g <- pa$g / 2
              }
              side <- -1L
            } else {
              pa <- ps
              if (side == 1L) {
                pb$g <- pb$g / 2
              }
              side <- 1L
            }
          }
        }

        # remember the curvature for the next sweep
        if (is.finite(pb$g) && is.finite(pa$g) && abs(pb$t - pa$t) > 0) {
          curv_pl[pl, g] <- (pb$g - pa$g) / (pb$t - pa$t)
        }
        } # secant plane search

        # 4. accept the best evaluated point (if it improves the objective)
        if (is.null(best) || !(best$q < q_current)) {
          lambda_list[[g]] <- mm_lambda
          next
        }
        lambda_list[[g]] <- best$mm_lambda
        m_list[[g]] <- plane_apply(m_list[[g]], c1, c2, best$t,
          phi12 = phi12)$mm_lambda
        q_current <- best$q
        grad_list <- best$grad
        n_accepted <- n_accepted + 1L
        if (!orthogonal) {
          phi_list[[g]] <- phi_update(phi_list[[g]], c1, c2, best$t, best$gam)
        }
      } # planes
    } # groups

    if (!use_grad) {
      q_current <- as.numeric(objf(lambda_list, grad = FALSE))
      neval <- neval + 1L
    }
    diff <- abs(q_old - q_current)

    if (lav_verbose()) {
      cat(
        "  iter = ", sprintf("%4d", iter),
        " Q = ", sprintf("%13.11f", q_current),
        " change = ", sprintf("%13.11f", diff),
        if (use_grad) paste0(" max|g| = ", sprintf("%10.3e", g_max)) else "",
        "\n"
      )
    }

    if (!is.finite(diff)) {
      # something went wrong ... bail out
      break
    }
    if (use_grad) {
      if (g_max < gpa_tol) {
        converged <- TRUE
        break
      }
      if (n_accepted == 0L) {
        # no plane could improve the objective: stop here
        break
      }
    } else if (diff < tol) {
      converged <- TRUE
      break
    }

    # sweep-level extrapolation (gradient path only)
    if (use_grad && accelerate) {
      if (!is.null(m_prev)) {
        d_list <- lapply(seq_len(ngroups), function(g) {
          m_list[[g]] - m_prev[[g]]
        })
        d_norm <- sqrt(sum(unlist(d_list)^2))
        ratio <- if (is.finite(d_prev_norm) && d_prev_norm > 0) {
          d_norm / d_prev_norm
        } else {
          NA_real_
        }
        if (is.finite(ratio) && ratio > 0.05 && ratio < 0.98) {
          s_fac <- min(ratio / (1 - ratio), 20)
          lambda_e <- lambda_list
          phi_e <- phi_list
          m_e <- m_list
          ok <- TRUE
          for (g in seq_len(ngroups)) {
            m_g <- m_list[[g]] + s_fac * d_list[[g]]
            if (orthogonal) {
              # polar projection onto the orthogonal matrices
              svd_out <- try(svd(m_g), silent = TRUE)
              if (inherits(svd_out, "try-error")) {
                ok <- FALSE
                break
              }
              m_g <- svd_out$u %*% t(svd_out$v)
              phi_g <- phi_0[[g]]
            } else {
              # rescale the factors to unit variance
              m_inv <- try(solve(m_g), silent = TRUE)
              if (inherits(m_inv, "try-error")) {
                ok <- FALSE
                break
              }
              phi_g <- m_inv %*% phi_0[[g]] %*% t(m_inv)
              d_g <- 1 / sqrt(diag(phi_g))
              if (any(!is.finite(d_g))) {
                ok <- FALSE
                break
              }
              m_g <- t(t(m_g) / d_g)
              phi_g <- t(t(phi_g * d_g) * d_g)
            }
            m_e[[g]] <- m_g
            phi_e[[g]] <- phi_g
            lambda_e[[g]] <- lambda_0[[g]] %*% m_g
          }
          if (ok) {
            q_e <- objf(lambda_e, grad = TRUE)
            neval <- neval + 1L
            if (is.finite(q_e) && as.numeric(q_e) < q_current &&
                !is.null(attr(q_e, "grad"))) {
              lambda_list <- lambda_e
              phi_list <- phi_e
              m_list <- m_e
              q_current <- as.numeric(q_e)
              grad_list <- attr(q_e, "grad")
              n_extrap <- n_extrap + 1L
              if (lav_verbose()) {
                cat(
                  "  iter = ", sprintf("%4d", iter), " (extrapolated, s = ",
                  sprintf("%5.2f", s_fac), ") Q = ",
                  sprintf("%13.11f", q_current), "\n", sep = ""
                )
              }
              # restart the displacement history from the new point
              m_prev <- NULL
              d_prev_norm <- NA_real_
              next
            }
          }
        }
        d_prev_norm <- d_norm
      }
      m_prev <- m_list
    }
  } # iter

  list(
    lambda_list = lambda_list, phi_list = phi_list,
    q = q_current, iter = iter, converged = converged, neval = neval,
    n_extrap = n_extrap
  )
}
