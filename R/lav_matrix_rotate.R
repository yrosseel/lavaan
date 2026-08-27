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
            tol = tol,
            max_iter = max_iter
          )
        }
        info <- attr(rot, "info")
        attr(rot, "info") <- NULL
        c(info$method.value, lav_mat_vec(rot))
      }, error = function(e) c(Inf, rep(NA_real_, m * m)))

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

  # 4.b reorder the columns
  order_idx <- lav_efa_order_idx(mm_lambda, order_lv_by)

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

lav_mat_rotate_pairwise <- function(a = NULL, # original matrix
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
  # note: 'a' is the *current* loading matrix for this plane sweep (passed in by
  # optimize() below), not the original unrotated matrix of the same name in the
  # enclosing function.
  objf_orth <- function(theta = 0, a = NULL, col1 = 0L, col2 = 0L) {
    # construct ROT
    rot <- lav_mat_givens_orth(m, col1, col2, theta)

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
    rot <- lav_mat_givens_obliq(m, col1, col2, delta, phi12)

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
          a = mm_lambda, col1 = col1, col2 = col2,
          maximum = FALSE, tol = .Machine$double.eps^0.25
        )
        # best rotation - for this plane
        theta <- out$minimum

        # construct ROT
        rot <- lav_mat_givens_orth(m, col1, col2, theta)
      } else {
        phi12 <- phi[col1, col2]
        out <- optimize(
          f = objf_obliq, interval = c(-1, +1),
          a = mm_lambda, col1 = col1, col2 = col2,
          phi12 = phi12,
          maximum = FALSE, tol = .Machine$double.eps^0.25
        )

        # best rotation - for this plane
        delta <- out$minimum

        # construct ROT (rot[col1, col1] == gamma, the PHI scaling factor)
        rot <- lav_mat_givens_obliq(m, col1, col2, delta, phi12)
        gamma <- rot[col1, col1]
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
