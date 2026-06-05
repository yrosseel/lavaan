# utility functions
#
# initial version: YR 25/03/2009

# multivariate normal random number generation
# replacement for MASS::mvrnorm for better cross-machine reproducibility
# see: https://blog.djnavarro.net/posts/2025-05-18_multivariate-normal-sampling-floating-point/  # nolint
#
# the issue with MASS::mvrnorm is that it uses eigendecomposition with a
# transformation matrix (sqrt(Lambda) Q') that is NOT invariant to eigenvector
# sign flips. this leads to irreproducible results across machines even with
# the same random seed. this implementation follows mvtnorm::rmvnorm which uses:
# - eigen: Q sqrt(Lambda) Q' which IS invariant to sign flips (default)
# - svd: similar to eigen, also invariant
# - chol: unique factorization
# eigen is the default to match MASS::mvrnorm's approach while fixing the
# sign-flip issue with the invariant formula
lav_mvrnorm <- function(n = 1, mu, sigma_1, tol = 1e-06, empirical = FALSE,
                      method = "eigen", check_symmetry = TRUE, byrow = FALSE) {
  p <- length(mu)

  # check symmetry (like mvtnorm)
  if (check_symmetry && !isSymmetric(sigma_1,
    tol = sqrt(.Machine$double.eps),
    check.attributes = FALSE
  )) {
    lav_msg_stop(gettext("'Sigma' must be a symmetric matrix in lav_mvrnorm"))
  }

  # check dimensions (like mvtnorm)
  if (p != nrow(sigma_1)) {
    lav_msg_stop(gettext("incompatible arguments in lav_mvrnorm"))
  }

  # compute transformation matrix R based on method (following mvtnorm exactly)
  method <- match.arg(method, c("chol", "eigen", "svd"))

  r <- if (method == "eigen") {
    ev <- eigen(sigma_1, symmetric = TRUE)
    if (!all(ev$values >= -tol * abs(ev$values[1L]))) {
      warning("sigma is numerically not positive semidefinite")
    }
    # Q sqrt(Lambda) Q' - invariant to sign flips
    t(ev$vectors %*% (t(ev$vectors) * sqrt(pmax(ev$values, 0))))
  } else if (method == "svd") {
    s <- svd(sigma_1)
    if (!all(s$d >= -tol * abs(s$d[1L]))) {
      warning("sigma is numerically not positive semidefinite")
    }
    t(s$v %*% (t(s$u) * sqrt(pmax(s$d, 0))))
  } else if (method == "chol") {
    r <- chol(sigma_1, pivot = TRUE)
    r[, order(attr(r, "pivot"))]
  }

  # for names (fallback to dimnames of Sigma if mu has no names)
  nm <- names(mu)
  if (is.null(nm) && !is.null(dn <- dimnames(sigma_1))) {
    nm <- dn[[1L]]
  }

  if (empirical) {
    # generate standard normal, then apply empirical transformation
    x <- matrix(stats::rnorm(p * n), n, p, byrow = byrow)
    x <- scale(x, center = TRUE, scale = FALSE) # center
    x <- x %*% svd(x, nu = 0)$v # orthogonalize
    x <- scale(x, center = FALSE, scale = TRUE) # unit variance

    # transform by R
    x <- sweep(x %*% r, 2, mu, "+")
    colnames(x) <- nm

    if (n == 1) drop(x) else x
  } else {
    # generate samples (following mvtnorm exactly)
    x <- matrix(stats::rnorm(n * p), nrow = n, byrow = byrow) %*% r
    x <- sweep(x, 2, mu, "+")
    colnames(x) <- nm

    if (n == 1) drop(x) else x
  }
}

# sd with trimming
lav_sample_trimmed_sd <- function(x, na_rm = TRUE, trim = 0) {
  if (isTRUE(na_rm)) {
    x <- x[!is.na(x)]
  }
  n <- length(x)
  if (trim > 0 && n) {
    if (is.complex(x)) {
      lav_msg_stop(gettext("trimmed means are not defined for complex data"))
    }
    if (anyNA(x)) {
      return(NA_real_)
    }
    if (trim >= 0.5) {
      return(stats::median(x, na.rm = FALSE))
    }
    lo <- floor(n * trim) + 1
    hi <- n + 1 - lo
    x <- sort.int(x, partial = unique(c(lo, hi)))[lo:hi]
  }
  sd(x)
}

# convert correlation matrix + standard deviations to covariance matrix
# based on cov2cor in package:stats
lav_cor2cov <- function(R, sds, names = NULL) {             # nolint
  p <- (d <- dim(R))[1L]
  if (!is.numeric(R) || length(d) != 2L || p != d[2L]) {
    lav_msg_stop(gettext("'V' is not a square numeric matrix"))
  }

  if (any(!is.finite(sds))) {
    lav_msg_warn(gettext(
      "sds had 0 or NA entries; non-finite result is doubtful"
    ))
  }

  # if(sum(diag(R)) != p)
  #    stop("The diagonal of a correlation matrix should be all ones.")

  if (p != length(sds)) {
    lav_msg_stop(gettext("The standard deviation vector and correlation matrix
                         have a different number of variables"))
  }

  s <- R
  s[] <- sds * R * rep(sds, each = p)

  # optionally, add names
  if (!is.null(names)) {
    stopifnot(length(names) == p)
    rownames(s) <- colnames(s) <- names
  }

  s
}

# convert characters within single quotes to numeric vector
# eg. s <- '3 4.3 8e-3 2.0'
#     x <- lav_char2num(s)
lav_char2num <- function(s = "") {
  # first, strip all ',' or ';'
  s_1 <- gsub(",", " ", s)
  s_1 <- gsub(";", " ", s_1)
  tc <- textConnection(s_1)
  x <- scan(tc, quiet = TRUE)
  close(tc)
  x
}

# create full matrix based on lower.tri or upper.tri elements; add names
# always ROW-WISE!!
lav_getcov <- function(x, lower = TRUE, diagonal = TRUE, sds = NULL,
                       names = paste("V", 1:nvar, sep = "")) {
  # check x and sds
  if (is.character(x)) x <- lav_char2num(x)
  if (is.character(sds)) sds <- lav_char2num(sds)

  if (lower) {
    cov_1 <- lav_mat_lower2full(x, diagonal = diagonal)
  } else {
    cov_1 <- lav_mat_upper2full(x, diagonal = diagonal)
  }
  nvar <- ncol(cov_1)

  # if diagonal is false, assume unit diagonal
  if (!diagonal) diag(cov_1) <- 1

  # check if we have a sds argument
  if (!is.null(sds)) {
    stopifnot(length(sds) == nvar)
    cov_1 <- lav_cor2cov(cov_1, sds)
  }

  # names
  stopifnot(length(names) == nvar)
  rownames(cov_1) <- colnames(cov_1) <- names

  cov_1
}

lav_char2hash <- function(s = "") {
  stopifnot(is.character(s))
  nums <- utf8ToInt(paste(s, collapse = "\n"))
  rval <- 0x7EDCBA98L
  for (i in nums) {
    positions <- 1L + (bitwAnd(i, 15L))
    bits_r <- bitwShiftL(1L, positions) - 1L
    wrap1 <- bitwShiftL(bitwAnd(bits_r, rval), 31L - positions)
    wrap2 <- bitwShiftR(bitwAnd(bitwNot(bits_r), rval), positions)
    saw <- bitwOr(wrap1, wrap2)
    rval <- bitwXor(saw, i)
  }
  as.hexmode(rval)
}

# vectorize all (h0 or h1) sample statistics, in the same order
# as Gamma
lav_implied_to_vec <- function(implied = NULL, lavmodel = NULL,
                               drop_list = TRUE) {
  ngroups <- lavmodel@ngroups

  wls_obs <- vector("list", ngroups)
  for (g in seq_len(ngroups)) {

    var <- NULL
    res_var <- NULL
    if (lavmodel@conditional.x) {
      res_var <- diag(implied$res.cov[[g]])
      res_var <- res_var[res_var != 1]
    } else {
      var <- diag(implied$cov[[g]])
      var <- var[var != 1]
    }

    wls_obs[[g]] <- lav_samp_wls_obs(
      # plain
      mean_g = implied$mean[[g]],
      cov_g = implied$cov[[g]],
      var_g = var,
      th_g = implied$th[[g]],
      th_idx_g = lavmodel@th.idx[[g]],

      # conditional.x
      res_int_g = implied$res.int[[g]],
      res_cov_g = implied$res.cov[[g]],
      res_var_g = res_var,
      res_th_g = implied$res.th[[g]],
      res_slopes_g = implied$res.slopes[[g]],
      group_w_g = implied$group.w[[g]],

      # flags
      categorical = lavmodel@categorical,
      conditional_x = lavmodel@conditional.x,
      meanstructure = lavmodel@meanstructure,
      correlation = lavmodel@correlation,
      slopestructure = lavmodel@conditional.x,
      group_w_free = lavmodel@group.w.free
    )
  }

  if (drop_list) {
    out <- unlist(wls_obs)
  } else {
    out <- wls_obs
  }

  out
}

# given a vector of sample statistics, reconstruct implied
# entries (cov, mean, th, res.cov, ...)
lav_vec_to_implied <- function(x = NULL, lavmodel) {

  ngroups <- lavmodel@ngroups
  implied <- list()

  for (g in seq_len(ngroups)) {

    # number of variables
    nvar <- lavmodel@nvar[g]

    # if group.w.free, always comes first
    if (lavmodel@group.w.free) {
      idx <- 1L
      group_w <- x[idx]
      x <- x[-idx]
    } else {
      group_w <- 1
    }

    if (lavmodel@categorical) {
      if (lavmodel@conditional.x) {
        # TODO
        cat("\n NOT READY YET! \n")
      }

      # th
      nth <- length(lavmodel@th.idx[[g]])
      idx <- seq_len(nth)
      th_g <- x[idx]
      x <- x[-idx]

      # mean/var ov.num
      mean_g <- rep(0, nvar)
      var_g <- rep(1, nvar)
      if (any(lavmodel@th.idx[[g]] == 0)) {
        num_idx <- seq_len(nvar)
        num_idx <- num_idx[!num_idx %in% lavmodel@th.idx[[g]]]
        idx <- seq_along(num_idx)
        var_g[num_idx] <- x[idx]
        x <- x[-idx]
        # FIXME: change sign?
        mean_g[num_idx] <- th_g[lavmodel@th.idx[[g]] == 0]
      }

      # cov - lower only
      idx <- seq_len(nvar * (nvar - 1) / 2)
      cov_g <- lav_mat_vech_rev(x[idx], diagonal = FALSE)
      x <- x[-idx]
      diag(cov_g) <- var_g

      # fill in
      implied$cov[[g]] <- cov_g
      implied$mean[[g]] <- mean_g
      implied$th[[g]] <- th_g

    } else {
      # diag flag
      diag_flag <- TRUE
      if (lavmodel@correlation) {
        diag_flag <- FALSE
      }
      if (lavmodel@conditional.x) {
        # TODO
        cat("\n NOT READY YET! \n")
      }
      if (lavmodel@meanstructure) {
        idx <- seq_len(nvar)
        mean_g <- x[idx]
        x <- x[-idx]
      } else {
        mean_g <- numeric(nvar)
      }
      if (diag_flag) {
        idx <- seq_len(nvar * (nvar + 1) / 2)
      } else {
        idx <- seq_len(nvar * (nvar - 1) / 2)
      }
      cov_g <- lav_mat_vech_rev(x[idx], diagonal = diag_flag)
      x <- x[-idx]

      # fill in
      implied$cov[[g]] <- cov_g
      implied$mean[[g]] <- mean_g
      implied$th[g] <- list(NULL)
    }

    implied$group.w[[g]] <- group_w
  }

  implied
}


# find the solution of the linear equation:
#
# theta = solve(t(Delta) %*% W %*% Delta) %*% t(Delta) %*% W %*% svec
#
# where Delta is the Jacobian of Sigma wrt the free parameters,
# W is lav_mat_bdiag(S.inv, S22_inv)
# where S22_inv = 0.5 * lav_mat_dup_pre_post(S.inv %x% S.inv)
#  - S can be I (ULS) (then W = 0.5 * t(D) %*% D)
#  - S can be sample.cov (GLS)
#  - S can be Sigma (2RLS or RLS)
#
# svec is the vector with sample statistics (first means, then vech(S))
#
# strategy:
# - first we compute Q = W %*% Delta, using the fact that W only contains S.inv
#   elements (without explicitly constructing the kronecker product)
#   note that for the cov block:
#     W %*% Delta = w * vech(S.inv * unvech(v) * S.inv)
#   where w is 1/2 for diagonal elements, and 1 elsewhere
# - compute A = t(Delta) %*% Q, and b = Q %*% svec
# - compute solve(A,b) using cholesky back/forward, as A is symmetric and pd
#
# optionally, return H = solve(t(Delta) %*% W %*% Delta) %*% t(Delta) %*% W
# as an attribute
lav_utils_wls_linearization <- function(delta = NULL, s = NULL,
                                        meanstructure = FALSE,
                                        categorical = FALSE,
                                        #fixed.x = FALSE, # FIXME: needed?
                                        #x.idx = integer(0L),
                                        svec = NULL,
                                        return_h = FALSE) {
  nvar <- nrow(s)
  if (categorical) {
    meanstructure <- FALSE
    # all ordinal (for now)
    pstar <- nvar * (nvar - 1L) / 2L
  } else {
    pstar <- nvar * (nvar + 1L) / 2L
  }

  jac_cov <- delta
  nc <- ncol(jac_cov)
  s_cov <- svec

  # vech
  w <- rep(1.0, pstar)

  # split Delta/svec
  if (categorical) {
    cov_idx  <- seq(nrow(delta) - pstar + 1L, nrow(delta))
    jac_cov  <- delta[cov_idx, , drop = FALSE]
    s_cov    <- svec[cov_idx]
  } else if (!categorical && meanstructure) {
    mean_idx <- seq_len(nvar)
    jac_mean <- delta[mean_idx, , drop = FALSE]
    jac_cov <- delta[-mean_idx, , drop = FALSE]
    s_mean <- svec[mean_idx]
    s_cov <- svec[-mean_idx]
    w[lav_mat_diagh_idx(nvar)] <- 0.5
  } else if (nrow(delta) != pstar) {
    lav_msg_stop(gettext("nrow(Delta) != pstar"))
  }

  # special case: S = I (ULS)
  uls_flag <- isTRUE(all.equal(s, diag(nvar), check.attributes = FALSE))
  if (uls_flag || categorical) {
    # cov part
    q_cov <- w * jac_cov
    a <- crossprod(jac_cov, q_cov)
    b <- drop(crossprod(q_cov, s_cov))

    # mean part
    if (!categorical && meanstructure) {
      a <- a + crossprod(jac_mean)
      b <- b + drop(crossprod(jac_mean, s_mean))
    }

    # tQ only needed if H is requested
    if (return_h) {
      t_q <- t(q_cov)
      if (meanstructure) t_q <- cbind(t(jac_mean), t_q)
    }
  } else {
    # S inverse
    c_s <- chol(s)
    si <- chol2inv(c_s)

    # cov part
    q_cov <- matrix(0, pstar, nc)
    for (j in seq_len(nc)) {
      vmat <- lav_mat_vech_rev(jac_cov[, j])
      m_w <- si %*% vmat %*% si
      q_cov[, j] <- w * lav_mat_vech(m_w)
    }

    a <- crossprod(jac_cov, q_cov)
    b <- drop(crossprod(q_cov, s_cov))

    # mean part
    if (meanstructure) {
      q_mean <- backsolve(c_s, forwardsolve(t(c_s), jac_mean))
      a <- a + crossprod(jac_mean, q_mean)
      b <- b + drop(crossprod(q_mean, s_mean))
    }

    if (return_h) {
      t_q <- t(q_cov)
      if (meanstructure) t_q <- cbind(t(q_mean), t_q)
    }
  }

  r <- chol(a)
  if (!return_h) {
    out <- backsolve(r, forwardsolve(t(r), b))
  } else {
    ainv  <- chol2inv(r)
    out <- drop(ainv %*% b)
    h <- ainv %*% t_q
    attr(out, "H") <- h
  }

  out
}

# function to transform names of variables to snake_case
# this function is used mainly to rename function arguments given in a list
#     where 'old' names are still accepted to avoid breaking other packages
lav_snake_case <- function(old_names) {
  curval <- c("B", "C", "D", "E", "K", "W", "PI",
             "TAU", "DELTA", "NU", "LAMBDA", "eXo",
              "WMAT", "THETA", "ALPHA", "BETA", "GAMMA", "PSI",
              "SminTheta")
  newval <- c("m_b", "m_c", "m_d", "m_e", "m_k", "m_w", "pi0",
             "mm_tau", "mm_delta", "mm_nu", "mm_lambda", "exo",
             "mm_wmat", "mm_theta", "mm_alpha", "mm_beta", "mm_gamma", "mm_psi",
             "s_min_theta")
  # transform dot.case and CamelCase to snake_case
  varnames_new <- tolower(chartr(".", "_",
                   gsub("([a-z])([A-Z])", "\\1_\\2", old_names)))
  # apply standard modifications
  mtch <- match(old_names, curval)
  for (j in seq_along(old_names)) {
    if (!is.na(mtch[j])) varnames_new[j] <- newval[mtch[j]]
  }
  # remove trailing underscores in new names
  varnames_new <- gsub("_$", "", varnames_new)
  # check no doubles in new names
  tocheck <- varnames_new[varnames_new != ""]
  doubles <- anyDuplicated(tocheck)
  if (doubles) {
    lav_msg_stop(gettextf("At least one snake_cased name (%s) is duplicated!",
                   tocheck[doubles]))
  }
  varnames_new
}

# function to put arguments with old names (in ...) in the new named argument
# the function to adapt must have a ... argument and call this function in the
# beginning as follows :
#   dotdotdot <- list(...)
#   lav_adapt_func(environment(), dotdotdot, TRUE/FALSE)
# The argument include_ddd if set to TRUE transforms the names of dotdotdot also.
lav_adapt_func <- function(envir, dotdotdot, include_ddd = TRUE) {
  lijst <- ls(pos=envir)
  if (length(dotdotdot) > 0L) {
    dddnames <- names(dotdotdot)
    dddnewnames <- lav_snake_case(dddnames)
    newargs <- dddnewnames %in% lijst
    for (j in which(newargs)) {
      assign(dddnewnames[j], dotdotdot[[dddnames[j]]], envir)
      dotdotdot[[dddnames[j]]] <- NULL
    }
    if (include_ddd && length(dotdotdot) > 0L) {
      names(dotdotdot) <- lav_snake_case(names(dotdotdot))
    }
    assign("dotdotdot", dotdotdot, envir)
  }
}
