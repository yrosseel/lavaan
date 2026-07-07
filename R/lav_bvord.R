# the weighted bivariate ordinal model
# YR 19 Feb 2020 (replacing the old lav_polychor.R routines)
#
# - polychoric (and tetrachoric) correlations
# - bivariate ordinal regression
# - using sampling weights wt
#
# YR/LDW 2026: refactored (see lav_uvbv_common.R for the shared cache/minfns
#              design); no change in behavior
#
#  info concerning lintr package, Luc DW May 7, 2026
#  there is a known problem in codetools, transferred to lintr, that
#  variables assigned in a with() are marked as 'may not be used';
#  to this end the code in with() statements are excluded from linting!
#

# two-way frequency table
# only works if Y = 1,2,3,...
lav_bvord_freq <- function(y1, y2, wt = NULL) {
  max_y1 <- max(y1, na.rm = TRUE)
  max_y2 <- max(y2, na.rm = TRUE)

  bin <- y1 - 1L
  bin <- bin + max_y1 * (y2 - 1L)
  bin <- bin + 1L

  if (is.null(wt)) {
    bin <- bin[!is.na(bin)]
    out <- array(tabulate(bin, nbins = max_y1 * max_y2),
      dim = c(max_y1, max_y2)
    )
  } else {
    if (anyNA(y1) || anyNA(y2)) {
      wt[is.na(y1) | is.na(y2)] <- 0
      bin[is.na(bin)] <- 0
    }
    y_ncat <- max_y1 * max_y2
    y_freq <- numeric(y_ncat)
    for (cat in seq_len(y_ncat)) {
      y_freq[cat] <- sum(wt[bin == cat])
    }
    out <- array(y_freq, dim = c(max_y1, max_y2))
  }

  out
}

# polychoric correlation
#
# zero.add is a vector: first element is for 2x2 tables only, second element
#                       for general tables
# zero.keep.margins is only used for 2x2 tables
#
lav_bvord_cor_twostep_fit <- function(y1, y2, exo = NULL, wt = NULL,
                                      fit_y1 = NULL, fit_y2 = NULL,
                                      freq = NULL,
                                      zero_add = c(0.5, 0.0),
                                      zero_keep_margins = TRUE,
                                      zero_cell_warn = FALSE,
                                      zero_cell_flag = FALSE,
                                      optim_method = "nlminb2",
                                      optim_scale = 1.0,
                                      init_theta = NULL,
                                      control = list(step.min = 0.1), # 0.6-7
                                      y1_name = NULL, y2_name = NULL) {
  # create cache environment
  cache <- lav_bvord_cache_from_args(
    y1 = y1, y2 = y2, exo = exo, wt = wt,
    fit_y1 = fit_y1, fit_y2 = fit_y2
  )

  # empty cells or not
  empty_cells <- FALSE

  # check for zero cells (if not exo), and catch some special cases
  if (cache$nexo == 0L) {
    freq <- cache$freq
    nr <- nrow(freq)
    nc <- ncol(freq)

    # check for empty cells
    if (any(freq == 0L)) {
      empty_cells <- TRUE
      if (zero_cell_warn) {
        if (!is.null(y1_name) && !is.null(y2_name)) {
          lav_msg_warn(gettextf(
            "empty cell(s) in bivariate table of %1$s x %2$s",
            y1_name, y2_name))
        } else {
          lav_msg_warn(gettext("empty cell(s) in bivariate table"))
        }
      }
    }

    # treat 2x2 tables
    if (nr == 2L && nc == 2L) {
      idx <- which(freq == 0L)
      # catch 2 empty cells: perfect correlation!
      if (length(idx) == 2L) {
        lav_msg_warn(gettext("two empty cells in 2x2 table"))
        if (freq[1, 1] > 0L) {
          rho <- 1.0
          if (zero_cell_flag) {
            attr(rho, "zero.cell.flag") <- empty_cells
          }
          return(rho)
        } else {
          rho <- -1.0
          if (zero_cell_flag) {
            attr(rho, "zero.cell.flag") <- empty_cells
          }
          return(rho)
        }
      } else if (length(idx) == 1L && zero_add[1] > 0.0) {
        if (zero_keep_margins) {
          # add + compensate to preserve margins
          if (idx == 1L || idx == 4L) { # main diagonal
            freq[1, 1] <- freq[1, 1] + zero_add[1]
            freq[2, 2] <- freq[2, 2] + zero_add[1]
            freq[2, 1] <- freq[2, 1] - zero_add[1]
            freq[1, 2] <- freq[1, 2] - zero_add[1]
          } else {
            freq[1, 1] <- freq[1, 1] - zero_add[1]
            freq[2, 2] <- freq[2, 2] - zero_add[1]
            freq[2, 1] <- freq[2, 1] + zero_add[1]
            freq[1, 2] <- freq[1, 2] + zero_add[1]
          }
        } else {
          freq[idx] <- freq[idx] + zero_add[1]
        }
      }
      # general table
    } else {
      if (any(freq == 0L) && zero_add[2] > 0.0) {
        # general table: just add zero.add to the empty cell(s)
        freq[freq == 0] <- zero_add[2]
      }
    }

    # update (possibly change) freq table
    cache$freq <- freq

    # catch special cases for 2x2 tables
    if (nr == 2L && nc == 2L) {
      # 1. a*d == c*d
      storage.mode(freq) <- "numeric" # to avoid integer overflow
      if (freq[1, 1] * freq[2, 2] == freq[1, 2] * freq[2, 1]) {
        rho <- 0.0
        if (zero_cell_flag) {
          attr(rho, "zero.cell.flag") <- empty_cells
        }
        return(rho)
      }
      # 2. equal margins (th1 = th2 = 0)
      if (cache$th_y1[1] == 0 && cache$th_y2[1] == 0) {
        # see eg Brown & Benedetti 1977 eq 2
        rho <- -cos(2 * pi * freq[1, 1] / sum(freq))
        if (zero_cell_flag) {
          attr(rho, "zero.cell.flag") <- empty_cells
        }
        return(rho)
      }
    }
  } # non-exo

  # optim.method
  minfns <- lav_uvbv_optim_fns(optim_method, lav_bvord_min_fns())

  # optimize
  if (is.null(control$trace)) {
    control$trace <- ifelse(lav_verbose(), 1, 0)
  }

  # init theta?
  if (!is.null(init_theta)) {
    start_x <- init_theta
  } else {
    start_x <- cache$theta
  }

  # try 1
  optim <- nlminb(
    start = start_x, objective = minfns$objective,
    gradient = minfns$gradient, hessian = minfns$hessian,
    control = control,
    scale = optim_scale, lower = -0.999, upper = +0.999,
    cache = cache
  )

  # try 2
  if (optim$convergence != 0L) {
    # try again, with different starting value
    optim <- nlminb(
      start = 0, objective = minfns$objective,
      gradient = NULL, hessian = NULL,
      control = control,
      scale = optim_scale, lower = -0.995, upper = +0.995,
      cache = cache
    )
  }

  # check convergence
  if (optim$convergence != 0L) {
    if (!is.null(y1_name) && !is.null(y2_name)) {
      lav_msg_warn(gettextf(
        "estimation polychoric correlation did not converge for
                    variables %1$s and %2$s", y1_name, y2_name))
    } else {
      lav_msg_warn(gettext(
        "estimation polychoric correlation(s) did not always converge"))
    }
    rho <- start_x
  } else {
    rho <- optim$par
  }

  # zero.cell.flag
  if (zero_cell_flag) {
    attr(rho, "zero.cell.flag") <- empty_cells
  }

  rho
}

# build the cache from user-level arguments: fit the univariate models if
# needed, apply parameter overrides, and initialize the cache environment
# (shared by lav_bvord_cor_twostep_fit/_cor_sc/_logl/_lik)
lav_bvord_cache_from_args <- function(y1, y2, exo = NULL, wt = NULL,
                                      rho = NULL,
                                      fit_y1 = NULL, fit_y2 = NULL,
                                      th_y1 = NULL, th_y2 = NULL,
                                      sl_y1 = NULL, sl_y2 = NULL,
                                      scores = FALSE) {
  if (is.null(fit_y1)) {
    fit_y1 <- lav_uvord_fit(y = y1, x = exo, wt = wt)
  }
  if (is.null(fit_y2)) {
    fit_y2 <- lav_uvord_fit(y = y2, x = exo, wt = wt)
  }

  # update z1/z2 if needed (used in lav_pml_dploglik_dimplied() in
  #                                     lav_model_gradient_pml.R)
  fit_y1 <- lav_uvord_update_fit(
    fit_y = fit_y1,
    th_new = th_y1, sl_new = sl_y1
  )
  fit_y2 <- lav_uvord_update_fit(
    fit_y = fit_y2,
    th_new = th_y2, sl_new = sl_y2
  )

  cache <- lav_bvord_init_cache(
    fit_y1 = fit_y1, fit_y2 = fit_y2, wt = wt,
    scores = scores
  )
  if (!is.null(rho)) {
    cache$theta <- rho
  }

  cache
}

# prepare cache environment
lav_bvord_init_cache <- function(fit_y1 = NULL,
                                 fit_y2 = NULL,
                                 wt = NULL,
                                 scores = FALSE,
                                 parent = parent.frame()) {
  # data
  y1 <- fit_y1$y
  y2 <- fit_y2$y
  exo <- fit_y1$x

  # exo?
  if (is.null(exo)) {
    nexo <- 0L
    freq <- lav_bvord_freq(y1 = y1, y2 = y2, wt = wt)
    th_y1 <- fit_y1$theta[fit_y1$th_idx]
    th_y2 <- fit_y2$theta[fit_y2$th_idx]
    nth_y1 <- length(th_y1)
    nth_y2 <- length(th_y2)
    pth_y1 <- pnorm(th_y1)
    pth_y2 <- pnorm(th_y2)
  } else {
    nexo <- ncol(exo)
    freq <- NULL
    fit_y1_z1 <- fit_y1$z1
    fit_y2_z1 <- fit_y2$z1
    fit_y1_z2 <- fit_y1$z2
    fit_y2_z2 <- fit_y2$z2

    # take care of missing values
    if (length(fit_y1$missing_idx) > 0L || length(fit_y2$missing_idx) > 0L) {
      missing_idx <- unique(c(fit_y1$missing_idx, fit_y2$missing_idx))
      fit_y1_z1[missing_idx] <- 0
      fit_y2_z1[missing_idx] <- 0
      fit_y1_z2[missing_idx] <- 0
      fit_y2_z2[missing_idx] <- 0
    } else {
      missing_idx <- integer(0L)
    }
  }

  # nobs
  if (is.null(wt)) {
    n <- length(y1)
  } else {
    n <- sum(wt)
  }

  # starting value (for both exo and not-exo)
  # if (is.null(wt)) {
  if (sd(y1, na.rm = TRUE) == 0 || sd(y2, na.rm = TRUE) == 0) {
    rho_init <- 0.0
  } else {
    rho_init <- cor(y1, y2, use = "pairwise.complete.obs")
  }
  # }
  # cov.wt does not handle missing values...
  # rho.init <- cov.wt(cbind(Y1, Y2), wt = wt, cor = TRUE)$cor[2,1]
  if (is.na(rho_init) || abs(rho_init) >= 1.0) {
    rho_init <- 0.0
  }

  # parameter vector
  theta <- rho_init # only, for now

  # cache environment: common part + exo/noexo specific parts
  # (+ the case x th indicator matrices when scores are needed)
  cache_list <- list(
    nexo = nexo, theta = theta, wt = wt, n = n,
    y1 = y1, y2 = y2, freq = freq
  )
  if (nexo == 0L) {
    cache_list <- c(cache_list, list(
      th_y1 = th_y1, th_y2 = th_y2,
      nth_y1 = nth_y1, nth_y2 = nth_y2,
      pth_y1 = pth_y1, pth_y2 = pth_y2
    ))
    if (scores) {
      cache_list <- c(cache_list, list(
        fit_y1_z1 = fit_y1$z1, fit_y1_z2 = fit_y1$z2,
        fit_y2_z1 = fit_y2$z1, fit_y2_z2 = fit_y2$z2
      ))
    }
  } else {
    cache_list <- c(cache_list, list(
      fit_y1_z1 = fit_y1_z1, fit_y1_z2 = fit_y1_z2,
      fit_y2_z1 = fit_y2_z1, fit_y2_z2 = fit_y2_z2,
      missing_idx = missing_idx
    ))
    if (scores) {
      cache_list <- c(cache_list, list(exo = exo))
    }
  }
  if (scores) {
    cache_list <- c(cache_list, list(
      y1_y1 = fit_y1$y1, y1_y2 = fit_y1$y2,
      y2_y1 = fit_y2$y1, y2_y2 = fit_y2$y2
    ))
  }

  list2env(cache_list, parent = parent)
}

# probabilities for each cell, given rho, th_y1 and th_y2
lav_bvord_noexo_pi_cache <- function(cache = NULL) {
  with(cache, {                     # nolint start
    pi0 <- lav_bvord_noexo_pi(
      rho = theta[1L], th_y1 = th_y1, th_y2 = th_y2
    )
    return(pi0)
  })                                         # nolint end
}

# the four-corner accumulation: given the (nth_y1 x nth_y2) matrix of values
# at all threshold crossings, compute the (nth_y1+1 x nth_y2+1) per-cell
# rectangle combination x[i,j] - x[i-1,j] - x[i,j-1] + x[i-1,j-1]
lav_bvord_corner_sum <- function(x = NULL, nth_y1 = NULL, nth_y2 = NULL) {
  p1 <- p2 <- p3 <- p4 <- matrix(0, nth_y1 + 1L, nth_y2 + 1L)
  t1_idx <- seq_len(nth_y1)
  t2_idx <- seq_len(nth_y2)

  # p1 is left-upper corner
  p1[t1_idx, t2_idx] <- x
  # p2 is left-lower corner
  p2[t1_idx + 1L, t2_idx] <- x
  # p3 is right-upper corner
  p3[t1_idx, t2_idx + 1L] <- x
  # p4 is right-lower corner
  p4[t1_idx + 1L, t2_idx + 1L] <- x

  p1 - p2 - p3 + p4
}

# partial derivative of CDF(th_y1, th_y2, rho) with respect to rho
lav_bvord_noexo_phi_cache <- function(cache = NULL) {
  with(cache, {                  #  nolint start
    rho <- theta[1L]

    # compute lav_dbinorm for all possible combinations
    t1 <- rep(th_y1, times = nth_y2)
    t2 <- rep(th_y2, each = nth_y1)
    dbi_norm <- matrix(lav_dbinorm(t1, t2, rho),
      nrow = nth_y1, ncol = nth_y2
    )

    phi <- lav_bvord_corner_sum(dbi_norm, nth_y1, nth_y2)
    return(phi)
  })                                          # nolint end
}

# derivative of phi_2(y1,y2;rho) wrt rho equals
# phi_2(y1,y2;rho) * guv(y1,y2;rho) -- see lav_dbinorm_guv() in lav_bvreg.R
lav_bvord_noexo_gnorm_cache <- function(cache = NULL) {
  with(cache, {                         # nolint start
    rho <- theta[1L]

    # compute gnorm for all possible combinations
    gnorm_1 <- dbi_norm * matrix(lav_dbinorm_guv(t1, t2, rho), nth_y1, nth_y2)

    gnorm <- lav_bvord_corner_sum(gnorm_1, nth_y1, nth_y2)
    return(gnorm)
  })                                          # nolint end
}


# casewise likelihoods, unweighted!
lav_bvord_lik_cache <- function(cache = NULL) {
  with(cache, {                     # nolint start
    rho <- theta[1L]

    # no exo
    if (nexo == 0L) {
      pi0 <- lav_bvord_noexo_pi_cache(cache)
      lik <- pi0[cbind(y1, y2)]

      # exo
    } else {
      lik <- pbinorm(
        upper_x = fit_y1_z1, upper_y = fit_y2_z1,
        lower_x = fit_y1_z2, lower_y = fit_y2_z2, rho = rho
      )
      if (length(missing_idx) > 0L) {
        lik[missing_idx] <- NA
      }
      # catch very small values
      lik_toosmall_idx <- which(lik < .Machine$double.eps ^ (2 / 3))
      lik[lik_toosmall_idx] <- as.numeric(NA)
    }

    return(lik)
  })                                     # nolint end
}

lav_bvord_logl_cache <- function(cache = NULL) {
  with(cache, {                  # nolint start
    rho <- theta[1L]

    # no exo
    if (nexo == 0L) {
      pi0 <- lav_bvord_noexo_pi_cache(cache)
      logl <- sum(freq * log(pi0), na.rm = TRUE)

      # exo
    } else {
      lik <- lav_bvord_lik_cache(cache) # unweighted!
      if (!is.null(wt)) {
        logl <- sum(wt * log(lik), na.rm = TRUE)
      } else {
        logl <- sum(log(lik), na.rm = TRUE)
      }
    }

    return(logl)
  })                                     # nolint end
}

lav_bvord_grad_cache <- function(cache = NULL) {
  with(cache, {              # nolint start
    rho <- theta[1L]

    # no exo
    if (nexo == 0L) {
      phi <- lav_bvord_noexo_phi_cache(cache)
      bad_idx <- which(pi0 < .Machine$double.eps^(2/3))
      if (length(bad_idx) > 0L) {
        pi0[bad_idx] <- as.numeric(NA)
      }
      dx_rho <- sum((freq * phi) / pi0, na.rm = TRUE)

      # exo
    } else {
      d1 <- lav_dbinorm(fit_y1_z1, fit_y2_z1, rho)
      d2 <- lav_dbinorm(fit_y1_z2, fit_y2_z1, rho)
      d3 <- lav_dbinorm(fit_y1_z1, fit_y2_z2, rho)
      d4 <- lav_dbinorm(fit_y1_z2, fit_y2_z2, rho)
      phi <- (d1 - d2 - d3 + d4)

      # avoid dividing by very tiny numbers (new in 0.6-6)
      # -> done automatically: lik == NA in this case

      dx2 <- phi / lik

      if (is.null(wt)) {
        dx_rho <- sum(dx2, na.rm = TRUE)
      } else {
        dx_rho <- sum(wt * dx2, na.rm = TRUE)
      }
    }

    return(dx_rho)
  })                                     # nolint end
}

lav_bvord_hessian_cache <- function(cache = NULL) {
  with(cache, {                    # nolint start
    rho <- theta[1L]

    # no exo
    if (nexo == 0L) {
      bad_idx <- which(pi0 < .Machine$double.eps^(2/3))
      if (length(bad_idx) > 0L) {
        pi0[bad_idx] <- as.numeric(NA)
      }
      gnorm <- lav_bvord_noexo_gnorm_cache(cache)
      # H <- sum( freq * (gnorm/PI - (phi*phi)/(PI*PI)), na.rm = TRUE)
      h <- (sum((freq * gnorm) / pi0, na.rm = TRUE) -
        sum((freq * phi * phi) / (pi0 * pi0), na.rm = TRUE))
      dim(h) <- c(1L, 1L)

      # exo
    } else {
      gnorm <- ((d1 * lav_dbinorm_guv(fit_y1_z1, fit_y2_z1, rho)) -
        (d2 * lav_dbinorm_guv(fit_y1_z2, fit_y2_z1, rho)) -
        (d3 * lav_dbinorm_guv(fit_y1_z1, fit_y2_z2, rho)) +
        (d4 * lav_dbinorm_guv(fit_y1_z2, fit_y2_z2, rho)))

      if (is.null(wt)) {
        h <- sum(gnorm / lik - (phi * phi) / (lik * lik), na.rm = TRUE)
      } else {
        h <- sum(wt * (gnorm / lik - (phi * phi) / (lik * lik)), na.rm = TRUE)
      }

      dim(h) <- c(1L, 1L)
    }

    return(h)
  })                                              # nolint end
}

# nlminb objective/gradient/hessian (see lav_uvbv_common.R)
lav_bvord_min_fns <- function() {
  lav_uvbv_min_fns(
    logl_fun = lav_bvord_logl_cache,
    grad_fun = lav_bvord_grad_cache,
    hessian_fun = lav_bvord_hessian_cache
  )
}


# casewise scores
lav_bvord_cor_sc_cache <- function(cache = NULL, na_zero = FALSE,
                                       use_weights = TRUE) {
  # make the arguments visible inside with(cache, ...): symbols there
  # resolve through the cache environment (and its parents), NOT through
  # this function's frame
  cache$na_zero <- na_zero
  cache$use_weights <- use_weights
  with(cache, {                 # nolint start
    rho <- theta[1L]
    r <- sqrt(1 - rho * rho)

    # lik
    lik <- lav_bvord_lik_cache(cache = cache)
    bad_idx <- which(lik < .Machine$double.eps^(2/3))
    if (length(bad_idx) > 0L) {
      lik[bad_idx] <- as.numeric(NA)
    }

    d_y1_z1 <- dnorm(fit_y1_z1)
    d_y1_z2 <- dnorm(fit_y1_z2)
    d_y2_z1 <- dnorm(fit_y2_z1)
    d_y2_z2 <- dnorm(fit_y2_z2)

    # th_y1
    if (identical(r, 0.0)) {
      y1_z1 <- d_y1_z1 * 0.5
      y1_z2 <- d_y1_z2 * 0.5
    } else {
      y1_z1 <- (d_y1_z1 * pnorm((fit_y2_z1 - rho * fit_y1_z1) / r) -
        d_y1_z1 * pnorm((fit_y2_z2 - rho * fit_y1_z1) / r))
      y1_z2 <- (d_y1_z2 * pnorm((fit_y2_z1 - rho * fit_y1_z2) / r) -
        d_y1_z2 * pnorm((fit_y2_z2 - rho * fit_y1_z2) / r))
    }
    dx_th_y1 <- (y1_y1 * y1_z1 - y1_y2 * y1_z2) / lik
    if (na_zero) {
      dx_th_y1[is.na(dx_th_y1)] <- 0
    }

    # th_y2
    if (identical(r, 0.0)) {
      y2_z1 <- d_y2_z1 * 0.5
      y2_z2 <- d_y2_z2 * 0.5
    } else {
      y2_z1 <- (d_y2_z1 * pnorm((fit_y1_z1 - rho * fit_y2_z1) / r) -
        d_y2_z1 * pnorm((fit_y1_z2 - rho * fit_y2_z1) / r))
      y2_z2 <- (d_y2_z2 * pnorm((fit_y1_z1 - rho * fit_y2_z2) / r) -
        d_y2_z2 * pnorm((fit_y1_z2 - rho * fit_y2_z2) / r))
    }
    dx_th_y2 <- (y2_y1 * y2_z1 - y2_y2 * y2_z2) / lik
    if (na_zero) {
      dx_th_y2[is.na(dx_th_y2)] <- 0
    }

    # slopes
    dx_sl_y1 <- dx_sl_y2 <- NULL
    if (nexo > 0L) {
      # sl.y1
      dx_sl_y1 <- (y1_z2 - y1_z1) * exo / lik
      if (na_zero) {
        dx_sl_y1[is.na(dx_sl_y1)] <- 0
      }

      # sl.y2
      dx_sl_y2 <- (y2_z2 - y2_z1) * exo / lik
      if (na_zero) {
        dx_sl_y2[is.na(dx_sl_y2)] <- 0
      }
    }

    # rho
    if (nexo == 0L) {
      phi <- lav_bvord_noexo_phi_cache(cache)
      dx <- phi[cbind(y1, y2)]
    } else {
      dx <- (lav_dbinorm(fit_y1_z1, fit_y2_z1, rho) -
        lav_dbinorm(fit_y1_z2, fit_y2_z1, rho) -
        lav_dbinorm(fit_y1_z1, fit_y2_z2, rho) +
        lav_dbinorm(fit_y1_z2, fit_y2_z2, rho))
    }
    dx_rho <- dx / lik
    if (na_zero) {
      dx_rho[is.na(dx_rho)] <- 0
    }

    if (!is.null(wt) && use_weights) {
      dx_th_y1 <- dx_th_y1 * wt
      dx_th_y2 <- dx_th_y2 * wt
      if (nexo > 0L) {
        dx_sl_y1 <- dx_sl_y1 * wt
        dx_sl_y2 <- dx_sl_y2 * wt
      }
      dx_rho <- dx_rho * wt
    }

    out <- list(
      dx_th_y1 = dx_th_y1, dx_th_y2 = dx_th_y2,
      dx_sl_y1 = dx_sl_y1, dx_sl_y2 = dx_sl_y2, dx_rho = dx_rho
    )
    return(out)
  })                                                   # nolint end
}


# casewise scores - no cache
lav_bvord_cor_sc <- function(y1, y2, exo = NULL, wt = NULL,
                                 rho = NULL,
                                 fit_y1 = NULL, fit_y2 = NULL,
                                 th_y1 = NULL, th_y2 = NULL,
                                 sl_y1 = NULL, sl_y2 = NULL,
                                 na_zero = FALSE, use_weights = TRUE) {
  cache <- lav_bvord_cache_from_args(
    y1 = y1, y2 = y2, exo = exo, wt = wt, rho = rho,
    fit_y1 = fit_y1, fit_y2 = fit_y2,
    th_y1 = th_y1, th_y2 = th_y2,
    sl_y1 = sl_y1, sl_y2 = sl_y2,
    scores = TRUE
  )

  lav_bvord_cor_sc_cache(
    cache = cache, na_zero = na_zero,
    use_weights = use_weights
  )
}

# logl - no cache
lav_bvord_logl <- function(y1, y2, exo = NULL, wt = NULL,
                           rho = NULL,
                           fit_y1 = NULL, fit_y2 = NULL,
                           th_y1 = NULL, th_y2 = NULL,
                           sl_y1 = NULL, sl_y2 = NULL) {
  cache <- lav_bvord_cache_from_args(
    y1 = y1, y2 = y2, exo = exo, wt = wt, rho = rho,
    fit_y1 = fit_y1, fit_y2 = fit_y2,
    th_y1 = th_y1, th_y2 = th_y2,
    sl_y1 = sl_y1, sl_y2 = sl_y2
  )

  lav_bvord_logl_cache(cache = cache)
}

# lik - no cache
lav_bvord_lik <- function(y1, y2, exo = NULL, wt = NULL,
                          rho = NULL,
                          fit_y1 = NULL, fit_y2 = NULL,
                          th_y1 = NULL, th_y2 = NULL,
                          sl_y1 = NULL, sl_y2 = NULL,
                          .log = FALSE) {
  cache <- lav_bvord_cache_from_args(
    y1 = y1, y2 = y2, exo = exo, wt = wt, rho = rho,
    fit_y1 = fit_y1, fit_y2 = fit_y2,
    th_y1 = th_y1, th_y2 = th_y2,
    sl_y1 = sl_y1, sl_y2 = sl_y2
  )

  lik <- lav_bvord_lik_cache(cache = cache) # unweighted
  if (.log) {
    lik <- log(lik)
  }

  if (!is.null(wt)) {
    if (.log) {
      lik <- wt * lik
    } else {
      tmp <- wt * log(lik)
      lik <- exp(tmp)
    }
  }

  lik
}

# probabilities for each cell, given rho, th_y1 and th_y2
# (also used directly by lav_tables.R and the PML objective)
lav_bvord_noexo_pi <- function(rho = NULL, th_y1 = NULL, th_y2 = NULL) {
  nth_y1 <- length(th_y1)
  nth_y2 <- length(th_y2)
  pth_y1 <- pnorm(th_y1)
  pth_y2 <- pnorm(th_y2)

  # catch special case: rho = 0.0
  if (rho == 0.0) {
    row_pi <- base::diff(c(0, pth_y1, 1))
    col_pi <- base::diff(c(0, pth_y2, 1))
    pi_ij <- base::outer(row_pi, col_pi)
    return(pi_ij)
  }

  # prepare for a single call to pbivnorm
  upper_y <- rep(th_y2, times = rep.int(nth_y1, nth_y2))
  upper_x <- rep(th_y1, times = ceiling(length(upper_y)) / nth_y1)

  bi <- pbivnorm::pbivnorm(x = upper_x, y = upper_y, rho = rho)
  dim(bi) <- c(nth_y1, nth_y2)
  bi <- rbind(0, bi, pth_y2, deparse.level = 0L)
  bi <- cbind(0, bi, c(0, pth_y1, 1), deparse.level = 0L)

  # get probabilities
  nr <- nrow(bi)
  nc <- ncol(bi)
  pi0 <- bi[-1L, -1L] - bi[-1L, -nc] - bi[-nr, -1L] + bi[-nr, -nc]

  # all elements should be strictly positive
  pi0[pi0 < .Machine$double.eps ^ (2 / 3)] <-
                   .Machine$double.eps ^ (2 / 3)
  pi0
}
