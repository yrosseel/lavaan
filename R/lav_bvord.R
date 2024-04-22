# the weighted bivariate ordinal model
# YR 19 Feb 2020 (replacing the old lav_polychor.R routines)
#
# - polychoric (and tetrachoric) correlations
# - bivariate ordinal regression
# - using sampling weights wt

# two-way frequency table
# only works if Y = 1,2,3,...
lav_bvord_freq <- function(Y1, Y2, wt = NULL) {
  max.y1 <- max(Y1, na.rm = TRUE)
  max.y2 <- max(Y2, na.rm = TRUE)

  bin <- Y1 - 1L
  bin <- bin + max.y1 * (Y2 - 1L)
  bin <- bin + 1L

  if (is.null(wt)) {
    bin <- bin[!is.na(bin)]
    out <- array(tabulate(bin, nbins = max.y1 * max.y2),
      dim = c(max.y1, max.y2)
    )
  } else {
    if (anyNA(Y1) || anyNA(Y2)) {
      wt[is.na(Y1) | is.na(Y2)] <- 0
      bin[is.na(bin)] <- 0
    }
    y.ncat <- max.y1 * max.y2
    y.freq <- numeric(y.ncat)
    for (cat in seq_len(y.ncat)) {
      y.freq[cat] <- sum(wt[bin == cat])
    }
    out <- array(y.freq, dim = c(max.y1, max.y2))
  }

  out
}

# polychoric correlation
#
# zero.add is a vector: first element is for 2x2 tables only, second element
#                       for general tables
# zero.keep.margins is only used for 2x2 tables
#
lav_bvord_cor_twostep_fit <- function(Y1, Y2, eXo = NULL, wt = NULL,
                                      fit.y1 = NULL, fit.y2 = NULL,
                                      freq = NULL,
                                      zero.add = c(0.5, 0.0),
                                      zero.keep.margins = TRUE,
                                      zero.cell.warn = FALSE,
                                      zero.cell.flag = FALSE,
                                      verbose = FALSE,
                                      optim.method = "nlminb2",
                                      optim.scale = 1.0,
                                      init.theta = NULL,
                                      control = list(step.min = 0.1), # 0.6-7
                                      Y1.name = NULL, Y2.name = NULL) {
  if (is.null(fit.y1)) {
    fit.y1 <- lav_uvord_fit(y = Y1, X = eXo, wt = wt)
  }
  if (is.null(fit.y2)) {
    fit.y2 <- lav_uvord_fit(y = Y2, X = eXo, wt = wt)
  }

  # create cache environment
  cache <- lav_bvord_init_cache(fit.y1 = fit.y1, fit.y2 = fit.y2, wt = wt)

  # empty cells or not
  empty.cells <- FALSE

  # check for zero cells (if not exo), and catch some special cases
  if (cache$nexo == 0L) {
    freq <- cache$freq
    nr <- nrow(freq)
    nc <- ncol(freq)

    # check for empty cells
    if (any(freq == 0L)) {
      empty.cells <- TRUE
      if (zero.cell.warn) {
        if (!is.null(Y1.name) && !is.null(Y2.name)) {
          lav_msg_warn(gettextf(
            "empty cell(s) in bivariate table of %1$s x %2$s",
            Y1.name, Y2.name))
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
          if (zero.cell.flag) {
            attr(rho, "zero.cell.flag") <- empty.cells
          }
          return(rho)
        } else {
          rho <- -1.0
          if (zero.cell.flag) {
            attr(rho, "zero.cell.flag") <- empty.cells
          }
          return(rho)
        }
      } else if (length(idx) == 1L && zero.add[1] > 0.0) {
        if (zero.keep.margins) {
          # add + compensate to preserve margins
          if (idx == 1L || idx == 4L) { # main diagonal
            freq[1, 1] <- freq[1, 1] + zero.add[1]
            freq[2, 2] <- freq[2, 2] + zero.add[1]
            freq[2, 1] <- freq[2, 1] - zero.add[1]
            freq[1, 2] <- freq[1, 2] - zero.add[1]
          } else {
            freq[1, 1] <- freq[1, 1] - zero.add[1]
            freq[2, 2] <- freq[2, 2] - zero.add[1]
            freq[2, 1] <- freq[2, 1] + zero.add[1]
            freq[1, 2] <- freq[1, 2] + zero.add[1]
          }
        } else {
          freq[idx] <- freq[idx] + zero.add[1]
        }
      }
      # general table
    } else {
      if (any(freq == 0L) && zero.add[2] > 0.0) {
        # general table: just add zero.add to the empty cell(s)
        freq[freq == 0] <- zero.add[2]
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
        if (zero.cell.flag) {
          attr(rho, "zero.cell.flag") <- empty.cells
        }
        return(rho)
      }
      # 2. equal margins (th1 = th2 = 0)
      if (cache$th.y1[1] == 0 && cache$th.y2[1] == 0) {
        # see eg Brown & Benedetti 1977 eq 2
        rho <- -cos(2 * pi * freq[1, 1] / sum(freq))
        if (zero.cell.flag) {
          attr(rho, "zero.cell.flag") <- empty.cells
        }
        return(rho)
      }
    }
  } # non-exo

  # optim.method
  minObjective <- lav_bvord_min_objective
  minGradient <- lav_bvord_min_gradient
  minHessian <- lav_bvord_min_hessian
  if (optim.method == "nlminb" || optim.method == "nlminb2") {
    # nothing to do
  } else if (optim.method == "nlminb0") {
    minGradient <- minHessian <- NULL
  } else if (optim.method == "nlminb1") {
    minHessian <- NULL
  }

  # optimize
  if (is.null(control$trace)) {
    control$trace <- ifelse(verbose, 1, 0)
  }

  # init theta?
  if (!is.null(init.theta)) {
    start.x <- init.theta
  } else {
    start.x <- cache$theta
  }

  # try 1
  optim <- nlminb(
    start = start.x, objective = minObjective,
    gradient = minGradient, hessian = minHessian,
    control = control,
    scale = optim.scale, lower = -0.999, upper = +0.999,
    cache = cache
  )

  # try 2
  if (optim$convergence != 0L) {
    # try again, with different starting value
    optim <- nlminb(
      start = 0, objective = minObjective,
      gradient = NULL, hessian = NULL,
      control = control,
      scale = optim.scale, lower = -0.995, upper = +0.995,
      cache = cache
    )
  }

  # check convergence
  if (optim$convergence != 0L) {
    if (!is.null(Y1.name) && !is.null(Y2.name)) {
      lav_msg_warn(gettextf(
        "estimation polychoric correlation did not converge for
                    variables %1$s and %2$s", Y1.name, Y2.name))
    } else {
      lav_msg_warn(gettext(
        "estimation polychoric correlation(s) did not always converge"))
    }
    rho <- start.x
  } else {
    rho <- optim$par
  }

  # zero.cell.flag
  if (zero.cell.flag) {
    attr(rho, "zero.cell.flag") <- empty.cells
  }

  rho
}


# prepare cache environment
lav_bvord_init_cache <- function(fit.y1 = NULL,
                                 fit.y2 = NULL,
                                 wt = NULL,
                                 scores = FALSE,
                                 parent = parent.frame()) {
  # data
  Y1 <- fit.y1$y
  Y2 <- fit.y2$y
  eXo <- fit.y1$X

  # exo?
  if (is.null(eXo)) {
    nexo <- 0L
    freq <- lav_bvord_freq(Y1 = Y1, Y2 = Y2, wt = wt)
    th.y1 <- fit.y1$theta[fit.y1$th.idx]
    th.y2 <- fit.y2$theta[fit.y2$th.idx]
    nth.y1 <- length(th.y1)
    nth.y2 <- length(th.y2)
    pth.y1 <- pnorm(th.y1)
    pth.y2 <- pnorm(th.y2)
    upper.y <- rep(th.y2, times = rep.int(nth.y1, nth.y2))
    upper.x <- rep(th.y1, times = ceiling(length(upper.y)) / nth.y1)
  } else {
    nexo <- ncol(eXo)
    freq <- NULL
    fit.y1.z1 <- fit.y1$z1
    fit.y2.z1 <- fit.y2$z1
    fit.y1.z2 <- fit.y1$z2
    fit.y2.z2 <- fit.y2$z2

    # take care of missing values
    if (length(fit.y1$missing.idx) > 0L || length(fit.y2$missing.idx) > 0L) {
      missing.idx <- unique(c(fit.y1$missing.idx, fit.y2$missing.idx))
      fit.y1.z1[missing.idx] <- 0
      fit.y2.z1[missing.idx] <- 0
      fit.y1.z2[missing.idx] <- 0
      fit.y2.z2[missing.idx] <- 0
    } else {
      missing.idx <- integer(0L)
    }
  }

  # nobs
  if (is.null(wt)) {
    N <- length(Y1)
  } else {
    N <- sum(wt)
  }

  # starting value (for both exo and not-exo)
  # if(is.null(wt)) {
  rho.init <- cor(Y1, Y2, use = "pairwise.complete.obs")
  # }
  # cov.wt does not handle missing values...
  # rho.init <- cov.wt(cbind(Y1, Y2), wt = wt, cor = TRUE)$cor[2,1]
  if (is.na(rho.init) || abs(rho.init) >= 1.0) {
    rho.init <- 0.0
  }

  # parameter vector
  theta <- rho.init # only, for now

  # different cache if exo or not
  if (nexo == 0L) {
    if (scores) {
      out <- list2env(
        list(
          nexo = nexo, theta = theta, N = N,
          fit.y1.z1 = fit.y1$z1, fit.y1.z2 = fit.y1$z2,
          fit.y2.z1 = fit.y2$z1, fit.y2.z2 = fit.y2$z2,
          y1.Y1 = fit.y1$Y1, y1.Y2 = fit.y1$Y2,
          y2.Y1 = fit.y2$Y1, y2.Y2 = fit.y2$Y2,
          Y1 = Y1, Y2 = Y2, freq = freq,
          th.y1 = th.y1, th.y2 = th.y2,
          nth.y1 = nth.y1, nth.y2 = nth.y2,
          pth.y1 = pth.y1, pth.y2 = pth.y2,
          upper.y = upper.y, upper.x = upper.x
        ),
        parent = parent
      )
    } else {
      out <- list2env(
        list(
          nexo = nexo, theta = theta, N = N,
          Y1 = Y1, Y2 = Y2, freq = freq,
          th.y1 = th.y1, th.y2 = th.y2,
          nth.y1 = nth.y1, nth.y2 = nth.y2,
          pth.y1 = pth.y1, pth.y2 = pth.y2,
          upper.y = upper.y, upper.x = upper.x
        ),
        parent = parent
      )
    }
  } else {
    if (scores) {
      out <- list2env(
        list(
          nexo = nexo, theta = theta, wt = wt, N = N,
          eXo = eXo,
          y1.Y1 = fit.y1$Y1, y1.Y2 = fit.y1$Y2,
          y2.Y1 = fit.y2$Y1, y2.Y2 = fit.y2$Y2,
          fit.y1.z1 = fit.y1.z1, fit.y1.z2 = fit.y1.z2,
          fit.y2.z1 = fit.y2.z1, fit.y2.z2 = fit.y2.z2,
          missing.idx = missing.idx
        ),
        parent = parent
      )
    } else {
      out <- list2env(
        list(
          nexo = nexo, theta = theta, wt = wt, N = N,
          fit.y1.z1 = fit.y1.z1, fit.y1.z2 = fit.y1.z2,
          fit.y2.z1 = fit.y2.z1, fit.y2.z2 = fit.y2.z2,
          missing.idx = missing.idx
        ),
        parent = parent
      )
    }
  }

  out
}

# probabilities for each cell, given rho, th.y1 and th.y2
lav_bvord_noexo_pi_cache <- function(cache = NULL) {
  with(cache, {
    rho <- theta[1L]

    # catch special case: rho = 0.0
    if (rho == 0.0) {
      rowPI <- base::diff(c(0, pth.y1, 1))
      colPI <- base::diff(c(0, pth.y2, 1))
      PI.ij <- base::outer(rowPI, colPI)
      return(PI.ij)
    }

    BI <- pbivnorm::pbivnorm(x = upper.x, y = upper.y, rho = rho)
    dim(BI) <- c(nth.y1, nth.y2)
    BI <- rbind(0, BI, pth.y2, deparse.level = 0L)
    BI <- cbind(0, BI, c(0, pth.y1, 1), deparse.level = 0L)

    # get probabilities
    nr <- nrow(BI)
    nc <- ncol(BI)
    PI <- BI[-1L, -1L] - BI[-1L, -nc] - BI[-nr, -1L] + BI[-nr, -nc]

    # all elements should be strictly positive
    PI[PI < sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)

    return(PI)
  })
}

# partial derivative of CDF(th.y1, th.y2, rho) with respect to rho
lav_bvord_noexo_phi_cache <- function(cache = NULL) {
  with(cache, {
    rho <- theta[1L]

    # compute dbinorm for all possible combinations
    t1 <- rep(th.y1, times = nth.y2)
    t2 <- rep(th.y2, each = nth.y1)
    dbiNorm <- matrix(dbinorm(t1, t2, rho),
      nrow = nth.y1, ncol = nth.y2
    )

    p1 <- p2 <- p3 <- p4 <- matrix(0, nth.y1 + 1L, nth.y2 + 1L)
    t1.idx <- seq_len(nth.y1)
    t2.idx <- seq_len(nth.y2)

    # p1 is left-upper corner
    p1[t1.idx, t2.idx] <- dbiNorm
    # p2 is left-lower corner
    p2[t1.idx + 1L, t2.idx] <- dbiNorm
    # p3 is right-upper corner
    p3[t1.idx, t2.idx + 1L] <- dbiNorm
    # p3 is right-lower corner
    p4[t1.idx + 1L, t2.idx + 1L] <- dbiNorm

    phi <- p1 - p2 - p3 + p4
    return(phi)
  })
}

# Olsson 1979 A2
lav_bvord_noexo_gnorm_cache <- function(cache = NULL) {
  with(cache, {
    rho <- theta[1L]

    # note: Olsson 1979 A2 contains an error!!
    # derivative of phi_2(y1,y2;rho) wrt to rho equals
    # phi_2(y1,y2;rho) * guv(y1,y2;rho), where guv() is defined below:
    guv <- function(u, v, rho) {
      R <- (1 - rho * rho)
      (u * v * R - rho * ((u * u) - 2 * rho * u * v + (v * v)) + rho * R) / (R * R)
    }

    # compute gnorm for all possible combinations
    Gnorm <- dbiNorm * matrix(guv(t1, t2, rho), nth.y1, nth.y2)

    p1 <- p2 <- p3 <- p4 <- matrix(0, nth.y1 + 1L, nth.y2 + 1L)
    t1.idx <- seq_len(nth.y1)
    t2.idx <- seq_len(nth.y2)

    # p1 is left-upper corner
    p1[t1.idx, t2.idx] <- Gnorm
    # p2 is left-lower corner
    p2[t1.idx + 1L, t2.idx] <- Gnorm
    # p3 is right-upper corner
    p3[t1.idx, t2.idx + 1L] <- Gnorm
    # p3 is right-lower corner
    p4[t1.idx + 1L, t2.idx + 1L] <- Gnorm

    gnorm <- p1 - p2 - p3 + p4
    return(gnorm)
  })
}


# casewise likelihoods, unweighted!
lav_bvord_lik_cache <- function(cache = NULL) {
  with(cache, {
    rho <- theta[1L]

    # no exo
    if (nexo == 0L) {
      PI <- lav_bvord_noexo_pi_cache(cache)
      lik <- PI[cbind(Y1, Y2)]

      # exo
    } else {
      lik <- pbinorm(
        upper.x = fit.y1.z1, upper.y = fit.y2.z1,
        lower.x = fit.y1.z2, lower.y = fit.y2.z2, rho = rho
      )
      if (length(missing.idx) > 0L) {
        lik[missing.idx] <- NA
      }
      # catch very small values
      lik.toosmall.idx <- which(lik < sqrt(.Machine$double.eps))
      lik[lik.toosmall.idx] <- as.numeric(NA)
    }

    return(lik)
  })
}

lav_bvord_logl_cache <- function(cache = NULL) {
  with(cache, {
    rho <- theta[1L]

    # no exo
    if (nexo == 0L) {
      PI <- lav_bvord_noexo_pi_cache(cache)
      logl <- sum(freq * log(PI), na.rm = TRUE)

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
  })
}

lav_bvord_gradient_cache <- function(cache = NULL) {
  with(cache, {
    rho <- theta[1L]

    # no exo
    if (nexo == 0L) {
      phi <- lav_bvord_noexo_phi_cache(cache)
      bad.idx <- which(PI <= sqrt(.Machine$double.eps))
      if (length(bad.idx) > 0L) {
        PI[bad.idx] <- as.numeric(NA)
      }
      dx.rho <- sum((freq * phi) / PI, na.rm = TRUE)

      # exo
    } else {
      d1 <- dbinorm(fit.y1.z1, fit.y2.z1, rho)
      d2 <- dbinorm(fit.y1.z2, fit.y2.z1, rho)
      d3 <- dbinorm(fit.y1.z1, fit.y2.z2, rho)
      d4 <- dbinorm(fit.y1.z2, fit.y2.z2, rho)
      phi <- (d1 - d2 - d3 + d4)

      # avoid dividing by very tine numbers (new in 0.6-6)
      # -> done automatically: lik == NA in this case
      # bad.idx <- which(lik <= sqrt(.Machine$double.eps))
      # if(length(bad.idx) > 0L) {
      #    lik[bad.idx] <- as.numeric(NA)
      # }

      dx2 <- phi / lik

      if (is.null(wt)) {
        dx.rho <- sum(dx2, na.rm = TRUE)
      } else {
        dx.rho <- sum(wt * dx2, na.rm = TRUE)
      }
    }

    return(dx.rho)
  })
}

lav_bvord_hessian_cache <- function(cache = NULL) {
  with(cache, {
    rho <- theta[1L]

    # no exo
    if (nexo == 0L) {
      bad.idx <- which(PI <= sqrt(.Machine$double.eps))
      if (length(bad.idx) > 0L) {
        PI[bad.idx] <- as.numeric(NA)
      }
      gnorm <- lav_bvord_noexo_gnorm_cache(cache)
      # H <- sum( freq * (gnorm/PI - (phi*phi)/(PI*PI)), na.rm = TRUE)
      H <- (sum((freq * gnorm) / PI, na.rm = TRUE) -
        sum((freq * phi * phi) / (PI * PI), na.rm = TRUE))
      dim(H) <- c(1L, 1L)

      # exo
    } else {
      guv <- function(u, v, rho) {
        R <- (1 - rho * rho)
        (u * v * R - rho * ((u * u) - 2 * rho * u * v + (v * v)) + rho * R) / (R * R)
      }

      gnorm <- ((d1 * guv(fit.y1.z1, fit.y2.z1, rho)) -
        (d2 * guv(fit.y1.z2, fit.y2.z1, rho)) -
        (d3 * guv(fit.y1.z1, fit.y2.z2, rho)) +
        (d4 * guv(fit.y1.z2, fit.y2.z2, rho)))

      if (is.null(wt)) {
        H <- sum(gnorm / lik - (phi * phi) / (lik * lik), na.rm = TRUE)
      } else {
        H <- sum(wt * (gnorm / lik - (phi * phi) / (lik * lik)), na.rm = TRUE)
      }

      dim(H) <- c(1L, 1L)
    }

    return(H)
  })
}



# compute total (log)likelihood, for specific 'x' (nlminb)
lav_bvord_min_objective <- function(x, cache = NULL) {
  cache$theta <- x
  -1 * lav_bvord_logl_cache(cache = cache) / cache$N
}

# compute gradient, for specific 'x' (nlminb)
lav_bvord_min_gradient <- function(x, cache = NULL) {
  # check if x has changed
  if (!all(x == cache$theta)) {
    cache$theta <- x
    tmp <- lav_bvord_logl_cache(cache = cache)
  }
  -1 * lav_bvord_gradient_cache(cache = cache) / cache$N
}

# compute hessian, for specific 'x' (nlminb)
lav_bvord_min_hessian <- function(x, cache = NULL) {
  # check if x has changed
  if (!all(x == cache$theta)) {
    cache$theta <- x
    tmp <- lav_bvord_logl_cache(cache = cache)
    tmp <- lav_bvord_gradient_cache(cache = cache)
  }
  -1 * lav_bvord_hessian_cache(cache = cache) / cache$N
}




# casewise scores
lav_bvord_cor_scores_cache <- function(cache = NULL, na.zero = FALSE,
                                       use.weights = TRUE) {
  with(cache, {
    rho <- theta[1L]
    R <- sqrt(1 - rho * rho)

    # lik
    lik <- lav_bvord_lik_cache(cache = cache)
    bad.idx <- which(lik <= sqrt(.Machine$double.eps))
    if (length(bad.idx) > 0L) {
      lik[bad.idx] <- as.numeric(NA)
    }

    d.y1.z1 <- dnorm(fit.y1.z1)
    d.y1.z2 <- dnorm(fit.y1.z2)
    d.y2.z1 <- dnorm(fit.y2.z1)
    d.y2.z2 <- dnorm(fit.y2.z2)

    # th.y1
    if (identical(R, 0.0)) {
      y1.Z1 <- d.y1.z1 * 0.5
      y1.Z2 <- d.y1.z2 * 0.5
    } else {
      y1.Z1 <- (d.y1.z1 * pnorm((fit.y2.z1 - rho * fit.y1.z1) / R) -
        d.y1.z1 * pnorm((fit.y2.z2 - rho * fit.y1.z1) / R))
      y1.Z2 <- (d.y1.z2 * pnorm((fit.y2.z1 - rho * fit.y1.z2) / R) -
        d.y1.z2 * pnorm((fit.y2.z2 - rho * fit.y1.z2) / R))
    }
    dx.th.y1 <- (y1.Y1 * y1.Z1 - y1.Y2 * y1.Z2) / lik
    if (na.zero) {
      dx.th.y1[is.na(dx.th.y1)] <- 0
    }

    # th.y2
    if (identical(R, 0.0)) {
      y2.Z1 <- d.y2.z1 * 0.5
      y2.Z2 <- d.y2.z2 * 0.5
    } else {
      y2.Z1 <- (d.y2.z1 * pnorm((fit.y1.z1 - rho * fit.y2.z1) / R) -
        d.y2.z1 * pnorm((fit.y1.z2 - rho * fit.y2.z1) / R))
      y2.Z2 <- (d.y2.z2 * pnorm((fit.y1.z1 - rho * fit.y2.z2) / R) -
        d.y2.z2 * pnorm((fit.y1.z2 - rho * fit.y2.z2) / R))
    }
    dx.th.y2 <- (y2.Y1 * y2.Z1 - y2.Y2 * y2.Z2) / lik
    if (na.zero) {
      dx.th.y2[is.na(dx.th.y2)] <- 0
    }

    # slopes
    dx.sl.y1 <- dx.sl.y2 <- NULL
    if (nexo > 0L) {
      # sl.y1
      dx.sl.y1 <- (y1.Z2 - y1.Z1) * eXo / lik
      if (na.zero) {
        dx.sl.y1[is.na(dx.sl.y1)] <- 0
      }

      # sl.y2
      dx.sl.y2 <- (y2.Z2 - y2.Z1) * eXo / lik
      if (na.zero) {
        dx.sl.y2[is.na(dx.sl.y2)] <- 0
      }
    }

    # rho
    if (nexo == 0L) {
      phi <- lav_bvord_noexo_phi_cache(cache)
      dx <- phi[cbind(Y1, Y2)]
    } else {
      dx <- (dbinorm(fit.y1.z1, fit.y2.z1, rho) -
        dbinorm(fit.y1.z2, fit.y2.z1, rho) -
        dbinorm(fit.y1.z1, fit.y2.z2, rho) +
        dbinorm(fit.y1.z2, fit.y2.z2, rho))
    }
    dx.rho <- dx / lik
    if (na.zero) {
      dx.rho[is.na(dx.rho)] <- 0
    }

    if (!is.null(wt) && use.weights) {
      dx.th.y1 <- dx.th.y1 * wt
      dx.th.y2 <- dx.th.y2 * wt
      if (nexo > 0L) {
        dx.sl.y1 <- dx.sl.y1 * wt
        dx.sl.y2 <- dx.sl.y2 * wt
      }
      dx.rho <- dx.rho * wt
    }

    out <- list(
      dx.th.y1 = dx.th.y1, dx.th.y2 = dx.th.y2,
      dx.sl.y1 = dx.sl.y1, dx.sl.y2 = dx.sl.y2, dx.rho = dx.rho
    )
    return(out)
  })
}


# casewise scores - no cache
lav_bvord_cor_scores <- function(Y1, Y2, eXo = NULL, wt = NULL,
                                 rho = NULL,
                                 fit.y1 = NULL, fit.y2 = NULL,
                                 th.y1 = NULL, th.y2 = NULL,
                                 sl.y1 = NULL, sl.y2 = NULL,
                                 na.zero = FALSE, use.weights = TRUE) {
  if (is.null(fit.y1)) {
    fit.y1 <- lav_uvord_fit(y = Y1, X = eXo, wt = wt)
  }
  if (is.null(fit.y2)) {
    fit.y2 <- lav_uvord_fit(y = Y2, X = eXo, wt = wt)
  }

  # update z1/z2 if needed (used in pml_deriv1() in lav_model_gradient_pml.R)
  fit.y1 <- lav_uvord_update_fit(
    fit.y = fit.y1,
    th.new = th.y1, sl.new = sl.y1
  )
  fit.y2 <- lav_uvord_update_fit(
    fit.y = fit.y2,
    th.new = th.y2, sl.new = sl.y2
  )

  # create cache environment
  cache <- lav_bvord_init_cache(
    fit.y1 = fit.y1, fit.y2 = fit.y2, wt = wt,
    scores = TRUE
  )
  cache$theta <- rho

  SC <- lav_bvord_cor_scores_cache(
    cache = cache, na.zero = na.zero,
    use.weights = use.weights
  )

  SC
}

# logl - no cache
lav_bvord_logl <- function(Y1, Y2, eXo = NULL, wt = NULL,
                           rho = NULL,
                           fit.y1 = NULL, fit.y2 = NULL,
                           th.y1 = NULL, th.y2 = NULL,
                           sl.y1 = NULL, sl.y2 = NULL) {
  if (is.null(fit.y1)) {
    fit.y1 <- lav_uvord_fit(y = Y1, X = eXo, wt = wt)
  }
  if (is.null(fit.y2)) {
    fit.y2 <- lav_uvord_fit(y = Y2, X = eXo, wt = wt)
  }

  # update z1/z2 if needed (used in pml_deriv1() in lav_model_gradient_pml.R)
  fit.y1 <- lav_uvord_update_fit(
    fit.y = fit.y1,
    th.new = th.y1, sl.new = sl.y1
  )
  fit.y2 <- lav_uvord_update_fit(
    fit.y = fit.y2,
    th.new = th.y2, sl.new = sl.y2
  )

  # create cache environment
  cache <- lav_bvord_init_cache(fit.y1 = fit.y1, fit.y2 = fit.y2, wt = wt)
  cache$theta <- rho

  lav_bvord_logl_cache(cache = cache)
}

# lik - no cache
lav_bvord_lik <- function(Y1, Y2, eXo = NULL, wt = NULL,
                          rho = NULL,
                          fit.y1 = NULL, fit.y2 = NULL,
                          th.y1 = NULL, th.y2 = NULL,
                          sl.y1 = NULL, sl.y2 = NULL,
                          .log = FALSE) {
  if (is.null(fit.y1)) {
    fit.y1 <- lav_uvord_fit(y = Y1, X = eXo, wt = wt)
  }
  if (is.null(fit.y2)) {
    fit.y2 <- lav_uvord_fit(y = Y2, X = eXo, wt = wt)
  }

  # update fit.y1/fit.y2
  fit.y1 <- lav_uvord_update_fit(
    fit.y = fit.y1,
    th.new = th.y1, sl.new = sl.y1
  )
  fit.y2 <- lav_uvord_update_fit(
    fit.y = fit.y2,
    th.new = th.y2, sl.new = sl.y2
  )

  # create cache environment
  cache <- lav_bvord_init_cache(fit.y1 = fit.y1, fit.y2 = fit.y2, wt = wt)
  cache$theta <- rho

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

# noexo_pi - for backwards compatibility
lav_bvord_noexo_pi <- function(rho = NULL, th.y1 = NULL, th.y2 = NULL) {
  nth.y1 <- length(th.y1)
  nth.y2 <- length(th.y2)
  pth.y1 <- pnorm(th.y1)
  pth.y2 <- pnorm(th.y2)

  # catch special case: rho = 0.0
  if (rho == 0.0) {
    rowPI <- base::diff(c(0, pth.y1, 1))
    colPI <- base::diff(c(0, pth.y2, 1))
    PI.ij <- base::outer(rowPI, colPI)
    return(PI.ij)
  }

  # prepare for a single call to pbinorm
  upper.y <- rep(th.y2, times = rep.int(nth.y1, nth.y2))
  upper.x <- rep(th.y1, times = ceiling(length(upper.y)) / nth.y1)
  # rho <- rep(rho, length(upper.x)) # only one rho here

  BI <- pbivnorm::pbivnorm(x = upper.x, y = upper.y, rho = rho)
  dim(BI) <- c(nth.y1, nth.y2)
  BI <- rbind(0, BI, pth.y2, deparse.level = 0L)
  BI <- cbind(0, BI, c(0, pth.y1, 1), deparse.level = 0L)

  # get probabilities
  nr <- nrow(BI)
  nc <- ncol(BI)
  PI <- BI[-1L, -1L] - BI[-1L, -nc] - BI[-nr, -1L] + BI[-nr, -nc]

  # all elements should be strictly positive
  PI[PI < sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  PI
}
