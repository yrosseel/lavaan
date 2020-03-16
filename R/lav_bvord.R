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

    if(is.null(wt)) {
        bin <- bin[!is.na(bin)]
        out <- array(tabulate(bin, nbins = max.y1 * max.y2),
                     dim = c(max.y1, max.y2))
    } else {
        if(anyNA(Y1) || anyNA(Y2)) {
            wt[is.na(Y1) | is.na(Y2)] <- 0
            bin[is.na(bin)] <- 0
        }
        y.ncat <- max.y1 * max.y2
        y.freq <- numeric(y.ncat)
        for(cat in seq_len(y.ncat)) {
            y.freq[cat] <- sum(wt[bin == cat])
        }
        out <- array(y.freq, dim = c(max.y1, max.y2))
    }

    out
}

# probabilities for each cell
# given rho, th.y1 and th.y2
lav_bvord_noexo_pi <- function(rho, th.y1, th.y2) {

    nth.y1 <- length(th.y1); nth.y2 <- length(th.y2)
    pth.y1 <- pnorm(th.y1);  pth.y2 <- pnorm(th.y2)

    # catch special case: rho = 0.0
    if(rho == 0.0) {
        rowPI <- diff(c(0,pth.y1,1))
        colPI <- diff(c(0,pth.y2,1))
        PI.ij <- outer(rowPI, colPI)
        return(PI.ij)
    }

    # prepare for a single call to pbinorm
    upper.y <- rep(th.y2, times = rep.int(nth.y1, nth.y2))
    upper.x <- rep(th.y1, times = ceiling(length(upper.y))/nth.y1)
    #rho <- rep(rho, length(upper.x)) # only one rho here

    BI <- pbivnorm::pbivnorm(x = upper.x, y = upper.y, rho = rho)
    dim(BI) <- c(nth.y1, nth.y2)
    BI <- rbind(0, BI, pth.y2, deparse.level = 0)
    BI <- cbind(0, BI, c(0, pth.y1, 1), deparse.level = 0)

    # get probabilities
    nr <- nrow(BI); nc <- ncol(BI)
    PI <- BI[-1L,-1L] - BI[-1L,-nc] - BI[-nr,-1L] + BI[-nr,-nc]

    # all elements should be strictly positive
    PI[PI < .Machine$double.eps] <- .Machine$double.eps

    PI
}

# partial derivative of CDF(th.y1, th.y2, rho) with respect
# to rho
lav_bvord_noexo_phi <- function(rho, th.y1, th.y2) {

    n.th1 <- length(th.y1); t1.idx <- seq_len(n.th1)
    n.th2 <- length(th.y2); t2.idx <- seq_len(n.th2)

    # compute dbinorm for all possible combinations
    dbiNorm <- matrix(dbinorm(rep(th.y1, times = n.th2),
                              rep(th.y2, each = n.th1), rho),
                      nrow = n.th1, ncol = n.th2)

    # p1 is left-upper corner
    p1 <- p2 <- p3 <- p4 <- matrix(0, n.th1 + 1L, n.th2 + 1L)

    # p1 is left-upper corner
    p1[t1.idx     , t2.idx     ] <- dbiNorm
    # p2 is left-lower corner
    p2[t1.idx + 1L, t2.idx     ] <- dbiNorm
    # p3 is right-upper corner
    p3[t1.idx     , t2.idx + 1L] <- dbiNorm
    # p3 is right-lower corner
    p4[t1.idx + 1L, t2.idx + 1L] <- dbiNorm

    p1 - p2 - p3 + p4
}

# Olsson 1979 A2
lav_bvord_noexo_gnorm <- function(rho, th.y1, th.y2) {

    # note: Olsson 1979 A2 contains an error!!
    guv <- function(u, v, rho) {
        R <- (1 - rho*rho)
        ( u*v*R - rho*((u*u) - 2*rho*u*v + (v*v)) + rho*R ) / (R*R)
    }

    n.th1 <- length(th.y1); t1.idx <- seq_len(n.th1)
    n.th2 <- length(th.y2); t2.idx <- seq_len(n.th2)

    # compute gnorm for all possible combinations
    t1 <- rep(th.y1, times = n.th2)
    t2 <- rep(th.y2, each = n.th1)
    Gnorm <- matrix( dbinorm(t1, t2, rho) * guv(t1, t2, rho),
                 nrow = n.th1, ncol = n.th2)

    # p1 is left-upper corner
    p1 <- p2 <- p3 <- p4 <- matrix(0, n.th1 + 1L, n.th2 + 1L)

    # p1 is left-upper corner
    p1[t1.idx     , t2.idx     ] <- Gnorm
    # p2 is left-lower corner
    p2[t1.idx + 1L, t2.idx     ] <- Gnorm
    # p3 is right-upper corner
    p3[t1.idx     , t2.idx + 1L] <- Gnorm
    # p3 is right-lower corner
    p4[t1.idx + 1L, t2.idx + 1L] <- Gnorm

    p1 - p2 - p3 + p4
}

# (summed) loglikelihood
lav_bvord_logl <- function(Y1, Y2, wt = NULL, eXo = NULL,
                           rho = NULL,
                           fit.y1 = NULL, fit.y2 = NULL,
                           freq = NULL) {
    stopifnot(!is.null(rho))

    if(is.null(fit.y1)) {
        fit.y1 <- lav_uvord_fit(y = Y1, X = eXo, wt = wt)
    }
    if(is.null(fit.y2)) {
        fit.y2 <- lav_uvord_fit(y = Y2, X = eXo, wt = wt)
    }

    # if no eXo, use shortcut (grouped)
    if(length(fit.y1$slope.idx) == 0L) {
        if(is.null(freq)) {
            freq <- lav_bvord_freq(Y1 = fit.y1$y, Y2 = fit.y2$y, wt = wt)
        }

        # grouped lik
        PI <- lav_bvord_noexo_pi(rho, th.y1 = fit.y1$theta[fit.y1$th.idx],
                                      th.y2 = fit.y2$theta[fit.y2$th.idx])
        # PIs are always strictly positive, no check needed
        logl <- sum( freq * log(PI) )

    } else {
        loglik <- lav_bvord_lik(Y1 = Y1, Y2 = Y2, rho = rho, wt = wt,
                                fit.y1 = fit.y1, fit.y2 = fit.y2, .log = TRUE)
        logl <- sum(loglik, na.rm = TRUE)
        # logl will be -Inf if any element in '(log)lik' is -Inf
    }

    logl
}

# casewise likelihoods
# (loglikelihoods if .log = TRUE)
lav_bvord_lik <- function(Y1, Y2, wt = NULL, eXo = NULL, rho = NULL,
                          fit.y1 = NULL, fit.y2 = NULL, .log = FALSE) {

    stopifnot(!is.null(rho))

    if(is.null(fit.y1)) {
        fit.y1 <- lav_uvord_fit(y = Y1, X = eXo, wt = wt)
    }
    if(is.null(fit.y2)) {
        fit.y2 <- lav_uvord_fit(y = Y2, X = eXo, wt = wt)
    }

    # if no eXo, use shortcut (grouped)
    if(length(fit.y1$slope.idx) == 0L) {
        # probability per cell
        PI <- lav_bvord_noexo_pi(rho, th.y1 = fit.y1$theta[fit.y1$th.idx],
                                      th.y2 = fit.y2$theta[fit.y2$th.idx])
        if(.log) {
            # PI elements are always strictly positive
            logPI <- log(PI)
            lik <- logPI[ cbind(fit.y1$y, fit.y2$y) ]
        } else {
            lik <- PI[ cbind(fit.y1$y, fit.y2$y) ]
        }
    } else {
        # take care of missing values
        if(length(fit.y1$missing.idx) > 0L ||
           length(fit.y2$missing.idx) > 0L) {
            missing.idx <- unique(c(fit.y1$missing.idx,
                                    fit.y2$missing.idx))
            fit.y1.z1 <- fit.y1$z1; fit.y1.z1[missing.idx] <- 0
            fit.y2.z1 <- fit.y2$z1; fit.y2.z1[missing.idx] <- 0
            fit.y1.z2 <- fit.y1$z2; fit.y1.z2[missing.idx] <- 0
            fit.y2.z2 <- fit.y2$z2; fit.y2.z2[missing.idx] <- 0

            lik <- pbinorm(upper.x = fit.y1.z1, upper.y = fit.y2.z1,
                           lower.x = fit.y1.z2, lower.y = fit.y2.z2, rho = rho)
            lik[missing.idx] <- NA
        } else {
            lik <- pbinorm(upper.x = fit.y1$z1, upper.y = fit.y2$z1,
                           lower.x = fit.y1$z2, lower.y = fit.y2$z2, rho = rho)
        }

        if(.log) {
            # lik[lik < .Machine$double.eps] <- .Machine$double.eps
            lik[lik < .Machine$double.eps] <- -Inf
            lik <- log(lik)
        }
    }

    # take care of observation weights wt
    if(!is.null(wt)) {
        if(.log) {
            lik <- wt * lik
        } else {
            tmp <- wt * log(lik)
            lik <- exp(tmp)
        }
    }

    lik
}


# polychoric correlation
#
# zero.add is a vector: first element is for 2x2 tables only, second element
#                       for general tables
# zero.keep.margins is only used for 2x2 tables
lav_bvord_cor_twostep <- function(Y1, Y2, eXo = NULL, wt = NULL,
                                  fit.y1 = NULL, fit.y2 = NULL, freq = NULL,
                                  method = "nlminb",  control = list(),
                                  zero.add = c(0.5, 0.0),
                                  zero.keep.margins = TRUE,
                                  zero.cell.warn = FALSE,
                                  zero.cell.flag = FALSE,
                                  verbose = FALSE,
                                  Y1.name = NULL, Y2.name = NULL) {

    if(is.null(fit.y1)) {
        fit.y1 <- lav_uvord_fit(y = Y1, X = eXo, wt = wt)
    } else {
        Y1 <- fit.y1$y
        eXo <- fit.y1$X
    }
    if(is.null(fit.y2)) {
        fit.y2 <- lav_uvord_fit(y = Y2, X = eXo, wt = wt)
    } else {
        Y2 <- fit.y2$y
    }

    stopifnot(method %in% c("nlminb", "BFGS", "optimize"))

    # empty cells or not
    empty.cells <- FALSE

    # exo or not?
    exo <- ifelse(length(fit.y1$slope.idx) > 0L, TRUE, FALSE)

    # thresholds
    th.y1 <- fit.y1$theta[fit.y1$th.idx]
    th.y2 <- fit.y2$theta[fit.y2$th.idx]

    # check for zero cells (if not exo), and catch some special cases
    if(!exo) {
        if(is.null(freq)) {
            freq <- lav_bvord_freq(Y1 = fit.y1$y, Y2 = fit.y2$y, wt = wt)
        }
        nr <- nrow(freq); nc <- ncol(freq)

        # check for empty cells
        if(any(freq == 0L)) {
            empty.cells <- TRUE
            if(zero.cell.warn) {
                if(!is.null(Y1.name) && !is.null(Y2.name)) {
                    warning("lavaan WARNING: ",
                            "empty cell(s) in bivariate table of ",
                            Y1.name, " x ", Y2.name)
                } else {
                    warning("lavaan WARNING: empty cell(s) in bivariate table")
                }
            }
        }

        # treat 2x2 tables
        if(nr == 2L && nc == 2L) {
            idx <- which(freq == 0L)
            # catch 2 empty cells: perfect correlation!
            if(length(idx) == 2L) {
                warning("lavaan WARNING: two empty cells in 2x2 table")
                if(freq[1,1] > 0L) {
                    rho <- 1.0
                    if(zero.cell.flag) {
                        attr(rho, "zero.cell.flag") <- empty.cells
                    }
                    return(rho)
                } else {
                    rho <- -1.0
                    if(zero.cell.flag) {
                        attr(rho, "zero.cell.flag") <- empty.cells
                    }
                    return(rho)
                }
            } else if(length(idx) == 1L && zero.add[1] > 0.0) {
                if(zero.keep.margins) {
                    # add + compensate to preserve margins
                    if(idx == 1L || idx == 4L) { # main diagonal
                        freq[1,1] <- freq[1,1] + zero.add[1]
                        freq[2,2] <- freq[2,2] + zero.add[1]
                        freq[2,1] <- freq[2,1] - zero.add[1]
                        freq[1,2] <- freq[1,2] - zero.add[1]
                    } else {
                        freq[1,1] <- freq[1,1] - zero.add[1]
                        freq[2,2] <- freq[2,2] - zero.add[1]
                        freq[2,1] <- freq[2,1] + zero.add[1]
                        freq[1,2] <- freq[1,2] + zero.add[1]
                    }
                } else {
                    freq[idx] <- freq[idx] + zero.add[1]
                }
            }
        # general table
        } else {
            if(any(freq == 0L) && zero.add[2] > 0.0) {
                # general table: just add zero.add to the empty cell(s)
                freq[freq == 0] <- zero.add[2]
            }
        }

        # catch special cases for 2x2 tables
        if(nr == 2L && nc == 2L) {
            # 1. a*d == c*d
            storage.mode(freq) <- "numeric" # to avoid integer overflow
            if(freq[1,1]*freq[2,2] == freq[1,2]*freq[2,1]) {
                rho <- 0.0
                if(zero.cell.flag) {
                    attr(rho, "zero.cell.flag") <- empty.cells
                }
                return(rho)
            }
            # 2. equal margins (th1 = th2 = 0)
            if(th.y1[1] == 0 && th.y2[1] == 0) {
                # see eg Brown & Benedetti 1977 eq 2
                rho <- - cos( 2*pi*freq[1,1]/sum(freq) )
                if(zero.cell.flag) {
                    attr(rho, "zero.cell.flag") <- empty.cells
                }
                return(rho)
            }
        }
    } # !exo

    objectiveFunction <- function(x) {
        logl <- lav_bvord_logl(rho = tanh(x[1L]),
                               fit.y1 = fit.y1, fit.y2 = fit.y2, wt = wt,
                               freq = freq)
        -logl # to minimize!
    }

    gradientFunction <- function(x) {
        rho <- tanh(x[1L])
        if(!exo) {
            PI  <- lav_bvord_noexo_pi( rho, th.y1, th.y2)
            phi <- lav_bvord_noexo_phi(rho, th.y1, th.y2)
            dx.rho <- sum(freq/PI * phi)
        } else {
            lik <- lav_bvord_lik(rho = rho,
                                 fit.y1 = fit.y1, fit.y2 = fit.y2, wt = NULL,
                                 .log = FALSE)
            dx <- ( dbinorm(fit.y1$z1, fit.y2$z1, rho) -
                    dbinorm(fit.y1$z2, fit.y2$z1, rho) -
                    dbinorm(fit.y1$z1, fit.y2$z2, rho) +
                    dbinorm(fit.y1$z2, fit.y2$z2, rho) ) / lik
            if(is.null(wt)) {
                dx.rho <- sum(dx, na.rm = TRUE)
            } else {
                dx.rho <- sum(wt * dx, na.rm = TRUE)
            }
        }
        # dF/drho * drho/dx, dtanh = 1/cosh(x)^2
        -dx.rho * 1/(cosh(x)*cosh(x))
    }

    # starting value
    rho.init <- cor(Y1, Y2, use = "pairwise.complete.obs")

    # if rho.init is missing, set starting value to zero
    if(is.na(rho.init)) {
      rho.init <- 0.0
    }
    # check range of rho.init is within [-1,+1]
    if(abs(rho.init) >= 1.0) {
        rho.init <- 0.0
    }

    # default values
    control.nlminb <- list(eval.max = 20000L,
                           iter.max = 10000L,
                           trace = ifelse(verbose, 1L, 0L),
                           #abs.tol = 1e-20, ### important!! fx never negative
                           abs.tol = (.Machine$double.eps * 10),
                           rel.tol = ifelse(method == "nlminb", 1e-10, 1e-7),
                           x.tol = 1.5e-8,
                           xf.tol = 2.2e-14)
    control.nlminb <- modifyList(control.nlminb, control)
    control <- control.nlminb[c("eval.max", "iter.max", "trace",
                                "abs.tol", "rel.tol", "x.tol", "xf.tol")]

    if(method == "nlminb") {
        out <- nlminb(start = atanh(rho.init), objective = objectiveFunction,
                      gradient = gradientFunction,
                      scale = 10,
                      control = control)
        if(out$convergence != 0L) {
            warning("no convergence")
        }
        rho <- tanh(out$par)
    } else if(method == "BFGS") {
        # NOTE: known to fail if rho.init is too far from final value
        # seems to be better with parscale = 0.1??
        out <- optim(par = atanh(rho.init), fn = objectiveFunction,
                     gr = gradientFunction,
                     control = list(parscale = 0.1, reltol = 1e-10,
                                    trace = ifelse(verbose, 1L, 0L),
                                    REPORT = 1L,
                                    abstol=(.Machine$double.eps * 10)),
                     method = "BFGS")
        if(out$convergence != 0L) {
            warning("no convergence")
        }
        rho <- tanh(out$par)
    } else if(method == "optimize") {
        # not atanh/tanh transform
        objectiveFunction2 <- function(x) {
            logl <- lav_bvord_logl(rho = x[1L],
                                   fit.y1 = fit.y1, fit.y2 = fit.y2,
                                   wt = wt, freq=freq)
            -logl # to minimize!
        }
        out <- optimize(f = objectiveFunction2, interval = c(-0.9999,0.9999))
        rho <- out$minimum
    }

    if(zero.cell.flag) {
        attr(rho, "zero.cell.flag") <- empty.cells
    }

    rho
}

# casewise scores
lav_bvord_cor_scores <- function(Y1, Y2, eXo = NULL, wt = NULL,
                                 rho = NULL,
                                 fit.y1 = NULL, fit.y2 = NULL,
                                 th.y1 = NULL, th.y2 = NULL,
                                 sl.y1 = NULL, sl.y2 = NULL,
                                 na.zero = FALSE) {
    R <- sqrt(1 - rho*rho)
    if(is.null(fit.y1)) {
        fit.y1 <- lav_uvord_fit(y = Y1, X = eXo, wt = wt)
    }
    if(is.null(fit.y2)) {
        fit.y2 <- lav_uvord_fit(y = Y2, X = eXo, wt = wt)
    }

    # update z1/z2 if needed (used in pml_deriv1() in lav_model_gradient_pml.R)
    fit.y1 <- lav_uvord_update_fit(fit.y = fit.y1,
                                   th.new = th.y1, sl.new = sl.y1)
    fit.y2 <- lav_uvord_update_fit(fit.y = fit.y2,
                                   th.new = th.y2, sl.new = sl.y2)

    # lik
    lik <- lav_bvord_lik(rho = rho, fit.y1 = fit.y1, fit.y2 = fit.y2, wt = NULL)

    # th.y1
    if(identical(R, 0.0)) {
        y1.Z1 <- dnorm(fit.y1$z1) * 0.5
        y1.Z2 <- dnorm(fit.y1$z2) * 0.5
    } else {
        y1.Z1 <- ( dnorm(fit.y1$z1) * pnorm( (fit.y2$z1-rho*fit.y1$z1) / R) -
                   dnorm(fit.y1$z1) * pnorm( (fit.y2$z2-rho*fit.y1$z1) / R) )
        y1.Z2 <- ( dnorm(fit.y1$z2) * pnorm( (fit.y2$z1-rho*fit.y1$z2) / R) -
                   dnorm(fit.y1$z2) * pnorm( (fit.y2$z2-rho*fit.y1$z2) / R) )
    }
    dx.th.y1 <- (fit.y1$Y1*y1.Z1 - fit.y1$Y2*y1.Z2) / lik
    if(!is.null(wt)) {
        dx.th.y1 <- dx.th.y1 * wt
    }
    if(na.zero) {
        dx.th.y1[is.na(dx.th.y1)] <- 0
    }

    # th.y2
    if(identical(R, 0.0)) {
        y2.Z1  <- dnorm(fit.y2$z1) * 0.5
        y2.Z2  <- dnorm(fit.y2$z2) * 0.5
    } else {
        y2.Z1  <- ( dnorm(fit.y2$z1) * pnorm( (fit.y1$z1-rho*fit.y2$z1) / R) -
                    dnorm(fit.y2$z1) * pnorm( (fit.y1$z2-rho*fit.y2$z1) / R) )
        y2.Z2  <- ( dnorm(fit.y2$z2) * pnorm( (fit.y1$z1-rho*fit.y2$z2) / R) -
                    dnorm(fit.y2$z2) * pnorm( (fit.y1$z2-rho*fit.y2$z2) / R) )
    }
    dx.th.y2 <- (fit.y2$Y1*y2.Z1 - fit.y2$Y2*y2.Z2) / lik
    if(!is.null(wt)) {
        dx.th.y2 <- dx.th.y2 * wt
    }
    if(na.zero) {
        dx.th.y2[is.na(dx.th.y2)] <- 0
    }

    # slopes
    dx.sl.y1 <- dx.sl.y2 <- NULL
    if(length(fit.y1$slope.idx) > 0L) {

        # sl.y1
        dx.sl.y1 <- (y1.Z2 - y1.Z1) * fit.y1$X / lik
        if(!is.null(wt)) {
            dx.sl.y1 <- dx.sl.y1 * wt
        }
        if(na.zero) {
            dx.sl.y1[is.na(dx.sl.y1)] <- 0
        }

        # sl.y2
        dx.sl.y2 <- (y2.Z2 - y2.Z1) * fit.y1$X / lik
        if(!is.null(wt)) {
            dx.sl.y2 <- dx.sl.y2 * wt
        }
        if(na.zero) {
            dx.sl.y2[is.na(dx.sl.y2)] <- 0
        }
    }

    # rho
    if(length(fit.y1$slope.idx) == 0L) {
        phi <- lav_bvord_noexo_phi(rho, th.y1 = fit.y1$theta[fit.y1$th.idx],
                                        th.y2 = fit.y2$theta[fit.y2$th.idx])
        dx <- phi[cbind(fit.y1$y, fit.y2$y)]
    } else {
        dx <- ( dbinorm(fit.y1$z1, fit.y2$z1, rho) -
                dbinorm(fit.y1$z2, fit.y2$z1, rho) -
                dbinorm(fit.y1$z1, fit.y2$z2, rho) +
                dbinorm(fit.y1$z2, fit.y2$z2, rho) )
    }
    dx.rho <- dx / lik
    if(!is.null(wt)) {
        dx.rho <- dx.rho * wt
    }
    if(na.zero) {
        dx.rho[is.na(dx.rho)] <- 0
    }

    list(dx.th.y1 = dx.th.y1, dx.th.y2 = dx.th.y2,
         dx.sl.y1 = dx.sl.y1, dx.sl.y2 = dx.sl.y2, dx.rho = dx.rho)
}



