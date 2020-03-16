# the weighted bivariate ordinal/linear model
# YR 08 March 2020 (replacing the old lav_polyserial.R routines)
#
# - polyserial (and biserial) correlations
# - bivariate ordinal/linear regression
# - using sampling weights wt

# (summed) loglikelihood
lav_bvmix_logl <- function(Y1 = NULL, Y2 = NULL, wt = NULL, eXo = NULL,
                           rho = NULL,
                           fit.y1 = NULL, fit.y2 = NULL,
                           freq = NULL) {

    loglik <- lav_bvmix_lik(Y1 = Y1, Y2 = Y2, eXo = eXo, rho = rho, wt = wt,
                            fit.y1 = fit.y1, fit.y2 = fit.y2, .log = TRUE)
    logl <- sum(loglik, na.rm = TRUE)
    # logl will be -Inf if any element in '(log)lik' is -Inf

    logl
}

# casewise likelihoods
# (loglikelihoods if .log = TRUE)
# Y1 = linear
# Y2 = ordinal
lav_bvmix_lik <- function(Y1 = NULL, Y2 = NULL, wt = NULL, eXo = NULL,
                          rho = NULL,
                          evar.y1 = NULL, beta.y1 = NULL,
                          th.y2 = NULL, sl.y2 = NULL,
                          fit.y1 = NULL, fit.y2 = NULL, .log = FALSE) {

    stopifnot(!is.null(rho))

    if(is.null(fit.y1)) {
        fit.y1 <- lav_uvreg_fit(y = Y1, X = eXo, wt = wt)
    } else {
        Y1 <- fit.y1$y
    }
    if(is.null(fit.y2)) {
        fit.y2 <- lav_uvord_fit(y = Y2, X = eXo, wt = wt)
    }

    # user specified parameters
    if(!is.null(evar.y1) || !is.null(beta.y1)) {
        fit.y1 <- lav_uvreg_update_fit(fit.y = fit.y1,
                                       evar.new = evar.y1, beta.new = beta.y1)
    }
    if(!is.null(th.y2) || !is.null(sl.y2)) {
        fit.y2 <- lav_uvord_update_fit(fit.y = fit.y2,
                                       th.new = th.y2, sl.new = sl.y2)
    }

    R <- sqrt(1 - rho*rho)
    y1.SD  <- sqrt(fit.y1$theta[fit.y1$var.idx])
    y1.ETA <- fit.y1$yhat
    Z <- (Y1 - y1.ETA) / y1.SD

    # p(Y2|Y1)
    tauj.star  <- (fit.y2$z1 - rho*Z)/R
    tauj1.star <- (fit.y2$z2 - rho*Z)/R
    py2y1 <- pnorm(tauj.star) - pnorm(tauj1.star)
    py2y1[py2y1 < .Machine$double.eps] <- .Machine$double.eps

    # p(Y1)
    py1 <- dnorm(Y1, mean = y1.ETA, sd = y1.SD)

    # lik
    lik <- py1 * py2y1

    # log?
    if(.log) {
        # lik[lik < .Machine$double.eps] <- .Machine$double.eps
        lik[lik < .Machine$double.eps] <- -Inf
        lik <- log(lik)
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


# polyserial correlation
#
# Y1 = linear
# Y2 = ordinal
lav_bvmix_cor_twostep <- function(Y1, Y2, eXo = NULL, wt = NULL,
                                  fit.y1 = NULL, fit.y2 = NULL,
                                  method = "nlminb",  control = list(),
                                  verbose = FALSE) {

    if(is.null(fit.y1)) {
        fit.y1 <- lav_uvreg_fit(y = Y1, X = eXo, wt = wt)
    } else {
        Y1 <- fit.y1$y
        eXo <- fit.y1$X
    }
    if(is.null(fit.y2)) {
        fit.y2 <- lav_uvord_fit(y = Y2, X = eXo, wt = wt)
    } else {
        Y2 <- fit.y2$y
    }

    y1.VAR <- fit.y1$theta[fit.y1$var.idx]; y1.SD  <- sqrt(y1.VAR)
    y1.ETA <- fit.y1$yhat
    Z <- (Y1 - y1.ETA) / y1.SD


    objectiveFunction <- function(x) {
        logl <- lav_bvmix_logl(rho = tanh(x[1L]),
                               fit.y1 = fit.y1, fit.y2 = fit.y2, wt = wt)
        -logl
    }

    gradientFunction <- function(x) {
        rho <- tanh(x[1L])
        R <- sqrt(1 - rho*rho)
        tauj.star  <- (fit.y2$z1 - rho*Z)/R
        tauj1.star <- (fit.y2$z2 - rho*Z)/R
        y.Z1  <- dnorm(tauj.star); y.Z2 <- dnorm(tauj1.star)

        # p(y|x)
        py2y1 <- pnorm(tauj.star) - pnorm(tauj1.star)
        py2y1[py2y1 < .Machine$double.eps] <- .Machine$double.eps
        pyx.inv <- 1/py2y1

        # rho
        TAUj  <- y.Z1 * (fit.y2$z1*rho - Z)
        TAUj1 <- y.Z2 * (fit.y2$z2*rho - Z)
        dx <- pyx.inv * 1/(R*R*R) * (TAUj - TAUj1)
        if(is.null(wt)) {
            dx.rho <- sum(dx, na.rm = TRUE)
        } else {
            dx.rho <- sum(wt * dx, na.rm = TRUE)
        }

        # tanh + minimize
        -dx.rho * 1/(cosh(x)*cosh(x))
    }

    # starting value -- Olsson 1982 eq 38
    if(length(fit.y2$slope.idx) > 0L) {
        # exo
        if(is.null(wt)) {
            COR <- cor(Z, Y2, use = "pairwise.complete.obs")
            SD <- sd(Y2, na.rm = TRUE)
        } else {
            tmp <- na.omit(cbind(Z, Y2, wt))
            COR <- cov.wt(x = tmp[,1:2], wt = tmp[,3], cor = TRUE)$cor[2,1]
            SD <- sqrt(lav_matrix_var_wt(tmp[,2], wt = tmp[,3]))
        }
        rho.init <- ( COR * SD / sum(dnorm(fit.y2$theta[fit.y2$th.idx])) )
    } else {
        # no exo
        if(is.null(wt)) {
            COR <- cor(Y1, Y2, use = "pairwise.complete.obs")
            SD <- sd(Y2, na.rm = TRUE)
        } else {
            tmp <- na.omit(cbind(Y1, Y2, wt))
            COR <- cov.wt(x = tmp[,1:2], wt = tmp[,3], cor = TRUE)$cor[2,1]
            SD <- sqrt(lav_matrix_var_wt(tmp[,2], wt = tmp[,3]))
        }
        rho.init <- ( COR * SD / sum(dnorm(fit.y2$theta[fit.y2$th.idx])) )
    }

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
        out <- optim(par = atanh(rho.init), fn = objectiveFunction,
                     gr = gradientFunction,
                     control = list(parscale = 1, reltol = 1e-10,
                                    trace = ifelse(verbose, 1L, 0L),
                                    abstol=(.Machine$double.eps * 10)),
                     method = "BFGS")
        if(out$convergence != 0L) {
            warning("no convergence")
        }
        rho <- tanh(out$par)
    }

    rho
}

# casewise scores
#
# Y1 = linear
# Y2 = ordinal
lav_bvmix_cor_scores <- function(Y1, Y2, eXo = NULL, wt = NULL,
                                 rho = NULL,
                                 fit.y1 = NULL, fit.y2 = NULL,
                                 evar.y1 = NULL, beta.y1 = NULL,
                                 th.y2 = NULL, sl.y2 = NULL,
                                 sigma.correction = FALSE,
                                 na.zero = FALSE) {
    if(is.null(fit.y1)) {
        fit.y1 <- lav_uvreg_fit(y = Y1, X = eXo, wt = wt)
    } else {
        Y1 <- fit.y1$y
        eXo <- fit.y1$X
    }
    if(is.null(fit.y2)) {
        fit.y2 <- lav_uvord_fit(y = Y2, X = eXo, wt = wt)
    } else {
        Y2 <- fit.y2$y
    }

    # update z1/z2 if needed (used in pml_deriv1() in lav_model_gradient_pml.R)
    fit.y1 <- lav_uvreg_update_fit(fit.y = fit.y1, evar.new = evar.y1,
                                   beta.new = beta.y1)
    fit.y2 <- lav_uvord_update_fit(fit.y = fit.y2, th.new = th.y2,
                                   sl.new = sl.y2)

    R <- sqrt(1 - rho*rho)
    y1.VAR <- fit.y1$theta[fit.y1$var.idx]; y1.SD  <- sqrt(y1.VAR)
    y1.ETA <- fit.y1$yhat
    Z <- (Y1 - y1.ETA) / y1.SD

    tauj.star  <- (fit.y2$z1 - rho*Z)/R
    tauj1.star <- (fit.y2$z2 - rho*Z)/R
    y.Z1  <- dnorm(tauj.star); y.Z2 <- dnorm(tauj1.star)

    # p(Y2|Y1)
    py2y1 <- pnorm(tauj.star) - pnorm(tauj1.star)
    py2y1[py2y1 < .Machine$double.eps] <- .Machine$double.eps
    pyx.inv <- 1/py2y1

    # mu.y1
    y.Z1.y.Z2 <- y.Z1-y.Z2
    dx.mu.y1 <- 1/y1.SD * (Z + (pyx.inv * (rho/R) * y.Z1.y.Z2))
    if(!is.null(wt)) {
        dx.mu.y1 <- wt * dx.mu.y1
    }

    # var.y1
    dx.var.y1 <-  1/(2*y1.VAR) * ( ((Z*Z)-1) + (pyx.inv*rho*Z/R)*(y.Z1.y.Z2) )
    if(!is.null(wt)) {
        dx.var.y1 <- wt * dx.var.y1
    }

    # th.y2
    dx.th.y2 <- (fit.y2$Y1*y.Z1 - fit.y2$Y2*y.Z2) * 1/R * pyx.inv
    if(!is.null(wt)) {
        dx.th.y2 <- wt * dx.th.y2
    }

    # sl.y1
    dx.sl.y1 <- NULL
    if(length(fit.y1$slope.idx) > 0L) {
        dx.sl.y1 <- dx.mu.y1 * eXo
        #if(!is.null(wt)) {
            # dx.mu.y1 had already been weighted
        #}
    }

    # sl.y2
    dx.sl.y2 <- NULL
    if(length(fit.y2$slope.idx) > 0L) {
        dx.sl.y2 <- (y.Z2 - y.Z1) * eXo * 1/R * pyx.inv
        if(!is.null(wt)) {
            dx.sl.y2 <- wt * dx.sl.y2
        }
    }

    # rho
    TAUj  <- y.Z1 * (fit.y2$z1*rho - Z)
    TAUj1 <- y.Z2 * (fit.y2$z2*rho - Z)
    dx.rho <- pyx.inv * 1/(R*R*R) * (TAUj - TAUj1)
    if(!is.null(wt)) {
        dx.rho <- wt * dx.rho
    }

    # FIXME: only tested for non_exo!
    # used by pml_deriv1()
    if(sigma.correction) {
        dx.rho.orig <- dx.rho
        dx.var.y1.orig <- dx.var.y1

        # sigma
        dx.rho <- dx.rho.orig / y1.SD

        # var
        COV <- rho * y1.SD
        dx.var.y1 <- ( dx.var.y1.orig -
                       1/2 * COV/y1.VAR * 1/y1.SD * dx.rho.orig )
    }


    list(dx.mu.y1 = dx.mu.y1, dx.var.y1 = dx.var.y1, dx.th.y2 = dx.th.y2,
         dx.sl.y1 = dx.sl.y1, dx.sl.y2 = dx.sl.y2, dx.rho = dx.rho)
}



