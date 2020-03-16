# the weighted bivariate linear regression model
# YR 14 March 2020 ((replacing the old lav_pearson.R routines)
#
# - pearson correlation
# - bivariate linear regression
# - using sampling weights wt

# (summed) loglikelihood
lav_bvreg_logl <- function(Y1 = NULL, Y2 = NULL, eXo = NULL, wt = NULL,
                           rho = NULL, fit.y1 = NULL, fit.y2 = NULL) {

    loglik <- lav_bvreg_lik(Y1 = Y1, Y2 = Y2, eXo = eXo, wt = wt, rho = rho,
                            fit.y1 = fit.y1, fit.y2 = fit.y2, .log = TRUE)
    logl <- sum(loglik, na.rm = TRUE)
    # logl will be -Inf if any element in '(log)lik' is -Inf
    logl
}

# individual likelihoods
lav_bvreg_lik <- function(Y1 = NULL, Y2 = NULL, eXo = NULL, wt = NULL,
                          rho = NULL,
                          evar.y1 = NULL, beta.y1 = NULL,
                          evar.y2 = NULL, beta.y2 = NULL,
                          fit.y1 = NULL, fit.y2 = NULL, .log = FALSE) {

    stopifnot(!is.null(rho))

    if(is.null(fit.y1)) {
        fit.y1 <- lav_uvreg_fit(y = Y1, X = eXo, wt = wt)
    } else {
        Y1 <- fit.y1$y
    }
    if(is.null(fit.y2)) {
        fit.y2 <- lav_uvreg_fit(y = Y2, X = eXo, wt = wt)
    } else {
        Y2 <- fit.y2$y
    }

    # user specified parameters
    if(!is.null(evar.y1) || !is.null(beta.y1)) {
        fit.y1 <- lav_uvreg_update_fit(fit.y = fit.y1,
                                       evar.new = evar.y1, beta.new = beta.y1)
    }
    if(!is.null(evar.y2) || !is.null(beta.y2)) {
        fit.y2 <- lav_uvreg_update_fit(fit.y = fit.y2,
                                       evar.new = evar.y2, beta.new = beta.y2)
    }

    # lik
    evar.y1 <- fit.y1$theta[fit.y1$var.idx]
    evar.y2 <- fit.y2$theta[fit.y1$var.idx]
    cov.y12 <- rho * sqrt(evar.y1) * sqrt(evar.y2)
    sigma <- matrix(c(evar.y1, cov.y12, cov.y12, evar.y2), 2L, 2L)
    eta.y1 <- fit.y1$yhat
    eta.y2 <- fit.y2$yhat
    lik <- dmnorm( cbind(Y1, Y2), mean = cbind(eta.y1, eta.y2), varcov = sigma)

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

# pearson correlation
# if no missing, solution is just cor(Y1,Y2) or cor(e1,e2)
# but if missing, two-step solution is NOT the same as cor(Y1,Y2) or cor(e1,e2)
lav_bvreg_cor_twostep <- function(Y1, Y2, eXo = NULL, wt = NULL,
                                  fit.y1 = NULL, fit.y2 = NULL,
                                  method = "nlminb", verbose = FALSE) {

    if(is.null(fit.y1)) {
        fit.y1 <- lav_uvreg_fit(y = Y1, X = eXo, wt = wt)
    } else {
        Y1 <- fit.y1$y
        eXo <- fit.y1$X
    }
    if(is.null(fit.y2)) {
        fit.y2 <- lav_uvreg_fit(y = Y2, X = eXo, wt = wt)
    } else {
        Y2 <- fit.y2$y
    }

    # the complete case
    if(!anyNA(Y1) && !anyNA(Y2)) {
        if(fit.y1$nexo > 0L) {
            E1 <- Y1 - fit.y1$yhat
            E2 <- Y2 - fit.y2$yhat
            if(is.null(wt)) {
                rho <- cor(E1, E2)
            } else {
                rho <- cov.wt(cbind(E1, E2), wt = wt, cor = TRUE)$cor[2,1]
            }
        } else {
            if(is.null(wt)) {
                rho <- cor(Y1, Y2)
            } else {
                rho <- cov.wt(cbind(Y1, Y2), wt = wt, cor = TRUE)$cor[2,1]
            }
        }
        return(rho)
    }

    Y1c <- Y1 - fit.y1$yhat
    Y2c <- Y2 - fit.y2$yhat
    var.y1 <- fit.y1$theta[fit.y1$var.idx]; sd.y1 <- sqrt(var.y1)
    var.y2 <- fit.y2$theta[fit.y2$var.idx]; sd.y2 <- sqrt(var.y2)

    objectiveFunction <- function(x) {
        rho = tanh(x[1L])
        logl <- lav_bvreg_logl(rho = rho, fit.y1 = fit.y1, fit.y2 = fit.y2,
                               wt = wt)
        -logl
    }

    gradientFunction <- function(x) {
        rho = tanh(x[1L])
        R <- (1 - rho*rho)
        z <- (Y1c*Y1c)/var.y1 - 2*rho*Y1c*Y2c/(sd.y1*sd.y2) + (Y2c*Y2c)/var.y2
        dx <- rho/R + (Y1c*Y2c/(sd.y1*sd.y2*R) - z*rho/(R*R))
        if(is.null(wt)) {
            dx.rho <- sum(dx, na.rm = TRUE)
        } else {
            dx.rho <- sum(wt * dx, na.rm = TRUE)
        }

        # tanh + minimize
        -dx.rho * 1/(cosh(x)*cosh(x))
    }

    # starting value (just to see if the gradient works)
    rho.init <- 0

    # minimize
    out <- nlminb(start = atanh(rho.init), objective = objectiveFunction,
                  gradient = gradientFunction, scale = 10,
                  control = list(trace = ifelse(verbose,1L,0L), rel.tol=1e-10))
    if(out$convergence != 0L) {
        warning("no convergence")
    }
    rho <- tanh(out$par)

    rho
}


# casewise scores
lav_bvreg_cor_scores <- function(Y1, Y2, eXo = NULL, wt = NULL,
                                 rho = NULL,
                                 fit.y1 = NULL, fit.y2 = NULL) {

    if(is.null(fit.y1)) {
        fit.y1 <- lav_uvreg_fit(y = Y1, X = eXo, wt = wt)
    } else {
        Y1 <- fit.y1$y
        eXo <- fit.y1$X
    }
    if(is.null(fit.y2)) {
        fit.y2 <- lav_uvreg_fit(y = Y2, X = eXo, wt = wt)
    } else {
        Y2 <- fit.y2$y
    }

    R <- (1 - rho*rho)

    Y1c <- Y1 - fit.y1$yhat
    Y2c <- Y2 - fit.y2$yhat
    var.y1 <- fit.y1$theta[fit.y1$var.idx]; sd.y1 <- sqrt(var.y1)
    var.y2 <- fit.y2$theta[fit.y2$var.idx]; sd.y2 <- sqrt(var.y2)

    # mu.y1
    dx.mu.y1 <- (2*Y1c/var.y1 - 2*rho*Y2c/(sd.y1*sd.y2))/(2*R)
    if(!is.null(wt)) {
        dx.mu.y1 <- wt * dx.mu.y1
    }

    # mu.y2
    dx.mu.y2 <- - (2*rho*Y1c/(sd.y1*sd.y2) - 2*Y2c/var.y2)/(2*R)
    if(!is.null(wt)) {
        dx.mu.y2 <- wt * dx.mu.y2
    }

    # var.y1
    dx.var.y1 <- - ( 0.5/var.y1 - ( (Y1c*Y1c)/(var.y1*var.y1) -
                       rho*Y1c*Y2c/(var.y1*sd.y1*sd.y2) ) / (2*R) )
    if(!is.null(wt)) {
        dx.var.y1 <- wt * dx.var.y1
    }

    # var.y2
    dx.var.y2 <- - ( 0.5/var.y2 + ( rho*Y1c*Y2c/(var.y2*sd.y1*sd.y2) -
                     (Y2c*Y2c)/(var.y2*var.y2) ) / (2*R) )
    if(!is.null(wt)) {
        dx.var.y2 <- wt * dx.var.y2
    }

    # sl.y1
    dx.sl.y1 <- NULL
    if(length(fit.y1$slope.idx) > 0L) {
        dx.sl.y1 <- dx.mu.y1 * eXo
        # weights already included in dx.mu.y1
    }

    # sl.y2
    dx.sl.y2 <- NULL
    if(length(fit.y2$slope.idx) > 0L) {
        dx.sl.y2 <- dx.mu.y2 * eXo
        # weights already included in dx.mu.y2
    }

    # rho
    z <- (Y1c*Y1c)/var.y1 - 2*rho*Y1c*Y2c/(sd.y1*sd.y2) + (Y2c*Y2c)/var.y2
    dx.rho <- rho/R + (Y1c*Y2c/(sd.y1*sd.y2*R) - z*rho/(R*R))
    if(!is.null(wt)) {
        dx.rho <- wt * dx.rho
    }

    list(dx.mu.y1 = dx.mu.y1, dx.var.y1 = dx.var.y1,
         dx.mu.y2 = dx.mu.y2, dx.var.y2 = dx.var.y2,
         dx.sl.y1 = dx.sl.y1, dx.sl.y2 = dx.sl.y2,
         dx.rho = dx.rho)
}

