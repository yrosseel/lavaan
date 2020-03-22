# the weighted bivariate linear regression model
# YR 14 March 2020 ((replacing the old lav_pearson.R + lav_binorm.R routines)
#
# - bivariate standard normal
# - pearson correlation
# - bivariate linear regression
# - using sampling weights wt


# density of a bivariate __standard__ normal
lav_dbinorm <- dbinorm <- function(u, v, rho, force.zero = FALSE) {
    # dirty hack to handle extreme large values for rho
    # note that u, v, and rho are vectorized!
    RHO.limit <- 0.9999
    abs.rho <- abs(rho); idx <- which(abs.rho > RHO.limit)
    if(length(idx) > 0L) { 
        rho[idx] <- sign(rho[idx]) * RHO.limit
    }

    R <- 1 - rho*rho
    out <- 1/(2*pi*sqrt(R)) * exp( - 0.5*(u*u - 2*rho*u*v + v*v)/R )

    # if abs(u) or abs(v) are very large (say, >10), set result equal 
    # to exactly zero
    idx <- which( abs(u) > 10 | abs(v) > 10)
    if(length(idx) > 0L && force.zero) {
        out[idx] <- 0    
    }

    out
}

# partial derivative - rho
lav_dbinorm_drho <- function(u, v, rho) {
    R <- 1 - rho*rho
    dbinorm(u,v,rho) * (u*v*R -rho*(u*u - 2*rho*u*v + v*v) + rho*R )/(R*R)
}

# partial derivative - u
lav_dbinorm_du <- function(u, v, rho) {
    R <- 1 - rho*rho
    -dbinorm(u,v,rho) * (u - rho*v)/R
}

# partial derivative - v
lav_dbinorm_dv <- function(u, v, rho) {
    R <- 1 - rho*rho
    -dbinorm(u,v,rho) * (v - rho*u)/R
}


# CDF of bivariate standard normal
# function pbinorm(upper.x, upper.y, rho)

# partial derivative pbinorm - upper.x
lav_pbinorm_dupperx <- function(upper.x, upper.y, rho=0.0) {
    R <- 1 - rho*rho
    dnorm(upper.x) * pnorm( (upper.y - rho*upper.x)/sqrt(R) )
}

lav_pbinorm_duppery <- function(upper.x, upper.y, rho=0.0) {
    R <- 1 - rho*rho
    dnorm(upper.y) * pnorm( (upper.x - rho*upper.y)/sqrt(R) )
}

lav_pbinorm_drho <- function(upper.x, upper.y, rho=0.0) {
    dbinorm(upper.x, upper.y, rho)
}


# switch between pbivnorm, mnormt, ...
pbinorm <- function(upper.x=NULL, upper.y=NULL, rho=0.0,
                    lower.x=-Inf, lower.y=-Inf, check=FALSE) {

    pbinorm2(upper.x=upper.x, upper.y=upper.y, rho=rho,
             lower.x=lower.x, lower.y=lower.y, check=check)

}

# using vectorized version (a la pbivnorm)
pbinorm2 <- function(upper.x=NULL, upper.y=NULL, rho=0.0,
                     lower.x=-Inf, lower.y=-Inf, check=FALSE) {

    N <- length(upper.x)
    stopifnot(length(upper.y) == N)
    if(N > 1L) {
        if(length(rho) == 1L)
            rho <- rep(rho, N)
        if(length(lower.x) == 1L)
            lower.x <- rep(lower.x, N)
        if(length(lower.y) == 1L)
            lower.y <- rep(lower.y, N)
    }

   upper.only <- all(lower.x == -Inf & lower.y == -Inf)
   if(upper.only) {
        upper.x[upper.x == +Inf] <-  exp(10) # better pnorm?
        upper.y[upper.y == +Inf] <-  exp(10)
        upper.x[upper.x == -Inf] <- -exp(10)
        upper.y[upper.y == -Inf] <- -exp(10)
        res <- pbivnorm(upper.x, upper.y, rho=rho)
    } else {
        # pbivnorm does not handle -Inf well...
        lower.x[lower.x == -Inf] <- -exp(10)
        lower.y[lower.y == -Inf] <- -exp(10)
        res <- pbivnorm(upper.x, upper.y, rho=rho) -
               pbivnorm(lower.x, upper.y, rho=rho) -
               pbivnorm(upper.x, lower.y, rho=rho) +
               pbivnorm(lower.x, lower.y, rho=rho)
    }

    res
}




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

