# pearson product moment correlation: both Y1 and Y2 are NUMERIC

# (summed) loglikelihood
pp_logl <- function(Y1, Y2, eXo=NULL, rho=NULL, fit.y1=NULL, fit.y2=NULL) {

    lik <- pp_lik(Y1=Y1, Y2=Y2, eXo=eXo, rho=rho, fit.y1=fit.y1, fit.y2=fit.y2)
    if(all(lik > 0, na.rm = TRUE))
        logl <- sum(log(lik), na.rm = TRUE)
    else
        logl <- -Inf
    logl
}

# individual likelihoods
pp_lik <- function(Y1, Y2, eXo=NULL,
                   eta.y1 = NULL, eta.y2 = NULL,
                   var.y1 = NULL, var.y2 = NULL,
                   rho=NULL, fit.y1=NULL, fit.y2=NULL) {

    stopifnot(!is.null(rho))
    if(is.null(fit.y1)) fit.y1 <- lavOLS(Y1, X=eXo)
    if(is.null(fit.y2)) fit.y2 <- lavOLS(Y2, X=eXo)
    if(missing(Y1)) Y1 <- fit.y1$y
    if(missing(Y2)) Y2 <- fit.y2$y

    if(is.null(var.y1)) {
        var.y1 <- fit.y1$theta[fit.y1$var.idx]
    }
    if(is.null(var.y2)) {
        var.y2 <- fit.y2$theta[fit.y2$var.idx]
    }
    if(is.null(eta.y1)) {
        eta.y1 <- fit.y1$yhat
    }
    if(is.null(eta.y2)) {
        eta.y2 <- fit.y2$yhat
    }

    # lik
    cov.y12 <- rho*sqrt(var.y1)*sqrt(var.y2)
    sigma <- matrix(c(var.y1,cov.y12,cov.y12,var.y2), 2L, 2L)
    #lik <- numeric(length(Y1))
    #for(i in 1:length(Y1))
    #    lik[i] <- dmvnorm(c(Y1[i],Y2[i]), mean=c(eta.y1[i], eta.y2[i]),
    #                      sigma=sigma)
    lik <- dmnorm( cbind(Y1,Y2), mean=cbind(eta.y1, eta.y2), varcov=sigma)

    lik
}

# loglikelihood (x-version)
pp_logl_x <- function(x, Y1, Y2, eXo=NULL) {

    nexo <- ifelse(is.null(eXo), 0L, ncol(eXo)); S <- seq_len
    stopifnot(length(x) == (5L + 2*nexo))

       rho = x[1L]
     mu.y1 = x[2L]
    var.y1 = x[3L]
     mu.y2 = x[4L]
    var.y2 = x[5L]
     sl.y1 = x[5L + S(nexo)]
     sl.y2 = x[5L +   nexo + S(nexo)]

    fit.y1 <- lavOLS(y=Y1, X=eXo)
    fit.y1$theta[fit.y1$int.idx] <- mu.y1
    fit.y1$theta[fit.y1$var.idx] <- var.y1
    fit.y1$theta[fit.y1$slope.idx] <- sl.y1
    fit.y1$lik()

    fit.y2 <- lavOLS(y=Y2, X=eXo)
    fit.y2$theta[fit.y2$int.idx] <- mu.y2
    fit.y2$theta[fit.y2$var.idx] <- var.y2
    fit.y2$theta[fit.y2$slope.idx] <- sl.y2
    fit.y2$lik()

    pp_logl(Y1=Y1, Y2=Y2, eXo=eXo, rho=rho, fit.y1=fit.y1, fit.y2=fit.y2)
}

# pearson correlation (just for fun); solution is cor(Y1,Y2) or cor(e1,e2)
pp_cor_TS <- function(Y1, Y2, eXo=NULL, fit.y1=NULL, fit.y2=NULL,
                      method="nlminb", verbose=FALSE) {

    stopifnot(method == "nlminb")

    if(is.null(fit.y1)) fit.y1 <- lavOLS(Y1, X=eXo)
    if(is.null(fit.y2)) fit.y2 <- lavOLS(Y2, X=eXo)
    if(missing(Y1)) Y1 <- fit.y1$y
    if(missing(Y2)) Y2 <- fit.y2$y
    if(missing(eXo) && length(fit.y1$slope.idx) > 0L) eXo <- fit.y2$X[,-1]

    Y1c <- Y1 - fit.y1$yhat
    Y2c <- Y2 - fit.y2$yhat
    var.y1 <- fit.y1$theta[fit.y1$var.idx]; sd.y1 <- sqrt(var.y1)
    var.y2 <- fit.y2$theta[fit.y2$var.idx]; sd.y2 <- sqrt(var.y2)

    objectiveFunction <- function(x) {
        rho = tanh(x[1L])
        logl <- pp_logl(rho=rho, fit.y1=fit.y1, fit.y2=fit.y2)
        -logl
    }

    gradientFunction <- function(x) {
        rho = tanh(x[1L])
        R <- (1 - rho*rho)
        z <- (Y1c*Y1c)/var.y1 - 2*rho*Y1c*Y2c/(sd.y1*sd.y2) + (Y2c*Y2c)/var.y2
        dx.rho <- sum( rho/R + (Y1c*Y2c/(sd.y1*sd.y2*R) - z*rho/(R*R)),
                       na.rm = TRUE )

        # tanh + minimize
        -dx.rho * 1/(cosh(x)*cosh(x))
    }

    # FIXME:: TODO!!!
    hessianFunction <- function(x) {
    }

    # starting value (just to see if the gradient works)
    rho.init <- 0

    # minimize
    out <- nlminb(start=atanh(rho.init), objective=objectiveFunction,
                  gradient=gradientFunction,
                  scale=10,
                  control=list(trace=ifelse(verbose,1L,0L), rel.tol=1e-10))
    if(out$convergence != 0L) warning("no convergence")
    rho <- tanh(out$par)

    rho
}

pp_cor_scores <- function(Y1, Y2, eXo=NULL, rho=NULL,
                          fit.y1=NULL, fit.y2=NULL) {

    stopifnot(!is.null(rho))
    if(is.null(fit.y1)) fit.y1 <- lavOLS(Y1, X=eXo)
    if(is.null(fit.y2)) fit.y2 <- lavOLS(Y2, X=eXo)
    if(missing(Y1)) Y1 <- fit.y1$y
    if(missing(Y2)) Y2 <- fit.y2$y
    if(missing(eXo) && length(fit.y1$slope.idx) > 0L) eXo <- fit.y2$X[,-1]

    R <- (1 - rho*rho)

    Y1c <- Y1 - fit.y1$yhat
    Y2c <- Y2 - fit.y2$yhat
    var.y1 <- fit.y1$theta[fit.y1$var.idx]; sd.y1 <- sqrt(var.y1)
    var.y2 <- fit.y2$theta[fit.y2$var.idx]; sd.y2 <- sqrt(var.y2)

    # mu.y1
    dx.mu.y1 <- (2*Y1c/var.y1 - 2*rho*Y2c/(sd.y1*sd.y2))/(2*R)

    # mu.y2
    dx.mu.y2 <- - (2*rho*Y1c/(sd.y1*sd.y2) - 2*Y2c/var.y2)/(2*R)

    # var.y1
    dx.var.y1 <- - (0.5/var.y1 -
                   ((Y1c*Y1c)/(var.y1*var.y1) - rho*Y1c*Y2c/(var.y1*sd.y1*sd.y2))/(2*R))

    # var.y2
    dx.var.y2 <- -(0.5/var.y2 +
                  (rho*Y1c*Y2c/(var.y2*sd.y1*sd.y2) - (Y2c*Y2c)/(var.y2*var.y2))/(2*R))

    # sl.y1
    dx.sl.y1 <- NULL
    if(length(fit.y1$slope.idx) > 0L)
        dx.sl.y1 <- dx.mu.y1 * eXo

    # sl.y2
    dx.sl.y2 <- NULL
    if(length(fit.y2$slope.idx) > 0L)
        dx.sl.y2 <- dx.mu.y2 * eXo

    # rho
    z <- (Y1c*Y1c)/var.y1 - 2*rho*Y1c*Y2c/(sd.y1*sd.y2) + (Y2c*Y2c)/var.y2
    dx.rho <- rho/R + (Y1c*Y2c/(sd.y1*sd.y2*R) - z*rho/(R*R))

    list(dx.mu.y1=dx.mu.y1, dx.var.y1=dx.var.y1,
         dx.mu.y2=dx.mu.y2, dx.var.y2=dx.var.y2,
         dx.sl.y1=dx.sl.y1, dx.sl.y2=dx.sl.y2,
         dx.rho=dx.rho)
}

pp_cor_scores_no_exo <- function(Y1, Y2,
                                 eta.y1 = NULL, var.y1 = NULL,
                                 eta.y2 = NULL, var.y2 = NULL,
                                 rho = NULL) {

    stopifnot(!is.null(rho))

    R <- (1 - rho*rho)

    Y1c <- Y1 - eta.y1
    Y2c <- Y2 - eta.y2
    sd.y1 <- sqrt(var.y1)
    sd.y2 <- sqrt(var.y2)

    # mu.y1
    dx.mu.y1 <- (2*Y1c/var.y1 - 2*rho*Y2c/(sd.y1*sd.y2))/(2*R)

    # mu.y2
    dx.mu.y2 <- - (2*rho*Y1c/(sd.y1*sd.y2) - 2*Y2c/var.y2)/(2*R)

    # var.y1
    dx.var.y1 <- - (0.5/var.y1 -
                   ((Y1c*Y1c)/(var.y1*var.y1) - rho*Y1c*Y2c/(var.y1*sd.y1*sd.y2))/(2*R))

    # var.y2
    dx.var.y2 <- -(0.5/var.y2 +
                  (rho*Y1c*Y2c/(var.y2*sd.y1*sd.y2) - (Y2c*Y2c)/(var.y2*var.y2))/(2*R))

    # rho
    z <- (Y1c*Y1c)/var.y1 - 2*rho*Y1c*Y2c/(sd.y1*sd.y2) + (Y2c*Y2c)/var.y2
    dx.rho <- rho/R + (Y1c*Y2c/(sd.y1*sd.y2*R) - z*rho/(R*R))

    list(dx.mu.y1=dx.mu.y1, dx.var.y1=dx.var.y1,
         dx.mu.y2=dx.mu.y2, dx.var.y2=dx.var.y2,
         dx.rho=dx.rho)
}

