# polyserial Y1 is numeric, Y2 is ordinal

# (summed) loglikelihood
ps_logl <- function(Y1, Y2, eXo=NULL, rho=NULL, fit.y1=NULL, fit.y2=NULL) {

    lik <- ps_lik(Y1=Y1, Y2=Y2, eXo=eXo, rho=rho, fit.y1=fit.y1, fit.y2=fit.y2)
    if(all(lik > 0)) 
        logl <- sum(log(lik))
    else
        logl <- -Inf
    logl
}

# individual likelihoods
ps_lik <- function(Y1, Y2, eXo=NULL, rho=NULL, fit.y1=NULL, fit.y2=NULL) {

    stopifnot(!is.null(rho))
    if(is.null(fit.y1)) fit.y1 <- lavOLS(Y1, X=eXo)
    if(is.null(fit.y2)) fit.y2 <- lavProbit(Y2, X=eXo)
    if(missing(Y1)) Y1 <- fit.y1$y

    R <- sqrt(1-rho^2)
    y1.SD  <- sqrt(fit.y1$theta[fit.y1$var.idx])
    y1.ETA <- fit.y1$yhat
    Z <- (Y1 - y1.ETA) / y1.SD
    
    # p(Y2|Y1)
    tauj.star  <- (fit.y2$z1 - rho*Z)/R
    tauj1.star <- (fit.y2$z2 - rho*Z)/R
    py2y1 <- pnorm(tauj.star) - pnorm(tauj1.star)
    py2y1[py2y1 < .Machine$double.eps] <- .Machine$double.eps

    # p(Y1)
    py1 <- dnorm(Y1, mean=y1.ETA, sd=y1.SD)

    # lik
    py1 * py2y1
}

# loglikelihood (x-version)
ps_logl_x <- function(x, Y1, Y2, eXo=NULL, nth.y2) {

    nexo <- ifelse(is.null(eXo), 0L, ncol(eXo)); S <- seq_len
    stopifnot(length(x) == (3L + nth.y2 + 2*nexo))

       rho = x[1L]
     mu.y1 = x[2L]
    var.y1 = x[3L]
     th.y2 = x[3L + S(nth.y2)]
     sl.y1 = x[3L + nth.y2 + S(nexo)]
     sl.y2 = x[3L + nth.y2 +   nexo + S(nexo)]

    fit.y1 <- lavOLS(y=Y1, X=eXo)
    fit.y1$theta[fit.y1$int.idx] <- mu.y1
    fit.y1$theta[fit.y1$var.idx] <- var.y1
    fit.y1$theta[fit.y1$slope.idx] <- sl.y1
    fit.y1$lik()

    fit.y2 <- lavProbit(y=Y2, X=eXo)
    fit.y2$theta[fit.y2$th.idx] <- th.y2
    fit.y2$theta[fit.y2$slope.idx] <- sl.y2
    fit.y2$lik()


    ps_logl(Y1=Y1, Y2=Y2, eXo=eXo, rho=rho, fit.y1=fit.y1, fit.y2=fit.y2)
}

# polyserial correlation
ps_cor_TS <- function(Y1, Y2, eXo=NULL, fit.y1=NULL, fit.y2=NULL,
                      method="nlminb", verbose=FALSE) {

    stopifnot(method == "nlminb")

    if(is.null(fit.y1)) fit.y1 <- lavOLS(Y1, X=eXo)
    if(is.null(fit.y2)) fit.y2 <- lavProbit(Y2, X=eXo)
    if(missing(Y1)) Y1 <- fit.y1$y 
    if(missing(Y2)) Y2 <- fit.y2$y else as.integer(Y2)
    if(missing(eXo) && length(fit.y2$slope.idx) > 0L) eXo <- fit.y2$X

    y1.VAR <- fit.y1$theta[fit.y1$var.idx]; y1.SD  <- sqrt(y1.VAR)
    y1.ETA <- fit.y1$yhat
    Z <- (Y1 - y1.ETA) / y1.SD

    objectiveFunction <- function(x) {
        rho = tanh(x[1L])
        logl <- ps_logl(rho=rho, fit.y1=fit.y1, fit.y2=fit.y2)
        -logl
    }

    gradientFunction <- function(x) {
        rho = tanh(x[1L])
        R <- sqrt(1-rho^2)
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
        dx.rho <- sum( pyx.inv * 1/R^3 * (TAUj - TAUj1) )

        # tanh + minimize
        -dx.rho * 1/cosh(x)^2
    }

    # FIXME::: TODO!!!
    hessianFunction <- function(x) {
    }

    # starting value -- Olsson 1982 eq 38
    if(length(fit.y2$slope.idx) > 0L) {
        # exo
        rho.init <- cor(Z,Y2)*sd(Y2) /  sum(dnorm(fit.y2$theta[fit.y2$th.idx]))
    } else {
        # no exo
        rho.init <- cor(Y1,Y2)*sd(Y2) / sum(dnorm(fit.y2$theta[fit.y2$th.idx]))
    }

    # check range of rho.init is within [-1,+1]
    if(abs(rho.init) > 1.0) {
        rho.init <- 0.0
    }

    # minimize
    out <- nlminb(start=atanh(rho.init), objective=objectiveFunction,
                  gradient=gradientFunction,
                  scale=10,
                  control=list(trace=ifelse(verbose,1L,0L), rel.tol=1e-10))
    if(out$convergence != 0L) warning("no convergence")
    rho <- tanh(out$par)

    rho
}

ps_cor_scores <- function(Y1, Y2, eXo=NULL, rho=NULL, 
                          fit.y1=NULL, fit.y2=NULL) {

    stopifnot(!is.null(rho))
    if(is.null(fit.y1)) fit.y1 <- lavOLS(Y1, X=eXo)
    if(is.null(fit.y2)) fit.y2 <- lavProbit(Y2, X=eXo)
    if(missing(Y1)) Y1 <- fit.y1$y
    if(missing(eXo) && length(fit.y2$slope.idx) > 0L) eXo <- fit.y2$X

    R <- sqrt(1-rho^2)
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

    # var.y1
    dx.var.y1 <-  1/(2*y1.VAR) * ( (Z^2-1) + (pyx.inv*rho*Z/R)*(y.Z1.y.Z2) )

    # th.y2
    dx.th.y2 <- (fit.y2$Y1*y.Z1 - fit.y2$Y2*y.Z2) * 1/R * pyx.inv

    # sl.y1
    dx.sl.y1 <- NULL
    if(length(fit.y1$slope.idx) > 0L)
        dx.sl.y1 <- dx.mu.y1 * eXo

    # sl.y2
    dx.sl.y2 <- NULL
    if(length(fit.y2$slope.idx) > 0L)
        dx.sl.y2 <- (y.Z2 - y.Z1) * eXo * 1/R * pyx.inv

    # rho
    TAUj  <- y.Z1 * (fit.y2$z1*rho - Z)
    TAUj1 <- y.Z2 * (fit.y2$z2*rho - Z)
    dx.rho <- pyx.inv * 1/R^3 * (TAUj - TAUj1)


    list(dx.mu.y1=dx.mu.y1, dx.var.y1=dx.var.y1, dx.th.y2=dx.th.y2,
         dx.sl.y1=dx.sl.y1, dx.sl.y2=dx.sl.y2, dx.rho=dx.rho)
}
