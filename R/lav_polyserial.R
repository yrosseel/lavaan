# polyserial Y1 is numeric, Y2 is ordinal

# (summed) loglikelihood
ps_logl <- function(Y1, Y2, eXo=NULL, rho=NULL, fit.y1=NULL, fit.y2=NULL) {

    lik <- ps_lik(Y1=Y1, Y2=Y2, eXo=eXo, rho=rho, fit.y1=fit.y1, fit.y2=fit.y2)
    if(all(lik > 0, na.rm = TRUE))
        logl <- sum(log(lik), na.rm = TRUE)
    else
        logl <- -Inf
    logl
}

# individual likelihoods
ps_lik <- function(Y1, Y2, eXo=NULL, var.y1 = NULL, eta.y1 = NULL,
                   rho=NULL, fit.y1=NULL, fit.y2=NULL) {

    stopifnot(!is.null(rho))
    if(is.null(fit.y1)) fit.y1 <- lavOLS(Y1, X=eXo)
    if(is.null(fit.y2)) fit.y2 <- lavProbit(Y2, X=eXo)
    if(missing(Y1)) Y1 <- fit.y1$y
    if(is.null(var.y1)) {
        var.y1 <- fit.y1$theta[fit.y1$var.idx]
    }
    if(is.null(eta.y1)) {
        eta.y1 <- fit.y1$yhat
    }

    R <- sqrt(1 - rho*rho)
    y1.SD  <- sqrt(var.y1)
    y1.ETA <- eta.y1
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

# individual loglikelihoods - no exo!
ps_loglik_no_exo <- function(Y1, Y2, var.y1 = NULL, eta.y1 = NULL, rho=NULL,
                             cov.y12 = NULL, th.y2 = NULL) {

    if(is.null(rho)) {
        rho <- cov.y12 / sqrt(var.y1)
    }

    #if(abs(rho) > 1) {
    #    # should never happen
    #    cat("rho = ", rho, "\n")
    #    cat("cov = ", cov.y12, "\n")
    #    cat("var.y1 = ", var.y1, "\n")
    #    stop("rho > 1")
    #}

    R <- sqrt(1 - rho*rho)
    Z <- (Y1 - eta.y1) / sqrt(var.y1)

    # p(Y2|Y1)
    TH <- c(-Inf, th.y2, +Inf)
    z1 <- pmin( 100, TH[Y2 + 1L])
    z2 <- pmax(-100, TH[Y2 + 1L - 1L])

    tauj.star  <- (z1 - rho*Z)/R
    tauj1.star <- (z2 - rho*Z)/R
    py2y1 <- pnorm(tauj.star) - pnorm(tauj1.star)
    py2y1[py2y1 < .Machine$double.eps] <- .Machine$double.eps
    py2y1.log <- log(py2y1)

    # p(Y1)
    py1.log <- dnorm(Y1, mean=eta.y1, sd=sqrt(var.y1), log = TRUE)

    # loglik
    loglik <- py1.log + py2y1.log

    loglik
}

# polyserial correlation
ps_cor_TS <- function(Y1, Y2, eXo=NULL, fit.y1=NULL, fit.y2=NULL,
                      method="nlminb", verbose=FALSE) {

    stopifnot(method %in% c("nlminb", "BFGS"))

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
        dx.rho <- sum( pyx.inv * 1/(R*R*R) * (TAUj - TAUj1), na.rm = TRUE )

        # tanh + minimize
        -dx.rho * 1/(cosh(x)*cosh(x))
    }

    # FIXME:: TODO!!!
    hessianFunction <- function(x) {
    }

    # starting value -- Olsson 1982 eq 38
    if(length(fit.y2$slope.idx) > 0L) {
        # exo
        rho.init <- ( cor(Z,Y2,use="pairwise.complete.obs")*sd(Y2,na.rm=TRUE) /
                      sum(dnorm(fit.y2$theta[fit.y2$th.idx])) )
    } else {
        # no exo
        rho.init <- ( cor(Y1,Y2,use="pairwise.complete.obs")*sd(Y2,na.rm=TRUE) /
                      sum(dnorm(fit.y2$theta[fit.y2$th.idx])) )
    }

    # if rho.init is missing, set starting value to zero
    if(is.na(rho.init)) {
      rho.init <- 0.0
    }
    # check range of rho.init is within [-1,+1]
    if(abs(rho.init) >= 1.0) {
        rho.init <- 0.0
    }

    # minimize
    if(method == "nlminb") {
        out <- nlminb(start=atanh(rho.init), objective=objectiveFunction,
                      gradient=gradientFunction,
                      scale=10,
                      control=list(trace=ifelse(verbose,1L,0L), rel.tol=1e-10))
    } else if(method == "BFGS") {
        out <- optim(par = atanh(rho.init), fn = objectiveFunction,
                     gr = gradientFunction,
                     control = list(parscale = 1, reltol = 1e-10,
                                    abstol=(.Machine$double.eps * 10)),
                     method = "BFGS")
    }
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

    # var.y1
    dx.var.y1 <-  1/(2*y1.VAR) * ( ((Z*Z)-1) + (pyx.inv*rho*Z/R)*(y.Z1.y.Z2) )

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
    dx.rho <- pyx.inv * 1/(R*R*R) * (TAUj - TAUj1)


    list(dx.mu.y1=dx.mu.y1, dx.var.y1=dx.var.y1, dx.th.y2=dx.th.y2,
         dx.sl.y1=dx.sl.y1, dx.sl.y2=dx.sl.y2, dx.rho=dx.rho)
}


# note: input must be rho, not cov_y1.y2, because the latter
#       depends on var.y1, and this complicates the gradient
ps_cor_scores_no_exo <- function(Y1, Y2,
                                 var.y1 = NULL, eta.y1 = NULL,
                                 th.y2 = NULL, rho = NULL,
                                 sigma.correction = FALSE) {

    R <- sqrt(1 - rho*rho)
    Z <- (Y1 - eta.y1) / sqrt(var.y1)
    y1.SD <- sqrt(var.y1)

    nth <- length(th.y2)
    nobs <- length(Y1)
    y2.Y1 <- matrix(1:nth, nobs, nth, byrow=TRUE) == Y2
    y2.Y2 <- matrix(1:nth, nobs, nth, byrow=TRUE) == Y2 - 1L


    # p(Y2|Y1)
    TH <- c(-Inf, th.y2, +Inf)
    z1 <- pmin( 100, TH[Y2 + 1L])
    z2 <- pmax(-100, TH[Y2 + 1L - 1L])

    tauj.star  <- (z1 - rho*Z)/R
    tauj1.star <- (z2 - rho*Z)/R
    y.Z1  <- dnorm(tauj.star); y.Z2 <- dnorm(tauj1.star)

    # p(Y2|Y1)
    py2y1 <- pnorm(tauj.star) - pnorm(tauj1.star)
    py2y1[py2y1 < .Machine$double.eps] <- .Machine$double.eps
    pyx.inv <- 1/py2y1

    # mu.y1
    y.Z1.y.Z2 <- y.Z1-y.Z2
    dx.mu.y1 <- 1/y1.SD * (Z + (pyx.inv * (rho/R) * y.Z1.y.Z2))

    # var.y1
    dx.var.y1 <-  1/(2*var.y1) * ( ((Z*Z)-1) + (pyx.inv*rho*Z/R)*(y.Z1.y.Z2) )

    # th.y2
    dx.th.y2 <- (y2.Y1*y.Z1 - y2.Y2*y.Z2) * 1/R * pyx.inv

    # rho
    TAUj  <- y.Z1 * (z1*rho - Z)
    TAUj1 <- y.Z2 * (z2*rho - Z)
    dx.rho <- pyx.inv * 1/(R*R*R) * (TAUj - TAUj1)

    if(sigma.correction) {
        dx.rho.orig <- dx.rho
        dx.var.y1.orig <- dx.var.y1

        # sigma
        dx.rho <- dx.rho.orig / y1.SD

        # var
        COV <- rho * y1.SD
        dx.var.y1 <- ( dx.var.y1.orig -
                       1/2 * COV/var.y1 * 1/y1.SD * dx.rho.orig )
    }

    list(dx.mu.y1=dx.mu.y1, dx.var.y1=dx.var.y1, dx.th.y2=dx.th.y2,
         dx.rho=dx.rho)
}

