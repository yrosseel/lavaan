# Note: this is just a stub to get us going
# we will revise this sooner or later! June 16 2012


# this is MASS::polr with some minor changes
# YR May 2010
# - y is assumed numeric (not a factor)
# - if y is binary, we do simple probit regression
# - no weights: wt = rep(1, n) 
# - no offsets: offset = rep(0, n)
# - always returns gradient
# - optionally returns hessian
# - if x=NULL, only intercepts/thresholds
uni.ordinal <- function(x, y, method="probit", hessian=FALSE)
{
    # covariates?
    if(is.null(x)) {
        x <- matrix(1, nrow=length(y), ncol=0)
    }

    nlev <- length(unique(y))
    logit <- function(p) log(p/(1 - p))
    n <- nrow(x)
    wt <- rep(1, n)

    if(nlev == 2) { # binary case only
        y <- y - 1
        X <- cbind(rep(1,n), x)
        fmin <- function(beta) {
            p <- pnorm(X %*% beta)
            -sum(2 * wt * ifelse(y, log(p), log(1-p)))
        }
        gmin <- function(beta) {
            eta <- X %*% beta; p <- pnorm(eta)
            pre <- -2 * (wt *dnorm(eta) * ifelse(y, 1/p, -1/(1-p)))
            out <- t(pre) %*% X 
            out
        }
        gmin2 <- function(beta) {
            eta <- X %*% beta; p <- pnorm(eta)
            pre <- -2 * (wt *dnorm(eta) * ifelse(y, 1/p, -1/(1-p)))
            ## FIXME !!!
            out <- pre * X
            out
        }
        s0 <- rep(0, ncol(X))
        res <- optim(s0, fmin, gmin, method="BFGS", hessian=hessian)
        beta <- res$par[-1]
        intercept <- -1 * res$par[1] 
        dx <- gmin2(res$par)
        ## WARNING: we may need to change signs in DX/HESSIAN???
        if(hessian) {
            out <- list(slopes=beta, intercepts=intercept, dx=dx, hessian=res$hessian)
        } else {
            out <- list(slopes=beta, intercepts=intercept, dx=dx)
        }
        return(out)
    }

    fmin <- function(beta) {
        theta <- beta[pc + 1:q]
        gamm <- c(-100, cumsum(c(theta[1], exp(theta[-1]))), 100)
        eta <- offset
        if (pc > 0)
            eta <- eta + drop(x %*% beta[1:pc])
        pr <- pfun(gamm[y + 1] - eta) - pfun(gamm[y] - eta)
        if (all(pr > 0))
            -sum(wt * log(pr))
        else Inf
    }

    gmin <- function(beta)
    {
        jacobian <- function(theta) { ## dgamma by dtheta matrix
            k <- length(theta)
            etheta <- exp(theta)
            mat <- matrix(0 , k, k)
            mat[, 1] <- rep(1, k)
            for (i in 2:k) mat[i:k, i] <- etheta[i]
            mat
        }
        theta <- beta[pc + 1:q]
        gamm <- c(-100, cumsum(c(theta[1], exp(theta[-1]))), 100)
        eta <- offset
        if(pc > 0) eta <- eta + drop(x %*% beta[1:pc])
        pr <- pfun(gamm[y+1] - eta) - pfun(gamm[y] - eta)
        p1 <- dfun(gamm[y+1] - eta)
        p2 <- dfun(gamm[y] - eta)
        g1 <- if(pc > 0) t(x) %*% (wt*(p1 - p2)/pr) else numeric(0)
        xx <- .polrY1*p1 - .polrY2*p2
        g2 <- - t(xx) %*% (wt/pr)
        g2 <- t(g2) %*% jacobian(theta)
        if(all(pr > 0)) c(g1, g2) else rep(NA, pc+q)
    }

    pfun <- switch(method, logistic = plogis, probit = pnorm,
                   cloglog = pgumbel, cauchit = pcauchy)
    dfun <- switch(method, logistic = dlogis, probit = dnorm,
                   cloglog = dgumbel, cauchit = dcauchy)
    #n <- nrow(x)
    pc <- ncol(x)
    #lev <- levels(y)
    #if(length(lev) <= 2) stop("response must have 3 or more levels")
    y <- unclass(y)
    q <- nlev - 1
    Y <- matrix(0, n, q)
    .polrY1 <- col(Y) == y
    .polrY2 <- col(Y) == y - 1
    #wt <- rep(1, n)
    offset <- rep(0, n)
    # always get starting values
        # try logistic/probit regression on 'middle' cut
        q1 <- nlev %/% 2
        y1 <- (y > q1)
        X <- cbind(Intercept = rep(1, n), x)
        fit <-
            switch(method,
                   "logistic"= glm.fit(X, y1, wt, family = binomial(), offset = offset),
                   "probit" = glm.fit(X, y1, wt, family = binomial("probit"), offset = offset),
                   ## this is deliberate, a better starting point
                   "cloglog" = glm.fit(X, y1, wt, family = binomial("probit"), offset = offset),
                   "cauchit" = glm.fit(X, y1, wt, family = binomial("cauchit"), offset = offset))
        if(!fit$converged)
            stop("attempt to find suitable starting values failed")
        coefs <- fit$coefficients
        if(any(is.na(coefs))) {
            warning("design appears to be rank-deficient, so dropping some coefs")
            keep <- names(coefs)[!is.na(coefs)]
            coefs <- coefs[keep]
            x <- x[, keep[-1], drop = FALSE]
            pc <- ncol(x)
        }
        spacing <- logit((1:q)/(q+1)) # just a guess
        if(method != "logistic") spacing <- spacing/1.7
        gammas <- -coefs[1] + spacing - spacing[q1]
        thetas <- c(gammas[1], log(diff(gammas)))
        s0 <- c(coefs[-1], thetas)
    # end starting values
    res <- optim(s0, fmin, gmin, method="BFGS", hessian=hessian)
    beta <- res$par[seq_len(pc)]
    theta <- res$par[pc + 1:q]
    zeta <- cumsum(c(theta[1],exp(theta[-1])))
    #names(zeta) <- paste(lev[-length(lev)], lev[-1], sep="|")
    if(pc > 0) {
        #names(beta) <- colnames(x)
        eta <- drop(x %*% beta)
    } else {
        eta <- rep(0, n)
    }
    dx <- gmin(res$par)
    if(hessian) {
        list(slopes=beta, intercepts=zeta, dx=dx, hessian=res$hessian)
    } else {
        list(slopes=beta, intercepts=zeta, dx=dx)
    }
}
