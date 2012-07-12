
numLogl <- function(X, mu.x, var.x) {
    logl <- sum(dnorm(X, mean=mu.x, sd=sqrt(var.x), log=TRUE))
    logl
}

numLogl_x <- function(x, X) {
    mu.x=x[1L]; var.x=x[2L]
    logl <- numLogl(X, mu.x=mu.x, var.x=var.x)
    logl
}

ml_mu <- function(X, verbose=FALSE) {
    N <- length(X); var.x <- var(X)*(N-1)/N
    objectiveFunction <- function(x) {
        mu.x=x[1L]; logl <- numLogl(X, mu.x=mu.x, var.x=var.x)
        -logl
    }

    gradientFunction <- function(x) {
        mu.x=x[1L]; dx.mu.x <- sum( 1/var.x * (X-mu.x) )
        -dx.mu.x
    }
    
    out <- nlminb(start=0, objective=objectiveFunction,
                  gradient=gradientFunction,  
                  control=list(trace=ifelse(verbose,1L,0L), rel.tol=1e-10))
    if(out$convergence != 0L) warning("no convergence")
    out$par
}

ml_var <- function(X, verbose=FALSE) {
    mu.x = mean(X)
    objectiveFunction <- function(x) {
        var.x=x[1L]
        logl <- numLogl(X, mu.x=mu.x, var.x=var.x)
        -logl
    }
    gradientFunction <- function(x) {
        var.x=x[1L]; dx.var.x <- sum( -1/(2*var.x) + (X-mu.x)^2/(2*var.x^2) )
        -dx.var.x
    }
   
    out <- nlminb(start=1, objective=objectiveFunction,
                  gradient=gradientFunction,
                  control=list(trace=ifelse(verbose,1L,0L), rel.tol=1e-10))
    if(out$convergence != 0L) warning("no convergence")
    out$par
}

ml_muvar <- function(X, verbose=FALSE) {
    objectiveFunction <- function(x) {
        mu.x=x[1L]; var.x=x[2L]
        logl <- numLogl(X, mu.x=mu.x, var.x=var.x)
        -logl
    }
    gradientFunction <- function(x) {
        mu.x=x[1L]; var.x=x[2L]
        dx.mu.x <- sum( 1/var.x * (X-mu.x) )
        dx.var.x <- sum( -1/(2*var.x) + (X-mu.x)^2/(2*var.x^2) )
        -c(dx.mu.x, dx.var.x)
    }
  
    out <- nlminb(start=c(0,1), objective=objectiveFunction,
                  gradient=gradientFunction,
                  control=list(trace=ifelse(verbose,1L,0L), rel.tol=1e-10))
    if(out$convergence != 0L) warning("no convergence")
    out$par
}

scores_mu <- function(X, mu.x=NULL, var.x=NULL) {
    if(is.null(mu.x)) mu.x <- mean(X)
    if(is.null(var.x)) { N <- length(X); var.x <- var(X)*(N-1)/N }
    1/var.x * (X-mu.x)
}

scores_var <- function(X, mu.x=NULL, var.x=NULL) {
    if(is.null(mu.x)) mu.x <- mean(X)
    if(is.null(var.x)) { N <- length(X); var.x <- var(X)*(N-1)/N }
     -1/(2*var.x) + (X-mu.x)^2/(2*var.x^2)
}

scores_cor <- function(X, Y, rho=NULL, 
                       mu.x=NULL, var.x=NULL, mu.y=NULL, var.y=NULL) {
    R <- (1-rho^2); sd.x <- sqrt(var.x); sd.y <- sqrt(var.y)
    Xc <- X-mu.x; Yc <- Y-mu.y

    # rho
    z <- Xc^2/var.x - 2*rho*Xc*Yc/(sd.x*sd.y) + Yc^2/var.y
    dx.rho <- rho/R + (Xc*Yc/(sd.x*sd.y*R) - z*rho/R^2)
    dx.rho
}

# pearson product moment correlation 
# cor(X,Y)
ppcor_TS <- function(X,Y, mu.x=NULL, var.x=NULL, mu.y=NULL, var.y=NULL,
                     verbose=TRUE) {
    if(is.null(mu.x)) mu.x <- mean(X); if(is.null(mu.y)) mu.y <- mean(Y)
    if(is.null(var.x)) { N <- length(X); var.x <- var(X)*(N-1)/N }
    if(is.null(var.y)) { N <- length(Y); var.y <- var(Y)*(N-1)/N }
    sd.x <- sqrt(var.x); sd.y <- sqrt(var.y)
    ZX <- (X - mu.x)/sd.x; ZY <- (Y - mu.y)/sd.y
    objectiveFunction <- function(x) {
        rho=x[1L]
        cat("rho = ", rho, "\n")
        logl <- ml_bi(ZX=ZX, ZY=ZY, rho=rho) 
        -logl
    }
    gradientFunction <- function(x) {
        rho=x[1L]; R <- (1-rho^2)
        dx.rho <- sum( rho/R*(1 - (ZX^2 - 2*rho*ZX*ZY + ZY^2)/R) + ZX*ZY/R )
        -dx.rho
    }
  
    out <- nlminb(start=0, objective=objectiveFunction,
                  gradient=gradientFunction,
                  scale=10,
                  control=list(trace=ifelse(verbose,1L,0L), rel.tol=1e-10))
    if(out$convergence != 0L) warning("no convergence")
    out$par
}

scores_cor_uni <- function(X, Y, rho=NULL, mu.x=NULL, var.x=NULL,
                                           mu.y=NULL, var.y=NULL) {
    R <- (1-rho^2); sd.x <- sqrt(var.x); sd.y <- sqrt(var.y)
    Xc <- X-mu.x; Yc <- Y-mu.y

    # mu.x
    dx.mu.x <- (2*Xc/var.x - 2*rho*Yc/(sd.x*sd.y))/(2*R)

    # mu.y
    dx.mu.y <- - (2*rho*Xc/(sd.x*sd.y) - 2*Yc/var.y)/(2*R)

    # var.x
    dx.var.x <- - (0.5/var.x -
                   (Xc^2/var.x^2 - rho*Xc*Yc/(var.x*sd.x*sd.y))/(2*R))

    # var.y
    dx.var.y <- -(0.5/var.y +
                  (rho*Xc*Yc/(var.y*sd.x*sd.y) - Yc^2/var.y^2)/(2*R))

    # rho
    #z <- Xc^2/var.x - 2*rho*Xc*Yc/(sd.x*sd.y) + Yc^2/var.y
    #dx.rho <- rho/R + (Xc*Yc/(sd.x*sd.y*R) - z*rho/R^2)

    #dx <- cbind(dx.mu.x, dx.mu.y, dx.var.x, dx.var.y, dx.rho)
    list(dx.mu.x=dx.mu.x, dx.var.x=dx.var.x, 
         dx.mu.y=dx.mu.y, dx.var.x=dx.var.x)
}

scores_corX_uni <- function(X, Y, rho=NULL, eXo=NULL,
                            eta.x=NULL, var.x=NULL, eta.y=NULL, var.y=NULL) {
    R <- (1-rho^2); sd.x <- sqrt(var.x); sd.y <- sqrt(var.y)
    Xc <- X-eta.x; Yc <- Y-eta.y

    # mu.x
    dx.mu.x <- (2*Xc/var.x - 2*rho*Yc/(sd.x*sd.y))/(2*R)

    # mu.y
    dx.mu.y <- - (2*rho*Xc/(sd.x*sd.y) - 2*Yc/var.y)/(2*R)

    # var.x
    dx.var.x <- - (0.5/var.x -
                   (Xc^2/var.x^2 - rho*Xc*Yc/(var.x*sd.x*sd.y))/(2*R))

    # var.y
    dx.var.y <- -(0.5/var.y +
                  (rho*Xc*Yc/(var.y*sd.x*sd.y) - Yc^2/var.y^2)/(2*R))

    # sl.x
    dx.sl.x <- dx.mu.x * eXo

    # sl.y
    dx.sl.y <- dx.mu.y * eXo

    # rho
    #z <- Xc^2/var.x - 2*rho*Xc*Yc/(sd.x*sd.y) + Yc^2/var.y
    #dx.rho <- rho/R + (Xc*Yc/(sd.x*sd.y*R) - z*rho/R^2)

    #dx <- cbind(dx.mu.x, dx.mu.y, dx.var.x, dx.var.y, dx.rho)
    list(dx.mu.x=dx.mu.x, dx.var.x=dx.var.x, dx.mu.y=dx.mu.y, dx.var.y=dx.var.y,
         dx.sl.x=dx.sl.x, dx.sl.y=dx.sl.y)
}


F_ij_ps <- function(x, X, Y, nth.y, verbose=FALSE) {
    rho = x[1L]
    mu.x = x[2L]
    var.x = x[3L]
    th.y = x[3L + 1:nth.y]

    Z <- (X - mu.x) / sqrt(var.x)
    logl <- psLogl(rho=rho, mu.x=mu.x, var.x=var.x, th.y=th.y, X=X, Y=Y,
                   Z=Z)
    logl
}
#numDeriv:::grad(F_ij_ps, x=x.par)
#scores <- scores_pscor_uni(X=X, Y=Y, rho=rho, mu.x=mu.x, var.x=var.x, th.y=th.y)
#apply(scores, 2, sum)

F_ij_cor <- function(x, X, Y, verbose=FALSE) {

    Lx <- function(X, Y, rho, mu.x, var.x, mu.y, var.y) {
        cov.xy <- rho*sqrt(var.x)*sqrt(var.y)
        sigma <- matrix(c(var.x,cov.xy,cov.xy,var.y), 2L, 2L)
        dx <- numeric(length(X))
        for(i in 1:length(X))
            dx[i] <- dmvnorm(c(X[i],Y[i]), mean=c(mu.x, mu.y), sigma=sigma)
        sum(log(dx))
    }

    rho = x[1L]
    mu.x = x[2L]
    var.x = x[3L]
    mu.y = x[4L]
    var.y = x[5L]

    logl <- Lx(X=X, Y=Y, rho=rho, 
               mu.x=mu.x, var.x=var.x, mu.y=mu.y, var.y=var.y)
    logl
}
#numDeriv:::grad(F_ij_corX, x=x.par)[-1]

F_ij_psX <- function(x, X, Y, eXo, nth.y, verbose=FALSE) {

    nexo <- ncol(eXo)

    rho = x[1L]
    mu.x = x[2L]
    var.x = x[3L]
    th.y = x[3L + 1:nth.y]
    sl.x = x[3L + nth.y + 1:nexo]
    sl.y = x[3L + nth.y + nexo + 1:nexo]

    eta.x <- drop( cbind(1,eXo) %*% c(mu.x,sl.x) )
    ZRESID <- (X - eta.x)/sqrt(var.x)

    TH.y <- c(-Inf, th.y, +Inf); eta.y <- drop(eXo %*% sl.y)
    y.z1 <- pmin( 100, TH.y[Y+1L   ] - eta.y)
    y.z2 <- pmax(-100, TH.y[Y+1L-1L] - eta.y)

    if(verbose) {
        cat("rho        = ", rho, "\n")
        cat("X[1:6]     = ", X[1:6], "\n")
        cat("eta.x[1:6] = ", eta.x[1:6], "\n")
        cat("var.x      = ", var.x, "\n")
        cat("y.z1[1:6]  = ", y.z1[1:6], "\n")
        cat("y.z2[1:6]  = ", y.z2[1:6], "\n")
    }

    logl <- psLoglx(rho=rho, X=X, eta.x=eta.x, var.x=var.x, y.z1=y.z1,
                    y.z2=y.z2, ZRESID=ZRESID)
    logl
}
#numDeriv:::grad(F_ij_psX, x=x.par)[-1]

F_ij_corX <- function(x, X, Y, eXo, verbose=FALSE) {

    Lx <- function(X, Y, rho, eta.x, var.x, eta.y, var.y) {
        cov.xy <- rho*sqrt(var.x)*sqrt(var.y)
        sigma <- matrix(c(var.x,cov.xy,cov.xy,var.y), 2L, 2L)
        dx <- numeric(length(X))
        for(i in 1:length(X))
            dx[i] <- dmvnorm(c(X[i],Y[i]), mean=c(eta.x[i], eta.y[i]), sigma=sigma)
        sum(log(dx))
    }

    nexo <- ncol(eXo)

    rho = x[1L]
    mu.x = x[2L]
    var.x = x[3L]
    mu.y = x[4L]
    var.y = x[5L]
    sl.x = x[5L + 1:nexo]
    sl.y = x[5L + nexo + 1:nexo]

    eta.x <- drop( cbind(1,eXo) %*% c(mu.x,sl.x) )
    eta.y <- drop( cbind(1,eXo) %*% c(mu.y,sl.y) )

    logl <- Lx(X=X, Y=Y, rho=rho, eta.x=eta.x, var.x=var.x,
               eta.y=eta.y, var.y=var.y)
    logl
}
#numDeriv:::grad(F_ij_corX, x=x.par)[-1]




