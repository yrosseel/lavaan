
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

# low-level function, no checking
bifreq <- function(x, y) {
    max.x <- max(x, na.rm=TRUE); max.y <- max(y, na.rm=TRUE)
    bin <- x - 1L; bin <- bin + max.x * (y - 1L); bin <- bin[!is.na(bin)]
    if (length(bin)) bin <- bin + 1L
    freq <- array(tabulate(bin, nbins = max.x*max.y), dim=c(max.x, max.y))
    freq
}


# low-level function, no checking
unithord <- function(X, freq=NULL, prop=NULL) {
    if(is.null(prop)) {
        if(is.null(freq)) freq <- tabulate(X)
        prop <- freq / sum(freq)
    }
    th <- qnorm(cumsum(prop))[-length(prop)]
    th
}

pccor_TS <- function(x, y, th.x=NULL, th.y=NULL, freq=NULL,
                  method="nlminb", verbose=FALSE,
                  acov=FALSE) {

    stopifnot(min(x) == 1L, min(y) == 1L,
              method %in% c("nlminb", "optimize", "nlminb.hessian"))
  
    # thresholds
    if(is.null(th.x)) th.x <- unithord(freq=tabulate(x))
    if(is.null(th.y)) th.y <- unithord(freq=tabulate(y))

    # create cross-table low-level (see table() function)
    if(is.null(freq)) freq <- bifreq(x,y)
    nr <- nrow(freq); nc <- ncol(freq)

    objectiveFunction <- function(x) {
        logl <- pcLogl_freq(freq, rho=tanh(x[1L]), th.x, th.y)
        #cat("rho = ", tanh(x[1L]), "\n")
        -logl
    }

    if(method == "optimize") {
        out <- optimize(f=objectiveFunction,
                        interval=c(-5, +5), # tanh scale!
                        tol=1e-8, maximum=FALSE)
        rho <- tanh(out$minimum)
    } else if(method == "nlminb" || method == "nlminb.hessian") {

        # catch tetrachoric case
        if(nr == 2L && nc == 2L) {
           # Divgi 1979 initial value
            h <- max(abs(th.x), abs(th.y)); k <- min(abs(th.x), abs(th.y))
            R <- (freq[1,1]*freq[2,2])/(freq[1,2]*freq[2,1])
            D <- k*(.79289 + 4.28981/(1+3.30231*h));D <- D*sign(th.x)*sign(th.y)
            C <- 0.07557*h + (h-k)^2 * (0.51141/(h+2.05793) - 0.07557/h)
            B <- 0.5/(1 + (h^2 + k^2)*(0.82281-1.03514*(k/sqrt(h^2+k^2))))
            A <- 0.5/(1 + (h^2 + k^2)*(0.12454-0.27102*(1-h/sqrt(h^2+k^2))))
            alpha <- A + B*(-1 + 1/(1 + C*(log(R)-D)^2))
            rho.init <- cos(pi/(1+R^alpha))
        } else {
            rho.init <- cov(x,y)
        }

        gradientFunction <- function(x) {
            rho <- tanh(x[1L])
            PI  <- pcComputePI(rho, th.x, th.y)
            phi <- pcComputephi(rho, th.x, th.y)
            dx.rho <- sum(freq/PI * phi)
            -dx.rho * 1/cosh(x)^2 # dF/drho * drho/dx, dtanh = 1/cosh(x)^2
        }

        # OLSSON 1979 A2 + A3
        hessianFunction <- function(x) {
            rho <- tanh(x[1L])
            PI  <- pcComputePI(rho, th.x, th.y)
            phi <- pcComputephi(rho, th.x, th.y)
            gnorm <- pcComputegnorm(rho, th.x, th.y)
            H <-  sum(freq/PI * gnorm) - sum(freq/PI^2 * phi^2)

            # to compensate for tanh
            # u=f(x), d^2y/dx^2 = d^2y/du^2 * (du/dx)^2 + dy/du * d^2u/dx^2
            # dtanh = 1/cosh(x)^2
            # dtanh_2 = 8*exp(2*x)*(1-exp(2*x))/(exp(2*x)+1)^3
            grad <- sum(freq/PI * phi)
            u1 <- 1/cosh(x)^2
            u2 <- 8*exp(2*x)*(1-exp(2*x))/(exp(2*x)+1)^3
            H <- H * u1^2 + grad * u2
            dim(H) <- c(1L,1L) # for nlminb

            -H
        }


        if(method == "nlminb") {
            out <- nlminb(start=atanh(rho.init), objective=objectiveFunction,
                          gradient=gradientFunction,
                          scale=10,
                          control=list(trace=ifelse(verbose,1L,0L), 
                                       rel.tol=1e-10))
        } else if(method == "nlminb.hessian") {
            # to mimic Mplus
            out <- nlminb(start=atanh(rho.init), objective=objectiveFunction,
                          gradient=gradientFunction,
                          hessian=hessianFunction,
                          scale=10, # not needed?
                          control=list(trace=ifelse(verbose,1L,0L), 
                                       rel.tol=1e-5))
        }
        if(out$convergence != 0L) warning("no convergence")
        rho <- tanh(out$par)
    }

    if(acov) {
        Gamma <- pcGamma(freq, rho=rho, th.x=th.x, th.y=th.y)
        omega <- attr(Gamma, "omega"); attr(Gamma, "omega") <- NULL
           PI <- attr(Gamma, "PI");    attr(Gamma, "PI") <- NULL
        attr(rho, "Gamma") <- Gamma
        attr(rho, "omega") <- omega
        attr(rho, "PI")    <- PI
    }

    rho
}

pccorX_TS <- function(x, y, eXo=NULL, 
                      x.z1=NULL, x.z2=NULL, y.z1=NULL, y.z2=NULL,
                      method="nlminb", verbose=FALSE, scores=FALSE) {

    stopifnot(min(x) == 1L, min(y) == 1L,
              method %in% c("nlminb", "optimize"))
  
    # thresholds
    if(is.null(x.z1) || is.null(x.z2)) {
        fit.x <- lavProbit(y=x, X=eXo)
        x.z1 <- fit.x$z1; x.z2 <- fit.x$z2
    }
    if(is.null(y.z1) || is.null(y.z2)) {
        fit.y <- lavProbit(y=y, X=eXo)
        y.z1 <- fit.y$z1; y.z2 <- fit.y$z2
    }

    objectiveFunction <- function(x) {
        rho = tanh(x[1L])
        logl <- pcLogl_i(x.z1, x.z2, y.z1, y.z2, rho)
        #cat("rho = ", tanh(x[1L]), "\n")
        -logl
    }

    if(method == "optimize") {
        out <- optimize(f=objectiveFunction,
                        interval=c(-5, +5), # tanh scale!
                        tol=1e-8, maximum=FALSE)
        rho <- tanh(out$minimum)
    } else if(method == "nlminb") {
        rho.init <- cov(x,y)

        gradientFunction <- function(x) {
            rho = tanh(x[1L])
            lik <- pcL_i(x.z1, x.z2, y.z1, y.z2, rho)
            dx <- ( dbinorm(x.z1, y.z1, rho) -
                    dbinorm(x.z2, y.z1, rho) -
                    dbinorm(x.z1, y.z2, rho) +
                    dbinorm(x.z2, y.z2, rho) ) / lik
            -sum(dx) * 1/cosh(x)^2
        }

        hessianFunction <- function(x) {
            ### FIXME: TODO
        }

        if(method == "nlminb") {
            out <- nlminb(start=atanh(rho.init), objective=objectiveFunction,
                          gradient=gradientFunction,
                          scale=10,
                          control=list(trace=ifelse(verbose,1L,0L), 
                                       rel.tol=1e-10))
        } 
        #else if(method == "nlminb.hessian") {
        #    # to mimic Mplus
        #    out <- nlminb(start=atanh(rho.init), objective=objectiveFunction,
        #                  gradient=gradientFunction,
        #                  hessian=hessianFunction,
        #                  scale=10, # not needed?
        #                  control=list(trace=ifelse(verbose,1L,0L), 
        #                               rel.tol=1e-5))
        #}
        if(out$convergence != 0L) warning("no convergence")
        rho <- tanh(out$par)
    }

    if(scores) {
        lik <- pcL_i(x.z1, x.z2, y.z1, y.z2, rho)
        dx <- ( dbinorm(x.z1, y.z1, rho) -
                dbinorm(x.z2, y.z1, rho) -
                dbinorm(x.z1, y.z2, rho) +
                dbinorm(x.z2, y.z2, rho) ) / lik    
        attr(rho, "scores") <- dx
    }

    rho
}

# polyserial x=numeric, y=ordinal -- TWO-STEP method
pscor_TS <- function(X, Y, th.y=NULL, verbose=FALSE, method="nlminb",
                     rescale.var=TRUE, acov=FALSE) {

    stopifnot(min(Y) == 1L, method %in% c("nlminb", "optimize"))

    # Y = ordinal
    if(is.null(th.y)) th.y <- unithord(Y)

    # X = numeric
    mu.x <- mean(X); var.x <- var(X);
    N <- length(X); if(rescale.var) var.x <- var.x * (N-1) / N
    Z <- (X - mu.x) / sqrt(var.x)

    # ML version to get Mplus results -- see Olsson 1982 eq 19 + 20
    objectiveFunction <- function(x) {
        rho = tanh(x[1L])
        logl <- psLogl(rho=rho, mu.x=mu.x, var.x=var.x, th.y=th.y, X=X, Y=Y,
                       Z=Z)
        -logl
    }

    if(method == "optimize") {
        out <- optimize(f=objectiveFunction,
                        interval=c(-5, +5),  # on tanh scale!
                        tol=1e-8, maximum=FALSE)
        rho <- out$minimum
    } else if(method == "nlminb") {
        rho.init <- cor(X,Y) * sd(Y) / sum(dnorm(th.y)) # Olsson 1982 eq 38

        gradientFunction <- function(x) {
            rho = tanh(x[1L]); TH.Y <- c(-Inf, th.y, Inf); R <- sqrt(1-rho^2)
            tauj.star  <- (TH.Y[Y+1L   ] - rho*Z)/R
            tauj1.star <- (TH.Y[Y+1L-1L] - rho*Z)/R
            d.tauj.star  <- dnorm(tauj.star); d.tauj1.star <- dnorm(tauj1.star)

            # p(y|x)
            pyx <- pspyx(rho=rho, th.y=th.y, Y=Y, Z=Z,
                         tauj.star=tauj.star, tauj1.star=tauj1.star)
            pyx.inv <- 1/pyx

            # rho
            TH.Y0 <- c(0,th.y,0) # if TH.Y = Inf or -Inf, dnorm() will be zero
            pyx.inv.c <- pyx.inv * 1/R^3
            TAUj  <- d.tauj.star  * (TH.Y0[Y+1L   ]*rho - Z)
            TAUj1 <- d.tauj1.star * (TH.Y0[Y+1L-1L]*rho - Z)
            dx.rho <- sum( pyx.inv.c * (TAUj - TAUj1) )

            dx.x <- dx.rho * 1/cosh(x)^2 # dF/drho * drho/dx 
            -dx.x
        }

        out <- nlminb(start=atanh(rho.init), objective=objectiveFunction,
                      gradient=gradientFunction,
                      scale=10,
                      control=list(trace=ifelse(verbose,1L,0L), rel.tol=1e-10))
        if(out$convergence != 0L) warning("no convergence")
        rho <- tanh(out$par)
    }

    rho
}

# X = residuals!
pscorX_TS <- function(X, Y, eXo=NULL, 
                      mu.x=NULL, var.x=NULL, eta.x=NULL,
                      y.z1=NULL, y.z2=NULL, 
                      verbose=FALSE, method="nlminb", scores=FALSE) {

    stopifnot(min(Y) == 1L, method %in% c("nlminb", "optimize"))

    # Y = ordinal
    if(is.null(y.z1) || is.null(y.z2)) {
        fit.y <- lavProbit(y=Y, X=eXo)
        y.z1 <- fit.y$z1; y.z2 <- fit.y$z2
    }

    # X = numeric
    if(is.null(mu.x) || is.null(var.x) || is.null(eta.x)) {
        fit.x <- lavOLS(y=X, X=eXo)
        var.x <- fit.x$theta[fit.x$npar]
        mu.x  <- fit.x$theta[1L]
        eta.x <- fit.x$yhat
    }
    ZRESID <- (X - eta.x)/sqrt(var.x)

    objectiveFunction <- function(x) {
        rho = tanh(x[1L])
        logl <- psLoglx(rho=rho, X=X, eta.x=eta.x, var.x=var.x, 
                        y.z1=y.z1, y.z2=y.z2, ZRESID=ZRESID)
        -logl
    }

    if(method == "optimize") {
        out <- optimize(f=objectiveFunction,
                        interval=c(-5, +5),  # on tanh scale!
                        tol=1e-8, maximum=FALSE)
        rho <- out$minimum
    } else if(method == "nlminb") {
        rho.init <- cor(X,Y)

        gradientFunction <- function(x) {
            rho = tanh(x[1L])
            R <- sqrt(1-rho^2)
            tauj.star  <- (y.z1 - rho*ZRESID)/R
            tauj1.star <- (y.z2 - rho*ZRESID)/R
            d.tauj.star  <- dnorm(tauj.star); d.tauj1.star <- dnorm(tauj1.star)
   
            # p(y|x)
            pyx <- pnorm(tauj.star) - pnorm(tauj1.star)
            pyx[pyx < .Machine$double.eps] <- .Machine$double.eps
            pyx.inv <- 1/pyx
     
            # rho
            pyx.inv.c <- pyx.inv * 1/R^3
            TAUj  <- d.tauj.star  * (y.z1*rho - ZRESID)
            TAUj1 <- d.tauj1.star * (y.z2*rho - ZRESID)
            dx.rho <- sum( pyx.inv.c * (TAUj - TAUj1) )
            -dx.rho * 1/cosh(x)^2
        }

        out <- nlminb(start=atanh(rho.init), objective=objectiveFunction,
                      gradient=gradientFunction,
                      scale=10,
                      control=list(trace=ifelse(verbose,1L,0L), rel.tol=1e-10))
        if(out$convergence != 0L) warning("no convergence")
        rho <- tanh(out$par)
    }

    if(scores) {
        R <- sqrt(1-rho^2)
        tauj.star  <- (y.z1 - rho*ZRESID)/R
        tauj1.star <- (y.z2 - rho*ZRESID)/R
        d.tauj.star  <- dnorm(tauj.star); d.tauj1.star <- dnorm(tauj1.star)

        # p(y|x)
        pyx <- pnorm(tauj.star) - pnorm(tauj1.star)
        pyx[pyx < .Machine$double.eps] <- .Machine$double.eps
        pyx.inv <- 1/pyx

        # rho
        pyx.inv.c <- pyx.inv * 1/R^3
        TAUj  <- d.tauj.star  * (y.z1*rho - ZRESID)
        TAUj1 <- d.tauj1.star * (y.z2*rho - ZRESID)
        dx.rho <-  pyx.inv.c * (TAUj - TAUj1)

        attr(rho, "scores") <- dx.rho
    }

    rho
}


# logl for X=numeric, Y=ordinal, (Z=scale(X), only for two.step)
psLogl <- function(rho, mu.x, var.x, th.y, X, Y, Z=NULL) {
    sd.x <- sqrt(var.x)
    if(is.null(Z)) Z <- (X - mu.x) / sd.x
    pyx <- pspyx(rho=rho, th.y=th.y, Y=Y, Z=Z)
    px <- dnorm(X, mean=mu.x, sd=sd.x)
    logl <- sum(log(px) + log(pyx))
    logl
}

psLoglx <- function(rho, X, eta.x, var.x, y.z1=NULL, y.z2=NULL, ZRESID=NULL) {
    R <- sqrt(1-rho^2)
    tauj.star  <- (y.z1 - rho*ZRESID)/R
    tauj1.star <- (y.z2 - rho*ZRESID)/R
    pyx <- pnorm(tauj.star) - pnorm(tauj1.star)
    pyx[pyx < .Machine$double.eps] <- .Machine$double.eps
    px <- dnorm(X, mean=eta.x, sd=sqrt(var.x))
    logl <- sum(log(px) + log(pyx))
    logl
}

# p(y|x)
pspyx <- function(rho, th.y, Y, Z, tauj.star=NULL, tauj1.star=NULL) {
    TH.Y <- c(-Inf, th.y, Inf); R <- sqrt(1-rho^2)
    if(is.null(tauj.star )) tauj.star  <- (TH.Y[Y+1L   ] - rho*Z)/R
    if(is.null(tauj1.star)) tauj1.star <- (TH.Y[Y+1L-1L] - rho*Z)/R
    pyx <- pnorm(tauj.star) - pnorm(tauj1.star)
    pyx[pyx < .Machine$double.eps] <- .Machine$double.eps
    pyx
}





unithord_logl_x <- function(x, X) {
    TH <- c(-Inf, x, Inf) 
    freq <- tabulate(X)
    PI <- pnorm(TH)[-1L] - pnorm(TH)[-length(TH)]
    logl <- sum( freq * log(PI) )
    logl
}

# testing only
# solution is: th = qnorm(cumsum(prop)[-length(prop)])
ml_th <- function(X, verbose=TRUE) {
    freq <- tabulate(X)
    objectiveFunction <- function(x) {
        th.x=x; logl <- unithord_logl_x(x, X)
        -logl
    }

    gradientFunction <- function(x) {
        TH <- c(-Inf, x, Inf)
        PI <- pnorm(TH)[-1L] - pnorm(TH)[-length(TH)]; D_pi.inv <- diag(1/PI)
        dth  <- dnorm(x)
        if(length(x) > 1L) {
            Delta <- matrix(0, length(x)+1L, length(x))
            diag(Delta) <- dth; diag(Delta[-1L,]) <- -dth
        } else {
            Delta <- matrix(c(dth, -dth), 2L, 1L)
        }
        dFdTH <- D_pi.inv %*% Delta
        # scores <- dFdTH[X,]; dx.th.x <- colSums(scores)
        dx.th.x <- freq %*% dFdTH
        -dx.th.x
    }

    freq <- tabulate(X)
    start <- seq(-1,+1,length.out=(length(freq)-1L))
    out <- nlminb(start=start, objective=objectiveFunction,
                  gradient=gradientFunction,
                  control=list(trace=ifelse(verbose,1L,0L), rel.tol=1e-10))
    if(out$convergence != 0L) warning("no convergence")
    out$par
}

unithord_logl_dx <- function(x, X) {
    TH <- c(-Inf, x, Inf)  
    freq <- tabulate(X)
    PI <- pnorm(TH)[-1L] - pnorm(TH)[-length(TH)]
    logl <- sum( freq * log(PI) )
    logl
}

scores_th <- function(X, th.x=NULL) {
    if(is.null(th.x)) th.x <- unithord(X)
    TH <- c(-Inf, th.x, Inf)
    PI <- pnorm(TH)[-1L] - pnorm(TH)[-length(TH)]; D_pi.inv <- diag(1/PI)
    dth  <- dnorm(th.x)
    if(length(th.x) > 1L) {
        Delta <- matrix(0, length(th.x)+1L, length(th.x))
        diag(Delta) <- dth; diag(Delta[-1L,]) <- -dth
    } else {
        Delta <- matrix(c(dth, -dth), 2L, 1L)
    }
    dFdTH <- D_pi.inv %*% Delta
    dFdTH[X,]   
}

scores_pscor <- function(X, Y, rho=NULL, mu.x=NULL, var.x=NULL, th.y=NULL) {
    sd.x=sqrt(var.x)
    Z <- (X - mu.x) / sd.x # ALWAYS recompute Z
    TH.Y <- c(-Inf, th.y, Inf); TH.Y0 <- c(0,th.y,0)

    R <- sqrt(1-rho^2)
    tauj.star  <- (TH.Y[Y+1L   ] - rho*Z)/R
    tauj1.star <- (TH.Y[Y+1L-1L] - rho*Z)/R
    d.tauj.star  <- dnorm(tauj.star); d.tauj1.star <- dnorm(tauj1.star)

    # p(y|x)
    pyx <- pspyx(rho=rho, th.y=th.y, Y=Y, Z=Z,
                 tauj.star=tauj.star, tauj1.star=tauj1.star)
    pyx.inv <- 1/pyx

    # rho
    pyx.inv.c <- pyx.inv * 1/R^3
    TAUj  <- d.tauj.star  * (TH.Y0[Y+1L   ]*rho - Z)
    TAUj1 <- d.tauj1.star * (TH.Y0[Y+1L-1L]*rho - Z)
    pyx.inv.c * (TAUj - TAUj1)
}

scores_pccor <- function(X, Y, rho=NULL, th.x=NULL, th.y=NULL) {
    PI  <- pcComputePI(rho, th.x, th.y)
    phi <- pcComputephi(rho, th.x, th.y)
    PP <- phi/PI
    PP[cbind(X,Y)]
}

scores_pccor_uni <- function(X, Y, rho, th.x, th.y) {

    PI  <- pcComputePI(rho, th.x, th.y); PI.XY.inv <- 1/PI[ cbind(X,Y) ]
    TH.X <- c(-Inf, th.x, Inf); TH.Y <- c(-Inf, th.y, Inf); R <- sqrt(1-rho^2)
 
    # th.x
    x.Y1 <- matrix(1:length(th.x),length(X),length(th.x), byrow=TRUE) == X
    x.Y2 <- matrix(1:length(th.x),length(X),length(th.x), byrow=TRUE) == (X-1L)
    x.z1 <- pmin( 100, TH.X[X+1L   ]); x.z2 <- pmax(-100, TH.X[X+1L-1L])
    x.Z1 <- ( dnorm(x.z1) * pnorm( (TH.Y[Y+1L   ]-rho*x.z1)/R) -
              dnorm(x.z1) * pnorm( (TH.Y[Y+1L-1L]-rho*x.z1)/R) )
    x.Z2 <- ( dnorm(x.z2) * pnorm( (TH.Y[Y+1L   ]-rho*x.z2)/R) -
              dnorm(x.z2) * pnorm( (TH.Y[Y+1L-1L]-rho*x.z2)/R) )
    dx.th.x <- (x.Y1*x.Z1 - x.Y2*x.Z2) * PI.XY.inv

    # th.y
    y.Y1 <- matrix(1:length(th.y),length(Y),length(th.y), byrow=TRUE) == Y
    y.Y2 <- matrix(1:length(th.y),length(Y),length(th.y), byrow=TRUE) == (Y-1L)
    y.z1 <- pmin( 100, TH.Y[Y+1L   ]); y.z2 <- pmax(-100, TH.Y[Y+1L-1L])
    y.Z1  <- ( dnorm(y.z1) * pnorm( (TH.X[X+1L   ]-rho*y.z1)/R) -
               dnorm(y.z1) * pnorm( (TH.X[X+1L-1L]-rho*y.z1)/R) )
    y.Z2  <- ( dnorm(y.z2) * pnorm( (TH.X[X+1L   ]-rho*y.z2)/R) -
               dnorm(y.z2) * pnorm( (TH.X[X+1L-1L]-rho*y.z2)/R) )
    dx.th.y <- (y.Y1*y.Z1 - y.Y2*y.Z2) * PI.XY.inv

    list(dx.th.x=dx.th.x, dx.th.y=dx.th.y)
}

# old version, not used anymore
scores_pccor_uni2 <- function(X, Y, rho, th.x, th.y) {

    PI  <- pcComputePI(rho, th.x, th.y); PI.XY.inv <- 1/PI[ cbind(X,Y) ]
    TH.X <- c(-Inf, th.x, Inf); TH.Y <- c(-Inf, th.y, Inf)
    dth.x <- dnorm(th.x); dth.y <- dnorm(th.y); R <- sqrt(1-rho^2) 

    dx.th.x <- matrix(0, length(X), length(th.x))
    ### FIXME: we can speed this up!
    for(k in 1:length(th.x)) {
        kj  <- dth.x[k] * pnorm((TH.Y[Y+1L   ]-rho*th.x[k])/R)
        kj1 <- dth.x[k] * pnorm((TH.Y[Y+1L-1L]-rho*th.x[k])/R)
        DpiDth <- ifelse(X == k, (kj - kj1), ifelse(X == k+1, (-kj + kj1), 0))
        dx.th.x[,k] <- PI.XY.inv * DpiDth
    }

    dx.th.y <- matrix(0, length(Y), length(th.y))
    for(m in 1:length(th.y)) {
        ki  <- dth.y[m] * pnorm((TH.X[X+1L   ]-rho*th.y[m])/R)
        ki1 <- dth.y[m] * pnorm((TH.X[X+1L-1L]-rho*th.y[m])/R)
        DpiDth <- ifelse(Y == m, (ki - ki1), ifelse(Y == m+1, (-ki + ki1), 0))
        dx.th.y[,m] <- PI.XY.inv * DpiDth
    }

    cbind(dx.th.x, dx.th.y)
}

scores_pccorX_uni <- function(X, Y, rho, th.x, th.y, sl.x, sl.y,
                              x.z1, x.z2, y.z1, y.z2, eXo=NULL) {

    lik <- pcL_i(x.z1=x.z1, x.z2=x.z2, y.z1=y.z1, y.z2=y.z2, rho=rho)
    R <- sqrt(1-rho^2)

    # th.x
    x.Y1 <- matrix(1:length(th.x),length(X),length(th.x), byrow=TRUE) == X
    x.Y2 <- matrix(1:length(th.x),length(X),length(th.x), byrow=TRUE) == (X-1L)
    x.Z1 <- ( dnorm(x.z1) * pnorm( (y.z1-rho*x.z1)/R) -
              dnorm(x.z1) * pnorm( (y.z2-rho*x.z1)/R) )
    x.Z2 <- ( dnorm(x.z2) * pnorm( (y.z1-rho*x.z2)/R) -
              dnorm(x.z2) * pnorm( (y.z2-rho*x.z2)/R) )
    dx.th.x <- (x.Y1*x.Z1 - x.Y2*x.Z2) / lik

    # th.y
    y.Y1 <- matrix(1:length(th.y),length(Y),length(th.y), byrow=TRUE) == Y
    y.Y2 <- matrix(1:length(th.y),length(Y),length(th.y), byrow=TRUE) == (Y-1L)
    y.Z1  <- ( dnorm(y.z1) * pnorm( (x.z1-rho*y.z1)/R) -
               dnorm(y.z1) * pnorm( (x.z2-rho*y.z1)/R) )
    y.Z2  <- ( dnorm(y.z2) * pnorm( (x.z1-rho*y.z2)/R) -
               dnorm(y.z2) * pnorm( (x.z2-rho*y.z2)/R) )
    dx.th.y <- (y.Y1*y.Z1 - y.Y2*y.Z2) / lik

    # sl.x
    dx.sl.x <- (x.Z2 - x.Z1) * eXo / lik

    # sl.y
    dx.sl.y <- (y.Z2 - y.Z1) * eXo / lik

    list(dx.th.x=dx.th.x, dx.th.y=dx.th.y, dx.sl.x=dx.sl.x, dx.sl.y=dx.sl.y)
}

scores_pscor_uni <- function(X, Y, rho=NULL, mu.x=NULL, var.x=NULL, th.y=NULL) {
    sd.x=sqrt(var.x); Z <- (X - mu.x) / sd.x # ALWAYS recompute Z
    TH.Y <- c(-Inf, th.y, Inf); R <- sqrt(1-rho^2)

    tauj.star  <- (TH.Y[Y+1L   ] - rho*Z)/R
    tauj1.star <- (TH.Y[Y+1L-1L] - rho*Z)/R
    y.Z1  <- dnorm(tauj.star); y.Z2 <- dnorm(tauj1.star)

    # p(y|x)
    pyx <- pnorm(tauj.star) - pnorm(tauj1.star)
    pyx[pyx < .Machine$double.eps] <- .Machine$double.eps
    pyx.inv <- 1/pyx

    # mu.x
    y.Z1.y.Z2 <- y.Z1-y.Z2
    dx.mu.x <- 1/sd.x * (Z + (pyx.inv * (rho/R) * y.Z1.y.Z2))

    # var.x
    dx.var.x <-  1/(2*var.x) * ( (Z^2-1) + (pyx.inv*rho*Z/R)*(y.Z1.y.Z2) )

    # th.y
    y.Y1 <- matrix(1:length(th.y),length(Y),length(th.y), byrow=TRUE) == Y
    y.Y2 <- matrix(1:length(th.y),length(Y),length(th.y), byrow=TRUE) == (Y-1L)
    dx.th.y <- (y.Y1*y.Z1 - y.Y2*y.Z2) * 1/R * pyx.inv 

    list(dx.mu.x=dx.mu.x, dx.var.x=dx.var.x, dx.th.y=dx.th.y)
}

# old version, not used anymore
scores_pscor_uni2 <- function(X, Y, rho=NULL, mu.x=NULL, var.x=NULL, th.y=NULL) {
    sd.x=sqrt(var.x)
    Z <- (X - mu.x) / sd.x # ALWAYS recompute Z
    TH.Y <- c(-Inf, th.y, Inf); TH.Y0 <- c(0,th.y,0)

    R <- sqrt(1-rho^2)
    tauj.star  <- (TH.Y[Y+1L   ] - rho*Z)/R
    tauj1.star <- (TH.Y[Y+1L-1L] - rho*Z)/R
    d.tauj.star  <- dnorm(tauj.star); d.tauj1.star <- dnorm(tauj1.star)

    # p(y|x)
    pyx <- pspyx(rho=rho, th.y=th.y, Y=Y, Z=Z,
                 tauj.star=tauj.star, tauj1.star=tauj1.star)
    pyx.inv <- 1/pyx

    # mu.x
    d.tauj.star.d.tauj1.star <- d.tauj.star-d.tauj1.star
    dx.mu.x <- 1/sd.x * (Z + (pyx.inv * (rho/R) * d.tauj.star.d.tauj1.star))

    # var.x
    dx.var.x <-  1/(2*var.x) * ( (Z^2-1) +
                           (pyx.inv*rho*Z/R)*(d.tauj.star.d.tauj1.star) )

    # th.y
    dx.th.y <- matrix(0, length(X), length(th.y))
    for(j in 1:length(th.y)) {
        dj  <- ifelse(Y == j,    d.tauj.star, 0)
        dj  <- ifelse(Y == j+1L, -d.tauj1.star, dj)
        dx.th.y[,j] <- 1/R * pyx.inv * dj 
    }

    cbind(dx.mu.x, dx.var.x, dx.th.y)
}

scores_pscorX_uni <- function(X, Y, eXo=NULL, rho=NULL,
                              eta.x=NULL, var.x=NULL, y.z1=NULL, y.z2=NULL) {

    ZRESID <- (X - eta.x)/sqrt(var.x)
    R <- sqrt(1-rho^2)
    tauj.star  <- (y.z1 - rho*ZRESID)/R
    tauj1.star <- (y.z2 - rho*ZRESID)/R
    y.Z1  <- dnorm(tauj.star); y.Z2 <- dnorm(tauj1.star)

    # p(y|x)
    pyx <- pnorm(tauj.star) - pnorm(tauj1.star)
    pyx[pyx < .Machine$double.eps] <- .Machine$double.eps
    pyx.inv <- 1/pyx

    # mu.x
    y.Z1.y.Z2 <- y.Z1-y.Z2
    dx.mu.x <- 1/sqrt(var.x) * (ZRESID + (pyx.inv * (rho/R) * y.Z1.y.Z2))

    # var.x
    dx.var.x <-  1/(2*var.x) * ( (ZRESID^2-1) +
                     (pyx.inv*rho*ZRESID/R)*(y.Z1.y.Z2) )

    # th.y
    nth <- NCOL(y.z1)
    y.Y1 <- matrix(1:nth, length(Y), nth, byrow=TRUE) == Y
    y.Y2 <- matrix(1:nth, length(Y), nth, byrow=TRUE) == (Y-1L)
    dx.th.y <- (y.Y1*y.Z1 - y.Y2*y.Z2) * 1/R * pyx.inv

    # sl.x
    
    dx.sl.x <-  1/var.x * eXo * (X - eta.x) # * p(X) * 1/p(X)

    # sl.y
    dx.sl.y <- (y.Z2 - y.Z1) * eXo * pyx.inv

    list(dx.mu.x=dx.mu.x, dx.var.x=dx.var.x, dx.th.y=dx.th.y, 
         dx.sl.x=dx.sl.x, dx.sl.y=dx.sl.y)
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
    list(dx.mu.x=dx.mu.x, dx.mu.y=dx.mu.y, dx.var.x=dx.var.x, dx.var.y=dx.var.y)
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


# faster!
pcLogl_freq <- function(freq, rho, th.x, th.y) {

    PI <- pcComputePI(rho, th.x, th.y)
    logl <- sum( freq * log(PI) )

    logl
}

pcComputePI <- function(rho, th.x, th.y) {

    TH.X <- c(-Inf, th.x, Inf); TH.Y <- c(-Inf, th.y, Inf)
    pTH.X <- pnorm(TH.X); rowPI <- pTH.X[-1L] - pTH.X[-length(TH.X)]
    pTH.Y <- pnorm(TH.Y); colPI <- pTH.Y[-1L] - pTH.Y[-length(TH.Y)]

    # catch special case: rho = 0.0
    if(rho == 0.0) {
        PI.ij <- outer(rowPI, colPI)
        return(PI.ij)
    }

    nr <- length(TH.X) - 1L; nc <- length(TH.Y) - 1L
    sigma <- matrix(c(1,rho,rho,1), 2L, 2L)
    PI <- matrix(0, nr, nc)
    for(i in seq_len(nr-1L)) {
        for(j in seq_len(nc-1L)) {
            PI[i,j] <- pmvnorm(lower=c(TH.X[i-1L+1L], TH.Y[j-1L+1L]),
                                  upper=c(TH.X[i   +1L], TH.Y[j   +1L]),
                                  mean=c(0,0), corr=sigma)
        }
    }
    # add last col (rowSums(PI) must correspond with TH.X)
    PI[,nc] <- rowPI - rowSums(PI[,1:(nc-1L),drop=FALSE])
    # PI[nr,nc] will be wrong at this point, but gets overridden
    # add last row (colSums(PI) must correspond with TH.Y)
    PI[nr,] <- colPI - colSums(PI[1:(nr-1L),,drop=FALSE])

    # all elements should be strictly positive
    PI[PI < .Machine$double.eps] <- .Machine$double.eps

    PI
}

pcComputephi <- function(rho, th.x, th.y) {

    TH.X <- c(-Inf, th.x, Inf); TH.Y <- c(-Inf, th.y, Inf)
    nr <- length(TH.X) - 1L; nc <- length(TH.Y) - 1L
    phi <- matrix(0, nr, nc)
    for(i in seq_len(nr)) {
        for(j in seq_len(nc)) {
            p1 <- p2 <- p3 <- p4 <- 0
            if(i < nr && j < nc)
                p1 <- dbinorm(TH.X[i   +1L], TH.Y[j   +1L], rho)
            if(i > 1L && j < nc)
                p2 <- dbinorm(TH.X[i-1L+1L], TH.Y[j   +1L], rho)
            if(i < nr && j > 1L)
                p3 <- dbinorm(TH.X[i   +1L], TH.Y[j-1L+1L], rho)
            if(i > 1L && j > 1L)
                p4 <- dbinorm(TH.X[i-1L+1L], TH.Y[j-1L+1L], rho)
            phi[i,j] <- (p1 - p2 - p3 + p4)
        }
    }

    phi
}

pcComputegnorm <- function(rho, th.x, th.y) {

    # note: Olsson 1979 A2 contains an error!!
    guv <- function(u, v, rho) {
        R <- (1-rho^2)
        ( u*v*R - rho*(u^2 - 2*rho*u*v + v^2) + rho*R ) / R^2
    }

    TH.X <- c(-Inf, th.x, Inf); TH.Y <- c(-Inf, th.y, Inf)
    nr <- length(TH.X) - 1L; nc <- length(TH.Y) - 1L
    gnorm <- matrix(0, nr, nc)
    for(i in seq_len(nr)) {
        for(j in seq_len(nc)) {
            g1 <- g2 <- g3 <- g4 <- 0
            if(i < nr && j < nc) {
                u <- TH.X[i   +1L]; v <- TH.Y[j   +1L]
                g1 <- dbinorm(u, v, rho) * guv(u,v,rho)
            }
            if(i > 1L && j < nc) {
                u <- TH.X[i-1L+1L]; v <- TH.Y[j   +1L]
                g2 <- dbinorm(u, v, rho) * guv(u,v,rho)
            }
            if(i < nr && j > 1L) {
                u <- TH.X[i   +1L]; v <- TH.Y[j-1L+1L]
                g3 <- dbinorm(u, v, rho) * guv(u,v,rho)
            }
            if(i > 1L && j > 1L) {
                u <- TH.X[i-1L+1L]; v <- TH.Y[j-1L+1L]
                g4 <- dbinorm(u, v, rho) * guv(u,v,rho)
            }
            gnorm[i,j] <- (g1 - g2 - g3 + g4)
        }
    }

    gnorm
}

# density of a bivariate standard normal 
dbinorm <- function(u, v, rho) {
    R <- 1-rho^2
    1/(2*pi*sqrt(R)) * exp( - 0.5*(u^2 - 2*rho*u*v + v^2)/R )
}

# partial derivative - rho
dbinorm_drho <- function(u, v, rho) {
    R <- 1 - rho^2
    dbinorm(u,v,rho) * (u*v*R -rho*(u^2 - 2*rho*u*v + v^2) + rho*R )/R^2
}

# partial derivative - u
dbinorm_du <- function(u, v, rho) {
    R <- 1 - rho^2
    -dbinorm(u,v,rho) * (u - rho*v)/R
}

# partial derivative - v
dbinorm_dv <- function(u, v, rho) {
    R <- 1 - rho^2
    -dbinorm(u,v,rho) * (v - rho*u)/R
}

# CDF of bivariate standard normal
# function pbinorm(upper.x, upper.y, rho)

# partial derivative pbinorm - upper.x
pbinorm_dupper.x <- function(upper.x, upper.y, rho=0.0) {
    R <- 1 - rho^2
    dnorm(upper.x) * pnorm( (upper.y - rho*upper.x)/R )
}

pbinorm_dupper.y <- function(upper.x, upper.y, rho=0.0) {
    R <- 1 - rho^2
    dnorm(upper.y) * pnorm( (upper.x - rho*upper.y)/R )
}

pbinorm_drho <- function(upper.x, upper.y, rho=0.0) {
    dbinorm(upper.x, upper.y, rho)    
}



pcL_i <- function(x.z1, x.z2, y.z1, y.z2, rho) {
   # using pbivnorm is MUCH faster (loop in fortran)
   #
   #lik <-  ( pbivnorm(x=x.z1, y=y.z1, rho=rho) -
   #          pbivnorm(x=x.z2, y=y.z1, rho=rho) -
   #          pbivnorm(x=x.z1, y=y.z2, rho=rho) +
   #          pbivnorm(x=x.z2, y=y.z2, rho=rho)  )

   # this uses mvtnorm
   lik <- pbinorm(upper.x=x.z1, upper.y=y.z1,
                  lower.x=x.z2, lower.y=y.z2, rho=rho)
   lik
}


pcLogl_i <- function(x.z1, x.z2, y.z1, y.z2, rho) {

    lik <- pcL_i(x.z1=x.z1, x.z2=x.z2, y.z1=y.z1, y.z2=y.z2, rho=rho)
    logl <- sum( log(lik) )

    logl
}

# testing only
F_ij_pc <- function(x) {
    rho = x[1L]
    th.x = x[1L + 1:length(th.x)]
    th.y = x[1L + length(th.x) + 1:length(th.y)]
    
    logl <- pcLogl_freq(freq, rho=tanh(x[1L]), th.x, th.y)
    logl
}
#numDeriv:::grad(F_ij_pc, x=x.par)
#scores <- scores_pccor_uni(X=X, Y=Y, rho=rho, th.x=th.x, th.y=th.y)
#apply(scores, 2, sum)

F_ij_ps <- function(x) {
    rho = x[1L]
    mu.x = x[2L]
    var.x = x[3L]
    th.y = x[3L + 1:length(th.y)]

    Z <- (X - mu.x) / sqrt(var.x)
    logl <- psLogl(rho=rho, mu.x=mu.x, var.x=var.x, th.y=th.y, X=X, Y=Y,
                   Z=Z)
    logl
}

#numDeriv:::grad(F_ij_ps, x=x.par)
#scores <- scores_pscor_uni(X=X, Y=Y, rho=rho, mu.x=mu.x, var.x=var.x, th.y=th.y)
#apply(scores, 2, sum)

F_ij_pcX <- function(x) {
    rho = x[1L]
mu.x = x[2L]
    var.x = x[3L]
    th.y = x[3L + 1:length(th.y)]



    th.x = x[1L + 1:length(th.x)]
    th.y = x[1L + length(th.x) + 1:length(th.y)]
    sl.x = x[1L + length(th.x) + length(th.y) + 1:length(sl.x)]
    sl.y = x[1L + length(th.x) + length(th.y) + length(sl.x) + 1:length(sl.y)]

    TH.x <- c(-Inf, th.x, +Inf); eta.x <- drop(eXo %*% sl.x)
    x.z1 <- pmin( 100, TH.x[X+1L   ] - eta.x)
    x.z2 <- pmax(-100, TH.x[X+1L-1L] - eta.x)
    TH.y <- c(-Inf, th.y, +Inf); eta.y <- drop(eXo %*% sl.y)
    y.z1 <- pmin( 100, TH.y[Y+1L   ] - eta.y)
    y.z2 <- pmax(-100, TH.y[Y+1L-1L] - eta.y)

    logl <- pcLogl_i(x.z1, x.z2, y.z1, y.z2, rho)
    logl
}
#numDeriv:::grad(F_ij_pcX, x=x.par)[-1]

F_ij_pcX <- function(x) {
    rho = x[1L]
    mu.x = x[2L]
    var.x = x[3L]
    th.y = x[3L + 1:length(th.y)]

    th.x = x[1L + 1:length(th.x)]
    th.y = x[1L + length(th.x) + 1:length(th.y)]
    sl.x = x[1L + length(th.x) + length(th.y) + 1:length(sl.x)]
    sl.y = x[1L + length(th.x) + length(th.y) + length(sl.x) + 1:length(sl.y)]

    TH.x <- c(-Inf, th.x, +Inf); eta.x <- drop(eXo %*% sl.x)
    x.z1 <- pmin( 100, TH.x[X+1L   ] - eta.x)
    x.z2 <- pmax(-100, TH.x[X+1L-1L] - eta.x)
    TH.y <- c(-Inf, th.y, +Inf); eta.y <- drop(eXo %*% sl.y)
    y.z1 <- pmin( 100, TH.y[Y+1L   ] - eta.y)
    y.z2 <- pmax(-100, TH.y[Y+1L-1L] - eta.y)

    logl <- pcLogl_i(x.z1, x.z2, y.z1, y.z2, rho)
    logl
}
#numDeriv:::grad(F_ij_pcX, x=x.par)[-1]


F_ij_psX <- function(x) {

    rho = x[1L]
    mu.x = x[2L]
    var.x = x[3L]
    th.y = x[3L + 1:length(th.y)]
    sl.x = x[3L + length(th.y) + 1:length(sl.x)]
    sl.y = x[3L + length(th.y) + length(sl.x) + 1:length(sl.y)]

    eta.x <- drop( cbind(1,eXo) %*% c(mu.x,sl.x) )
    ZRESID <- (fit.x$y - eta.x)/sqrt(var.x)

    TH.y <- c(-Inf, th.y, +Inf); eta.y <- drop(eXo %*% sl.y)
    y.z1 <- pmin( 100, TH.y[Y+1L   ] - eta.y)
    y.z2 <- pmax(-100, TH.y[Y+1L-1L] - eta.y)

    logl <- psLoglx(rho=rho, X=X, eta.x=eta.x, var.x=var.x, y.z1=y.z1,
                    y.z2=y.z2, ZRESID=ZRESID)
    logl
}
#numDeriv:::grad(F_ij_psX, x=x.par)[-1]

F_ij_corX <- function(x) {

    Lx <- function(X, Y, rho, eta.x, var.x, eta.y, var.y) {
        cov.xy <- rho*sqrt(var.x)*sqrt(var.y)
        sigma <- matrix(c(var.x,cov.xy,cov.xy,var.y), 2L, 2L)
        dx <- numeric(length(X))
        for(i in 1:length(X))
            dx[i] <- dmvnorm(c(X[i],Y[i]), mean=c(eta.x[i], eta.y[i]), sigma=sigma)
        sum(log(dx))
    }

    rho = x[1L]
    mu.x = x[2L]
    var.x = x[3L]
    mu.y = x[4L]
    var.y = x[5L]
    sl.x = x[5L + 1:length(sl.x)]
    sl.y = x[5L + length(sl.x) + 1:length(sl.y)]

    eta.x <- drop( cbind(1,eXo) %*% c(mu.x,sl.x) )
    eta.y <- drop( cbind(1,eXo) %*% c(mu.y,sl.y) )

    logl <- Lx(X=X, Y=Y, rho=rho, eta.x=eta.x, var.x=var.x,
               eta.y=eta.y, var.y=var.y)
    logl
}
#numDeriv:::grad(F_ij_corX, x=x.par)[-1]




