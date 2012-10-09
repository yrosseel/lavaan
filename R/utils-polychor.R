# polychoric Y1 and Y2 are both ORDINAL variables

# two-way frequency table
pc_freq <- function(Y1, Y2) {
    max.y1 <- max(Y1, na.rm=TRUE); max.y2 <- max(Y2, na.rm=TRUE)
    bin <- Y1 - 1L; bin <- bin + max.y1 * (Y2 - 1L); bin <- bin[!is.na(bin)]
    if (length(bin)) bin <- bin + 1L
    array(tabulate(bin, nbins = max.y1*max.y2), dim=c(max.y1, max.y2))
}

# compute thresholds (no eXo)
pc_th <- function(Y, freq=NULL, prop=NULL) {
    if(is.null(prop)) {
        if(is.null(freq)) freq <- tabulate(Y)
        prop <- freq / sum(freq)
    }
    prop1 <- prop[-length(prop)]
    qnorm(cumsum(prop1))
}

pc_PI <- function(rho, th.y1, th.y2) {

    TH.Y1 <- c(-Inf, th.y1, Inf); TH.Y2 <- c(-Inf, th.y2, Inf)
    pTH.Y1 <- pnorm(TH.Y1); rowPI <- pTH.Y1[-1L] - pTH.Y1[-length(TH.Y1)]
    pTH.Y2 <- pnorm(TH.Y2); colPI <- pTH.Y2[-1L] - pTH.Y2[-length(TH.Y2)]

    # catch special case: rho = 0.0
    if(rho == 0.0) {
        PI.ij <- outer(rowPI, colPI)
        return(PI.ij)
    }

    nr <- length(TH.Y1) - 1L; nc <- length(TH.Y2) - 1L
    PI <- matrix(0, nr, nc)
    for(i in seq_len(nr-1L)) {
        for(j in seq_len(nc-1L)) {
            PI[i,j] <- pbinorm(lower.x=TH.Y1[i-1L+1L], lower.y=TH.Y2[j-1L+1L],
                               upper.x=TH.Y1[i   +1L], upper.y=TH.Y2[j   +1L],
                               rho)
        }
    }
    # add last col (rowSums(PI) must correspond with TH.Y1)
    PI[,nc] <- rowPI - rowSums(PI[,1:(nc-1L),drop=FALSE])
    # PI[nr,nc] will be wrong at this point, but gets overridden
    # add last row (colSums(PI) must correspond with TH.Y2)
    PI[nr,] <- colPI - colSums(PI[1:(nr-1L),,drop=FALSE])

    # all elements should be strictly positive
    PI[PI < .Machine$double.eps] <- .Machine$double.eps

    PI
}

pc_PHI <- function(rho, th.y1, th.y2) {

    TH.Y1 <- c(-Inf, th.y1, Inf); TH.Y2 <- c(-Inf, th.y2, Inf)
    nr <- length(TH.Y1) - 1L; nc <- length(TH.Y2) - 1L
    phi <- matrix(0, nr, nc)
    for(i in seq_len(nr)) {
        for(j in seq_len(nc)) {
            p1 <- p2 <- p3 <- p4 <- 0
            if(i < nr && j < nc)
                p1 <- dbinorm(TH.Y1[i   +1L], TH.Y2[j   +1L], rho)
            if(i > 1L && j < nc)
                p2 <- dbinorm(TH.Y1[i-1L+1L], TH.Y2[j   +1L], rho)
            if(i < nr && j > 1L)
                p3 <- dbinorm(TH.Y1[i   +1L], TH.Y2[j-1L+1L], rho)
            if(i > 1L && j > 1L)
                p4 <- dbinorm(TH.Y1[i-1L+1L], TH.Y2[j-1L+1L], rho)
            phi[i,j] <- (p1 - p2 - p3 + p4)
        }
    }

    phi
}

pc_gnorm <- function(rho, th.y1, th.y2) {

    # note: Olsson 1979 A2 contains an error!!
    guv <- function(u, v, rho) {
        R <- (1-rho^2)
        ( u*v*R - rho*(u^2 - 2*rho*u*v + v^2) + rho*R ) / R^2
    }

    TH.Y1 <- c(-Inf, th.y1, Inf); TH.Y2 <- c(-Inf, th.y2, Inf)
    nr <- length(TH.Y1) - 1L; nc <- length(TH.Y2) - 1L
    gnorm <- matrix(0, nr, nc)
    for(i in seq_len(nr)) {
        for(j in seq_len(nc)) {
            g1 <- g2 <- g3 <- g4 <- 0
            if(i < nr && j < nc) {
                u <- TH.Y1[i   +1L]; v <- TH.Y2[j   +1L]
                g1 <- dbinorm(u, v, rho) * guv(u,v,rho)
            }
            if(i > 1L && j < nc) {
                u <- TH.Y1[i-1L+1L]; v <- TH.Y2[j   +1L]
                g2 <- dbinorm(u, v, rho) * guv(u,v,rho)
            }
            if(i < nr && j > 1L) {
                u <- TH.Y1[i   +1L]; v <- TH.Y2[j-1L+1L]
                g3 <- dbinorm(u, v, rho) * guv(u,v,rho)
            }
            if(i > 1L && j > 1L) {
                u <- TH.Y1[i-1L+1L]; v <- TH.Y2[j-1L+1L]
                g4 <- dbinorm(u, v, rho) * guv(u,v,rho)
            }
            gnorm[i,j] <- (g1 - g2 - g3 + g4)
        }
    }

    gnorm
}

# (summed) loglikelihood
pc_logl <- function(Y1, Y2, eXo=NULL, rho=NULL, fit.y1=NULL, fit.y2=NULL,
                    freq=NULL) {

    stopifnot(!is.null(rho))

    if(is.null(fit.y1)) fit.y1 <- lavProbit(y=Y1, X=eXo)
    if(is.null(fit.y2)) fit.y2 <- lavProbit(y=Y2, X=eXo)

    # if no eXo, use shortcut (grouped)
    if(length(fit.y1$slope.idx) == 0L) {
        if(is.null(freq)) freq <- pc_freq(fit.y1$y,fit.y2$y)
        # grouped lik
        PI <- pc_PI(rho, th.y1=fit.y1$theta[fit.y1$th.idx], 
                         th.y2=fit.y2$theta[fit.y2$th.idx])
        if(all(PI > 0)) 
            logl <- sum( freq * log(PI) )
        else logl <- -Inf
    } else {
        lik <- pc_lik(Y1=Y1, Y2=Y2, rho=rho, fit.y1=fit.y1, fit.y2=fit.y2)
        if(all(lik > 0))
            logl <- sum( log(lik) )
        else logl <- -Inf 
    }

    logl
}

# individual likelihoods
pc_lik <- function(Y1, Y2, eXo=NULL, rho=NULL, fit.y1=NULL, fit.y2=NULL) {

    stopifnot(!is.null(rho))

    if(is.null(fit.y1)) fit.y1 <- lavProbit(y=Y1, X=eXo)
    if(is.null(fit.y2)) fit.y2 <- lavProbit(y=Y2, X=eXo)
    
    # if no eXo, use shortcut (grouped)
    if(length(fit.y1$slope.idx) == 0L) {
        # probability per cell
        PI <- pc_PI(rho, th.y1=fit.y1$theta[fit.y1$th.idx],
                         th.y2=fit.y2$theta[fit.y2$th.idx])
        lik <- PI[ cbind(fit.y1$y, fit.y2$y) ]
    } else {
        # individual likelihoods
        # pbivnorm package is MUCH faster (loop in fortran)
        # lik <-  ( pbivnorm(x=fit.y1$z1, y=fit.y2$z1, rho=rho) -
        #           pbivnorm(x=fit.y1$z2, y=fit.y2$z1, rho=rho) -
        #           pbivnorm(x=fit.y1$z1, y=fit.y2$z2, rho=rho) +
        #           pbivnorm(x=fit.y1$z2, y=fit.y2$z2, rho=rho)  )

        # this uses mvtnorm
        lik <- pbinorm(upper.x=fit.y1$z1, upper.y=fit.y2$z1,
                       lower.x=fit.y1$z2, lower.y=fit.y2$z2, rho=rho)
    }

    lik
}

# loglikelihood (x-version)
pc_logl_x <- function(x, Y1, Y2, eXo=NULL, nth.y1, nth.y2, freq=NULL) {

    nexo <- ifelse(is.null(eXo), 0L, ncol(eXo)); S <- seq_len
    stopifnot(length(x) == (1L + nth.y1 + nth.y2 + 2*nexo))

      rho = x[1L]
    th.y1 = x[1L + S(nth.y1)]
    th.y2 = x[1L +   nth.y1 + S(nth.y2)]
    sl.y1 = x[1L +   nth.y1 +   nth.y2 + S(nexo)]
    sl.y2 = x[1L +   nth.y1 +   nth.y2 +   nexo + S(nexo)]

    fit.y1 <- lavProbit(y=Y1, X=eXo)
    fit.y1$theta[fit.y1$th.idx] <- th.y1
    fit.y1$theta[fit.y1$slope.idx] <- sl.y1
    fit.y1$lik()

    fit.y2 <- lavProbit(y=Y2, X=eXo)
    fit.y2$theta[fit.y2$th.idx] <- th.y2
    fit.y2$theta[fit.y2$slope.idx] <- sl.y2
    fit.y2$lik()

    pc_logl(Y1=Y1, Y2=Y2, eXo=eXo, rho=rho, fit.y1=fit.y1, fit.y2=fit.y2,
            freq=freq)
}

# polychoric correlation
pc_cor_TS <- function(Y1, Y2, eXo=NULL, fit.y1=NULL, fit.y2=NULL, freq=NULL,
                      method="nlminb", zerofreq=0.5, verbose=FALSE) {

    if(is.null(fit.y1)) fit.y1 <- lavProbit(y=Y1, X=eXo)
    if(is.null(fit.y2)) fit.y2 <- lavProbit(y=Y2, X=eXo)
    if(missing(Y1)) Y1 <- fit.y1$y else as.integer(Y1)
    if(missing(Y2)) Y2 <- fit.y2$y else as.integer(Y2)
    if(missing(eXo) && length(fit.y1$slope.idx) > 0L) eXo <- fit.y1$X

    stopifnot(min(Y1) == 1L, min(Y2) == 1L,
              method %in% c("nlminb", "nlminb.hessian"))

    # exo or not?
    exo <- ifelse(length(fit.y1$slope.idx) > 0L, TRUE, FALSE)

    # thresholds
    th.y1 <- fit.y1$theta[fit.y1$th.idx]
    th.y2 <- fit.y2$theta[fit.y2$th.idx]
    
    # freq
    if(!exo) {
        if(is.null(freq)) freq <- pc_freq(fit.y1$y,fit.y2$y)
        nr <- nrow(freq); nc <- ncol(freq)
        # check for empty cells -- FIXME: make this an option!
        if(any(freq == 0) && zerofreq > 0) {
            #freq[freq == 0] <- zerofreq/length(fit.y1$y)
            freq[freq == 0] <- zerofreq
        }
    }

    objectiveFunction <- function(x) {
        logl <- pc_logl(rho=tanh(x[1L]), fit.y1=fit.y1, fit.y2=fit.y2,
                        freq=freq)
        -logl # to minimize!
    }

    gradientFunction <- function(x) {
        rho <- tanh(x[1L])
        if(!exo) {
            PI  <- pc_PI(rho, th.y1, th.y2)
            phi <- pc_PHI(rho, th.y1, th.y2)
            dx.rho <- sum(freq/PI * phi)
        } else {
            lik <- pc_lik(rho=rho, fit.y1=fit.y1, fit.y2=fit.y2)
            dx <- ( dbinorm(fit.y1$z1, fit.y2$z1, rho) -
                    dbinorm(fit.y1$z2, fit.y2$z1, rho) -
                    dbinorm(fit.y1$z1, fit.y2$z2, rho) +
                    dbinorm(fit.y1$z2, fit.y2$z2, rho) ) / lik
            dx.rho <- sum(dx)
        }
        -dx.rho * 1/cosh(x)^2 # dF/drho * drho/dx, dtanh = 1/cosh(x)^2
    }

    hessianFunction2 <- function(x) {
        numDeriv:::hessian(func=objectiveFunction, x=x)
    }

    # OLSSON 1979 A2 + A3 (no EXO!!)
    hessianFunction <- function(x) {
        rho <- tanh(x[1L])
        PI  <- pc_PI(rho, th.y1, th.y2)
        phi <- pc_PHI(rho, th.y1, th.y2)
        gnorm <- pc_gnorm(rho, th.y1, th.y2)
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

    # starting value

    # catch tetrachoric case
    #if(!exo && (nr == 2L && nc == 2L) && !any(freq == 0)) {
        # Divgi 1979 initial value
    #    h <- max(abs(th.y1), abs(th.y2)); k <- min(abs(th.y1), abs(th.y2))
        # h can not be zero;
    #    if(h == 0) h <- 1e-5
    #    R <- (freq[1,1]*freq[2,2])/(freq[1,2]*freq[2,1])
    #    D <- k*(.79289 + 4.28981/(1+3.30231*h));D <- D*sign(th.y1)*sign(th.y2)
    #    C <- 0.07557*h + (h-k)^2 * (0.51141/(h+2.05793) - 0.07557/h)
    #    B <- 0.5/(1 + (h^2 + k^2)*(0.82281-1.03514*(k/sqrt(h^2+k^2))))
    #    A <- 0.5/(1 + (h^2 + k^2)*(0.12454-0.27102*(1-h/sqrt(h^2+k^2))))
    #    alpha <- A + B*(-1 + 1/(1 + C*(log(R)-D)^2))
    #    rho.init <- cos(pi/(1+R^alpha))
    #} else {
        rho.init <- cor(Y1,Y2)
    #}

    if(method == "nlminb") {
        out <- nlminb(start=atanh(rho.init), objective=objectiveFunction,
                      gradient=gradientFunction,
                      scale=10,
                      control=list(trace=ifelse(verbose,1L,0L),
                                   rel.tol=1e-10))
    } else if(method == "nlminb.hessian") {
        stopifnot(!exo)
        out <- nlminb(start=atanh(rho.init), objective=objectiveFunction,
                      gradient=gradientFunction,
                      hessian=hessianFunction2,
                      scale=100, # not needed?
                      control=list(trace=ifelse(verbose,1L,0L),
                                   rel.tol=1e-7))
    }
    if(out$convergence != 0L) warning("no convergence")
    rho <- tanh(out$par)

    rho
}

pc_cor_scores <- function(Y1, Y2, eXo=NULL, rho, fit.y1=NULL, fit.y2=NULL,
                          th.y1=NULL, th.y2=NULL,
                          sl.y1=NULL, sl.y2=NULL) {

    R <- sqrt(1-rho^2)
    if(is.null(fit.y1)) fit.y1 <- lavProbit(y=Y1, X=eXo)
    if(is.null(fit.y2)) fit.y2 <- lavProbit(y=Y2, X=eXo)
    y1.update <- y2.update <- FALSE
    if(!is.null(th.y1)) { # update thresholds fit.y1
        y1.update <- TRUE
        fit.y1$theta[fit.y1$th.idx] <- th.y1
    }
    if(!is.null(th.y2)) { # update thresholds fit.y1
        y2.update <- TRUE
        fit.y2$theta[fit.y2$th.idx] <- th.y2
    }
    if(!is.null(sl.y1)) { # update slopes
        y1.update <- TRUE
        fit.y1$theta[fit.y1$slope.idx] <- sl.y1
    }
    if(!is.null(sl.y2)) { # update slopes
        y2.update <- TRUE
        fit.y2$theta[fit.y2$slope.idx] <- sl.y2
    }
    if(y1.update) tmp <- fit.y1$lik()
    if(y2.update) tmp <- fit.y2$lik()
    if(missing(Y1)) Y1 <- fit.y1$y
    if(missing(Y2)) Y2 <- fit.y2$y
    if(missing(eXo) && length(fit.y1$slope.idx) > 0L) eXo <- fit.y1$X
   
    # lik
    lik <- pc_lik(rho=rho, fit.y1=fit.y1, fit.y2=fit.y2)

    # th.y1
    y1.Z1 <- ( dnorm(fit.y1$z1) * pnorm( (fit.y2$z1-rho*fit.y1$z1)/R) -
               dnorm(fit.y1$z1) * pnorm( (fit.y2$z2-rho*fit.y1$z1)/R) )
    y1.Z2 <- ( dnorm(fit.y1$z2) * pnorm( (fit.y2$z1-rho*fit.y1$z2)/R) -
               dnorm(fit.y1$z2) * pnorm( (fit.y2$z2-rho*fit.y1$z2)/R) )
    dx.th.y1 <- (fit.y1$Y1*y1.Z1 - fit.y1$Y2*y1.Z2) / lik

    # th.y2
    y2.Z1  <- ( dnorm(fit.y2$z1) * pnorm( (fit.y1$z1-rho*fit.y2$z1)/R) -
                dnorm(fit.y2$z1) * pnorm( (fit.y1$z2-rho*fit.y2$z1)/R) )
    y2.Z2  <- ( dnorm(fit.y2$z2) * pnorm( (fit.y1$z1-rho*fit.y2$z2)/R) -
                dnorm(fit.y2$z2) * pnorm( (fit.y1$z2-rho*fit.y2$z2)/R) )
    dx.th.y2 <- (fit.y2$Y1*y2.Z1 - fit.y2$Y2*y2.Z2) / lik

    dx.sl.y1 <- dx.sl.y2 <- NULL
    if(length(fit.y1$slope.idx) > 0L) {
        # sl.y1
        dx.sl.y1 <- (y1.Z2 - y1.Z1) * eXo / lik
        # sl.y2
        dx.sl.y2 <- (y2.Z2 - y2.Z1) * eXo / lik
    }

    # rho
    if(length(fit.y1$slope.idx) == 0L) {
        phi <- pc_PHI(rho, th.y1=fit.y1$theta[fit.y1$th.idx], 
                           th.y2=fit.y2$theta[fit.y2$th.idx])
        #PP <- phi/PI
        dx <- phi[cbind(Y1,Y2)]
    } else {
        dx <- ( dbinorm(fit.y1$z1, fit.y2$z1, rho) -
                dbinorm(fit.y1$z2, fit.y2$z1, rho) -
                dbinorm(fit.y1$z1, fit.y2$z2, rho) +
                dbinorm(fit.y1$z2, fit.y2$z2, rho) )
    }
    dx.rho <- dx / lik

    list(dx.th.y1=dx.th.y1, dx.th.y2=dx.th.y2, 
         dx.sl.y1=dx.sl.y1, dx.sl.y2=dx.sl.y2, dx.rho=dx.rho)
}
