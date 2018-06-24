# polychoric Y1 and Y2 are both ORDINAL variables

# two-way frequency table
pc_freq <- function(Y1, Y2) {
    # FIXME: for 2x2, we could use
    # array(tabulate((Y1-1) + (Y2-1)*2 + 1, nbins = 4L), dim=c(2L ,2L))
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

    nth.y1 <- length(th.y1); nth.y2 <- length(th.y2)
    pth.y1 <- pnorm(th.y1);  pth.y2 <- pnorm(th.y2)

    # catch special case: rho = 0.0
    if(rho == 0.0) {
        rowPI <- diff(c(0,pth.y1,1))
        colPI <- diff(c(0,pth.y2,1))
        PI.ij <- outer(rowPI, colPI)
        return(PI.ij)
    }

    # prepare for a single call to pbinorm
    upper.y <- rep(th.y2, times=rep.int(nth.y1, nth.y2))
    upper.x <- rep(th.y1, times=ceiling(length(upper.y))/nth.y1)
    #rho <- rep(rho, length(upper.x)) # only one rho here

    BI <- pbivnorm::pbivnorm(x=upper.x, y=upper.y, rho=rho)
    #BI <- pbinorm1(upper.x=upper.x, upper.y=upper.y, rho=rho)
    dim(BI) <- c(nth.y1, nth.y2)
    BI <- rbind(0, BI, pth.y2, deparse.level = 0)
    BI <- cbind(0, BI, c(0, pth.y1, 1), deparse.level = 0)



    # get probabilities
    nr <- nrow(BI); nc <- ncol(BI)
    PI <- BI[-1L,-1L] - BI[-1L,-nc] - BI[-nr,-1L] + BI[-nr,-nc]

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
        R <- (1 - rho*rho)
        ( u*v*R - rho*((u*u) - 2*rho*u*v + (v*v)) + rho*R ) / (R*R)
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
        if(all(lik > 0, na.rm = TRUE))
            logl <- sum( log(lik), na.rm = TRUE )
        else logl <- -Inf
    }

    logl
}

# pc_lik, with user-specified parameters
pc_lik2 <- function(Y1, Y2, eXo=NULL, rho, fit.y1=NULL, fit.y2=NULL,
                    th.y1=NULL, th.y2=NULL,
                    sl.y1=NULL, sl.y2=NULL) {

    R <- sqrt(1 - rho*rho)
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

    lik
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

        # take care of missing values
        if(fit.y1$missing.values || fit.y2$missing.values) {
            lik <- rep(as.numeric(NA), length(fit.y1$z1))
            missing.idx <- unique(c(fit.y1$missing.idx,
                                    fit.y2$missing.idx))
            fit.y1.z1 <- fit.y1$z1; fit.y1.z1[missing.idx] <- 0
            fit.y2.z1 <- fit.y2$z1; fit.y2.z1[missing.idx] <- 0
            fit.y1.z2 <- fit.y1$z2; fit.y1.z2[missing.idx] <- 0
            fit.y2.z2 <- fit.y2$z2; fit.y2.z2[missing.idx] <- 0

            lik <- pbinorm(upper.x=fit.y1.z1, upper.y=fit.y2.z1,
                           lower.x=fit.y1.z2, lower.y=fit.y2.z2, rho=rho)
            lik[missing.idx] <- NA
        } else {
            lik <- pbinorm(upper.x=fit.y1$z1, upper.y=fit.y2$z1,
                           lower.x=fit.y1$z2, lower.y=fit.y2$z2, rho=rho)
        }

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
#
# zero.add is a vector: first element is for 2x2 tables only, second element
#                       for general tables
# zero.keep.margins is only used for 2x2 tables
pc_cor_TS <- function(Y1, Y2, eXo=NULL, fit.y1=NULL, fit.y2=NULL, freq=NULL,
                      method="nlminb", zero.add = c(0.5, 0.0), control=list(),
                      zero.keep.margins = TRUE, zero.cell.warn = FALSE,
                      zero.cell.flag = FALSE,
                      verbose=FALSE, Y1.name=NULL, Y2.name=NULL) {

    # cat("DEBUG: method = ", method, "\n")

    if(is.null(fit.y1)) fit.y1 <- lavProbit(y=Y1, X=eXo)
    if(is.null(fit.y2)) fit.y2 <- lavProbit(y=Y2, X=eXo)
    if(missing(Y1)) Y1 <- fit.y1$y else as.integer(Y1)
    if(missing(Y2)) Y2 <- fit.y2$y else as.integer(Y2)
    if(missing(eXo) && length(fit.y1$slope.idx) > 0L) eXo <- fit.y1$X

    stopifnot(min(Y1, na.rm=TRUE) == 1L, min(Y2, na.rm=TRUE) == 1L,
              method %in% c("nlminb", "BFGS", "nlminb.hessian", "optimize"))

    # empty cells or not
    empty.cells <- FALSE

    # exo or not?
    exo <- ifelse(length(fit.y1$slope.idx) > 0L, TRUE, FALSE)

    # thresholds
    th.y1 <- fit.y1$theta[fit.y1$th.idx]
    th.y2 <- fit.y2$theta[fit.y2$th.idx]

    # freq
    if(!exo) {
        if(is.null(freq)) freq <- pc_freq(fit.y1$y,fit.y2$y)
        nr <- nrow(freq); nc <- ncol(freq)

        # check for empty cells
        if(any(freq == 0L)) {
            empty.cells <- TRUE
            if(zero.cell.warn) {
                if(!is.null(Y1.name) && !is.null(Y2.name)) {
                warning("lavaan WARNING: empty cell(s) in bivariate table of ",
                            Y1.name, " x ", Y2.name)
                } else {
                    warning("lavaan WARNING: empty cell(s) in bivariate table")
                }
            }
        }

        # treat 2x2 tables
        if(nr == 2L && nc == 2L) {
            idx <- which(freq == 0L)
            # catch 2 empty cells: perfect correlation!
            if(length(idx) == 2L) {
                warning("lavaan WARNING: two empty cells in 2x2 table")
                if(freq[1,1] > 0L) {
                    rho <- 1.0
                    if(zero.cell.flag) {
                        attr(rho, "zero.cell.flag") <- empty.cells
                    }
                    return(rho)
                } else {
                    rho <- -1.0
                    if(zero.cell.flag) {
                        attr(rho, "zero.cell.flag") <- empty.cells
                    }
                    return(rho)
                }
            } else if(length(idx) == 1L && zero.add[1] > 0.0) {
                if(zero.keep.margins) {
                    # add + compensate to preserve margins
                    if(idx == 1L || idx == 4L) { # main diagonal
                        freq[1,1] <- freq[1,1] + zero.add[1]
                        freq[2,2] <- freq[2,2] + zero.add[1]
                        freq[2,1] <- freq[2,1] - zero.add[1]
                        freq[1,2] <- freq[1,2] - zero.add[1]
                    } else {
                        freq[1,1] <- freq[1,1] - zero.add[1]
                        freq[2,2] <- freq[2,2] - zero.add[1]
                        freq[2,1] <- freq[2,1] + zero.add[1]
                        freq[1,2] <- freq[1,2] + zero.add[1]
                    }
                } else {
                    freq[idx] <- freq[idx] + zero.add[1]
                }
            }
        # general table
        } else {
            if(any(freq == 0L) && zero.add[2] > 0.0) {
                # general table: just add zero.add to the empty cell(s)
                freq[freq == 0] <- zero.add[2]
            }
        }

        # catch special cases for 2x2 tables
        if(nr == 2L && nc == 2L) {
            # 1. a*d == c*d
            storage.mode(freq) <- "numeric" # to avoid integer overflow
            if(freq[1,1]*freq[2,2] == freq[1,2]*freq[2,1]) {
                rho <- 0.0
                if(zero.cell.flag) {
                    attr(rho, "zero.cell.flag") <- empty.cells
                }
                return(rho)
            }
            # 2. equal margins (th1 = th2 = 0)
            if(th.y1[1] == 0L && th.y2[1] == 0L) {
                # see eg Brown & Benedetti 1977 eq 2
                rho <- - cos( 2*pi*freq[1,1]/sum(freq) )
                if(zero.cell.flag) {
                    attr(rho, "zero.cell.flag") <- empty.cells
                }
                return(rho)
            }
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
            dx.rho <- sum(dx, na.rm = TRUE)
        }
        -dx.rho * 1/(cosh(x)*cosh(x)) # dF/drho * drho/dx, dtanh = 1/cosh(x)^2
    }

    #hessianFunction2 <- function(x) {
    #    numDeriv::hessian(func=objectiveFunction, x=x)
    #}

    # OLSSON 1979 A2 + A3 (no EXO!!)
    hessianFunction <- function(x) {
        rho <- tanh(x[1L])
        PI  <- pc_PI(rho, th.y1, th.y2)
        phi <- pc_PHI(rho, th.y1, th.y2)
        gnorm <- pc_gnorm(rho, th.y1, th.y2)
        H <-  sum(freq/PI * gnorm) - sum(freq/(PI*PI) * (phi*phi))

        # to compensate for tanh
        # u=f(x), d^2y/dx^2 = d^2y/du^2 * (du/dx)^2 + dy/du * d^2u/dx^2
        # dtanh = 1/cosh(x)^2
        # dtanh_2 = 8*exp(2*x)*(1-exp(2*x))/(exp(2*x)+1)^3
        grad <- sum(freq/PI * phi)
        u1 <- 1/(cosh(x)*cosh(x))
        tmp3 <- (exp(2*x)+1) * (exp(2*x)+1) * (exp(2*x)+1)
        u2 <- 8*exp(2*x)*(1-exp(2*x))/tmp3
        H <- H * (u1*u1) + grad * u2
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
        rho.init <- cor(Y1,Y2, use="pairwise.complete.obs")
    #}

    # if rho.init is missing, set starting value to zero
    if(is.na(rho.init)) {
      rho.init <- 0.0
    }
    # check range of rho.init is within [-1,+1]
    if(abs(rho.init) >= 1.0) {
        rho.init <- 0.0
    }

    # default values
    control.nlminb <- list(eval.max=20000L,
                               iter.max=10000L,
                               trace=ifelse(verbose, 1L, 0L),
                               #abs.tol=1e-20, ### important!! fx never negative
                               abs.tol=(.Machine$double.eps * 10),
                               rel.tol=ifelse(method == "nlminb", 1e-10, 1e-7),
                               x.tol=1.5e-8,
                               xf.tol=2.2e-14)
    control.nlminb <- modifyList(control.nlminb, control)
    control <- control.nlminb[c("eval.max", "iter.max", "trace",
                                "abs.tol", "rel.tol", "x.tol", "xf.tol")]

    if(method == "nlminb") {
        out <- nlminb(start=atanh(rho.init), objective=objectiveFunction,
                      gradient=gradientFunction,
                      scale=10,
                      control=control)
        if(out$convergence != 0L) warning("no convergence")
        rho <- tanh(out$par)
    } else if(method == "BFGS") {
        # NOTE: known to fail if rho.init is too far from final value
        # seems to be better with parscale = 0.1??
        out <- optim(par = atanh(rho.init), fn = objectiveFunction,
                     gr = gradientFunction,
                     control = list(parscale = 0.1, reltol = 1e-10,
                                    trace = ifelse(verbose, 1L, 0L),
                                    REPORT = 1L,
                                    abstol=(.Machine$double.eps * 10)),
                     method = "BFGS")
        if(out$convergence != 0L) warning("no convergence")
        rho <- tanh(out$par)
    } else if(method == "optimize") {
        # not atanh/tanh transform
        objectiveFunction2 <- function(x) {
            logl <- pc_logl(rho=x[1L], fit.y1=fit.y1, fit.y2=fit.y2,
                            freq=freq)
            -logl # to minimize!
        }
        out <- optimize(f = objectiveFunction2, interval = c(-0.9999,0.9999))
        rho <- out$minimum
    } else if(method == "nlminb.hessian") {
        stopifnot(!exo)
        out <- nlminb(start=atanh(rho.init), objective=objectiveFunction,
                      gradient=gradientFunction,
                      hessian=hessianFunction,
                      scale=100, # not needed?
                      control=control)
        if(out$convergence != 0L) warning("no convergence")
        rho <- tanh(out$par)
    }

    if(zero.cell.flag) {
        attr(rho, "zero.cell.flag") <- empty.cells
    }

    rho
}

pc_cor_gradient_noexo <- function(Y1, Y2, rho, th.y1=NULL, th.y2=NULL,
                                  freq=NULL) {

    R <- sqrt(1- rho*rho)
    TH.Y1 <- c(-Inf, th.y1, Inf); TH.Y2 <- c(-Inf, th.y2, Inf)
    dth.y1 <- dnorm(th.y1); dth.y2 <- dnorm(th.y2)
    if(is.null(freq)) freq <- pc_freq(Y1, Y2)

    # rho
    PI  <- pc_PI(rho, th.y1, th.y2)
    phi <- pc_PHI(rho, th.y1, th.y2)
    dx.rho <- sum(freq/PI * phi)

    # th.y2
    PI.XY.inv <- 1/PI[ cbind(Y1,Y2) ]
    dx.th.y2 <- matrix(0, length(Y2), length(th.y2))
    for(m in 1:length(th.y2)) {
        ki  <- dth.y2[m] * pnorm((TH.Y1[Y1+1L   ]-rho*th.y2[m])/R)
        ki1 <- dth.y2[m] * pnorm((TH.Y1[Y1+1L-1L]-rho*th.y2[m])/R)
        DpiDth <- ifelse(Y2 == m, (ki - ki1), ifelse(Y2 == m+1, (-ki + ki1), 0))
        dx.th.y2[,m] <- PI.XY.inv * DpiDth
    }

}

pc_cor_scores <- function(Y1, Y2, eXo=NULL, rho, fit.y1=NULL, fit.y2=NULL,
                          th.y1=NULL, th.y2=NULL,
                          sl.y1=NULL, sl.y2=NULL,
                          na.zero=FALSE) {

    # check if rho >

    R <- sqrt(1 - rho*rho)
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
    if(identical(R, 0.0)) {
        y1.Z1 <- dnorm(fit.y1$z1) * 0.5
        y1.Z2 <- dnorm(fit.y1$z2) * 0.5
    } else {
        y1.Z1 <- ( dnorm(fit.y1$z1) * pnorm( (fit.y2$z1-rho*fit.y1$z1)/R) -
                   dnorm(fit.y1$z1) * pnorm( (fit.y2$z2-rho*fit.y1$z1)/R) )
        y1.Z2 <- ( dnorm(fit.y1$z2) * pnorm( (fit.y2$z1-rho*fit.y1$z2)/R) -
                   dnorm(fit.y1$z2) * pnorm( (fit.y2$z2-rho*fit.y1$z2)/R) )
    }
    dx.th.y1 <- (fit.y1$Y1*y1.Z1 - fit.y1$Y2*y1.Z2) / lik
    if(na.zero) {
        dx.th.y1[is.na(dx.th.y1)] <- 0
    }

    # th.y2
    if(identical(R, 0.0)) {
        y2.Z1  <- dnorm(fit.y2$z1) * 0.5
        y2.Z2  <- dnorm(fit.y2$z2) * 0.5
    } else {
        y2.Z1  <- ( dnorm(fit.y2$z1) * pnorm( (fit.y1$z1-rho*fit.y2$z1)/R) -
                    dnorm(fit.y2$z1) * pnorm( (fit.y1$z2-rho*fit.y2$z1)/R) )
        y2.Z2  <- ( dnorm(fit.y2$z2) * pnorm( (fit.y1$z1-rho*fit.y2$z2)/R) -
                    dnorm(fit.y2$z2) * pnorm( (fit.y1$z2-rho*fit.y2$z2)/R) )
    }
    dx.th.y2 <- (fit.y2$Y1*y2.Z1 - fit.y2$Y2*y2.Z2) / lik
    if(na.zero) {
        dx.th.y2[is.na(dx.th.y2)] <- 0
    }

    dx.sl.y1 <- dx.sl.y2 <- NULL
    if(length(fit.y1$slope.idx) > 0L) {
        # sl.y1
        dx.sl.y1 <- (y1.Z2 - y1.Z1) * eXo / lik
        if(na.zero) {
            dx.sl.y1[is.na(dx.sl.y1)] <- 0
        }
        # sl.y2
        dx.sl.y2 <- (y2.Z2 - y2.Z1) * eXo / lik
        if(na.zero) {
            dx.sl.y2[is.na(dx.sl.y2)] <- 0
        }
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
    if(na.zero) {
        dx.rho[is.na(dx.rho)] <- 0
    }

    list(dx.th.y1=dx.th.y1, dx.th.y2=dx.th.y2,
         dx.sl.y1=dx.sl.y1, dx.sl.y2=dx.sl.y2, dx.rho=dx.rho)
}
