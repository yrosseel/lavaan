muthen1984 <- function(Data, ov.names=NULL, ov.types=NULL) {

    require(mvtnorm)

    # This function was written in January 2012 -- Yves Rosseel
    # First success: Friday 20 Jan 2012: the standard errors for
    #                thresholds and polychoric correlations (in an 
    #                unrestricted/saturated model) are spot on!

    nvar <- ncol(Data); stopifnot(nvar > 1L); N <- nrow(Data)

    # variable types; default = numeric
    ord.idx <- which(ov.types == "ordered")
    num.idx <- which(ov.types == "numeric"); nnum <- length(num.idx)

    # means and thresholds
    TH <- vector("list", length=nvar)
    # slopes (only if fixed.x)
    SLOPES <- vector("list", length=nvar)
    # variances (for continuous variables only)
    VAR <- vector("list", length=nvar)
    # correlations
    COR <- diag(nvar); colnames(COR) <- rownames(COR) <- ov.names

    # stage one - TH/VAR only
    for(i in 1:nvar) {
        if(ov.types[i] == "numeric") {
            # compute mean and variance
            TH[[i]] <- mean(Data[,i], na.rm=TRUE)
            VAR[[i]] <- var(Data[,i], na.rm=TRUE) * (N-1)/N
        } else if(ov.types[i] == "ordered") {
            TH[[i]] <- unithord(X=Data[,i]); VAR[[i]] <- 1.0
        } else {
            stop("unknown ov.types:", ov.types[i])
        }
    }

    # stage one - SC
    SC.TH  <- matrix(0, N, length(unlist(TH)))
    colnames(SC.TH) <- unlist(lapply(as.list(1:nvar), 
        function(x) paste(ov.names[x],"|",1:length(TH[[x]]),sep="")))
    SC.SL  <- matrix(0, N, nvar)
    SC.VAR <- matrix(0, N, nvar); colnames(SC.VAR) <- ov.names
    nTH <- unlist(lapply(TH, length)) # number of TH per variable
    end.idx <- cumsum(nTH); start.idx <- end.idx - (nTH - 1L)
    for(i in 1:nvar) {
        th.idx <- start.idx[i]:end.idx[i]
        if(ov.types[i] == "numeric") {
            SC.TH[,th.idx] <-
                scores_mu(Data[,i], mu.x=TH[[i]], var.x=VAR[[i]])
            SC.VAR[,i] <- scores_var(Data[,i], mu.x=TH[[i]],var.x=VAR[[i]])
        } else if(ov.types[i] == "ordered") {
            SC.TH[,th.idx] <- scores_th(Data[,i], TH[[i]])
            SC.VAR[,i] <- rep(0, N)
        } else {
            stop("unknown ov.types:", ov.types[i])
        }
    }
    # rm VAR columns from ordinal variables
    SC.VAR <- SC.VAR[,-ord.idx] 


    # stage two
    pstar <- nvar*(nvar-1)/2
    SC.COR <- matrix(0, N, pstar)
    PSTAR <- matrix(0, nvar, nvar)
    COR.NAMES <- character(pstar)
    # LAVAAN style: col-wise!
    PSTAR[lavaan:::vech.idx(nvar, diag=FALSE)] <- 1:pstar
    # LISREL style: row-wise
    #PSTAR[lavaan:::vechr.idx(nvar, diag=FALSE)] <- 1:pstar
    for(j in 1:(nvar-1L)) {
        for(i in (j+1L):nvar) {
            pstar.idx <- PSTAR[i,j]
            COR.NAMES[pstar.idx] <- paste(ov.names[i],"~~",ov.names[j],sep="")
            if(ov.types[i] == "numeric" && ov.types[j] == "numeric") {
                COR[i,j] <- COR[j,i] <- cor(Data[,i], Data[,j], 
                                            use="pairwise.complete.obs")
                SC.COR[,pstar.idx] <-
                    scores_cor(X=Data[,i], Y=Data[,j], rho=COR[i,j],
                               mu.x=TH[[i]], var.x=VAR[[i]], 
                               mu.y=TH[[j]], var.y=VAR[[j]])
            } else if(ov.types[i] == "numeric" && ov.types[j] == "ordered") {
                # polyserial
                out <- pscor_TS(X=Data[,i], Y=Data[,j], th.y=TH[[j]])
                COR[i,j] <- COR[j,i] <- out
                SC.COR[,pstar.idx] <- 
                    scores_pscor(X=Data[,i], Y=Data[,j], rho=COR[i,j], 
                                 mu.x=TH[[i]], var.x=VAR[[i]], th.y=TH[[j]])
            } else if(ov.types[j] == "numeric" && ov.types[i] == "ordered") {
                # polyserial
                out <- pscor_TS(X=Data[,j], Y=Data[,i], th.y=TH[[i]])
                COR[i,j] <- COR[j,i] <- out
                SC.COR[,pstar.idx] <- 
                    scores_pscor(X=Data[,j], Y=Data[,i], rho=COR[i,j], 
                                 mu.x=TH[[j]], var.x=VAR[[j]], th.y=TH[[i]])
            } else if(ov.types[i] == "ordered" && ov.types[j] == "ordered") {
                # polychoric correlation
                out <- pccor_TS(Data[,i], Data[,j],
                                th.x = TH[[i]], th.y = TH[[j]])
                COR[i,j] <- COR[j,i] <- out
                SC.COR[,pstar.idx] <- 
                    scores_pccor(X=Data[,i], Y=Data[,j], rho=COR[i,j],
                                 th.x=TH[[i]], th.y=TH[[j]])
            }
        }
    }
    colnames(SC.COR) <- COR.NAMES

    # stage three
    SC <- cbind(SC.TH, SC.VAR, SC.COR)
    INNER <- crossprod(SC)

    # A11
    A11_block <- vector("list", length=(nvar + nnum))
    # A11 - TH
    for(i in 1:nvar) {
        th.idx <- start.idx[i]:end.idx[i]
        SC.TH_i <- SC.TH[,th.idx]
        A11_block[[i]] <- crossprod(SC.TH_i)
    }
    # A11 - VAR
    if(nnum > 0L) {
        for(i in 1:ncol(SC.VAR)) {
            A11_block[[ nvar + i ]] <- crossprod(SC.VAR[,i])
        }
    }
    A11 <- lavaan:::bdiag(A11_block)

    # A21
    A21 <- matrix(0, pstar, ncol(A11))
    # for this one, we need new scores: for each F_ij (cor), the
    # scores with respect to the TH, VAR, ...
    for(j in 1:(nvar-1L)) {
        for(i in (j+1L):nvar) {
            pstar.idx <- PSTAR[i,j]
            th.idx_i <- start.idx[i]:end.idx[i]
            th.idx_j <- start.idx[j]:end.idx[j]
            if(ov.types[i] == "numeric" && ov.types[j] == "numeric") {
                SC.COR.UNI <-
                    scores_cor_uni(X=Data[,i], Y=Data[,j], rho=COR[i,j],
                                   mu.x=TH[[i]], var.x=VAR[[i]], 
                                   mu.y=TH[[j]], var.y=VAR[[j]])
                A21[pstar.idx, th.idx_i] <-
                    crossprod(SC.COR[,pstar.idx],
                              SC.COR.UNI[,1L])
                A21[pstar.idx, th.idx_j] <-
                    crossprod(SC.COR[,pstar.idx],
                              SC.COR.UNI[,2L])
                A21[pstar.idx, ncol(SC.TH) + match(i, num.idx)] <-
                    crossprod(SC.COR[,pstar.idx],
                              SC.COR.UNI[,3L])
                A21[pstar.idx, ncol(SC.TH) + match(j, num.idx)] <-
                    crossprod(SC.COR[,pstar.idx],
                              SC.COR.UNI[,4L])
            } else if(ov.types[i] == "numeric" && ov.types[j] == "ordered") {
                # polyserial
                SC.COR.UNI <-
                    scores_pscor_uni(X=Data[,i], Y=Data[,j], rho=COR[i,j], 
                                     mu.x=TH[[i]], var.x=VAR[[i]], th.y=TH[[j]])
                A21[pstar.idx, th.idx_i] <-
                    crossprod(SC.COR[,pstar.idx],
                              SC.COR.UNI[,1L])
                A21[pstar.idx, th.idx_j] <-
                    crossprod(SC.COR[,pstar.idx],
                              SC.COR.UNI[,2L+(1:length(TH[[j]]))])
                A21[pstar.idx, ncol(SC.TH) + match(i, num.idx)] <-
                    crossprod(SC.COR[,pstar.idx],
                              SC.COR.UNI[,2L])
            } else if(ov.types[j] == "numeric" && ov.types[i] == "ordered") {
                # polyserial
                SC.COR.UNI <-
                    scores_pscor_uni(X=Data[,j], Y=Data[,i], rho=COR[i,j], 
                                     mu.x=TH[[j]], var.x=VAR[[j]], th.y=TH[[i]])
                A21[pstar.idx, th.idx_j] <-
                    crossprod(SC.COR[,pstar.idx],
                              SC.COR.UNI[,1L])
                A21[pstar.idx, th.idx_i] <-
                    crossprod(SC.COR[,pstar.idx],
                              SC.COR.UNI[,2L+(1:length(TH[[i]]))])
                A21[pstar.idx, ncol(SC.TH) + match(j, num.idx)] <-
                    crossprod(SC.COR[,pstar.idx],
                              SC.COR.UNI[,2L])
            } else if(ov.types[i] == "ordered" && ov.types[j] == "ordered") {
                # polychoric correlation
                SC.COR.UNI <-
                    scores_pccor_uni(X=Data[,i], Y=Data[,j], rho=COR[i,j],
                                     th.x=TH[[i]], th.y=TH[[j]])
                A21[pstar.idx, th.idx_i] <- 
                    crossprod(SC.COR[,pstar.idx], 
                              SC.COR.UNI[,1:length(TH[[i]])])
                A21[pstar.idx, th.idx_j] <- 
                    crossprod(SC.COR[,pstar.idx],
                              SC.COR.UNI[,length(TH[[i]])+1:length(TH[[j]])])
            }
        }
    }
    

    # A22
    A22 <- matrix(0, pstar, pstar)
    for(i in 1:pstar) {
        A22[i,i] <- sum( SC.COR[,i]^2 )
    }

    # A12
    A12 <- matrix(0, nrow(A11), ncol(A22))

    B <- rbind( cbind(A11,A12),
                cbind(A21,A22) )
    B.inv <- solve(B)

    ACOV <- (B.inv %*% INNER %*% t(B.inv)) * N

    # COV matrix
    COV <- cor2cov(R=COR, sds=sqrt(unlist(VAR)))

    out <- list(TH=TH, SLOPES=SLOPES, VAR=VAR, COR=COR, COV=COV,
                INNER=INNER, A11=A11, A12=A12, A21=A21, A22=A22,
                ACOV=ACOV)
    out
}


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
              method %in% c("nlminb", "optimize"))

    # thresholds
    if(is.null(th.x)) th.x <- unithord_freq(tabulate(x))
    if(is.null(th.y)) th.y <- unithord_freq(tabulate(y))

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
    } else if(method == "nlminb") {

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

        out <- nlminb(start=atanh(rho.init), objective=objectiveFunction,
                      gradient=gradientFunction,
                      scale=10,
                      control=list(trace=ifelse(verbose,1L,0L), rel.tol=1e-10))
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

# logl for X=numeric, Y=ordinal, (Z=scale(X), only for two.step)
psLogl <- function(rho, mu.x, var.x, th.y, X, Y, Z=NULL) {
    sd.x <- sqrt(var.x)
    if(is.null(Z)) Z <- (X - mu.x) / sd.x
    pyx <- pspyx(rho=rho, th.y=th.y, Y=Y, Z=Z)
    px <- dnorm(X, mean=mu.x, sd=sd.x)
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

scores_pscor_uni <- function(X, Y, rho=NULL, mu.x=NULL, var.x=NULL, th.y=NULL) {
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
    cbind(dx.mu.x, dx.mu.y, dx.var.x, dx.var.y)
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
    # add last col (rowSums(PI) must correspond witlh TH.X)
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
    sigma <- matrix(c(1,rho,rho,1), 2L, 2L)
    phi <- matrix(0, nr, nc)
    for(i in seq_len(nr)) {
        for(j in seq_len(nc)) {
            p1 <- p2 <- p3 <- p4 <- 0
            if(i < nr && j < nc)
                p1 <- dmvnorm(c(TH.X[i   +1L], TH.Y[j   +1L]), sigma=sigma)
            if(i > 1L && j < nc)
                p2 <- dmvnorm(c(TH.X[i-1L+1L], TH.Y[j   +1L]), sigma=sigma)
            if(i < nr && j > 1L)
                p3 <- dmvnorm(c(TH.X[i   +1L], TH.Y[j-1L+1L]), sigma=sigma)
            if(i > 1L && j > 1L)
                p4 <- dmvnorm(c(TH.X[i-1L+1L], TH.Y[j-1L+1L]), sigma=sigma)
            phi[i,j] <- (p1 - p2 - p3 + p4)
        }
    }

    phi
}

