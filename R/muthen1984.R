muthen1984 <- function(Data, ov.names=NULL, ov.types=NULL, ov.levels=NULL,
                       ov.names.x=character(0L), eXo=NULL, verbose=FALSE) {

    require(mvtnorm)

    # This function was written in January 2012 -- Yves Rosseel
    # First success: Friday 20 Jan 2012: the standard errors for
    #                thresholds and polychoric correlations (in an 
    #                unrestricted/saturated model) are spot on!
    # Second success: Saturday 9 June 2012: support for mixed (ordinal + metric)
    #                 variables; thanks to the delta method to get the ACOV 
    #                 right (see H matrix)
    # Third success: Monday 2 July 2012: support for fixed.x covariates

    nvar <- ncol(Data); N <- nrow(Data)
    nTH <- ov.levels - 1L; nTH[nTH == -1L] <- 1L
    nth <- sum(nTH)
    th.end.idx <- cumsum(nTH); th.start.idx <- th.end.idx - (nTH - 1L)

    # variable types; default = numeric
    ord.idx <- which(ov.types == "ordered")
    num.idx <- which(ov.types == "numeric"); nnum <- length(num.idx)
    nexo <- length(ov.names.x)
    if(nexo > 0L) stopifnot(ncol(eXo) == nexo)

    if(verbose) {
        cat("number of endo variables: ", nvar, "\n")
        cat("endo variable names:\n"); print(ov.names); cat("\n")
        cat("endo ov types:\n"); print(ov.types); cat("\n")
        cat("endo ov levels:\n "); print(ov.levels); cat("\n")
        cat("number of exo variables: ", nexo, "\n")
        cat(" exo variable names:\n"); print(ov.names.x); cat("\n")
    }

    # means and thresholds
    TH <- vector("list", length=nvar)
    TH.NAMES <- vector("list", length=nvar)
    TH.IDX <- vector("list", length=nvar)
    # slopes (only if fixed.x)
    SLOPES <- vector("list", length=nvar)
    # variances (for continuous variables only)
    VAR <- vector("list", length=nvar)
    # correlations
    COR <- diag(nvar); colnames(COR) <- rownames(COR) <- ov.names

    # SCORES
    SC.VAR <- matrix(0, N, nvar); colnames(SC.VAR) <- ov.names
    SC.SL  <- matrix(0, N, nvar*nexo)
    colnames(SC.SL) <- paste(rep(ov.names, each=nexo), 
                             rep(ov.names.x,nvar), sep="")
    SC.TH  <- matrix(0, N, nth)
    colnames(SC.TH) <- unlist(lapply(as.list(1:nvar),
        function(x) paste(ov.names[x],"|",1:nTH[x],sep="")))
    FIT <- vector("list", length=nvar)

    # stage one - TH/SLOPES/VAR only
    ov.num <- 0L
    for(i in 1:nvar) {
        #cat("var = ", ov.names[i], " type = ", ov.types[i], "\n")
        th.idx <- th.start.idx[i]:th.end.idx[i]
        if(ov.types[i] == "numeric") {
            fit <- lavOLS(y=Data[,i], X=eXo); scores <- fit$scores()
            FIT[[i]] <- fit
            ov.num <- ov.num + 1L
            # compute mean and variance
            TH[[i]] <- fit$theta[1L]
            VAR[[i]] <- fit$theta[fit$npar]
            TH.NAMES[[i]] <- ov.names[i]; TH.IDX[[i]] <- 0L
            #SC.TH[,th.idx] <-
            #    scores_mu(Data[,i], mu.x=TH[[i]], var.x=VAR[[i]])
            SC.TH[,th.idx] <- scores[,1L]
            #SC.VAR[,i] <- scores_var(Data[,i], mu.x=TH[[i]],var.x=VAR[[i]])
            SC.VAR[,i] <- scores[,fit$npar]
            if(nexo > 0L) {
                SLOPES[[i]] <- fit$theta[-c(1L, fit$npar)]
                sl.end.idx <- (i*nexo); sl.start.idx <- (i-1L)*nexo + 1L
                SC.SL[,sl.start.idx:sl.end.idx] <- 
                    scores[,-c(1L, fit$npar),drop=FALSE]
            }
        } else if(ov.types[i] == "ordered") {
            if(nexo == 0L) {
                # FIXME: merge with lavProbit...
                TH[[i]] <- unithord(X=Data[,i])
                SC.TH[,th.idx] <- scores_th(Data[,i], TH[[i]])
                SC.VAR[,i] <- rep(0, N)
            } else {
                fit <- lavProbit(y=Data[,i], X=eXo); scores <- fit$scores()
                FIT[[i]] <- fit
                TH[[i]] <- fit$theta[fit$th.idx]
                SC.TH[,th.idx] <- scores[,fit$th.idx,drop=FALSE]
                SLOPES[[i]] <- fit$theta[fit$beta.idx]
                sl.end.idx <- (i*nexo); sl.start.idx <- (i-1L)*nexo + 1L
                SC.SL[,sl.start.idx:sl.end.idx] <- 
                    scores[,fit$beta.idx,drop=FALSE]
            }
            VAR[[i]] <- 1.0
            TH.NAMES[[i]] <- paste(ov.names[i], "|t", 1:length(TH[[i]]), 
                                   sep="")
            TH.IDX[[i]] <- rep(i, length(TH[[i]]))
        } else {
            stop("unknown ov.types:", ov.types[i])
        }
    }

    # rm VAR columns from ordinal variables
    SC.VAR <- SC.VAR[,-ord.idx]

    if(verbose) {
        cat("STEP 1: univariate statistics\n")
        cat("Threshold + means:\n")
        print(unlist(TH))
        cat("Slopes (if any):\n")
        print(unlist(SLOPES))
        cat("Variances:\n")
        print(unlist(VAR))
    }

    # stage two

    
    pstar <- nvar*(nvar-1)/2
    SC.COR <- matrix(0, N, pstar)
    PSTAR <- matrix(0, nvar, nvar)
    COR.NAMES <- character(pstar)
    H22 <- diag(pstar) # for the delta rule
    # LAVAAN style: col-wise!
    PSTAR[lavaan:::vech.idx(nvar, diag=FALSE)] <- 1:pstar
    # LISREL style: row-wise
    #PSTAR[lavaan:::vechr.idx(nvar, diag=FALSE)] <- 1:pstar
    if(nvar > 1L) {
        for(j in 1:(nvar-1L)) {
            for(i in (j+1L):nvar) {
                #if(verbose) { cat(" i = ", i, " j = ", j, "\n") }
                pstar.idx <- PSTAR[i,j]
                COR.NAMES[pstar.idx] <- paste(ov.names[i],"~~",ov.names[j],sep="")
                if(ov.types[i] == "numeric" && ov.types[j] == "numeric") {
                    if(nexo > 0L) {
                        Y1 <- Data[,i]-FIT[[i]]$yhat; Y2 <- Data[,j]-FIT[[j]]$yhat
                        my1 <- FIT[[i]]$yhat; my2 <- FIT[[j]]$yhat
                    } else {
                        Y1 <- Data[,i]; Y2 <- Data[,j]
                        my1 <- TH[[i]]; my2 <- TH[[j]]
                    }
                    COR[i,j] <- COR[j,i] <- cor(Y1, Y2, use="pairwise.complete.obs")
                    
                    SC.COR[,pstar.idx] <-
                        scores_cor(X=Data[,i], Y=Data[,j], rho=COR[i,j],
                                   mu.x=my1, var.x=VAR[[i]], 
                                   mu.y=my2, var.y=VAR[[j]])
                    H22[pstar.idx,pstar.idx] <- sqrt(VAR[[i]]) * sqrt(VAR[[j]])
                } else if(ov.types[i] == "numeric" && ov.types[j] == "ordered") {
                    # polyserial
                    if(nexo > 0L) {
                        out <- pscorX_TS(X=Data[,i], Y=Data[,j], eXo=NULL,
                                         mu.x=TH[[i]], var.x=VAR[[i]], 
                                         eta.x=FIT[[i]]$yhat,
                                         y.z1=FIT[[j]]$z1, y.z2=FIT[[j]]$z2,
                                         scores=TRUE)
                        SC.COR[,pstar.idx] <- attr(out, "scores")
                    } else {
                        out <- pscor_TS(X=Data[,i], Y=Data[,j], th.y=TH[[j]])
                        COR[i,j] <- COR[j,i] <- out
                        SC.COR[,pstar.idx] <- 
                            scores_pscor(X=Data[,i], Y=Data[,j], rho=out,
                                         mu.x=TH[[i]], var.x=VAR[[i]], th.y=TH[[j]])
                    }
                    COR[i,j] <- COR[j,i] <- out
                    H22[pstar.idx,pstar.idx] <- sqrt(VAR[[i]])
                } else if(ov.types[j] == "numeric" && ov.types[i] == "ordered") {
                    # polyserial
                    if(nexo > 0L) {
                        out <- pscorX_TS(X=Data[,j], Y=Data[,i], eXo=NULL,
                                         mu.x=TH[[j]], var.x=VAR[[j]],
                                         eta.x=FIT[[j]]$yhat,
                                         y.z1=FIT[[i]]$z1, y.z2=FIT[[i]]$z2,
                                         scores=TRUE)
                        SC.COR[,pstar.idx] <- attr(out, "scores")
                    } else {
                        out <- pscor_TS(X=Data[,j], Y=Data[,i], th.y=TH[[i]])
                        SC.COR[,pstar.idx] <- 
                        scores_pscor(X=Data[,j], Y=Data[,i], rho=out,
                                     mu.x=TH[[j]], var.x=VAR[[j]], th.y=TH[[i]])
                    }
                    COR[i,j] <- COR[j,i] <- out
                    H22[pstar.idx,pstar.idx] <- sqrt(VAR[[j]])
                } else if(ov.types[i] == "ordered" && ov.types[j] == "ordered") {
                    # polychoric correlation
                    if(nexo > 0L) {
                        out <- pccorX_TS(x=Data[,i], y=Data[,j], eXo=NULL,
                                         x.z1=FIT[[i]]$z1, x.z2=FIT[[i]]$z2,
                                         y.z1=FIT[[j]]$z1, y.z2=FIT[[j]]$z2,
                                         scores=TRUE)
                        SC.COR[,pstar.idx] <- attr(out, "scores")
                    } else {
                        out <- pccor_TS(Data[,i], Data[,j],
                                        th.x = TH[[i]], th.y = TH[[j]])
                        SC.COR[,pstar.idx] <- 
                            scores_pccor(X=Data[,i], Y=Data[,j], rho=out,
                                         th.x=TH[[i]], th.y=TH[[j]])
                    }
                    COR[i,j] <- COR[j,i] <- out
                }
            }
        }
        colnames(SC.COR) <- COR.NAMES
    }

    if(verbose) {
        cat("\n\nSTEP 2: covariances/correlations:\n")
        print(COR)
    }

    # stage three
    SC <- cbind(SC.TH, SC.SL, SC.VAR, SC.COR)
    INNER <- crossprod(SC)

    # A11
    # new approach (2 June 2012): A11 is just a 'sparse' version of 
    # (the left upper block of) INNER
    A11.size <- ncol(SC.TH) + ncol(SC.SL) + ncol(SC.VAR)
    A11 <- matrix(0, A11.size, A11.size)
    for(i in 1:nvar) {
        th.idx <- th.start.idx[i]:th.end.idx[i]
        sl.idx <- integer(0L)
        var.idx <- integer(0L)
        if(nexo > 0L) {
            sl.end.idx <- (i*nexo); sl.start.idx <- (i-1L)*nexo + 1L
            sl.idx <- ncol(SC.TH) + c(sl.start.idx, sl.end.idx)
        } 
        if(ov.types[i] == "numeric") {
            var.idx <- ncol(SC.TH) + ncol(SC.SL) + match(i, num.idx)
        }
        var.idx <- c(th.idx, sl.idx, var.idx)
        A11[var.idx, var.idx] <- INNER[var.idx, var.idx]
    }

    # A21
    A21 <- matrix(0, pstar, ncol(A11))
    H21 <- matrix(0, pstar, ncol(A11))
    # for this one, we need new scores: for each F_ij (cor), the
    # scores with respect to the TH, VAR, ...
    if(nvar > 1L) {
        for(j in 1:(nvar-1L)) {
            for(i in (j+1L):nvar) {
                pstar.idx <- PSTAR[i,j]
                th.idx_i <- th.start.idx[i]:th.end.idx[i]
                th.idx_j <- th.start.idx[j]:th.end.idx[j]
                if(nexo > 0L) {
                    sl.end.idx <- (i*nexo); sl.start.idx <- (i-1L)*nexo + 1L
                    sl.idx_i <- ncol(SC.TH) + (sl.start.idx:sl.end.idx)
                    sl.end.idx <- (j*nexo); sl.start.idx <- (j-1L)*nexo + 1L
                    sl.idx_j <- ncol(SC.TH) + (sl.start.idx:sl.end.idx)

                    var.idx_i <- ncol(SC.TH) + ncol(SC.SL) + match(i, num.idx)
                    var.idx_j <- ncol(SC.TH) + ncol(SC.SL) + match(j, num.idx)
                } else {
                    var.idx_i <- ncol(SC.TH) + match(i, num.idx)
                    var.idx_j <- ncol(SC.TH) + match(j, num.idx)
                }
                if(ov.types[i] == "numeric" && ov.types[j] == "numeric") {
                    if(nexo > 0) {
                        SC.COR.UNI <-
                            scores_corX_uni(X=Data[,i], Y=Data[,j], eXo=eXo,
                                            rho=COR[i,j],
                                            eta.x=FIT[[i]]$yhat, var.x=VAR[[i]],
                                            eta.y=FIT[[j]]$yhat, var.y=VAR[[j]])
                    } else {
                        SC.COR.UNI <-
                            scores_cor_uni(X=Data[,i], Y=Data[,j], rho=COR[i,j],
                                           mu.x=TH[[i]], var.x=VAR[[i]], 
                                           mu.y=TH[[j]], var.y=VAR[[j]])
                    }
                    # TH
                    A21[pstar.idx, th.idx_i] <-
                        crossprod(SC.COR[,pstar.idx], SC.COR.UNI$dx.mu.x)
                    A21[pstar.idx, th.idx_j] <-
                        crossprod(SC.COR[,pstar.idx], SC.COR.UNI$dx.mu.y)
                    # VAR
                    A21[pstar.idx, var.idx_i] <-
                        crossprod(SC.COR[,pstar.idx], SC.COR.UNI$dx.var.x)
                    A21[pstar.idx, var.idx_j] <-
                        crossprod(SC.COR[,pstar.idx], SC.COR.UNI$dx.var.y)
                    # H21 only needed for VAR
                    H21[pstar.idx, var.idx_i] <-
                        (sqrt(VAR[[j]]) * COR[i,j]) / (2*sqrt(VAR[[i]]))
                    H21[pstar.idx, var.idx_j] <-
                        (sqrt(VAR[[i]]) * COR[i,j]) / (2*sqrt(VAR[[j]]))
                } else if(ov.types[i] == "numeric" && ov.types[j] == "ordered") {
                    # polyserial
                    if(nexo > 0L) {
                        SC.COR.UNI <-
                            scores_pscorX_uni(X=Data[,i], Y=Data[,j], eXo=eXo,
                                              rho=COR[i,j],
                                              eta.x=FIT[[i]]$yhat, var.x=VAR[[i]], 
                                              y.z1=FIT[[j]]$z1, y.z2=FIT[[j]]$z2)
                    } else {
                        SC.COR.UNI <-
                            scores_pscor_uni(X=Data[,i], Y=Data[,j], rho=COR[i,j], 
                                             mu.x=TH[[i]], var.x=VAR[[i]], 
                                             th.y=TH[[j]])
                    }
                    # TH
                    A21[pstar.idx, th.idx_i] <-
                        crossprod(SC.COR[,pstar.idx], SC.COR.UNI$dx.mu.x)
                    A21[pstar.idx, th.idx_j] <-
                        crossprod(SC.COR[,pstar.idx], SC.COR.UNI$dx.th.y)
                    # SL
                    if(nexo > 0L) {
                        A21[pstar.idx, sl.idx_i] <-
                            crossprod(SC.COR[,pstar.idx], SC.COR.UNI$dx.sl.x)
                        A21[pstar.idx, sl.idx_j] <-
                            crossprod(SC.COR[,pstar.idx], SC.COR.UNI$dx.sl.y)
                    }
                    # VAR
                    A21[pstar.idx, var.idx_i] <-
                        crossprod(SC.COR[,pstar.idx], SC.COR.UNI$dx.var.x)
                    # H21 only need for VAR
                    H21[pstar.idx,  var.idx_i] <- COR[i,j] / (2*sqrt(VAR[[i]]))
                } else if(ov.types[j] == "numeric" && ov.types[i] == "ordered") {
                    # polyserial
                    if(nexo > 0L) {
                        SC.COR.UNI <-
                            scores_pscorX_uni(X=Data[,j], Y=Data[,i], eXo=eXo,
                                              rho=COR[i,j],
                                              eta.x=FIT[[j]]$yhat, var.x=VAR[[j]], 
                                              y.z1=FIT[[i]]$z1, y.z2=FIT[[i]]$z2)
                    } else {
                        SC.COR.UNI <-
                            scores_pscor_uni(X=Data[,j], Y=Data[,i], rho=COR[i,j], 
                                             mu.x=TH[[j]], var.x=VAR[[j]],
                                             th.y=TH[[i]])
                    }
                    # TH
                    A21[pstar.idx, th.idx_j] <-
                        crossprod(SC.COR[,pstar.idx], SC.COR.UNI$dx.mu.x)
                    A21[pstar.idx, th.idx_i] <-
                        crossprod(SC.COR[,pstar.idx], SC.COR.UNI$dx.th.y)
                    # SL
                    if(nexo > 0L) {
                        A21[pstar.idx, sl.idx_j] <-
                            crossprod(SC.COR[,pstar.idx], SC.COR.UNI$dx.sl.x)
                        A21[pstar.idx, sl.idx_i] <-
                            crossprod(SC.COR[,pstar.idx], SC.COR.UNI$dx.sl.y)
                    }
                    # VAR
                    A21[pstar.idx, var.idx_j] <-
                        crossprod(SC.COR[,pstar.idx], SC.COR.UNI$dx.var.x)
                    # H21 only for VAR
                    H21[pstar.idx, var.idx_j] <- COR[i,j] / (2*sqrt(VAR[[j]]))
                } else if(ov.types[i] == "ordered" && ov.types[j] == "ordered") {
                    # polychoric correlation
                    if(nexo > 0L) {
                        SC.COR.UNI <-
                            scores_pccorX_uni(X=Data[,i], Y=Data[,j], rho=COR[i,j],
                                              th.x=TH[[i]], th.y=TH[[j]],
                                              sl.x=SLOPES[[i]], sl.y=SLOPES[[j]],
                                              x.z1=FIT[[i]]$z1, x.z2=FIT[[i]]$z2,
                                              y.z1=FIT[[j]]$z1, y.z2=FIT[[j]]$z2,
                                              eXo=eXo)
                    } else {
                        SC.COR.UNI <-
                            scores_pccor_uni(X=Data[,i], Y=Data[,j], rho=COR[i,j],
                                             th.x=TH[[i]], th.y=TH[[j]])
                    }
                    # TH
                   A21[pstar.idx, th.idx_i] <- 
                        crossprod(SC.COR[,pstar.idx], SC.COR.UNI$dx.th.x)
                    A21[pstar.idx, th.idx_j] <- 
                        crossprod(SC.COR[,pstar.idx], SC.COR.UNI$dx.th.y)
                    # SL
                    if(nexo > 0L) {
                    A21[pstar.idx, sl.idx_i] <-
                        crossprod(SC.COR[,pstar.idx], SC.COR.UNI$dx.sl.x)
                    A21[pstar.idx, sl.idx_j] <-
                        crossprod(SC.COR[,pstar.idx], SC.COR.UNI$dx.sl.y)
                    }
                    # NO VAR
                }
            }
        }
    }
       

    # A22
    A22 <- matrix(0, pstar, pstar)
    for(i in seq_len(pstar)) {
        A22[i,i] <- sum( SC.COR[,i]^2 )
    }

    # A12
    A12 <- matrix(0, nrow(A11), ncol(A22))

    B <- rbind( cbind(A11,A12),
                cbind(A21,A22) )
    B.inv <- solve(B)

    NACOV <- (B.inv %*% INNER %*% t(B.inv)) * N

    # COV matrix?
    if(any("numeric" %in% ov.types)) {
        COV <- cor2cov(R=COR, sds=sqrt(unlist(VAR)))

        # construct H matrix to apply delta rule (for the tranformation
        # of rho_ij to cov_ij)
        H11 <- diag(nrow(A11))
        H12 <- matrix(0, nrow(A11), ncol(A22))
        # H22 and H21 already filled in
        H <- rbind( cbind(H11,H12),
                    cbind(H21,H22) )

        NACOV <- H %*% NACOV %*% t(H)
    } else {
        COV <- COR
          H <- diag(ncol(NACOV))
    }
    

    out <- list(TH=TH, SLOPES=SLOPES, VAR=VAR, COR=COR, COV=COV,
                SC=SC,
                INNER=INNER, A11=A11, A12=A12, A21=A21, A22=A22,
                NACOV=NACOV, H=H, TH.NAMES=TH.NAMES, TH.IDX=TH.IDX)
    out
}

