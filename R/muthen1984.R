muthen1984 <- function(Data, ov.names=NULL, ov.types=NULL, ov.levels=NULL,
                       ov.names.x=NULL, eXo=NULL) {

    require(mvtnorm)

    # This function was written in January 2012 -- Yves Rosseel
    # First success: Friday 20 Jan 2012: the standard errors for
    #                thresholds and polychoric correlations (in an 
    #                unrestricted/saturated model) are spot on!
    # Second success: Saturday 9 June 2012: support for mixed (ordinal + metric)
    #                 variables; thanks to the delta method to get the ACOV 
    #                 right (see H matrix)

    nvar <- ncol(Data); stopifnot(nvar > 1L); N <- nrow(Data)
    nTH <- ov.levels - 1L; nTH[nTH == -1L] <- 1L
    nth <- sum(nTH)
    th.end.idx <- cumsum(nTH); th.start.idx <- th.end.idx - (nTH - 1L)

    # variable types; default = numeric
    ord.idx <- which(ov.types == "ordered")
    num.idx <- which(ov.types == "numeric"); nnum <- length(num.idx)
    nexo <- length(ov.names.x)

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
        cat("var = ", ov.names[i], " type = ", ov.types[i], "\n")
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


    # stage two

    if(nexo > 0L) {
        # regress eXo out of numeric y's first
        exo.fit <- lm.fit(x=cbind(1,eXo), y=Data[,num.idx])
        Data[,num.idx] <- exo.fit$residuals
    }

    pstar <- nvar*(nvar-1)/2
    SC.COR <- matrix(0, N, pstar)
    PSTAR <- matrix(0, nvar, nvar)
    COR.NAMES <- character(pstar)
    H22 <- diag(pstar) # for the delta rule
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
                H22[pstar.idx,pstar.idx] <- sqrt(VAR[[i]]) * sqrt(VAR[[j]])
            } else if(ov.types[i] == "numeric" && ov.types[j] == "ordered") {
                # polyserial
                if(nexo > 0L) {
                    out <- pscorX_TS(X=Data[,i], Y=Data[,j], eXo=NULL,
                                     mu.x=TH[[i]], var.x=VAR[[i]], 
                                     XRESID=Data[,i], 
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
                                     XRESID=Data[,j], 
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

    # stage three
    SC <- cbind(SC.TH, SC.SL, SC.VAR, SC.COR)
    INNER <- crossprod(SC)

    # A11
    if(nexo > 0L) {
        A11_block <- vector("list", length=(nvar + nvar + nnum))
    } else {
        A11_block <- vector("list", length=(nvar + nnum))
    }
    # A11 - TH
    for(i in 1:nvar) {
        th.idx <- th.start.idx[i]:th.end.idx[i]
        SC.TH_i <- SC.TH[,th.idx]
        A11_block[[i]] <- crossprod(SC.TH_i)
    }
    # A11 - SLOPES
    if(nexo > 0L) {
        for(i in 1:nvar) {
            sl.end.idx <- (i*nexo)
            sl.start.idx <- (i-1L)*nexo + 1L
            A11_block[[ nvar + i ]] <- 
                crossprod(SC.SL[,sl.start.idx:sl.end.idx])
        }

    }
    # A11 - VAR
    if(nnum > 0L) {
        offset <- 0L
        if(nexo > 0L) 
            offset <- nvar
        for(i in 1:ncol(SC.VAR)) {
            A11_block[[ nvar + offset + i ]] <- crossprod(SC.VAR[,i])
        }
    }
    #print(A11_block)
    A11 <- lavaan:::bdiag(A11_block)

    # A21
    A21 <- matrix(0, pstar, ncol(A11))
    H21 <- matrix(0, pstar, ncol(A11))
    # for this one, we need new scores: for each F_ij (cor), the
    # scores with respect to the TH, VAR, ...
    for(j in 1:(nvar-1L)) {
        for(i in (j+1L):nvar) {
            pstar.idx <- PSTAR[i,j]
            th.idx_i <- th.start.idx[i]:th.end.idx[i]
            th.idx_j <- th.start.idx[j]:th.end.idx[j]
            if(ov.types[i] == "numeric" && ov.types[j] == "numeric") {
                SC.COR.UNI <-
                    scores_cor_uni(X=Data[,i], Y=Data[,j], rho=COR[i,j],
                                   mu.x=TH[[i]], var.x=VAR[[i]], 
                                   mu.y=TH[[j]], var.y=VAR[[j]])
                # TH
                A21[pstar.idx, th.idx_i] <-
                    crossprod(SC.COR[,pstar.idx],
                              SC.COR.UNI[,1L])
                A21[pstar.idx, th.idx_j] <-
                    crossprod(SC.COR[,pstar.idx],
                              SC.COR.UNI[,2L])
                # VAR
                A21[pstar.idx, ncol(SC.TH) + match(i, num.idx)] <-
                    crossprod(SC.COR[,pstar.idx],
                              SC.COR.UNI[,3L])
                A21[pstar.idx, ncol(SC.TH) + match(j, num.idx)] <-
                    crossprod(SC.COR[,pstar.idx],
                              SC.COR.UNI[,4L])
                # H21 only needed for VAR
                H21[pstar.idx, ncol(SC.TH) + match(i, num.idx)] <-
                    (sqrt(VAR[[j]]) * COR[i,j]) / (2*sqrt(VAR[[i]]))
                H21[pstar.idx,  ncol(SC.TH) + match(j, num.idx)] <- 
                    (sqrt(VAR[[i]]) * COR[i,j]) / (2*sqrt(VAR[[j]]))
            } else if(ov.types[i] == "numeric" && ov.types[j] == "ordered") {
                # polyserial
                SC.COR.UNI <-
                    scores_pscor_uni(X=Data[,i], Y=Data[,j], rho=COR[i,j], 
                                     mu.x=TH[[i]], var.x=VAR[[i]], th.y=TH[[j]])
                # TH
                A21[pstar.idx, th.idx_i] <-
                    crossprod(SC.COR[,pstar.idx],
                              SC.COR.UNI[,1L])
                A21[pstar.idx, th.idx_j] <-
                    crossprod(SC.COR[,pstar.idx],
                              SC.COR.UNI[,2L+(1:length(TH[[j]]))])
                # VAR
                A21[pstar.idx, ncol(SC.TH) + match(i, num.idx)] <-
                    crossprod(SC.COR[,pstar.idx],
                              SC.COR.UNI[,2L])
                # H21 only need for VAR
                H21[pstar.idx, ncol(SC.TH) + match(i, num.idx)] <- 
                    COR[i,j] / (2*sqrt(VAR[[i]]))
            } else if(ov.types[j] == "numeric" && ov.types[i] == "ordered") {
                # polyserial
                SC.COR.UNI <-
                    scores_pscor_uni(X=Data[,j], Y=Data[,i], rho=COR[i,j], 
                                     mu.x=TH[[j]], var.x=VAR[[j]], th.y=TH[[i]])
                # TH
                A21[pstar.idx, th.idx_j] <-
                    crossprod(SC.COR[,pstar.idx],
                              SC.COR.UNI[,1L])
                A21[pstar.idx, th.idx_i] <-
                    crossprod(SC.COR[,pstar.idx],
                              SC.COR.UNI[,2L+(1:length(TH[[i]]))])
                # VAR
                A21[pstar.idx, ncol(SC.TH) + match(j, num.idx)] <-
                    crossprod(SC.COR[,pstar.idx],
                              SC.COR.UNI[,2L])
                # H21 only for VAR
                H21[pstar.idx, ncol(SC.TH) + match(j, num.idx)] <- 
                    COR[i,j] / (2*sqrt(VAR[[j]]))
            } else if(ov.types[i] == "ordered" && ov.types[j] == "ordered") {
                # polychoric correlation
                SC.COR.UNI <-
                    scores_pccor_uni(X=Data[,i], Y=Data[,j], rho=COR[i,j],
                                     th.x=TH[[i]], th.y=TH[[j]])
                # TH
                A21[pstar.idx, th.idx_i] <- 
                    crossprod(SC.COR[,pstar.idx], 
                              SC.COR.UNI[,1:length(TH[[i]])])
                A21[pstar.idx, th.idx_j] <- 
                    crossprod(SC.COR[,pstar.idx],
                              SC.COR.UNI[,length(TH[[i]])+1:length(TH[[j]])])
                # NO VAR
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
    # DEBUG ONLY
    A21 <- matrix(0, nrow(A22), ncol(A11))

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
                INNER=INNER, A11=A11, A12=A12, A21=A21, A22=A22,
                NACOV=NACOV, H=H, TH.NAMES=TH.NAMES, TH.IDX=TH.IDX)
    out
}

