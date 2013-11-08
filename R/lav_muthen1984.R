muthen1984 <- function(Data, ov.names=NULL, ov.types=NULL, ov.levels=NULL,
                       ov.names.x=character(0L), eXo=NULL, verbose=FALSE,
                       missing="listwise",
                       WLS.W=TRUE, # do we need asymptotic variance of stats?
                       zero.add = c(0.5, 0.0),
                       zero.keep.margins = TRUE,
                       group=1L) { # group only for error messages

    # internal function lav_crossprod2
    if(missing == "listwise") {
        lav_crossprod2 <- base::crossprod
    } else { # pairwise, we can have missing values
        lav_crossprod2 <- function(x, y) sum(x * y, na.rm = TRUE)
    }

    # pairwise version
    # FIXME: surely a much better/faster solution is possible??
    lav_crossprod_matrix <- function(A) {
        ndim <- ncol(A)
        # off-diagonal
        upper <- apply(combn(ncol(A),2),2,
                       function(x) sum(A[,x[1]] * A[,x[2]], na.rm=TRUE))
        tmp <- diag(apply(A, 2, function(x) sum(x*x, na.rm=TRUE)))
        tmp[ vechru.idx(ndim, diagonal=FALSE) ] <- upper
        tmp[ vech.idx(  ndim, diagonal=FALSE) ] <- upper
        tmp
    }
    

    #require(mvtnorm)

    # This function was written in January 2012 -- Yves Rosseel
    # First success: Friday 20 Jan 2012: the standard errors for
    #                thresholds and polychoric correlations (in an 
    #                unrestricted/saturated model) are spot on!
    # Second success: Saturday 9 June 2012: support for mixed (ordinal + metric)
    #                 variables; thanks to the delta method to get the ACOV 
    #                 right (see H matrix)
    # Third success: Monday 2 July 2012: support for fixed.x covariates
    # 
    # Friday 13 July: merge exo + non-exo code
    # Monday 16 July: fixed sign numeric in WLS.W; I think we got it right now
    nvar <- ncol(Data); N <- nrow(Data)
    nTH <- ov.levels - 1L; nTH[nTH == -1L] <- 1L
    nth <- sum(nTH)
    th.end.idx <- cumsum(nTH); th.start.idx <- th.end.idx - (nTH - 1L)

    # variable types; default = numeric
    ord.idx <- which(ov.types == "ordered")
    num.idx <- which(ov.types == "numeric"); nnum <- length(num.idx)
    nexo <- length(ov.names.x)
    if(nexo > 0L) stopifnot(ncol(eXo) == nexo)
    pstar <- nvar*(nvar-1)/2

    if(verbose) {
        cat("\nPreparing for WLS estimation -- STEP 1 + 2\n")
        cat("Number of endogenous variables: ", nvar, "\n")
        cat("Endogenous variable names:\n"); print(ov.names); cat("\n")
        cat("Endogenous ov types:\n"); print(ov.types); cat("\n")
        cat("Endogenous ov levels:\n "); print(ov.levels); cat("\n")
        cat("Number of exogenous variables: ", nexo, "\n")
        cat("Exogenous variable names:\n"); print(ov.names.x); cat("\n")
    }

    # means and thresholds
    TH       <- vector("list", length=nvar)
    TH.NOX   <- vector("list", length=nvar)
    TH.NAMES <- vector("list", length=nvar)
    TH.IDX   <- vector("list", length=nvar)
    # slopes (only if fixed.x)
    SLOPES <- matrix(as.numeric(NA), nrow=nvar, ncol=nexo)
    # variances (for continuous variables only)
    VAR <- numeric(length=nvar)
    # correlations
    COR <- diag(nvar); colnames(COR) <- rownames(COR) <- ov.names

    # SCORES
    SC.VAR <- matrix(0, N, nvar); colnames(SC.VAR) <- ov.names
    SC.SL  <- matrix(0, N, nvar*nexo)
    #colnames(SC.SL) <- paste(rep(ov.names, times=nexo), 
    #                         rep(ov.names.x,nvar), sep="")
    SC.TH  <- matrix(0, N, nth)
    SC.COR <- matrix(0, N, pstar)
    COR.NAMES <- character(pstar)
    colnames(SC.TH) <- unlist(lapply(as.list(1:nvar),
        function(x) paste(ov.names[x],"|",1:nTH[x],sep="")))
    FIT <- vector("list", length=nvar)

    # stage one - TH/SLOPES/VAR only
    ov.num <- 0L
    for(i in 1:nvar) {
        th.idx <- th.start.idx[i]:th.end.idx[i]
        sl.idx <- seq(i, by=nvar, length.out=nexo)
        if(ov.types[i] == "numeric") {
            fit <- lavOLS(y=Data[,i], X=eXo); scores <- fit$scores()
            FIT[[i]] <- fit
            ov.num <- ov.num + 1L
            # compute mean and variance
            TH[[i]] <- TH.NOX[[i]] <- fit$theta[1L]
            VAR[i] <- fit$theta[fit$npar]
            TH.NAMES[[i]] <- ov.names[i]; TH.IDX[[i]] <- 0L
            if(WLS.W) {
                SC.TH[,th.idx] <- scores[,1L]
                SC.VAR[,i] <- scores[,fit$npar]
            }
            if(nexo > 0L) {
                SLOPES[i,] <- fit$theta[-c(1L, fit$npar)]
                if(WLS.W) {
                    SC.SL[,sl.idx] <- scores[,-c(1L, fit$npar),drop=FALSE]
                }
                TH.NOX[[i]] <- mean(Data[,i], na.rm=TRUE)
            }
        } else if(ov.types[i] == "ordered") {
            # check if we have enough categories in this group
            # FIXME: should we more tolerant here???
            y.freq <- tabulate(Data[,i], nbins=ov.levels[i])
            if(length(y.freq) != ov.levels[i])
                stop("lavaan ERROR: variable ", ov.names[i], " has fewer categories (", length(y.freq), ") than expected (", ov.levels[i], ") in group ", group)
            if(any(y.freq == 0L))
                stop("lavaan ERROR: some categories of variable `", ov.names[i], "' are empty in group ", group, "; frequencies are [", paste(y.freq, collapse=" "), "]")
            fit <- lavProbit(y=Data[,i], X=eXo); scores <- fit$scores()
            FIT[[i]] <- fit
            TH[[i]] <- fit$theta[fit$th.idx]
            TH.NOX[[i]] <- pc_th(Y=Data[,i])
            if(WLS.W) {
                SC.TH[,th.idx] <- scores[,fit$th.idx,drop=FALSE]
            }
            SLOPES[i,] <- fit$theta[fit$slope.idx]
            if(WLS.W) {
                SC.SL[,sl.idx] <- scores[,fit$slope.idx,drop=FALSE]
            }
            VAR[i] <- 1.0
            TH.NAMES[[i]] <- paste(ov.names[i], "|t", 1:length(TH[[i]]), 
                                   sep="")
            TH.IDX[[i]] <- rep(i, length(TH[[i]]))
        } else {
            stop("unknown ov.types:", ov.types[i])
        }
    }

    # rm VAR columns from ordinal variables
    if(WLS.W) {
        SC.VAR <- SC.VAR[,-ord.idx, drop=FALSE]
    }

    if(verbose) {
        cat("STEP 1: univariate statistics\n")
        cat("Threshold + means:\n")
        TTHH <- unlist(TH)
        names(TTHH) <- unlist(TH.NAMES)
        print(TTHH)
        cat("Slopes (if any):\n")
        colnames(SLOPES) <- ov.names.x
        rownames(SLOPES) <- ov.names
        print(SLOPES)
        cat("Variances:\n")
        names(VAR) <- ov.names
        print(unlist(VAR))
    }

    # stage two
 
    if(verbose) cat("\n\nSTEP 2: covariances/correlations:\n")

    # LAVAAN style: col-wise! (LISREL style: row-wise using vechr.idx)
    PSTAR <- matrix(0, nvar, nvar)
    PSTAR[vech.idx(nvar, diagonal=FALSE)] <- 1:pstar
    for(j in seq_len(nvar-1L)) {
        for(i in (j+1L):nvar) {
            if(verbose) { cat(" i = ", i, " j = ", j, 
                              "[",ov.names[i], "-", ov.names[j], "] ",
                              "(",ov.types[i], "-", ov.types[j], ")\n") }
            pstar.idx <- PSTAR[i,j]
            COR.NAMES[pstar.idx] <- paste(ov.names[i],"~~",ov.names[j],sep="")
            if(ov.types[i] == "numeric" && ov.types[j] == "numeric") {
                if(nexo > 0L) {
                    Y1 <- Data[,i]-FIT[[i]]$yhat; Y2 <- Data[,j]-FIT[[j]]$yhat
                } else {
                    Y1 <- Data[,i]; Y2 <- Data[,j]
                }
                COR[i,j] <- COR[j,i] <- cor(Y1, Y2, use="pairwise.complete.obs")
            } else if(ov.types[i] == "numeric" && ov.types[j] == "ordered") {
                # polyserial
                out <- ps_cor_TS(fit.y1=FIT[[i]], fit.y2=FIT[[j]])
                COR[i,j] <- COR[j,i] <- out
            } else if(ov.types[j] == "numeric" && ov.types[i] == "ordered") {
                # polyserial
                out <- ps_cor_TS(fit.y1=FIT[[j]], fit.y2=FIT[[i]])
                COR[i,j] <- COR[j,i] <- out
            } else if(ov.types[i] == "ordered" && ov.types[j] == "ordered") {
                # polychoric correlation
                out <- pc_cor_TS(fit.y1=FIT[[i]], fit.y2=FIT[[j]],
                                 zero.add = zero.add, 
                                 zero.keep.margins = zero.keep.margins)
                COR[i,j] <- COR[j,i] <- out
            }
            # check for near 1.0 correlations
            if(abs(COR[i,j]) > 0.99) {
                warning("lavaan WARNING: correlation between variables ", ov.names[i], " and ", ov.names[j], " is (nearly) 1.0")
            }
        }
    }

    if(WLS.W) {
        colnames(SC.COR) <- COR.NAMES
    }

    if(verbose) {
        print(COR)
    }

    if(!WLS.W) { # we do not need the asymptotic variance matrix
        if(any("numeric" %in% ov.types)) {
            COV <- cor2cov(R=COR, sds=sqrt(unlist(VAR)))
        } else {
            COV <- COR
        }
        out <- list(TH=TH, SLOPES=SLOPES, VAR=VAR, COR=COR, COV=COV,
                SC=NULL, TH.NOX=TH.NOX,TH.NAMES=TH.NAMES, TH.IDX=TH.IDX,
                INNER=NULL, A11=NULL, A12=NULL, A21=NULL, A22=NULL,
                WLS.W=NULL, H=NULL)
        return(out)
    }

    A11.size <- ncol(SC.TH) + ncol(SC.SL) + ncol(SC.VAR)

    # A21
    A21 <- matrix(0, pstar, A11.size)
    H22 <- diag(pstar) # for the delta rule
    H21 <- matrix(0, pstar, A11.size)
    # for this one, we need new scores: for each F_ij (cor), the
    # scores with respect to the TH, VAR, ...
    for(j in seq_len(nvar-1L)) {
        for(i in (j+1L):nvar) {
            pstar.idx <- PSTAR[i,j]
            th.idx_i <- th.start.idx[i]:th.end.idx[i]
            th.idx_j <- th.start.idx[j]:th.end.idx[j]
            if(nexo > 0L) {
                sl.idx_i <- ncol(SC.TH) + seq(i, by=nvar, length.out=nexo)
                sl.idx_j <- ncol(SC.TH) + seq(j, by=nvar, length.out=nexo)

                var.idx_i <- ncol(SC.TH) + ncol(SC.SL) + match(i, num.idx)
                var.idx_j <- ncol(SC.TH) + ncol(SC.SL) + match(j, num.idx)
            } else {
                var.idx_i <- ncol(SC.TH) + match(i, num.idx)
                var.idx_j <- ncol(SC.TH) + match(j, num.idx)
            }
            if(ov.types[i] == "numeric" && ov.types[j] == "numeric") {
                SC.COR.UNI <- pp_cor_scores(rho=COR[i,j],
                                            fit.y1=FIT[[i]],
                                            fit.y2=FIT[[j]])

                # RHO
                SC.COR[,pstar.idx] <- SC.COR.UNI$dx.rho

                # TH
                A21[pstar.idx, th.idx_i] <-
                    lav_crossprod2(SC.COR[,pstar.idx], SC.COR.UNI$dx.mu.y1)
                A21[pstar.idx, th.idx_j] <-
                    lav_crossprod2(SC.COR[,pstar.idx], SC.COR.UNI$dx.mu.y2)
                # SL
                if(nexo > 0L) {
                    A21[pstar.idx, sl.idx_i] <-
                        lav_crossprod2(SC.COR[,pstar.idx], SC.COR.UNI$dx.sl.y1)
                    A21[pstar.idx, sl.idx_j] <-
                        lav_crossprod2(SC.COR[,pstar.idx], SC.COR.UNI$dx.sl.y2)
                }
                # VAR
                A21[pstar.idx, var.idx_i] <-
                    lav_crossprod2(SC.COR[,pstar.idx], SC.COR.UNI$dx.var.y1)
                A21[pstar.idx, var.idx_j] <-
                    lav_crossprod2(SC.COR[,pstar.idx], SC.COR.UNI$dx.var.y2)
                # H21 only needed for VAR
                H21[pstar.idx, var.idx_i] <-
                    (sqrt(VAR[j]) * COR[i,j]) / (2*sqrt(VAR[i]))
                H21[pstar.idx, var.idx_j] <-
                    (sqrt(VAR[i]) * COR[i,j]) / (2*sqrt(VAR[j]))
                H22[pstar.idx, pstar.idx] <- sqrt(VAR[i]) * sqrt(VAR[j])
            } else if(ov.types[i] == "numeric" && ov.types[j] == "ordered") {
                SC.COR.UNI <- ps_cor_scores(rho=COR[i,j],
                                            fit.y1=FIT[[i]], 
                                            fit.y2=FIT[[j]])
                # RHO
                SC.COR[,pstar.idx] <- SC.COR.UNI$dx.rho

                # TH
                A21[pstar.idx, th.idx_i] <-
                    lav_crossprod2(SC.COR[,pstar.idx], SC.COR.UNI$dx.mu.y1)
                A21[pstar.idx, th.idx_j] <-
                    lav_crossprod2(SC.COR[,pstar.idx], SC.COR.UNI$dx.th.y2)
                # SL
                if(nexo > 0L) {
                    A21[pstar.idx, sl.idx_i] <-
                        lav_crossprod2(SC.COR[,pstar.idx], SC.COR.UNI$dx.sl.y1)
                    A21[pstar.idx, sl.idx_j] <-
                        lav_crossprod2(SC.COR[,pstar.idx], SC.COR.UNI$dx.sl.y2)
                }
                # VAR
                A21[pstar.idx, var.idx_i] <-
                    lav_crossprod2(SC.COR[,pstar.idx], SC.COR.UNI$dx.var.y1)
                # H21 only need for VAR
                H21[pstar.idx,  var.idx_i] <- COR[i,j] / (2*sqrt(VAR[i]))
                H22[pstar.idx, pstar.idx] <- sqrt(VAR[i])
            } else if(ov.types[j] == "numeric" && ov.types[i] == "ordered") {
                SC.COR.UNI <- ps_cor_scores(rho=COR[i,j],
                                            fit.y1=FIT[[j]], 
                                            fit.y2=FIT[[i]])
                # RHO
                SC.COR[,pstar.idx] <- SC.COR.UNI$dx.rho

                # TH
                A21[pstar.idx, th.idx_j] <-
                    lav_crossprod2(SC.COR[,pstar.idx], SC.COR.UNI$dx.mu.y1)
                A21[pstar.idx, th.idx_i] <-
                    lav_crossprod2(SC.COR[,pstar.idx], SC.COR.UNI$dx.th.y2)
                # SL
                if(nexo > 0L) {
                    A21[pstar.idx, sl.idx_j] <-
                        lav_crossprod2(SC.COR[,pstar.idx], SC.COR.UNI$dx.sl.y1)
                    A21[pstar.idx, sl.idx_i] <-
                        lav_crossprod2(SC.COR[,pstar.idx], SC.COR.UNI$dx.sl.y2)
                }
                # VAR
                A21[pstar.idx, var.idx_j] <-
                    lav_crossprod2(SC.COR[,pstar.idx], SC.COR.UNI$dx.var.y1)
                # H21 only for VAR
                H21[pstar.idx, var.idx_j] <- COR[i,j] / (2*sqrt(VAR[j]))
                H22[pstar.idx, pstar.idx] <- sqrt(VAR[j])
            } else if(ov.types[i] == "ordered" && ov.types[j] == "ordered") {
                # polychoric correlation
                SC.COR.UNI <- pc_cor_scores(rho=COR[i,j], 
                                            fit.y1=FIT[[i]], 
                                            fit.y2=FIT[[j]])
                # RHO
                SC.COR[,pstar.idx] <- SC.COR.UNI$dx.rho

                # TH
                A21[pstar.idx, th.idx_i] <- 
                    lav_crossprod2(SC.COR[,pstar.idx], SC.COR.UNI$dx.th.y1)
                A21[pstar.idx, th.idx_j] <- 
                    lav_crossprod2(SC.COR[,pstar.idx], SC.COR.UNI$dx.th.y2)
                # SL
                if(nexo > 0L) {
                    A21[pstar.idx, sl.idx_i] <-
                        lav_crossprod2(SC.COR[,pstar.idx], SC.COR.UNI$dx.sl.y1)
                    A21[pstar.idx, sl.idx_j] <-
                        lav_crossprod2(SC.COR[,pstar.idx], SC.COR.UNI$dx.sl.y2)
                }
                # NO VAR
            }
        }
    }

    # stage three
    SC <- cbind(SC.TH, SC.SL, SC.VAR, SC.COR)
    if(missing == "listwise") {
        INNER <- crossprod(SC)
    } else {
        INNER <- lav_crossprod_matrix(SC)
    }

    # A11
    # new approach (2 June 2012): A11 is just a 'sparse' version of 
    # (the left upper block of) INNER
    A11 <- matrix(0, A11.size, A11.size)
    for(i in 1:nvar) {
        th.idx <- th.start.idx[i]:th.end.idx[i]
        sl.idx <- integer(0L)
        var.idx <- integer(0L)
        if(nexo > 0L) {
            sl.idx <- ncol(SC.TH) + seq(i, by=nvar, length.out=nexo)
            #sl.end.idx <- (i*nexo); sl.start.idx <- (i-1L)*nexo + 1L
            #sl.idx <- ncol(SC.TH) + (sl.start.idx:sl.end.idx)
        }
        if(ov.types[i] == "numeric") {
            var.idx <- ncol(SC.TH) + ncol(SC.SL) + match(i, num.idx)
        }
        a11.idx <- c(th.idx, sl.idx, var.idx)
        A11[a11.idx, a11.idx] <- INNER[a11.idx, a11.idx]
    }

    # A22
    A22 <- matrix(0, pstar, pstar)
    for(i in seq_len(pstar)) {
        A22[i,i] <- sum( SC.COR[,i]^2, na.rm=TRUE )
    }

    # A12
    A12 <- matrix(0, nrow(A11), ncol(A22))

    B <- rbind( cbind(A11,A12),
                cbind(A21,A22) )
    B.inv <- solve(B)
    ## FIXME: we need to invert B as a partioned matrix

    #  weight matrix (correlation metric)
    WLS.W <- B.inv %*% INNER %*% t(B.inv)

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

        WLS.W <- H %*% WLS.W %*% t(H)
    } else {
        COV <- COR
          H <- diag(ncol(WLS.W))
    }

    # reverse sign numeric TH (because we provide -mu in WLS.obs)
    # (WOW, it took me a LOOONGGG time to realize this!)
    # YR 16 July 2012
    if(length(num.idx) > 0L) {
        WLS.W[num.idx,] <- -WLS.W[num.idx,]
        WLS.W[,num.idx] <- -WLS.W[,num.idx]
    }
    
    out <- list(TH=TH, SLOPES=SLOPES, VAR=VAR, COR=COR, COV=COV,
                SC=SC, TH.NOX=TH.NOX,TH.NAMES=TH.NAMES, TH.IDX=TH.IDX,
                INNER=INNER, A11=A11, A12=A12, A21=A21, A22=A22,
                WLS.W=WLS.W, H=H)
    out
}

