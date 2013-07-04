# matrix representation of latent variable models
# and matrix-representation specific functions:
# - computeSigmaHat
# - computeMuHat
# - derivative.F

# initital version: YR 2011-01-21: LISREL stuff
# updates:          YR 2011-12-01: group specific extraction
#                   YR 2012-05-17: thresholds

representation.LISREL <- function(partable=NULL, target=NULL, 
                                  extra=FALSE) {

    # prepare target list
    if(is.null(target)) target <- partable 

    # prepare output
    N <- length(target$lhs)
    tmp.mat <- character(N); tmp.row <- integer(N); tmp.col <- integer(N)

    # global settings
    meanstructure <- any(partable$op == "~1")
    categorical   <- any(partable$op == "|")
    gamma <- categorical

    # number of groups
    ngroups <- max(partable$group)
    ov.dummy.names <- vector("list", ngroups)
    if(extra) {
        REP.mmNames     <- vector("list", ngroups)
        REP.mmNumber    <- vector("list", ngroups)
        REP.mmRows      <- vector("list", ngroups)
        REP.mmCols      <- vector("list", ngroups)
        REP.mmDimNames  <- vector("list", ngroups)
        REP.mmSymmetric <- vector("list", ngroups)
    }

    for(g in 1:ngroups) {

        # info from user model per group
        if(gamma) {
            ov.names <- vnames(partable, "ov.nox",  group=g)
        } else {
            ov.names <- vnames(partable, "ov",  group=g)
        }
        nvar <- length(ov.names)
        lv.names   <- vnames(partable, "lv",  group=g); nfac <- length(lv.names)
        ov.th      <- vnames(partable, "th",  group=g); nth  <- length(ov.th)
        ov.names.x <- vnames(partable, "ov.x",group=g); nexo <- length(ov.names.x)

        # in this representation, we need to create 'phantom/dummy' latent 
        # variables for all `x' and `y' variables not in lv.names
        # (only y if categorical)

        # regression dummys
        if(categorical) {
            tmp.names <-
                unique( partable$lhs[(partable$op == "~" | 
                                        partable$op == "<~") &
                                        partable$group == g] )
        } else {
            tmp.names <- 
                unique( c(partable$lhs[(partable$op == "~" | 
                                        partable$op == "<~") & 
                                        partable$group == g],
                          partable$rhs[(partable$op == "~" | 
                                        partable$op == "<~") & 
                                        partable$group == g]) )
        }
        dummy.names1 <- tmp.names[ !tmp.names %in% lv.names ]

        # covariances involving dummys
        dummy.cov.idx <- which(partable$op == "~~" & partable$group == g &
                               (partable$lhs %in% dummy.names1 |
                                partable$rhs %in% dummy.names1))
        dummy.names2 <- unique( c(partable$lhs[dummy.cov.idx],
                                  partable$rhs[dummy.cov.idx]) )
        # collect all dummy variables
        dummy.names <- unique(c(dummy.names1, dummy.names2))


        if(length(dummy.names)) {
            # make sure order is the same as ov.names
            ov.dummy.names[[g]] <- ov.names[ ov.names %in% dummy.names ]

            # extend lv.names
            lv.names <- c(lv.names, ov.dummy.names[[g]])
            nfac <- length(lv.names)

            # add 'dummy' =~ entries
            dummy.mat <- rep("lambda", length(dummy.names))
        } else {
            ov.dummy.names[[g]] <- character(0)
        }

        # 1a. "=~" regular indicators
        idx <- which(target$group == g &
                     target$op == "=~" & !(target$rhs %in% lv.names))
        tmp.mat[idx] <- "lambda"
        tmp.row[idx] <- match(target$rhs[idx], ov.names)
        tmp.col[idx] <- match(target$lhs[idx], lv.names)

        # 1b. "=~" regular higher-order lv indicators
        idx <- which(target$group == g &
                     target$op == "=~" & !(target$rhs %in% ov.names))
        tmp.mat[idx] <- "beta"
        tmp.row[idx] <- match(target$rhs[idx], lv.names)
        tmp.col[idx] <- match(target$lhs[idx], lv.names)
    
        # 1c. "=~" indicators that are both in ov and lv
        idx <- which(target$group == g &
                     target$op == "=~" & target$rhs %in% ov.names
                                       & target$rhs %in% lv.names)
        tmp.mat[idx] <- "beta"
        tmp.row[idx] <- match(target$rhs[idx], lv.names)
        tmp.col[idx] <- match(target$lhs[idx], lv.names)
    
        # 2. "~" regressions
        if(categorical) {
            # gamma
            idx <- which(target$rhs %in% ov.names.x &
                         target$group == g & (target$op == "~" |
                                              target$op == "<~") )
            tmp.mat[idx] <- "gamma"
            tmp.row[idx] <- match(target$lhs[idx], lv.names)
            tmp.col[idx] <- match(target$rhs[idx], ov.names.x)

            # beta
            idx <- which(!target$rhs %in% ov.names.x &
                         target$group == g & (target$op == "~" |
                                              target$op == "<~") )
            tmp.mat[idx] <- "beta"
            tmp.row[idx] <- match(target$lhs[idx], lv.names)
            tmp.col[idx] <- match(target$rhs[idx], lv.names)
        } else {
            idx <- which(target$group == g & (target$op == "~" |
                                              target$op == "<~") )
            tmp.mat[idx] <- "beta"
            tmp.row[idx] <- match(target$lhs[idx], lv.names)
            tmp.col[idx] <- match(target$rhs[idx], lv.names)
        }
  
        # 3a. "~~" ov
        idx <- which(target$group == g &
                     target$op == "~~" & !(target$lhs %in% lv.names))
        tmp.mat[idx] <- "theta"
        tmp.row[idx] <- match(target$lhs[idx], ov.names)
        tmp.col[idx] <- match(target$rhs[idx], ov.names)
    
        # 3b. "~~" lv
        idx <- which(target$group == g &
                     target$op == "~~" & target$rhs %in% lv.names)
        tmp.mat[idx] <- "psi"
        tmp.row[idx] <- match(target$lhs[idx], lv.names)
        tmp.col[idx] <- match(target$rhs[idx], lv.names)
  
        # 4a. "~1" ov
        idx <- which(target$group == g &
                     target$op == "~1" & !(target$lhs %in% lv.names))
        tmp.mat[idx] <- "nu"
        tmp.row[idx] <- match(target$lhs[idx], ov.names)
        tmp.col[idx] <- 1L
    
        # 4b. "~1" lv
        idx <- which(target$group == g &
                     target$op == "~1" & target$lhs %in% lv.names)
        tmp.mat[idx] <- "alpha"
        tmp.row[idx] <- match(target$lhs[idx], lv.names)
        tmp.col[idx] <- 1L

        # 5. "|" th
        LABEL <- paste(target$lhs, target$op, target$rhs, sep="")
        idx <-  which(target$group == g & 
                      target$op == "|" & LABEL %in% ov.th)
        TH <- paste(target$lhs[idx], "|", target$rhs[idx], sep="")
        tmp.mat[idx] <- "tau"
        tmp.row[idx] <- match(TH, ov.th)
        tmp.col[idx] <- 1L

        # 6. "~*~" scales
        idx <- which(target$group == g &
                     target$op == "~*~")
        tmp.mat[idx] <- "delta"
        tmp.row[idx] <- match(target$lhs[idx], ov.names)
        tmp.col[idx] <- 1L

        # new 0.5-12: catch lower-elements in theta/psi
        idx.lower <- which(tmp.mat %in% c("theta","psi") & tmp.row > tmp.col)
        if(length(idx.lower) > 0L) {
            tmp <- tmp.row[idx.lower]
            tmp.row[idx.lower] <- tmp.col[idx.lower]
            tmp.col[idx.lower] <- tmp
        }

        if(extra) {
            # mRows
            mmRows <- list(tau    = nth,
                           delta  = nvar,
                           nu     = nvar,
                           lambda = nvar,
                           theta  = nvar,
                           alpha  = nfac,
                           beta   = nfac,
                           gamma  = nfac,
                           psi    = nfac)

            # mCols
            mmCols <- list(tau    = 1L,
                           delta  = 1L,
                           nu     = 1L,
                           lambda = nfac,
                           theta  = nvar,
                           alpha  = 1L,
                           beta   = nfac,
                           gamma  = nexo,
                           psi    = nfac)

            # dimNames for LISREL model matrices
            mmDimNames <- list(tau    = list( ov.th,    "threshold"),
                               delta  = list( ov.names,    "scales"),
                               nu     = list( ov.names, "intercept"),
                               lambda = list( ov.names,    lv.names),
                               theta  = list( ov.names,    ov.names),
                               alpha  = list( lv.names, "intercept"),
                               beta   = list( lv.names,    lv.names),
                               gamma  = list( lv.names,  ov.names.x),
                               psi    = list( lv.names,    lv.names))
    
            # isSymmetric
            mmSymmetric <- list(tau    = FALSE,
                                delta  = FALSE,
                                nu     = FALSE,
                                lambda = FALSE,
                                theta  = TRUE,
                                alpha  = FALSE,
                                beta   = FALSE,
                                gamma  = FALSE,
                                psi    = TRUE)
    
            # which mm's do we need? (always include lambda, theta and psi)
            mmNames <- c("lambda", "theta", "psi")
            if("beta" %in% tmp.mat) mmNames <- c(mmNames, "beta")
            if(meanstructure) mmNames <- c(mmNames, "nu", "alpha")
            if("tau" %in% tmp.mat) mmNames <- c(mmNames, "tau")
            if("delta" %in% tmp.mat) mmNames <- c(mmNames, "delta")
            if("gamma" %in% tmp.mat) mmNames <- c(mmNames, "gamma")

            REP.mmNames[[g]]     <- mmNames
            REP.mmNumber[[g]]    <- length(mmNames)
            REP.mmRows[[g]]      <- unlist(mmRows[ mmNames ])
            REP.mmCols[[g]]      <- unlist(mmCols[ mmNames ])
            REP.mmDimNames[[g]]  <- mmDimNames[ mmNames ]
            REP.mmSymmetric[[g]] <- unlist(mmSymmetric[ mmNames ])
        } # extra
    } # ngroups

    REP <- list(mat = tmp.mat,
                row = tmp.row,
                col = tmp.col)

    # always add 'ov.dummy.names' attribute
    attr(REP, "ov.dummy.names") <- ov.dummy.names

    if(extra) {
        attr(REP, "mmNames")     <- REP.mmNames
        attr(REP, "mmNumber")    <- REP.mmNumber
        attr(REP, "mmRows")      <- REP.mmRows
        attr(REP, "mmCols")      <- REP.mmCols
        attr(REP, "mmDimNames")  <- REP.mmDimNames
        attr(REP, "mmSymmetric") <- REP.mmSymmetric
    }

    REP
}


# compute SigmaHat for a single group
computeSigmaHat.LISREL <- function(MLIST=NULL, delta=TRUE) {

    LAMBDA <- MLIST$lambda; nvar <- nrow(LAMBDA)
    PSI    <- MLIST$psi
    THETA  <- MLIST$theta
    BETA   <- MLIST$beta

    # beta?
    if(is.null(BETA)) {
        LAMBDA..IB.inv <- LAMBDA
    } else {
        tmp <- -BETA; nr <- nrow(BETA); i <- seq_len(nr);
        tmp[cbind(i, i)] <- 1
        IB.inv <- solve(tmp)
        LAMBDA..IB.inv <- LAMBDA %*% IB.inv
    }

    # compute Sigma Hat
    Sigma.hat <- tcrossprod(LAMBDA..IB.inv %*% PSI, LAMBDA..IB.inv) + THETA

    # if delta, scale
    if(delta && !is.null(MLIST$delta)) {
        DELTA <- diag(MLIST$delta[,1L], nrow=nvar, ncol=nvar)
        Sigma.hat <- DELTA %*% Sigma.hat %*% DELTA
    }

    Sigma.hat
}

# compute MuHat for a single group
computeMuHat.LISREL <- function(MLIST=NULL) {

    NU     <- MLIST$nu
    ALPHA  <- MLIST$alpha
    LAMBDA <- MLIST$lambda
    BETA   <- MLIST$beta

    # shortcut
    if(is.null(ALPHA) || is.null(NU)) return(matrix(0, nrow(LAMBDA), 1L))

    # beta?
    if(is.null(BETA)) {
        LAMBDA..IB.inv <- LAMBDA
    } else {
        tmp <- -BETA; nr <- nrow(BETA); i <- seq_len(nr);
        tmp[cbind(i, i)] <- 1
        IB.inv <- solve(tmp)
        LAMBDA..IB.inv <- LAMBDA %*% IB.inv
    }
    
    # compute Mu Hat
    Mu.hat <- NU + LAMBDA..IB.inv %*% ALPHA

    Mu.hat
}

# compute TH for a single group
computeTH.LISREL <- function(MLIST=NULL, th.idx=NULL) {

    LAMBDA <- MLIST$lambda; nvar <- nrow(LAMBDA); nfac <- ncol(LAMBDA)
    BETA   <- MLIST$beta
    TAU    <- MLIST$tau; nth <- nrow(TAU)

    # missing alpha
    if(is.null(MLIST$alpha))
        ALPHA <- matrix(0, nfac, 1L)
    else
        ALPHA  <- MLIST$alpha

    # missing nu
    if(is.null(MLIST$nu))
        NU <- matrix(0, nvar, 1L)
    else
        NU <- MLIST$nu

    if(is.null(th.idx)) {
        th.idx <- 1:nth
        nlev <- rep(1L, nvar)
        K_nu <- diag(nvar)
    } else {
        nlev <- tabulate(th.idx, nbins=nvar); nlev[nlev == 0L] <- 1L
        K_nu <- matrix(0, sum(nlev), nvar)
        K_nu[ cbind(1:sum(nlev), rep(1:nvar, times=nlev)) ] <- 1.0
    }

    # shortcut
    if(is.null(TAU)) return(matrix(0, length(th.idx), 1L))

    # beta?
    if(is.null(BETA)) {
        LAMBDA..IB.inv <- LAMBDA
    } else {
        tmp <- -BETA; nr <- nrow(BETA); i <- seq_len(nr);
        tmp[cbind(i, i)] <- 1
        IB.inv <- solve(tmp)
        LAMBDA..IB.inv <- LAMBDA %*% IB.inv
    }
   
    # compute pi0
    pi0 <- NU + LAMBDA..IB.inv %*% ALPHA

    # interleave th's with zeros where we have numeric variables
    th <- numeric( length(th.idx) )
    th[ th.idx > 0L ] <- TAU[,1L]

    # compute TH
    TH <- th - (K_nu %*% pi0)

    # if delta, scale
    if(!is.null(MLIST$delta)) {
        DELTA.diag <- MLIST$delta[,1L]
        DELTA.star.diag <- rep(DELTA.diag, times=nlev)
        TH <- TH * DELTA.star.diag
    }

    as.vector(TH)
}

# compute PI for a single group
computePI.LISREL <- function(MLIST=NULL) {

    LAMBDA <- MLIST$lambda
    BETA   <- MLIST$beta
    GAMMA  <- MLIST$gamma

    # shortcut
    if(is.null(GAMMA)) return(matrix(0, nrow(LAMBDA), 0L))

    # beta?
    if(is.null(BETA)) {
        LAMBDA..IB.inv <- LAMBDA
    } else {
        tmp <- -BETA; nr <- nrow(BETA); i <- seq_len(nr);
        tmp[cbind(i, i)] <- 1
        IB.inv <- solve(tmp)
        LAMBDA..IB.inv <- LAMBDA %*% IB.inv
    }

    # compute PI
    PI <- LAMBDA..IB.inv %*% GAMMA

    # if delta, scale
    if(!is.null(MLIST$delta)) {
        DELTA.diag <- MLIST$delta[,1L]
        PI <- PI * DELTA.diag
    }

    PI
}

# compute V(ETA): variances/covariances of latents
computeVETA.LISREL <- function(MLIST=NULL, cov.x=NULL) {

    LAMBDA <- MLIST$lambda; nvar <- nrow(LAMBDA)
    PSI    <- MLIST$psi
    THETA  <- MLIST$theta
    BETA   <- MLIST$beta

    # beta?
    if(is.null(BETA)) {
        SY <- PSI
    } else {
        tmp <- -BETA; nr <- nrow(BETA); i <- seq_len(nr);
        tmp[cbind(i, i)] <- 1
        IB.inv <- solve(tmp)
        SY <- tcrossprod(IB.inv %*% PSI, IB.inv)
    }

    # if GAMMA, also x part
    GAMMA <- MLIST$gamma
    if(!is.null(GAMMA)) {
        stopifnot(!is.null(cov.x))
        if(is.null(BETA)) {
            SX <- tcrossprod(GAMMA %*% cov.x, GAMMA)
        } else {
            IB.inv..GAMMA <- IB.inv %*% GAMMA
            SX <- tcrossprod(IB.inv..GAMMA %*% cov.x, IB.inv..GAMMA)
        }
        SYX <- SX + SY
    } else {
        SYX <- SY
    }

    SYX
}

# compute E(ETA): expected value of latents
computeEETA.LISREL <- function(MLIST=NULL, x=NULL) {

    LAMBDA <- MLIST$lambda; nvar <- nrow(LAMBDA); nfac <- ncol(LAMBDA)
    ALPHA <- MLIST$alpha
    if(is.null(ALPHA)) {
        eeta <- numeric(nfac)    
    } else {
        BETA   <- MLIST$beta
        # beta?
        if(is.null(BETA)) {
            eeta <- ALPHA
        } else {
            tmp <- -BETA; nr <- nrow(BETA); i <- seq_len(nr);
            tmp[cbind(i, i)] <- 1
            IB.inv <- solve(tmp)
            eeta <- IB.inv %*% ALPHA
        }
    }

    # if GAMMA, also x part, and we return a matrix of (nfac x N)
    GAMMA <- MLIST$gamma
    if(!is.null(GAMMA)) {
        stopifnot(!is.null(x))
        if(is.null(BETA)) {
            EX <- x %*% t(GAMMA)
        } else {
            IB.inv..GAMMA <- IB.inv %*% GAMMA
            EX <- x %*% t(IB.inv..GAMMA)
        }
        EETA <- sweep(EX, MARGIN=2, STATS=eeta, FUN="+")
    } else {
        EETA <- eeta
    }

    EETA
}

# compute Sigma/ETA: variances/covariances of BOTH observed and latent variables
computeCOV.LISREL <- function(MLIST=NULL, cov.x=NULL, delta=TRUE) {

    LAMBDA <- MLIST$lambda; nvar <- nrow(LAMBDA)
    PSI    <- MLIST$psi;    nlat <- nrow(PSI)
    THETA  <- MLIST$theta
    BETA   <- MLIST$beta

    # 'extend' matrices
    LAMBDA2 <- rbind(LAMBDA, diag(nlat))
    THETA2  <- bdiag(THETA,  matrix(0,nlat,nlat))


    # beta?
    if(is.null(BETA)) {
        LAMBDA..IB.inv <- LAMBDA2
    } else {
        tmp <- -BETA; nr <- nrow(BETA); i <- seq_len(nr);
        tmp[cbind(i, i)] <- 1
        IB.inv <- solve(tmp)
        LAMBDA..IB.inv <- LAMBDA2 %*% IB.inv
    }

    # compute augment COV matrix
    COV <- tcrossprod(LAMBDA..IB.inv %*% PSI, LAMBDA..IB.inv) + THETA2

    # if delta, scale
    if(delta && !is.null(MLIST$delta)) {
        DELTA <- diag(MLIST$delta[,1L], nrow=nvar, ncol=nvar)
        COV[1:nvar,1:nvar] <- DELTA %*% COV[1:nvar,1:nvar] %*% DELTA
    }


    # if GAMMA, also x part
    GAMMA <- MLIST$gamma
    if(!is.null(GAMMA)) {
        stopifnot(!is.null(cov.x))
        if(is.null(BETA)) {
            SX <- tcrossprod(GAMMA %*% cov.x, GAMMA)
        } else {
            IB.inv..GAMMA <- IB.inv %*% GAMMA
            SX <- tcrossprod(IB.inv..GAMMA %*% cov.x, IB.inv..GAMMA)
        }
        COV[(nvar+1):(nvar+nlat),(nvar+1):(nvar+nlat)] <- 
            COV[(nvar+1):(nvar+nlat),(nvar+1):(nvar+nlat)] + SX
    }

    COV
}

# compute the *un*conditional variance of y: V(Y) or V(Y*)
computeVY.LISREL <- function(MLIST=NULL, cov.x=NULL, num.idx=NULL) {

    LAMBDA <- MLIST$lambda; nvar <- nrow(LAMBDA)
    PSI    <- MLIST$psi
    THETA  <- MLIST$theta
    BETA   <- MLIST$beta

    # beta?
    if(is.null(BETA)) {
        LAMBDA..IB.inv <- LAMBDA
    } else {
        tmp <- -BETA; nr <- nrow(BETA); i <- seq_len(nr);
        tmp[cbind(i, i)] <- 1
        IB.inv <- solve(tmp)
        LAMBDA..IB.inv <- LAMBDA %*% IB.inv
    }

    SY1 <- tcrossprod(LAMBDA..IB.inv %*% PSI, LAMBDA..IB.inv)

    # if TAU, we need to adjust the diagonal of THETA
    TAU <- MLIST$tau
    if(!is.null(TAU)) {
        if(!is.null(MLIST$delta)) {
            DELTA.inv2 <- 1/(MLIST$delta[,1L]*MLIST$delta[,1L]) 
        } else {
            DELTA.inv2 <- rep(1, nvar)
        }    
        THETA.diag <- DELTA.inv2 - diag(SY1)
        # but not for continuous 
        if(length(num.idx) > 0L)
            THETA.diag[num.idx] <- diag(THETA)[num.idx]
        # replace diagonal of THETA
        diag(THETA) <- THETA.diag
    }

    # compute Sigma Hat
    SY <- SY1 + THETA

    # if GAMMA, also x part
    GAMMA <- MLIST$gamma
    if(!is.null(GAMMA)) {
        stopifnot(!is.null(cov.x))
        LAMBDA..IB.inv..GAMMA <- LAMBDA..IB.inv %*% GAMMA
        SX <- tcrossprod(LAMBDA..IB.inv..GAMMA %*% cov.x, LAMBDA..IB.inv..GAMMA)
        SYX <- SX + SY
    } else {
        SYX <- SY
    }

    # variances only
    diag(SYX)
}

# derivative of the objective function
derivative.F.LISREL <- function(MLIST=NULL, Omega=NULL, Omega.mu=NULL) {

    LAMBDA <- MLIST$lambda
    PSI    <- MLIST$psi
    BETA   <- MLIST$beta
    ALPHA  <- MLIST$alpha 

    # beta?
    if(is.null(BETA)) {
        LAMBDA..IB.inv <- LAMBDA
    } else {
        tmp <- -BETA; nr <- nrow(BETA); i <- seq_len(nr);
        tmp[cbind(i, i)] <- 1
        IB.inv <- solve(tmp)
        LAMBDA..IB.inv <- LAMBDA %*% IB.inv
    }

    # meanstructure?
    meanstructure <- FALSE; if(!is.null(Omega.mu)) meanstructure <- TRUE

    # pre-compute some values
    tLAMBDA..IB.inv <- t(LAMBDA..IB.inv)
    if(!is.null(BETA)) {
        Omega..LAMBDA..IB.inv..PSI..tIB.inv <- 
            ( Omega %*% LAMBDA..IB.inv %*% PSI %*% t(IB.inv) )
    } else {
        Omega..LAMBDA <- Omega %*% LAMBDA
    }
    
    # 1. LAMBDA
    if(!is.null(BETA)) {
        if(meanstructure) {
            LAMBDA.deriv <- -1.0 * ( Omega.mu %*% t(ALPHA) %*% t(IB.inv) +
                                     Omega..LAMBDA..IB.inv..PSI..tIB.inv )
        } else {
            LAMBDA.deriv <- -1.0 * Omega..LAMBDA..IB.inv..PSI..tIB.inv
        }
    } else {
        # no BETA
        if(meanstructure) {
            LAMBDA.deriv <- -1.0 * ( Omega.mu %*% t(ALPHA) + 
                                     Omega..LAMBDA %*% PSI )
        } else {
            LAMBDA.deriv <- -1.0 * (Omega..LAMBDA %*% PSI)
        }
    }

    # 2. BETA
    if(!is.null(BETA)) {
        if(meanstructure) {
            BETA.deriv <- -1.0*(( t(IB.inv) %*% 
                                    (t(LAMBDA) %*% Omega.mu %*% t(ALPHA)) %*% 
                                  t(IB.inv)) +
                                (tLAMBDA..IB.inv %*% 
                                 Omega..LAMBDA..IB.inv..PSI..tIB.inv))
        } else {
            BETA.deriv <- -1.0 * ( tLAMBDA..IB.inv %*%
                                   Omega..LAMBDA..IB.inv..PSI..tIB.inv )
        }
    } else {
        BETA.deriv <- NULL
    }
   
    # 3. PSI 
    PSI.deriv <- -1.0 * ( tLAMBDA..IB.inv %*% Omega %*% LAMBDA..IB.inv )
    diag(PSI.deriv) <- 0.5 * diag(PSI.deriv)

    # 4. THETA
    THETA.deriv <- -1.0 * Omega
    diag(THETA.deriv) <- 0.5 * diag(THETA.deriv)

    if(meanstructure) {
        # 5. NU
        NU.deriv <- -1.0 * Omega.mu

        # 6. ALPHA
        ALPHA.deriv <- -1.0 * t( t(Omega.mu) %*% LAMBDA..IB.inv )
    } else {
        NU.deriv <- NULL
        ALPHA.deriv <- NULL
    }

    list(lambda = LAMBDA.deriv,
         beta   = BETA.deriv,
         theta  = THETA.deriv,
         psi    = PSI.deriv,
         nu     = NU.deriv,
         alpha  = ALPHA.deriv)
}

# dSigma/dx -- per model matrix
# note:
# we avoid using the duplication and elimination matrices
# for now (perhaps until we'll use the Matrix package)
derivative.sigma.LISREL <- function(m="lambda", 
                                    # all model matrix elements, or only a few?
                                    # NOTE: for symmetric matrices, 
                                    # we assume that the have full size 
                                    # (nvar*nvar) (but already correct for 
                                    # symmetry)
                                    idx=1:length(MLIST[[m]]),
                                    MLIST=NULL) {

    LAMBDA <- MLIST$lambda; nvar <- nrow(LAMBDA); nfac <- ncol(LAMBDA)
    PSI    <- MLIST$psi
 
    # only lower.tri part of sigma (not same order as elimination matrix?)
    v.idx <- lavaan:::vech.idx( nvar  ); pstar <- nvar*(nvar+1)/2

    # shortcut for gamma, nu, alpha and tau: empty matrix
    if(m == "nu" || m == "alpha" || m == "tau" || m == "gamma") {
        return( matrix(0.0, nrow=pstar, ncol=length(idx)) )
    }

    # Delta?
    delta.flag <- FALSE
    if(!is.null(MLIST$delta)) {
        DELTA <- MLIST$delta
        delta.flag <- TRUE
    }

    # beta?
    if(!is.null(MLIST$ibeta.inv)) {
        IB.inv <- MLIST$ibeta.inv
    } else if(!is.null(MLIST$beta)) {
        tmp <- -1.0 * MLIST$beta; diag(tmp) <- 1.0
        IB.inv <- solve(tmp)
    } else {
        IB.inv <- diag(nfac)
    }

    # pre
    if(m == "lambda" || m == "beta" || m == "delta") 
        IK <- diag(nvar^2) + commutationMatrix(nvar, nvar)
    if(m == "lambda" || m == "beta") {
        IB.inv..PSI..tIB.inv..tLAMBDA <-
            IB.inv %*% PSI %*% t(IB.inv) %*% t(LAMBDA)
    }
    if(m == "beta" || m == "psi") {
        LAMBDA..IB.inv <- LAMBDA %*% IB.inv
    }

    # here we go:
    if(m == "lambda") {
        DX <- IK %*% t(IB.inv..PSI..tIB.inv..tLAMBDA %x% diag(nvar))
        if(delta.flag)
             DX <- DX * as.numeric(DELTA %x% DELTA)
    } else if(m == "beta") {
        DX <- IK %*% ( t(IB.inv..PSI..tIB.inv..tLAMBDA) %x% LAMBDA..IB.inv )
        # this is not really needed (because we select idx=m.el.idx)
        DX[,lavaan:::diag.idx(nfac)] <- 0.0
        if(delta.flag) 
             DX <- DX * as.numeric(DELTA %x% DELTA)
    } else if(m == "psi") {
        DX <- (LAMBDA..IB.inv %x% LAMBDA..IB.inv) 
        # symmetry correction, but keeping all duplicated elements
        # since we depend on idx=m.el.idx
        # otherwise, we could simply postmultiply with the duplicationMatrix

        # we sum up lower.tri + upper.tri (but not the diagonal elements!)
        #imatrix <- matrix(1:nfac^2,nfac,nfac)
        #lower.idx <- imatrix[lower.tri(imatrix, diag=FALSE)]
        #upper.idx <- imatrix[upper.tri(imatrix, diag=FALSE)]
        lower.idx <- lavaan:::vech.idx(nfac, diag=FALSE)
        upper.idx <- lavaan:::vechru.idx(nfac, diag=FALSE)
        # NOTE YR: upper.idx (see 3 lines up) is wrong in MH patch!
        # fixed again 13/06/2012 after bug report of Mijke Rhemtulla.

        offdiagSum <- DX[,lower.idx] + DX[,upper.idx]
        DX[,c(lower.idx, upper.idx)] <- cbind(offdiagSum, offdiagSum)
        if(delta.flag)
            DX <- DX * as.numeric(DELTA %x% DELTA)
    } else if(m == "theta") {
        DX <- diag(nvar^2) # very sparse...
        # symmetry correction not needed, since all off-diagonal elements
        # are zero?
        if(delta.flag)
            DX <- DX * as.numeric(DELTA %x% DELTA)
    } else if(m == "delta") {
        Omega <- computeSigmaHat.LISREL(MLIST, delta=FALSE)
        DD <- diag(DELTA[,1], nvar, nvar)
        DD.Omega <- (DD %*% Omega)
        A <- DD.Omega %x% diag(nvar); B <- diag(nvar) %x% DD.Omega
        DX <- A[,lavaan:::diag.idx(nvar),drop=FALSE] + 
              B[,lavaan:::diag.idx(nvar),drop=FALSE]
    } else {
        stop("wrong model matrix names: ", m, "\n")
    }

    DX <- DX[v.idx, idx, drop=FALSE]
    DX
}

# dMu/dx -- per model matrix
derivative.mu.LISREL <- function(m="alpha", 
                                 # all model matrix elements, or only a few?
                                 idx=1:length(MLIST[[m]]), 
                                 MLIST=NULL) {


    LAMBDA <- MLIST$lambda; nvar <- nrow(LAMBDA); nfac <- ncol(LAMBDA)

    # shortcut for empty matrices
    if(m == "gamma" || m == "psi" || m == "theta" || 
       m == "tau" || m == "delta") {
        return( matrix(0.0, nrow=nvar, ncol=length(idx) ) )
    }

    # missing alpha
    if(is.null(MLIST$alpha)) 
        ALPHA <- matrix(0, nfac, 1L)
    else
        ALPHA  <- MLIST$alpha
 

    # beta?
    if(!is.null(MLIST$ibeta.inv)) {
        IB.inv <- MLIST$ibeta.inv
    } else if(!is.null(MLIST$beta)) {
        tmp <- -1.0 * MLIST$beta; diag(tmp) <- 1.0
        IB.inv <- solve(tmp)
    } else {
        IB.inv <- diag(nfac)
    }

    if(m == "nu") {
        DX <- diag(nvar)
    } else if(m == "lambda") {
        DX <- t(IB.inv %*% ALPHA) %x% diag(nvar)
    } else if(m == "beta") {
        DX <- t(IB.inv %*% ALPHA) %x% (LAMBDA %*% IB.inv)
        # this is not really needed (because we select idx=m.el.idx)
        DX[,lavaan:::diag.idx(nfac)] <- 0.0
    } else if(m == "alpha") {
        DX <- LAMBDA %*% IB.inv
    } else {
        stop("wrong model matrix names: ", m, "\n")
    }

    DX <- DX[, idx, drop=FALSE]
    DX
}

# dTh/dx -- per model matrix
derivative.th.LISREL <- function(m="tau",
                                 # all model matrix elements, or only a few?
                                 idx=1:length(MLIST[[m]]),
                                 th.idx=NULL,
                                 MLIST=NULL) {


    LAMBDA <- MLIST$lambda; nvar <- nrow(LAMBDA); nfac <- ncol(LAMBDA)
    TAU <- MLIST$tau; nth <- nrow(TAU)

    # missing alpha
    if(is.null(MLIST$alpha))
        ALPHA <- matrix(0, nfac, 1L)
    else
        ALPHA  <- MLIST$alpha

    # missing nu
    if(is.null(MLIST$nu))
        NU <- matrix(0, nvar, 1L)
    else
        NU <- MLIST$nu

    # Delta?
    delta.flag <- FALSE
    if(!is.null(MLIST$delta)) {
        DELTA <- MLIST$delta
        delta.flag <- TRUE
    }

    if(is.null(th.idx)) {
        th.idx <- 1:nth
        nlev <- rep(1L, nvar)
        K_nu <- diag(nvar)
    } else {
        nlev <- tabulate(th.idx, nbins=nvar); nlev[nlev == 0L] <- 1L
        K_nu <- matrix(0, sum(nlev), nvar)
        K_nu[ cbind(1:sum(nlev), rep(1:nvar, times=nlev)) ] <- 1.0
    }

    # shortcut for empty matrices
    if(m == "gamma" || m == "psi" || m == "theta") {
        return( matrix(0.0, nrow=length(th.idx), ncol=length(idx) ) )
    }

    # beta?
    if(!is.null(MLIST$ibeta.inv)) {
        IB.inv <- MLIST$ibeta.inv
    } else if(!is.null(MLIST$beta)) {
        tmp <- -1.0 * MLIST$beta; diag(tmp) <- 1.0
        IB.inv <- solve(tmp)
    } else {
        IB.inv <- diag(nfac)
    }

    if(m == "tau") {
        DX <- matrix(0, nrow=length(th.idx), ncol=nth)
        DX[ th.idx > 0L, ] <-  diag(nth)
        if(delta.flag)
            DX <- DX * as.numeric(K_nu %*% DELTA)
    } else if(m == "nu") {
        DX <- (-1) * K_nu
        if(delta.flag)
            DX <- DX * as.numeric(K_nu %*% DELTA)
    } else if(m == "lambda") {
        DX <- (-1) * t(IB.inv %*% ALPHA) %x% diag(nvar)
        DX <- K_nu %*% DX
        if(delta.flag)
            DX <- DX * as.numeric(K_nu %*% DELTA)
    } else if(m == "beta") {
        DX <- (-1) * t(IB.inv %*% ALPHA) %x% (LAMBDA %*% IB.inv)
        # this is not really needed (because we select idx=m.el.idx)
        DX[,lavaan:::diag.idx(nfac)] <- 0.0
        DX <- K_nu %*% DX
        if(delta.flag)
            DX <- DX * as.numeric(K_nu %*% DELTA)
    } else if(m == "alpha") {
        DX <- (-1) * LAMBDA %*% IB.inv
        DX <- K_nu %*% DX
        if(delta.flag)
            DX <- DX * as.numeric(K_nu %*% DELTA)
    } else if(m == "delta") {
        DX1 <- matrix(0, nrow=length(th.idx), ncol=1)
        DX1[ th.idx > 0L, ] <-  TAU
        DX2 <- NU + LAMBDA %*% IB.inv %*% ALPHA
        DX2 <- K_nu %*% DX2
        DX <- K_nu * as.numeric(DX1 - DX2)
    } else {
        stop("wrong model matrix names: ", m, "\n")
    }

    DX <- DX[, idx, drop=FALSE]
    DX
}

# dPi/dx -- per model matrix
derivative.pi.LISREL <- function(m="lambda",
                                 # all model matrix elements, or only a few?
                                 idx=1:length(MLIST[[m]]),
                                 MLIST=NULL) {

    LAMBDA <- MLIST$lambda; nvar <- nrow(LAMBDA); nfac <- ncol(LAMBDA)
    GAMMA <- MLIST$gamma; nexo <- ncol(GAMMA)

    # Delta?
    delta.flag <- FALSE
    if(!is.null(MLIST$delta)) {
        DELTA.diag <- MLIST$delta[,1L]
        delta.flag <- TRUE
    }

    # shortcut for empty matrices
    if(m == "tau" || m == "nu" || m == "alpha" || m == "psi" || m == "theta") {
        return( matrix(0.0, nrow=nvar*nexo, ncol=length(idx) ) )
    }

    # beta?
    if(!is.null(MLIST$ibeta.inv)) {
        IB.inv <- MLIST$ibeta.inv
    } else if(!is.null(MLIST$beta)) {
        tmp <- -1.0 * MLIST$beta; diag(tmp) <- 1.0
        IB.inv <- solve(tmp)
    } else {
        IB.inv <- diag(nfac)
    }

    if(m == "lambda") {
        DX <- t(IB.inv %*% GAMMA) %x% diag(nvar)
        if(delta.flag)
            DX <- DX * DELTA.diag
    } else if(m == "beta") {
        DX <- t(IB.inv %*% GAMMA) %x% (LAMBDA %*% IB.inv)
        # this is not really needed (because we select idx=m.el.idx)
        DX[,lavaan:::diag.idx(nfac)] <- 0.0
        if(delta.flag)
            DX <- DX * DELTA.diag
    } else if(m == "gamma") {
        DX <- diag(nexo) %x% (LAMBDA %*% IB.inv)
        if(delta.flag)
            DX <- DX * DELTA.diag
    } else if(m == "delta") {
        PRE <- rep(1, nexo) %x% diag(nvar)
        DX <- PRE * as.numeric(LAMBDA %*% IB.inv %*% GAMMA)
    } else {
        stop("wrong model matrix names: ", m, "\n")
    }

    DX <- DX[, idx, drop=FALSE]
    DX
}

TESTING_derivatives.LISREL <- function(MLIST = NULL, meanstructure=TRUE,
                                       th=FALSE, delta=FALSE, pi=FALSE) {

    if(is.null(MLIST)) {
        # create artificial matrices, compare 'numerical' vs 'analytical' 
        # derivatives
        #nvar <- 12; nfac <- 3; nexo <- 4 # this combination is special?
        nvar <- 20; nfac <- 6; nexo <- 5
        th.idx <- rep(seq_len(nvar), times=sample(c(1,1,2,6),nvar,replace=TRUE))
        # add some numeric
        num.idx <- sample(seq_len(nvar), ceiling(nvar/2))
        th.idx[ th.idx %in% num.idx & 
               !th.idx %in% th.idx[duplicated(th.idx)] ] <- 0
        nth <- sum(th.idx > 0L)

        MLIST <- list()
        MLIST$lambda <- matrix(0,nvar,nfac) 
        MLIST$beta   <- matrix(0,nfac,nfac)
        MLIST$theta  <- matrix(0,nvar,nvar)
        MLIST$psi    <- matrix(0,nfac,nfac)
        if(meanstructure) {
            MLIST$alpha  <- matrix(0,nfac,1L)
            MLIST$nu     <- matrix(0,nvar,1L)
        }
        if(th) MLIST$tau    <- matrix(0,nth,1L)
        if(delta) MLIST$delta  <- matrix(0,nvar,1L)
        MLIST$gamma <- matrix(0,nfac,nexo)

        # feed random numbers
        MLIST <- lapply(MLIST, function(x) {x[,] <- rnorm(length(x)); x})
        # fix
        diag(MLIST$beta) <- 0.0
        MLIST$psi[ lavaan:::vechru.idx(nfac) ] <-  
            MLIST$psi[ lavaan:::vech.idx(nfac) ]
        MLIST$theta[ lavaan:::vechru.idx(nvar) ] <-  
            MLIST$theta[ lavaan:::vech.idx(nvar) ]
        if(delta) MLIST$delta[,] <- abs(MLIST$delta)*10
    }

    compute.sigma <- function(x, mm="lambda", MLIST=NULL) {
        mlist <- MLIST
        if(mm %in% c("psi", "theta")) {
            mlist[[mm]] <- vech.reverse(x)
        } else {
            mlist[[mm]][,] <- x
        }
        vech(computeSigmaHat.LISREL(mlist))
    }

    compute.mu <- function(x, mm="lambda", MLIST=NULL) {
        mlist <- MLIST
        if(mm %in% c("psi", "theta")) {
            mlist[[mm]] <- vech.reverse(x)
        } else {
            mlist[[mm]][,] <- x
        }
        computeMuHat.LISREL(mlist)
    }

    compute.th2 <- function(x, mm="tau", MLIST=NULL, th.idx) {
        mlist <- MLIST
        if(mm %in% c("psi", "theta")) {
            mlist[[mm]] <- vech.reverse(x)
        } else {
            mlist[[mm]][,] <- x
        }
        computeTH.LISREL(mlist, th.idx=th.idx)
    }

    compute.pi <- function(x, mm="lambda", MLIST=NULL) {
        mlist <- MLIST
        if(mm %in% c("psi", "theta")) {
            mlist[[mm]] <- vech.reverse(x)
        } else {
            mlist[[mm]][,] <- x
        }
        computePI.LISREL(mlist)
    } 

    for(mm in names(MLIST)) {
        if(mm %in% c("psi", "theta")) {
            x <- lavaan:::vech(MLIST[[mm]])
        } else {
            x <- lavaan:::vec(MLIST[[mm]])
        }

        # 1. sigma
        DX1 <- lavJacobianC(func=compute.sigma, x=x, mm=mm, MLIST=MLIST)
        DX2 <- derivative.sigma.LISREL(m=mm, idx=1:length(MLIST[[mm]]),
                                       MLIST=MLIST)
        if(mm %in% c("psi","theta")) {
            # remove duplicated columns of symmetric matrices 
            idx <- lavaan:::vechru.idx(sqrt(ncol(DX2)), diag=FALSE)
            DX2 <- DX2[,-idx]
        }
        cat("[SIGMA] mm = ", sprintf("%-8s:", mm), "sum delta = ", 
            sprintf("%12.9f", sum(DX1-DX2)), "  max delta = ",
            sprintf("%12.9f", max(DX1-DX2)), "\n")

        # 2. mu
        DX1 <- lavJacobianC(func=compute.mu, x=x, mm=mm, MLIST=MLIST)
        DX2 <- derivative.mu.LISREL(m=mm, idx=1:length(MLIST[[mm]]),
                                       MLIST=MLIST)
        if(mm %in% c("psi","theta")) {
            # remove duplicated columns of symmetric matrices 
            idx <- lavaan:::vechru.idx(sqrt(ncol(DX2)), diag=FALSE)
            DX2 <- DX2[,-idx]
        }
        cat("[MU   ] mm = ", sprintf("%-8s:", mm), "sum delta = ",
            sprintf("%12.9f", sum(DX1-DX2)), "  max delta = ",
            sprintf("%12.9f", max(DX1-DX2)), "\n")

        # 3. th
        if(th) {
        DX1 <- lavJacobianC(func=compute.th2, x=x, mm=mm, MLIST=MLIST, 
                            th.idx=th.idx)
        DX2 <- derivative.th.LISREL(m=mm, idx=1:length(MLIST[[mm]]),
                                    MLIST=MLIST, th.idx=th.idx)
        if(mm %in% c("psi","theta")) {
            # remove duplicated columns of symmetric matrices 
            idx <- lavaan:::vechru.idx(sqrt(ncol(DX2)), diag=FALSE)
            DX2 <- DX2[,-idx]
        }
        cat("[TH   ] mm = ", sprintf("%-8s:", mm), "sum delta = ",
            sprintf("%12.9f", sum(DX1-DX2)), "  max delta = ",
            sprintf("%12.9f", max(DX1-DX2)), "\n")
        }

        # 4. pi
        if(pi) {
        DX1 <- lavJacobianC(func=compute.pi, x=x, mm=mm, MLIST=MLIST)
        DX2 <- derivative.pi.LISREL(m=mm, idx=1:length(MLIST[[mm]]),
                                    MLIST=MLIST)
        if(mm %in% c("psi","theta")) {
            # remove duplicated columns of symmetric matrices 
            idx <- lavaan:::vechru.idx(sqrt(ncol(DX2)), diag=FALSE)
            DX2 <- DX2[,-idx]
        }
        cat("[PI   ] mm = ", sprintf("%-8s:", mm), "sum delta = ",
            sprintf("%12.9f", sum(DX1-DX2)), "  max delta = ",
            sprintf("%12.9f", max(DX1-DX2)), "\n")
        }
    }

    MLIST
}


