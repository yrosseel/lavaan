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
    group.w.free  <- any(partable$lhs == "group" & partable$op == "%")
    gamma <- categorical

    # number of groups
    ngroups <- max(partable$group)
    ov.dummy.names.nox <- vector("list", ngroups)
    ov.dummy.names.x   <- vector("list", ngroups)
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
        ov.names.nox <- vnames(partable, "ov.nox",group=g)

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
            ov.dummy.names.nox[[g]] <-
                ov.names.nox[ ov.names.nox %in% dummy.names ]
            ov.dummy.names.x[[g]]   <- 
                ov.names.x[     ov.names.x %in% dummy.names ]

            # combine them, make sure order is identical to ov.names
            tmp <-  ov.names[ ov.names %in% dummy.names ]

            # extend lv.names
            lv.names <- c(lv.names, tmp)
            nfac <- length(lv.names)

            # add 'dummy' =~ entries
            dummy.mat <- rep("lambda", length(dummy.names))
        } else {
            ov.dummy.names.nox[[g]] <- character(0)
            ov.dummy.names.x[[g]]   <- character(0)
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

        # new 0.5-16: group weights
        idx <- which(target$group == g & target$lhs == "group" &
                     target$op == "%")
        tmp.mat[idx] <- "gw"
        tmp.row[idx] <- 1L
        tmp.col[idx] <- 1L

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
                           gw     = 1L,
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
                           gw     = 1L,
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
                               gw     = list( "group",     "weight"),
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
                                gw     = FALSE,
                                psi    = TRUE)
    
            # which mm's do we need? (always include lambda, theta and psi)
            mmNames <- c("lambda", "theta", "psi")
            if("beta" %in% tmp.mat) mmNames <- c(mmNames, "beta")
            if(meanstructure) mmNames <- c(mmNames, "nu", "alpha")
            if("tau" %in% tmp.mat) mmNames <- c(mmNames, "tau")
            if("delta" %in% tmp.mat) mmNames <- c(mmNames, "delta")
            if("gamma" %in% tmp.mat) mmNames <- c(mmNames, "gamma")
            if("gw" %in% tmp.mat) mmNames <- c(mmNames, "gw")

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

    # always add 'ov.dummy.*.names' attributes
    attr(REP, "ov.dummy.names.nox") <- ov.dummy.names.nox
    attr(REP, "ov.dummy.names.x")   <- ov.dummy.names.x

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

computeLAMBDA.LISREL <- function(MLIST=NULL) {
    return(MLIST$lambda)
}

# compute V(ETA): variances/covariances of latent variables
# - if no eXo (and GAMMA)
#     V(ETA) = (I-B)^-1 PSI (I-B)^-T
# - if eXo and GAMMA: (cfr lisrel submodel 3a with ksi=x)
#     V(ETA) = (I-B)^-1 [ GAMMA  cov.x t(GAMMA) + PSI] (I-B)^-T
computeVETA.LISREL <- function(MLIST=NULL, cov.x=NULL) {

    LAMBDA <- MLIST$lambda; nvar <- nrow(LAMBDA)
    PSI    <- MLIST$psi
    THETA  <- MLIST$theta
    BETA   <- MLIST$beta
    GAMMA  <- MLIST$gamma

    if(!is.null(GAMMA)) {
        stopifnot(!is.null(cov.x))
        # we treat 'x' as 'ksi' in the LISREL model; cov.x is PHI
        PSI <- tcrossprod(GAMMA %*% cov.x, GAMMA) + PSI
    }

    # beta?
    if(is.null(BETA)) {
        VETA <- PSI
    } else {
        tmp <- -BETA; nr <- nrow(BETA); i <- seq_len(nr);
        tmp[cbind(i, i)] <- 1
        IB.inv <- solve(tmp)
        VETA <- tcrossprod(IB.inv %*% PSI, IB.inv)
    }

    VETA
}

# compute V(ETA|x_i): variances/covariances of latent variables
#     V(ETA) = (I-B)^-1 PSI (I-B)^-T  + remove dummies
computeVETAx.LISREL <- function(MLIST=NULL, lv.dummy.idx=NULL) {

    PSI    <- MLIST$psi
    BETA   <- MLIST$beta

    if(!is.null(lv.dummy.idx)) {
        PSI  <-  PSI[-lv.dummy.idx, -lv.dummy.idx, drop=FALSE]
        BETA <- BETA[-lv.dummy.idx, -lv.dummy.idx, drop=FALSE]
        if(ncol(BETA) == 0L)
            BETA <- NULL
    }

    # beta?
    if(is.null(BETA)) {
        VETA <- PSI
    } else {
        tmp <- -BETA; nr <- nrow(BETA); i <- seq_len(nr);
        tmp[cbind(i, i)] <- 1
        IB.inv <- solve(tmp)
        VETA <- tcrossprod(IB.inv %*% PSI, IB.inv)
    }

    VETA
}

# compute E(ETA): expected value of latent variables
# - if no eXo (and GAMMA): 
#     E(ETA) = (I-B)^-1 ALPHA 
# - if eXo and GAMMA:
#     E(ETA) = (I-B)^-1 ALPHA + (I-B)^-1 GAMMA mean.x
computeEETA.LISREL <- function(MLIST=NULL, mean.x=NULL, 
                               sample.mean=NULL,
                               ov.y.dummy.ov.idx=NULL,
                               ov.x.dummy.ov.idx=NULL,
                               ov.y.dummy.lv.idx=NULL,
                               ov.x.dummy.lv.idx=NULL) {

    LAMBDA <- MLIST$lambda; nvar <- nrow(LAMBDA); nfac <- ncol(LAMBDA)
    BETA <- MLIST$beta; ALPHA <- MLIST$alpha; GAMMA <- MLIST$gamma

    ov.dummy.idx = c(ov.y.dummy.ov.idx, ov.x.dummy.ov.idx)
    lv.dummy.idx = c(ov.y.dummy.lv.idx, ov.x.dummy.lv.idx)

    if(is.null(ALPHA)) {
        if(length(ov.dummy.idx) > 0L) {
            eeta <- matrix(0, nfac, 1)
            # fill in exo values
            eeta[lv.dummy.idx] <- sample.mean[ov.dummy.idx]
 
            # Note: instead of sample.mean, we need 'intercepts'
            # sample.mean = NU + LAMBDA..IB.inv %*% ALPHA
            # so,
            # solve(LAMBDA..IB.inv) %*% (sample.mean - NU) = ALPHA
            # where
            # - LAMBDA..IB.inv only contains 'dummy' variables, and is square
            # - NU elements are not needed (since not in ov.dummy.idx)
            tmp <- -BETA; nr <- nrow(BETA); i <- seq_len(nr);
            tmp[cbind(i, i)] <- 1; IB.inv <- solve(tmp)
            LAMBDA..IB.inv <- LAMBDA %*% IB.inv
            LAMBDA..IB.inv.dummy <- LAMBDA..IB.inv[ov.dummy.idx, lv.dummy.idx]
            eeta[lv.dummy.idx] <- ( solve(LAMBDA..IB.inv.dummy) %*% 
                                    sample.mean[ov.dummy.idx] )
        } else { 
            eeta <- matrix(0, nfac, 1)
        }
    } else {
        eeta <- ALPHA
    }

    # IB.inv
    if(!is.null(BETA)) {
        tmp <- -BETA; nr <- nrow(BETA); i <- seq_len(nr);
        tmp[cbind(i, i)] <- 1
        IB.inv <- solve(tmp)
 
        eeta <- as.numeric(IB.inv %*% eeta)
        if(!is.null(GAMMA))
            eeta <- eeta + as.numeric( IB.inv %*% GAMMA %*% mean.x )
    } else {
        if(!is.null(GAMMA))
            eeta <- eeta + as.numeric( GAMMA %*% mean.x )
    }
        
    eeta
}

# compute E(ETA|x_i): conditional expected value of latent variable,
#                     given specific value of x_i
# - if no eXo (and GAMMA): 
#     E(ETA) = (I-B)^-1 ALPHA
#     we return a matrix of size [nobs x nfac] replicating E(ETA)
# - if eXo and GAMMA:
#     E(ETA|x_i) = (I-B)^-1 ALPHA + (I-B)^-1 GAMMA x_i
#     we return  a matrix of size [nobs x nfac]
#
# usually, exo = mean.x
# but if we change it, and sample.mean/ALPHA contains mean.x, we
# need to adapt
computeEETAx.LISREL <- function(MLIST=NULL, eXo=NULL, 
                                sample.mean=NULL,
                                ov.y.dummy.ov.idx=NULL,
                                ov.x.dummy.ov.idx=NULL,
                                ov.y.dummy.lv.idx=NULL,
                                ov.x.dummy.lv.idx=NULL) {

    LAMBDA <- MLIST$lambda; nvar <- nrow(LAMBDA); nfac <- ncol(LAMBDA)
    BETA <- MLIST$beta; ALPHA <- MLIST$alpha; GAMMA <- MLIST$gamma
    N <- nrow(eXo)

    ov.dummy.idx = c(ov.y.dummy.ov.idx, ov.x.dummy.ov.idx)
    lv.dummy.idx = c(ov.y.dummy.lv.idx, ov.x.dummy.lv.idx)

    if(is.null(ALPHA)) {
        if(length(ov.dummy.idx) > 0L) {
            eeta <- matrix(0, nfac, 1)
            eeta[lv.dummy.idx] <- sample.mean[ov.dummy.idx]

            # Note: instead of sample.mean, we need 'intercepts'
            # sample.mean = NU + LAMBDA..IB.inv %*% ALPHA
            # so,
            # solve(LAMBDA..IB.inv) %*% (sample.mean - NU) = ALPHA
            # where
            # - LAMBDA..IB.inv only contains 'dummy' variables, and is square
            # - NU elements are not needed (since not in ov.dummy.idx)
            tmp <- -BETA; nr <- nrow(BETA); i <- seq_len(nr);
            tmp[cbind(i, i)] <- 1; IB.inv <- solve(tmp)
            LAMBDA..IB.inv <- LAMBDA %*% IB.inv
            LAMBDA..IB.inv.dummy <- LAMBDA..IB.inv[ov.dummy.idx, lv.dummy.idx]
            eeta[lv.dummy.idx] <- ( solve(LAMBDA..IB.inv.dummy) %*%
                                    sample.mean[ov.dummy.idx] )
        } else {
            eeta <- matrix(0, nfac, 1)
        }
        EETA <- matrix(eeta, N, nfac, byrow=TRUE)
        if(length(ov.x.dummy.lv.idx) > 0L) {
            EETA[,ov.x.dummy.lv.idx] <- eXo
        }
    } else {
        EETA <- matrix(ALPHA, N, nfac, byrow=TRUE)
        if(length(ov.x.dummy.lv.idx) > 0L) {
            EETA[,ov.x.dummy.lv.idx] <- eXo
        }
    }

    # IB.inv
    if(!is.null(BETA)) {
        tmp <- -BETA; nr <- nrow(BETA); i <- seq_len(nr);
        tmp[cbind(i, i)] <- 1
        IB.inv <- solve(tmp)
 
        EETA <- EETA %*% t(IB.inv)
        if(!is.null(GAMMA))
            EETA <- EETA + eXo %*% t(IB.inv %*% GAMMA)
    } else {
        if(!is.null(GAMMA))
            EETA <- EETA + eXo %*% t(GAMMA)
    }
        
    EETA
}

# compute E(Y|x_i): conditional expected value of observed variable
#                     given specific value of x_i
#
# y*_i = nu + lambda eta_i + K x_i + epsilon_i
# 
# where eta_i = predict(fit) = factor scores
#
computeYHATx.LISREL <- function(MLIST=NULL, eXo=NULL, ETA=NULL,
                                sample.mean=NULL,
                                ov.y.dummy.ov.idx=NULL,
                                ov.x.dummy.ov.idx=NULL,
                                ov.y.dummy.lv.idx=NULL,
                                ov.x.dummy.lv.idx=NULL) {

    LAMBDA <- MLIST$lambda; nvar <- nrow(LAMBDA); nfac <- ncol(LAMBDA)
    BETA <- MLIST$beta; ALPHA <- MLIST$alpha; GAMMA <- MLIST$gamma
    DELTA <- MLIST$delta
    N <- nrow(ETA)

    ov.dummy.idx = c(ov.y.dummy.ov.idx, ov.x.dummy.ov.idx)
    lv.dummy.idx = c(ov.y.dummy.lv.idx, ov.x.dummy.lv.idx)

    # exogenous variables?
    if(is.null(eXo)) {
        nexo <- 0L
    } else {
        nexo <- ncol(eXo)
    }
    nvar <- nrow(MLIST$lambda)

    # fix NU
    NU <- MLIST$nu
    if(!is.null(NU)) {
        if(length(lv.dummy.idx) > 0L) {
            NU[ov.dummy.idx, 1L] <- MLIST$alpha[lv.dummy.idx, 1L]
        }
    } else {
        # if nexo == 0L, fill in unrestricted mean
        NU <- sample.mean
        # if nexo > 0, substract lambda %*% EETA
        if(nexo > 0L) {
            EETA <- computeEETA.LISREL(MLIST, mean.x=NULL,
                sample.mean=sample.mean,
                ov.y.dummy.ov.idx=ov.y.dummy.ov.idx,
                ov.x.dummy.ov.idx=ov.x.dummy.ov.idx,
                ov.y.dummy.lv.idx=ov.y.dummy.lv.idx,
                ov.x.dummy.lv.idx=ov.x.dummy.lv.idx)

            LAMBDA.X <- MLIST$lambda
            if(length(ov.y.dummy.ov.idx) > 0L) {
                LAMBDA.X[ov.y.dummy.ov.idx,] <-
                    MLIST$beta[ov.y.dummy.lv.idx,
                               ,drop=FALSE]
            }
            # 'regress' NU on X
            NU <- NU - LAMBDA.X %*% EETA
            NU[ov.x.dummy.ov.idx] <- sample.mean[ov.x.dummy.ov.idx]
        }
    }

    # fix LAMBDA
    LAMBDA <- MLIST$lambda
    if(length(lv.dummy.idx) > 0L) {
        LAMBDA <- LAMBDA[, -lv.dummy.idx, drop=FALSE]
        nfac <- ncol(LAMBDA)
        LAMBDA[ov.y.dummy.ov.idx,] <-
            MLIST$beta[ov.y.dummy.lv.idx,
                       1:nfac, drop=FALSE]
    }

    # compute YHAT
    YHAT <- sweep(ETA %*% t(LAMBDA), MARGIN=2, NU, "+")

    # Kappa + eXo?
    # note: Kappa elements are either in Gamma or in Beta
    if(nexo > 0L) {
        KAPPA <- matrix(0, nvar, nexo)
        if(!is.null(MLIST$gamma)) {
            KAPPA[ov.y.dummy.ov.idx,] <-
                MLIST$gamma[ov.y.dummy.lv.idx,,drop=FALSE]
        } else if(length(ov.x.dummy.ov.idx) > 0L) {
            KAPPA[ov.y.dummy.ov.idx,] <-
                MLIST$beta[ov.y.dummy.lv.idx,
                           ov.x.dummy.lv.idx, drop=FALSE]
        }

        # add fixed part
        YHAT <- YHAT + (eXo %*% t(KAPPA))

        # put back eXo
        if(length(ov.x.dummy.ov.idx) > 0L) {
            YHAT[, ov.x.dummy.ov.idx] <- eXo
        }
    }

    # delta?
    # FIXME: not used here?
    #if(!is.null(DELTA)) {
    #    YHAT <- sweep(YHAT, MARGIN=2, DELTA, "*")
    #}

    YHAT
}

# compute E(Y): expected value of observed
# E(Y) = nu + lambda * E(eta) 
# if delta -> E(Y) = delta * E(Y)
computeEY.LISREL <- function(MLIST=NULL, mean.x=NULL, sample.mean=NULL) {

    # E(Y) = nu + lambda * E(eta) 
    # if delta -> E(Y) = delta * E(Y)
    LAMBDA <- MLIST$lambda; nvar <- nrow(LAMBDA); nfac <- ncol(LAMBDA)
    NU <- MLIST$nu; ALPHA <- MLIST$alpha

    # compute E(ETA)
    EETA <- computeEETA.LISREL(MLIST = MLIST, mean.x = mean.x)

    # EY
    EY <- as.numeric(LAMBDA %*% EETA)

    # nu?
    if(!is.null(NU)) {
        EY <- EY + as.numeric(NU)
    } else {
        # use sample mean
        EY <- EY + sample.mean
    }

    # if delta, scale
    if(!is.null(MLIST$delta)) {
        EY <- EY * as.numeric(MLIST$delta)
    }

    EY
}


# deal with 'dummy' OV.X latent variables 
# create additional matrices (eg GAMMA), and resize
# remove all ov.x related entries
MLIST2MLISTX <- function(MLIST=NULL,
                         ov.x.dummy.ov.idx = NULL,
                         ov.x.dummy.lv.idx = NULL) {

    lv.idx <- ov.x.dummy.lv.idx
    ov.idx <- ov.x.dummy.ov.idx
    if(length(lv.idx) == 0L) return(MLIST)
    if(!is.null(MLIST$gamma)) {
        nexo <- ncol(MLIST$gamma)
    } else {
        nexo <- length(ov.x.dummy.ov.idx)
    }
    nvar <- nrow(MLIST$lambda)
    nfac <- ncol(MLIST$lambda) - length(lv.idx)

    # copy
    MLISTX <- MLIST

    # fix LAMBDA: 
    # - remove all ov.x related columns/rows
    MLISTX$lambda <- MLIST$lambda[-ov.idx, -lv.idx,drop=FALSE]

    # fix THETA:
    # - remove ov.x related columns/rows
    MLISTX$theta <- MLIST$theta[-ov.idx, -ov.idx, drop=FALSE]

    # fix PSI:
    # - remove ov.x related columns/rows
    MLISTX$psi <- MLIST$psi[-lv.idx, -lv.idx, drop=FALSE]

    # create GAMMA
    if(length(ov.x.dummy.lv.idx) > 0L) {
        MLISTX$gamma <- MLIST$beta[-lv.idx, lv.idx, drop=FALSE]
    }

    # fix BETA (remove if empty)
    if(!is.null(MLIST$beta)) {
        MLISTX$beta <- MLIST$beta[-lv.idx, -lv.idx, drop=FALSE]
        if(ncol(MLISTX$beta) == 0L) MLISTX$beta <- NULL
    }

    # fix NU
    if(!is.null(MLIST$nu)) {
        MLISTX$nu <- MLIST$nu[-ov.idx, 1L, drop=FALSE]
    }
    
    # fix ALPHA
    if(!is.null(MLIST$alpha)) {
        MLISTX$alpha <- MLIST$alpha[-lv.idx, 1L, drop=FALSE]
    }

    MLISTX
}


# create MLIST from MLISTX
MLISTX2MLIST <- function(MLISTX=NULL,
                         ov.x.dummy.ov.idx = NULL,
                         ov.x.dummy.lv.idx = NULL,
                         mean.x=NULL,
                         cov.x=NULL) {

    lv.idx <- ov.x.dummy.lv.idx; ndum <- length(lv.idx)
    ov.idx <- ov.x.dummy.ov.idx
    if(length(lv.idx) == 0L) return(MLISTX)
    stopifnot(!is.null(cov.x), !is.null(mean.x))
    nvar <- nrow(MLISTX$lambda); nfac <- ncol(MLISTX$lambda)

    # copy
    MLIST <- MLISTX

    # resize matrices
    MLIST$lambda <- rbind(cbind(MLISTX$lambda, matrix(0, nvar, ndum)),
                          matrix(0, ndum, nfac+ndum))
    MLIST$psi <- rbind(cbind(MLISTX$psi, matrix(0, nfac, ndum)),
                       matrix(0, ndum, nfac+ndum))
    MLIST$theta <- rbind(cbind(MLISTX$theta, matrix(0, nvar, ndum)),
                         matrix(0, ndum, nvar+ndum))
    if(!is.null(MLISTX$beta)) {
        MLIST$beta <- rbind(cbind(MLISTX$beta, matrix(0, nfac, ndum)),
                            matrix(0, ndum, nfac+ndum))
    }
    if(!is.null(MLISTX$alpha)) {
        MLIST$alpha <- rbind(MLISTX$alpha, matrix(0, ndum, 1))
    }
    if(!is.null(MLISTX$nu)) {
        MLIST$nu <- rbind(MLISTX$nu, matrix(0, ndum, 1))
    }

    # fix LAMBDA: 
    # - add columns for all dummy latent variables
    MLIST$lambda[ cbind(ov.idx, lv.idx) ] <- 1

    # fix PSI
    # - move cov.x elements to PSI
    MLIST$psi[lv.idx, lv.idx] <- cov.x

    # move (ov.x.dummy elements of) GAMMA to BETA
    MLIST$beta[1:nfac, ov.x.dummy.lv.idx] <- MLISTX$gamma
    MLIST$gamma <- NULL

    # fix ALPHA
    if(!is.null(MLIST$alpha)) {
        MLIST$alpha[lv.idx] <- mean.x
    }
    
    MLIST
}

# if DELTA parameterization, compute residual elements (in theta, or psi)
# of observed categorical variables, as a function of other model parameters
setResidualElements.LISREL <- function(MLIST=NULL,
                                       num.idx=NULL,
                                       ov.y.dummy.ov.idx=NULL,
                                       ov.y.dummy.lv.idx=NULL) {

    # remove num.idx from ov.y.dummy.*
    if(length(num.idx) > 0L && length(ov.y.dummy.ov.idx) > 0L) {
        n.idx <- which(ov.y.dummy.ov.idx %in% num.idx)
        if(length(n.idx) > 0L) {
            ov.y.dummy.ov.idx <- ov.y.dummy.ov.idx[-n.idx]
            ov.y.dummy.lv.idx <- ov.y.dummy.lv.idx[-n.idx]
        }
    }

    # force non-numeric theta elements to be zero
    if(length(num.idx) > 0L) {
        diag(MLIST$theta)[-num.idx] <- 0.0
    } else {
        diag(MLIST$theta) <- 0.0
    }
    if(length(ov.y.dummy.ov.idx) > 0L) {
        MLIST$psi[ov.y.dummy.lv.idx, ov.y.dummy.lv.idx] <- 0.0
    }

    # special case: PSI=0, and lambda=I (eg ex3.12)
    if(ncol(MLIST$psi) > 0L && 
       sum(diag(MLIST$psi)) == 0.0 && all(diag(MLIST$lambda) == 1)) {
        ### FIXME: more elegant/general solution??
        diag(MLIST$psi) <- 1
        Sigma.hat <- computeSigmaHat.LISREL(MLIST = MLIST, delta=FALSE)
        diag.Sigma <- diag(Sigma.hat) - 1.0
    } else if(ncol(MLIST$psi) == 0L) {
        diag.Sigma <- rep(0, ncol(MLIST$theta))
    } else {
        Sigma.hat <- computeSigmaHat.LISREL(MLIST = MLIST, delta=FALSE)
        diag.Sigma <- diag(Sigma.hat)
    }

    if(is.null(MLIST$delta)) {
        delta <- rep(1, length(diag.Sigma))
    } else {
        delta <- MLIST$delta
    }
    # theta = DELTA^(-1/2) - diag( LAMBDA (I-B)^-1 PSI (I-B)^-T t(LAMBDA) )
    RESIDUAL <- as.numeric(1/delta^2 - diag.Sigma)
    if(length(num.idx) > 0L) {
        diag(MLIST$theta)[-num.idx] <- RESIDUAL[-num.idx]
    } else {
        diag(MLIST$theta) <- RESIDUAL
    }

    # move ov.y.dummy elements from THETA to PSI
    if(length(ov.y.dummy.ov.idx) > 0L) {
        MLIST$psi[ov.y.dummy.lv.idx, ov.y.dummy.lv.idx] <- 
            MLIST$theta[ov.y.dummy.ov.idx, ov.y.dummy.ov.idx]
        MLIST$theta[ov.y.dummy.ov.idx, ov.y.dummy.ov.idx] <- 0.0
    }

    MLIST
}

# if THETA parameterization, compute delta elements 
# of observed categorical variables, as a function of other model parameters
setDeltaElements.LISREL <- function(MLIST=NULL, num.idx=NULL) {

    Sigma.hat <- computeSigmaHat.LISREL(MLIST = MLIST, delta=FALSE)
    diag.Sigma <- diag(Sigma.hat)

    # (1/delta^2) = diag( LAMBDA (I-B)^-1 PSI (I-B)^-T t(LAMBDA) ) + THETA
    #tmp <- diag.Sigma + THETA
    tmp <- diag.Sigma
    MLIST$delta[, 1L] <- sqrt(1/tmp)

    # numeric delta's stay 1.0
    if(length(num.idx) > 0L) {
        MLIST$delta[num.idx] <- 1.0
    }

    MLIST
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

    # group weight?
    group.w.free <- FALSE; if(!is.null(MLIST$gw)) group.w.free <- TRUE

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

    if(group.w.free) {
        GROUP.W.deriv <- 0.0
    } else {
        GROUP.W.deriv <- NULL
    }

    list(lambda = LAMBDA.deriv,
         beta   = BETA.deriv,
         theta  = THETA.deriv,
         psi    = PSI.deriv,
         nu     = NU.deriv,
         alpha  = ALPHA.deriv,
         gw     = GROUP.W.deriv)
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
                                    MLIST=NULL,
                                    delta = TRUE) {

    LAMBDA <- MLIST$lambda; nvar <- nrow(LAMBDA); nfac <- ncol(LAMBDA)
    PSI    <- MLIST$psi
 
    # only lower.tri part of sigma (not same order as elimination matrix?)
    v.idx <- vech.idx( nvar  ); pstar <- nvar*(nvar+1)/2

    # shortcut for gamma, nu, alpha and tau: empty matrix
    if(m == "nu" || m == "alpha" || m == "tau" || m == "gamma" || m == "gw") {
        return( matrix(0.0, nrow=pstar, ncol=length(idx)) )
    }

    # Delta?
    delta.flag <- FALSE
    if(delta && !is.null(MLIST$delta)) {
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
        DX[,diag.idx(nfac)] <- 0.0
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
        lower.idx <- vech.idx(nfac, diagonal=FALSE)
        upper.idx <- vechru.idx(nfac, diagonal=FALSE)
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
        DX <- A[,diag.idx(nvar),drop=FALSE] + 
              B[,diag.idx(nvar),drop=FALSE]
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
       m == "tau" || m == "delta"|| m == "gw") {
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
        DX[,diag.idx(nfac)] <- 0.0
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
                                 MLIST=NULL,
                                 delta = TRUE) {


    LAMBDA <- MLIST$lambda; nvar <- nrow(LAMBDA); nfac <- ncol(LAMBDA)
    TAU <- MLIST$tau; nth <- nrow(TAU)

    # missing alpha
    if(is.null(MLIST$alpha)) {
        ALPHA <- matrix(0, nfac, 1L)
    } else {
        ALPHA  <- MLIST$alpha
    }

    # missing nu
    if(is.null(MLIST$nu)) {
        NU <- matrix(0, nvar, 1L)
    } else {
        NU <- MLIST$nu
    }

    # Delta?
    delta.flag <- FALSE
    if(delta && !is.null(MLIST$delta)) {
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
    if(m == "gamma" || m == "psi" || m == "theta" || m == "gw") {
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
        DX[,diag.idx(nfac)] <- 0.0
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
    if(m == "tau" || m == "nu" || m == "alpha" || m == "psi" || 
       m == "theta" || m == "gw") {
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
        DX[,diag.idx(nfac)] <- 0.0
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

# dGW/dx -- per model matrix
derivative.gw.LISREL <- function(m="gw", 
                                 # all model matrix elements, or only a few?
                                 idx=1:length(MLIST[[m]]), 
                                 MLIST=NULL) {

    # shortcut for empty matrices
    if(m != "gw") {
        return( matrix(0.0, nrow=1L, ncol=length(idx) ) )
    } else {
        # m == "gw"
        DX <- matrix(1.0, 1, 1)
    }

    DX <- DX[, idx, drop=FALSE]
    DX
}


# MLIST = NULL; meanstructure=TRUE; th=TRUE; delta=TRUE; pi=TRUE; gw=FALSE
# vech.idx <- lavaan:::vech.idx; vechru.idx <- lavaan:::vechru.idx
# vec <- lavaan:::vec; lavJacobianC <- lavaan:::lavJacobianC
# computeSigmaHat.LISREL <- lavaan:::computeSigmaHat.LISREL
# setDeltaElements.LISREL <- lavaan:::setDeltaElements.LISREL
TESTING_derivatives.LISREL <- function(MLIST = NULL,
                                       nvar = NULL, nfac = NULL, nexo = NULL,
                                       th.idx = NULL, num.idx = NULL,
                                       meanstructure = TRUE,
                                       th = TRUE, delta = TRUE, pi = TRUE,
                                       gw = FALSE, theta = FALSE,
                                       debug = FALSE) {

    if(is.null(MLIST)) {
        # create artificial matrices, compare 'numerical' vs 'analytical' 
        # derivatives
        #nvar <- 12; nfac <- 3; nexo <- 4 # this combination is special?
        if(is.null(nvar)) {
            nvar <- 20
        }
        if(is.null(nfac)) {
            nfac <- 6
        }
        if(is.null(nexo)) {
            nexo <- 5
        }
        if(is.null(num.idx)) {
            num.idx <- sort(sample(seq_len(nvar), ceiling(nvar/2)))
        }
        if(is.null(th.idx)) {
            th.idx <- integer(0L)
            for(i in 1:nvar) {
                if(i %in% num.idx) {
                    th.idx <- c(th.idx, 0)
                } else {
                    th.idx <- c(th.idx, rep(i, sample(c(1,1,2,6), 1L)))
                }
            }
        }
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
        if(gw) MLIST$gw <- matrix(0, 1L, 1L)

        # feed random numbers
        MLIST <- lapply(MLIST, function(x) {x[,] <- rnorm(length(x)); x})
        # fix
        diag(MLIST$beta) <- 0.0
        diag(MLIST$theta) <- diag(MLIST$theta)^2 * 10 
        diag(MLIST$psi)   <- diag(MLIST$psi)^2 * 10
        MLIST$psi[ vechru.idx(nfac) ] <-  
            MLIST$psi[ vech.idx(nfac) ]
        MLIST$theta[ vechru.idx(nvar) ] <-  
            MLIST$theta[ vech.idx(nvar) ]
        if(delta) MLIST$delta[,] <- abs(MLIST$delta)*10
    } else {
        nvar <- nrow(MLIST$lambda)
    }

    compute.sigma <- function(x, mm="lambda", MLIST=NULL) {
        mlist <- MLIST
        if(mm %in% c("psi", "theta")) {
            mlist[[mm]] <- vech.reverse(x)
        } else {
            mlist[[mm]][,] <- x
        }
        if(theta) {
            mlist <- setDeltaElements.LISREL(MLIST = mlist, num.idx = num.idx)
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
        if(theta) {
            mlist <- setDeltaElements.LISREL(MLIST = mlist, num.idx = num.idx)
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
        if(theta) {
            mlist <- setDeltaElements.LISREL(MLIST = mlist, num.idx = num.idx)
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
        if(theta) {
            mlist <- setDeltaElements.LISREL(MLIST = mlist, num.idx = num.idx)
        }
        computePI.LISREL(mlist)
    }

    compute.gw <- function(x, mm="gw", MLIST=NULL) {
        mlist <- MLIST
        if(mm %in% c("psi", "theta")) {
            mlist[[mm]] <- vech.reverse(x)
        } else {
            mlist[[mm]][,] <- x
        }
        if(theta) {
            mlist <- setDeltaElements.LISREL(MLIST = mlist, num.idx = num.idx)
        }
        mlist$gw[1,1]
    }

    # if theta, set MLIST$delta
    if(theta) {
        MLIST <- setDeltaElements.LISREL(MLIST = MLIST, num.idx = num.idx)
    }

    for(mm in names(MLIST)) {
        if(mm %in% c("psi", "theta")) {
            x <- vech(MLIST[[mm]])
        } else {
            x <- vec(MLIST[[mm]])
        }
        if(mm == "delta" && theta) next
        if(debug) {
            cat("### mm = ", mm, "\n")
        }

        # 1. sigma
        DX1 <- lavJacobianC(func=compute.sigma, x=x, mm=mm, MLIST=MLIST)
        DX2 <- derivative.sigma.LISREL(m=mm, idx=1:length(MLIST[[mm]]),
                                       MLIST=MLIST, delta = !theta)
        if(mm %in% c("psi","theta")) {
            # remove duplicated columns of symmetric matrices 
            idx <- vechru.idx(sqrt(ncol(DX2)), diagonal=FALSE)
            if(length(idx) > 0L) DX2 <- DX2[,-idx]
        }
        if(theta) {
            sigma.hat <- computeSigmaHat.LISREL(MLIST=MLIST, delta=FALSE)
            R <- lav_deriv_cov2cor(sigma.hat, num.idx = num.idx)
            
            DX3 <- DX2
            DX2 <- R %*% DX2
        }
        if(debug) {
            cat("[SIGMA] mm = ", sprintf("%-8s:", mm), "DX1 (numerical):\n");
            print(zapsmall(DX1)); cat("\n")
            cat("[SIGMA] mm = ", sprintf("%-8s:", mm), "DX2 (analytical):\n");
            print(DX2); cat("\n")
            cat("[SIGMA] mm = ", sprintf("%-8s:", mm), "DX3 (analytical):\n");
            print(DX3); cat("\n")
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
            idx <- vechru.idx(sqrt(ncol(DX2)), diagonal=FALSE)
            if(length(idx) > 0L) DX2 <- DX2[,-idx]
        }
        cat("[MU   ] mm = ", sprintf("%-8s:", mm), "sum delta = ",
            sprintf("%12.9f", sum(DX1-DX2)), "  max delta = ",
            sprintf("%12.9f", max(DX1-DX2)), "\n")
        if(debug) {
            cat("[MU   ] mm = ", sprintf("%-8s:", mm), "DX1 (numerical):\n");
            print(zapsmall(DX1)); cat("\n")
            cat("[MU   ] mm = ", sprintf("%-8s:", mm), "DX2 (analytical):\n");
            print(DX2); cat("\n")
        }

        # 3. th
        if(th) {
            DX1 <- lavJacobianC(func=compute.th2, x=x, mm=mm, MLIST=MLIST, 
                                th.idx=th.idx)
            DX2 <- derivative.th.LISREL(m=mm, idx=1:length(MLIST[[mm]]),
                                        MLIST=MLIST, th.idx=th.idx,
                                        delta=TRUE)
            if(theta) {
                # 1. compute dDelta.dx
                dxSigma <- 
                    derivative.sigma.LISREL(m=mm, idx=1:length(MLIST[[mm]]),
                                            MLIST=MLIST, delta = !theta)
                var.idx <- which(!lavaan:::vech.idx(nvar) %in% 
                                  lavaan:::vech.idx(nvar, diag=FALSE))
                sigma.hat <- computeSigmaHat.LISREL(MLIST=MLIST, delta=FALSE)
                dsigma <- diag(sigma.hat)
                # dy/ddsigma = -0.5/(ddsigma*sqrt(ddsigma))
                dDelta.dx <- dxSigma[var.idx,] * -0.5 / (dsigma*sqrt(dsigma))

                # 2. compute dth.dDelta
                dth.dDelta <- 
                    derivative.th.LISREL(m="delta", 
                                         idx=1:length(MLIST[["delta"]]),
                                         MLIST=MLIST, th.idx=th.idx)

                # 3. add dth.dDelta %*% dDelta.dx
                no.num.idx <- which(th.idx > 0)
                DX2[no.num.idx,] <- DX2[no.num.idx,,drop=FALSE] + 
                                    (dth.dDelta %*% dDelta.dx)[no.num.idx,,drop=FALSE]
                #DX2 <- DX2 + dth.dDelta %*% dDelta.dx
            }
            if(mm %in% c("psi","theta")) {
                # remove duplicated columns of symmetric matrices 
                idx <- vechru.idx(sqrt(ncol(DX2)), diagonal=FALSE)
                if(length(idx) > 0L) DX2 <- DX2[,-idx]
            }
            cat("[TH   ] mm = ", sprintf("%-8s:", mm), "sum delta = ",
                sprintf("%12.9f", sum(DX1-DX2)), "  max delta = ",
                sprintf("%12.9f", max(DX1-DX2)), "\n")
            if(debug) {
                cat("[TH   ] mm = ",sprintf("%-8s:", mm),"DX1 (numerical):\n")
                print(zapsmall(DX1)); cat("\n")
                cat("[TH   ] mm = ",sprintf("%-8s:", mm),"DX2 (analytical):\n")
                print(DX2); cat("\n")
            }
        }

        # 4. pi
        if(pi) {
            DX1 <- lavJacobianC(func=compute.pi, x=x, mm=mm, MLIST=MLIST)
            DX2 <- derivative.pi.LISREL(m=mm, idx=1:length(MLIST[[mm]]),
                                        MLIST=MLIST)
            if(mm %in% c("psi","theta")) {
                # remove duplicated columns of symmetric matrices 
                idx <- vechru.idx(sqrt(ncol(DX2)), diagonal=FALSE)
                if(length(idx) > 0L) DX2 <- DX2[,-idx]
            }
            if(theta) {
                # 1. compute dDelta.dx
                dxSigma <- 
                    derivative.sigma.LISREL(m=mm, idx=1:length(MLIST[[mm]]),
                                            MLIST=MLIST, delta = !theta)
                if(mm %in% c("psi","theta")) {
                    # remove duplicated columns of symmetric matrices 
                    idx <- vechru.idx(sqrt(ncol(dxSigma)), diagonal=FALSE)
                    if(length(idx) > 0L) dxSigma <- dxSigma[,-idx]
                }
                var.idx <- which(!lavaan:::vech.idx(nvar) %in% 
                                  lavaan:::vech.idx(nvar, diag=FALSE))
                sigma.hat <- computeSigmaHat.LISREL(MLIST=MLIST, delta=FALSE)
                dsigma <- diag(sigma.hat)
                # dy/ddsigma = -0.5/(ddsigma*sqrt(ddsigma))
                dDelta.dx <- dxSigma[var.idx,] * -0.5 / (dsigma*sqrt(dsigma))

                # 2. compute dpi.dDelta
                dpi.dDelta <- 
                    derivative.pi.LISREL(m="delta", 
                                         idx=1:length(MLIST[["delta"]]),
                                         MLIST=MLIST)

                # 3. add dpi.dDelta %*% dDelta.dx
                no.num.idx <- which(! seq.int(1L, nvar) %in% num.idx )
                no.num.idx <- rep(seq.int(0,nexo-1) * nvar, 
                                  each=length(no.num.idx)) + no.num.idx
                DX2[no.num.idx,] <- DX2[no.num.idx,,drop=FALSE] + 
                                    (dpi.dDelta %*% dDelta.dx)[no.num.idx,,drop=FALSE]
            }
            cat("[PI   ] mm = ", sprintf("%-8s:", mm), "sum delta = ",
                sprintf("%12.9f", sum(DX1-DX2)), "  max delta = ",
                sprintf("%12.9f", max(DX1-DX2)), "\n")
            if(debug) {
                cat("[PI   ] mm = ",sprintf("%-8s:", mm),"DX1 (numerical):\n")
                print(zapsmall(DX1)); cat("\n")
                cat("[PI   ] mm = ",sprintf("%-8s:", mm),"DX2 (analytical):\n")
                print(DX2); cat("\n")
            }
        }

        # 5. gw
        if(gw) {
            DX1 <- lavJacobianC(func=compute.gw, x=x, mm=mm, MLIST=MLIST)
            DX2 <- derivative.gw.LISREL(m=mm, idx=1:length(MLIST[[mm]]),
                                    MLIST=MLIST)
            if(mm %in% c("psi","theta")) {
                # remove duplicated columns of symmetric matrices 
                idx <- vechru.idx(sqrt(ncol(DX2)), diagonal=FALSE)
                if(length(idx) > 0L) DX2 <- DX2[,-idx]
            }
            cat("[GW   ] mm = ", sprintf("%-8s:", mm), "sum delta = ",
                sprintf("%12.9f", sum(DX1-DX2)), "  max delta = ",
                sprintf("%12.9f", max(DX1-DX2)), "\n")
            if(debug) {
                cat("[GW   ] mm = ",sprintf("%-8s:", mm),"DX1 (numerical):\n")
                print(DX1); cat("\n\n")
                cat("[GW   ] mm = ",sprintf("%-8s:", mm),"DX2 (analytical):\n")
                print(DX2); cat("\n\n")
            }
        }
    }

    MLIST$th.idx <- th.idx
    MLIST$num.idx <- num.idx

    MLIST
}


