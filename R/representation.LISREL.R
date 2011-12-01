# matrix representation of latent variable models
# and matrix-representation specific functions:
# - computeSigmaHat
# - computeMuHat
# - derivative.F

# initital version: YR 2011-01-21: LISREL stuff
# updates:          YR 2011-12-01: group specific extraction

representation.LISREL <- function(user=NULL, target=NULL, 
                                  extra=FALSE) {

    # prepare target list
    if(is.null(target)) target <- user

    # prepare output
    N <- length(target$lhs)
    tmp.mat <- character(N); tmp.row <- integer(N); tmp.col <- integer(N)

    # global settings
    meanstructure <- any(user$op == "~1")

    # number of groups
    ngroups <- max(user$group)
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
        ov.names <- vnames(user, "ov", group=g); nvar <- length(ov.names)
        lv.names <- vnames(user, "lv", group=g); nfac <- length(lv.names)

        # in this representation, we need to create 'phantom/dummy' latent 
        # variables for all `x' and `y' variables not in lv.names
        tmp.names <- 
            unique( c(user$lhs[user$op == "~" & user$group == g],
                      user$rhs[user$op == "~" & user$group == g]) )
        dummy.names <- tmp.names[ !tmp.names %in% lv.names ]
        if(length(dummy.names)) {
            # make sure order is the same as ov.names
            ov.dummy.names[[g]] <- ov.names[ ov.names %in% dummy.names ]

            # extend lv.names
            lv.names <- c(lv.names, ov.dummy.names)
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
        idx <- which(target$group == g & target$op == "~")
        tmp.mat[idx] <- "beta"
        tmp.row[idx] <- match(target$lhs[idx], lv.names)
        tmp.col[idx] <- match(target$rhs[idx], lv.names)

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

        if(extra) {
            # mRows
            mmRows <- list(nu     = nvar,
                           lambda = nvar,
                           theta  = nvar,
                           alpha  = nfac,
                           beta   = nfac,
                           psi    = nfac)

            # mCols
            mmCols <- list(nu     = 1L,
                           lambda = nfac,
                           theta  = nvar,
                           alpha  = 1L,
                           beta   = nfac,
                           psi    = nfac)

            # dimNames for LISREL model matrices
            mmDimNames <- list(nu     = list( ov.names, "intercept"),
                               lambda = list( ov.names,    lv.names),
                               theta  = list( ov.names,    ov.names),
                               alpha  = list( lv.names, "intercept"),
                               beta   = list( lv.names,    lv.names),
                               psi    = list( lv.names,    lv.names))
    
            # isSymmetric
            mmSymmetric <- list(nu     = FALSE,
                                lambda = FALSE,
                                theta  = TRUE,
                                alpha  = FALSE,
                                beta   = FALSE,
                                psi    = TRUE)
    
            # which mm's do we need? (always include lambda, theta and psi)
            mmNames <- c("lambda", "theta", "psi")
            if("beta" %in% tmp.mat) mmNames <- c(mmNames, "beta")
            if(meanstructure) mmNames <- c(mmNames, "nu", "alpha")

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
computeSigmaHat.LISREL <- function(MLIST=NULL) {

    LAMBDA <- MLIST$lambda
    PSI    <- MLIST$psi
    THETA  <- MLIST$theta
    BETA   <- MLIST$beta

    # beta?
    if(is.null(BETA)) {
        LAMBDA..IB.inv <- LAMBDA
    } else {
        tmp <- -1.0 * BETA; diag(tmp) <- 1.0
        IB.inv <- solve(tmp)
        LAMBDA..IB.inv <- LAMBDA %*% IB.inv
    }

    # compute Sigma Hat
    Sigma.hat <- ( LAMBDA..IB.inv %*% PSI %*% t(LAMBDA..IB.inv) + THETA )

    Sigma.hat
}

# compute MuHat for a single group
computeMuHat.LISREL <- function(MLIST=NULL) {

    NU     <- MLIST$nu
    ALPHA  <- MLIST$alpha
    LAMBDA <- MLIST$lambda
    BETA   <- MLIST$beta

    # beta?
    if(is.null(BETA)) {
        LAMBDA..IB.inv <- LAMBDA
    } else {
        tmp <- -1.0 * BETA; diag(tmp) <- 1.0
        IB.inv <- solve(tmp)
        LAMBDA..IB.inv <- LAMBDA %*% IB.inv
    }
    
    # compute Mu Hat
    Mu.hat <- NU + LAMBDA..IB.inv %*% ALPHA

    Mu.hat
}


# derivative of the objective function
derivative.F.LISREL <- function(MLIST=NULL, Omega=NULL, Omega.mu=NULL) {

    LAMBDA <- MLIST$lambda
    PSI    <- MLIST$psi
    THETA  <- MLIST$theta
    BETA   <- MLIST$beta
    NU     <- MLIST$nu
    ALPHA  <- MLIST$alpha 

    # beta?
    if(is.null(BETA)) {
        LAMBDA..IB.inv <- LAMBDA
    } else {
        tmp <- -1.0 * BETA; diag(tmp) <- 1.0
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
derivative.sigma.LISREL <- function(m="lambda", idx=NULL, MLIST=NULL) {

    LAMBDA <- MLIST$lambda; nvar <- nrow(LAMBDA); nfac <- ncol(LAMBDA)
    PSI    <- MLIST$psi
    THETA  <- MLIST$theta
 
    # only lower.tri part of sigma (not same order as elimination matrix?)
    v.idx <- vech.idx( nvar  ); pstar <- nvar*(nvar+1)/2

    # all model matrix elements, or only a few?
    # NOTE: for symmetric matrices, we assume that the have full 
    #       size (nvar*nvar) (but already correct for symmetry)
    if(is.null(idx)) idx <- 1:length(MLIST[[m]])

    # shortcut for nu and alpha: empty matrix
    if(m == "nu" || m == "alpha") {
        return( matrix(0.0, nrow=pstar, ncol=length(idx)) )
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
    if(m == "lambda" || m == "beta") {
        IK <- diag(nvar^2) + commutationMatrix(nvar, nvar)
        IB.inv..PSI..tIB.inv..tLAMBDA <-
            IB.inv %*% PSI %*% t(IB.inv) %*% t(LAMBDA)
    }
    if(m == "beta" || m == "psi") {
        LAMBDA..IB.inv <- LAMBDA %*% IB.inv
    }

    # here we go:
    if(m == "lambda") {
        DX <- IK %*% t(IB.inv..PSI..tIB.inv..tLAMBDA %x% diag(nvar))
    } else if(m == "beta") {
        DX <- IK %*% ( t(IB.inv..PSI..tIB.inv..tLAMBDA) %x% LAMBDA..IB.inv )
    } else if(m == "psi") {
        DX <- (LAMBDA..IB.inv %x% LAMBDA..IB.inv) 
        # symmetry correction, but keeping all duplicated elements
        # since we depend on idx=m.el.idx
        # otherwise, we could simply postmultiply with the duplicationMatrix

        ## FIXME!! find a more elegant solution here...
        # we sum up lower.tri + upper.tri (but not the diagonal elements!)
        imatrix <- matrix(1:nfac^2,nfac,nfac)
        lower.idx <- imatrix[lower.tri(t(imatrix), diag=FALSE)]
        upper.idx <- t(imatrix)[lower.tri(t(imatrix), diag=FALSE)]
        tmp <- DX[,lower.idx]
        DX[,lower.idx] <- DX[,lower.idx] + DX[,upper.idx]
        DX[,upper.idx] <- DX[,upper.idx] + tmp
    } else if(m == "theta") {
        DX <- diag(nvar^2) # very sparse...
        # symmetry correction not needed, since all off-diagonal elements
        # are zero?
    } else {
        stop("wrong model matrix names: ", m, "\n")
    }

    DX <- DX[v.idx, idx, drop=FALSE]
    DX
}

# dMu/dx -- per model matrix
derivative.mu.LISREL <- function(m="alpha", idx=NULL, MLIST=NULL) {

    LAMBDA <- MLIST$lambda; nvar <- nrow(LAMBDA); nfac <- ncol(LAMBDA)

    # all model matrix elements, or only a few?
    if(is.null(idx)) idx <- 1:length(MLIST[[m]])

    # shortcut for empty matrices
    if(m == "psi" || m == "theta") {
        return( matrix(0.0, nrow=nvar, ncol=length(idx) ) )
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
    NU     <- MLIST$nu
    ALPHA  <- MLIST$alpha

    if(m == "nu") {
        DX <- diag(nvar)
    } else if(m == "lambda") {
        DX <- t(IB.inv %*% ALPHA) %x% diag(nvar)
    } else if(m == "beta") {
        DX <- t(IB.inv %*% ALPHA) %x% (LAMBDA %*% IB.inv)
    } else if(m == "alpha") {
        DX <- LAMBDA %*% IB.inv
    } else {
        stop("wrong model matrix names: ", m, "\n")
    }

    DX <- DX[, idx, drop=FALSE]
    DX
}



