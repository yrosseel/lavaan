# Model-methods

#
# initial version: YR 25/03/2009: `methods' for the Model class

getModelParameters <- function(object, GLIST=NULL, type="free",
                               extra=TRUE) {

    # type == "free": only non-redundant free parameters (x)
    # type == "unco": all free parameters (including constrained ones)
    # type == "user": all parameters listed in User model

    # state or final?
    if(is.null(GLIST)) GLIST <- object@GLIST

    if(type == "free") {
        N <- object@nx.free
    } else if(type == "unco") {
        N <- object@nx.unco
    } else if(type == "user") {
        N <- object@nx.user
    }
    x <- numeric(N)

    for(mm in 1:length(object@GLIST)) {
        if(type == "free") {
            m.idx <- object@m.free.idx[[mm]]
            x.idx <- object@x.free.idx[[mm]]
        } else if(type == "unco") { 
            m.idx <- object@m.unco.idx[[mm]]
            x.idx <- object@x.unco.idx[[mm]]
        } else if(type == "user") {
            m.idx <- object@m.user.idx[[mm]]
            x.idx <- object@x.user.idx[[mm]]
        }
        x[x.idx] <- GLIST[[mm]][m.idx]
    }

    if(type == "user" && extra && sum(object@x.def.idx,
                                      object@x.ceq.idx, 
                                      object@x.cin.idx) > 0L) {
        # we need 'free' x
        x.free <- getModelParameters(object, GLIST=GLIST, type="free")
        if(length(object@x.def.idx) > 0L) {
            x[object@x.def.idx] <- object@def.function(x.free)
        }
        if(length(object@x.ceq.idx) > 0L) {
            x[object@x.ceq.idx] <- object@ceq.function(x.free)
        }
        if(length(object@x.cin.idx) > 0L) {
            x[object@x.cin.idx] <- object@cin.function(x.free)
        }
    }

    x
}

# warning: this will make a copy of object
setModelParameters <- function(object, x=NULL) {

    tmp <- object@GLIST
    for(mm in 1:length(object@GLIST)) {
        m.free.idx <- object@m.free.idx[[mm]]
        x.free.idx <- object@x.free.idx[[mm]]
        tmp[[mm]][m.free.idx] <- x[x.free.idx]
    }

    # categorical? set theta elements (if any)
    if(object@categorical) {
        nmat <- object@nmat
        if(object@representation == "LISREL") {
            theta.idx <- which(names(tmp) == "theta")
            Sigma.hat <- computeSigmaHat(object, GLIST=tmp)
            for(g in 1:object@ngroups) {
                # which mm belong to group g?
                mm.in.group <- 1:nmat[g] + cumsum(c(0L,nmat))[g]
                num.idx <- object@num.idx[[g]]
                if(length(num.idx) > 0L) {
                    diag(tmp[[theta.idx[g]]])[-num.idx] <- 0.0
                } else {
                    diag(tmp[[theta.idx[g]]]) <- 0.0
                }
                MLIST <- tmp[mm.in.group]
                Sigma.hat <- computeSigmaHat.LISREL(MLIST = MLIST)
                if(length(num.idx) > 0L) {
                    diag(tmp[[theta.idx[g]]])[-num.idx] <-
                        (1 - diag(Sigma.hat)[-num.idx])
                } else {
                    diag(tmp[[theta.idx[g]]]) <- (1 - diag(Sigma.hat))
                }
            }
        } else {
            cat("FIXME: deal with theta elements in the categorical case")
        }
    }

    object@GLIST <- tmp

    object
}

# create a standalone GLIST, filled with (new) x values
# (avoiding a copy of object)
x2GLIST <- function(object, x=NULL, type="free") {

    GLIST <- object@GLIST
    for(mm in 1:length(GLIST)) {
        if(type == "free") {
            m.el.idx <- object@m.free.idx[[mm]]
            x.el.idx <- object@x.free.idx[[mm]]
        } else if(type == "full") {
            if(object@isSymmetric[mm]) {
                N <- ncol(GLIST[[mm]])
                m.el.idx <- vech.idx(N)
            } else {
                m.el.idx <- 1:length(GLIST[[mm]])
            }
            x.el.idx <- 1:length(m.el.idx)
            if(mm > 1) x.el.idx <- x.el.idx + sum(object@mmSize[1:(mm-1)])
        }

        # assign
        GLIST[[mm]][m.el.idx] <- x[x.el.idx]

        # make symmetric (if full)
        if(type == "full" && object@isSymmetric[mm]) {
            T <- t(GLIST[[mm]])
            GLIST[[mm]][upper.tri(GLIST[[mm]])] <- T[upper.tri(T)]
        }
    }

    GLIST
}

computeSigmaHat <- function(object, GLIST=NULL, extra=FALSE) {

    # state or final?
    if(is.null(GLIST)) GLIST <- object@GLIST

    nmat           <- object@nmat
    nvar           <- object@nvar
    ngroups        <- object@ngroups
    representation <- object@representation
    categorical    <- object@categorical

    # return a list
    Sigma.hat <- vector("list", length=ngroups)

    for(g in 1:ngroups) {

        # which mm belong to group g?
        mm.in.group <- 1:nmat[g] + cumsum(c(0,nmat))[g]
        MLIST <- GLIST[mm.in.group]

        if(representation == "LISREL") {
            Sigma.hat[[g]] <- computeSigmaHat.LISREL(MLIST = MLIST)
        } else {
            stop("only representation LISREL has been implemented for now")
        }

        if(extra) {
            # check if matrix is positive definite
            ev <- eigen(Sigma.hat[[g]], symmetric=TRUE, only.values=TRUE)$values
            if(any(ev < 0)) {
                attr(Sigma.hat[[g]], "po") <- FALSE
            } else {
                ## FIXME
                ## since we already do an 'eigen' decomposition, we should
                ## 'reuse' that information, instead of doing a new cholesky
                Sigma.hat.inv <-  inv.chol(Sigma.hat[[g]], logdet=TRUE)
                Sigma.hat.log.det <- attr(Sigma.hat.inv, "logdet")
                attr(Sigma.hat[[g]], "po") <- TRUE
                attr(Sigma.hat[[g]], "inv") <- Sigma.hat.inv
                attr(Sigma.hat[[g]], "log.det") <- Sigma.hat.log.det
            }
        }
    } # ngroups

    Sigma.hat
}

computeMuHat <- function(object, GLIST=NULL) {

    # state or final?
    if(is.null(GLIST)) GLIST <- object@GLIST

    nmat           <- object@nmat
    ngroups        <- object@ngroups
    representation <- object@representation
    meanstructure  <- object@meanstructure

    # return a list
    Mu.hat <- vector("list", length=ngroups)

    for(g in 1:ngroups) {

        # which mm belong to group g?
        mm.in.group <- 1:nmat[g] + cumsum(c(0,nmat))[g]

        if(!meanstructure) {
            Mu.hat[[g]] <- numeric( object@nvar[g] )
        } else
        if(representation == "LISREL") {
            Mu.hat[[g]] <- computeMuHat.LISREL(MLIST = GLIST[ mm.in.group ])
        } else {
            stop("only representation LISREL has been implemented for now")
        }
    } # ngroups

    Mu.hat
}

# TH.star = DELTA.star * (th.star - pi0.star)
# see Muthen 1984 eq 11
computeTH <- function(object, GLIST=NULL) {

    # state or final?
    if(is.null(GLIST)) GLIST <- object@GLIST

    ngroups <- object@ngroups
    nmat    <- object@nmat
    representation <- object@representation
    th.idx  <- object@th.idx

    # return a list
    TH <- vector("list", length=ngroups)

    # compute TH for each group
    for(g in 1:ngroups) {
        # which mm belong to group g?
        mm.in.group <- 1:nmat[g] + cumsum(c(0,nmat))[g]

        if(representation == "LISREL") {
            TH[[g]] <- computeTH.LISREL(MLIST = GLIST[ mm.in.group ],
                                        th.idx=th.idx[[g]])
        } else {
            stop("only representation LISREL has been implemented for now")
        }
    }

    TH
}

# PI = slope structure
# see Muthen 1984 eq 12
computePI <- function(object, GLIST=NULL) {

    # state or final?
    if(is.null(GLIST)) GLIST <- object@GLIST

    ngroups        <- object@ngroups
    nmat           <- object@nmat
    representation <- object@representation
    fixed.x        <- object@fixed.x

    # return a list
    PI <- vector("list", length=ngroups)

    # compute TH for each group
    for(g in 1:ngroups) {
        # which mm belong to group g?
        mm.in.group <- 1:nmat[g] + cumsum(c(0,nmat))[g]
        MLIST <- GLIST[ mm.in.group ]

        if(!fixed.x) {
            PI.g <- numeric( object@nvar[g] )
        } else
        if(representation == "LISREL") {
            PI.g <- computePI.LISREL(MLIST = MLIST)
        } else {
            stop("only representation LISREL has been implemented for now")
        }

        PI[[g]] <- PI.g
    }

    PI
}

# unconditional variances of Y
#  - same as diag(Sigma.hat) if all Y are continuous)
#  - 1.0 (or delta^2) if categorical
#  - if also Gamma, cov.x is used (only if categorical)
computeVY <- function(object, GLIST=NULL, samplestats=NULL) {
    # state or final?
    if(is.null(GLIST)) GLIST <- object@GLIST

    ngroups        <- object@ngroups
    nmat           <- object@nmat
    representation <- object@representation

    # return a list
    VY <- vector("list", length=ngroups)

    # compute TH for each group
    for(g in 1:ngroups) {
        # which mm belong to group g?
        mm.in.group <- 1:nmat[g] + cumsum(c(0,nmat))[g]
        MLIST <- GLIST[ mm.in.group ]

        cov.x <- samplestats@cov.x[[g]]
        num.idx <- object@num.idx[[g]]

        if(representation == "LISREL") {
            VY.g <- computeVY.LISREL(MLIST = MLIST, cov.x = cov.x,
                                     num.idx = num.idx)
        } else {
            stop("only representation LISREL has been implemented for now")
        }

        VY[[g]] <- VY.g
    }

    VY
}

# ETA: latent variances variances/covariances
computeETA <- function(object, GLIST=NULL, samplestats=NULL) {
    # state or final?
    if(is.null(GLIST)) GLIST <- object@GLIST

    ngroups        <- object@ngroups
    nmat           <- object@nmat
    representation <- object@representation

    # return a list
    ETA <- vector("list", length=ngroups)

    # compute TH for each group
    for(g in 1:ngroups) {
        # which mm belong to group g?
        mm.in.group <- 1:nmat[g] + cumsum(c(0,nmat))[g]
        MLIST <- GLIST[ mm.in.group ]

        cov.x <- samplestats@cov.x[[g]]
        num.idx <- object@num.idx[[g]]

        if(representation == "LISREL") {
            ETA.g <- computeETA.LISREL(MLIST = MLIST, cov.x = cov.x,
                                       num.idx = num.idx)
        } else {
            stop("only representation LISREL has been implemented for now")
        }

        ETA[[g]] <- ETA.g
    }

    ETA
}

computeObjective <- function(object, GLIST=NULL, 
                             samplestats=NULL, X = NULL,
                             estimator="ML", verbose=FALSE, forcePD=TRUE) {

    # state or final?
    if(is.null(GLIST)) GLIST <- object@GLIST

    # shortcut for data.type == "none"
    if(length(samplestats@cov) == 0L) {
        fx <- as.numeric(NA)
        attr(fx, "fx.group") <- rep(as.numeric(NA), samplestats@ngroups)
        return(fx)
    }

    meanstructure <- object@meanstructure
    categorical   <- object@categorical
    fixed.x       <- object@fixed.x
    num.idx       <- object@num.idx
    th.idx        <- object@th.idx

    # compute moments for all groups
    Sigma.hat <- computeSigmaHat(object, GLIST=GLIST, extra=(estimator=="ML"))
    if(meanstructure && !categorical) {
        Mu.hat <- computeMuHat(object, GLIST=GLIST)
    } else if(categorical) {
        TH <- computeTH(object, GLIST=GLIST)
        if(fixed.x) 
            PI <- computePI(object, GLIST=GLIST)
    }
 
    fx <- 0.0
    fx.group <- numeric( samplestats@ngroups )
    for(g in 1:samplestats@ngroups) {

        # incomplete data and fiml?
        if(samplestats@missing.flag) {
            if(estimator == "ML") {
                # FIML

                # catch non-pd Sigma.hat
                # 0.4-5: instead of using force.pd, we return Inf
                if(!attr(Sigma.hat[[g]], "po")) return(Inf)

                group.fx <- estimator.FIML(Sigma.hat=Sigma.hat[[g]],
                                           Mu.hat=Mu.hat[[g]],
                                           M=samplestats@missing[[g]],
                                           h1=samplestats@missing.h1[[g]]$h1)
            } else {
                stop("this estimator: `", estimator, 
                     "' can not be used with incomplete data and the missing=\"ml\" option")
            }
        } else if(estimator == "ML") {
        # complete data
            # ML and friends
            group.fx <- estimator.ML(Sigma.hat=Sigma.hat[[g]], 
                                     Mu.hat=Mu.hat[[g]],
                                     data.cov=samplestats@cov[[g]], 
                                     data.mean=samplestats@mean[[g]], 
                                     data.cov.log.det=samplestats@cov.log.det[[g]],
                                     meanstructure=meanstructure)
        } else if(estimator == "GLS"  || 
                  estimator == "WLS"  || 
                  estimator == "DWLS" ||
                  estimator == "ULS") {
            if(categorical) {
                # order of elements is important here:
                # 1. thresholds + means (interleaved)
                # 2. slopes (if any, columnwise per exo)
                # 3. variances (if any)
                # 4. correlations (no diagonal!)
                if(fixed.x) {
                    WLS.est <- c(TH[[g]],vec(PI[[g]]), 
                                 diag(Sigma.hat[[g]])[num.idx[[g]]],
                                 vech(Sigma.hat[[g]], diagonal=FALSE))
                } else {
                    WLS.est <- c(TH[[g]],diag(Sigma.hat[[g]])[num.idx[[g]]],
                                 vech(Sigma.hat[[g]], diagonal=FALSE))
                }
                #cat("WLS.obs = \n")
                #print(samplestats@WLS.obs[[g]])
                #cat("WLS.est = \n")
                #print(WLS.est)
                #cat("TH = \n"); print(TH[[g]])
                #cat("PI = \n"); print(PI[[g]])
                #cat("VAR = \n"); print(diag(Sigma.hat[[g]])[num.idx[[g]]])
                #cat("COR = \n"); print(vech(Sigma.hat[[g]], diag=FALSE))
            } else if(meanstructure) {
                WLS.est <- c(Mu.hat[[g]], vech(Sigma.hat[[g]]))
            } else {
                WLS.est <- vech(Sigma.hat[[g]])
            }
            group.fx <- estimator.WLS(WLS.est = WLS.est,
                                      WLS.obs = samplestats@WLS.obs[[g]], 
                                      WLS.V=samplestats@WLS.V[[g]])  
        } else if(estimator == "PML") {
            # Pairwise maximum likelihood
            group.fx <- estimator.PML(Sigma.hat = Sigma.hat[[g]],
                                      TH        = TH[[g]],
                                      th.idx    = th.idx[[g]],
                                      num.idx   = num.idx[[g]],
                                      X         = X[[g]])
        } else {
            stop("unsupported estimator: ", estimator)
        }

        if(estimator == "ML") {
            group.fx <- 0.5 * group.fx
        } else if(estimator == "PML") {
            # do nothing
        } else {
            group.fx <- 0.5 * (samplestats@nobs[[g]]-1)/samplestats@nobs[[g]] * group.fx
        }

        fx.group[g] <- group.fx
    } # g

    if(samplestats@ngroups > 1) {
        nobs <- unlist(samplestats@nobs)
        fx <- weighted.mean(fx.group, w=nobs)
    } else { # single group
        fx <- fx.group[1]
    }

    attr(fx, "fx.group") <- fx.group

    fx
}

# for testing purposes only
#computeDeltaNumerical <- function(object, GLIST=NULL, g=1) {
#
#    # state or final?
#   if(is.null(GLIST)) GLIST <- object@GLIST
#   
#   compute.moments <- function(x) {
#       GLIST <- x2GLIST(object, x=x, type="free")
#       Sigma.hat <- computeSigmaHat(object, GLIST=GLIST)
#        S.vec <- vech(Sigma.hat[[g]])
#        if(object@meanstructure) {
#            Mu.hat <- computeMuHat(object, GLIST=GLIST)
#            out <- c(Mu.hat[[g]], S.vec)
#        } else {   
#            out <- S.vec
#        }
#        out
#    }
#
#    x <- getModelParameters(object, GLIST=GLIST, type="free")
#    Delta <- numDeriv::jacobian(func=compute.moments, x=x)
#
#    Delta
#}


### FIXME: should we here also:
###        - weight for groups? (no, for now)
###        - handle equality constraints? (yes, for now)
computeDelta <- function(object, GLIST.=NULL, m.el.idx.=NULL, x.el.idx.=NULL) {

    representation <- object@representation
    categorical    <- object@categorical
    nmat           <- object@nmat
    ngroups        <- object@ngroups
    nvar           <- object@nvar
    num.idx        <- object@num.idx
    th.idx         <- object@th.idx
    nexo           <- object@nexo

    # number of thresholds per group (if any)
    nth <- sapply(th.idx, function(x) sum(x > 0L))

    # state or final?
    if(is.null(GLIST.)) 
        GLIST <- object@GLIST
    else
        GLIST <- GLIST.

    # type = "free" or something else?
    type <- "nonfree"
    m.el.idx <- m.el.idx.; x.el.idx <- x.el.idx.
    if(is.null(m.el.idx) && is.null(x.el.idx)) 
         type <- "free"

    # number of rows in DELTA.group
    pstar <- integer(ngroups)
    for(g in 1:ngroups) {
        pstar[g] <- as.integer(nvar[g] * (nvar[g] + 1) / 2)
        if(object@meanstructure) {
            pstar[g] <- nvar[g] + pstar[g]  # first the means, then sigma
        }
        if(categorical) {
            pstar[g] <- pstar[g] - nvar[g] # remove variances
            pstar[g] <- pstar[g] - nvar[g] # remove means

            pstar[g] <- pstar[g] + nth[g]  # add thresholds
            pstar[g] <- pstar[g] + length(num.idx[[g]]) # add num means
            pstar[g] <- pstar[g] + length(num.idx[[g]]) # add num vars
            if(nexo[g] > 0L)
                pstar[g] <- pstar[g] + (nvar[g] * nexo[g]) # add slopes
        }
    }

    # number of columns in DELTA + m.el.idx/x.el.idx
    if(type == "free") {
        NCOL <- object@nx.unco
        m.el.idx <- x.el.idx <- vector("list", length=length(GLIST))
        for(mm in 1:length(GLIST)) {
            m.el.idx[[mm]] <- object@m.unco.idx[[mm]]
            x.el.idx[[mm]] <- object@x.unco.idx[[mm]]
            # handle symmetric matrices
            if(object@isSymmetric[mm]) {
                # since we use 'x.unco.idx', only symmetric elements
                # are duplicated (not the equal ones, only in x.free.unco)
                dix <- duplicated(x.el.idx[[mm]])
                if(any(dix)) {
                    m.el.idx[[mm]] <- m.el.idx[[mm]][!dix]
                    x.el.idx[[mm]] <- x.el.idx[[mm]][!dix]
                }
            }    
        }
    } else {
        ## FIXME: this does *not* take into account symmetric
        ##        matrices; hence NCOL will be too large, and empty
        ##        columns will be added 
        ##        this is ugly, but it doesn't hurt
        ## alternative could be:
        ## NCOL <- sum(unlist(lapply(x.el.idx, function(x) length(unique(x)))))
        #NCOL <- sum(unlist(lapply(m.el.idx, length)))
        NCOL <- sum(unlist(lapply(x.el.idx, function(x) length(unique(x)))))
        # sanity check
        #nx <- sum(unlist(lapply(x.el.idx, length)))
        #stopifnot(NCOL == nx)
    }


    # compute Delta
    Delta <- vector("list", length=ngroups)
    for(g in 1:ngroups) {
        Delta.group <- matrix(0, nrow=pstar[g], ncol=NCOL)

        # which mm belong to group g?
        mm.in.group <- 1:nmat[g] + cumsum(c(0,nmat))[g]

        for(mm in mm.in.group) {
            mname <- names(object@GLIST)[mm]

            # skip empty ones
            if(!length(m.el.idx[[mm]])) next

            # get Delta columns for this model matrix
            if(representation == "LISREL") {
                DELTA <- derivative.sigma.LISREL(m=mname,
                                                 idx=m.el.idx[[mm]],
                                                 MLIST=GLIST[ mm.in.group ])
                if(categorical) {
                    # reorder: first variances (of numeric), then covariances
                    cov.idx  <- lavaan:::vech.idx(nvar[g])
                    covd.idx <- lavaan:::vech.idx(nvar[g], diag=FALSE)

                    var.idx <- which(is.na(match(cov.idx, 
                                                 covd.idx)))[num.idx[[g]]]
                    cor.idx <- match(covd.idx, cov.idx)
 
                    DELTA <- rbind(DELTA[var.idx,,drop=FALSE], 
                                   DELTA[cor.idx,,drop=FALSE])
                }
                if(object@meanstructure && !categorical) {
                    DELTA.mu <- derivative.mu.LISREL(m=mname,
                                                     idx=m.el.idx[[mm]],
                                                     MLIST=GLIST[ mm.in.group ])
                    DELTA <- rbind(DELTA.mu, DELTA)
                } else if(categorical) {
                    DELTA.th <- derivative.th.LISREL(m=mname,
                                                     idx=m.el.idx[[mm]],
                                                     th.idx=th.idx[[g]],
                                                     MLIST=GLIST[ mm.in.group ])
                    if(object@nexo[g] > 0L) {
                        DELTA.pi <- derivative.pi.LISREL(m=mname,
                                                         idx=m.el.idx[[mm]],
                                                         MLIST=GLIST[ mm.in.group ])
                        DELTA <- rbind(DELTA.th, DELTA.pi, DELTA)
                    } else {
                        DELTA <- rbind(DELTA.th, DELTA)
                    }
                }
            } else {
                stop("representation", representation, "not implemented yet")
            }

            Delta.group[ ,x.el.idx[[mm]]] <- DELTA
        } # mm

        # if type == "free" take care of equality constraints
        if(type == "free" && object@eq.constraints) {
            Delta.group <- Delta.group %*% object@eq.constraints.K
        }

        Delta[[g]] <- Delta.group

    } # g

    Delta
}

computeOmega <- function(Sigma.hat=NULL, Mu.hat=NULL,  
                         samplestats=NULL, estimator="ML", meanstructure=FALSE) {

    Omega    <- vector("list", length=samplestats@ngroups)
    Omega.mu <- vector("list", length=samplestats@ngroups)

    for(g in 1:samplestats@ngroups) {

        # ML
        if(estimator == "ML") {

            if(attr(Sigma.hat[[g]], "po") == FALSE) {
                # FIXME: WHAT IS THE BEST THING TO DO HERE??
                # CURRENTLY: stop
                stop("computeGradient: Sigma.hat is not positive definite\n")
                #Sigma.hat[[g]] <- force.pd(Sigma.hat[[g]])
                #Sigma.hat.inv <- inv.chol(Sigma.hat[[g]], logdet=TRUE)
                #Sigma.hat.log.det <- attr(Sigma.hat.inv, "logdet")
            } else {
                Sigma.hat.inv <-  attr(Sigma.hat[[g]], "inv")
                Sigma.hat.log.det <- attr(Sigma.hat[[g]], "log.det")
            }

            if(!samplestats@missing.flag) { # complete data
                if(meanstructure) {
                    diff <- samplestats@mean[[g]] - Mu.hat[[g]]
                    W.tilde <- samplestats@cov[[g]] + tcrossprod(diff)
                    # Browne 1995 eq 4.55
                    Omega.mu[[g]] <- t(t(diff) %*% Sigma.hat.inv)
                    Omega[[g]] <- 
                        ( Sigma.hat.inv %*% (W.tilde - Sigma.hat[[g]]) %*%
                          Sigma.hat.inv )
                } else {
                    W.tilde <- samplestats@cov[[g]]
                    Omega[[g]] <- 
                        ( Sigma.hat.inv %*% (W.tilde - Sigma.hat[[g]]) %*%
                          Sigma.hat.inv )
                }
            } else { # missing data
                M <- samplestats@missing[[g]]

                nvar <- ncol(samplestats@cov[[g]])
                OMEGA    <- matrix(0, nvar, nvar)
                OMEGA.MU <- matrix(0, nvar, 1)

                for(p in 1:length(M)) {
                    SX <- M[[p]][["SX"]]
                    MX <- M[[p]][["MX"]]
                    nobs <- M[[p]][["nobs"]]
                    var.idx <- M[[p]][["var.idx"]]

                    Sigma.inv <- inv.chol(Sigma.hat[[g]][var.idx, var.idx],
                                          logdet=FALSE)
                    Mu <- Mu.hat[[g]][var.idx]
                    W.tilde <- SX + tcrossprod(MX - Mu)

                    OMEGA.MU[var.idx, 1] <-
                        ( OMEGA.MU[var.idx, 1] + nobs/samplestats@ntotal *
                          t(t(MX - Mu) %*% Sigma.inv) )

                    OMEGA[var.idx, var.idx] <-
                        ( OMEGA[var.idx, var.idx] + nobs/samplestats@ntotal *
                          (Sigma.inv %*%
                           (W.tilde - Sigma.hat[[g]][var.idx,var.idx]) %*%
                           Sigma.inv ) )
                }
                Omega.mu[[g]] <- OMEGA.MU
                Omega[[g]]    <- OMEGA
            } # missing

        # GLS
        } else if(estimator == "GLS") {
            W.inv <- samplestats@icov[[g]]
            W     <- samplestats@cov[[g]]
            Omega[[g]] <- (samplestats@nobs[[g]]-1)/samplestats@nobs[[g]] *
                              (W.inv %*% (W - Sigma.hat[[g]]) %*% W.inv)
            if(meanstructure) {
                diff <- as.matrix(samplestats@mean[[g]] - Mu.hat[[g]])
                Omega.mu[[g]] <- t( t(diff) %*% W.inv )
            }
        }

    } # g

    if(meanstructure) attr(Omega, "mu") <- Omega.mu

    Omega
}


computeGradient <- function(object, GLIST=NULL, samplestats=NULL, 
                            X=NULL, type="free", 
                            estimator="ML", verbose=FALSE, forcePD=TRUE, 
                            group.weight=TRUE, constraints=TRUE,
                            Delta=NULL) {

    nmat           <- object@nmat
    representation <- object@representation
    meanstructure  <- object@meanstructure
    categorical    <- object@categorical
    fixed.x        <- object@fixed.x
    num.idx        <- object@num.idx
    nx.unco        <- object@nx.unco

    # state or final?
    if(is.null(GLIST)) GLIST <- object@GLIST

    # group.w
    if(group.weight) {
        group.w <- (unlist(samplestats@nobs)/samplestats@ntotal)
    } else {
        group.w <- rep(1.0, samplestats@ngroups)
    }

    # Sigma.hat + Mu.hat
    Sigma.hat <- computeSigmaHat(object, GLIST=GLIST, extra=(estimator == "ML"))
    if(meanstructure && !categorical) {
        Mu.hat <- computeMuHat(object, GLIST=GLIST)
    } else if(categorical) {
        TH <- computeTH(object, GLIST=GLIST)
        if(fixed.x)
            PI <- computePI(object, GLIST=GLIST)
    }
    

    # two approaches:
    # - ML/GLS approach: using Omega (and Omega.mu)
    # - WLS: using Delta

    # 1. ML/GLS approach
    if(estimator == "ML" || estimator == "GLS") {
        if(meanstructure) {
            Omega <- computeOmega(Sigma.hat=Sigma.hat, Mu.hat=Mu.hat,
                                  samplestats=samplestats, estimator=estimator, 
                                  meanstructure=TRUE)
            Omega.mu <- attr(Omega, "mu")
        } else {
            Omega <- computeOmega(Sigma.hat=Sigma.hat, Mu.hat=NULL,
                                  samplestats=samplestats, estimator=estimator,
                                  meanstructure=FALSE)
            Omega.mu <- vector("list", length=samplestats@ngroups)
        }

        # compute DX (for all elements in every model matrix)
        DX <- vector("list", length=length(GLIST))
      
        for(g in 1:samplestats@ngroups) {
            # which mm belong to group g?
            mm.in.group <- 1:nmat[g] + cumsum(c(0,nmat))[g]
            mm.names <- names( GLIST[mm.in.group] )

            if(representation == "LISREL") {
                DX.group <- derivative.F.LISREL(GLIST[mm.in.group], 
                                                Omega[[g]],
                                                Omega.mu[[g]])

                # only save what we need
                DX[mm.in.group] <- DX.group[ mm.names ] 
            } else {
                stop("only representation LISREL has been implemented for now")
            }

            # weight by group
            if(samplestats@ngroups > 1L) {
                for(mm in mm.in.group) {
                    DX[[mm]] <- group.w[g] * DX[[mm]]
                }
            }
        }

        # extract free parameters + weight by group
        if(type == "free") {
            dx <- numeric( nx.unco )
            for(g in 1:samplestats@ngroups) {
                mm.in.group <- 1:nmat[g] + cumsum(c(0,nmat))[g]
                for(mm in mm.in.group) {
                      m.unco.idx  <- object@m.unco.idx[[mm]]
                         x.unco.idx  <- object@x.unco.idx[[mm]]
                      dx[x.unco.idx] <- DX[[mm]][m.unco.idx]
                }
            }

            # handle equality constraints
            if(object@eq.constraints && constraints) {
                dx <- as.numeric( t(object@eq.constraints.K) %*% dx )
            }
        } else {
            dx <- DX
            # handle equality constraints
            ### FIXME!!!! TODO!!!!
        } 

    } else # ML/GLS

    # 2. WLS approach
    if(estimator == "WLS" || estimator == "DWLS" || estimator == "ULS") {

        if(type != "free") {
            if(is.null(Delta))
                stop("FIXME: Delta should be given if type != free")
            #stop("FIXME: WLS gradient with type != free needs fixing!")
        } else {
            Delta <- computeDelta(object, GLIST.=GLIST)
        }

        for(g in 1:samplestats@ngroups) {
            if(categorical) {                                    
                # order of elements is important here:
                # 1. thresholds + means (interleaved)
                # 2. slopes (if any)
                # 3. variances (if any)
                # 4. correlations (no diagonal!)
                if(fixed.x) {
                    WLS.est <- c(TH[[g]],PI[[g]],
                                 diag(Sigma.hat[[g]])[num.idx[[g]]],
                                 vech(Sigma.hat[[g]], diagonal=FALSE))
                } else {
                    WLS.est <- c(TH[[g]],diag(Sigma.hat[[g]])[num.idx[[g]]],
                                 vech(Sigma.hat[[g]], diagonal=FALSE))
                }
            } else if(meanstructure) {                             
                WLS.est <- c(Mu.hat[[g]], vech(Sigma.hat[[g]]))    
            } else {                                           
                WLS.est <- vech(Sigma.hat[[g]])                
            }
            WLS.obs <- samplestats@WLS.obs[[g]]
            diff <- as.matrix(WLS.obs - WLS.est)
            group.dx <- -1 * ( t(Delta[[g]]) %*% samplestats@WLS.V[[g]] %*% diff)
            group.dx <- group.w[g] * group.dx

            if(g == 1) {
                dx <- group.dx
            } else {
                dx <- dx + group.dx
            }
        } # g

        if(type == "free") {
            # nothing to do
        } else {
            # make a GLIST
            dx <- x2GLIST(object, x=dx, type="full")
        }

    } # WLS

    else if(estimator == "PML") {

        if(type != "free")
            stop("FIXME: type != free in computeGradient for estimator PML")

        for(g in 1:samplestats@ngroups) {
            group.dx <- numeric( nx.unco )

            if(g == 1) {
                dx <- group.dx
            } else {
                dx <- dx + group.dx
            }
        } # g
    } 

    else {
        stop("lavaan ERROR: no analytical gradient available for estimator ",
             estimator)
    }

    dx
}


estimateModel <- function(object, samplestats=NULL, X=NULL, do.fit=TRUE, 
                          options=NULL, control=list()) {

    estimator     <- options$estimator
    verbose       <- options$verbose
    debug         <- options$debug

    if(samplestats@missing.flag) { 
        group.weight <- FALSE
    } else {
        group.weight <- TRUE
    }

    # function to be minimized
    minimize.this.function <- function(x, verbose=FALSE, infToMax=FALSE) {
      
        #cat("DEBUG: x = ", x, "\n")

        # current strategy: forcePD is by default FALSE, except
        # if missing patterns are used
        #if(any(samplestats@missing.flag)) {
        #    forcePD <- TRUE
        #} else {
            forcePD <- FALSE
        #}

        # transform variances back
        #x[object@x.free.var.idx] <- tan(x[object@x.free.var.idx])

        # update GLIST (change `state') and make a COPY!
        GLIST <- x2GLIST(object, x=x)

        fx <- computeObjective(object, GLIST=GLIST, 
                               samplestats=samplestats, X=X,
                               estimator=estimator, verbose=verbose,
                               forcePD=forcePD)	
        if(debug || verbose) { 
            cat("Objective function  = ", sprintf("%12.10f", fx), "\n", sep="") 
        }
        if(debug) {
            cat("Current unconstrained parameter values =\n")
            tmp.x <- getModelParameters(object, GLIST=GLIST, type="unco")
            print(tmp.x); cat("\n")
            cat("Current free parameter values =\n"); print(x); cat("\n")
        }

        # for L-BFGS-B
        if(infToMax && is.infinite(fx)) fx <- 1e20

        fx
    }

    first.derivative.param <- function(x, verbose=FALSE, infToMax=FALSE) {

        # transform variances back
        #x[object@x.free.var.idx] <- tan(x[object@x.free.var.idx])

        # update GLIST (change `state') and make a COPY!
        GLIST <- x2GLIST(object, x=x)

        dx <- computeGradient(object, GLIST=GLIST, samplestats,
                              type="free", 
                              group.weight=group.weight, ### check me!!
                              estimator=estimator,
                              verbose=verbose, forcePD=TRUE)

        if(debug) {
            cat("Gradient function (analytical) =\n"); print(dx); cat("\n")
        }

        dx
    } 

    first.derivative.param.numerical <- function(x, verbose=FALSE) {

        # transform variances back
        #x[object@x.free.var.idx] <- tan(x[object@x.free.var.idx])

        # numerical approximation using the Richardson method
        npar <- length(x)
        h <- 10e-6
        dx <- numeric( npar )
 
        ## FIXME: call computeObjective directly!!
        for(i in 1:npar) {
            x.left <- x.left2 <- x.right <- x.right2 <- x
            x.left[i]  <- x[i] - h; x.left2[i]  <- x[i] - 2*h
            x.right[i] <- x[i] + h; x.right2[i] <- x[i] + 2*h
            fx.left   <- minimize.this.function(x.left)
            fx.left2  <- minimize.this.function(x.left2)
            fx.right  <- minimize.this.function(x.right)
            fx.right2 <- minimize.this.function(x.right2)
            dx[i] <- (fx.left2 - 8*fx.left + 8*fx.right - fx.right2)/(12*h)
        }

        #dx <- lavGradientC(func=minimize.this.function, x=x)
        # does not work if pnorm is involved... (eg PML)

        if(debug) {
            cat("Gradient function (numerical) =\n"); print(dx); cat("\n")
        }        

        dx
    }
 
    # starting values
    start.x <- getModelParameters(object)
    if(debug) {
        cat("start.unco = ", getModelParameters(object, type="unco"), "\n")
        cat("start.x = ", start.x, "\n")
    }

    # scaling factors
    # FIXME: what is the best way to set the scale??
    # current strategy: if startx > 1.0, we rescale by using
    # 1/startx
    SCALE <- rep(1.0, length(start.x))
    #idx <- which(abs(start.x) > 10.0)
    idx <- which(abs(start.x) > 1.0)
    if(length(idx) > 0L) SCALE[idx] <- abs(1.0/start.x[idx])
    #idx <- which(abs(start.x) < 1.0 & start.x != 0.0)
    #if(length(idx) > 0L) SCALE[idx] <- abs(1.0/start.x[idx])
    if(debug) {
        cat("SCALE = ", SCALE, "\n")
    }

    # transforming variances using atan (or another sigmoid function?)
    # FIXME: better approach?
    #start.x[object@x.free.var.idx] <- atan(start.x[object@x.free.var.idx])


    # first some nelder mead steps? (default = FALSE)
    if(is.null(control$init_nelder_mead)) {
        INIT_NELDER_MEAD <- FALSE
    } else {
        INIT_NELDER_MEAD <- control$init_nelder_mead
    }

    # optimizer
    if(is.null(body(object@ceq.function)) && 
       is.null(body(object@cin.function)) ) {
        if(is.null(control$optim.method)) {
            OPTIMIZER <- "NLMINB"
            #OPTIMIZER <- "BFGS"  # slightly slower, no bounds; better scaling!
            #OPTIMIZER <- "L-BFGS-B"  # trouble with Inf values for fx!
        } else {
            OPTIMIZER <- toupper(control$optim.method)
            stopifnot(OPTIMIZER %in% c("NLMINB", "BFGS", "L-BFGS-B"))
        }
    } else {
        if(is.null(control$optim.method)) {
            OPTIMIZER <- "NLMINB.CONSTR"
            #OPTIMIZER <- "ALABAMA"
        } else {
            OPTIMIZER <- toupper(control$optim.method)
            stopifnot(OPTIMIZER %in% c("NLMINB.CONSTR", "ALABAMA"))
        }
    }

    if(INIT_NELDER_MEAD) {
        if(verbose) cat("Initial Nelder-Mead step:\n")
        trace <- 0L; if(verbose) trace <- 1L
        optim.out <- optim(par=start.x,
                           fn=minimize.this.function,
                           method="Nelder-Mead",
                           #control=list(maxit=10L, 
                           #             parscale=SCALE,
                           #             trace=trace),
                           hessian=FALSE,
                           verbose=verbose)
        cat("\n")
        start.x <- optim.out$par
    }

    if(OPTIMIZER == "NLMINB") {
        if(verbose) cat("Quasi-Newton steps using NLMINB:\n")
        #if(debug) control$trace <- 1L;
        control.nlminb <- list(eval.max=20000L,
                               iter.max=10000L,
                               trace=0L,
                               #abs.tol=1e-20, ### important!! fx never negative
                               abs.tol=(.Machine$double.eps * 10),
                               rel.tol=1e-10,
                               x.tol=1.5e-8,
                               step.min=2.2e-14)
        control.nlminb <- modifyList(control.nlminb, control)
        control <- control.nlminb[c("eval.max", "iter.max", "trace",
                                    "abs.tol", "rel.tol", "x.tol",
                                    "step.min")]
        #cat("DEBUG: control = ", unlist(control.nlminb), "\n")
        optim.out <- nlminb(start=start.x,
                            objective=minimize.this.function,
                            #gradient=first.derivative.param,
                            gradient=first.derivative.param.numerical,
                            control=control,
                            scale=SCALE,
                            verbose=verbose) 
        if(verbose) {
            cat("convergence status (0=ok): ", optim.out$convergence, "\n")
            cat("nlminb message says: ", optim.out$message, "\n")
            cat("number of iterations: ", optim.out$iterations, "\n")
            cat("number of function evaluations [objective, gradient]: ",
                optim.out$evaluations, "\n")
        }

        iterations <- optim.out$iterations
        x          <- optim.out$par
        if(optim.out$convergence == 0) {
            converged <- TRUE
        } else {
            converged <- FALSE
        }
    } else if(OPTIMIZER == "BFGS") {

        # warning: Bollen example with estimator=GLS does NOT converge!
        # (but WLS works!)
        # - BB.ML works too

        control.bfgs <- list(trace=0L, fnscale=1, 
                             parscale=SCALE, ## or not?
                             ndeps=1e-3,
                             maxit=10000,
                             abstol=1e-20,
                             reltol=1e-10,
                             REPORT=1L)
        control.bfgs <- modifyList(control.bfgs, control)
        control <- control.bfgs[c("trace", "fnscale", "parscale", "ndeps",
                                  "maxit", "abstol", "reltol", "REPORT")]
        #trace <- 0L; if(verbose) trace <- 1L
        optim.out <- optim(par=start.x,
                           fn=minimize.this.function,
                           gr=first.derivative.param,
                           method="BFGS",
                           control=control,
                           hessian=FALSE,
                           verbose=verbose)
        if(verbose) {
            cat("convergence status (0=ok): ", optim.out$convergence, "\n")
            cat("optim BFGS message says: ", optim.out$message, "\n")
            #cat("number of iterations: ", optim.out$iterations, "\n")
            cat("number of function evaluations [objective, gradient]: ",
                optim.out$counts, "\n")
        }

        #iterations <- optim.out$iterations
        iterations <- optim.out$counts[1]
        x          <- optim.out$par
        if(optim.out$convergence == 0L) {
            converged <- TRUE
        } else {
            converged <- FALSE
        }
    } else if(OPTIMIZER == "L-BFGS-B") {

        # warning, does not cope with Inf values!!

        control.lbfgsb <- list(trace=0L, fnscale=1,
                               parscale=SCALE, ## or not?
                               ndeps=1e-3,
                               maxit=10000,
                               REPORT=1L,
                               lmm=5L,
                               factr=1e7,
                               pgtol=0)
        control.lbfgsb <- modifyList(control.lbfgsb, control)
        control <- control.lbfgsb[c("trace", "fnscale", "parscale", 
                                    "ndeps", "maxit", "REPORT", "lmm", 
                                    "factr", "pgtol")]
        optim.out <- optim(par=start.x,
                           fn=minimize.this.function,
                           gr=first.derivative.param,
                           method="L-BFGS-B",
                           control=control,
                           hessian=FALSE,
                           verbose=verbose,
                           infToMax=TRUE)
        if(verbose) {
            cat("convergence status (0=ok): ", optim.out$convergence, "\n")
            cat("optim L-BFGS-B message says: ", optim.out$message, "\n")
            #cat("number of iterations: ", optim.out$iterations, "\n")
            cat("number of function evaluations [objective, gradient]: ",
                optim.out$counts, "\n")
        }

        #iterations <- optim.out$iterations
        iterations <- optim.out$counts[1]
        x          <- optim.out$par
        if(optim.out$convergence == 0L) {
            converged <- TRUE
        } else {
            converged <- FALSE
        }
    } else if(OPTIMIZER == "NLMINB.CONSTR") {

        control.nlminb <- list(eval.max=20000L,
                               iter.max=10000L,
                               trace=0L,
                               #abs.tol=1e-20,
                               abs.tol=(.Machine$double.eps * 10),
                               rel.tol=1e-9, # 1e-10 seems 'too strict'
                               x.tol=1.5e-8,
                               step.min=2.2e-14)
        control.nlminb <- modifyList(control.nlminb, control)
        control <- control.nlminb[c("eval.max", "iter.max", "trace",
                                    "abs.tol", "rel.tol", "x.tol",
                                    "step.min")]
        cin <- cin.jac <- ceq <- ceq.jac <- NULL
        if(!is.null(body(object@cin.function))) cin     <- object@cin.function
        if(!is.null(body(object@cin.jacobian))) cin.jac <- object@cin.jacobian
        if(!is.null(body(object@ceq.function))) ceq     <- object@ceq.function
        if(!is.null(body(object@ceq.jacobian))) ceq.jac <- object@ceq.jacobian
        trace <- FALSE; if(verbose) trace <- TRUE
        optim.out <- nlminb.constr(start = start.x,
                                   objective=minimize.this.function,
                                   gradient=first.derivative.param,
                                   #gradient=first.derivative.param.numerical,
                                   control=control,
                                   scale=SCALE,
                                   verbose=verbose,
                                   cin = cin, cin.jac = cin.jac,
                                   ceq = ceq, ceq.jac = ceq.jac,
                                   control.outer = list(verbose=verbose)
                                  )
        if(verbose) {
            cat("convergence status (0=ok): ", optim.out$convergence, "\n")
            cat("nlminb.constr message says: ", optim.out$message, "\n")
            cat("number of outer iterations: ", optim.out$outer.iterations, "\n")
            cat("number of inner iterations: ", optim.out$iterations, "\n")
            cat("number of function evaluations [objective, gradient]: ",
                optim.out$evaluations, "\n")
        }

        iterations <- optim.out$iterations
        x          <- optim.out$par
        if(optim.out$convergence == 0) {
            converged <- TRUE
        } else {
            converged <- FALSE
        }
    }

    fx <- minimize.this.function(x) # to get "fx.group" attribute

    # transform variances back
    #x[object@x.free.var.idx] <- tan(x[object@x.free.var.idx])

    attr(x, "converged")  <- converged
    attr(x, "iterations") <- iterations
    attr(x, "control")    <- control
    attr(x, "fx")         <- fx
    if(!is.null(optim.out$con.jac)) attr(x, "con.jac")    <- optim.out$con.jac
    if(!is.null(optim.out$lambda))  attr(x, "con.lambda") <- optim.out$lambda

    x

}


