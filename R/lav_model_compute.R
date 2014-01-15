## FIXME: allow for GAMMA in continuous case
computeSigmaHat <- function(lavmodel = NULL, GLIST = NULL, extra = FALSE, 
                            delta = TRUE, debug = FALSE) {

    # state or final?
    if(is.null(GLIST)) GLIST <- lavmodel@GLIST

    nmat           <- lavmodel@nmat
    nvar           <- lavmodel@nvar
    ngroups        <- lavmodel@ngroups
    representation <- lavmodel@representation
    categorical    <- lavmodel@categorical

    # return a list
    Sigma.hat <- vector("list", length=ngroups)

    for(g in 1:ngroups) {

        # which mm belong to group g?
        mm.in.group <- 1:nmat[g] + cumsum(c(0,nmat))[g]
        MLIST <- GLIST[mm.in.group]

        if(representation == "LISREL") {
            Sigma.hat[[g]] <- computeSigmaHat.LISREL(MLIST = MLIST, 
                                                     delta = delta)
        } else {
            stop("only representation LISREL has been implemented for now")
        }
        if(debug) print(Sigma.hat[[g]])

        if(extra) {
            # check if matrix is positive definite
            ev <- eigen(Sigma.hat[[g]], symmetric=TRUE, only.values=TRUE)$values
            if(any(ev < .Machine$double.eps) || sum(ev) == 0) {
                Sigma.hat.inv <-  MASS::ginv(Sigma.hat[[g]])
                Sigma.hat.log.det <- log(.Machine$double.eps)
                attr(Sigma.hat[[g]], "po") <- FALSE
                attr(Sigma.hat[[g]], "inv") <- Sigma.hat.inv
                attr(Sigma.hat[[g]], "log.det") <- Sigma.hat.log.det
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

### FIXME:: replace by computeEY (to allow for GAMMA in continuous case)?
computeMuHat <- function(lavmodel = NULL, GLIST = NULL) {

    # state or final?
    if(is.null(GLIST)) GLIST <- lavmodel@GLIST

    nmat           <- lavmodel@nmat
    ngroups        <- lavmodel@ngroups
    representation <- lavmodel@representation
    meanstructure  <- lavmodel@meanstructure

    # return a list
    Mu.hat <- vector("list", length=ngroups)

    for(g in 1:ngroups) {

        # which mm belong to group g?
        mm.in.group <- 1:nmat[g] + cumsum(c(0,nmat))[g]

        if(!meanstructure) {
            Mu.hat[[g]] <- numeric( lavmodel@nvar[g] )
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
computeTH <- function(lavmodel = NULL, GLIST = NULL) {

    # state or final?
    if(is.null(GLIST)) GLIST <- lavmodel@GLIST

    ngroups <- lavmodel@ngroups
    nmat    <- lavmodel@nmat
    representation <- lavmodel@representation
    th.idx  <- lavmodel@th.idx

    # return a list
    TH <- vector("list", length=ngroups)

    # compute TH for each group
    for(g in 1:ngroups) {

        if(length(th.idx[[g]]) == 0) {
            TH[[g]] <- numeric(0L)
            next
        }
   
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
computePI <- function(lavmodel = NULL, GLIST = NULL) {

    # state or final?
    if(is.null(GLIST)) GLIST <- lavmodel@GLIST

    ngroups        <- lavmodel@ngroups
    nmat           <- lavmodel@nmat
    representation <- lavmodel@representation
    fixed.x        <- lavmodel@fixed.x

    # return a list
    PI <- vector("list", length=ngroups)

    # compute TH for each group
    for(g in 1:ngroups) {
        # which mm belong to group g?
        mm.in.group <- 1:nmat[g] + cumsum(c(0,nmat))[g]
        MLIST <- GLIST[ mm.in.group ]

        if(!fixed.x) {
            PI.g <- numeric( lavmodel@nvar[g] )
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


# GW = group weight
computeGW <- function(lavmodel = NULL, GLIST=NULL) {

    # state or final?
    if(is.null(GLIST)) GLIST <- lavmodel@GLIST

    ngroups        <- lavmodel@ngroups
    nmat           <- lavmodel@nmat
    representation <- lavmodel@representation
    group.w.free   <- lavmodel@group.w.free

    # return a list
    GW <- vector("list", length=ngroups)

    # compute GW for each group
    for(g in 1:ngroups) {
        # which mm belong to group g?
        mm.in.group <- 1:nmat[g] + cumsum(c(0,nmat))[g]
        MLIST <- GLIST[ mm.in.group ]

        if(!group.w.free) {
            GW.g <- 0.0 # FIXME
        } else
        if(representation == "LISREL") {
            GW.g <- as.numeric(MLIST$gw[1,1])
        } else {
            stop("only representation LISREL has been implemented for now")
        }

        GW[[g]] <- GW.g
    }

    # transform to proportions
    #gw <- unlist(GW)
    #gw <- exp(gw) / sum(exp(gw))
    #for(g in 1:ngroups) {
    #    GW[[g]] <- gw[g]
    #}

    GW
}

# unconditional variances of Y
#  - same as diag(Sigma.hat) if all Y are continuous)
#  - 1.0 (or delta^2) if categorical
#  - if also Gamma, cov.x is used (only if categorical)

# only for semTools compatibility
computeVY <- function(lavmodel = NULL, GLIST = NULL, lavsamplestats = NULL,
                      samplestats = NULL) {
    # state or final?
    if(is.null(GLIST)) GLIST <- lavmodel@GLIST

    # backwards compatibility (semTools!)
    if(!is.null(samplestats)) {
        lavsamplestats <- samplestats
    }

    ngroups        <- lavmodel@ngroups
    nmat           <- lavmodel@nmat
    representation <- lavmodel@representation

    # return a list
    VY <- vector("list", length=ngroups)

    # compute TH for each group
    for(g in 1:ngroups) {
        # which mm belong to group g?
        mm.in.group <- 1:nmat[g] + cumsum(c(0,nmat))[g]
        MLIST <- GLIST[ mm.in.group ]

        cov.x <- lavsamplestats@cov.x[[g]]
        num.idx <- lavmodel@num.idx[[g]]

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

# THETA: observed (residual) variances
computeTHETA <- function(lavmodel = NULL, GLIST = NULL) {
    # state or final?
    if(is.null(GLIST)) GLIST <- lavmodel@GLIST

    ngroups        <- lavmodel@ngroups
    nmat           <- lavmodel@nmat
    representation <- lavmodel@representation

    # return a list
    THETA <- vector("list", length=ngroups)

    # compute THETA for each group
    for(g in 1:ngroups) {
        # which mm belong to group g?
        mm.in.group <- 1:nmat[g] + cumsum(c(0,nmat))[g]
        MLIST <- GLIST[ mm.in.group ]

        if(representation == "LISREL") {
            THETA.g <- MLIST$theta
            # fix theta
            lv.idx <- c(lavmodel@ov.y.dummy.lv.idx[[g]],
                        lavmodel@ov.x.dummy.lv.idx[[g]])
            ov.idx <- c(lavmodel@ov.y.dummy.ov.idx[[g]],
                        lavmodel@ov.x.dummy.ov.idx[[g]])
            if(length(ov.idx) > 0L) {
                THETA.g[ov.idx, ov.idx] <- MLIST$psi[lv.idx, lv.idx]
            }

        } else {
            stop("only representation LISREL has been implemented for now")
        }

        THETA[[g]] <- THETA.g
    }

    THETA
}

# V(ETA): latent variances variances/covariances
# - if conditional=TRUE, ignore all ov.dummy. * and GAMMA
computeVETA <- function(lavmodel = NULL, GLIST = NULL, 
                        remove.dummy.lv = FALSE, lavsamplestats = NULL) {
    # state or final?
    if(is.null(GLIST)) GLIST <- lavmodel@GLIST

    ngroups        <- lavmodel@ngroups
    nmat           <- lavmodel@nmat
    representation <- lavmodel@representation

    # return a list
    ETA <- vector("list", length=ngroups)

    # compute ETA for each group
    for(g in 1:ngroups) {
        # which mm belong to group g?
        mm.in.group <- 1:nmat[g] + cumsum(c(0,nmat))[g]
        MLIST <- GLIST[ mm.in.group ]

        cov.x <- NULL
        if(!is.null(lavsamplestats)) {
            cov.x <- lavsamplestats@cov.x[[g]]
        }
       
        if(representation == "LISREL") {
            ETA.g <- computeVETA.LISREL(MLIST = MLIST, cov.x = cov.x)

            if(remove.dummy.lv) {
                # remove all dummy latent variables
                lv.idx <- c(lavmodel@ov.y.dummy.lv.idx[[g]],
                            lavmodel@ov.x.dummy.lv.idx[[g]])
                if(!is.null(lv.idx)) {
                    ETA.g <- ETA.g[-lv.idx, -lv.idx, drop=FALSE]
                }
            }
        } else {
            stop("only representation LISREL has been implemented for now")
        }

        ETA[[g]] <- ETA.g
    }

    ETA
}

# V(ETA|x_i): latent variances variances/covariances, conditional on x_
# - this is always (I-B)^-1 PSI (I-B)^-T, after REMOVING lv dummies
computeVETAx <- function(lavmodel = NULL, GLIST = NULL) {
    # state or final?
    if(is.null(GLIST)) GLIST <- lavmodel@GLIST

    ngroups        <- lavmodel@ngroups
    nmat           <- lavmodel@nmat
    representation <- lavmodel@representation

    # return a list
    ETA <- vector("list", length=ngroups)

    # compute ETA for each group
    for(g in 1:ngroups) {
        # which mm belong to group g?
        mm.in.group <- 1:nmat[g] + cumsum(c(0,nmat))[g]
        MLIST <- GLIST[ mm.in.group ]

        if(representation == "LISREL") {
            lv.idx <- c(lavmodel@ov.y.dummy.lv.idx[[g]],
                        lavmodel@ov.x.dummy.lv.idx[[g]])
            ETA.g <- computeVETAx.LISREL(MLIST = MLIST, 
                                         lv.dummy.idx = lv.idx)
        } else {
            stop("only representation LISREL has been implemented for now")
        }

        ETA[[g]] <- ETA.g
    }

    ETA
}

# COV: observed+latent variances variances/covariances
computeCOV <- function(lavmodel = NULL, GLIST = NULL, lavsamplestats = NULL) {
    # state or final?
    if(is.null(GLIST)) GLIST <- lavmodel@GLIST

    ngroups        <- lavmodel@ngroups
    nmat           <- lavmodel@nmat
    representation <- lavmodel@representation

    # return a list
    COV <- vector("list", length=ngroups)

    # compute COV for each group
    for(g in 1:ngroups) {
        # which mm belong to group g?
        mm.in.group <- 1:nmat[g] + cumsum(c(0,nmat))[g]
        MLIST <- GLIST[ mm.in.group ]

        cov.x <- lavsamplestats@cov.x[[g]]

        if(representation == "LISREL") {
            COV.g <- computeCOV.LISREL(MLIST = MLIST, cov.x = cov.x)
        } else {
            stop("only representation LISREL has been implemented for now")
        }

        COV[[g]] <- COV.g
    }

    COV
}


# E(ETA): expectation (means) of latent variables (return vector)
computeEETA <- function(lavmodel = NULL, GLIST = NULL, lavsamplestats = NULL, 
                        remove.dummy.lv = FALSE) {
  
    # state or final?
    if(is.null(GLIST)) GLIST <- lavmodel@GLIST

    ngroups        <- lavmodel@ngroups
    nmat           <- lavmodel@nmat
    representation <- lavmodel@representation

    # return a list
    EETA <- vector("list", length=ngroups)

    # compute E(ETA) for each group
    for(g in 1:ngroups) {
        # which mm belong to group g?
        mm.in.group <- 1:nmat[g] + cumsum(c(0,nmat))[g]
        MLIST <- GLIST[ mm.in.group ]

        if(representation == "LISREL") {
            EETA.g <- computeEETA.LISREL(MLIST, 
                mean.x=lavsamplestats@mean.x[[g]],
                sample.mean=lavsamplestats@mean[[g]],
                ov.y.dummy.lv.idx=lavmodel@ov.y.dummy.lv.idx[[g]],
                ov.x.dummy.lv.idx=lavmodel@ov.x.dummy.lv.idx[[g]],
                ov.y.dummy.ov.idx=lavmodel@ov.y.dummy.ov.idx[[g]],
                ov.x.dummy.ov.idx=lavmodel@ov.x.dummy.ov.idx[[g]])
            if(remove.dummy.lv) {
                # remove dummy
                lv.dummy.idx <- c(lavmodel@ov.y.dummy.lv.idx[[g]],
                                  lavmodel@ov.x.dummy.lv.idx[[g]])
                if(length(lv.dummy.idx) > 0L) {
                    EETA.g <- EETA.g[-lv.dummy.idx]
                }
            }
        } else {
            stop("only representation LISREL has been implemented for now")
        }

        EETA[[g]] <- EETA.g
    }

    EETA
}

# E(ETA|x_i): conditional expectation (means) of latent variables
# for a given value of x_i (instead of E(x_i))
computeEETAx <- function(lavmodel = NULL, GLIST = NULL, lavsamplestats = NULL, 
                         eXo = NULL, remove.dummy.lv = FALSE) {
  
    # state or final?
    if(is.null(GLIST)) GLIST <- lavmodel@GLIST

    ngroups        <- lavmodel@ngroups
    nmat           <- lavmodel@nmat
    representation <- lavmodel@representation

    # return a list
    EETAx <- vector("list", length=ngroups)

    # compute E(ETA) for each group
    for(g in 1:ngroups) {
        # which mm belong to group g?
        mm.in.group <- 1:nmat[g] + cumsum(c(0,nmat))[g]
        MLIST <- GLIST[ mm.in.group ]

        EXO <- eXo[[g]]
        if(is.null(EXO)) {
            # create empty matrix
            EXO <- matrix(0, lavsamplestats@nobs[[g]], 0L)
        }

        if(representation == "LISREL") {
            EETAx.g <- computeEETAx.LISREL(MLIST, 
                eXo=EXO,
                sample.mean=lavsamplestats@mean[[g]],
                ov.y.dummy.lv.idx=lavmodel@ov.y.dummy.lv.idx[[g]],
                ov.x.dummy.lv.idx=lavmodel@ov.x.dummy.lv.idx[[g]],
                ov.y.dummy.ov.idx=lavmodel@ov.y.dummy.ov.idx[[g]],
                ov.x.dummy.ov.idx=lavmodel@ov.x.dummy.ov.idx[[g]])

            if(remove.dummy.lv) {
                # remove dummy
                lv.dummy.idx <- c(lavmodel@ov.y.dummy.lv.idx[[g]],
                                  lavmodel@ov.x.dummy.lv.idx[[g]])
                if(length(lv.dummy.idx) > 0L) {
                    EETAx.g <- EETAx.g[ ,-lv.dummy.idx, drop=FALSE]
                }
            }
        } else {
            stop("only representation LISREL has been implemented for now")
        }

        EETAx[[g]] <- EETAx.g
    }

    EETAx
}

# return 'regular' LAMBDA
computeLAMBDA <- function(lavmodel = NULL, GLIST = NULL, 
                          remove.dummy.lv = FALSE) {
  
    # state or final?
    if(is.null(GLIST)) GLIST <- lavmodel@GLIST

    ngroups        <- lavmodel@ngroups
    nmat           <- lavmodel@nmat
    representation <- lavmodel@representation

    # return a list
    LAMBDA <- vector("list", length=ngroups)

    # compute LAMBDA for each group
    for(g in 1:ngroups) {
        # which mm belong to group g?
        mm.in.group <- 1:nmat[g] + cumsum(c(0,nmat))[g]
        MLIST <- GLIST[ mm.in.group ]

        if(representation == "LISREL") {
            LAMBDA.g <- computeLAMBDA.LISREL(MLIST)
            if(remove.dummy.lv) {
                # remove dummy
                lv.dummy.idx <- c(lavmodel@ov.y.dummy.lv.idx[[g]],
                                  lavmodel@ov.x.dummy.lv.idx[[g]])
                if(length(lv.dummy.idx) > 0L) {
                    LAMBDA.g <- LAMBDA.g[,-lv.dummy.idx,drop=FALSE]
                }
            }
        } else {
            stop("only representation LISREL has been implemented for now")
        }

        LAMBDA[[g]] <- LAMBDA.g
    }

    LAMBDA
}

# E(Y): expectation (mean) of observed variables
# returns vector 1 x nvar
computeEY <- function(lavmodel = NULL, GLIST = NULL, lavsamplestats = NULL) {
  
    # state or final?
    if(is.null(GLIST)) GLIST <- lavmodel@GLIST

    ngroups        <- lavmodel@ngroups
    nmat           <- lavmodel@nmat
    representation <- lavmodel@representation

    # return a list
    EY <- vector("list", length=ngroups)

    # compute E(Y) for each group
    for(g in 1:ngroups) {
        # which mm belong to group g?
        mm.in.group <- 1:nmat[g] + cumsum(c(0,nmat))[g]
        MLIST <- GLIST[ mm.in.group ]

        if(representation == "LISREL") {
            EY.g <- computeEY.LISREL(MLIST = MLIST, 
                                     mean.x=lavsamplestats@mean.x[[g]],
                                     sample.mean=lavsamplestats@mean[[g]])
        } else {
            stop("only representation LISREL has been implemented for now")
        }

        EY[[g]] <- EY.g
    }

    EY
}


# E(Y | ETA, x_i): conditional expectation (means) of observed variables
# for a given value of x_i AND eta_i
computeYHAT <- function(lavmodel = NULL, GLIST = NULL, lavsamplestats = NULL, 
                        eXo = NULL, ETA = NULL, duplicate = FALSE) {
  
    # state or final?
    if(is.null(GLIST)) GLIST <- lavmodel@GLIST

    # return a list
    YHAT <- vector("list", length=lavsamplestats@ngroups)

    # compute YHAT for each group
    for(g in group) {
        # which mm belong to group g?
        mm.in.group <- 1:lavmodel@nmat[g] + cumsum(c(0L,lavmodel@nmat))[g]
        MLIST <- GLIST[ mm.in.group ]

        if(is.null(eXo[[g]]) && duplicate) {
            Nobs <- lavsamplestats@nobs[[g]]
        } else {
            Nobs <- 1L
        }

        if(lavmodel@representation == "LISREL") {
            YHAT[[g]] <- computeYHATx.LISREL(MLIST = MLIST,
                          eXo = eXo[[g]], ETA = ETA[[g]],
                          sample.mean = lavsamplestats@mean[[g]],
                          ov.y.dummy.ov.idx = lavmodel@ov.y.dummy.ov.idx[[g]],
                          ov.x.dummy.ov.idx = lavmodel@ov.x.dummy.ov.idx[[g]],
                          ov.y.dummy.lv.idx = lavmodel@ov.y.dummy.lv.idx[[g]],
                          ov.x.dummy.lv.idx = lavmodel@ov.x.dummy.lv.idx[[g]],
                          Nobs = Nobs)
        } else {
            stop("lavaan ERROR: representation ", lavmodel@representation,
                 " not supported yet.")
        }
    }

    YHAT
}

