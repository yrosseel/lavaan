# model gradient

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
#    Delta <- lavJacobianC(func=compute.moments, x = x)
#
#    Delta
#}


### FIXME: should we here also:
###        - weight for groups? (no, for now)
###        - handle equality constraints? (yes, for now)
computeDelta <- function(object, GLIST.=NULL, m.el.idx.=NULL, x.el.idx.=NULL) {

    representation   <- object@representation
    categorical      <- object@categorical
    group.w.free     <- object@group.w.free
    nmat             <- object@nmat
    ngroups          <- object@ngroups
    nvar             <- object@nvar
    num.idx          <- object@num.idx
    th.idx           <- object@th.idx
    nexo             <- object@nexo
    parameterization <- object@parameterization

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

    ################################################## 
    ##################################################
    # special treatment: parameterization == "theta"
    # we do it numerically (for now) -- FIXME!!!
#    if(object@parameterization == "theta") {
#        compute.moments <- function(x, g=1L) {

#            # which mm belong to group g?
#            mm.in.group <- 1:nmat[g] + cumsum(c(0,nmat))[g]
#            GLIST <- x2GLIST(object, x=x, type="free")
#            MLIST <- GLIST[mm.in.group]
    
   
#            # 1. TH
#            out <- computeTH.LISREL(MLIST = MLIST, th.idx = th.idx[[g]])
     
#            # 2. PI
#            if(object@nexo[g] > 0L) {
#                PI <- computePI.LISREL(MLIST = MLIST)
#                out <- c(out, as.numeric(PI))
#            }
    
#           # 3. Sigma.hat
#            Sigma.hat <- computeSigmaHat.LISREL(MLIST = MLIST, delta = TRUE)
#            # reorder: first variances (of numeric), then covariances
#            var.num <- diag(Sigma.hat)[num.idx[[g]]]    
#            COR <- Sigma.hat[ vech.idx(nvar[g], diagonal=FALSE) ]
#            out <- c(out, var.num, COR)
#    
#           out
#       }
# 
#        Delta <- vector("list", length=ngroups)    
#        for(g in 1:ngroups) {
#            x <- getModelParameters(object, GLIST=GLIST, type="free")
#            Delta[[g]] <- lavJacobianC(func=compute.moments, x=x, g=g)
#        }        
#     
#        return(Delta)
#    }
    ################################################## 
    ##################################################



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
        if(group.w.free) {
            pstar[g] <- pstar[g] + 1L # add group weight
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

        # if theta, do some preparation
        if(parameterization == "theta") {
            sigma.hat <- computeSigmaHat.LISREL(MLIST=GLIST[mm.in.group], 
                                                delta=FALSE)
            dsigma <- diag(sigma.hat)
            # dcor/dcov for sigma
            R <- lav_deriv_cov2cor(sigma.hat, num.idx = object@num.idx[[g]])

            theta.var.idx <- which(!vech.idx(nvar[g]) %in%
                                    vech.idx(nvar[g], diagonal=FALSE))
        }

        for(mm in mm.in.group) {
            mname <- names(object@GLIST)[mm]

            # skip empty ones
            if(!length(m.el.idx[[mm]])) next

            # get Delta columns for this model matrix
            if(representation == "LISREL") {

                # Sigma
                DELTA <- dxSigma <-
                    derivative.sigma.LISREL(m = mname,
                                            idx = m.el.idx[[mm]],
                                            MLIST = GLIST[ mm.in.group ],
                                            delta = parameterization == "delta")
                if(categorical && parameterization == "theta") {
                    DELTA <- R %*% DELTA
                }
                
                if(categorical) {
                    # reorder: first variances (of numeric), then covariances
                    cov.idx  <- vech.idx(nvar[g])
                    covd.idx <- vech.idx(nvar[g], diagonal=FALSE)

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
                                                     MLIST=GLIST[ mm.in.group ],
                                                     delta = TRUE)
                    if(parameterization == "theta") {
                        # dy/ddsigma = -0.5/(ddsigma*sqrt(ddsigma))
                        dDelta.dx <- 
                            ( dxSigma[theta.var.idx,,drop=FALSE] * 
                              -0.5 / (dsigma*sqrt(dsigma)) )
                        dth.dDelta <- 
                            derivative.th.LISREL(m = "delta",
                                                 idx = 1:nvar[g],
                                                 MLIST = GLIST[ mm.in.group ],
                                                 th.idx = th.idx[[g]])
                        # add dth.dDelta %*% dDelta.dx
                        no.num.idx <- which(th.idx[[g]] > 0)
                        DELTA.th[no.num.idx,] <- 
                            DELTA.th[no.num.idx,,drop=FALSE] +
                            (dth.dDelta %*% dDelta.dx)[no.num.idx,,drop=FALSE]
                    }
                    if(object@nexo[g] > 0L) {
                        DELTA.pi <- 
                            derivative.pi.LISREL(m=mname,
                                                 idx=m.el.idx[[mm]],
                                                 MLIST=GLIST[ mm.in.group ])
                        if(parameterization == "theta") {
                            dpi.dDelta <-
                                derivative.pi.LISREL(m = "delta",
                                    idx = 1:nvar[g],
                                    MLIST = GLIST[ mm.in.group ])
                            # add dpi.dDelta %*% dDelta.dx
                            no.num.idx <- which(!seq.int(1L,nvar) %in% num.idx)
                            no.num.idx <- rep(seq.int(0,nexo-1) * nvar,
                                          each=length(no.num.idx)) + no.num.idx
                            DELTA.pi[no.num.idx,] <- 
                                DELTA.pi[no.num.idx,,drop=FALSE] +
                                (dpi.dDelta %*% dDelta.dx)[no.num.idx,,drop=FALSE]
                        }
                        DELTA <- rbind(DELTA.th, DELTA.pi, DELTA)
                    } else {
                        DELTA <- rbind(DELTA.th, DELTA)
                    }
                }
                if(group.w.free) {
                    DELTA.gw <- derivative.gw.LISREL(m=mname,
                                                     idx=m.el.idx[[mm]],
                                                     MLIST=GLIST[ mm.in.group ])
                    DELTA <- rbind(DELTA.gw, DELTA)
                }
            } else {
                stop("representation", representation, "not implemented yet")
            }

            Delta.group[ ,x.el.idx[[mm]]] <- DELTA
        } # mm

        # save(Delta.group, file=paste0("delta_NO_EQ",g,".Rdata"))

        # if type == "free" take care of equality constraints
        if(type == "free" && object@eq.constraints) {
            Delta.group <- Delta.group %*% object@eq.constraints.K
        }
  
        #Delta.eq <- Delta.group
        # save(Delta.eq, file=paste0("delta_NO_EQ",g,".Rdata"))

        Delta[[g]] <- Delta.group

    } # g

    Delta
}

computeDeltaDx <- function(object, GLIST=NULL, target="lambda") {

    # state or final?
    if(is.null(GLIST)) GLIST <- object@GLIST

    representation   <- object@representation
    nmat             <- object@nmat
    ngroups          <- object@ngroups
    num.idx          <- object@num.idx
    th.idx           <- object@th.idx

    # number of columns in DELTA + m.el.idx/x.el.idx
    type <- "free"
    #if(type == "free") {
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
    #} else {
    #    NCOL <- sum(unlist(lapply(x.el.idx, function(x) length(unique(x)))))
    #}

    # compute Delta per group
    Delta <- vector("list", length=ngroups)
    for(g in 1:ngroups) {
        mm.in.group <- 1:nmat[g] + cumsum(c(0,nmat))[g]
        Delta.group <- NULL
        for(mm in mm.in.group) {
            mname <- names(object@GLIST)[mm]

            # skip empty ones
            if(!length(m.el.idx[[mm]])) next

            # get Delta columns for this model matrix
            if(representation == "LISREL") {
                if(target == "lambda") {
                    DELTA <- derivative.lambda.LISREL(m=mname,
                               idx=m.el.idx[[mm]], MLIST=GLIST[ mm.in.group ])
                } else if(target == "th") {
                    DELTA <- derivative.th.LISREL(m=mname, th.idx = th.idx[[g]],
                               idx=m.el.idx[[mm]], MLIST=GLIST[ mm.in.group ],
                               delta=TRUE)
                } else if(target == "mu") {
                    DELTA <- derivative.mu.LISREL(m=mname,
                               idx=m.el.idx[[mm]], MLIST=GLIST[ mm.in.group ])
                } else if(target == "theta") {
                    DELTA <- derivative.theta.LISREL(m=mname, 
                               idx=m.el.idx[[mm]], MLIST=GLIST[ mm.in.group ])
                } else if(target == "sigma") {
                    DELTA <- derivative.sigma.LISREL(m=mname,
                               idx=m.el.idx[[mm]], MLIST=GLIST[ mm.in.group ],
                               delta=TRUE)
                } else {
                    stop("lavaan ERROR: target ", target, " not implemented yet")
                }

                # initialize?
                if(is.null(Delta.group)) {
                    Delta.group <- matrix(0, nrow=nrow(DELTA), ncol=NCOL)
                }
                Delta.group[ ,x.el.idx[[mm]]] <- DELTA
            }
        } # mm

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
                warning("computeGradient: Sigma.hat is not positive definite\n")
                #Sigma.hat[[g]] <- force.pd(Sigma.hat[[g]])
                Sigma.hat.inv <- MASS::ginv(Sigma.hat[[g]])
                Sigma.hat.log.det <- log(.Machine$double.eps)
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
                            X=NULL, cache=NULL, type="free", 
                            estimator="ML", verbose=FALSE, forcePD=TRUE, 
                            group.weight=TRUE, constraints=TRUE,
                            Delta=NULL) {

    nmat           <- object@nmat
    representation <- object@representation
    meanstructure  <- object@meanstructure
    categorical    <- object@categorical
    group.w.free   <- object@group.w.free
    fixed.x        <- object@fixed.x
    num.idx        <- object@num.idx
    th.idx         <- object@th.idx
    nx.unco        <- object@nx.unco

    # state or final?
    if(is.null(GLIST)) GLIST <- object@GLIST

    # group.weight
    if(group.weight) {
        if(estimator %in% c("ML","PML","FML","MML")) {
            group.w <- (unlist(samplestats@nobs)/samplestats@ntotal)
        } else {
            # FIXME: double check!
            group.w <- ((unlist(samplestats@nobs)-1)/samplestats@ntotal)
        }
    } else {
        group.w <- rep(1.0, samplestats@ngroups)
    }

    # do we need WLS.est?
    if(estimator == "WLS"  || estimator == "DWLS" || estimator == "ULS") {
        WLS.est <- lav_model_wls_est(object = object, GLIST = GLIST)
    } else {
        # compute moments for all groups
        Sigma.hat <- computeSigmaHat(object, GLIST=GLIST,
                                     extra=(estimator=="ML"))
        if(meanstructure && !categorical) {
            Mu.hat <- computeMuHat(object, GLIST=GLIST)
        } else if(categorical) {
            TH <- computeTH(object, GLIST=GLIST)
            if(fixed.x)
                 PI <- computePI(object, GLIST=GLIST)
        }
        if(group.w.free) {
            GW <- computeGW(object, GLIST=GLIST)
        }
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

        # extract free parameters
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
            diff <- as.matrix(samplestats@WLS.obs[[g]]  - WLS.est[[g]])
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

    else if(estimator == "PML" || estimator == "FML" ||
            estimator == "MML") {

        if(type != "free") {
            stop("FIXME: type != free in computeGradient for estimator PML")
        } else {
            Delta <- computeDelta(object, GLIST.=GLIST)
        }

        for(g in 1:samplestats@ngroups) {

            #print(GLIST)
            #print(getModelParameters(object, GLIST=GLIST))
            #print(Sigma.hat[[g]])
            #print(TH[[g]])
            #cat("*****\n")

            # compute partial derivative of logLik with respect to 
            # thresholds/means, slopes, variances, correlations
            if(estimator == "PML") {
                d1 <- pml_deriv1(Sigma.hat = Sigma.hat[[g]],
                                 TH        = TH[[g]],
                                 th.idx    = th.idx[[g]],
                                 num.idx   = num.idx[[g]],
                                 X         = X[[g]],
                                 cache     = cache[[g]])
            } else {
                d1 <- fml_deriv1(Sigma.hat = Sigma.hat[[g]],
                                 TH        = TH[[g]],
                                 th.idx    = th.idx[[g]],
                                 num.idx   = num.idx[[g]],
                                 X         = X[[g]],
                                 cache     = cache[[g]])
            }

            # chain rule (logLik)
            #group.dx <- as.numeric(t(d1) %*% Delta[[g]])

            # chain rule (fmin)
            ### FIXME why -1L ???
            group.dx <- as.numeric(t(d1) %*% Delta[[g]])/samplestats@nobs[[g]]

            # group weights (if any)
            group.dx <- group.w[g] * group.dx
            if(g == 1) {
                dx <- group.dx
            } else {
                dx <- dx + group.dx
            }
        } # g
    } else {
        stop("lavaan ERROR: no analytical gradient available for estimator ",
             estimator)
    }


    # group.w.free for ML
    if(object@group.w.free && estimator %in% c("ML","MML","FML","PML")) {
        #est.prop <- unlist( computeGW(object, GLIST=GLIST) )
        #obs.prop <- unlist(samplestats@group.w)
        # FIXME: G2 based -- ML and friends only!!
        #dx.GW <- - (obs.prop - est.prop)

        # poisson version
        est.freq <- exp(unlist(computeGW(object, GLIST=GLIST)))
        obs.freq <- unlist(samplestats@group.w) * samplestats@ntotal
        dx.GW <- - (obs.freq - est.freq)
        # divide by N (to be consistent with the rest of lavaan)
        dx.GW <- dx.GW / samplestats@ntotal

        # remove last element (fixed LAST group to zero)
        # dx.GW <- dx.GW[-length(dx.GW)]
        
        # fill in in dx
        gw.mat.idx <- which(names(object@GLIST) == "gw")
        gw.x.idx <- unlist( object@x.free.idx[gw.mat.idx] )
        dx[gw.x.idx] <- dx.GW
    }

    dx
}

