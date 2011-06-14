# Model-methods

#
# initial version: YR 25/03/2009: `methods' for the Model class

getModelParameters <- function(object, GLIST=NULL, type="free") {

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

    x
}

# warning: this will make a copy of object
setModelParameters <- function(object, x=NULL) {

    names.GLIST <- names(object@GLIST)

    tmp <- object@GLIST
    for(mm in 1:length(object@GLIST)) {
        m.free.idx <- object@m.free.idx[[mm]]
        x.free.idx <- object@x.free.idx[[mm]]
        tmp[[mm]][m.free.idx] <- x[x.free.idx]
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
                m.el.idx <- vecs.idx(N)
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

setMethod("computeSigmaHat", "Model",
function(object, GLIST=NULL, extra=FALSE) {

    # state or final?
    if(is.null(GLIST)) GLIST <- object@GLIST

    nmat           <- object@nmat
    ngroups        <- object@ngroups
    representation <- object@representation

    # return a list
    Sigma.hat <- vector("list", length=ngroups)

    for(g in 1:ngroups) {

        # which mm belong to group g?
        mm.in.group <- nmat * (g - 1L) + 1:nmat

        if(representation == "LISREL") {
            Sigma.hat[[g]] <- computeSigmaHat.LISREL(MLIST = GLIST[mm.in.group])
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
})

setMethod("computeMuHat", "Model",
function(object, GLIST=NULL) {

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
        mm.in.group <- nmat * (g - 1L) + 1:nmat

        if(!meanstructure) {
            Mu.hat[[g]] <- numeric( object@nvar )
        } else
        if(representation == "LISREL") {
            Mu.hat[[g]] <- computeMuHat.LISREL(MLIST = GLIST[ mm.in.group ])
        } else {
            stop("only representation LISREL has been implemented for now")
        }
    } # ngroups

    Mu.hat
})

setMethod("computeObjective", "Model",
function(object, GLIST=NULL, sample, estimator="ML", 
         verbose=FALSE, forcePD=TRUE) {

    meanstructure <- object@meanstructure

    # compute moments for all groups
    Sigma.hat <- computeSigmaHat(object, GLIST=GLIST, extra=(estimator=="ML"))
    if(meanstructure) Mu.hat <- computeMuHat(object, GLIST=GLIST)
    
    fx <- 0.0
    fx.group <- numeric( sample@ngroups )
    for(g in 1:sample@ngroups) {

        # incomplete data and missing.flag[g]=TRUE?
        if(sample@missing.flag[g]) {
            if(estimator == "ML") {
                # FIML

                # catch non-pd Sigma.hat
                # 0.4-5: instead of using force.pd, we return Inf
                if(!attr(Sigma.hat[[g]], "po")) return(Inf)

                group.fx <- estimator.FIML(Sigma.hat=Sigma.hat[[g]],
                                           Mu.hat=Mu.hat[[g]],
                                           M=sample@missing[[g]])
            } else {
                stop("this estimator: `", estimator, 
                     "' can not be used with incomplete data and the missing=\"ml\" option")
            }
        } else if(estimator == "ML") {
        # complete data
            # ML and friends
            group.fx <- estimator.ML(Sigma.hat=Sigma.hat[[g]], 
                                     Mu.hat=Mu.hat[[g]],
                                     data.cov=sample@cov[[g]], 
                                     data.mean=sample@mean[[g]], 
                                     data.cov.log.det=sample@cov.log.det[[g]],
                                     meanstructure=meanstructure)
        } else if(estimator == "GLS" || estimator == "WLS") {
            group.fx <- estimator.WLS(Sigma.hat=Sigma.hat[[g]], 
                                      Mu.hat=Mu.hat[[g]],
                                      w.vecs=sample@cov.vecs[[g]], 
                                      data.mean=sample@mean[[g]],
                                      WLS.V=sample@WLS.V[[g]],  
                                      meanstructure=meanstructure)
        } else {
            stop("unsupported estimator: ", estimator)
        }

        if(estimator == "ML") {
            group.fx <- 0.5 * group.fx
        } else {
            group.fx <- 0.5 * (sample@nobs[[g]]-1)/sample@nobs[[g]] * group.fx
        }

        fx.group[g] <- group.fx
    } # g

    if(sample@ngroups > 1) {
        nobs <- unlist(sample@nobs)
        fx <- weighted.mean(fx.group, w=nobs)
    } else { # single group
        fx <- fx.group[1]
    }

    attr(fx, "fx.group") <- fx.group

    fx
})

# for testing purposes only
#computeDeltaNumerical <- function(object, GLIST=NULL, g=1) {
#
#    # state or final?
#    if(is.null(GLIST)) GLIST <- object@GLIST
#    
#    compute.moments <- function(x) {
#        GLIST <- x2GLIST(object, x=x, type="free")
#        Sigma.hat <- computeSigmaHat(object, GLIST=GLIST)
#        S.vec <- vecs(Sigma.hat[[g]])
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
computeDelta <- function(object, GLIST=NULL, m.el.idx=NULL, x.el.idx=NULL) {

    representation <- object@representation
    nmat           <- object@nmat
    ngroups        <- object@ngroups
    nvar           <- object@nvar

    # state or final?
    if(is.null(GLIST)) GLIST <- object@GLIST

    # type = "free" or something else?
    type <- "nonfree"
    if(is.null(m.el.idx) && is.null(x.el.idx)) type <- "free"

    # number of rows in DELTA
    pstar <- nvar * (nvar + 1) / 2
    if(object@meanstructure) {
        pstar <- nvar + pstar  # first the means, then sigma
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
        Delta.group <- matrix(0, nrow=pstar, ncol=NCOL)

        # which mm belong to group g?
        mm.in.group <- nmat * (g - 1L) + 1:nmat

        for(mm in mm.in.group) {
            mname <- names(object@GLIST)[mm]

            # skip empty ones
            if(!length(m.el.idx[[mm]])) next

            # get Delta columns for this model matrix
            if(representation == "LISREL") {
                DELTA <- derivative.sigma.LISREL(m=mname,
                                                 idx=m.el.idx[[mm]],
                                                 MLIST=GLIST[ mm.in.group ])
                if(object@meanstructure) {
                    DELTA.mu <- derivative.mu.LISREL(m=mname,
                                                     idx=m.el.idx[[mm]],
                                                     MLIST=GLIST[ mm.in.group ])
                    DELTA <- rbind(DELTA.mu, DELTA)
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
                         sample=NULL, estimator="ML", meanstructure=FALSE) {

    Omega    <- vector("list", length=sample@ngroups)
    Omega.mu <- vector("list", length=sample@ngroups)

    for(g in 1:sample@ngroups) {

        # ML
        if(estimator == "ML") {

            if(attr(Sigma.hat[[g]], "po") == FALSE) {
                # FIXME: WHAT IS THE BEST THING TO DO HERE??
                # CURRENTLY: force matrix to be POS DEFINITE if requested
                # but with a warning
                warning("computeGradient: Sigma.hat is not positive definite\n")
                Sigma.hat[[g]] <- force.pd(Sigma.hat[[g]])
                Sigma.hat.inv <- inv.chol(Sigma.hat[[g]], logdet=TRUE)
                Sigma.hat.log.det <- attr(Sigma.hat.inv, "logdet")
            } else {
                Sigma.hat.inv <-  attr(Sigma.hat[[g]], "inv")
                Sigma.hat.log.det <- attr(Sigma.hat[[g]], "log.det")
            }

            if(!sample@missing.flag[g]) { # complete data
                if(meanstructure) {
                    diff <- sample@mean[[g]] - Mu.hat[[g]]
                    W.tilde <- sample@cov[[g]] + tcrossprod(diff)
                    # Browne 1995 eq 4.55
                    Omega.mu[[g]] <- t(t(diff) %*% Sigma.hat.inv)
                    Omega[[g]] <- 
                        ( Sigma.hat.inv %*% (W.tilde - Sigma.hat[[g]]) %*%
                          Sigma.hat.inv )
                } else {
                    W.tilde <- sample@cov[[g]]
                    Omega[[g]] <- 
                        ( Sigma.hat.inv %*% (W.tilde - Sigma.hat[[g]]) %*%
                          Sigma.hat.inv )
                }
            } else { # missing data
                M <- sample@missing[[g]]

                OMEGA    <- matrix(0, sample@nvar, sample@nvar)
                OMEGA.MU <- matrix(0, sample@nvar, 1)

                for(p in 1:M$npatterns) {
                    SX <- M$data[[p]][["SX"]]
                    MX <- M$data[[p]][["MX"]]
                    nobs <- M$data[[p]][["nobs"]]
                    var.idx <- M$data[[p]][["var.idx"]]

                    Sigma.inv <- inv.chol(Sigma.hat[[g]][var.idx, var.idx],
                                          logdet=FALSE)
                    Mu <- Mu.hat[[g]][var.idx]
                    W.tilde <- SX + tcrossprod(MX - Mu)

                    OMEGA.MU[var.idx, 1] <-
                        ( OMEGA.MU[var.idx, 1] + nobs/sample@ntotal *
                          t(t(MX - Mu) %*% Sigma.inv) )

                    OMEGA[var.idx, var.idx] <-
                        ( OMEGA[var.idx, var.idx] + nobs/sample@ntotal *
                          (Sigma.inv %*%
                           (W.tilde - Sigma.hat[[g]][var.idx,var.idx]) %*%
                           Sigma.inv ) )
                }
                Omega.mu[[g]] <- OMEGA.MU
                Omega[[g]]    <- OMEGA
            } # missing

        # GLS
        } else if(estimator == "GLS") {
            W.inv <- sample@icov[[g]]
            W     <- sample@cov[[g]]
            Omega[[g]] <- (sample@nobs[[g]]-1)/sample@nobs[[g]] *
                              (W.inv %*% (W - Sigma.hat[[g]]) %*% W.inv)
            if(meanstructure) {
                diff <- as.matrix(sample@mean[[g]] - Mu.hat[[g]])
                Omega.mu[[g]] <- t( t(diff) %*% W.inv )
            }
        }

    } # g

    if(meanstructure) attr(Omega, "mu") <- Omega.mu

    Omega
}


setMethod("computeGradient", "Model",
function(object, GLIST=NULL, sample=NULL, type="free", 
         estimator="ML", verbose=FALSE, forcePD=TRUE, 
         group.weight=TRUE, constraints=TRUE) {

    nmat           <- object@nmat
    representation <- object@representation
    meanstructure  <- object@meanstructure
    nx.unco  <- object@nx.unco

    # state or final?
    if(is.null(GLIST)) GLIST <- object@GLIST

    # group.w
    if(group.weight) {
        group.w <- (unlist(sample@nobs)/sample@ntotal)
    } else {
        group.w <- rep(1.0, sample@ngroups)
    }

    # Sigma.hat + Mu.hat
    Sigma.hat <- computeSigmaHat(object, GLIST=GLIST, extra=(estimator == "ML"))
    if(meanstructure) Mu.hat <- computeMuHat(object, GLIST=GLIST)
    

    # two approaches:
    # - ML/GLS approach: using Omega (and Omega.mu)
    # - WLS: using Delta

    # 1. ML/GLS approach
    if(estimator == "ML" || estimator == "GLS") {
        if(meanstructure) {
            Omega <- computeOmega(Sigma.hat=Sigma.hat, Mu.hat=Mu.hat,
                                  sample=sample, estimator=estimator, 
                                  meanstructure=TRUE)
            Omega.mu <- attr(Omega, "mu")
        } else {
            Omega <- computeOmega(Sigma.hat=Sigma.hat, Mu.hat=NULL,
                                  sample=sample, estimator=estimator,
                                  meanstructure=FALSE)
            Omega.mu <- vector("list", length=sample@ngroups)
        }

        # compute DX (for all elements in every model matrix)
        DX <- vector("list", length=length(GLIST))
      
        for(g in 1:sample@ngroups) {
            # which mm belong to group g?
            mm.in.group <- nmat * (g - 1L) + 1:nmat
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
            if(sample@ngroups > 1L) {
                for(mm in mm.in.group) {
                    DX[[mm]] <- group.w[g] * DX[[mm]]
                }
            }
        }

        # extract free parameters + weight by group
        if(type == "free") {
            dx <- numeric( nx.unco )
            for(g in 1:sample@ngroups) {
                mm.in.group <- nmat * (g - 1L) + 1:nmat
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
    if(estimator == "WLS") {

        if(type != "free") {
            stop("FIXME: WLS gradient with type != free needs fixing!")
        }
        Delta <- computeDelta(object, GLIST=GLIST)

        for(g in 1:sample@ngroups) {
            # Browne & Arminger 1995 eq 4.49
            if(!meanstructure) {
                obs <- sample@cov.vecs[[g]]
                est <- vecs(Sigma.hat[[g]])
            } else {
                obs <- c(sample@mean[[g]], sample@cov.vecs[[g]])
                est <- c(Mu.hat[[g]], vecs(Sigma.hat[[g]]))
            }
            diff <- as.matrix(obs - est)
            group.dx <- -1 * ( t(Delta[[g]]) %*% sample@WLS.V[[g]] %*% diff)
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

    dx
})


setMethod("estimateModel", "Model",
function(object, sample, do.fit=TRUE, options=NULL) {

    estimator     <- options$estimator
    verbose       <- options$verbose
    debug         <- options$debug


    # function to be minimized
    minimize.this.function <- function(x, verbose=FALSE) {
      
        #cat("DEBUG: x = ", x, "\n")

        # current strategy: forcePD is by default FALSE, except
        # if missing patterns are used
        if(any(sample@missing.flag)) {
            forcePD <- TRUE
        } else {
            forcePD <- FALSE
        }

        # update GLIST (change `state') and make a COPY!
        GLIST <- x2GLIST(object, x=x)

        fx <- computeObjective(object, GLIST=GLIST, sample,
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

        fx
    }

    first.derivative.param <- function(x, verbose=FALSE) {

        # update GLIST (change `state') and make a COPY!
        GLIST <- x2GLIST(object, x=x)

        dx <- computeGradient(object, GLIST=GLIST, sample, 
                              type="free", 
                              estimator=estimator,
                              verbose=verbose, forcePD=TRUE)

        if(debug) {
            cat("Gradient function (analytical) =\n"); print(dx); cat("\n")
        }

        dx
    } 

    first.derivative.param.numerical <- function(x, verbose=FALSE) {
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

    if(do.fit) {
        iter.max <- 10000
    } else {
        iter.max <- 0
    }

    INIT_NELDER_MEAD <- FALSE
    if(is.null(body(object@ceq.function)) && 
       is.null(body(object@cin.function)) ) {
        OPTIMIZER <- "NLMINB"
        #OPTIMIZER <- "BFGS"      # slightly slower, no bounds; better scaling!
        #OPTIMIZER <- "L-BFGS-B"  # trouble with Inf values for fx!
    } else {
        OPTIMIZER <- "NLMINB.CONSTR"
        #OPTIMIZER <- "ALABAMA"
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
        trace <- 0L; if(debug) trace <- 1L;
        optim.out <- nlminb(start=start.x,
                            objective=minimize.this.function,
                            gradient=first.derivative.param,
                            #gradient=first.derivative.param.numerical,
                            control=list(iter.max=iter.max,
                                         eval.max=iter.max*2, 
                                         trace=trace),
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

        trace <- 0L; if(verbose) trace <- 1L
        optim.out <- optim(par=start.x,
                           fn=minimize.this.function,
                           gr=first.derivative.param,
                           method="BFGS",
                           control=list(maxit=iter.max, REPORT=1L, 
                                        parscale=SCALE,
                                        trace=trace, reltol=1e-12),
                           hessian=FALSE)
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

        trace <- 0L; if(verbose) trace <- 1L
        optim.out <- optim(par=start.x,
                           fn=minimize.this.function,
                           gr=first.derivative.param,
                           method="L-BFGS-B",
                           control=list(maxit=iter.max, REPORT=1L, 
                                        trace=trace, factr=1e-12),
                           hessian=FALSE)
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
    } else if(OPTIMIZER == "BBoptim") {

        trace <- FALSE; if(verbose) trace <- TRUE
        optim.out <- BBoptim(par=start.x,
                             fn=minimize.this.function,
                             gr=first.derivative.param,
                             control=list(maxit=iter.max, triter=1L,
                                          trace=trace),
                             quiet=TRUE
                            )
        if(verbose) {
            cat("convergence status (0=ok): ", optim.out$convergence, "\n")
            cat("optim BBoptim message says: ", optim.out$message, "\n")
            cat("number of iterations: ", optim.out$iterations, "\n")
            cat("number of function evaluations [objective]: ",
                optim.out$feval, "\n")
        }

        iterations <- as.integer(optim.out$iter)
        x          <- optim.out$par
        if(optim.out$convergence == 0L) {
            converged <- TRUE
        } else {
            converged <- FALSE
        }
    } else if(OPTIMIZER == "ALABAMA") {
        #require(alabama)
        hin <- hin.jac <- heq <- heq.jac <- NULL
        if(!is.null(body(object@cin.function))) hin     <- object@cin.function
        if(!is.null(body(object@cin.jacobian))) hin.jac <- object@cin.jacobian
        if(!is.null(body(object@ceq.function))) heq     <- object@ceq.function
        if(!is.null(body(object@ceq.jacobian))) heq.jac <- object@ceq.jacobian
        trace <- FALSE; if(verbose) trace <- TRUE
        optim.out <- auglag(par=start.x,
                            fn=minimize.this.function,
                            gr=first.derivative.param,
                            #hin=hin,
                            #hin.jac=NULL,
                            heq=heq,
                            #heq.jac=heq.jac,
                            control.outer=list(trace=trace, method="BFGS"),
                            control.optim=list(maxit=iter.max, REPORT=1L,
                                               parscale=SCALE,
                                               trace=trace, reltol=1e-12),
                           )
        if(verbose) {
            cat("convergence status (0=ok): ", optim.out$convergence, "\n")
            cat("optim auglag message says: ", optim.out$message, "\n")
            cat("number of outer iterations: ", 
                optim.out$outer.iterations, "\n")
            cat("number of function evaluations [objective gradient]: ",
                optim.out$counts, "\n")
        }

        iterations <- as.integer(optim.out$counts[1])
        x          <- optim.out$par
        if(optim.out$convergence == 0L) {
            converged <- TRUE
        } else {
            converged <- FALSE
        }
    } else if(OPTIMIZER == "NLMINB.CONSTR") {
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
                                   control=list(iter.max=iter.max,
                                                eval.max=iter.max*2,
                                                trace=trace),
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

    # adjust fx if FIML (using h1 value)
    if(any(sample@missing.flag)) {
        fx.group <- attr(fx,"fx.group") 
        h1 <- lapply(sample@missing, "[[", "h1")
        # in case we have a complete group
        h1[unlist(lapply(h1, is.null))] <- 0.0
        fx.group <- fx.group - unlist(h1)/2.0
        fx <- sum(fx.group)
        attr(fx, "fx.group") <- fx.group
    }

    attr(x, "converged")  <- converged
    attr(x, "iterations") <- iterations
    attr(x, "fx")         <- fx
    if(!is.null(optim.out$con.jac)) attr(x, "con.jac") <- optim.out$con.jac

    x

})


