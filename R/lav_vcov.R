Nvcov.standard <- function(lavmodel = lavmodel, lavsamplestats = NULL, 
                           lavdata = NULL, estimator = "ML", 
                           lavcache = NULL, information = "observed") {

    # compute information matrix
    if(information == "observed") {
        if(lavsamplestats@missing.flag) {
            group.weight <- FALSE
        } else {
            group.weight <- TRUE
        }
        E <- computeObservedInformation(lavmodel = lavmodel, 
            lavsamplestats = lavsamplestats, lavdata = lavdata, 
            lavcache = lavcache, estimator = estimator, type = "free",
            group.weight = group.weight)
    } else {
        E <- computeExpectedInformation(lavmodel = lavmodel, 
            lavsamplestats = lavsamplestats, lavdata = lavdata, 
            lavcache = lavcache, estimator = estimator)
    }

    # handle constraints
    if(nrow(lavmodel@con.jac) > 0L) {
        H <- lavmodel@con.jac
        inactive.idx <- attr(H, "inactive.idx")
        lambda <- lavmodel@con.lambda # lagrangean coefs
        if(length(inactive.idx) > 0L) {
            H <- H[-inactive.idx,,drop=FALSE]
            lambda <- lambda[-inactive.idx]
        }
        if(nrow(H) > 0L) {
            H0 <- matrix(0,nrow(H),nrow(H))
            H10 <- matrix(0, ncol(E), nrow(H))
            DL <- 2*diag(lambda, nrow(H), nrow(H))
            E3 <- rbind( cbind(     E,  H10, t(H)),
                         cbind(t(H10),   DL,  H0),
                         cbind(     H,   H0,  H0)  )
            NVarCov <- MASS::ginv(E3)[1:ncol(E), 1:ncol(E)]
            # FIXME: better include inactive + slacks??
        } else {
            NVarCov <- solve(E)
        }
    #} else if(lavmodel@group.w.free) {
        ## TESTING ONLY!!  should be combined with con.jac above...
        #gw.mat.idx <- which(names(lavmodel@GLIST) == "gw")
        #gw.x.idx <- unlist( lavmodel@x.free.idx[gw.mat.idx] )

        # divide by N (FIXME: why: the group weight?)
        #E[gw.x.idx, gw.x.idx] <- E[gw.x.idx, gw.x.idx] / lavsamplestats@ntotal
        #prop <- diag( E[gw.x.idx, gw.x.idx] )
        #x.prop <- numeric( ncol(E) ); x.prop[gw.x.idx] <- prop
        #p11 <- matrix(1,1,1)
        #p12 <- matrix(x.prop, 1, ncol(E))
        #p21 <- matrix(x.prop, nrow(E), 1)
        #Eplus <- rbind( cbind( p11, p12),
        #                cbind( p21, E)    )
        #Eplus.inv <- solve(Eplus)
        #NVarCov <- Eplus.inv[-1, -1]
    #    NVarCov <- solve(E)
    } else {
        NVarCov <- solve(E)
    }

    NVarCov
}

Nvcov.bootstrap <- function(lavmodel = NULL, lavsamplestats = NULL, 
                            lavoptions = NULL, lavdata = NULL,
                            lavcache = NULL, lavpartable = NULL, 
                            control=list()) {

    # number of bootstrap draws
    if(!is.null(lavoptions$bootstrap)) {
        R <- lavoptions$bootstrap
    } else {
        R <- 1000L
    }
  
    boot.type <- "ordinary"
    if(lavoptions$test == "bollen.stine") boot.type <- "bollen.stine"

    TEST <- NULL
    COEF <- bootstrap.internal(object          = NULL,
                               lavmodel.       = lavmodel, 
                               lavsamplestats. = lavsamplestats, 
                               lavpartable.    = lavpartable, 
                               lavoptions.     = lavoptions, 
                               lavdata.        = lavdata,
                               R               = R, 
                               verbose         = lavoptions$verbose,
                               type            = boot.type,
                               FUN  = ifelse(boot.type == "bollen.stine",
                                          "coeftest", "coef"),
                               warn            = -1L,
                               parallel        = control$parallel,
                               ncpus           = control$ncpus,
                               cl              = control$cl)
    if(boot.type == "bollen.stine") {
        nc <- ncol(COEF)
        TEST <- COEF[,nc]
        COEF <- COEF[,-nc]
    }

    # FIXME: cov rescale? Yes for now
    nboot <- nrow(COEF)
    NVarCov <- lavsamplestats@ntotal * (cov(COEF) * (nboot-1)/nboot )

    # save COEF and TEST (if any)
    attr(NVarCov, "BOOT.COEF") <- COEF
    attr(NVarCov, "BOOT.TEST") <- TEST
 
    NVarCov
}

Nvcov.first.order <- function(lavmodel = NULL, lavsamplestats = NULL, 
                              lavdata = NULL, lavcache = NULL, 
                              estimator = "ML") {

    B0.group <- vector("list", lavsamplestats@ngroups)

    if(estimator == "PML") {
        Delta <- computeDelta(lavmodel = lavmodel)
        Sigma.hat <- computeSigmaHat(lavmodel = lavmodel)
        TH <- computeTH(lavmodel = lavmodel)
    } else {
        Sigma.hat <- computeSigmaHat(lavmodel = lavmodel)
        Mu.hat <- computeMuHat(lavmodel = lavmodel)
        Delta <- computeDelta(lavmodel = lavmodel)
    }

    for(g in 1:lavsamplestats@ngroups) {
        if(estimator == "PML") {
            # slow approach: compute outer product of case-wise scores
            SC <- pml_deriv1(Sigma.hat = Sigma.hat[[g]],
                             TH        = TH[[g]],
                             th.idx    = lavmodel@th.idx[[g]],
                             num.idx   = lavmodel@num.idx[[g]],
                             X         = lavdata@X[[g]],
                             lavcache  = lavcache[[g]],
                             scores    = TRUE,
                             negative  = FALSE)

            # chain rule
            group.SC <- SC %*% Delta[[g]]

            # outer product
            B0.group[[g]] <- crossprod(group.SC)

            # to get NACOV instead of ACOV
            #B0.group[[g]] <- B0.group[[g]] / lavsamplestats@ntotal

            # group weights (if any)
            #group.dx <- group.w[g] * group.dx 
        } else {
            B1 <- compute.Bbeta(Sigma.hat=Sigma.hat[[g]], 
                                Mu.hat=Mu.hat[[g]],
                                lavsamplestats=lavsamplestats, 
                                lavdata=lavdata, group=g)

            B0.group[[g]] <- t(Delta[[g]]) %*% B1 %*% Delta[[g]] 

            #if(type=="free.fixed.x" && lavmodel@fixed.x && length(lavmodel@x.idx)>0) {
            #    # for resid(., type="standardized")
            #    idx.all <- which(object$free > 0 | (object$fixed.x > 0 &
            #                     object$row >= object$col))
            #    idx.free <- which(object$free > 0)
            #    idx.x <- which(object$fixed.x > 0 & 
            #                   object$row >= object$col)
            #        
            #    info.free.idx  <- idx.all %in% idx.free
            #    info.fixed.idx <- idx.all %in% idx.x
            #    B0.group[[g]][info.fixed.idx, info.free.idx] <- 0
            #    B0.group[[g]][info.free.idx, info.fixed.idx] <- 0
            #}
        }
    } # g

    if(lavsamplestats@ngroups > 1L) {
        # groups weights
        B0 <- (lavsamplestats@nobs[[1]]/lavsamplestats@ntotal) * B0.group[[1]]
        for(g in 2:lavsamplestats@ngroups) {
            B0 <- B0 + (lavsamplestats@nobs[[g]]/lavsamplestats@ntotal) * B0.group[[g]]
        }
    } else {
        B0 <- B0.group[[1]]
    }

    # handle constraints
    E <- B0
    if(nrow(lavmodel@con.jac) > 0L) {
        H <- lavmodel@con.jac
        inactive.idx <- attr(H, "inactive.idx")
        lambda <- lavmodel@con.lambda # lagrangean coefs
        if(length(inactive.idx) > 0L) {
            H <- H[-inactive.idx,,drop=FALSE]
            lambda <- lambda[-inactive.idx]
        }
        if(nrow(H) > 0L) {
            H0 <- matrix(0,nrow(H),nrow(H))
            H10 <- matrix(0, ncol(E), nrow(H))
            DL <- 2*diag(lambda, nrow(H), nrow(H))
            E3 <- rbind( cbind(     E,  H10, t(H)),
                         cbind(t(H10),   DL,  H0),
                         cbind(     H,   H0,  H0)  )
            NVarCov <- MASS::ginv(E3)[1:ncol(E), 1:ncol(E)]
            # FIXME: better include inactive + slacks??
        } else {
            # check if E is pd 
            eigvals <- eigen(E, symmetric = TRUE, only.values = TRUE)$values
            if(any(eigvals < -1 * .Machine$double.eps^(3/4))) {
                warning("lavaan WARNING: matrix based on first order outer product of the derivatives is not positive definite; the standard errors may not be trustworthy")
            }
            NVarCov <- MASS::ginv(E) ## FIXME: should we allow this?
        }
    } else {
        # NVarCov <- solve(E)
        # check if E is pd 
        eigvals <- eigen(E, symmetric = TRUE, only.values = TRUE)$values
        if(any(eigvals < -1 * .Machine$double.eps^(3/4))) {
            warning("lavaan WARNING: matrix based on first order outer product of the derivatives is not positive definite; the standard errors may not be trustworthy")
        }

        NVarCov <- MASS::ginv(E) ## FIXME: should we allow this?
    }

    attr(NVarCov, "B0") <- B0
    attr(NVarCov, "B0.group") <- B0.group

    NVarCov
}

Nvcov.robust.sem <- function(lavmodel = NULL, lavsamplestats = NULL, 
                             lavdata = NULL, lavcache = NULL,
                             estimator = "ML", mimic = "lavaan") {

    # compute information matrix
    if(estimator == "ML" && mimic == "Mplus") {
        # YR - 11 aug 2010 - what Mplus seems to do is (see Muthen apx 4 eq102)
        # - WLS.V is not based on Sigma.hat and Mu.hat (as it
        #   should be?), but on lavsamplestats@cov and lavsamplestats@mean...
        # - Gamma is not identical to what is used for WLS; closer to EQS
        # - N/N-1 bug in G11 for NVarCov (but not test statistic)
        # - we divide by N-1! (just like EQS)
        E <- computeExpectedInformationMLM(lavmodel = lavmodel, 
                                           lavsamplestats = lavsamplestats)
        Gamma <- lavsamplestats@NACOV
        # 'fix' G11 part of Gamma (NOTE: this is NOT needed for SB test 
        # statistic
        for(g in 1:lavsamplestats@ngroups) {
            gg1 <- (lavsamplestats@nobs[[g]]-1)/lavsamplestats@nobs[[g]]
            G11 <- gg1 * lavsamplestats@cov[[g]]
            Gamma[[g]][1:nrow(G11), 1:nrow(G11)] <- G11
        } # g
    } else {
        E <- computeExpectedInformation(lavmodel = lavmodel, 
                                        lavsamplestats = lavsamplestats, 
                                        lavdata = lavdata, 
                                        estimator = estimator, 
                                        extra = TRUE) 
        Gamma <- lavsamplestats@NACOV
    }

    # handle constraints
    if(nrow(lavmodel@con.jac) > 0L) {
        H <- lavmodel@con.jac
        inactive.idx <- attr(H, "inactive.idx")
        lambda <- lavmodel@con.lambda # lagrangean coefs
        if(length(inactive.idx) > 0L) {
            H <- H[-inactive.idx,,drop=FALSE]
            lambda <- lambda[-inactive.idx]
        }
        if(nrow(H) > 0L) {
            H0 <- matrix(0,nrow(H),nrow(H))
            H10 <- matrix(0, ncol(E), nrow(H))
            DL <- 2*diag(lambda, nrow(H), nrow(H))
            E3 <- rbind( cbind(     E,  H10, t(H)),
                         cbind(t(H10),   DL,  H0),
                         cbind(     H,   H0,  H0)  )
            E.inv <- MASS::ginv(E3)[1:ncol(E), 1:ncol(E)]
            # FIXME: better include inactive + slacks??
        } else {
            E.inv <- solve(E)
        }
    } else {
        E.inv <- solve(E)
    }

    Delta <- attr(E, "Delta")
    WLS.V <- attr(E, "WLS.V")

    tDVGVD <- matrix(0, ncol=ncol(E), nrow=nrow(E))
    for(g in 1:lavsamplestats@ngroups) {
        fg  <-  lavsamplestats@nobs[[g]]   /lavsamplestats@ntotal
        fg1 <- (lavsamplestats@nobs[[g]]-1)/lavsamplestats@ntotal
        # fg twice for WLS.V, 1/fg1 once for GaMMA
        # if fg==fg1, there would be only one fg, as in Satorra 1999 p.8
        # t(Delta) * WLS.V %*% Gamma %*% WLS.V %*% Delta
        if(estimator == "DWLS" || estimator == "ULS") {
            # diagonal weight matrix
            WD <- WLS.V[[g]] * Delta[[g]]
        } else {
            # full weight matrix
            WD <- WLS.V[[g]] %*% Delta[[g]]
        }
        tDVGVD <- tDVGVD + fg*fg/fg1 * crossprod(WD, Gamma[[g]] %*% WD)
    } # g
    NVarCov <- (E.inv %*% tDVGVD %*% E.inv)

    # to be reused by lavaanTest
    attr(NVarCov, "E.inv") <- E.inv
    attr(NVarCov, "Delta") <- Delta
    attr(NVarCov, "WLS.V") <- WLS.V

    NVarCov
}

Nvcov.robust.huber.white <- function(lavmodel = lavmodel, 
                                     lavsamplestats = NULL, 
                                     lavdata = NULL,
                                     information = "observed", 
                                     lavcache = NULL, 
                                     estimator="ML") {

    # compute standard Nvcov
    E.inv <- Nvcov.standard(lavmodel       = lavmodel,
                            lavsamplestats = lavsamplestats,
                            lavdata        = lavdata,
                            estimator      = estimator,
                            lavcache       = lavcache,
                            information    = information)

    # compute first.order Nvcov
    Nvcov <- Nvcov.first.order(lavmodel       = lavmodel,
                               lavsamplestats = lavsamplestats,
                               lavdata        = lavdata,
                               lavcache       = lavcache,
                               estimator      = estimator)
    B0 <- attr(Nvcov, "B0")
    B0.group <- attr(Nvcov, "B0.group")

    # compute sandwich
    NVarCov <- E.inv %*% B0 %*% E.inv

    attr(NVarCov, "B0.group") <- B0.group
    attr(NVarCov, "E.inv") <- E.inv

    NVarCov
}


lav_model_vcov <- function(lavmodel       = NULL, 
                           lavsamplestats = NULL, 
                           lavoptions     = NULL, 
                           lavdata        = NULL, 
                           lavpartable    = NULL, 
                           lavcache       = NULL, 
                           control=list()) {

    estimator   <- lavoptions$estimator
    likelihood  <- lavoptions$likelihood
    information <- lavoptions$information
    se          <- lavoptions$se
    verbose     <- lavoptions$verbose
    mimic       <- lavoptions$mimic
  
    if(se == "none") return(NULL)

    # some require meanstructure (for now)
    if(se %in% c("first.order", "robust.sem", "robust.huber.white") && 
       !lavoptions$meanstructure) {
        stop("se (", se, ") requires meanstructure (for now)")
    }

    if(se == "standard") {
        NVarCov <- try( Nvcov.standard(lavmodel       = lavmodel,
                                       lavsamplestats = lavsamplestats,
                                       lavdata        = lavdata,
                                       estimator      = estimator, 
                                       lavcache       = lavcache,
                                       information    = information),
                         silent=TRUE )

    } else if(se == "first.order") {
        NVarCov <- try( Nvcov.first.order(lavmodel       = lavmodel,
                                          lavsamplestats = lavsamplestats,
                                          lavdata        = lavdata,
                                          lavcache       = lavcache,
                                          estimator      = estimator),
                         silent=TRUE )

    } else if(se == "robust.sem") {
        NVarCov <- try( Nvcov.robust.sem(lavmodel       = lavmodel,
                                         lavsamplestats = lavsamplestats,
                                         estimator      = estimator,
                                         mimic          = mimic,
                                         lavcache       = lavcache,
                                         lavdata        = lavdata),
                        silent=TRUE )

    } else if(se == "robust.huber.white") {
        NVarCov <- try( Nvcov.robust.huber.white(lavmodel = lavmodel,
                                           lavsamplestats = lavsamplestats,
                                           lavdata        = lavdata,
                                           information    = information,
                                           lavcache       = lavcache,
                                           estimator      = estimator),
                        silent=TRUE )

    } else if(se == "bootstrap") {
        NVarCov <- try( Nvcov.bootstrap(lavmodel       = lavmodel,
                                        lavsamplestats = lavsamplestats,
                                        lavoptions     = lavoptions,
                                        lavdata        = lavdata,
                                        lavcache       = lavcache,
                                        lavpartable    = lavpartable,
                                        control        = control),
                        silent=TRUE )
    }

    if(! inherits(NVarCov, "try-error") ) {

        # denominator!
        if(estimator %in% c("ML","PML","FML") && likelihood == "normal") {
            N <- lavsamplestats@ntotal
        } else {
            N <- lavsamplestats@ntotal - lavsamplestats@ngroups
        }

        if(estimator %in% c("PML", "MML")) {
            VarCov <- NVarCov
        } else {
            VarCov <- 1/N * NVarCov
        }

    } else {
        warning("lavaan WARNING: could not compute standard errors!\n  lavaan NOTE: this may be a symptom that the model is not identified.\n")
        VarCov <- NULL
    }

    VarCov
}


