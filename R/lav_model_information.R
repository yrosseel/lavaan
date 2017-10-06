# here, we compute various versions of the `information' matrix
# NOTE:
# 1) we ALWAYS compute the UNIT information (not the total information)
# 
# 2) by default, we ignore the constraints (we deal with this when we
#    take the inverse later on)

lav_model_information <- function(lavmodel       = NULL,
                                  lavsamplestats = NULL,
                                  lavdata        = NULL,
                                  Delta          = NULL,
                                  lavcache       = NULL,
                                  lavoptions     = NULL,
                                  extra          = FALSE,
                                  augmented      = FALSE,
                                  inverted       = FALSE,
                                  use.ginv       = FALSE) {

    estimator   <- lavmodel@estimator
    information <- lavoptions$information

    # compute information matrix
    if(information == "observed") {
        if(lavsamplestats@missing.flag || lavdata@nlevels > 1L) {
            group.weight <- FALSE
        } else {
            group.weight <- TRUE
        }
        E <- lav_model_information_observed(lavmodel = lavmodel,
            lavsamplestats = lavsamplestats, lavdata = lavdata,
            lavcache = lavcache, group.weight = group.weight,
            lavoptions = lavoptions,
            augmented = augmented, inverted = inverted, use.ginv = use.ginv)
    } else if(information == "expected") {
        E <- lav_model_information_expected(lavmodel = lavmodel,
            lavsamplestats = lavsamplestats, lavdata = lavdata,
            lavcache = lavcache, lavoptions = lavoptions, extra = extra, 
            augmented = augmented, inverted = inverted, use.ginv = use.ginv)
    } else if(information == "first.order") {
        E <- lav_model_information_firstorder(lavmodel = lavmodel,
            lavsamplestats = lavsamplestats, lavdata = lavdata,
            lavcache = lavcache, lavoptions = lavoptions, #extra = extra,
            check.pd = FALSE,
            augmented = augmented, inverted = inverted, use.ginv = use.ginv)
    }

    # information, augmented information, or inverted information
    E
}

# fisher/expected information
#
# information = Delta' H Delta, where H is the unit information of
# the saturated model (evaluated either at the structured or unstructured
# estimates)
lav_model_information_expected <- function(lavmodel       = NULL,
                                           lavsamplestats = NULL,
                                           lavdata        = NULL,
                                           lavoptions     = NULL,
                                           Delta          = NULL,
                                           lavcache       = NULL,
                                           extra          = FALSE,
                                           augmented      = FALSE,
                                           inverted       = FALSE,
                                           use.ginv       = FALSE) {

    estimator <- lavmodel@estimator
    # structured of unstructured? (since 0.5-23)
    if(!is.null(lavoptions) &&
       !is.null(lavoptions$h1.information) &&
       lavoptions$h1.information == "unstructured") {
        structured <- FALSE
    } else {
        structured <- TRUE
    }

    if(inverted) {
        augmented <- TRUE
    }

    # compute DELTA
    if(is.null(Delta)) {
        Delta <- computeDelta(lavmodel = lavmodel)
    }

    # compute/get WLS.V == h1 information (A1)
    # if DWLS or ULS, this is the diagonal only! (since 0.5-17)
    WLS.V <- lav_model_wls_v(lavmodel       = lavmodel,
                             lavsamplestats = lavsamplestats,
                             structured     = structured,
                             lavdata        = lavdata)

    # compute Information per group
    Info.group  <- vector("list", length=lavsamplestats@ngroups)
    for(g in 1:lavsamplestats@ngroups) {
        # note LISREL documentation suggest (Ng - 1) instead of Ng...
        fg <- lavsamplestats@nobs[[g]]/lavsamplestats@ntotal
        # compute information for this group
        if(estimator %in% c("DWLS", "ULS")) {
            # diagonal weight matrix
            Delta2 <- sqrt(WLS.V[[g]]) * Delta[[g]]
            Info.group[[g]] <- fg * crossprod(Delta2)
        } else {
            # full weight matrix
            # Info.group[[g]] <- 
                # fg * (t(Delta[[g]]) %*% WLS.V[[g]] %*% Delta[[g]])
            Info.group[[g]] <- 
                fg * ( crossprod(Delta[[g]], WLS.V[[g]]) %*% Delta[[g]] )
        }
    }

    # assemble over groups
    Information <- Info.group[[1]]
    if(lavsamplestats@ngroups > 1) {
        for(g in 2:lavsamplestats@ngroups) {
            Information <- Information + Info.group[[g]]
        }
    }

    # augmented information?
    if(augmented) {
        Information <- 
            lav_model_information_augment_invert(lavmodel    = lavmodel,
                                                 information = Information,
                                                 inverted    = inverted,
                                                 use.ginv    = use.ginv)
    }

    if(extra) {
        attr(Information, "Delta") <- Delta
        attr(Information, "WLS.V") <- WLS.V # unweighted
    }

    # possibly augmented/inverted
    Information
}

# only for Mplus MLM
lav_model_information_expected_MLM <- function(lavmodel       = NULL, 
                                               lavsamplestats = NULL, 
                                               Delta          = NULL,
                                               extra          = FALSE,
                                               augmented      = FALSE,
                                               inverted       = FALSE,
                                               use.ginv       = FALSE) {

    if(inverted) {
        augmented <- TRUE
    }

    if(is.null(Delta)) {
        Delta = computeDelta(lavmodel = lavmodel)
    }

    # compute WLS.V 
    WLS.V <- vector("list", length=lavsamplestats@ngroups)
    if(lavmodel@group.w.free) {
        GW <- unlist(computeGW(lavmodel = lavmodel))
    }
    for(g in 1:lavsamplestats@ngroups) {
        WLS.V[[g]] <- lav_mvnorm_h1_information_expected(
                          sample.cov     = lavsamplestats@cov[[g]],
                          sample.cov.inv = lavsamplestats@icov[[g]])
        # the same as GLS... (except for the N/N-1 scaling)
        if(lavmodel@group.w.free) {
            # unweight!!
            a <- exp(GW[g]) / lavsamplestats@nobs[[g]]
            # a <- exp(GW[g]) * lavsamplestats@ntotal / lavsamplestats@nobs[[g]]
            WLS.V[[g]] <- lav_matrix_bdiag( matrix(a,1,1), WLS.V[[g]])
        }
    }

    # compute Information per group
    Info.group  <- vector("list", length=lavsamplestats@ngroups)
    for(g in 1:lavsamplestats@ngroups) {
        fg <- lavsamplestats@nobs[[g]]/lavsamplestats@ntotal
        # compute information for this group
        Info.group[[g]] <- fg * (t(Delta[[g]]) %*% WLS.V[[g]] %*% Delta[[g]])
    }

    # assemble over groups
    Information <- Info.group[[1]]
    if(lavsamplestats@ngroups > 1) {
        for(g in 2:lavsamplestats@ngroups) {
            Information <- Information + Info.group[[g]]
        }
    }

    # augmented information?
    if(augmented) {
        Information <-
            lav_model_information_augment_invert(lavmodel    = lavmodel,
                                                 information = Information,
                                                 inverted    = inverted,
                                                 use.ginv    = use.ginv)
    }

    if(extra) {
        attr(Information, "Delta") <- Delta
        attr(Information, "WLS.V") <- WLS.V # unweighted
    }

    Information
}

lav_model_information_observed <- function(lavmodel       = NULL,
                                           lavsamplestats = NULL,
                                           lavdata        = NULL,
                                           lavcache       = NULL,
                                           lavoptions     = NULL,
                                           group.weight   = TRUE,
                                           augmented      = FALSE,
                                           inverted       = FALSE,
                                           use.ginv       = FALSE) {
    estimator <- lavmodel@estimator

    if(inverted) {
        augmented <- TRUE
    }

    # observed.information:
    #     - "hessian": second derivative of objective function
    #     - "h1": observed information matrix of saturated (h1) model, 
    #             pre- and post-multiplied by the jacobian of the model
    #             parameters (Delta), usually evaluated at the structured
    #             sample statistics
    if(!is.null(lavoptions) &&
       !is.null(lavoptions$observed.information) &&
       lavoptions$observed.information == "h1") {
        observed.information <- "h1"
    } else {
        observed.information <- "hessian"
    }
    # structured?
    if(!is.null(lavoptions) &&
       !is.null(lavoptions$h1.information) &&
       lavoptions$h1.information == "unstructured") {
        structured <- FALSE
    } else {
        structured <- TRUE
    }
 

    if(observed.information == "hessian") {
        Hessian <- lav_model_hessian(lavmodel       = lavmodel,
                                     lavsamplestats = lavsamplestats,
                                     lavdata        = lavdata,
                                     lavoptions     = lavoptions,
                                     lavcache       = lavcache,
                                     group.weight   = group.weight)

        # NOTE! What is the relationship between the Hessian of the objective
        # function, and the `information' matrix (unit or total)
      
        # 1. in lavaan, we ALWAYS minimize, so the Hessian is already pos def
        # 2. currently, all estimators give unit information, except MML and PML
        Information <- Hessian

        # divide by 'N' for MML and PML
        if(estimator == "PML" || estimator == "MML") {
            Information <- Information / lavsamplestats@ntotal
        }


    # using 'observed h1 information'
    } else {

        # compute DELTA
        Delta <- computeDelta(lavmodel = lavmodel)
        # compute observed information h1
        if(lavmodel@estimator == "GLS"  || lavmodel@estimator == "WLS") {
            WLS.V <- lavsamplestats@WLS.V
        } else if(lavmodel@estimator == "DWLS" || lavmodel@estimator == "ULS") {
            # diagonal only!!
            WLS.V <- lavsamplestats@WLS.VD
        } else if(lavmodel@estimator == "ML") {
            WLS.V <- vector("list", length=lavsamplestats@ngroups)
            # four options: 
            #     - complete data, structured (default)
            #     - complete data, unstructured
            #     - incomplete data, structured (default)
            #     - incomplete data, unstructured
            if(structured) {
                SIGMA <- computeSigmaHat(lavmodel = lavmodel)
                if(lavmodel@meanstructure) {
                    MU <- computeMuHat(lavmodel = lavmodel)
                } else {
                    MU <- lavsamplestats@mean
                }
            } else {
                SIGMA <- lavsamplestats@cov
                MU    <- lavsamplestats@mean
            }

            # - if missing = two.stage, MU/SIGMA can be EM estimates
            #     if unstructured, or model-implied moments if structured
            for(g in 1:lavsamplestats@ngroups) {
                WLS.V[[g]] <- 
                    lav_mvnorm_information_observed_samplestats(
                        sample.mean = lavsamplestats@mean[[g]],
                        sample.cov  = lavsamplestats@cov[[g]],
                        Mu          = MU[[g]], 
                        Sigma       = SIGMA[[g]],
                        meanstructure = lavmodel@meanstructure)
            }
        } else {
            stop("lavaan ERROR: observed.information = ", 
                 dQuote(observed.information), " not supported for estimator ",
                 dQuote(lavmodel@estimator) )
        }

        # compute Information per group
        Info.group  <- vector("list", length=lavsamplestats@ngroups)
        for(g in 1:lavsamplestats@ngroups) {
            fg <- lavsamplestats@nobs[[g]]/lavsamplestats@ntotal
            # compute information for this group
            if(estimator %in% c("DWLS", "ULS")) {
                # diagonal weight matrix
                Delta2 <- sqrt(WLS.V[[g]]) * Delta[[g]]
                Info.group[[g]] <- fg * crossprod(Delta2)
            } else {
                # full weight matrix
                Info.group[[g]] <-
                    fg * ( crossprod(Delta[[g]], WLS.V[[g]]) %*% Delta[[g]] )
            }
        }

        # assemble over groups
        Information <- Info.group[[1]]
        if(lavsamplestats@ngroups > 1) {
            for(g in 2:lavsamplestats@ngroups) {
                Information <- Information + Info.group[[g]]
            }
        }
    }

    # augmented information?
    if(augmented) {
        Information <-
            lav_model_information_augment_invert(lavmodel    = lavmodel,
                                                 information = Information,
                                                 inverted    = inverted,
                                                 use.ginv    = use.ginv)
    }

    # for two.stage + observed.hession = "h1"
    if(observed.information != "hessian") {
        attr(Information, "Delta") <- Delta
        attr(Information, "WLS.V") <- WLS.V
    }

    Information
}

# outer product of the case-wise scores (gradients)
lav_model_information_firstorder <- function(lavmodel       = NULL,
                                             lavsamplestats = NULL,
                                             lavdata        = NULL,
                                             lavcache       = NULL,
                                             lavoptions     = NULL,
                                             check.pd       = FALSE,
                                             extra          = FALSE,
                                             augmented      = FALSE,
                                             inverted       = FALSE,
                                             use.ginv       = FALSE) {
    estimator <- lavmodel@estimator
    if(!is.null(lavoptions) &&
       !is.null(lavoptions$h1.information) &&
       lavoptions$h1.information == "unstructured") {
        structured <- FALSE
    } else {
        structured <- TRUE
    }

    if(inverted) {
        augmented <- TRUE
    }

    B0.group <- vector("list", lavsamplestats@ngroups)

    if(estimator == "PML") {
        Delta <- computeDelta(lavmodel = lavmodel)
        Sigma.hat <- computeSigmaHat(lavmodel = lavmodel)
        TH <- computeTH(lavmodel = lavmodel)
        if(lavmodel@conditional.x) {
            PI <- computePI(lavmodel = lavmodel)
        } else {
            PI <- vector("list", length = lavsamplestats@ngroups)
        }
    } else {
        Sigma.hat <- computeSigmaHat(lavmodel = lavmodel)
        Mu.hat <- computeMuHat(lavmodel = lavmodel)
        Delta <- computeDelta(lavmodel = lavmodel)
    }

    for(g in 1:lavsamplestats@ngroups) {
        if(estimator == "PML") {
            # slow approach: compute outer product of case-wise scores
            SC <- pml_deriv1(Sigma.hat  = Sigma.hat[[g]],
                             TH         = TH[[g]],
                             th.idx     = lavmodel@th.idx[[g]],
                             num.idx    = lavmodel@num.idx[[g]],
                             X          = lavdata@X[[g]],
                             eXo        = lavdata@eXo[[g]],
                             PI         = PI[[g]],
                             lavcache   = lavcache[[g]],
                             missing    = lavdata@missing,
                             scores     = TRUE,
                             negative   = FALSE)

            # chain rule
            group.SC <- SC %*% Delta[[g]]

            # outer product
            B0.group[[g]] <- crossprod(group.SC)

        } else if(estimator == "ML") {
            if(lavsamplestats@missing.flag) {
                B1 <- 
                  lav_mvnorm_missing_information_firstorder(Y = lavdata@X[[g]],
                      Mp = lavdata@Mp[[g]], wt = lavdata@weights[[g]],
                      Mu = Mu.hat[[g]], Sigma = Sigma.hat[[g]])
            } else {
                if(lavmodel@meanstructure) {
                    MEAN <- Mu.hat[[g]]
                } else {
                    # NOTE: the information matrix will be the same (minus
                    # the meanstructure block), but once INVERTED, the
                    # standard errors will be (slightly) smaller!!!
                    # This is only visibile when estimator = "MLF"
                    # (or information = "first.order")
                    MEAN <- lavsamplestats@mean[[g]] # saturated
                }

                if(structured) {
                    B1 <- lav_mvnorm_information_firstorder(
                              Y = lavdata@X[[g]],
                              Mu = MEAN, Sigma = Sigma.hat[[g]],
                              wt = lavdata@weights[[g]],
                              meanstructure = lavmodel@meanstructure)
                } else {
                    B1 <- lav_mvnorm_h1_information_firstorder(
                              Y = lavdata@X[[g]],
                              sample.cov.inv = lavsamplestats@icov[[g]],
                              Gamma = lavsamplestats@NACOV[[g]],
                              wt = lavdata@weights[[g]],
                              meanstructure = lavmodel@meanstructure)
                }
            }
            B0.group[[g]] <- t(Delta[[g]]) %*% B1 %*% Delta[[g]]
        } else {
            stop("lavaan ERROR: information = \"first.order\" not available for estimator ", sQuote(estimator))
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

    Information <- B0

    # NOTE: for MML and PML, we get 'total' information (instead of unit)
    # divide by 'N' for MML and PML
    if(estimator == "PML" || estimator == "MML") {
        Information <- Information / lavsamplestats@ntotal
        for(g in 1:lavsamplestats@ngroups) {
            B0.group[[g]] <- B0.group[[g]] / lavsamplestats@ntotal
        }
    }

    # augmented information?
    if(augmented) {
        Information <-
            lav_model_information_augment_invert(lavmodel    = lavmodel,
                                                 information = Information,
                                                 check.pd    = check.pd,
                                                 inverted    = inverted,
                                                 use.ginv    = use.ginv)
    }

    if(extra) {
        attr(Information, "B0.group") <- B0.group
    }

    Information
}


# create augmented information matrix (if needed), and take the inverse
# (if inverted = TRUE), returning only the [1:npar, 1:npar] elements
lav_model_information_augment_invert <- function(lavmodel    = NULL,
                                                 information = NULL,
                                                 inverted    = FALSE,
                                                 check.pd    = FALSE,
                                                 use.ginv    = FALSE) {

    npar <- nrow(information)
    is.augmented <- FALSE

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
            is.augmented <- TRUE   
            H0 <- matrix(0,nrow(H),nrow(H))
            H10 <- matrix(0, ncol(information), nrow(H))
            DL <- 2*diag(lambda, nrow(H), nrow(H))
            # FIXME: better include inactive + slacks??
            E3 <- rbind( cbind(     information,  H10, t(H)),
                         cbind(          t(H10),   DL,  H0),
                         cbind(               H,   H0,  H0)  )
            information <- E3   
        }
    }

    if(check.pd) {
        eigvals <- eigen(information, symmetric = TRUE, 
                         only.values = TRUE)$values
        if(any(eigvals < -1 * .Machine$double.eps^(3/4))) {
            warning("lavaan WARNING: matrix based on first order outer product of the derivatives is not positive definite; the model may not be identified")
        }
    }

    if(inverted) {
        if(is.augmented) {
            # note: default tol in MASS::ginv is sqrt(.Machine$double.eps)
            #       which seems a bit too conservative
            #       from 0.5-20, we changed this to .Machine$double.eps^(3/4)
            information <- 
                try( MASS::ginv(information, 
                                tol = .Machine$double.eps^(3/4))[1:npar, 
                                                                 1:npar, 
                                                                 drop = FALSE],
                     silent = TRUE )
        } else {
            if(use.ginv) {
                information <- try( MASS::ginv(information,
                                               tol = .Machine$double.eps^(3/4)),
                                    silent = TRUE )
            } else {
                information <- try( solve(information), silent = TRUE )
            }
        }
    }

    # augmented/inverted information
    information
}
