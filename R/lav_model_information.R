# here, we compute various versions of the `information' matrix
# NOTE:
# 1) we ALWAYS compute the UNIT information (not the total information)
# 
# 2) by default, we ignore the constraints (we deal with this when we
#    take the inverse later on)

lav_model_information <- function(lavmodel       = NULL,
                                  lavsamplestats = NULL,
                                  lavdata        = NULL,
                                  estimator      = "ML",
                                  Delta          = NULL,
                                  lavcache       = NULL,
                                  information    = "observed",
                                  extra          = FALSE,
                                  augmented      = FALSE,
                                  inverted       = FALSE) {

    # compute information matrix
    if(information == "observed") {
        if(lavsamplestats@missing.flag) {
            group.weight <- FALSE
        } else {
            group.weight <- TRUE
        }
        E <- lav_model_information_observed(lavmodel = lavmodel,
            lavsamplestats = lavsamplestats, lavdata = lavdata,
            lavcache = lavcache, estimator = estimator,
            group.weight = group.weight, 
            augmented = augmented, inverted = inverted)
    } else {
        E <- lav_model_information_expected(lavmodel = lavmodel,
            lavsamplestats = lavsamplestats, lavdata = lavdata,
            lavcache = lavcache, estimator = estimator,
            augmented = augmented, inverted = inverted)
    }

    # information, augmented information, or inverted information
    E
}

# fisher/expected information
lav_model_information_expected <- function(lavmodel       = NULL,
                                           lavsamplestats = NULL,
                                           lavdata        = NULL,
                                           estimator      = "ML",
                                           Delta          = NULL,
                                           lavcache       = NULL,
                                           extra          = FALSE,
                                           augmented      = FALSE,
                                           inverted       = FALSE) {

    if(inverted) {
        augmented <- TRUE
    }

    # compute DELTA
    if(is.null(Delta)) {
        Delta <- computeDelta(lavmodel = lavmodel)
    }

    # compute/get WLS.V
    # if DWLS or ULS, this is the diagonal only! (since 0.5-17)
    WLS.V <- lav_model_wls_v(lavmodel       = lavmodel, 
                             lavsamplestats = lavsamplestats,
                             estimator      = estimator,
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
                                                 inverted    = inverted)
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
                                                     inverted       = FALSE) {

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
        WLS.V[[g]] <- compute.A1.sample(lavsamplestats=lavsamplestats, group=g,
                                        meanstructure=TRUE, 
                                        information="expected")
        # the same as GLS... (except for the N/N-1 scaling)
        if(lavmodel@group.w.free) {
            # unweight!!
            a <- exp(GW[g]) / lavsamplestats@nobs[[g]]
            # a <- exp(GW[g]) * lavsamplestats@ntotal / lavsamplestats@nobs[[g]]
            WLS.V[[g]] <- bdiag( matrix(a,1,1), WLS.V[[g]])
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
                                                 inverted    = inverted)
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
                                                 estimator      = "ML",
                                                 lavcache       = NULL,
                                                 group.weight   = TRUE,
                                                 augmented      = FALSE,
                                                 inverted       = FALSE) {

    if(inverted) {
        augmented <- TRUE
    }

    Hessian <- lav_model_hessian(lavmodel       = lavmodel,
                                 lavsamplestats = lavsamplestats,
                                 lavdata        = lavdata,
                                 lavcache       = lavcache,
                                 estimator      = estimator,
                                 group.weight   = group.weight)

    # NOTE! What is the relationship between the Hessian of the objective
    # function, and the `information' matrix (unit or total)
    
    # 1. in lavaan, we ALWAYS minimize, so the Hessian is already pos def
    # 2. currently, all estimators give unit information, except MML and PML
    Information <- Hessian

    # TODO: divide by 'N' for MML and PML?

    # augmented information?
    if(augmented) {
        Information <-
            lav_model_information_augment_invert(lavmodel    = lavmodel,
                                                 information = Information,
                                                 inverted    = inverted)
    }

    Information
}

# outer product of the case-wise scores (gradients)
lav_model_information_firstorder <- function(lavmodel       = NULL,
                                                   lavsamplestats = NULL,
                                                   lavdata        = NULL,
                                                   estimator      = "ML",
                                                   lavcache       = NULL,
                                                   extra          = FALSE,
                                                   check.pd       = FALSE,
                                                   augmented      = FALSE,
                                                   inverted       = FALSE) {
    if(inverted) {
        augmented <- TRUE
    }

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

        } else {
            B1 <- compute.Bbeta(Sigma.hat=Sigma.hat[[g]],
                                Mu.hat=Mu.hat[[g]],
                                lavsamplestats=lavsamplestats,
                                lavdata=lavdata, group=g)

            B0.group[[g]] <- t(Delta[[g]]) %*% B1 %*% Delta[[g]]
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

    # NOTE: for MML and PML, we get 'total' information (instead of unit)?

    # augmented information?
    if(augmented) {
        Information <-
            lav_model_information_augment_invert(lavmodel    = lavmodel,
                                                 information = Information,
                                                 check.pd    = check.pd,
                                                 inverted    = inverted)
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
                                                 check.pd    = FALSE) {

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
            information <- 
                try( MASS::ginv(information)[1:npar, 1:npar, drop = FALSE],
                     silent = TRUE )
        } else {
            information <- try( solve(information), silent = TRUE )
        }
    }

    # augmented/inverted information
    information
}
