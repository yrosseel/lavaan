# constructor for the 'lavSampleStats' class
#
# initial version: YR 25/03/2009
# major revision: YR 5/11/2011: separate data.obs and sample statistics

lavSampleStatsFromData <- function(Data          = NULL,
                                   DataX         = NULL,
                                   DataOvnames   = NULL,
                                   DataOv        = NULL,
                                   missing       = "listwise",
                                   rescale       = FALSE,
                                   missing.h1    = TRUE,
                                   estimator     = "ML",
                                   mimic         = "lavaan",
                                   meanstructure = FALSE,
                                   WLS.V         = NULL,
                                   NACOV         = NULL,
                                   ridge         = 1e-5,
                                   debug         = FALSE,
                                   verbose       = FALSE) {

    # ridge default
    ridge.eps <- 0.0

    # get X and Mp
    if(!is.null(Data)) {
        X <- Data@X; Mp <- Data@Mp
        ngroups <- Data@ngroups
        nobs <- Data@nobs
        ov.names <- Data@ov.names
        DataOv <- Data@ov
    } else if(!is.null(DataX)) {
        stopifnot(is.list(DataX), is.matrix(DataX[[1L]]))
        X <- DataX
        ngroups <- length(X)
        Mp <- vector("list", length=ngroups)
        nobs <- vector("list", length=ngroups)
        for(g in 1:ngroups) {
            if(missing != "listwise") {
                Mp[[g]] <- getMissingPatterns(X[[g]])
            }
            nobs[[g]] <- nrow(X[[g]])
        }
        ov.names <- DataOvnames
    } else {
        stop("both Data and DataX argument are NULL")
    }

    # sample statistics per group
    cov         <- vector("list", length=ngroups)
    var         <- vector("list", length=ngroups)
    mean        <- vector("list", length=ngroups)
    th          <- vector("list", length=ngroups)
    th.nox      <- vector("list", length=ngroups)
    th.idx      <- vector("list", length=ngroups)
    th.names    <- vector("list", length=ngroups)
    slopes      <- vector("list", length=ngroups)
    cov.x       <- vector("list", length=ngroups)
    bifreq      <- vector("list", length=ngroups)
    # extra sample statistics per group
    icov          <- vector("list", length=ngroups)
    cov.log.det   <- vector("list", length=ngroups)
    WLS.obs       <- vector("list", length=ngroups)
    missing.      <- vector("list", length=ngroups)
    missing.h1.   <- vector("list", length=ngroups)
    missing.flag. <- FALSE

    if(is.null(WLS.V)) {
        WLS.V      <- vector("list", length=ngroups)
        WLS.V.user <- FALSE   
    } else {
        if(!is.list(WLS.V)) {
            if(ngroups == 1L) {
                WLS.V <- list(WLS.V)
            } else {
                stop("lavaan ERROR: WLS.V argument should be a list of length ",
                     ngroups)
            }
        } else {
            if(length(WLS.V) != ngroups)
                stop("lavaan ERROR: WLS.V assumes ", length(WLS.V),
                     " groups; data contains ", ngroups, " groups")
        }
        WLS.V.user <- TRUE
        # FIXME: check dimension of WLS.V!!
    }

    if(is.null(NACOV)) {
        NACOV      <- vector("list", length=ngroups)
        NACOV.user <- FALSE
    } else {
        if(!is.list(NACOV)) {
            if(ngroups == 1L) {
                NACOV <- list(NACOV)
            } else {
                stop("lavaan ERROR: NACOV argument should be a list of length ",
                     ngroups)
            }
        } else {
            if(length(NACOV) != ngroups)
                stop("lavaan ERROR: NACOV assumes ", length(NACOV),
                     " groups; data contains ", ngroups, " groups")
        }
        NACOV.user <- TRUE
        # FIXME: check dimension of NACOV!!
    }

    for(g in 1:ngroups) {

        # check if we have categorical data in this group
        categorical <- FALSE
        ov.types  <- DataOv$type[ match(ov.names[[g]], DataOv$name) ]
        ov.levels <- DataOv$nlev[ match(ov.names[[g]], DataOv$name) ]
        CAT <- list()
        if(!is.null(Data) && "ordered" %in% ov.types) {
            categorical <- TRUE
            if(estimator %in% c("PML","ML")) {
                WLS.W <- FALSE
            } else {
                WLS.W <- TRUE
            }
            if(verbose) {
                cat("Estimating sample thresholds and correlations ... ")
            }
            CAT <- muthen1984(Data=X[[g]], 
                              ov.names=ov.names[[g]], 
                              ov.types=ov.types,
                              ov.levels=ov.levels,
                              ov.names.x=Data@ov.names.x[[g]],
                              eXo=Data@eXo[[g]], ## FIXME, will not work with bootstrap
                              group = g, # for error messages only
                              WLS.W = WLS.W,
                              verbose=debug)
            if(verbose) cat("done\n")
            # if (and only if) all variables are ordinal, store pairwise
            # tables
            #if(all(ov.types == "ordered") && estimator == "PML" &&
            #   length(ov.types) > 1L) {
            #    # pairwise tables, as a long vector
            #    PW <- pairwiseTables(data=X[[g]], no.x=ncol(X[[g]]))$pairTables
            #    # FIXME: handle zero cells here???
            #    bifreq[[g]] <- as.numeric(unlist(PW)) 
            #    #bifreq[[g]] <- PW
            #}
        }

        # fill in the other slots
        if(categorical) {
            var[[g]]  <- CAT$VAR
            cov[[g]]  <- unname(CAT$COV)
            mean[[g]] <- apply(X[[g]], 2, mean, na.rm=TRUE)
            notnum.idx <- which(ov.types != "numeric")
            mean[[g]][notnum.idx] <- 0.0

            # th also contains the means of numeric variables
            th[[g]] <- unlist(CAT$TH)
            th.nox[[g]] <- unlist(CAT$TH.NOX)
            th.idx[[g]] <- unlist(CAT$TH.IDX)
            th.names[[g]] <- unlist(CAT$TH.NAMES)

            slopes[[g]] <- CAT$SLOPES
            if(!is.null(Data@eXo[[g]])) {
                cov.x[[g]] <- cov(Data@eXo[[g]], use="pairwise")
                if(rescale) {
                    # we 'transform' the sample cov (divided by n-1) 
                    # to a sample cov divided by 'n'
                    cov.x[[g]] <- (nobs[[g]]-1)/nobs[[g]] * cov.x[[g]]
                }               
            }
        } else {
            cov[[g]]  <-   cov(X[[g]], use="pairwise") # must be pairwise
            var[[g]]  <-   diag(cov[[g]])
            # rescale cov by (N-1)/N? (only COV!)
            if(rescale) {
                # we 'transform' the sample cov (divided by n-1) 
                # to a sample cov divided by 'n'
                cov[[g]] <- (nobs[[g]]-1)/nobs[[g]] * cov[[g]]
            }
            mean[[g]] <- apply(X[[g]], 2, mean, na.rm=TRUE)
        
            # icov and cov.log.det (but not if missing)
            if(is.null(Mp[[g]])) {
                tmp <- try(inv.chol(cov[[g]], logdet=TRUE), silent=TRUE)
                if(inherits(tmp, "try-error")) {
                    if(ngroups > 1) {
                        warning("lavaan WARNING sample covariance can not be inverted in group: ", g)
                    } else {
                        warning("lavaan WARNING: sample covariance can not be inverted")
                    }
                    # ok, try ridging for exogenous x only
                    ## FIXME -- only x (but all for now)
                    ridge.eps <- ridge
                    diag(cov[[g]]) <- diag(cov[[g]]) + ridge.eps
                    tmp <- try(inv.chol(cov[[g]], logdet=TRUE), silent=TRUE)
                    if(inherits(tmp, "try-error")) {
                        # emergency values
                        icov[[g]] <- MASS:::ginv(cov[[g]])
                        cov.log.det[[g]] <- log(.Machine$double.eps)
                    } else {
                        cov.log.det[[g]] <- attr(tmp, "logdet")
                        attr(tmp, "logdet") <- NULL
                        icov[[g]]        <- tmp
                    }
                } else {
                    cov.log.det[[g]] <- attr(tmp, "logdet")
                    attr(tmp, "logdet") <- NULL
                    icov[[g]]        <- tmp
                }
            }
        }

        # WLS.obs
        if(categorical) {
            # order of elements is important here:
            # 1. thresholds + (negative) means (interleaved)
            # 2. slopes (if any)
            # 3. variances (if any)
            # 4. covariance matrix (no diagonal!)
            TH <- th[[g]]
            TH[ov.types == "numeric"] <- -1*TH[ov.types == "numeric"]
            WLS.obs[[g]] <- c(TH,
                              vec(CAT$SLOPES), # FIXME
                              unlist(CAT$VAR[ov.types == "numeric"]),
                              vech(CAT$COV, diagonal=FALSE))
        } else if(!categorical && meanstructure) {
            WLS.obs[[g]] <- c(mean[[g]], vech(cov[[g]]))
        } else {
            WLS.obs[[g]] <- vech(cov[[g]])
        }

        # if missing = "fiml", sample statistics per pattern
        if(!is.null(Mp[[g]])) {
            missing.flag. <- TRUE
            missing.[[g]] <- 
                getMissingPatternStats(X  = X[[g]],
                                       Mp = Mp[[g]])

            #cat("missing.h1 = "); print(missing.h1); cat("\n")
            if(missing.h1) {
                # estimate moments unrestricted model
                out <- estimate.moments.EM(X=X[[g]], M=missing.[[g]],
                                           verbose=verbose)
                missing.h1.[[g]]$sigma <- out$sigma
                missing.h1.[[g]]$mu    <- out$mu
                missing.h1.[[g]]$h1    <- out$fx
            }
        }

        # NACOV (=GAMMA)
        if(!NACOV.user) {
            if(estimator == "ML") {
                NACOV[[g]] <- compute.Gamma(X[[g]], meanstructure=meanstructure)
            } else if(estimator %in% c("WLS","DWLS","ULS")) {
                if(!categorical) {
                    # sample size large enough?
                    nvar <- ncol(X[[g]]); pstar <- nvar*(nvar+1)/2
                    if(meanstructure) pstar <- pstar + nvar
                    if(nrow(X[[g]]) < pstar) {
                        if(ngroups > 1L) {
                            txt <- cat(" in group: ", g, "\n", sep="")
                        } else {
                            txt <- "\n"
                        }
                        warning("lavaan WARNING: number of observations (", 
                                nrow(X[[g]]), ") too small to compute Gamma", 
                                txt)
                    }
                    NACOV[[g]] <- compute.Gamma(X[[g]], 
                                                meanstructure=meanstructure,
                                                Mplus.WLS=(mimic=="Mplus"))
                } else { # categorical case
                    NACOV[[g]]  <- CAT$WLS.W  * (nobs[[g]] - 1L)
                }
            } else if(estimator == "PML") {
                # no NACOV ... for now
            }
        }

        # WLS.V
        if(!WLS.V.user) {
            if(estimator == "GLS") {
                if(meanstructure) {
                    V11 <- icov[[g]]
                    if(mimic == "Mplus") { # is this a bug in Mplus?
                        V11 <- V11 * nobs[[g]]/(nobs[[g]]-1)
                    }
                    V22 <- 0.5 * D.pre.post(icov[[g]] %x% icov[[g]])
                    WLS.V[[g]] <- bdiag(V11,V22)
                } else {
                    WLS.V[[g]] <-
                        0.5 * D.pre.post(icov[[g]] %x% icov[[g]])
                }
            } else if(estimator == "ML") {
                # no WLS.V here, since function of model-implied moments
            } else if(estimator %in% c("WLS","DWLS","ULS")) {
                if(!categorical) {
                    if(estimator == "WLS") {
                        # Gamma should be po before we invert
                        ev <- eigen(NACOV[[g]], symmetric=FALSE, 
                                    only.values=TRUE)$values
                        if(is.complex(ev) || any(Re(ev) < 0)) {
                           stop("lavaan ERROR: Gamma (NACOV) matrix is not positive-definite")
                        }
                        WLS.V[[g]] <- inv.chol(NACOV[[g]])
                    } else if(estimator == "DWLS") {
                        dacov <- diag(NACOV[[g]])
                        WLS.V[[g]] <- diag(1/dacov, nrow=NROW(NACOV[[g]]), 
                                                    ncol=NCOL(NACOV[[g]]))
                    } else if(estimator == "ULS") {
                        WLS.V[[g]] <- diag(length(WLS.obs[[g]]))
                    }
                } else {
                    if(estimator == "WLS") {
                        WLS.V[[g]] <- inv.chol(CAT$WLS.W * nobs[[g]])
                    } else if(estimator == "DWLS") {
                        dacov <- diag(CAT$WLS.W * nobs[[g]])
                        WLS.V[[g]] <- diag(1/dacov, nrow=NROW(CAT$WLS.W),
                                                    ncol=NCOL(CAT$WLS.W))
                    } else if(estimator == "ULS") {
                        DWLS <- diag(NROW(CAT$WLS.W))
                        WLS.V[[g]] <- DWLS
                    }
                }
            } else if(estimator == "PML") {
                # no WLS.V here
            }
        }
    } # ngroups

    # construct SampleStats object
    lavSampleStats <- new("lavSampleStats",
                       CAT = CAT, # debug only
                       # sample moments
                       th           = th,
                       th.nox       = th.nox,
                       th.idx       = th.idx,
                       th.names     = th.names,
                       mean         = mean,
                       cov          = cov,
                       var          = var,
                       slopes       = slopes,
                       cov.x        = cov.x,
                       bifreq       = bifreq,
 
                       # convenience
                       nobs         = nobs,
                       ntotal       = sum(unlist(nobs)),
                       ngroups      = ngroups,

                       # extra sample statistics
                       icov         = icov,
                       cov.log.det  = cov.log.det,
                       ridge        = ridge.eps,
                       WLS.obs      = WLS.obs,
                       WLS.V        = WLS.V,                     
                       NACOV        = NACOV,

                       # missingness
                       missing.flag = missing.flag.,
                       missing      = missing.,
                       missing.h1   = missing.h1.
                      )

    lavSampleStats
}


lavSampleStatsFromMoments <- function(sample.cov    = NULL,
                                      sample.mean   = NULL,
                                      sample.nobs   = NULL,
                                      rescale       = FALSE,
                                      ov.names      = NULL,
                                      estimator     = "ML",
                                      mimic         = "lavaan",
                                      WLS.V         = NULL,
                                      NACOV         = NULL,
                                      ridge         = 1e-5,
                                      meanstructure = FALSE) {

    # ridge default
    ridge.eps <- 0.0

    # matrix -> list
    if(!is.list(sample.cov)) sample.cov  <- list(sample.cov)
    if(!is.null(sample.mean) && !is.list(sample.mean)) {
        # check if sample.mean is string (between single quotes)
        if(is.character(sample.mean)) {
            sample.mean <- char2num(sample.mean)
        }
        sample.mean <- list(sample.mean)
    }

    # number of groups
    ngroups <- length(sample.cov)


    # sample statistics per group
    cov         <- vector("list", length=ngroups)
    var         <- vector("list", length=ngroups)
    mean        <- vector("list", length=ngroups)
    th          <- vector("list", length=ngroups)
    th.nox      <- vector("list", length=ngroups)
    th.idx      <- vector("list", length=ngroups)
    th.names    <- vector("list", length=ngroups)
    slopes      <- vector("list", length=ngroups)
    cov.x       <- vector("list", length=ngroups)
    bifreq      <- vector("list", length=ngroups)
    # extra sample statistics per group
    icov          <- vector("list", length=ngroups)
    cov.log.det   <- vector("list", length=ngroups)
    WLS.obs       <- vector("list", length=ngroups)
    missing.      <- vector("list", length=ngroups)
    missing.h1.   <- vector("list", length=ngroups)
    missing.flag. <- FALSE

    if(is.null(WLS.V)) {
        WLS.V      <- vector("list", length=ngroups)
        WLS.V.user <- FALSE
    } else {
        if(!is.list(WLS.V)) {
            if(ngroups == 1L) {
                WLS.V <- list(WLS.V)
            } else {
                stop("lavaan ERROR: WLS.V argument should be a list of length ",
                     ngroups)
            }
        } else {
            if(length(WLS.V) != ngroups)
                stop("lavaan ERROR: WLS.V assumes ", length(WLS.V),
                     " groups; data contains ", ngroups, " groups")
        }
        WLS.V.user <- TRUE
        # FIXME: check dimension of WLS.V!!
    }

    if(is.null(NACOV)) {
        NACOV      <- vector("list", length=ngroups)
        NACOV.user <- FALSE
    } else {
        if(!is.list(NACOV)) {
            if(ngroups == 1L) {
                NACOV <- list(NACOV)
            } else {
                stop("lavaan ERROR: NACOV argument should be a list of length ",
                     ngroups)
            }
        } else {
            if(length(NACOV) != ngroups)
                stop("lavaan ERROR: NACOV assumes ", length(NACOV),
                     " groups; data contains ", ngroups, " groups")
        }
        NACOV.user <- TRUE
        # FIXME: check dimension of NACOV!!
    }

    nobs    <- as.list(as.integer(sample.nobs))


    for(g in 1:ngroups) {

        tmp.cov <- sample.cov[[g]]

        # make sure that the matrix is fully symmetric (NEEDED?)
        T <- t(tmp.cov)
        tmp.cov[upper.tri(tmp.cov)] <- T[upper.tri(T)]

        # check dimnames
        if(!is.null(rownames(tmp.cov))) {
            cov.names <- rownames(tmp.cov)
        } else if(!is.null(colnames(tmp.cov))) {
            cov.names <- colnames(tmp.cov)
        } else {
            stop("lavaan ERROR: please provide row/col names ",
                "for the covariance matrix!\n")
        }

        # extract only the part we need (using ov.names)
        idx <- match(ov.names[[g]], cov.names)
        if(any(is.na(idx))) {
            cat("found: ", cov.names, "\n")
            cat("expected: ", ov.names[[g]], "\n")
            stop("lavaan ERROR: rownames of covariance matrix do not match ",
                 "the model!\n", 
                 "  found: ", paste(cov.names, collapse=" "), "\n",
                 "  expected: ", paste(ov.names[[g]], collapse=" "), "\n")
        } else {
            tmp.cov <- tmp.cov[idx,idx]
        }

        # strip dimnames
        dimnames(tmp.cov) <- NULL

        if(is.null(sample.mean)) {
            # assume zero mean vector
            tmp.mean <- numeric(ncol(tmp.cov))
        } else {
            # extract only the part we need
            tmp.mean <- as.numeric(sample.mean[[g]][idx])
        }

        cov[[g]]  <- tmp.cov
        var[[g]]  <- diag(tmp.cov)
        mean[[g]] <- tmp.mean

        # rescale cov by (N-1)/N?
        if(rescale) {
            # we 'transform' the sample cov (divided by n-1) 
            # to a sample cov divided by 'n'
            cov[[g]] <- (nobs[[g]]-1)/nobs[[g]] * cov[[g]]
        }

        # icov and cov.log.det
        tmp <- try(inv.chol(cov[[g]], logdet=TRUE), silent=TRUE)
        if(inherits(tmp, "try-error")) {
            if(ngroups > 1) {
                warning("lavaan WARNING sample covariance can not be inverted in group: ", g)
            } else {
                warning("lavaan WARNING: sample covariance can not be inverted")
            }
            # ok, try ridging for exogenous x only
            ## FIXME -- only x (but all for now)
            ridge.eps <- ridge
            diag(cov[[g]]) <- diag(cov[[g]]) + ridge.eps
            tmp <- try(inv.chol(cov[[g]], logdet=TRUE), silent=TRUE)
            if(inherits(tmp, "try-error")) {
                # emergency values
                icov[[g]] <- MASS:::ginv(cov[[g]])
                cov.log.det[[g]] <- log(.Machine$double.eps)
            } else {
                cov.log.det[[g]] <- attr(tmp, "logdet")
                attr(tmp, "logdet") <- NULL
                icov[[g]]        <- tmp
            }
        } else {
            cov.log.det[[g]] <- attr(tmp, "logdet")
            attr(tmp, "logdet") <- NULL
            icov[[g]]        <- tmp
        }

        # WLS.obs
        if(meanstructure) 
            WLS.obs[[g]] <- c(mean[[g]], vech(cov[[g]]))
        else
            WLS.obs[[g]] <- vech(cov[[g]])

        # NACOV
        # NACOV (=GAMMA)
        #if(!NACOV.user) {
            # nothing to do here; only used if provided by user
        #}

        # WLS.V
        if(!WLS.V.user) {
            if(estimator == "GLS") {
                if(meanstructure) {
                    V11 <- icov[[g]]
                    if(mimic == "Mplus") { # is this a bug in Mplus?
                        V11 <- V11 * nobs[[g]]/(nobs[[g]]-1)
                    }
                    V22 <- 0.5 * D.pre.post(icov[[g]] %x% icov[[g]])
                    WLS.V[[g]] <- bdiag(V11,V22)
                } else {
                    WLS.V[[g]] <-
                        0.5 * D.pre.post(icov[[g]] %x% icov[[g]])
                }
            } else if(estimator == "ULS") {
                WLS.V[[g]] <- diag(length(WLS.obs[[g]]))
            } else if(estimator == "WLS" || estimator == "DWLS") {
                if(is.null(WLS.V[[g]]))
                    stop("lavaan ERROR: the (D)WLS estimator is only available with full data or with a user-provided WLS.V")
            }
        }

    } # ngroups

    # construct SampleStats object
    lavSampleStats <- new("lavSampleStats",

                       # sample moments
                       th           = th,
                       th.nox       = th.nox,
                       th.idx       = th.idx,
                       th.names     = th.names,
                       mean         = mean,
                       cov          = cov,
                       var          = var,
                       slopes       = slopes,
                       cov.x        = cov.x,
                       bifreq       = bifreq,

                       # convenience
                       nobs         = nobs,
                       ntotal       = sum(unlist(nobs)),
                       ngroups      = ngroups,

                       # extra sample statistics
                       icov         = icov,
                       cov.log.det  = cov.log.det,
                       ridge        = ridge.eps,
                       WLS.obs      = WLS.obs,
                       WLS.V        = WLS.V,
                       NACOV        = NACOV,

                       # missingness
                       missing.flag = missing.flag.,
                       missing      = missing.,
                       missing.h1   = missing.h1.

                      )

    lavSampleStats
}

