# constructor for the 'lavSampleStats' class
#
# initial version: YR 25/03/2009
# major revision: YR 5/11/2011: separate data.obs and sample statistics
# YR 5/01/2016: add rescov, resvar, ... if conditional.x = TRUE

lav_samplestats_from_data <- function(lavdata           = NULL,
                                      missing           = "listwise",
                                      rescale           = FALSE,
                                      missing.h1        = TRUE,
                                      estimator         = "ML",
                                      mimic             = "lavaan",
                                      meanstructure     = FALSE,
                                      conditional.x     = FALSE,
                                      fixed.x           = FALSE,
                                      group.w.free      = FALSE,
                                      WLS.V             = NULL,
                                      NACOV             = NULL,
                                      gamma.n.minus.one = FALSE,
                                      se                = "standard",
                                      information       = "expected",
                                      ridge             = 1e-5,
                                      optim.method      = "nlminb",
                                      zero.add          = c(0.5, 0.0),
                                      zero.keep.margins = TRUE,
                                      zero.cell.warn    = TRUE,
                                      debug             = FALSE,
                                      verbose           = FALSE) {

    # ridge default
    ridge.eps <- 0.0

    # check lavdata
    stopifnot(!is.null(lavdata))

    # lavdata slots (FIXME: keep lavdata@ names)
    X <- lavdata@X; Mp <- lavdata@Mp
    ngroups <- lavdata@ngroups
    nlevels <- lavdata@nlevels
    nobs <- lavdata@nobs
    ov.names <- lavdata@ov.names
    ov.names.x <- lavdata@ov.names.x
    DataOv <- lavdata@ov
    eXo <- lavdata@eXo
    WT  <- lavdata@weights

    # sample statistics per group

    # joint (y,x)
    cov         <- vector("list", length = ngroups)
    var         <- vector("list", length = ngroups)
    mean        <- vector("list", length = ngroups)
    th          <- vector("list", length = ngroups)
    th.idx      <- vector("list", length = ngroups)
    th.names    <- vector("list", length = ngroups)

    # residual (y | x)
    res.cov     <- vector("list", length = ngroups)
    res.var     <- vector("list", length = ngroups)
    res.th      <- vector("list", length = ngroups)
    res.th.nox  <- vector("list", length = ngroups)
    res.slopes  <- vector("list", length = ngroups)
    res.int     <- vector("list", length = ngroups)

    # fixed.x
    mean.x      <- vector("list", length = ngroups)
    cov.x       <- vector("list", length = ngroups)

    # binary/ordinal
    bifreq      <- vector("list", length = ngroups)

    # extra sample statistics per group
    icov             <- vector("list", length = ngroups)
    cov.log.det      <- vector("list", length = ngroups)
    res.icov         <- vector("list", length = ngroups)
    res.cov.log.det  <- vector("list", length = ngroups)
    WLS.obs          <- vector("list", length = ngroups)
    missing.         <- vector("list", length = ngroups)
    missing.h1.      <- vector("list", length = ngroups)
    missing.flag.    <- FALSE
    zero.cell.tables <- vector("list", length = ngroups)
    YLp              <- vector("list", length = ngroups)

    # group weights
    group.w       <- vector("list", length = ngroups)

    # convenience? # FIXME!
    x.idx         <- vector("list", length = ngroups)


    WLS.VD <- vector("list", length = ngroups)
    if(is.null(WLS.V)) {
        WLS.V      <- vector("list", length = ngroups)
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

        # is WLS.V full? check first
        if(is.null(dim(WLS.V[[1]]))) {
            # we will assume it is the diagonal only
            WLS.VD <- WLS.V
            WLS.V  <- lapply(WLS.VD, diag)
        } else {
            # create WLS.VD
            WLS.VD <- lapply(WLS.V, diag)
        }

        WLS.V.user <- TRUE
        # FIXME: check dimension of WLS.V!!
    }

    NACOV.compute <- TRUE
    if(is.null(NACOV)) {
        NACOV      <- vector("list", length = ngroups)
        NACOV.user <- FALSE
    } else if(is.logical(NACOV)) {
        if(!NACOV) {
            NACOV.compute <- FALSE
        } else {
            NACOV.compute <- TRUE
        }
        NACOV.user <- FALSE
        NACOV      <- vector("list", length = ngroups)
    } else {
        NACOV.compute <- FALSE
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

    # compute some sample statistics per group
    for(g in 1:ngroups) {

        # check nobs
        if(nobs[[g]] < 2L) {
            if(nobs[[g]] == 0L) {
                stop("lavaan ERROR: data contains no observations",
                     ifelse(ngroups > 1L, paste(" in group ", g, sep=""), ""))
            } else {
                stop("lavaan ERROR: data contains only a single observation",
                     ifelse(ngroups > 1L, paste(" in group ", g, sep=""), ""))
            }
        }

        # exogenous x?
        nexo <- length(ov.names.x[[g]])
        if(nexo) {
            stopifnot( nexo == NCOL(eXo[[g]]) )

            # two cases: ov.names contains 'x' variables, or not
            if(conditional.x) {
                # ov.names.x are NOT in ov.names
                x.idx[[g]] <- length(ov.names[[g]]) + seq_len(nexo)
            } else {
                if(fixed.x) {
                    # ov.names.x are a subset of ov.names
                    x.idx[[g]] <- match(ov.names.x[[g]], ov.names[[g]])
                    stopifnot( !anyNA(x.idx[[g]]) )
                } else {
                    x.idx[[g]] <- integer(0L)
                }
            }
        } else {
            x.idx[[g]] <- integer(0L)
            conditional.x <- FALSE
            fixed.x <- FALSE
        }

        # group weight
        group.w[[g]] <- nobs[[g]] / sum(unlist(nobs))

        # check if we have categorical data in this group
        categorical <- FALSE
        ov.types  <- DataOv$type[ match(ov.names[[g]], DataOv$name) ]
        ov.levels <- DataOv$nlev[ match(ov.names[[g]], DataOv$name) ]
        CAT <- list()
        if("ordered" %in% ov.types) {
            categorical <- TRUE
            if(nlevels > 1L) {
                warning("lavaan ERROR: multilevel + categorical not supported yet.")
            }
        }

        if(categorical) {
            if(estimator %in% c("ML","REML","PML","FML","MML","none","ULS")) {
                WLS.W <- FALSE
                if(estimator == "ULS" && se == "robust.sem") { #||
                        #test %in% c("satorra.bentler", "scaled.shifted",
                        #            "mean.var.adjusted"))) {
                    WLS.W <- TRUE
                }
            } else {
                WLS.W <- TRUE
            }
            if(verbose) {
                cat("Estimating sample thresholds and correlations ... ")
            }

            if(conditional.x) {
                CAT <- muthen1984(Data=X[[g]],
                                  ov.names=ov.names[[g]],
                                  ov.types=ov.types,
                                  ov.levels=ov.levels,
                                  ov.names.x=ov.names.x[[g]],
                                  eXo=eXo[[g]],
                                  group = g, # for error messages only
                                  missing = missing, # listwise or pairwise?
                                  WLS.W = WLS.W,
                                  optim.method = optim.method,
                                  zero.add = zero.add,
                                  zero.keep.margins = zero.keep.margins,
                                  zero.cell.warn = FALSE,
                                  zero.cell.tables = TRUE,
                                  verbose=debug)
            } else {
                CAT <- muthen1984(Data=X[[g]],
                                  ov.names=ov.names[[g]],
                                  ov.types=ov.types,
                                  ov.levels=ov.levels,
                                  ov.names.x=NULL,
                                  eXo=NULL,
                                  group = g, # for error messages only
                                  missing = missing, # listwise or pairwise?
                                  WLS.W = WLS.W,
                                  optim.method = optim.method,
                                  zero.add = zero.add,
                                  zero.keep.margins = zero.keep.margins,
                                  zero.cell.warn = FALSE,
                                  zero.cell.tables = TRUE,
                                  verbose=debug)
            }
            # empty cell tables
            zero.cell.tables[[g]] <- CAT$zero.cell.tables
            if(verbose) cat("done\n")
        }

        if(categorical) {

            # convenience
            th.idx[[g]] <- unlist(CAT$TH.IDX)
            th.names[[g]] <- unlist(CAT$TH.NAMES)

            if(conditional.x) {
                # residual var/cov
                res.var[[g]]  <- unlist(CAT$VAR)
                res.cov[[g]]  <- unname(CAT$COV)

                # th also contains the means of numeric variables
                res.th[[g]]     <- unlist(CAT$TH)
                res.th.nox[[g]] <- unlist(CAT$TH.NOX)

                # for convenience, we store the intercept of numeric
                # variables in res.int
                NVAR <- NCOL(res.cov[[g]])
                mean[[g]] <- res.int[[g]] <- numeric(NVAR)
                num.idx <- which(!seq_len(NVAR) %in% th.idx[[g]])
                if(length(num.idx) > 0L) {
                    NUM.idx <- which(th.idx[[g]] == 0L)
                    mean[[g]][num.idx] <- res.th.nox[[g]][NUM.idx]
                    res.int[[g]][num.idx]  <- res.th[[g]][NUM.idx]
                }

                # slopes
                res.slopes[[g]] <- CAT$SLOPES
            } else {
                # var/cov
                var[[g]]  <- unlist(CAT$VAR)
                cov[[g]]  <- unname(CAT$COV)

                # th also contains the means of numeric variables
                th[[g]] <- unlist(CAT$TH)

                # mean (numeric only)
                NVAR <- NCOL(cov[[g]])
                mean[[g]] <- numeric(NVAR)
                num.idx <- which(!seq_len(NVAR) %in% th.idx[[g]])
                if(length(num.idx) > 0L) {
                    NUM.idx <- which(th.idx[[g]] == 0L)
                    mean[[g]][num.idx] <- th[[g]][NUM.idx]
                }

            }

        } else if(nlevels > 1L) { # continuous, multilevel setting

            # level-based sample statistics
            YLp[[g]] <- lav_samplestats_cluster_patterns(Y  = X[[g]],
                                                         Lp = lavdata@Lp[[g]])

            # FIXME: needed?
            cov[[g]] <- unname( stats::cov(X[[g]], use="pairwise"))
            mean[[g]] <- unname( colMeans(X[[g]], na.rm=TRUE) )
            var[[g]] <- diag(cov[[g]])

        } else { # continuous, single-level case

            if(conditional.x) {

                # FIXME!
                # no handling of missing data yet....
                if(missing %in% c("ml", "ml.x",
                                  "two.stage", "robust.two.stage")) {
                    stop("lavaan ERROR: missing = ", missing, " + conditional.x not supported yet")
                }

                # residual covariances!

                   Y <- cbind(X[[g]], eXo[[g]])
                 COV <- unname( stats::cov(Y, use="pairwise"))
                MEAN <- unname( colMeans(Y, na.rm=TRUE) )
                var[[g]] <- diag(COV)
                cov[[g]] <- COV
                # rescale cov by (N-1)/N? (only COV!)
                if(rescale) {
                    # we 'transform' the sample cov (divided by n-1)
                    # to a sample cov divided by 'n'
                    COV <- (nobs[[g]]-1)/nobs[[g]] * COV
                }
                cov[[g]] <- COV
                mean[[g]] <- MEAN

                A <- COV[-x.idx[[g]], -x.idx[[g]], drop=FALSE]
                B <- COV[-x.idx[[g]],  x.idx[[g]], drop=FALSE]
                C <- COV[ x.idx[[g]],  x.idx[[g]], drop=FALSE]
                # FIXME: make robust against singular C!!!
                res.cov[[g]] <- A - B %*% solve(C) %*% t(B)
                res.var[[g]] <- diag( cov[[g]] )


                MY <- MEAN[-x.idx[[g]]]; MX <- MEAN[x.idx[[g]]]
                C3 <- rbind(c(1,MX),
                            cbind(MX, C + tcrossprod(MX)))
                B3 <- cbind(MY, B + tcrossprod(MY,MX))
                COEF <- unname(solve(C3, t(B3)))

                res.int[[g]]    <- COEF[1,]                  # intercepts
                res.slopes[[g]] <- t(COEF[-1,,drop = FALSE]) # slopes

            } else if(missing == "two.stage" ||
                      missing == "robust.two.stage") {
                missing.flag. <- FALSE #!!! just use sample statistics
                missing.[[g]] <-
                    lav_samplestats_missing_patterns(Y  = X[[g]],
                                                     Mp = Mp[[g]],
                                                     wt = WT[[g]])
                out <- lav_mvnorm_missing_h1_estimate_moments(Y = X[[g]],
                          wt = WT[[g]],
                          Mp = Mp[[g]], Yp = missing.[[g]], verbose = verbose)
                missing.h1.[[g]]$sigma <- out$Sigma
                missing.h1.[[g]]$mu    <- out$Mu
                missing.h1.[[g]]$h1    <- out$fx

                # here, sample statistics == EM estimates
                cov[[g]]  <- missing.h1.[[g]]$sigma
                var[[g]]  <- diag(cov[[g]])
                mean[[g]] <- missing.h1.[[g]]$mu
            } else if(missing %in% c("ml", "ml.x")) {
                missing.flag. <- TRUE
                missing.[[g]] <-
                    lav_samplestats_missing_patterns(Y  = X[[g]],
                                                     Mp = Mp[[g]],
                                                     wt = WT[[g]])

                if(missing.h1 && nlevels == 1L) {
                    # estimate moments unrestricted model
                    out <- lav_mvnorm_missing_h1_estimate_moments(Y = X[[g]],
                              wt = WT[[g]],
                              Mp = Mp[[g]], Yp = missing.[[g]], verbose = verbose)
                    missing.h1.[[g]]$sigma <- out$Sigma
                    missing.h1.[[g]]$mu    <- out$Mu
                    missing.h1.[[g]]$h1    <- out$fx
                }

                if(!is.null(WT[[g]])) {
                    # here, sample statistics == EM estimates
                    cov[[g]]  <- missing.h1.[[g]]$sigma
                    var[[g]]  <- diag(cov[[g]])
                    mean[[g]] <- missing.h1.[[g]]$mu
                } else {
                    # NEEDED? why not just EM-based?
                    cov[[g]]  <-   stats::cov(X[[g]], use = "pairwise")
                    var[[g]]  <-   diag(cov[[g]])
                    # rescale cov by (N-1)/N? (only COV!)
                    if(rescale) {
                        # we 'transform' the sample cov (divided by n-1)
                        # to a sample cov divided by 'n'
                        cov[[g]] <- (nobs[[g]]-1)/nobs[[g]] * cov[[g]]
                    }
                    mean[[g]] <- colMeans(X[[g]], na.rm=TRUE)
                }
            } else {
                # LISTWISE
                if(!is.null(WT[[g]])) {
                    out <- stats::cov.wt(X[[g]], wt = WT[[g]],
                                        method = "ML")
                    cov[[g]]  <- out$cov
                    var[[g]]  <- diag(cov[[g]])
                    mean[[g]] <- out$center
                } else {
                    cov[[g]]  <-   stats::cov(X[[g]], use = "pairwise")
                    var[[g]]  <-   diag(cov[[g]])
                    # rescale cov by (N-1)/N? (only COV!)
                    if(rescale) {
                        # we 'transform' the sample cov (divided by n-1)
                        # to a sample cov divided by 'n'
                        cov[[g]] <- (nobs[[g]]-1)/nobs[[g]] * cov[[g]]
                    }
                    mean[[g]] <- colMeans(X[[g]], na.rm=TRUE)
                }
            }


            # icov and cov.log.det (but not if missing)
            if(!missing %in% c("ml", "ml.x")) {
                out <- lav_samplestats_icov(COV = cov[[g]], ridge = ridge,
                           x.idx = x.idx[[g]],
                           ngroups = ngroups, g = g, warn = TRUE)
                icov[[g]] <- out$icov
                cov.log.det[[g]] <- out$cov.log.det

                # the same for res.cov if conditional.x = TRUE
                if(conditional.x) {
                    out <- lav_samplestats_icov(COV = res.cov[[g]], ridge=ridge,
                               x.idx = x.idx[[g]],
                               ngroups = ngroups, g = g, warn = TRUE)
                    res.icov[[g]] <- out$icov
                    res.cov.log.det[[g]] <- out$cov.log.det
                }
            }
        }


        # WLS.obs
        if(nlevels == 1L) {
            WLS.obs[[g]] <- lav_samplestats_wls_obs(mean.g = mean[[g]],
                cov.g = cov[[g]], var.g = var[[g]], th.g = th[[g]],
                th.idx.g = th.idx[[g]], res.int.g = res.int[[g]],
                res.cov.g = res.cov[[g]], res.var.g = res.var[[g]],
                res.th.g = res.th[[g]],  res.slopes.g = res.slopes[[g]],
                group.w.g = group.w[[g]],
                categorical = categorical, conditional.x = conditional.x,
                meanstructure = meanstructure, slopestructure = conditional.x,
                group.w.free = group.w.free)

        }

        # fill in the other slots
        if(!is.null(eXo[[g]])) {
            if(!is.null(WT[[g]])) {
                if(missing != "listwise") {
                    cov.x[[g]]  <- missing.h1.[[g]]$sigma[ x.idx[[g]],
                                                           x.idx[[g]],
                                                           drop = FALSE ]
                    mean.x[[g]] <- missing.h1.[[g]]$mu[  x.idx[[g]] ]
                } else {
                    out <- stats::cov.wt(eXo[[g]], wt = WT[[g]],
                                         method = "ML")
                    cov.x[[g]]  <- out$cov
                    mean.x[[g]] <- out$center
                }
            } else {
                cov.x[[g]] <- cov(eXo[[g]], use="pairwise")
                if(rescale) {
                    # we 'transform' the sample cov (divided by n-1)
                    # to a sample cov divided by 'n'
                    cov.x[[g]] <- (nobs[[g]]-1)/nobs[[g]] * cov.x[[g]]
                }
                mean.x[[g]] <- colMeans(eXo[[g]])
            }
        }

        # NACOV (=GAMMA)
        if(!NACOV.user && nlevels == 1L) {
            if(estimator == "ML" && !missing.flag. && NACOV.compute) {
                if(conditional.x) {
                     Y <- Y
                } else {
                    Y <- X[[g]]
                }

                if(length(lavdata@cluster) > 0L) {
                    cluster.idx <- lavdata@Lp[[g]]$cluster.idx[[2]]
                } else {
                    cluster.idx <- NULL
                }

                NACOV[[g]] <-
                    lav_samplestats_Gamma(Y              = Y,
                                          x.idx          = x.idx[[g]],
                                          cluster.idx    = cluster.idx,
                                          fixed.x        = fixed.x,
                                          conditional.x  = conditional.x,
                                          meanstructure  = meanstructure,
                                          slopestructure = conditional.x,
                                          gamma.n.minus.one = gamma.n.minus.one,
                                          Mplus.WLS      = FALSE)
            } else if(estimator %in% c("WLS","DWLS","ULS")) {
                if(!categorical) {
                    # sample size large enough?
                    nvar <- ncol(X[[g]])
                    #if(conditional.x && nexo > 0L) {
                    #    nvar <- nvar - nexo
                    #}
                    pstar <- nvar*(nvar+1)/2
                    if(meanstructure) pstar <- pstar + nvar
                    if(conditional.x && nexo > 0L) {
                        pstar <- pstar + (nvar * nexo)
                    }
                    if(nrow(X[[g]]) < pstar) {
                        if(ngroups > 1L) {
                            txt <- paste(" in group: ", g, "\n", sep="")
                        } else {
                            txt <- "\n"
                        }
                        warning("lavaan WARNING: number of observations (",
                                nrow(X[[g]]), ") too small to compute Gamma",
                                txt)
                    }
                    if(conditional.x) {
                        Y <- Y
                    } else {
                        Y <- X[[g]]
                    }

                    if(length(lavdata@cluster) > 0L) {
                        cluster.idx <- lavdata@Lp[[g]]$cluster.idx[[2]]
                    } else {
                        cluster.idx <- NULL
                    }
                    NACOV[[g]] <-
                        lav_samplestats_Gamma(Y             = Y,
                                              x.idx         = x.idx[[g]],
                                              cluster.idx   = cluster.idx,
                                              fixed.x       = fixed.x,
                                              conditional.x = conditional.x,
                                              meanstructure = meanstructure,
                                              slopestructure = conditional.x,
                                              gamma.n.minus.one = gamma.n.minus.one,
                                              Mplus.WLS = (mimic=="Mplus"))
                } else { # categorical case
                    NACOV[[g]]  <- CAT$WLS.W  * nobs[[g]]
                    if(gamma.n.minus.one) {
                        NACOV[[g]] <- NACOV[[g]] * nobs[[g]] / (nobs[[g]] - 1L)
                    }
                }
            } else if(estimator == "PML") {
                # no NACOV ... for now
            }
        }

        # WLS.V
        if(!WLS.V.user && nlevels == 1L) {
            if(estimator == "GLS") {
                # Note: we need the 'original' COV/MEAN/ICOV
                #        sample statistics; not the 'residual' version
                WLS.V[[g]] <- lav_samplestats_Gamma_inverse_NT(
                    ICOV           = icov[[g]],
                    COV            = cov[[g]],
                    MEAN           = mean[[g]],
                    rescale        = FALSE,
                    x.idx          = x.idx[[g]],
                    fixed.x        = fixed.x,
                    conditional.x  = conditional.x,
                    meanstructure  = meanstructure,
                    slopestructure = conditional.x)
                if(mimic == "Mplus" && !conditional.x && meanstructure) {
                    # bug in Mplus? V11 rescaled by nobs[[g]]/(nobs[[g]]-1)
                    nvar <- NCOL(cov[[g]])
                    WLS.V[[g]][1:nvar, 1:nvar] <- WLS.V[[g]][1:nvar, 1:nvar,
                                        drop = FALSE] * nobs[[g]]/(nobs[[g]]-1)
                }


            } else if(estimator == "ML") {
                # no WLS.V here, since function of model-implied moments


            } else if(estimator %in% c("WLS","DWLS","ULS")) {
                if(!categorical) {
                    if(estimator == "WLS") {
                        if(!fixed.x) {
                            # Gamma should be po before we invert
                            ev <- eigen(NACOV[[g]], # symmetric=FALSE,
                                        only.values=TRUE)$values
                            if(is.complex(ev)) {
                               stop("lavaan ERROR: Gamma (NACOV) matrix is not positive-definite")
                            }
                            if(any(Re(ev) < 0)) {
                                stop("lavaan ERROR: Gamma (NACOV) matrix is not positive-definite")
                            }
                            WLS.V[[g]] <- lav_matrix_symmetric_inverse(NACOV[[g]])
                        } else {
                            # fixed.x: we have zero cols/rows
                            # ginv does the trick, but perhaps this is overkill
                            # just removing the zero rows/cols, invert, and
                            # fill back in the zero rows/cols would do it
                            WLS.V[[g]] <- MASS::ginv(NACOV[[g]])
                        }
                    } else if(estimator == "DWLS") {
                        dacov <- diag(NACOV[[g]])
                        if(!all(is.finite(dacov))) {
                            stop("lavaan ERROR: diagonal of Gamma (NACOV) contains non finite values")
                        }
                        if(fixed.x) {
                            # structural zeroes!
                            zero.idx <- which(dacov == 0.0)
                            idacov <- 1/dacov
                            idacov[zero.idx] <- 0.0
                        } else {
                            idacov <- 1/dacov
                        }
                        WLS.V[[g]] <- diag(idacov, nrow=NROW(NACOV[[g]]),
                                                   ncol=NCOL(NACOV[[g]]))
                        WLS.VD[[g]] <- idacov
                    } else if(estimator == "ULS") {
                        #WLS.V[[g]] <- diag(length(WLS.obs[[g]]))
                        WLS.VD[[g]] <- rep(1, length(WLS.obs[[g]]))
                    }
                } else {
                    if(estimator == "WLS") {
                        WLS.V[[g]] <- inv.chol(CAT$WLS.W * nobs[[g]])
                    } else if(estimator == "DWLS") {
                        dacov <- diag(CAT$WLS.W * nobs[[g]])
                        #WLS.V[[g]] <- diag(1/dacov, nrow=NROW(CAT$WLS.W),
                        #                            ncol=NCOL(CAT$WLS.W))
                        WLS.VD[[g]] <- 1/dacov
                    } else if(estimator == "ULS") {
                        #WLS.V[[g]] <- diag(length(WLS.obs[[g]]))
                        WLS.VD[[g]] <- rep(1, length(WLS.obs[[g]]))
                    }
                }
            } else if(estimator == "PML" || estimator == "FML") {
                # no WLS.V here
            }

            # group.w.free
            if(!is.null(WLS.V[[g]]) && group.w.free) {
                # unweight!!
                #a <- group.w[[g]] * sum(unlist(nobs)) / nobs[[g]]
                # always 1!!!
                ### FIXME: this is consistent with expected information
                ###        but why not group.w[[g]] * (1 - group.w[[g]])?
                #a <- 1
                a <- group.w[[g]] * (1 - group.w[[g]]) * sum(unlist(nobs)) / nobs[[g]]
                # invert
                a <- 1/a
                WLS.V[[g]] <- lav_matrix_bdiag( matrix(a, 1, 1), WLS.V[[g]] )
            }

        }
    } # ngroups

    # remove 'CAT', unless debug -- this is to save memory
    if(!debug) {
        CAT <- list()
    }

    # construct SampleStats object
    lavSampleStats <- new("lavSampleStats",
                       # sample moments
                       th           = th,
                       th.idx       = th.idx,
                       th.names     = th.names,
                       mean         = mean,
                       cov          = cov,
                       var          = var,

                       # residual (y | x)
                       res.cov      = res.cov,
                       res.var      = res.var,
                       res.th       = res.th,
                       res.th.nox   = res.th.nox,
                       res.slopes   = res.slopes,
                       res.int      = res.int,

                       mean.x       = mean.x,
                       cov.x        = cov.x,
                       bifreq       = bifreq,
                       group.w      = group.w,

                       # convenience
                       nobs         = nobs,
                       ntotal       = sum(unlist(nobs)),
                       ngroups      = ngroups,
                       x.idx        = x.idx,

                       # extra sample statistics
                       icov            = icov,
                       cov.log.det     = cov.log.det,
                       res.icov        = res.icov,
                       res.cov.log.det = res.cov.log.det,
                       ridge           = ridge.eps,
                       WLS.obs         = WLS.obs,
                       WLS.V           = WLS.V,
                       WLS.VD          = WLS.VD,
                       NACOV           = NACOV,
                       NACOV.user      = NACOV.user,

                       # cluster/levels
                       YLp             = YLp,

                       # missingness
                       missing.flag    = missing.flag.,
                       missing         = missing.,
                       missing.h1      = missing.h1.,
                       zero.cell.tables = zero.cell.tables
                      )

    # just a SINGLE warning if we have empty cells
    if(categorical && zero.cell.warn &&
       any(sapply(zero.cell.tables, nrow) > 0L)) {
        nempty <- sum(sapply(zero.cell.tables, nrow))
        warning("lavaan WARNING: ", nempty,
                " bivariate tables have empty cells; to see them, use:\n",
                "                  lavInspect(fit, \"zero.cell.tables\")")
    }

    lavSampleStats
}


lav_samplestats_from_moments <- function(sample.cov    = NULL,
                                         sample.mean   = NULL,
                                         sample.th     = NULL,
                                         sample.nobs   = NULL,
                                         rescale       = FALSE,
                                         ov.names      = NULL, # including x
                                         ov.names.x    = NULL,
                                         estimator     = "ML",
                                         mimic         = "lavaan",
                                         WLS.V         = NULL,
                                         NACOV         = NULL,
                                         ridge         = 1e-5,
                                         group.w.free  = FALSE) {

    # no multilevel yet
    nlevels <- 1L

    # ridge default
    ridge.eps <- 0.0

    # new in 0.6-3:
    # check if sample.cov has attributes if conditional.x = TRUE
    sample.res.slopes <- attr(sample.cov, "res.slopes")
    sample.cov.x      <- attr(sample.cov, "cov.x")
    sample.mean.x     <- attr(sample.cov, "mean.x")
    if(!is.null(sample.res.slopes)) {
        conditional.x = TRUE
        # strip attributes
        attr(sample.cov, "res.slopes") <- NULL
        attr(sample.cov, "cov.x") <- NULL
        attr(sample.cov, "mean.x") <- NULL
        # make list
        if(!is.list(sample.res.slopes)) {
            sample.res.slopes <- list(sample.res.slopes)
        }
        if(!is.list(sample.cov.x)) {
            sample.cov.x <- list(sample.cov.x)
        }
        if(!is.list(sample.mean.x)) {
            sample.mean.x <- list(sample.mean.x)
        }
    } else if(!is.null(sample.cov.x)) {
        conditional.x <- FALSE
        fixed.x <- TRUE

        # strip attributes
        attr(sample.cov, "cov.x") <- NULL
        attr(sample.cov, "mean.x") <- NULL
        # make list
        if(!is.list(sample.cov.x)) {
            sample.cov.x <- list(sample.cov.x)
        }
        if(!is.list(sample.mean.x)) {
            sample.mean.x <- list(sample.mean.x)
        }
    } else {
        conditional.x <- FALSE
        fixed.x <- FALSE
    }

    # matrix -> list
    if(!is.list(sample.cov)) {
        sample.cov  <- list(sample.cov)
    }

    # number of groups
    ngroups <- length(sample.cov)

    # ov.names
    if(!is.list(ov.names)) {
        ov.names <- rep(list(ov.names), ngroups)
    }
    if(!is.list(ov.names.x)) {
        ov.names.x <- rep(list(ov.names.x), ngroups)
    }

    if(!is.null(sample.mean)) {
        meanstructure <- TRUE
        if(!is.list(sample.mean)) {
            # check if sample.mean is string (between single quotes)
            if(is.character(sample.mean)) {
                sample.mean <- char2num(sample.mean)
            }
            sample.mean <- list(unname(sample.mean))
        } else {
            sample.mean <- lapply(lapply(sample.mean, unname), unclass)
        }
    } else {
        meanstructure <- FALSE
    }

    if(!is.null(sample.th)) {
        th.idx <- attr(sample.th, "th.idx")
        attr(sample.th, "th.idx") <- NULL
        if(is.null(th.idx)) {
            stop("lavaan ERROR: sample.th should have a th.idx attribute")
        } else {
            if(is.list(th.idx)) {
                th.names <- lapply(th.idx, names)
                th.idx <- lapply(lapply(th.idx, unname), unclass)
            } else {
                th.names <- list(names(th.idx))
                th.idx   <- list(unclass(unname(th.idx)))
            }
        }
        if(is.list(sample.th)) {
            # strip names and lavaan.vector class
            sample.th <- lapply(lapply(sample.th, unname), unclass)
        } else {
            # strip names and lavaan.vector class, make list
            sample.th <- list(unclass(unname(sample.th)))
        }
    } else {
        th.idx <- vector("list", length = ngroups)
        th.names <- vector("list", length = ngroups)
    }

    # sample statistics per group
    cov         <- vector("list", length = ngroups)
    var         <- vector("list", length = ngroups)
    mean        <- vector("list", length = ngroups)
    th          <- vector("list", length = ngroups)
    #th.idx      <- vector("list", length = ngroups)
    #th.names    <- vector("list", length = ngroups)

    # residual (y | x)
    res.cov     <- vector("list", length = ngroups)
    res.var     <- vector("list", length = ngroups)
    res.slopes  <- vector("list", length = ngroups)
    res.int     <- vector("list", length = ngroups)
    res.th      <- vector("list", length = ngroups)
    res.th.nox  <- vector("list", length = ngroups)

    # fixed.x / conditional.x
    mean.x      <- vector("list", length = ngroups)
    cov.x       <- vector("list", length = ngroups)

    bifreq      <- vector("list", length = ngroups)

    # extra sample statistics per group
    icov          <- vector("list", length = ngroups)
    cov.log.det   <- vector("list", length = ngroups)
    res.icov        <- vector("list", length = ngroups)
    res.cov.log.det <- vector("list", length = ngroups)
    WLS.obs         <- vector("list", length = ngroups)
    missing.        <- vector("list", length = ngroups)
    missing.h1.     <- vector("list", length = ngroups)
    missing.flag.   <- FALSE
    zero.cell.tables <- vector("list", length = ngroups)
    YLp              <- vector("list", length = ngroups)

    # group weights
    group.w         <- vector("list", length = ngroups)
    x.idx           <- vector("list", length = ngroups)

    categorical <- FALSE
    if(!is.null(sample.th)) {
        categorical <- TRUE
    }

    WLS.VD <- vector("list", length = ngroups)
    if(is.null(WLS.V)) {
        WLS.V      <- vector("list", length = ngroups)
        WLS.V.user <- FALSE
    } else {
        if(!is.list(WLS.V)) {
            if(ngroups == 1L) {
                WLS.V <- list(unclass(WLS.V))
            } else {
                stop("lavaan ERROR: WLS.V argument should be a list of length ",
                     ngroups)
            }
        } else {
            if(length(WLS.V) != ngroups) {
                stop("lavaan ERROR: WLS.V assumes ", length(WLS.V),
                     " groups; data contains ", ngroups, " groups")
            }
            WLS.V <- lapply(WLS.V, unclass)
        }

        # is WLS.V full? check first
        if(is.null(dim(WLS.V[[1]]))) {
            # we will assume it is the diagonal only
            WLS.VD <- WLS.V
            WLS.V  <- lapply(WLS.VD, diag)
        } else {
            # create WLS.VD
            WLS.VD <- lapply(WLS.V, diag)
            # we could remove WLS.V to save space...
        }

        WLS.V.user <- TRUE
        # FIXME: check dimension of WLS.V!!
    }

    if(is.null(NACOV)) {
        NACOV      <- vector("list", length = ngroups)
        NACOV.user <- FALSE
    } else {
        if(!is.list(NACOV)) {
            if(ngroups == 1L) {
                NACOV <- list(unclass(NACOV))
            } else {
                stop("lavaan ERROR: NACOV argument should be a list of length ",
                     ngroups)
            }
        } else {
            if(length(NACOV) != ngroups) {
                stop("lavaan ERROR: NACOV assumes ", length(NACOV),
                     " groups; data contains ", ngroups, " groups")
            }
            NACOV <- lapply(NACOV, unclass)
        }
        NACOV.user <- TRUE
        # FIXME: check dimension of NACOV!!
    }

    nobs    <- as.list(as.integer(sample.nobs))


    for(g in 1:ngroups) {

        # exogenous x?
        nexo <- length(ov.names.x[[g]])
        if(nexo) {
            # two cases: ov.names contains 'x' variables, or not
            if(conditional.x) {
                # ov.names.x are NOT in ov.names
                x.idx[[g]] <- which(ov.names[[g]] %in% ov.names.x[[g]])
            } else {
                if(fixed.x) {
                    # ov.names.x are a subset of ov.names
                    x.idx[[g]] <- match(ov.names.x[[g]], ov.names[[g]])
                    stopifnot( !anyNA(x.idx[[g]]) )
                } else {
                    x.idx[[g]] <- integer(0L)
                }
            }
        } else {
            x.idx[[g]] <- integer(0L)
            conditional.x <- FALSE
            fixed.x <- FALSE
        }


        # group weight
        group.w[[g]] <- nobs[[g]] / sum(unlist(nobs))

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
        if(conditional.x) {
            idx <- match(ov.names[[g]][-x.idx[[g]]], cov.names)
        } else {
            idx <- match(ov.names[[g]], cov.names)
        }
        if(any(is.na(idx))) {
            cat("found: ", cov.names, "\n")
            cat("expected: ", ov.names[[g]], "\n")
            stop("lavaan ERROR: rownames of covariance matrix do not match ",
                 "the model!\n",
                 "  found: ", paste(cov.names, collapse=" "), "\n",
                 "  expected: ", paste(ov.names[[g]], collapse=" "), "\n")
        } else {
            tmp.cov <- tmp.cov[idx,idx,drop=FALSE]
        }

        # strip dimnames
        dimnames(tmp.cov) <- NULL

        if(is.null(sample.mean)) {
            # assume zero mean vector
            tmp.mean <- numeric(ncol(tmp.cov))
        } else {
            # extract only the part we need
            tmp.mean <- unclass(sample.mean[[g]][idx])
        }



        if(categorical) {

            # categorical + conditional.x = TRUE
            if(conditional.x) {

                th.g <- numeric( length(th.idx[[g]]) )
                ord.idx <- which(th.idx[[g]] > 0)
                num.idx <- which(th.idx[[g]] == 0)
                if(length(ord.idx) > 0L) {
                    th.g[ord.idx] <- sample.th[[g]]
                }
                if(length(num.idx) > 0L) {
                    ord.var.idx <- unique(th.idx[[g]][th.idx[[g]] > 0])
                    th.g[num.idx] <- -1 * sample.mean[[g]][ -ord.var.idx ]
                }
                res.th[[g]] <- th.g
                res.th.nox[[g]] <- sample.th[[g]]

                res.cov[[g]] <- tmp.cov
                res.var[[g]] <- diag(tmp.cov)
                res.int[[g]] <- tmp.mean

                res.slopes[[g]] <- unclass(unname(sample.res.slopes[[g]]))
                cov.x[[g]]      <- unclass(unname(sample.cov.x[[g]]))
                mean.x[[g]]     <- unclass(unname(sample.mean.x[[g]]))

                # th.idx and th.names are already ok

            # categorical + conditional.x = FALSE
            } else {

                th.g <- numeric( length(th.idx[[g]]) )
                ord.idx <- which(th.idx[[g]] > 0)
                num.idx <- which(th.idx[[g]] == 0)
                if(length(ord.idx) > 0L) {
                    th.g[ord.idx] <- sample.th[[g]]
                }
                if(length(num.idx) > 0L) {
                    ord.var.idx <- unique(th.idx[[g]][th.idx[[g]] > 0])
                    th.g[num.idx] <- -1 * sample.mean[[g]][ -ord.var.idx ]
                }
                th[[g]] <- th.g

                cov[[g]]  <- tmp.cov
                var[[g]]  <- diag(tmp.cov)
                mean[[g]] <- tmp.mean

                # fixed.x? (needed?)
                if(fixed.x) {
                    cov.x[[g]]  <- unclass(unname(sample.cov.x[[g]]))
                    mean.x[[g]] <- unclass(unname(sample.mean.x[[g]]))
                }

                # th, th.idx and th.names are already ok
            }

        # multilevel
        } else if(nlevels > 1L) {
            stop("lavaan ERROR: multilevel + sample stats not ready yet")


        # single level
        } else {

            # single-level + continuous + conditional.x = TRUE
            if(conditional.x) {
                res.cov[[g]] <- tmp.cov
                res.var[[g]] <- diag(tmp.cov)
                res.int[[g]] <- tmp.mean
                res.slopes[[g]] <- unclass(unname(sample.res.slopes[[g]]))
                cov.x[[g]]      <- unclass(unname(sample.cov.x[[g]]))
                mean.x[[g]]     <- unclass(unname(sample.mean.x[[g]]))

                # no rescale!

                # icov and cov.log.det
                out <- lav_samplestats_icov(COV = res.cov[[g]], ridge = ridge,
                       x.idx = x.idx[[g]],
                       ngroups = ngroups, g = g, warn = TRUE)
                res.icov[[g]] <- out$icov
                res.cov.log.det[[g]] <- out$cov.log.det

            # continuous + conditional.x = FALSE
            } else {
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
                out <- lav_samplestats_icov(COV = cov[[g]], ridge = ridge,
                           x.idx = x.idx[[g]], ngroups = ngroups, g = g,
                           warn = TRUE)
                icov[[g]] <- out$icov; cov.log.det[[g]] <- out$cov.log.det

                # fixed.x?
                if(fixed.x) {
                    cov.x[[g]]  <- unclass(unname(sample.cov.x[[g]]))
                    mean.x[[g]] <- unclass(unname(sample.mean.x[[g]]))
                }
            }
        }

        # WLS.obs
        WLS.obs[[g]] <- lav_samplestats_wls_obs(mean.g = mean[[g]],
            cov.g = cov[[g]], var.g = var[[g]], th.g = th[[g]],
            th.idx.g = th.idx[[g]], res.int.g = res.int[[g]],
            res.cov.g = res.cov[[g]], res.var.g = res.var[[g]],
            res.th.g = res.th[[g]], res.slopes.g = res.slopes[[g]],
            group.w.g = group.w[[g]],
            categorical = categorical, conditional.x = conditional.x,
            meanstructure = meanstructure, slopestructure = conditional.x,
            group.w.free = group.w.free)

        # WLS.V
        if(!WLS.V.user) {
            if(estimator == "GLS") {
                # FIXME: in <0.5-21, we had
                #V11 <- icov[[g]]
                #    if(mimic == "Mplus") { # is this a bug in Mplus?
                #        V11 <- V11 * nobs[[g]]/(nobs[[g]]-1)
                #    }
                WLS.V[[g]] <- lav_samplestats_Gamma_inverse_NT(ICOV = icov[[g]],
                                                 COV  = cov[[g]],
                                                 MEAN = mean[[g]],
                                                 rescale = FALSE,
                                                 x.idx = x.idx[[g]],
                                                 fixed.x = fixed.x,
                                                 conditional.x = conditional.x,
                                                 meanstructure = meanstructure,
                                                 slopestructure = conditional.x)
            } else if(estimator == "ULS") {
                WLS.V[[g]] <- diag(length(WLS.obs[[g]]))
                WLS.VD[[g]] <- rep(1, length(WLS.obs[[g]]))
            } else if(estimator == "WLS" || estimator == "DWLS") {
                if(is.null(WLS.V[[g]]))
                    stop("lavaan ERROR: the (D)WLS estimator is only available with full data or with a user-provided WLS.V")
            }

            # group.w.free
            if(!is.null(WLS.V[[g]]) && group.w.free) {
                # FIXME!!!
                WLS.V[[g]] <- lav_matrix_bdiag( matrix(1, 1, 1), WLS.V[[g]] )
            }

        }

    } # ngroups

    # construct SampleStats object
    lavSampleStats <- new("lavSampleStats",
                       # sample moments
                       th           = th,
                       th.idx       = th.idx,
                       th.names     = th.names,
                       mean         = mean,
                       cov          = cov,
                       var          = var,

                       # residual (y | x)
                       res.cov      = res.cov,
                       res.var      = res.var,
                       res.th       = res.th,
                       res.th.nox   = res.th.nox,
                       res.slopes   = res.slopes,
                       res.int      = res.int,

                       # fixed.x
                       mean.x       = mean.x,
                       cov.x        = cov.x,

                       # other
                       bifreq       = bifreq,
                       group.w      = group.w,

                       # convenience
                       nobs         = nobs,
                       ntotal       = sum(unlist(nobs)),
                       ngroups      = ngroups,
                       x.idx        = x.idx,

                       # extra sample statistics
                       icov         = icov,
                       cov.log.det  = cov.log.det,
                       res.icov         = res.icov,
                       res.cov.log.det  = res.cov.log.det,
                       ridge            = ridge.eps,
                       WLS.obs          = WLS.obs,
                       WLS.V            = WLS.V,
                       WLS.VD           = WLS.VD,
                       NACOV            = NACOV,
                       NACOV.user       = NACOV.user,

                       # cluster/level
                       YLp              = YLp,

                       # missingness
                       missing.flag     = missing.flag.,
                       missing          = missing.,
                       missing.h1       = missing.h1.,
                       zero.cell.tables = zero.cell.tables
                      )

    lavSampleStats
}

# compute sample statistics, per missing pattern
lav_samplestats_missing_patterns <- function(Y = NULL, Mp = NULL, wt = NULL) {

    # coerce Y to matrix
    Y <- as.matrix(Y)

    if(is.null(Mp)) {
        Mp <- lav_data_missing_patterns(Y, sort.freq = FALSE, coverage = FALSE)
    }

    Yp <- vector("list", length = Mp$npatterns)

    # fill in pattern statistics
    for(p in seq_len(Mp$npatterns)) {

        # extract raw data for these cases
        RAW <- Y[Mp$case.idx[[p]], Mp$pat[p, ], drop = FALSE]

        # more than one case
        if (Mp$freq[p] > 1L) {
            if(!is.null(wt)) {
                out <- stats::cov.wt(RAW, wt = wt[Mp$case.idx[[p]]],
                                     method = "ML")
                SY <- out$cov
                MY <- out$center
            } else {
                MY <- colMeans(RAW)
                SY <- crossprod(RAW)/Mp$freq[p] - tcrossprod(MY)
            }
        }
        # only a single observation (no need to weight!)
        else {
            SY <- 0
            MY <- as.numeric(RAW)
        }

        if(!is.null(wt)) {
            FREQ <- sum( wt[Mp$case.idx[[p]]] )
        } else {
            FREQ <- Mp$freq[p]
        }

        # store sample statistics, var.idx and freq
        Yp[[p]] <- list(SY = SY, MY = MY, var.idx = Mp$pat[p,],
                        freq = FREQ)
    }

    Yp
}

# compute sample statistics, per cluster
lav_samplestats_cluster_patterns <- function(Y = NULL, Lp = NULL) {

    # coerce Y to matrix
    Y1 <- as.matrix(Y); N <- NROW(Y1); P <- NCOL(Y1)

    if(is.null(Lp)) {
        stop("lavaan ERROR: Lp is NULL")
    }

    # how many levels?
    nlevels <- length(Lp$cluster) + 1L

    # compute some sample statistics per level
    YLp     <- vector("list", length = nlevels)
    for(l in 2:nlevels) {
        ncluster.sizes  <- Lp$ncluster.sizes[[l]]
        cluster.size    <- Lp$cluster.size[[l]]
        cluster.sizes   <- Lp$cluster.sizes[[l]]
        nclusters       <- Lp$nclusters[[l]]
        both.idx        <- Lp$both.idx[[l]]
        within.idx      <- Lp$within.idx[[l]]
        between.idx     <- Lp$between.idx[[l]]
        cluster.idx     <- Lp$cluster.idx[[l]]
        cluster.size.ns <- Lp$cluster.size.ns[[l]]
        ov.idx.1        <- sort.int(c(both.idx, within.idx))
        ov.idx.2        <- sort.int(c(both.idx, between.idx))

        #s <- (N^2 - sum(cluster.size^2)) / (N*(nclusters - 1L))
        # same as
        s <- (N - sum(cluster.size^2)/N)/(nclusters - 1)
        # NOTE: must be (nclusters - 1), otherwise, s is not average cluster
        # size even in the balanced case

        Y1.means <- colMeans(Y1, na.rm = TRUE)
        Y1Y1 <- crossprod(Y1)
        both.idx <- all.idx <- seq_len(P)
        if(length(within.idx) > 0L ||
           length(between.idx) > 0L) {
            both.idx <- all.idx[-c(within.idx, between.idx)]
        }

        # cluster-means
        # WARNING: aggregate() converts to FACTOR (changing the ORDER!)
        Y2 <- unname(as.matrix(aggregate(Y1, by = list(cluster.idx),
                               FUN = mean, na.rm = TRUE)[,-1]))
        Y2c <- t( t(Y2) - Y1.means )

        # compute S.w
        Y1a <- Y1 - Y2[cluster.idx, , drop = FALSE]

        # NOTE: we divide here by 'N - nclusters'
        # - this is just arbitrary, we could as well divide by 'N'
        # - this matters when we compute the objective function,
        #   where we need to 'multiply' again with the same constant
        # - a slight advantage of the 'N - nclusters' is that the same
        #   constant is needed for multiplying sigma.w.logdet, so we can
        #   can combine them
        S.w <- lav_matrix_crossprod(Y1a) / (N - nclusters)
        #S.w <- lav_matrix_crossprod(Y1a) / N


        # S.b
        # three parts: within/within, between/between, between/within
        S.b <- lav_matrix_crossprod(Y2c * cluster.size, Y2c) / nclusters

        # what if (nj*S.b - (nj-s)*S.w)/s is not-pd?
        #NJ <- max(cluster.size)
        #Sigma.j.max <- (NJ*S.b - (NJ-s)*S.w)/s
        #EV <- eigen(Sigma.j.max, symmetric = TRUE, only.values = TRUE)$values
        #if(any(EV < 0)) {
        #    # 1. spit out warning
        #    warning("lavaan WARNING: Sigma.j.max is not positive-definite.")
        #}


        S <- cov(Y1, use = "pairwise.complete.obs") * (N - 1L)/N
        S.PW.start <- S.w
        if(length(within.idx) > 0L) {
             S.PW.start[within.idx, within.idx] <-
                      S[within.idx, within.idx, drop = FALSE]
        }

        if(length(between.idx) > 0L) {
            S.w[between.idx,] <- 0
            S.w[,between.idx] <- 0
            S.PW.start[between.idx,] <- 0
            S.PW.start[,between.idx] <- 0
        }

        if(length(between.idx) > 0L) {
            # this is what is needed for MUML:
            S.b[, between.idx] <-
                (s * nclusters/N) * S.b[, between.idx, drop = FALSE]
            S.b[between.idx, ] <-
                (s * nclusters/N) * S.b[between.idx, , drop = FALSE]
            S.b[between.idx, between.idx] <-
                ( s * lav_matrix_crossprod(Y2c[, between.idx, drop = FALSE],
                         Y2c[, between.idx, drop = FALSE]) / nclusters  )
        }

        Sigma.B <- (S.b - S.w)/s
        Sigma.B[within.idx,] <- 0
        Sigma.B[,within.idx] <- 0


        Mu.W <- numeric( P )
        Mu.W[within.idx] <- Y1.means[within.idx]

        Mu.B <- Y1.means
        Mu.B[within.idx] <- 0
        if(length(between.idx) > 0L) {
            # replace between.idx by cov(Y2)[,] elements...
            Mu.B[between.idx] <- colMeans(Y2[,between.idx,drop = FALSE],
                                          na.rm = TRUE)

            S2 <- ( cov(Y2, use = "pairwise.complete.obs") *
                     (nclusters - 1L) / nclusters )

            Sigma.B[  between.idx, between.idx] <-
              S2[between.idx, between.idx, drop = FALSE]
        }

        # FIXME: Mu.B not quite ok for (fixed.x) x variables if they
        # occur both at level 1 AND level 2
        Mu.B.start <- Mu.B
        #Mu.B.start[both.idx] <- Mu.B.start[both.idx] - colMeans(Y2c[,both.idx])

        # per cluster-size
        cov.d  <- vector("list", length = ncluster.sizes)
        mean.d <- vector("list", length = ncluster.sizes)

        for(clz in seq_len(ncluster.sizes)) {
            nj <- cluster.sizes[clz]
            # select clusters with this size
            d.idx <- which(cluster.size == nj)
            ns <- length(d.idx)
            # NOTE:!!!!
            # reorder columns
            # to match A.inv and m.k later on in objective!!!
            tmp2 <- Y2[d.idx,
                       c(between.idx, sort.int(c(both.idx, within.idx))),
                       drop = FALSE]
            mean.d[[clz]] <- colMeans(tmp2, na.rm = TRUE)
            if(length(d.idx) > 1L) {
                cov.d[[clz]] <- ( cov(tmp2, use = "pairwise.complete.obs") *
                                      (ns-1) / ns )
            } else {
                cov.d[[clz]] <- 0
            }
        }

        YLp[[l]] <- list(Y1Y1 = Y1Y1,
                         Y2 = Y2, s = s, S.b = S.b, S.PW.start = S.PW.start,
                         Sigma.W = S.w, Mu.W = Mu.W,
                         Sigma.B = Sigma.B, Mu.B = Mu.B,
                         Mu.B.start = Mu.B.start,
                         mean.d = mean.d, cov.d = cov.d)
    } # l

    YLp
}
