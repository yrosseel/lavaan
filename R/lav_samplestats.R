# constructor for the 'lavSampleStats' class
#
# initial version: YR 25/03/2009
# major revision: YR 5/11/2011: separate data.obs and sample statistics
# YR 5/01/2016: add rescov, resvar, ... if conditional.x = TRUE

lav_samplestats_from_data <- function(lavdata           = NULL,
                                      DataX             = NULL,
                                      DataeXo           = NULL,
                                      DataOvnames       = NULL,
                                      DataOvnamesx      = NULL,
                                      DataOv            = NULL,
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

    # get X and Mp
    if(!is.null(lavdata)) {
        X <- lavdata@X; Mp <- lavdata@Mp
        ngroups <- lavdata@ngroups
        nlevels <- lavdata@nlevels
        nobs <- lavdata@nobs
        ov.names <- lavdata@ov.names
        ov.names.x <- lavdata@ov.names.x
        DataOv <- lavdata@ov
        eXo <- lavdata@eXo
    } else if(!is.null(DataX)) {
        stopifnot(is.list(DataX), is.matrix(DataX[[1L]]))
        X <- DataX
        eXo <- DataeXo
        ngroups <- length(X)
        nlevels <- 1L # for now
        Mp   <- vector("list", length = ngroups)
        nobs <- vector("list", length = ngroups)
        for(g in 1:ngroups) {
            if(missing != "listwise") {
                Mp[[g]] <- lav_data_missing_patterns(X[[g]], sort.freq = FALSE,
                               coverage = FALSE)
            }
            nobs[[g]] <- nrow(X[[g]])
        }
        ov.names   <- DataOvnames
        ov.names.x <- DataOvnamesx
    } else {
        stop("both lavdata and DataX argument are NULL")
    }

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
            if(nlevels > 1L) {
                stop("lavaan ERROR: multilevel + categorical not supported yet.")
            } else {
                categorical <- TRUE
            }
        }

        if(categorical) {
            if(estimator %in% c("ML","REML","PML","FML","MML","none")) {
                WLS.W <- FALSE
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

        # fill in the other slots
        if(!is.null(eXo[[g]])) {
            cov.x[[g]] <- cov(eXo[[g]], use="pairwise")
            if(rescale) {
                # we 'transform' the sample cov (divided by n-1) 
                # to a sample cov divided by 'n'
                cov.x[[g]] <- (nobs[[g]]-1)/nobs[[g]] * cov.x[[g]]
            }
            mean.x[[g]] <- colMeans(eXo[[g]])
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

            # 1. adapt lav_samplestats_cluster_patterns to use 
            #    within.idx, between.idx, ...
            # using S.PW for cov[[1]]
            cov[[g]] <- unname( stats::cov(X[[g]], use="pairwise"))
            #YLp[[g]][[2]]$S.PW
            mean[[g]] <- unname( colMeans(X[[g]], na.rm=TRUE) )
            var[[g]] <- diag(cov[[g]])

        } else { # continuous, single-level case
 
            if(conditional.x) {
                # residual covariances!

                # FIXME: how to handle missing data here?
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
                stopifnot(!conditional.x) # for now
                missing.flag. <- FALSE #!!! just use sample statistics
                missing.[[g]] <-
                    lav_samplestats_missing_patterns(Y  = X[[g]],
                                                     Mp = Mp[[g]])
                out <- lav_mvnorm_missing_h1_estimate_moments(Y = X[[g]],
                          Mp = Mp[[g]], Yp = missing.[[g]], verbose = verbose)
                missing.h1.[[g]]$sigma <- out$Sigma
                missing.h1.[[g]]$mu    <- out$Mu
                missing.h1.[[g]]$h1    <- out$fx

                # here, sample statistics == EM estimates
                cov[[g]]  <- missing.h1.[[g]]$sigma
                var[[g]]  <- diag(cov[[g]])
                mean[[g]] <- missing.h1.[[g]]$mu
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


            # icov and cov.log.det (but not if missing)
            if(missing != "ml") {
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

        # if missing = "fiml", sample statistics per pattern
        if(missing == "ml") {
            #if(nlevels > 1L) {
            #    stop("lavaan ERROR: multilevel + fiml not supported yet")
            #}
            if(conditional.x) {
                stop("lavaan ERROR: multilevel + conditional.x not supported yet")
            }
            stopifnot(!conditional.x) # for now
            missing.flag. <- TRUE
            missing.[[g]] <- 
                lav_samplestats_missing_patterns(Y  = X[[g]],
                                                 Mp = Mp[[g]])

            #cat("missing.h1 = "); print(missing.h1); cat("\n")
            if(missing.h1 && nlevels == 1L) {
                # estimate moments unrestricted model
                out <- lav_mvnorm_missing_h1_estimate_moments(Y = X[[g]],
                          Mp = Mp[[g]], Yp = missing.[[g]], verbose = verbose)
                missing.h1.[[g]]$sigma <- out$Sigma
                missing.h1.[[g]]$mu    <- out$Mu
                missing.h1.[[g]]$h1    <- out$fx
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
                NACOV[[g]] <- 
                    lav_samplestats_Gamma(Y              = Y,
                                          x.idx          = x.idx[[g]],
                                          # FALSE for now, until we use to
                                          # compute SB
                                          #fixed.x       = fixed.x,
                                          conditional.x  = conditional.x,
                                          meanstructure  = meanstructure,
                                          slopestructure = conditional.x,
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
                    NACOV[[g]] <-
                        lav_samplestats_Gamma(Y             = Y,
                                              x.idx         = x.idx[[g]],
                                              fixed.x       = fixed.x,
                                              conditional.x = conditional.x,
                                              meanstructure = meanstructure,
                                              slopestructure = conditional.x,
                                              Mplus.WLS = (mimic=="Mplus"))
                } else { # categorical case
                    NACOV[[g]]  <- CAT$WLS.W  * (nobs[[g]] - 1L)
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
                            if(is.complex(ev) || any(Re(ev) < 0)) {
                               stop("lavaan ERROR: Gamma (NACOV) matrix is not positive-definite")
                            }
                            WLS.V[[g]] <- inv.chol(NACOV[[g]])
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
                                         sample.nobs   = NULL,
                                         rescale       = FALSE,
                                         ov.names      = NULL,
                                         estimator     = "ML",
                                         mimic         = "lavaan",
                                         WLS.V         = NULL,
                                         NACOV         = NULL,
                                         ridge         = 1e-5,
                                         meanstructure = FALSE,
                                         group.w.free  = FALSE) {

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

    # for now, we do NOT support categorical data (using moments only),
    # fixed.x, and conditional.x
    categorical <- FALSE
    fixed.x <- FALSE
    conditional.x <- FALSE


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

    if(is.null(NACOV)) {
        NACOV      <- vector("list", length = ngroups)
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

        # FIXME: if the user provides x.idx, we could use this!!!
        # exogenous x?
        x.idx[[g]] <- integer(0L)
        conditional.x <- FALSE
        fixed.x <- FALSE

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
        idx <- match(ov.names[[g]], cov.names)
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

        # FIXME: create res.cov if conditional.x = TRUE!!!
        stopifnot(!conditional.x)

        # icov and cov.log.det
        out <- lav_samplestats_icov(COV = cov[[g]], ridge = ridge,
                   x.idx = x.idx[[g]],
                   ngroups = ngroups, g = g, warn = TRUE)
        icov[[g]] <- out$icov; cov.log.det[[g]] <- out$cov.log.det

        # the same for res.cov if conditional.x = TRUE
        if(conditional.x) {
            out <- lav_samplestats_icov(COV = res.cov[[g]], ridge = ridge,
                       x.idx = x.idx[[g]],
                       ngroups = ngroups, g = g, warn = TRUE)
            res.icov[[g]] <- out$icov
            res.cov.log.det[[g]] <- out$cov.log.det
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

        # NACOV
        # NACOV (=GAMMA)
        #if(!NACOV.user) {
            # nothing to do here; only used if provided by user
        #}

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
lav_samplestats_missing_patterns <- function(Y = NULL, Mp = NULL) {

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
            MY <- colMeans(RAW)
            SY <- crossprod(RAW)/Mp$freq[p] - tcrossprod(MY)
        }
        # only a single observation
        else {
            SY <- 0
            MY <- as.numeric(RAW)
        }

        # store sample statistics, var.idx and freq
        Yp[[p]] <- list(SY = SY, MY = MY, var.idx = Mp$pat[p,],
                        freq = Mp$freq[p])
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

        Y1.means <- colMeans(Y1, na.rm = TRUE)
        S <- cov(Y1, use = "pairwise.complete.obs") * (N - 1L)/N
        both.idx <- all.idx <- seq_len(P)
        if(length(within.idx) > 0L ||
           length(between.idx) > 0L) {
            both.idx <- all.idx[-c(within.idx, between.idx)]
        }

        # WARNING: aggregate() converts to FACTOR (changing the ORDER!)
        Y2 <- unname(as.matrix(aggregate(Y1, by = list(cluster.idx), 
                               FUN = mean, na.rm = TRUE)[,-1]))
        Y2c <- t( t(Y2) - Y1.means )
        if(length(within.idx) > 0L) {
        #    Y2[, within.idx] <- 0
        #    Y2c[,within.idx] <- 0
        }

        # compute S.w
        S.w <- matrix(0, P, P)
        Y1a <- Y1 - Y2[cluster.idx, , drop = FALSE]
        S.w <- lav_matrix_crossprod(Y1a) / (N - nclusters)
        S.PW.start <- S.w

        if(length(within.idx) > 0L) {
            #bw.idx <- c(both.idx, within.idx)
            #S.w[bw.idx, within.idx] <- S[bw.idx, within.idx]
            #S.w[within.idx, bw.idx] <- S[within.idx, bw.idx]

            #S.w[within.idx, within.idx] <- 
            #  S[within.idx, within.idx, drop = FALSE]
             S.PW.start[within.idx, within.idx] <- 
                      S[within.idx, within.idx, drop = FALSE]
        }
        if(length(between.idx) > 0L) {
            S.w[between.idx,] <- 0
            S.w[,between.idx] <- 0
            S.PW.start[between.idx,] <- 0
            S.PW.start[,between.idx] <- 0
        }
        #s <- (N^2 - sum(cluster.size^2)) / (N*(nclusters - 1L))
        # same as
        s <- (N - sum(cluster.size^2)/N)/(nclusters - 1)

        # S.b
        # three parts: within/within, between/between, between/within
        #S.b <- crossprod(Y2c * cluster.size, Y2c) / (nclusters - 1L)
        S.b <- lav_matrix_crossprod(Y2c * cluster.size, Y2c) / nclusters

        if(length(between.idx) > 0L) {
            # this is what is needed for MUML:
            S.b[, between.idx] <- 
                (s * nclusters/N) * S.b[, between.idx, drop = FALSE]
            S.b[between.idx, ] <- 
                (s * nclusters/N) * S.b[between.idx, , drop = FALSE]
            S.b[between.idx, between.idx] <-
                ( s * lav_matrix_crossprod(Y2c[, between.idx, drop = FALSE],
              #          Y2c[, between.idx, drop = FALSE]) / (nclusters - 1L))
                         Y2c[, between.idx, drop = FALSE]) / nclusters  )
        }


        V2 <- (S.b - S.w)/s
        V2[within.idx,] <- 0
        V2[,within.idx] <- 0

        Mu.W <- numeric( P )
        Mu.W[within.idx] <- Y1.means[within.idx]

        Mu.B <- colMeans(Y1, na.rm = TRUE)
        Mu.B[within.idx] <- 0
        if(length(between.idx) > 0L) {
            # replace between.idx by cov(Y2)[,] elements...
            Mu.B[between.idx] <- colMeans(Y2[,between.idx,drop = FALSE], 
                                          na.rm = TRUE)

            S2 <- ( cov(Y2, use = "pairwise.complete.obs") * 
                     (nclusters - 1L) / nclusters )
            #S2 <- cov(Y2)

            #bb.idx <- c(both.idx, between.idx)
            #V2[bb.idx, between.idx] <- S2[bb.idx, between.idx]
            #V2[between.idx, bb.idx] <- S2[between.idx, bb.idx]

            V2[  between.idx, between.idx] <- 
              S2[between.idx, between.idx, drop = FALSE]
        }

        # FIXME: Mu.B not quite ok for (fixed.x) x variables if they 
        # occur both at level 1 AND level 2
        Mu.B.start <- Mu.B
        #Mu.B.start[both.idx] <- Mu.B.start[both.idx] - colMeans(Y2c[,both.idx])

        # per cluster-size
        cov.d  <- vector("list", length = ncluster.sizes)
        mean.d <- vector("list", length = ncluster.sizes)

        ### FIXME  ####
        # WARNING: aggregate() converts to FACTOR (changing the ORDER!)
        Y2 <- unname(as.matrix(aggregate(Y1, by = list(cluster.idx), 
                                         FUN = mean, na.rm = TRUE)[,-1]))
        ### NO zero columns for within.idx columns #####

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

        YLp[[l]] <- list(Y2 = Y2, s = s, S.b = S.b, S.PW.start = S.PW.start,
                         Sigma.W = S.w, Mu.W = Mu.W,
                         Sigma.B = V2, Mu.B = Mu.B, Mu.B.start = Mu.B.start,
                         mean.d = mean.d, cov.d = cov.d)
    } # l

    YLp
}
