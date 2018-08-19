# lavPredict() contains a collection of `predict' methods
# the unifying theme is that they all rely on the (unknown, to be estimated)
# or (known, apriori specified) values for the latent variables
#
# lv: lavtent variables (aka `factor scores')
# ov: predict linear part of y_i
#
# - YR 11 June 2013: first version, in order to get factor scores for the
#                    categorical case
# - YR 12 Jan 2014: refactoring + lav_predict_fy (to be used by estimator MML)
#

# overload standard R function `predict'
setMethod("predict", "lavaan",
function(object, newdata = NULL) {
    lavPredict(object = object, newdata = newdata, type="lv", method="EBM",
               fsm = FALSE,
               optim.method = "bfgs")
})

# main function
lavPredict <- function(object, type = "lv", newdata = NULL, method = "EBM",
                       se = "none", label = TRUE, fsm = FALSE, level = 1L,
                       optim.method = "bfgs", ETA = NULL) {

    stopifnot(inherits(object, "lavaan"))
    lavmodel       <- object@Model
    lavdata        <- object@Data
    lavsamplestats <- object@SampleStats
    lavimplied     <- object@implied
    lavpta         <- object@pta

    # type
    type <- tolower(type)
    if(type %in% c("latent", "lv", "factor", "factor.score", "factorscore"))
        type <- "lv"
    if(type %in% c("ov","yhat"))
        type <- "yhat"

    # se?
    if(se != "none") {
        if(is.logical(se) && se) {
            se <- "standard"
        }
        if(type != "lv") {
            stop("lavaan ERROR: standard errors only available if type = \"lv\"")
        }
    }

    # need full data set supplied
    if(is.null(newdata)) {
        # use internal copy:
        if(lavdata@data.type != "full") {
            stop("lavaan ERROR: sample statistics were used for fitting and newdata is empty")
        } else if(is.null(lavdata@X[[1]])) {
            stop("lavaan ERROR: no local copy of data; FIXME!")
        } else {
            data.obs <- lavdata@X
        }
        eXo <- lavdata@eXo
    } else {
        OV <- lavdata@ov
        newData <- lavData(data        = newdata,
                           group       = lavdata@group,
                           ov.names    = lavdata@ov.names,
                           ov.names.x  = lavdata@ov.names.x,
                           ordered     = OV$name[ OV$type == "ordered" ],
                           lavoptions  = list(std.ov = lavdata@std.ov,
                                              group.label = lavdata@group.label,
                                              missing = lavdata@missing,
                                              warn = FALSE),
                           allow.single.case = TRUE)
        data.obs <- newData@X
        eXo <- newData@eXo
    }

    if(type == "lv") {
        if(!is.null(ETA)) {
            warning("lavaan WARNING: lvs will be predicted here; supplying ETA has no effect")
        }

        # post fit check (lv pd?)
        ok <- lav_object_post_check(object)
        #if(!ok) {
        #    stop("lavaan ERROR: lavInspect(,\"post.check\") is not TRUE; factor scores can not be computed. See the WARNING message.")
        #}

        out <- lav_predict_eta(lavobject = NULL, lavmodel = lavmodel,
                   lavdata = lavdata, lavsamplestats = lavsamplestats,
                   lavimplied = lavimplied, se = se, level = level,
                   data.obs = data.obs, eXo = eXo, method = method,
                   fsm = fsm, optim.method = optim.method)

        # extract fsm here
        if(fsm) {
            FSM <- attr(out, "fsm")
        }

        # extract se here
        if(se != "none") {
            SE <- attr(out, "se")
        }

        # remove dummy lv? (removes attr!)
        out <- lapply(seq_len(lavdata@ngroups), function(g) {
                   # determine block
                   if(lavdata@nlevels == 1L) {
                        bb <- g
                   } else {
                        bb <- (g - 1)*lavdata@nlevels + level
                   }
                   lv.idx <- c(lavmodel@ov.y.dummy.lv.idx[[bb]],
                               lavmodel@ov.x.dummy.lv.idx[[bb]])
                   ret <- out[[g]]
                   if(length(lv.idx) > 0L) {
                       ret <- out[[g]][, -lv.idx, drop=FALSE]
                   }
                   ret
               })

        if(fsm) {
            #FSM <- attr(out, "fsm")
            FSM <- lapply(seq_len(lavdata@ngroups), function(g) {
                       # determine block
                       if(lavdata@nlevels == 1L) {
                           bb <- g
                       } else {
                           bb <- (g - 1)*lavdata@nlevels + level
                       }
                       lv.idx <- c(lavmodel@ov.y.dummy.lv.idx[[bb]],
                                   lavmodel@ov.x.dummy.lv.idx[[bb]])
                       ov.idx <- c(lavmodel@ov.y.dummy.ov.idx[[bb]],
                                   lavmodel@ov.x.dummy.ov.idx[[bb]])
                       ret <- FSM[[g]]
                       if(length(lv.idx) > 0L) {
                           ret <- FSM[[g]][-lv.idx, -ov.idx, drop=FALSE]
                       }
                       ret
                   })
        }

        if(se != "none") {
            SE <- lapply(seq_len(lavdata@ngroups), function(g) {
                       # determine block
                       if(lavdata@nlevels == 1L) {
                           bb <- g
                       } else {
                           bb <- (g - 1)*lavdata@nlevels + level
                       }
                       lv.idx <- c(lavmodel@ov.y.dummy.lv.idx[[bb]],
                                   lavmodel@ov.x.dummy.lv.idx[[bb]])
                       ret <- SE[[g]]
                       if(length(lv.idx) > 0L) {
                           ret <- SE[[g]][, -lv.idx, drop=FALSE]
                       }
                       ret
                   })
        }

        # label?
        if(label) {
            for(g in seq_len(lavdata@ngroups)) {
                if(lavdata@nlevels > 1L && level == 2L) {
                    gg <- (g - 1)*lavdata@nlevels + 2L
                } else {
                    gg <- g
                }
                colnames(out[[g]]) <- lavpta$vnames$lv[[gg]]
            }
            if(se != "none") {
                if(lavdata@nlevels > 1L && level == 2L) {
                    gg <- (g - 1)*lavdata@nlevels + 2L
                } else {
                    gg <- g
                }
                colnames(SE[[g]]) <- lavpta$vnames$lv[[gg]]
            }
        }

    # estimated value for the observed indicators, given (estimated)
    # factor scores
    } else if(type == "yhat") {
        out <- lav_predict_yhat(lavobject = NULL, lavmodel = lavmodel,
                   lavdata = lavdata, lavsamplestats = lavsamplestats,
                   lavimplied = lavimplied,
                   data.obs = data.obs, eXo = eXo,
                   ETA = ETA, method = method, optim.method = optim.method)

        # label?
        if(label) {
            for(g in seq_len(lavdata@ngroups)) {
                colnames(out[[g]]) <- lavpta$vnames$ov[[g]]
            }
        }

    # density for each observed item, given (estimated) factor scores
    } else if(type == "fy") {
        out <- lav_predict_fy(lavobject = NULL, lavmodel = lavmodel,
                   lavdata = lavdata, lavsamplestats = lavsamplestats,
                   lavimplied = lavimplied,
                   data.obs = data.obs, eXo = eXo,
                   ETA = ETA, method = method, optim.method = optim.method)

        # label?
        if(label) {
            for(g in seq_len(lavdata@ngroups)) {
                colnames(out[[g]]) <- lavpta$vnames$ov[[g]]
            }
        }

    } else {
        stop("lavaan ERROR: type must be one of: lv yhat fy")
    }

    # lavaan.matrix
    out <- lapply(out, "class<-", c("lavaan.matrix", "matrix"))

    if(lavdata@ngroups == 1L) {
        out <- out[[1L]]
    } else {
        out
    }

    if(fsm) {
        attr(out, "fsm") <- FSM
    }

    if(se != "none") {
        attr(out, "se") <- SE
    }

    out
}

# internal function
lav_predict_eta <- function(lavobject = NULL,  # for convenience
                            # sub objects
                            lavmodel = NULL, lavdata = NULL,
                            lavsamplestats = NULL,
                            lavimplied = NULL,
                            # new data
                            data.obs = NULL, eXo = NULL,
                            # options
                            method = "EBM",
                            fsm = FALSE,
                            se = "none",
                            level = 1L,
                            optim.method = "bfgs") {

    # full object?
    if(inherits(lavobject, "lavaan")) {
        lavdata <- lavobject@Data
    } else {
        stopifnot(!is.null(lavdata))
    }

    # method
    method <- tolower(method)

    # alias
    if(method == "regression") {
        method <- "ebm"
    } else if(method == "bartlett" || method == "bartlet") {
        method <- "ml"
    }

    # normal case?
    if(all(lavdata@ov$type == "numeric")) {
        if(method == "ebm") {
            out <- lav_predict_eta_normal(lavobject = lavobject,
                       lavmodel = lavmodel, lavdata = lavdata,
                       lavimplied = lavimplied, se = se, level = level,
                       lavsamplestats = lavsamplestats,
                       data.obs = data.obs, eXo = eXo, fsm = fsm)
        } else if(method == "ml") {
            out <- lav_predict_eta_bartlett(lavobject = lavobject,
                       lavmodel = lavmodel, lavdata = lavdata,
                       lavimplied = lavimplied, se = se, level = level,
                       lavsamplestats = lavsamplestats,
                       data.obs = data.obs, eXo = eXo, fsm = fsm)
        } else {
            stop("lavaan ERROR: unkown method: ", method)
        }
    } else {
        if(method == "ebm") {
            out <- lav_predict_eta_ebm_ml(lavobject = lavobject,
                       lavmodel = lavmodel, lavdata = lavdata,
                       lavsamplestats = lavsamplestats, se = se,
                       level = level, data.obs = data.obs, eXo = eXo,
                       ML = FALSE, optim.method = optim.method)
        } else if(method == "ml") {
            out <- lav_predict_eta_ebm_ml(lavobject = lavobject,
                       lavmodel = lavmodel, lavdata = lavdata,
                       lavsamplestats = lavsamplestats, se = se,
                       level = level, data.obs = data.obs, eXo = eXo,
                       ML = TRUE, optim.method = optim.method)
        } else {
            stop("lavaan ERROR: unkown method: ", method)
        }
    }

    out
}


# factor scores - normal case
# NOTE: this is the classic 'regression' method; for the linear/continuous
#       case, this is equivalent to both EB and EBM
lav_predict_eta_normal <- function(lavobject = NULL,  # for convenience
                                   # sub objects
                                   lavmodel = NULL, lavdata = NULL,
                                   lavsamplestats = NULL,
                                   lavimplied = NULL,
                                   # optional new data
                                   data.obs = NULL, eXo = NULL,
                                   se = "none", level = 1L,
                                   fsm = FALSE) {

    # full object?
    if(inherits(lavobject, "lavaan")) {
        lavmodel       <- lavobject@Model
        lavdata        <- lavobject@Data
        lavsamplestats <- lavobject@SampleStats
        lavimplied     <- lavobject@implied
    } else {
        stopifnot(!is.null(lavmodel), !is.null(lavdata),
                  !is.null(lavsamplestats), !is.null(lavimplied))
    }

    if(is.null(data.obs)) {
        data.obs <- lavdata@X
        newdata.flag <- FALSE
    } else {
        newdata.flag <- TRUE
    }
    # eXo not needed

    # missings? and missing = "ml"?
    # impute values under the normal
    if(lavdata@missing %in% c("ml", "ml.x")) {
        for(g in seq_len(lavdata@ngroups)) {
            if(newdata.flag) {
                DATA <- data.obs[[g]]
                MP   <- lav_data_missing_patterns(data.obs[[g]])
            } else {
                DATA <- lavdata@X[[g]]
                MP   <- lavdata@Mp[[g]]
            }
            data.obs[[g]] <-
                lav_mvnorm_missing_impute_pattern(Y = DATA,
                                                  Mp = MP,
                                                  Mu = lavimplied$mean[[g]],
                                                  Sigma = lavimplied$cov[[g]])
        }
    }

    LAMBDA <- computeLAMBDA(lavmodel = lavmodel, remove.dummy.lv = FALSE)
    Sigma.hat <- lavimplied$cov
    Sigma.hat.inv <- lapply(Sigma.hat, MASS::ginv)
    VETA   <- computeVETA(lavmodel = lavmodel)
    EETA   <- computeEETA(lavmodel = lavmodel, lavsamplestats = lavsamplestats)
    EY     <- computeEY(  lavmodel = lavmodel, lavsamplestats = lavsamplestats)

    FS <- vector("list", length = lavdata@ngroups)
    if(fsm) {
        FSM <- vector("list", length = lavdata@ngroups)
    }
    if(se != "none") {
        SE <- vector("list", length = lavdata@ngroups)
    }

    for(g in 1:lavdata@ngroups) {

        if(lavdata@nlevels > 1L) {
            Lp <- lavdata@Lp[[g]]
            YLp <- lavsamplestats@YLp[[g]]

            # implied for this group
            group.idx <- (g - 1)*lavdata@nlevels + seq_len(lavdata@nlevels)
            implied.group <- lapply(lavimplied, function(x) x[group.idx])

            # random effects (=random intercepts or cluster means)
            out <- lav_mvnorm_cluster_implied22l(Lp = Lp,
                                                 implied = implied.group)
            MB.j <- lav_mvnorm_cluster_em_estep_ranef(YLp = YLp, Lp = Lp,
                        sigma.w = out$sigma.w, sigma.b = out$sigma.b,
                        sigma.zz = out$sigma.zz, sigma.yz = out$sigma.yz,
                        mu.z = out$mu.z, mu.w = out$mu.w, mu.b = out$mu.b,
                        se = FALSE)

            ov.idx <- Lp$ov.idx

            if(level == 1L) {
                data.W <- data.obs[[g]][, ov.idx[[1]] ]
                data.B <- MB.j[ Lp$cluster.idx[[2]], , drop = FALSE]

                # center
                data.obs.g <- data.W - data.B
            } else if(level == 2L) {
                Data.B <- matrix(0, nrow = nrow(MB.j),
                                    ncol = ncol(data.obs[[g]]))
                Data.B[, ov.idx[[1]] ] <- MB.j
                data.obs.g <- Data.B[, ov.idx[[2]] ]
            } else {
                stop("lavaan ERROR: only 2 levels are supported")
            }

            gg <- (g-1)*lavdata@nlevels + level
            VETA.g     <- VETA[[gg]]
            EETA.g     <- EETA[[gg]]
            LAMBDA.g   <- LAMBDA[[gg]]
            EY.g       <- EY[[gg]]
            Sigma.hat.inv.g <- Sigma.hat.inv[[gg]]

        } else {
            data.obs.g <- data.obs[[g]]
            VETA.g     <- VETA[[g]]
            EETA.g     <- EETA[[g]]
            LAMBDA.g   <- LAMBDA[[g]]
            EY.g       <- EY[[g]]
            Sigma.hat.inv.g <- Sigma.hat.inv[[g]]
        }


        nfac <- ncol(VETA[[g]])
        if(nfac == 0L) {
            FS[[g]] <- matrix(0, lavdata@nobs[[g]], nfac)
            next
        }

        # factor score coefficient matrix 'C'
        FSC <- VETA.g %*% t(LAMBDA.g) %*% Sigma.hat.inv.g
        if(fsm) {
            FSM[[g]] <- FSC
        }

        # standard error
        if(se == "standard") {
            tmp <- (VETA.g -
             VETA.g %*% t(LAMBDA.g) %*% Sigma.hat.inv.g %*% LAMBDA.g %*% VETA.g)
            tmp.d <- diag(tmp)
            tmp.d[ tmp.d < 1e-05 ] <- as.numeric(NA)
            SE[[g]] <- matrix(sqrt(tmp.d), nrow = 1L)
        }

        RES  <- sweep(data.obs.g,  MARGIN = 2L, STATS = EY.g,   FUN = "-")
        FS.g <- sweep(RES %*% t(FSC), MARGIN = 2L, STATS = EETA.g, FUN = "+")

        # replace values in dummy lv's by their observed counterpart
        if(length(lavmodel@ov.y.dummy.lv.idx[[g]]) > 0L && level == 1L) {
            FS.g[, lavmodel@ov.y.dummy.lv.idx[[g]] ] <-
                data.obs.g[, lavmodel@ov.y.dummy.ov.idx[[g]], drop = FALSE]
        }
        if(length(lavmodel@ov.x.dummy.lv.idx[[g]]) > 0L && level == 1L) {
            FS.g[, lavmodel@ov.x.dummy.lv.idx[[g]] ] <-
                data.obs.g[, lavmodel@ov.x.dummy.ov.idx[[g]], drop = FALSE]
        }

        FS[[g]] <- FS.g
    }

    if(fsm) {
        attr(FS, "fsm") <- FSM
    }
    if(se != "none") {
        attr(FS, "se") <- SE
    }

    FS
}

# factor scores - normal case - Bartlett method
# NOTES: 1) this is the classic 'Bartlett' method; for the linear/continuous
#           case, this is equivalent to 'ML'
#        2) the usual formula is:
#               FSC = solve(lambda' theta.inv lambda) (lambda' theta.inv)
#           BUT to deal with zero or negative variances, we use the
#           'GLS' version instead:
#               FSC = solve(lambda' sigma.inv lambda) (lambda' sigma.inv)
#           Reference: Bentler & Yuan (1997) 'Optimal Conditionally Unbiased
#                      Equivariant Factor Score Estimators'
#                      in Berkane (Ed) 'Latent variable modeling with
#                      applications to causality' (Springer-Verlag)
#        3) instead of solve(), we use MASS::ginv, for special settings where
#           -by construction- (lambda' sigma.inv lambda) is singular
lav_predict_eta_bartlett <- function(lavobject = NULL, # for convenience
                                     # sub objects
                                     lavmodel = NULL, lavdata = NULL,
                                     lavsamplestats = NULL,
                                     lavimplied = NULL,
                                     # optional new data
                                     data.obs = NULL, eXo = NULL,
                                     se = "none", level = 1L,
                                     fsm = FALSE) {

    # full object?
    if(inherits(lavobject, "lavaan")) {
        lavmodel       <- lavobject@Model
        lavdata        <- lavobject@Data
        lavsamplestats <- lavobject@SampleStats
        lavimplied     <- lavobject@implied
    } else {
        stopifnot(!is.null(lavmodel), !is.null(lavdata),
                  !is.null(lavsamplestats), !is.null(lavimplied))
    }

    if(is.null(data.obs)) {
        data.obs <- lavdata@X
        newdata.flag <- FALSE
    } else {
        newdata.flag <- TRUE
    }
    # eXo not needed

    # missings? and missing = "ml"?
    # impute values under the normal
    #
    # FIXME: THIS IS NOT CORRECT; we should use FIML!!!
    #
    if(lavdata@missing %in% c("ml", "ml.x")) {
        for(g in seq_len(lavdata@ngroups)) {
            if(newdata.flag) {
                DATA <- data.obs[[g]]
                MP   <- lav_data_missing_patterns(data.obs[[g]])
            } else {
                DATA <- lavdata@X[[g]]
                MP   <- lavdata@Mp[[g]]
            }
            data.obs[[g]] <-
                lav_mvnorm_missing_impute_pattern(Y = DATA,
                                                  Mp = MP,
                                                  Mu = lavimplied$mean[[g]],
                                                  Sigma = lavimplied$cov[[g]])
        }
    }


    LAMBDA <- computeLAMBDA(lavmodel = lavmodel, remove.dummy.lv = FALSE)
    Sigma.hat.inv <- lapply(lavimplied$cov, MASS::ginv)
    VETA   <- computeVETA(lavmodel = lavmodel) # for se only
    EETA   <- computeEETA(lavmodel = lavmodel, lavsamplestats = lavsamplestats)
    EY     <- computeEY(  lavmodel = lavmodel, lavsamplestats = lavsamplestats)

    FS <- vector("list", length = lavdata@ngroups)
    if(fsm) {
        FSM <- vector("list", length = lavdata@ngroups)
    }
    if(se != "none") {
        SE <- vector("list", length = lavdata@ngroups)
    }

    for(g in 1:lavdata@ngroups) {

        if(lavdata@nlevels > 1L) {
            Lp <- lavdata@Lp[[g]]
            YLp <- lavsamplestats@YLp[[g]]

            # implied for this group
            group.idx <- (g - 1)*lavdata@nlevels + seq_len(lavdata@nlevels)
            implied.group <- lapply(lavimplied, function(x) x[group.idx])

            # random effects (=random intercepts or cluster means)
            out <- lav_mvnorm_cluster_implied22l(Lp = Lp,
                                                 implied = implied.group)
            MB.j <- lav_mvnorm_cluster_em_estep_ranef(YLp = YLp, Lp = Lp,
                        sigma.w = out$sigma.w, sigma.b = out$sigma.b,
                        sigma.zz = out$sigma.zz, sigma.yz = out$sigma.yz,
                        mu.z = out$mu.z, mu.w = out$mu.w, mu.b = out$mu.b,
                        se = FALSE)

            ov.idx <- Lp$ov.idx

            if(level == 1L) {
                data.W <- data.obs[[g]][, ov.idx[[1]] ]
                data.B <- MB.j[ Lp$cluster.idx[[2]], , drop = FALSE]

                # center
                data.obs.g <- data.W - data.B
            } else if(level == 2L) {
                Data.B <- matrix(0, nrow = nrow(MB.j),
                                    ncol = ncol(data.obs[[g]]))
                Data.B[, ov.idx[[1]] ] <- MB.j
                data.obs.g <- Data.B[, ov.idx[[2]] ]
            } else {
                stop("lavaan ERROR: only 2 levels are supported")
            }

            gg <- (g-1)*lavdata@nlevels + level
            VETA.g     <- VETA[[gg]]
            EETA.g     <- EETA[[gg]]
            LAMBDA.g   <- LAMBDA[[gg]]
            EY.g       <- EY[[gg]]
            Sigma.hat.inv.g <- Sigma.hat.inv[[gg]]

        } else {
            data.obs.g <- data.obs[[g]]
            VETA.g     <- VETA[[g]]
            EETA.g     <- EETA[[g]]
            LAMBDA.g   <- LAMBDA[[g]]
            EY.g       <- EY[[g]]
            Sigma.hat.inv.g <- Sigma.hat.inv[[g]]
        }

        nfac <- length(EETA[[g]])
        if(nfac == 0L) {
            FS[[g]] <- matrix(0, lavdata@nobs[[g]], nfac)
            next
        }

        # factor score coefficient matrix 'C'
        FSC <- ( MASS::ginv(t(LAMBDA.g) %*% Sigma.hat.inv.g %*% LAMBDA.g)
                %*% t(LAMBDA.g) %*% Sigma.hat.inv.g )

        if(fsm) {
            FSM[[g]] <- FSC
        }

        # standard error
        if(se == "standard") {
            # the traditional formula is:
            #     solve(t(lambda) %*% solve(theta) %*% lambda)
            # but we replace it by
            #     solve( t(lambda) %*% solve(sigma) %*% lambda ) - psi
            # to handle negative variances
            # in addition, we use ginv
            tmp <- ( MASS::ginv(t(LAMBDA.g) %*%
                                  Sigma.hat.inv.g %*% LAMBDA.g)
                     - VETA.g )
            tmp.d <- diag(tmp)
            tmp.d[ tmp.d < 1e-05 ] <- as.numeric(NA)
            SE[[g]] <- matrix(sqrt(tmp.d), nrow = 1L)
        }

        RES  <- sweep(data.obs.g,  MARGIN = 2L, STATS = EY.g,   FUN = "-")
        FS.g <- sweep(RES %*% t(FSC), MARGIN = 2L, STATS = EETA.g, FUN = "+")

        # remove dummy lv's
        #lv.dummy.idx <- c(lavmodel@ov.y.dummy.lv.idx[[g]],
        #                  lavmodel@ov.x.dummy.lv.idx[[g]])
        #if(length(lv.dummy.idx) > 0L) {
        #    FS.g <- FS.g[,-lv.dummy.idx,drop=FALSE]
        #}

        # replace values in dummy lv's by their observed counterpart
        if(length(lavmodel@ov.y.dummy.lv.idx[[g]]) > 0L && level == 1L) {
            FS.g[, lavmodel@ov.y.dummy.lv.idx[[g]] ] <-
                data.obs[[g]][, lavmodel@ov.y.dummy.ov.idx[[g]], drop = FALSE]
        }
        if(length(lavmodel@ov.x.dummy.lv.idx[[g]]) > 0L && level == 1L) {
            FS.g[, lavmodel@ov.x.dummy.lv.idx[[g]] ] <-
                data.obs[[g]][, lavmodel@ov.x.dummy.ov.idx[[g]], drop = FALSE]
        }

        FS[[g]] <- FS.g
    }

    if(fsm) {
        attr(FS, "fsm") <- FSM
    }
    if(se != "none") {
        attr(FS, "se") <- SE
    }

    FS
}

# factor scores - EBM or ML
lav_predict_eta_ebm_ml <- function(lavobject = NULL,  # for convenience
                                   # sub objects
                                   lavmodel = NULL, lavdata = NULL,
                                   lavsamplestats = NULL,
                                   # optional new data
                                   data.obs = NULL, eXo = NULL,
                                   se = "none", level = 1L,
                                   ML = FALSE,
                                   optim.method = "bfgs") {

    optim.method <- tolower(optim.method)

    stopifnot(optim.method %in% c("nlminb", "bfgs"))

    ### FIXME: if all indicators of a factor are normal, can we not
    ###        just use the `classic' regression method??
    ###        (perhaps after whitening, to get uncorrelated factors...)

    # full object?
    if(inherits(lavobject, "lavaan")) {
        lavmodel       <- lavobject@Model
        lavdata        <- lavobject@Data
        lavsamplestats <- lavobject@SampleStats
    } else {
        stopifnot(!is.null(lavmodel), !is.null(lavdata),
                  !is.null(lavsamplestats))
    }

    # new data?
    if(is.null(data.obs)) {
        data.obs <- lavdata@X
    }
    if(is.null(eXo)) {
        eXo <- lavdata@eXo
    }

    # se?
    if(se != "none") {
        warning("lavaan WARNING: standard errors are not available (yet) for the non-normal case")
    }

    VETAx <- computeVETAx(lavmodel = lavmodel)
    VETAx.inv <- VETAx
    for(g in seq_len(lavdata@ngroups)) {
        if(nrow(VETAx[[g]]) > 0L) {
            VETAx.inv[[g]] <- solve(VETAx[[g]])
        }
    }
    EETAx <- computeEETAx(lavmodel = lavmodel, lavsamplestats = lavsamplestats,
                          eXo = eXo, nobs = lavdata@norig,
                          remove.dummy.lv = TRUE) ## FIXME?
    TH    <- computeTH(   lavmodel = lavmodel)
    THETA <- computeTHETA(lavmodel = lavmodel)

    # local objective function: x = lv values
    f.eta.i <- function(x, y.i, x.i, mu.i) {

        # add 'dummy' values (if any) for ov.y
        if(length(lavmodel@ov.y.dummy.lv.idx[[g]]) > 0L) {
            x2 <- c(x, data.obs[[g]][i,
                           lavmodel@ov.y.dummy.ov.idx[[g]], drop = FALSE])
        } else {
            x2 <- x
        }

        # conditional density of y, given eta.i(=x)
        log.fy <- lav_predict_fy_eta.i(lavmodel       = lavmodel,
                                       lavdata        = lavdata,
                                       lavsamplestats = lavsamplestats,
                                       y.i            = y.i,
                                       x.i            = x.i,
                                       eta.i          = matrix(x2, nrow=1L), # <---- eta!
                                       theta.sd       = theta.sd,
                                       th             = th,
                                       th.idx         = th.idx,
                                       log            = TRUE)

        if(ML) {
            # NOTE: 'true' ML is simply  -1*sum(log.fy)
            #  - but there is no upper/lower bound for the extrema:
            #    a pattern of all (in)correct drives the 'theta' parameter
            #    towards +/- Inf
            # - therefore, we add a vague prior, just to stabilize
            #
            diff <- t(x) - mu.i
            V <- diag( length(x) ) * 1e-05
            tmp <- as.numeric(0.5 * diff %*% V %*% t(diff))
            out <- 1 + tmp - sum(log.fy, na.rm=TRUE)
        } else {
            diff <- t(x) - mu.i
            V <- VETAx.inv[[g]]
            tmp <- as.numeric(0.5 * diff %*% V %*% t(diff))
            out <- tmp - sum(log.fy, na.rm=TRUE)
        }
        out
    }

    FS <- vector("list", length=lavdata@ngroups)
    for(g in seq_len(lavdata@ngroups)) {
        nfac <- ncol(VETAx[[g]])
        nfac2 <- nfac
        if(length(lavmodel@ov.y.dummy.lv.idx[[g]]) > 0L) {
            nfac2 <- nfac2 + length(lavmodel@ov.y.dummy.lv.idx[[g]])
        }
        FS[[g]] <- matrix(as.numeric(NA), nrow(data.obs[[g]]), nfac2)

        # special case: no regular lv's
        if(nfac == 0) {
            # impute dummy ov.y (if any)
            FS[[g]][, lavmodel@ov.y.dummy.ov.idx[[g]] ] <-
                data.obs[[g]][, lavmodel@ov.y.dummy.ov.idx[[g]], drop = FALSE]
            next
        }

        ## FIXME: factor scores not identical (but close) to Mplus
        #         if delta elements not equal to 1??
        mm.in.group <- 1:lavmodel@nmat[g] + cumsum(c(0,lavmodel@nmat))[g]
        MLIST     <- lavmodel@GLIST[ mm.in.group ]

        # check for negative values
        neg.var.idx <- which(diag(THETA[[g]]) < 0)
        if(length(neg.var.idx) > 0) {
            warning("lavaan WARNING: factor scores could not be computed due to at least one negative (residual) variance")
            next
        }

        # common values
        theta.sd <- sqrt(diag(THETA[[g]]))
        th       <- TH[[g]]
        th.idx   <- lavmodel@th.idx[[g]]

        # casewise for now
        N <- nrow(data.obs[[g]])
        for(i in 1:N) {

            # eXo?
            if(!is.null(eXo[[g]])) {
                 x.i <- eXo[[g]][i,,drop=FALSE]
            } else {
                 x.i <- NULL
            }
            mu.i <- EETAx[[g]][i,,drop=FALSE]
            y.i <- data.obs[[g]][i,,drop=FALSE]

            START <- numeric(nfac) # initial values for eta

            if(!all(is.na(y.i))) {
                # find best values for eta.i
                if(optim.method == "nlminb") {
                    out <- nlminb(start=START, objective=f.eta.i,
                                  gradient=NULL, # for now
                                  control=list(rel.tol=1e-8),
                                  y.i=y.i, x.i=x.i, mu.i=mu.i)
                } else if(optim.method == "bfgs") {
                    out <- optim(par = START, fn = f.eta.i,
                                 gr = NULL,
                                 control = list(reltol = 1e-8, fnscale = 1.1),
                                 method = "BFGS",
                                 y.i = y.i, x.i = x.i, mu.i = mu.i)
                }
                if(out$convergence == 0L) {
                    eta.i <- out$par
                } else {
                    eta.i <- rep(as.numeric(NA), nfac)
                }
            } else {
                eta.i <- rep(as.numeric(NA), nfac)
            }

            # add dummy ov.y lv values
            if(length(lavmodel@ov.y.dummy.lv.idx[[g]]) > 0L) {
                eta.i <- c(eta.i, data.obs[[g]][i,
                       lavmodel@ov.y.dummy.ov.idx[[g]], drop = FALSE])
            }

            FS[[g]][i,] <- eta.i
        }
    }

    FS
}

# predicted value for response y*_i, conditional on the predicted latent
# variable scores
# `measurement part':
#     y*_i = nu + lambda eta_i + K x_i + epsilon_i
#
#    where eta_i = latent variable value for i (either given or from predict)
#
# Two types: 1) nrow(ETA) = nrow(X) (factor scores)
#            2) nrow(ETA) = 1L (given values)
#
# in both cases, we return [nobs x nvar] matrix per group
lav_predict_yhat <- function(lavobject = NULL, # for convience
                             # sub objects
                             lavmodel = NULL, lavdata = NULL,
                             lavsamplestats = NULL,
                             lavimplied = NULL,
                             # new data
                             data.obs = NULL, eXo = NULL,
                             # ETA values
                             ETA = NULL,
                             # options
                             method = "EBM",
                             duplicate = FALSE,
                             optim.method = "bfgs") {

    # full object?
    if(inherits(lavobject, "lavaan")) {
        lavmodel       <- lavobject@Model
        lavdata        <- lavobject@Data
        lavsamplestats <- lavobject@SampleStats
        lavimplied     <- lavobject@implied
    } else {
        stopifnot(!is.null(lavmodel), !is.null(lavdata),
                  !is.null(lavsamplestats), !is.null(lavimplied))
    }

    # new data?
    if(is.null(data.obs)) {
        data.obs <- lavdata@X
    }
    if(is.null(eXo)) {
        eXo <- lavdata@eXo
    }

    # do we get values for ETA? If not, use `predict' to get plausible values
    if(is.null(ETA)) {
        ETA <- lav_predict_eta(lavobject = NULL, lavmodel = lavmodel,
                   lavdata = lavdata, lavsamplestats = lavsamplestats,
                   lavimplied = lavimplied,
                   data.obs = data.obs, eXo = eXo, method = method,
                   optim.method = optim.method)
    } else {
        # matrix
        if(is.matrix(ETA)) { # user-specified?
            if(nrow(ETA) == 1L) {
                tmp <- matrix(ETA, lavsamplestats@ntotal, length(ETA),
                              byrow = TRUE)
            } else if(nrow(ETA) != lavsamplestats@ntotal) {
                stop("lavaan ERROR: nrow(ETA) != lavsamplestats@ntotal")
            } else {
                tmp <- ETA
            }
            ETA <- lapply(1:lavdata@ngroups, function(i) tmp[lavdata@case.idx[[i]],])
        # vector: just 1 row of factor-scores
        } else if(is.numeric(ETA)) {
            # convert to matrix
            tmp <- matrix(ETA, lavsamplestats@ntotal, length(ETA), byrow = TRUE)
            ETA <- lapply(1:lavdata@ngroups, function(i) tmp[lavdata@case.idx[[i]],])
        } else if(is.list(ETA)) {
            stopifnot(lavdata@ngroups == length(ETA))
        }
    }

    YHAT <- computeYHAT(lavmodel = lavmodel, GLIST = NULL,
                        lavsamplestats = lavsamplestats, eXo = eXo,
                        nobs = lavdata@norig,
                        ETA = ETA, duplicate = duplicate)

    # if conditional.x, paste eXo
    if(lavmodel@categorical && !is.null(eXo)) {
         YHAT <- lapply(seq_len(lavdata@ngroups), function(g) {
                   ret <- cbind(YHAT[[g]], eXo[[g]])
                   ret })
    }

    YHAT
}

# conditional density y -- assuming independence!!
# f(y_i | eta_i, x_i) for EACH item
#
lav_predict_fy <- function(lavobject = NULL, # for convience
                           # sub objects
                           lavmodel = NULL, lavdata = NULL,
                           lavsamplestats = NULL,
                           lavimplied = NULL,
                           # new data
                           data.obs = NULL, eXo = NULL,
                           # ETA values
                           ETA = NULL,
                           # options
                           method = "EBM",
                           log. = FALSE,
                           optim.method = "bfgs") {

    # full object?
    if(inherits(lavobject, "lavaan")) {
        lavmodel       <- lavobject@Model
        lavdata        <- lavobject@Data
        lavsamplestats <- lavobject@SampleStats
        lavimplied     <- lavobject@implied
    } else {
        stopifnot(!is.null(lavmodel), !is.null(lavdata),
                  !is.null(lavsamplestats), !is.null(lavimplied))
    }

    # new data?
    if(is.null(data.obs)) {
        data.obs <- lavdata@X
    }
    if(is.null(eXo)) {
        eXo <- lavdata@eXo
    }

    # we need the YHATs (per group)
    YHAT <- lav_predict_yhat(lavobject = NULL, lavmodel = lavmodel,
                lavdata = lavdata, lavsamplestats = lavsamplestats,
                lavimplied = lavimplied,
                data.obs = data.obs, eXo = eXo, ETA = ETA, method = method,
                duplicate = FALSE, optim.method = optim.method)

    THETA <- computeTHETA(lavmodel = lavmodel)
    TH    <- computeTH(   lavmodel = lavmodel)

    # all normal?
    NORMAL <- all(lavdata@ov$type == "numeric")

    FY <- vector("list", length=lavdata@ngroups)
    for(g in seq_len(lavdata@ngroups)) {
        FY[[g]] <- lav_predict_fy_internal(X = data.obs[[g]], yhat = YHAT[[g]],
                       TH = TH[[g]], THETA = THETA[[g]],
                       num.idx = lavmodel@num.idx[[g]],
                       th.idx  = lavmodel@th.idx[[g]],
                       link    = lavmodel@link, log. = log.)
    }

    FY
}


# single group, internal function
lav_predict_fy_internal <- function(X = NULL, yhat = NULL,
                                    TH = NULL, THETA = NULL,
                                    num.idx = NULL, th.idx = NULL,
                                    link = NULL, log. = FALSE) {


    # shortcuts
    theta.var <- diag(THETA)

    # check size YHAT (either 1L or Nobs rows)
    if(! (nrow(yhat) == 1L || nrow(yhat) == nrow(X)) ) {
        stop("lavaan ERROR: nrow(YHAT[[g]]) not 1L and not nrow(X))")
    }

    FY.group <- matrix(0, nrow(X), ncol(X))
    #if(NORMAL) {
    #    if(nrow(yhat) == nrow(X)) {
    #        tmp <- (X - yhat)^2
    #    } else {
    #        tmp <- sweep(X, MARGIN=2, STATS=yhat, FUN="-")^2
    #    }
    #    tmp1 <- sweep(tmp, MARGIN=2, theta.var, "/")
    #    tmp2 <- exp( -0.5 * tmp1 )
    #    tmp3 <- sweep(tmp2, MARGIN=2, sqrt(2*pi*theta.var), "/")
    #    if(log.) {
    #        FY.group <- log(tmp3)
    #    } else {
    #        FY.group <- tmp3
    #    }
    #} else {
        # mixed items

    ord.idx <- unique( th.idx[th.idx > 0L] )

    # first, NUMERIC variables
    if(length(num.idx) > 0L) {
        # multivariate
        # FY.group[,num.idx] <-
        #    dmnorm(X[,num.idx],
        #           mean = yhat[n,num.idx],
            #           varcov = THETA[[g]][num.idx, num.idx], log = log.)
        for(v in num.idx) {
            FY.group[,v] <- dnorm(X[,v],
                                  # YHAT may change or not per case
                                  mean = yhat[,v],
                                  sd   = sqrt(theta.var[v]),
                                  log  = log.)
        }
    }

    # second, ORDERED variables
    for(v in ord.idx) {
        th.y <- TH[ th.idx == v ]; TH.Y <- c(-Inf, th.y, Inf)
        ncat <- length(th.y) + 1L
        fy <- numeric(ncat)
        theta.v <- sqrt(theta.var[v])
        yhat.v  <- yhat[,v]

        # two cases: yhat.v is a scalar, or has length = nobs
        fy <- matrix(0, nrow=length(yhat.v), ncol=ncat)

        # for each category
        for(k in seq_len(ncat)) {
            if(link == "probit") {
                fy[,k] = pnorm(  (TH.Y[k+1] - yhat.v) / theta.v) -
                         pnorm(  (TH.Y[k  ] - yhat.v) / theta.v)
            } else if(link == "logit") {
                fy[,k] = plogis( (TH.Y[k+1] - yhat.v) / theta.v) -
                         plogis( (TH.Y[k  ] - yhat.v) / theta.v)
            } else {
                stop("lavaan ERROR: link must be probit or logit")
            }
        }

        # underflow
        idx <- which(fy < .Machine$double.eps)
        if(length(idx) > 0L) {
            fy[idx] <- .Machine$double.eps
        }

        # log?
        if(log.) {
            fy <- log(fy)
        }

        # case-wise expansion/selection
        if(length(yhat.v) == 1L) {
            # expand category probabilities for all observations
            FY.group[,v] <- fy[1L, X[,v]]
        } else {
            # select correct category probability per observation
            FY.group[,v] <- fy[ cbind(seq_len(nrow(fy)), X[,v]) ]
        }
    } # ord

    FY.group
}



# conditional density y -- assuming independence!!
# f(y_i | eta_i, x_i)
#
# but for a SINGLE observation y_i (and x_i), for given values of eta_i
#
lav_predict_fy_eta.i <- function(lavmodel = NULL, lavdata = NULL,
                                 lavsamplestats = NULL,
                                 y.i = NULL, x.i = NULL,
                                 eta.i = NULL, theta.sd = NULL, g = 1L,
                                 th = NULL, th.idx = NULL, log = TRUE) {

    mm.in.group <- 1:lavmodel@nmat[g] + cumsum(c(0,lavmodel@nmat))[g]
    MLIST <- lavmodel@GLIST[ mm.in.group ]

    # linear predictor for all items
    YHAT <-
        computeEYetax.LISREL(MLIST             = MLIST,
                             eXo               = x.i,
                             ETA               = eta.i,
                             sample.mean       = lavsamplestats@mean[[g]],
                            ov.y.dummy.ov.idx = lavmodel@ov.y.dummy.ov.idx[[g]],
                            ov.x.dummy.ov.idx = lavmodel@ov.x.dummy.ov.idx[[g]],
                            ov.y.dummy.lv.idx = lavmodel@ov.y.dummy.lv.idx[[g]],
                            ov.x.dummy.lv.idx = lavmodel@ov.x.dummy.lv.idx[[g]])

    # P(y_i | eta_i, x_i) for all items
    if(all(lavdata@ov$type == "numeric")) {
        # NORMAL case
        FY <- dnorm(y.i, mean = YHAT, sd = theta.sd, log = log)
    } else {
        FY <- numeric(lavmodel@nvar[g])
        for(v in seq_len(lavmodel@nvar[g])) {
            if(lavdata@ov$type[v] == "numeric") {
                ### FIXME!!! we can do all numeric vars at once!!
                FY[v] <- dnorm(y.i[v], mean = YHAT[v], sd = theta.sd[v],
                               log = log)
            } else if(lavdata@ov$type[v] == "ordered") {
                # handle missing value
                if(is.na(y.i[v])) {
                    FY[v] <- as.numeric(NA)
                } else {
                    th.y <- th[ th.idx == v ]; TH.Y <-  c(-Inf, th.y, Inf)
                    k <- y.i[v]
                    p1 <- pnorm( (TH.Y[ k + 1 ] - YHAT[v])/theta.sd[v] )
                    p2 <- pnorm( (TH.Y[ k     ] - YHAT[v])/theta.sd[v] )
                    prob <- (p1 - p2)
                    if(prob < .Machine$double.eps) {
                       prob <- .Machine$double.eps
                    }
                    if(log) {
                        FY[v] <- log(prob)
                    } else {
                        FY[v] <- prob
                    }
                }
            } else {
                stop("lavaan ERROR: unknown type: `",
                      lavdata@ov$type[v], "' for variable: ",
                      lavdata@ov$name[v])
            }
        }
    }

    FY
}

