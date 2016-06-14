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
               optim.method = "nlminb")
})

# main function
lavPredict <- function(object, type = "lv", newdata = NULL, method = "EBM",
                       se.fit = FALSE, label = TRUE, fsm = FALSE,
                       optim.method = "nlminb") {

    stopifnot(inherits(object, "lavaan"))
    lavmodel       <- object@Model
    lavdata        <- object@Data
    lavsamplestats <- object@SampleStats
    lavpta         <- object@pta

    # type
    type <- tolower(type)
    if(type %in% c("latent", "lv", "factor", "factor.score", "factorscore"))
        type <- "lv"
    if(type %in% c("ov","yhat")) 
        type <- "yhat"

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
                           group.label = lavdata@group.label,
                           ov.names    = lavdata@ov.names,
                           ordered     = OV$name[ OV$type == "ordered" ],
                           ov.names.x  = lavdata@ov.names.x,
                           std.ov      = lavdata@std.ov,
                           missing     = lavdata@missing,
                           # warn      = FALSE,
                           allow.single.case = TRUE)
        data.obs <- newData@X
        eXo <- newData@eXo
    }

    if(type == "lv") {

        # post fit check (lv pd?)
        ok <- lav_object_post_check(object)
        #if(!ok) {
        #    stop("lavaan ERROR: lavInspect(,\"post.check\") is not TRUE; factor scores can not be computed. See the WARNING message.")
        #}

        out <- lav_predict_eta(lavobject = NULL, lavmodel = lavmodel,
                   lavdata = lavdata, lavsamplestats = lavsamplestats,
                   data.obs = data.obs, eXo = eXo, method = method,
                   fsm = fsm, optim.method = optim.method)

        # remove dummy lv?
        if(fsm) {
            FSM <- attr(out, "fsm")
        }
        out <- lapply(seq_len(lavdata@ngroups), function(g) {
                   lv.idx <- c(lavmodel@ov.y.dummy.lv.idx[[g]],
                               lavmodel@ov.x.dummy.lv.idx[[g]])
                   ret <- out[[g]]
                   if(length(lv.idx) > 0L) {
                       ret <- out[[g]][, -lv.idx, drop=FALSE]
                   }
                   ret 
               })

        # label?
        if(label) {
            for(g in seq_len(lavdata@ngroups)) {
                colnames(out[[g]]) <- lavpta$vnames$lv[[g]]
            }
        }

    # estimated value for the observed indicators, given (estimated)
    # factor scores
    } else if(type == "yhat") {
        out <- lav_predict_yhat(lavobject = NULL, lavmodel = lavmodel,
                   lavdata = lavdata, lavsamplestats = lavsamplestats,
                   data.obs = data.obs, eXo = eXo,
                   ETA = NULL, method = method, optim.method = optim.method)

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
                   data.obs = data.obs, eXo = eXo,
                   ETA = NULL, method = method, optim.method = optim.method)

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

    out
}

# internal function
lav_predict_eta <- function(lavobject = NULL,  # for convenience
                            # sub objects
                            lavmodel = NULL, lavdata = NULL,
                            lavsamplestats = NULL,
                            # new data
                            data.obs = NULL, eXo = NULL,
                            # options
                            method = "EBM",
                            fsm = FALSE,
                            optim.method = "nlminb") {

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
    }

    # normal case?
    if(all(lavdata@ov$type == "numeric")) {
        if(method == "ebm") {
            out <- lav_predict_eta_normal(lavobject = lavobject,
                       lavmodel = lavmodel, lavdata = lavdata, 
                       lavsamplestats = lavsamplestats,
                       data.obs = data.obs, eXo = eXo, fsm = fsm)
        } else if(method == "bartlett" || method == "bartlet") {
            out <- lav_predict_eta_bartlett(lavobject = lavobject,
                       lavmodel = lavmodel, lavdata = lavdata,
                       lavsamplestats = lavsamplestats,
                       data.obs = data.obs, eXo = eXo, fsm = fsm)
        } else {
            stop("lavaan ERROR: unkown method: ", method)
        }
    } else {
        out <- lav_predict_eta_ebm(lavobject = lavobject,
                   lavmodel = lavmodel, lavdata = lavdata,
                   lavsamplestats = lavsamplestats,
                   data.obs = data.obs, eXo = eXo,
                   optim.method = optim.method)
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
                                   # optional new data
                                   data.obs = NULL, eXo = NULL,
                                   fsm = FALSE) { 

    # full object?
    if(inherits(lavobject, "lavaan")) {
        lavmodel       <- lavobject@Model
        lavdata        <- lavobject@Data
        lavsamplestats <- lavobject@SampleStats
    } else {
        stopifnot(!is.null(lavmodel), !is.null(lavdata),
                  !is.null(lavsamplestats))
    }

    if(is.null(data.obs)) {
        data.obs <- lavdata@X
    }
    # eXo not needed

    Sigma.hat <- computeSigmaHat(lavmodel = lavmodel)
    Sigma.hat.inv <- lapply(Sigma.hat, solve)
    VETA   <- computeVETA(lavmodel = lavmodel, lavsamplestats = lavsamplestats)
    EETA   <- computeEETA(lavmodel = lavmodel, lavsamplestats = lavsamplestats)
    EY     <- computeEY(  lavmodel = lavmodel, lavsamplestats = lavsamplestats)
    LAMBDA <- computeLAMBDA(lavmodel = lavmodel, remove.dummy.lv = FALSE)
     
    FS <- vector("list", length = lavdata@ngroups)
    if(fsm) {
        FSM <- vector("list", length = lavdata@ngroups)
    }
    for(g in 1:lavdata@ngroups) {
        nfac <- ncol(VETA[[g]])
        if(nfac == 0L) {
            FS[[g]] <- matrix(0, lavdata@nobs[[g]], nfac)
            next
        }

        # factor score coefficient matrix 'C'
        FSC <- VETA[[g]] %*% t(LAMBDA[[g]]) %*% Sigma.hat.inv[[g]]
        if(fsm) {
            FSM[[g]] <- FSC
        }

        RES  <- sweep(data.obs[[g]],  MARGIN = 2L, STATS = EY[[g]],   FUN = "-")
        FS.g <- sweep(RES %*% t(FSC), MARGIN = 2L, STATS = EETA[[g]], FUN = "+")

        # remove dummy lv's
        #lv.dummy.idx <- c(lavmodel@ov.y.dummy.lv.idx[[g]],
        #                  lavmodel@ov.x.dummy.lv.idx[[g]])
        #if(length(lv.dummy.idx) > 0L) {
        #    FS.g <- FS.g[,-lv.dummy.idx,drop=FALSE]
        #}

        # replace values in dummy lv's by their observed counterpart
        if(length(lavmodel@ov.y.dummy.lv.idx[[g]]) > 0L) {
            FS.g[, lavmodel@ov.y.dummy.lv.idx[[g]] ] <-
                data.obs[[g]][, lavmodel@ov.y.dummy.ov.idx[[g]], drop = FALSE]
        }
        if(length(lavmodel@ov.x.dummy.lv.idx[[g]]) > 0L) {
            FS.g[, lavmodel@ov.x.dummy.lv.idx[[g]] ] <-
                data.obs[[g]][, lavmodel@ov.x.dummy.ov.idx[[g]], drop = FALSE]
        }

        FS[[g]] <- FS.g
    }

    if(fsm) {
        attr(FS, "fsm") <- FSM
    }

    FS
}

# factor scores - normal case - Bartlett method
# NOTE: this is the classic 'Bartlett' method; for the linear/continuous 
#       case, this is equivalent to 'ML'
lav_predict_eta_bartlett <- function(lavobject = NULL, # for convenience
                                     # sub objects
                                     lavmodel = NULL, lavdata = NULL, 
                                     lavsamplestats = NULL,
                                     # optional new data
                                     data.obs = NULL, eXo = NULL,
                                     fsm = FALSE) { 

    # full object?
    if(inherits(lavobject, "lavaan")) {
        lavmodel       <- lavobject@Model
        lavdata        <- lavobject@Data
        lavsamplestats <- lavobject@SampleStats
    } else {
        stopifnot(!is.null(lavmodel), !is.null(lavdata),
                  !is.null(lavsamplestats))
    }

    if(is.null(data.obs)) {
        data.obs <- lavdata@X
    }
    # eXo not needed

    LAMBDA <- computeLAMBDA(lavmodel = lavmodel, remove.dummy.lv = FALSE)
    THETA  <- computeTHETA(lavmodel = lavmodel)
    THETA.inv <- lapply(THETA, solve)

    EETA   <- computeEETA(lavmodel = lavmodel, lavsamplestats = lavsamplestats)
    EY     <- computeEY(  lavmodel = lavmodel, lavsamplestats = lavsamplestats)
     
    FS <- vector("list", length = lavdata@ngroups)
    if(fsm) {
        FSM <- vector("list", length = lavdata@ngroups)
    }
    for(g in 1:lavdata@ngroups) {
        nfac <- length(EETA[[g]])
        if(nfac == 0L) {
            FS[[g]] <- matrix(0, lavdata@nobs[[g]], nfac)
            next
        }

        # factor score coefficient matrix 'C'
        FSC = ( solve(t(LAMBDA[[g]]) %*% THETA.inv[[g]] %*% LAMBDA[[g]]) %*%
                t(LAMBDA[[g]]) %*% THETA.inv[[g]] )

        if(fsm) {
            FSM[[g]] <- FSC
        }

        RES  <- sweep(data.obs[[g]],  MARGIN = 2L, STATS = EY[[g]],   FUN = "-")
        FS.g <- sweep(RES %*% t(FSC), MARGIN = 2L, STATS = EETA[[g]], FUN = "+")

        # remove dummy lv's
        #lv.dummy.idx <- c(lavmodel@ov.y.dummy.lv.idx[[g]],
        #                  lavmodel@ov.x.dummy.lv.idx[[g]])
        #if(length(lv.dummy.idx) > 0L) {
        #    FS.g <- FS.g[,-lv.dummy.idx,drop=FALSE]
        #}

        # replace values in dummy lv's by their observed counterpart
        if(length(lavmodel@ov.y.dummy.lv.idx[[g]]) > 0L) {
            FS.g[, lavmodel@ov.y.dummy.lv.idx[[g]] ] <-
                data.obs[[g]][, lavmodel@ov.y.dummy.ov.idx[[g]], drop = FALSE]
        }
        if(length(lavmodel@ov.x.dummy.lv.idx[[g]]) > 0L) {
            FS.g[, lavmodel@ov.x.dummy.lv.idx[[g]] ] <-
                data.obs[[g]][, lavmodel@ov.x.dummy.ov.idx[[g]], drop = FALSE]
        }

        FS[[g]] <- FS.g
    }

    if(fsm) {
        attr(FS, "fsm") <- FSM
    }

    FS
}

# factor scores - EBM
lav_predict_eta_ebm <- function(lavobject = NULL,  # for convenience
                                # sub objects
                                lavmodel = NULL, lavdata = NULL,
                                lavsamplestats = NULL,
                                # optional new data
                                data.obs = NULL, eXo = NULL,
                                optim.method = "nlminb") {

    stopifnot(optim.method %in% c("nlminb", "BFGS"))

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

    VETAx <- computeVETAx(lavmodel = lavmodel)
    VETAx.inv <- VETAx
    for(g in seq_len(lavdata@ngroups)) {
        if(nrow(VETAx[[g]]) > 0L) {
            VETAx.inv[[g]] <- solve(VETAx[[g]])
        }
    }
    EETAx <- computeEETAx(lavmodel = lavmodel, lavsamplestats = lavsamplestats,
                          eXo = eXo, remove.dummy.lv = TRUE) ## FIXME? 
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

        diff <- t(x) - mu.i
        V <- VETAx.inv[[g]]
        # handle missing values: we skip them
        tmp <- as.numeric(0.5 * diff %*% V %*% t(diff))
        # handle missing values: we skip them: FIXME!!!
        out <- tmp - sum(log.fy, na.rm=TRUE)
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

            # find best values for eta.i
            if(optim.method == "nlminb") {
                out <- nlminb(start=START, objective=f.eta.i,
                              gradient=NULL, # for now
                              control=list(rel.tol=1e-8),
                              y.i=y.i, x.i=x.i, mu.i=mu.i)
            } else if(optim.method == "BFGS") {
                out <- optim(par = START, fn = f.eta.i,
                             gr = NULL,
                             control = list(reltol = 1e-8),
                             method = "BFGS",
                             y.i = y.i, x.i = x.i, mu.i = mu.i)
            }
            if(out$convergence == 0L) {
                eta.i <- out$par
            } else {
                eta.i <- rep(as.numeric(NA), nfac2)
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
                             # new data
                             data.obs = NULL, eXo = NULL,
                             # ETA values
                             ETA = NULL,
                             # options
                             method = "EBM", 
                             duplicate = FALSE,
                             optim.method = "nlminb") {

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

    # do we get values for ETA? If not, use `predict' to get plausible values
    if(is.null(ETA)) {
        ETA <- lav_predict_eta(lavobject = NULL, lavmodel = lavmodel,
                   lavdata = lavdata, lavsamplestats = lavsamplestats,
                   data.obs = data.obs, eXo = eXo, method = method,
                   optim.method = optim.method)
    } else {
        # list
        if(is.matrix(ETA)) { # user-specified?
            tmp <- ETA; ETA <- vector("list", length=lavdata@ngroups)
            ETA[seq_len(lavdata@ngroups)] <- list(tmp)
        } else if(is.list(ETA)) {
            stopifnot(lavdata@ngroups == length(ETA))
        }
    }

    YHAT <- computeYHAT(lavmodel = lavmodel, GLIST = NULL,
                        lavsamplestats = lavsamplestats, eXo = eXo,
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
                           # new data
                           data.obs = NULL, eXo = NULL,
                           # ETA values
                           ETA = NULL,
                           # options
                           method = "EBM",
                           log. = FALSE,
                           optim.method = "nlminb") {

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

    # we need the YHATs (per group)
    YHAT <- lav_predict_yhat(lavobject = NULL, lavmodel = lavmodel,
                lavdata = lavdata, lavsamplestats = lavsamplestats,
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
                stop("lavaan ERROR: unknown type: ", 
                      lavdata@ov$type[v], " for variable", 
                      lavdata@ov$name[v])
            }
        }
    }

    FY
}

