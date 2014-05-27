# inspect a fitted lavaan object

setMethod("inspect", "lavaan",
function(object, what="free") {

    if(length(what) > 1) {
        stop("`what' arguments contains multiple arguments; only one is allowed")
    }

    # be case insensitive
    what <- tolower(what)

    if(what == "free") {
        free.parameters(object, type="free")
    } else if(what == "unco") {
        free.parameters(object, type="unco")
    } else if(what == "user") {
        free.parameters(object, type="user")
    } else if(what == "se" ||
              what == "std.err" ||
              what == "standard.errors") {
        standard.errors(object)
    } else if(what == "start" ||
              what == "starting.values") {
        starting.values(object)
    } else if(what == "coef" ||
              what == "coefficients" ||
              what == "parameters" ||
              what == "parameter.estimates" ||
              what == "parameter.values" ||
              what == "estimates" ||
              what == "x" ||
              what == "est") {
        parameter.values(object)
    } else if(what == "list") {
        parameter.list(object)
    } else if(what == "dx" ||
              what == "gradient" ||
              what == "derivatives") {
        derivatives(object)
    } else if(what == "fit" ||
              what == "fitmeasures" ||
              what == "fit.measures" ||
              what == "fit.indices") {
        fitMeasures(object)
    } else if(what == "std" ||
              what == "std.coef" ||
              what == "standardized" ||
              what == "standardizedsolution" ||
              what == "standardized.solution") {
        standardizedSolution(object)
    } else if(what == "mi" ||
              what == "modindices" ||
              what == "modification" ||
              what == "modificationindices" ||
              what == "modification.indices") {
         modificationIndices(object)
    } else if(what == "samp" ||
              what == "sample" ||
              what == "samplestatistics" ||
              what == "sampstat") {
        sampStat(object)
    } else if(what == "data") {
        lav_inspect_data(object)
    } else if(what == "rsquare" || 
              what == "r-square" ||
              what == "r2") {
         rsquare(object)
    } else if(what == "wls.v") {
        getWLS.V(object, drop.list.single.group=TRUE)
    } else if(what == "nacov") {
        getSampleStatsNACOV(object)    
    } else if(what == "modelcovlv"  ||
              what == "modelcov.lv" ||
              what == "cov.lv") {
        getModelCovLV(object, correlation.metric=FALSE, labels=TRUE)
    } else if(what == "modelcorlv"  ||
              what == "modelcor.lv" ||
              what == "cor.lv") {
        getModelCorLV(object, labels=TRUE)
    } else if(what == "modelcovall"  ||
              what == "modelcov.all" ||
              what == "cov.all") {
        getModelCov(object, correlation.metric=FALSE, labels=TRUE)
    } else if(what == "modelcorall"  ||
              what == "modelcor.all" ||
              what == "cor.all") {
        getModelCor(object, labels=TRUE)
    } else if(what == "theta") {
        getModelTheta(object, labels=TRUE)
    } else if(what == "mean.lv") {
        getModelMeanLV(object, labels=TRUE)
    } else if(what == "coverage") {
        getDataMissingCoverage(object, labels=TRUE)
    } else if(what %in% c("patterns", "pattern")) {
        getDataMissingPatterns(object, labels=TRUE)
    } else if(what == "converged") {
        object@Fit@converged
    } else {
        stop("unknown `what' argument in inspect function: `", what, "'")
    }

})

free.parameters <- function(object, type="free") {

    GLIST <- object@Model@GLIST

    for(mm in 1:length(GLIST)) {
        # erase everything
        GLIST[[mm]][,] <- 0.0

        dimnames(GLIST[[mm]]) <- object@Model@dimNames[[mm]]

        # fill in free/unco parameter counts
        if(type == "free") {
            m.el.idx <- object@Model@m.free.idx[[mm]]
            x.el.idx <- object@Model@x.free.idx[[mm]]
        } else if(type == "unco") {
            m.el.idx <- object@Model@m.unco.idx[[mm]]
            x.el.idx <- object@Model@x.unco.idx[[mm]]
        } else if(type == "user") {
            m.el.idx <- object@Model@m.user.idx[[mm]]
            x.el.idx <- object@Model@x.user.idx[[mm]]
        } else {
            stop("lavaan ERROR: unknown type argument:", type, )
        }
        GLIST[[mm]][m.el.idx] <- x.el.idx

        # class
        if(object@Model@isSymmetric[mm]) {
            class(GLIST[[mm]]) <- c("lavaan.matrix.symmetric", "matrix")
        } else {
            class(GLIST[[mm]]) <- c("lavaan.matrix", "matrix")
        }
    }

    GLIST
}


standard.errors <- function(object) {

    GLIST <- object@Model@GLIST

    for(mm in 1:length(GLIST)) {
        # labels
        dimnames(GLIST[[mm]]) <- object@Model@dimNames[[mm]]

        # fill in starting values
        m.user.idx <- object@Model@m.user.idx[[mm]]
        x.user.idx <- object@Model@x.user.idx[[mm]]
        GLIST[[mm]][m.user.idx] <- object@Fit@se[x.user.idx]

        # class
        if(object@Model@isSymmetric[mm]) {
            class(GLIST[[mm]]) <- c("lavaan.matrix.symmetric", "matrix")
        } else {
            class(GLIST[[mm]]) <- c("lavaan.matrix", "matrix")
        }
    }

    GLIST
}

starting.values <- function(object) {

    GLIST <- object@Model@GLIST

    for(mm in 1:length(GLIST)) {
        # labels
        dimnames(GLIST[[mm]]) <- object@Model@dimNames[[mm]]

        # fill in starting values
        m.user.idx <- object@Model@m.user.idx[[mm]]
        x.user.idx <- object@Model@x.user.idx[[mm]]
        GLIST[[mm]][m.user.idx] <- object@Fit@start[x.user.idx]

        # class
        if(object@Model@isSymmetric[mm]) {
            class(GLIST[[mm]]) <- c("lavaan.matrix.symmetric", "matrix")
        } else {
            class(GLIST[[mm]]) <- c("lavaan.matrix", "matrix")
        }
    }

    GLIST
}


parameter.values <- function(object) {

    GLIST <- object@Model@GLIST

    for(mm in 1:length(GLIST)) {
        # labels
        dimnames(GLIST[[mm]]) <- object@Model@dimNames[[mm]]

        # fill in starting values
        m.user.idx <- object@Model@m.user.idx[[mm]]
        x.user.idx <- object@Model@x.user.idx[[mm]]
        GLIST[[mm]][m.user.idx] <- object@Fit@est[x.user.idx]

        # class
        if(object@Model@isSymmetric[mm]) {
            class(GLIST[[mm]]) <- c("lavaan.matrix.symmetric", "matrix")
        } else {
            class(GLIST[[mm]]) <- c("lavaan.matrix", "matrix")
        }
    }

    GLIST
}


parameter.list <- function(object) {

    LIST <- as.data.frame(object@ParTable, stringsAsFactors = FALSE)
    LIST
}


derivatives <- function(object) {

    GLIST <- lav_model_gradient(lavmodel       = object@Model,
                                GLIST          = NULL,
                                lavsamplestats = object@SampleStats,
                                lavdata        = object@Data,
                                lavcache       = object@Cache,
                                type           = "allofthem",
                                estimator      = object@Options$estimator,
                                verbose        = FALSE,
                                forcePD        = TRUE,
                                group.weight   = TRUE,
                                constraints    = FALSE,
                                Delta          = NULL)
    names(GLIST) <- names(object@Model@GLIST)

    for(mm in 1:length(GLIST)) {
        # labels
        dimnames(GLIST[[mm]]) <- object@Model@dimNames[[mm]]

        # class
        if(object@Model@isSymmetric[mm]) {
            class(GLIST[[mm]]) <- c("lavaan.matrix.symmetric", "matrix")
        } else {
            class(GLIST[[mm]]) <- c("lavaan.matrix", "matrix")
        }
    }

    GLIST
}



# fixme, should we export this function?
sampStat <- function(object, labels=TRUE) {

    G <- object@Data@ngroups
    ov.names <- object@Data@ov.names

    OUT <- vector("list", length=G)
    for(g in 1:G) {
        OUT[[g]]$cov  <- object@SampleStats@cov[[g]]
        if(labels) 
            rownames(OUT[[g]]$cov) <- colnames(OUT[[g]]$cov) <- ov.names[[g]]
        class(OUT[[g]]$cov) <- c("lavaan.matrix.symmetric", "matrix")

        #if(object@Model@meanstructure) {
            OUT[[g]]$mean <- as.numeric(object@SampleStats@mean[[g]])
            if(labels) names(OUT[[g]]$mean) <- ov.names[[g]]
            class(OUT[[g]]$mean) <- c("lavaan.vector", "numeric")
        #}

        if(object@Model@categorical) {
            OUT[[g]]$th <- as.numeric(object@SampleStats@th[[g]])
            if(length(object@Model@num.idx[[g]]) > 0L) {
                OUT[[g]]$th <- OUT[[g]]$th[-object@Model@num.idx[[g]]]
            }
            if(labels) {
                names(OUT[[g]]$th) <- 
                    vnames(object@ParTable, type="th", group=g)
            }
            class(OUT[[g]]$th) <- c("lavaan.vector", "numeric")
        }

        if(object@Model@categorical &&
           object@Model@nexo > 0L) {
            OUT[[g]]$slopes  <- object@SampleStats@slopes[[g]]
            if(labels) {
                rownames(OUT[[g]]$slopes) <- ov.names[[g]]
                colnames(OUT[[g]]$slopes) <- 
                    vnames(object@ParTable, type="ov.x", group=g)
                class(OUT[[g]]$slopes) <- c("lavaan.matrix", "matrix")
            }
        }

        if(object@Model@group.w.free) {
            OUT[[g]]$group.w <- object@SampleStats@group.w[[g]]
            if(labels) {
                names(OUT[[g]]$group.w) <- "w"
                class(OUT[[g]]$group.w) <- c("lavaan.vector", "numeric")
            }
        }
    }

    if(G == 1) {
        OUT <- OUT[[1]]
    } else {
        names(OUT) <- unlist(object@Data@group.label)
    }

    OUT
}


rsquare <- function(object, est.std.all=NULL) {

    if(is.null(est.std.all)) est.std.all <- standardize.est.all(object)
    ngroups <- object@Data@ngroups
    partable <- object@ParTable
    partable$rsquare <- 1.0 - est.std.all
    # no values > 1.0
    partable$rsquare[partable$rsquare > 1.0] <- as.numeric(NA)
    r2 <- vector("list", length=ngroups)

    for(g in 1:ngroups) {
        ind.names   <- partable$rhs[ which(partable$op == "=~" & partable$group == g) ]
        eqs.y.names <- partable$lhs[ which(partable$op == "~"  & partable$group == g) ]
        y.names <- unique( c(ind.names, eqs.y.names) )

        idx <- which(partable$op == "~~" & partable$lhs %in% y.names & 
                     partable$rhs == partable$lhs & partable$group == g)
        tmp <- partable$rsquare[idx]; names(tmp) <- partable$lhs[idx]
        r2[[g]] <- tmp
    }

    if(ngroups == 1) {
        r2 <- r2[[1]]
    } else {
        names(r2) <- unlist(object@Data@group.label)
    }

    r2
}

getWLS.est <- function(object, drop.list.single.group=FALSE) {

    G <- object@Data@ngroups
    OUT <- lav_model_wls_est(object@Model)

    if(G == 1 && drop.list.single.group) {
        OUT <- OUT[[1]]
    } else {
        if(G > 1) names(OUT) <- unlist(object@Data@group.label)
    }

    OUT
}

getWLS.V <- function(object, Delta=computeDelta(lavmodel = object@Model), 
                     drop.list.single.group=FALSE) {

    # shortcuts
    lavsamplestats = object@SampleStats
    estimator   = object@Options$estimator
    G           = object@Data@ngroups

    WLS.V       <- vector("list", length=lavsamplestats@ngroups)
    if(estimator == "GLS"  ||
       estimator == "WLS"  ||
       estimator == "DWLS" ||
       estimator == "ULS") {
        # for GLS, the WLS.V22 part is: 0.5 * t(D) %*% [S.inv %x% S.inv] %*% D
        # for WLS, the WLS.V22 part is: Gamma
        WLS.V <- lavsamplestats@WLS.V
    } else if(estimator == "ML") {
        Sigma.hat <- computeSigmaHat(lavmodel = object@Model)
        if(object@Model@meanstructure) {
            Mu.hat <- computeMuHat(lavmodel = object@Model)
        }
        for(g in seq_len(lavsamplestats@ngroups)) {
            if(lavsamplestats@missing.flag) {
                WLS.V[[g]] <- compute.Abeta(Sigma.hat=Sigma.hat[[g]],
                                            Mu.hat=Mu.hat[[g]],
                                            lavsamplestats=lavsamplestats,
                                            lavdata=object@Data, group=g,
                                            information="expected")
            } else {
                # WLS.V22 = 0.5*t(D) %*% [Sigma.hat.inv %x% Sigma.hat.inv]%*% D
                WLS.V[[g]] <-
                    compute.Abeta.complete(Sigma.hat=Sigma.hat[[g]],
                                           meanstructure=object@Model@meanstructure)
            }
        }
    } # ML


    OUT <- WLS.V
    if(G == 1 && drop.list.single.group) {
        OUT <- OUT[[1]]
    } else {
        if(G > 1) names(OUT) <- unlist(object@Data@group.label)
    }

    OUT
}

getSampleStatsNACOV <- function(object) {

    if(object@Options$se == "robust.mlr")
        stop("not done yet; FIX THIS!")

    # shortcuts
    lavsamplestats = object@SampleStats
    estimator   = object@Options$estimator

    NACOV       <- vector("list", length=lavsamplestats@ngroups)

    if(estimator == "GLS"  ||
       estimator == "WLS"  ||
       estimator == "DWLS" ||
       estimator == "ULS") {
        NACOV <- lavsamplestats@NACOV
    } else if(estimator == "ML") {
        for(g in 1:lavsamplestats@ngroups) {
            NACOV[[g]] <- 
                compute.Gamma(object@Data@X[[g]], 
                              meanstructure=object@Options$meanstructure)
            if(object@Options$mimic == "Mplus") {
                G11 <- ( lavsamplestats@cov[[g]] * (lavsamplestats@nobs[[g]]-1)/lavsamplestats@nobs[[g]] )
                NACOV[[g]][1:nrow(G11), 1:nrow(G11)] <- G11
            }
        }
    }

    NACOV
}

getHessian <- function(object) {
    # lazy approach: take -1 the observed information
    E <- computeObservedInformation(lavmodel       = object@Model, 
                                    lavsamplestats = object@SampleStats,
                                    lavdata        = object@Data,
                                    lavcache       = object@Cache,
                                    type           = "free",
                                    estimator      = object@Options$estimator,
                                    group.weight   = TRUE)

    -E
}

getVariability <- function(object) {
    # lazy approach: get it from Nvcov.first.order
    NACOV <- Nvcov.first.order(lavmodel       = object@Model,
                               lavsamplestats = object@SampleStats,
                               lavdata        = object@Data,
                               estimator      = object@Options$estimator)

    B0 <- attr(NACOV, "B0")

    if(object@Options$estimator == "PML") {
        B0 <- B0 * object@SampleStats@ntotal
    }
    
    B0
}

getModelCorLV <- function(object, labels=TRUE) {
    getModelCovLV(object, correlation.metric=TRUE, labels=labels)
}

getModelCovLV <- function(object, correlation.metric=FALSE, labels=TRUE) {

    G <- object@Data@ngroups

    # compute lv covar
    OUT <- computeVETA(lavmodel       = object@Model, 
                       lavsamplestats = object@SampleStats)

    # correlation?
    if(correlation.metric) {
        OUT <- lapply(OUT, cov2cor)
    }
    
    # we need psi matrix for labels
    psi.group <- which(names(object@Model@GLIST) == "psi")
    if(labels) {
        for(g in 1:G) {
            psi.idx <- psi.group[g]
            NAMES <- object@Model@dimNames[[psi.idx]][[1L]]
            # remove all dummy latent variables
            lv.idx <- c(object@Model@ov.y.dummy.lv.idx[[g]],
                        object@Model@ov.x.dummy.lv.idx[[g]])
            if(length(lv.idx) > 0L) {
                NAMES <- NAMES[-lv.idx]
                OUT[[g]] <- OUT[[g]][-lv.idx, -lv.idx, drop=FALSE]
            }
            colnames(OUT[[g]]) <- rownames(OUT[[g]]) <- NAMES
            class(OUT[[g]]) <- c("lavaan.matrix.symmetric", "matrix")
        }
    }

    if(G == 1) {
        OUT <- OUT[[1]]
    } else {
        names(OUT) <- unlist(object@Data@group.label)
    }

    OUT
}

getModelCor <- function(object, labels=TRUE) {
    getModelCov(object, correlation.metric=TRUE, labels=labels)
}

getModelMeanLV <- function(object, labels=TRUE) {

    G <- object@Data@ngroups

    # compute lv means
    OUT <- computeEETA(lavmodel       = object@Model, 
                       lavsamplestats = object@SampleStats,
                       remove.dummy.lv = TRUE)

    OUT <- lapply(OUT, as.numeric)

    # we need psi matrix for labels
    psi.group <- which(names(object@Model@GLIST) == "psi")
    if(labels) {
        for(g in 1:G) {
            psi.idx <- psi.group[g]
            NAMES <- object@Model@dimNames[[psi.idx]][[1L]]
            # remove all dummy latent variables
            lv.idx <- c(object@Model@ov.y.dummy.lv.idx[[g]],
                        object@Model@ov.x.dummy.lv.idx[[g]])
            if(length(lv.idx) > 0L) {
                NAMES <- NAMES[-lv.idx]
                OUT[[g]] <- OUT[[g]][-lv.idx]
            }
            names(OUT[[g]]) <- NAMES
            class(OUT[[g]]) <- c("lavaan.vector", "numeric")
        }
    }

    if(G == 1) {
        OUT <- OUT[[1]]
    } else {
        names(OUT) <- unlist(object@Data@group.label)
    }

    OUT
}

getModelCov <- function(object, correlation.metric=FALSE, labels=TRUE) {

    G <- object@Data@ngroups

    # compute extended model implied covariance matrix (both ov and lv)
    OUT <- computeCOV(lavmodel = object@Model, 
                      lavsamplestats = object@SampleStats)

    # correlation?
    if(correlation.metric) {
        OUT <- lapply(OUT, cov2cor)
    }

    # we need lambda + psi matrix for labels
    lambda.group <- which(names(object@Model@GLIST) == "lambda")
    psi.group <- which(names(object@Model@GLIST) == "psi")
    if(labels) {
        for(g in 1:G) {
            lambda.idx <- lambda.group[g]
            NAMES.ov <- object@Model@dimNames[[lambda.idx]][[1L]]

            psi.idx <- psi.group[g]
            NAMES.lv <- object@Model@dimNames[[psi.idx]][[1L]]
            # remove all dummy latent variables
            lv.idx <- c(object@Model@ov.y.dummy.lv.idx[[g]],
                        object@Model@ov.x.dummy.lv.idx[[g]])
            if(length(lv.idx) > 0L) {
                NAMES.lv <- NAMES.lv[-lv.idx]
                lv.idx <- lv.idx + length(NAMES.ov)
                OUT[[g]] <- OUT[[g]][-lv.idx, -lv.idx, drop=FALSE]
            }

            NAMES <- c(NAMES.ov, NAMES.lv)
            colnames(OUT[[g]]) <- rownames(OUT[[g]]) <- NAMES
            class(OUT[[g]]) <- c("lavaan.matrix.symmetric", "matrix")
        }
    }

    if(G == 1) {
        OUT <- OUT[[1]]
    } else {
        names(OUT) <- unlist(object@Data@group.label)
    }

    OUT
}


getModelTheta <- function(object, correlation.metric=FALSE, labels=TRUE) {

    G <- object@Data@ngroups

    # get residual covariances
    OUT <- computeTHETA(lavmodel = object@Model)

    # correlation?
    if(correlation.metric) {
        OUT <- lapply(OUT, cov2cor)
    }

    # we need lambda for labels
    lambda.group <- which(names(object@Model@GLIST) == "lambda")
    if(labels) {
        for(g in 1:G) {
            lambda.idx <- lambda.group[g]
            NAMES <- object@Model@dimNames[[lambda.idx]][[1L]]
            colnames(OUT[[g]]) <- rownames(OUT[[g]]) <- NAMES
            class(OUT[[g]]) <- c("lavaan.matrix.symmetric", "matrix")
        }
    }

    if(G == 1) {
        OUT <- OUT[[1]]
    } else {
        names(OUT) <- unlist(object@Data@group.label)
    }

    OUT
}

getDataMissingCoverage <- function(object, labels=TRUE) {

    G <- object@Data@ngroups

    # get missing covarage
    OUT <- vector("list", G)
   
    for(g in 1:G) {
        if(!is.null(object@Data@Mp[[g]])) {
            OUT[[g]] <- object@Data@Mp[[g]]$coverage
        } else {
            nvar <- length(object@Data@ov.names[[g]])
            OUT[[g]] <- matrix(1.0, nvar, nvar)
        }

        if(labels) {
            NAMES <- object@Data@ov.names[[g]]
            colnames(OUT[[g]]) <- rownames(OUT[[g]]) <- NAMES
            class(OUT[[g]]) <- c("lavaan.matrix.symmetric", "matrix")
        }
    }

    if(G == 1) {
        OUT <- OUT[[1]]
    } else {
        names(OUT) <- unlist(object@Data@group.label)
    }

    OUT
}

getDataMissingPatterns <- function(object, labels = TRUE) {

    G <- object@Data@ngroups

    # get missing covarage
    OUT <- vector("list", G)
   
    for(g in 1:G) {
        if(!is.null(object@Data@Mp[[g]])) {
            OUT[[g]] <- object@Data@Mp[[g]]$pat
        } else {
            nvar <- length(object@Data@ov.names[[g]])
            OUT[[g]] <- matrix(TRUE, 1L, nvar)
            rownames(OUT[[g]]) <- object@Data@nobs[[g]]
        }

        if(labels) {
            NAMES <- object@Data@ov.names[[g]]
            colnames(OUT[[g]]) <- NAMES
            class(OUT[[g]]) <- c("lavaan.matrix", "matrix")
        }
    }

    if(G == 1) {
        OUT <- OUT[[1]]
    } else {
        names(OUT) <- unlist(object@Data@group.label)
    }

    OUT
}


lav_inspect_data <- function(object, labels = TRUE) {

    G <- object@Data@ngroups
    OUT <- object@Data@X

    if(labels) {
        for(g in 1:G) {
            colnames(OUT[[g]]) <- object@Data@ov.names[[g]]
        }
    }

    if(G == 1) {
        OUT <- OUT[[1]]
    } else {
        names(OUT) <- unlist(object@Data@group.label)
    }

    OUT
}

