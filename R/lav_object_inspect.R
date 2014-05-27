# inspect a fitted lavaan object

# backward compatibility -- wrapper around lavInspect
inspect <- function(object, what = "free") {
    lavInspect(object = object,
               what = what,
               add.labels = TRUE,
               add.class = TRUE,
               drop.list.single.group = TRUE)
}

# the `tech' version: no labels, full matrices, ... for further processing
lavTech <- function(object, 
                    what = "free", 
                    add.labels = FALSE, 
                    add.class = FALSE,
                    drop.list.single.group = TRUE) { # or FALSE?
    lavInspect(object = object, what = what,
               add.labels = add.labels, add.class = add.class,
                drop.list.single.group =  drop.list.single.group)
}

# the `user' version: for display only
lavInspect <- function(object,
                       what = "free",
                       add.labels = TRUE,
                       add.class = TRUE,
                       drop.list.single.group = TRUE) {

    # object must inherit from class lavaan
    stopifnot(inherits(object, "lavaan"))

    # only a single argument
    if(length(what) > 1) {
        stop("`what' arguments contains multiple arguments; only one is allowed")
    }

    # be case insensitive
    what <- tolower(what)


    #### model matrices, with different contents ####
    if(what == "free") {
        lav_object_inspect_modelmatrices(object, what = "free", type = "free",
            add.labels = add.labels, add.class = add.class)
    } else if(what == "unco") {
        lav_object_inspect_modelmatrices(object, what = "free", type="unco",
            add.labels = add.labels, add.class = add.class)
    } else if(what == "user") {
        lav_object_inspect_modelmatrices(object, what = "free", type="user",
            add.labels = add.labels, add.class = add.class)
    } else if(what == "se" ||
              what == "std.err" ||
              what == "standard.errors") {
        lav_object_inspect_modelmatrices(object, what = "se",
            add.labels = add.labels, add.class = add.class)
    } else if(what == "start" || what == "starting.values") {
        lav_object_inspect_modelmatrices(object, what = "start",
            add.labels = add.labels, add.class = add.class)
    } else if(what == "est"  || what == "estimates" ||
              what == "coef" || what == "coefficients" ||
              what == "x") {
        lav_object_inspect_modelmatrices(object, what = "est",
            add.labels = add.labels, add.class = add.class)
    } else if(what == "dx.free") {
        lav_object_inspect_modelmatrices(object, what = "dx.free",
            add.labels = add.labels, add.class = add.class)
    } else if(what == "dx" || what == "gradient" || what == "derivatives") {
        lav_object_inspect_modelmatrices(object, what = "dx",
            add.labels = add.labels, add.class = add.class)
    } else if(what == "std" || what == "std.all") {
        lav_object_inspect_modelmatrices(object, what = "std.all",
           add.labels = add.labels, add.class = add.class)
    } else if(what == "std.lv") {
        lav_object_inspect_modelmatrices(object, what = "std.lv",
           add.labels = add.labels, add.class = add.class)
    } else if(what == "std.nox") {
        lav_object_inspect_modelmatrices(object, what = "std.nox",
           add.labels = add.labels, add.class = add.class)


    #### parameter table ####
    } else if(what == "list") {
        parTable(object)

    #### fit indices ####
    } else if(what == "fit" ||
              what == "fitmeasures" ||
              what == "fit.measures" ||
              what == "fit.indices") {
        fitMeasures(object)


    #### modification indices ####
    } else if(what == "mi" ||
              what == "modindices" ||
              what == "modification.indices") {
        modificationIndices(object)


    #### sample statistics #####
    } else if(what == "sampstat" ||
              what == "samp" ||
              what == "sample" ||
              what == "samplestatistics") {
        lav_object_inspect_sampstat(object, 
            add.labels = add.labels, add.class = add.class, 
            drop.list.single.group = drop.list.single.group)


    #### data ####
    } else if(what == "data") {
        lav_object_inspect_data(object, 
            drop.list.single.group = drop.list.single.group)


    #### rsquare ####
    } else if(what == "rsquare" || what == "r-square" || what == "r2") {
         lav_object_inspect_rsquare(object, 
             add.labels = add.labels, add.class = add.class,
             drop.list.single.group = drop.list.single.group)


    #### model-implied sample statistics ####
    } else if(what == "cov.lv") {
        lav_object_inspect_cov_lv(object,
            correlation.metric = FALSE,
            add.labels = add.labels, add.class = add.class,
            drop.list.single.group = drop.list.single.group)
    } else if(what == "cor.lv") {
        lav_object_inspect_cov_lv(object,
            correlation.metric = TRUE,
            add.labels = add.labels, add.class = add.class,
            drop.list.single.group = drop.list.single.group)
    } else if(what == "mean.lv") {
        lav_object_inspect_mean_lv(object,
            add.labels = add.labels, add.class = add.class,
            drop.list.single.group = drop.list.single.group)
    } else if(what == "cov.all") {
        lav_object_inspect_cov_all(object,
            correlation.metric = FALSE,
            add.labels = add.labels, add.class = add.class,
            drop.list.single.group = drop.list.single.group)
    } else if(what == "cor.all") {
        lav_object_inspect_cov_all(object,
            correlation.metric = TRUE,
            add.labels = add.labels, add.class = add.class,
            drop.list.single.group = drop.list.single.group)
    } else if(what == "cov.ov") {
        lav_object_inspect_cov_ov(object,
            correlation.metric = FALSE,
            add.labels = add.labels, add.class = add.class,
            drop.list.single.group = drop.list.single.group)
    } else if(what == "cor.ov") {
        lav_object_inspect_cov_lv(object,
            correlation.metric = TRUE,
            add.labels = add.labels, add.class = add.class,
            drop.list.single.group = drop.list.single.group)
    } else if(what == "mean.ov") {
        lav_object_inspect_mean_ov(object,
            add.labels = add.labels, add.class = add.class,
            drop.list.single.group = drop.list.single.group)
    } else if(what == "th") {
        lav_object_inspect_th(object,
            add.labels = add.labels, add.class = add.class,
            drop.list.single.group = drop.list.single.group)


    #### wls.v - nacov####
    } else if(what == "wls.v") {
        getWLS.V(object, drop.list.single.group = drop.list.single.group)
    } else if(what == "nacov") {
        getSampleStatsNACOV(object)



    #### specific model matrices? ####
    } else if(what == "theta") {
        getModelTheta(object, labels=TRUE,
                      drop.list.single.group = drop.list.single.group)


    #### missingness ####
    } else if(what == "coverage") {
        getDataMissingCoverage(object, labels=TRUE,
                               drop.list.single.group = drop.list.single.group)
    } else if(what %in% c("patterns", "pattern")) {
        getDataMissingPatterns(object, labels=TRUE,
                               drop.list.single.group = drop.list.single.group)


    #### convergence ####
    } else if(what == "converged") {
        object@Fit@converged

    #### not found ####
    } else {
        stop("unknown `what' argument in inspect function: `", what, "'")
    }

}


lav_object_inspect_modelmatrices <- function(object, what = "free",
    type = "free", add.labels = FALSE, add.class = FALSE) {

    GLIST <- object@Model@GLIST

    if(what == "dx.free") {
        DX <- lav_model_gradient(lavmodel       = object@Model,
                                 GLIST          = NULL,
                                 lavsamplestats = object@SampleStats,
                                 lavdata        = object@Data,
                                 lavcache       = object@Cache,
                                 type           = "free",
                                 estimator      = object@Options$estimator,
                                 verbose        = FALSE,
                                 forcePD        = TRUE,
                                 group.weight   = TRUE,
                                 constraints    = FALSE,
                                 Delta          = NULL)
    } else if(what == "dx") {
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
    } else if(what == "std.all") {
        STD <- standardize.est.all(object)
    } else if(what == "std.lv") {
        STD <- standardize.est.lv(object)
    } else if(what == "std.nox") {
        STD <- standardize.est.all.nox(object)
    }

    for(mm in 1:length(GLIST)) {

        if(what != "dx") {
            # erase everything
            GLIST[[mm]][,] <- 0.0
        }

        if(add.labels) {
            dimnames(GLIST[[mm]]) <- object@Model@dimNames[[mm]]
        }

        if(what == "free") {
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
        } else if(what == "se") {
            # fill in standard errors
            m.user.idx <- object@Model@m.user.idx[[mm]]
            x.user.idx <- object@Model@x.user.idx[[mm]]
            GLIST[[mm]][m.user.idx] <- object@Fit@se[x.user.idx]
        } else if(what == "start") {
            # fill in starting values
            m.user.idx <- object@Model@m.user.idx[[mm]]
            x.user.idx <- object@Model@x.user.idx[[mm]]
            GLIST[[mm]][m.user.idx] <- object@Fit@start[x.user.idx]
        } else if(what == "est") {
            # fill in estimated parameter values
            m.user.idx <- object@Model@m.user.idx[[mm]]
            x.user.idx <- object@Model@x.user.idx[[mm]]
            GLIST[[mm]][m.user.idx] <- object@Fit@est[x.user.idx]
        } else if(what == "dx.free") {
            # fill in derivatives free parameters
            m.el.idx <- object@Model@m.free.idx[[mm]]
            x.el.idx <- object@Model@x.free.idx[[mm]]
            GLIST[[mm]][m.el.idx] <- DX[x.el.idx]
        } else if(what %in% c("std.all", "std.lv", "std.nox")) {
            m.user.idx <- object@Model@m.user.idx[[mm]]
            x.user.idx <- object@Model@x.user.idx[[mm]]
            GLIST[[mm]][m.user.idx] <- STD[x.user.idx]
        } 

        # class
        if(add.class) {
            if(object@Model@isSymmetric[mm]) {
                class(GLIST[[mm]]) <- c("lavaan.matrix.symmetric", "matrix")
            } else {
                class(GLIST[[mm]]) <- c("lavaan.matrix", "matrix")
            }
        }
    }

    GLIST
}




# fixme, should we export this function?
lav_object_inspect_sampstat <- function(object, 
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {

    G <- object@Data@ngroups
    ov.names <- object@Data@ov.names

    OUT <- vector("list", length=G)
    for(g in 1:G) {
        OUT[[g]]$cov  <- object@SampleStats@cov[[g]]
        if(add.labels) {
            rownames(OUT[[g]]$cov) <- colnames(OUT[[g]]$cov) <- ov.names[[g]]
        }
        if(add.class) {
            class(OUT[[g]]$cov) <- c("lavaan.matrix.symmetric", "matrix")
        }

        #if(object@Model@meanstructure) {
            OUT[[g]]$mean <- as.numeric(object@SampleStats@mean[[g]])
            if(add.labels) {
                names(OUT[[g]]$mean) <- ov.names[[g]]
            }
            if(add.class) {
                class(OUT[[g]]$mean) <- c("lavaan.vector", "numeric")
            }
        #}

        if(object@Model@categorical) {
            OUT[[g]]$th <- as.numeric(object@SampleStats@th[[g]])
            if(length(object@Model@num.idx[[g]]) > 0L) {
                OUT[[g]]$th <- OUT[[g]]$th[-object@Model@num.idx[[g]]]
            }
            if(add.labels) {
                names(OUT[[g]]$th) <- 
                    vnames(object@ParTable, type="th", group=g)
            }
            if(add.class) {
                class(OUT[[g]]$th) <- c("lavaan.vector", "numeric")
            }
        }

        if(object@Model@categorical &&
           object@Model@nexo > 0L) {
            OUT[[g]]$slopes  <- object@SampleStats@slopes[[g]]
            if(add.labels) {
                rownames(OUT[[g]]$slopes) <- ov.names[[g]]
                colnames(OUT[[g]]$slopes) <- 
                    vnames(object@ParTable, type="ov.x", group=g)
            }
            if(add.class) {
                class(OUT[[g]]$slopes) <- c("lavaan.matrix", "matrix")
            }
        }

        if(object@Model@group.w.free) {
            OUT[[g]]$group.w <- object@SampleStats@group.w[[g]]
            if(add.labels) {
                names(OUT[[g]]$group.w) <- "w"
            }
            if(add.class) {
                class(OUT[[g]]$group.w) <- c("lavaan.vector", "numeric")
            }
        }
    }

    if(G == 1L && drop.list.single.group) {
        OUT <- OUT[[1]]
    } else {
        if(length(object@Data@group.label) > 0L) {
            names(OUT) <- unlist(object@Data@group.label)
        }
    }

    OUT
}


lav_object_inspect_data <- function(object, add.labels = FALSE,
                                    drop.list.single.group = FALSE) {

    G <- object@Data@ngroups
    OUT <- object@Data@X

    if(add.labels) {
        for(g in 1:G) {
            colnames(OUT[[g]]) <- object@Data@ov.names[[g]]
        }
    }

    if(G == 1L && drop.list.single.group) {
        OUT <- OUT[[1]]
    } else {
        if(length(object@Data@group.label) > 0L) {
            names(OUT) <- unlist(object@Data@group.label)
        }
    }

    OUT
}


lav_object_inspect_rsquare <- function(object, est.std.all=NULL,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {

    G <- object@Data@ngroups
    OUT <- vector("list", length=G)

    if(is.null(est.std.all)) {
        est.std.all <- standardize.est.all(object)
    }

    partable <- object@ParTable
    partable$rsquare <- 1.0 - est.std.all
    # no values > 1.0
    partable$rsquare[partable$rsquare > 1.0] <- as.numeric(NA)

    for(g in 1:G) {
        ind.names <- partable$rhs[ which(partable$op == "=~" & 
                                         partable$group == g) ]
        eqs.y.names <- partable$lhs[ which(partable$op == "~"  & 
                                           partable$group == g) ]
        y.names <- unique( c(ind.names, eqs.y.names) )

        idx <- which(partable$op == "~~" & partable$lhs %in% y.names & 
                     partable$rhs == partable$lhs & partable$group == g)
        tmp <- partable$rsquare[idx]

        if(add.labels) {
            names(tmp) <- partable$lhs[idx]
        }
        if(add.class) {
            class(tmp) <- c("lavaan.vector", "numeric")
        }

        OUT[[g]] <- tmp
    }

    if(G == 1L && drop.list.single.group) {
        OUT <- OUT[[1]]
    } else {
        if(length(object@Data@group.label) > 0L) {
            names(OUT) <- unlist(object@Data@group.label)
        }
    }

    OUT 
}


lav_object_inspect_cov_lv <- function(object, correlation.metric = FALSE,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {

    G <- object@Data@ngroups

    # compute lv covar
    OUT <- computeVETA(lavmodel       = object@Model, 
                       lavsamplestats = object@SampleStats,
                       remove.dummy.lv = TRUE)

    # cor + labels + class
    for(g in 1:G) {

        if(correlation.metric && nrow(OUT[[g]]) > 1L) {
            # note: cov2cor fails if matrix is empty!
            OUT[[g]] <- cov2cor(OUT[[g]])
        }

        if(add.labels) {
            colnames(OUT[[g]]) <- rownames(OUT[[g]]) <- 
                object@pta$vnames$lv.regular[[g]]
        }

        if(add.class) {
            class(OUT[[g]]) <- c("lavaan.matrix.symmetric", "matrix")
        }
    }

    if(G == 1L && drop.list.single.group) {
        OUT <- OUT[[1]]
    } else {
        if(length(object@Data@group.label) > 0L) {
            names(OUT) <- unlist(object@Data@group.label)
        }
    }

    OUT
}

lav_object_inspect_mean_lv <- function(object,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {

    G <- object@Data@ngroups

    # compute lv means
    OUT <- computeEETA(lavmodel       = object@Model, 
                       lavsamplestats = object@SampleStats,
                       remove.dummy.lv = TRUE)
    OUT <- lapply(OUT, as.numeric)

    # labels + class
    for(g in 1:G) {
        if(add.labels) {
            names(OUT[[g]]) <- object@pta$vnames$lv.regular[[g]]
        }
        if(add.class) {
            class(OUT[[g]]) <- c("lavaan.vector", "numeric")
        }
    }

    if(G == 1L && drop.list.single.group) {
        OUT <- OUT[[1]]
    } else {
        if(length(object@Data@group.label) > 0L) {
            names(OUT) <- unlist(object@Data@group.label)
        }
    }

    OUT
}

lav_object_inspect_cov_all <- function(object, correlation.metric = FALSE,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {

    G <- object@Data@ngroups

    # compute extended model implied covariance matrix (both ov and lv)
    OUT <- computeCOV(lavmodel = object@Model, 
                      lavsamplestats = object@SampleStats)

    # cor + labels + class
    for(g in 1:G) {

        if(correlation.metric && nrow(OUT[[g]]) > 1L) {
            # note: cov2cor fails if matrix is empty!
            OUT[[g]] <- cov2cor(OUT[[g]])
        }

        if(add.labels) {
            NAMES <- c(object@pta$vnames$ov[[g]],
                       object@pta$vnames$lv.regular[[g]])
            colnames(OUT[[g]]) <- rownames(OUT[[g]]) <- NAMES
        }
        if(add.class) {
            class(OUT[[g]]) <- c("lavaan.matrix.symmetric", "matrix")
        }
    }

    if(G == 1L && drop.list.single.group) {
        OUT <- OUT[[1]]
    } else {
        if(length(object@Data@group.label) > 0L) {
            names(OUT) <- unlist(object@Data@group.label)
        }
    }

    OUT
}


lav_object_inspect_cov_ov <- function(object, correlation.metric = FALSE,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {

    G <- object@Data@ngroups

    # get model-implied covariance matrix observed
    OUT <- object@Fit@Sigma.hat

    # cor + labels + class
    for(g in 1:G) {
 
        if(correlation.metric && nrow(OUT[[g]]) > 1L) {
            # note: cov2cor fails if matrix is empty!
            OUT[[g]] <- cov2cor(OUT[[g]])
        }

        if(add.labels) {
            colnames(OUT[[g]]) <- rownames(OUT[[g]]) <- 
                object@pta$vnames$ov[[g]]
        }
        if(add.class) {
            class(OUT[[g]]) <- c("lavaan.matrix.symmetric", "matrix")
        }
    }

    if(G == 1L && drop.list.single.group) {
        OUT <- OUT[[1]]
    } else {
        if(length(object@Data@group.label) > 0L) {
            names(OUT) <- unlist(object@Data@group.label)
        }
    }

    OUT
}

lav_object_inspect_mean_ov <- function(object,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {

    G <- object@Data@ngroups

    # compute lv means
    OUT <- object@Fit@Mu.hat
    OUT <- lapply(OUT, as.numeric)

    # labels + class
    for(g in 1:G) {
        if(add.labels) {
            names(OUT[[g]]) <- object@pta$vnames$ov[[g]]
        }
        if(add.class) {
            class(OUT[[g]]) <- c("lavaan.vector", "numeric")
        }
    }

    if(G == 1L && drop.list.single.group) {
        OUT <- OUT[[1]]
    } else {
        if(length(object@Data@group.label) > 0L) {
            names(OUT) <- unlist(object@Data@group.label)
        }
    }

    OUT
}

lav_object_inspect_th <- function(object,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {

    G <- object@Data@ngroups

    # thresholds
    OUT <- object@Fit@TH
    OUT <- lapply(OUT, as.numeric)

    # labels + class
    for(g in 1:G) {
        if(length(object@Model@num.idx[[g]]) > 0L) {
            OUT[[g]] <- OUT[[g]][ -object@Model@num.idx[[g]] ]
        }
        if(add.labels) {
            names(OUT[[g]]) <- object@pta$vnames$th[[g]]
        }
        if(add.class) {
            class(OUT[[g]]) <- c("lavaan.vector", "numeric")
        }
    }

    if(G == 1L && drop.list.single.group) {
        OUT <- OUT[[1]]
    } else {
        if(length(object@Data@group.label) > 0L) {
            names(OUT) <- unlist(object@Data@group.label)
        }
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

