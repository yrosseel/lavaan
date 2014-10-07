# inspect a fitted lavaan object

# backward compatibility -- wrapper around lavInspect
setMethod("inspect", "lavaan",
function(object, what = "free") {
    lavInspect(lavobject              = object,
               what                   = what,
               add.labels             = TRUE,
               add.class              = TRUE,
               drop.list.single.group = TRUE)
})

# the `tech' version: no labels, full matrices, ... for further processing
lavTech <- function(lavobject, 
                    what                   = "free",
                    add.labels             = FALSE,
                    add.class              = FALSE,
                    drop.list.single.group = FALSE) {

    lavInspect(lavobject = lavobject, what = what,
               add.labels = add.labels, add.class = add.class,
               drop.list.single.group =  drop.list.single.group)
}

# the `user' version: with defaults for display only
lavInspect <- function(lavobject,
                       what                   = "free",
                       add.labels             = TRUE,
                       add.class              = TRUE,
                       drop.list.single.group = TRUE) {

    # lavobject must inherit from class lavaan
    stopifnot(inherits(lavobject, "lavaan"))

    # only a single argument
    if(length(what) > 1) {
        stop("`what' arguments contains multiple arguments; only one is allowed")
    }

    # be case insensitive
    what <- tolower(what)


    #### model matrices, with different contents ####
    if(what == "free") {
        lav_object_inspect_modelmatrices(lavobject, what = "free", type = "free",
            add.labels = add.labels, add.class = add.class)
    } else if(what == "unco" || what == "unconstrained") {
        lav_object_inspect_modelmatrices(lavobject, what = "free", type="unco",
            add.labels = add.labels, add.class = add.class)
    } else if(what == "partable" || what == "user") {
        lav_object_inspect_modelmatrices(lavobject, what = "free", type="partable",
            add.labels = add.labels, add.class = add.class)
    } else if(what == "se" ||
              what == "std.err" ||
              what == "standard.errors") {
        lav_object_inspect_modelmatrices(lavobject, what = "se",
            add.labels = add.labels, add.class = add.class)
    } else if(what == "start" || what == "starting.values") {
        lav_object_inspect_modelmatrices(lavobject, what = "start",
            add.labels = add.labels, add.class = add.class)
    } else if(what == "est"  || what == "estimates" ||
              what == "coef" || what == "coefficients" ||
              what == "x") {
        lav_object_inspect_modelmatrices(lavobject, what = "est",
            add.labels = add.labels, add.class = add.class)
    } else if(what == "dx.free") {
        lav_object_inspect_modelmatrices(lavobject, what = "dx.free",
            add.labels = add.labels, add.class = add.class)
    } else if(what == "dx.all") {
        lav_object_inspect_modelmatrices(lavobject, what = "dx.all",
            add.labels = add.labels, add.class = add.class)
    } else if(what == "std" || what == "std.all" || what == "standardized") {
        lav_object_inspect_modelmatrices(lavobject, what = "std.all",
           add.labels = add.labels, add.class = add.class)
    } else if(what == "std.lv") {
        lav_object_inspect_modelmatrices(lavobject, what = "std.lv",
           add.labels = add.labels, add.class = add.class)
    } else if(what == "std.nox") {
        lav_object_inspect_modelmatrices(lavobject, what = "std.nox",
           add.labels = add.labels, add.class = add.class)


    #### parameter table ####
    } else if(what == "list") {
        parTable(lavobject)

    #### fit indices ####
    } else if(what == "fit" ||
              what == "fitmeasures" ||
              what == "fit.measures" ||
              what == "fit.indices") {
        fitMeasures(lavobject)


    #### modification indices ####
    } else if(what == "mi" ||
              what == "modindices" ||
              what == "modification.indices") {
        modificationIndices(lavobject)


    #### sample statistics #####
    } else if(what == "sampstat" ||
              what == "samp" ||
              what == "sample" ||
              what == "samplestatistics") {
        lav_object_inspect_sampstat(lavobject, 
            add.labels = add.labels, add.class = add.class, 
            drop.list.single.group = drop.list.single.group)


    #### data ####
    } else if(what == "data") {
        lav_object_inspect_data(lavobject, 
            drop.list.single.group = drop.list.single.group)
    } else if(what == "case.idx") {
        lav_object_inspect_case_idx(lavobject,
            drop.list.single.group = drop.list.single.group) 


    #### rsquare ####
    } else if(what == "rsquare" || what == "r-square" || what == "r2") {
         lav_object_inspect_rsquare(lavobject, 
             add.labels = add.labels, add.class = add.class,
             drop.list.single.group = drop.list.single.group)


    #### model-implied sample statistics ####
    } else if(what == "cov.lv" || what == "veta") {
        lav_object_inspect_cov_lv(lavobject,
            correlation.metric = FALSE,
            add.labels = add.labels, add.class = add.class,
            drop.list.single.group = drop.list.single.group)
    } else if(what == "cor.lv") {
        lav_object_inspect_cov_lv(lavobject,
            correlation.metric = TRUE,
            add.labels = add.labels, add.class = add.class,
            drop.list.single.group = drop.list.single.group)
    } else if(what == "mean.lv" || what == "eeta") {
        lav_object_inspect_mean_lv(lavobject,
            add.labels = add.labels, add.class = add.class,
            drop.list.single.group = drop.list.single.group)
    } else if(what == "cov.all") {
        lav_object_inspect_cov_all(lavobject,
            correlation.metric = FALSE,
            add.labels = add.labels, add.class = add.class,
            drop.list.single.group = drop.list.single.group)
    } else if(what == "cor.all") {
        lav_object_inspect_cov_all(lavobject,
            correlation.metric = TRUE,
            add.labels = add.labels, add.class = add.class,
            drop.list.single.group = drop.list.single.group)
    } else if(what == "cov.ov" || what == "sigma" || what == "sigma.hat") {
        lav_object_inspect_cov_ov(lavobject,
            correlation.metric = FALSE,
            add.labels = add.labels, add.class = add.class,
            drop.list.single.group = drop.list.single.group)
    } else if(what == "cor.ov") {
        lav_object_inspect_cov_ov(lavobject,
            correlation.metric = TRUE,
            add.labels = add.labels, add.class = add.class,
            drop.list.single.group = drop.list.single.group)
    } else if(what == "mean.ov" || what == "mu" || what == "mu.hat") {
        lav_object_inspect_mean_ov(lavobject,
            add.labels = add.labels, add.class = add.class,
            drop.list.single.group = drop.list.single.group)
    } else if(what == "th" || what == "thresholds") {
        lav_object_inspect_th(lavobject,
            add.labels = add.labels, add.class = add.class,
            drop.list.single.group = drop.list.single.group)
    } else if(what == "vy") {
        lav_object_inspect_vy(lavobject,
            add.labels = add.labels, add.class = add.class,
            drop.list.single.group = drop.list.single.group)


    #### specific model matrices? ####
    } else if(what == "theta" || what == "theta.cov") {
        lav_object_inspect_theta(lavobject,  correlation.metric = FALSE,
            add.labels = add.labels, add.class = add.class,
            drop.list.single.group = drop.list.single.group)
    } else if(what == "theta.cor") {
        lav_object_inspect_theta(lavobject,  correlation.metric = TRUE,
            add.labels = add.labels, add.class = add.class,
            drop.list.single.group = drop.list.single.group)
 


    #### missingness ####
    } else if(what == "coverage") {
        lav_object_inspect_missing_coverage(lavobject,
            add.labels = add.labels, add.class = add.class,
            drop.list.single.group = drop.list.single.group)
    } else if(what %in% c("patterns", "pattern")) {
        lav_object_inspect_missing_patterns(lavobject,
            add.labels = add.labels, add.class = add.class,
            drop.list.single.group = drop.list.single.group)
    } else if(what == "empty.idx") {
        lav_object_inspect_empty_idx(lavobject,
            drop.list.single.group = drop.list.single.group)

    #### convergence, meanstructure, categorical ####
    } else if(what == "converged") {
        lavobject@Fit@converged
    } else if(what == "meanstructure") {
        lavobject@Model@meanstructure
    } else if(what == "categorical") {
        lavobject@Model@categorical
    


    #### wls.est - wls.obs - wls.v ####
    } else if(what == "wls.est") {
        lav_object_inspect_wls_est(lavobject,
            add.labels = add.labels, add.class = add.class,
            drop.list.single.group = drop.list.single.group)
    } else if(what == "wls.obs") {
        lav_object_inspect_wls_obs(lavobject,
            add.labels = add.labels, add.class = add.class,
            drop.list.single.group = drop.list.single.group)
    } else if(what == "wls.v") {
        lav_object_inspect_wls_v(lavobject,
            add.labels = add.labels, add.class = add.class,
            drop.list.single.group = drop.list.single.group)


    #### NACOV samplestats ####
    } else if(what == "gamma" || what == "sampstat.nacov") {
        lav_object_inspect_sampstat_nacov(lavobject,
            add.labels = add.labels, add.class = add.class,
            drop.list.single.group = drop.list.single.group)


    #### Hessian, first.order ####
    } else if(what == "hessian") {
        lav_object_inspect_hessian(lavobject,
            add.labels = add.labels, add.class = add.class)
    } else if(what == "first.order") {
        lav_object_inspect_firstorder(lavobject,
            add.labels = add.labels, add.class = add.class)

    




    #### not found ####
    } else {
        stop("unknown `what' argument in inspect function: `", what, "'")
    }

}


lav_object_inspect_modelmatrices <- function(lavobject, what = "free",
    type = "free", add.labels = FALSE, add.class = FALSE) {

    GLIST <- lavobject@Model@GLIST

    if(what == "dx.free") {
        DX <- lav_model_gradient(lavmodel       = lavobject@Model,
                                 GLIST          = NULL,
                                 lavsamplestats = lavobject@SampleStats,
                                 lavdata        = lavobject@Data,
                                 lavcache       = lavobject@Cache,
                                 type           = "free",
                                 estimator      = lavobject@Options$estimator,
                                 verbose        = FALSE,
                                 forcePD        = TRUE,
                                 group.weight   = TRUE,
                                 constraints    = FALSE,
                                 Delta          = NULL)
    } else if(what == "dx.all") {
        GLIST <- lav_model_gradient(lavmodel   = lavobject@Model,
                                GLIST          = NULL,
                                lavsamplestats = lavobject@SampleStats,
                                lavdata        = lavobject@Data,
                                lavcache       = lavobject@Cache,
                                type           = "allofthem",
                                estimator      = lavobject@Options$estimator,
                                verbose        = FALSE,
                                forcePD        = TRUE,
                                group.weight   = TRUE,
                                constraints    = FALSE,
                                Delta          = NULL)
        names(GLIST) <- names(lavobject@Model@GLIST)
    } else if(what == "std.all") {
        STD <- standardize.est.all(lavobject)
    } else if(what == "std.lv") {
        STD <- standardize.est.lv(lavobject)
    } else if(what == "std.nox") {
        STD <- standardize.est.all.nox(lavobject)
    }

    for(mm in 1:length(GLIST)) {

        if(what != "dx.all") {
            # erase everything
            GLIST[[mm]][,] <- 0.0
        }

        if(add.labels) {
            dimnames(GLIST[[mm]]) <- lavobject@Model@dimNames[[mm]]
        }

        if(what == "free") {
            # fill in free/unco parameter counts
            if(type == "free") {
                m.el.idx <- lavobject@Model@m.free.idx[[mm]]
                x.el.idx <- lavobject@Model@x.free.idx[[mm]]
            } else if(type == "unco") {
                m.el.idx <- lavobject@Model@m.unco.idx[[mm]]
                x.el.idx <- lavobject@Model@x.unco.idx[[mm]]
            } else if(type == "partable") {
                m.el.idx <- lavobject@Model@m.user.idx[[mm]]
                x.el.idx <- lavobject@Model@x.user.idx[[mm]]
            } else {
                stop("lavaan ERROR: unknown type argument:", type, )
            }
            GLIST[[mm]][m.el.idx] <- x.el.idx
        } else if(what == "se") {
            # fill in standard errors
            m.user.idx <- lavobject@Model@m.user.idx[[mm]]
            x.user.idx <- lavobject@Model@x.user.idx[[mm]]
            GLIST[[mm]][m.user.idx] <- lavobject@Fit@se[x.user.idx]
        } else if(what == "start") {
            # fill in starting values
            m.user.idx <- lavobject@Model@m.user.idx[[mm]]
            x.user.idx <- lavobject@Model@x.user.idx[[mm]]
            GLIST[[mm]][m.user.idx] <- lavobject@Fit@start[x.user.idx]
        } else if(what == "est") {
            # fill in estimated parameter values
            m.user.idx <- lavobject@Model@m.user.idx[[mm]]
            x.user.idx <- lavobject@Model@x.user.idx[[mm]]
            GLIST[[mm]][m.user.idx] <- lavobject@Fit@est[x.user.idx]
        } else if(what == "dx.free") {
            # fill in derivatives free parameters
            m.el.idx <- lavobject@Model@m.free.idx[[mm]]
            x.el.idx <- lavobject@Model@x.free.idx[[mm]]
            GLIST[[mm]][m.el.idx] <- DX[x.el.idx]
        } else if(what %in% c("std.all", "std.lv", "std.nox")) {
            m.user.idx <- lavobject@Model@m.user.idx[[mm]]
            x.user.idx <- lavobject@Model@x.user.idx[[mm]]
            GLIST[[mm]][m.user.idx] <- STD[x.user.idx]
        } 

        # class
        if(add.class) {
            if(lavobject@Model@isSymmetric[mm]) {
                class(GLIST[[mm]]) <- c("lavaan.matrix.symmetric", "matrix")
            } else {
                class(GLIST[[mm]]) <- c("lavaan.matrix", "matrix")
            }
        }
    }

    GLIST
}




# fixme, should we export this function?
lav_object_inspect_sampstat <- function(lavobject, 
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {

    G <- lavobject@Data@ngroups
    ov.names <- lavobject@Data@ov.names

    OUT <- vector("list", length=G)
    for(g in 1:G) {
        OUT[[g]]$cov  <- lavobject@SampleStats@cov[[g]]
        if(add.labels) {
            rownames(OUT[[g]]$cov) <- colnames(OUT[[g]]$cov) <- ov.names[[g]]
        }
        if(add.class) {
            class(OUT[[g]]$cov) <- c("lavaan.matrix.symmetric", "matrix")
        }

        #if(lavobject@Model@meanstructure) {
            OUT[[g]]$mean <- as.numeric(lavobject@SampleStats@mean[[g]])
            if(add.labels) {
                names(OUT[[g]]$mean) <- ov.names[[g]]
            }
            if(add.class) {
                class(OUT[[g]]$mean) <- c("lavaan.vector", "numeric")
            }
        #}

        if(lavobject@Model@categorical) {
            OUT[[g]]$th <- as.numeric(lavobject@SampleStats@th[[g]])
            if(length(lavobject@Model@num.idx[[g]]) > 0L) {
                OUT[[g]]$th <- OUT[[g]]$th[-lavobject@Model@num.idx[[g]]]
            }
            if(add.labels) {
                names(OUT[[g]]$th) <- 
                    vnames(lavobject@ParTable, type="th", group=g)
            }
            if(add.class) {
                class(OUT[[g]]$th) <- c("lavaan.vector", "numeric")
            }
        }

        if(lavobject@Model@categorical &&
           lavobject@Model@nexo > 0L) {
            OUT[[g]]$slopes  <- lavobject@SampleStats@slopes[[g]]
            if(add.labels) {
                rownames(OUT[[g]]$slopes) <- ov.names[[g]]
                colnames(OUT[[g]]$slopes) <- 
                    vnames(lavobject@ParTable, type="ov.x", group=g)
            }
            if(add.class) {
                class(OUT[[g]]$slopes) <- c("lavaan.matrix", "matrix")
            }
        }

        if(lavobject@Model@group.w.free) {
            OUT[[g]]$group.w <- lavobject@SampleStats@group.w[[g]]
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
        if(length(lavobject@Data@group.label) > 0L) {
            names(OUT) <- unlist(lavobject@Data@group.label)
        }
    }

    OUT
}


lav_object_inspect_data <- function(lavobject, add.labels = FALSE,
                                    drop.list.single.group = FALSE) {

    G <- lavobject@Data@ngroups
    OUT <- lavobject@Data@X

    if(add.labels) {
        for(g in 1:G) {
            colnames(OUT[[g]]) <- lavobject@Data@ov.names[[g]]
        }
    }

    if(G == 1L && drop.list.single.group) {
        OUT <- OUT[[1]]
    } else {
        if(length(lavobject@Data@group.label) > 0L) {
            names(OUT) <- unlist(lavobject@Data@group.label)
        }
    }

    OUT
}

lav_object_inspect_case_idx <- function(lavobject,
                                        drop.list.single.group = FALSE) {

    G <- lavobject@Data@ngroups
    OUT <- lavobject@Data@case.idx

    if(G == 1L && drop.list.single.group) {
        OUT <- OUT[[1]]
    } else {
        if(length(lavobject@Data@group.label) > 0L) {
            names(OUT) <- unlist(lavobject@Data@group.label)
        }
    }

    OUT
}



lav_object_inspect_rsquare <- function(lavobject, est.std.all=NULL,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {

    G <- lavobject@Data@ngroups
    OUT <- vector("list", length=G)

    if(is.null(est.std.all)) {
        est.std.all <- standardize.est.all(lavobject)
    }

    partable <- lavobject@ParTable
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

        if(add.labels && length(tmp) > 0L) {
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
        if(length(lavobject@Data@group.label) > 0L) {
            names(OUT) <- unlist(lavobject@Data@group.label)
        }
    }

    OUT 
}


lav_object_inspect_cov_lv <- function(lavobject, correlation.metric = FALSE,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {

    G <- lavobject@Data@ngroups

    # compute lv covar
    OUT <- computeVETA(lavmodel       = lavobject@Model, 
                       lavsamplestats = lavobject@SampleStats,
                       remove.dummy.lv = TRUE)

    # cor + labels + class
    for(g in 1:G) {

        if(correlation.metric && nrow(OUT[[g]]) > 1L) {
            # note: cov2cor fails if matrix is empty!
            OUT[[g]] <- cov2cor(OUT[[g]])
        }

        if(add.labels) {
            colnames(OUT[[g]]) <- rownames(OUT[[g]]) <- 
                lavobject@pta$vnames$lv.regular[[g]]
        }

        if(add.class) {
            class(OUT[[g]]) <- c("lavaan.matrix.symmetric", "matrix")
        }
    }

    if(G == 1L && drop.list.single.group) {
        OUT <- OUT[[1]]
    } else {
        if(length(lavobject@Data@group.label) > 0L) {
            names(OUT) <- unlist(lavobject@Data@group.label)
        }
    }

    OUT
}

lav_object_inspect_mean_lv <- function(lavobject,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {

    G <- lavobject@Data@ngroups

    # compute lv means
    OUT <- computeEETA(lavmodel       = lavobject@Model, 
                       lavsamplestats = lavobject@SampleStats,
                       remove.dummy.lv = TRUE)
    OUT <- lapply(OUT, as.numeric)

    # labels + class
    for(g in 1:G) {
        if(add.labels && length(OUT[[g]]) > 0L) {
            names(OUT[[g]]) <- lavobject@pta$vnames$lv.regular[[g]]
        }
        if(add.class) {
            class(OUT[[g]]) <- c("lavaan.vector", "numeric")
        }
    }

    if(G == 1L && drop.list.single.group) {
        OUT <- OUT[[1]]
    } else {
        if(length(lavobject@Data@group.label) > 0L) {
            names(OUT) <- unlist(lavobject@Data@group.label)
        }
    }

    OUT
}

lav_object_inspect_cov_all <- function(lavobject, correlation.metric = FALSE,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {

    G <- lavobject@Data@ngroups

    # compute extended model implied covariance matrix (both ov and lv)
    OUT <- computeCOV(lavmodel = lavobject@Model, 
                      lavsamplestats = lavobject@SampleStats,
                      remove.dummy.lv = TRUE)

    # cor + labels + class
    for(g in 1:G) {

        if(correlation.metric && nrow(OUT[[g]]) > 1L) {
            # note: cov2cor fails if matrix is empty!
            OUT[[g]] <- cov2cor(OUT[[g]])
        }

        if(add.labels) {
            NAMES <- c(lavobject@pta$vnames$ov.model[[g]],
                       lavobject@pta$vnames$lv.regular[[g]])
            colnames(OUT[[g]]) <- rownames(OUT[[g]]) <- NAMES
        }
        if(add.class) {
            class(OUT[[g]]) <- c("lavaan.matrix.symmetric", "matrix")
        }
    }

    if(G == 1L && drop.list.single.group) {
        OUT <- OUT[[1]]
    } else {
        if(length(lavobject@Data@group.label) > 0L) {
            names(OUT) <- unlist(lavobject@Data@group.label)
        }
    }

    OUT
}


lav_object_inspect_cov_ov <- function(lavobject, correlation.metric = FALSE,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {

    G <- lavobject@Data@ngroups

    # get model-implied covariance matrix observed
    OUT <- lavobject@Fit@Sigma.hat

    # cor + labels + class
    for(g in 1:G) {
 
        if(correlation.metric && nrow(OUT[[g]]) > 1L) {
            # note: cov2cor fails if matrix is empty!
            OUT[[g]] <- cov2cor(OUT[[g]])
        }

        if(add.labels) {
            colnames(OUT[[g]]) <- rownames(OUT[[g]]) <- 
                lavobject@pta$vnames$ov.model[[g]]
        }
        if(add.class) {
            class(OUT[[g]]) <- c("lavaan.matrix.symmetric", "matrix")
        }
    }

    if(G == 1L && drop.list.single.group) {
        OUT <- OUT[[1]]
    } else {
        if(length(lavobject@Data@group.label) > 0L) {
            names(OUT) <- unlist(lavobject@Data@group.label)
        }
    }

    OUT
}

lav_object_inspect_mean_ov <- function(lavobject,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {

    G <- lavobject@Data@ngroups

    # compute lv means
    OUT <- lavobject@Fit@Mu.hat
    OUT <- lapply(OUT, as.numeric)

    # labels + class
    for(g in 1:G) {
        if(add.labels && length(OUT[[g]]) > 0L) {
            names(OUT[[g]]) <- lavobject@pta$vnames$ov.model[[g]]
        }
        if(add.class) {
            class(OUT[[g]]) <- c("lavaan.vector", "numeric")
        }
    }

    if(G == 1L && drop.list.single.group) {
        OUT <- OUT[[1]]
    } else {
        if(length(lavobject@Data@group.label) > 0L) {
            names(OUT) <- unlist(lavobject@Data@group.label)
        }
    }

    OUT
}

lav_object_inspect_th <- function(lavobject,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {

    G <- lavobject@Data@ngroups

    # thresholds
    OUT <- lavobject@Fit@TH
    OUT <- lapply(OUT, as.numeric)

    # labels + class
    for(g in 1:G) {
        if(length(lavobject@Model@num.idx[[g]]) > 0L) {
            OUT[[g]] <- OUT[[g]][ -lavobject@Model@num.idx[[g]] ]
        }
        if(add.labels && length(OUT[[g]]) > 0L) {
            names(OUT[[g]]) <- lavobject@pta$vnames$th[[g]]
        }
        if(add.class) {
            class(OUT[[g]]) <- c("lavaan.vector", "numeric")
        }
    }

    if(G == 1L && drop.list.single.group) {
        OUT <- OUT[[1]]
    } else {
        if(length(lavobject@Data@group.label) > 0L) {
            names(OUT) <- unlist(lavobject@Data@group.label)
        }
    }

    OUT
}

lav_object_inspect_vy <- function(lavobject,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {

    G <- lavobject@Data@ngroups

    # 'unconditional' model-implied variances
    #  - same as diag(Sigma.hat) if all Y are continuous)
    #  - 1.0 (or delta^2) if categorical
    #  - if also Gamma, cov.x is used (only if categorical)

    OUT <- computeVY(lavmodel = lavobject@Model, GLIST = NULL, 
                     lavsamplestats = lavobject@SampleStats)
                     

    # labels + class
    for(g in 1:G) {
        if(add.labels && length(OUT[[g]]) > 0L) {
            if(lavobject@Model@categorical) {
                 names(OUT[[g]]) <- lavobject@pta$vnames$ov.nox[[g]]
            } else {
                 names(OUT[[g]]) <- lavobject@pta$vnames$ov[[g]]
            }
        }
        if(add.class) {
            class(OUT[[g]]) <- c("lavaan.vector", "numeric")
        }
    }

    if(G == 1L && drop.list.single.group) {
        OUT <- OUT[[1]]
    } else {
        if(length(lavobject@Data@group.label) > 0L) {
            names(OUT) <- unlist(lavobject@Data@group.label)
        }
    }

    OUT
}


lav_object_inspect_theta <- function(lavobject, correlation.metric = FALSE,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {

    G <- lavobject@Data@ngroups

    # get residual covariances
    OUT <- computeTHETA(lavmodel = lavobject@Model)

    # labels + class
    for(g in 1:G) {
        
        if(correlation.metric && nrow(OUT[[g]]) > 0L) {
            OUT[[g]] <- cov2cor(OUT[[g]])
        }

        if(add.labels && length(OUT[[g]]) > 0L) {
            colnames(OUT[[g]]) <- rownames(OUT[[g]]) <- 
                lavobject@pta$vnames$ov.model[[g]]
        }

        if(add.class) {
            class(OUT[[g]]) <- c("lavaan.vector", "numeric")
        }
    }

    if(G == 1L && drop.list.single.group) {
        OUT <- OUT[[1]]
    } else {
        if(length(lavobject@Data@group.label) > 0L) {
            names(OUT) <- unlist(lavobject@Data@group.label)
        }
    }

    OUT
}


lav_object_inspect_missing_coverage <- function(lavobject,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {

    G <- lavobject@Data@ngroups

    # get missing covarage
    OUT <- vector("list", G)
   
    for(g in 1:G) {
        if(!is.null(lavobject@Data@Mp[[g]])) {
            OUT[[g]] <- lavobject@Data@Mp[[g]]$coverage
        } else {
            nvar <- length(lavobject@Data@ov.names[[g]])
            OUT[[g]] <- matrix(1.0, nvar, nvar)
        }

        if(add.labels && length(OUT[[g]]) > 0L) {
            colnames(OUT[[g]]) <- rownames(OUT[[g]]) <-
                lavobject@pta$vnames$ov.model[[g]]
        }

        if(add.class) {
            class(OUT[[g]]) <- c("lavaan.matrix.symmetric", "matrix")
        }
    }

    if(G == 1L && drop.list.single.group) {
        OUT <- OUT[[1]]
    } else {
        if(length(lavobject@Data@group.label) > 0L) {
            names(OUT) <- unlist(lavobject@Data@group.label)
        }
    }

    OUT
}

lav_object_inspect_missing_patterns <- function(lavobject,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {

    G <- lavobject@Data@ngroups

    # get missing covarage
    OUT <- vector("list", G)
   
    for(g in 1:G) {
        if(!is.null(lavobject@Data@Mp[[g]])) {
            OUT[[g]] <- lavobject@Data@Mp[[g]]$pat
        } else {
            nvar <- length(lavobject@Data@ov.names[[g]])
            OUT[[g]] <- matrix(TRUE, 1L, nvar)
            rownames(OUT[[g]]) <- lavobject@Data@nobs[[g]]
        }

        if(add.labels && length(OUT[[g]]) > 0L) {
            colnames(OUT[[g]]) <- lavobject@pta$vnames$ov.model[[g]]
        }

        if(add.class) {
            class(OUT[[g]]) <- c("lavaan.matrix", "matrix")
        }
    }

    if(G == 1L && drop.list.single.group) {
        OUT <- OUT[[1]]
    } else {
        if(length(lavobject@Data@group.label) > 0L) {
            names(OUT) <- unlist(lavobject@Data@group.label)
        }
    }

    OUT
}

lav_object_inspect_empty_idx <- function(lavobject,
                                         drop.list.single.group = FALSE) {
   
    G <- lavobject@Data@ngroups

    # get empty idx
    OUT <- vector("list", G)

    for(g in 1:G) {
        if(!is.null(lavobject@Data@Mp[[g]])) {
            OUT[[g]] <- lavobject@Data@Mp[[g]]$empty.idx
        } else {
            OUT[[g]] <- integer(0L)
        }
    }
    
    if(G == 1L && drop.list.single.group) {
        OUT <- OUT[[1]]
    } else {
        if(length(lavobject@Data@group.label) > 0L) {
            names(OUT) <- unlist(lavobject@Data@group.label)
        }
    }

    OUT
}


lav_object_inspect_wls_est <- function(lavobject,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {

    G <- lavobject@Data@ngroups
    OUT <- lav_model_wls_est(lavobject@Model)

    for(g in 1:G) {
        if(add.labels && length(OUT[[g]]) > 0L) {
            #FIXME!!!!
            #names(OUT[[g]]) <- ??
        }

        if(add.class) {
            class(OUT[[g]]) <- c("lavaan.vector", "numeric")
        }
    }

    if(G == 1L && drop.list.single.group) {
        OUT <- OUT[[1]]
    } else {
        if(length(lavobject@Data@group.label) > 0L) {
            names(OUT) <- unlist(lavobject@Data@group.label)
        }
    }

    OUT
}

lav_object_inspect_wls_obs <- function(lavobject,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {

    G <- lavobject@Data@ngroups
    OUT <- lavobject@SampleStats@WLS.obs

    for(g in 1:G) {
        if(add.labels && length(OUT[[g]]) > 0L) {
            #FIXME!!!!
            #names(OUT[[g]]) <- ??
        }

        if(add.class) {
            class(OUT[[g]]) <- c("lavaan.vector", "numeric")
        }
    }

    if(G == 1L && drop.list.single.group) {
        OUT <- OUT[[1]]
    } else {
        if(length(lavobject@Data@group.label) > 0L) {
            names(OUT) <- unlist(lavobject@Data@group.label)
        }
    }

    OUT
}

lav_object_inspect_wls_v <- function(lavobject,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {

    # shortcuts
    G <- lavobject@Data@ngroups

    OUT <- lav_model_wls_v(lavmodel       = lavobject@Model,
                           lavsamplestats = lavobject@SampleStats,
                           estimator      = lavobject@Options$estimator,
                           lavdata        = lavobject@Data)

    # if estimator == "DWLS" or "ULS", we only stored the diagonal
    # hence, we create a full matrix here
    if(lavobject@Options$estimator %in% c("DWLS", "ULS")) {
        OUT <- lapply(OUT, 
            function(x) { nr = NROW(x); diag(x, nrow=nr, ncol=nr) })
    }

    # label + class
    for(g in 1:G) {
        if(add.labels && nrow(OUT[[g]]) > 0L) {
            #FIXME!!!!
            #names(OUT[[g]]) <- ??
        }

        if(add.class) {
            class(OUT[[g]]) <- c("lavaan.matrix", "matrix")
        }
    }

    if(G == 1L && drop.list.single.group) {
        OUT <- OUT[[1]]
    } else {
        if(length(lavobject@Data@group.label) > 0L) {
            names(OUT) <- unlist(lavobject@Data@group.label)
        }
    }

    OUT
}


lav_object_inspect_sampstat_nacov <- function(lavobject,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {

    if(lavobject@Options$se == "robust.mlr")
        stop("not done yet; FIX THIS!")

    # shortcuts
    lavsamplestats = lavobject@SampleStats
    estimator   = lavobject@Options$estimator
    G <- lavobject@Data@ngroups

    OUT <- vector("list", length=lavsamplestats@ngroups)

    if(estimator == "GLS"  ||
       estimator == "WLS"  ||
       estimator == "DWLS" ||
       estimator == "ULS") {
        OUT <- lavsamplestats@NACOV
    } else if(estimator == "ML" && !lavobject@SampleStats@missing.flag) {
        for(g in 1:lavsamplestats@ngroups) {
            OUT[[g]] <- 
                compute.Gamma(lavobject@Data@X[[g]], 
                              meanstructure=lavobject@Options$meanstructure)
            if(lavobject@Options$mimic == "Mplus") {
                G11 <- ( lavsamplestats@cov[[g]] * (lavsamplestats@nobs[[g]]-1)/lavsamplestats@nobs[[g]] )
                OUT[[g]][1:nrow(G11), 1:nrow(G11)] <- G11
            }
        }
    }

    if(G == 1L && drop.list.single.group) {
        OUT <- OUT[[1]]
    } else {
        if(length(lavobject@Data@group.label) > 0L) {
            names(OUT) <- unlist(lavobject@Data@group.label)
        }
    }

    OUT
}




lav_object_inspect_hessian <- function(lavobject,
    add.labels = FALSE, add.class = FALSE) {

    # lazy approach: take -1 the observed information
    E <- computeObservedInformation(lavmodel       = lavobject@Model, 
                                    lavsamplestats = lavobject@SampleStats,
                                    lavdata        = lavobject@Data,
                                    lavcache       = lavobject@Cache,
                                    type           = "free",
                                    estimator      = lavobject@Options$estimator,
                                    group.weight   = TRUE)

    OUT <- -E

    # labels
    if(add.labels) {
        # FIXME!
    }

    # class
    if(add.class) {
        class(OUT) <- c("lavaan.matrix.symmetric", "matrix")
    }

    OUT
}

lav_object_inspect_firstorder <- function(lavobject,
    add.labels = FALSE, add.class = FALSE) {

    # lazy approach: get it from Nvcov.first.order
    NACOV <- Nvcov.first.order(lavmodel       = lavobject@Model,
                               lavsamplestats = lavobject@SampleStats,
                               lavdata        = lavobject@Data,
                               estimator      = lavobject@Options$estimator)

    B0 <- attr(NACOV, "B0")

    if(lavobject@Options$estimator == "PML") {
        B0 <- B0 * lavobject@SampleStats@ntotal
    }
    
    OUT <- B0

    # labels
    if(add.labels) {
        # FIXME!
    }

    # class
    if(add.class) {
        class(OUT) <- c("lavaan.matrix.symmetric", "matrix")
    }

    OUT
}

