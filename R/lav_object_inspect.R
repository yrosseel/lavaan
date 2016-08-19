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
                    list.by.group          = FALSE,  
                    drop.list.single.group = FALSE) {

    lavInspect(lavobject = lavobject, what = what,
               add.labels = add.labels, add.class = add.class,
               list.by.group = list.by.group,
               drop.list.single.group =  drop.list.single.group)
}

# the `user' version: with defaults for display only
lavInspect <- function(lavobject,
                       what                   = "free",
                       add.labels             = TRUE,
                       add.class              = TRUE,
                       list.by.group          = TRUE,
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
        lav_object_inspect_modelmatrices(lavobject, what = "free", 
            type = "free", add.labels = add.labels, add.class = add.class,
            list.by.group = list.by.group, 
            drop.list.single.group = drop.list.single.group)
    } else if(what == "partable" || what == "user") {
        lav_object_inspect_modelmatrices(lavobject, what = "free", 
            type="partable", add.labels = add.labels, add.class = add.class,
            list.by.group = list.by.group, 
            drop.list.single.group = drop.list.single.group)
    } else if(what == "se" ||
              what == "std.err" ||
              what == "standard.errors") {
        lav_object_inspect_modelmatrices(lavobject, what = "se",
            add.labels = add.labels, add.class = add.class,
            list.by.group = list.by.group, 
            drop.list.single.group = drop.list.single.group)
    } else if(what == "start" || what == "starting.values") {
        lav_object_inspect_modelmatrices(lavobject, what = "start",
            add.labels = add.labels, add.class = add.class,
            list.by.group = list.by.group, 
            drop.list.single.group = drop.list.single.group)
    } else if(what == "est"  || what == "estimates" ||
              what == "coef" || what == "coefficients" ||
              what == "x") {
        lav_object_inspect_modelmatrices(lavobject, what = "est",
            add.labels = add.labels, add.class = add.class,
            list.by.group = list.by.group, 
            #list.by.group = FALSE, for semTools only
            drop.list.single.group = drop.list.single.group)
    } else if(what == "dx.free") {
        lav_object_inspect_modelmatrices(lavobject, what = "dx.free",
            add.labels = add.labels, add.class = add.class,
            list.by.group = list.by.group, 
            drop.list.single.group = drop.list.single.group)
    } else if(what == "dx.all") {
        lav_object_inspect_modelmatrices(lavobject, what = "dx.all",
            add.labels = add.labels, add.class = add.class,
            list.by.group = list.by.group, 
            drop.list.single.group = drop.list.single.group)
    } else if(what == "std" || what == "std.all" || what == "standardized") {
        lav_object_inspect_modelmatrices(lavobject, what = "std.all",
           add.labels = add.labels, add.class = add.class,
           list.by.group = list.by.group, 
           drop.list.single.group = drop.list.single.group)
    } else if(what == "std.lv") {
        lav_object_inspect_modelmatrices(lavobject, what = "std.lv",
           add.labels = add.labels, add.class = add.class,
           list.by.group = list.by.group, 
           drop.list.single.group = drop.list.single.group)
    } else if(what == "std.nox") {
        lav_object_inspect_modelmatrices(lavobject, what = "std.nox",
           add.labels = add.labels, add.class = add.class,
           list.by.group = list.by.group, 
           drop.list.single.group = drop.list.single.group)


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
              what == "sampstats" ||
              what == "samplestats" ||
              what == "samp" ||
              what == "sample" ||
              what == "samplestatistics") {
        lav_object_inspect_sampstat(lavobject, h1 = FALSE,
            add.labels = add.labels, add.class = add.class, 
            drop.list.single.group = drop.list.single.group)
    } else if(what == "h1" || what == "missing.h1" || what == "sampstat.h1") {
        lav_object_inspect_sampstat(lavobject, h1 = TRUE,
            add.labels = add.labels, add.class = add.class,
            drop.list.single.group = drop.list.single.group)

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



    #### data + missingness ####
    } else if(what == "data") {
        lav_object_inspect_data(lavobject, 
            drop.list.single.group = drop.list.single.group)
    } else if(what == "case.idx") {
        lav_object_inspect_case_idx(lavobject,
            drop.list.single.group = drop.list.single.group) 
    } else if(what == "ngroups") {
        lavobject@Data@ngroups
    } else if(what == "group") {
        lavobject@Data@group
    } else if(what == "group.label") {
        lavobject@Data@group.label
    } else if(what == "nobs") {
        unlist( lavobject@Data@nobs )
    } else if(what == "norig") {
        unlist( lavobject@Data@norig )
    } else if(what == "ntotal") {
        sum(unlist( lavobject@Data@nobs ))
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
 


    #### convergence, meanstructure, categorical ####
    } else if(what == "converged") {
        lavobject@optim$converged
    } else if(what == "iterations" ||
              what == "iter" ||
              what == "niter") {
        lavobject@optim$iterations
    } else if(what == "meanstructure") {
        lavobject@Model@meanstructure
    } else if(what == "categorical") {
        lavobject@Model@categorical
    } else if(what == "fixed.x") {
        lavobject@Model@fixed.x
    } else if(what == "parameterization") {
        lavobject@Model@parameterization
    


    #### NACOV samplestats ####
    } else if(what == "gamma") {
        lav_object_inspect_sampstat_gamma(lavobject,
            add.labels = add.labels, add.class = add.class,
            drop.list.single.group = drop.list.single.group)


    #### gradient, Hessian, information, first.order, vcov ####
    } else if(what == "gradient") {
        lav_object_inspect_gradient(lavobject,
            add.labels = add.labels, add.class = add.class)
    } else if(what == "hessian") {
        lav_object_inspect_hessian(lavobject,
            add.labels = add.labels, add.class = add.class)

    } else if(what == "information") {
        lav_object_inspect_information(lavobject, information = "default",
            augmented = FALSE, inverted = FALSE,
            add.labels = add.labels, add.class = add.class)
    } else if(what == "information.expected") {
        lav_object_inspect_information(lavobject, information = "expected",
            augmented = FALSE, inverted = FALSE,
            add.labels = add.labels, add.class = add.class)
    } else if(what == "information.observed") {
        lav_object_inspect_information(lavobject, information = "observed",
            augmented = FALSE, inverted = FALSE,
            add.labels = add.labels, add.class = add.class)
    } else if(what == "information.first.order" || what == "first.order") {
        lav_object_inspect_information(lavobject, information = "first.order",
            augmented = FALSE, inverted = FALSE,
            add.labels = add.labels, add.class = add.class)

    } else if(what == "augmented.information") {
        lav_object_inspect_information(lavobject, information = "default",
            augmented = TRUE, inverted = FALSE,
            add.labels = add.labels, add.class = add.class)
    } else if(what == "augmented.information.expected") {
        lav_object_inspect_information(lavobject, information = "expected",
            augmented = TRUE, inverted = FALSE,
            add.labels = add.labels, add.class = add.class)
    } else if(what == "augmented.information.observed") {
        lav_object_inspect_information(lavobject, information = "observed",
            augmented = TRUE, inverted = FALSE,
            add.labels = add.labels, add.class = add.class)
    } else if(what == "augmented.information.first.order" ||
              what == "augmented.first.order") {
        lav_object_inspect_information(lavobject, information = "first.order",
            augmented = TRUE, inverted = FALSE,
            add.labels = add.labels, add.class = add.class)

    } else if(what == "inverted.information") {
        lav_object_inspect_information(lavobject, information = "default",
            augmented = TRUE, inverted = TRUE,
            add.labels = add.labels, add.class = add.class)
    } else if(what == "inverted.information.expected") {
        lav_object_inspect_information(lavobject, information = "expected",
            augmented = TRUE, inverted = TRUE,
            add.labels = add.labels, add.class = add.class)
    } else if(what == "inverted.information.observed") {
        lav_object_inspect_information(lavobject, information = "observed",
            augmented = TRUE, inverted = TRUE,
            add.labels = add.labels, add.class = add.class)
    } else if(what == "inverted.information.first.order" || 
              what == "inverted.first.order") {
        lav_object_inspect_information(lavobject, information = "first.order",
            augmented = TRUE, inverted = TRUE,
            add.labels = add.labels, add.class = add.class)

    } else if(what == "vcov") {
        lav_object_inspect_vcov(lavobject,
            standardized = FALSE,
            add.labels = add.labels, add.class = add.class)
    } else if(what == "vcov.std.all" || what == "vcov.standardized" ||
              what == "vcov.std") {
        lav_object_inspect_vcov(lavobject,
            standardized = TRUE, type = "std.all",
            add.labels = add.labels, add.class = add.class)
    } else if(what == "vcov.std.lv") {
        lav_object_inspect_vcov(lavobject,
            standardized = TRUE, type = "std.lv",
            add.labels = add.labels, add.class = add.class)
    } else if(what == "vcov.std.nox") {
        lav_object_inspect_vcov(lavobject,
            standardized = TRUE, type = "std.nox",
            add.labels = add.labels, add.class = add.class)
    } else if(what == "vcov.def") {
        lav_object_inspect_vcov_def(lavobject,
            standardized = FALSE,
            add.labels = add.labels, add.class = add.class)
    } else if(what == "vcov.def.std.all" || what == "vcov.def.standardized" ||
              what == "vcov.def.std") {
        lav_object_inspect_vcov_def(lavobject,
            standardized = TRUE, type = "std.all",
            add.labels = add.labels, add.class = add.class)
    } else if(what == "vcov.def.std.lv") {
        lav_object_inspect_vcov_def(lavobject,
            standardized = TRUE, type = "std.lv",
            add.labels = add.labels, add.class = add.class)
    } else if(what == "vcov.def.std.nox") {
        lav_object_inspect_vcov_def(lavobject,
            standardized = TRUE, type = "std.nox",
            add.labels = add.labels, add.class = add.class)
    } else if(what == "ugamma" || what == "ug" || what == "u.gamma") {
        lav_object_inspect_UGamma(lavobject,
            add.labels = add.labels, add.class = add.class)

    ### jacobians ####
    } else if(what == "delta") {
        lav_object_inspect_delta(lavobject,
            add.labels = add.labels, add.class = add.class,
            drop.list.single.group = drop.list.single.group)

    # post-checking
    } else if(what == "post.check" || what == "post") {
        lav_object_post_check(lavobject)

    # options
    } else if(what == "options" || what == "lavoptions") {
        lavobject@Options

    # call
    } else if(what == "call") {
        as.list( lavobject@call )

    # timing
    } else if(what == "timing") {
        lavobject@timing

    # optim
    } else if(what == "optim") {
        lavobject@optim

    # test
    } else if(what == "test") {
        lavobject@test

    #### not found ####
    } else {
        stop("unknown `what' argument in inspect function: `", what, "'")
    }

}


# helper functions (mostly to deal with older 'object' that may have
# been save somewhere)
lav_object_inspect_est <- function(lavobject) {
    
    # from 0.5-19, they are in the partable
    if(!is.null(lavobject@ParTable$est)) {
        OUT <- lavobject@ParTable$est
    } else if("Fit" %in% slotNames(lavobject)) {
        # in < 0.5-19, we should look in @Fit@est
        OUT <- lavobject@Fit@est
    } else {
        PT <- parTable(lavobject)
        OUT <- rep(as.numeric(NA), length(PT$lhs))
    }

    OUT
}

lav_object_inspect_se <- function(lavobject) {
    
    # from 0.5-19, they are in the partable
    if(!is.null(lavobject@ParTable$se)) {
        OUT <- lavobject@ParTable$se
    } else if("Fit" %in% slotNames(lavobject)) {
        # in < 0.5-19, we should look in @Fit@se
        OUT <- lavobject@Fit@se
    } else {
        PT <- parTable(lavobject)
        OUT <- rep(as.numeric(NA), length(PT$lhs))
    }

    OUT
}

lav_object_inspect_start <- function(lavobject) {

    # from 0.5-19, they are in the partable
    if(!is.null(lavobject@ParTable$start)) {
        OUT <- lavobject@ParTable$start
    } else {
        # in < 0.5-19, we should look in @Fit@start
        OUT <- lavobject@Fit@start
    }

    OUT
}

lav_object_inspect_boot <- function(lavobject) {

    # from 0.5-19. they are in a separate slot
    tmp <- try(slot(lavobject,"boot"), silent = TRUE)
    if(inherits(tmp, "try-error")) {
        # older version of object?
        est <- lav_object_inspect_est(lavobject)
        BOOT <- attr(est, "BOOT.COEF")
    } else {
        # 0.5-19 way
        BOOT <- lavobject@boot$coef
    }

    BOOT
}


lav_object_inspect_modelmatrices <- function(lavobject, what = "free",
    type = "free", add.labels = FALSE, add.class = FALSE,
    list.by.group = FALSE,
    drop.list.single.group = FALSE) {

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

        if(add.labels) {
            dimnames(GLIST[[mm]]) <- lavobject@Model@dimNames[[mm]]
        }

        if(what == "free") {
            # fill in free parameter counts
            if(type == "free") {
                m.el.idx <- lavobject@Model@m.free.idx[[mm]]
                x.el.idx <- lavobject@Model@x.free.idx[[mm]]
            #} else if(type == "unco") {
            #    m.el.idx <- lavobject@Model@m.unco.idx[[mm]]
            #    x.el.idx <- lavobject@Model@x.unco.idx[[mm]]
            } else if(type == "partable") {
                m.el.idx <- lavobject@Model@m.user.idx[[mm]]
                x.el.idx <- lavobject@Model@x.user.idx[[mm]]
            } else {
                stop("lavaan ERROR: unknown type argument:", type, )
            }
            # erase everything
            GLIST[[mm]][,] <- 0.0
            GLIST[[mm]][m.el.idx] <- x.el.idx
        } else if(what == "se") {
            # fill in standard errors
            m.user.idx <- lavobject@Model@m.user.idx[[mm]]
            x.user.idx <- lavobject@Model@x.user.idx[[mm]]
            SE <- lav_object_inspect_se(lavobject)
            # erase everything
            GLIST[[mm]][,] <- 0.0
            GLIST[[mm]][m.user.idx] <- SE[x.user.idx]
        } else if(what == "start") {
            # fill in starting values
            m.user.idx <- lavobject@Model@m.user.idx[[mm]]
            x.user.idx <- lavobject@Model@x.user.idx[[mm]]
            START <- lav_object_inspect_start(lavobject)
            GLIST[[mm]][m.user.idx] <- START[x.user.idx]
        } else if(what == "est") {
            # fill in estimated parameter values
            m.user.idx <- lavobject@Model@m.user.idx[[mm]]
            x.user.idx <- lavobject@Model@x.user.idx[[mm]]
            EST <- lav_object_inspect_est(lavobject)
            GLIST[[mm]][m.user.idx] <- EST[x.user.idx]
        } else if(what == "dx.free") {
            # fill in derivatives free parameters
            m.el.idx <- lavobject@Model@m.free.idx[[mm]]
            x.el.idx <- lavobject@Model@x.free.idx[[mm]]
            # erase everything
            GLIST[[mm]][,] <- 0.0
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

    # try to reflect `equality constraints'
    con.flag <- FALSE
    if(what == "free" && lavobject@Model@eq.constraints) {
        # extract constraints from parameter table
        PT <- parTable(lavobject)
        CON <-  PT[PT$op %in% c("==","<",">") ,c("lhs","op","rhs")]
        rownames(CON) <- NULL

        # replace 'labels' by parameter numbers
        ID <- lav_partable_constraints_label_id(PT)
        LABEL <- names(ID)
        for(con in 1:nrow(CON)) {
            # lhs
            LHS.labels <- all.vars(as.formula(paste("~",CON[con,"lhs"]))) 

            if(length(LHS.labels) > 0L) {
                # par id
                LHS.freeid <- ID[match(LHS.labels, LABEL)]

                # substitute
                tmp <- CON[con,"lhs"]
                for(pat in 1:length(LHS.labels)) {
                    tmp <- sub(LHS.labels[pat], LHS.freeid[pat], tmp)
                }
                CON[con,"lhs"] <- tmp
            }

            # rhs
            RHS.labels <- all.vars(as.formula(paste("~",CON[con,"rhs"])))

            if(length(RHS.labels) > 0L) {
                # par id
                RHS.freeid <- ID[match(RHS.labels, LABEL)]
                # substitute
                tmp <- CON[con,"rhs"]
                for(pat in 1:length(RHS.labels)) {
                    tmp <- sub(RHS.labels[pat], RHS.freeid[pat], tmp)
                }
                CON[con,"rhs"] <- tmp
            }
        } # con

        # add this info at the top
        #GLIST <- c(constraints = list(CON), GLIST)
        #no, not a good idea, it does not work with list.by.group
  
        # add it as a 'header' attribute?
        attr(CON, "header") <- "Note: model contains equality constraints:"
        con.flag <- TRUE
    }

    # should we group them per group?
    if(list.by.group) {
        lavsamplestats <- lavobject@SampleStats
        lavmodel       <- lavobject@Model
        nmat           <- lavmodel@nmat

        OUT <- vector("list", length = lavsamplestats@ngroups)
        for(g in 1:lavsamplestats@ngroups) {
            # which mm belong to group g?
            mm.in.group <- 1:nmat[g] + cumsum(c(0,nmat))[g]
            mm.names <- names( GLIST[mm.in.group] )

            OUT[[g]] <- GLIST[mm.in.group]
        }

        if(lavsamplestats@ngroups == 1L && drop.list.single.group) {
            OUT <- OUT[[1]]
        } else {
            if(length(lavobject@Data@group.label) > 0L) {
                names(OUT) <- unlist(lavobject@Data@group.label)
            }
        }
    } else {
        OUT <- GLIST
    }

    # header
    if(con.flag) {
        attr(OUT, "header") <- CON
    }

    # lavaan.list
    if(add.class) {
        class(OUT) <- c("lavaan.list", "list")
    }

    OUT
}




# - fixme, should we export this function?
# - since 0.5-21, conditional.x = TRUE returns residual sample statistics
#    for ML, we have both joint and residual cov/var/...; but for 
#    categorical = TRUE, we only have residual cov/var...; so, we
#    only return residual in both cases, whenever residual
lav_object_inspect_sampstat <- function(lavobject, h1 = FALSE,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {

    G <- lavobject@Data@ngroups
    ov.names <- lavobject@pta$vnames$ov
    ov.names.res <- lavobject@pta$vnames$ov.nox
    ov.names.x   <- lavobject@pta$vnames$ov.x
    lavsamplestats <- lavobject@SampleStats

    OUT <- vector("list", length=G)
    for(g in 1:G) {

        if(!lavobject@Model@conditional.x) {

            # covariance matrix
            if(h1 && !is.null(lavsamplestats@missing.h1[[g]])) {
                OUT[[g]]$cov  <- lavsamplestats@missing.h1[[g]]$sigma
            } else {
                OUT[[g]]$cov  <- lavsamplestats@cov[[g]]
            }
            if(add.labels && !is.null(OUT[[g]]$cov)) {
                rownames(OUT[[g]]$cov) <- colnames(OUT[[g]]$cov) <- 
                    ov.names[[g]]
            }
            if(add.class) {
                class(OUT[[g]]$cov) <- c("lavaan.matrix.symmetric", "matrix")
            }

            # mean vector
            if(h1 && !is.null(lavsamplestats@missing.h1[[g]])) {
                OUT[[g]]$mean <- lavsamplestats@missing.h1[[g]]$mu
            } else {
                OUT[[g]]$mean <- as.numeric(lavsamplestats@mean[[g]])
            }
            if(add.labels) {
                names(OUT[[g]]$mean) <- ov.names[[g]]
            }
            if(add.class) {
                class(OUT[[g]]$mean) <- c("lavaan.vector", "numeric")
            }

            # thresholds
            if(lavobject@Model@categorical) {
                OUT[[g]]$th <- as.numeric(lavsamplestats@th[[g]])
                if(length(lavobject@Model@num.idx[[g]]) > 0L) {
                    NUM.idx <- which(lavobject@Model@th.idx[[g]] == 0)
                    OUT[[g]]$th <- OUT[[g]]$th[ -NUM.idx ]
                }
                if(add.labels) {
                    names(OUT[[g]]$th) <- lavobject@pta$vnames$th[[g]]
                }
                if(add.class) {
                    class(OUT[[g]]$th) <- c("lavaan.vector", "numeric")
                }
            }
        } # !conditional.x

        else { # if conditional.x = TRUE

            # residual covariance matrix
            OUT[[g]]$res.cov  <- lavsamplestats@res.cov[[g]]
            if(add.labels) {
                rownames(OUT[[g]]$res.cov) <- colnames(OUT[[g]]$res.cov) <- 
                    ov.names.res[[g]]
            }
            if(add.class) {
                class(OUT[[g]]$res.cov) <- 
                    c("lavaan.matrix.symmetric", "matrix")
            }
   
            # intercepts
            if(lavobject@Model@conditional.x) {
                OUT[[g]]$res.int <- as.numeric(lavsamplestats@res.int[[g]])
                if(add.labels) {
                    names(OUT[[g]]$res.int) <- ov.names.res[[g]]
                }
                if(add.class) {
                    class(OUT[[g]]$res.int) <- c("lavaan.vector", "numeric")
                }
            }

            # thresholds
            if(lavobject@Model@categorical) {
                OUT[[g]]$res.th <- as.numeric(lavsamplestats@res.th[[g]])
                if(length(lavobject@Model@num.idx[[g]]) > 0L) {
                    NUM.idx <- which(lavobject@Model@th.idx[[g]] == 0)
                    OUT[[g]]$res.th <- OUT[[g]]$res.th[ -NUM.idx ]
                }
                if(add.labels) {
                    names(OUT[[g]]$res.th) <- lavobject@pta$vnames$th[[g]]
                }
                if(add.class) {
                    class(OUT[[g]]$res.th) <- c("lavaan.vector", "numeric")
                }
            }

            # slopes
            if(lavobject@Model@nexo > 0L) {
                OUT[[g]]$res.slopes  <- lavsamplestats@res.slopes[[g]]
                if(add.labels) {
                    rownames(OUT[[g]]$res.slopes) <- ov.names.res[[g]]
                    colnames(OUT[[g]]$res.slopes) <- ov.names.x[[g]]
                }
                if(add.class) {
                    class(OUT[[g]]$res.slopes) <- c("lavaan.matrix", "matrix")
                }
            }

            # cov.x
            if(lavobject@Model@nexo > 0L) {
                OUT[[g]]$cov.x  <- lavsamplestats@cov.x[[g]]
                if(add.labels) {
                    rownames(OUT[[g]]$cov.x) <- ov.names.x[[g]]
                    colnames(OUT[[g]]$cov.x) <- ov.names.x[[g]]
                }
                if(add.class) {
                    class(OUT[[g]]$cov.x) <- 
                        c("lavaan.matrix.symmetric", "matrix")
                }
            }

        } # conditional.x

        # stochastic weights
        if(lavobject@Model@group.w.free) {
            OUT[[g]]$group.w <- lavsamplestats@group.w[[g]]
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
                lavobject@pta$vnames$lv[[g]]
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
    OUT <- lavobject@implied$cov

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
    OUT <- lavobject@implied$mean
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
    OUT <- lavobject@implied$th
    OUT <- lapply(OUT, as.numeric)

    # labels + class
    for(g in 1:G) {
        if(length(lavobject@Model@num.idx[[g]]) > 0L) {
            NUM.idx <- which(lavobject@Model@th.idx[[g]] == 0)
            OUT[[g]] <- OUT[[g]][ -NUM.idx ]
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
                     lavsamplestats = lavobject@SampleStats,
                     diagonal.only = TRUE)
                     

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
    OUT <- lav_model_wls_est(lavobject@Model) #,
                             #cov.x = lavobject@SampleStats@cov.x)

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


lav_object_inspect_sampstat_gamma <- function(lavobject,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {

    # shortcuts
    G <- lavobject@Data@ngroups

    if(!is.null(lavobject@SampleStats@NACOV[[1]])) {
        OUT <- lavobject@SampleStats@NACOV
    } else {
        OUT <- lavGamma(lavobject)
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


lav_object_inspect_gradient <- function(lavobject,
    add.labels = FALSE, add.class = FALSE) {

    if(lavobject@SampleStats@missing.flag ||
       lavobject@Options$estimator == "PML") {
        group.weight <- FALSE
    } else {
        group.weight <- TRUE
    }

    OUT <- lav_model_gradient(lavmodel       = lavobject@Model,
                              GLIST          = NULL,
                              lavsamplestats = lavobject@SampleStats,
                              lavdata        = lavobject@Data,
                              lavcache       = lavobject@Cache,
                              type           = "free",
                              estimator      = lavobject@Options$estimator,
                              verbose        = FALSE,
                              group.weight   = group.weight)

    # labels
    if(add.labels) {
        names(OUT) <- lav_partable_labels(lavobject@ParTable, type="free")
    }

    # class
    if(add.class) {
        class(OUT) <- c("lavaan.vector", "numeric")
    }

    OUT
}

lav_object_inspect_hessian <- function(lavobject,
    add.labels = FALSE, add.class = FALSE) {

    OUT <- lav_model_hessian(lavmodel       = lavobject@Model, 
                             lavsamplestats = lavobject@SampleStats,
                             lavdata        = lavobject@Data,
                             lavcache       = lavobject@Cache,
                             estimator      = lavobject@Options$estimator,
                             group.weight   = TRUE)

    # labels
    if(add.labels) {
        colnames(OUT) <- rownames(OUT) <-
            lav_partable_labels(lavobject@ParTable, type="free")
    }

    # class
    if(add.class) {
        class(OUT) <- c("lavaan.matrix.symmetric", "matrix")
    }

    OUT
}

lav_object_inspect_information <- function(lavobject, 
    information = "default", augmented = FALSE, inverted = FALSE,
    add.labels = FALSE, add.class = FALSE) {

    if(information == "default") {
        information <- lavobject@Options$information
    } 

    if(information == "expected" || information == "observed") {
        OUT <- lav_model_information(lavmodel =  lavobject@Model,
                  lavsamplestats = lavobject@SampleStats,
                  lavdata        = lavobject@Data,
                  estimator      = lavobject@Options$estimator,
                  lavcache       = lavobject@Cache,
                  information    = information,
                  augmented      = augmented,
                  inverted       = inverted)
    } else if(information == "first.order") {
        B0 <- lav_model_information_firstorder(lavmodel =  lavobject@Model,
              lavsamplestats = lavobject@SampleStats,
              lavdata        = lavobject@Data,
              estimator      = lavobject@Options$estimator,
              lavcache       = lavobject@Cache,
              check.pd       = FALSE,
              augmented      = augmented,
              inverted       = inverted)
        attr(B0, "B0.group") <- NULL
        OUT <- B0
    }

    # labels
    if(add.labels) {
        NAMES <- lav_partable_labels(lavobject@ParTable, type="free")
        if(augmented) {
            nExtra <- nrow(OUT) - length(NAMES)
            if(nExtra > 0L) {
                NAMES <- c(NAMES, paste("aug", 1:nExtra, sep=""))
            }
        }
        colnames(OUT) <- rownames(OUT) <- NAMES
    }

    # class
    if(add.class) {
        class(OUT) <- c("lavaan.matrix.symmetric", "matrix")
    }

    OUT
}

# only to provide a direct function to the old 'getVariability()' function
lav_object_inspect_firstorder <- function(lavobject, 
    add.labels = FALSE, add.class = FALSE) {

     B0 <- lav_model_information_firstorder(lavmodel =  lavobject@Model,
              lavsamplestats = lavobject@SampleStats,
              lavdata        = lavobject@Data,
              estimator      = lavobject@Options$estimator,
              lavcache       = lavobject@Cache,
              check.pd       = FALSE,
              augmented      = FALSE,
              inverted       = FALSE)
    attr(B0, "B0.group") <- NULL
    OUT <- B0

    # labels
    if(add.labels) {
        colnames(OUT) <- rownames(OUT) <-
            lav_partable_labels(lavobject@ParTable, type="free")
    }

    # class
    if(add.class) {
        class(OUT) <- c("lavaan.matrix.symmetric", "matrix")
    }

    OUT
}

lav_object_inspect_vcov <- function(lavobject, standardized = FALSE,
    type = "std.all", free.only = TRUE,
    add.labels = FALSE, add.class = FALSE, remove.duplicated = FALSE) {

    npar <- max(lavobject@ParTable$free)
    if(lavobject@optim$npar == 0) {
        OUT <- matrix(0,0,0)
    } else {
        # check if we already have it
        tmp <- try(slot(lavobject, "vcov"), silent = TRUE)
        if(!inherits(tmp, "try-error") && !is.null(lavobject@vcov$vcov)) {
            OUT <- lavobject@vcov$vcov
        } else {
        # compute it again
            OUT <- lav_model_vcov(lavmodel       = lavobject@Model,
                                  lavsamplestats = lavobject@SampleStats,
                                  lavoptions     = lavobject@Options,
                                  lavdata        = lavobject@Data,
                                  lavcache       = lavobject@Cache
                                 )
        }
    }
   
    # strip attributes
    attr(OUT, "E.inv") <- NULL
    attr(OUT, "B0") <- NULL
    attr(OUT, "B0.group") <- NULL
    attr(OUT, "Delta") <- NULL
    attr(OUT, "WLS.V") <- NULL
    attr(OUT, "BOOT.COEF") <- NULL
    attr(OUT, "BOOT.TEST") <- NULL

    # standardized?
    if(standardized) {
        if(type == "std.lv") {
            JAC <- try(lav_func_jacobian_complex(func = standardize.est.lv.x,
                x = lavobject@optim$x, lavobject = lavobject), silent = TRUE)
            if(inherits(JAC, "try-error")) { # eg. pnorm()
                JAC <- lav_func_jacobian_simple(func = standardize.est.lv.x,
                    x = lavobject@optim$x, lavobject=lavobject)
            }
        } else if(type == "std.all") {
            JAC <- try(lav_func_jacobian_complex(func = standardize.est.all.x,
                x = lavobject@optim$x, lavobject = lavobject), silent = TRUE)
            if(inherits(JAC, "try-error")) { # eg. pnorm()
                JAC <- lav_func_jacobian_simple(func = standardize.est.all.x,
                    x = lavobject@optim$x, lavobject=lavobject)
            }
        } else if(type == "std.nox") {
            JAC <- 
                try(lav_func_jacobian_complex(func = standardize.est.all.nox.x,
                    x = lavobject@optim$x, lavobject = lavobject), silent = TRUE)
            if(inherits(JAC, "try-error")) { # eg. pnorm()
                JAC <- 
                    lav_func_jacobian_simple(func = standardize.est.all.nox.x,
                        x = lavobject@optim$x, lavobject=lavobject)
            }
        }

        # JAC contains *all* parameters in the parameter table
        if(free.only) { 
            free.idx <- which(lavobject@ParTable$free > 0L)
            JAC <- JAC[free.idx,, drop = FALSE]
        }
        OUT <- JAC %*% OUT %*% t(JAC)
    }

    # labels
    if(add.labels) {
        colnames(OUT) <- rownames(OUT) <-
            lav_partable_labels(lavobject@ParTable, type="free")
    }

    # alias?
    if(remove.duplicated && lavobject@Model@eq.constraints) {
        simple.flag <- lav_constraints_check_simple(lavobject@Model)
        if(simple.flag) {
            LAB <- lav_partable_labels(lavobject@ParTable, type="free")
            dup.flag <- duplicated(LAB)
            OUT <- OUT[!dup.flag, !dup.flag, drop = FALSE]
        } else {
            warning("lavaan WARNING: alias is TRUE, but equality constraints do not appear to be simple; returning full vcov")
        }
    }

    # class
    if(add.class) {
        class(OUT) <- c("lavaan.matrix.symmetric", "matrix")
    }

    OUT
}

lav_object_inspect_vcov_def <- function(lavobject, standardized = FALSE,
    type = "std.all", add.labels = FALSE, add.class = FALSE) {

    lavmodel    <- lavobject@Model
    lavpartable <- lavobject@ParTable
    def.idx <- which(lavpartable$op == ":=")

    if(length(def.idx) == 0L) {
        return( matrix(0,0,0) )
    }

    if(standardized) {
        # compute VCOV for "free" parameters only
        VCOV <- lav_object_inspect_vcov(lavobject = lavobject, 
                                        standardized = TRUE,
                                        type = type, free.only = FALSE,
                                        add.labels = FALSE, add.class = FALSE)
        OUT <- VCOV[def.idx, def.idx, drop = FALSE]
    } else {

        # get free parameters
        x <- lav_model_get_parameters(lavmodel, type = "free")

        # bootstrap or not?
        if(!is.null(lavobject@boot$coef)) {
            BOOT <- lavobject@boot$coef
            BOOT.def <- apply(BOOT, 1L, lavmodel@def.function)
            if(length(def.idx) == 1L) {
                BOOT.def <- as.matrix(BOOT.def)
            } else {
                BOOT.def <- t(BOOT.def)
            }
            OUT <- cov(BOOT.def)
        } else {
            # VCOV
            VCOV <- lav_object_inspect_vcov(lavobject = lavobject,
                                            standardized = FALSE,
                                            type = type, free.only = TRUE,
                                            add.labels = FALSE, 
                                            add.class = FALSE)
    
            # regular delta method
            JAC <- try(lav_func_jacobian_complex(func = lavmodel@def.function,
                       x = x), silent=TRUE)
            if(inherits(JAC, "try-error")) { # eg. pnorm()
                JAC <- lav_func_jacobian_simple(func = lavmodel@def.function,
                                                x = x)
            }
            OUT <- JAC %*% VCOV %*% t(JAC)
        }
    }

    # labels
    if(add.labels) {
        LHS.names <- lavpartable$lhs[def.idx]
        colnames(OUT) <- rownames(OUT) <- LHS.names
    }

    # class
    if(add.class) {
        class(OUT) <- c("lavaan.matrix.symmetric", "matrix")
    }

    OUT
}

lav_object_inspect_UGamma <- function(lavobject,
    add.labels = FALSE, add.class = FALSE) {

    out <- lav_test_satorra_bentler(lavobject     = lavobject,
                                    return.ugamma = TRUE)
    OUT <- out$UGamma

    # labels
    if(add.labels) {
       # colnames(OUT) <- rownames(OUT) <-
    }

    # class
    if(add.class) {
        class(OUT) <- c("lavaan.matrix.symmetric", "matrix")
    }

    OUT
}

# Delta (jacobian: d samplestats / d free_parameters)
lav_object_inspect_delta <- function(lavobject,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {

    OUT <- computeDelta(lavobject@Model)

    # labels
    lavmodel <- lavobject@Model
    categorical    <- lavmodel@categorical
    conditional.x  <- lavmodel@conditional.x
    group.w.free   <- lavmodel@group.w.free
    nvar           <- lavmodel@nvar
    num.idx        <- lavmodel@num.idx
    th.idx         <- lavmodel@th.idx
    nexo           <- lavmodel@nexo
    ngroups        <- lavmodel@ngroups

    if(add.labels) {
        PNAMES <- lav_partable_labels(lavobject@ParTable, type="free")

        for(g in 1:ngroups) {
            colnames(OUT[[g]]) <- PNAMES

            if(conditional.x) {
                ov.names <- lavobject@pta$vnames$ov.nox[[g]]
            } else {
                ov.names <- lavobject@pta$vnames$ov[[g]]
            }
            ov.names.x <- lavobject@pta$vnames$ov.x[[g]]
            nvar <- length(ov.names)


            names.cov <- names.cor <- names.var <- character(0L)
            names.mu <- names.pi <- names.th <- character(0L)
            names.gw <- character(0L)
 
            # Sigma
            # - if continuous: vech(Sigma)
            # - if categorical: first numeric variances, then
            tmp <- apply(expand.grid(ov.names, ov.names), 1L, 
                         paste, collapse = "~~")
            if(categorical) {
                names.cor <- tmp[lav_matrix_vech_idx(nvar, diagonal = FALSE)]
                names.var <- tmp[lav_matrix_diag_idx(nvar)[num.idx[[g]]]]
            } else {
                names.cov <- tmp[lav_matrix_vech_idx(nvar, diagonal = TRUE)]
            }

            # Mu
            if(!categorical && lavmodel@meanstructure) {
                names.mu <- paste(ov.names, "~1", sep = "")
            }

            # Pi
            if(conditional.x && lavmodel@nexo[g] > 0L) {
               names.pi <- apply(expand.grid(ov.names, ov.names.x), 1L, 
                                 paste, collapse = "~")
            }

            # th
            if(categorical) {
                names.th <- lavobject@pta$vnames$th[[g]]
                # interweave numeric intercepts, if any
                if(length(num.idx[[g]]) > 0L) {
                    tmp <- character( length(th.idx[[g]]) )
                    tmp[ th.idx[[g]] > 0 ] <- names.th
                    tmp[ th.idx[[g]] == 0 ] <- paste(ov.names[ num.idx[[g]] ],
                                                     "~1", sep = "")
                    names.th <- tmp
                }
            }

            # gw
            if(group.w.free) {
                names.gw <- "w"
            }

            rownames(OUT[[g]]) <- c(names.gw,
                                    names.th, names.mu, 
                                    names.pi, 
                                    names.cov, names.var, names.cor)

            # class
            if(add.class) {
                class(OUT[[g]]) <- c("lavaan.matrix", "matrix")
            }

        } # g
    } # labels

    if(ngroups == 1L && drop.list.single.group) {
        OUT <- OUT[[1]]
    } else {
        if(length(lavobject@Data@group.label) > 0L) {
            names(OUT) <- unlist(lavobject@Data@group.label)
        }
    }

    OUT
}


