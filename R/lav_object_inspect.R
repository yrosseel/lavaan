# inspect a fitted lavaan object

# backward compatibility -- wrapper around lavInspect
inspect.lavaan <- function(object, what = "free", ...) {
    lavInspect.lavaan(object              = object,
                      what                   = what,
                      add.labels             = TRUE,
                      add.class              = TRUE,
                      drop.list.single.group = TRUE)
}

# the `tech' version: no labels, full matrices, ... for further processing
lavTech.lavaan <- function(object,
                           what                   = "free",
                           add.labels             = FALSE,
                           add.class              = FALSE,
                           list.by.group          = FALSE,
                           drop.list.single.group = FALSE) {

    lavInspect.lavaan(object, what = what,
                      add.labels = add.labels, add.class = add.class,
                      list.by.group = list.by.group,
                      drop.list.single.group =  drop.list.single.group)
}

# the `user' version: with defaults for display only
lavInspect.lavaan <- function(object,
                              what                   = "free",
                              add.labels             = TRUE,
                              add.class              = TRUE,
                              list.by.group          = TRUE,
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
        lav_object_inspect_modelmatrices(object, what = "free",
            type = "free", add.labels = add.labels, add.class = add.class,
            list.by.group = list.by.group,
            drop.list.single.group = drop.list.single.group)
    } else if(what == "impute" ||
              what == "imputed") { # just to ease the transition for semTools!
        object@imputed
    } else if(what == "partable" || what == "user") {
        lav_object_inspect_modelmatrices(object, what = "free",
            type="partable", add.labels = add.labels, add.class = add.class,
            list.by.group = list.by.group,
            drop.list.single.group = drop.list.single.group)
    } else if(what == "se" ||
              what == "std.err" ||
              what == "standard.errors") {
        lav_object_inspect_modelmatrices(object, what = "se",
            add.labels = add.labels, add.class = add.class,
            list.by.group = list.by.group,
            drop.list.single.group = drop.list.single.group)
    } else if(what == "start" || what == "starting.values") {
        lav_object_inspect_modelmatrices(object, what = "start",
            add.labels = add.labels, add.class = add.class,
            list.by.group = list.by.group,
            drop.list.single.group = drop.list.single.group)
    } else if(what == "est"  || what == "estimates" ||
              what == "x") {
        lav_object_inspect_modelmatrices(object, what = "est",
            add.labels = add.labels, add.class = add.class,
            list.by.group = list.by.group,
            #list.by.group = FALSE, for semTools only
            drop.list.single.group = drop.list.single.group)
    } else if(what == "dx.free") {
        lav_object_inspect_modelmatrices(object, what = "dx.free",
            add.labels = add.labels, add.class = add.class,
            list.by.group = list.by.group,
            drop.list.single.group = drop.list.single.group)
    } else if(what == "dx.all") {
        lav_object_inspect_modelmatrices(object, what = "dx.all",
            add.labels = add.labels, add.class = add.class,
            list.by.group = list.by.group,
            drop.list.single.group = drop.list.single.group)
    } else if(what == "std" || what == "std.all" || what == "standardized") {
        lav_object_inspect_modelmatrices(object, what = "std.all",
           add.labels = add.labels, add.class = add.class,
           list.by.group = list.by.group,
           drop.list.single.group = drop.list.single.group)
    } else if(what == "std.lv") {
        lav_object_inspect_modelmatrices(object, what = "std.lv",
           add.labels = add.labels, add.class = add.class,
           list.by.group = list.by.group,
           drop.list.single.group = drop.list.single.group)
    } else if(what == "std.nox") {
        lav_object_inspect_modelmatrices(object, what = "std.nox",
           add.labels = add.labels, add.class = add.class,
           list.by.group = list.by.group,
           drop.list.single.group = drop.list.single.group)


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
    } else if(what == "obs" ||
              what == "observed" ||
              what == "sampstat" ||
              what == "sampstats" ||
              what == "samplestats" ||
              what == "samp" ||
              what == "sample" ||
              what == "samplestatistics") {
        # new in 0.6-3: always use h1 = TRUE!!!
        lav_object_inspect_sampstat(object, h1 = TRUE, std = FALSE,
            add.labels = add.labels, add.class = add.class,
            drop.list.single.group = drop.list.single.group)
    } else if(what == "obs.std" ||
              what == "observed.std" ||
              what == "sampstat.std" ||
              what == "sampstats.std" ||
              what == "samplestats.std" ||
              what == "samp.std" ||
              what == "sample.std" ||
              what == "samplestatistics.std") {
        lav_object_inspect_sampstat(object, h1 = TRUE, std = TRUE,
            add.labels = add.labels, add.class = add.class,
            drop.list.single.group = drop.list.single.group)
    } else if(what == "h1" || what == "missing.h1" || what == "sampstat.h1") {
        lav_object_inspect_sampstat(object, h1 = TRUE,
            add.labels = add.labels, add.class = add.class,
            drop.list.single.group = drop.list.single.group)

    #### wls.est - wls.obs - wls.v ####
    } else if(what == "wls.est") {
        lav_object_inspect_wls_est(object,
            add.labels = add.labels, add.class = add.class,
            drop.list.single.group = drop.list.single.group)
    } else if(what == "wls.obs") {
        lav_object_inspect_wls_obs(object,
            add.labels = add.labels, add.class = add.class,
            drop.list.single.group = drop.list.single.group)
    } else if(what == "wls.v") {
        lav_object_inspect_wls_v(object,
            add.labels = add.labels, add.class = add.class,
            drop.list.single.group = drop.list.single.group)



    #### data + missingness ####
    } else if(what == "data") {
        lav_object_inspect_data(object, add.labels = add.labels,
            drop.list.single.group = drop.list.single.group)
    } else if(what == "case.idx") {
        lav_object_inspect_case_idx(object,
            drop.list.single.group = drop.list.single.group)
    } else if(what == "ngroups") {
        object@Data@ngroups
    } else if(what == "group") {
        object@Data@group
    } else if(what == "cluster") {
      object@Data@cluster
    } else if(what == "nlevels") {
      object@Data@nlevels
    } else if(what == "ordered") {
        object@Data@ordered
    } else if(what == "group.label") {
        object@Data@group.label
    } else if(what == "nobs") {
        unlist( object@Data@nobs )
    } else if(what == "norig") {
        unlist( object@Data@norig )
    } else if(what == "ntotal") {
        sum(unlist( object@Data@nobs ))
    } else if(what == "coverage") {
        lav_object_inspect_missing_coverage(object,
            add.labels = add.labels, add.class = add.class,
            drop.list.single.group = drop.list.single.group)
    } else if(what %in% c("patterns", "pattern")) {
        lav_object_inspect_missing_patterns(object,
            add.labels = add.labels, add.class = add.class,
            drop.list.single.group = drop.list.single.group)
    } else if(what == "empty.idx") {
        lav_object_inspect_empty_idx(object,
            drop.list.single.group = drop.list.single.group)


    #### rsquare ####
    } else if(what == "rsquare" || what == "r-square" || what == "r2") {
         lav_object_inspect_rsquare(object,
             add.labels = add.labels, add.class = add.class,
             drop.list.single.group = drop.list.single.group)


    #### model-implied sample statistics ####
    } else if(what == "implied" || what == "fitted" ||
              what == "expected" || what == "exp") {
        lav_object_inspect_implied(object,
            add.labels = add.labels, add.class = add.class,
            drop.list.single.group = drop.list.single.group)
    } else if(what == "resid" || what == "res" || what == "residual" ||
              what == "residuals") {
        lav_object_inspect_residuals(object, h1 = TRUE,
            add.labels = add.labels, add.class = add.class,
            drop.list.single.group = drop.list.single.group)
    } else if(what == "cov.lv" || what == "veta") {
        lav_object_inspect_cov_lv(object,
            correlation.metric = FALSE,
            add.labels = add.labels, add.class = add.class,
            drop.list.single.group = drop.list.single.group)
    } else if(what == "cor.lv") {
        lav_object_inspect_cov_lv(object,
            correlation.metric = TRUE,
            add.labels = add.labels, add.class = add.class,
            drop.list.single.group = drop.list.single.group)
    } else if(what == "mean.lv" || what == "eeta") {
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
    } else if(what == "cov.ov" || what == "sigma" || what == "sigma.hat") {
        lav_object_inspect_cov_ov(object,
            correlation.metric = FALSE,
            add.labels = add.labels, add.class = add.class,
            drop.list.single.group = drop.list.single.group)
    } else if(what == "cor.ov") {
        lav_object_inspect_cov_ov(object,
            correlation.metric = TRUE,
            add.labels = add.labels, add.class = add.class,
            drop.list.single.group = drop.list.single.group)
    } else if(what == "mean.ov" || what == "mu" || what == "mu.hat") {
        lav_object_inspect_mean_ov(object,
            add.labels = add.labels, add.class = add.class,
            drop.list.single.group = drop.list.single.group)
    } else if(what == "th" || what == "thresholds") {
        lav_object_inspect_th(object,
            add.labels = add.labels, add.class = add.class,
            drop.list.single.group = drop.list.single.group)
    } else if(what == "th.idx") {
        lav_object_inspect_th_idx(object,
            add.labels = add.labels, add.class = add.class,
            drop.list.single.group = drop.list.single.group)
    } else if(what == "vy") {
        lav_object_inspect_vy(object,
            add.labels = add.labels, add.class = add.class,
            drop.list.single.group = drop.list.single.group)


    #### specific model matrices? ####
    } else if(what == "theta" || what == "theta.cov") {
        lav_object_inspect_theta(object,  correlation.metric = FALSE,
            add.labels = add.labels, add.class = add.class,
            drop.list.single.group = drop.list.single.group)
    } else if(what == "theta.cor") {
        lav_object_inspect_theta(object,  correlation.metric = TRUE,
            add.labels = add.labels, add.class = add.class,
            drop.list.single.group = drop.list.single.group)



    #### convergence, meanstructure, categorical ####
    } else if(what == "converged") {
        object@optim$converged
    } else if(what == "iterations" ||
              what == "iter" ||
              what == "niter") {
        object@optim$iterations
    } else if(what == "meanstructure") {
        object@Model@meanstructure
    } else if(what == "categorical") {
        object@Model@categorical
    } else if(what == "fixed.x") {
        object@Model@fixed.x
    } else if(what == "parameterization") {
        object@Model@parameterization
    } else if(what == "npar") {
        lav_object_inspect_npar(object, type = "free")
    } else if(what == "coef") {
        # this breaks simsem and semTools -- 0.6-1
        #lav_object_inspect_coef(object, type = "free",
        #                        add.labels = add.labels, add.class = add.class)
        lav_object_inspect_modelmatrices(object, what = "est",
            type = "free", add.labels = add.labels, add.class = add.class,
            list.by.group = list.by.group,
            drop.list.single.group = drop.list.single.group)


    #### NACOV samplestats ####
    } else if(what == "gamma") {
        lav_object_inspect_sampstat_gamma(object,
            add.labels = add.labels, add.class = add.class,
            drop.list.single.group = drop.list.single.group)


    #### gradient, Hessian, information, first.order, vcov ####
    } else if(what == "gradient") {
        lav_object_inspect_gradient(object,
            add.labels = add.labels, add.class = add.class, logl = FALSE)
    } else if(what == "gradient.logl") {
        lav_object_inspect_gradient(object,
            add.labels = add.labels, add.class = add.class, logl = TRUE)
    } else if(what == "optim.gradient") {
        lav_object_inspect_gradient(object,
            add.labels = add.labels, add.class = add.class, optim = TRUE)
    } else if(what == "hessian") {
        lav_object_inspect_hessian(object,
            add.labels = add.labels, add.class = add.class)

    } else if(what == "information") {
        lav_object_inspect_information(object, information = "default",
            augmented = FALSE, inverted = FALSE,
            add.labels = add.labels, add.class = add.class)
    } else if(what == "information.expected") {
        lav_object_inspect_information(object, information = "expected",
            augmented = FALSE, inverted = FALSE,
            add.labels = add.labels, add.class = add.class)
    } else if(what == "information.observed") {
        lav_object_inspect_information(object, information = "observed",
            augmented = FALSE, inverted = FALSE,
            add.labels = add.labels, add.class = add.class)
    } else if(what == "information.first.order" ||
              what == "information.firstorder"  ||
              what == "first.order") {
        lav_object_inspect_information(object, information = "first.order",
            augmented = FALSE, inverted = FALSE,
            add.labels = add.labels, add.class = add.class)

    } else if(what == "augmented.information") {
        lav_object_inspect_information(object, information = "default",
            augmented = TRUE, inverted = FALSE,
            add.labels = add.labels, add.class = add.class)
    } else if(what == "augmented.information.expected") {
        lav_object_inspect_information(object, information = "expected",
            augmented = TRUE, inverted = FALSE,
            add.labels = add.labels, add.class = add.class)
    } else if(what == "augmented.information.observed") {
        lav_object_inspect_information(object, information = "observed",
            augmented = TRUE, inverted = FALSE,
            add.labels = add.labels, add.class = add.class)
    } else if(what == "augmented.information.first.order" ||
              what == "augmented.first.order") {
        lav_object_inspect_information(object, information = "first.order",
            augmented = TRUE, inverted = FALSE,
            add.labels = add.labels, add.class = add.class)

    } else if(what == "inverted.information") {
        lav_object_inspect_information(object, information = "default",
            augmented = TRUE, inverted = TRUE,
            add.labels = add.labels, add.class = add.class)
    } else if(what == "inverted.information.expected") {
        lav_object_inspect_information(object, information = "expected",
            augmented = TRUE, inverted = TRUE,
            add.labels = add.labels, add.class = add.class)
    } else if(what == "inverted.information.observed") {
        lav_object_inspect_information(object, information = "observed",
            augmented = TRUE, inverted = TRUE,
            add.labels = add.labels, add.class = add.class)
    } else if(what == "inverted.information.first.order" ||
              what == "inverted.first.order") {
        lav_object_inspect_information(object, information = "first.order",
            augmented = TRUE, inverted = TRUE,
            add.labels = add.labels, add.class = add.class)

    } else if(what == "vcov") {
        lav_object_inspect_vcov(object,
            standardized = FALSE,
            add.labels = add.labels, add.class = add.class)
    } else if(what == "vcov.std.all" || what == "vcov.standardized" ||
              what == "vcov.std") {
        lav_object_inspect_vcov(object,
            standardized = TRUE, type = "std.all",
            add.labels = add.labels, add.class = add.class)
    } else if(what == "vcov.std.lv") {
        lav_object_inspect_vcov(object,
            standardized = TRUE, type = "std.lv",
            add.labels = add.labels, add.class = add.class)
    } else if(what == "vcov.std.nox") {
        lav_object_inspect_vcov(object,
            standardized = TRUE, type = "std.nox",
            add.labels = add.labels, add.class = add.class)

    } else if(what == "vcov.def") {
        lav_object_inspect_vcov_def(object, joint = FALSE,
            standardized = FALSE,
            add.labels = add.labels, add.class = add.class)
    } else if(what == "vcov.def.std.all" || what == "vcov.def.standardized" ||
              what == "vcov.def.std") {
        lav_object_inspect_vcov_def(object, joint = FALSE,
            standardized = TRUE, type = "std.all",
            add.labels = add.labels, add.class = add.class)
    } else if(what == "vcov.def.std.lv") {
        lav_object_inspect_vcov_def(object, joint = FALSE,
            standardized = TRUE, type = "std.lv",
            add.labels = add.labels, add.class = add.class)
    } else if(what == "vcov.def.std.nox") {
        lav_object_inspect_vcov_def(object, joint = FALSE,
            standardized = TRUE, type = "std.nox",
            add.labels = add.labels, add.class = add.class)

    } else if(what == "vcov.def.joint") {
        lav_object_inspect_vcov_def(object, joint = TRUE,
            standardized = FALSE,
            add.labels = add.labels, add.class = add.class)
    } else if(what == "vcov.def.joint.std.all" ||
              what == "vcov.def.joint.standardized" ||
              what == "vcov.def.joint.std") {
        lav_object_inspect_vcov_def(object, joint = TRUE,
            standardized = TRUE, type = "std.all",
            add.labels = add.labels, add.class = add.class)
    } else if(what == "vcov.def.joint.std.lv") {
        lav_object_inspect_vcov_def(object, joint = TRUE,
            standardized = TRUE, type = "std.lv",
            add.labels = add.labels, add.class = add.class)
    } else if(what == "vcov.def.joint.std.nox") {
        lav_object_inspect_vcov_def(object, joint = TRUE,
            standardized = TRUE, type = "std.nox",
            add.labels = add.labels, add.class = add.class)

    } else if(what == "ugamma" || what == "ug" || what == "u.gamma") {
        lav_object_inspect_UGamma(object,
            add.labels = add.labels, add.class = add.class)
    } else if(what == "ufromugamma" || what == "u") {
        lav_object_inspect_UfromUGamma(object,
            add.labels = add.labels, add.class = add.class,
            drop.list.single.group = drop.list.single.group)

    ### jacobians ####
    } else if(what == "delta") {
        lav_object_inspect_delta(object,
            add.labels = add.labels, add.class = add.class,
            drop.list.single.group = drop.list.single.group)

    # multilevel #
    } else if(what == "icc") {
        lav_object_inspect_icc(object,
            add.labels = add.labels, add.class = add.class,
            drop.list.single.group = drop.list.single.group)
    } else if(what == "ranef") {
        lav_object_inspect_ranef(object,
            add.labels = add.labels, add.class = add.class,
            drop.list.single.group = drop.list.single.group)

    # post-checking
    } else if(what == "post.check" || what == "post") {
        lav_object_post_check(object)

    # options
    } else if(what == "options" || what == "lavoptions") {
        object@Options

    # version
    } else if(what == "version") {
        object@version

    # call
    } else if(what == "call") {
        as.list( object@call )

    # timing
    } else if(what == "timing") {
        object@timing

    # optim
    } else if(what == "optim") {
        object@optim

    # test
    } else if(what == "test") {
        object@test

    # zero cell tables
    } else if(what == "zero.cell.tables") {
        lav_object_inspect_zero_cell_tables(object,
            add.labels = add.labels, add.class = add.class,
            drop.list.single.group = drop.list.single.group)

    #### not found ####
    } else {
        stop("unknown `what' argument in inspect function: `", what, "'")
    }

}


# helper functions (mostly to deal with older 'object' that may have
# been saved somewhere)
lav_object_inspect_est <- function(object) {

    if(inherits(object, "lavaan")) {
        # from 0.5-19, they are in the partable
        if(!is.null(object@ParTable$est)) {
            OUT <- object@ParTable$est
        } else if(.hasSlot(object, "Fit")) {
            # in < 0.5-19, we should look in @Fit@est
            OUT <- object@Fit@est
        } else {
            PT <- parTable(object)
            OUT <- rep(as.numeric(NA), length(PT$lhs))
        }
    } else {
        # try generic coef()
        OUT <- coef(object, type = "user")
        if(is.matrix(OUT)) {
            # lavaanList?
            OUT <- rowMeans(OUT)
        }
    }

    OUT
}

lav_object_inspect_se <- function(object) {

    # from 0.5-19, they are in the partable
    if(!is.null(object@ParTable$se)) {
        OUT <- object@ParTable$se
    } else if(.hasSlot(object, "Fit")) {
        # in < 0.5-19, we should look in @Fit@se
        OUT <- object@Fit@se
    } else {
        PT <- parTable(object)
        OUT <- rep(as.numeric(NA), length(PT$lhs))
    }

    OUT
}

lav_object_inspect_start <- function(object) {

    # from 0.5-19, they are in the partable
    if(!is.null(object@ParTable$start)) {
        OUT <- object@ParTable$start
    } else {
        # in < 0.5-19, we should look in @Fit@start
        OUT <- object@Fit@start
    }

    OUT
}

lav_object_inspect_boot <- function(object) {

    # from 0.5-19. they are in a separate slot
    tmp <- try(slot(object,"boot"), silent = TRUE)
    if(inherits(tmp, "try-error")) {
        # older version of object?
        est <- lav_object_inspect_est(object)
        BOOT <- attr(est, "BOOT.COEF")
    } else {
        # 0.5-19 way
        BOOT <- object@boot$coef
    }

    BOOT
}


lav_object_inspect_modelmatrices <- function(object, what = "free",
    type = "free", add.labels = FALSE, add.class = FALSE,
    list.by.group = FALSE,
    drop.list.single.group = FALSE) {

    GLIST <- object@Model@GLIST

    if(what == "dx.free") {
        DX <- lav_model_gradient(lavmodel       = object@Model,
                                 GLIST          = NULL,
                                 lavsamplestats = object@SampleStats,
                                 lavdata        = object@Data,
                                 lavcache       = object@Cache,
                                 type           = "free",
                                 verbose        = FALSE,
                                 forcePD        = TRUE,
                                 group.weight   = TRUE,
                                 Delta          = NULL)
    } else if (what == "dx.all") {
        GLIST <- lav_model_gradient(lavmodel   = object@Model,
                                GLIST          = NULL,
                                lavsamplestats = object@SampleStats,
                                lavdata        = object@Data,
                                lavcache       = object@Cache,
                                type           = "allofthem",
                                verbose        = FALSE,
                                forcePD        = TRUE,
                                group.weight   = TRUE,
                                Delta          = NULL)
        names(GLIST) <- names(object@Model@GLIST)
    } else if (what == "std.all") {
        STD <- lav_standardize_all(object)
    } else if (what == "std.lv") {
        STD <- lav_standardize_lv(object)
    } else if (what == "std.nox") {
        STD <- lav_standardize_all_nox(object)
    } else if(what == "se") {
        SE <- lav_object_inspect_se(object)
    } else if (what == "start") {
        START <- lav_object_inspect_start(object)
    } else if (what == "est") {
        EST <- lav_object_inspect_est(object)
    }

    for(mm in 1:length(GLIST)) {

        if(add.labels) {
            dimnames(GLIST[[mm]]) <- object@Model@dimNames[[mm]]
        }

        if(what == "free") {
            # fill in free parameter counts
            if(type == "free") {
                m.el.idx <- object@Model@m.free.idx[[mm]]
                x.el.idx <- object@Model@x.free.idx[[mm]]
            #} else if(type == "unco") {
            #    m.el.idx <- object@Model@m.unco.idx[[mm]]
            #    x.el.idx <- object@Model@x.unco.idx[[mm]]
            } else if(type == "partable") {
                m.el.idx <- object@Model@m.user.idx[[mm]]
                x.el.idx <- object@Model@x.user.idx[[mm]]
            } else {
                stop("lavaan ERROR: unknown type argument:", type, )
            }
            # erase everything
            GLIST[[mm]][,] <- 0.0
            GLIST[[mm]][m.el.idx] <- x.el.idx
        } else if(what == "se") {
            # fill in standard errors
            m.user.idx <- object@Model@m.user.idx[[mm]]
            x.user.idx <- object@Model@x.user.idx[[mm]]
            # erase everything
            GLIST[[mm]][,] <- 0.0
            GLIST[[mm]][m.user.idx] <- SE[x.user.idx]
        } else if(what == "start") {
            # fill in starting values
            m.user.idx <- object@Model@m.user.idx[[mm]]
            x.user.idx <- object@Model@x.user.idx[[mm]]
            GLIST[[mm]][m.user.idx] <- START[x.user.idx]
        } else if(what == "est") {
            # fill in estimated parameter values
            m.user.idx <- object@Model@m.user.idx[[mm]]
            x.user.idx <- object@Model@x.user.idx[[mm]]
            GLIST[[mm]][m.user.idx] <- EST[x.user.idx]
        } else if(what == "dx.free") {
            # fill in derivatives free parameters
            m.el.idx <- object@Model@m.free.idx[[mm]]
            x.el.idx <- object@Model@x.free.idx[[mm]]
            # erase everything
            GLIST[[mm]][,] <- 0.0
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

    # try to reflect `equality constraints'
    con.flag <- FALSE
    if(what == "free" && object@Model@eq.constraints) {
        # extract constraints from parameter table
        PT <- parTable(object)
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

    # should we group them per block?
    if(list.by.group) {
        lavsamplestats <- object@SampleStats
        lavmodel       <- object@Model
        nmat           <- lavmodel@nmat

        OUT <- vector("list", length = lavmodel@nblocks)
        for(b in seq_len(lavmodel@nblocks)) {
            # which mm belong to this block?
            mm.in.group <- 1:nmat[b] + cumsum(c(0,nmat))[b]
            mm.names <- names( GLIST[mm.in.group] )

            OUT[[b]] <- GLIST[mm.in.group]
        }

        if(lavmodel@nblocks == 1L && drop.list.single.group) {
            OUT <- OUT[[1]]
        } else {
            if(object@Data@nlevels == 1L &&
               length(object@Data@group.label) > 0L) {
                names(OUT) <- unlist(object@Data@group.label)
            } else if(object@Data@nlevels > 1L &&
                      length(object@Data@group.label) == 0L) {
                names(OUT) <- object@Data@level.label
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
# - since 0.6-3, we always extract the values from the @h1 slot (if present)
#   if meanstructure = FALSE, do NOT include $mean elements any longer
lav_object_inspect_sampstat <- function(object, h1 = TRUE,
    std = FALSE,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {

    nblocks <- object@Model@nblocks
    ov.names <- object@pta$vnames$ov
    ov.names.res <- object@pta$vnames$ov.nox
    ov.names.x   <- object@pta$vnames$ov.x

    # slots
    lavsamplestats <- object@SampleStats
    lavmodel <- object@Model

    # if nlevels, override h1 to be TRUE
    if(object@Data@nlevels > 1L) {
        h1 <- TRUE
    }

    # check if we have a non-empty @h1 slot
    if(!.hasSlot(object, "h1")) {
        h1 <- FALSE
    } else if(length(object@h1) == 0L) {
        h1 <- FALSE
    } else {
        H1 <- object@h1$implied
    }

    # if h1 = FALSE and nlevels > 1L, nothing can show...
    if(!h1 && object@Data@nlevels > 1L) {
        stop("lavaan ERROR: sample statistics not available; refit with option h1 = TRUE")
    }

    OUT <- vector("list", length = nblocks)
    for(b in seq_len(nblocks)) {

        if(!lavmodel@conditional.x) {

            # covariance matrix
            if(h1) {
                OUT[[b]]$cov  <- H1$cov[[b]]
            } else {
                OUT[[b]]$cov  <- lavsamplestats@cov[[b]]
            }
            if(std) {
                diag.orig <- diag(OUT[[b]]$cov)
                OUT[[b]]$cov <- cov2cor(OUT[[b]]$cov)
            }
            if(add.labels && !is.null(OUT[[b]]$cov)) {
                rownames(OUT[[b]]$cov) <- colnames(OUT[[b]]$cov) <-
                    ov.names[[b]]
            }
            if(add.class) {
                class(OUT[[b]]$cov) <- c("lavaan.matrix.symmetric", "matrix")
            }

            # mean vector
            if(lavmodel@meanstructure) {
                if(h1) {
                    OUT[[b]]$mean <- as.numeric(H1$mean[[b]])
                } else {
                    OUT[[b]]$mean <- as.numeric(lavsamplestats@mean[[b]])
                }
                if(std) {
                    diag.orig[ diag.orig < .Machine$double.eps ] <- NA
                    OUT[[b]]$mean <- OUT[[b]]$mean/sqrt(diag.orig)
                }
                if(add.labels) {
                    names(OUT[[b]]$mean) <- ov.names[[b]]
                }
                if(add.class) {
                    class(OUT[[b]]$mean) <- c("lavaan.vector", "numeric")
                }
            }

            # thresholds
            if(lavmodel@categorical) {
                if(h1) {
                    OUT[[b]]$th <- as.numeric(H1$th[[b]])
                } else {
                    OUT[[b]]$th <- as.numeric(lavsamplestats@th[[b]])
                }
                if(length(lavmodel@num.idx[[b]]) > 0L) {
                    NUM.idx <- which(lavmodel@th.idx[[b]] == 0)
                    OUT[[b]]$th <- OUT[[b]]$th[ -NUM.idx ]
                }
                # FIXME: what to do if std = TRUE (depends on delta/theta)
                if(add.labels) {
                    names(OUT[[b]]$th) <- object@pta$vnames$th[[b]]
                }
                if(add.class) {
                    class(OUT[[b]]$th) <- c("lavaan.vector", "numeric")
                }
            }
        } # !conditional.x

        else { # if conditional.x = TRUE

            # residual covariance matrix
            if(h1) {
                OUT[[b]]$res.cov  <- H1$res.cov[[b]]
            } else {
                OUT[[b]]$res.cov  <- lavsamplestats@res.cov[[b]]
            }
            if(std) {
                diag.orig <- diag(OUT[[b]]$res.cov)
                OUT[[b]]$res.cov <- cov2cor(OUT[[b]]$res.cov)
            }
            if(add.labels) {
                rownames(OUT[[b]]$res.cov) <- colnames(OUT[[b]]$res.cov) <-
                    ov.names.res[[b]]
            }
            if(add.class) {
                class(OUT[[b]]$res.cov) <-
                    c("lavaan.matrix.symmetric", "matrix")
            }

            # intercepts
            if(lavmodel@meanstructure) {
                if(h1) {
                    OUT[[b]]$res.int <- as.numeric(H1$res.int[[b]])
                } else {
                    OUT[[b]]$res.int <- as.numeric(lavsamplestats@res.int[[b]])
                }
                if(std) {
                    diag.orig[ diag.orig < .Machine$double.eps ] <- NA
                    OUT[[b]]$res.int <- OUT[[b]]$res.int/sqrt(diag.orig)
                }
                if(add.labels) {
                    names(OUT[[b]]$res.int) <- ov.names.res[[b]]
                }
                if(add.class) {
                    class(OUT[[b]]$res.int) <- c("lavaan.vector", "numeric")
                }
            }

            # thresholds
            if(lavmodel@categorical) {
                if(h1) {
                    OUT[[b]]$res.th <- as.numeric(H1$res.th[[b]])
                } else {
                    OUT[[b]]$res.th <- as.numeric(lavsamplestats@res.th[[b]])
                }
                if(length(lavmodel@num.idx[[b]]) > 0L) {
                    NUM.idx <- which(lavmodel@th.idx[[b]] == 0)
                    OUT[[b]]$res.th <- OUT[[b]]$res.th[ -NUM.idx ]
                }
                # FIXME: if std: what to do?
                if(add.labels) {
                    names(OUT[[b]]$res.th) <- object@pta$vnames$th[[b]]
                }
                if(add.class) {
                    class(OUT[[b]]$res.th) <- c("lavaan.vector", "numeric")
                }
            }

            # slopes
            if(lavmodel@nexo[b] > 0L) {
                if(h1) {
                    OUT[[b]]$res.slopes  <- H1$res.slopes[[b]]
                } else {
                    OUT[[b]]$res.slopes  <- lavsamplestats@res.slopes[[b]]
                }
                # FIXME: if std: what to do? (here: b.z = b * s.x /s.y)
                if(std) {
                    tmp.y <- matrix(sqrt(diag.orig),
                                    nrow(OUT[[b]]$res.slopes),
                                    ncol(OUT[[b]]$res.slopes))
                    tmp.x <- matrix(sqrt(diag(lavsamplestats@cov.x[[b]])),
                                    nrow(OUT[[b]]$res.slopes),
                                    ncol(OUT[[b]]$res.slopes), byrow = TRUE)
                    OUT[[b]]$res.slopes <- OUT[[b]]$res.slopes / tmp.y * tmp.x
                }
                if(add.labels) {
                    rownames(OUT[[b]]$res.slopes) <- ov.names.res[[b]]
                    colnames(OUT[[b]]$res.slopes) <- ov.names.x[[b]]
                }
                if(add.class) {
                    class(OUT[[b]]$res.slopes) <- c("lavaan.matrix", "matrix")
                }
            }

            # cov.x
            if(lavmodel@nexo[b] > 0L) {
                OUT[[b]]$cov.x  <- lavsamplestats@cov.x[[b]]
                if(std) {
                    diag.orig <- diag(OUT[[b]]$cov.x)
                    OUT[[b]]$cov.x <- cov2cor(OUT[[b]]$cov.x)
                }
                if(add.labels) {
                    rownames(OUT[[b]]$cov.x) <- ov.names.x[[b]]
                    colnames(OUT[[b]]$cov.x) <- ov.names.x[[b]]
                }
                if(add.class) {
                    class(OUT[[b]]$cov.x) <-
                        c("lavaan.matrix.symmetric", "matrix")
                }
            }

            # mean.x
            if(lavmodel@nexo[b] > 0L) {
                OUT[[b]]$mean.x <- as.numeric(object@SampleStats@mean.x[[b]])
                if(std) {
                    diag.orig[ diag.orig < .Machine$double.eps ] <- NA
                    OUT[[b]]$mean.x <- OUT[[b]]$mean.x/sqrt(diag.orig)
                }
                if(add.labels) {
                    names(OUT[[b]]$mean.x) <- ov.names.x[[b]]
                }
                if(add.class) {
                    class(OUT[[b]]$mean.x) <- c("lavaan.vector", "numeric")
                }
            }

        } # conditional.x

        # stochastic weights
        if(lavmodel@group.w.free) {
            # to be consistent with the 'implied' values,
            # transform so group.w is the 'log(group.freq)'
            OUT[[b]]$group.w <-
                log(lavsamplestats@group.w[[b]] * lavsamplestats@ntotal)
            if(add.labels) {
                names(OUT[[b]]$group.w) <- "w"
            }
            if(add.class) {
                class(OUT[[b]]$group.w) <- c("lavaan.vector", "numeric")
            }
        }
    }

    if(nblocks == 1L && drop.list.single.group) {
        OUT <- OUT[[1]]
    } else {
        if(object@Data@nlevels == 1L &&
           length(object@Data@group.label) > 0L) {
            names(OUT) <- unlist(object@Data@group.label)
        } else if(object@Data@nlevels > 1L &&
                  length(object@Data@group.label) == 0L) {
            names(OUT) <- object@Data@level.label
        }
    }

    OUT
}


lav_object_inspect_data <- function(object, add.labels = FALSE,
                                    drop.list.single.group = FALSE) {

    G <- object@Data@ngroups
    if(object@Model@conditional.x) {
        OUT <- vector("list", length = G)
        for(g in 1:G) {
            OUT[[g]] <- cbind(object@Data@X[[g]],
                              object@Data@eXo[[g]])
        }
    } else {
        OUT <- object@Data@X
    }

    if(add.labels) {
        for(g in 1:G) {
            if(object@Model@conditional.x) {
                colnames(OUT[[g]]) <- c(object@Data@ov.names[[g]],
                                        object@Data@ov.names.x[[g]])
            } else {
                colnames(OUT[[g]]) <- object@Data@ov.names[[g]]
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

lav_object_inspect_case_idx <- function(object,
                                        drop.list.single.group = FALSE) {

    G <- object@Data@ngroups
    OUT <- object@Data@case.idx

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

    nblocks <- object@Model@nblocks
    OUT <- vector("list", length = nblocks)

    if(is.null(est.std.all)) {
        est.std.all <- lav_standardize_all(object)
    }

    partable <- object@ParTable
    partable$rsquare <- 1.0 - est.std.all
    # no values > 1.0
    partable$rsquare[partable$rsquare > 1.0] <- as.numeric(NA)

    for(b in seq_len(nblocks)) {
        ind.names <- partable$rhs[ which(partable$op == "=~" &
                                         partable$block == b) ]
        eqs.y.names <- partable$lhs[ which(partable$op == "~"  &
                                           partable$block == b) ]
        y.names <- unique( c(ind.names, eqs.y.names) )

        idx <- which(partable$op == "~~" & partable$lhs %in% y.names &
                     partable$rhs == partable$lhs & partable$block == b)
        tmp <- partable$rsquare[idx]

        if(add.labels && length(tmp) > 0L) {
            names(tmp) <- partable$lhs[idx]
        }
        if(add.class) {
            class(tmp) <- c("lavaan.vector", "numeric")
        }

        OUT[[b]] <- tmp
    }

    if(nblocks == 1L && drop.list.single.group) {
        OUT <- OUT[[1]]
    } else {
        if(object@Data@nlevels == 1L &&
           length(object@Data@group.label) > 0L) {
            names(OUT) <- unlist(object@Data@group.label)
        } else if(object@Data@nlevels > 1L &&
                  length(object@Data@group.label) == 0L) {
            names(OUT) <- object@Data@level.label
        }
    }

    OUT
}

# model implied sample stats
lav_object_inspect_implied <- function(object,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {

    nblocks <- object@Model@nblocks
    ov.names <- object@pta$vnames$ov
    ov.names.res <- object@pta$vnames$ov.nox
    ov.names.x   <- object@pta$vnames$ov.x

    # slots
    lavimplied <- object@implied
    lavmodel   <- object@Model

    OUT <- vector("list", length = nblocks)
    for(b in seq_len(nblocks)) {

        if(!lavmodel@conditional.x) {

            # covariance matrix
            OUT[[b]]$cov  <- lavimplied$cov[[b]]
            if(add.labels && !is.null(OUT[[b]]$cov)) {
                rownames(OUT[[b]]$cov) <- colnames(OUT[[b]]$cov) <-
                    ov.names[[b]]
            }
            if(add.class) {
                class(OUT[[b]]$cov) <- c("lavaan.matrix.symmetric", "matrix")
            }

            # mean vector
            if(lavmodel@meanstructure) {
                OUT[[b]]$mean <- as.numeric(lavimplied$mean[[b]])
                if(add.labels) {
                    names(OUT[[b]]$mean) <- ov.names[[b]]
                }
                if(add.class) {
                    class(OUT[[b]]$mean) <- c("lavaan.vector", "numeric")
                }
            }

            # thresholds
            if(lavmodel@categorical) {
                OUT[[b]]$th <- as.numeric(lavimplied$th[[b]])
                if(length(object@Model@num.idx[[b]]) > 0L) {
                    NUM.idx <- which(object@Model@th.idx[[b]] == 0)
                    OUT[[b]]$th <- OUT[[b]]$th[ -NUM.idx ]
                }
                if(add.labels) {
                    names(OUT[[b]]$th) <- object@pta$vnames$th[[b]]
                }
                if(add.class) {
                    class(OUT[[b]]$th) <- c("lavaan.vector", "numeric")
                }
            }
        } # !conditional.x

       else { # if conditional.x = TRUE

            # residual covariance matrix
            OUT[[b]]$res.cov  <- lavimplied$res.cov[[b]]
            if(add.labels) {
                rownames(OUT[[b]]$res.cov) <- colnames(OUT[[b]]$res.cov) <-
                    ov.names.res[[b]]
            }
            if(add.class) {
                class(OUT[[b]]$res.cov) <-
                    c("lavaan.matrix.symmetric", "matrix")
            }

            # intercepts
            if(lavmodel@meanstructure) {
                OUT[[b]]$res.int <- as.numeric(lavimplied$res.int[[b]])
                if(add.labels) {
                    names(OUT[[b]]$res.int) <- ov.names.res[[b]]
                }
                if(add.class) {
                    class(OUT[[b]]$res.int) <- c("lavaan.vector", "numeric")
                }
            }

            # thresholds
            if(lavmodel@categorical) {
                OUT[[b]]$res.th <- as.numeric(lavimplied$res.th[[b]])
                if(length(object@Model@num.idx[[b]]) > 0L) {
                    NUM.idx <- which(object@Model@th.idx[[b]] == 0)
                    OUT[[b]]$res.th <- OUT[[b]]$res.th[ -NUM.idx ]
                }
                if(add.labels) {
                    names(OUT[[b]]$res.th) <- object@pta$vnames$th[[b]]
                }
                if(add.class) {
                    class(OUT[[b]]$res.th) <- c("lavaan.vector", "numeric")
                }
            }

            # slopes
            if(lavmodel@nexo[b] > 0L) {
                OUT[[b]]$res.slopes  <- lavimplied$res.slopes[[b]]
                if(add.labels) {
                    rownames(OUT[[b]]$res.slopes) <- ov.names.res[[b]]
                    colnames(OUT[[b]]$res.slopes) <- ov.names.x[[b]]
                }
                if(add.class) {
                    class(OUT[[b]]$res.slopes) <- c("lavaan.matrix", "matrix")
                }
            }

            # cov.x
            if(lavmodel@nexo[b] > 0L) {
                OUT[[b]]$cov.x  <- object@SampleStats@cov.x[[b]]
                if(add.labels) {
                    rownames(OUT[[b]]$cov.x) <- ov.names.x[[b]]
                    colnames(OUT[[b]]$cov.x) <- ov.names.x[[b]]
                }
                if(add.class) {
                    class(OUT[[b]]$cov.x) <-
                        c("lavaan.matrix.symmetric", "matrix")
                }
            }

            # mean.x
            if(lavmodel@nexo[b] > 0L) {
                OUT[[b]]$mean.x  <- as.numeric(object@SampleStats@mean.x[[b]])
                if(add.labels) {
                    names(OUT[[b]]$mean.x) <- ov.names.x[[b]]
                }
                if(add.class) {
                    class(OUT[[b]]$mean.x) <- c("lavaan.vector", "numeric")
                }
            }


        } # conditional.x

        # stochastic weights
        if(lavmodel@group.w.free) {
            OUT[[b]]$group.w <- lavimplied$group.w[[b]]
            if(add.labels) {
                names(OUT[[b]]$group.w) <- "w" # somewhat redundant
            }
            if(add.class) {
                class(OUT[[b]]$group.w) <- c("lavaan.vector", "numeric")
            }
        }
    }

    if(nblocks == 1L && drop.list.single.group) {
        OUT <- OUT[[1]]
    } else {
        if(object@Data@nlevels == 1L &&
           length(object@Data@group.label) > 0L) {
            names(OUT) <- unlist(object@Data@group.label)
        } else if(object@Data@nlevels > 1L &&
                  length(object@Data@group.label) == 0L) {
            names(OUT) <- object@Data@level.label
        }
    }

    OUT
}


# residuals: _inspect_sampstat - _inspect_implied
lav_object_inspect_residuals <- function(object, h1 = TRUE,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {

    lav_residuals(object, type = "raw", h1 = h1,
                  add.labels = add.labels, add.class = add.class,
                  drop.list.single.group = drop.list.single.group)
}


lav_object_inspect_cov_lv <- function(object, correlation.metric = FALSE,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {

    # compute lv covar
    OUT <- computeVETA(lavmodel = object@Model, remove.dummy.lv = TRUE)

    # nblocks
    nblocks <- length(OUT)

    # cor + labels + class
    for(b in seq_len(nblocks)) {

        if(correlation.metric && nrow(OUT[[b]]) > 1L) {
            # note: cov2cor fails if matrix is empty!
            OUT[[b]] <- cov2cor(OUT[[b]])
        }

        if(add.labels) {
            colnames(OUT[[b]]) <- rownames(OUT[[b]]) <-
                object@pta$vnames$lv[[b]]
        }

        if(add.class) {
            class(OUT[[b]]) <- c("lavaan.matrix.symmetric", "matrix")
        }
    }

    if(nblocks == 1L && drop.list.single.group) {
        OUT <- OUT[[1]]
    } else {
          if(object@Data@nlevels == 1L &&
           length(object@Data@group.label) > 0L) {
            names(OUT) <- unlist(object@Data@group.label)
        } else if(object@Data@nlevels > 1L &&
                  length(object@Data@group.label) == 0L) {
            names(OUT) <- object@Data@level.label
        }
    }

    OUT
}

lav_object_inspect_mean_lv <- function(object,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {

    # compute lv means
    OUT <- computeEETA(lavmodel       = object@Model,
                       lavsamplestats = object@SampleStats,
                       remove.dummy.lv = TRUE)

    # nblocks
    nblocks <- length(OUT)

    # ensure numeric
    OUT <- lapply(OUT, as.numeric)

    # labels + class
    for(b in seq_len(nblocks)) {
        if(add.labels && length(OUT[[b]]) > 0L) {
            names(OUT[[b]]) <- object@pta$vnames$lv.regular[[b]]
        }
        if(add.class) {
            class(OUT[[b]]) <- c("lavaan.vector", "numeric")
        }
    }

    if(nblocks == 1L && drop.list.single.group) {
        OUT <- OUT[[1]]
    } else {
        if(object@Data@nlevels == 1L &&
           length(object@Data@group.label) > 0L) {
            names(OUT) <- unlist(object@Data@group.label)
        } else if(object@Data@nlevels > 1L &&
                  length(object@Data@group.label) == 0L) {
            names(OUT) <- object@Data@level.label
        }
    }

    OUT
}

lav_object_inspect_cov_all <- function(object, correlation.metric = FALSE,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {

    # compute extended model implied covariance matrix (both ov and lv)
    OUT <- computeCOV(lavmodel = object@Model,
                      remove.dummy.lv = TRUE)

    # nblocks
    nblocks <- length(OUT)

    # cor + labels + class
    for(b in seq_len(nblocks)) {

        if(correlation.metric && nrow(OUT[[b]]) > 1L) {
            # note: cov2cor fails if matrix is empty!
            OUT[[b]] <- cov2cor(OUT[[b]])
        }

        if(add.labels) {
            NAMES <- c(object@pta$vnames$ov.model[[b]],
                       object@pta$vnames$lv.regular[[b]])
            colnames(OUT[[b]]) <- rownames(OUT[[b]]) <- NAMES
        }
        if(add.class) {
            class(OUT[[b]]) <- c("lavaan.matrix.symmetric", "matrix")
        }
    }

    if(nblocks == 1L && drop.list.single.group) {
        OUT <- OUT[[1]]
    } else {
        if(object@Data@nlevels == 1L &&
           length(object@Data@group.label) > 0L) {
            names(OUT) <- unlist(object@Data@group.label)
        } else if(object@Data@nlevels > 1L &&
                  length(object@Data@group.label) == 0L) {
            names(OUT) <- object@Data@level.label
        }
    }
    OUT
}


lav_object_inspect_cov_ov <- function(object, correlation.metric = FALSE,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {

    # get model-implied covariance matrix observed
    if(object@Model@conditional.x) {
        OUT <- object@implied$res.cov
    } else {
        OUT <- object@implied$cov
    }

    # nblocks
    nblocks <- length(OUT)

    # cor + labels + class
    for(b in seq_len(nblocks)) {

        if(correlation.metric && nrow(OUT[[b]]) > 1L) {
            # note: cov2cor fails if matrix is empty!
            OUT[[b]] <- cov2cor(OUT[[b]])
        }

        if(add.labels) {
            colnames(OUT[[b]]) <- rownames(OUT[[b]]) <-
                object@pta$vnames$ov.model[[b]]
        }
        if(add.class) {
            class(OUT[[b]]) <- c("lavaan.matrix.symmetric", "matrix")
        }
    }

    if(nblocks == 1L && drop.list.single.group) {
        OUT <- OUT[[1]]
    } else {
        if(object@Data@nlevels == 1L &&
           length(object@Data@group.label) > 0L) {
            names(OUT) <- unlist(object@Data@group.label)
        } else if(object@Data@nlevels > 1L &&
                  length(object@Data@group.label) == 0L) {
            names(OUT) <- object@Data@level.label
        }
    }

    OUT
}

lav_object_inspect_mean_ov <- function(object,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {

    # compute ov means
    if(object@Model@conditional.x) {
        OUT <- object@implied$res.int
    } else {
        OUT <- object@implied$mean
    }

    # nblocks
    nblocks <- length(OUT)

    # make numeric
    OUT <- lapply(OUT, as.numeric)

    # labels + class
    for(b in seq_len(nblocks)) {
        if(add.labels && length(OUT[[b]]) > 0L) {
            names(OUT[[b]]) <- object@pta$vnames$ov.model[[b]]
        }
        if(add.class) {
            class(OUT[[b]]) <- c("lavaan.vector", "numeric")
        }
    }

    if(nblocks == 1L && drop.list.single.group) {
        OUT <- OUT[[1]]
    } else {
        if(object@Data@nlevels == 1L &&
           length(object@Data@group.label) > 0L) {
            names(OUT) <- unlist(object@Data@group.label)
        } else if(object@Data@nlevels > 1L &&
                  length(object@Data@group.label) == 0L) {
            names(OUT) <- object@Data@level.label
        }
    }

    OUT
}

lav_object_inspect_th <- function(object,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {

    # thresholds
    if(object@Model@conditional.x) {
        OUT <- object@implied$res.th
    } else {
        OUT <- object@implied$th
    }

    # nblocks
    nblocks <- length(OUT)

    # make numeric
    OUT <- lapply(OUT, as.numeric)

    # labels + class
    for(b in seq_len(nblocks)) {
        if(length(object@Model@num.idx[[b]]) > 0L) {
            NUM.idx <- which(object@Model@th.idx[[b]] == 0)
            OUT[[b]] <- OUT[[b]][ -NUM.idx ]
        }
        if(add.labels && length(OUT[[b]]) > 0L) {
            names(OUT[[b]]) <- object@pta$vnames$th[[b]]
        }
        if(add.class) {
            class(OUT[[b]]) <- c("lavaan.vector", "numeric")
        }
    }

    if(nblocks == 1L && drop.list.single.group) {
        OUT <- OUT[[1]]
    } else {
        if(object@Data@nlevels == 1L &&
           length(object@Data@group.label) > 0L) {
            names(OUT) <- unlist(object@Data@group.label)
        } else if(object@Data@nlevels > 1L &&
                  length(object@Data@group.label) == 0L) {
            names(OUT) <- object@Data@level.label
        }
    }

    OUT
}

lav_object_inspect_th_idx <- function(object,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {

    # thresholds idx
    OUT <- object@SampleStats@th.idx

    # nblocks
    nblocks <- length(OUT)

    # labels + class
    for(b in seq_len(nblocks)) {
        if(add.labels && length(OUT[[b]]) > 0L) {
            names(OUT[[b]]) <- object@SampleStats@th.names[[b]]
        }
        if(add.class) {
            class(OUT[[b]]) <- c("lavaan.vector", "numeric")
        }
    }

    if(nblocks == 1L && drop.list.single.group) {
        OUT <- OUT[[1]]
    } else {
        if(object@Data@nlevels == 1L &&
           length(object@Data@group.label) > 0L) {
            names(OUT) <- unlist(object@Data@group.label)
        } else if(object@Data@nlevels > 1L &&
                  length(object@Data@group.label) == 0L) {
            names(OUT) <- object@Data@level.label
        }
    }

    OUT
}

lav_object_inspect_vy <- function(object,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {

    # 'unconditional' model-implied variances
    #  - same as diag(Sigma.hat) if all Y are continuous)
    #  - 1.0 (or delta^2) if categorical
    #  - if also Gamma, cov.x is used (only if categorical)

    OUT <- computeVY(lavmodel = object@Model, GLIST = NULL,
                     diagonal.only = TRUE)

    # nblocks
    nblocks <- length(OUT)

    # labels + class
    for(b in seq_len(nblocks)) {
        if(add.labels && length(OUT[[b]]) > 0L) {
            if(object@Model@categorical) {
                 names(OUT[[b]]) <- object@pta$vnames$ov.nox[[b]]
            } else {
                 names(OUT[[b]]) <- object@pta$vnames$ov[[b]]
            }
        }
        if(add.class) {
            class(OUT[[b]]) <- c("lavaan.vector", "numeric")
        }
    }

    if(nblocks == 1L && drop.list.single.group) {
        OUT <- OUT[[1]]
    } else {
        if(object@Data@nlevels == 1L &&
           length(object@Data@group.label) > 0L) {
            names(OUT) <- unlist(object@Data@group.label)
        } else if(object@Data@nlevels > 1L &&
                  length(object@Data@group.label) == 0L) {
            names(OUT) <- object@Data@level.label
        }
    }

    OUT
}


lav_object_inspect_theta <- function(object, correlation.metric = FALSE,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {

    # get residual covariances
    OUT <- computeTHETA(lavmodel = object@Model)

    # nblocks
    nblocks <- length(OUT)

    # labels + class
    for(b in seq_len(nblocks)) {

        if(correlation.metric && nrow(OUT[[b]]) > 0L) {
            if(all(OUT[[b]] == 0)) {
                OUT[[b]] <- OUT[[b]]
            } else {
                OUT[[b]] <- cov2cor(OUT[[b]])
            }
        }

        if(add.labels && length(OUT[[b]]) > 0L) {
            colnames(OUT[[b]]) <- rownames(OUT[[b]]) <-
                object@pta$vnames$ov.model[[b]]
        }

        if(add.class) {
            class(OUT[[b]]) <- c("lavaan.matrix.symmetric", "matrix")
        }
    }

    if(nblocks == 1L && drop.list.single.group) {
        OUT <- OUT[[1]]
    } else {
        if(object@Data@nlevels == 1L &&
           length(object@Data@group.label) > 0L) {
            names(OUT) <- unlist(object@Data@group.label)
        } else if(object@Data@nlevels > 1L &&
                  length(object@Data@group.label) == 0L) {
            names(OUT) <- object@Data@level.label
        }
    }

    OUT
}


lav_object_inspect_missing_coverage <- function(object,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {

    G <- object@Data@ngroups
    OUT <- vector("list", G)

    for(g in 1:G) {
        if(!is.null(object@Data@Mp[[g]])) {
            OUT[[g]] <- object@Data@Mp[[g]]$coverage
        } else {
            nvar <- length(object@Data@ov.names[[g]])
            OUT[[g]] <- matrix(1.0, nvar, nvar)
        }

        if(add.labels && length(OUT[[g]]) > 0L) {
            colnames(OUT[[g]]) <- rownames(OUT[[g]]) <-
                object@pta$vnames$ov.model[[g]]
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

lav_object_inspect_missing_patterns <- function(object,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {

    G <- object@Data@ngroups
    OUT <- vector("list", G)

    for(g in 1:G) {
        if(!is.null(object@Data@Mp[[g]])) {
            OUT[[g]] <- object@Data@Mp[[g]]$pat
        } else {
            nvar <- length(object@Data@ov.names[[g]])
            OUT[[g]] <- matrix(TRUE, 1L, nvar)
            rownames(OUT[[g]]) <- object@Data@nobs[[g]]
        }

        if(add.labels && length(OUT[[g]]) > 0L) {
            colnames(OUT[[g]]) <- object@pta$vnames$ov.model[[g]]
        }

        if(add.class) {
            class(OUT[[g]]) <- c("lavaan.matrix", "matrix")
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

lav_object_inspect_empty_idx <- function(object,
                                         drop.list.single.group = FALSE) {

    G <- object@Data@ngroups

    # get empty idx
    OUT <- vector("list", G)

    for(g in 1:G) {
        if(!is.null(object@Data@Mp[[g]])) {
            OUT[[g]] <- object@Data@Mp[[g]]$empty.idx
        } else {
            OUT[[g]] <- integer(0L)
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


lav_object_inspect_wls_est <- function(object,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {

    OUT <- lav_model_wls_est(object@Model)

    # nblocks
    nblocks <- length(OUT)

    for(b in seq_len(nblocks)) {
        if(add.labels && length(OUT[[b]]) > 0L) {
            #FIXME!!!!
            #names(OUT[[b]]) <- ??
        }

        if(add.class) {
            class(OUT[[b]]) <- c("lavaan.vector", "numeric")
        }
    }

    if(nblocks == 1L && drop.list.single.group) {
        OUT <- OUT[[1]]
    } else {
        if(object@Data@nlevels == 1L &&
           length(object@Data@group.label) > 0L) {
            names(OUT) <- unlist(object@Data@group.label)
        } else if(object@Data@nlevels > 1L &&
                  length(object@Data@group.label) == 0L) {
            names(OUT) <- object@Data@level.label
        }
    }

    OUT
}

lav_object_inspect_wls_obs <- function(object,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {

    OUT <- object@SampleStats@WLS.obs ### FIXME: should be in @h1??

    # nblocks
    nblocks <- length(OUT)

    for(b in seq_len(nblocks)) {
        if(add.labels && length(OUT[[b]]) > 0L) {
            #FIXME!!!!
            #names(OUT[[b]]) <- ??
        }

        if(add.class) {
            class(OUT[[b]]) <- c("lavaan.vector", "numeric")
        }
    }

    if(nblocks == 1L && drop.list.single.group) {
        OUT <- OUT[[1]]
    } else {
        if(object@Data@nlevels == 1L &&
           length(object@Data@group.label) > 0L) {
            names(OUT) <- unlist(object@Data@group.label)
        } else if(object@Data@nlevels > 1L &&
                  length(object@Data@group.label) == 0L) {
            names(OUT) <- object@Data@level.label
        }
    }

    OUT
}

lav_object_inspect_wls_v <- function(object,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {

    #OUT <- lav_model_wls_v(lavmodel       = object@Model,
    #                       lavsamplestats = object@SampleStats,
    #                       structured     = TRUE,
    #                       lavdata        = object@Data)
    # WLS.V == (traditionally) h1 expected information
    OUT <- lav_model_h1_information_expected(lavobject = object)
    # this affects fit measures gfi, agfi, pgfi


    # nblocks
    nblocks <- length(OUT)

    # if estimator == "DWLS" or "ULS", we only stored the diagonal
    # hence, we create a full matrix here
    if(object@Options$estimator %in% c("DWLS", "ULS")) {
        OUT <- lapply(OUT,
            function(x) { nr = NROW(x); diag(x, nrow=nr, ncol=nr) })
    }

    # label + class
    for(b in seq_len(nblocks)) {
        if(add.labels && nrow(OUT[[b]]) > 0L) {
            #FIXME!!!!
            #names(OUT[[b]]) <- ??
        }

        if(add.class) {
            class(OUT[[b]]) <- c("lavaan.matrix", "matrix")
        }
    }

    if(nblocks == 1L && drop.list.single.group) {
        OUT <- OUT[[1]]
    } else {
        if(object@Data@nlevels == 1L &&
           length(object@Data@group.label) > 0L) {
            names(OUT) <- unlist(object@Data@group.label)
        } else if(object@Data@nlevels > 1L &&
                  length(object@Data@group.label) == 0L) {
            names(OUT) <- object@Data@level.label
        }
    }

    OUT
}


lav_object_inspect_sampstat_gamma <- function(object,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {

    if(!is.null(object@SampleStats@NACOV[[1]])) {
        OUT <- object@SampleStats@NACOV
    } else {
        OUT <- lavGamma(object)
    }

    # nblocks
    nblocks <- length(OUT)

    if(nblocks == 1L && drop.list.single.group) {
        OUT <- OUT[[1]]
        if(add.class) {
            class(OUT) <- c("lavaan.matrix.symmetric", "matrix")
        }
    } else {
        if(object@Data@nlevels == 1L && length(object@Data@group.label) > 0L) {
            names(OUT) <- unlist(object@Data@group.label)
            if(add.class) {
                for(g in seq_len(object@Data@ngroups)) {
                    class(OUT[[g]]) <- c("lavaan.matrix.symmetric", "matrix")
                }
            }
        } else if(object@Data@nlevels > 1L &&
                  length(object@Data@group.label) == 0L) {
            names(OUT) <- object@Data@level.label
        }
    }

    OUT
}


lav_object_inspect_gradient <- function(object,
    add.labels = FALSE, add.class = FALSE, logl = FALSE,
    optim = FALSE) {

    lavmodel       <- object@Model
    lavdata        <- object@Data
    lavsamplestats <- object@SampleStats

    if(optim) {
        logl <- FALSE
    }

    if(lavsamplestats@missing.flag ||
       object@Options$estimator == "PML") {
        group.weight <- FALSE
    } else {
        group.weight <- TRUE
    }

    dx <- lav_model_gradient(lavmodel       = lavmodel,
                              GLIST          = NULL,
                              lavsamplestats = lavsamplestats,
                              lavdata        = object@Data,
                              lavcache       = object@Cache,
                              type           = "free",
                              verbose        = FALSE,
                              group.weight   = group.weight)

    # if logl, rescale to get gradient wrt the loglikelihood
    if(logl) {
        if(lavmodel@estimator %in% c("ML")) {
            if(lavdata@nlevels == 1L) {
                # currently, this is just a sign switch
                dx <- -1 * dx
            } else {
                lavpartable <- object@ParTable
                # gradient.log = gradient.obj * (2 * N) / nclusters

                if(lavdata@ngroups == 1L) {
                    N <- lavdata@Lp[[1]]$nclusters[[1]]
                    nclusters <- lavdata@Lp[[1]]$nclusters[[2]]
                    dx <- dx * (2 * N) / nclusters
                } else {
                    for(g in seq_len(lavdata@ngroups)) {
                        N <- lavdata@Lp[[g]]$nclusters[[1]]
                        nclusters <- lavdata@Lp[[g]]$nclusters[[2]]
                        g.idx <-
                          which((lavpartable$group == g)[lavpartable$free > 0L])
                        dx[g.idx] <- dx[g.idx] * (2 * N) / nclusters
                    }
                }
            }
        } else {
            # do nothing (for now)
        }
    }

    # optim?
    if(optim) {
        # 1. scale (note: divide, not multiply!)
        if(!is.null(object@optim$parscale)) {
            dx <- dx / object@optim$parscale
        }

        # 2. pack
        if(lavmodel@eq.constraints) {
            dx <- as.numeric( dx %*% lavmodel@eq.constraints.K )
        }
        # only for PML: divide by N (to speed up convergence)
        if(lavmodel@estimator == "PML") {
            dx <- dx / lavsamplestats@ntotal
        }
    }

    # labels
    if(add.labels) {
        if(optim && lavmodel@eq.constraints) {
            NAMES.all <- lav_partable_labels(object@ParTable, type="free")
            SEQ <- seq_len( length(NAMES.all) )
            pack.SEQ <- as.numeric( (SEQ - lavmodel@eq.constraints.k0) %*%
+                                   lavmodel@eq.constraints.K )
            ok.idx <- which(pack.SEQ %in% SEQ)
            NAMES <- rep("(eq.con)", length(pack.SEQ))
            NAMES[ok.idx] <- NAMES.all[ pack.SEQ[ok.idx] ]
            names(dx) <- NAMES
        } else {
            names(dx) <- lav_partable_labels(object@ParTable, type="free")
        }
    }

    # class
    if(add.class) {
        class(dx) <- c("lavaan.vector", "numeric")
    }

    dx
}

lav_object_inspect_hessian <- function(object,
    add.labels = FALSE, add.class = FALSE) {

    OUT <- lav_model_hessian(lavmodel       = object@Model,
                             lavsamplestats = object@SampleStats,
                             lavdata        = object@Data,
                             lavcache       = object@Cache,
                             lavoptions     = object@Options,
                             group.weight   = TRUE)

    # labels
    if(add.labels) {
        colnames(OUT) <- rownames(OUT) <-
            lav_partable_labels(object@ParTable, type="free")
    }

    # class
    if(add.class) {
        class(OUT) <- c("lavaan.matrix.symmetric", "matrix")
    }

    OUT
}

lav_object_inspect_information <- function(object,
    information = "default", augmented = FALSE, inverted = FALSE,
    add.labels = FALSE, add.class = FALSE) {

    if(information != "default") {
        # override option
        object@Options$information <- information
    }

    OUT <- lav_model_information(lavmodel =  object@Model,
                  lavsamplestats = object@SampleStats,
                  lavdata        = object@Data,
                  lavcache       = object@Cache,
                  lavimplied     = object@implied,
                  lavoptions     = object@Options,
                  extra          = FALSE,
                  augmented      = augmented,
                  inverted       = inverted)

    # labels
    if(add.labels) {
        NAMES <- lav_partable_labels(object@ParTable, type="free")
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
lav_object_inspect_firstorder <- function(object,
    add.labels = FALSE, add.class = FALSE) {

     B0 <- lav_model_information_firstorder(lavmodel =  object@Model,
              lavsamplestats = object@SampleStats,
              lavdata        = object@Data,
              lavcache       = object@Cache,
              lavoptions     = object@Options,
              check.pd       = FALSE,
              augmented      = FALSE,
              inverted       = FALSE)
    attr(B0, "B0.group") <- NULL
    OUT <- B0

    # labels
    if(add.labels) {
        colnames(OUT) <- rownames(OUT) <-
            lav_partable_labels(object@ParTable, type="free")
    }

    # class
    if(add.class) {
        class(OUT) <- c("lavaan.matrix.symmetric", "matrix")
    }

    OUT
}

lav_object_inspect_vcov <- function(object, standardized = FALSE,
    type = "std.all", free.only = TRUE,
    add.labels = FALSE, add.class = FALSE, remove.duplicated = FALSE) {

    npar <- max(object@ParTable$free)
    if(object@optim$npar == 0) {
        OUT <- matrix(0,0,0)
    } else {
        # check if we already have it
        tmp <- try(slot(object, "vcov"), silent = TRUE)
        if(!inherits(tmp, "try-error") && !is.null(object@vcov$vcov)) {
            OUT <- object@vcov$vcov
        } else {
        # compute it again
            OUT <- lav_model_vcov(lavmodel       = object@Model,
                                  lavsamplestats = object@SampleStats,
                                  lavoptions     = object@Options,
                                  lavdata        = object@Data,
                                  lavcache       = object@Cache,
                                  lavimplied     = object@implied,
                                  lavh1          = object@h1
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
            JAC <- try(lav_func_jacobian_complex(func = lav_standardize_lv_x,
                x = object@optim$x, lavobject = object), silent = TRUE)
            if(inherits(JAC, "try-error")) { # eg. pnorm()
                JAC <- lav_func_jacobian_simple(func = lav_standardize_lv_x,
                    x = object@optim$x, lavobject = object)
            }
        } else if(type == "std.all") {
            JAC <- try(lav_func_jacobian_complex(func = lav_standardize_all_x,
                x = object@optim$x, object = object), silent = TRUE)
            if(inherits(JAC, "try-error")) { # eg. pnorm()
                JAC <- lav_func_jacobian_simple(func = lav_standardize_all_x,
                    x = object@optim$x, lavobject = object)
            }
        } else if(type == "std.nox") {
            JAC <-
                try(lav_func_jacobian_complex(func = lav_standardize_all_nox_x,
                    x = object@optim$x, lavobject = object), silent = TRUE)
            if(inherits(JAC, "try-error")) { # eg. pnorm()
                JAC <-
                    lav_func_jacobian_simple(func = lav_standardize_all_nox_x,
                        x = object@optim$x, lavobject = object)
            }
        }

        # JAC contains *all* parameters in the parameter table
        if(free.only) {
            free.idx <- which(object@ParTable$free > 0L)
            JAC <- JAC[free.idx,, drop = FALSE]
        }
        OUT <- JAC %*% OUT %*% t(JAC)
    }

    # labels
    if(add.labels) {
        colnames(OUT) <- rownames(OUT) <-
            lav_partable_labels(object@ParTable, type="free")
    }

    # alias?
    if(remove.duplicated && object@Model@eq.constraints) {
        simple.flag <- lav_constraints_check_simple(object@Model)
        if(simple.flag) {
            LAB <- lav_partable_labels(object@ParTable, type="free")
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

lav_object_inspect_vcov_def <- function(object, joint = FALSE,
    standardized = FALSE, type = "std.all",
    add.labels = FALSE, add.class = FALSE) {

    lavmodel    <- object@Model
    lavpartable <- object@ParTable
    free.idx <- which(lavpartable$free > 0L)
    def.idx <- which(lavpartable$op == ":=")
    joint.idx <- c(free.idx, def.idx)

    if(!joint && length(def.idx) == 0L) {
        return( matrix(0,0,0) )
    } else if(joint && length(joint.idx) == 0L) {
        return( matrix(0,0,0) )
    }

    if(standardized) {
        # compute VCOV for "free" parameters only
        VCOV <- lav_object_inspect_vcov(object,
                                        standardized = TRUE,
                                        type = type, free.only = FALSE,
                                        add.labels = FALSE, add.class = FALSE)
        if(joint) {
            OUT <- VCOV[joint.idx, joint.idx, drop = FALSE]
        } else {
            OUT <- VCOV[def.idx, def.idx, drop = FALSE]
        }
    } else {

        # get free parameters
        x <- lav_model_get_parameters(lavmodel, type = "free")

        # bootstrap or not?
        if(!is.null(object@boot$coef)) {
            BOOT <- object@boot$coef
            BOOT.def <- apply(BOOT, 1L, lavmodel@def.function)
            if(length(def.idx) == 1L) {
                BOOT.def <- as.matrix(BOOT.def)
            } else {
                BOOT.def <- t(BOOT.def)
            }
            OUT <- cov(BOOT.def)
        } else {
            # VCOV
            VCOV <- lav_object_inspect_vcov(object,
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
            if(joint) {
                JAC2 <- rbind( diag(nrow = ncol(JAC)), JAC )
                OUT <- JAC2 %*% VCOV %*% t(JAC2)
            } else {
                OUT <- JAC %*% VCOV %*% t(JAC)
            }
        }
    }

    # labels
    if(add.labels) {
        if(joint) {
            LHS.names <- lavpartable$lhs[joint.idx]
        } else {
            LHS.names <- lavpartable$lhs[def.idx]
        }
        colnames(OUT) <- rownames(OUT) <- LHS.names
    }

    # class
    if(add.class) {
        class(OUT) <- c("lavaan.matrix.symmetric", "matrix")
    }

    OUT
}

lav_object_inspect_UGamma <- function(object,
    add.labels = FALSE, add.class = FALSE) {

    out <- lav_test_satorra_bentler(lavobject     = object,
                                    method        = "original",
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

lav_object_inspect_UfromUGamma <- function(object,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {

    out <- lav_test_satorra_bentler(lavobject     = object,
                                    method        = "original",
                                    return.u      = TRUE)
    OUT <- out$UfromUGamma

    # labels
    #if(add.labels) {
       # colnames(OUT) <- rownames(OUT) <-
    #}

    lavmodel <- object@Model
    lavdata  <- object@Data

    if(lavmodel@ngroups == 1L && drop.list.single.group) {
        OUT <- OUT[[1]]
        # class
        if(add.class) {
            class(OUT) <- c("lavaan.matrix.symmetric", "matrix")
        }
    } else {
        if(length(lavdata@group.label) > 0L) {
            names(OUT) <- unlist(lavdata@group.label)
        }
        if(add.class) {
            for(g in seq_len(lavmodel@ngroups)) {
                class(OUT[[g]]) <- c("lavaan.matrix.symmetric", "matrix")
            }
        }
    }

    OUT
}


# Delta (jacobian: d samplestats / d free_parameters)
lav_object_inspect_delta <- function(object,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {

    lavmodel    <- object@Model
    lavdata     <- object@Data
    lavpartable <- object@ParTable
    lavpta      <- object@pta

    OUT <- lav_object_inspect_delta_internal(lavmodel = lavmodel,
        lavdata = lavdata, lavpartable = lavpartable, lavpta = lavpta,
        add.labels = add.labels, add.class = add.class,
        drop.list.single.group = drop.list.single.group)

    OUT
}

lav_object_inspect_delta_internal <- function(lavmodel = NULL, lavdata  = NULL,
    lavpartable = NULL, lavpta = NULL,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {

    OUT <- computeDelta(lavmodel)

    categorical    <- lavmodel@categorical
    conditional.x  <- lavmodel@conditional.x
    group.w.free   <- lavmodel@group.w.free
    nvar           <- lavmodel@nvar
    num.idx        <- lavmodel@num.idx
    th.idx         <- lavmodel@th.idx
    nexo           <- lavmodel@nexo
    nblocks        <- lavmodel@nblocks

    if(add.labels) {
        PNAMES <- lav_partable_labels(lavpartable, type="free")

        if(is.null(lavpta)) {
            lavpta <- lav_partable_attributes(lavpartable)
        }

        if(lavdata@nlevels > 1L) {
            # store names per block, rbind later
            NAMES <- vector("list", length = nblocks)
        }

        for(g in 1:nblocks) {

            if(lavdata@nlevels == 1L) {
                colnames(OUT[[g]]) <- PNAMES
            }

            if(conditional.x) {
                ov.names <- lavpta$vnames$ov.nox[[g]]
            } else {
                ov.names <- lavpta$vnames$ov[[g]]
            }
            ov.names.x <- lavpta$vnames$ov.x[[g]]
            nvar <- length(ov.names)


            names.cov <- names.cor <- names.var <- character(0L)
            names.mu <- names.pi <- names.th <- character(0L)
            names.gw <- character(0L)

            # Sigma
            # - if continuous: vech(Sigma)
            # - if categorical: first numeric variances, then
            #tmp <- apply(expand.grid(ov.names, ov.names), 1L,
            #             paste, collapse = "~~")

            #if(categorical) {
            #    names.cor <- tmp[lav_matrix_vech_idx(nvar, diagonal = FALSE)]
            #    names.var <- tmp[lav_matrix_diag_idx(nvar)[num.idx[[g]]]]
            #} else {
            #    names.cov <- tmp[lav_matrix_vech_idx(nvar, diagonal = TRUE)]
            #}

            # NOTE: in 0.6-1, we use the same order, but 'label' in row-wise
            # format (eg x1 ~~ x2 instead of x2 ~~ x1)
            tmp <- matrix(apply(expand.grid(ov.names, ov.names), 1L,
                          paste, collapse = "~~"), nrow = nvar)
            if(categorical) {
                names.cor <- lav_matrix_vechru(tmp, diagonal = FALSE)
                names.var <- diag(tmp)[num.idx[[g]]]
            } else {
                names.cov <- lav_matrix_vechru(tmp, diagonal = TRUE)
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
                names.th <- lavpta$vnames$th[[g]]
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

            if(lavdata@nlevels == 1L) {
                rownames(OUT[[g]]) <- c(names.gw,
                                        names.th, names.mu,
                                        names.pi,
                                        names.cov, names.var, names.cor)
                # class
                if(add.class) {
                    class(OUT[[g]]) <- c("lavaan.matrix", "matrix")
                }
            } else {
                NAMES[[g]] <- c(names.gw,
                                names.th, names.mu,
                                names.pi,
                                names.cov, names.var, names.cor)
            }

        } # g
    } # labels

    # multilevel
    if(lavmodel@multilevel) {
        for(g in 1:lavmodel@ngroups) {
            if(add.labels) {
                colnames(OUT[[g]]) <- PNAMES
                rownames(OUT[[g]]) <- c(NAMES[[(g-1)*2 + 1]],
                                        NAMES[[(g-1)*2 + 2]] )
            }
            if(add.class) {
                class(OUT[[g]]) <- c("lavaan.matrix", "matrix")
            }
        }
    }

    if(lavmodel@ngroups == 1L && drop.list.single.group) {
        OUT <- OUT[[1]]
    } else {
        if(length(lavdata@group.label) > 0L) {
            names(OUT) <- unlist(lavdata@group.label)
        }
    }

    OUT
}

lav_object_inspect_zero_cell_tables <- function(object,
            add.labels = FALSE, add.class = FALSE,
            drop.list.single.group = FALSE) {

    # categorical?
    if(!object@Model@categorical) {
        warning("lavaan WARNING: no categorical variables in fitted model")
        return(invisible(list()))
    }

    lavdata <- object@Data

    # create 2-way tables
    TABLE <- lavTables(object, dimension = 2L, output = "data.frame",
                       statistic = NULL)

    # select tables with empty cells
    empty.id <- TABLE$id[which(TABLE$obs.freq == 0)]


    if(length(empty.id) == 0L) {
        # only when lavInspect() is used, give message
        if(add.class) {
            cat("(There are no tables with empty cells for this fitted model)\n")
        }
        return(invisible(list()))
    } else {
        OUT <- lav_tables_cells_format(TABLE[TABLE$id %in% empty.id,],
                   lavdata = lavdata,
                   drop.list.single.group = drop.list.single.group)
    }

    OUT
}

lav_object_inspect_coef <- function(object, type = "free",
                                    add.labels = FALSE, add.class = FALSE) {

    if(type == "user" || type == "all") {
        type <- "user"
        idx <- 1:length( object@ParTable$lhs )
    } else if(type == "free") {
        idx <- which(object@ParTable$free > 0L & !duplicated(object@ParTable$free))
    } else {
        stop("lavaan ERROR: argument `type' must be one of free or user")
    }
    EST <- lav_object_inspect_est(object)
    cof <- EST[idx]

    # labels?
    if(add.labels) {
        names(cof) <- lav_partable_labels(object@ParTable, type = type)
    }

    # class
    if(add.class) {
        class(cof) <- c("lavaan.vector", "numeric")
    }

    cof
}

lav_object_inspect_npar <- function(object, type = "free") {

    if(type == "free") {
        npar <- sum(object@ParTable$free > 0L &
                    !duplicated(object@ParTable$free))
    } else {
        npar <- length(object@ParTable$lhs)
    }

    npar
}

lav_object_inspect_icc <- function(object, add.labels = FALSE,
                                   add.class = FALSE,
                                   drop.list.single.group = FALSE) {

    lavdata <- object@Data
    G <- lavdata@ngroups
    OUT <- vector("list", G)

    # multilevel?
    if(lavdata@nlevels == 1L) {
        stop("lavaan ERROR: intraclass correlation only available for clustered data")
    }

    if(length(object@h1) == 0L) {
        stop("lavaan ERROR: h1 slot is of available; refit with h1 = TRUE")
    }

    # implied statistics
    implied <- object@h1$implied

    for(g in 1:G) {
        Sigma.W <- implied$cov[[  (g-1)*lavdata@nlevels + 1 ]]
        Sigma.B <- implied$cov[[  (g-1)*lavdata@nlevels + 2 ]]

        W.diag <- diag(Sigma.W)
        B.diag <- diag(Sigma.B)

        OUT[[g]] <- numeric(length(W.diag))

        ov.names.l <- lavdata@ov.names.l[[g]]
        w.idx <- which(ov.names.l[[1]] %in% ov.names.l[[2]])
        w.names <- ov.names.l[[1]][w.idx]
        b.idx <- match(w.names, ov.names.l[[2]])

        OUT[[g]][w.idx] <- B.diag[b.idx]/(W.diag[w.idx] + B.diag[b.idx])

        # label
        if(add.labels) {
            names(OUT[[g]]) <- ov.names.l[[1]]
        }

        # class
        if(add.class) {
            class(OUT[[g]]) <- c("lavaan.vector", "numeric")
        }
    } # g

    if(G == 1L && drop.list.single.group) {
        OUT <- OUT[[1]]
    } else {
        if(length(object@Data@group.label) > 0L) {
            names(OUT) <- unlist(object@Data@group.label)
        }
    }

    OUT
}

lav_object_inspect_ranef <- function(object, add.labels = FALSE,
                                     add.class = FALSE,
                                     drop.list.single.group = FALSE) {

    lavdata <- object@Data
    lavsamplestats <- object@SampleStats

    G <- lavdata@ngroups
    OUT <- vector("list", G)

    # multilevel?
    if(lavdata@nlevels == 1L) {
        stop("lavaan ERROR: random effects only available for clustered data (in the long format)")
    }

    # implied statistics
    lavimplied <- object@implied

    for(g in 1:G) {

        Lp <- lavdata@Lp[[g]]
        YLp <- lavsamplestats@YLp[[g]]

        # implied for this group
        group.idx <- (g - 1)*lavdata@nlevels + seq_len(lavdata@nlevels)
        implied.group <- lapply(lavimplied, function(x) x[group.idx])

        # random effects (=random intercepts or cluster means)
        out <- lav_mvnorm_cluster_implied22l(Lp = Lp, implied = implied.group)
        MB.j <- lav_mvnorm_cluster_em_estep_ranef(YLp = YLp, Lp = Lp,
                        sigma.w = out$sigma.w, sigma.b = out$sigma.b,
                        sigma.zz = out$sigma.zz, sigma.yz = out$sigma.yz,
                        mu.z = out$mu.z, mu.w = out$mu.w, mu.b = out$mu.b,
                        se = FALSE)
        OUT[[g]] <- MB.j

        ov.names.l <- lavdata@ov.names.l[[g]]

        # label
        if(add.labels) {
            colnames(OUT[[g]]) <- ov.names.l[[1]]
        }

        # class
        if(add.class) {
            class(OUT[[g]]) <- c("lavaan.matrix", "matrix")

        }
    } # g

    if(G == 1L && drop.list.single.group) {
        OUT <- OUT[[1]]
    } else {
        if(length(object@Data@group.label) > 0L) {
            names(OUT) <- unlist(object@Data@group.label)
        }
    }

    OUT
}
