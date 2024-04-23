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
lavInspect.lavaan <- function(object,                                # nolint
                              what                   = "free",
                              add.labels             = TRUE,
                              add.class              = TRUE,
                              list.by.group          = TRUE,
                              drop.list.single.group = TRUE) {
  # object must inherit from class lavaan
  stopifnot(inherits(object, "lavaan"))

  # store partable with pta in object to use cache in called functions
  object@ParTable <- lav_partable_set_cache(object@ParTable, object@pta)

  # old (<0.6) object?
  if (!.hasSlot(object@Data, "block.label")) {
    object@Data@block.label <- object@Data@group.label
  }

  # only a single argument
  if (length(what) > 1) {
    lav_msg_stop(gettextf("argument %s cannot have more than one element",
                          "what"))
  }

  # be case insensitive
  what <- tolower(what)

  #### model matrices, with different contents ####
  if (what == "free") {
    lav_object_inspect_modelmatrices(object, what = "free",
      type = "free", add.labels = add.labels, add.class = add.class,
      list.by.group = list.by.group,
      drop.list.single.group = drop.list.single.group)
  } else if (what == "impute" ||
    what == "imputed") { # just to ease the transition for semTools!
    object@imputed
  } else if (what == "partable" || what == "user") {
    lav_object_inspect_modelmatrices(object, what = "free",
      type = "partable", add.labels = add.labels, add.class = add.class,
      list.by.group = list.by.group,
      drop.list.single.group = drop.list.single.group)
  } else if (what == "se" ||
    what == "std.err" ||
    what == "standard.errors") {
    lav_object_inspect_modelmatrices(object, what = "se",
      add.labels = add.labels, add.class = add.class,
      list.by.group = list.by.group,
      drop.list.single.group = drop.list.single.group)
  } else if (what == "se.std" ||
    what == "std.se") {
    lav_object_inspect_modelmatrices(object, what = "std.se",
      add.labels = add.labels, add.class = add.class,
      list.by.group = list.by.group,
      drop.list.single.group = drop.list.single.group)
  } else if (what == "start" || what == "starting.values") {
    lav_object_inspect_modelmatrices(object, what = "start",
      add.labels = add.labels, add.class = add.class,
      list.by.group = list.by.group,
      drop.list.single.group = drop.list.single.group)
  } else if (what == "est"  || what == "estimates" ||
    what == "x") {
    lav_object_inspect_modelmatrices(object, what = "est",
      add.labels = add.labels, add.class = add.class,
      list.by.group = list.by.group,
      drop.list.single.group = drop.list.single.group)
  } else if (what == "est.unrotated") {
    lav_object_inspect_modelmatrices(object, what = "est.unrotated",
      add.labels = add.labels, add.class = add.class,
      list.by.group = list.by.group,
      drop.list.single.group = drop.list.single.group)
  } else if (what == "dx.free") {
    lav_object_inspect_modelmatrices(object, what = "dx.free",
      add.labels = add.labels, add.class = add.class,
      list.by.group = list.by.group,
      drop.list.single.group = drop.list.single.group)
  } else if (what == "dx.all") {
    lav_object_inspect_modelmatrices(object, what = "dx.all",
      add.labels = add.labels, add.class = add.class,
      list.by.group = list.by.group,
      drop.list.single.group = drop.list.single.group)
  } else if (what == "std" || what == "std.all" ||
    what == "est.std" || what == "std.est" ||
    what == "standardized") {
    lav_object_inspect_modelmatrices(object, what = "std.all",
      add.labels = add.labels, add.class = add.class,
      list.by.group = list.by.group,
      drop.list.single.group = drop.list.single.group)
  } else if (what == "std.lv") {
    lav_object_inspect_modelmatrices(object, what = "std.lv",
      add.labels = add.labels, add.class = add.class,
      list.by.group = list.by.group,
      drop.list.single.group = drop.list.single.group)
  } else if (what == "std.nox") {
    lav_object_inspect_modelmatrices(object, what = "std.nox",
      add.labels = add.labels, add.class = add.class,
      list.by.group = list.by.group,
      drop.list.single.group = drop.list.single.group)


    #### parameter table ####
  } else if (what == "list") {
    parTable(object)

    #### bootstrap coef ####
  } else if (what %in% c("boot", "bootstrap", "boot.coef", "coef.boot")) {
    lav_object_inspect_boot(object, add.labels = add.labels,
      add.class = add.class)

    #### fit indices ####
  } else if (what == "fit" ||
    what == "fitmeasures" ||
    what == "fit.measures" ||
    what == "fit.indices") {
    fitMeasures(object)

    #### baseline model ####
  } else if (what == "baseline.partable") {
    out <- as.data.frame(object@baseline$partable, stringsAsFactors = FALSE)
    if (add.class) {
      class(out) <- c("lavaan.data.frame", "data.frame")
    }
    return(out)
  } else if (what == "baseline.test") {
    object@baseline$test

    #### modification indices ####
  } else if (what == "mi" ||
    what == "modindices" ||
    what == "modification.indices") {
    modificationIndices(object)


    #### sample statistics #####
  } else if (what == "obs" ||
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
  } else if (what == "obs.std" ||
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
  } else if (what == "h1" || what == "missing.h1" || what == "sampstat.h1") {
    lav_object_inspect_sampstat(object, h1 = TRUE,
      add.labels = add.labels, add.class = add.class,
      drop.list.single.group = drop.list.single.group)

    #### wls.est - wls.obs - wls.v ####
  } else if (what == "wls.est") {
    lav_object_inspect_wls_est(object,
      add.labels = add.labels, add.class = add.class,
      drop.list.single.group = drop.list.single.group)
  } else if (what == "wls.obs") {
    lav_object_inspect_wls_obs(object,
      add.labels = add.labels, add.class = add.class,
      drop.list.single.group = drop.list.single.group)
  } else if (what == "wls.v") {
    lav_object_inspect_wls_v(object,
      add.labels = add.labels, add.class = add.class,
      drop.list.single.group = drop.list.single.group)



    #### data + missingness ####
  } else if (what == "data") {
    lav_object_inspect_data(object, add.labels = add.labels,
      drop.list.single.group = drop.list.single.group)
  } else if (what == "case.idx") {
    lav_object_inspect_case_idx(object,
      drop.list.single.group = drop.list.single.group)
  } else if (what == "ngroups") {
    object@Data@ngroups
  } else if (what == "group") {
    object@Data@group
  } else if (what == "cluster") {
    object@Data@cluster
  } else if (what == "nlevels") {
    object@Data@nlevels
  } else if (what == "nclusters") {
    lav_object_inspect_cluster_info(object, level = 2L,
      what = "nclusters",
      drop.list.single.group = drop.list.single.group)
  } else if (what == "ncluster.size") {
    lav_object_inspect_cluster_info(object, level = 2L,
      what = "ncluster.size",
      drop.list.single.group = drop.list.single.group)
  } else if (what == "cluster.size") {
    lav_object_inspect_cluster_info(object, level = 2L,
      what = "cluster.size",
      drop.list.single.group = drop.list.single.group)
  } else if (what == "cluster.id") {
    lav_object_inspect_cluster_info(object, level = 2L,
      what = "cluster.id",
      drop.list.single.group = drop.list.single.group)
  } else if (what == "cluster.idx") {
    lav_object_inspect_cluster_info(object, level = 2L,
      what = "cluster.idx",
      drop.list.single.group = drop.list.single.group)
  } else if (what == "cluster.label") {
    lav_object_inspect_cluster_info(object, level = 2L,
      what = "cluster.label",
      drop.list.single.group = drop.list.single.group)
  } else if (what == "cluster.sizes") {
    lav_object_inspect_cluster_info(object, level = 2L,
      what = "cluster.sizes",
      drop.list.single.group = drop.list.single.group)
  } else if (what == "average.cluster.size") {
    lav_object_inspect_cluster_info(object, level = 2L,
      what = "average.cluster.size",
      drop.list.single.group = drop.list.single.group)
  } else if (what == "ordered") {
    object@Data@ordered
  } else if (what == "group.label") {
    object@Data@group.label
  } else if (what == "level.label") {
    object@Data@level.label
  } else if (what == "nobs") {
    unlist(object@Data@nobs)
  } else if (what == "norig") {
    unlist(object@Data@norig)
  } else if (what == "ntotal") {
    sum(unlist(object@Data@nobs))
  } else if (what == "coverage") {
    lav_object_inspect_missing_coverage(object,
      add.labels = add.labels, add.class = add.class,
      drop.list.single.group = drop.list.single.group)
  } else if (what %in% c("patterns", "pattern")) {
    lav_object_inspect_missing_patterns(object,
      add.labels = add.labels, add.class = add.class,
      drop.list.single.group = drop.list.single.group)
  } else if (what == "empty.idx") {
    lav_object_inspect_empty_idx(object,
      drop.list.single.group = drop.list.single.group)


    #### rsquare ####
  } else if (what == "rsquare" || what == "r-square" || what == "r2") {
    lav_object_inspect_rsquare(object,
      add.labels = add.labels, add.class = add.class,
      drop.list.single.group = drop.list.single.group)


    #### model-implied sample statistics ####
  } else if (what == "implied" || what == "fitted" ||
    what == "expected" || what == "exp") {
    lav_object_inspect_implied(object,
      add.labels = add.labels, add.class = add.class,
      drop.list.single.group = drop.list.single.group)
  } else if (what == "resid" || what == "res" || what == "residual" ||
    what == "residuals") {
    lav_object_inspect_residuals(object, h1 = TRUE,
      add.labels = add.labels, add.class = add.class,
      drop.list.single.group = drop.list.single.group)
  } else if (what == "cov.lv" || what == "veta") {
    lav_object_inspect_cov_lv(object,
      correlation.metric = FALSE,
      add.labels = add.labels, add.class = add.class,
      drop.list.single.group = drop.list.single.group)
  } else if (what == "cor.lv") {
    lav_object_inspect_cov_lv(object,
      correlation.metric = TRUE,
      add.labels = add.labels, add.class = add.class,
      drop.list.single.group = drop.list.single.group)
  } else if (what == "mean.lv" || what == "eeta") {
    lav_object_inspect_mean_lv(object,
      add.labels = add.labels, add.class = add.class,
      drop.list.single.group = drop.list.single.group)
  } else if (what == "cov.all") {
    lav_object_inspect_cov_all(object,
      correlation.metric = FALSE,
      add.labels = add.labels, add.class = add.class,
      drop.list.single.group = drop.list.single.group)
  } else if (what == "cor.all") {
    lav_object_inspect_cov_all(object,
      correlation.metric = TRUE,
      add.labels = add.labels, add.class = add.class,
      drop.list.single.group = drop.list.single.group)
  } else if (what == "cov.ov" || what == "sigma" || what == "sigma.hat") {
    lav_object_inspect_cov_ov(object,
      correlation.metric = FALSE,
      add.labels = add.labels, add.class = add.class,
      drop.list.single.group = drop.list.single.group)
  } else if (what == "cor.ov") {
    lav_object_inspect_cov_ov(object,
      correlation.metric = TRUE,
      add.labels = add.labels, add.class = add.class,
      drop.list.single.group = drop.list.single.group)
  } else if (what == "mean.ov" || what == "mu" || what == "mu.hat") {
    lav_object_inspect_mean_ov(object,
      add.labels = add.labels, add.class = add.class,
      drop.list.single.group = drop.list.single.group)
  } else if (what == "th" || what == "thresholds") {
    lav_object_inspect_th(object,
      add.labels = add.labels, add.class = add.class,
      drop.list.single.group = drop.list.single.group)
  } else if (what == "th.idx") {
    lav_object_inspect_th_idx(object,
      add.labels = add.labels, add.class = add.class,
      drop.list.single.group = drop.list.single.group)
  } else if (what == "vy") {
    lav_object_inspect_vy(object,
      add.labels = add.labels, add.class = add.class,
      drop.list.single.group = drop.list.single.group)


    #### specific model matrices? ####
  } else if (what == "theta" || what == "theta.cov") {
    lav_object_inspect_theta(object,  correlation.metric = FALSE,
      add.labels = add.labels, add.class = add.class,
      drop.list.single.group = drop.list.single.group)
  } else if (what == "theta.cor") {
    lav_object_inspect_theta(object,  correlation.metric = TRUE,
      add.labels = add.labels, add.class = add.class,
      drop.list.single.group = drop.list.single.group)

    #### (squared) Mahalanobis distances ####
  } else if (what == "mdist2.fs") {
    lav_object_inspect_mdist2(object, type = "lv", squared = TRUE,
      add.labels = add.labels, add.class = add.class,
      drop.list.single.group = drop.list.single.group)
  } else if (what == "mdist2.resid") {
    lav_object_inspect_mdist2(object, type = "resid", squared = TRUE,
      add.labels = add.labels, add.class = add.class,
      drop.list.single.group = drop.list.single.group)
  } else if (what == "mdist.fs") {
    lav_object_inspect_mdist2(object, type = "lv", squared = FALSE,
      add.labels = add.labels, add.class = add.class,
      drop.list.single.group = drop.list.single.group)
  } else if (what == "mdist.resid") {
    lav_object_inspect_mdist2(object, type = "resid", squared = FALSE,
      add.labels = add.labels, add.class = add.class,
      drop.list.single.group = drop.list.single.group)

    #### convergence, meanstructure, categorical ####
  } else if (what == "converged") {
    object@optim$converged
  } else if (what == "iterations" ||
    what == "iter" ||
    what == "niter") {
    object@optim$iterations
  } else if (what == "meanstructure") {
    object@Model@meanstructure
  } else if (what == "categorical") {
    object@Model@categorical
  } else if (what == "fixed.x") {
    object@Model@fixed.x
  } else if (what == "parameterization") {
    object@Model@parameterization
  } else if (what == "npar") {
    lav_object_inspect_npar(object, type = "free")
  } else if (what == "coef") {
    # this breaks simsem and semTools -- 0.6-1
    # lav_object_inspect_coef(object, type = "free",
    #                        add.labels = add.labels, add.class = add.class)
    lav_object_inspect_modelmatrices(object, what = "est",
      type = "free", add.labels = add.labels, add.class = add.class,
      list.by.group = list.by.group,
      drop.list.single.group = drop.list.single.group)


    #### NACOV samplestats ####
  } else if (what == "gamma") {
    lav_object_inspect_sampstat_gamma(object,
      add.labels = add.labels, add.class = add.class,
      drop.list.single.group = drop.list.single.group)


    #### gradient, Hessian, information, first.order, vcov ####
  } else if (what == "gradient") {
    lav_object_inspect_gradient(object,
      add.labels = add.labels, add.class = add.class, logl = FALSE)
  } else if (what == "gradient.logl") {
    lav_object_inspect_gradient(object,
      add.labels = add.labels, add.class = add.class, logl = TRUE)
  } else if (what == "optim.gradient") {
    lav_object_inspect_gradient(object,
      add.labels = add.labels, add.class = add.class, optim = TRUE)
  } else if (what == "hessian") {
    lav_object_inspect_hessian(object,
      add.labels = add.labels, add.class = add.class)

  } else if (what == "information") {
    lav_object_inspect_information(object, information = "default",
      augmented = FALSE, inverted = FALSE,
      add.labels = add.labels, add.class = add.class)
  } else if (what == "information.expected") {
    lav_object_inspect_information(object, information = "expected",
      augmented = FALSE, inverted = FALSE,
      add.labels = add.labels, add.class = add.class)
  } else if (what == "information.observed") {
    lav_object_inspect_information(object, information = "observed",
      augmented = FALSE, inverted = FALSE,
      add.labels = add.labels, add.class = add.class)
  } else if (what == "information.first.order" ||
    what == "information.firstorder"  ||
    what == "first.order") {
    lav_object_inspect_information(object, information = "first.order",
      augmented = FALSE, inverted = FALSE,
      add.labels = add.labels, add.class = add.class)

  } else if (what == "augmented.information") {
    lav_object_inspect_information(object, information = "default",
      augmented = TRUE, inverted = FALSE,
      add.labels = add.labels, add.class = add.class)
  } else if (what == "augmented.information.expected") {
    lav_object_inspect_information(object, information = "expected",
      augmented = TRUE, inverted = FALSE,
      add.labels = add.labels, add.class = add.class)
  } else if (what == "augmented.information.observed") {
    lav_object_inspect_information(object, information = "observed",
      augmented = TRUE, inverted = FALSE,
      add.labels = add.labels, add.class = add.class)
  } else if (what == "augmented.information.first.order" ||
    what == "augmented.first.order") {
    lav_object_inspect_information(object, information = "first.order",
      augmented = TRUE, inverted = FALSE,
      add.labels = add.labels, add.class = add.class)

  } else if (what == "inverted.information") {
    lav_object_inspect_information(object, information = "default",
      augmented = TRUE, inverted = TRUE,
      add.labels = add.labels, add.class = add.class)
  } else if (what == "inverted.information.expected") {
    lav_object_inspect_information(object, information = "expected",
      augmented = TRUE, inverted = TRUE,
      add.labels = add.labels, add.class = add.class)
  } else if (what == "inverted.information.observed") {
    lav_object_inspect_information(object, information = "observed",
      augmented = TRUE, inverted = TRUE,
      add.labels = add.labels, add.class = add.class)
  } else if (what == "inverted.information.first.order" ||
    what == "inverted.first.order") {
    lav_object_inspect_information(object, information = "first.order",
      augmented = TRUE, inverted = TRUE,
      add.labels = add.labels, add.class = add.class)

  } else if (what == "h1.information") {
    lav_object_inspect_h1_information(object, information = "default",
      h1.information = "default", inverted = FALSE,
      add.labels = add.labels, add.class = add.class,
      drop.list.single.group = drop.list.single.group)
  } else if (what == "h1.information.expected") {
    lav_object_inspect_h1_information(object, information = "expected",
      h1.information = "default", inverted = FALSE,
      add.labels = add.labels, add.class = add.class,
      drop.list.single.group = drop.list.single.group)
  } else if (what == "h1.information.observed") {
    lav_object_inspect_h1_information(object, information = "observed",
      h1.information = "default", inverted = FALSE,
      add.labels = add.labels, add.class = add.class,
      drop.list.single.group = drop.list.single.group)
  } else if (what == "h1.information.first.order" ||
    what == "h1.information.firstorder"  ||
    what == "h1.first.order") {
    lav_object_inspect_h1_information(object,
      information = "first.order", h1.information = "default",
      inverted = FALSE, add.labels = add.labels, add.class = add.class,
      drop.list.single.group = drop.list.single.group)

  } else if (what == "vcov") {
    lav_object_inspect_vcov(object,
      standardized = FALSE,
      add.labels = add.labels, add.class = add.class)
  } else if (what == "vcov.std.all" || what == "vcov.standardized" ||
    what == "vcov.std") {
    lav_object_inspect_vcov(object,
      standardized = TRUE, type = "std.all",
      add.labels = add.labels, add.class = add.class)
  } else if (what == "vcov.std.lv") {
    lav_object_inspect_vcov(object,
      standardized = TRUE, type = "std.lv",
      add.labels = add.labels, add.class = add.class)
  } else if (what == "vcov.std.nox") {
    lav_object_inspect_vcov(object,
      standardized = TRUE, type = "std.nox",
      add.labels = add.labels, add.class = add.class)

  } else if (what == "vcov.def") {
    lav_object_inspect_vcov_def(object, joint = FALSE,
      standardized = FALSE,
      add.labels = add.labels, add.class = add.class)
  } else if (what == "vcov.def.std.all" || what == "vcov.def.standardized" ||
    what == "vcov.def.std") {
    lav_object_inspect_vcov_def(object, joint = FALSE,
      standardized = TRUE, type = "std.all",
      add.labels = add.labels, add.class = add.class)
  } else if (what == "vcov.def.std.lv") {
    lav_object_inspect_vcov_def(object, joint = FALSE,
      standardized = TRUE, type = "std.lv",
      add.labels = add.labels, add.class = add.class)
  } else if (what == "vcov.def.std.nox") {
    lav_object_inspect_vcov_def(object, joint = FALSE,
      standardized = TRUE, type = "std.nox",
      add.labels = add.labels, add.class = add.class)

  } else if (what == "vcov.def.joint") {
    lav_object_inspect_vcov_def(object, joint = TRUE,
      standardized = FALSE,
      add.labels = add.labels, add.class = add.class)
  } else if (what == "vcov.def.joint.std.all" ||
    what == "vcov.def.joint.standardized" ||
    what == "vcov.def.joint.std") {
    lav_object_inspect_vcov_def(object, joint = TRUE,
      standardized = TRUE, type = "std.all",
      add.labels = add.labels, add.class = add.class)
  } else if (what == "vcov.def.joint.std.lv") {
    lav_object_inspect_vcov_def(object, joint = TRUE,
      standardized = TRUE, type = "std.lv",
      add.labels = add.labels, add.class = add.class)
  } else if (what == "vcov.def.joint.std.nox") {
    lav_object_inspect_vcov_def(object, joint = TRUE,
      standardized = TRUE, type = "std.nox",
      add.labels = add.labels, add.class = add.class)

  } else if (what == "ugamma" || what == "ug" || what == "u.gamma") {
    lav_object_inspect_UGamma(object,
      add.labels = add.labels, add.class = add.class)
  } else if (what == "ufromugamma" || what == "u") {
    lav_object_inspect_UfromUGamma(object,
      add.labels = add.labels, add.class = add.class,
      drop.list.single.group = drop.list.single.group)

    ### jacobians ####
  } else if (what == "delta") {
    lav_object_inspect_delta(object,
      add.labels = add.labels, add.class = add.class,
      drop.list.single.group = drop.list.single.group)
  } else if (what == "delta.rownames") {
    lav_object_inspect_delta_rownames(object,
      drop.list.single.group = drop.list.single.group)

    ### casewise loglikehoods ###
  } else if (what == "loglik.casewise") {
    lav_object_inspect_loglik_casewise(object, log. = TRUE,
      add.labels = add.labels, add.class = add.class,
      drop.list.single.group = drop.list.single.group)
  } else if (what == "lik.casewise") {
    lav_object_inspect_loglik_casewise(object, log. = FALSE,
      add.labels = add.labels, add.class = add.class,
      drop.list.single.group = drop.list.single.group)

    # multilevel #
  } else if (what == "icc") {
    lav_object_inspect_icc(object,
      add.labels = add.labels, add.class = add.class,
      drop.list.single.group = drop.list.single.group)
  } else if (what == "ranef") {
    lav_object_inspect_ranef(object,
      add.labels = add.labels, add.class = add.class,
      drop.list.single.group = drop.list.single.group)

    # post-checking
  } else if (what == "post.check" || what == "post") {
    lav_object_post_check(object)

    # options
  } else if (what == "options" || what == "lavoptions") {
    object@Options

    # version
  } else if (what == "version") {
    object@version

    # call
  } else if (what == "call") {
    as.list(object@call)

    # timing
  } else if (what == "timing") {
    object@timing

    # optim
  } else if (what == "optim") {
    object@optim

    # test
  } else if (what == "test") {
    object@test

    # zero cell tables
  } else if (what == "zero.cell.tables") {
    lav_object_inspect_zero_cell_tables(object,
      add.labels = add.labels, add.class = add.class,
      drop.list.single.group = drop.list.single.group)

    #### not found ####
  } else {
    lav_msg_stop(gettextf(
      "%1$s argument unknown: %2$s",
      "what", lav_msg_view(what)
    ))
  }

}


# helper functions (mostly to deal with older 'object' that may have
# been saved somewhere)
lav_object_inspect_est <- function(object, unrotated = FALSE) {

  if (inherits(object, "lavaan")) {
    # from 0.5-19, they are in the partable
    if (!is.null(object@ParTable$est)) {
      if (unrotated) {
        return.value <- object@ParTable$est.unrotated
      } else {
        return.value <- object@ParTable$est # if this changes, tag @TDJorgensen in commit message
      }
    } else if (.hasSlot(object, "Fit")) {
      # in < 0.5-19, we should look in @Fit@est
      return.value <- object@Fit@est
    } else {
      partable <- parTable(object)
      return.value <- rep(as.numeric(NA), length(partable$lhs))
    }
  } else {
    # try generic coef()
    return.value <- coef(object, type = "user")
    if (is.matrix(return.value)) {
      # lavaanList?
      return.value <- rowMeans(return.value)
    }
  }

  return.value
}

lav_object_inspect_se <- function(object) {
  # from 0.5-19, they are in the partable
  if (!is.null(object@ParTable$se)) {
    return.value <- object@ParTable$se
  } else if (.hasSlot(object, "Fit")) {
    # in < 0.5-19, we should look in @Fit@se
    return.value <- object@Fit@se
  } else {
    partable <- parTable(object)
    return.value <- rep(as.numeric(NA), length(partable$lhs))
  }

  return.value
}

lav_object_inspect_std_se <- function(object) {

  if (!is.null(object@ParTable$se.std)) {
    return.value <- object@ParTable$se.std
  } else {
    tmp.std <- standardizedSolution(object)
    return.value <- tmp.std$se
  }

  return.value
}

lav_object_inspect_start <- function(object) {
  # from 0.5-19, they are in the partable
  if (!is.null(object@ParTable$start)) {
    return.value <- object@ParTable$start
  } else {
    # in < 0.5-19, we should look in @Fit@start
    return.value <- object@Fit@start
  }

  return.value
}

lav_object_inspect_boot <- function(object, add.labels = FALSE,
                                    add.class = FALSE) {

  if (object@Options$se   != "bootstrap" &&
    !any(c("bootstrap", "bollen.stine") %in% object@Options$test)) {
    lav_msg_stop(gettext("bootstrap was not used."))
  }

  # from 0.5-19. they are in a separate slot
  tmp <- try(slot(object, "boot"), silent = TRUE)
  if (inherits(tmp, "try-error")) {
    # older version of object?
    est <- lav_object_inspect_est(object)
    tmp.boot <- attr(est, "tmp.boot.COEF")
  } else {
    # 0.5-19 way
    tmp.boot <- object@boot$coef
  }

  # add coef names
  if (add.labels) {
    colnames(tmp.boot) <- names(coef(object))
  }

  # add class
  if (add.class) {
    class(tmp.boot) <- c("lavaan.matrix", "matrix")
  }

  tmp.boot
}


lav_object_inspect_modelmatrices <- function(object, what = "free", # nolint
    type = "free", add.labels = FALSE, add.class = FALSE,
    list.by.group = FALSE,
    drop.list.single.group = FALSE) {

  glist <- object@Model@GLIST

  if (what == "dx.free") {
    tmp.dx <- lav_model_gradient(
      lavmodel       = object@Model,
      GLIST          = NULL,
      lavsamplestats = object@SampleStats,
      lavdata        = object@Data,
      lavcache       = object@Cache,
      type           = "free",
      verbose        = FALSE,
      group.weight   = TRUE,
      ceq.simple     = TRUE,
      Delta          = NULL)
  } else if (what == "dx.all") {
    glist <- lav_model_gradient(lavmodel   = object@Model,
      GLIST          = NULL,
      lavsamplestats = object@SampleStats,
      lavdata        = object@Data,
      lavcache       = object@Cache,
      type           = "allofthem",
      verbose        = FALSE,
      group.weight   = TRUE,
      ceq.simple     = FALSE,
      Delta          = NULL)
    names(glist) <- names(object@Model@GLIST)
  } else if (what == "std.all") {
    tmp.std <- lav_standardize_all(object)
  } else if (what == "std.lv") {
    tmp.std <- lav_standardize_lv(object)
  } else if (what == "std.nox") {
    tmp.std <- lav_standardize_all_nox(object)
  } else if (what == "se") {
    tmp.se <- lav_object_inspect_se(object)
  } else if (what == "std.se") {
    tmp.se <- lav_object_inspect_std_se(object)
  } else if (what == "start") {
    tmp.start <- lav_object_inspect_start(object)
  } else if (what == "est") {
    tmp.est <- lav_object_inspect_est(object)
  } else if (what == "est.unrotated") {
    if (!is.null(object@Options$rotation) &&
      object@Options$rotation == "none") {
      tmp.est <- lav_object_inspect_est(object, unrotated = FALSE)
    } else {
      tmp.est <- lav_object_inspect_est(object, unrotated = TRUE)
    }
  }

  for (mm in seq_along(glist)) {

    if (add.labels) {
      dimnames(glist[[mm]]) <- object@Model@dimNames[[mm]]
    }

    if (what == "free") {
      # fill in free parameter counts
      if (type == "free") {
        m.el.idx <- object@Model@m.free.idx[[mm]]
        x.el.idx <- object@Model@x.free.idx[[mm]]
        # } else if(type == "unco") {
        #    m.el.idx <- object@Model@m.unco.idx[[mm]]
        #    x.el.idx <- object@Model@x.unco.idx[[mm]]
      } else if (type == "partable") {
        m.el.idx <- object@Model@m.user.idx[[mm]]
        x.el.idx <- object@Model@x.user.idx[[mm]]
      } else {
        lav_msg_stop(gettextf(
          "%1$s argument unknown: %2$s",
          "type", lav_msg_view(type)
        ))
      }
      # erase everything
      glist[[mm]][, ] <- 0.0
      glist[[mm]][m.el.idx] <- x.el.idx
    } else if (what == "se" || what == "std.se") {
      # fill in standard errors
      m.user.idx <- object@Model@m.user.idx[[mm]]
      x.user.idx <- object@Model@x.user.idx[[mm]]
      # erase everything
      glist[[mm]][, ] <- 0.0
      glist[[mm]][m.user.idx] <- tmp.se[x.user.idx]
    } else if (what == "start") {
      # fill in starting values
      m.user.idx <- object@Model@m.user.idx[[mm]]
      x.user.idx <- object@Model@x.user.idx[[mm]]
      glist[[mm]][m.user.idx] <- tmp.start[x.user.idx]
    } else if (what %in% c("est", "est.unrotated")) {
      # fill in estimated parameter values
      m.user.idx <- object@Model@m.user.idx[[mm]]
      x.user.idx <- object@Model@x.user.idx[[mm]]
      glist[[mm]][m.user.idx] <- tmp.est[x.user.idx]
    } else if (what == "dx.free") {
      # fill in derivatives free parameters
      m.el.idx <- object@Model@m.free.idx[[mm]]
      x.el.idx <- object@Model@x.free.idx[[mm]]
      # erase everything
      glist[[mm]][, ] <- 0.0
      glist[[mm]][m.el.idx] <- tmp.dx[x.el.idx]
    } else if (what %in% c("std.all", "std.lv", "std.nox")) {
      m.user.idx <- object@Model@m.user.idx[[mm]]
      x.user.idx <- object@Model@x.user.idx[[mm]]
      glist[[mm]][m.user.idx] <- tmp.std[x.user.idx]
    }

    # class
    if (add.class) {
      if (object@Model@isSymmetric[mm]) {
        class(glist[[mm]]) <- c("lavaan.matrix.symmetric", "matrix")
      } else {
        class(glist[[mm]]) <- c("lavaan.matrix", "matrix")
      }
    }
  }

  # try to reflect `equality constraints'
  con.flag <- FALSE
  if (what == "free" && object@Model@eq.constraints) {
    # extract constraints from parameter table
    partable <- parTable(object)
    tmp.con <-  partable[partable$op %in% c("==", "<", ">"),
                         c("lhs", "op", "rhs")]
    rownames(tmp.con) <- NULL

    # replace 'labels' by parameter numbers
    tmp.id <- lav_partable_constraints_label_id(partable)
    tmp.label <- names(tmp.id)
    for (con in seq_len(nrow(tmp.con))) {
      # lhs
      lhs.labels <- all.vars(as.formula(paste("~", tmp.con[con, "lhs"])))

      if (length(lhs.labels) > 0L) {
        # par id
        lhs.freeid <- tmp.id[match(lhs.labels, tmp.label)]

        # substitute
        tmp <- tmp.con[con, "lhs"]
        for (pat in seq_along(lhs.labels)) {
          tmp <- sub(lhs.labels[pat], lhs.freeid[pat], tmp)
        }
        tmp.con[con, "lhs"] <- tmp
      }

      # rhs
      rhs.labels <- all.vars(as.formula(paste("~", tmp.con[con, "rhs"])))

      if (length(rhs.labels) > 0L) {
        # par id
        rhs.freeid <- tmp.id[match(rhs.labels, tmp.label)]
        # substitute
        tmp <- tmp.con[con, "rhs"]
        for (pat in seq_along(rhs.labels)) {
          tmp <- sub(rhs.labels[pat], rhs.freeid[pat], tmp)
        }
        tmp.con[con, "rhs"] <- tmp
      }
    } # con

    # add this info at the top
    # glist <- c(constraints = list(tmp.con), glist)
    # no, not a good idea, it does not work with list.by.group

    # add it as a 'header' attribute?
    attr(tmp.con, "header") <- "Note: model contains equality constraints:"
    con.flag <- TRUE
  }

  # should we group them per block?
  if (list.by.group) {
    lavmodel       <- object@Model
    nmat           <- lavmodel@nmat

    return.value <- vector("list", length = lavmodel@nblocks)
    for (b in seq_len(lavmodel@nblocks)) {
      # which mm belong to this block?
      mm.in.group <- 1:nmat[b] + cumsum(c(0, nmat))[b]

      return.value[[b]] <- glist[mm.in.group]
    }

    if (lavmodel@nblocks == 1L && drop.list.single.group) {
      return.value <- return.value[[1]]
    } else if (lavmodel@nblocks > 1L) {
      names(return.value) <- object@Data@block.label
    }
  } else {
    return.value <- glist
  }

  # header
  if (con.flag) {
    attr(return.value, "header") <- tmp.con
  }

  # lavaan.list
  if (add.class) {
    class(return.value) <- c("lavaan.list", "list")
  }

  return.value
}




# - fixme, should we export this function?
# - since 0.5-21, conditional.x = TRUE returns residual sample statistics
#    for ML, we have both joint and residual cov/var/...; but for
#    categorical = TRUE, we only have residual cov/var...; so, we
#    only return residual in both cases, whenever residual
# - since 0.6-3, we always extract the values from the @h1 slot (if present)
#   if meanstructure = FALSE, do NOT include $mean elements any longer
lav_object_inspect_sampstat <- function(object, h1 = TRUE,        # nolint
    std = FALSE,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {
  # old (<0.6) object?
  if (!.hasSlot(object@Data, "block.label")) {
    object@Data@block.label <- object@Data@group.label
  }

  nblocks <- object@Model@nblocks
  ov.names <- object@pta$vnames$ov
  ov.names.res <- object@pta$vnames$ov.nox
  ov.names.x   <- object@pta$vnames$ov.x

  # slots
  lavsamplestats <- object@SampleStats
  lavmodel <- object@Model

  # if nlevels, override h1 to be TRUE, and set conditional.x = FALSE
  if (object@Data@nlevels > 1L) {
    h1 <- TRUE
    conditional.x <- FALSE # for now (0.6-12)
  } else {
    conditional.x <- lavmodel@conditional.x
  }

  # check if we have a non-empty @h1 slot
  if (!.hasSlot(object, "h1")) {
    h1 <- FALSE
  } else if (length(object@h1) == 0L) {
    h1 <- FALSE
  } else {
    h1.implied <- object@h1$implied
  }

  # if h1 = FALSE and nlevels > 1L, nothing can show...
  if (!h1 && object@Data@nlevels > 1L) {
    lav_msg_stop(gettext(
      "sample statistics not available; refit with option h1 = TRUE"))
  }

  return.value <- vector("list", length = nblocks)
  for (b in seq_len(nblocks)) {

    if (!conditional.x) {
      # covariance matrix
      if (h1) {
        return.value[[b]]$cov  <- h1.implied$cov[[b]]
      } else {
        return.value[[b]]$cov  <- lavsamplestats@cov[[b]]
      }
      if (std) {
        diag.orig <- diag(return.value[[b]]$cov)
        return.value[[b]]$cov <- cov2cor(return.value[[b]]$cov)
      }
      if (add.labels && !is.null(return.value[[b]]$cov)) {
        rownames(return.value[[b]]$cov) <- colnames(return.value[[b]]$cov) <-
          ov.names[[b]]
      }
      if (add.class) {
        class(return.value[[b]]$cov) <- c("lavaan.matrix.symmetric", "matrix")
      }

      # mean vector
      if (lavmodel@meanstructure) {
        if (h1) {
          return.value[[b]]$mean <- as.numeric(h1.implied$mean[[b]])
        } else {
          return.value[[b]]$mean <- as.numeric(lavsamplestats@mean[[b]])
        }
        if (std) {
          diag.orig[diag.orig < .Machine$double.eps] <- NA
          return.value[[b]]$mean <- return.value[[b]]$mean / sqrt(diag.orig)
        }
        if (add.labels) {
          names(return.value[[b]]$mean) <- ov.names[[b]]
        }
        if (add.class) {
          class(return.value[[b]]$mean) <- c("lavaan.vector", "numeric")
        }
      }

      # thresholds
      if (lavmodel@categorical) {
        if (h1) {
          return.value[[b]]$th <- as.numeric(h1.implied$th[[b]])
        } else {
          return.value[[b]]$th <- as.numeric(lavsamplestats@th[[b]])
        }
        if (length(lavmodel@num.idx[[b]]) > 0L) {
          num.idx <- which(lavmodel@th.idx[[b]] == 0)
          return.value[[b]]$th <- return.value[[b]]$th[-num.idx]
        }
        # FIXME: what to do if std = TRUE (depends on delta/theta)
        if (add.labels) {
          names(return.value[[b]]$th) <- object@pta$vnames$th[[b]]
        }
        if (add.class) {
          class(return.value[[b]]$th) <- c("lavaan.vector", "numeric")
        }
      }
      # !conditional.x
    } else {
      # if conditional.x = TRUE

      # residual covariance matrix
      if (h1) {
        return.value[[b]]$res.cov  <- h1.implied$res.cov[[b]]
      } else {
        return.value[[b]]$res.cov  <- lavsamplestats@res.cov[[b]]
      }
      if (std) {
        diag.orig <- diag(return.value[[b]]$res.cov)
        return.value[[b]]$res.cov <- cov2cor(return.value[[b]]$res.cov)
      }
      if (add.labels) {
        rownames(return.value[[b]]$res.cov) <-
          colnames(return.value[[b]]$res.cov) <-
          ov.names.res[[b]]
      }
      if (add.class) {
        class(return.value[[b]]$res.cov) <-
          c("lavaan.matrix.symmetric", "matrix")
      }

      # intercepts
      if (lavmodel@meanstructure) {
        if (h1) {
          return.value[[b]]$res.int <- as.numeric(h1.implied$res.int[[b]])
        } else {
          return.value[[b]]$res.int <- as.numeric(lavsamplestats@res.int[[b]])
        }
        if (std) {
          diag.orig[diag.orig < .Machine$double.eps] <- NA
          return.value[[b]]$res.int <- return.value[[b]]$res.int /
            sqrt(diag.orig)
        }
        if (add.labels) {
          names(return.value[[b]]$res.int) <- ov.names.res[[b]]
        }
        if (add.class) {
          class(return.value[[b]]$res.int) <- c("lavaan.vector", "numeric")
        }
      }

      # thresholds
      if (lavmodel@categorical) {
        if (h1) {
          return.value[[b]]$res.th <- as.numeric(h1.implied$res.th[[b]])
        } else {
          return.value[[b]]$res.th <- as.numeric(lavsamplestats@res.th[[b]])
        }
        if (length(lavmodel@num.idx[[b]]) > 0L) {
          num.idx <- which(lavmodel@th.idx[[b]] == 0)
          return.value[[b]]$res.th <- return.value[[b]]$res.th[-num.idx]
        }
        # FIXME: if std: what to do?
        if (add.labels) {
          names(return.value[[b]]$res.th) <- object@pta$vnames$th[[b]]
        }
        if (add.class) {
          class(return.value[[b]]$res.th) <- c("lavaan.vector", "numeric")
        }
      }

      # slopes
      if (lavmodel@nexo[b] > 0L) {
        if (h1) {
          return.value[[b]]$res.slopes  <- h1.implied$res.slopes[[b]]
        } else {
          return.value[[b]]$res.slopes  <- lavsamplestats@res.slopes[[b]]
        }
        # FIXME: if std: what to do? (here: b.z = b * s.x /s.y)
        if (std) {
          tmp.y <- matrix(sqrt(diag.orig),
            nrow(return.value[[b]]$res.slopes),
            ncol(return.value[[b]]$res.slopes))
          tmp.x <- matrix(sqrt(diag(lavsamplestats@cov.x[[b]])),
            nrow(return.value[[b]]$res.slopes),
            ncol(return.value[[b]]$res.slopes), byrow = TRUE)
          return.value[[b]]$res.slopes <- return.value[[b]]$res.slopes /
                                          tmp.y * tmp.x
        }
        if (add.labels) {
          rownames(return.value[[b]]$res.slopes) <- ov.names.res[[b]]
          colnames(return.value[[b]]$res.slopes) <- ov.names.x[[b]]
        }
        if (add.class) {
          class(return.value[[b]]$res.slopes) <- c("lavaan.matrix", "matrix")
        }
      }

      # cov.x
      if (lavmodel@nexo[b] > 0L) {
        return.value[[b]]$cov.x  <- lavsamplestats@cov.x[[b]]
        if (std) {
          diag.orig <- diag(return.value[[b]]$cov.x)
          return.value[[b]]$cov.x <- cov2cor(return.value[[b]]$cov.x)
        }
        if (add.labels) {
          rownames(return.value[[b]]$cov.x) <- ov.names.x[[b]]
          colnames(return.value[[b]]$cov.x) <- ov.names.x[[b]]
        }
        if (add.class) {
          class(return.value[[b]]$cov.x) <-
            c("lavaan.matrix.symmetric", "matrix")
        }
      }

      # mean.x
      if (lavmodel@nexo[b] > 0L) {
        return.value[[b]]$mean.x <- as.numeric(object@SampleStats@mean.x[[b]])
        if (std) {
          diag.orig[diag.orig < .Machine$double.eps] <- NA
          return.value[[b]]$mean.x <- return.value[[b]]$mean.x / sqrt(diag.orig)
        }
        if (add.labels) {
          names(return.value[[b]]$mean.x) <- ov.names.x[[b]]
        }
        if (add.class) {
          class(return.value[[b]]$mean.x) <- c("lavaan.vector", "numeric")
        }
      }

    } # conditional.x

    # stochastic weights
    if (lavmodel@group.w.free) {
      # to be consistent with the 'implied' values,
      # transform so group.w is the 'log(group.freq)'
      return.value[[b]]$group.w <-
        log(lavsamplestats@group.w[[b]] * lavsamplestats@ntotal)
      if (add.labels) {
        names(return.value[[b]]$group.w) <- "w"
      }
      if (add.class) {
        class(return.value[[b]]$group.w) <- c("lavaan.vector", "numeric")
      }
    }
  }

  if (nblocks == 1L && drop.list.single.group) {
    return.value <- return.value[[1]]
  } else if (nblocks > 1L) {
    names(return.value) <- object@Data@block.label
  }

  return.value
}


lav_object_inspect_data <- function(object, add.labels = FALSE,
                                    drop.list.single.group = FALSE) {

  n.g <- object@Data@ngroups
  if (object@Model@conditional.x) {
    return.value <- vector("list", length = n.g)
    for (g in 1:n.g) {
      return.value[[g]] <- cbind(object@Data@X[[g]],
        object@Data@eXo[[g]])
    }
  } else {
    return.value <- object@Data@X
  }

  if (add.labels) {
    for (g in 1:n.g) {
      if (object@Model@conditional.x) {
        colnames(return.value[[g]]) <- c(object@Data@ov.names[[g]],
          object@Data@ov.names.x[[g]])
      } else {
        colnames(return.value[[g]]) <- object@Data@ov.names[[g]]
      }
    }
  }

  if (n.g == 1L && drop.list.single.group) {
    return.value <- return.value[[1]]
  } else {
    if (length(object@Data@group.label) > 0L) {
      names(return.value) <- unlist(object@Data@group.label)
    }
  }

  return.value
}

lav_object_inspect_case_idx <- function(object,
                                        drop.list.single.group = FALSE) {
  n.g <- object@Data@ngroups

  return.value <- object@Data@case.idx

  if (n.g == 1L && drop.list.single.group) {
    return.value <- return.value[[1]]
  } else {
    if (length(object@Data@group.label) > 0L) {
      names(return.value) <- unlist(object@Data@group.label)
    }
  }

  return.value
}

# lav_object_inspect_case_idx <- function(object, level = 1L,
#                                        drop.list.single.group = FALSE) {
#    #FIXME: if lavaan ever allows 3-level or cross-classifed models,
#    #      "level=" should be a character indicating the clustering variable
#
#    n.g <- object@Data@ngroups
#    nlevels <- object@Data@nlevels
#    if (nlevels == 1L) level <- 1L
#                         # if what="cluster.idx" for single-level model
#
#    if (level == 2L) {
#        # level-2 (cluster) IDs
#        return.value <- lapply(object@Data@Lp, function(gg)
#                           gg$cluster.id[[2]][ gg$cluster.idx[[2]] ])
#        #FIXME: update if lavaan ever accepts 3-level or
#                                           cross-classified models
#
#    } else return.value <- object@Data@case.idx # level-1 (casewise) IDs
#
#    if(n.g == 1L && drop.list.single.group) {
#        return.value <- return.value[[1]]
#    } else {
#        if(length(object@Data@group.label) > 0L) {
#            names(return.value) <- unlist(object@Data@group.label)
#        }
#    }
#    return.value
# }
#

# cluster info
lav_object_inspect_cluster_info <- function(                # nolint
    object,
    what = "cluster.size",
    level = 2L,
    drop.list.single.group = FALSE) {

  n.g <- object@Data@ngroups
  nlevels <- object@Data@nlevels

  # just in case we have no clusters
  if (nlevels == 1L) {
    if (what %in% c("nclusters", "ncluster.size", "cluster.id")) {
      return.value <- as.list(rep(1L, n.g))
    } else if (what %in% c("cluster.size", "cluster.sizes")) {
      return.value <- object@Data@nobs
    } else if (what %in% c("cluster.idx", "cluster.label")) {
      # everybody belongs to cluster 1
      return.value <- lapply(seq_len(n.g),
        function(gg) rep(1L, object@Data@nobs[[gg]]))
    }
  }

  # if we do have clusters
  if (nlevels > 1L) {
    return.value <- vector("list", length = n.g)

    for (g in seq_len(n.g)) {
      tmp.lp <- object@Data@Lp[[g]]
      if (what == "nclusters") {
        return.value[[g]] <- tmp.lp$nclusters[[level]]
      } else if (what == "ncluster.size") {
        return.value[[g]] <- tmp.lp$ncluster.size[[level]]
      } else if (what == "cluster.size") {
        return.value[[g]] <- tmp.lp$cluster.size[[level]]
      } else if (what == "cluster.id") {
        return.value[[g]] <- tmp.lp$cluster.id[[level]]
      } else if (what == "cluster.idx") {
        return.value[[g]] <- tmp.lp$cluster.idx[[level]]
      } else if (what == "cluster.label") {
        return.value[[g]] <-
          tmp.lp$cluster.id[[level]][tmp.lp$cluster.idx[[level]]]
      } else if (what == "cluster.sizes") {
        return.value[[g]] <- tmp.lp$cluster.sizes[[level]]
      } else if (what == "average.cluster.size") {
        nn.g <- object@Data@nobs[[g]]
        cluster.size <- tmp.lp$cluster.size[[level]]
        nclusters <- tmp.lp$nclusters[[level]]
        return.value[[g]] <- (nn.g^2 - sum(cluster.size^2)) /
          (nn.g * (nclusters - 1L))
      }
    } # g
  } # nlevels > 1L

  if (n.g == 1L && drop.list.single.group) {
    return.value <- return.value[[1]]
  } else {
    if (length(object@Data@group.label) > 0L) {
      names(return.value) <- unlist(object@Data@group.label)
    }
  }

  return.value
}


# count the number of clusters, or obtain tmp.n within each cluster
# lav_object_inspect_ncluster <- function(object, sizes = FALSE, #level = 2L,
#                                        drop.list.single.group = FALSE) {
#    n.g <- object@Data@ngroups
#    nlevels <- object@Data@nlevels
#
#    if (nlevels == 1L) {
#      # single-level model, return sample size(s) or count 1 cluster per group
#      return.value <- if (sizes) unlist(object@Data@nobs) else rep(1L, n.g)
#
#    } else if (sizes) {
#      # for each group, a vector of cluster sizes
#      return.value <- lapply(object@Data@Lp, function(gg) gg$cluster.size[[2]])
#      #FIXME: update if lavaan ever accepts 3-level or cross-classified models
#
#      if (n.g == 1L && drop.list.single.group)
#                       return.value <- return.value[[1]]
#
#    } else {
#      # number of clusters in each group
#      return.value <- sapply(object@Data@Lp, function(gg) gg$nclusters[[2]])
#      #FIXME: update if lavaan ever accepts 3-level or cross-classified models
#    }
#
#    # assign group names, if applicable
#    if (n.g > 1L) names(return.value) <- unlist(object@Data@group.label)
#    return.value
# }

lav_object_inspect_rsquare <- function(object, est.std.all = NULL,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {

  nblocks <- object@Model@nblocks
  return.value <- vector("list", length = nblocks)

  if (is.null(est.std.all)) {
    est.std.all <- lav_standardize_all(object)
  }

  partable <- object@ParTable
  partable$rsquare <- 1.0 - est.std.all
  # no values > 1.0
  partable$rsquare[partable$rsquare > 1.0] <- as.numeric(NA)

  for (b in seq_len(nblocks)) {
    ind.names <- partable$rhs[which(partable$op == "=~" &
      partable$block == b)]
    eqs.y.names <- partable$lhs[which(partable$op == "~"  &
      partable$block == b)]
    y.names <- unique(c(ind.names, eqs.y.names))

    idx <- which(partable$op == "~~" & partable$lhs %in% y.names &
      partable$rhs == partable$lhs & partable$block == b)
    tmp <- partable$rsquare[idx]

    if (add.labels && length(tmp) > 0L) {
      names(tmp) <- partable$lhs[idx]
    }
    if (add.class) {
      class(tmp) <- c("lavaan.vector", "numeric")
    }

    return.value[[b]] <- tmp
  }

  if (nblocks == 1L && drop.list.single.group) {
    return.value <- return.value[[1]]
  } else if (nblocks > 1L) {
    names(return.value) <- object@Data@block.label
  }

  return.value
}

# model implied sample stats
lav_object_inspect_implied <- function(object,                    # nolint
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {
  # old (<0.6) object?
  if (!.hasSlot(object@Data, "block.label")) {
    object@Data@block.label <- object@Data@group.label
  }


  nblocks <- object@Model@nblocks
  ov.names <- object@pta$vnames$ov
  ov.names.res <- object@pta$vnames$ov.nox
  ov.names.x   <- object@pta$vnames$ov.x

  # slots
  lavimplied <- object@implied
  lavmodel   <- object@Model

  # if nlevels, always set conditional.x = FALSE
  if (object@Data@nlevels > 1L) {
    lavimplied <- lav_model_implied_cond2uncond(lavimplied)
    conditional.x <- FALSE # for now (0.6-12)
  } else {
    conditional.x <- lavmodel@conditional.x
  }

  return.value <- vector("list", length = nblocks)
  for (b in seq_len(nblocks)) {

    if (!conditional.x) {
      # covariance matrix
      return.value[[b]]$cov  <- lavimplied$cov[[b]]
      if (add.labels && !is.null(return.value[[b]]$cov)) {
        rownames(return.value[[b]]$cov) <- colnames(return.value[[b]]$cov) <-
          ov.names[[b]]
      }
      if (add.class) {
        class(return.value[[b]]$cov) <- c("lavaan.matrix.symmetric", "matrix")
      }

      # mean vector
      if (lavmodel@meanstructure) {
        return.value[[b]]$mean <- as.numeric(lavimplied$mean[[b]])
        if (add.labels) {
          names(return.value[[b]]$mean) <- ov.names[[b]]
        }
        if (add.class) {
          class(return.value[[b]]$mean) <- c("lavaan.vector", "numeric")
        }
      }

      # thresholds
      if (lavmodel@categorical) {
        return.value[[b]]$th <- as.numeric(lavimplied$th[[b]])
        if (length(object@Model@num.idx[[b]]) > 0L) {
          num.idx <- which(object@Model@th.idx[[b]] == 0)
          return.value[[b]]$th <- return.value[[b]]$th[-num.idx]
        }
        if (add.labels) {
          names(return.value[[b]]$th) <- object@pta$vnames$th[[b]]
        }
        if (add.class) {
          class(return.value[[b]]$th) <- c("lavaan.vector", "numeric")
        }
      }
      # !conditional.x
    } else {
      # if conditional.x = TRUE
      # residual covariance matrix
      return.value[[b]]$res.cov  <- lavimplied$res.cov[[b]]
      if (add.labels) {
        rownames(return.value[[b]]$res.cov) <-
          colnames(return.value[[b]]$res.cov) <-
          ov.names.res[[b]]
      }
      if (add.class) {
        class(return.value[[b]]$res.cov) <-
          c("lavaan.matrix.symmetric", "matrix")
      }

      # intercepts
      if (lavmodel@meanstructure) {
        return.value[[b]]$res.int <- as.numeric(lavimplied$res.int[[b]])
        if (add.labels) {
          names(return.value[[b]]$res.int) <- ov.names.res[[b]]
        }
        if (add.class) {
          class(return.value[[b]]$res.int) <- c("lavaan.vector", "numeric")
        }
      }

      # thresholds
      if (lavmodel@categorical) {
        return.value[[b]]$res.th <- as.numeric(lavimplied$res.th[[b]])
        if (length(object@Model@num.idx[[b]]) > 0L) {
          num.idx <- which(object@Model@th.idx[[b]] == 0)
          return.value[[b]]$res.th <- return.value[[b]]$res.th[-num.idx]
        }
        if (add.labels) {
          names(return.value[[b]]$res.th) <- object@pta$vnames$th[[b]]
        }
        if (add.class) {
          class(return.value[[b]]$res.th) <- c("lavaan.vector", "numeric")
        }
      }

      # slopes
      if (lavmodel@nexo[b] > 0L) {
        return.value[[b]]$res.slopes  <- lavimplied$res.slopes[[b]]
        if (add.labels) {
          rownames(return.value[[b]]$res.slopes) <- ov.names.res[[b]]
          colnames(return.value[[b]]$res.slopes) <- ov.names.x[[b]]
        }
        if (add.class) {
          class(return.value[[b]]$res.slopes) <- c("lavaan.matrix", "matrix")
        }
      }

      # cov.x
      if (lavmodel@nexo[b] > 0L) {
        return.value[[b]]$cov.x  <- lavimplied$cov.x[[b]]
        if (add.labels) {
          rownames(return.value[[b]]$cov.x) <- ov.names.x[[b]]
          colnames(return.value[[b]]$cov.x) <- ov.names.x[[b]]
        }
        if (add.class) {
          class(return.value[[b]]$cov.x) <-
            c("lavaan.matrix.symmetric", "matrix")
        }
      }

      # mean.x
      if (lavmodel@nexo[b] > 0L) {
        return.value[[b]]$mean.x  <- as.numeric(lavimplied$mean.x[[b]])
        if (add.labels) {
          names(return.value[[b]]$mean.x) <- ov.names.x[[b]]
        }
        if (add.class) {
          class(return.value[[b]]$mean.x) <- c("lavaan.vector", "numeric")
        }
      }


    } # conditional.x

    # stochastic weights
    if (lavmodel@group.w.free) {
      return.value[[b]]$group.w <- lavimplied$group.w[[b]]
      if (add.labels) {
        names(return.value[[b]]$group.w) <- "w" # somewhat redundant
      }
      if (add.class) {
        class(return.value[[b]]$group.w) <- c("lavaan.vector", "numeric")
      }
    }
  }

  if (nblocks == 1L && drop.list.single.group) {
    return.value <- return.value[[1]]
  } else if (nblocks > 1L) {
    names(return.value) <- object@Data@block.label
  }

  return.value
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
  return.value <- computeVETA(lavmodel = object@Model, remove.dummy.lv = TRUE)

  # nblocks
  nblocks <- length(return.value)

  # cor + labels + class
  for (b in seq_len(nblocks)) {

    if (correlation.metric && nrow(return.value[[b]]) > 1L) {
      # note: cov2cor fails if matrix is empty!
      return.value[[b]] <- cov2cor(return.value[[b]])
    }

    if (add.labels) {
      colnames(return.value[[b]]) <- rownames(return.value[[b]]) <-
        object@pta$vnames$lv[[b]]
    }

    if (add.class) {
      class(return.value[[b]]) <- c("lavaan.matrix.symmetric", "matrix")
    }
  }

  if (nblocks == 1L && drop.list.single.group) {
    return.value <- return.value[[1]]
  } else if (nblocks > 1L) {
    names(return.value) <- object@Data@block.label
  }

  return.value
}

lav_object_inspect_mean_lv <- function(object,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {
  # compute lv means
  return.value <- computeEETA(lavmodel       = object@Model,
    lavsamplestats = object@SampleStats,
    remove.dummy.lv = TRUE)

  # nblocks
  nblocks <- length(return.value)

  # ensure numeric
  return.value <- lapply(return.value, as.numeric)

  # labels + class
  for (b in seq_len(nblocks)) {
    if (add.labels && length(return.value[[b]]) > 0L) {
      names(return.value[[b]]) <- object@pta$vnames$lv[[b]]
    }
    if (add.class) {
      class(return.value[[b]]) <- c("lavaan.vector", "numeric")
    }
  }

  if (nblocks == 1L && drop.list.single.group) {
    return.value <- return.value[[1]]
  } else if (nblocks > 1L) {
    names(return.value) <- object@Data@block.label
  }

  return.value
}

lav_object_inspect_cov_all <- function(object, correlation.metric = FALSE,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {
  # compute extended model implied covariance matrix (both ov and lv)
  return.value <- computeCOV(lavmodel = object@Model,
    remove.dummy.lv = TRUE)

  # nblocks
  nblocks <- length(return.value)

  # cor + labels + class
  for (b in seq_len(nblocks)) {

    if (correlation.metric && nrow(return.value[[b]]) > 1L) {
      # note: cov2cor fails if matrix is empty!
      return.value[[b]] <- cov2cor(return.value[[b]])
    }

    if (add.labels) {
      tmp.names <- c(object@pta$vnames$ov.model[[b]],
        object@pta$vnames$lv[[b]])
      colnames(return.value[[b]]) <- rownames(return.value[[b]]) <- tmp.names
    }
    if (add.class) {
      class(return.value[[b]]) <- c("lavaan.matrix.symmetric", "matrix")
    }
  }

  if (nblocks == 1L && drop.list.single.group) {
    return.value <- return.value[[1]]
  } else if (nblocks > 1L) {
    names(return.value) <- object@Data@block.label
  }

  return.value
}


lav_object_inspect_cov_ov <- function(object, correlation.metric = FALSE,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {
  # get model-implied covariance matrix observed
  if (object@Model@conditional.x) {
    return.value <- object@implied$res.cov
  } else {
    return.value <- object@implied$cov
  }

  # nblocks
  nblocks <- length(return.value)

  # cor + labels + class
  for (b in seq_len(nblocks)) {

    if (correlation.metric && nrow(return.value[[b]]) > 1L) {
      # note: cov2cor fails if matrix is empty!
      return.value[[b]] <- cov2cor(return.value[[b]])
    }

    if (add.labels) {
      colnames(return.value[[b]]) <- rownames(return.value[[b]]) <-
        object@pta$vnames$ov.model[[b]]
    }
    if (add.class) {
      class(return.value[[b]]) <- c("lavaan.matrix.symmetric", "matrix")
    }
  }

  if (nblocks == 1L && drop.list.single.group) {
    return.value <- return.value[[1]]
  } else {
    names(return.value) <- object@Data@block.label
  }

  return.value
}

lav_object_inspect_mean_ov <- function(object,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {
  # compute ov means
  if (object@Model@conditional.x) {
    return.value <- object@implied$res.int
  } else {
    return.value <- object@implied$mean
  }

  # nblocks
  nblocks <- length(return.value)

  # make numeric
  return.value <- lapply(return.value, as.numeric)

  # labels + class
  for (b in seq_len(nblocks)) {
    if (add.labels && length(return.value[[b]]) > 0L) {
      names(return.value[[b]]) <- object@pta$vnames$ov.model[[b]]
    }
    if (add.class) {
      class(return.value[[b]]) <- c("lavaan.vector", "numeric")
    }
  }

  if (nblocks == 1L && drop.list.single.group) {
    return.value <- return.value[[1]]
  } else if (nblocks > 1L) {
    names(return.value) <- object@Data@block.label
  }

  return.value
}

lav_object_inspect_th <- function(object,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {
  # thresholds
  if (object@Model@conditional.x) {
    return.value <- object@implied$res.th
  } else {
    return.value <- object@implied$th
  }

  # nblocks
  nblocks <- length(return.value)

  # make numeric
  return.value <- lapply(return.value, as.numeric)

  # labels + class
  for (b in seq_len(nblocks)) {
    if (length(object@Model@num.idx[[b]]) > 0L) {
      num.idx <- which(object@Model@th.idx[[b]] == 0)
      return.value[[b]] <- return.value[[b]][-num.idx]
    }
    if (add.labels && length(return.value[[b]]) > 0L) {
      names(return.value[[b]]) <- object@pta$vnames$th[[b]]
    }
    if (add.class) {
      class(return.value[[b]]) <- c("lavaan.vector", "numeric")
    }
  }

  if (nblocks == 1L && drop.list.single.group) {
    return.value <- return.value[[1]]
  } else if (nblocks > 1L) {
    names(return.value) <- object@Data@block.label
  }

  return.value
}

lav_object_inspect_th_idx <- function(object,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {
  # thresholds idx
  return.value <- object@SampleStats@th.idx

  # nblocks
  nblocks <- length(return.value)

  # labels + class
  for (b in seq_len(nblocks)) {
    if (add.labels && length(return.value[[b]]) > 0L) {
      names(return.value[[b]]) <- object@SampleStats@th.names[[b]]
    }
    if (add.class && !is.null(return.value[[b]])) {
      class(return.value[[b]]) <- c("lavaan.vector", "numeric")
    }
  }

  if (nblocks == 1L && drop.list.single.group) {
    return.value <- return.value[[1]]
  } else if (nblocks > 1L) {
    names(return.value) <- object@Data@block.label
  }

  return.value
}

lav_object_inspect_vy <- function(object,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {
  # 'unconditional' model-implied variances
  #  - same as diag(Sigma.hat) if all Y are continuous)
  #  - 1.0 (or delta^2) if categorical
  #  - if also Gamma, cov.x is used (only if categorical)

  return.value <- computeVY(lavmodel = object@Model, GLIST = NULL,
    diagonal.only = TRUE)

  # nblocks
  nblocks <- length(return.value)

  # labels + class
  for (b in seq_len(nblocks)) {
    if (add.labels && length(return.value[[b]]) > 0L) {
      if (object@Model@categorical) {
        names(return.value[[b]]) <- object@pta$vnames$ov.nox[[b]]
      } else {
        names(return.value[[b]]) <- object@pta$vnames$ov[[b]]
      }
    }
    if (add.class) {
      class(return.value[[b]]) <- c("lavaan.vector", "numeric")
    }
  }

  if (nblocks == 1L && drop.list.single.group) {
    return.value <- return.value[[1]]
  } else if (nblocks > 1L) {
    names(return.value) <- object@Data@block.label
  }

  return.value
}


lav_object_inspect_theta <- function(object, correlation.metric = FALSE,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {
  # get residual covariances
  return.value <- computeTHETA(lavmodel = object@Model)

  # nblocks
  nblocks <- length(return.value)

  # labels + class
  for (b in seq_len(nblocks)) {

    if (correlation.metric && nrow(return.value[[b]]) > 0L) {
      if (all(return.value[[b]] == 0)) {
        return.value[[b]] <- return.value[[b]]
      } else {
        return.value[[b]] <- cov2cor(return.value[[b]])
      }
    }

    if (add.labels && length(return.value[[b]]) > 0L) {
      colnames(return.value[[b]]) <- rownames(return.value[[b]]) <-
        object@pta$vnames$ov.model[[b]]
    }

    if (add.class) {
      class(return.value[[b]]) <- c("lavaan.matrix.symmetric", "matrix")
    }
  }

  if (nblocks == 1L && drop.list.single.group) {
    return.value <- return.value[[1]]
  } else if (nblocks > 1L) {
    names(return.value) <- object@Data@block.label
  }

  return.value
}


lav_object_inspect_missing_coverage <- function(object,          # nolint
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {

  n.g <- object@Data@ngroups
  return.value <- vector("list", n.g)

  for (g in 1:n.g) {
    if (!is.null(object@Data@Mp[[g]])) {
      return.value[[g]] <- object@Data@Mp[[g]]$coverage
    } else {
      nvar <- length(object@Data@ov.names[[g]])
      return.value[[g]] <- matrix(1.0, nvar, nvar)
    }

    if (add.labels && length(return.value[[g]]) > 0L) {
      colnames(return.value[[g]]) <- rownames(return.value[[g]]) <-
        object@pta$vnames$ov.model[[g]]
    }

    if (add.class) {
      class(return.value[[g]]) <- c("lavaan.matrix.symmetric", "matrix")
    }
  }

  if (n.g == 1L && drop.list.single.group) {
    return.value <- return.value[[1]]
  } else {
    if (length(object@Data@group.label) > 0L) {
      names(return.value) <- unlist(object@Data@group.label)
    }
  }

  return.value
}

lav_object_inspect_missing_patterns <- function(object,           # nolint
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {

  n.g <- object@Data@ngroups
  return.value <- vector("list", n.g)

  for (g in 1:n.g) {
    if (!is.null(object@Data@Mp[[g]])) {
      return.value[[g]] <- object@Data@Mp[[g]]$pat
    } else {
      nvar <- length(object@Data@ov.names[[g]])
      return.value[[g]] <- matrix(TRUE, 1L, nvar)
      rownames(return.value[[g]]) <- object@Data@nobs[[g]]
    }

    if (add.labels && length(return.value[[g]]) > 0L) {
      colnames(return.value[[g]]) <- object@pta$vnames$ov.model[[g]]
    }

    if (add.class) {
      class(return.value[[g]]) <- c("lavaan.matrix", "matrix")
    }
  }

  if (n.g == 1L && drop.list.single.group) {
    return.value <- return.value[[1]]
  } else {
    if (length(object@Data@group.label) > 0L) {
      names(return.value) <- unlist(object@Data@group.label)
    }
  }

  return.value
}

lav_object_inspect_empty_idx <- function(object,
                                         drop.list.single.group = FALSE) {

  n.g <- object@Data@ngroups

  # get empty idx
  return.value <- vector("list", n.g)

  for (g in 1:n.g) {
    if (!is.null(object@Data@Mp[[g]])) {
      return.value[[g]] <- object@Data@Mp[[g]]$empty.idx
    } else {
      return.value[[g]] <- integer(0L)
    }
  }

  if (n.g == 1L && drop.list.single.group) {
    return.value <- return.value[[1]]
  } else {
    if (length(object@Data@group.label) > 0L) {
      names(return.value) <- unlist(object@Data@group.label)
    }
  }

  return.value
}


lav_object_inspect_wls_est <- function(object,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {
  # old (<0.6) object?
  if (!.hasSlot(object@Data, "block.label")) {
    object@Data@block.label <- object@Data@group.label
  }


  return.value <- lav_model_wls_est(object@Model)

  if (add.labels) {
    tmp.names <- lav_object_inspect_delta_rownames(object,
      drop.list.single.group = FALSE)
  }

  # nblocks
  nblocks <- length(return.value)

  for (b in seq_len(nblocks)) {
    if (add.labels && length(return.value[[b]]) > 0L &&
        object@Data@nlevels == 1L) {
      names(return.value[[b]]) <- tmp.names[[b]]
    }

    if (add.class) {
      class(return.value[[b]]) <- c("lavaan.vector", "numeric")
    }
  }

  if (nblocks == 1L && drop.list.single.group) {
    return.value <- return.value[[1]]
  } else if (nblocks > 1L) {
    names(return.value) <- object@Data@block.label
  }

  return.value
}

lav_object_inspect_wls_obs <- function(object,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {
  # old (<0.6) object?
  if (!.hasSlot(object@Data, "block.label")) {
    object@Data@block.label <- object@Data@group.label
  }


  return.value <- object@SampleStats@WLS.obs ### FIXME: should be in @h1??

  if (add.labels) {
    tmp.names <- lav_object_inspect_delta_rownames(object,
      drop.list.single.group = FALSE)
  }

  # nblocks
  nblocks <- length(return.value)

  for (b in seq_len(nblocks)) {
    if (add.labels && length(return.value[[b]]) > 0L &&
        object@Data@nlevels == 1L) {
      names(return.value[[b]]) <- tmp.names[[b]]
    }

    if (add.class) {
      class(return.value[[b]]) <- c("lavaan.vector", "numeric")
    }
  }

  if (nblocks == 1L && drop.list.single.group) {
    return.value <- return.value[[1]]
  } else if (nblocks > 1L) {
    names(return.value) <- object@Data@block.label
  }

  return.value
}

lav_object_inspect_wls_v <- function(object,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {
  # old (<0.6) object?
  if (!.hasSlot(object@Data, "block.label")) {
    object@Data@block.label <- object@Data@group.label
  }


  # return.value <- lav_model_wls_v(lavmodel       = object@Model,
  #                       lavsamplestats = object@SampleStats,
  #                       structured     = TRUE,
  #                       lavdata        = object@Data)
  # WLS.V == (traditionally) h1 expected information
  return.value <- lav_model_h1_information_expected(lavobject = object)
  # this affects fit measures gfi, agfi, pgfi


  # nblocks
  nblocks <- length(return.value)

  # if estimator == "DWLS" or "ULS", we only stored the diagonal
  # hence, we create a full matrix here
  if (object@Options$estimator %in% c("DWLS", "ULS")) {
    return.value <- lapply(return.value,
      function(x) {
        nr <- NROW(x)
        diag(x, nrow = nr, ncol = nr)
      })
  }

  if (add.labels) {
    tmp.names <- lav_object_inspect_delta_rownames(object,
      drop.list.single.group = FALSE)
  }

  # label + class
  for (b in seq_len(nblocks)) {
    # labels
    if (add.labels && nrow(return.value[[b]]) > 0L &&
        object@Data@nlevels == 1L) {
      colnames(return.value[[b]]) <-
        rownames(return.value[[b]]) <- tmp.names[[b]]
    }

    # class
    if (add.class) {
      class(return.value[[b]]) <- c("lavaan.matrix", "matrix")
    }
  }

  if (nblocks == 1L && drop.list.single.group) {
    return.value <- return.value[[1]]
  } else if (nblocks > 1L) {
    names(return.value) <- object@Data@block.label
  }

  return.value
}


lav_object_inspect_sampstat_gamma <- function(object,         # nolint
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {

  if (!is.null(object@SampleStats@NACOV[[1]])) {
    return.value <- object@SampleStats@NACOV
  } else {
    return.value <- lav_object_gamma(object)
  }

  if (add.labels) {
    tmp.names <- lav_object_inspect_delta_rownames(object,
      drop.list.single.group = FALSE)
  }

  # nblocks
  nblocks <- length(return.value)

  if (nblocks == 1L && drop.list.single.group) {
    return.value <- return.value[[1]]

    # labels
    if (add.labels) {
      colnames(return.value) <- rownames(return.value) <- tmp.names[[1]]
    }

    # class
    if (add.class) {
      class(return.value) <- c("lavaan.matrix.symmetric", "matrix")
    }
  } else {
    if (object@Data@nlevels == 1L && length(object@Data@group.label) > 0L) {
      names(return.value) <- unlist(object@Data@group.label)

      # labels
      if (add.labels) {
        for (g in seq_len(object@Data@ngroups)) {
          colnames(return.value[[g]]) <-
            rownames(return.value[[g]]) <- tmp.names[[g]]
        }
      }

      # class
      if (add.class) {
        for (g in seq_len(object@Data@ngroups)) {
          class(return.value[[g]]) <- c("lavaan.matrix.symmetric", "matrix")
        }
      }
    } else if (object@Data@nlevels > 1L &&
      length(object@Data@group.label) == 0L) {
      names(return.value) <- object@Data@level.label
    }
  }

  return.value
}


lav_object_inspect_gradient <- function(object,
    add.labels = FALSE, add.class = FALSE, logl = FALSE,
    optim = FALSE) {

  lavmodel       <- object@Model
  lavdata        <- object@Data
  lavsamplestats <- object@SampleStats

  if (optim) {
    logl <- FALSE
  }

  if (lavsamplestats@missing.flag ||
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
  if (logl) {
    if (lavmodel@estimator %in% c("ML")) {
      if (lavdata@nlevels == 1L) {
        # currently, this is just a sign switch
        dx <- -1 * dx
      } else {
        lavpartable <- object@ParTable
        # gradient.log = gradient.obj * (2 * tmp.n) / nclusters

        if (lavdata@ngroups == 1L) {
          tmp.n <- lavdata@Lp[[1]]$nclusters[[1]]
          nclusters <- lavdata@Lp[[1]]$nclusters[[2]]
          dx <- dx * (2 * tmp.n) / nclusters
        } else {
          group.values <- lav_partable_group_values(lavpartable)
          for (g in seq_len(lavdata@ngroups)) {
            tmp.n <- lavdata@Lp[[g]]$nclusters[[1]]
            nclusters <- lavdata@Lp[[g]]$nclusters[[2]]
            g.idx <-
              which((lavpartable$group ==
                group.values[g])[lavpartable$free > 0L])
            dx[g.idx] <- dx[g.idx] * (2 * tmp.n) / nclusters
          }
        }
      }
    } else {
      # FIXME:
      # non-likelihood: what to do? just switch the sign for now.
      # Note: this is used in lavTestScore()
      dx <- -1 * dx
    }
  }

  # optim?
  if (optim) {
    # 1. scale (note: divide, not multiply!)
    if (!is.null(object@optim$parscale)) {
      dx <- dx / object@optim$parscale
    }

    # 2. pack
    if (lavmodel@eq.constraints) {
      dx <- as.numeric(dx %*% lavmodel@eq.constraints.K)
    }
    # only for PML: divide by tmp.n (to speed up convergence)
    if (lavmodel@estimator == "PML") {
      dx <- dx / lavsamplestats@ntotal
    }
  }

  # labels
  if (add.labels) {
    if (optim && lavmodel@eq.constraints) {
      tmp.names.all <- lav_partable_labels(object@ParTable, type = "free")
      tmp.seq <- seq_len(length(tmp.names.all))
      pack.seq <- as.numeric((tmp.seq - lavmodel@eq.constraints.k0) %*%
        +lavmodel@eq.constraints.K)
      ok.idx <- which(pack.seq %in% tmp.seq)
      tmp.names <- rep("(eq.con)", length(pack.seq))
      tmp.names[ok.idx] <- tmp.names.all[pack.seq[ok.idx]]
      names(dx) <- tmp.names
    } else {
      names(dx) <- lav_partable_labels(object@ParTable, type = "free")
    }
  }

  # class
  if (add.class) {
    class(dx) <- c("lavaan.vector", "numeric")
  }

  dx
}

lav_object_inspect_hessian <- function(object,
    add.labels = FALSE, add.class = FALSE) {

  return.value <- lav_model_hessian(lavmodel       = object@Model,
    lavsamplestats = object@SampleStats,
    lavdata        = object@Data,
    lavcache       = object@Cache,
    lavoptions     = object@Options,
    group.weight   = TRUE)

  # labels
  if (add.labels) {
    colnames(return.value) <- rownames(return.value) <-
      lav_partable_labels(object@ParTable, type = "free")
  }

  # class
  if (add.class) {
    class(return.value) <- c("lavaan.matrix.symmetric", "matrix")
  }

  return.value
}

lav_object_inspect_information <- function(object,
    information = "default", augmented = FALSE, inverted = FALSE,
    add.labels = FALSE, add.class = FALSE) {

  if (information != "default") {
    # override option
    object@Options$information <- information
  }

  # backward compatibility
  if (.hasSlot(object, "h1")) {
    lavh1 <- object@h1
  } else {
    lavh1 <- lav_h1_implied_logl(lavdata = object@Data,
      lavsamplestats = object@SampleStats,
      lavoptions = object@Options)
  }

  return.value <- lav_model_information(lavmodel =  object@Model,
    lavsamplestats = object@SampleStats,
    lavdata        = object@Data,
    lavcache       = object@Cache,
    lavimplied     = object@implied,
    lavh1          = lavh1,
    lavoptions     = object@Options,
    extra          = FALSE,
    augmented      = augmented,
    inverted       = inverted)

  # labels
  if (add.labels) {
    tmp.names <- lav_partable_labels(object@ParTable, type = "free")
    if (augmented) {
      n.extra <- nrow(return.value) - length(tmp.names)
      if (n.extra > 0L) {
        tmp.names <- c(tmp.names, paste("aug", 1:n.extra, sep = ""))
      }
    }
    colnames(return.value) <- rownames(return.value) <- tmp.names
  }

  # class
  if (add.class) {
    class(return.value) <- c("lavaan.matrix.symmetric", "matrix")
  }

  return.value
}

lav_object_inspect_h1_information <- function(object,            # nolint
    information = "default", h1.information = "default", inverted = FALSE,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {

  if (information != "default") {
    # override option
    object@Options$information <- information
  }
  if (h1.information != "default") {
    # override option
    object@Options$h1.information <- h1.information
  }

  lavmodel <- object@Model
  lavdata  <- object@Data

  # list!
  return.value <- lav_model_h1_information(lavmodel = lavmodel,
    lavsamplestats = object@SampleStats,
    lavdata        = lavdata,
    lavcache       = object@Cache,
    lavimplied     = object@implied,
    lavh1          = object@h1,
    lavoptions     = object@Options)

  # inverted? (NOT USED)
  # if(inverted) {
  #    return.value <- lapply(return.value, solve) # FIXME: handle errors...
  # }

  if (add.labels) {
    tmp.names <- lav_object_inspect_delta_rownames(object,
      drop.list.single.group = FALSE)
  }

  # labels/class per group
  for (g in seq_len(lavmodel@ngroups)) {
    # labels
    if (add.labels) {
      colnames(return.value[[g]]) <-
        rownames(return.value[[g]]) <- tmp.names[[g]]
    }

    # class
    if (add.class) {
      class(return.value[[g]]) <- c("lavaan.matrix.symmetric", "matrix")
    }
  }

  # drop list?
  if (lavmodel@ngroups == 1L && drop.list.single.group) {
    return.value <- return.value[[1]]
  } else if (!is.null(lavdata)) {
    if (length(lavdata@group.label) > 0L) {
      names(return.value) <- unlist(lavdata@group.label)
    }
  }

  return.value
}

# only to provide a direct function to the old 'getVariability()' function
lav_object_inspect_firstorder <- function(object,
    add.labels = FALSE, add.class = FALSE) {

  tmp.b0 <- lav_model_information_firstorder(lavmodel =  object@Model,
    lavsamplestats = object@SampleStats,
    lavdata        = object@Data,
    lavcache       = object@Cache,
    lavoptions     = object@Options,
    check.pd       = FALSE,
    augmented      = FALSE,
    inverted       = FALSE)
  attr(tmp.b0, "B0.group") <- NULL
  return.value <- tmp.b0

  # labels
  if (add.labels) {
    colnames(return.value) <- rownames(return.value) <-
      lav_partable_labels(object@ParTable, type = "free")
  }

  # class
  if (add.class) {
    class(return.value) <- c("lavaan.matrix.symmetric", "matrix")
  }

  return.value
}

lav_object_inspect_vcov <- function(object, standardized = FALSE,
    type = "std.all", free.only = TRUE,
    add.labels = FALSE, add.class = FALSE, remove.duplicated = FALSE) {

  lavmodel <- object@Model
  lavoptions <- object@Options

  # store partable with pta in object to use cache in called functions
  if (is.null(attr(object, "vnames"))) {
    object@ParTable <- lav_partable_set_cache(object@ParTable, object@pta)
  }
  # rotation?
  # if( .hasSlot(lavmodel, "nefa") && (lavmodel@nefa > 0L) &&
  #    lavoptions$rotation != "none" #&& lavoptions$rotation.se == "delta"
  #  ) {
  #    rotation <- TRUE
  # } else {
  #    rotation <- FALSE
  # }

  if (object@optim$npar == 0) {
    return.value <- matrix(0, 0, 0)
  } else {
    # check if we already have it
    # tmp <- try(slot(object, "vcov"), silent = TRUE)
    # if( !inherits(tmp, "try-error") && !is.null(object@vcov$vcov)
    #   && !(rotation && standardized)) {
    if (.hasSlot(object, "vcov") && !is.null(object@vcov$vcov)) {
      return.value <- object@vcov$vcov # if this changes, tag @TDJorgensen in commit message
    } else {
      # compute it again
      # if(rotation && standardized) {
      #    lavmodel <- lav_model_set_parameters(lavmodel,
      #                    x = object@optim$x)
      #    lavoptions <- object@Options
      #    lavoptions$rotation.se <- "delta"
      # }
      return.value <- lav_model_vcov(lavmodel       = lavmodel,
        lavsamplestats = object@SampleStats,
        lavoptions     = lavoptions,
        lavdata        = object@Data,
        lavcache       = object@Cache,
        lavimplied     = object@implied,
        lavh1          = object@h1
      )
      # if(rotation && !standardized) {
      #    # fixme: compute tmp.vcov.rot manually...
      #    stop("lavaan ERROR: rerun with store.vcov = TRUE")
      # }

      if (is.null(return.value)) {
        return(return.value)
      }
    }
  }

  # strip attributes
  attr(return.value, "E.inv") <- NULL
  attr(return.value, "B0") <- NULL
  attr(return.value, "B0.group") <- NULL
  attr(return.value, "Delta") <- NULL
  attr(return.value, "WLS.V") <- NULL
  attr(return.value, "tmp.boot.COEF") <- NULL
  attr(return.value, "tmp.boot.TEST") <- NULL

  # standardized?
  if (standardized) {
    if (type == "std.lv") {
      tmp.fun <- lav_standardize_lv_x
    } else if (type == "std.all") {
      tmp.fun <- lav_standardize_all_x
    } else if (type == "std.nox") {
      tmp.fun <- lav_standardize_all_nox_x
    }



    # if(rotation) {
    #    if(.hasSlot(object@Model, "ceq.simple.only") &&
    #       object@Model@ceq.simple.only) {
    #        x.vec <- drop(object@optim$x %*% t(object@Model@ceq.simple.K))
    #    } else {
    #        x.vec <- object@optim$x
    #    }
    #    tmp.jac <- numDeriv::jacobian(func = tmp.fun, x = x.vec,
    #               method = "simple",
    #               method.args = list(eps = 1e-03), # default is 1e-04
    #               lavobject = object, rotation = rotation)
    # } else {
    # if(.hasSlot(object@Model, "ceq.simple.only") &&
    #   object@Model@ceq.simple.only) {
    #    x <- lav_model_get_parameters(lavmodel)
    #    x.vec <- drop(x %*% t(object@Model@ceq.simple.K))
    # } else {
    x.vec <- lav_model_get_parameters(lavmodel)
    # }
    tmp.jac <- try(lav_func_jacobian_complex(func = tmp.fun, x = x.vec,
      lavobject = object),
    silent = TRUE)
    if (inherits(tmp.jac, "try-error")) { # eg. pnorm()
      tmp.jac <- lav_func_jacobian_simple(func = tmp.fun, x = x.vec,
        lavobject = object)
    }
    # }

    # tmp.jac contains *all* parameters in the parameter table
    if (free.only) {
      if (.hasSlot(object@Model, "ceq.simple.only") &&
        object@Model@ceq.simple.only) {
        free.idx <- which(object@ParTable$free > 0L &&
          !duplicated(object@ParTable$free))
      } else {
        free.idx <- which(object@ParTable$free > 0L)
      }
      tmp.jac <- tmp.jac[free.idx, , drop = FALSE]
    }

    return.value <- tmp.jac %*% return.value %*% t(tmp.jac)

    # force return.value to be symmetric and pd
    return.value <- (return.value + t(return.value)) / 2
    # return.value <- lav_matrix_symmetric_force_pd(return.value,
    #                                     tol = 1e-09) # was 1e-06 < 0.6-9
  }

  # labels
  if (add.labels) {
    # if(rotation && !free.only) {
    #    # todo
    # } else {
    colnames(return.value) <- rownames(return.value) <-
      lav_partable_labels(object@ParTable, type = "free")
    # }
  }

  # alias?
  if (remove.duplicated && lavmodel@eq.constraints) {
    simple.flag <- lav_constraints_check_simple(lavmodel)
    if (simple.flag) {
      tmp.lab <- lav_partable_labels(object@ParTable, type = "free")
      dup.flag <- duplicated(tmp.lab)
      return.value <- return.value[!dup.flag, !dup.flag, drop = FALSE]
    } else {
      lav_msg_warn(
        gettext("alias is TRUE, but equality constraints do not appear
                to be simple; returning full vcov"))
    }
  }

  # class
  if (add.class) {
    class(return.value) <- c("lavaan.matrix.symmetric", "matrix")
  }

  return.value
}

lav_object_inspect_vcov_def <- function(object, joint = FALSE,
    standardized = FALSE, type = "std.all",
    add.labels = FALSE, add.class = FALSE) {

  lavmodel    <- object@Model
  lavpartable <- object@ParTable
  free.idx <- which(lavpartable$free > 0L)
  def.idx <- which(lavpartable$op == ":=")
  joint.idx <- c(free.idx, def.idx)

  if (!joint && length(def.idx) == 0L) {
    return(matrix(0, 0, 0))
  } else if (joint && length(joint.idx) == 0L) {
    return(matrix(0, 0, 0))
  }

  if (standardized) {
    # compute tmp.vcov for "free" parameters only
    tmp.vcov <- lav_object_inspect_vcov(object,
      standardized = TRUE,
      type = type, free.only = FALSE,
      add.labels = FALSE, add.class = FALSE)
    if (joint) {
      return.value <- tmp.vcov[joint.idx, joint.idx, drop = FALSE]
    } else {
      return.value <- tmp.vcov[def.idx, def.idx, drop = FALSE]
    }
  } else {
    # get free parameters
    x <- lav_model_get_parameters(lavmodel, type = "free")

    # bootstrap or not?
    if (!is.null(object@boot$coef)) {
      tmp.boot <- object@boot$coef
      # remove NA rows
      error.idx <- attr(tmp.boot, "error.idx")
      if (length(error.idx) > 0L) {
        tmp.boot <- tmp.boot[-error.idx, , drop = FALSE] # drops attributes
      }
      tmp.boot.def <- apply(tmp.boot, 1L, lavmodel@def.function)
      if (length(def.idx) == 1L) {
        tmp.boot.def <- as.matrix(tmp.boot.def)
      } else {
        tmp.boot.def <- t(tmp.boot.def)
      }
      return.value <- cov(tmp.boot.def)
    } else {
      # tmp.vcov
      tmp.vcov <- lav_object_inspect_vcov(object,
        standardized = FALSE,
        type = type, free.only = TRUE,
        add.labels = FALSE,
        add.class = FALSE)

      # regular delta method
      tmp.jac <- try(lav_func_jacobian_complex(func = lavmodel@def.function,
        x = x), silent = TRUE)
      if (inherits(tmp.jac, "try-error")) { # eg. pnorm()
        tmp.jac <- lav_func_jacobian_simple(func = lavmodel@def.function,
          x = x)
      }
      if (joint) {
        tmp.jac2 <- rbind(diag(nrow = ncol(tmp.jac)), tmp.jac)
        return.value <- tmp.jac2 %*% tmp.vcov %*% t(tmp.jac2)
      } else {
        return.value <- tmp.jac %*% tmp.vcov %*% t(tmp.jac)
      }
    }
  }

  # labels
  if (add.labels) {
    if (joint) {
      lhs.names <- lavpartable$lhs[joint.idx]
    } else {
      lhs.names <- lavpartable$lhs[def.idx]
    }
    colnames(return.value) <- rownames(return.value) <- lhs.names
  }

  # class
  if (add.class) {
    class(return.value) <- c("lavaan.matrix.symmetric", "matrix")
  }

  return.value
}

lav_object_inspect_UGamma <- function(object,                  # nolint
    add.labels = FALSE, add.class = FALSE) {

  out <- lav_test_satorra_bentler(lavobject     = object,
    method        = "original",
    return.ugamma = TRUE)
  return.value <- out$UGamma

  # labels
  # if(add.labels) {
  # colnames(return.value) <- rownames(return.value) <-
  # }

  # class
  if (add.class) {
    class(return.value) <- c("lavaan.matrix.symmetric", "matrix")
  }

  return.value
}

lav_object_inspect_UfromUGamma <- function(object,            # nolint
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {

  out <- lav_test_satorra_bentler(lavobject     = object,
    method        = "original",
    return.u      = TRUE)
  return.value <- out$UfromUGamma

  # labels
  # if(add.labels) {
  # colnames(return.value) <- rownames(return.value) <-
  # }

  # class
  if (add.class) {
    class(return.value) <- c("lavaan.matrix.symmetric", "matrix")
  }

  return.value
}


# Delta (jacobian: d samplestats / d free_parameters)
lav_object_inspect_delta <- function(object,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {

  lavmodel    <- object@Model
  lavdata     <- object@Data
  lavpartable <- object@ParTable
  lavpta      <- object@pta

  return.value <- lav_object_inspect_delta_internal(lavmodel = lavmodel,
    lavdata = lavdata, lavpartable = lavpartable,
    add.labels = add.labels, add.class = add.class,
    drop.list.single.group = drop.list.single.group)

  return.value
}

lav_object_inspect_delta_rownames <- function(                # nolint
    object,
    lavmodel = NULL,
    lavpartable = NULL,
    drop.list.single.group = FALSE) {
  if (!is.null(object)) {
    lavmodel    <- object@Model
    lavpartable <- object@ParTable
    lavpta      <- object@pta
    lavdata     <- object@Data
  } else {
    lavdata     <- NULL
    lavpta <- lav_partable_attributes(lavpartable)
  }

  categorical    <- lavmodel@categorical
  correlation    <- FALSE
  if (.hasSlot(lavmodel, "correlation")) {
    correlation    <- lavmodel@correlation
  }
  conditional.x  <- lavmodel@conditional.x
  group.w.free   <- lavmodel@group.w.free
  nvar           <- lavmodel@nvar
  num.idx        <- lavmodel@num.idx
  th.idx         <- lavmodel@th.idx
  nblocks        <- lavmodel@nblocks

  # store names per block, rbind later
  tmp.names <- vector("list", length = nblocks)

  # output is per group
  return.value <- vector("list", lavmodel@ngroups)

  for (g in 1:nblocks) {

    if (conditional.x) {
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
    # tmp <- apply(expand.grid(ov.names, ov.names), 1L,
    #             paste, collapse = "~~")

    # if(categorical) {
    #    names.cor <- tmp[lav_matrix_vech_idx(nvar, diagonal = FALSE)]
    #    names.var <- tmp[lav_matrix_diag_idx(nvar)[num.idx[[g]]]]
    # } else {
    #    names.cov <- tmp[lav_matrix_vech_idx(nvar, diagonal = TRUE)]
    # }

    # NOTE: in 0.6-1, we use the same order, but 'label' in row-wise
    # format (eg x1 ~~ x2 instead of x2 ~~ x1)
    tmp <- matrix(apply(expand.grid(ov.names, ov.names), 1L,
      paste, collapse = "~~"), nrow = nvar)
    if (categorical) {
      names.cor <- lav_matrix_vechru(tmp, diagonal = FALSE)
      names.var <- diag(tmp)[num.idx[[g]]]
    } else if (correlation) {
      names.cor <- lav_matrix_vechru(tmp, diagonal = FALSE)
    } else {
      names.cov <- lav_matrix_vechru(tmp, diagonal = TRUE)
    }


    # Mu
    if (!categorical && lavmodel@meanstructure) {
      names.mu <- paste(ov.names, "~1", sep = "")
    }

    # Pi
    if (conditional.x && lavmodel@nexo[g] > 0L) {
      names.pi <- apply(expand.grid(ov.names, ov.names.x), 1L,
        paste, collapse = "~")
    }

    # th
    if (categorical) {
      names.th <- lavpta$vnames$th[[g]]
      # interweave numeric intercepts, if any
      if (length(num.idx[[g]]) > 0L) {
        tmp <- character(length(th.idx[[g]]))
        tmp[th.idx[[g]] > 0] <- names.th
        tmp[th.idx[[g]] == 0] <- paste(ov.names[num.idx[[g]]],
          "~1", sep = "")
        names.th <- tmp
      }
    }

    # gw
    if (group.w.free) {
      names.gw <- "w"
    }

    tmp.names[[g]] <- c(names.gw,
      names.th, names.mu,
      names.pi,
      names.cov, names.var, names.cor)

  } # blocks

  # multilevel?
  if (.hasSlot(lavmodel, "multilevel") && lavmodel@multilevel) {
    for (g in 1:lavmodel@ngroups) {
      return.value[[g]] <- c(tmp.names[[(g - 1) * 2 + 1]],
        tmp.names[[(g - 1) * 2 + 2]])
    }
  } else {
    return.value <- tmp.names
  }

  # drop list?
  if (lavmodel@ngroups == 1L && drop.list.single.group) {
    return.value <- return.value[[1]]
  } else if (!is.null(lavdata)) {
    if (length(lavdata@group.label) > 0L) {
      names(return.value) <- unlist(lavdata@group.label)
    }
  }

  return.value
}

lav_object_inspect_delta_internal <- function(                   # nolint
    lavmodel = NULL, lavdata  = NULL,
    lavpartable = NULL,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {

  return.value <- computeDelta(lavmodel)

  if (add.labels) {
    tmp.pnames <- lav_partable_labels(lavpartable, type = "free")
    tmp.rownames  <- lav_object_inspect_delta_rownames(object = NULL,
      lavmodel = lavmodel, lavpartable = lavpartable,
      drop.list.single.group = FALSE)
  }

  for (g in seq_len(lavmodel@ngroups)) {
    # add labels
    if (add.labels) {
      colnames(return.value[[g]]) <- tmp.pnames
      rownames(return.value[[g]]) <- tmp.rownames[[g]]
    }

    # add class
    if (add.class) {
      class(return.value[[g]]) <- c("lavaan.matrix", "matrix")
    }

  } # ngroups

  # drop list?
  if (lavmodel@ngroups == 1L && drop.list.single.group) {
    return.value <- return.value[[1]]
  } else {
    if (length(lavdata@group.label) > 0L) {
      names(return.value) <- unlist(lavdata@group.label)
    }
  }

  return.value
}

lav_object_inspect_zero_cell_tables <-                          # nolint
  function(object, add.labels = FALSE, add.class = FALSE,
            drop.list.single.group = FALSE) {
  # categorical?
  if (!object@Model@categorical) {
    lav_msg_warn(gettext("no categorical variables in fitted model"))
    return(invisible(list()))
  }

  lavdata <- object@Data

  # create 2-way tables
  tmp.table <- lavTables(object, dimension = 2L, output = "data.frame",
    statistic = NULL)

  # select tables with empty cells
  empty.id <- tmp.table$id[which(tmp.table$obs.freq == 0)]


  if (length(empty.id) == 0L) {
    # only when lavInspect() is used, give message
    if (add.class) {
      cat("(There are no tables with empty cells for this fitted model)\n")
    }
    return(invisible(list()))
  } else {
    return.value <- lav_tables_cells_format(
      tmp.table[tmp.table$id %in% empty.id, ],
      lavdata = lavdata,
      drop.list.single.group = drop.list.single.group)
  }

  return.value
}

lav_object_inspect_coef <- function(object, type = "free",
                                    add.labels = FALSE, add.class = FALSE) {

  if (type == "user" || type == "all") {
    type <- "user"
    idx <- seq_along(object@ParTable$lhs)
  } else if (type == "free") {
    #idx <- which(object@ParTable$free > 0L & !duplicated(object@ParTable$free))
    idx <- which(object@ParTable$free > 0L)
  } else {
    lav_msg_stop(gettextf(
      "%1$s argument must be either %2$s or %3$s",
      "type", "free", "user"))
  }
  tmp.est <- lav_object_inspect_est(object)
  cof <- tmp.est[idx]

  # labels?
  if (add.labels) {
    names(cof) <- lav_partable_labels(object@ParTable, type = type)
  }

  # class
  if (add.class) {
    class(cof) <- c("lavaan.vector", "numeric")
  }

  cof
}

lav_object_inspect_npar <- function(object, type = "free") {

  if (type == "free") {
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
  n.g <- lavdata@ngroups
  return.value <- vector("list", n.g)

  # multilevel?
  if (lavdata@nlevels == 1L) {
    lav_msg_stop(gettext(
      "intraclass correlation only available for clustered data"))
  }

  if (length(object@h1) == 0L) {
    lav_msg_stop(gettext("h1 slot is not available; refit with h1 = TRUE"))
  }

  # implied statistics
  implied <- object@h1$implied

  for (g in 1:n.g) {
    sigma.w <- implied$cov[[(g - 1) * lavdata@nlevels + 1]]
    sigma.b <- implied$cov[[(g - 1) * lavdata@nlevels + 2]]

    w.diag <- diag(sigma.w)
    b.diag <- diag(sigma.b)

    return.value[[g]] <- numeric(length(w.diag))

    ov.names.l <- lavdata@ov.names.l[[g]]
    w.idx <- which(ov.names.l[[1]] %in% ov.names.l[[2]])
    w.names <- ov.names.l[[1]][w.idx]
    b.idx <- match(w.names, ov.names.l[[2]])

    return.value[[g]][w.idx] <- b.diag[b.idx] / (w.diag[w.idx] + b.diag[b.idx])

    # label
    if (add.labels) {
      names(return.value[[g]]) <- ov.names.l[[1]]
    }

    # class
    if (add.class) {
      class(return.value[[g]]) <- c("lavaan.vector", "numeric")
    }
  } # g

  if (n.g == 1L && drop.list.single.group) {
    return.value <- return.value[[1]]
  } else {
    if (length(object@Data@group.label) > 0L) {
      names(return.value) <- unlist(object@Data@group.label)
    }
  }

  return.value
}

lav_object_inspect_ranef <- function(object, add.labels = FALSE,
                                     add.class = FALSE,
                                     drop.list.single.group = FALSE) {

  lavdata <- object@Data
  lavsamplestats <- object@SampleStats

  n.g <- lavdata@ngroups
  return.value <- vector("list", n.g)

  # multilevel?
  if (lavdata@nlevels == 1L) {
    lav_msg_stop(gettext(
      "random effects only available for clustered data (in the long format)"))
  }

  # implied statistics
  lavimplied <- object@implied

  for (g in 1:n.g) {

    tmp.lp <- lavdata@Lp[[g]]
    tmp.ylp <- lavsamplestats@YLp[[g]]

    # implied for this group
    group.idx <- (g - 1) * lavdata@nlevels + seq_len(lavdata@nlevels)
    implied.group <- lapply(lavimplied, function(x) x[group.idx])

    # random effects (=random intercepts or cluster means)
    out <- lav_mvnorm_cluster_implied22l(Lp = tmp.lp, implied = implied.group)
    mb.j <- lav_mvnorm_cluster_em_estep_ranef(YLp = tmp.ylp, Lp = tmp.lp,
      sigma.w = out$sigma.w, sigma.b = out$sigma.b,
      sigma.zz = out$sigma.zz, sigma.yz = out$sigma.yz,
      mu.z = out$mu.z, mu.w = out$mu.w, mu.b = out$mu.b,
      se = FALSE)
    return.value[[g]] <- mb.j

    ov.names.l <- lavdata@ov.names.l[[g]]

    # label
    if (add.labels) {
      colnames(return.value[[g]]) <- ov.names.l[[1]]
    }

    # class
    if (add.class) {
      class(return.value[[g]]) <- c("lavaan.matrix", "matrix")

    }
  } # g

  if (n.g == 1L && drop.list.single.group) {
    return.value <- return.value[[1]]
  } else {
    if (length(object@Data@group.label) > 0L) {
      names(return.value) <- unlist(object@Data@group.label)
    }
  }

  return.value
}

# casewise loglikelihood contributions
lav_object_inspect_loglik_casewise <- function(object, log. = TRUE, # nolint
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {

  lavdata <- object@Data
  lavsamplestats <- object@SampleStats
  lavimplied <- object@implied
  lavoptions <- object@Options

  n.g <- lavdata@ngroups
  return.value <- vector("list", n.g)

  # multilevel?
  if (lavdata@nlevels > 1L) {
    lav_msg_stop(gettext(
  "casewise (log)likeloods contributions not yet available for clustered data"))
  }

  # estimator ML?
  if (object@Options$estimator != "ML") {
    lav_msg_stop(gettextf(
      "casewise (log)likeloods contributions only available for estimator = %s",
      dQuote("ML")))
  }

  for (g in 1:n.g) {

    if (lavsamplestats@missing.flag) {
      return.value[[g]] <-
        lav_mvnorm_missing_llik_casewise(Y = lavdata@X[[g]],
          wt = lavdata@weights[[g]],
          Mu = lavimplied$mean[[g]],
          Sigma = lavimplied$cov[[g]],
          x.idx  = lavsamplestats@x.idx[[g]])
    } else { # single-level, complete data
      if (lavoptions$conditional.x) {
        if (!is.null(lavdata@weights[[g]])) {
          lav_msg_stop(gettext("no support (yet) if weights are used."))
        }
        return.value[[g]] <- lav_mvreg_loglik_data(
          Y          = lavdata@X[[g]],
          eXo        = lavdata@eXo[[g]],
          res.int    = lavimplied$res.int[[g]],
          res.slopes = lavimplied$res.slopes[[g]],
          res.cov    = lavimplied$res.cov[[g]],
          casewise   = TRUE)

      } else {

        if (object@Model@meanstructure) {
          tmp.mean <- lavimplied$mean[[g]]
        } else {
          tmp.mean <- lavsamplestats@mean[[g]]
        }
        return.value[[g]] <-
          lav_mvnorm_loglik_data(Y  = lavdata@X[[g]],
            wt = lavdata@weights[[g]],
            Mu = tmp.mean,
            Sigma = lavimplied$cov[[g]],
            x.idx = lavsamplestats@x.idx[[g]],
            casewise = TRUE)
      }
    } # single-level, complete data

    # log. = FALSE?
    if (!log.) {
      return.value[[g]] <- exp(return.value[[g]])
    }

    # label
    # if(add.labels) {
    # }

    # class
    if (add.class) {
      class(return.value[[g]]) <- c("lavaan.vector", "numeric")

    }

  } # g

  if (n.g == 1L && drop.list.single.group) {
    return.value <- return.value[[1]]
  } else {
    if (length(object@Data@group.label) > 0L) {
      names(return.value) <- unlist(object@Data@group.label)
    }
  }

  return.value
}

# Mahalanobis distances for factor scores or casewise residuals
# type = "lv" -> factor scores
# type = "resid" -> casewise residuals
#
# we always use Bartlett factor scores (see Yuan & Hayashi 2010)
# (this has no impact on the m-distances for the factor scores,
#  and only a very slight impact on the m-distances for the casewise
#  residuals; but asymptotically, only when we use Bartlett factor
#  scores are the 'true scores' (=LAMBDA %*% FS) orthogonal to the
#  casewise residuals)
lav_object_inspect_mdist2 <- function(object, type = "resid", squared = TRUE,
                                      add.labels = FALSE, add.class = FALSE,
                                      drop.list.single.group = FALSE) {

  lavdata <- object@Data
  n.g <- lavdata@ngroups

  # lavPredict()
  out <- lavPredict(object, type = type, method = "ML", # = Bartlett
    label = FALSE, fsm = TRUE, mdist = TRUE,
    se = "none", acov = "none")
  return.value <- attr(out, "mdist")

  for (g in seq_len(n.g)) {
    # squared?
    if (!squared) {
      return.value[[g]] <- sqrt(return.value[[g]])
    }

    # labels?
    # if(add.labels) {
    # }

    # class
    if (add.class) {
      class(return.value[[g]]) <- c("lavaan.vector", "numeric")
    }
  } # g

  if (n.g == 1L && drop.list.single.group) {
    return.value <- return.value[[1]]
  } else {
    if (length(object@Data@group.label) > 0L) {
      names(return.value) <- unlist(object@Data@group.label)
    }
  }

  return.value
}
