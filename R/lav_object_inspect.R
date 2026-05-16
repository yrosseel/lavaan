# inspect a fitted lavaan object

# backward compatibility -- wrapper around lavInspect
lav_lavaan_inspect <- function(object, what = "free", ...) {
  dotdotdot <- list(...)
  if (length(dotdotdot) > 0L) {
    for (j in seq_along(dotdotdot)) {
      lav_msg_warn(gettextf(
        "Unknown argument %s for %s", sQuote(names(dotdotdot)[j]),
        sQuote("inspect"))
      )
    }
  }
  lav_lavaan_lavinspect(object              = object,
    what                   = what,
    add.labels             = TRUE,
    add.class              = TRUE,
    drop.list.single.group = TRUE)
}

# the `tech' version: no labels, full matrices, ... for further processing
lav_lavaan_lavtech <- function(object,                                  # nolint start
                           what                   = "free",
                           add.labels             = FALSE,
                           add.class              = FALSE,
                           list.by.group          = FALSE,
                           drop.list.single.group = FALSE) {            # nolint end

  lav_lavaan_lavinspect(object, what = what,
    add.labels = add.labels, add.class = add.class,
    list.by.group =  list.by.group,
    drop.list.single.group =  drop.list.single.group)
}

# the `user' version: with defaults for display only
lav_lavaan_lavinspect <- function(object,                                # nolint start
                              what                   = "free",
                              add.labels             = TRUE,
                              add.class              = TRUE,
                              list.by.group          = TRUE,
                              drop.list.single.group = TRUE) {           # nolint end
  # object must inherit from class lavaan
  stopifnot(inherits(object, "lavaan"))

  # check object
  object <- lav_object_check_version(object)

  # store partable with pta in object to use cache in called functions
  object@ParTable <- lav_partable_set_cache(object@ParTable, object@pta)

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
      type = "free", add_labels = add.labels, add_class = add.class,
      list_by_group =  list.by.group,
      drop_list_single_group = drop.list.single.group)
  } else if (what == "impute" ||
    what == "imputed") { # just to ease the transition for semTools!
    object@imputed
  } else if (what == "partable" || what == "user") {
    lav_object_inspect_modelmatrices(object, what = "free",
      type = "partable", add_labels = add.labels, add_class = add.class,
      list_by_group =  list.by.group,
      drop_list_single_group = drop.list.single.group)
  } else if (what == "se" ||
    what == "std.err" ||
    what == "standard.errors") {
    lav_object_inspect_modelmatrices(object, what = "se",
      add_labels = add.labels, add_class = add.class,
      list_by_group =  list.by.group,
      drop_list_single_group = drop.list.single.group)
  } else if (what == "se.std" ||
    what == "std.se") {
    lav_object_inspect_modelmatrices(object, what = "std.se",
      add_labels = add.labels, add_class = add.class,
      list_by_group =  list.by.group,
      drop_list_single_group = drop.list.single.group)
  } else if (what == "start" || what == "starting.values") {
    lav_object_inspect_modelmatrices(object, what = "start",
      add_labels = add.labels, add_class = add.class,
      list_by_group =  list.by.group,
      drop_list_single_group = drop.list.single.group)
  } else if (what == "est"  || what == "estimates" ||
    what == "x") {
    lav_object_inspect_modelmatrices(object, what = "est",
      add_labels = add.labels, add_class = add.class,
      list_by_group =  list.by.group,
      drop_list_single_group = drop.list.single.group)
  } else if (what == "est.unrotated") {
    lav_object_inspect_modelmatrices(object, what = "est.unrotated",
      add_labels = add.labels, add_class = add.class,
      list_by_group =  list.by.group,
      drop_list_single_group = drop.list.single.group)
  } else if (what == "dx.free") {
    lav_object_inspect_modelmatrices(object, what = "dx.free",
      add_labels = add.labels, add_class = add.class,
      list_by_group =  list.by.group,
      drop_list_single_group = drop.list.single.group)
  } else if (what == "dx.all") {
    lav_object_inspect_modelmatrices(object, what = "dx.all",
      add_labels = add.labels, add_class = add.class,
      list_by_group =  list.by.group,
      drop_list_single_group = drop.list.single.group)
  } else if (what == "std" || what == "std.all" ||
    what == "est.std" || what == "std.est" ||
    what == "standardized") {
    lav_object_inspect_modelmatrices(object, what = "std.all",
      add_labels = add.labels, add_class = add.class,
      list_by_group =  list.by.group,
      drop_list_single_group = drop.list.single.group)
  } else if (what == "std.lv") {
    lav_object_inspect_modelmatrices(object, what = "std.lv",
      add_labels = add.labels, add_class = add.class,
      list_by_group =  list.by.group,
      drop_list_single_group = drop.list.single.group)
  } else if (what == "std.nox") {
    lav_object_inspect_modelmatrices(object, what = "std.nox",
      add_labels = add.labels, add_class = add.class,
      list_by_group =  list.by.group,
      drop_list_single_group = drop.list.single.group)


    #### parameter table ####
  } else if (what == "list") {
    parTable(object)

    #### bootstrap coef ####
  } else if (what %in% c("boot", "bootstrap", "boot.coef", "coef.boot")) {
    lav_object_inspect_boot(object, add_labels = add.labels,
      add_class = add.class)

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
    out
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
      add_labels = add.labels, add_class = add.class,
      drop_list_single_group = drop.list.single.group)
  } else if (what == "obs.std" ||
    what == "observed.std" ||
    what == "sampstat.std" ||
    what == "sampstats.std" ||
    what == "samplestats.std" ||
    what == "samp.std" ||
    what == "sample.std" ||
    what == "samplestatistics.std") {
    lav_object_inspect_sampstat(object, h1 = TRUE, std = TRUE,
      add_labels = add.labels, add_class = add.class,
      drop_list_single_group = drop.list.single.group)
  } else if (what == "h1" || what == "missing.h1" || what == "sampstat.h1") {
    lav_object_inspect_sampstat(object, h1 = TRUE,
      add_labels = add.labels, add_class = add.class,
      drop_list_single_group = drop.list.single.group)

    #### wls.est - wls.obs - wls.v ####
  } else if (what == "wls.est") {
    lav_object_inspect_wls_est(object,
      add_labels = add.labels, add_class = add.class,
      drop_list_single_group = drop.list.single.group)
  } else if (what == "wls.obs") {
    lav_object_inspect_wls_obs(object,
      add_labels = add.labels, add_class = add.class,
      drop_list_single_group = drop.list.single.group)
  } else if (what == "wls.v") {
    lav_object_inspect_wls_v(object,
      add_labels = add.labels, add_class = add.class,
      drop_list_single_group = drop.list.single.group)



    #### data + missingness ####
  } else if (what == "data") {
    lav_object_inspect_data(object, add_labels = add.labels,
      drop_list_single_group = drop.list.single.group)
  } else if (what == "case.idx") {
    lav_object_inspect_case_idx(object,
      drop_list_single_group = drop.list.single.group)
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
      drop_list_single_group = drop.list.single.group)
  } else if (what == "ncluster.size") {
    lav_object_inspect_cluster_info(object, level = 2L,
      what = "ncluster.size",
      drop_list_single_group = drop.list.single.group)
  } else if (what == "cluster.size") {
    lav_object_inspect_cluster_info(object, level = 2L,
      what = "cluster.size",
      drop_list_single_group = drop.list.single.group)
  } else if (what == "cluster.id") {
    lav_object_inspect_cluster_info(object, level = 2L,
      what = "cluster.id",
      drop_list_single_group = drop.list.single.group)
  } else if (what == "cluster.idx") {
    lav_object_inspect_cluster_info(object, level = 2L,
      what = "cluster.idx",
      drop_list_single_group = drop.list.single.group)
  } else if (what == "cluster.label") {
    lav_object_inspect_cluster_info(object, level = 2L,
      what = "cluster.label",
      drop_list_single_group = drop.list.single.group)
  } else if (what == "cluster.sizes") {
    lav_object_inspect_cluster_info(object, level = 2L,
      what = "cluster.sizes",
      drop_list_single_group = drop.list.single.group)
  } else if (what == "average.cluster.size") {
    lav_object_inspect_cluster_info(object, level = 2L,
      what = "average.cluster.size",
      drop_list_single_group = drop.list.single.group)
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
      add_labels = add.labels, add_class = add.class,
      drop_list_single_group = drop.list.single.group)
  } else if (what %in% c("patterns", "pattern")) {
    lav_object_inspect_missing_patterns(object,
      add_labels = add.labels, add_class = add.class,
      drop_list_single_group = drop.list.single.group)
  } else if (what == "empty.idx") {
    lav_object_inspect_empty_idx(object,
      drop_list_single_group = drop.list.single.group)


    #### rsquare ####
  } else if (what == "rsquare" || what == "r-square" || what == "r2") {
    lav_object_inspect_rsquare(object,
      add_labels = add.labels, add_class = add.class,
      drop_list_single_group = drop.list.single.group)


    #### model-implied sample statistics ####
  } else if (what == "implied" || what == "fitted" ||
    what == "expected" || what == "exp") {
    lav_object_inspect_implied(object,
      add_labels = add.labels, add_class = add.class,
      drop_list_single_group = drop.list.single.group)
  } else if (what == "resid" || what == "res" || what == "residual" ||
    what == "residuals") {
    lav_object_inspect_residuals(object, h1 = TRUE,
      add_labels = add.labels, add_class = add.class,
      drop_list_single_group = drop.list.single.group)
  } else if (what == "cov.lv" || what == "veta") {
    lav_object_inspect_cov_lv(object,
      correlation_metric = FALSE,
      add_labels = add.labels, add_class = add.class,
      drop_list_single_group = drop.list.single.group)
  } else if (what == "cor.lv") {
    lav_object_inspect_cov_lv(object,
      correlation_metric = TRUE,
      add_labels = add.labels, add_class = add.class,
      drop_list_single_group = drop.list.single.group)
  } else if (what == "mean.lv" || what == "eeta") {
    lav_object_inspect_mean_lv(object,
      add_labels = add.labels, add_class = add.class,
      drop_list_single_group = drop.list.single.group)
  } else if (what == "cov.all") {
    lav_object_inspect_cov_all(object,
      correlation_metric = FALSE,
      add_labels = add.labels, add_class = add.class,
      drop_list_single_group = drop.list.single.group)
  } else if (what == "cor.all") {
    lav_object_inspect_cov_all(object,
      correlation_metric = TRUE,
      add_labels = add.labels, add_class = add.class,
      drop_list_single_group = drop.list.single.group)
  } else if (what == "cov.ov" || what == "sigma" || what == "sigma.hat") {
    lav_object_inspect_cov_ov(object,
      correlation_metric = FALSE,
      add_labels = add.labels, add_class = add.class,
      drop_list_single_group = drop.list.single.group)
  } else if (what == "cor.ov") {
    lav_object_inspect_cov_ov(object,
      correlation_metric = TRUE,
      add_labels = add.labels, add_class = add.class,
      drop_list_single_group = drop.list.single.group)
  } else if (what == "mean.ov" || what == "mu" || what == "mu.hat") {
    lav_object_inspect_mean_ov(object,
      add_labels = add.labels, add_class = add.class,
      drop_list_single_group = drop.list.single.group)
  } else if (what == "th" || what == "thresholds") {
    lav_object_inspect_th(object,
      add_labels = add.labels, add_class = add.class,
      drop_list_single_group = drop.list.single.group)
  } else if (what == "th.idx") {
    lav_object_inspect_th_idx(object,
      add_labels = add.labels, add_class = add.class,
      drop_list_single_group = drop.list.single.group)
  } else if (what == "vy") {
    lav_object_inspect_vy(object,
      add_labels = add.labels, add_class = add.class,
      drop_list_single_group = drop.list.single.group)
  } else if (what %in% c("fs.reliability", "fs.rel", "fs.reliabilities")) {
    lav_object_inspect_fs_determinacy(object, squared = TRUE,
      fs_method = "regression",
      add_labels = add.labels, add_class = add.class,
      drop_list_single_group = drop.list.single.group)
  } else if (what %in% c("fs.determinacy", "fs.det", "fs.determin",
                         "fs.determinacies")) {
    lav_object_inspect_fs_determinacy(object, squared = FALSE,
      fs_method = "regression",
      add_labels = add.labels, add_class = add.class,
      drop_list_single_group = drop.list.single.group)
  } else if (what %in% c("fs.reliability.bartlett",
                         "fs.reliability.Bartlett",
                         "fs.rel.bartlett", "fs.rel.Bartlett",
                         "fs.reliabilities.bartlett",
                         "fs.reliabilities.Bartlett")) {
    lav_object_inspect_fs_determinacy(object, squared = TRUE,
      fs_method = "Bartlett",
      add_labels = add.labels, add_class = add.class,
      drop_list_single_group = drop.list.single.group)
  } else if (what %in% c("fs.determinacy.bartlett",
                         "fs.determinacy.Bartlett",
                         "fs.det.bartlett", "fs.det.Bartlett",
                         "fs.determin.bartlett", "fs.determin.Bartlett",
                         "fs.determinacies.bartlett",
                         "fs.determinacies.Bartlett")) {
    lav_object_inspect_fs_determinacy(object, squared = FALSE,
      fs_method = "Bartlett",
      add_labels = add.labels, add_class = add.class,
      drop_list_single_group = drop.list.single.group)



    #### specific model matrices? ####
  } else if (what == "theta" || what == "theta.cov") {
    lav_object_inspect_theta(object,  correlation_metric = FALSE,
      add_labels = add.labels, add_class = add.class,
      drop_list_single_group = drop.list.single.group)
  } else if (what == "theta.cor") {
    lav_object_inspect_theta(object,  correlation_metric = TRUE,
      add_labels = add.labels, add_class = add.class,
      drop_list_single_group = drop.list.single.group)

    #### (squared) Mahalanobis distances ####
  } else if (what == "mdist2.fs") {
    lav_object_inspect_mdist2(object, type = "lv", squared = TRUE,
      add_labels = add.labels, add_class = add.class,
      drop_list_single_group = drop.list.single.group)
  } else if (what == "mdist2.resid") {
    lav_object_inspect_mdist2(object, type = "resid", squared = TRUE,
      add_labels = add.labels, add_class = add.class,
      drop_list_single_group = drop.list.single.group)
  } else if (what == "mdist.fs") {
    lav_object_inspect_mdist2(object, type = "lv", squared = FALSE,
      add_labels = add.labels, add_class = add.class,
      drop_list_single_group = drop.list.single.group)
  } else if (what == "mdist.resid") {
    lav_object_inspect_mdist2(object, type = "resid", squared = FALSE,
      add_labels = add.labels, add_class = add.class,
      drop_list_single_group = drop.list.single.group)

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
    lav_object_inspect_npar(object, ceq = FALSE) # ignore equality constraints
  } else if (what == "coef") {
    # this breaks simsem and semTools -- 0.6-1
    # lav_object_inspect_coef(object, type = "free",
    #                        add_labels = add.labels, add_class = add.class)
    lav_object_inspect_modelmatrices(object, what = "est",
      type = "free", add_labels = add.labels, add_class = add.class,
      list_by_group =  list.by.group,
      drop_list_single_group = drop.list.single.group)


    #### NACOV samplestats ####
  } else if (what == "gamma") {
    lav_object_inspect_sampstat_gamma(object,
      add_labels = add.labels, add_class = add.class,
      drop_list_single_group = drop.list.single.group)


    #### gradient, Hessian, information, first.order, vcov ####
  } else if (what == "gradient") {
    lav_object_inspect_gradient(object,
      add_labels = add.labels, add_class = add.class, logl = FALSE)
  } else if (what == "gradient.logl") {
    lav_object_inspect_gradient(object,
      add_labels = add.labels, add_class = add.class, logl = TRUE)
  } else if (what == "optim.gradient") {
    lav_object_inspect_gradient(object,
      add_labels = add.labels, add_class = add.class, optim = TRUE)
  } else if (what == "hessian") {
    lav_object_inspect_hessian(object,
      add_labels = add.labels, add_class = add.class)

  } else if (what == "information") {
    lav_object_inspect_information(object, information = "default",
      augmented = FALSE, inverted = FALSE,
      add_labels = add.labels, add_class = add.class)
  } else if (what == "information.expected") {
    lav_object_inspect_information(object, information = "expected",
      augmented = FALSE, inverted = FALSE,
      add_labels = add.labels, add_class = add.class)
  } else if (what == "information.observed") {
    lav_object_inspect_information(object, information = "observed",
      augmented = FALSE, inverted = FALSE,
      add_labels = add.labels, add_class = add.class)
  } else if (what == "information.first.order" ||
    what == "information.firstorder"  ||
    what == "first.order") {
    lav_object_inspect_information(object, information = "first.order",
      augmented = FALSE, inverted = FALSE,
      add_labels = add.labels, add_class = add.class)

  } else if (what == "augmented.information") {
    lav_object_inspect_information(object, information = "default",
      augmented = TRUE, inverted = FALSE,
      add_labels = add.labels, add_class = add.class)
  } else if (what == "augmented.information.expected") {
    lav_object_inspect_information(object, information = "expected",
      augmented = TRUE, inverted = FALSE,
      add_labels = add.labels, add_class = add.class)
  } else if (what == "augmented.information.observed") {
    lav_object_inspect_information(object, information = "observed",
      augmented = TRUE, inverted = FALSE,
      add_labels = add.labels, add_class = add.class)
  } else if (what == "augmented.information.first.order" ||
    what == "augmented.first.order") {
    lav_object_inspect_information(object, information = "first.order",
      augmented = TRUE, inverted = FALSE,
      add_labels = add.labels, add_class = add.class)

  } else if (what == "inverted.information") {
    lav_object_inspect_information(object, information = "default",
      augmented = TRUE, inverted = TRUE,
      add_labels = add.labels, add_class = add.class)
  } else if (what == "inverted.information.expected") {
    lav_object_inspect_information(object, information = "expected",
      augmented = TRUE, inverted = TRUE,
      add_labels = add.labels, add_class = add.class)
  } else if (what == "inverted.information.observed") {
    lav_object_inspect_information(object, information = "observed",
      augmented = TRUE, inverted = TRUE,
      add_labels = add.labels, add_class = add.class)
  } else if (what == "inverted.information.first.order" ||
    what == "inverted.first.order") {
    lav_object_inspect_information(object, information = "first.order",
      augmented = TRUE, inverted = TRUE,
      add_labels = add.labels, add_class = add.class)

  } else if (what == "h1.information") {
    lav_object_inspect_h1_information(object, information = "default",
      h1_information = "default", inverted = FALSE,
      add_labels = add.labels, add_class = add.class,
      drop_list_single_group = drop.list.single.group)
  } else if (what == "h1.information.expected") {
    lav_object_inspect_h1_information(object, information = "expected",
      h1_information = "default", inverted = FALSE,
      add_labels = add.labels, add_class = add.class,
      drop_list_single_group = drop.list.single.group)
  } else if (what == "h1.information.observed") {
    lav_object_inspect_h1_information(object, information = "observed",
      h1_information = "default", inverted = FALSE,
      add_labels = add.labels, add_class = add.class,
      drop_list_single_group = drop.list.single.group)
  } else if (what == "h1.information.first.order" ||
    what == "h1.information.firstorder"  ||
    what == "h1.first.order") {
    lav_object_inspect_h1_information(object,
      information = "first.order", h1_information = "default",
      inverted = FALSE, add_labels = add.labels, add_class = add.class,
      drop_list_single_group = drop.list.single.group)

  } else if (what == "vcov") {
    lav_object_inspect_vcov(object,
      standardized = FALSE,
      add_labels = add.labels, add_class = add.class)
  } else if (what == "vcov.std.all" || what == "vcov.standardized" ||
    what == "vcov.std") {
    lav_object_inspect_vcov(object,
      standardized = TRUE, type = "std.all",
      add_labels = add.labels, add_class = add.class)
  } else if (what == "vcov.std.lv") {
    lav_object_inspect_vcov(object,
      standardized = TRUE, type = "std.lv",
      add_labels = add.labels, add_class = add.class)
  } else if (what == "vcov.std.nox") {
    lav_object_inspect_vcov(object,
      standardized = TRUE, type = "std.nox",
      add_labels = add.labels, add_class = add.class)

  } else if (what == "vcov.def") {
    lav_object_inspect_vcov_def(object, joint = FALSE,
      standardized = FALSE,
      add_labels = add.labels, add_class = add.class)
  } else if (what == "vcov.def.std.all" || what == "vcov.def.standardized" ||
    what == "vcov.def.std") {
    lav_object_inspect_vcov_def(object, joint = FALSE,
      standardized = TRUE, type = "std.all",
      add_labels = add.labels, add_class = add.class)
  } else if (what == "vcov.def.std.lv") {
    lav_object_inspect_vcov_def(object, joint = FALSE,
      standardized = TRUE, type = "std.lv",
      add_labels = add.labels, add_class = add.class)
  } else if (what == "vcov.def.std.nox") {
    lav_object_inspect_vcov_def(object, joint = FALSE,
      standardized = TRUE, type = "std.nox",
      add_labels = add.labels, add_class = add.class)

  } else if (what == "vcov.def.joint") {
    lav_object_inspect_vcov_def(object, joint = TRUE,
      standardized = FALSE,
      add_labels = add.labels, add_class = add.class)
  } else if (what == "vcov.def.joint.std.all" ||
    what == "vcov.def.joint.standardized" ||
    what == "vcov.def.joint.std") {
    lav_object_inspect_vcov_def(object, joint = TRUE,
      standardized = TRUE, type = "std.all",
      add_labels = add.labels, add_class = add.class)
  } else if (what == "vcov.def.joint.std.lv") {
    lav_object_inspect_vcov_def(object, joint = TRUE,
      standardized = TRUE, type = "std.lv",
      add_labels = add.labels, add_class = add.class)
  } else if (what == "vcov.def.joint.std.nox") {
    lav_object_inspect_vcov_def(object, joint = TRUE,
      standardized = TRUE, type = "std.nox",
      add_labels = add.labels, add_class = add.class)

  } else if (what == "ugamma" || what == "ug" || what == "u.gamma") {
    lav_object_inspect_ugamma(object,
      add_labels = add.labels, add_class = add.class)
  } else if (what == "ufromugamma" || what == "u") {
    lav_object_inspect_u_from_ugamma(object,
      add_labels = add.labels, add_class = add.class,
      drop_list_single_group = drop.list.single.group)

    ### jacobians ####
  } else if (what == "delta") {
    lav_object_inspect_delta(object,
      add_labels = add.labels, add_class = add.class,
      drop_list_single_group = drop.list.single.group)
  } else if (what == "delta.rownames") {
    lav_object_inspect_delta_rownames(object,
      drop_list_single_group = drop.list.single.group)

    ### casewise loglikelihoods ###
  } else if (what == "loglik.casewise") {
    lav_object_inspect_loglik_casewise(object, log_1 = TRUE,
      add_labels = add.labels, add_class = add.class,
      drop_list_single_group = drop.list.single.group)
  } else if (what == "lik.casewise") {
    lav_object_inspect_loglik_casewise(object, log_1 = FALSE,
      add_labels = add.labels, add_class = add.class,
      drop_list_single_group = drop.list.single.group)

    # multilevel #
  } else if (what == "icc") {
    lav_object_inspect_icc(object,
      add_labels = add.labels, add_class = add.class,
      drop_list_single_group = drop.list.single.group)
  } else if (what == "ranef") {
    lav_object_inspect_ranef(object,
      add_labels = add.labels, add_class = add.class,
      drop_list_single_group = drop.list.single.group)

    # instrumental variables
  } else if (what %in% c("iv", "ivs", "miiv", "miivs", "instr", "instruments")) {
    lav_object_inspect_iv(object,
      drop_list_single_group = drop.list.single.group)
  } else if (what %in% c("eqs")) {
    lav_object_inspect_eqs(object,
      drop_list_single_group = drop.list.single.group)
  } else if (what %in% c("sargan")) {
    lav_object_inspect_sargan(object,
      drop_list_single_group = drop.list.single.group)

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
      add_labels = add.labels, add_class = add.class,
      drop_list_single_group = drop.list.single.group)

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
        return_value <- object@ParTable$est.unrotated
      } else {
        return_value <- object@ParTable$est 
                  # if this changes, tag @TDJorgensen in commit message
      }
    } else {
      partable <- parTable(object)
      return_value <- rep(as.numeric(NA), length(partable$lhs))
    }
  } else {
    # try generic coef()
    return_value <- coef(object, type = "user")
    if (is.matrix(return_value)) {
      # lavaanList?
      return_value <- rowMeans(return_value)
    }
  }

  return_value
}

lav_object_inspect_se <- function(object) {
  # from 0.5-19, they are in the partable
  if (!is.null(object@ParTable$se)) {
    return_value <- object@ParTable$se
  } else {
    partable <- parTable(object)
    return_value <- rep(as.numeric(NA), length(partable$lhs))
  }

  return_value
}

lav_object_inspect_std_se <- function(object) {

  if (!is.null(object@ParTable$se.std)) {
    return_value <- object@ParTable$se.std
  } else {
    tmp_std <- standardizedSolution(object)
    return_value <- tmp_std$se
  }

  return_value
}

lav_object_inspect_start <- function(object) {
  # from 0.5-19, they are in the partable
  if (!is.null(object@ParTable$start)) {
    return_value <- object@ParTable$start
  } else {
    # in < 0.5-19, we should look in @Fit@start
    return_value <- object@Fit@start
  }

  return_value
}

lav_object_inspect_boot <- function(object, add_labels = FALSE,
                                    add_class = FALSE) {

  if (object@Options$se   != "bootstrap" &&
    !any(c("bootstrap", "bollen.stine") %in% object@Options$test)) {
    lav_msg_stop(gettext("bootstrap was not used."))
  }

  # from 0.5-19. they are in a separate slot
  tmp <- try(slot(object, "boot"), silent = TRUE)
  if (inherits(tmp, "try-error")) {
    # older version of object?
    est <- lav_object_inspect_est(object)
    tmp_boot <- attr(est, "tmp.boot.COEF")
  } else {
    # 0.5-19 way
    tmp_boot <- object@boot$coef
  }

  # add coef names
  if (add_labels) {
    colnames(tmp_boot) <- names(coef(object))
  }

  # add class
  if (add_class) {
    class(tmp_boot) <- c("lavaan.matrix", "matrix")
  }

  tmp_boot
}


lav_object_inspect_modelmatrices <- function(object, what = "free", # nolint
    type = "free", add_labels = FALSE, add_class = FALSE,
    list_by_group = FALSE,
    drop_list_single_group = FALSE) {

  glist <- object@Model@GLIST

  current_verbose <- lav_verbose()
  if (what == "dx.free") {
    if (lav_verbose(FALSE)) on.exit(lav_verbose(current_verbose), TRUE)
    tmp_dx <- lav_model_gradient(
      lavmodel       = object@Model,
      glist          = NULL,
      lavsamplestats = object@SampleStats,
      lavdata        = object@Data,
      lavcache       = object@Cache,
      type           = "free",
      group_weight   = TRUE,
      ceq_simple     = TRUE,
      delta          = NULL)
  } else if (what == "dx.all") {
    if (lav_verbose(FALSE)) on.exit(lav_verbose(current_verbose), TRUE)
    glist <- lav_model_gradient(
      lavmodel   = object@Model,
      glist          = NULL,
      lavsamplestats = object@SampleStats,
      lavdata        = object@Data,
      lavcache       = object@Cache,
      type           = "allofthem",
      group_weight   = TRUE,
      ceq_simple     = FALSE,
      delta          = NULL)
    names(glist) <- names(object@Model@GLIST)
  } else if (what == "std.all") {
    tmp_std <- lav_standardize_all(object)
  } else if (what == "std.lv") {
    tmp_std <- lav_standardize_lv(object)
  } else if (what == "std.nox") {
    tmp_std <- lav_standardize_all_nox(object)
  } else if (what == "se") {
    tmp_se <- lav_object_inspect_se(object)
  } else if (what == "std.se") {
    tmp_se <- lav_object_inspect_std_se(object)
  } else if (what == "start") {
    tmp_start <- lav_object_inspect_start(object)
  } else if (what == "est") {
    tmp_est <- lav_object_inspect_est(object)
  } else if (what == "est.unrotated") {
    if (!is.null(object@Options$rotation) &&
      object@Options$rotation == "none") {
      tmp_est <- lav_object_inspect_est(object, unrotated = FALSE)
    } else {
      tmp_est <- lav_object_inspect_est(object, unrotated = TRUE)
    }
  }

  for (mm in seq_along(glist)) {

    if (add_labels) {
      dimnames(glist[[mm]]) <- object@Model@dimNames[[mm]]
    }

    if (what == "free") {
      # fill in free parameter counts
      if (type == "free") {
        m_el_idx <- object@Model@m.free.idx[[mm]]
        x_el_idx <- object@Model@x.free.idx[[mm]]
        # } else if(type == "unco") {
        #    m.el.idx <- object@Model@m.unco.idx[[mm]]
        #    x.el.idx <- object@Model@x.unco.idx[[mm]]
      } else if (type == "partable") {
        m_el_idx <- object@Model@m.user.idx[[mm]]
        x_el_idx <- object@Model@x.user.idx[[mm]]
      } else {
        lav_msg_stop(gettextf(
          "%1$s argument unknown: %2$s",
          "type", lav_msg_view(type)
        ))
      }
      # erase everything
      glist[[mm]][, ] <- 0.0
      glist[[mm]][m_el_idx] <- x_el_idx
    } else if (what == "se" || what == "std.se") {
      # fill in standard errors
      m_user_idx <- object@Model@m.user.idx[[mm]]
      x_user_idx <- object@Model@x.user.idx[[mm]]
      # erase everything
      glist[[mm]][, ] <- 0.0
      glist[[mm]][m_user_idx] <- tmp_se[x_user_idx]
    } else if (what == "start") {
      # fill in starting values
      m_user_idx <- object@Model@m.user.idx[[mm]]
      x_user_idx <- object@Model@x.user.idx[[mm]]
      glist[[mm]][m_user_idx] <- tmp_start[x_user_idx]
    } else if (what %in% c("est", "est.unrotated")) {
      # fill in estimated parameter values
      m_user_idx <- object@Model@m.user.idx[[mm]]
      x_user_idx <- object@Model@x.user.idx[[mm]]
      glist[[mm]][m_user_idx] <- tmp_est[x_user_idx]
    } else if (what == "dx.free") {
      # fill in derivatives free parameters
      m_el_idx <- object@Model@m.free.idx[[mm]]
      x_el_idx <- object@Model@x.free.idx[[mm]]
      # erase everything
      glist[[mm]][, ] <- 0.0
      glist[[mm]][m_el_idx] <- tmp_dx[x_el_idx]
    } else if (what %in% c("std.all", "std.lv", "std.nox")) {
      m_user_idx <- object@Model@m.user.idx[[mm]]
      x_user_idx <- object@Model@x.user.idx[[mm]]
      glist[[mm]][m_user_idx] <- tmp_std[x_user_idx]
    }

    # class
    if (add_class) {
      if (object@Model@isSymmetric[mm]) {
        class(glist[[mm]]) <- c("lavaan.matrix.symmetric", "matrix")
      } else {
        class(glist[[mm]]) <- c("lavaan.matrix", "matrix")
      }
    }
  }

  # try to reflect `equality constraints'
  con_flag <- FALSE
  if (what == "free" && object@Model@eq.constraints) {
    # extract constraints from parameter table
    partable <- parTable(object)
    tmp_con <-  partable[partable$op %in% c("==", "<", ">"),
                         c("lhs", "op", "rhs")]
    rownames(tmp_con) <- NULL

    # replace 'labels' by parameter numbers
    tmp_id <- lav_partable_constraints_label_id(partable)
    tmp_label <- names(tmp_id)
    for (con in seq_len(nrow(tmp_con))) {
      # lhs
      lhs_labels <- all.vars(as.formula(paste("~", tmp_con[con, "lhs"])))

      if (length(lhs_labels) > 0L) {
        # par id
        lhs_freeid <- tmp_id[match(lhs_labels, tmp_label)]

        # substitute
        tmp <- tmp_con[con, "lhs"]
        for (pat in seq_along(lhs_labels)) {
          tmp <- sub(lhs_labels[pat], lhs_freeid[pat], tmp)
        }
        tmp_con[con, "lhs"] <- tmp
      }

      # rhs
      rhs_labels <- all.vars(as.formula(paste("~", tmp_con[con, "rhs"])))

      if (length(rhs_labels) > 0L) {
        # par id
        rhs_freeid <- tmp_id[match(rhs_labels, tmp_label)]
        # substitute
        tmp <- tmp_con[con, "rhs"]
        for (pat in seq_along(rhs_labels)) {
          tmp <- sub(rhs_labels[pat], rhs_freeid[pat], tmp)
        }
        tmp_con[con, "rhs"] <- tmp
      }
    } # con

    # add this info at the top
    # glist <- c(constraints = list(tmp.con), glist)
    # no, not a good idea, it does not work with list.by.group

    # add it as a 'header' attribute?
    attr(tmp_con, "header") <- "Note: model contains equality constraints:"
    con_flag <- TRUE
  }

  # should we group them per block?
  if (list_by_group) {
    lavmodel       <- object@Model
    nmat           <- lavmodel@nmat

    return_value <- vector("list", length = lavmodel@nblocks)
    for (b in seq_len(lavmodel@nblocks)) {
      # which mm belong to this block?
      mm_in_group <- 1:nmat[b] + cumsum(c(0, nmat))[b]

      return_value[[b]] <- glist[mm_in_group]
    }

    if (lavmodel@nblocks == 1L && drop_list_single_group) {
      return_value <- return_value[[1]]
    } else if (lavmodel@nblocks > 1L) {
      names(return_value) <- object@Data@block.label
    }
  } else {
    return_value <- glist
  }

  # header
  if (con_flag) {
    attr(return_value, "header") <- tmp_con
  }

  # lavaan.list
  if (add_class) {
    class(return_value) <- c("lavaan.list", "list")
  }

  return_value
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
    add_labels = FALSE, add_class = FALSE, drop_list_single_group = FALSE) {

  nblocks <- object@Model@nblocks
  ov_names <- object@pta$vnames$ov
  ov_names_res <- object@pta$vnames$ov.nox
  ov_names_x   <- object@pta$vnames$ov.x

  # slots
  lavsamplestats <- object@SampleStats
  lavmodel <- object@Model

  # if nlevels, override h1 to be TRUE, and set conditional.x = FALSE
  if (object@Data@nlevels > 1L) {
    h1 <- TRUE
    conditional_x <- FALSE # for now (0.6-12)
  } else {
    conditional_x <- lavmodel@conditional.x
  }

  # check if we have a non-empty @h1 slot
  if (length(object@h1) == 0L) {
    h1 <- FALSE
  } else {
    h1_implied <- object@h1$implied
  }

  # if h1 = FALSE and nlevels > 1L, nothing to show...
  if (!h1 && object@Data@nlevels > 1L) {
    lav_msg_stop(gettext(
      "sample statistics not available; refit with option h1 = TRUE"))
  }

  return_value <- vector("list", length = nblocks)
  for (b in seq_len(nblocks)) {

    if (!conditional_x) {
      # covariance matrix
      if (h1) {
        return_value[[b]]$cov  <- h1_implied$cov[[b]]
      } else {
        return_value[[b]]$cov  <- lavsamplestats@cov[[b]]
      }
      if (std) {
        diag_orig <- diag(return_value[[b]]$cov)
        return_value[[b]]$cov <- cov2cor(return_value[[b]]$cov)
      }
      if (add_labels && !is.null(return_value[[b]]$cov)) {
        rownames(return_value[[b]]$cov) <- colnames(return_value[[b]]$cov) <-
          ov_names[[b]]
      }
      if (add_class) {
        class(return_value[[b]]$cov) <- c("lavaan.matrix.symmetric", "matrix")
      }

      # mean vector
      if (lavmodel@meanstructure) {
        if (h1) {
          return_value[[b]]$mean <- as.numeric(h1_implied$mean[[b]])
        } else {
          return_value[[b]]$mean <- as.numeric(lavsamplestats@mean[[b]])
        }
        if (std) {
          diag_orig[diag_orig < .Machine$double.eps] <- NA
          return_value[[b]]$mean <- return_value[[b]]$mean / sqrt(diag_orig)
        }
        if (add_labels) {
          names(return_value[[b]]$mean) <- ov_names[[b]]
        }
        if (add_class) {
          class(return_value[[b]]$mean) <- c("lavaan.vector", "numeric")
        }
      }

      # thresholds
      if (lavmodel@categorical) {
        if (h1) {
          return_value[[b]]$th <- as.numeric(h1_implied$th[[b]])
        } else {
          return_value[[b]]$th <- as.numeric(lavsamplestats@th[[b]])
        }
        if (length(lavmodel@num.idx[[b]]) > 0L) {
          num_idx <- which(lavmodel@th.idx[[b]] == 0)
          return_value[[b]]$th <- return_value[[b]]$th[-num_idx]
        }
        # FIXME: what to do if std = TRUE (depends on delta/theta)
        if (add_labels) {
          names(return_value[[b]]$th) <- object@pta$vnames$th[[b]]
        }
        if (add_class) {
          class(return_value[[b]]$th) <- c("lavaan.vector", "numeric")
        }
      }
      # !conditional.x
    } else {
      # if conditional.x = TRUE

      # residual covariance matrix
      if (h1) {
        return_value[[b]]$res.cov  <- h1_implied$res.cov[[b]]
      } else {
        return_value[[b]]$res.cov  <- lavsamplestats@res.cov[[b]]
      }
      if (std) {
        diag_orig <- diag(return_value[[b]]$res.cov)
        return_value[[b]]$res.cov <- cov2cor(return_value[[b]]$res.cov)
      }
      if (add_labels) {
        rownames(return_value[[b]]$res.cov) <-
          colnames(return_value[[b]]$res.cov) <-
          ov_names_res[[b]]
      }
      if (add_class) {
        class(return_value[[b]]$res.cov) <-
          c("lavaan.matrix.symmetric", "matrix")
      }

      # intercepts
      if (lavmodel@meanstructure) {
        if (h1) {
          return_value[[b]]$res.int <- as.numeric(h1_implied$res.int[[b]])
        } else {
          return_value[[b]]$res.int <- as.numeric(lavsamplestats@res.int[[b]])
        }
        if (std) {
          diag_orig[diag_orig < .Machine$double.eps] <- NA
          return_value[[b]]$res.int <- return_value[[b]]$res.int /
            sqrt(diag_orig)
        }
        if (add_labels) {
          names(return_value[[b]]$res.int) <- ov_names_res[[b]]
        }
        if (add_class) {
          class(return_value[[b]]$res.int) <- c("lavaan.vector", "numeric")
        }
      }

      # thresholds
      if (lavmodel@categorical) {
        if (h1) {
          return_value[[b]]$res.th <- as.numeric(h1_implied$res.th[[b]])
        } else {
          return_value[[b]]$res.th <- as.numeric(lavsamplestats@res.th[[b]])
        }
        if (length(lavmodel@num.idx[[b]]) > 0L) {
          num_idx <- which(lavmodel@th.idx[[b]] == 0)
          return_value[[b]]$res.th <- return_value[[b]]$res.th[-num_idx]
        }
        # FIXME: if std: what to do?
        if (add_labels) {
          names(return_value[[b]]$res.th) <- object@pta$vnames$th[[b]]
        }
        if (add_class) {
          class(return_value[[b]]$res.th) <- c("lavaan.vector", "numeric")
        }
      }

      # slopes
      if (lavmodel@nexo[b] > 0L) {
        if (h1) {
          return_value[[b]]$res.slopes  <- h1_implied$res.slopes[[b]]
        } else {
          return_value[[b]]$res.slopes  <- lavsamplestats@res.slopes[[b]]
        }
        # FIXME: if std: what to do? (here: b.z = b * s.x /s.y)
        if (std) {
          tmp_y <- matrix(sqrt(diag_orig),
            nrow(return_value[[b]]$res.slopes),
            ncol(return_value[[b]]$res.slopes))
          tmp_x <- matrix(sqrt(diag(lavsamplestats@cov.x[[b]])),
            nrow(return_value[[b]]$res.slopes),
            ncol(return_value[[b]]$res.slopes), byrow = TRUE)
          return_value[[b]]$res.slopes <- return_value[[b]]$res.slopes /
                                          tmp_y * tmp_x
        }
        if (add_labels) {
          rownames(return_value[[b]]$res.slopes) <- ov_names_res[[b]]
          colnames(return_value[[b]]$res.slopes) <- ov_names_x[[b]]
        }
        if (add_class) {
          class(return_value[[b]]$res.slopes) <- c("lavaan.matrix", "matrix")
        }
      }

      # cov.x
      if (lavmodel@nexo[b] > 0L) {
        return_value[[b]]$cov.x  <- lavsamplestats@cov.x[[b]]
        if (std) {
          diag_orig <- diag(return_value[[b]]$cov.x)
          return_value[[b]]$cov.x <- cov2cor(return_value[[b]]$cov.x)
        }
        if (add_labels) {
          rownames(return_value[[b]]$cov.x) <- ov_names_x[[b]]
          colnames(return_value[[b]]$cov.x) <- ov_names_x[[b]]
        }
        if (add_class) {
          class(return_value[[b]]$cov.x) <-
            c("lavaan.matrix.symmetric", "matrix")
        }
      }

      # mean.x
      if (lavmodel@nexo[b] > 0L) {
        return_value[[b]]$mean.x <- as.numeric(object@SampleStats@mean.x[[b]])
        if (std) {
          diag_orig[diag_orig < .Machine$double.eps] <- NA
          return_value[[b]]$mean.x <- return_value[[b]]$mean.x / sqrt(diag_orig)
        }
        if (add_labels) {
          names(return_value[[b]]$mean.x) <- ov_names_x[[b]]
        }
        if (add_class) {
          class(return_value[[b]]$mean.x) <- c("lavaan.vector", "numeric")
        }
      }

    } # conditional.x

    # stochastic weights
    if (lavmodel@group.w.free) {
      # to be consistent with the 'implied' values,
      # transform so group.w is the 'log(group.freq)'
      return_value[[b]]$group.w <-
        log(lavsamplestats@group.w[[b]] * lavsamplestats@ntotal)
      if (add_labels) {
        names(return_value[[b]]$group.w) <- "w"
      }
      if (add_class) {
        class(return_value[[b]]$group.w) <- c("lavaan.vector", "numeric")
      }
    }
  }

  if (nblocks == 1L && drop_list_single_group) {
    return_value <- return_value[[1]]
  } else if (nblocks > 1L) {
    names(return_value) <- object@Data@block.label
  }

  return_value
}


lav_object_inspect_data <- function(object, add_labels = FALSE,
                                    drop_list_single_group = FALSE) {

  n_g <- object@Data@ngroups
  if (object@Model@conditional.x) {
    return_value <- vector("list", length = n_g)
    for (g in 1:n_g) {
      return_value[[g]] <- cbind(object@Data@X[[g]],
        object@Data@eXo[[g]])
    }
  } else {
    return_value <- object@Data@X
  }

  if (add_labels) {
    for (g in 1:n_g) {
      if (object@Model@conditional.x) {
        colnames(return_value[[g]]) <- c(object@Data@ov.names[[g]],
          object@Data@ov.names.x[[g]])
      } else {
        colnames(return_value[[g]]) <- object@Data@ov.names[[g]]
      }
    }
  }

  if (n_g == 1L && drop_list_single_group) {
    return_value <- return_value[[1]]
  } else {
    if (length(object@Data@group.label) > 0L) {
      names(return_value) <- unlist(object@Data@group.label)
    }
  }

  return_value
}

lav_object_inspect_case_idx <- function(object,
                                        drop_list_single_group = FALSE) {
  n_g <- object@Data@ngroups

  return_value <- object@Data@case.idx

  if (n_g == 1L && drop_list_single_group) {
    return_value <- return_value[[1]]
  } else {
    if (length(object@Data@group.label) > 0L) {
      names(return_value) <- unlist(object@Data@group.label)
    }
  }

  return_value
}

# lav_object_inspect_case_idx <- function(object, level = 1L,
#                                        drop_list_single_group = FALSE) {
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
    drop_list_single_group = FALSE) {

  n_g <- object@Data@ngroups
  nlevels <- object@Data@nlevels

  # just in case we have no clusters
  if (nlevels == 1L) {
    if (what %in% c("nclusters", "ncluster.size", "cluster.id")) {
      return_value <- as.list(rep(1L, n_g))
    } else if (what %in% c("cluster.size", "cluster.sizes")) {
      return_value <- object@Data@nobs
    } else if (what %in% c("cluster.idx", "cluster.label")) {
      # everybody belongs to cluster 1
      return_value <- lapply(seq_len(n_g),
        function(gg) rep(1L, object@Data@nobs[[gg]]))
    }
  }

  # if we do have clusters
  if (nlevels > 1L) {
    return_value <- vector("list", length = n_g)

    for (g in seq_len(n_g)) {
      tmp_lp <- object@Data@Lp[[g]]
      if (what == "nclusters") {
        return_value[[g]] <- tmp_lp$nclusters[[level]]
      } else if (what == "ncluster.size") {
        return_value[[g]] <- tmp_lp$ncluster.size[[level]]
      } else if (what == "cluster.size") {
        return_value[[g]] <- tmp_lp$cluster.size[[level]]
      } else if (what == "cluster.id") {
        return_value[[g]] <- tmp_lp$cluster.id[[level]]
      } else if (what == "cluster.idx") {
        return_value[[g]] <- tmp_lp$cluster.idx[[level]]
      } else if (what == "cluster.label") {
        return_value[[g]] <-
          tmp_lp$cluster.id[[level]][tmp_lp$cluster.idx[[level]]]
      } else if (what == "cluster.sizes") {
        return_value[[g]] <- tmp_lp$cluster.sizes[[level]]
      } else if (what == "average.cluster.size") {
        nn_g <- object@Data@nobs[[g]]
        cluster_size <- tmp_lp$cluster.size[[level]]
        nclusters <- tmp_lp$nclusters[[level]]
        return_value[[g]] <- (nn_g^2 - sum(cluster_size^2)) /
          (nn_g * (nclusters - 1L))
      }
    } # g
  } # nlevels > 1L

  if (n_g == 1L && drop_list_single_group) {
    return_value <- return_value[[1]]
  } else {
    if (length(object@Data@group.label) > 0L) {
      names(return_value) <- unlist(object@Data@group.label)
    }
  }

  return_value
}


# count the number of clusters, or obtain tmp.n within each cluster
# lav_object_inspect_ncluster <- function(object, sizes = FALSE, #level = 2L,
#                                        drop_list_single_group = FALSE) {
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

lav_object_inspect_rsquare <- function(object, est_std_all = NULL,
    add_labels = FALSE, add_class = FALSE, drop_list_single_group = FALSE) {

  nblocks <- object@Model@nblocks
  return_value <- vector("list", length = nblocks)

  if (is.null(est_std_all)) {
    est_std_all <- lav_standardize_all(object)
  }

  partable <- object@ParTable
  partable$rsquare <- 1.0 - est_std_all
  # no values > 1.0
  partable$rsquare[partable$rsquare > 1.0] <- as.numeric(NA)

  for (b in seq_len(nblocks)) {
    ind_names <- partable$rhs[which(partable$op == "=~" &
      partable$block == b)]
    eqs_y_names <- partable$lhs[which(partable$op == "~"  &
      partable$block == b)]
    y_names <- unique(c(ind_names, eqs_y_names))

    idx <- which(partable$op == "~~" & partable$lhs %in% y_names &
      partable$rhs == partable$lhs & partable$block == b)
    tmp <- partable$rsquare[idx]

    if (add_labels && length(tmp) > 0L) {
      names(tmp) <- partable$lhs[idx]
    }
    if (add_class) {
      class(tmp) <- c("lavaan.vector", "numeric")
    }

    return_value[[b]] <- tmp
  }

  if (nblocks == 1L && drop_list_single_group) {
    return_value <- return_value[[1]]
  } else if (nblocks > 1L) {
    names(return_value) <- object@Data@block.label
  }

  return_value
}

# model implied sample stats
lav_object_inspect_implied <- function(object,                    # nolint
    add_labels = FALSE, add_class = FALSE, drop_list_single_group = FALSE) {

  nblocks <- object@Model@nblocks
  ov_names <- object@pta$vnames$ov
  ov_names_res <- object@pta$vnames$ov.nox
  ov_names_x   <- object@pta$vnames$ov.x

  # slots
  lavimplied <- object@implied
  lavmodel   <- object@Model

  # if nlevels, always set conditional.x = FALSE
  if (object@Data@nlevels > 1L) {
    lavimplied <- lav_model_implied_cond2uncond(lavimplied)
    conditional_x <- FALSE # for now (0.6-12)
  } else {
    conditional_x <- lavmodel@conditional.x
  }

  return_value <- vector("list", length = nblocks)
  for (b in seq_len(nblocks)) {

    if (!conditional_x) {
      # covariance matrix
      return_value[[b]]$cov  <- lavimplied$cov[[b]]
      if (add_labels && !is.null(return_value[[b]]$cov)) {
        rownames(return_value[[b]]$cov) <- colnames(return_value[[b]]$cov) <-
          ov_names[[b]]
      }
      if (add_class) {
        class(return_value[[b]]$cov) <- c("lavaan.matrix.symmetric", "matrix")
      }

      # mean vector
      if (lavmodel@meanstructure) {
        return_value[[b]]$mean <- as.numeric(lavimplied$mean[[b]])
        if (add_labels) {
          names(return_value[[b]]$mean) <- ov_names[[b]]
        }
        if (add_class) {
          class(return_value[[b]]$mean) <- c("lavaan.vector", "numeric")
        }
      }

      # thresholds
      if (lavmodel@categorical) {
        return_value[[b]]$th <- as.numeric(lavimplied$th[[b]])
        if (length(object@Model@num.idx[[b]]) > 0L) {
          num_idx <- which(object@Model@th.idx[[b]] == 0)
          return_value[[b]]$th <- return_value[[b]]$th[-num_idx]
        }
        if (add_labels) {
          names(return_value[[b]]$th) <- object@pta$vnames$th[[b]]
        }
        if (add_class) {
          class(return_value[[b]]$th) <- c("lavaan.vector", "numeric")
        }
      }
      # !conditional.x
    } else {
      # if conditional.x = TRUE
      # residual covariance matrix
      return_value[[b]]$res.cov  <- lavimplied$res.cov[[b]]
      if (add_labels) {
        rownames(return_value[[b]]$res.cov) <-
          colnames(return_value[[b]]$res.cov) <-
          ov_names_res[[b]]
      }
      if (add_class) {
        class(return_value[[b]]$res.cov) <-
          c("lavaan.matrix.symmetric", "matrix")
      }

      # intercepts
      if (lavmodel@meanstructure) {
        return_value[[b]]$res.int <- as.numeric(lavimplied$res.int[[b]])
        if (add_labels) {
          names(return_value[[b]]$res.int) <- ov_names_res[[b]]
        }
        if (add_class) {
          class(return_value[[b]]$res.int) <- c("lavaan.vector", "numeric")
        }
      }

      # thresholds
      if (lavmodel@categorical) {
        return_value[[b]]$res.th <- as.numeric(lavimplied$res.th[[b]])
        if (length(object@Model@num.idx[[b]]) > 0L) {
          num_idx <- which(object@Model@th.idx[[b]] == 0)
          return_value[[b]]$res.th <- return_value[[b]]$res.th[-num_idx]
        }
        if (add_labels) {
          names(return_value[[b]]$res.th) <- object@pta$vnames$th[[b]]
        }
        if (add_class) {
          class(return_value[[b]]$res.th) <- c("lavaan.vector", "numeric")
        }
      }

      # slopes
      if (lavmodel@nexo[b] > 0L) {
        return_value[[b]]$res.slopes  <- lavimplied$res.slopes[[b]]
        if (add_labels) {
          rownames(return_value[[b]]$res.slopes) <- ov_names_res[[b]]
          colnames(return_value[[b]]$res.slopes) <- ov_names_x[[b]]
        }
        if (add_class) {
          class(return_value[[b]]$res.slopes) <- c("lavaan.matrix", "matrix")
        }
      }

      # cov.x
      if (lavmodel@nexo[b] > 0L) {
        return_value[[b]]$cov.x  <- lavimplied$cov.x[[b]]
        if (add_labels) {
          rownames(return_value[[b]]$cov.x) <- ov_names_x[[b]]
          colnames(return_value[[b]]$cov.x) <- ov_names_x[[b]]
        }
        if (add_class) {
          class(return_value[[b]]$cov.x) <-
            c("lavaan.matrix.symmetric", "matrix")
        }
      }

      # mean.x
      if (lavmodel@nexo[b] > 0L) {
        return_value[[b]]$mean.x  <- as.numeric(lavimplied$mean.x[[b]])
        if (add_labels) {
          names(return_value[[b]]$mean.x) <- ov_names_x[[b]]
        }
        if (add_class) {
          class(return_value[[b]]$mean.x) <- c("lavaan.vector", "numeric")
        }
      }


    } # conditional.x

    # stochastic weights
    if (lavmodel@group.w.free) {
      return_value[[b]]$group.w <- lavimplied$group.w[[b]]
      if (add_labels) {
        names(return_value[[b]]$group.w) <- "w" # somewhat redundant
      }
      if (add_class) {
        class(return_value[[b]]$group.w) <- c("lavaan.vector", "numeric")
      }
    }
  }

  if (nblocks == 1L && drop_list_single_group) {
    return_value <- return_value[[1]]
  } else if (nblocks > 1L) {
    names(return_value) <- object@Data@block.label
  }

  return_value
}


# residuals: _inspect_sampstat - _inspect_implied
lav_object_inspect_residuals <- function(object, h1 = TRUE,
    add_labels = FALSE, add_class = FALSE, drop_list_single_group = FALSE) {

  lav_residuals(object, type = "raw", h1 = h1,
    add.labels = add_labels, add.class = add_class,
    drop.list.single.group = drop_list_single_group)
}


lav_object_inspect_cov_lv <- function(object, correlation_metric = FALSE,
    add_labels = FALSE, add_class = FALSE, drop_list_single_group = FALSE) {
  # compute lv covar
  return_value <- lav_model_veta(lavmodel = object@Model,
                                 remove_dummy_lv = TRUE)

  # nblocks
  nblocks <- length(return_value)

  # cor + labels + class
  for (b in seq_len(nblocks)) {

    if (correlation_metric && nrow(return_value[[b]]) > 1L) {
      # note: cov2cor fails if matrix is empty!
      return_value[[b]] <- cov2cor(return_value[[b]])
    }

    if (add_labels) {
      colnames(return_value[[b]]) <- rownames(return_value[[b]]) <-
        object@pta$vnames$lv[[b]]
    }

    if (add_class) {
      class(return_value[[b]]) <- c("lavaan.matrix.symmetric", "matrix")
    }
  }

  if (nblocks == 1L && drop_list_single_group) {
    return_value <- return_value[[1]]
  } else if (nblocks > 1L) {
    names(return_value) <- object@Data@block.label
  }

  return_value
}

lav_object_inspect_mean_lv <- function(object,
    add_labels = FALSE, add_class = FALSE, drop_list_single_group = FALSE) {
  # compute lv means
  return_value <- lav_model_eeta(lavmodel = object@Model,
    lavsamplestats = object@SampleStats,
    remove_dummy_lv = TRUE)

  # nblocks
  nblocks <- length(return_value)

  # ensure numeric
  return_value <- lapply(return_value, as.numeric)

  # labels + class
  for (b in seq_len(nblocks)) {
    if (add_labels && length(return_value[[b]]) > 0L) {
      names(return_value[[b]]) <- object@pta$vnames$lv[[b]]
    }
    if (add_class) {
      class(return_value[[b]]) <- c("lavaan.vector", "numeric")
    }
  }

  if (nblocks == 1L && drop_list_single_group) {
    return_value <- return_value[[1]]
  } else if (nblocks > 1L) {
    names(return_value) <- object@Data@block.label
  }

  return_value
}

lav_object_inspect_cov_all <- function(object, correlation_metric = FALSE,
    add_labels = FALSE, add_class = FALSE, drop_list_single_group = FALSE) {
  # compute extended model implied covariance matrix (both ov and lv)
  return_value <- lav_model_cov_both(lavmodel = object@Model,
    remove_dummy_lv = TRUE)

  # nblocks
  nblocks <- length(return_value)

  # cor + labels + class
  for (b in seq_len(nblocks)) {

    if (correlation_metric && nrow(return_value[[b]]) > 1L) {
      # note: cov2cor fails if matrix is empty!
      return_value[[b]] <- cov2cor(return_value[[b]])
    }

    if (add_labels) {
      tmp_names <- c(object@pta$vnames$ov.model[[b]],
        object@pta$vnames$lv[[b]])
      colnames(return_value[[b]]) <- rownames(return_value[[b]]) <- tmp_names
    }
    if (add_class) {
      class(return_value[[b]]) <- c("lavaan.matrix.symmetric", "matrix")
    }
  }

  if (nblocks == 1L && drop_list_single_group) {
    return_value <- return_value[[1]]
  } else if (nblocks > 1L) {
    names(return_value) <- object@Data@block.label
  }

  return_value
}


lav_object_inspect_cov_ov <- function(object, correlation_metric = FALSE,
    add_labels = FALSE, add_class = FALSE, drop_list_single_group = FALSE) {
  # get model-implied covariance matrix observed
  if (object@Model@conditional.x) {
    return_value <- object@implied$res.cov
  } else {
    return_value <- object@implied$cov
  }

  # nblocks
  nblocks <- length(return_value)

  # cor + labels + class
  for (b in seq_len(nblocks)) {

    if (correlation_metric && nrow(return_value[[b]]) > 1L) {
      # note: cov2cor fails if matrix is empty!
      return_value[[b]] <- cov2cor(return_value[[b]])
    }

    if (add_labels) {
      colnames(return_value[[b]]) <- rownames(return_value[[b]]) <-
        object@pta$vnames$ov.model[[b]]
    }
    if (add_class) {
      class(return_value[[b]]) <- c("lavaan.matrix.symmetric", "matrix")
    }
  }

  if (nblocks == 1L && drop_list_single_group) {
    return_value <- return_value[[1]]
  } else {
    names(return_value) <- object@Data@block.label
  }

  return_value
}

lav_object_inspect_mean_ov <- function(object,
    add_labels = FALSE, add_class = FALSE, drop_list_single_group = FALSE) {
  # compute ov means
  if (object@Model@conditional.x) {
    return_value <- object@implied$res.int
  } else {
    return_value <- object@implied$mean
  }

  # nblocks
  nblocks <- length(return_value)

  # make numeric
  return_value <- lapply(return_value, as.numeric)

  # labels + class
  for (b in seq_len(nblocks)) {
    if (add_labels && length(return_value[[b]]) > 0L) {
      names(return_value[[b]]) <- object@pta$vnames$ov.model[[b]]
    }
    if (add_class) {
      class(return_value[[b]]) <- c("lavaan.vector", "numeric")
    }
  }

  if (nblocks == 1L && drop_list_single_group) {
    return_value <- return_value[[1]]
  } else if (nblocks > 1L) {
    names(return_value) <- object@Data@block.label
  }

  return_value
}

lav_object_inspect_th <- function(object,
    add_labels = FALSE, add_class = FALSE, drop_list_single_group = FALSE) {
  # thresholds
  if (object@Model@conditional.x) {
    return_value <- object@implied$res.th
  } else {
    return_value <- object@implied$th
  }

  # nblocks
  nblocks <- length(return_value)

  # make numeric
  return_value <- lapply(return_value, as.numeric)

  # labels + class
  for (b in seq_len(nblocks)) {
    if (length(object@Model@num.idx[[b]]) > 0L) {
      num_idx <- which(object@Model@th.idx[[b]] == 0)
      return_value[[b]] <- return_value[[b]][-num_idx]
    }
    if (add_labels && length(return_value[[b]]) > 0L) {
      names(return_value[[b]]) <- object@pta$vnames$th[[b]]
    }
    if (add_class) {
      class(return_value[[b]]) <- c("lavaan.vector", "numeric")
    }
  }

  if (nblocks == 1L && drop_list_single_group) {
    return_value <- return_value[[1]]
  } else if (nblocks > 1L) {
    names(return_value) <- object@Data@block.label
  }

  return_value
}

lav_object_inspect_th_idx <- function(object,
    add_labels = FALSE, add_class = FALSE, drop_list_single_group = FALSE) {
  # thresholds idx
  return_value <- object@SampleStats@th.idx

  # nblocks
  nblocks <- length(return_value)

  # labels + class
  for (b in seq_len(nblocks)) {
    if (add_labels && length(return_value[[b]]) > 0L) {
      names(return_value[[b]]) <- object@SampleStats@th.names[[b]]
    }
    if (add_class && !is.null(return_value[[b]])) {
      class(return_value[[b]]) <- c("lavaan.vector", "numeric")
    }
  }

  if (nblocks == 1L && drop_list_single_group) {
    return_value <- return_value[[1]]
  } else if (nblocks > 1L) {
    names(return_value) <- object@Data@block.label
  }

  return_value
}

lav_object_inspect_vy <- function(object,
    add_labels = FALSE, add_class = FALSE, drop_list_single_group = FALSE) {
  # 'unconditional' model-implied variances
  #  - same as diag(Sigma.hat) if all Y are continuous)
  #  - 1.0 (or delta^2) if categorical
  #  - if also Gamma, cov.x is used (only if categorical)

  return_value <- lav_model_vy(lavmodel = object@Model, glist = NULL,
    diagonal_only = TRUE)

  # nblocks
  nblocks <- length(return_value)

  # labels + class
  for (b in seq_len(nblocks)) {
    if (add_labels && length(return_value[[b]]) > 0L) {
      if (object@Model@categorical) {
        names(return_value[[b]]) <- object@pta$vnames$ov.nox[[b]]
      } else {
        names(return_value[[b]]) <- object@pta$vnames$ov[[b]]
      }
    }
    if (add_class) {
      class(return_value[[b]]) <- c("lavaan.vector", "numeric")
    }
  }

  if (nblocks == 1L && drop_list_single_group) {
    return_value <- return_value[[1]]
  } else if (nblocks > 1L) {
    names(return_value) <- object@Data@block.label
  }

  return_value
}

lav_object_inspect_fs_determinacy <- function(object,  # nolint
  squared = TRUE, fs_method = "regression",
  add_labels = FALSE, add_class = FALSE, drop_list_single_group = FALSE) {

  fs <- lavPredict(object, type = "lv", method = fs_method, rel = TRUE)
  return_value <- attr(fs, "rel")

  # determinacies or reliabilities?
  if (!squared) {
    return_value <- lapply(return_value, sqrt)
  }

  # nblocks
  nblocks <- length(return_value)

  # ensure numeric
  return_value <- lapply(return_value, as.numeric)

  # labels + class
  for (b in seq_len(nblocks)) {
    if (add_labels && length(return_value[[b]]) > 0L) {
      names(return_value[[b]]) <- object@pta$vnames$lv[[b]]
    }
    if (add_class) {
      class(return_value[[b]]) <- c("lavaan.vector", "numeric")
    }
  }

  if (nblocks == 1L && drop_list_single_group) {
    return_value <- return_value[[1]]
  } else if (nblocks > 1L) {
    names(return_value) <- object@Data@block.label
  }

  return_value
}


lav_object_inspect_theta <- function(object, correlation_metric = FALSE,
    add_labels = FALSE, add_class = FALSE, drop_list_single_group = FALSE) {
  # get residual covariances
  return_value <- lav_model_theta(lavmodel = object@Model)

  # nblocks
  nblocks <- length(return_value)

  # labels + class
  for (b in seq_len(nblocks)) {

    if (correlation_metric && nrow(return_value[[b]]) > 0L) {
      if (all(return_value[[b]] == 0)) {
        return_value[[b]] <- return_value[[b]]
      } else {
        return_value[[b]] <- cov2cor(return_value[[b]])
      }
    }

    if (add_labels && length(return_value[[b]]) > 0L) {
      colnames(return_value[[b]]) <- rownames(return_value[[b]]) <-
        object@pta$vnames$ov.model[[b]]
    }

    if (add_class) {
      class(return_value[[b]]) <- c("lavaan.matrix.symmetric", "matrix")
    }
  }

  if (nblocks == 1L && drop_list_single_group) {
    return_value <- return_value[[1]]
  } else if (nblocks > 1L) {
    names(return_value) <- object@Data@block.label
  }

  return_value
}


lav_object_inspect_missing_coverage <- function(object,          # nolint
    add_labels = FALSE, add_class = FALSE, drop_list_single_group = FALSE) {

  n_g <- object@Data@ngroups
  return_value <- vector("list", n_g)

  for (g in 1:n_g) {
    if (!is.null(object@Data@Mp[[g]])) {
      return_value[[g]] <- object@Data@Mp[[g]]$coverage
    } else {
      nvar <- length(object@Data@ov.names[[g]])
      return_value[[g]] <- matrix(1.0, nvar, nvar)
    }

    if (add_labels && length(return_value[[g]]) > 0L) {
      colnames(return_value[[g]]) <- rownames(return_value[[g]]) <-
        object@pta$vnames$ov.model[[g]]
    }

    if (add_class) {
      class(return_value[[g]]) <- c("lavaan.matrix.symmetric", "matrix")
    }
  }

  if (n_g == 1L && drop_list_single_group) {
    return_value <- return_value[[1]]
  } else {
    if (length(object@Data@group.label) > 0L) {
      names(return_value) <- unlist(object@Data@group.label)
    }
  }

  return_value
}

lav_object_inspect_missing_patterns <- function(object,           # nolint
    add_labels = FALSE, add_class = FALSE, drop_list_single_group = FALSE) {

  n_g <- object@Data@ngroups
  return_value <- vector("list", n_g)

  for (g in 1:n_g) {
    if (!is.null(object@Data@Mp[[g]])) {
      return_value[[g]] <- object@Data@Mp[[g]]$pat
    } else {
      nvar <- length(object@Data@ov.names[[g]])
      return_value[[g]] <- matrix(TRUE, 1L, nvar)
      rownames(return_value[[g]]) <- object@Data@nobs[[g]]
    }

    if (add_labels && length(return_value[[g]]) > 0L) {
      colnames(return_value[[g]]) <- object@pta$vnames$ov.model[[g]]
    }

    if (add_class) {
      class(return_value[[g]]) <- c("lavaan.matrix", "matrix")
    }
  }

  if (n_g == 1L && drop_list_single_group) {
    return_value <- return_value[[1]]
  } else {
    if (length(object@Data@group.label) > 0L) {
      names(return_value) <- unlist(object@Data@group.label)
    }
  }

  return_value
}

lav_object_inspect_empty_idx <- function(object,
                                         drop_list_single_group = FALSE) {

  n_g <- object@Data@ngroups

  # get empty idx
  return_value <- vector("list", n_g)

  for (g in 1:n_g) {
    if (!is.null(object@Data@Mp[[g]])) {
      return_value[[g]] <- object@Data@Mp[[g]]$empty.idx
    } else {
      return_value[[g]] <- integer(0L)
    }
  }

  if (n_g == 1L && drop_list_single_group) {
    return_value <- return_value[[1]]
  } else {
    if (length(object@Data@group.label) > 0L) {
      names(return_value) <- unlist(object@Data@group.label)
    }
  }

  return_value
}


lav_object_inspect_wls_est <- function(object,
    add_labels = FALSE, add_class = FALSE, drop_list_single_group = FALSE) {

  return_value <- lav_model_wls_est(object@Model)

  if (add_labels) {
    tmp_names <- lav_object_inspect_delta_rownames(object,
      drop_list_single_group = FALSE)
  }

  # nblocks
  nblocks <- length(return_value)

  for (b in seq_len(nblocks)) {
    if (add_labels && length(return_value[[b]]) > 0L &&
        object@Data@nlevels == 1L) {
      names(return_value[[b]]) <- tmp_names[[b]]
    }

    if (add_class) {
      class(return_value[[b]]) <- c("lavaan.vector", "numeric")
    }
  }

  if (nblocks == 1L && drop_list_single_group) {
    return_value <- return_value[[1]]
  } else if (nblocks > 1L) {
    names(return_value) <- object@Data@block.label
  }

  return_value
}

lav_object_inspect_wls_obs <- function(object,
    add_labels = FALSE, add_class = FALSE, drop_list_single_group = FALSE) {

  return_value <- object@SampleStats@WLS.obs ### FIXME: should be in @h1??

  if (add_labels) {
    tmp_names <- lav_object_inspect_delta_rownames(object,
      drop_list_single_group = FALSE)
  }

  # nblocks
  nblocks <- length(return_value)

  for (b in seq_len(nblocks)) {
    if (add_labels && length(return_value[[b]]) > 0L &&
        object@Data@nlevels == 1L) {
      names(return_value[[b]]) <- tmp_names[[b]]
    }

    if (add_class) {
      class(return_value[[b]]) <- c("lavaan.vector", "numeric")
    }
  }

  if (nblocks == 1L && drop_list_single_group) {
    return_value <- return_value[[1]]
  } else if (nblocks > 1L) {
    names(return_value) <- object@Data@block.label
  }

  return_value
}

lav_object_inspect_wls_v <- function(object,
    add_labels = FALSE, add_class = FALSE, drop_list_single_group = FALSE) {

  # return.value <- lav_model_wls_v(lavmodel       = object@Model,
  #                       lavsamplestats = object@SampleStats,
  #                       structured     = TRUE,
  #                       lavdata        = object@Data)
  # WLS.V == (traditionally) h1 expected information
  return_value <- lav_model_h1_information_expected(lavobject = object)
  # this affects fit measures gfi, agfi, pgfi


  # nblocks
  nblocks <- length(return_value)

  # if estimator == "DWLS" or "ULS", we only stored the diagonal
  # hence, we create a full matrix here
  if (object@Options$estimator %in% c("DWLS", "ULS")) {
    return_value <- lapply(return_value,
      function(x) {
        nr <- NROW(x)
        diag(x, nrow = nr, ncol = nr)
      })
  }

  if (add_labels) {
    tmp_names <- lav_object_inspect_delta_rownames(object,
      drop_list_single_group = FALSE)
  }

  # label + class
  for (b in seq_len(nblocks)) {
    # labels
    if (add_labels && nrow(return_value[[b]]) > 0L &&
        object@Data@nlevels == 1L) {
      colnames(return_value[[b]]) <-
        rownames(return_value[[b]]) <- tmp_names[[b]]
    }

    # class
    if (add_class) {
      class(return_value[[b]]) <- c("lavaan.matrix", "matrix")
    }
  }

  if (nblocks == 1L && drop_list_single_group) {
    return_value <- return_value[[1]]
  } else if (nblocks > 1L) {
    names(return_value) <- object@Data@block.label
  }

  return_value
}


lav_object_inspect_sampstat_gamma <- function(object,         # nolint
    add_labels = FALSE, add_class = FALSE, drop_list_single_group = FALSE) {

  if (!is.null(object@SampleStats@NACOV[[1]])) {
    return_value <- object@SampleStats@NACOV
  } else {
    return_value <- lav_object_gamma(object)
  }

  if (add_labels) {
    tmp_names <- lav_object_inspect_delta_rownames(object,
      drop_list_single_group = FALSE)
  }

  # nblocks
  nblocks <- length(return_value)

  if (nblocks == 1L && drop_list_single_group) {
    return_value <- return_value[[1]]

    # labels
    if (add_labels) {
      colnames(return_value) <- rownames(return_value) <- tmp_names[[1]]
    }

    # class
    if (add_class) {
      class(return_value) <- c("lavaan.matrix.symmetric", "matrix")
    }
  } else {
    if (object@Data@nlevels == 1L && length(object@Data@group.label) > 0L) {
      names(return_value) <- unlist(object@Data@group.label)

      # labels
      if (add_labels) {
        for (g in seq_len(object@Data@ngroups)) {
          colnames(return_value[[g]]) <-
            rownames(return_value[[g]]) <- tmp_names[[g]]
        }
      }

      # class
      if (add_class) {
        for (g in seq_len(object@Data@ngroups)) {
          class(return_value[[g]]) <- c("lavaan.matrix.symmetric", "matrix")
        }
      }
    } else if (object@Data@nlevels > 1L &&
      length(object@Data@group.label) == 0L) {
      names(return_value) <- object@Data@level.label
    }
  }

  return_value
}


lav_object_inspect_gradient <- function(object,
    add_labels = FALSE, add_class = FALSE, logl = FALSE,
    optim = FALSE) {

  lavmodel       <- object@Model
  lavdata        <- object@Data
  lavsamplestats <- object@SampleStats

  if (optim) {
    logl <- FALSE
  }

  if (lavsamplestats@missing.flag ||
    object@Options$estimator == "PML") {
    group_weight <- FALSE
  } else {
    group_weight <- TRUE
  }
  current_verbose <- lav_verbose()
  if (lav_verbose(FALSE)) on.exit(lav_verbose(current_verbose), TRUE)
  dx <- lav_model_gradient(
    lavmodel       = lavmodel,
    glist          = NULL,
    lavsamplestats = lavsamplestats,
    lavdata        = object@Data,
    lavcache       = object@Cache,
    type           = "free",
    group_weight   = group_weight)
  lav_verbose(current_verbose)
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
          tmp_n <- lavdata@Lp[[1]]$nclusters[[1]]
          nclusters <- lavdata@Lp[[1]]$nclusters[[2]]
          dx <- dx * (2 * tmp_n) / nclusters
        } else {
          group_values <- lav_partable_group_values(lavpartable)
          for (g in seq_len(lavdata@ngroups)) {
            tmp_n <- lavdata@Lp[[g]]$nclusters[[1]]
            nclusters <- lavdata@Lp[[g]]$nclusters[[2]]
            g_idx <-
              which((lavpartable$group ==
                group_values[g])[lavpartable$free > 0L])
            dx[g_idx] <- dx[g_idx] * (2 * tmp_n) / nclusters
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
  if (add_labels) {
    if (optim && lavmodel@eq.constraints) {
      tmp_names_all <- lav_partable_labels(object@ParTable, type = "free")
      tmp_seq <- seq_along(tmp_names_all)
      pack_seq <- as.numeric((tmp_seq - lavmodel@eq.constraints.k0) %*%
        +lavmodel@eq.constraints.K)
      ok_idx <- which(pack_seq %in% tmp_seq)
      tmp_names <- rep("(eq.con)", length(pack_seq))
      tmp_names[ok_idx] <- tmp_names_all[pack_seq[ok_idx]]
      names(dx) <- tmp_names
    } else {
      names(dx) <- lav_partable_labels(object@ParTable, type = "free")
    }
  }

  # class
  if (add_class) {
    class(dx) <- c("lavaan.vector", "numeric")
  }

  dx
}

lav_object_inspect_hessian <- function(object,
    add_labels = FALSE, add_class = FALSE) {

  return_value <- lav_model_hessian(
    lavmodel       = object@Model,
    lavsamplestats = object@SampleStats,
    lavdata        = object@Data,
    lavcache       = object@Cache,
    lavoptions     = object@Options,
    group_weight   = TRUE)

  # labels
  if (add_labels) {
    colnames(return_value) <- rownames(return_value) <-
      lav_partable_labels(object@ParTable, type = "free")
  }

  # class
  if (add_class) {
    class(return_value) <- c("lavaan.matrix.symmetric", "matrix")
  }

  return_value
}

lav_object_inspect_information <- function(object,
    information = "default", augmented = FALSE, inverted = FALSE,
    add_labels = FALSE, add_class = FALSE) {

  if (information != "default") {
    # override option
    object@Options$information <- information
  }

  return_value <- lav_model_information(lavmodel =  object@Model,
    lavsamplestats = object@SampleStats,
    lavdata        = object@Data,
    lavcache       = object@Cache,
    lavimplied     = object@implied,
    lavh1          = object@h1,
    lavoptions     = object@Options,
    extra          = FALSE,
    augmented      = augmented,
    inverted       = inverted)

  # labels
  if (add_labels) {
    tmp_names <- lav_partable_labels(object@ParTable, type = "free")
    if (augmented) {
      n_extra <- nrow(return_value) - length(tmp_names)
      if (n_extra > 0L) {
        tmp_names <- c(tmp_names, paste("aug", 1:n_extra, sep = ""))
      }
    }
    colnames(return_value) <- rownames(return_value) <- tmp_names
  }

  # class
  if (add_class) {
    class(return_value) <- c("lavaan.matrix.symmetric", "matrix")
  }

  return_value
}

lav_object_inspect_h1_information <- function(object,            # nolint
    information = "default", h1_information = "default", inverted = FALSE,
    add_labels = FALSE, add_class = FALSE, drop_list_single_group = FALSE) {

  if (information != "default") {
    # override option
    object@Options$information <- information
  }
  if (h1_information != "default") {
    # override option
    object@Options$h1.information <- h1_information
  }

  lavmodel <- object@Model
  lavdata  <- object@Data

  # list!
  return_value <- lav_model_h1_information(lavmodel = lavmodel,
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

  if (add_labels) {
    tmp_names <- lav_object_inspect_delta_rownames(object,
      drop_list_single_group = FALSE)
  }

  # labels/class per group
  for (g in seq_len(lavmodel@ngroups)) {
    # labels
    if (add_labels) {
      colnames(return_value[[g]]) <-
        rownames(return_value[[g]]) <- tmp_names[[g]]
    }

    # class
    if (add_class) {
      class(return_value[[g]]) <- c("lavaan.matrix.symmetric", "matrix")
    }
  }

  # drop list?
  if (lavmodel@ngroups == 1L && drop_list_single_group) {
    return_value <- return_value[[1]]
  } else if (!is.null(lavdata)) {
    if (length(lavdata@group.label) > 0L) {
      names(return_value) <- unlist(lavdata@group.label)
    }
  }

  return_value
}

lav_object_inspect_vcov <- function(object, standardized = FALSE,
    type = "std.all", free_only = TRUE,
    add_labels = FALSE, add_class = FALSE, remove_duplicated = FALSE) {

  lavmodel <- object@Model
  lavoptions <- object@Options

  # store partable with pta in object to use cache in called functions
  if (is.null(attr(object, "vnames"))) {
    object@ParTable <- lav_partable_set_cache(object@ParTable, object@pta)
  }

  if (object@optim$npar == 0) {
    return_value <- matrix(0, 0, 0)
  } else {
    # check if we already have it
    # tmp <- try(slot(object, "vcov"), silent = TRUE)
    # if( !inherits(tmp, "try-error") && !is.null(object@vcov$vcov)
    #   && !(rotation && standardized)) {
    if (!is.null(object@vcov$vcov) &&
        !(standardized && lavmodel@ceq.simple.only)) {
      return_value <- object@vcov$vcov 
                 # if this changes, tag @TDJorgensen in commit message
    } else {
      # compute it again
      # if(rotation && standardized) {
      #    lavmodel <- lav_model_set_parameters(lavmodel,
      #                    x = object@optim$x)
      #    lavoptions <- object@Options
      #    lavoptions$rotation.se <- "delta"
      # }

      # special case: ceq.simple.only + standardized
      if (standardized && lavmodel@ceq.simple.only) {
        # create 'compact' vcov (compatible with tmp.jac later)
        lavmodel@ceq.simple.only <- FALSE
        lavmodel@con.jac <- matrix(0, 0L, 0L)
      }

      return_value <- lav_model_vcov(
        lavmodel       = lavmodel,
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

      if (is.null(return_value)) {
        return(return_value)
      }
    }
  }

  # strip attributes
  attr(return_value, "E.inv") <- NULL
  attr(return_value, "B0") <- NULL
  attr(return_value, "B0.group") <- NULL
  attr(return_value, "Delta") <- NULL
  attr(return_value, "WLS.V") <- NULL
  attr(return_value, "tmp.boot.COEF") <- NULL
  attr(return_value, "tmp.boot.TEST") <- NULL

  # standardized?
  if (standardized) {
    if (type == "std.lv") {
      tmp_fun <- lav_standardize_lv_x
    } else if (type == "std.all") {
      tmp_fun <- lav_standardize_all_x
    } else if (type == "std.nox") {
      tmp_fun <- lav_standardize_all_nox_x
    }

    x_vec <- lav_model_get_parameters(lavmodel)
    tmp_jac <- try(lav_func_jacobian_complex(func = tmp_fun, x = x_vec,
      lavobject = object),
    silent = TRUE)
    if (inherits(tmp_jac, "try-error")) { # eg. pnorm()
      tmp_jac <- lav_func_jacobian_simple(func = tmp_fun, x = x_vec,
        lavobject = object)
    }
    # }

    # handle ceq.simple.only
    #if (object@Model@ceq.simple.only) {
    #  tmp.jac <- tmp.jac %*% t(object@Model@ceq.simple.K)
    #}

    # rows of tmp.jac contain *all* the parameters in the parameter table
    if (free_only) {
      if (object@Model@ceq.simple.only) {
        pt_free_idx <- which(object@ParTable$free > 0L &
          !duplicated(object@ParTable$free))
      } else {
        pt_free_idx <- which(object@ParTable$free > 0L)
      }
      tmp_jac <- tmp_jac[pt_free_idx, , drop = FALSE]
    }

    # apply Delta method
    return_value <- tmp_jac %*% return_value %*% t(tmp_jac)

    # force return.value to be symmetric and pd
    return_value <- (return_value + t(return_value)) / 2
    # return.value <- lav_matrix_symmetric_force_pd(return.value,
    #                                     tol = 1e-09) # was 1e-06 < 0.6-9
  }

  # labels
  if (add_labels) {
    # if(rotation && !free.only) {
    #    # todo
    # } else {
    colnames(return_value) <- rownames(return_value) <-
      lav_partable_labels(object@ParTable,
                          ## add "user" labels?
                          type = ifelse(standardized && !free_only,
                                        "user", "free"))
    # }
  }

  # alias?
  if (remove_duplicated && lavmodel@eq.constraints) {
    simple_flag <- lav_constraints_check_simple(lavmodel)
    if (simple_flag) {
      tmp_lab <- lav_partable_labels(object@ParTable, type = "free")
      dup_flag <- duplicated(tmp_lab)
      return_value <- return_value[!dup_flag, !dup_flag, drop = FALSE]
    } else {
      lav_msg_warn(
        gettext("alias is TRUE, but equality constraints do not appear
                to be simple; returning full vcov"))
    }
  }

  # class
  if (add_class) {
    class(return_value) <- c("lavaan.matrix.symmetric", "matrix")
  }

  return_value
}

lav_object_inspect_vcov_def <- function(object, joint = FALSE,
    standardized = FALSE, type = "std.all",
    add_labels = FALSE, add_class = FALSE) {

  lavmodel    <- object@Model
  lavpartable <- object@ParTable
  free_idx <- which(lavpartable$free > 0L)
  def_idx <- which(lavpartable$op == ":=")
  joint_idx <- c(free_idx, def_idx)

  if (!joint && length(def_idx) == 0L) {
    return(matrix(0, 0, 0))
  } else if (joint && length(joint_idx) == 0L) {
    return(matrix(0, 0, 0))
  }

  if (standardized) {
    # compute tmp.vcov for "free" parameters only
    tmp_vcov <- lav_object_inspect_vcov(object,
      standardized = TRUE,
      type = type, free_only = FALSE,
      add_labels = FALSE, add_class = FALSE)
    if (joint) {
      return_value <- tmp_vcov[joint_idx, joint_idx, drop = FALSE]
    } else {
      return_value <- tmp_vcov[def_idx, def_idx, drop = FALSE]
    }
  } else {
    # get free parameters
    x <- lav_model_get_parameters(lavmodel, type = "free")

    # bootstrap or not?
    if (!is.null(object@boot$coef)) {
      tmp_boot <- object@boot$coef
      # remove NA rows
      error_idx <- attr(tmp_boot, "error.idx")
      if (length(error_idx) > 0L) {
        tmp_boot <- tmp_boot[-error_idx, , drop = FALSE] # drops attributes
      }
      tmp_boot_def <- apply(tmp_boot, 1L, lavmodel@def.function)
      if (length(def_idx) == 1L) {
        tmp_boot_def <- as.matrix(tmp_boot_def)
      } else {
        tmp_boot_def <- t(tmp_boot_def)
      }
      if (joint) {
        tmp_boot_def <- cbind(tmp_boot, tmp_boot_def)
      }
      return_value <- cov(tmp_boot_def)
    } else {
      # tmp.vcov
      tmp_vcov <- lav_object_inspect_vcov(object,
        standardized = FALSE,
        type = type, free_only = TRUE,
        add_labels = FALSE,
        add_class = FALSE)

      # regular delta method
      tmp_jac <- try(lav_func_jacobian_complex(func = lavmodel@def.function,
        x = x), silent = TRUE)
      if (inherits(tmp_jac, "try-error")) { # eg. pnorm()
        tmp_jac <- lav_func_jacobian_simple(func = lavmodel@def.function,
          x = x)
      }
      if (joint) {
        tmp_jac2 <- rbind(diag(nrow = ncol(tmp_jac)), tmp_jac)
        return_value <- tmp_jac2 %*% tmp_vcov %*% t(tmp_jac2)
      } else {
        return_value <- tmp_jac %*% tmp_vcov %*% t(tmp_jac)
      }
    }
  }

  # labels
  if (add_labels) {
    if (joint) {
      lhs_names <- lavpartable$lhs[joint_idx]
    } else {
      lhs_names <- lavpartable$lhs[def_idx]
    }
    colnames(return_value) <- rownames(return_value) <- lhs_names
  }

  # class
  if (add_class) {
    class(return_value) <- c("lavaan.matrix.symmetric", "matrix")
  }

  return_value
}

lav_object_inspect_ugamma <- function(object,                  # nolint
    add_labels = FALSE, add_class = FALSE) {

  out <- lav_test_satorra_bentler(lavobject = object,
    method = "original",
    return_ugamma = TRUE)
  return_value <- out$UGamma

  # labels
  # if (add_labels) {
  # colnames(return.value) <- rownames(return.value) <-
  # }

  # class
  if (add_class) {
    class(return_value) <- c("lavaan.matrix.symmetric", "matrix")
  }

  return_value
}

lav_object_inspect_u_from_ugamma <- function(object,            # nolint
    add_labels = FALSE, add_class = FALSE, drop_list_single_group = FALSE) {

  out <- lav_test_satorra_bentler(lavobject     = object,
    method        = "original",
    return_u      = TRUE)
  return_value <- out$UfromUGamma

  # labels
  # if (add_labels) {
  # colnames(return.value) <- rownames(return.value) <-
  # }

  # class
  if (add_class) {
    class(return_value) <- c("lavaan.matrix.symmetric", "matrix")
  }

  return_value
}


# Delta (jacobian: d samplestats / d free_parameters)
lav_object_inspect_delta <- function(object,
    add_labels = FALSE, add_class = FALSE, drop_list_single_group = FALSE) {

  lavmodel    <- object@Model
  lavdata     <- object@Data
  lavpartable <- object@ParTable
  # lavpta      <- object@pta

  return_value <- lav_object_inspect_delta_internal(lavmodel = lavmodel,
    lavdata = lavdata, lavpartable = lavpartable,
    add_labels = add_labels, add_class = add_class,
    drop_list_single_group = drop_list_single_group)

  return_value
}

lav_object_inspect_delta_rownames <- function(                # nolint
    object,
    lavmodel = NULL,
    lavpartable = NULL,
    drop_list_single_group = FALSE) {
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
  correlation    <- lavmodel@correlation
  conditional_x  <- lavmodel@conditional.x
  group_w_free   <- lavmodel@group.w.free
  nvar           <- lavmodel@nvar
  num_idx        <- lavmodel@num.idx
  th_idx         <- lavmodel@th.idx
  nblocks        <- lavmodel@nblocks

  # store names per block, rbind later
  tmp_names <- vector("list", length = nblocks)

  # output is per group
  return_value <- vector("list", lavmodel@ngroups)

  for (g in 1:nblocks) {

    if (conditional_x) {
      ov_names <- lavpta$vnames$ov.nox[[g]]
    } else {
      ov_names <- lavpta$vnames$ov[[g]]
    }
    ov_names_x <- lavpta$vnames$ov.x[[g]]
    nvar <- length(ov_names)


    names_cov <- names_cor <- names_var <- character(0L)
    names_mu <- names_pi <- names_th <- character(0L)
    names_gw <- character(0L)

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
    tmp <- matrix(apply(expand.grid(ov_names, ov_names), 1L,
      paste, collapse = "~~"), nrow = nvar)
    if (categorical) {
      names_cor <- lav_matrix_vechru(tmp, diagonal = FALSE)
      names_var <- diag(tmp)[num_idx[[g]]]
    } else if (correlation) {
      names_cor <- lav_matrix_vechru(tmp, diagonal = FALSE)
    } else {
      names_cov <- lav_matrix_vechru(tmp, diagonal = TRUE)
    }


    # Mu
    if (!categorical && lavmodel@meanstructure) {
      names_mu <- paste(ov_names, "~1", sep = "")
    }

    # Pi
    if (conditional_x && lavmodel@nexo[g] > 0L) {
      names_pi <- apply(expand.grid(ov_names, ov_names_x), 1L,
        paste, collapse = "~")
    }

    # th
    if (categorical) {
      names_th <- lavpta$vnames$th[[g]]
      # interweave numeric intercepts, if any
      if (length(num_idx[[g]]) > 0L) {
        tmp <- character(length(th_idx[[g]]))
        tmp[th_idx[[g]] > 0] <- names_th
        tmp[th_idx[[g]] == 0] <- paste(ov_names[num_idx[[g]]],
          "~1", sep = "")
        names_th <- tmp
      }
    }

    # gw
    if (group_w_free) {
      names_gw <- "w"
    }

    tmp_names[[g]] <- c(names_gw,
      names_th, names_mu,
      names_pi,
      names_cov, names_var, names_cor)

  } # blocks

  # multilevel?
  if (lavmodel@multilevel) {
    for (g in 1:lavmodel@ngroups) {
      return_value[[g]] <- c(tmp_names[[(g - 1) * 2 + 1]],
        tmp_names[[(g - 1) * 2 + 2]])
    }
  } else {
    return_value <- tmp_names
  }

  # drop list?
  if (lavmodel@ngroups == 1L && drop_list_single_group) {
    return_value <- return_value[[1]]
  } else if (!is.null(lavdata)) {
    if (length(lavdata@group.label) > 0L) {
      names(return_value) <- unlist(lavdata@group.label)
    }
  }

  return_value
}

lav_object_inspect_delta_internal <- function(                   # nolint
    lavmodel = NULL, lavdata  = NULL,
    lavpartable = NULL,
    add_labels = FALSE, add_class = FALSE, drop_list_single_group = FALSE) {

  return_value <- lav_model_delta(lavmodel)

  if (add_labels) {
    tmp_pnames <- lav_partable_labels(lavpartable, type = "free")
    tmp_rownames  <- lav_object_inspect_delta_rownames(object = NULL,
      lavmodel = lavmodel, lavpartable = lavpartable,
      drop_list_single_group = FALSE)
  }

  for (g in seq_len(lavmodel@ngroups)) {
    # add labels
    if (add_labels) {
      colnames(return_value[[g]]) <- tmp_pnames
      rownames(return_value[[g]]) <- tmp_rownames[[g]]
    }

    # add class
    if (add_class) {
      class(return_value[[g]]) <- c("lavaan.matrix", "matrix")
    }

  } # ngroups

  # drop list?
  if (lavmodel@ngroups == 1L && drop_list_single_group) {
    return_value <- return_value[[1]]
  } else {
    if (length(lavdata@group.label) > 0L) {
      names(return_value) <- unlist(lavdata@group.label)
    }
  }

  return_value
}

lav_object_inspect_zero_cell_tables <-                          # nolint
  function(object, add_labels = FALSE, add_class = FALSE,
            drop_list_single_group = FALSE) {
  # categorical?
  if (!object@Model@categorical) {
    lav_msg_warn(gettext("no categorical variables in fitted model"))
    return(invisible(list()))
  }

  lavdata <- object@Data

  # create 2-way tables
  tmp_table <- lavTables(object, dimension = 2L, output = "data.frame",
    statistic = NULL)

  # select tables with empty cells
  empty_id <- tmp_table$id[which(tmp_table$obs.freq == 0)]


  if (length(empty_id) == 0L) {
    # only when lavInspect() is used, give message
    if (add_class) {
      cat("(There are no tables with empty cells for this fitted model)\n")
    }
    return(invisible(list()))
  }
  lav_tables_cells_format(
      tmp_table[tmp_table$id %in% empty_id, ],
      lavdata = lavdata,
      drop.list.single.group = drop_list_single_group)
}

lav_object_inspect_coef <- function(object, type = "free",
                                    add_labels = FALSE, add_class = FALSE) {

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
  tmp_est <- lav_object_inspect_est(object)
  cof <- tmp_est[idx]

  # labels?
  if (add_labels) {
    names(cof) <- lav_partable_labels(object@ParTable, type = type)
  }

  # class
  if (add_class) {
    class(cof) <- c("lavaan.vector", "numeric")
  }

  cof
}

lav_object_inspect_npar <- function(object, ceq = FALSE) {

  # free parameters (going to the optimizer)
  npar <- object@Model@nx.free

  # account for equality constraints?
  if (ceq && nrow(object@Model@con.jac) > 0L) {
    ceq_idx <- attr(object@Model@con.jac, "ceq.idx")
    if (length(ceq_idx) > 0L) {
      neq <- qr(object@Model@con.jac[ceq_idx, , drop = FALSE])$rank
      npar <- npar - neq
    }
  }

  npar
}

lav_object_inspect_icc <- function(object, add_labels = FALSE,
                                   add_class = FALSE,
                                   drop_list_single_group = FALSE) {

  lavdata <- object@Data
  n_g <- lavdata@ngroups
  return_value <- vector("list", n_g)

  # clustered data?
  if (length(lavdata@cluster) == 0L) {
    lav_msg_stop(gettext(
      "intraclass correlation only available for clustered data"))
  } else if (lavdata@nlevels == 1L) {
    lav_msg_stop(gettext(
      "intraclass correlation only available if the model syntax
       contains levels"))
  }

  if (length(object@h1) == 0L) {
    lav_msg_stop(gettext("h1 slot is not available; refit with h1 = TRUE"))
  }

  # implied statistics
  implied <- object@h1$implied

  for (g in 1:n_g) {
    sigma_w <- implied$cov[[(g - 1) * lavdata@nlevels + 1]]
    sigma_b <- implied$cov[[(g - 1) * lavdata@nlevels + 2]]

    w_diag <- diag(sigma_w)
    b_diag <- diag(sigma_b)

    return_value[[g]] <- numeric(length(w_diag))

    ov_names_l <- lavdata@ov.names.l[[g]]
    w_idx <- which(ov_names_l[[1]] %in% ov_names_l[[2]])
    w_names <- ov_names_l[[1]][w_idx]
    b_idx <- match(w_names, ov_names_l[[2]])

    return_value[[g]][w_idx] <- b_diag[b_idx] / (w_diag[w_idx] + b_diag[b_idx])

    # label
    if (add_labels) {
      names(return_value[[g]]) <- ov_names_l[[1]]
    }

    # class
    if (add_class) {
      class(return_value[[g]]) <- c("lavaan.vector", "numeric")
    }
  } # g

  if (n_g == 1L && drop_list_single_group) {
    return_value <- return_value[[1]]
  } else {
    if (length(object@Data@group.label) > 0L) {
      names(return_value) <- unlist(object@Data@group.label)
    }
  }

  return_value
}

lav_object_inspect_ranef <- function(object, add_labels = FALSE,
                                     add_class = FALSE,
                                     drop_list_single_group = FALSE) {

  lavdata <- object@Data
  lavsamplestats <- object@SampleStats

  n_g <- lavdata@ngroups
  return_value <- vector("list", n_g)

  # multilevel?
  if (lavdata@nlevels == 1L) {
    lav_msg_stop(gettext(
      "random effects only available for clustered data (in the long format)"))
  }

  # implied statistics
  lavimplied <- object@implied

  for (g in 1:n_g) {

    tmp_lp <- lavdata@Lp[[g]]
    tmp_ylp <- lavsamplestats@YLp[[g]]

    # implied for this group
    group_idx <- (g - 1) * lavdata@nlevels + seq_len(lavdata@nlevels)
    implied_group <- lapply(lavimplied, function(x) x[group_idx])

    # random effects (=random intercepts or cluster means)
    out <- lav_mvnorm_cluster_implied22l(lp = tmp_lp, implied = implied_group)
    mb_j <- lav_mvnorm_cluster_em_estep_ranef(ylp = tmp_ylp, lp = tmp_lp,
      sigma_w = out$sigma.w, sigma_b = out$sigma.b,
      sigma_zz = out$sigma.zz, sigma_yz = out$sigma.yz,
      mu_z = out$mu.z, mu_w = out$mu.w, mu_b = out$mu.b,
      se = FALSE)
    return_value[[g]] <- mb_j

    ov_names_l <- lavdata@ov.names.l[[g]]

    # label
    if (add_labels) {
      colnames(return_value[[g]]) <- ov_names_l[[1]]
    }

    # class
    if (add_class) {
      class(return_value[[g]]) <- c("lavaan.matrix", "matrix")

    }
  } # g

  if (n_g == 1L && drop_list_single_group) {
    return_value <- return_value[[1]]
  } else {
    if (length(object@Data@group.label) > 0L) {
      names(return_value) <- unlist(object@Data@group.label)
    }
  }

  return_value
}

# casewise loglikelihood contributions
lav_object_inspect_loglik_casewise <- function(object, log_1 = TRUE, # nolint
    add_labels = FALSE, add_class = FALSE, drop_list_single_group = FALSE) {

  lavdata <- object@Data
  lavsamplestats <- object@SampleStats
  lavimplied <- object@implied
  lavoptions <- object@Options

  n_g <- lavdata@ngroups
  return_value <- vector("list", n_g)

  # multilevel?
  if (lavdata@nlevels > 1L) {
    lav_msg_stop(gettext("casewise (log)likelihoods contributions
     not yet available for clustered data"))
  }

  # estimator ML?
  if (object@Options$estimator != "ML") {
    lav_msg_stop(gettextf("casewise (log)likelihoods contributions
     only available for estimator = %s",
    dQuote("ML")))
  }

  for (g in 1:n_g) {

    if (lavsamplestats@missing.flag) {
      return_value[[g]] <-
        lav_mvnorm_missing_llik_casewise(y = lavdata@X[[g]],
          wt = lavdata@weights[[g]],
          mu = lavimplied$mean[[g]],
          sigma_1 = lavimplied$cov[[g]],
          x_idx  = lavsamplestats@x.idx[[g]])
    } else { # single-level, complete data
      if (lavoptions$conditional.x) {
        if (!is.null(lavdata@weights[[g]])) {
          lav_msg_stop(gettext("no support (yet) if weights are used."))
        }
        return_value[[g]] <- lav_mvreg_loglik_data(
          y          = lavdata@X[[g]],
          exo        = lavdata@eXo[[g]],
          res_int    = lavimplied$res.int[[g]],
          res_slopes = lavimplied$res.slopes[[g]],
          res_cov    = lavimplied$res.cov[[g]],
          casewise   = TRUE)

      } else {

        if (object@Model@meanstructure) {
          tmp_mean <- lavimplied$mean[[g]]
        } else {
          tmp_mean <- lavsamplestats@mean[[g]]
        }
        return_value[[g]] <-
          lav_mvnorm_loglik_data(y = lavdata@X[[g]],
            wt = lavdata@weights[[g]],
            mu = tmp_mean,
            sigma_1 = lavimplied$cov[[g]],
            x_idx = lavsamplestats@x.idx[[g]],
            casewise = TRUE)
      }
    } # single-level, complete data

    # log. = FALSE?
    if (!log_1) {
      return_value[[g]] <- exp(return_value[[g]])
    }

    # label
    # if(add.labels) {
    # }

    # class
    if (add_class) {
      class(return_value[[g]]) <- c("lavaan.vector", "numeric")

    }

  } # g

  if (n_g == 1L && drop_list_single_group) {
    return_value <- return_value[[1]]
  } else {
    if (length(object@Data@group.label) > 0L) {
      names(return_value) <- unlist(object@Data@group.label)
    }
  }

  return_value
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
                                      add_labels = FALSE, add_class = FALSE,
                                      drop_list_single_group = FALSE) {

  lavdata <- object@Data
  n_g <- lavdata@ngroups

  # lavPredict()
  out <- lavPredict(object, type = type, method = "ML", # = Bartlett
    label = FALSE, fsm = TRUE, mdist = TRUE,
    se = "none", acov = "none")
  return_value <- attr(out, "mdist")

  for (g in seq_len(n_g)) {
    # squared?
    if (!squared) {
      return_value[[g]] <- sqrt(return_value[[g]])
    }

    # labels?
    # if(add.labels) {
    # }

    # class
    if (add_class) {
      class(return_value[[g]]) <- c("lavaan.vector", "numeric")
    }
  } # g

  if (n_g == 1L && drop_list_single_group) {
    return_value <- return_value[[1]]
  } else {
    if (length(object@Data@group.label) > 0L) {
      names(return_value) <- unlist(object@Data@group.label)
    }
  }

  return_value
}

# N versus N-1 (or N versus N-G in the multiple group setting)
# Changed 0.5-15: suggestion by Mark Seeto
lav_object_inspect_ntotal <- function(object) {
  if (object@Options$estimator %in% c("ML", "PML", "FML", "catML") &&
    object@Options$likelihood %in% c("default", "normal")) {
    n <- object@SampleStats@ntotal
  } else {
    n <- object@SampleStats@ntotal - object@SampleStats@ngroups
  }

  n
}

lav_object_inspect_iv <- function(object, drop_list_single_group = FALSE) {

  if (is.null(object@internal$eqs)) {
    lav_msg_stop(gettext("no equations/ivs found"))
  }
  lavmodel <- object@Model
  lavdata <- object@Data

  # grab equations
  iv_list <- object@internal$eqs

  # nblocks
  nblocks <- object@pta$nblocks

  table <- vector("list", length = nblocks)
  for (b in seq_len(nblocks)) {
    eqs <- iv_list[[b]]
    lhs <- sapply(eqs, "[[", "lhs")
    rhs <- sapply(lapply(eqs, "[[", "rhs"), paste, collapse = " + ")
    lhs_new <- sapply(eqs, "[[", "lhs_new")
    rhs_new <- sapply(lapply(eqs, "[[", "rhs_new"), paste, collapse = " + ")
    iv <- sapply(lapply(eqs, "[[", "iv"), paste, collapse = ", ")
    type <- sapply(eqs, "[[", "iv_type")
    table[[b]] <- data.frame(
      lhs = lhs, rhs = rhs,
      lhs.new = lhs_new, rhs.new = rhs_new, type = type, instruments = iv
    )
    class(table[[b]]) <- c("lavaan.data.frame", "data.frame")
  }

  # return value
  return_value <- table

  # drop list?
  if (lavmodel@ngroups == 1L && drop_list_single_group) {
    return_value <- return_value[[1]]
  } else if (!is.null(lavdata)) {
    if (length(lavdata@group.label) > 0L) {
      names(return_value) <- unlist(lavdata@group.label)
    }
  }

  return_value
}

lav_object_inspect_eqs <- function(object, drop_list_single_group = FALSE) {

  if (is.null(object@internal$eqs)) {
    lav_msg_stop(gettext("no equations/ivs found"))
  }
  lavmodel <- object@Model
  lavdata <- object@Data

  # grab equations
  eqs <- object@internal$eqs

  # return value
  return_value <- eqs

  # drop list?
  if (lavmodel@ngroups == 1L && drop_list_single_group) {
    return_value <- return_value[[1]]
  } else if (!is.null(lavdata)) {
    if (length(lavdata@group.label) > 0L) {
      names(return_value) <- unlist(lavdata@group.label)
    }
  }

  return_value
}

lav_object_inspect_sargan <- function(object, drop_list_single_group = FALSE) {

  if (is.null(object@internal$eqs)) {
    lav_msg_stop(gettext("no equations/ivs found"))
  }
  lavmodel <- object@Model
  lavdata <- object@Data

  # grab equations
  iv_list <- object@internal$eqs

  # nblocks
  nblocks <- object@pta$nblocks

  table <- vector("list", length = nblocks)
  for (b in seq_len(nblocks)) {
    eqs <- iv_list[[b]]
    lhs <- sapply(eqs, "[[", "lhs")
    rhs <- sapply(lapply(eqs, "[[", "rhs"), paste, collapse = " + ")
    miiv <- sapply(lapply(eqs, "[[", "miiv"), paste, collapse = ", ")
    sargan_stat <- sapply(seq_along(eqs),
      function(x) eqs[[x]][["sargan"]]["stat"])
    sargan_df <- sapply(seq_along(eqs),
      function(x) eqs[[x]][["sargan"]]["df"])
    sargan_pvalue <- sapply(seq_along(eqs),
      function(x) eqs[[x]][["sargan"]]["pvalue"])
    table[[b]] <- data.frame(
      lhs = lhs, rhs = rhs, instruments = miiv,
      sargan.stat = sargan_stat, df = sargan_df, pvalue = sargan_pvalue
    )

    # remove rows for which the Sargan statistic is NA
    na_idx <- which(is.na(sargan_stat))
    if (length(na_idx) > 0L) {
      table[[b]] <- table[[b]][-na_idx, , drop = FALSE]
    }

    class(table[[b]]) <- c("lavaan.data.frame", "data.frame")
  }

  # return value
  return_value <- table

  # drop list?
  if (lavmodel@ngroups == 1L && drop_list_single_group) {
    return_value <- return_value[[1]]
  } else if (!is.null(lavdata)) {
    if (length(lavdata@group.label) > 0L) {
      names(return_value) <- unlist(lavdata@group.label)
    }
  }

  return_value
}
