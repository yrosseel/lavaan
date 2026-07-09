# lav_gamma_recipe() -- single source of truth for "which Gamma (= NACOV,
# N times the asymptotic variance matrix of the sample statistics) does this
# analysis use, and with which settings?"
#
# Both the fit-time orchestrator (lav_samp_from_data) and the post-fit
# accessor (lav_object_gamma) consult this recipe, so that a Gamma that is
# recomputed after fitting cannot diverge from the one used at fit time.
#
# The recipe describes the DEFAULT policy implied by lavoptions + lavdata;
# caller-level overrides (a user-supplied NACOV= argument, or nacov =
# TRUE/FALSE) are handled by the callers themselves.

# YR 07 July 2026

lav_gamma_recipe <- function(lavoptions = NULL, lavdata = NULL,
                             nacov_compute = NULL) {
  stopifnot(!is.null(lavoptions), !is.null(lavdata))

  missing <- lavoptions$missing
  estimator <- lavoptions$estimator
  se <- lavoptions$se
  test <- lavoptions$test
  correlation <- isTRUE(lavoptions$correlation) # absent in old objects

  # data properties
  categorical <- any(lavdata@ov$type == "ordered")
  multilevel <- lavdata@nlevels > 1L
  clustered <- !multilevel && length(lavdata@cluster) > 0L
  wt <- length(lavdata@weights) > 0L && !is.null(lavdata@weights[[1]])

  # design vs frequency weighting of the (categorical/continuous) Gamma
  swt_type <- if (!is.null(lavoptions$sampling.weights.type)) {
    lavoptions$sampling.weights.type
  } else {
    "design"
  }

  # partial correlation structure: only these observed variables are scaled
  # to unit variance (empty -> all of them, the usual correlation structure)
  correlation_ov <- lavoptions$.correlation.ov
  if (is.null(correlation_ov)) {
    correlation_ov <- character(0L)
  }

  # 1. should NACOV be computed at fit time? (default policy)
  if (is.null(nacov_compute)) {
    nacov_compute <- FALSE
    if (se %in% c("robust.sem", "robust.sem.nt", "robust.cluster.sem") &&
        missing == "listwise") {
      nacov_compute <- TRUE
    }
    # note: test can be a vector...
    if (missing == "listwise" && any(test %in% c(
      "satorra.bentler",
      "mean.var.adjusted",
      "scaled.shifted"
    ))) {
      nacov_compute <- TRUE
    }
    if (missing == "listwise" &&
        any(vapply(test, lav_test_fmg_is_fmg, logical(1L)))) {
      nacov_compute <- TRUE
    }
    if (estimator == "IV" &&
        lavoptions$estimator.args$iv_vcov_stage1 == "gamma") {
      nacov_compute <- TRUE
    }
    if (any(missing == c("two.stage", "robust.two.stage")) &&
        estimator %in% c("ULS", "GLS", "WLS", "DLS")) {
      nacov_compute <- TRUE
    }
  }

  # 2. which flavor? (mirrors the branch order of the NACOV block in
  #    lav_samp_from_data, and the two-level builders in xxx_lavaan)
  #    - "adf"          ADF Gamma (lav_samp_gamma)
  #    - "nt"           normal-theory Gamma (lav_samp_gamma_nt); note: only
  #                     reachable for the (D)WLS-family estimators -- for
  #                     estimator = "ML", se = "robust.sem.nt" still stores
  #                     the ADF Gamma
  #    - "adf.cor"      (partial-)correlation ADF Gamma via the delta method
  #    - "missing.two.stage" / "missing.robust.two.stage"
  #                     NACOV of the saturated EM moments
  #    - "cat"          categorical: WLS.W (muthen1984), rescaled by N
  #    - "twolevel" / "twolevel.cat"
  #                     two-level (D)WLS builders (cluster sandwich)
  #    - "none"         no NACOV computed at fit time
  wls_family <- estimator %in% c("WLS", "DWLS", "ULS", "DLS", "IV", "catML")
  flavor <- "none"
  if (multilevel) {
    if (estimator %in% c("WLS", "DWLS", "ULS")) {
      flavor <- if (categorical) "twolevel.cat" else "twolevel"
    }
  } else if (nacov_compute && !categorical &&
             any(missing == c("two.stage", "robust.two.stage")) &&
             estimator %in% c("ULS", "GLS", "WLS", "DLS")) {
    flavor <- paste0("missing.", missing)
  } else if (estimator %in% c("ML", "GLS") && missing == "listwise" &&
             nacov_compute) {
    flavor <- if (correlation) "adf.cor" else "adf"
  } else if (wls_family) {
    if (categorical) {
      flavor <- "cat"
    } else if (correlation) {
      flavor <- "adf.cor"
    } else if (estimator %in% c("ULS", "DWLS") &&
               "robust.sem.nt" %in% se) {
      # only for ULS/DWLS (their continuous default se): for WLS/DLS
      # the (full) weight matrix is derived from the NACOV, which must
      # remain ADF
      flavor <- "nt"
    } else {
      flavor <- "adf"
    }
  }

  # 3. the knobs, as passed to the constructors
  #    note the asymmetry (as at fit time): gamma.wls.mplus only applies to
  #    the plain ADF Gamma of the (D)WLS-family estimators; the ML/GLS branch
  #    always uses mplus_wls = FALSE
  mplus_wls <- FALSE
  if (flavor == "adf" && wls_family && isTRUE(lavoptions$gamma.wls.mplus)) {
    mplus_wls <- TRUE
  }

  list(
    # decisions
    compute = nacov_compute,
    flavor = flavor,
    # data properties
    categorical = categorical,
    multilevel = multilevel,
    clustered = clustered,
    wt = wt,
    swt.type = swt_type,
    # model/options knobs
    meanstructure = lavoptions$meanstructure,
    fixed.x = lavoptions$fixed.x,
    conditional.x = lavoptions$conditional.x,
    correlation = correlation,
    correlation.ov = correlation_ov,
    n.minus.one = lavoptions$gamma.n.minus.one,
    unbiased = lavoptions$gamma.unbiased,
    mplus.wls = mplus_wls,
    group.w.free = lavoptions$group.w.free,
    missing = missing
  )
}

# which variables are scaled to unit variance under a (partial) correlation
# structure? (all of them, unless a subset was requested via
# correlation = c(...))
lav_gamma_recipe_cor_idx <- function(recipe, ov_names_g) {
  if (length(recipe$correlation.ov) > 0L) {
    which(ov_names_g %in% recipe$correlation.ov)
  } else {
    seq_along(ov_names_g)
  }
}
