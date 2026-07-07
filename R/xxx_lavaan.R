# # # # # # # # # # # # # # #
# # lavaan main function  # #
# # # # # # # # # # # # # # #
#
# main user-visible cfa/sem/growth functions
#
# initial version: YR 25/03/2009
# added lavoptions YR 02/08/2010
# major revision: YR 9/12/2010: - new workflow (since 0.4-5)
#                               - merge cfa/sem/growth functions
# YR 25/02/2012: changed data slot (from list() to S4); data@X contains data

# YR 26 Jan 2017: use '...' to capture the never-ending list of options
# YR 07 Feb 2023: add ov.order= argument
# HJ 18 Oct 2023: extend PML to allow sampling weights
# LDW 26 Feb 2024: split lavaan in smaller steps
# LDW June 2026: argument names in snake_case

lavaan <- function(
    # user specified model: can be syntax, parameter Table, ...
    model = NULL,
    # data (second argument, most used)
    data = NULL,
    # variable information
    ordered = NULL,
    # auxiliary observed variables (not part of the model)
    aux = NULL,
    # sampling weights
    sampling_weights = NULL,
    # summary data
    sample_cov = NULL,
    sample_mean = NULL,
    sample_th = NULL,
    sample_nobs = NULL,
    # multiple groups?
    group = NULL,
    # multiple levels?
    cluster = NULL,
    # constraints
    constraints = "",
    # user-specified variance matrices
    wls_v = NULL,
    nacov = NULL,
    # internal order of ov.names
    ov_order = "model",
    # full slots from previous fits
    slot_options = NULL,
    slot_par_table = NULL,
    slot_sample_stats = NULL,
    slot_data = NULL,
    slot_model = NULL,
    slot_cache = NULL,
    sloth1 = NULL,
    # options (dotdotdot)
    ...) {
  dotdotdot <- list(...)
  lav_adapt_func(environment(), dotdotdot, FALSE)
  # start timer
  start_time0 <- proc.time()[3]
  timing <- list()
  timing$start_time <- start_time0

  # ------------- adapt parameters -----------------
  mc <- match.call(expand.dots = TRUE)
  temp <- lav_step00_parameters(
    matchcall = mc,
    syscall   = sys.call(), # to get main arguments without partial matching
    dotdotdot = dotdotdot
  )
  lavmc <- temp$mc
  dotdotdot <- temp$dotdotdot
  cluster <- lavmc$cluster
  rm(mc)

  # store current random seed (if any)
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    temp_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    temp_seed <- NULL
  }

  # ------------- handling of warn/debug/verbose switches ----------
  if (!is.null(dotdotdot$debug)) {
    current_debug <- lav_debug()
    if (lav_debug(dotdotdot$debug))
      on.exit(lav_debug(current_debug), TRUE)
    dotdotdot$debug <- NULL
    if (lav_debug()) {
      dotdotdot$warn <- TRUE       # force warnings if debug
      dotdotdot$verbose <- TRUE    # force verbose if debug
    }
  }
  if (!is.null(dotdotdot$warn)) {
    current_warn <- lav_warn()
    if (lav_warn(dotdotdot$warn))
      on.exit(lav_warn(current_warn), TRUE)
    dotdotdot$warn <- NULL
  }
  if (!is.null(dotdotdot$verbose)) {
    current_verbose <- lav_verbose()
    if (lav_verbose(dotdotdot$verbose))
      on.exit(lav_verbose(current_verbose), TRUE)
    dotdotdot$verbose <- NULL
  }

  # ------------ check data ------------------------
  temp <- lav_step00_checkdata(
    data        = data,
    dotdotdot   = dotdotdot,
    sample_cov  = sample_cov,
    sample_nobs = sample_nobs,
    sample_mean = sample_mean,
    sample_th   = sample_th,
    nacov       = nacov,
    wls_v       = wls_v,
    ov_order    = ov_order
  )
  data <- temp$data
  dotdotdot <- temp$dotdotdot
  sample_cov <- temp$sample_cov
  sample_nobs <- temp$sample_nobs
  sample_mean <- temp$sample_mean
  sample_th <- temp$sample_th
  nacov <- temp$nacov
  wls_v <- temp$wls_v
  ov_order <- temp$ov_order

  timing <- lav_add_timing(timing, "init")

  # ------------ ov.names 1 -----  initial flat model --------------------
  # if parser not specified, take default one
  if (is.null(dotdotdot$parser)) {
    opt_default <- lav_options_default()
    useparser <- opt_default$parser
  } else {
    useparser <- dotdotdot$parser
  }

  flat_model <- lav_step01_ovnames_initflat(
    slot_par_table   = slot_par_table,
    model            = model,
    dotdotdot_parser = useparser
  )

  # ------------ ov.names 1b ----- handle 'old way' for composites -------
  if (!is.null(dotdotdot[["composites", exact = TRUE]]) &&
     !dotdotdot[["composites", exact = TRUE]] &&
     any(flat_model$op == "<~")) {
    flat_model <- lav_step01_ovnames_composites(flat_model)
  }

  # ------------ ov.names 2 ------ handle ov_order -----------------------
  flat_model <- lav_step01_ovnames_ovorder(
    flat_model = flat_model,
    ov_order   = ov_order,
    data       = data,
    sample_cov = sample_cov,
    slot_data  = slot_data
  )

  # ------------ ov.names 3 ------- group blocks ------------------
  ngroups <- 1L # default value
  temp <- lav_step01_ovnames_group(
    flat_model = flat_model,
    ngroups    = ngroups
  )
  flat_model <- temp$flat_model
  ov_names <- temp$ov_names
  ov_names_x <- temp$ov_names_x
  ov_names_y <- temp$ov_names_y
  lv_names <- temp$lv_names
  group_values <- temp$group_values
  ngroups <- temp$ngroups

  # ------------ ov.names 4 ------ sanity checks ------------------
  lav_step01_ovnames_checklv(
    lv_names     = lv_names,
    ov_names     = ov_names,
    data         = data,
    sample_cov   = sample_cov,
    dotdotdot    = dotdotdot,
    slot_options = slot_options
  )

  # ------------ ov.names 5 ------ handle ov.names.l --------------
  temp <- lav_step01_ovnames_namesl(
    data         = data,
    cluster      = cluster,
    flat_model   = flat_model,
    group_values = group_values,
    ngroups      = ngroups
  )
  flat_model <- temp$flat_model
  ov_names_l <- temp$ov_names_l

  # ------------ ov.names 6 ------ sanity check ordered --------------
  ordered_orig <- ordered
  ordered <- lav_step01_ovnames_ordered(
    ordered    = ordered,
    flat_model = flat_model,
    data       = data
  )
  timing <- lav_add_timing(timing, "ov.names")

  # ------------ auxiliary variables (FIML saturated correlates) ----------
  # if aux= is supplied and the user requests full-information ML, add the
  # auxiliary variables to the model as a block of 'saturated correlates'
  # (Graham, 2003): they become ordinary (endogenous) observed variables that
  # are freely correlated with each other and with every model variable, and
  # have a free mean. The (enlarged) model is then estimated with casewise
  # FIML, so the auxiliary variables affect the estimates and standard errors.
  # The auxiliary parameters are hidden from the default output. The two-stage
  # missing-data options use a different mechanism (see lav_data_aux_check()).
  aux_satcor         <- FALSE
  aux_fiml_attempted <- FALSE
  aux_names          <- character(0L)
  if (!is.null(aux)) {
    user_missing <- tolower(
      if (!is.null(dotdotdot$missing)) dotdotdot$missing else "")
    user_estimator <- toupper(
      if (!is.null(dotdotdot$estimator)) dotdotdot$estimator else "ML")
    is_fiml <- user_missing %in% c("ml", "fiml", "direct", "direct.ml")
    is_ml_estimator <- user_estimator %in%
      c("", "DEFAULT", "ML", "MLR", "MLM", "MLMV", "MLMVS", "MLF")
    if (is_fiml && is_ml_estimator) {
      aux_fiml_attempted <- TRUE
      aux_names <- lav_aux_clean_names(
        aux = aux, data = data, ov_names = ov_names, ordered = ordered
      )
      if (length(aux_names) > 0L) {
        aux_satcor <- TRUE
        # add the saturated-correlates block to the (flat) model
        flat_model <- lav_aux_satcor_flat(flat_model, aux_names, ov_names)
        # add the auxiliary variables to the observed variable names, so they
        # are read from the data and enter the FIML likelihood (endogenous)
        if (is.list(ov_names)) {
          ov_names   <- lapply(ov_names,   function(x) unique(c(x, aux_names)))
          ov_names_y <- lapply(ov_names_y, function(x) unique(c(x, aux_names)))
        } else {
          ov_names   <- unique(c(ov_names,   aux_names))
          ov_names_y <- unique(c(ov_names_y, aux_names))
        }
        # auxiliary means require a mean structure; auxiliary covariances with
        # exogenous covariates require fixed.x = FALSE
        dotdotdot$meanstructure <- TRUE
        dotdotdot$fixed.x       <- FALSE
        lav_msg_note(gettextf(
          "auxiliary (aux) variable(s) added to the model as saturated
           correlates (estimated with FIML and hidden from the output): %s",
          paste(aux_names, collapse = " ")))
      }
    }
  }

  # ------------ external instruments (estimator = "IV") ----------------
  # instruments supplied with the |~ operator that are NOT otherwise part of
  # the model are 'external' instruments. Like auxiliary variables, they must
  # be read from the data (so that their covariances with the model variables
  # are available to the 2SLS estimator), but they must NOT enter the
  # model-implied summary statistics. They are routed through the generic
  # extra-observed-variable channel (ov_names_aux -> lavData@aux); the IV
  # estimator then augments the covariance matrix with these columns (see
  # lav_sem_miiv.R).
  ext_iv_names <- character(0L)
  if (!aux_satcor) {
    user_estimator2 <- toupper(
      if (!is.null(dotdotdot$estimator)) {
        as.character(dotdotdot$estimator)[1L]
      } else "")
    if (identical(user_estimator2, "IV")) {
      ext_iv_names <- lav_iv_external_names(
        flat_model = flat_model, ov_names = ov_names, data = data)
    }
  }

  # ------------ lavoptions --------------------
  lavoptions <- lav_step02_options(
    slot_options      = slot_options,
    slot_data         = slot_data,
    flat_model       = flat_model,
    ordered          = ordered,
    ordered_orig     = ordered_orig,
    sample_cov       = sample_cov,
    sample_mean      = sample_mean,
    sample_th        = sample_th,
    sample_nobs      = sample_nobs,
    ov_names_l       = ov_names_l,
    sampling_weights = sampling_weights,
    constraints      = constraints,
    group            = group,
    ov_names_x       = ov_names_x,
    ov_names_y       = ov_names_y,
    dotdotdot        = dotdotdot,
    cluster          = cluster,
    data             = data
  )
  # fixed.x = FALSE? set ov_names_x = character(0L)
  # new in 0.6-1
  if (!lavoptions$fixed.x) {
    ov_names_x <- character(0L)
  }
  timing <- lav_add_timing(timing, "Options")

  # ------------ auxiliary variables --------------
  # if the auxiliary variables were already added as saturated correlates
  # (FIML, see above), they are part of the model now; otherwise, validate
  # them for the EM-based (two-stage) channel (see lav_data_aux_check())
  if (aux_satcor) {
    ov_names_aux       <- character(0L)
    lavoptions$aux     <- aux_names
  } else if (aux_fiml_attempted) {
    # FIML saturated-correlates path was attempted, but no usable auxiliary
    # variables remained (already validated/warned above); nothing to do
    ov_names_aux <- character(0L)
  } else if (length(ext_iv_names) > 0L) {
    # external instruments (estimator = "IV"): route through the extra-ov
    # channel so they are read from the data, but keep them out of the model
    ov_names_aux <- ext_iv_names
  } else {
    ov_names_aux <- lav_data_aux_check(
      aux        = aux,
      data       = data,
      ov_names   = ov_names,
      ordered    = ordered,
      lavoptions = lavoptions
    )
  }

  # ------------ lavdata ------------------------
  temp <- lav_step03_data(
    slot_data         = slot_data,
    lavoptions       = lavoptions,
    ov_names         = ov_names,
    ov_names_y       = ov_names_y,
    group            = group,
    data             = data,
    cluster          = cluster,
    ov_names_x       = ov_names_x,
    ov_names_l       = ov_names_l,
    ov_names_aux     = ov_names_aux,
    ordered          = ordered,
    sampling_weights = sampling_weights,
    sample_cov       = sample_cov,
    sample_mean      = sample_mean,
    sample_th        = sample_th,
    sample_nobs      = sample_nobs,
    slot_par_table   = slot_par_table,
    ngroups          = ngroups,
    dotdotdot        = dotdotdot,
    flat_model       = flat_model,
    model            = model, # in case model is a lavaan object
    nacov            = nacov,
    wls_v            = wls_v
  )
  lavdata <- temp$lavdata
  lavoptions <- temp$lavoptions

  timing <- lav_add_timing(timing, "Data")

  # ------------ lavpartable -------------------
  temp <- lav_step04_pt(
    slot_par_table = slot_par_table,
    model          = model,
    flat_model     = flat_model,
    lavoptions     = lavoptions,
    lavdata        = lavdata,
    constraints    = constraints
  )
  lavoptions <- temp$lavoptions
  lavpartable <- temp$lavpartable
  timing <- lav_add_timing(timing, "ParTable")

  # ------------ lavpta ------------------------
  # lavpta <- lav_lavaan_step04_pta(
  #   lavpartable = lavpartable,
  #   lavoptions  = lavoptions
  # )
  # timing <- lav_add_timing(timing, "lavpta")

  # ------------ lavsamplestats ---------------
  lavsamplestats <- lav_step05_samp(
    slot_sample_stats = slot_sample_stats,
    lavdata           = lavdata,
    lavoptions        = lavoptions,
    wls_v             = wls_v,
    nacov             = nacov,
    sample_cov        = sample_cov,
    sample_mean       = sample_mean,
    sample_th         = sample_th,
    sample_nobs       = sample_nobs,
    ov_names          = ov_names,
    ov_names_x        = ov_names_x,
    lavpartable       = lavpartable
  )
  timing <- lav_add_timing(timing, "SampleStats")

  # ------------ lavh1 ------------------------
  lavh1 <- lav_step06_h1(
    sloth1         = sloth1,
    lavoptions     = lavoptions,
    lavsamplestats = lavsamplestats,
    lavdata        = lavdata,
    lavpartable    = lavpartable
  )
  timing <- lav_add_timing(timing, "h1")

  # ------------ two-level WLS statistics -----
  # for two-level (D)WLS/ULS, the sample 'statistics' are the h1 estimates,
  # and their acov (Gamma) is the cluster-based sandwich; both can only be
  # computed once the h1 model has been estimated
  if (lavdata@nlevels > 1L &&
    lavoptions$estimator %in% c("WLS", "DWLS", "ULS")) {
    lavsamplestats <- lav_samp_wls_2l(
      lavsamplestats = lavsamplestats,
      lavh1          = lavh1,
      lavdata        = lavdata,
      lavoptions     = lavoptions
    )
    timing <- lav_add_timing(timing, "wls2l")
  }

  # ------------ adapt marker ------------------
  # if the first indicator (the default marker) of a latent variable turns
  # out to be a poor item, switch to a better marker (and warn); this avoids
  # convergence problems caused by a (very) poor marker item. We keep the
  # original parameter table around: if the switched model fails to converge,
  # we retry with the original marker(s) further below, so that the
  # bad.marker.crit mechanism can never make convergence worse.
  lavpartable_orig <- lavpartable
  marker_switched  <- FALSE
  # composites use marker scaling (first weight fixed to 1) under BOTH
  # auto.fix.first and std.lv (see lav_partable_flat()), so the check must
  # also run when std.lv = TRUE and the model contains composites; the
  # 'clean default' guards in lav_pt_marker_adapt() keep reflective factors
  # (whose loadings are all free under std.lv) out of the check
  if ((isTRUE(lavoptions$auto.fix.first) ||
       (isTRUE(lavoptions$std.lv) && any(lavpartable$op == "<~"))) &&
      lavoptions$bad.marker.crit > 0) {
    adapt <- lav_pt_marker_adapt(
      lavpartable = lavpartable,
      lavh1       = lavh1,
      lavoptions  = lavoptions,
      lavdata     = lavdata,
      threshold   = lavoptions$bad.marker.crit
    )
    if (!is.null(adapt)) {
      info_fac <- adapt$info[adapt$info$type == "factor", , drop = FALSE]
      if (nrow(info_fac) > 0L) {
        lav_msg_warn(gettextf(
          "the first indicator of the following latent variable(s) is a poor
           item; switching to another marker item (to set the metric) to avoid
           convergence problems; use bad.marker.crit = 0 to switch off
           this behavior: %s",
          paste0(info_fac$lv, " (", info_fac$old, " -> ",
                 info_fac$new, ")", collapse = ", ")))
      }
      info_comp <- adapt$info[adapt$info$type == "composite", , drop = FALSE]
      if (nrow(info_comp) > 0L) {
        lav_msg_warn(gettextf(
          "the first indicator of the following composite(s) has a (near)
           zero implied weight; switching to another marker indicator (to set
           the metric) to avoid convergence problems; use bad.marker.crit = 0
           to switch off this behavior: %s",
          paste0(info_comp$lv, " (", info_comp$old, " -> ",
                 info_comp$new, ")", collapse = ", ")))
      }
      # rebuild the parameter table using the new marker(s)
      temp <- lav_step04_pt(
        slot_par_table = slot_par_table,
        model          = model,
        flat_model     = flat_model,
        lavoptions     = lavoptions,
        lavdata        = lavdata,
        constraints    = constraints,
        marker         = adapt$marker
      )
      lavpartable <- temp$lavpartable
      marker_switched <- TRUE
    }
  }

  # ------------ bounds + start + model + cache + optim ------------
  # steps 07-11 are wrapped in a local function so they can be re-run with the
  # original marker(s) as a fail-safe (see the marker switch above and the
  # retry just below)
  lav_bounds_to_optim <- function(lavpartable, timing) {
    # ------------ bounds ------------------------
    lavpartable <- lav_step07_bounds(
      lavoptions     = lavoptions,
      lavh1          = lavh1,
      lavdata        = lavdata,
      lavsamplestats = lavsamplestats,
      lavpartable    = lavpartable
    )
    timing <- lav_add_timing(timing, "bounds")

    # ------------ lavstart ----------------------
    lavpartable <- lav_step08_start(
      slot_model      = slot_model,
      lavoptions      = lavoptions,
      lavpartable     = lavpartable,
      lavsamplestats  = lavsamplestats,
      lavh1           = lavh1
    )
    timing <- lav_add_timing(timing, "start")

    # ------------ model -------------------------
    temp <- lav_step09_model(
      slot_model     = slot_model,
      lavoptions     = lavoptions,
      lavpartable    = lavpartable,
      lavsamplestats = lavsamplestats,
      lavdata        = lavdata
    )
    lavpartable <- temp$lavpartable
    lavmodel <- temp$lavmodel
    timing <- lav_add_timing(timing, "Model")

    # -------- lavcache ----------------------------------
    lavcache <- lav_step10_cache(
      slot_cache       = slot_cache,
      lavdata          = lavdata,
      lavmodel         = lavmodel,
      lavpartable      = lavpartable,
      lavoptions       = lavoptions,
      sampling_weights = sampling_weights
    )
    timing <- lav_add_timing(timing, "cache")

    # -------- est + lavoptim ----------------------------
    temp <- lav_step11_estoptim(
      lavdata        = lavdata,
      lavmodel       = lavmodel,
      lavcache       = lavcache,
      lavsamplestats = lavsamplestats,
      lavh1          = lavh1,
      lavoptions     = lavoptions,
      lavpartable    = lavpartable
    )
    timing <- lav_add_timing(timing, "optim")

    list(lavpartable = temp$lavpartable,
         lavmodel    = temp$lavmodel,
         lavcache    = lavcache,
         lavoptim    = temp$lavoptim,
         x           = temp$x,
         lavoptions  = temp$lavoptions,
         timing      = timing)
  }

  temp <- lav_bounds_to_optim(lavpartable, timing)

  # fail-safe: if we switched marker(s) but the model did not converge, retry
  # with the original marker(s). bad.marker.crit should never make convergence
  # worse: if the original (un-switched) model converges we keep that; if it
  # does not converge either, we still revert (so the result matches what
  # bad.marker.crit = 0 would have produced).
  if (marker_switched && !isTRUE(attr(temp$x, "converged"))) {
    temp_orig <- lav_bounds_to_optim(lavpartable_orig, temp$timing)
    if (isTRUE(attr(temp_orig$x, "converged"))) {
      lav_msg_note(gettext(
        "the model with the switched marker item(s) did not converge;
         reverting to the original marker item(s) (which did converge)."))
    } else {
      lav_msg_note(gettext(
        "the model with the switched marker item(s) did not converge;
         reverting to the original marker item(s)."))
    }
    temp <- temp_orig
  }

  lavpartable <- temp$lavpartable
  lavmodel    <- temp$lavmodel
  lavcache    <- temp$lavcache
  lavoptim    <- temp$lavoptim
  x           <- temp$x
  timing      <- temp$timing
  if (!is.null(temp$lavoptions)) {
    # e.g., the (quiet) em -> nlminb fallback: record the optimizer
    # that was actually used
    lavoptions <- temp$lavoptions
  }

  # store eqs if present in x
  laveqs <- list()
  if (!is.null(attr(x, "eqs"))) {
    laveqs <- attr(x, "eqs")
  }

  # -------- lavimplied + lavloglik --------------------
  lavimplied <- lav_step12_implied(
    lavoptions = lavoptions,
    lavmodel   = lavmodel
  )
  timing <- lav_add_timing(timing, "implied")

  lavloglik <- lav_step12_loglik(
    lavcache = lavcache,
    lavoptions     = lavoptions,
    lavdata        = lavdata,
    lavsamplestats = lavsamplestats,
    lavh1          = lavh1,
    lavimplied     = lavimplied,
    lavmodel       = lavmodel
  )
  timing <- lav_add_timing(timing, "loglik")

  # ----------- lavvcov + lavboot -------------------
  temp <- lav_step13_vcov_boot(
    lavoptions     = lavoptions,
    lavmodel       = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata        = lavdata,
    lavpartable    = lavpartable,
    lavcache       = lavcache,
    lavimplied     = lavimplied,
    lavh1          = lavh1,
    x              = x
  )
  lavpartable <- temp$lavpartable
  lavvcov <- temp$lavvcov
  vcov <- temp$VCOV
  lavmodel <- temp$lavmodel
  lavboot <- temp$lavboot
  lav_monte_carlo <- temp$lav_monte_carlo

  timing <- lav_add_timing(timing, "vcov")

  # ----------- lavtest ----------
  lavtest <- lav_step14_test(
    lavoptions     = lavoptions,
    lavmodel       = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata        = lavdata,
    lavpartable    = lavpartable,
    lavcache       = lavcache,
    lavimplied     = lavimplied,
    lavh1          = lavh1,
    x              = x,
    vcov           = vcov,
    lavloglik      = lavloglik
  )
  timing <- lav_add_timing(timing, "test")

  # ----------- lavfit ----------
  lavfit <- lav_step14_fit(
    lavpartable = lavpartable,
    lavmodel    = lavmodel,
    lavimplied  = lavimplied,
    x           = x,
    vcov        = vcov,
    lavtest     = lavtest
  )
  timing <- lav_add_timing(timing, "Fit")

  # ----------- baseline ----------------------------
  lavbaseline <- lav_step15_baseline(
    lavoptions     = lavoptions,
    lavsamplestats = lavsamplestats,
    lavdata        = lavdata,
    lavcache       = lavcache,
    lavh1          = lavh1,
    lavpartable    = lavpartable
  )
  timing <- lav_add_timing(timing, "baseline")

  # ----------- rotation ---------------------------
  temp <- lav_step16_rotation(
    lavoptions     = lavoptions,
    lavmodel       = lavmodel,
    lavpartable    = lavpartable,
    lavh1          = lavh1,
    lavdata        = lavdata,
    x              = x,
    lavvcov        = lavvcov,
    vcov           = vcov,
    lavcache       = lavcache,
    lavimplied     = lavimplied,
    lavsamplestats = lavsamplestats
  )
  lavpartable <- temp$lavpartable
  lavmodel <- temp$lavmodel
  lavvcov <- temp$lavvcov
  # rotated bootstrap (se = "bootstrap" + rotation) is computed in step 16,
  # after the rotated parameter table exists; keep step 13's lavboot otherwise
  if (!is.null(temp$lavboot)) {
    lavboot <- temp$lavboot
  }

  timing <- lav_add_timing(timing, "rotation")

  # ------ lavaan result  ----------------
  out <- lav_step17_lavaan(
    lavmc          = lavmc,
    timing         = timing,
    lavoptions     = lavoptions,
    lavpartable    = lavpartable,
    lavdata        = lavdata,
    lavsamplestats = lavsamplestats,
    lavmodel       = lavmodel,
    lavcache       = lavcache,
    lavfit         = lavfit,
    lavboot        = lavboot,
    lavoptim       = lavoptim,
    lavimplied     = lavimplied,
    lavloglik      = lavloglik,
    lavvcov        = lavvcov,
    lavtest        = lavtest,
    lavh1          = lavh1,
    lavbaseline    = lavbaseline,
    laveqs         = laveqs,
    start_time0    = start_time0,
    lav_monte_carlo = lav_monte_carlo
  )

  # restore random seed
  if (!is.null(temp_seed)) {
    assign(".Random.seed", temp_seed, envir = .GlobalEnv)             # nolint
  } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    # initially there was no .Random.seed, but we created one along the way
    # clean up
    remove(".Random.seed", envir = .GlobalEnv)
  }

  out
}
# # # # # #
# # cfa # #
# # # # # #
cfa <- function(
  model = NULL,
  data = NULL,
  ordered = NULL,
  aux = NULL,
  sampling_weights = NULL,
  sample_cov = NULL,
  sample_mean = NULL,
  sample_th = NULL,
  sample_nobs = NULL,
  group = NULL,
  cluster = NULL,
  constraints = "",
  wls_v = NULL,
  nacov = NULL,
  ov_order = "model",
  ...) {
  sc <- sys.call()
  sc[["model.type"]] <- quote("cfa")
  sc[["cmd"]] <- quote("cfa")
  # call mother function
  sc[[1L]] <- quote(lavaan::lavaan)
  eval(sc, parent.frame())
}
# # # # # #
# # sem # #
# # # # # #
sem <- function(
    model = NULL,
    data = NULL,
    ordered = NULL,
    aux = NULL,
    sampling_weights = NULL,
    sample_cov = NULL,
    sample_mean = NULL,
    sample_th = NULL,
    sample_nobs = NULL,
    group = NULL,
    cluster = NULL,
    constraints = "",
    wls_v = NULL,
    nacov = NULL,
    ov_order = "model",
    ...) {
  sc <- sys.call()
  sc[["model.type"]] <- quote("sem")
  sc[["cmd"]] <- quote("sem")
  # call mother function
  sc[[1L]] <- quote(lavaan::lavaan)
  eval(sc, parent.frame())
}
# # # # # # # #
# # growth  # #
# # # # # # # #
growth <- function(
    model = NULL,
    data = NULL,
    ordered = NULL,
    aux = NULL,
    sampling_weights = NULL,
    sample_cov = NULL,
    sample_mean = NULL,
    sample_th = NULL,
    sample_nobs = NULL,
    group = NULL,
    cluster = NULL,
    constraints = "",
    wls_v = NULL,
    nacov = NULL,
    ov_order = "model",
    ...) {
  sc <- sys.call()
  sc[["model.type"]] <- quote("growth")
  sc[["cmd"]] <- quote("growth")
  # call mother function
  sc[[1L]] <- quote(lavaan::lavaan)
  eval(sc, parent.frame())
}

# # # # # # # # # # # # # # # # # # #
# # help function lav_add_timing  # #
# # # # # # # # # # # # # # # # # # #
lav_add_timing <- function(timing, part) {
  # timing is a list with element start_time
  # this function adds an element with name as specified in parameter part
  # and the duration of the interval from start_time up to now
  # thereafter the element start_time is set to now (prepare for next call)
  # the adapted list is returned
  timenow <- proc.time()[3]
  timing[[part]] <- (timenow - timing$start_time)
  timing$start_time <- timenow

  timing
}
