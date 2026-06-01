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
  lav_adapt_func(environment(), dotdotdot)
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
  if (!is.null(dotdotdot$composites) && !dotdotdot$composites &&
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
  lavoptim <- temp$lavoptim
  lavmodel <- temp$lavmodel
  lavpartable <- temp$lavpartable
  x <- temp$x
  # store eqs if present in x
  laveqs <- list()
  if (!is.null(attr(x, "eqs"))) {
    laveqs <- attr(x, "eqs")
  }

  timing <- lav_add_timing(timing, "optim")

  # -------- lavimplied + lavloglik --------------------
  lavimplied <- lav_step12_implied(
    lavoptions = lavoptions,
    lavmodel   = lavmodel
  )
  timing <- lav_add_timing(timing, "implied")

  lavloglik <- lav_step12_loglik(
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
