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
#
lavaan <- function(
    # user specified model: can be syntax, parameter Table, ...
    model = NULL,
    # data (second argument, most used)
    data = NULL,
    # variable information
    ordered = NULL,
    # sampling weights
    sampling.weights = NULL,
    # summary data
    sample.cov = NULL,
    sample.mean = NULL,
    sample.th = NULL,
    sample.nobs = NULL,
    # multiple groups?
    group = NULL,
    # multiple levels?
    cluster = NULL,
    # constraints
    constraints = "",
    # user-specified variance matrices
    WLS.V = NULL, # nolint
    NACOV = NULL, # nolint
    # internal order of ov.names
    ov.order = "model",
    # full slots from previous fits
    slotOptions = NULL, # nolint
    slotParTable = NULL, # nolint
    slotSampleStats = NULL, # nolint
    slotData = NULL, # nolint
    slotModel = NULL, # nolint
    slotCache = NULL, # nolint
    sloth1 = NULL,
    # options (dotdotdot)
    ...) {
  # start timer
  start.time0 <- proc.time()[3]
  timing <- list()
  timing$start.time <- start.time0

  # ------------- adapt parameters -----------------
  mc <- match.call(expand.dots = TRUE)
  temp <- lav_lavaan_step00_parameters(
    matchcall = mc,
    syscall   = sys.call(), # to get main arguments without partial matching
    dotdotdot = list(...)
  )
  lavmc <- temp$mc
  dotdotdot <- temp$dotdotdot
  cluster <- lavmc$cluster
  rm(mc)

  # ------------ check data ----------------------
  temp <- lav_lavaan_step00_checkdata(
    data        = data,
    dotdotdot   = dotdotdot,
    sample.cov  = sample.cov,
    sample.nobs = sample.nobs,
    sample.mean = sample.mean,
    sample.th   = sample.th,
    NACOV       = NACOV,
    WLS.V       = WLS.V,
    ov.order    = ov.order
  )
  data <- temp$data
  dotdotdot <- temp$dotdotdot
  sample.cov <- temp$sample.cov
  sample.nobs <- temp$sample.nobs
  sample.mean <- temp$sample.mean
  sample.th <- temp$sample.th
  NACOV <- temp$NACOV # nolint
  WLS.V <- temp$WLS.V # nolint
  ov.order <- temp$ov.order

  timing <- ldw_add_timing(timing, "init")

  # ------------ ov.names 1 -----  initial flat model --------------------
  flat.model <- lav_lavaan_step01_ovnames_initflat(
    slotParTable     = slotParTable,
    model            = model,
    dotdotdot.parser = dotdotdot$parser
  )

  # ------------ ov.names 2 ------ handle ov.order -----------------------
  flat.model <- lav_lavaan_step01_ovnames_ovorder(
    flat.model = flat.model,
    ov.order   = ov.order,
    data       = data,
    sample.cov = sample.cov,
    slotData   = slotData
  )

  # ------------ ov.names 3 ------- group blocks ------------------
  ngroups <- 1L # default value
  temp <- lav_lavaan_step01_ovnames_group(
    flat.model = flat.model,
    ngroups    = ngroups
  )
  flat.model <- temp$flat.model
  ov.names <- temp$ov.names
  ov.names.x <- temp$ov.names.x
  ov.names.y <- temp$ov.names.y
  lv.names <- temp$lv.names
  group.values <- temp$group.values
  ngroups <- temp$ngroups

  # ------------ ov.names 4 ------ sanity checks ------------------
  lav_lavaan_step01_ovnames_checklv(
    lv.names    = lv.names,
    data        = data,
    sample.cov  = sample.cov,
    dotdotdot   = dotdotdot,
    slotOptions = slotOptions
  )

  # ------------ ov.names 5 ------ handle ov.names.l --------------
  temp <- lav_lavaan_step01_ovnames_namesl(
    data         = data,
    cluster      = cluster,
    flat.model   = flat.model,
    group.values = group.values,
    ngroups      = ngroups
  )
  flat.model <- temp$flat.model
  ov.names.l <- temp$ov.names.l

  # ------------ ov.names 6 ------ sanity check ordered --------------
  ordered <- lav_lavaan_step01_ovnames_ordered(
    ordered    = ordered,
    flat.model = flat.model,
    data       = data
  )
  timing <- ldw_add_timing(timing, "ov.names")

  # ------------ lavoptions --------------------
  lavoptions <- lav_lavaan_step02_options(
    slotOptions      = slotOptions,
    slotData         = slotData,
    flat.model       = flat.model,
    ordered          = ordered,
    sample.cov       = sample.cov,
    sample.mean      = sample.mean,
    sample.th        = sample.th,
    sample.nobs      = sample.nobs,
    ov.names.l       = ov.names.l,
    sampling.weights = sampling.weights,
    constraints      = constraints,
    group            = group,
    ov.names.x       = ov.names.x,
    ov.names.y       = ov.names.y,
    dotdotdot        = dotdotdot,
    cluster          = cluster,
    data             = data
  )
  # fixed.x = FALSE? set ov.names.x = character(0L)
  # new in 0.6-1
  if (!lavoptions$fixed.x) {
    ov.names.x <- character(0L)
  }

  timing <- ldw_add_timing(timing, "Options")

  # ------------ lavdata ------------------------
  temp <- lav_lavaan_step03_data(
    slotData         = slotData,
    lavoptions       = lavoptions,
    ov.names         = ov.names,
    ov.names.y       = ov.names.y,
    group            = group,
    data             = data,
    cluster          = cluster,
    ov.names.x       = ov.names.x,
    ov.names.l       = ov.names.l,
    ordered          = ordered,
    sampling.weights = sampling.weights,
    sample.cov       = sample.cov,
    sample.mean      = sample.mean,
    sample.th        = sample.th,
    sample.nobs      = sample.nobs,
    slotParTable     = slotParTable,
    ngroups          = ngroups,
    dotdotdot        = dotdotdot,
    flat.model       = flat.model,
    model            = model, # in case model is a lavaan object
    NACOV            = NACOV,
    WLS.V            = WLS.V
  )
  lavdata <- temp$lavdata
  lavoptions <- temp$lavoptions

  timing <- ldw_add_timing(timing, "Data")

  # ------------ lavpartable -------------------
  temp <- lav_lavaan_step04_partable(
    slotParTable = slotParTable,
    model        = model,
    flat.model   = flat.model,
    lavoptions   = lavoptions,
    lavdata      = lavdata,
    constraints  = constraints
  )
  lavoptions <- temp$lavoptions
  lavpartable <- temp$lavpartable
  timing <- ldw_add_timing(timing, "ParTable")

  # ------------ lavpta ------------------------
  # lavpta <- lav_lavaan_step04_pta(
  #   lavpartable = lavpartable,
  #   lavoptions  = lavoptions
  # )
  # timing <- ldw_add_timing(timing, "lavpta")

  # ------------ lavsamplestats ---------------
  lavsamplestats <- lav_lavaan_step05_samplestats(
    slotSampleStats = slotSampleStats,
    lavdata         = lavdata,
    lavoptions      = lavoptions,
    WLS.V           = WLS.V,
    NACOV           = NACOV,
    sample.cov      = sample.cov,
    sample.mean     = sample.mean,
    sample.th       = sample.th,
    sample.nobs     = sample.nobs,
    ov.names        = ov.names,
    ov.names.x      = ov.names.x,
    lavpartable     = lavpartable
  )
  timing <- ldw_add_timing(timing, "SampleStats")

  # ------------ lavh1 ------------------------
  lavh1 <- lav_lavaan_step06_h1(
    sloth1         = sloth1,
    lavoptions     = lavoptions,
    lavsamplestats = lavsamplestats,
    lavdata        = lavdata,
    lavpartable    = lavpartable
  )
  timing <- ldw_add_timing(timing, "h1")

  # ------------ bounds ------------------------
  lavpartable <- lav_lavaan_step07_bounds(
    lavoptions     = lavoptions,
    lavh1          = lavh1,
    lavdata        = lavdata,
    lavsamplestats = lavsamplestats,
    lavpartable    = lavpartable
  )
  timing <- ldw_add_timing(timing, "bounds")

  # ------------ lavstart ----------------------
  lavpartable <- lav_lavaan_step08_start(
    slotModel      = slotModel,
    lavoptions     = lavoptions,
    lavpartable    = lavpartable,
    lavsamplestats = lavsamplestats,
    lavh1          = lavh1
  )

  timing <- ldw_add_timing(timing, "start")

  # ------------ model -------------------------
  temp <- lav_lavaan_step09_model(
    slotModel      = slotModel,
    lavoptions     = lavoptions,
    lavpartable    = lavpartable,
    lavsamplestats = lavsamplestats,
    lavdata        = lavdata
  )
  lavpartable <- temp$lavpartable
  lavmodel <- temp$lavmodel

  timing <- ldw_add_timing(timing, "Model")

  # -------- lavcache ----------------------------------
  lavcache <- lav_lavaan_step10_cache(
    slotCache        = slotCache,
    lavdata          = lavdata,
    lavmodel         = lavmodel,
    lavpartable      = lavpartable,
    lavoptions       = lavoptions,
    sampling.weights = sampling.weights
  )
  timing <- ldw_add_timing(timing, "cache")

  # -------- est + lavoptim ----------------------------
  temp <- lav_lavaan_step11_estoptim(
    lavdata        = lavdata,
    lavmodel       = lavmodel,
    lavcache       = lavcache,
    lavsamplestats = lavsamplestats,
    lavoptions     = lavoptions,
    lavpartable    = lavpartable
  )
  lavoptim <- temp$lavoptim
  lavmodel <- temp$lavmodel
  lavpartable <- temp$lavpartable
  x <- temp$x

  timing <- ldw_add_timing(timing, "optim")

  # -------- lavimplied + lavloglik --------------------
  lavimplied <- lav_lavaan_step12_implied(
    lavoptions = lavoptions,
    lavmodel   = lavmodel
  )
  timing <- ldw_add_timing(timing, "implied")

  lavloglik <- lav_lavaan_step12_loglik(
    lavoptions     = lavoptions,
    lavdata        = lavdata,
    lavsamplestats = lavsamplestats,
    lavimplied     = lavimplied,
    lavmodel       = lavmodel
  )
  timing <- ldw_add_timing(timing, "loglik")

  # ----------- lavvcov + lavboot -------------------
  temp <- lav_lavaan_step13_vcov_boot(
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
  VCOV <- temp$VCOV # nolint
  lavmodel <- temp$lavmodel
  lavboot <- temp$lavboot

  timing <- ldw_add_timing(timing, "vcov")

  # ----------- lavtest ----------
  lavtest <- lav_lavaan_step14_test(
    lavoptions     = lavoptions,
    lavmodel       = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata        = lavdata,
    lavpartable    = lavpartable,
    lavcache       = lavcache,
    lavimplied     = lavimplied,
    lavh1          = lavh1,
    x              = x,
    VCOV           = VCOV,
    lavloglik      = lavloglik
  )
  timing <- ldw_add_timing(timing, "test")

  # ----------- lavfit ----------
  lavfit <- lav_lavaan_step14_fit(
    lavpartable = lavpartable,
    lavmodel    = lavmodel,
    lavimplied  = lavimplied,
    x           = x,
    VCOV        = VCOV,
    lavtest     = lavtest
  )
  timing <- ldw_add_timing(timing, "Fit")

  # ----------- baseline ----------------------------
  lavbaseline <- lav_lavaan_step15_baseline(
    lavoptions = lavoptions,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavcache = lavcache,
    lavh1 = lavh1,
    lavpartable = lavpartable
  )
  timing <- ldw_add_timing(timing, "baseline")

  # ----------- rotation ---------------------------
  temp <- lav_lavaan_step16_rotation(
    lavoptions     = lavoptions,
    lavmodel       = lavmodel,
    lavpartable    = lavpartable,
    lavh1          = lavh1,
    lavdata        = lavdata,
    x              = x,
    lavvcov        = lavvcov,
    VCOV           = VCOV,
    lavcache       = lavcache,
    lavimplied     = lavimplied,
    lavsamplestats = lavsamplestats
  )
  lavpartable <- temp$lavpartable
  lavmodel <- temp$lavmodel
  lavvcov <- temp$lavvcov

  timing <- ldw_add_timing(timing, "rotation")

  # ------ lavaan result  ----------------
  out <- lav_lavaan_step17_lavaan(
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
    start.time0    = start.time0
  )

  out
}
# # # # # #
# # cfa # #
# # # # # #
cfa <- function(
  model = NULL,
  data = NULL,
  ordered = NULL,
  sampling.weights = NULL,
  sample.cov = NULL,
  sample.mean = NULL,
  sample.th = NULL,
  sample.nobs = NULL,
  group = NULL,
  cluster = NULL,
  constraints = "",
  WLS.V = NULL, # nolint
  NACOV = NULL, # nolint
  ov.order = "model",
  ...) {
  sc <- sys.call()
  sc[["model.type"]] <- quote("cfa")
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
    sampling.weights = NULL,
    sample.cov = NULL,
    sample.mean = NULL,
    sample.th = NULL,
    sample.nobs = NULL,
    group = NULL,
    cluster = NULL,
    constraints = "",
    WLS.V = NULL, # nolint
    NACOV = NULL, # nolint
    ov.order = "model",
    ...) {
  sc <- sys.call()
  sc[["model.type"]] <- quote("sem")
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
    sampling.weights = NULL,
    sample.cov = NULL,
    sample.mean = NULL,
    sample.th = NULL,
    sample.nobs = NULL,
    group = NULL,
    cluster = NULL,
    constraints = "",
    WLS.V = NULL, # nolint
    NACOV = NULL, # nolint
    ov.order = "model",
    ...) {
  sc <- sys.call()
  sc[["model.type"]] <- quote("growth")
  # call mother function
  sc[[1L]] <- quote(lavaan::lavaan)
  eval(sc, parent.frame())
}

# # # # # # # # # # # # # # # # # # #
# # help function ldw_add_timing  # #
# # # # # # # # # # # # # # # # # # #
ldw_add_timing <- function(timing, part) {
  # timing is a list with element start.time
  # this function adds an element with name as specified in parameter part
  # and the duration of the interval from start.time upto now
  # thereafter the element start.time is set to now (prepare for next call)
  # the adapted list is returned
  timenow <- proc.time()[3]
  timing[[part]] <- (timenow - timing$start.time)
  timing$start.time <- timenow

  timing
}
