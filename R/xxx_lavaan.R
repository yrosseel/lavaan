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
# LDW 26 Feb 2024: split lavaan into smaller steps
# 
lavaan <- function(
  # user-specified model: can be syntax, parameter Table, ...
  model              = NULL,                  # not changed
  # data (second argument, most used)
  data               = NULL,                  # changed in step0

  # variable information
  ordered            = NULL,                  # changed in step1.6

  # sampling weights
  sampling.weights   = NULL,                  # not changed

  # summary data
  sample.cov         = NULL,                  # changed in step0
  sample.mean        = NULL,                  # changed in step0
  sample.th          = NULL,                  # changed in step0
  sample.nobs        = NULL,                  # changed in step0

  # multiple groups?
  group              = NULL,                  # not changed

  # multiple levels?
  cluster            = NULL,                  # changed in ldw_adapt_match_call

  # constraints
  constraints        = "",                    # not changed

  # user-specified variance matrices
  WLS.V              = NULL,                  # changed in step0
  NACOV              = NULL,                  # changed in step0

  # internal order of ov.names
  ov.order           = "model",               # changed in step0

  # full slots from previous fits
  slotOptions        = NULL,                  # not changed
  slotParTable       = NULL,                  # not changed
  slotSampleStats    = NULL,                  # not changed
  slotData           = NULL,                  # not changed
  slotModel          = NULL,                  # not changed
  slotCache          = NULL,                  # not changed
  sloth1             = NULL,                  # not changed

  # options (dotdotdot)
  ...
) {
  # start timer
  start.time0 <- proc.time()[3]
  timing <- list()
  timing$start.time <- start.time0
  mc <- match.call(expand.dots = TRUE)
  temp <- ldw_adapt_match_call(matchcall = mc,
                               defaults = NULL,
                               syscall = sys.call(), # to get main arguments without partial matching
                               dotdotdot = list(...))
  mc <- temp[[1]]
  dotdotdot <- temp[[2]]
  rm(temp)
  cluster <- mc$cluster
  # ------------ check data ----------------------
  step0 <- lav_lavaan_step00_checkdata(data, dotdotdot, sample.cov, sample.nobs,
                                     sample.mean, sample.th, NACOV, WLS.V, ov.order)
  data <- step0$data
  dotdotdot <- step0$dotdotdot
  sample.cov <- step0$sample.cov
  sample.nobs <- step0$sample.nobs
  sample.mean <- step0$sample.mean
  sample.th <- step0$sample.th
  NACOV <- step0$NACOV
  WLS.V <- step0$WLS.V
  ov.order <- step0$ov.order
  timing <- ldw_add_timing(timing, "init")
  # ------------ ov.names 1 -----  initial flat model --------------------
  flat.model <- lav_lavaan_step01_ovnames_1(slotParTable, model, dotdotdot$parser)
  # ------------ ov.names 2 ------ handle ov.order -----------------------
  flat.model <- lav_lavaan_step01_ovnames_2(flat.model, ov.order, data, sample.cov, slotData)
  # ------------ ov.names 3 ------- group blocks ------------------
  ngroups <- 1L # default value
  step1 <- lav_lavaan_step01_ovnames_3(flat.model, ov.names, ngroups)
  flat.model <- step1$flat.model
  ov.names <- step1$ov.names
  ov.names.x <- step1$ov.names.x
  ov.names.y <- step1$ov.names.y
  lv.names <- step1$lv.names
  group.values <- step1$group.values
  ngroups <- step1$ngroups
  # ------------ ov.names 4 ------ sanity checks ------------------
  lav_lavaan_step01_ovnames_4(lv.names, data, sample.cov, dotdotdot, slotOptions)
  # ------------ ov.names 5 ------ handle ov.names.l --------------
  step1 <- lav_lavaan_step01_ovnames_5(data, cluster, flat.model, group.values, ngroups)
  flat.model <- step1$flat.model
  ov.names.l <- step1$ov.names.l
  # ------------ ov.names 6 ------ sanity check ordered --------------
  ordered <- lav_lavaan_step01_ovnames_6(ordered, flat.model, data)
  timing <- ldw_add_timing(timing, "ov.names")
  # ------------ lavoptions --------------------
  lavoptions <- lav_lavaan_step02_options(slotOptions, slotData, flat.model, ordered,
                                          sample.cov, sample.mean, sample.th,
                                          sample.nobs, ov.names.l, sampling.weights,
                                          constraints, group, ov.names.x, ov.names.y,
                                          dotdotdot, cluster, data)
  timing <- ldw_add_timing(timing, "Options")
  # fixed.x = FALSE? set ov.names.x = character(0L)
  # new in 0.6-1
  if (!lavoptions$fixed.x) {
    ov.names.x <- character(0L)
  }
  # re-order ov.names.* if requested (new in 0.6-7)
  # ------------ lavdata ------------------------
  step3 <- lav_lavaan_step03_data(slotData, lavoptions, ov.names, ov.names.y, group,
                                  data, cluster, ov.names.x, ov.names.l, ordered,
                                  sampling.weights, sample.cov, sample.mean, sample.th, sample.nobs,
                                  slotParTable, ngroups, dotdotdot, flat.model, model, NACOV, WLS.V)
  lavdata <- step3$lavdata
  lavoptions <- step3$lavoptions
  timing <- ldw_add_timing(timing, "Data")
  # ------------ lavpartable -------------------
  step4 <- lav_lavaan_step04_partable(slotParTable, model, flat.model, lavoptions, lavdata,
                                       constraints)
  lavoptions <- step4$lavoptions
  lavpartable <- step4$lavpartable
  timing <- ldw_add_timing(timing, "ParTable")
  # ------------ lavpta ------------------------
  lavpta <- lav_lavaan_step04_pta(lavpartable, lavoptions) 
  timing <- ldw_add_timing(timing, "lavpta")
  # ------------ lavsamplestats ---------------
  lavsamplestats <- lav_lavaan_step05_samplestats(slotSampleStats, lavdata, lavoptions, WLS.V, NACOV,
                                                  sample.cov, sample.mean, sample.th, sample.nobs,
                                                  ov.names, ov.names.x, lavpta)
  timing <- ldw_add_timing(timing, "SampleStats")
  # ------------ lavh1 ------------------------
  lavh1 <- lav_lavaan_step06_h1(sloth1, lavoptions, lavsamplestats, lavdata, lavpta)
  timing <- ldw_add_timing(timing, "h1")
  # ------------ bounds ------------------------
  lavpartable <- lav_lavaan_step07_bounds(lavoptions, lavh1, lavdata, lavsamplestats, lavpartable)
  timing <- ldw_add_timing(timing, "bounds")
  # ------------ lavstart ----------------------
  lavpartable <- lav_lavaan_step08_start(slotModel, lavoptions, lavpartable, lavsamplestats, lavh1)
  timing <- ldw_add_timing(timing, "start")
  # ------------ model -------------------------
  step9 <- lav_lavaan_step09_model(slotModel, lavoptions, lavpartable, lavsamplestats, lavpta, lavdata)
  lavpartable <- step9$lavpartable
  lavmodel <- step9$lavmodel
  timing <- ldw_add_timing(timing, "Model")
  # -------- lavcache ----------------------------------
  lavcache <- lav_lavaan_step10_cache(slotCache, lavdata, lavmodel, lavpta, lavoptions, sampling.weights)
  timing <- ldw_add_timing(timing, "cache")
  # -------- est + lavoptim ----------------------------
  step11 <- lav_lavaan_step11_estoptim(lavdata, lavmodel, lavpta, lavcache,
                                     lavsamplestats, lavoptions, lavpartable)
  lavoptim <- step11$lavoptim
  lavmodel <- step11$lavmodel
  lavpartable <- step11$lavpartable
  x <- step11$x
  timing <- ldw_add_timing(timing, "optim")
  # -------- lavimplied + lavloglik --------------------
  lavimplied <- lav_lavaan_step12_implied(lavoptions, lavmodel)
  timing <- ldw_add_timing(timing, "implied")
  lavloglik <- lav_lavaan_step12_loglik(lavoptions, lavdata, lavsamplestats, lavimplied, lavmodel)
  timing <- ldw_add_timing(timing, "loglik")
  # ----------- lavvcov + lavboot -------------------
  step13 <- lav_lavaan_step13_vcov_boot(lavoptions, lavmodel, lavsamplestats, lavdata, lavpartable,
                                     lavcache, lavimplied, lavh1, x)
  lavpartable <- step13$lavpartable
  lavvcov <- step13$lavvcov
  VCOV <- step13$VCOV
  lavmodel <- step13$lavmodel
  lavboot <- step13$lavboot
  timing <- ldw_add_timing(timing, "vcov")
  # ----------- lavtest ----------
  lavtest <- lav_lavaan_step14_test(lavoptions, lavmodel, lavsamplestats, lavdata, lavpartable,
                                  lavcache, lavimplied, lavh1, lavpta, x, VCOV, lavloglik)
  timing <- ldw_add_timing(timing, "test")
  # ----------- lavfit ---------- 
  lavfit <- lav_lavaan_step14_fit(lavpartable, lavmodel, lavimplied, x, VCOV, lavtest)
  timing <- ldw_add_timing(timing, "Fit")
  # ----------- baseline ----------------------------
  lavbaseline <- lav_lavaan_step15_baseline(lavoptions, lavsamplestats, lavdata,
                                          lavcache, lavh1, lavpta)
  timing <- ldw_add_timing(timing, "baseline")
  # ----------- rotation ---------------------------
  step16 <- lav_lavaan_step16_rotation(lavoptions, lavmodel, lavpartable, lavh1, lavdata,
                                     x, lavvcov, VCOV, lavcache, lavimplied, lavsamplestats)
  lavpartable <- step16$lavpartable
  lavmodel <- step16$lavmodel
  timing <- ldw_add_timing(timing, "rotation")
  # ------ lavaan result  ----------------
  return(lav_lavaan_step17_lavaan(mc, timing, lavoptions, lavpartable,
                                lavpta, lavdata, lavsamplestats, lavmodel, lavcache, lavfit,
                                lavboot, lavoptim, lavimplied, lavloglik, lavvcov, lavtest,
                                lavh1, lavbaseline, start.time0))
}

# # # # # # # # #
# # cfa + sem # # 
# # # # # # # # #
cfa <- sem <- function(# user-specified model: can be syntax, parameter Table
  model              = NULL,
  # data (second argument, most used)
  data               = NULL,
  
  # variable information
  ordered            = NULL,
  
  # sampling weights
  sampling.weights   = NULL,
  
  # summary data
  sample.cov         = NULL,
  sample.mean        = NULL,
  sample.th          = NULL,
  sample.nobs        = NULL,
  
  # multiple groups?
  group              = NULL,
  
  # multiple levels?
  cluster            = NULL,
  
  # constraints
  constraints        = "",
  
  # user-specified variance matrices
  WLS.V              = NULL,
  NACOV              = NULL,
  
  # internal order of ov.names
  ov.order           = "model",
  
  # options (dotdotdot)
  ...) {
  
  # default options for sem/cfa call
  defaults <- list(
    int.ov.free = TRUE,
    int.lv.free = FALSE,
    auto.fix.first = TRUE, # (re)set in lav_options_set
    auto.fix.single = TRUE,
    auto.var = TRUE,
    auto.cov.lv.x = TRUE,
    auto.cov.y = TRUE,
    auto.th = TRUE,
    auto.delta = TRUE,
    auto.efa = TRUE
  )
  
  # set model.type
  mc <- match.call(expand.dots = TRUE)
  temp <- ldw_adapt_match_call(matchcall = mc,
                               defaults = defaults,
                               syscall = sys.call(), # to get main arguments without partial matching
                               dotdotdot = list(...))
  mc <- temp[[1]]
  # set model.type (cfa or sem)
  mc$model.type <- as.character(mc[[1L]])
  if(length(mc$model.type) == 3L) {
    mc$model.type <- mc$model.type[3L]
  }
  
  # call mother function
  mc[[1L]] <- quote(lavaan::lavaan)
  eval(mc, parent.frame())
}

# # # # # # # # # # # # # # #
# # simple growth models  # # 
# # # # # # # # # # # # # # #
growth <- function(# user-specified model: can be syntax, parameter Table
  model              = NULL,
  # data (second argument, most used)
  data               = NULL,
  
  # variable information
  ordered            = NULL,
  
  # sampling weights
  sampling.weights   = NULL,
  
  # summary data
  sample.cov         = NULL,
  sample.mean        = NULL,
  sample.th          = NULL,
  sample.nobs        = NULL,
  
  # multiple groups?
  group              = NULL,
  
  # multiple levels?
  cluster            = NULL,
  
  # constraints
  constraints        = "",
  
  # user-specified variance matrices
  WLS.V              = NULL,
  NACOV              = NULL,
  
  # internal order of ov.names
  ov.order           = "model",
  
  # options (dotdotdot)
  ...) {
  # default options for growth call
  defaults <- list(
    int.ov.free = FALSE,
    int.lv.free = TRUE,
    auto.fix.first = TRUE, # (re)set in lav_options_set
    auto.fix.single = TRUE,
    auto.var = TRUE,
    auto.cov.lv.x = TRUE,
    auto.cov.y = TRUE,
    auto.th = TRUE,
    auto.delta = TRUE,
    auto.efa = TRUE
  )
  
  mc <- match.call(expand.dots = TRUE)
  temp <- ldw_adapt_match_call(matchcall = mc,
                               defaults = defaults,
                               syscall = sys.call(), # to get main arguments without partial matching
                               dotdotdot = list(...))
  mc <- temp[[1]]
  # set model.type to growth
  mc$model.type <- "growth"
  
  # call mother function
  mc[[1L]] <- quote(lavaan::lavaan)
  eval(mc, parent.frame())
}
# # # # # # # # # # # # # # # # # # #
# # help function ldw_match_call  # #
# # # # # # # # # # # # # # # # # # #
ldw_adapt_match_call <- function(matchcall, defaults, syscall, dotdotdot) {
  # 1. to resolve problem where parameter 'cl' is matched to 'cluster' and shouldn't
  # 2. to apply defaut options for cfa/sem/growth functions
  mc <- matchcall
  sc <- syscall
  ddd <- dotdotdot
  # catch partial matching of 'cl' (expanded to cluster)
  if (!is.null(sc[["cl"]]) && is.null(sc[["cluster"]]) &&
      !is.null(mc[["cluster"]])) {
    mc[["cl"]] <- mc[["cluster"]]
    mc[["cluster"]] <- NULL
    ddd$cl <- sc[["cl"]]
  }
  if (!is.null(mc$cluster)) mc$cluster <- eval(mc$cluster, parent.frame(2))
  # default options
  if (!is.null(defaults)) {
    for (dflt.i in seq_along(defaults)) {
      argname <- names(defaults)[dflt.i]
      if (is.null(mc[[argname]])) mc[[argname]] <- defaults[[dflt.i]]
    }
  }
  return(list(mc, ddd))
}
# # # # # # # # # # # # # # # # # # #
# # help function ldw_add_timing  # #
# # # # # # # # # # # # # # # # # # #
ldw_add_timing <- function(timing, part) {
  # timing is a list with element start.time
  # this function adds an element with name as specified in parameter part and the duration of the interval
  #   from start.time upto now
  # thereafter the element start.time is set to now (prepare for next call)
  # the adapted list is returned
  timenow <- proc.time()[3]
  timing[[part]] <- (timenow - timing$start.time)
  timing$start.time <- timenow
  return(timing)
}
