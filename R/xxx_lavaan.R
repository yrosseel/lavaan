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
  # full slots from previous fits
  slotOptions        = NULL,
  slotParTable       = NULL,
  slotSampleStats    = NULL,
  slotData           = NULL,
  slotModel          = NULL,
  slotCache          = NULL,
  sloth1             = NULL,
  # options (dotdotdot)
  ...
) {
  # start timer
  start.time0 <- proc.time()[3]
  timing <- list()
  timing$start.time <- start.time0
  # ------------- adapt parameters -----------------
  mc <- match.call(expand.dots = TRUE)
  temp <- lav_lavaan_step00_parameters(matchcall = mc,
                               defaults = NULL,
                               syscall = sys.call(), # to get main arguments without partial matching
                               dotdotdot = list(...))
  lavmc <- temp$mc
  dotdotdot <- temp$dotdotdot
  cluster <- lavmc$cluster
  rm(mc)
  # ------------ check data ----------------------
  temp <- lav_lavaan_step00_checkdata(data, dotdotdot, sample.cov, sample.nobs,
                                     sample.mean, sample.th, NACOV, WLS.V, ov.order)
  data <- temp$data
  dotdotdot <- temp$dotdotdot
  sample.cov <- temp$sample.cov
  sample.nobs <- temp$sample.nobs
  sample.mean <- temp$sample.mean
  sample.th <- temp$sample.th
  NACOV <- temp$NACOV
  WLS.V <- temp$WLS.V
  ov.order <- temp$ov.order
  timing <- ldw_add_timing(timing, "init")
  # ------------ ov.names 1 -----  initial flat model --------------------
  flat.model <- lav_lavaan_step01_ovnames_initflat(slotParTable, model, dotdotdot$parser)
  # ------------ ov.names 2 ------ handle ov.order -----------------------
  flat.model <- lav_lavaan_step01_ovnames_ovorder(flat.model, ov.order, data, sample.cov, slotData)
  # ------------ ov.names 3 ------- group blocks ------------------
  ngroups <- 1L # default value
  temp <- lav_lavaan_step01_ovnames_group(flat.model, ov.names, ngroups)
  flat.model <- temp$flat.model
  ov.names <- temp$ov.names
  ov.names.x <- temp$ov.names.x
  ov.names.y <- temp$ov.names.y
  lv.names <- temp$lv.names
  group.values <- temp$group.values
  ngroups <- temp$ngroups
  # ------------ ov.names 4 ------ sanity checks ------------------
  lav_lavaan_step01_ovnames_checklv(lv.names, data, sample.cov, dotdotdot, slotOptions)
  # ------------ ov.names 5 ------ handle ov.names.l --------------
  temp <- lav_lavaan_step01_ovnames_namesl(data, cluster, flat.model, group.values, ngroups)
  flat.model <- temp$flat.model
  ov.names.l <- temp$ov.names.l
  # ------------ ov.names 6 ------ sanity check ordered --------------
  ordered <- lav_lavaan_step01_ovnames_ordered(ordered, flat.model, data)
  timing <- ldw_add_timing(timing, "ov.names")
  # ------------ lavoptions --------------------
  lavoptions <- lav_lavaan_step02_options(slotOptions, slotData, flat.model, ordered,
                                          sample.cov, sample.mean, sample.th,
                                          sample.nobs, ov.names.l, sampling.weights,
                                          constraints, group, ov.names.x, ov.names.y,
                                          dotdotdot, cluster, data)
  timing <- ldw_add_timing(timing, "Options")
  # ------------ lavdata ------------------------
  temp <- lav_lavaan_step03_data(slotData, lavoptions, ov.names, ov.names.y, group,
                                  data, cluster, ov.names.x, ov.names.l, ordered,
                                  sampling.weights, sample.cov, sample.mean, sample.th, sample.nobs,
                                  slotParTable, ngroups, dotdotdot, flat.model, model, NACOV, WLS.V)
  lavdata <- temp$lavdata
  lavoptions <- temp$lavoptions
  timing <- ldw_add_timing(timing, "Data")
  # ------------ lavpartable -------------------
  temp <- lav_lavaan_step04_partable(slotParTable, model, flat.model, lavoptions, lavdata,
                                       constraints)
  lavoptions <- temp$lavoptions
  lavpartable <- temp$lavpartable
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
  temp <- lav_lavaan_step09_model(slotModel, lavoptions, lavpartable, lavsamplestats, lavpta, lavdata)
  lavpartable <- temp$lavpartable
  lavmodel <- temp$lavmodel
  timing <- ldw_add_timing(timing, "Model")
  # -------- lavcache ----------------------------------
  lavcache <- lav_lavaan_step10_cache(slotCache, lavdata, lavmodel, lavpta, lavoptions, sampling.weights)
  timing <- ldw_add_timing(timing, "cache")
  # -------- est + lavoptim ----------------------------
  temp <- lav_lavaan_step11_estoptim(lavdata, lavmodel, lavpta, lavcache,
                                     lavsamplestats, lavoptions, lavpartable)
  lavoptim <- temp$lavoptim
  lavmodel <- temp$lavmodel
  lavpartable <- temp$lavpartable
  x <- temp$x
  timing <- ldw_add_timing(timing, "optim")
  # -------- lavimplied + lavloglik --------------------
  lavimplied <- lav_lavaan_step12_implied(lavoptions, lavmodel)
  timing <- ldw_add_timing(timing, "implied")
  lavloglik <- lav_lavaan_step12_loglik(lavoptions, lavdata, lavsamplestats, lavimplied, lavmodel)
  timing <- ldw_add_timing(timing, "loglik")
  # ----------- lavvcov + lavboot -------------------
  temp <- lav_lavaan_step13_vcov_boot(lavoptions, lavmodel, lavsamplestats, lavdata, lavpartable,
                                     lavcache, lavimplied, lavh1, x)
  lavpartable <- temp$lavpartable
  lavvcov <- temp$lavvcov
  VCOV <- temp$VCOV
  lavmodel <- temp$lavmodel
  lavboot <- temp$lavboot
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
  temp <- lav_lavaan_step16_rotation(lavoptions, lavmodel, lavpartable, lavh1, lavdata,
                                     x, lavvcov, VCOV, lavcache, lavimplied, lavsamplestats)
  lavpartable <- temp$lavpartable
  lavmodel <- temp$lavmodel
  timing <- ldw_add_timing(timing, "rotation")
  # ------ lavaan result  ----------------
  return(lav_lavaan_step17_lavaan(lavmc, timing, lavoptions, lavpartable,
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
  temp <- lav_lavaan_step00_parameters(matchcall = mc,
                               defaults = defaults,
                               syscall = sys.call(), # to get main arguments without partial matching
                               dotdotdot = list(...))
  mc <- temp$mc
  # set model.type (cfa or sem)
  mc$model.type <- as.character(mc[[1L]])
  if (length(mc$model.type) == 3L) {
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
  temp <- lav_lavaan_step00_parameters(matchcall = mc,
                               defaults = defaults,
                               syscall = sys.call(), # to get main arguments without partial matching
                               dotdotdot = list(...))
  mc <- temp$mc
  # set model.type to growth
  mc$model.type <- "growth"
  # call mother function
  mc[[1L]] <- quote(lavaan::lavaan)
  eval(mc, parent.frame())
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
