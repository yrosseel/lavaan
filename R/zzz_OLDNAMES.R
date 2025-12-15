# keep 'old' names for some function names that have been used
# (or are still being used) by external packages
# rsem
computeExpectedInformation <- lav_model_information_expected          # nolint start
# only for simsem ....
getParameterLabels <- lav_partable_labels

# standardize function names in lav_utils.R / 31 Oct 2025
getCov <- function(x, lower = TRUE, diagonal = TRUE, sds = NULL,
                   names = paste("V", 1:nvar, sep = "")) {
  lav_deprecated("lav_get_cov")
  if (diagonal) {
    nvar <- (sqrt(1 + 8 * length(x)) - 1)/2
  }
  else {
    nvar <- (sqrt(1 + 8 * length(x)) + 1)/2
  }
  sc <- sys.call()
  sc[[1L]] <- quote(lavaan::lav_get_cov)
  eval(sc, parent.frame())
}
char2num <- function(s = "") {
    lav_deprecated("lav_char2num")
    sc <- sys.call()
    sc[[1L]] <- quote(lavaan::lav_char2num)
    eval(sc, parent.frame())
}
cor2cov <- function(R, sds, names = NULL) {
    lav_deprecated("lav_cor2cov")
    sc <- sys.call()
    sc[[1L]] <- quote(lavaan::lav_cor2cov)
    eval(sc, parent.frame())
}

# standardize function names in lav_mplus_lavaan / 7 Nov 2025
mplus2lavaan.modelSyntax <- function(syntax) {
    lav_deprecated("lav_mplus_syntax_model")
    sc <- sys.call()
    sc[[1L]] <- quote(lavaan::lav_mplus_syntax_model)
    eval(sc, parent.frame())
}
mplus2lavaan  <- function(inpfile, run = TRUE)  {
    lav_deprecated("lav_mplus_lavaan")
    sc <- sys.call()
    sc[[1L]] <- quote(lavaan::lav_mplus_lavaan)
    eval(sc, parent.frame())
}

# standardize function names in lav_partable_vnames.R / 9 December 2025
lavaanNames <- function(object, type = "ov", ...) {
    lav_deprecated("lavNames", times = 0L) #--> for now no warning
    sc <- sys.call()
    sc[[1L]] <- quote(lavaan::lavNames)
    eval(sc, parent.frame())
}

# standardize function names in lav_simulate_old.R / 9 December 2025
simulateData <- function(
                         model = NULL,
                         model.type = "sem",
                         meanstructure = FALSE,
                         int.ov.free = TRUE,
                         int.lv.free = FALSE,
                         marker.int.zero = FALSE,
                         conditional.x = FALSE,
                         composites = TRUE,
                         fixed.x = FALSE,
                         orthogonal = FALSE,
                         std.lv = TRUE,
                         auto.fix.first = FALSE,
                         auto.fix.single = FALSE,
                         auto.var = TRUE,
                         auto.cov.lv.x = TRUE,
                         auto.cov.y = TRUE,
                         ...,
                         sample.nobs = 500L,
                         ov.var = NULL,
                         group.label = paste("G", 1:ngroups, sep = ""),
                         skewness = NULL,
                         kurtosis = NULL,
                         seed = NULL,
                         empirical = FALSE,
                         return.type = "data.frame",
                         return.fit = FALSE,
                         debug = FALSE,
                         standardized = FALSE) {
  lav_deprecated("lavSimulateData", times = 0L) #--> for now no warning
  if (is.list(model)) {
    ngroups <- lav_partable_ngroups(model)
  } else {
    ngroups <- sample.nobs
  }
  sc <- sys.call()
  # call new function
  sc[[1L]] <- quote(lavaan::lavSimulateData)
  eval(sc, parent.frame())
}
# standardize function names in lav_bootstrap.R / 9 December 2025
bootstrapLavaan <- function(object,
                            R = 1000L,
                            type = "ordinary",
                            verbose = FALSE,
                            FUN = "coef",
                            keep.idx = FALSE,
                            parallel = c("no", "multicore", "snow"),
                            ncpus = max(1L, parallel::detectCores() - 2L),
                            cl = NULL,
                            iseed = NULL,
                            h0.rmsea = NULL,
                            ...) {
    lav_deprecated("lavBootstrap", times = 0L) #--> for now no warning
    sc <- sys.call()
    # call new function
    sc[[1L]] <- quote(lavaan::lavBootstrap)
    eval(sc, parent.frame())
}

inspect <- function(object, what = "free", ...) {
    lav_deprecated("lavInspect", times = 0L) #--> for now no warning
    UseMethod("inspect", object)
}
# standardize exported functions  / 13 December 2025
inspectSampleCov <- function(model, data, ...) {
    lav_deprecated("lavInspectSampleCov", times = 0L) #--> for now no warning
    sc <- sys.call()
    sc[[1L]] <- quote(lavaan::lavInspectSampleCov)
    eval(sc, parent.frame())
}                                                                  # nolint end
parameterEstimates <- function(object,
                   se = TRUE, zstat = TRUE, pvalue = TRUE, ci = TRUE,
                   standardized = FALSE,
                   fmi = FALSE, plabel = FALSE,
                   level = 0.95, boot.ci.type = "perc",
                   cov.std = TRUE, fmi.options = list(),
                   rsquare = FALSE,
                   remove.system.eq = TRUE, remove.eq = TRUE,
                   remove.ineq = TRUE, remove.def = FALSE,
                   remove.nonfree = FALSE, remove.step1 = TRUE,
                   remove.unused = FALSE, add.attributes = FALSE,
                   output = "data.frame", header = FALSE) {
  lav_deprecated("lavParameterEstimates", times = 0L) #--> for now no warning
  sc <- sys.call()
  sc[[1L]] <- quote(lavaan::lavParameterEstimates)
  eval(sc, parent.frame())
}
parameterestimates <- function(object,
                               se = TRUE, zstat = TRUE, pvalue = TRUE, ci = TRUE,
                               standardized = FALSE,
                               fmi = FALSE, plabel = FALSE,
                               level = 0.95, boot.ci.type = "perc",
                               cov.std = TRUE, fmi.options = list(),
                               rsquare = FALSE,
                               remove.system.eq = TRUE, remove.eq = TRUE,
                               remove.ineq = TRUE, remove.def = FALSE,
                               remove.nonfree = FALSE, remove.step1 = TRUE,
                               remove.unused = FALSE, add.attributes = FALSE,
                               output = "data.frame", header = FALSE) {
  lav_deprecated("lavParameterEstimates", times = 0L) #--> for now no warning
  sc <- sys.call()
  sc[[1L]] <- quote(lavaan::lavParameterEstimates)
  eval(sc, parent.frame())
}
