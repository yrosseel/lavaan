# keep 'old' names for some function names that have been used
# (or are still being used) by external packages
# rsem ==> tools::testinstalledpackages successful without next line
# computeExpectedInformation <- lav_model_information_expected
# only for simsem ==> tools::testinstalledpackages successful without next line

# standardize function names in lav_utils.R / 31 Oct 2025
getCov <- function(x, lower = TRUE, diagonal = TRUE, sds = NULL,  # nolint start
                   names = paste("V", 1:nvar, sep = "")) {
  lav_deprecated("lav_getcov")
  if (diagonal) {
    nvar <- (sqrt(1 + 8 * length(x)) - 1)/2
  }
  else {
    nvar <- (sqrt(1 + 8 * length(x)) + 1)/2
  }
  sc <- sys.call()
  sc[[1L]] <- quote(lavaan::lav_getcov)
  eval(sc, parent.frame())
}
char2num <- function(s = "") {
    lav_deprecated("lav_char2num")
    sc <- sys.call()
    sc[[1L]] <- quote(lavaan::lav_char2num)
    eval(sc, parent.frame())
}
cor2cov <- function(r, sds, names = NULL, ...) {
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

# standardize function names in lav_data_simulate.R / 9 December 2025
simulateData <- function(
                         model = NULL,
                         model_type = "sem",
                         meanstructure = FALSE,
                         int_ov_free = TRUE,
                         int_lv_free = FALSE,
                         marker_int_zero = FALSE,
                         conditional_x = FALSE,
                         composites = TRUE,
                         fixed_x = FALSE,
                         orthogonal = FALSE,
                         std_lv = TRUE,
                         auto_fix_first = FALSE,
                         auto_fix_single = FALSE,
                         auto_var = TRUE,
                         auto_cov_lv_x = TRUE,
                         auto_cov_y = TRUE,
                         ...,
                         sample_nobs = 500L,
                         ov_var = NULL,
                         group_label = NULL,
                         skewness = NULL,
                         kurtosis = NULL,
                         cluster_idx = NULL,
                         seed = NULL,
                         empirical = FALSE,
                         mass = FALSE,
                         ordered_center = TRUE,
                         return_type = "data.frame",
                         return_fit = FALSE,
                         debug = FALSE,
                         standardized = FALSE) {
  lav_deprecated("lavSimulateData", times = 0L)   #--> for now no warning
  sc <- sys.call()
  # call new function
  sc[[1L]] <- quote(lavaan::lavSimulateData)
  eval(sc, parent.frame())
}
# standardize function names in lav_bootstrap.R / 9 December 2025
bootstrapLavaan <- function(object,
                            r = 1000L,
                            type = "ordinary",
                            verbose = FALSE,
                            fun = "coef",
                            keep_idx = FALSE,
                            parallel = c("no", "multicore", "snow"),
                            ncpus = max(1L, parallel::detectCores() - 2L, na.rm = TRUE),
                            cl = NULL,
                            iseed = NULL,
                            h0_rmsea = NULL,
                            ...) {
    lav_deprecated("lavBootstrap", times = 0L) #--> for now no warning
    sc <- sys.call()
    # call new function
    sc[[1L]] <- quote(lavaan::lavBootstrap)
    eval(sc, parent.frame())
}

# standardize exported functions  / 13 December 2025
inspectSampleCov <- function(model, data, ...) {
    lav_deprecated("lavInspectSampleCov") #--> for now no warning
    sc <- sys.call()
    sc[[1L]] <- quote(lavaan::lavInspectSampleCov)
    eval(sc, parent.frame())
}                                                                  # nolint end

# parameterestimates <- function(object,
#                                se = TRUE, zstat = TRUE, pvalue = TRUE,
#                                ci = TRUE, standardized = FALSE,
#                                fmi = FALSE, plabel = FALSE,
#                                level = 0.95, boot_ci_type = "perc",
#                                cov_std = TRUE, fmi_options = list(),
#                                rsquare = FALSE,
#                                remove_system_eq = TRUE, remove_eq = TRUE,
#                                remove_ineq = TRUE, remove_def = FALSE,
#                                remove_nonfree = FALSE, remove_step1 = TRUE,
#                                remove_unused = FALSE, add_attributes = FALSE,
#                                output = "data.frame", header = FALSE) {
#   lav_deprecated("lavParameterEstimates", times = 1L) #--> for now no warning
#   sc <- sys.call()
#   sc[[1L]] <- quote(lavaan::lavParameterEstimates)
#   eval(sc, parent.frame())
# }
parameterestimates <- lav_alias("lavParameterEstimates") # alias

# for tidySEM in 0.6-21 only
# vnames <- lav_partable_vnames
