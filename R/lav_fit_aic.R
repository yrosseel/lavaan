# functions related to AIC and other information criteria

# lower-level functions:
# - lav_fit_aic
# - lav_fit_bic
# - lav_fit_sabic

# higher-level functions:
# - lav_fit_aic_lavobject

# Y.R. 21 July 2022

lav_fit_aic <- function(logl = NULL, npar = NULL) {
  aic <- (-2 * logl) + (2 * npar)
  aic
}

lav_fit_bic <- function(logl = NULL, npar = NULL, n = NULL) {
  bic <- (-2 * logl) + (npar * log(n))
  bic
}

lav_fit_sabic <- function(logl = NULL, npar = NULL, n = NULL) {
  n_star <- (n + 2) / 24
  sabic <- (-2 * logl) + (npar * log(n_star))
  sabic
}

lav_fit_aic_lavobject <- function(lavobject = NULL, fit_measures = "aic",
                                  standard_test = "standard",
                                  scaled_test = "none",
                                  estimator = "ML") {
  # check lavobject
  stopifnot(inherits(lavobject, "lavaan"))
  # check object
  lavobject <- lav_object_check_version(lavobject)

  # tests
  test <- lavobject@test
  test_names <- sapply(lavobject@test, "[[", "test")
  if (test_names[1] == "none" || standard_test == "none") {
    return(list())
  }
  test_idx <- which(test_names == standard_test)[1]
  if (length(test_idx) == 0L) {
    return(list())
  }

  scaled_flag <- FALSE
  if (!scaled_test %in% c("none", "standard", "default")) {
    scaled_idx <- which(test_names == scaled_test)
    if (length(scaled_idx) > 0L) {
      scaled_idx <- scaled_idx[1] # only the first one
      scaled_flag <- TRUE
    }
  }

  # estimator?
  if (missing(estimator)) {
    estimator <- lavobject@Options$estimator
  }

  # supported fit measures in this function
  if (estimator == "MML") {
    fit_logl <- c("logl", "aic", "bic", "ntotal", "bic2")
  } else {
    fit_logl <- c(
      "logl", "unrestricted.logl", "aic", "bic",
      "ntotal", "bic2"
    )
  }
  if (scaled_flag &&
    scaled_test %in% c("yuan.bentler", "yuan.bentler.mplus")) {
    fit_logl <- c(fit_logl, "scaling.factor.h1", "scaling.factor.h0")
  }

  # which one do we need?
  if (missing(fit_measures)) {
    # default set
    fit_measures <- fit_logl
  } else {
    # remove any not-CFI related index from fit.measures
    rm_idx <- which(!fit_measures %in% fit_logl)
    if (length(rm_idx) > 0L) {
      fit_measures <- fit_measures[-rm_idx]
    }
    if (length(fit_measures) == 0L) {
      return(list())
    }
  }

  # output container
  indices <- list()

  # non-ML values
  indices["logl"] <- as.numeric(NA)
  indices["unrestricted.logl"] <- as.numeric(NA)
  indices["aic"] <- as.numeric(NA)
  indices["bic"] <- as.numeric(NA)
  indices["ntotal"] <- lavobject@SampleStats@ntotal
  indices["bic2"] <- as.numeric(NA)

  if (estimator %in% c("ML", "MML")) {
    # do we have a @h1 slot?
    if (length(lavobject@h1) > 0L) {
      indices["unrestricted.logl"] <- lavobject@h1$logl$loglik
    } else {
      lavh1 <- lav_h1_implied_logl(
        lavdata = lavobject@Data,
        lavsamplestats = lavobject@SampleStats,
        lavoptions = lavobject@Options
      )
      indices["unrestricted.logl"] <- lavh1$logl$loglik
    }

    # logl H0
    loglik <- lavobject@loglik
    indices["logl"] <- loglik$loglik
    indices["aic"] <- loglik$AIC
    indices["bic"] <- loglik$BIC
    indices["ntotal"] <- loglik$ntotal
    indices["bic2"] <- loglik$BIC2

    # scaling factor for MLR
    if (scaled_test %in% c("yuan.bentler", "yuan.bentler.mplus")) {
      indices["scaling.factor.h1"] <- test[[scaled_idx]]$scaling.factor.h1
      indices["scaling.factor.h0"] <- test[[scaled_idx]]$scaling.factor.h0
    }
  } # ML

  # return only those that were requested
  indices[fit_measures]
}
