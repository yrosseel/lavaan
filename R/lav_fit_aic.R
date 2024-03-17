# functions related to AIC and other information criteria

# lower-level functions:
# - lav_fit_aic
# - lav_fit_bic
# - lav_fit_sabic

# higher-level functions:
# - lav_fit_aic_lavobject

# Y.R. 21 July 2022

lav_fit_aic <- function(logl = NULL, npar = NULL) {
  AIC <- (-2 * logl) + (2 * npar)
  AIC
}

lav_fit_bic <- function(logl = NULL, npar = NULL, N = NULL) {
  BIC <- (-2 * logl) + (npar * log(N))
  BIC
}

lav_fit_sabic <- function(logl = NULL, npar = NULL, N = NULL) {
  N.star <- (N + 2) / 24
  SABIC <- (-2 * logl) + (npar * log(N.star))
  SABIC
}

lav_fit_aic_lavobject <- function(lavobject = NULL, fit.measures = "aic",
                                  standard.test = "standard",
                                  scaled.test = "none",
                                  estimator = "ML") {
  # check lavobject
  stopifnot(inherits(lavobject, "lavaan"))

  # tests
  TEST <- lavobject@test
  test.names <- sapply(lavobject@test, "[[", "test")
  if (test.names[1] == "none" || standard.test == "none") {
    return(list())
  }
  test.idx <- which(test.names == standard.test)[1]
  if (length(test.idx) == 0L) {
    return(list())
  }

  scaled.flag <- FALSE
  if (!scaled.test %in% c("none", "standard", "default")) {
    scaled.idx <- which(test.names == scaled.test)
    if (length(scaled.idx) > 0L) {
      scaled.idx <- scaled.idx[1] # only the first one
      scaled.flag <- TRUE
    }
  }

  # estimator?
  if (missing(estimator)) {
    estimator <- lavobject@Options$estimator
  }

  # supported fit measures in this function
  if (estimator == "MML") {
    fit.logl <- c("logl", "aic", "bic", "ntotal", "bic2")
  } else {
    fit.logl <- c(
      "logl", "unrestricted.logl", "aic", "bic",
      "ntotal", "bic2"
    )
  }
  if (scaled.flag &&
    scaled.test %in% c("yuan.bentler", "yuan.bentler.mplus")) {
    fit.logl <- c(fit.logl, "scaling.factor.h1", "scaling.factor.h0")
  }

  # which one do we need?
  if (missing(fit.measures)) {
    # default set
    fit.measures <- fit.logl
  } else {
    # remove any not-CFI related index from fit.measures
    rm.idx <- which(!fit.measures %in% fit.logl)
    if (length(rm.idx) > 0L) {
      fit.measures <- fit.measures[-rm.idx]
    }
    if (length(fit.measures) == 0L) {
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
    if (.hasSlot(lavobject, "h1") && length(lavobject@h1) > 0L) {
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
    if (.hasSlot(lavobject, "loglik")) {
      loglik <- lavobject@loglik
    } else {
      loglik <- lav_model_loglik(
        lavdata = lavobject@Data,
        lavsamplestats = lavobject@SampleStats,
        lavimplied = lavobject@implied,
        lavmodel = lavobject@Model,
        lavoptions = lavobject@Options
      )
    }
    indices["logl"] <- loglik$loglik
    indices["aic"] <- loglik$AIC
    indices["bic"] <- loglik$BIC
    indices["ntotal"] <- loglik$ntotal
    indices["bic2"] <- loglik$BIC2

    # scaling factor for MLR
    if (scaled.test %in% c("yuan.bentler", "yuan.bentler.mplus")) {
      indices["scaling.factor.h1"] <- TEST[[scaled.idx]]$scaling.factor.h1
      indices["scaling.factor.h0"] <- TEST[[scaled.idx]]$scaling.factor.h0
    }
  } # ML

  # return only those that were requested
  indices[fit.measures]
}
