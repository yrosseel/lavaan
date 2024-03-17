# functions related to GFI and other 'absolute' fit indices

# lower-level functions:
# - lav_fit_gfi
# - lav_fit_agfi
# - lav_fit_pgfi

# higher-level functions:
# - lav_fit_gfi_lavobject

# Y.R. 21 July 2022

# original formulas were given in Joreskog and Sorbom (1984) user's guide
# for LISREL VI (one for ML, and another for ULS)

# here we use the more 'general' formulas
# (generalized to allow for meanstructures etc)

# References:

# Mulaik, S. A., James, L. R., Van Alstine, J., Bennett, N., Lind, S., &
# Stilwell, C. D. (1989). Evaluation of goodness-of-fit indices for structural
# equation models. Psychological bulletin, 105(3), 430.

# Tanaka, J. S., & Huba, G. J. (1985). A fit index for covariance structure
# models under arbitrary GLS estimation. British Journal of Mathematical and
# Statistical Psychology, 38,197-201.

lav_fit_gfi <- function(WLS.obs = NULL, WLS.est = NULL, WLS.V = NULL,
                        NOBS = NULL) {
  # number of groups
  G <- length(WLS.obs)

  # compute gfi per group
  gfi.group <- numeric(G)
  for (g in 1:G) {
    wls.obs <- WLS.obs[[g]]
    wls.est <- WLS.est[[g]]
    wls.v <- WLS.V[[g]]

    if (is.null(wls.v)) {
      gfi.group[g] <- as.numeric(NA)
    } else {
      wls.diff <- wls.obs - wls.est
      if (is.matrix(wls.v)) {
        # full weight matrix
        t1 <- crossprod(wls.diff, wls.v) %*% wls.diff
        t2 <- crossprod(wls.obs, wls.v) %*% wls.obs
      } else {
        # diagonal weight matrix
        t1 <- as.numeric(crossprod(wls.diff^2, wls.v))
        t2 <- as.numeric(crossprod(wls.obs^2, wls.v))
      }
      gfi.group[g] <- 1 - t1 / t2
    }
  }

  if (G > 1) {
    ## CHECKME: get the scaling right
    NOBS <- unlist(NOBS)
    GFI <- as.numeric((NOBS %*% gfi.group) / sum(NOBS))
  } else {
    GFI <- gfi.group[1L]
  }

  GFI
}

# 'adjusted' GFI (adjusted for degrees of freedom)
lav_fit_agfi <- function(GFI = NULL, nel = NULL, df = NULL) {
  if (!is.finite(GFI) || !is.finite(nel) || !is.finite(df)) {
    AGFI <- as.numeric(NA)
  } else if (df > 0) {
    AGFI <- 1 - (nel / df) * (1 - GFI)
  } else {
    AGFI <- 1
  }

  AGFI
}

# PGFI: parsimony goodness-of-fit index

# Mulaik, S. A., James, L. R., Van Alstine, J., Bennett, N., Lind, S., &
# Stilwell, C. D. (1989). Evaluation of goodness-of-fit indices for structural
# equation models. Psychological bulletin, 105(3), 430.

# LISREL formula (Simplis book 2002, p. 126)
lav_fit_pgfi <- function(GFI = NULL, nel = NULL, df = NULL) {
  if (!is.finite(GFI) || !is.finite(nel) || !is.finite(df)) {
    PGFI <- as.numeric(NA)
  } else if (nel == 0) {
    PGFI <- as.numeric(NA)
  } else {
    PGFI <- (df / nel) * GFI
  }

  PGFI
}


lav_fit_gfi_lavobject <- function(lavobject = NULL, fit.measures = "gfi") {
  # check lavobject
  stopifnot(inherits(lavobject, "lavaan"))

  # possible fit measures
  fit.gfi <- c("gfi", "agfi", "pgfi")

  # which one do we need?
  if (missing(fit.measures)) {
    # default set
    fit.measures <- fit.gfi
  } else {
    # remove any not-GFI related index from fit.measures
    rm.idx <- which(!fit.measures %in% fit.gfi)
    if (length(rm.idx) > 0L) {
      fit.measures <- fit.measures[-rm.idx]
    }
    if (length(fit.measures) == 0L) {
      return(list())
    }
  }

  # extract ingredients
  WLS.obs <- lav_object_inspect_wls_obs(lavobject)
  WLS.est <- lav_object_inspect_wls_est(lavobject)
  WLS.V <- lav_object_inspect_wls_v(lavobject)
  NOBS <- lavobject@SampleStats@nobs

  # compute GFI
  GFI <- lav_fit_gfi(
    WLS.obs = WLS.obs, WLS.est = WLS.est,
    WLS.V = WLS.V, NOBS = NOBS
  )

  # total number of modeled sample stats
  nel <- length(unlist(WLS.obs))

  # degrees of freedom
  df <- lavobject@test[[1]]$df

  # container
  indices <- list()

  indices["gfi"] <- GFI
  indices["agfi"] <- lav_fit_agfi(GFI = GFI, nel = nel, df = df)
  indices["pgfi"] <- lav_fit_pgfi(GFI = GFI, nel = nel, df = df)

  # return only those that were requested
  indices[fit.measures]
}
