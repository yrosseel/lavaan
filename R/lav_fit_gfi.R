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

lav_fit_gfi <- function(wls_obs = NULL, wls_est = NULL, wls_v = NULL,
                        nobs_1 = NULL) {
  # number of groups
  g_1 <- length(wls_obs)

  # compute gfi per group
  gfi_group <- numeric(g_1)
  for (g in 1:g_1) {
    wls_obs_1 <- wls_obs[[g]]
    wls_est_1 <- wls_est[[g]]
    wls_v_1 <- wls_v[[g]]

    if (is.null(wls_v_1)) {
      gfi_group[g] <- as.numeric(NA)
    } else {
      wls_diff <- wls_obs_1 - wls_est_1
      if (is.matrix(wls_v_1)) {
        # full weight matrix
        t1 <- crossprod(wls_diff, wls_v_1) %*% wls_diff
        t2 <- crossprod(wls_obs_1, wls_v_1) %*% wls_obs_1
      } else {
        # diagonal weight matrix
        t1 <- as.numeric(crossprod(wls_diff^2, wls_v_1))
        t2 <- as.numeric(crossprod(wls_obs_1^2, wls_v_1))
      }
      gfi_group[g] <- 1 - t1 / t2
    }
  }

  if (g_1 > 1) {
    ## CHECKME: get the scaling right
    nobs_1 <- unlist(nobs_1)
    gfi <- as.numeric((nobs_1 %*% gfi_group) / sum(nobs_1))
  } else {
    gfi <- gfi_group[1L]
  }

  gfi
}

# 'adjusted' GFI (adjusted for degrees of freedom)
lav_fit_agfi <- function(gfi = NULL, nel = NULL, df = NULL) {
  if (!is.finite(gfi) || !is.finite(nel) || !is.finite(df)) {
    agfi <- as.numeric(NA)
  } else if (df > 0) {
    agfi <- 1 - (nel / df) * (1 - gfi)
  } else {
    agfi <- 1
  }

  agfi
}

# PGFI: parsimony goodness-of-fit index

# Mulaik, S. A., James, L. R., Van Alstine, J., Bennett, N., Lind, S., &
# Stilwell, C. D. (1989). Evaluation of goodness-of-fit indices for structural
# equation models. Psychological bulletin, 105(3), 430.

# LISREL formula (Simplis book 2002, p. 126)
lav_fit_pgfi <- function(gfi = NULL, nel = NULL, df = NULL) {
  if (!is.finite(gfi) || !is.finite(nel) || !is.finite(df)) {
    pgfi <- as.numeric(NA)
  } else if (nel == 0) {
    pgfi <- as.numeric(NA)
  } else {
    pgfi <- (df / nel) * gfi
  }

  pgfi
}


lav_fit_gfi_lavobject <- function(lavobject = NULL, fit_measures = "gfi") {
  # check lavobject
  stopifnot(inherits(lavobject, "lavaan"))

  # possible fit measures
  fit_gfi <- c("gfi", "agfi", "pgfi")

  # which one do we need?
  if (missing(fit_measures)) {
    # default set
    fit_measures <- fit_gfi
  } else {
    # remove any not-GFI related index from fit.measures
    rm_idx <- which(!fit_measures %in% fit_gfi)
    if (length(rm_idx) > 0L) {
      fit_measures <- fit_measures[-rm_idx]
    }
    if (length(fit_measures) == 0L) {
      return(list())
    }
  }

  # extract ingredients
  wls_obs <- lav_object_inspect_wls_obs(lavobject)
  wls_est <- lav_object_inspect_wls_est(lavobject)
  wls_v <- lav_object_inspect_wls_v(lavobject)
  nobs_1 <- lavobject@SampleStats@nobs

  # compute GFI
  gfi <- lav_fit_gfi(
    wls_obs = wls_obs, wls_est = wls_est,
    wls_v = wls_v, nobs_1 = nobs_1
  )

  # total number of modeled sample stats
  nel <- length(unlist(wls_obs))

  # degrees of freedom
  df <- lavobject@test[[1]]$df

  # container
  indices <- list()

  indices["gfi"] <- gfi
  indices["agfi"] <- lav_fit_agfi(gfi = gfi, nel = nel, df = df)
  indices["pgfi"] <- lav_fit_pgfi(gfi = gfi, nel = nel, df = df)

  # return only those that were requested
  indices[fit_measures]
}
