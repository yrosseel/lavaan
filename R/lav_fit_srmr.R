# functions related to the SRMR fit measures (single level only)

# lower-level functions:
# - lav_fit_srmr_mplus
# - lav_fit_srmr_twolevel

# higher-level functions:
# - lav_fit_srmr_lavobject

# Y.R. 22 July 2022

# Note: for rmrm/srmr/crmr, we use lav_residuals_summmary()

# SRMR for continuous data only
# see https://www.statmodel.com/download/SRMR.pdf
lav_fit_srmr_mplus <- function(lavobject) {

  lavsamplestats <- lavobject@SampleStats
  lavh1 <- lavobject@h1

  # ngroups
  g_1 <- lavobject@Data@ngroups

  # container per group
  srmr_mplus_group <- numeric(g_1)
  srmr_mplus_nomean_group <- numeric(g_1)

  # If you change how any of the observed/estimated moments below are retrieved,
  # please tag @TDJorgensen at the end of the commit message.
  for (g in 1:g_1) {
    # observed
    if (!lavsamplestats@missing.flag) {
      if (lavobject@Model@conditional.x) {
        s <- lavsamplestats@res.cov[[g]]
        m <- lavsamplestats@res.int[[g]]
      } else {
        s <- lavsamplestats@cov[[g]]
        m <- lavsamplestats@mean[[g]]
      }
    } else {
      # EM estimates
      if (!is.null(lavh1$implied$cov[[g]])) {
        s <- lavh1$implied$cov[[g]]
      } else {
        s <- lavsamplestats@missing.h1[[g]]$sigma
      }
      if (!is.null(lavh1$implied$mean[[g]])) {
        m <- lavh1$implied$mean[[g]]
      } else {
        m <- lavsamplestats@missing.h1[[g]]$mu
      }
    }
    nvar <- ncol(s)

    # estimated
    implied <- lavobject@implied
    lavmodel <- lavobject@Model
    sigma_hat <- if (lavmodel@conditional.x) {
      implied$res.cov[[g]]
    } else {
      implied$cov[[g]]
    }
    mu_hat <- if (lavmodel@conditional.x) {
      implied$res.int[[g]]
    } else {
      implied$mean[[g]]
    }

    # Bollen approach: simply using cov2cor ('correlation residuals')
    s_cor <- cov2cor(s)
    sigma_cor <- cov2cor(sigma_hat)
    r_cor <- (s_cor - sigma_cor)

    # meanstructure
    if (lavobject@Model@meanstructure) {
      # standardized residual mean vector
      r_cor_mean <- m / sqrt(diag(s)) - mu_hat / sqrt(diag(sigma_hat))

      e <- nvar * (nvar + 1) / 2 + nvar
      srmr_mplus_group[g] <-
        sqrt((sum(r_cor[lower.tri(r_cor, diag = FALSE)]^2) +
          sum(r_cor_mean^2) +
          sum(((diag(s) - diag(sigma_hat)) / diag(s))^2)) / e)

      e <- nvar * (nvar + 1) / 2
      srmr_mplus_nomean_group[g] <-
        sqrt((sum(r_cor[lower.tri(r_cor, diag = FALSE)]^2) +
          sum(((diag(s) - diag(sigma_hat)) / diag(s))^2)) / e)
    } else {
      e <- nvar * (nvar + 1) / 2
      srmr_mplus_nomean_group[g] <- srmr_mplus_group[g] <-
        sqrt((sum(r_cor[lower.tri(r_cor, diag = FALSE)]^2) +
          sum(((diag(s) - diag(sigma_hat)) / diag(s))^2)) / e)
    }
  } # G

  attr(srmr_mplus_group, "nomean") <- srmr_mplus_nomean_group
  srmr_mplus_group
}

lav_fit_srmr_twolevel <- function(lavobject = NULL) {
  nlevels <- lavobject@Data@nlevels
  g_1 <- lavobject@Data@ngroups

  srmr_within <- numeric(g_1)
  srmr_between <- numeric(g_1)
  for (g in 1:g_1) {
    b_within <- (g - 1L) * nlevels + 1L
    b_between <- (g - 1L) * nlevels + 2L

    # OBSERVED        # if these change, tag @TDJorgensen in commit message
    s_within <- lavobject@h1$implied$cov[[b_within]]
    # m_within <- lavobject@h1$implied$mean[[b_within]]
    s_between <- lavobject@h1$implied$cov[[b_between]]
    # m_between <- lavobject@h1$implied$mean[[b_between]]

    # ESTIMATED       # if these change, tag @TDJorgensen in commit message
    implied <- lav_model_implied_cond2uncond(lavobject@implied)
    sigma_within <- implied$cov[[b_within]]
    # mu_within <- implied$mean[[b_within]]
    sigma_between <- implied$cov[[b_between]]
    # mu_between <- implied$mean[[b_between]]

    # force pd for between
    #    S.between <- lav_mat_sym_force_pd(S.between)
    sigma_between <- lav_mat_sym_force_pd(sigma_between)

    # Bollen approach: simply using cov2cor ('residual correlations')
    s_within_cor <- cov2cor(s_within)
    s_between_cor <- cov2cor(s_between)
    sigma_within_cor <- cov2cor(sigma_within)
    if (all(diag(sigma_between) > 0)) {
      sigma_between_cor <- cov2cor(sigma_between)
    } else {
      sigma_between_cor <- matrix(as.numeric(NA),
        nrow = nrow(sigma_between),
        ncol = ncol(sigma_between)
      )
    }
    r_within_cor <- (s_within_cor - sigma_within_cor)
    r_between_cor <- (s_between_cor - sigma_between_cor)

    nvar_within <- NCOL(s_within)
    nvar_between <- NCOL(s_between)
    pstar_within <- nvar_within * (nvar_within + 1) / 2
    pstar_between <- nvar_between * (nvar_between + 1) / 2

    # SRMR
    srmr_within[g] <- sqrt(sum(lav_mat_vech(r_within_cor)^2) /
      pstar_within)
    srmr_between[g] <- sqrt(sum(lav_mat_vech(r_between_cor)^2) /
      pstar_between)
  }

  # adjust for group sizes
  ng <- unlist(lavobject@SampleStats@nobs)
              # if this changes, tag @TDJorgensen in commit message
  ntotal <- lavobject@SampleStats@ntotal
              # if this changes, tag @TDJorgensen in commit message
  srmr_within_1 <- sum(ng / ntotal * srmr_within)
  srmr_between_1 <- sum(ng / ntotal * srmr_between)
  srmr_total <- srmr_within_1 + srmr_between_1

  c(srmr_total, srmr_within_1, srmr_between_1)
}

lav_fit_srmr_lavobject <- function(lavobject = NULL, fit_measures = "rmsea") {
  # check lavobject
  stopifnot(inherits(lavobject, "lavaan"))

  # categorical?
  categorical <- lavobject@Model@categorical

  # supported fit measures in this function
  if (categorical) {
    fit_srmr <- c("srmr")
    fit_srmr2 <- c(
      "rmr", "rmr_nomean",
      "srmr", # per default equal to srmr_bentler_nomean
      "srmr_bentler", "srmr_bentler_nomean",
      "crmr", "crmr_nomean",
      "srmr_mplus", "srmr_mplus_nomean"
    )
  } else {
    if (lavobject@Data@nlevels > 1L) {
      fit_srmr <- c("srmr", "srmr_within", "srmr_between")
      fit_srmr2 <- c("srmr", "srmr_within", "srmr_between")
    } else {
      fit_srmr <- c("srmr")
      fit_srmr2 <- c(
        "rmr", "rmr_nomean",
        "srmr", # the default
        "srmr_bentler", "srmr_bentler_nomean",
        "crmr", "crmr_nomean",
        "srmr_mplus", "srmr_mplus_nomean"
      )
    }
  }

  # which one do we need?
  if (missing(fit_measures)) {
    # default set
    fit_measures <- fit_srmr
  } else {
    # remove any not-SRMR related index from fit.measures
    rm_idx <- which(!fit_measures %in% fit_srmr2)
    if (length(rm_idx) > 0L) {
      fit_measures <- fit_measures[-rm_idx]
    }
    if (length(fit_measures) == 0L) {
      return(list())
    }
  }

  # output container
  indices <- list()

  # 1. single level
  if (lavobject@Data@nlevels == 1L) {
    # RMR/SRMR/CRMR: we get it from lav_residuals_summary()
    out <- lav_residuals_summary(lavobject, se = FALSE, unbiased = FALSE)

    # name of the covariance-residual column in the summary tables
    cov_cor <- "cov"
    if (categorical) {
      # conditional.x categorical summary uses res.cov / res.th / res.slopes
      cov_cor <- if (lavobject@Model@conditional.x) "res.cov" else "cor"
    } else if (lavobject@Model@conditional.x) {
      # continuous conditional.x summary uses res.cov / res.int / res.slopes
      cov_cor <- "res.cov"
    }

    # only cov
    rmr_nomean_group <- sapply(lapply(out, "[[", "rmr"), "[[", cov_cor)
    srmr_nomean_group <- sapply(lapply(out, "[[", "srmr"), "[[", cov_cor)
    crmr_nomean_group <- sapply(lapply(out, "[[", "crmr"), "[[", cov_cor)

    # total
    if (lavobject@Model@meanstructure) {
      rmr_group <- sapply(lapply(out, "[[", "rmr"), "[[", "total")
      srmr_group <- sapply(lapply(out, "[[", "srmr"), "[[", "total")
      crmr_group <- sapply(lapply(out, "[[", "crmr"), "[[", "total")
    } else {
      # no 'total', only 'cov'
      rmr_group <- rmr_nomean_group
      srmr_group <- srmr_nomean_group
      crmr_group <- crmr_nomean_group
    }

    # the Mplus versions
    srmr_mplus_group <- lav_fit_srmr_mplus(lavobject = lavobject)
    srmr_mplus_nomean_group <- attr(srmr_mplus_group, "nomean")
    attr(srmr_mplus_group, "nomean") <- NULL

    # adjust for group sizes
    ng <- unlist(lavobject@SampleStats@nobs)
              # if this changes, tag @TDJorgensen in commit message
    ntotal <- lavobject@SampleStats@ntotal
              # if this changes, tag @TDJorgensen in commit message
    rmr <- sum(ng / ntotal * rmr_group)
    rmr_nomean <- sum(ng / ntotal * rmr_nomean_group)
    srmr_bentler <- sum(ng / ntotal * srmr_group)
    srmr_bentler_nomean <- sum(ng / ntotal * srmr_nomean_group)
    crmr <- sum(ng / ntotal * crmr_group)
    crmr_nomean <- sum(ng / ntotal * crmr_nomean_group)
    srmr_mplus <- sum(ng / ntotal * srmr_mplus_group)
    srmr_mplus_nomean <- sum(ng / ntotal * srmr_mplus_nomean_group)

    # srmr
    if (lavobject@Options$mimic %in% c("lavaan", "EQS")) {
      if (categorical) {
        indices["srmr"] <- srmr_bentler_nomean
      } else {
        indices["srmr"] <- srmr_bentler
      }
    } else if (lavobject@Options$mimic == "Mplus") {
      if (lavobject@Options$information[1] == "expected") {
        if (categorical) {
          indices["srmr"] <- srmr_bentler_nomean
        } else {
          indices["srmr"] <- srmr_bentler
        }
      } else {
        if (categorical) {
          indices["srmr"] <- srmr_mplus_nomean
        } else {
          indices["srmr"] <- srmr_mplus
        }
      }
    } # Mplus only

    # the others
    indices["srmr_bentler"] <- srmr_bentler
    indices["srmr_bentler_nomean"] <- srmr_bentler_nomean
    indices["crmr"] <- crmr
    indices["crmr_nomean"] <- crmr_nomean

    # only correct for non-categorical:
    if (lavobject@Model@categorical) {
      # FIXME! Compute Mplus 8.1 way to compute SRMR in the
      #        categorical setting
      #        See 'SRMR in Mplus (2018)' document on Mplus website
      indices["srmr_mplus"] <- as.numeric(NA)
      indices["srmr_mplus_nomean"] <- as.numeric(NA)
    } else {
      indices["srmr_mplus"] <- srmr_mplus
      indices["srmr_mplus_nomean"] <- srmr_mplus_nomean
    }
    indices["rmr"] <- rmr
    indices["rmr_nomean"] <- rmr_nomean
  } else {
    # 2. twolevel setting
    out <- lav_fit_srmr_twolevel(lavobject = lavobject)
    indices["srmr"] <- out[1]
    indices["srmr_within"] <- out[2]
    indices["srmr_between"] <- out[3]
  } # twolevel

  # return only those that were requested
  indices[fit_measures]
}
