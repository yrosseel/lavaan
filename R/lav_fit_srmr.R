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
  # ngroups
  G <- lavobject@Data@ngroups

  # container per group
  srmr_mplus.group <- numeric(G)
  srmr_mplus_nomean.group <- numeric(G)

  # If you change how any of the observed/estimated moments below are retrieved,
  # please tag @TDJorgensen at the end of the commit message.
  for (g in 1:G) {
    # observed
    if (!lavobject@SampleStats@missing.flag) {
      if (lavobject@Model@conditional.x) {
        S <- lavobject@SampleStats@res.cov[[g]]
        M <- lavobject@SampleStats@res.int[[g]]
      } else {
        S <- lavobject@SampleStats@cov[[g]]
        M <- lavobject@SampleStats@mean[[g]]
      }
    } else {
      # EM estimates
      S <- lavobject@SampleStats@missing.h1[[g]]$sigma
      M <- lavobject@SampleStats@missing.h1[[g]]$mu
    }
    nvar <- ncol(S)

    # estimated
    implied <- lavobject@implied
    lavmodel <- lavobject@Model
    Sigma.hat <- if (lavmodel@conditional.x) {
      implied$res.cov[[g]]
    } else {
      implied$cov[[g]]
    }
    Mu.hat <- if (lavmodel@conditional.x) {
      implied$res.int[[g]]
    } else {
      implied$mean[[g]]
    }

    # Bollen approach: simply using cov2cor ('correlation residuals')
    S.cor <- cov2cor(S)
    Sigma.cor <- cov2cor(Sigma.hat)
    R.cor <- (S.cor - Sigma.cor)

    # meanstructure
    if (lavobject@Model@meanstructure) {
      # standardized residual mean vector
      R.cor.mean <- M / sqrt(diag(S)) - Mu.hat / sqrt(diag(Sigma.hat))

      e <- nvar * (nvar + 1) / 2 + nvar
      srmr_mplus.group[g] <-
        sqrt((sum(R.cor[lower.tri(R.cor, diag = FALSE)]^2) +
          sum(R.cor.mean^2) +
          sum(((diag(S) - diag(Sigma.hat)) / diag(S))^2)) / e)

      e <- nvar * (nvar + 1) / 2
      srmr_mplus_nomean.group[g] <-
        sqrt((sum(R.cor[lower.tri(R.cor, diag = FALSE)]^2) +
          sum(((diag(S) - diag(Sigma.hat)) / diag(S))^2)) / e)
    } else {
      e <- nvar * (nvar + 1) / 2
      srmr_mplus_nomean.group[g] <- srmr_mplus.group[g] <-
        sqrt((sum(R.cor[lower.tri(R.cor, diag = FALSE)]^2) +
          sum(((diag(S) - diag(Sigma.hat)) / diag(S))^2)) / e)
    }
  } # G

  attr(srmr_mplus.group, "nomean") <- srmr_mplus_nomean.group
  srmr_mplus.group
}

lav_fit_srmr_twolevel <- function(lavobject = NULL) {
  nlevels <- lavobject@Data@nlevels
  G <- lavobject@Data@ngroups

  SRMR.within <- numeric(G)
  SRMR.between <- numeric(G)
  for (g in 1:G) {
    b.within <- (g - 1L) * nlevels + 1L
    b.between <- (g - 1L) * nlevels + 2L

    # OBSERVED        # if these change, tag @TDJorgensen in commit message
    S.within <- lavobject@h1$implied$cov[[b.within]]
    M.within <- lavobject@h1$implied$mean[[b.within]]
    S.between <- lavobject@h1$implied$cov[[b.between]]
    M.between <- lavobject@h1$implied$mean[[b.between]]

    # ESTIMATED       # if these change, tag @TDJorgensen in commit message
    implied <- lav_model_implied_cond2uncond(lavobject@implied)
    Sigma.within <- implied$cov[[b.within]]
    Mu.within <- implied$mean[[b.within]]
    Sigma.between <- implied$cov[[b.between]]
    Mu.between <- implied$mean[[b.between]]

    # force pd for between
    #    S.between <- lav_matrix_symmetric_force_pd(S.between)
    Sigma.between <- lav_matrix_symmetric_force_pd(Sigma.between)

    # Bollen approach: simply using cov2cor ('residual correlations')
    S.within.cor <- cov2cor(S.within)
    S.between.cor <- cov2cor(S.between)
    Sigma.within.cor <- cov2cor(Sigma.within)
    if (all(diag(Sigma.between) > 0)) {
      Sigma.between.cor <- cov2cor(Sigma.between)
    } else {
      Sigma.between.cor <- matrix(as.numeric(NA),
        nrow = nrow(Sigma.between),
        ncol = ncol(Sigma.between)
      )
    }
    R.within.cor <- (S.within.cor - Sigma.within.cor)
    R.between.cor <- (S.between.cor - Sigma.between.cor)

    nvar.within <- NCOL(S.within)
    nvar.between <- NCOL(S.between)
    pstar.within <- nvar.within * (nvar.within + 1) / 2
    pstar.between <- nvar.between * (nvar.between + 1) / 2

    # SRMR
    SRMR.within[g] <- sqrt(sum(lav_matrix_vech(R.within.cor)^2) /
      pstar.within)
    SRMR.between[g] <- sqrt(sum(lav_matrix_vech(R.between.cor)^2) /
      pstar.between)
  }

  # adjust for group sizes
  ng <- unlist(lavobject@SampleStats@nobs)  # if this changes, tag @TDJorgensen in commit message
  ntotal <- lavobject@SampleStats@ntotal    # if this changes, tag @TDJorgensen in commit message
  SRMR_WITHIN <- sum(ng / ntotal * SRMR.within)
  SRMR_BETWEEN <- sum(ng / ntotal * SRMR.between)
  SRMR_TOTAL <- SRMR_WITHIN + SRMR_BETWEEN

  c(SRMR_TOTAL, SRMR_WITHIN, SRMR_BETWEEN)
}

lav_fit_srmr_lavobject <- function(lavobject = NULL, fit.measures = "rmsea") {
  # check lavobject
  stopifnot(inherits(lavobject, "lavaan"))

  # categorical?
  categorical <- lavobject@Model@categorical

  # supported fit measures in this function
  if (categorical) {
    fit.srmr <- c("srmr")
    fit.srmr2 <- c(
      "rmr", "rmr_nomean",
      "srmr", # per default equal to srmr_bentler_nomean
      "srmr_bentler", "srmr_bentler_nomean",
      "crmr", "crmr_nomean",
      "srmr_mplus", "srmr_mplus_nomean"
    )
  } else {
    if (lavobject@Data@nlevels > 1L) {
      fit.srmr <- c("srmr", "srmr_within", "srmr_between")
      fit.srmr2 <- c("srmr", "srmr_within", "srmr_between")
    } else {
      fit.srmr <- c("srmr")
      fit.srmr2 <- c(
        "rmr", "rmr_nomean",
        "srmr", # the default
        "srmr_bentler", "srmr_bentler_nomean",
        "crmr", "crmr_nomean",
        "srmr_mplus", "srmr_mplus_nomean"
      )
    }
  }

  # which one do we need?
  if (missing(fit.measures)) {
    # default set
    fit.measures <- fit.srmr
  } else {
    # remove any not-SRMR related index from fit.measures
    rm.idx <- which(!fit.measures %in% fit.srmr2)
    if (length(rm.idx) > 0L) {
      fit.measures <- fit.measures[-rm.idx]
    }
    if (length(fit.measures) == 0L) {
      return(list())
    }
  }

  # output container
  indices <- list()

  # 1. single level
  if (lavobject@Data@nlevels == 1L) {
    # RMR/SRMR/CRMR: we get it from lav_residuals_summary()
    out <- lav_residuals_summary(lavobject, se = FALSE, unbiased = FALSE)

    cov.cor <- "cov"
    if (categorical) {
      cov.cor <- "cor"
    }

    # only cov
    rmr_nomean.group <- sapply(lapply(out, "[[", "rmr"), "[[", cov.cor)
    srmr_nomean.group <- sapply(lapply(out, "[[", "srmr"), "[[", cov.cor)
    crmr_nomean.group <- sapply(lapply(out, "[[", "crmr"), "[[", cov.cor)

    # total
    if (lavobject@Model@meanstructure) {
      rmr.group <- sapply(lapply(out, "[[", "rmr"), "[[", "total")
      srmr.group <- sapply(lapply(out, "[[", "srmr"), "[[", "total")
      crmr.group <- sapply(lapply(out, "[[", "crmr"), "[[", "total")
    } else {
      # no 'total', only 'cov'
      rmr.group <- rmr_nomean.group
      srmr.group <- srmr_nomean.group
      crmr.group <- crmr_nomean.group
    }

    # the Mplus versions
    srmr_mplus.group <- lav_fit_srmr_mplus(lavobject = lavobject)
    srmr_mplus_nomean.group <- attr(srmr_mplus.group, "nomean")
    attr(srmr_mplus.group, "nomean") <- NULL

    # adjust for group sizes
    ng <- unlist(lavobject@SampleStats@nobs)  # if this changes, tag @TDJorgensen in commit message
    ntotal <- lavobject@SampleStats@ntotal    # if this changes, tag @TDJorgensen in commit message
    RMR <- sum(ng / ntotal * rmr.group)
    RMR_NOMEAN <- sum(ng / ntotal * rmr_nomean.group)
    SRMR_BENTLER <- sum(ng / ntotal * srmr.group)
    SRMR_BENTLER_NOMEAN <- sum(ng / ntotal * srmr_nomean.group)
    CRMR <- sum(ng / ntotal * crmr.group)
    CRMR_NOMEAN <- sum(ng / ntotal * crmr_nomean.group)
    SRMR_MPLUS <- sum(ng / ntotal * srmr_mplus.group)
    SRMR_MPLUS_NOMEAN <- sum(ng / ntotal * srmr_mplus_nomean.group)

    # srmr
    if (lavobject@Options$mimic %in% c("lavaan", "EQS")) {
      if (categorical) {
        indices["srmr"] <- SRMR_BENTLER_NOMEAN
      } else {
        indices["srmr"] <- SRMR_BENTLER
      }
    } else if (lavobject@Options$mimic == "Mplus") {
      if (lavobject@Options$information[1] == "expected") {
        if (categorical) {
          indices["srmr"] <- SRMR_BENTLER_NOMEAN
        } else {
          indices["srmr"] <- SRMR_BENTLER
        }
      } else {
        if (categorical) {
          indices["srmr"] <- SRMR_MPLUS_NOMEAN
        } else {
          indices["srmr"] <- SRMR_MPLUS
        }
      }
    } # Mplus only

    # the others
    indices["srmr_bentler"] <- SRMR_BENTLER
    indices["srmr_bentler_nomean"] <- SRMR_BENTLER_NOMEAN
    indices["crmr"] <- CRMR
    indices["crmr_nomean"] <- CRMR_NOMEAN

    # only correct for non-categorical:
    if (lavobject@Model@categorical) {
      # FIXME! Compute Mplus 8.1 way to compute SRMR in the
      #        categorical setting
      #        See 'SRMR in Mplus (2018)' document on Mplus website
      indices["srmr_mplus"] <- as.numeric(NA)
      indices["srmr_mplus_nomean"] <- as.numeric(NA)
    } else {
      indices["srmr_mplus"] <- SRMR_MPLUS
      indices["srmr_mplus_nomean"] <- SRMR_MPLUS_NOMEAN
    }
    indices["rmr"] <- RMR
    indices["rmr_nomean"] <- RMR_NOMEAN
  } else {
    # 2. twolevel setting
    out <- lav_fit_srmr_twolevel(lavobject = lavobject)
    indices["srmr"] <- out[1]
    indices["srmr_within"] <- out[2]
    indices["srmr_between"] <- out[3]
  } # twolevel

  # return only those that were requested
  indices[fit.measures]
}
