# Browne's residual test statistic
# see Browne (1984) eq 2.20a

# T.B = (N-1) * t(RES) %*% Delta.c %*%
#                          solve(t(Delta.c) %*% Gamma %*% Delta.c) %*%
#                          t(Delta.c) %*% RES
#
#     = (N-1) * t(RES) %*% (Gamma.inv -
#                           Gamma.inv %*% Delta %*%
#                                  solve(t(Delta) %*% Gamma.inv %*% Delta) %*%
#                           t(Delta) %*% Gamma.inv) %*% RES

# Note: if Gamma == solve(Weight matrix), then:
# t(Delta) %*% solve(Gamma) %*% RES == 0-vector!
#
# Therefore:
# - if estimator = "WLS", X2 == Browne's residual ADF statistic
# - if estimator = "GLS", X2 == Browne's residual NT statistic
#
# - if estimator = "NTRLS", X2 == Browne's residual NT statistic (model-based)
#   also known as the RLS test statistic
#   ... except in multigroup + equality constraints, where
#   t(Delta) %*% solve(Gamma) %*% RES not zero everywhere!?

# YR 26 July 2022: add alternative slots, if lavobject = NULL
# YR 22 Jan  2023: allow for model-based 'structured' Sigma

# TODo: - allow for non-linear equality constraints
#         (see Browne, 1982, eq 1.7.19; although we may face singular matrices)

lav_test_browne <- function(lavobject = NULL,
                            # or
                            lavdata = NULL,
                            lavsamplestats = NULL, # WLS.obs, NACOV
                            lavmodel = NULL,
                            lavpartable = NULL, # DF
                            lavoptions = NULL,
                            lavh1 = NULL,
                            lavimplied = NULL,
                            # further options:
                            n.minus.one = "default",
                            ADF = TRUE,
                            model.based = FALSE) {
  if (!is.null(lavobject)) {
    # check input
    if (!inherits(lavobject, "lavaan")) {
      lav_msg_stop(gettext("object is not a lavaan object."))
    }

    # slots
    lavdata <- lavobject@Data
    lavsamplestats <- lavobject@SampleStats
    lavmodel <- lavobject@Model
    lavpartable <- lavobject@ParTable
    lavoptions <- lavobject@Options
    lavh1 <- lavobject@h1
    lavimplied <- lavobject@implied
  }

  if (!ADF && lavmodel@categorical) {
    lav_msg_stop(gettext("normal theory version not available in the categorical
                         setting."))
  }
  if (lavdata@missing != "listwise" && !model.based) {
    lav_msg_stop(gettext("Browne's test is not available when data is missing"))
  }
  if (lavdata@nlevels > 1L) {
    lav_msg_stop(gettext("Browne's test is not available when data is
                         multilevel."))
  }
  if (length(lavmodel@ceq.nonlinear.idx) > 0L) {
    lav_msg_stop(gettext("Browne's test is not available (yet) when nonlinear
                         equality constraints are involved."))
  }

  if (!is.logical(n.minus.one)) {
    if (lavoptions$estimator == "ML" &&
      lavoptions$likelihood == "normal") {
      n.minus.one <- FALSE
    } else {
      n.minus.one <- TRUE
    }
  }

  # ingredients
  Delta <- computeDelta(lavmodel)
  if (ADF) {
    # ADF version
    if (!is.null(lavsamplestats@NACOV[[1]])) {
      Gamma <- lavsamplestats@NACOV
    } else {
      if (!is.null(lavobject)) {
        if (lavobject@Data@data.type != "full") {
          lav_msg_stop(gettext("ADF version not available without full data or
                               user-provided Gamma/NACOV matrix"))
        }
        Gamma <- lav_object_gamma(lavobject,
          ADF = TRUE,
          model.based = model.based
        )
      } else {
        if (lavdata@data.type != "full") {
          lav_msg_stop(gettext("ADF version not available without full data or
                               user-provided Gamma/NACOV matrix"))
        }
        Gamma <- lav_object_gamma(
          lavobject = NULL,
          lavdata = lavdata,
          lavoptions = lavoptions,
          lavsamplestats = lavsamplestats,
          lavh1 = lavh1,
          lavimplied = lavimplied,
          ADF = TRUE,
          model.based = model.based
        )
      }
    }
  } else {
    # NT version
    if (!is.null(lavobject)) {
      Gamma <- lav_object_gamma(lavobject,
        ADF = FALSE,
        model.based = model.based
      )
    } else {
      Gamma <- lav_object_gamma(
        lavobject = NULL,
        lavdata = lavdata,
        lavoptions = lavoptions,
        lavsamplestats = lavsamplestats,
        lavh1 = lavh1,
        lavimplied = lavimplied,
        ADF = FALSE,
        model.based = model.based
      )
    }
  }
  WLS.obs <- lavsamplestats@WLS.obs
  WLS.est <- lav_model_wls_est(lavmodel)
  nobs <- lavsamplestats@nobs
  ntotal <- lavsamplestats@ntotal

  # linear equality constraints?
  lineq.flag <- FALSE
  if (lavmodel@eq.constraints) {
    lineq.flag <- TRUE
  } else if (.hasSlot(lavmodel, "ceq.simple.only") &&
    lavmodel@ceq.simple.only) {
    lineq.flag <- TRUE
  }

  # compute T.B per group
  ngroups <- length(WLS.obs)
  stat.group <- numeric(ngroups)

  # 1. standard setting: no equality constraints
  if (!lineq.flag) {
    for (g in seq_len(ngroups)) {
      RES <- WLS.obs[[g]] - WLS.est[[g]]
      Delta.g <- Delta[[g]]
      Delta.c <- lav_matrix_orthogonal_complement(Delta.g)
      tDGD <- crossprod(Delta.c, Gamma[[g]]) %*% Delta.c
      # if fixed.x = TRUE, Gamma[[g]] may contain zero col/rows
      tDGD.inv <- lav_matrix_symmetric_inverse(tDGD)
      if (n.minus.one) {
        Ng <- nobs[[g]] - 1L
      } else {
        Ng <- nobs[[g]]
      }
      tResDelta.c <- crossprod(RES, Delta.c)
      stat.group[g] <-
        Ng * drop(tResDelta.c %*% tDGD.inv %*% t(tResDelta.c))
    }
    STAT <- sum(stat.group)

    # 2. linear equality constraint
  } else if (lineq.flag) {
    RES.all <- do.call("c", WLS.obs) - do.call("c", WLS.est)
    Delta.all <- do.call("rbind", Delta)
    if (lavmodel@eq.constraints) {
      Delta.g <- Delta.all %*% lavmodel@eq.constraints.K
    } else if (.hasSlot(lavmodel, "ceq.simple.only") &&
      lavmodel@ceq.simple.only) {
      Delta.g <- Delta.all %*% lavmodel@ceq.simple.K
    }
    Gamma.inv.weighted <- vector("list", ngroups)
    for (g in seq_len(ngroups)) {
      if (n.minus.one) {
        Ng <- nobs[[g]] - 1L
      } else {
        Ng <- nobs[[g]]
      }
      Gamma.inv.temp <- try(solve(Gamma[[g]]), silent = TRUE)
      if (inherits(Gamma.inv.temp, "try-error")) {
        # TDJ: This will happen whenever an (otherwise) unrestricted
        #      covariance matrix has a structure to it, such as equal
        #      variances (and certain covariances) for 2 members of an
        #      indistinguishable dyad (represented as 2 columns).  In
        #      such cases, their (N)ACOV elements are also identical.
        Gamma.inv.temp <- MASS::ginv(Gamma[[g]])
      }
      Gamma.inv.weighted[[g]] <- Gamma.inv.temp * Ng / ntotal
    }
    GI <- lav_matrix_bdiag(Gamma.inv.weighted)
    tDGiD <- t(Delta.g) %*% GI %*% Delta.g
    tDGiD.inv <- MASS::ginv(tDGiD) # GI may be rank-deficient
    q1 <- drop(t(RES.all) %*% GI %*% RES.all)
    q2 <- drop(t(RES.all) %*%
      GI %*% Delta.g %*% tDGiD.inv %*% t(Delta.g) %*% GI %*%
      RES.all)
    STAT <- ntotal * (q1 - q2)
    stat.group <- STAT * unlist(nobs) / ntotal # proxy only

    # 3. nonlinear equality constraints
  } else {
    # TODO
  }


  # DF
  if (!is.null(lavobject)) {
    DF <- lavobject@test[[1]]$df
  } else {
    # same approach as in lav_test.R
    df <- lav_partable_df(lavpartable)
    if (nrow(lavmodel@con.jac) > 0L) {
      ceq.idx <- attr(lavmodel@con.jac, "ceq.idx")
      if (length(ceq.idx) > 0L) {
        neq <- qr(lavmodel@con.jac[ceq.idx, , drop = FALSE])$rank
        df <- df + neq
      }
    } else if (lavmodel@ceq.simple.only) {
      # needed??
      ndat <- lav_partable_ndat(lavpartable)
      npar <- max(lavpartable$free)
      df <- ndat - npar
    }
    DF <- df
  }

  if (ADF) {
    if (model.based) {
      # using model-based Gamma
      NAME <- "browne.residual.adf.model"
      LABEL <- "Browne's residual (ADF model-based) test"
    } else {
      # regular one
      NAME <- "browne.residual.adf"
      LABEL <- "Browne's residual-based (ADF) test"
    }
  } else {
    if (model.based) {
      # using model-implied Sigma (instead of S)
      # also called the 'reweighted least-squares (RLS)' version
      NAME <- "browne.residual.nt.model"
      LABEL <- "Browne's residual (NT model-based) test"
    } else {
      # regular one
      NAME <- "browne.residual.nt"
      LABEL <- "Browne's residual-based (NT) test"
    }
  }
  out <- list(
    test = NAME,
    stat = STAT,
    stat.group = stat.group,
    df = DF,
    refdistr = "chisq",
    pvalue = 1 - pchisq(STAT, DF),
    label = LABEL
  )
  out
}
