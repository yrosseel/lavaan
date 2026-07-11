# Browne's residual test statistic
# see Browne (1984) eq 2.20a

# T_B = (N-1) * t(RES) %*% Delta.c %*%
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
                            n_minus_one = "default",
                            adf = TRUE,
                            model_based = FALSE) {
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

  if (lavmodel@categorical) {
    if (!adf) {
      lav_msg_warn(gettext("normal theory version of Browne's residual test
        not available in the categorical setting; switching to ADF."))
      adf <- TRUE
    }
    if (model_based) {
      lav_msg_stop(gettext("model-based version of Browne's residual test not
        available in the categorical setting; switching to sample-based."))
      model_based <- FALSE
    }
  }
  if (lavdata@missing != "listwise" && !model_based) {
    lav_msg_stop(gettext("sample-based version of Browne's test is not
      available when data is missing; switching to model-based."))
    model_based <- TRUE
  }
  if (lavdata@nlevels > 1L) {
    lav_msg_stop(gettext("Browne's test is not available when data is
                         multilevel."))
  }
  if (length(lavmodel@ceq.nonlinear.idx) > 0L) {
    lav_msg_stop(gettext("Browne's test is not available (yet) when nonlinear
                         equality constraints are involved."))
  }

  if (!is.logical(n_minus_one)) {
    if (lavoptions$estimator == "ML" &&
      lavoptions$likelihood == "normal") {
      n_minus_one <- FALSE
    } else {
      n_minus_one <- TRUE
    }
  }

  # linear equality constraints? NOTE: do not test the packing flags
  # (@eq.constraints / @ceq.simple.only) here -- they are both FALSE when
  # equality constraints coexist with inequality constraints or bounds,
  # and the constraints would silently be dropped from the projection
  eq_basis <- lav_con_eq_basis(lavmodel)
  lineq_flag <- !is.null(eq_basis)

  # can we use the fast version?
  fast_flag <- FALSE
  if (!adf && !lineq_flag && !lavmodel@conditional.x &&
      !lavmodel@group.w.free) {
    fast_flag <- TRUE
    if (model_based) {
      implied <- lavimplied
    } else {
      implied <- lavh1$implied
    }
  }

  # ingredients
  delta <- lav_model_delta(lavmodel)
  if (adf) {
    # ADF version
    # note: the stored NACOV can only be used for the sample-based
    # version; the model-based version must recompute Gamma at the
    # model-implied moments (before 0.7-2, browne.residual.adf.model
    # silently returned the sample-based statistic whenever a NACOV
    # was stored, e.g. for the (D)WLS-family estimators)
    if (!model_based &&
        !is.null(lavsamplestats@NACOV[[1]]) &&
        !lavoptions$se == "robust.sem.nt") { # eg, estimator = "ULS"/"DWLS"
      gamma_1 <- lavsamplestats@NACOV
    } else {
      if (!is.null(lavobject)) {
        if (lavobject@Data@data.type != "full") {
          lav_msg_stop(gettext("ADF version not available without full data or
                               user-provided Gamma/NACOV matrix"))
        }
        gamma_1 <- lav_object_gamma(lavobject,
          adf = TRUE,
          model_based = model_based
        )
      } else {
        if (lavdata@data.type != "full") {
          lav_msg_stop(gettext("ADF version not available without full data or
                               user-provided Gamma/NACOV matrix"))
        }
        gamma_1 <- lav_object_gamma(
          lavobject = NULL,
          lavdata = lavdata,
          lavoptions = lavoptions,
          lavsamplestats = lavsamplestats,
          lavh1 = lavh1,
          lavimplied = lavimplied,
          adf = TRUE,
          model_based = model_based
        )
      }
    }
  } else if (!fast_flag) {
    # NT version
    if (!is.null(lavobject)) {
      gamma_1 <- lav_object_gamma(lavobject,
        adf = FALSE,
        model_based = model_based
      )
    } else {
      gamma_1 <- lav_object_gamma(
        lavobject = NULL,
        lavdata = lavdata,
        lavoptions = lavoptions,
        lavsamplestats = lavsamplestats,
        lavh1 = lavh1,
        lavimplied = lavimplied,
        adf = FALSE,
        model_based = model_based
      )
    }
  }
  wls_obs <- lavsamplestats@WLS.obs
  wls_est <- lav_model_wls_est(lavmodel)
  nobs <- lavsamplestats@nobs
  ntotal <- lavsamplestats@ntotal

  # compute T.B per group
  ngroups <- length(wls_obs)
  stat_group <- numeric(ngroups)

  # 1. standard setting: no equality constraints
  if (!lineq_flag) {
    for (g in seq_len(ngroups)) {
      res <- wls_obs[[g]] - wls_est[[g]]
      delta_g <- delta[[g]]
      if (n_minus_one) {
        ng <- nobs[[g]] - 1L
      } else {
        ng <- nobs[[g]]
      }
      if (fast_flag) {
        stat_group[g] <- lav_test_browne_nt_fast(
          res = res, delta = delta_g,
          sample_cov = implied$cov[[g]], sample_nobs = ng,
          meanstructure = lavmodel@meanstructure,
          x_idx = lavsamplestats@x.idx[[g]]
        )
      } else {
        # naive formula base computation (slow!)
        delta_c <- lav_mat_ortho_complement(delta_g)
        t_dgd <- crossprod(delta_c, gamma_1[[g]]) %*% delta_c
        # if fixed.x = TRUE, gamma_1[[g]] may contain zero col/rows
        t_dgd_inv <- lav_mat_sym_inverse(t_dgd)
        t_res_delta_c <- crossprod(res, delta_c)
        stat_group[g] <-
          ng * drop(t_res_delta_c %*% t_dgd_inv %*% t(t_res_delta_c))
      }
    }
    stat <- sum(stat_group)

    # 2. linear equality constraint
  } else if (lineq_flag) {
    res_all <- do.call("c", wls_obs) - do.call("c", wls_est)
    delta_all <- do.call("rbind", delta)
    delta_g <- delta_all %*% eq_basis
    gamma_inv_weighted <- vector("list", ngroups)
    for (g in seq_len(ngroups)) {
      if (n_minus_one) {
        ng <- nobs[[g]] - 1L
      } else {
        ng <- nobs[[g]]
      }
      gamma_inv_temp <- try(solve(gamma_1[[g]]), silent = TRUE)
      if (inherits(gamma_inv_temp, "try-error")) {
        # TDJ: This will happen whenever an (otherwise) unrestricted
        #      covariance matrix has a structure to it, such as equal
        #      variances (and certain covariances) for 2 members of an
        #      indistinguishable dyad (represented as 2 columns).  In
        #      such cases, their (N)ACOV elements are also identical.
        gamma_inv_temp <- MASS::ginv(gamma_1[[g]])
      }
      gamma_inv_weighted[[g]] <- gamma_inv_temp * ng / ntotal
    }
    gi <- lav_mat_bdiag(gamma_inv_weighted)
    t_dgi_d <- t(delta_g) %*% gi %*% delta_g
    t_dgi_d_inv <- MASS::ginv(t_dgi_d) # GI may be rank-deficient
    q1 <- drop(t(res_all) %*% gi %*% res_all)
    q2 <- drop(t(res_all) %*%
      gi %*% delta_g %*% t_dgi_d_inv %*% t(delta_g) %*% gi %*%
      res_all)
    stat <- ntotal * (q1 - q2)
    stat_group <- stat * unlist(nobs) / ntotal # proxy only
  }
  # (nonlinear equality constraints were caught at the top)


  # DF
  if (!is.null(lavobject)) {
    df_1 <- lavobject@test[[1]]$df
  } else {
    # same approach as in lav_test.R
    df_1 <- lav_test_df(lavpartable = lavpartable, lavmodel = lavmodel)
  }

  if (adf) {
    if (model_based) {
      # using model-based Gamma
      name <- "browne.residual.adf.model"
      label <- "Browne's residual (ADF model-based) test"
    } else {
      # regular one
      name <- "browne.residual.adf"
      label <- "Browne's residual-based (ADF) test"
    }
  } else {
    if (model_based) {
      # using model-implied Sigma (instead of S)
      # also called the 'reweighted least-squares (RLS)' version
      name <- "browne.residual.nt.model"
      label <- "Browne's residual (NT model-based) test"
    } else {
      # regular one
      name <- "browne.residual.nt"
      label <- "Browne's residual-based (NT) test"
    }
  }
  out <- list(
    test = name,
    stat = stat,
    stat.group = stat_group,
    df = df_1,
    refdistr = "chisq",
    pvalue = 1 - pchisq(stat, df_1),
    label = label
  )
  out
}


# faster version for the NT setting with no constraints
lav_test_browne_nt_fast <- function(res = NULL, delta = NULL,
                                    sample_cov = NULL, sample_nobs = NULL,
                                    meanstructure = FALSE,
                                    x_idx = integer(0L)) {
  nvar <- nrow(sample_cov)
  pstar <- nvar * (nvar + 1L) / 2L
  q <- ncol(delta)
  nrow_expected <- if (meanstructure) nvar + pstar else pstar
  stopifnot(length(res) == nrow_expected, nrow(delta) == nrow_expected)

  # only once
  s_inv <- solve(sample_cov)

  # handle fixed.x = TRUE
  if (length(x_idx) > 0L) {
    s_inv[x_idx, ] <- 0.0
    s_inv[, x_idx] <- 0.0
  }

  # vech
  diag_idx <- lav_mat_diagh_idx(nvar)

  # applies gamma_cov_inv = 0.5 * t(D) %*% (S.inv %x% S.inv) %*% D
  # to a pstar-vector 'x'
  apply_gamma_inv <- function(x) {
    x_cov <- x
    if (meanstructure) {
      x_mean <- x[seq_len(nvar)]
      x_cov <- x[nvar + seq_len(pstar)]
      out_mean <- drop(s_inv %*% x_mean)
    }

    m_w <- lav_mat_vech_rev(x_cov)
    z <- s_inv %*% m_w %*% s_inv
    out_cov <- lav_mat_vech(z)
    out_cov[diag_idx] <- out_cov[diag_idx] / 2.0
    out <- out_cov
    if (meanstructure) {
      out <- c(out_mean, out_cov)
    }
    out
  }

  u <- apply_gamma_inv(res)
  gamma_inv_delta <- matrix(0.0, nrow_expected, q)
  for (j in seq_len(q)) {
    gamma_inv_delta[, j] <- apply_gamma_inv(delta[, j])
  }
  a <- crossprod(delta, gamma_inv_delta)
  b <- crossprod(delta, u)
  ab_inv_b <- tryCatch(
    {
      r_chol <- chol(a)
      backsolve(r_chol, forwardsolve(t(r_chol), b))
    },
    error = function(e) {
      MASS::ginv(a) %*% b
    }
  )

  term1 <- as.numeric(crossprod(res, u)) # t(res) Ginv res
  term2 <- as.numeric(crossprod(b, ab_inv_b)) # t(b) A^{-1} b

  stat <- sample_nobs * (term1 - term2)
  stat
}
