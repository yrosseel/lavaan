# the information matrix of the unrestricted (H1) model
# taking into account:
#   - the estimator (ML or (D)WLS/ULS)
#   - missing or not
#   - fixed.x = TRUE or FALSE
#   - conditional.x = TRUE or FALSE
#   - h1.information is "structured" or "unstructured"
#
# Note: this replaces the (old) lav_model_wls_v() function
#
# - YR 22 Okt 2017: initial version
# - YR 03 Dec 2017: add lavh1, implied is either lavimplied or lavh1
#                   add support for clustered data: first.order
# - YR 03 Jan 2018: add support for clustered data: expected
# - YR 23 Aug 2018: lav_model_h1_acov (0.6-3)


## For the lavaan.mi package, TDJ provides pooled versions of all the
## sample moments called in these functions.  If any updates to these functions
## require NEW information (from @SampleStats or @implied or @h1),
## PLEASE ADD A TAG     @TDJorgensen
## at the end of the commit message on GitHub, so TDJ can check whether
## lavaan.mi::lavResiduals.mi() needs to be updated accordingly.



lav_model_h1_information <- function(lavobject = NULL,
                                     lavmodel = NULL,
                                     lavsamplestats = NULL,
                                     lavdata = NULL,
                                     lavimplied = NULL,
                                     lavh1 = NULL,
                                     lavcache = NULL,
                                     lavoptions = NULL) {
  if (!is.null(lavobject) && inherits(lavobject, "lavaan")) {
    lavmodel <- lavobject@Model
    lavsamplestats <- lavobject@SampleStats
    lavdata <- lavobject@Data
    lavimplied <- lavobject@implied
    lavh1 <- lavobject@h1
    lavcache <- lavobject@Cache
    lavoptions <- lavobject@Options
  }

  # sanity check
  if (length(lavh1) == 0L) {
    lavh1 <- lav_h1_implied_logl(
      lavdata = lavdata,
      lavsamplestats = lavsamplestats,
      lavoptions = lavoptions
    )
  }
  if (length(lavimplied) == 0L) {
    lavimplied <- lav_model_implied(lavmodel = lavmodel)
  }

  # information
  information <- lavoptions$information[1] # ALWAYS take the first one
  # the caller must control it


  # compute information matrix
  if (information == "observed") {
    i1 <- lav_model_h1_information_observed(
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats, lavdata = lavdata,
      lavimplied = lavimplied, lavh1 = lavh1,
      lavcache = lavcache, lavoptions = lavoptions
    )
  } else if (information == "expected") {
    i1 <- lav_model_h1_information_expected(
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats, lavdata = lavdata,
      lavimplied = lavimplied, lavh1 = lavh1,
      lavcache = lavcache, lavoptions = lavoptions
    )
  } else if (information == "first.order") {
    i1 <- lav_model_h1_information_firstorder(
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats, lavdata = lavdata,
      lavimplied = lavimplied, lavh1 = lavh1,
      lavcache = lavcache, lavoptions = lavoptions
    )
  }

  # I1 information, as a list per group
  i1
}

# fisher/expected information of H1
lav_model_h1_information_expected <- function(lavobject = NULL, # nolint
                                              lavmodel = NULL,
                                              lavsamplestats = NULL,
                                              lavdata = NULL,
                                              lavoptions = NULL,
                                              lavimplied = NULL,
                                              lavh1 = NULL,
                                              lavcache = NULL) {
  if (!is.null(lavobject) && inherits(lavobject, "lavaan")) {
    lavmodel <- lavobject@Model
    lavsamplestats <- lavobject@SampleStats
    lavdata <- lavobject@Data
    lavimplied <- lavobject@implied
    lavh1 <- lavobject@h1
    lavcache <- lavobject@Cache
    lavoptions <- lavobject@Options
  }

  # sanity check
  if (length(lavh1) == 0L) {
    lavh1 <- lav_h1_implied_logl(
      lavdata = lavdata,
      lavsamplestats = lavsamplestats,
      lavoptions = lavoptions
    )
  }
  if (length(lavimplied) == 0L) {
    lavimplied <- lav_model_implied(lavmodel = lavmodel)
  }

  # estimator <- lavmodel@estimator

  # structured of unstructured? (since 0.5-23)
  if (!is.null(lavoptions) &&
    !is.null(lavoptions$h1.information[1]) &&
    lavoptions$h1.information[1] == "unstructured") {
    structured <- FALSE
  } else {
    structured <- TRUE
  }

  # 1. WLS.V (=A1) for GLS/WLS
  if (lavmodel@estimator == "GLS" || lavmodel@estimator == "WLS") {
    a1 <- lavsamplestats@WLS.V

  # 1b.
  } else if (lavmodel@estimator == "DLS") {
    if (lavmodel@estimator.args$dls.GammaNT == "sample") {
      a1 <- lavsamplestats@WLS.V
    } else {
      a1 <- vector("list", length = lavsamplestats@ngroups)
      for (g in seq_len(lavsamplestats@ngroups)) {
        dls_a <- lavmodel@estimator.args$dls.a
        gamma_nt <- lav_samplestats_gamma_nt(
          m_cov          = lavimplied$cov[[g]],
          m_mean         = lavimplied$mean[[g]],
          x_idx          = lavsamplestats@x.idx[[g]],
          fixed_x        = lavmodel@fixed.x,
          conditional_x  = lavmodel@conditional.x,
          meanstructure  = lavmodel@meanstructure,
          slopestructure = lavmodel@conditional.x
        )
        w_dls <- (1 - dls_a) * lavsamplestats@NACOV[[g]] + dls_a * gamma_nt
        a1[[g]] <- lav_matrix_symmetric_inverse(w_dls)
      }
    }

  # 2. DWLS/ULS diagonal @WLS.VD slot
  } else if (lavmodel@estimator == "DWLS" || lavmodel@estimator == "ULS") {
    # diagonal only!!
    a1 <- lavsamplestats@WLS.VD

  # 3a. ML single level
  } else if (lavmodel@estimator %in% c("ML", "NTRLS", "DLS", "catML", "IV") &&
    lavdata@nlevels == 1L) {
    a1 <- vector("list", length = lavsamplestats@ngroups)

    # structured? compute model-implied statistics
    if (structured && length(lavimplied) == 0L) {
      lavimplied <- lav_model_implied(lavmodel)
    }

    for (g in 1:lavsamplestats@ngroups) {
      wt <- lavdata@weights[[g]]

      if (lavsamplestats@missing.flag) {
        # mvnorm
        # FIXME: allow for meanstructure = FALSE
        # FIXME: allow for conditional.x = TRUE
        if (lavmodel@meanstructure && structured) {
          mean_1 <- lavimplied$mean[[g]]
        } else {
          #MEAN <- lavsamplestats@missing.h1[[g]]$mu
          mean_1 <- lavh1$implied$mean[[g]]
        }

        if (structured) {
          a1[[g]] <-
            lav_mvnorm_missing_information_expected(
              y = lavdata@X[[g]],
              mp = lavdata@Mp[[g]],
              wt = wt,
              mu = mean_1,
              # meanstructure = lavmodel@meanstructure,
              sigma_1 = lavimplied$cov[[g]],
              x_idx = lavsamplestats@x.idx[[g]]
            )
        } else {
          a1[[g]] <-
            lav_mvnorm_missing_information_expected(
              y = lavdata@X[[g]],
              mp = lavdata@Mp[[g]],
              wt = wt,
              mu = mean_1,
              # meanstructure = lavmodel@meanstructure,
              #Sigma = lavsamplestats@missing.h1[[g]]$sigma,
              sigma_1 = lavh1$implied$cov[[g]],
              x_idx = lavsamplestats@x.idx[[g]]
            )
        }
      } else {
        if (lavmodel@conditional.x) {
          # mvreg
          if (lavmodel@meanstructure && structured) {
            res_int <- lavimplied$res.int[[g]]
            res_slopes <- lavimplied$res.slopes[[g]]
          } else {
            res_int <- lavsamplestats@res.int[[g]]
            res_slopes <- lavsamplestats@res.slopes[[g]]
          }

          if (structured) {
            a1[[g]] <- lav_mvreg_information_expected(
              sample.mean.x     = lavsamplestats@mean.x[[g]],
              sample.cov.x      = lavsamplestats@cov.x[[g]],
              sample.nobs       = lavsamplestats@nobs[[g]],
              res.int           = res_int,
              res.slopes        = res_slopes,
              # wt               = WT,
              # meanstructure    = lavmodel@meanstructure,
              res.cov           = lavimplied$res.cov[[g]]
            )
          } else {
            a1[[g]] <- lav_mvreg_information_expected(
              sample.mean.x     = lavsamplestats@mean.x[[g]],
              sample.cov.x      = lavsamplestats@cov.x[[g]],
              sample.nobs       = lavsamplestats@nobs[[g]],
              res.int           = lavsamplestats@res.int[[g]],
              res.slopes        = lavsamplestats@res.slopes[[g]],
              # wt               = WT,
              # meanstructure    = lavmodel@meanstructure,
              res.cov           = lavsamplestats@res.cov[[g]]
            )
          }
        } else {
          # conditional.x = FALSE
          # mvnorm
          if (lavmodel@meanstructure && structured) {
            mean_1 <- lavimplied$mean[[g]]
          } else {
            mean_1 <- lavsamplestats@mean[[g]]
          }

          correlation_flag <- lavmodel@correlation

          if (structured) {
            a1[[g]] <- lav_mvnorm_information_expected(
              sigma_1         = lavimplied$cov[[g]],
              # wt = WT, # not needed
              x_idx         = lavsamplestats@x.idx[[g]],
              meanstructure = lavmodel@meanstructure,
              correlation   = correlation_flag
            )
          } else {
            a1[[g]] <- lav_mvnorm_h1_information_expected(
              sample_cov_inv = lavsamplestats@icov[[g]],
              # wt = WT, not needed
              x_idx          = lavsamplestats@x.idx[[g]],
              meanstructure  = lavmodel@meanstructure,
              correlation    = correlation_flag
            )
          }
        } # conditional.x
      } # missing

      # stochastic group weight
      if (lavmodel@group.w.free) {
        # unweight!! (as otherwise, we would 'weight' again)
        a <- exp(lavimplied$group.w[[g]]) / lavsamplestats@nobs[[g]]
        a1[[g]] <- lav_matrix_bdiag(matrix(a, 1L, 1L), a1[[g]])
      }
    } # g
  # ML

  # 3b. ML + multilevel
  } else if (lavmodel@estimator == "ML" && lavdata@nlevels > 1L) {
    a1 <- vector("list", length = lavsamplestats@ngroups)

    # structured? compute model-implied statistics
    if (structured && length(lavimplied) == 0L) {
      lavimplied <- lav_model_implied(lavmodel)
    }

    # structured? lavimplied vs lavh1
    if (structured) {
      implied <- lavimplied
    } else {
      implied <- lavh1$implied
    }

    for (g in 1:lavsamplestats@ngroups) {
      mu_w <- implied$mean[[(g - 1) * lavdata@nlevels + 1L]]
      mu_b <- implied$mean[[(g - 1) * lavdata@nlevels + 2L]]
      sigma_w <- implied$cov[[(g - 1) * lavdata@nlevels + 1L]]
      sigma_b <- implied$cov[[(g - 1) * lavdata@nlevels + 2L]]

      # clustered data
      a1[[g]] <- lav_mvnorm_cluster_information_expected(
        lp           = lavdata@Lp[[g]],
        mu_w         = mu_w,
        sigma_w      = sigma_w,
        mu_b         = mu_b,
        sigma_b      = sigma_b,
        x_idx        = lavsamplestats@x.idx[[g]]
      )
    } # g
  } # ML + multilevel


  a1
}

lav_model_h1_information_observed <- function(lavobject = NULL, # nolint
                                              lavmodel = NULL,
                                              lavsamplestats = NULL,
                                              lavdata = NULL,
                                              lavimplied = NULL,
                                              lavh1 = NULL,
                                              lavcache = NULL,
                                              lavoptions = NULL) {
  if (!is.null(lavobject) && inherits(lavobject, "lavaan")) {
    lavmodel <- lavobject@Model
    lavsamplestats <- lavobject@SampleStats
    lavdata <- lavobject@Data
    lavimplied <- lavobject@implied
    lavh1 <- lavobject@h1
    lavcache <- lavobject@Cache
    lavoptions <- lavobject@Options
  }

  # sanity check
  if (length(lavh1) == 0L) {
    lavh1 <- lav_h1_implied_logl(
      lavdata = lavdata,
      lavsamplestats = lavsamplestats,
      lavoptions = lavoptions
    )
  }
  if (length(lavimplied) == 0L) {
    lavimplied <- lav_model_implied(lavmodel = lavmodel)
  }

  # estimator <- lavmodel@estimator

  # structured?
  if (!is.null(lavoptions) &&
    !is.null(lavoptions$h1.information[1]) &&
    lavoptions$h1.information[1] == "unstructured") {
    structured <- FALSE
  } else {
    structured <- TRUE
  }

  # 1. WLS.V (=A1) for GLS/WLS
  if (lavmodel@estimator == "GLS" || lavmodel@estimator == "WLS" ||
    lavmodel@estimator == "DLS") {
    a1 <- lavsamplestats@WLS.V
  # 1b.
  } else if (lavmodel@estimator == "DLS") {
    if (lavmodel@estimator.args$dls.GammaNT == "sample") {
      a1 <- lavsamplestats@WLS.V
    } else {
      a1 <- vector("list", length = lavsamplestats@ngroups)
      for (g in seq_len(lavsamplestats@ngroups)) {
        dls_a <- lavmodel@estimator.args$dls.a
        gamma_nt <- lav_samplestats_gamma_nt(
          m_cov          = lavimplied$cov[[g]],
          m_mean         = lavimplied$mean[[g]],
          x_idx          = lavsamplestats@x.idx[[g]],
          fixed_x        = lavmodel@fixed.x,
          conditional_x  = lavmodel@conditional.x,
          meanstructure  = lavmodel@meanstructure,
          slopestructure = lavmodel@conditional.x
        )
        w_dls <- (1 - dls_a) * lavsamplestats@NACOV[[g]] + dls_a * gamma_nt
        a1[[g]] <- lav_matrix_symmetric_inverse(w_dls)
      }
    }
  # 2. DWLS/ULS diagonal @WLS.VD slot
  } else if (lavmodel@estimator == "DWLS" || lavmodel@estimator == "ULS") {
    # diagonal only!!
    a1 <- lavsamplestats@WLS.VD
  # 3a. ML single level
  } else if (lavmodel@estimator == "ML" && lavdata@nlevels == 1L) {
    a1 <- vector("list", length = lavsamplestats@ngroups)

    # structured? compute model-implied statistics
    if (structured && length(lavimplied) == 0L) {
      lavimplied <- lav_model_implied(lavmodel)
    }

    for (g in 1:lavsamplestats@ngroups) {
      if (lavsamplestats@missing.flag) {
        # mvnorm
        # FIXME: allow for meanstructure = FALSE
        # FIXME: allow for conditional.x = TRUE
        if (lavmodel@meanstructure && structured) {
          mean_1 <- lavimplied$mean[[g]]
        } else {
          #MEAN <- lavsamplestats@missing.h1[[g]]$mu
          mean_1 <- lavh1$implied$mean[[g]]
        }

        if (structured) {
          a1[[g]] <-
            lav_mvnorm_missing_information_observed_samplestats(
              yp = lavsamplestats@missing[[g]],
              # wt not needed
              mu = mean_1,
              # meanstructure = lavmodel@meanstructure,
              sigma_1 = lavimplied$cov[[g]],
              x_idx = lavsamplestats@x.idx[[g]]
            )
        } else {
          a1[[g]] <-
            lav_mvnorm_missing_information_observed_samplestats(
              yp = lavsamplestats@missing[[g]],
              # wt not needed
              mu = mean_1,
              # meanstructure = lavmodel@meanstructure,
              #Sigma = lavsamplestats@missing.h1[[g]]$sigma,
              sigma_1 = lavh1$implied$cov[[g]],
              x_idx = lavsamplestats@x.idx[[g]]
            )
        }
      } else {
        if (lavmodel@conditional.x) {
          # mvreg
          if (lavmodel@meanstructure && structured) {
            res_int <- lavimplied$res.int[[g]]
            res_slopes <- lavimplied$res.slopes[[g]]
          } else {
            res_int <- lavsamplestats@res.int[[g]]
            res_slopes <- lavsamplestats@res.slopes[[g]]
          }

          if (structured) {
            a1[[g]] <- lav_mvreg_information_observed_samplestats(
              sample.res.int    = lavsamplestats@res.int[[g]],
              sample.res.slopes = lavsamplestats@res.slopes[[g]],
              sample.res.cov    = lavsamplestats@res.cov[[g]],
              sample.mean.x     = lavsamplestats@mean.x[[g]],
              sample.cov.x      = lavsamplestats@cov.x[[g]],
              res.int           = res_int,
              res.slopes        = res_slopes,
              # wt               = WT,
              # meanstructure    = lavmodel@meanstructure,
              res.cov           = lavimplied$res.cov[[g]]
            )
          } else {
            a1[[g]] <- lav_mvreg_information_observed_samplestats(
              sample.res.int    = lavsamplestats@res.int[[g]],
              sample.res.slopes = lavsamplestats@res.slopes[[g]],
              sample.res.cov    = lavsamplestats@res.cov[[g]],
              sample.mean.x     = lavsamplestats@mean.x[[g]],
              sample.cov.x      = lavsamplestats@cov.x[[g]],
              res.int           = lavsamplestats@res.int[[g]],
              res.slopes        = lavsamplestats@res.slopes[[g]],
              # wt               = WT,
              # meanstructure    = lavmodel@meanstructure,
              res.cov           = lavsamplestats@res.cov[[g]]
            )
          }
        } else {
          # conditional.x = FALSE
          # mvnorm
          if (lavmodel@meanstructure && structured) {
            mean_1 <- lavimplied$mean[[g]]
          } else {
            mean_1 <- lavsamplestats@mean[[g]]
          }

          if (structured) {
            a1[[g]] <- lav_mvnorm_information_observed_samplestats(
              sample_mean   = lavsamplestats@mean[[g]],
              sample_cov    = lavsamplestats@cov[[g]],
              mu            = mean_1,
              sigma_1       = lavimplied$cov[[g]],
              # wt           = WT, # not needed
              x_idx         = lavsamplestats@x.idx[[g]],
              meanstructure = lavmodel@meanstructure
            )
          } else {
            a1[[g]] <- lav_mvnorm_h1_information_observed_samplestats(
              sample_mean    = lavsamplestats@mean[[g]],
              sample_cov     = lavsamplestats@cov[[g]],
              sample_cov_inv = lavsamplestats@icov[[g]],
              # wt            = WT, not needed
              x_idx          = lavsamplestats@x.idx[[g]],
              meanstructure  = lavmodel@meanstructure
            )
          }
        } # conditional.x
      } # missing

      # stochastic group weight
      if (lavmodel@group.w.free) {
        # unweight!!
        a <- exp(lavimplied$group.w[[g]]) / lavsamplestats@nobs[[g]]
        a1[[g]] <- lav_matrix_bdiag(matrix(a, 1, 1), a1[[g]])
      }
    } # g
  # ML
  # 3b. ML + multilevel
  } else if (lavmodel@estimator == "ML" && lavdata@nlevels > 1L) {
    a1 <- vector("list", length = lavsamplestats@ngroups)

    # structured? compute model-implied statistics
    if (structured && length(lavimplied) == 0L) {
      lavimplied <- lav_model_implied(lavmodel)
    }

    # structured? lavimplied vs lavh1
    if (structured) {
      implied <- lavimplied
    } else {
      implied <- lavh1$implied
    }

    for (g in 1:lavsamplestats@ngroups) {
      mu_w <- implied$mean[[(g - 1) * lavdata@nlevels + 1L]]
      mu_b <- implied$mean[[(g - 1) * lavdata@nlevels + 2L]]
      sigma_w <- implied$cov[[(g - 1) * lavdata@nlevels + 1L]]
      sigma_b <- implied$cov[[(g - 1) * lavdata@nlevels + 2L]]

      if (lavdata@missing == "ml") {
          a1[[g]] <- lav_mvnorm_cluster_missing_information_observed(
            y1            = lavdata@X[[g]],
            y2            = lavsamplestats@YLp[[g]][[2]]$Y2,
            lp            = lavdata@Lp[[g]],
            mp            = lavdata@Mp[[g]],
            mu_w          = mu_w,
            sigma_w       = sigma_w,
            mu_b          = mu_b,
            sigma_b       = sigma_b,
            x_idx         = lavsamplestats@x.idx[[g]]
          )
      } else {
        a1[[g]] <- lav_mvnorm_cluster_information_observed(
          lp           = lavdata@Lp[[g]],
          ylp          = lavsamplestats@YLp[[g]],
          mu_w         = mu_w,
          sigma_w      = sigma_w,
          mu_b         = mu_b,
          sigma_b      = sigma_b,
          x_idx        = lavsamplestats@x.idx[[g]]
        )
      }
    } # g
  } # ML + multilevel

  a1
}

# outer product of the case-wise scores (gradients)
# HJ 18/10/2023: Adjust J matrix correctly using weights. Note: H matrix is
# based on lav_model_hessian so no changes required.
lav_model_h1_information_firstorder <- function(lavobject = NULL,     # nolint
                                                lavmodel = NULL,
                                                lavsamplestats = NULL,
                                                lavdata = NULL,
                                                lavimplied = NULL,
                                                lavh1 = NULL,
                                                lavcache = NULL,
                                                lavoptions = NULL) {
  if (!is.null(lavobject) && inherits(lavobject, "lavaan")) {
    lavmodel <- lavobject@Model
    lavsamplestats <- lavobject@SampleStats
    lavdata <- lavobject@Data
    lavimplied <- lavobject@implied
    lavh1 <- lavobject@h1
    lavcache <- lavobject@Cache
    lavoptions <- lavobject@Options
  }

  # sanity check
  if (length(lavh1) == 0L) {
    lavh1 <- lav_h1_implied_logl(
      lavdata = lavdata,
      lavsamplestats = lavsamplestats,
      lavoptions = lavoptions
    )
  }
  if (length(lavimplied) == 0L) {
    lavimplied <- lav_model_implied(lavmodel = lavmodel)
  }

  estimator <- lavmodel@estimator
  if (!estimator %in% c("ML", "PML")) {
    lav_msg_stop(gettext(
      "information = \"first.order\" not available for estimator"),
      sQuote(estimator))
  }

  # structured?
  if (!is.null(lavoptions) &&
    !is.null(lavoptions$h1.information[1]) &&
    lavoptions$h1.information[1] == "unstructured") {
    structured <- FALSE
  } else {
    structured <- TRUE
  }

  # clustered?
  if (!is.null(lavoptions) &&
    !is.null(lavoptions$.clustered) &&
    lavoptions$.clustered) {
    clustered <- TRUE
    if (is.null(lavdata@Lp[[1]])) {
      lav_msg_stop(gettext("lavdata@Lp is empty, while clustered = TRUE"))
    }
    # if (estimator == "PML") {
    #   lav_msg_stop(gettext(
    #     "clustered information is not (yet) available when estimator = 'PML'"))
    # }
    # if(lavsamplestats@missing.flag) {
    #    stop("lavaan ERROR: clustered information is not (yet) available
    #                               when missing = \"ML\"")
    # }
    # if(lavmodel@conditional.x) {
    #    stop("lavaan ERROR: clustered information is not (yet) available
    #                               when conditional.x = TRUE")
    # }
    # if(!structured) {
    #    stop("lavaan ERROR: clustered information is not (yet) available
    #                               when h1.information = \"unstructured\"")
    # }
  } else {
    clustered <- FALSE
  }

  # structured? compute model-implied statistics
  if (estimator == "PML" || structured) {
    if (length(lavimplied) == 0L) {
      lavimplied <- lav_model_implied(lavmodel)
    }
  }

  # structured? lavimplied vs lavh1
  if (structured) {
    implied <- lavimplied
  } else {
    implied <- lavh1$implied
  }

  b1 <- vector("list", length = lavsamplestats@ngroups)
  for (g in 1:lavdata@ngroups) {
    wt <- lavdata@weights[[g]]

    if (estimator == "PML") {
      # slow approach: compute outer product of case-wise scores

      if (lavmodel@conditional.x) {
        sigma_1 <- implied$res.cov[[g]]
        mu <- implied$res.mean[[g]]
        th <- implied$res.th[[g]]
        pi0 <- implied$res.slopes[[g]]
        exo <- lavdata@eXo[[g]]
      } else {
        sigma_1 <- implied$cov[[g]]
        mu <- implied$mean[[g]]
        th <- implied$th[[g]]
        pi0 <- NULL
        exo <- NULL
      }
      sc <- lav_pml_dploglik_dimplied(
        sigma_hat = sigma_1,
        mu_hat = mu,
        th = th,
        th_idx = lavmodel@th.idx[[g]],
        num_idx = lavmodel@num.idx[[g]],
        x = lavdata@X[[g]],
        exo = exo,
        wt = NULL,
        pi0 = pi0,
        lavcache = lavcache[[g]],
        missing = lavdata@missing,
        scores = TRUE,
        negative = FALSE
      )

      # >>>>>>>> HJ/MK PML CODE >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


      # information H1

      if (isTRUE(clustered)) {
        # For clustered data, need to compute (centred) crossprod within each
        # cluster and sum them all up.
        clusters     <- lavdata@Lp[[g]]$cluster.id[[2]]  # why list of 2?
        clusters_idx <- lavdata@Lp[[g]]$cluster.idx[[2]]
        nclust       <- length(clusters)
        zb <- list()

        if (is.null(wt)) wt <- rep(1, length(clusters_idx))

        for (b in seq_along(clusters)) {
          sc_b <- sc[clusters_idx == b, ]
          wt_b <- wt[clusters_idx == b]
          zb[[b]] <- apply(sc_b * wt_b, 2, sum)
        }
        zbar <- apply(do.call(cbind, zb), 1, mean)
        b1c <- Reduce(f = `+`, lapply(zb, function(z) tcrossprod(z - zbar)))
        b1[[g]] <- nclust / (nclust - 1) * b1c

      } else {
        if (is.null(wt)) {
          b1[[g]] <- lav_matrix_crossprod(sc)
        } else {
          b1[[g]] <- crossprod(wt * sc)
        }
      }


    # ML
    } else if (estimator == "ML" && lavdata@nlevels > 1L) {
      # if not-structured, we use lavh1, and that is always
      # 'unconditional' (for now)
      if (lavmodel@conditional.x && structured) {
      if (lavdata@missing == "ml") {
      lav_msg_stop(gettext("firstorder information matrix not available
                                (yet) if conditional.x + fiml"))
    }
        res_sigma_w <- implied$res.cov[[(g - 1) * lavdata@nlevels + 1L]]
        res_int_w <- implied$res.int[[(g - 1) * lavdata@nlevels + 1L]]
        res_pi_w <- implied$res.slopes[[(g - 1) * lavdata@nlevels + 1L]]

        res_sigma_b <- implied$res.cov[[(g - 1) * lavdata@nlevels + 2L]]
        res_int_b <- implied$res.int[[(g - 1) * lavdata@nlevels + 2L]]
        res_pi_b <- implied$res.slopes[[(g - 1) * lavdata@nlevels + 2L]]
        b1[[g]] <- lav_mvreg_cluster_information_firstorder(
          Y1            = lavdata@X[[g]],
          YLp           = lavsamplestats@YLp[[g]],
          Lp            = lavdata@Lp[[g]],
          Res.Sigma.W   = res_sigma_w,
          Res.Int.W     = res_int_w,
          Res.Pi.W      = res_pi_w,
          Res.Sigma.B   = res_sigma_b,
          Res.Int.B     = res_int_b,
          Res.Pi.B      = res_pi_b,
          divide.by.two = TRUE
        )
      } else {
        mu_w <- implied$mean[[(g - 1) * lavdata@nlevels + 1L]]
        mu_b <- implied$mean[[(g - 1) * lavdata@nlevels + 2L]]
        sigma_w <- implied$cov[[(g - 1) * lavdata@nlevels + 1L]]
        sigma_b <- implied$cov[[(g - 1) * lavdata@nlevels + 2L]]
    if (lavdata@missing == "ml") {
          b1[[g]] <- lav_mvnorm_cluster_missing_information_firstorder(
            y1            = lavdata@X[[g]],
            y2            = lavsamplestats@YLp[[g]][[2]]$Y2,
            lp            = lavdata@Lp[[g]],
            mp            = lavdata@Mp[[g]],
            mu_w          = mu_w,
            sigma_w       = sigma_w,
            mu_b          = mu_b,
            sigma_b       = sigma_b,
            x_idx         = lavsamplestats@x.idx[[g]],
            divide_by_two = TRUE
          )
    } else {
      # no missing values
          b1[[g]] <- lav_mvnorm_cluster_information_firstorder(
            y1            = lavdata@X[[g]],
            ylp           = lavsamplestats@YLp[[g]],
            lp            = lavdata@Lp[[g]],
            mu_w          = mu_w,
            sigma_w       = sigma_w,
            mu_b          = mu_b,
            sigma_b       = sigma_b,
            x_idx         = lavsamplestats@x.idx[[g]],
            divide_by_two = TRUE
          )
    }
      }
    } else if (estimator == "ML" && lavdata@nlevels == 1L) {
      if (length(lavdata@cluster) > 0L) {
        cluster_idx <- lavdata@Lp[[g]]$cluster.idx[[2]]
      } else {
        cluster_idx <- NULL
      }

      if (lavsamplestats@missing.flag) {
        # mvnorm
        # FIXME: allow for meanstructure = FALSE
        # FIXME: allow for conditional.x = TRUE
        if (lavmodel@meanstructure && structured) {
          mean_1 <- lavimplied$mean[[g]]
        } else {
          #MEAN <- lavsamplestats@missing.h1[[g]]$mu
          mean_1 <- lavh1$implied$mean[[g]]
        }

        b1[[g]] <- lav_mvnorm_missing_information_firstorder(
          y = lavdata@X[[g]],
          mp = lavdata@Mp[[g]], wt = wt,
          cluster_idx = cluster_idx,
          mu = mean_1,
          # meanstructure = lavmodel@meanstructure,
          sigma_1 = implied$cov[[g]],
          x_idx = lavsamplestats@x.idx[[g]]
        )
      } else {
        if (lavmodel@conditional.x) {
          # mvreg
          if (lavmodel@meanstructure && structured) {
            res_int <- lavimplied$res.int[[g]]
            res_slopes <- lavimplied$res.slopes[[g]]
          } else {
            res_int <- lavsamplestats@res.int[[g]]
            res_slopes <- lavsamplestats@res.slopes[[g]]
          }

          b1[[g]] <- lav_mvreg_information_firstorder(
            Y              = lavdata@X[[g]],
            eXo            = lavdata@eXo[[g]],
            res.int        = res_int,
            res.slopes     = res_slopes,
            # wt            = WT,
            # meanstructure = lavmodel@meanstructure,
            res.cov        = implied$res.cov[[g]]
          )
        } else {
          # conditional.x = FALSE
          # mvnorm
          if (lavmodel@meanstructure && structured) {
            mean_1 <- lavimplied$mean[[g]]
          } else {
            # NOTE: the information matrix will be the same (minus
            # the meanstructure block), but once INVERTED, the
            # standard errors will be (slightly) smaller!!!
            # This is only visible when estimator = "MLF"
            # (or information = "first.order")
            mean_1 <- lavsamplestats@mean[[g]] # saturated
          }

          if (structured) {
            b1[[g]] <- lav_mvnorm_information_firstorder(
              y = lavdata@X[[g]],
              mu = mean_1, sigma_1 = lavimplied$cov[[g]],
              wt = wt,
              cluster_idx = cluster_idx,
              x_idx = lavsamplestats@x.idx[[g]],
              meanstructure = lavmodel@meanstructure
            )
          } else {
            b1[[g]] <- lav_mvnorm_h1_information_firstorder(
              y = lavdata@X[[g]],
              sample_cov_inv = lavsamplestats@icov[[g]],
              gamma_1 = lavsamplestats@NACOV[[g]],
              wt = wt,
              cluster_idx = cluster_idx, # only if wt
              x_idx = lavsamplestats@x.idx[[g]],
              meanstructure = lavmodel@meanstructure
            )
          }
        } # mvnorm
      } # missing
    } # ML

    # stochastic group weight
    if (lavmodel@group.w.free) {
      # unweight!!
      a <- exp(lavimplied$group.w[[g]]) / lavsamplestats@nobs[[g]]
      b1[[g]] <- lav_matrix_bdiag(matrix(a, 1, 1), b1[[g]])
    }
  } # g

  b1
}


# asymptotic variance matrix (=Gamma/N) of the unrestricted (H1)
# sample statistics
#
# FIXME: make this work for categorical/GLS/WLS/...
#
lav_model_h1_acov <- function(lavobject = NULL,
                              lavmodel = NULL,
                              lavsamplestats = NULL,
                              lavdata = NULL,
                              lavoptions = NULL,
                              lavimplied = NULL,
                              lavh1 = NULL,
                              lavcache = NULL,
                              meanstructure = NULL, # if specified, use it
                              h1_information = NULL, # if specified, use it
                              se = NULL) { # if specified, use it

  if (!is.null(lavobject) && inherits(lavobject, "lavaan")) {
    lavmodel <- lavobject@Model
    lavsamplestats <- lavobject@SampleStats
    lavdata <- lavobject@Data
    lavimplied <- lavobject@implied
    lavh1 <- lavobject@h1
    lavcache <- lavobject@Cache
    lavoptions <- lavobject@Options
  }

  # sanity check
  if (length(lavh1) == 0L) {
    lavh1 <- lav_h1_implied_logl(
      lavdata = lavdata,
      lavsamplestats = lavsamplestats,
      lavoptions = lavoptions
    )
  }
  if (length(lavimplied) == 0L) {
    lavimplied <- lav_model_implied(lavmodel = lavmodel)
  }

  # override
  if (!is.null(meanstructure)) {
    lavoptions$meanstructure <- meanstructure
  }
  if (!is.null(h1_information)) {
    lavoptions$h1.information[1] <- h1_information
  }
  if (!is.null(se)) {
    lavoptions$se <- se
  }


  # information
  information <- lavoptions$information[1] # ALWAYS used the first

  # compute information matrix
  if (information == "observed") {
    i1 <- lav_model_h1_information_observed(
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats, lavdata = lavdata,
      lavimplied = lavimplied, lavh1 = lavh1,
      lavcache = lavcache, lavoptions = lavoptions
    )
  } else if (information == "expected") {
    i1 <- lav_model_h1_information_expected(
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats, lavdata = lavdata,
      lavimplied = lavimplied, lavh1 = lavh1,
      lavcache = lavcache, lavoptions = lavoptions
    )
  } else if (information == "first.order") {
    i1 <- lav_model_h1_information_firstorder(
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats, lavdata = lavdata,
      lavimplied = lavimplied, lavh1 = lavh1,
      lavcache = lavcache, lavoptions = lavoptions
    )
  }

  if (lavoptions$se %in% c("robust.huber.white", "robust.sem")) {
    j1 <- lav_model_h1_information_firstorder(
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats, lavdata = lavdata,
      lavimplied = lavimplied, lavh1 = lavh1,
      lavcache = lavcache, lavoptions = lavoptions
    )
  }

  # compute ACOV per group
  acov <- vector("list", length = lavdata@ngroups)
  for (g in 1:lavdata@ngroups) {
    # denominator
    if (lavdata@nlevels == 1L) {
      ng <- lavsamplestats@nobs[[g]]
    } else {
      ng <- lavdata@Lp[[g]]$nclusters[[2]]
    }

    # invert information
    i1_g_inv <- try(lav_matrix_symmetric_inverse(i1[[g]]), silent = TRUE)
    if (inherits(i1_g_inv, "try-error")) {
      lav_msg_stop(gettext(
        "could not invert h1 information matrix in group"), g)
    }

    # which type of se?
    if (lavoptions$se %in% c("standard", "none")) {
      acov[[g]] <- 1 / ng * i1_g_inv
    } else if (lavoptions$se %in% c("robust.huber.white", "robust.sem")) {
      acov[[g]] <- 1 / ng * (i1_g_inv %*% j1[[g]] %*% i1_g_inv)
    }
  }

  acov
}
