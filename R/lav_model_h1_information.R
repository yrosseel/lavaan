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
    if (.hasSlot(lavobject, "h1")) {
      lavh1 <- lavobject@h1
    } else {
      lavh1 <- lav_h1_implied_logl(
        lavdata = lavobject@Data,
        lavsamplestats = lavobject@SampleStats,
        lavoptions = lavobject@Options
      )
    }
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
    I1 <- lav_model_h1_information_observed(
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats, lavdata = lavdata,
      lavimplied = lavimplied, lavh1 = lavh1,
      lavcache = lavcache, lavoptions = lavoptions
    )
  } else if (information == "expected") {
    I1 <- lav_model_h1_information_expected(
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats, lavdata = lavdata,
      lavimplied = lavimplied, lavh1 = lavh1,
      lavcache = lavcache, lavoptions = lavoptions
    )
  } else if (information == "first.order") {
    I1 <- lav_model_h1_information_firstorder(
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats, lavdata = lavdata,
      lavimplied = lavimplied, lavh1 = lavh1,
      lavcache = lavcache, lavoptions = lavoptions
    )
  }

  # I1 information, as a list per group
  I1
}

# fisher/expected information of H1
lav_model_h1_information_expected <- function(lavobject = NULL,
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
    if (.hasSlot(lavobject, "h1")) {
      lavh1 <- lavobject@h1
    } else {
      lavh1 <- lav_h1_implied_logl(
        lavdata = lavobject@Data,
        lavsamplestats = lavobject@SampleStats,
        lavoptions = lavobject@Options
      )
    }
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
    A1 <- lavsamplestats@WLS.V
  }

  # 1b.
  else if (lavmodel@estimator == "DLS") {
    if (lavmodel@estimator.args$dls.GammaNT == "sample") {
      A1 <- lavsamplestats@WLS.V
    } else {
      A1 <- vector("list", length = lavsamplestats@ngroups)
      for (g in seq_len(lavsamplestats@ngroups)) {
        dls.a <- lavmodel@estimator.args$dls.a
        GammaNT <- lav_samplestats_Gamma_NT(
          COV            = lavimplied$cov[[g]],
          MEAN           = lavimplied$mean[[g]],
          rescale        = FALSE,
          x.idx          = lavsamplestats@x.idx[[g]],
          fixed.x        = lavmodel@fixed.x,
          conditional.x  = lavmodel@conditional.x,
          meanstructure  = lavmodel@meanstructure,
          slopestructure = lavmodel@conditional.x
        )
        W.DLS <- (1 - dls.a) * lavsamplestats@NACOV[[g]] + dls.a * GammaNT
        A1[[g]] <- lav_matrix_symmetric_inverse(W.DLS)
      }
    }
  }

  # 2. DWLS/ULS diagonal @WLS.VD slot
  else if (lavmodel@estimator == "DWLS" || lavmodel@estimator == "ULS") {
    # diagonal only!!
    A1 <- lavsamplestats@WLS.VD
  }

  # 3a. ML single level
  else if (lavmodel@estimator %in% c("ML", "NTRLS", "DLS", "catML") &&
    lavdata@nlevels == 1L) {
    A1 <- vector("list", length = lavsamplestats@ngroups)

    # structured? compute model-implied statistics
    if (structured && length(lavimplied) == 0L) {
      lavimplied <- lav_model_implied(lavmodel)
    }

    for (g in 1:lavsamplestats@ngroups) {
      if (.hasSlot(lavdata, "weights")) {
        WT <- lavdata@weights[[g]]
      } else {
        WT <- NULL
      }

      if (lavsamplestats@missing.flag) {
        # mvnorm
        # FIXME: allow for meanstructure = FALSE
        # FIXME: allow for conditional.x = TRUE
        if (lavmodel@meanstructure && structured) {
          MEAN <- lavimplied$mean[[g]]
        } else {
          MEAN <- lavsamplestats@missing.h1[[g]]$mu
        }

        if (structured) {
          A1[[g]] <-
            lav_mvnorm_missing_information_expected(
              Y = lavdata@X[[g]],
              Mp = lavdata@Mp[[g]],
              wt = WT,
              Mu = MEAN,
              # meanstructure = lavmodel@meanstructure,
              Sigma = lavimplied$cov[[g]],
              x.idx = lavsamplestats@x.idx[[g]]
            )
        } else {
          A1[[g]] <-
            lav_mvnorm_missing_information_expected(
              Y = lavdata@X[[g]],
              Mp = lavdata@Mp[[g]],
              wt = WT,
              Mu = MEAN,
              # meanstructure = lavmodel@meanstructure,
              Sigma = lavsamplestats@missing.h1[[g]]$sigma,
              x.idx = lavsamplestats@x.idx[[g]]
            )
        }
      } else {
        if (lavmodel@conditional.x) {
          # mvreg
          if (lavmodel@meanstructure && structured) {
            RES.INT <- lavimplied$res.int[[g]]
            RES.SLOPES <- lavimplied$res.slopes[[g]]
          } else {
            RES.INT <- lavsamplestats@res.int[[g]]
            RES.SLOPES <- lavsamplestats@res.slopes[[g]]
          }

          if (structured) {
            A1[[g]] <- lav_mvreg_information_expected(
              sample.mean.x     = lavsamplestats@mean.x[[g]],
              sample.cov.x      = lavsamplestats@cov.x[[g]],
              sample.nobs       = lavsamplestats@nobs[[g]],
              res.int           = RES.INT,
              res.slopes        = RES.SLOPES,
              # wt               = WT,
              # meanstructure    = lavmodel@meanstructure,
              res.cov           = lavimplied$res.cov[[g]]
            )
          } else {
            A1[[g]] <- lav_mvreg_information_expected(
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
            MEAN <- lavimplied$mean[[g]]
          } else {
            MEAN <- lavsamplestats@mean[[g]]
          }

          correlation.flag <- FALSE
          if (.hasSlot(lavmodel, "correlation")) {
            correlation.flag <- lavmodel@correlation
          }
          if (structured) {
            A1[[g]] <- lav_mvnorm_information_expected(
              Sigma         = lavimplied$cov[[g]],
              # wt = WT, # not needed
              x.idx         = lavsamplestats@x.idx[[g]],
              meanstructure = lavmodel@meanstructure,
              correlation   = correlation.flag
            )
          } else {
            A1[[g]] <- lav_mvnorm_h1_information_expected(
              sample.cov.inv = lavsamplestats@icov[[g]],
              # wt = WT, not needed
              x.idx          = lavsamplestats@x.idx[[g]],
              meanstructure  = lavmodel@meanstructure,
              correlation    = correlation.flag
            )
          }
        } # conditional.x
      } # missing

      # stochastic group weight
      if (lavmodel@group.w.free) {
        # unweight!! (as otherwise, we would 'weight' again)
        a <- exp(lavimplied$group.w[[g]]) / lavsamplestats@nobs[[g]]
        A1[[g]] <- lav_matrix_bdiag(matrix(a, 1L, 1L), A1[[g]])
      }
    } # g
  } # ML

  # 3b. ML + multilevel
  else if (lavmodel@estimator == "ML" && lavdata@nlevels > 1L) {
    A1 <- vector("list", length = lavsamplestats@ngroups)

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
      MU.W <- implied$mean[[(g - 1) * lavdata@nlevels + 1L]]
      MU.B <- implied$mean[[(g - 1) * lavdata@nlevels + 2L]]
      SIGMA.W <- implied$cov[[(g - 1) * lavdata@nlevels + 1L]]
      SIGMA.B <- implied$cov[[(g - 1) * lavdata@nlevels + 2L]]

      # clustered data
      A1[[g]] <- lav_mvnorm_cluster_information_expected(
        Lp           = lavdata@Lp[[g]],
        Mu.W         = MU.W,
        Sigma.W      = SIGMA.W,
        Mu.B         = MU.B,
        Sigma.B      = SIGMA.B,
        x.idx        = lavsamplestats@x.idx[[g]]
      )
    } # g
  } # ML + multilevel


  A1
}

lav_model_h1_information_observed <- function(lavobject = NULL,
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
    if (.hasSlot(lavobject, "h1")) {
      lavh1 <- lavobject@h1
    } else {
      lavh1 <- lav_h1_implied_logl(
        lavdata = lavobject@Data,
        lavsamplestats = lavobject@SampleStats,
        lavoptions = lavobject@Options
      )
    }
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
    A1 <- lavsamplestats@WLS.V
  }

  # 1b.
  else if (lavmodel@estimator == "DLS") {
    if (lavmodel@estimator.args$dls.GammaNT == "sample") {
      A1 <- lavsamplestats@WLS.V
    } else {
      A1 <- vector("list", length = lavsamplestats@ngroups)
      for (g in seq_len(lavsamplestats@ngroups)) {
        dls.a <- lavmodel@estimator.args$dls.a
        GammaNT <- lav_samplestats_Gamma_NT(
          COV            = lavimplied$cov[[g]],
          MEAN           = lavimplied$mean[[g]],
          rescale        = FALSE,
          x.idx          = lavsamplestats@x.idx[[g]],
          fixed.x        = lavmodel@fixed.x,
          conditional.x  = lavmodel@conditional.x,
          meanstructure  = lavmodel@meanstructure,
          slopestructure = lavmodel@conditional.x
        )
        W.DLS <- (1 - dls.a) * lavsamplestats@NACOV[[g]] + dls.a * GammaNT
        A1[[g]] <- lav_matrix_symmetric_inverse(W.DLS)
      }
    }
  }

  # 2. DWLS/ULS diagonal @WLS.VD slot
  else if (lavmodel@estimator == "DWLS" || lavmodel@estimator == "ULS") {
    # diagonal only!!
    A1 <- lavsamplestats@WLS.VD
  }

  # 3a. ML single level
  else if (lavmodel@estimator == "ML" && lavdata@nlevels == 1L) {
    A1 <- vector("list", length = lavsamplestats@ngroups)

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
          MEAN <- lavimplied$mean[[g]]
        } else {
          MEAN <- lavsamplestats@missing.h1[[g]]$mu
        }

        if (structured) {
          A1[[g]] <-
            lav_mvnorm_missing_information_observed_samplestats(
              Yp = lavsamplestats@missing[[g]],
              # wt not needed
              Mu = MEAN,
              # meanstructure = lavmodel@meanstructure,
              Sigma = lavimplied$cov[[g]],
              x.idx = lavsamplestats@x.idx[[g]]
            )
        } else {
          A1[[g]] <-
            lav_mvnorm_missing_information_observed_samplestats(
              Yp = lavsamplestats@missing[[g]],
              # wt not needed
              Mu = MEAN,
              # meanstructure = lavmodel@meanstructure,
              Sigma = lavsamplestats@missing.h1[[g]]$sigma,
              x.idx = lavsamplestats@x.idx[[g]]
            )
        }
      } else {
        if (lavmodel@conditional.x) {
          # mvreg
          if (lavmodel@meanstructure && structured) {
            RES.INT <- lavimplied$res.int[[g]]
            RES.SLOPES <- lavimplied$res.slopes[[g]]
          } else {
            RES.INT <- lavsamplestats@res.int[[g]]
            RES.SLOPES <- lavsamplestats@res.slopes[[g]]
          }

          if (structured) {
            A1[[g]] <- lav_mvreg_information_observed_samplestats(
              sample.res.int    = lavsamplestats@res.int[[g]],
              sample.res.slopes = lavsamplestats@res.slopes[[g]],
              sample.res.cov    = lavsamplestats@res.cov[[g]],
              sample.mean.x     = lavsamplestats@mean.x[[g]],
              sample.cov.x      = lavsamplestats@cov.x[[g]],
              res.int           = RES.INT,
              res.slopes        = RES.SLOPES,
              # wt               = WT,
              # meanstructure    = lavmodel@meanstructure,
              res.cov           = lavimplied$res.cov[[g]]
            )
          } else {
            A1[[g]] <- lav_mvreg_information_observed_samplestats(
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
            MEAN <- lavimplied$mean[[g]]
          } else {
            MEAN <- lavsamplestats@mean[[g]]
          }

          if (structured) {
            A1[[g]] <- lav_mvnorm_information_observed_samplestats(
              sample.mean   = lavsamplestats@mean[[g]],
              sample.cov    = lavsamplestats@cov[[g]],
              Mu            = MEAN,
              Sigma         = lavimplied$cov[[g]],
              # wt           = WT, # not needed
              x.idx         = lavsamplestats@x.idx[[g]],
              meanstructure = lavmodel@meanstructure
            )
          } else {
            A1[[g]] <- lav_mvnorm_h1_information_observed_samplestats(
              sample.mean    = lavsamplestats@mean[[g]],
              sample.cov     = lavsamplestats@cov[[g]],
              sample.cov.inv = lavsamplestats@icov[[g]],
              # wt            = WT, not needed
              x.idx          = lavsamplestats@x.idx[[g]],
              meanstructure  = lavmodel@meanstructure
            )
          }
        } # conditional.x
      } # missing

      # stochastic group weight
      if (lavmodel@group.w.free) {
        # unweight!!
        a <- exp(lavimplied$group.w[[g]]) / lavsamplestats@nobs[[g]]
        A1[[g]] <- lav_matrix_bdiag(matrix(a, 1, 1), A1[[g]])
      }
    } # g
  } # ML

  # 3b. ML + multilevel
  else if (lavmodel@estimator == "ML" && lavdata@nlevels > 1L) {
    A1 <- vector("list", length = lavsamplestats@ngroups)

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
      MU.W <- implied$mean[[(g - 1) * lavdata@nlevels + 1L]]
      MU.B <- implied$mean[[(g - 1) * lavdata@nlevels + 2L]]
      SIGMA.W <- implied$cov[[(g - 1) * lavdata@nlevels + 1L]]
      SIGMA.B <- implied$cov[[(g - 1) * lavdata@nlevels + 2L]]

      if (lavdata@missing == "ml") {
          A1[[g]] <- lav_mvnorm_cluster_missing_information_observed(
            Y1            = lavdata@X[[g]],
            Y2            = lavsamplestats@YLp[[g]][[2]]$Y2,
            Lp            = lavdata@Lp[[g]],
            Mp            = lavdata@Mp[[g]],
            Mu.W          = MU.W,
            Sigma.W       = SIGMA.W,
            Mu.B          = MU.B,
            Sigma.B       = SIGMA.B,
            x.idx         = lavsamplestats@x.idx[[g]]
          )
      } else {
        A1[[g]] <- lav_mvnorm_cluster_information_observed(
          Lp           = lavdata@Lp[[g]],
          YLp          = lavsamplestats@YLp[[g]],
          Mu.W         = MU.W,
          Sigma.W      = SIGMA.W,
          Mu.B         = MU.B,
          Sigma.B      = SIGMA.B,
          x.idx        = lavsamplestats@x.idx[[g]]
        )
      }
    } # g
  } # ML + multilevel

  A1
}

# outer product of the case-wise scores (gradients)
# HJ 18/10/2023: Adjust J matrix correctly using weights. Note: H matrix is
# based on lav_model_hessian so no changes required.
lav_model_h1_information_firstorder <- function(lavobject = NULL,
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
    if (.hasSlot(lavobject, "h1")) {
      lavh1 <- lavobject@h1
    } else {
      lavh1 <- lav_h1_implied_logl(
        lavdata = lavobject@Data,
        lavsamplestats = lavobject@SampleStats,
        lavoptions = lavobject@Options
      )
    }
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
    #    stop("lavaan ERROR: clustered information is not (yet) available when missing = \"ML\"")
    # }
    # if(lavmodel@conditional.x) {
    #    stop("lavaan ERROR: clustered information is not (yet) available when conditional.x = TRUE")
    # }
    # if(!structured) {
    #    stop("lavaan ERROR: clustered information is not (yet) available when h1.information = \"unstructured\"")
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

  B1 <- vector("list", length = lavsamplestats@ngroups)
  for (g in 1:lavdata@ngroups) {
    if (.hasSlot(lavdata, "weights")) {
      WT <- lavdata@weights[[g]]
    } else {
      WT <- NULL
    }

    if (estimator == "PML") {
      # slow approach: compute outer product of case-wise scores

      if (lavmodel@conditional.x) {
        SIGMA <- implied$res.cov[[g]]
        MU <- implied$res.mean[[g]]
        TH <- implied$res.th[[g]]
        PI <- implied$res.slopes[[g]]
        EXO <- lavdata@eXo[[g]]
      } else {
        SIGMA <- implied$cov[[g]]
        MU <- implied$mean[[g]]
        TH <- implied$th[[g]]
        PI <- NULL
        EXO <- NULL
      }
      SC <- pml_deriv1(
        Sigma.hat = SIGMA,
        Mu.hat = MU,
        TH = TH,
        th.idx = lavmodel@th.idx[[g]],
        num.idx = lavmodel@num.idx[[g]],
        X = lavdata@X[[g]],
        eXo = EXO,
        wt = NULL,
        PI = PI,
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
        clusters.idx <- lavdata@Lp[[g]]$cluster.idx[[2]]
        nclust       <- length(clusters)
        zb <- list()

        if (is.null(WT)) WT <- rep(1, length(clusters.idx))

        for (b in seq_along(clusters)) {
          SC_b <- SC[clusters.idx == b, ]
          WT_b <- WT[clusters.idx == b]
          zb[[b]] <- apply(SC_b * WT_b, 2, sum)
        }
        zbar <- apply(do.call(cbind, zb), 1, mean)
        B1c <-
          lapply(zb, \(z) tcrossprod(z - zbar)) |>
          Reduce(f = `+`)
        B1[[g]] <- nclust / (nclust - 1) * B1c

      } else {
        if (is.null(WT)) {
          B1[[g]] <- lav_matrix_crossprod(SC)
        } else {
          B1[[g]] <- crossprod(WT * SC)
        }
      }
    } else if (estimator == "ML" && lavdata@nlevels > 1L) {
      # if not-structured, we use lavh1, and that is always
      # 'unconditional' (for now)
      if (lavmodel@conditional.x && structured) {
	    if (lavdata@missing == "ml") {
		  lav_msg_stop(gettext("firstorder information matrix not available ",
		                       "(yet) if conditional.x + fiml"))
		}
        Res.Sigma.W <- implied$res.cov[[(g - 1) * lavdata@nlevels + 1L]]
        Res.Int.W <- implied$res.int[[(g - 1) * lavdata@nlevels + 1L]]
        Res.Pi.W <- implied$res.slopes[[(g - 1) * lavdata@nlevels + 1L]]

        Res.Sigma.B <- implied$res.cov[[(g - 1) * lavdata@nlevels + 2L]]
        Res.Int.B <- implied$res.int[[(g - 1) * lavdata@nlevels + 2L]]
        Res.Pi.B <- implied$res.slopes[[(g - 1) * lavdata@nlevels + 2L]]
        B1[[g]] <- lav_mvreg_cluster_information_firstorder(
          Y1            = lavdata@X[[g]],
          YLp           = lavsamplestats@YLp[[g]],
          Lp            = lavdata@Lp[[g]],
          Res.Sigma.W   = Res.Sigma.W,
          Res.Int.W     = Res.Int.W,
          Res.Pi.W      = Res.Pi.W,
          Res.Sigma.B   = Res.Sigma.B,
          Res.Int.B     = Res.Int.B,
          Res.Pi.B      = Res.Pi.B,
          divide.by.two = TRUE
        )
      } else {
        MU.W <- implied$mean[[(g - 1) * lavdata@nlevels + 1L]]
        MU.B <- implied$mean[[(g - 1) * lavdata@nlevels + 2L]]
        SIGMA.W <- implied$cov[[(g - 1) * lavdata@nlevels + 1L]]
        SIGMA.B <- implied$cov[[(g - 1) * lavdata@nlevels + 2L]]
		if (lavdata@missing == "ml") {
          B1[[g]] <- lav_mvnorm_cluster_missing_information_firstorder(
            Y1            = lavdata@X[[g]],
            Y2            = lavsamplestats@YLp[[g]][[2]]$Y2,
            Lp            = lavdata@Lp[[g]],
            Mp            = lavdata@Mp[[g]],
            Mu.W          = MU.W,
            Sigma.W       = SIGMA.W,
            Mu.B          = MU.B,
            Sigma.B       = SIGMA.B,
            x.idx         = lavsamplestats@x.idx[[g]],
            divide.by.two = TRUE
          )
		} else {
		  # no missing values
          B1[[g]] <- lav_mvnorm_cluster_information_firstorder(
            Y1            = lavdata@X[[g]],
            YLp           = lavsamplestats@YLp[[g]],
            Lp            = lavdata@Lp[[g]],
            Mu.W          = MU.W,
            Sigma.W       = SIGMA.W,
            Mu.B          = MU.B,
            Sigma.B       = SIGMA.B,
            x.idx         = lavsamplestats@x.idx[[g]],
            divide.by.two = TRUE
          )
		}
      }
    } else if (estimator == "ML" && lavdata@nlevels == 1L) {
      if (length(lavdata@cluster) > 0L) {
        cluster.idx <- lavdata@Lp[[g]]$cluster.idx[[2]]
      } else {
        cluster.idx <- NULL
      }

      if (lavsamplestats@missing.flag) {
        # mvnorm
        # FIXME: allow for meanstructure = FALSE
        # FIXME: allow for conditional.x = TRUE
        if (lavmodel@meanstructure && structured) {
          MEAN <- lavimplied$mean[[g]]
        } else {
          MEAN <- lavsamplestats@missing.h1[[g]]$mu
        }

        B1[[g]] <- lav_mvnorm_missing_information_firstorder(
          Y = lavdata@X[[g]],
          Mp = lavdata@Mp[[g]], wt = WT,
          cluster.idx = cluster.idx,
          Mu = MEAN,
          # meanstructure = lavmodel@meanstructure,
          Sigma = implied$cov[[g]],
          x.idx = lavsamplestats@x.idx[[g]]
        )
      } else {
        if (lavmodel@conditional.x) {
          # mvreg
          if (lavmodel@meanstructure && structured) {
            RES.INT <- lavimplied$res.int[[g]]
            RES.SLOPES <- lavimplied$res.slopes[[g]]
          } else {
            RES.INT <- lavsamplestats@res.int[[g]]
            RES.SLOPES <- lavsamplestats@res.slopes[[g]]
          }

          B1[[g]] <- lav_mvreg_information_firstorder(
            Y              = lavdata@X[[g]],
            eXo            = lavdata@eXo[[g]],
            res.int        = RES.INT,
            res.slopes     = RES.SLOPES,
            # wt            = WT,
            # meanstructure = lavmodel@meanstructure,
            res.cov        = implied$res.cov[[g]]
          )
        } else {
          # conditional.x = FALSE
          # mvnorm
          if (lavmodel@meanstructure && structured) {
            MEAN <- lavimplied$mean[[g]]
          } else {
            # NOTE: the information matrix will be the same (minus
            # the meanstructure block), but once INVERTED, the
            # standard errors will be (slightly) smaller!!!
            # This is only visibile when estimator = "MLF"
            # (or information = "first.order")
            MEAN <- lavsamplestats@mean[[g]] # saturated
          }

          if (structured) {
            B1[[g]] <- lav_mvnorm_information_firstorder(
              Y = lavdata@X[[g]],
              Mu = MEAN, Sigma = lavimplied$cov[[g]],
              wt = WT,
              cluster.idx = cluster.idx,
              x.idx = lavsamplestats@x.idx[[g]],
              meanstructure = lavmodel@meanstructure
            )
          } else {
            B1[[g]] <- lav_mvnorm_h1_information_firstorder(
              Y = lavdata@X[[g]],
              sample.cov.inv = lavsamplestats@icov[[g]],
              Gamma = lavsamplestats@NACOV[[g]],
              wt = WT,
              cluster.idx = cluster.idx, # only if wt
              x.idx = lavsamplestats@x.idx[[g]],
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
      B1[[g]] <- lav_matrix_bdiag(matrix(a, 1, 1), B1[[g]])
    }
  } # g

  B1
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
                              h1.information = NULL, # if specified, use it
                              se = NULL) { # if specified, use it

  if (!is.null(lavobject) && inherits(lavobject, "lavaan")) {
    lavmodel <- lavobject@Model
    lavsamplestats <- lavobject@SampleStats
    lavdata <- lavobject@Data
    lavimplied <- lavobject@implied
    if (.hasSlot(lavobject, "h1")) {
      lavh1 <- lavobject@h1
    } else {
      lavh1 <- lav_h1_implied_logl(
        lavdata = lavobject@Data,
        lavsamplestats = lavobject@SampleStats,
        lavoptions = lavobject@Options
      )
    }
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
  if (!is.null(h1.information)) {
    lavoptions$h1.information[1] <- h1.information
  }
  if (!is.null(se)) {
    lavoptions$se <- se
  }



  # information
  information <- lavoptions$information[1] # ALWAYS used the first

  # compute information matrix
  if (information == "observed") {
    I1 <- lav_model_h1_information_observed(
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats, lavdata = lavdata,
      lavimplied = lavimplied, lavh1 = lavh1,
      lavcache = lavcache, lavoptions = lavoptions
    )
  } else if (information == "expected") {
    I1 <- lav_model_h1_information_expected(
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats, lavdata = lavdata,
      lavimplied = lavimplied, lavh1 = lavh1,
      lavcache = lavcache, lavoptions = lavoptions
    )
  } else if (information == "first.order") {
    I1 <- lav_model_h1_information_firstorder(
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats, lavdata = lavdata,
      lavimplied = lavimplied, lavh1 = lavh1,
      lavcache = lavcache, lavoptions = lavoptions
    )
  }

  if (lavoptions$se %in% c("robust.huber.white", "robust.sem")) {
    J1 <- lav_model_h1_information_firstorder(
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats, lavdata = lavdata,
      lavimplied = lavimplied, lavh1 = lavh1,
      lavcache = lavcache, lavoptions = lavoptions
    )
  }

  # compute ACOV per group
  ACOV <- vector("list", length = lavdata@ngroups)
  for (g in 1:lavdata@ngroups) {
    # denominator
    if (lavdata@nlevels == 1L) {
      Ng <- lavsamplestats@nobs[[g]]
    } else {
      Ng <- lavdata@Lp[[g]]$nclusters[[2]]
    }

    # invert information
    I1.g.inv <- try(lav_matrix_symmetric_inverse(I1[[g]]), silent = TRUE)
    if (inherits(I1.g.inv, "try-error")) {
      lav_msg_stop(gettext(
        "could not invert h1 information matrix in group"), g)
    }

    # which type of se?
    if (lavoptions$se %in% c("standard", "none")) {
      ACOV[[g]] <- 1 / Ng * I1.g.inv
    } else if (lavoptions$se %in% c("robust.huber.white", "robust.sem")) {
      ACOV[[g]] <- 1 / Ng * (I1.g.inv %*% J1[[g]] %*% I1.g.inv)
    }
  }

  ACOV
}
