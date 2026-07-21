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
# - YR 08 Jul 2026: merge the (almost fully parallel) expected/observed
#                   scaffolding into a single worker (lav_model_h1_info_ed)


## For the lavaan.mi package, TDJ provides pooled versions of all the
## sample moments called in these functions.  If any updates to these functions
## require NEW information (from @SampleStats or @implied or @h1),
## PLEASE ADD A TAG     @TDJorgensen
## at the end of the commit message on GitHub, so TDJ can check whether
## lavaan.mi::lavResiduals.mi() needs to be updated accordingly.


# all lav_model_h1_* functions accept either a fitted lavaan object
# (lavobject) or the individual slots; this helper unpacks the object into
# the slot variables of the CALLER's frame, and (optionally) fills in
# lavh1/lavimplied when they were not supplied
lav_model_h1_unpack <- function(lavobject = NULL, env = parent.frame(),
                                need_h1 = FALSE, need_implied = FALSE) {
  if (!is.null(lavobject) && inherits(lavobject, "lavaan")) {
    env$lavmodel <- lavobject@Model
    env$lavsamplestats <- lavobject@SampleStats
    env$lavdata <- lavobject@Data
    env$lavimplied <- lavobject@implied
    env$lavh1 <- lavobject@h1
    env$lavcache <- lavobject@Cache
    env$lavoptions <- lavobject@Options
  }

  # sanity check
  if (need_h1 && length(env$lavh1) == 0L) {
    env$lavh1 <- lav_h1_implied_logl(
      lavdata = env$lavdata,
      lavsamplestats = env$lavsamplestats,
      lavoptions = env$lavoptions
    )
  }
  if (need_implied && length(env$lavimplied) == 0L) {
    env$lavimplied <- lav_model_implied(lavmodel = env$lavmodel)
  }

  invisible(NULL)
}


lav_model_h1_info <- function(lavobject = NULL,
                                     lavmodel = NULL,
                                     lavsamplestats = NULL,
                                     lavdata = NULL,
                                     lavimplied = NULL,
                                     lavh1 = NULL,
                                     lavcache = NULL,
                                     lavoptions = NULL) {
  lav_model_h1_unpack(lavobject, need_h1 = TRUE, need_implied = TRUE)

  # information
  information <- lavoptions$information[1] # ALWAYS take the first one
  # the caller must control it


  # compute information matrix
  if (information == "observed") {
    i1 <- lav_model_h1_info_observed(
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats, lavdata = lavdata,
      lavimplied = lavimplied, lavh1 = lavh1,
      lavcache = lavcache, lavoptions = lavoptions
    )
  } else if (information == "expected") {
    i1 <- lav_model_h1_info_expected(
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats, lavdata = lavdata,
      lavimplied = lavimplied, lavh1 = lavh1,
      lavcache = lavcache, lavoptions = lavoptions
    )
  } else if (information == "first.order") {
    i1 <- lav_model_h1_info_firstorder(
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats, lavdata = lavdata,
      lavimplied = lavimplied, lavh1 = lavh1,
      lavcache = lavcache, lavoptions = lavoptions
    )
  }

  # I1 information, as a list per group
  i1
}

# shared worker for the expected and observed h1 information: the two are
# structurally identical (same estimator dispatch, same structured /
# unstructured moment selection, same group.w.free handling) -- only the
# innermost kernel calls differ. Returns NULL for unsupported estimators
# (the named wrappers below raise the error, so the message carries their
# name).
lav_model_h1_info_ed <- function(what = "expected",
                                 lavmodel = NULL,
                                 lavsamplestats = NULL,
                                 lavdata = NULL,
                                 lavoptions = NULL,
                                 lavimplied = NULL,
                                 lavh1 = NULL) {
  observed <- (what == "observed")

  # sanity check
  # (not needed for the estimators that read the precomputed WLS.V/WLS.VD
  #  slots -- and for two-level categorical WLS there is no h1 object)
  if (length(lavh1) == 0L &&
    !lavmodel@estimator %in% c("GLS", "WLS", "DLS", "DWLS", "ULS")) {
    lavh1 <- lav_h1_implied_logl(
      lavdata = lavdata,
      lavsamplestats = lavsamplestats,
      lavoptions = lavoptions
    )
  }
  if (length(lavimplied) == 0L) {
    lavimplied <- lav_model_implied(lavmodel = lavmodel)
  }

  # structured or unstructured? (since 0.5-23)
  if (!is.null(lavoptions) &&
    !is.null(lavoptions$h1.information[1]) &&
    lavoptions$h1.information[1] == "unstructured") {
    structured <- FALSE
  } else {
    structured <- TRUE
  }

  # which single-level ML-family estimators are handled by branch 3a below?
  # the observed h1 information is only available for "ML" itself
  # (estimator "DLS" is caught by the dedicated branch 1b)
  if (observed) {
    ml_family <- "ML"
  } else {
    ml_family <- c("ML", "NTRLS", "catML", "IV", "JS", "JSA", "MGM")
  }

  # 1. WLS.V (=A1) for GLS/WLS
  if (lavmodel@estimator == "GLS" || lavmodel@estimator == "WLS") {
    a1 <- lavsamplestats@WLS.V
    # GLS: since 0.7-1 the NT weight matrix is no longer pre-computed at
    # the samplestats stage (estimation does not need it); build it here
    # on demand. Only reached for the 'plain' single-level variant: the
    # correlation/conditional.x/group.w.free/two.stage variants still
    # store WLS.V eagerly (see lav_samp_from_data).
    if (lavmodel@estimator == "GLS" && is.null(a1[[1]])) {
      for (g in seq_len(lavsamplestats@ngroups)) {
        a1[[g]] <- lav_samp_wls_v_nt_g(
          m_cov         = lavsamplestats@cov[[g]],
          m_mean        = lavsamplestats@mean[[g]],
          m_icov        = lavsamplestats@icov[[g]],
          correlation   = FALSE,
          x_idx         = lavsamplestats@x.idx[[g]],
          fixed_x       = lavmodel@fixed.x,
          conditional_x = FALSE,
          meanstructure = lavmodel@meanstructure
        )
        # mimic Mplus: V11 rescale (only for full-data input; the
        # moments-based lav_samp_from_moments never applied it)
        if (isTRUE(lavmodel@estimator.args$gls.v11.mplus) &&
          lavmodel@meanstructure &&
          !is.null(lavdata) && lavdata@data.type == "full") {
          nvar <- NCOL(lavsamplestats@cov[[g]])
          a1[[g]][1:nvar, 1:nvar] <- a1[[g]][1:nvar, 1:nvar, drop = FALSE] *
            (lavsamplestats@nobs[[g]] / (lavsamplestats@nobs[[g]] - 1))
        }
      }
    }

  # 1b. DLS: for dls.GammaNT = "model", the weight matrix is a function of
  #     the model-implied moments and must be recomputed here (the @WLS.V
  #     slot only holds the plain inv(NACOV) fallback)
  } else if (lavmodel@estimator == "DLS") {
    if (lavmodel@estimator.args$dls.GammaNT == "sample") {
      a1 <- lavsamplestats@WLS.V
    } else {
      a1 <- vector("list", length = lavsamplestats@ngroups)
      for (g in seq_len(lavsamplestats@ngroups)) {
        a1[[g]] <- lav_dls_wls_v_g(
          m_cov         = lavimplied$cov[[g]],
          m_mean        = lavimplied$mean[[g]],
          nacov_g       = lavsamplestats@NACOV[[g]],
          dls_a         = lavmodel@estimator.args$dls.a,
          x_idx         = lavsamplestats@x.idx[[g]],
          fixed_x       = lavmodel@fixed.x,
          conditional_x = lavmodel@conditional.x,
          meanstructure = lavmodel@meanstructure
        )
      }
    }

  # 2. DWLS/ULS diagonal @WLS.VD slot
  } else if (lavmodel@estimator == "DWLS" || lavmodel@estimator == "ULS") {
    # diagonal only!!
    a1 <- lavsamplestats@WLS.VD

  # 3a. ML single level
  } else if (lavmodel@estimator %in% ml_family && lavdata@nlevels == 1L) {
    a1 <- vector("list", length = lavsamplestats@ngroups)

    for (g in 1:lavsamplestats@ngroups) {
      if (lavsamplestats@missing.flag && lavmodel@conditional.x) {
        # mvreg, missing data
        if (lavmodel@meanstructure && structured) {
          res_int <- lavimplied$res.int[[g]]
          res_slopes <- lavimplied$res.slopes[[g]]
        } else {
          res_int <- lavh1$implied$res.int[[g]]
          res_slopes <- lavh1$implied$res.slopes[[g]]
        }
        if (structured) {
          res_cov <- lavimplied$res.cov[[g]]
        } else {
          res_cov <- lavh1$implied$res.cov[[g]]
        }

        if (observed) {
          a1[[g]] <- lav_mvreg_mi_information_observed_samplestats(
            yp = lavsamplestats@missing[[g]],
            res_int = res_int,
            res_slopes = res_slopes,
            res_cov = res_cov
          )
        } else {
          a1[[g]] <- lav_mvreg_mi_info_expected(
            yp = lavsamplestats@missing[[g]],
            res_cov = res_cov
          )
        }
      } else if (lavsamplestats@missing.flag) {
        # mvnorm, missing data
        # FIXME: allow for meanstructure = FALSE
        if (lavmodel@meanstructure && structured) {
          mean_1 <- lavimplied$mean[[g]]
        } else {
          mean_1 <- lavh1$implied$mean[[g]]
        }
        if (structured) {
          sigma_1 <- lavimplied$cov[[g]]
        } else {
          sigma_1 <- lavh1$implied$cov[[g]]
        }

        if (observed) {
          a1[[g]] <-
            lav_mvnorm_missing_information_observed_samplestats(
              yp = lavsamplestats@missing[[g]],
              # wt not needed
              mu = mean_1,
              sigma_1 = sigma_1,
              x_idx = lavsamplestats@x.idx[[g]]
            )
        } else {
          a1[[g]] <-
            lav_mvn_mi_info_expected(
              y = lavdata@X[[g]],
              mp = lavdata@Mp[[g]],
              wt = lavdata@weights[[g]],
              mu = mean_1,
              sigma_1 = sigma_1,
              x_idx = lavsamplestats@x.idx[[g]]
            )
        }
      } else if (lavmodel@conditional.x) {
        # mvreg
        if (lavmodel@meanstructure && structured) {
          res_int <- lavimplied$res.int[[g]]
          res_slopes <- lavimplied$res.slopes[[g]]
        } else {
          res_int <- lavsamplestats@res.int[[g]]
          res_slopes <- lavsamplestats@res.slopes[[g]]
        }
        if (structured) {
          res_cov <- lavimplied$res.cov[[g]]
        } else {
          res_cov <- lavsamplestats@res.cov[[g]]
        }

        if (observed) {
          a1[[g]] <- lav_mvreg_information_observed_samplestats(
            sample_res_int    = lavsamplestats@res.int[[g]],
            sample_res_slopes = lavsamplestats@res.slopes[[g]],
            sample_res_cov    = lavsamplestats@res.cov[[g]],
            sample_mean_x     = lavsamplestats@mean.x[[g]],
            sample_cov_x      = lavsamplestats@cov.x[[g]],
            res_int           = res_int,
            res_slopes        = res_slopes,
            res_cov           = res_cov
          )
        } else {
          a1[[g]] <- lav_mvreg_info_expected(
            sample_mean_x     = lavsamplestats@mean.x[[g]],
            sample_cov_x      = lavsamplestats@cov.x[[g]],
            sample_nobs       = lavsamplestats@nobs[[g]],
            res_int           = res_int,
            res_slopes        = res_slopes,
            res_cov           = res_cov
          )
        }
      } else {
        # complete data, conditional.x = FALSE
        # mvnorm
        if (observed) {
          if (lavmodel@meanstructure && structured) {
            mean_1 <- lavimplied$mean[[g]]
          } else {
            mean_1 <- lavsamplestats@mean[[g]]
          }

          if (structured) {
            a1[[g]] <- lav_mvn_info_observed_samp(
              sample_mean   = lavsamplestats@mean[[g]],
              sample_cov    = lavsamplestats@cov[[g]],
              mu            = mean_1,
              sigma_1       = lavimplied$cov[[g]],
              x_idx         = lavsamplestats@x.idx[[g]],
              meanstructure = lavmodel@meanstructure
            )
          } else {
            a1[[g]] <- lav_mvn_h1_info_observed_samp(
              sample_mean    = lavsamplestats@mean[[g]],
              sample_cov     = lavsamplestats@cov[[g]],
              sample_cov_inv = lavsamplestats@icov[[g]],
              x_idx          = lavsamplestats@x.idx[[g]],
              meanstructure  = lavmodel@meanstructure
            )
          }
        } else {
          # note: the correlation= argument selects the reduced
          # (correlation-metric) information (catML); the D-augmented ML
          # mode (free ~*~ scales) lives in the FULL moment space
          cor_reduced <- lavmodel@correlation &&
            !lav_model_delta_free(lavmodel)
          if (structured) {
            a1[[g]] <- lav_mvn_info_expected(
              sigma_1       = lavimplied$cov[[g]],
              x_idx         = lavsamplestats@x.idx[[g]],
              meanstructure = lavmodel@meanstructure,
              correlation   = cor_reduced
            )
          } else {
            a1[[g]] <- lav_mvn_h1_info_expected(
              sample_cov_inv = lavsamplestats@icov[[g]],
              x_idx          = lavsamplestats@x.idx[[g]],
              meanstructure  = lavmodel@meanstructure,
              correlation    = cor_reduced
            )
          }
        }
      } # complete/missing, conditional.x

      # stochastic group weight
      if (lavmodel@group.w.free) {
        # unweight!! (as otherwise, we would 'weight' again)
        a <- exp(lavimplied$group.w[[g]]) / lavsamplestats@nobs[[g]]
        a1[[g]] <- lav_mat_bdiag(matrix(a, 1L, 1L), a1[[g]])
      }
    } # g

  # 3b. ML + multilevel
  } else if (lavmodel@estimator == "ML" && lavdata@nlevels > 1L) {
    a1 <- vector("list", length = lavsamplestats@ngroups)

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

      if (observed && lavdata@missing %in% c("ml", "ml.x")) {
        a1[[g]] <- lav_mvn_cl_mi_info_observed(
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
      } else if (observed) {
        a1[[g]] <- lav_mvn_cl_info_observed(
          lp           = lavdata@Lp[[g]],
          ylp          = lavsamplestats@YLp[[g]],
          mu_w         = mu_w,
          sigma_w      = sigma_w,
          mu_b         = mu_b,
          sigma_b      = sigma_b,
          x_idx        = lavsamplestats@x.idx[[g]]
        )
      } else {
        a1[[g]] <- lav_mvn_cl_info_expected(
          lp           = lavdata@Lp[[g]],
          mu_w         = mu_w,
          sigma_w      = sigma_w,
          mu_b         = mu_b,
          sigma_b      = sigma_b,
          x_idx        = lavsamplestats@x.idx[[g]]
        )
      }
    } # g
  } else {
    # no branch matched (e.g. PML): the wrappers raise the error
    return(NULL)
  }

  a1
}

# fisher/expected information of H1
lav_model_h1_info_expected <- function(lavobject = NULL,
                                              lavmodel = NULL,
                                              lavsamplestats = NULL,
                                              lavdata = NULL,
                                              lavoptions = NULL,
                                              lavimplied = NULL,
                                              lavh1 = NULL,
                                              lavcache = NULL) {
  lav_model_h1_unpack(lavobject)

  a1 <- lav_model_h1_info_ed(
    what = "expected",
    lavmodel = lavmodel, lavsamplestats = lavsamplestats,
    lavdata = lavdata, lavoptions = lavoptions,
    lavimplied = lavimplied, lavh1 = lavh1
  )
  if (is.null(a1)) {
    # unsupported estimator (e.g. PML): fail with a clear message instead
    # of an obscure "object 'a1' not found"
    lav_msg_stop(gettextf(
      "the (expected) h1 information matrix is not available for
       estimator %s.", lavmodel@estimator))
  }

  a1
}

lav_model_h1_info_observed <- function(lavobject = NULL,
                                              lavmodel = NULL,
                                              lavsamplestats = NULL,
                                              lavdata = NULL,
                                              lavimplied = NULL,
                                              lavh1 = NULL,
                                              lavcache = NULL,
                                              lavoptions = NULL) {
  lav_model_h1_unpack(lavobject)

  a1 <- lav_model_h1_info_ed(
    what = "observed",
    lavmodel = lavmodel, lavsamplestats = lavsamplestats,
    lavdata = lavdata, lavoptions = lavoptions,
    lavimplied = lavimplied, lavh1 = lavh1
  )
  if (is.null(a1)) {
    # unsupported estimator (e.g. PML): fail with a clear message instead
    # of an obscure "object 'a1' not found"
    lav_msg_stop(gettextf(
      "the (observed) h1 information matrix is not available for
       estimator %s.", lavmodel@estimator))
  }

  a1
}

# outer product of the case-wise scores (gradients)
# HJ 18/10/2023: Adjust J matrix correctly using weights. Note: H matrix is
# based on lav_model_hessian so no changes required.
lav_model_h1_info_firstorder <- function(lavobject = NULL,
                                                lavmodel = NULL,
                                                lavsamplestats = NULL,
                                                lavdata = NULL,
                                                lavimplied = NULL,
                                                lavh1 = NULL,
                                                lavcache = NULL,
                                                lavoptions = NULL) {
  lav_model_h1_unpack(lavobject, need_h1 = TRUE, need_implied = TRUE)

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
          # drop = FALSE: a single-observation cluster must stay a matrix
          sc_b <- sc[clusters_idx == b, , drop = FALSE]
          wt_b <- wt[clusters_idx == b]
          zb[[b]] <- colSums(sc_b * wt_b)
        }
        zbar <- apply(do.call(cbind, zb), 1, mean)
        b1c <- Reduce(f = `+`, lapply(zb, function(z) tcrossprod(z - zbar)))
        b1[[g]] <- nclust / (nclust - 1) * b1c

      } else {
        if (is.null(wt)) {
          b1[[g]] <- lav_mat_crossprod(sc)
        } else {
          # sc holds the unit (unweighted) casewise scores. How the sampling
          # weights enter the first-order information (the meat of the robust
          # sandwich / the scaled test) is controlled by sampling.weights.type:
          #  - "design"    : crossprod(wt * sc), weighted by sum(wt^2), the
          #                  design-based version (matches Mplus).
          #  - "frequency" : crossprod(sqrt(wt) * sc), weighted by sum(wt),
          #                  i.e. treats the weights as frequencies, so the
          #                  result matches a fit on the row-replicated data.
          swt_type <- if (!is.null(lavoptions$sampling.weights.type)) {
            lavoptions$sampling.weights.type
          } else {
            "design"
          }
          b1[[g]] <- lav_samp_wt_crossprod(sc, wt, swt_type)
        }
      }


    # ML
    } else if (estimator == "ML" && lavdata@nlevels > 1L) {
      # if not-structured, we use lavh1, and that is always
      # 'unconditional' (for now)
      if (lavmodel@conditional.x && structured) {
      if (lavdata@missing %in% c("ml", "ml.x")) {
      lav_msg_stop(gettext("firstorder information matrix not available
                                (yet) if conditional.x + fiml"))
    }
        res_sigma_w <- implied$res.cov[[(g - 1) * lavdata@nlevels + 1L]]
        res_int_w <- implied$res.int[[(g - 1) * lavdata@nlevels + 1L]]
        res_pi_w <- implied$res.slopes[[(g - 1) * lavdata@nlevels + 1L]]

        res_sigma_b <- implied$res.cov[[(g - 1) * lavdata@nlevels + 2L]]
        res_int_b <- implied$res.int[[(g - 1) * lavdata@nlevels + 2L]]
        res_pi_b <- implied$res.slopes[[(g - 1) * lavdata@nlevels + 2L]]
        b1[[g]] <- lav_mvreg_cl_info_firstorder(
          y1            = lavdata@X[[g]],
          ylp           = lavsamplestats@YLp[[g]],
          lp            = lavdata@Lp[[g]],
          res_sigma_w   = res_sigma_w,
          res_int_w     = res_int_w,
          res_pi_w      = res_pi_w,
          res_sigma_b   = res_sigma_b,
          res_int_b     = res_int_b,
          res_pi_b      = res_pi_b,
          divide_by_two = TRUE
        )
        # reorder from the kernel statistic order to the Delta row order
        perm <- lav_mvreg_cl_stat_perm(
          res_int_w = res_int_w, res_pi_w = res_pi_w,
          res_int_b = res_int_b, res_pi_b = res_pi_b
        )
        b1[[g]] <- b1[[g]][perm, perm, drop = FALSE]
      } else {
        # note: if conditional.x, the (unconditional) unstructured b1
        # may only be combined with matching unconditional-space factors
        # (as in the yuan.bentler.mplus trace), never with the
        # conditional.x Delta; user-requested h1.information =
        # "unstructured" is blocked in lav_options_set()
        mu_w <- implied$mean[[(g - 1) * lavdata@nlevels + 1L]]
        mu_b <- implied$mean[[(g - 1) * lavdata@nlevels + 2L]]
        sigma_w <- implied$cov[[(g - 1) * lavdata@nlevels + 1L]]
        sigma_b <- implied$cov[[(g - 1) * lavdata@nlevels + 2L]]
    if (lavdata@missing %in% c("ml", "ml.x")) {
          b1[[g]] <- lav_mvn_cl_mi_info_firstorder(
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
          b1[[g]] <- lav_mvn_cl_info_firstorder(
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

      if (lavsamplestats@missing.flag && lavmodel@conditional.x) {
        # mvreg, missing data
        if (lavmodel@meanstructure && structured) {
          res_int <- lavimplied$res.int[[g]]
          res_slopes <- lavimplied$res.slopes[[g]]
        } else {
          res_int <- lavh1$implied$res.int[[g]]
          res_slopes <- lavh1$implied$res.slopes[[g]]
        }

        b1[[g]] <- lav_mvreg_mi_info_firstorder(
          y = lavdata@X[[g]],
          exo = lavdata@eXo[[g]],
          mp = lavdata@Mp[[g]], wt = wt,
          cluster_idx = cluster_idx,
          res_int = res_int,
          res_slopes = res_slopes,
          res_cov = implied$res.cov[[g]]
        )
      } else if (lavsamplestats@missing.flag) {
        # mvnorm
        # FIXME: allow for meanstructure = FALSE
        if (lavmodel@meanstructure && structured) {
          mean_1 <- lavimplied$mean[[g]]
        } else {
          #MEAN <- lavsamplestats@missing.h1[[g]]$mu
          mean_1 <- lavh1$implied$mean[[g]]
        }

        b1[[g]] <- lav_mvn_mi_info_firstorder(
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

          b1[[g]] <- lav_mvreg_info_firstorder(
            y              = lavdata@X[[g]],
            exo            = lavdata@eXo[[g]],
            res_int        = res_int,
            res_slopes     = res_slopes,
            # wt            = WT,
            # meanstructure = lavmodel@meanstructure,
            res_cov        = implied$res.cov[[g]]
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
            b1[[g]] <- lav_mvn_info_firstorder(
              y = lavdata@X[[g]],
              mu = mean_1, sigma_1 = lavimplied$cov[[g]],
              wt = wt,
              cluster_idx = cluster_idx,
              x_idx = lavsamplestats@x.idx[[g]],
              meanstructure = lavmodel@meanstructure
            )
          } else {
            b1[[g]] <- lav_mvn_h1_info_firstorder(
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
      b1[[g]] <- lav_mat_bdiag(matrix(a, 1, 1), b1[[g]])
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
  lav_model_h1_unpack(lavobject, need_h1 = TRUE, need_implied = TRUE)

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


  # compute information matrix (dispatch on lavoptions$information[1])
  i1 <- lav_model_h1_info(
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats, lavdata = lavdata,
    lavimplied = lavimplied, lavh1 = lavh1,
    lavcache = lavcache, lavoptions = lavoptions
  )

  if (lavoptions$se %in% c("robust.huber.white", "robust.sem",
                           "robust.cluster", "robust.cluster.sem")) {
    # for clustered data (robust.cluster[.sem]), the first-order kernels
    # aggregate the casewise scores within clusters (they receive
    # cluster_idx whenever lavdata@cluster is set), so the same sandwich
    # below yields the cluster-robust ACOV
    j1 <- lav_model_h1_info_firstorder(
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
    #
    # some h1 sample statistics can be fixed/degenerate, with zero information
    # (e.g. the within-level means in a multilevel model, which are fixed at
    # zero). Those rows/columns make the full information matrix singular, so
    # we drop them before inverting and reinsert zeros afterwards (a zero
    # asymptotic variance for fixed statistics). For single-level models there
    # are no such columns and this reduces to the plain inverse.
    diag_i1 <- diag(i1[[g]])
    keep <- which(abs(diag_i1) > sqrt(.Machine$double.eps) * max(abs(diag_i1)))
    if (length(keep) < nrow(i1[[g]])) {
      i1_g_inv <- matrix(0, nrow(i1[[g]]), ncol(i1[[g]]))
      inv_keep <- try(
        lav_mat_sym_inverse(i1[[g]][keep, keep, drop = FALSE]),
        silent = TRUE
      )
      if (inherits(inv_keep, "try-error")) {
        lav_msg_stop(gettext(
          "could not invert h1 information matrix in group"), g)
      }
      i1_g_inv[keep, keep] <- inv_keep
    } else {
      i1_g_inv <- try(lav_mat_sym_inverse(i1[[g]]), silent = TRUE)
      if (inherits(i1_g_inv, "try-error")) {
        lav_msg_stop(gettext(
          "could not invert h1 information matrix in group"), g)
      }
    }

    # which type of se?
    if (lavoptions$se %in% c("standard", "none")) {
      acov[[g]] <- 1 / ng * i1_g_inv
    } else if (lavoptions$se %in% c("robust.huber.white", "robust.sem",
                                    "robust.cluster", "robust.cluster.sem")) {
      acov[[g]] <- 1 / ng * (i1_g_inv %*% j1[[g]] %*% i1_g_inv)
    } else {
      lav_msg_stop(gettextf(
        "h1 ACOV not available for se = %s.", dQuote(lavoptions$se)))
    }
  }

  acov
}
