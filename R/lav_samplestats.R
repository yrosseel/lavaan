# constructor for the 'lavSampleStats' class
#
# initial version: YR 25/03/2009
# major revision: YR 5/11/2011: separate data.obs and sample statistics
# YR 5/01/2016: add rescov, resvar, ... if conditional.x = TRUE

# YR 18 Jan 2021: use lavoptions

lav_samp_from_data <- function(lavdata = NULL,        # nolint start
                                      lavoptions = NULL,
                                      wls_v = NULL,
                                      nacov = NULL) {        # nolint end
  # extra info from lavoptions
  stopifnot(!is.null(lavoptions))
  missing <- lavoptions$missing
  rescale <- lavoptions$sample.cov.rescale
  estimator <- lavoptions$estimator
  # mimic <- lavoptions$mimic
  meanstructure <- lavoptions$meanstructure
  correlation <- lavoptions$correlation
  conditional_x <- lavoptions$conditional.x
  fixed_x <- lavoptions$fixed.x
  group_w_free <- lavoptions$group.w.free
  se <- lavoptions$se
  test <- lavoptions$test
  ridge <- lavoptions$ridge
  zero_add <- lavoptions$zero.add
  zero_keep_margins <- lavoptions$zero.keep.margins
  zero_cell_warn <- lavoptions$zero.cell.warn
  allow_empty_cell <- lavoptions$allow.empty.cell
  dls_a <- lavoptions$estimator.args$dls.a
  dls_gamma_nt <- lavoptions$estimator.args$dls.GammaNT
  # design vs frequency weighting of the (categorical/continuous) Gamma
  swt_type <- if (!is.null(lavoptions$sampling.weights.type)) {
    lavoptions$sampling.weights.type
  } else {
    "design"
  }

  # sample.icov (new in 0.6-9; ensure it exists, for older objects)
  sample_icov <- TRUE
  if (!is.null(lavoptions$sample.icov)) {
    sample_icov <- lavoptions$sample.icov
  }

  # ridge default
  if (ridge) {
    if (is.numeric(lavoptions$ridge.constant)) {
      ridge_eps <- lavoptions$ridge.constant
    } else {
      ridge_eps <- 1e-5
    }
  } else {
    ridge_eps <- 0.0
  }

  # check lavdata
  stopifnot(!is.null(lavdata))

  # lavdata slots (FIXME: keep lavdata@ names)
  x <- lavdata@X
  mp <- lavdata@Mp
  ngroups <- lavdata@ngroups
  nlevels <- lavdata@nlevels
  nobs <- lavdata@nobs
  ov_names <- lavdata@ov.names
  ov_names_x <- lavdata@ov.names.x
  data_ov <- lavdata@ov
  exo <- lavdata@eXo
  wt <- lavdata@weights

  # new in 0.6-6
  # if sampling weights have been used, redefine nobs:
  # per group, we define nobs == sum(wt)
  for (g in seq_len(ngroups)) {
    if (!is.null(wt[[g]])) {
      nobs[[g]] <- sum(wt[[g]])
    }
  }

  # sample.cov.robust cannot be used if sampling weights are used
  if (lavoptions$sample.cov.robust) {
    if (!is.null(wt[[1]])) {
      lav_msg_stop(gettext(
        "sample.cov.robust = TRUE does not work (yet)
        if sampling weights are provided."))
    }
  }

  # sample statistics per group

  # joint (y,x)
  cov <- vector("list", length = ngroups)
  var <- vector("list", length = ngroups)
  mean <- vector("list", length = ngroups)
  th <- vector("list", length = ngroups)
  th_idx <- vector("list", length = ngroups)
  th_names <- vector("list", length = ngroups)

  # residual (y | x)
  res_cov <- vector("list", length = ngroups)
  res_var <- vector("list", length = ngroups)
  res_th <- vector("list", length = ngroups)
  res_th_nox <- vector("list", length = ngroups)
  res_slopes <- vector("list", length = ngroups)
  res_int <- vector("list", length = ngroups)

  # fixed.x
  mean_x <- vector("list", length = ngroups)
  cov_x <- vector("list", length = ngroups)

  # binary/ordinal
  bifreq <- vector("list", length = ngroups)

  # extra sample statistics per group
  icov <- vector("list", length = ngroups)
  cov_log_det <- vector("list", length = ngroups)
  res_icov <- vector("list", length = ngroups)
  res_cov_log_det <- vector("list", length = ngroups)
  wls_obs <- vector("list", length = ngroups)
  missing_1 <- vector("list", length = ngroups)
  missing_h1 <- vector("list", length = ngroups)
  missing_flag <- FALSE
  zero_cell_tables <- vector("list", length = ngroups)
  ylp <- vector("list", length = ngroups)

  # group weights
  group_w <- vector("list", length = ngroups)

  # convenience? # FIXME!
  x_idx <- vector("list", length = ngroups)


  wls_vd <- vector("list", length = ngroups)
  if (is.null(wls_v)) {
    wls_v <- vector("list", length = ngroups)
    wls_v_user <- FALSE
  } else {
    if (!is.list(wls_v)) {
      if (ngroups == 1L) {
        wls_v <- list(wls_v)
      } else {
        lav_msg_stop(gettextf(
          "wls_v argument should be a list of length %s", ngroups)
        )
      }
    } else {
      if (length(wls_v) != ngroups) {
        lav_msg_stop(gettextf(
          "wls_v assumes %1$s groups; data contains %2$s groups",
          length(wls_v), ngroups))
      }
    }

    # is wls_v full? check first
    if (is.null(dim(wls_v[[1]]))) {
      # we will assume it is the diagonal only
      wls_vd <- wls_v
      wls_v <- lapply(wls_vd, diag)
    } else {
      # create WLS.VD
      wls_vd <- lapply(wls_v, diag)
    }

    wls_v_user <- TRUE
    # FIXME: check dimension of wls_v!!
  }

  nacov_compute <- FALSE # since 0.6-6
  if (is.null(nacov)) {
    nacov <- vector("list", length = ngroups)
    nacov_user <- FALSE
    if (se %in% c("robust.sem", "robust.sem.nt", "robust.cluster.sem") &&
                                                    missing == "listwise") {
      nacov_compute <- TRUE
    }
    # note: test can be a vector...
    if (missing == "listwise" && any(test %in% c(
      "satorra.bentler",
      "mean.var.adjusted",
      "scaled.shifted"
    ))) {
      nacov_compute <- TRUE
    }
    if (missing == "listwise" &&
        any(vapply(test, lav_test_fmg_is_fmg, logical(1L)))) {
      nacov_compute <- TRUE
    }
    if (estimator == "IV" &&
        lavoptions$estimator.args$iv_vcov_stage1 == "gamma") {
      nacov_compute <- TRUE
    }
    # two-stage missing data for the (continuous) least-squares estimators:
    # the robust.sem SEs (and satorra.bentler test) need the two-stage NACOV
    # of the EM moments (computed below, see lav_mvn_mi_* functions)
    if (any(missing == c("two.stage", "robust.two.stage")) &&
        estimator %in% c("ULS", "GLS", "WLS", "DLS")) {
      nacov_compute <- TRUE
    }
  } else if (is.logical(nacov)) {
    if (!nacov) {
      nacov_compute <- FALSE
    } else {
      nacov_compute <- TRUE
    }
    nacov_user <- FALSE
    nacov <- vector("list", length = ngroups)
  } else {
    if (!is.list(nacov)) {
      if (ngroups == 1L) {
        nacov <- list(nacov)
      } else {
        lav_msg_stop(gettextf(
          "nacov= argument should be a list of length ", ngroups))
      }
    } else {
      if (length(nacov) != ngroups) {
        lav_msg_stop(gettextf(
          "nacov= assumes %1$s groups; data contains %2$s groups",
          length(nacov), ngroups))
      }
    }
    nacov_user <- TRUE
    # FIXME: check dimension of nacov!!
  }



  # compute some sample statistics per group
  for (g in 1:ngroups) {
    # switch off computing all sample statistics? (housekeeping only)
    if (!is.null(lavoptions$samplestats) && !lavoptions$samplestats) {
    next
  }
    # check nobs
    if (is.null(wt[[g]])) {
      if (nobs[[g]] < 2L) {
        if (nobs[[g]] == 0L) {
          if (ngroups > 1L) {
            lav_msg_stop(gettextf("data contains no observations in
                                   group %s", g))
          } else {
            lav_msg_stop(gettext("data contains no observations"))
          }
        } else if (!allow_empty_cell) {
          if (ngroups > 1L) {
            lav_msg_stop(gettextf("data contains only a single observation
                                   in group %s", g))
          } else {
            lav_msg_stop(gettext("data contains only a single observation"))
          }
        }
      }
    }

    # exogenous x?
    nexo <- length(ov_names_x[[g]])
    if (nexo) {
      stopifnot(nexo == NCOL(exo[[g]]))

      # two cases: ov.names contains 'x' variables, or not
      if (conditional_x) {
        # ov.names.x are NOT in ov.names
        x_idx[[g]] <- length(ov_names[[g]]) + seq_len(nexo)
      } else {
        if (fixed_x) {
          # ov.names.x are a subset of ov.names
          x_idx[[g]] <- match(ov_names_x[[g]], ov_names[[g]])
          stopifnot(!anyNA(x_idx[[g]]))
        } else {
          x_idx[[g]] <- integer(0L)
        }
      }
    } else {
      x_idx[[g]] <- integer(0L)
      conditional_x <- FALSE
      fixed_x <- FALSE
    }

    # group weight
    group_w[[g]] <- nobs[[g]] / sum(unlist(nobs))

    # check if we have categorical data in this group
    categorical <- FALSE
    ov_types <- data_ov$type[match(ov_names[[g]], data_ov$name)]
    ov_levels <- data_ov$nlev[match(ov_names[[g]], data_ov$name)]
    cat_1 <- list()
    if ("ordered" %in% ov_types) {
      categorical <- TRUE
      if (nlevels > 1L) {
        lav_msg_warn(gettext("multilevel + categorical not supported yet."))
      }
    }

    if (categorical) {
      # compute CAT

      if (estimator %in% c("ML", "REML", "PML", "FML", "MML", "none", "IV",
                           "ULS")) {
        wls_w <- FALSE
        if (estimator == "ULS" && se %in% c("robust.sem", "robust.sem.nt")) {
          wls_w <- TRUE
        } else if (estimator == "IV" &&
                  lavoptions$estimator.args$iv_vcov_stage1 == "gamma") {
          wls_w <- TRUE
        }
      } else {
        wls_w <- TRUE
      }
    # check cat.wls.w option (new in 0.6-18)
    if (!is.null(lavoptions$cat.wls.w) && !lavoptions$cat.wls.w) {
      wls_w <- FALSE # perhaps do.fit = FALSE? (eg sam())
    }
      if (lav_verbose()) {
        cat("Estimating sample thresholds and correlations ... ")
      }

      current_verbose <- lav_verbose()
      if (lav_verbose(lav_debug()))
        on.exit(lav_verbose(current_verbose), TRUE)
      if (conditional_x) {
        cat_1 <- muthen1984(
          data_1 = x[[g]],
          wt = wt[[g]],
          sampling_weights_type = swt_type,
          ov_names = ov_names[[g]],
          ov_types = ov_types,
          ov_levels = ov_levels,
          ov_names_x = ov_names_x[[g]],
          exo = exo[[g]],
          group = g, # for error messages only
          wls_w = wls_w,
          zero_add = zero_add,
          zero_keep_margins = zero_keep_margins,
          zero_cell_warn = FALSE,
          zero_cell_tables = TRUE,
          allow_empty_cell = allow_empty_cell
        )
      } else {
        cat_1 <- muthen1984(
          data_1 = x[[g]],
          wt = wt[[g]],
          sampling_weights_type = swt_type,
          ov_names = ov_names[[g]],
          ov_types = ov_types,
          ov_levels = ov_levels,
          ov_names_x = NULL,
          exo = NULL,
          group = g, # for error messages only
          wls_w = wls_w,
          zero_add = zero_add,
          zero_keep_margins = zero_keep_margins,
          zero_cell_warn = FALSE,
          zero_cell_tables = TRUE,
          allow_empty_cell = allow_empty_cell
        )
      }
      lav_verbose(current_verbose)
      # empty cell tables
      zero_cell_tables[[g]] <- cat_1$zero.cell.tables
      if (lav_verbose()) cat("done\n")
    }

    if (categorical) {
      # convenience
      th_idx[[g]] <- unlist(cat_1$TH.IDX)
      th_names[[g]] <- unlist(cat_1$TH.NAMES)

      if (conditional_x) {
        # residual var/cov
        res_var[[g]] <- unlist(cat_1$VAR)
        res_cov[[g]] <- unname(cat_1$COV)
        if (ridge) {
          diag(res_cov[[g]]) <- diag(res_cov[[g]]) + ridge_eps
          res_var[[g]] <- diag(res_cov[[g]])
        }

        # th also contains the means of numeric variables
        res_th[[g]] <- unlist(cat_1$TH)
        res_th_nox[[g]] <- unlist(cat_1$TH.NOX)

        # for convenience, we store the intercept of numeric
        # variables in res.int
        nvar_1 <- NCOL(res_cov[[g]])
        mean[[g]] <- res_int[[g]] <- numeric(nvar_1)
        num_idx <- which(!seq_len(nvar_1) %in% th_idx[[g]])
        if (length(num_idx) > 0L) {
          num_idx_1 <- which(th_idx[[g]] == 0L)
          mean[[g]][num_idx] <- res_th_nox[[g]][num_idx_1]
          res_int[[g]][num_idx] <- res_th[[g]][num_idx_1]
        }

        # slopes
        res_slopes[[g]] <- cat_1$SLOPES
      } else {
        # var/cov
        var[[g]] <- unlist(cat_1$VAR)
        cov[[g]] <- unname(cat_1$COV)
        if (ridge) {
          diag(cov[[g]]) <- diag(cov[[g]]) + ridge_eps
          var[[g]] <- diag(cov[[g]])
        }

        # th also contains the means of numeric variables
        th[[g]] <- unlist(cat_1$TH)

        # mean (numeric only)
        nvar_1 <- NCOL(cov[[g]])
        mean[[g]] <- numeric(nvar_1)
        num_idx <- which(!seq_len(nvar_1) %in% th_idx[[g]])
        if (length(num_idx) > 0L) {
          num_idx_1 <- which(th_idx[[g]] == 0L)
          mean[[g]][num_idx] <- th[[g]][num_idx_1]
        }
      }

      # only for catML
      if (estimator == "catML") {
        cov_1 <- cov2cor(lav_mat_sym_force_pd(cov[[g]],
          tol = 1e-04
        ))
        # overwrite
        cov[[g]] <- cov_1
        out <- lav_samp_icov(
          cov_1 = cov_1,
          x_idx = x_idx[[g]],
          ngroups = ngroups, g = g
        )
        icov[[g]] <- out$icov
        cov_log_det[[g]] <- out$cov.log.det

        # the same for res.cov if conditional.x = TRUE
        if (conditional_x) {
          res_cov_1 <-
            cov2cor(lav_mat_sym_force_pd(res_cov[[g]],
              tol = 1e-04
            ))
          # overwrite
          res_cov[[g]] <- res_cov_1
          out <- lav_samp_icov(
            cov_1 = res_cov_1,
            ridge = 1e-05,
            x_idx = x_idx[[g]],
            ngroups = ngroups, g = g
          )
          res_icov[[g]] <- out$icov
          res_cov_log_det[[g]] <- out$cov.log.det
        }
      }
    # categorical

    # continuous -- multilevel
    } else if (nlevels > 1L) {
      # level-based sample statistics
      ylp[[g]] <- lav_samp_cl_patterns(
        y = x[[g]],
        lp = lavdata@Lp[[g]],
        conditional_x = lavoptions$conditional.x
      )


      if (conditional_x) {
        # for starting values only
        # no handling of missing data yet....
        if (missing %in% c(
          "ml", "ml.x",
          "two.stage", "robust.two.stage"
        )) {
          lav_msg_stop(gettextf(
            "missing = %s + conditional.x + two.level not supported yet",
            missing))
        }

        # residual covariances!
        y <- x[[g]] # contains eXo
        cov_1 <- unname(stats::cov(y, use = "pairwise.complete.obs"))
    # if we have missing values (missing by design?), replace them by 0
    cov_1[is.na(cov_1)] <- 0
        mean_1 <- unname(colMeans(y, na.rm = TRUE))
        var[[g]] <- diag(cov_1)
        # rescale cov by (N-1)/N? (only COV!)
        if (rescale) {
          # we 'transform' the sample cov (divided by n-1)
          # to a sample cov divided by 'n'
          cov_1 <- ((nobs[[g]] - 1) / nobs[[g]]) * cov_1
        }
        cov[[g]] <- cov_1
        if (ridge) {
          diag(cov[[g]]) <- diag(cov[[g]]) + ridge_eps
          var[[g]] <- diag(cov[[g]])
        }
        mean[[g]] <- mean_1

        a_1 <- cov_1[-x_idx[[g]], -x_idx[[g]], drop = FALSE]
        m_b <- cov_1[-x_idx[[g]], x_idx[[g]], drop = FALSE]
        c_1 <- cov_1[x_idx[[g]], x_idx[[g]], drop = FALSE]
        c_inv <- try(solve(c_1), silent = TRUE)
        if (inherits(c_inv, "try-error")) {
          lav_msg_warn(gettext(
            "Exogenous covariance matrix is singular; using generalized inverse."
          ))
          c_inv <- MASS::ginv(c_1)
        }
        res_cov[[g]] <- a_1 - m_b %*% c_inv %*% t(m_b)
        res_var[[g]] <- diag(cov[[g]])

        my <- mean_1[-x_idx[[g]]]
        mx <- mean_1[x_idx[[g]]]
        c3 <- rbind(
          c(1, mx),
          cbind(mx, c_1 + tcrossprod(mx))
        )
        b3 <- cbind(my, m_b + tcrossprod(my, mx))
        coef_1 <- try(unname(solve(c3, t(b3))), silent = TRUE)
        if (inherits(coef_1, "try-error")) {
          lav_msg_warn(gettext(
            "Augmented exogenous matrix is singular; using generalized inverse."
          ))
          coef_1 <- unname(MASS::ginv(c3) %*% t(b3))
        }

        res_int[[g]] <- coef_1[1, ] # intercepts
        res_slopes[[g]] <- t(coef_1[-1, , drop = FALSE]) # slopes
      } else {
        # FIXME: needed?
    cov_1 <- unname(stats::cov(x[[g]], use = "pairwise.complete.obs"))
        # if we have missing values (missing by design?), replace them by 0
        cov_1[is.na(cov_1)] <- 0
        cov[[g]] <- cov_1
        mean[[g]] <- unname(colMeans(x[[g]], na.rm = TRUE))
        var[[g]] <- diag(cov[[g]])

        # missing patterns
        if (missing %in% c("ml", "ml.x")) {
          missing_flag <- TRUE
          missing_1[[g]] <-
            lav_samp_mi_patterns(
              y = x[[g]],
              mp = mp[[g]],
              wt = wt[[g]],
              lp = lavdata@Lp[[g]]
            )
        }
      }
    } # multilevel

    # continuous -- single-level
    else {
      if (conditional_x) {
        # FIXME!
        # no correlation structures yet
        if (correlation) {
          lav_msg_stop(gettext(
            "conditional.x = TRUE is not supported (yet) for
            correlation structures."))
        }

        # FIXME!
        # no handling of missing data yet....
        if (missing %in% c(
          "ml", "ml.x",
          "two.stage", "robust.two.stage"
        )) {
          lav_msg_stop(gettextf(
            "missing = %s + conditional.x not supported yet", missing))
        }

        # residual covariances!

        y <- cbind(x[[g]], exo[[g]])
    cov_1 <- unname(stats::cov(y, use = "pairwise.complete.obs"))
        # if we have missing values (missing by design?), replace them by 0
        cov_1[is.na(cov_1)] <- 0
        mean_1 <- unname(colMeans(y, na.rm = TRUE))
        # rescale cov by (N-1)/N? (only COV!)
        if (rescale) {
          # we 'transform' the sample cov (divided by n-1)
          # to a sample cov divided by 'n'
          cov_1 <- ((nobs[[g]] - 1) / nobs[[g]]) * cov_1
        }
        cov[[g]] <- cov_1
        var[[g]] <- diag(cov_1)
        if (ridge) {
          diag(cov[[g]]) <- diag(cov[[g]]) + ridge_eps
          var[[g]] <- diag(cov[[g]])
        }
        mean[[g]] <- mean_1

        a_1 <- cov_1[-x_idx[[g]], -x_idx[[g]], drop = FALSE]
        m_b <- cov_1[-x_idx[[g]], x_idx[[g]], drop = FALSE]
        c_1 <- cov_1[x_idx[[g]], x_idx[[g]], drop = FALSE]
        # FIXME: make robust against singular c_1!!!
        res_cov[[g]] <- a_1 - m_b %*% solve(c_1) %*% t(m_b)
        res_var[[g]] <- diag(cov[[g]])


        my <- mean_1[-x_idx[[g]]]
        mx <- mean_1[x_idx[[g]]]
        c3 <- rbind(
          c(1, mx),
          cbind(mx, c_1 + tcrossprod(mx))
        )
        b3 <- cbind(my, m_b + tcrossprod(my, mx))
        coef_1 <- unname(solve(c3, t(b3)))

        res_int[[g]] <- coef_1[1, ] # intercepts
        res_slopes[[g]] <- t(coef_1[-1, , drop = FALSE]) # slopes
      } else if (missing == "two.stage" ||
        missing == "robust.two.stage") {
        missing_flag <- FALSE # !!! just use sample statistics
        missing_1[[g]] <-
          lav_samp_mi_patterns(
            y = x[[g]],
            mp = mp[[g]],
            wt = wt[[g]]
          )
        current_warn <- lav_warn()
        if (lav_warn(lavoptions$em.h1.warn))
          on.exit(lav_warn(current_warn), TRUE)
        out <- lav_mvn_mi_h1_est_moments(
          y = x[[g]],
          wt = wt[[g]],
          mp = mp[[g]], yp = missing_1[[g]],
          max_iter = lavoptions$em.h1.iter.max,
          tol = lavoptions$em.h1.tol,
        )
        lav_warn(current_warn)
        missing_h1[[g]]$sigma <- out$Sigma
        missing_h1[[g]]$mu <- out$Mu
        missing_h1[[g]]$h1 <- out$fx

        # here, sample statistics == EM estimates
        cov[[g]] <- missing_h1[[g]]$sigma
        if (ridge) {
          diag(cov[[g]]) <- diag(cov[[g]]) + ridge_eps
        }
        var[[g]] <- diag(cov[[g]])
        mean[[g]] <- missing_h1[[g]]$mu
      } else if (missing %in% c("ml", "ml.x")) {
        missing_flag <- TRUE
        missing_1[[g]] <-
          lav_samp_mi_patterns(
            y = x[[g]],
            mp = mp[[g]],
            wt = wt[[g]]
          )

        if (nlevels == 1L) {
          # estimate moments unrestricted model
          current_warn <- lav_warn()
          if (lav_warn(lavoptions$em.h1.warn))
            on.exit(lav_warn(current_warn), TRUE)
          # zero coverage?
          if (any(lav_mat_vech(mp[[g]]$coverage, diagonal = FALSE) == 0)) {
            #out <- lav_mvn_mi_h1_est_moments_chol(
            #         lavdata = lavdata, lavoptions = lavoptions, group = g)
            #missing.h1.[[g]]$sigma <- out$Sigma
            #missing.h1.[[g]]$mu <- out$Mu
            #missing.h1.[[g]]$h1 <- out$fx
          } else if (!allow_empty_cell) {
            out <- lav_mvn_mi_h1_est_moments(
              y = x[[g]],
              wt = wt[[g]],
              mp = mp[[g]], yp = missing_1[[g]],
              max_iter = lavoptions$em.h1.iter.max,
              tol = lavoptions$em.h1.tol
            )
            missing_h1[[g]]$sigma <- out$Sigma
            missing_h1[[g]]$mu <- out$Mu
            missing_h1[[g]]$h1 <- out$fx
          }
          lav_warn(current_warn)
        }

        if (!is.null(wt[[g]])) {
          # here, sample statistics == EM estimates
          cov[[g]] <- missing_h1[[g]]$sigma
          if (ridge) {
            diag(cov[[g]]) <- diag(cov[[g]]) + ridge_eps
          }
          var[[g]] <- diag(cov[[g]])
          mean[[g]] <- missing_h1[[g]]$mu
        } else {
          # NEEDED? why not just EM-based?
          cov_1 <- unname(stats::cov(x[[g]], use = "pairwise.complete.obs"))
          # if we have missing values (missing by design?), replace them by 0
      cov_1[is.na(cov_1)] <- 0
          cov[[g]] <- cov_1
          # rescale cov by (N-1)/N? (only COV!)
          if (rescale) {
            # we 'transform' the sample cov (divided by n-1)
            # to a sample cov divided by 'n'
            cov[[g]] <- ((nobs[[g]] - 1) / nobs[[g]]) * cov[[g]]
          }
          if (ridge) {
            diag(cov[[g]]) <- diag(cov[[g]]) + ridge_eps
          }
          var[[g]] <- diag(cov[[g]])
          mean[[g]] <- colMeans(x[[g]], na.rm = TRUE)
        }
      } else {
        # LISTWISE
        if (!is.null(wt[[g]])) {
          out <- stats::cov.wt(x[[g]],
            wt = wt[[g]],
            method = "ML"
          )
      cov_1 <- out$cov
          # if we have missing values (missing by design?), replace them by 0
          cov_1[is.na(cov_1)] <- 0
          cov[[g]] <- cov_1
          # cov.wt(method = "ML") divides by sum(wt) (the 'N' version). If
          # rescale is FALSE (e.g. GLS/ULS/(D)WLS), the unweighted path uses
          # the unbiased 'N-1' version; mirror that here so that supplying
          # sampling weights does not change the covariance normalization
          # (nobs[[g]] == sum(wt[[g]]); see top of this function).
          if (!rescale) {
            cov[[g]] <- (nobs[[g]] / (nobs[[g]] - 1)) * cov[[g]]
          }
          if (ridge) {
            diag(cov[[g]]) <- diag(cov[[g]]) + ridge_eps
          }
          var[[g]] <- diag(cov[[g]])
          mean[[g]] <- out$center
        } else if (lavoptions$sample.cov.robust) {
          # fixme: allow prob/max.it to be options
          out <- lav_cov_huber(
            y = x[[g]], prob = 0.95,
            max_it = 200L, tol = 1e-07
          )
          cov[[g]] <- out$Sigma
          var[[g]] <- diag(cov[[g]])
          mean[[g]] <- out$Mu
        } else {
      cov_1 <- unname(stats::cov(x[[g]], use = "pairwise.complete.obs"))
          # if we have missing values (missing by design?), replace them by 0
          cov_1[is.na(cov_1)] <- 0
          cov[[g]] <- cov_1
          # rescale cov by (N-1)/N? (only COV!)
          if (rescale) {
            # we 'transform' the sample cov (divided by n-1)
            # to a sample cov divided by 'n'
            cov[[g]] <- ((nobs[[g]] - 1) / nobs[[g]]) * cov[[g]]
          }
          if (ridge) {
            diag(cov[[g]]) <- diag(cov[[g]]) + ridge_eps
          }
          var[[g]] <- diag(cov[[g]])
          mean[[g]] <- colMeans(x[[g]], na.rm = TRUE)
        }
      }

      # correlation structure?
      if (correlation) {
        cov[[g]] <- cov2cor(cov[[g]])
        var[[g]] <- rep(1, length(var[[g]]))
        if (conditional_x) {
          res_cov[[g]] <- cov2cor(res_cov[[g]])
          res_var[[g]] <- rep(1, length(res_var[[g]]))
          cov_x[[g]] <- cov2cor(cov_x[[g]])
          # FIXME: slopes? more?
        }
      }

      # icov and cov.log.det (but not if missing)
      if (sample_icov && !missing %in% c("ml", "ml.x")) {
        out <- lav_samp_icov(
          cov_1 = cov[[g]], ridge = 1e-05,
          x_idx = x_idx[[g]],
          ngroups = ngroups, g = g
        )
        icov[[g]] <- out$icov
        cov_log_det[[g]] <- out$cov.log.det

        # the same for res.cov if conditional.x = TRUE
        if (conditional_x) {
          out <- lav_samp_icov(
            cov_1 = res_cov[[g]],
            ridge = 1e-05,
            x_idx = x_idx[[g]],
            ngroups = ngroups, g = g
          )
          res_icov[[g]] <- out$icov
          res_cov_log_det[[g]] <- out$cov.log.det
        }
      }
    } # continuous - single level


    # WLS.obs
    if (nlevels == 1L) {
      if (estimator == "catML") {
        # correlations only (for now)
        tmp_categorical <- FALSE
        tmp_meanstructure <- FALSE
      } else {
        tmp_categorical <- categorical
        tmp_meanstructure <- meanstructure
      }
      wls_obs[[g]] <- lav_samp_wls_obs(
        mean_g = mean[[g]],
        cov_g = cov[[g]], var_g = var[[g]], th_g = th[[g]],
        th_idx_g = th_idx[[g]], res_int_g = res_int[[g]],
        res_cov_g = res_cov[[g]], res_var_g = res_var[[g]],
        res_th_g = res_th[[g]], res_slopes_g = res_slopes[[g]],
        group_w_g = log(nobs[[g]]),
        categorical = tmp_categorical, conditional_x = conditional_x,
        meanstructure = tmp_meanstructure, correlation = correlation,
        slopestructure = conditional_x,
        group_w_free = group_w_free
      )
    }

    # fill in the other slots
    if (!is.null(exo[[g]])) {
      if (!is.null(wt[[g]])) {
        if (missing != "listwise") {
          cov_x[[g]] <- missing_h1[[g]]$sigma[x_idx[[g]],
            x_idx[[g]],
            drop = FALSE
          ]
          mean_x[[g]] <- missing_h1[[g]]$mu[x_idx[[g]]]
        } else {
          out <- stats::cov.wt(exo[[g]],
            wt = wt[[g]],
            method = "ML"
          )
          cov_x[[g]] <- out$cov
          mean_x[[g]] <- out$center
        }
      } else {
        cov_x[[g]] <- cov(exo[[g]], use = "pairwise")
        if (rescale) {
          # we 'transform' the sample cov (divided by n-1)
          # to a sample cov divided by 'n'
          cov_x[[g]] <- ((nobs[[g]] - 1) / nobs[[g]]) * cov_x[[g]]
        }
        mean_x[[g]] <- colMeans(exo[[g]])
      }
    }

    # nacov (=GAMMA)
    if (!nacov_user && nlevels == 1L) {
      if (nacov_compute && !categorical &&
          any(missing == c("two.stage", "robust.two.stage")) &&
          estimator %in% c("ULS", "GLS", "WLS", "DLS")) {
        # two-stage missing data: the NACOV of the (saturated) EM moments,
        # in (mean, vech(cov)) order, to be used by the robust.sem sandwich
        # and the satorra.bentler test
        # - two.stage:        I_1^{-1}            (Savalei & Bentler, 2009)
        # - robust.two.stage: I_1^{-1} J_1 I_1^{-1} (Savalei & Falk, 2014)
        mu_g <- missing_h1[[g]]$mu
        sigma_g <- missing_h1[[g]]$sigma
        x_idx_g <- if (fixed_x) x_idx[[g]] else integer(0L)
        if (missing == "robust.two.stage") {
          nacov[[g]] <- lav_mvn_mi_h1_omega_sw(
            y = x[[g]], mp = mp[[g]],
            yp = missing_1[[g]], wt = wt[[g]],
            mu = mu_g, sigma_1 = sigma_g, x_idx = x_idx_g,
            information = "observed"
          )
        } else {
          i1 <- lav_mvnorm_missing_information_observed_samplestats(
            yp = missing_1[[g]],
            mu = mu_g, sigma_1 = sigma_g, x_idx = x_idx_g
          )
          nacov[[g]] <- lav_mat_sym_inverse(i1)
        }
      } else if (estimator %in% c("ML", "GLS") &&
                 !missing_flag && nacov_compute) {
        if (conditional_x) {
          y <- y
        } else {
          y <- x[[g]]
        }

        if (length(lavdata@cluster) > 0L) {
          cluster_idx <- lavdata@Lp[[g]]$cluster.idx[[2]]
        } else {
          cluster_idx <- NULL
        }

        if (correlation) {
      nacov[[g]] <- lav_samp_cor_gamma(
        m_y = y,
        meanstructure = meanstructure
      )
    } else {
          nacov[[g]] <-
            lav_samp_gamma(
              m_y = y,
              x_idx = x_idx[[g]],
              cluster_idx = cluster_idx,
              fixed_x = fixed_x,
              conditional_x = conditional_x,
              meanstructure = meanstructure,
              slopestructure = conditional_x,
              gamma_n_minus_one =
                lavoptions$gamma.n.minus.one,
              unbiased = lavoptions$gamma.unbiased,
              mplus_wls = FALSE,
              wt = wt[[g]],
              sampling_weights_type = swt_type
            )
        }
      } else if (estimator %in% c("WLS", "DWLS", "ULS", "DLS", "IV", "catML")) {
        if (!categorical) {
          # sample size large enough?
          nvar <- ncol(x[[g]])
          # if(conditional.x && nexo > 0L) {
          #    nvar <- nvar - nexo
          # }
          pstar <- nvar * (nvar + 1) / 2
          if (meanstructure) pstar <- pstar + nvar
          if (conditional_x && nexo > 0L) {
            pstar <- pstar + (nvar * nexo)
          }
          if (nrow(x[[g]]) < pstar && estimator != "IV") {
            if (ngroups > 1L) {
              lav_msg_warn(gettextf(
              "number of observations (%s) too small to compute Gamma",
              nrow(x[[g]]), " in group %s", g))
            } else {
              lav_msg_warn(gettextf(
              "number of observations (%s) too small to compute Gamma",
              nrow(x[[g]])))
            }
          }
          if (conditional_x) {
            y <- y
          } else {
            y <- x[[g]]
          }

          if (length(lavdata@cluster) > 0L) {
            cluster_idx <- lavdata@Lp[[g]]$cluster.idx[[2]]
          } else {
            cluster_idx <- NULL
          }
      if (correlation) {
            nacov[[g]] <- lav_samp_cor_gamma(
              m_y = y,
              meanstructure = meanstructure
            )
      } else {
            if ("robust.sem.nt" %in% lavoptions$se ||
                "browne.residual.nt" %in% lavoptions$test) {
              nacov[[g]] <-
                lav_samp_gamma_nt(
                  m_y = y,
                  m_cov = cov[[g]],
                  m_mean = mean[[g]],
                  x_idx = x_idx[[g]],
                  # cluster.idx = cluster.idx, # not available
                  fixed_x = fixed_x,
                  conditional_x = conditional_x,
                  meanstructure = meanstructure,
                  slopestructure = conditional_x
                )
            } else {
              nacov[[g]] <-
                lav_samp_gamma(
                  m_y = y,
                  x_idx = x_idx[[g]],
                  cluster_idx = cluster_idx,
                  fixed_x = fixed_x,
                  conditional_x = conditional_x,
                  meanstructure = meanstructure,
                  slopestructure = conditional_x,
                  gamma_n_minus_one =
                    lavoptions$gamma.n.minus.one,
                  unbiased =
                    lavoptions$gamma.unbiased,
                  mplus_wls = lavoptions$gamma.wls.mplus,
                  wt = wt[[g]],
                  sampling_weights_type = swt_type
                )
            }
          }
        } else { # categorical case
          if (!is.null(cat_1$WLS.W)) {
            nacov[[g]] <- cat_1$WLS.W * nobs[[g]]
            if (lavoptions$gamma.n.minus.one) {
              nacov[[g]] <- nacov[[g]] * (nobs[[g]] / (nobs[[g]] - 1L))
            }
          }
          if (estimator == "catML") {
            # remove all but the correlation part
            ntotal <- nrow(nacov[[g]])
            pstar <- nrow(cat_1$A22)
            nocor <- ntotal - pstar
            if (length(nocor) > 0L) {
              nacov[[g]] <- nacov[[g]][
                -seq_len(nocor),
                -seq_len(nocor)
              ]
            }
          }
        }
      } else if (estimator == "PML") {
        # no nacov ... for now
      }

      # group.w.free
      if (!is.null(nacov[[g]]) && group_w_free) {
        # unweight!!
        a <- group_w[[g]] * sum(unlist(nobs)) / nobs[[g]]
        # always 1!!!
        nacov[[g]] <- lav_mat_bdiag(matrix(a, 1, 1), nacov[[g]])
      }
    }

    # wls_v
    if (!wls_v_user && nlevels == 1L) {
      if (estimator == "DLS" && dls_gamma_nt == "sample" && dls_a < 1.0) {
        # compute GammaNT here
        if (correlation) {
          gamma_nt <- lav_samp_cor_gamma_nt(
            m_cov          = cov[[g]],
            m_mean         = mean[[g]],
            rescale        = FALSE,
            x_idx          = x_idx[[g]],     # not used yet
            fixed_x        = fixed_x,        # not used yet
            conditional_x  = conditional_x,  # not used yet
            meanstructure  = meanstructure,  # not used yet
            slopestructure = conditional_x   # not used yet
          )
        } else {
          gamma_nt <- lav_samp_gamma_nt(
            m_cov          = cov[[g]],
            m_mean         = mean[[g]],
            rescale        = FALSE,
            x_idx          = x_idx[[g]],
            fixed_x        = fixed_x,
            conditional_x  = conditional_x,
            meanstructure  = meanstructure,
            slopestructure = conditional_x
          )
        }
      }

      if (estimator == "GLS" ||
        (estimator == "DLS" && dls_gamma_nt == "sample" &&
          dls_a == 1.0)) {
        # Note: we need the 'original' COV/MEAN/ICOV
        #        sample statistics; not the 'residual' version
        if (correlation) {
          gamma_nt <- lav_samp_cor_gamma_nt(
            m_cov          = cov[[g]],
            m_mean         = mean[[g]],
            #rescale        = FALSE,
            x_idx          = x_idx[[g]],     # not used yet
            fixed_x        = fixed_x,        # not used yet
            conditional_x  = conditional_x,  # not used yet
            meanstructure  = meanstructure,  # not used yet
            slopestructure = conditional_x   # not used yet
          )
          wls_v[[g]] <- lav_mat_sym_inverse(gamma_nt)
        } else {
          wls_v[[g]] <- lav_samp_gamma_inverse_nt(
            m_icov         = icov[[g]],
            m_cov          = cov[[g]],
            m_mean         = mean[[g]],
            rescale        = FALSE,
            x_idx          = x_idx[[g]],
            fixed_x        = fixed_x,
            conditional_x  = conditional_x,
            meanstructure  = meanstructure,
            slopestructure = conditional_x
          )
        }
        if (lavoptions$gls.v11.mplus && !conditional_x && meanstructure) {
          # bug in Mplus? V11 rescaled by nobs[[g]]/(nobs[[g]]-1)
          nvar <- NCOL(cov[[g]])
          wls_v[[g]][1:nvar, 1:nvar] <-
            wls_v[[g]][1:nvar, 1:nvar,
              drop = FALSE
            ] * (nobs[[g]] / (nobs[[g]] - 1))
        }
      } else if (estimator == "ML") {
        # no wls_v here, since function of model-implied moments
      } else if (estimator %in% c("WLS", "DWLS", "ULS", "DLS", "IV")) {
        if (!categorical) {
          if (estimator == "WLS" || estimator == "DLS") {
            if (!fixed_x) {
              if (estimator != "DLS") {
                # Gamma should be po before we invert
                ev <- eigen(nacov[[g]], # symmetric=FALSE,
                  only.values = TRUE
                )$values
                if (is.complex(ev)) {
                  lav_msg_stop(gettext(
                    "Gamma (NACOV) matrix is not positive-definite"))
                }
                if (any(Re(ev) < 0)) {
                  lav_msg_stop(gettext(
                    "Gamma (NACOV) matrix is not positive-definite"))
                }
              }
              if (estimator == "DLS" && dls_gamma_nt == "sample") {
                if (dls_a == 1.0) {
                  # nothing to do, use GLS version
                } else {
                  w_dls <-
                    (1 - dls_a) * nacov[[g]] + dls_a * gamma_nt
                  wls_v[[g]] <-
                    lav_mat_sym_inverse(w_dls)
                }
              } else { # WLS
                wls_v[[g]] <-
                  lav_mat_sym_inverse(nacov[[g]])
              }
            } else {
              # fixed.x: we have zero cols/rows
              # ginv does the trick, but perhaps this is overkill
              # just removing the zero rows/cols, invert, and
              # fill back in the zero rows/cols would do it
              # wls_v[[g]] <- MASS::ginv(nacov[[g]])
              if (estimator == "DLS" && dls_gamma_nt == "sample") {
                w_dls <- (1 - dls_a) * nacov[[g]] + dls_a * gamma_nt
                wls_v[[g]] <-
                  lav_mat_sym_inverse(w_dls)
              } else { # WLS
                wls_v[[g]] <-
                  lav_mat_sym_inverse(nacov[[g]])
              }
            }
          } else if (estimator == "DWLS") {
            dacov <- diag(nacov[[g]])
            if (!all(is.finite(dacov))) {
              lav_msg_stop(gettext(
                "diagonal of Gamma (NACOV) contains non finite values"))
            }
            if (fixed_x) {
              # structural zeroes!
              zero_idx <- which(dacov == 0.0)
              idacov <- 1 / dacov
              idacov[zero_idx] <- 0.0
            } else {
              idacov <- 1 / dacov
            }
            wls_v[[g]] <- diag(idacov,
              nrow = NROW(nacov[[g]]),
              ncol = NCOL(nacov[[g]])
            )
            wls_vd[[g]] <- idacov
          } else if (estimator == "ULS") {
            # wls_v[[g]] <- diag(length(WLS.obs[[g]]))
            wls_vd[[g]] <- rep(1, length(wls_obs[[g]]))
          }
        } else {
          if (estimator == "WLS") {
            c_s <- tryCatch(chol(cat_1$WLS.W * nobs[[g]]),
                           error = function(e) NULL)
            if (is.null(c_s)) {
              lav_msg_stop(gettext("could not invert CAT$WLS.W"))
            }
            wls_v[[g]] <- chol2inv(c_s)
          } else if (estimator == "DWLS") {
            dacov <- diag(cat_1$WLS.W * nobs[[g]])
            # wls_v[[g]] <- diag(1/dacov, nrow=NROW(CAT$WLS.W),
            #                            ncol=NCOL(CAT$WLS.W))
            wls_vd[[g]] <- 1 / dacov
          } else if (estimator == "ULS" || estimator == "IV") {
            # wls_v[[g]] <- diag(length(WLS.obs[[g]]))
            wls_vd[[g]] <- rep(1, length(wls_obs[[g]]))
          }
        }
      } else if (estimator == "PML" || estimator == "FML") {
        # no wls_v here
      }

      # group.w.free (only if categorical)
      if (group_w_free && categorical) {
        if (!is.null(wls_v[[g]])) {
          # unweight!!
          a <- group_w[[g]] * sum(unlist(nobs)) / nobs[[g]]
          # always 1!!!
          # invert
          a <- 1 / a
          wls_v[[g]] <- lav_mat_bdiag(matrix(a, 1, 1), wls_v[[g]])
        }
        if (!is.null(wls_vd[[g]])) {
          # unweight!!
          a <- group_w[[g]] * sum(unlist(nobs)) / nobs[[g]]
          # always 1!!!
          # invert
          a <- 1 / a
          wls_vd[[g]] <- c(a, wls_vd[[g]])
        }
      }
    }
  } # ngroups

  # remove 'CAT', unless debug -- this is to save memory
  if (!lav_debug()) {
    cat_1 <- list()
  }

  # construct SampleStats object
  lav_sample_stats <- new("lavSampleStats",
    # sample moments
    th = th,
    th.idx = th_idx,
    th.names = th_names,
    mean = mean,
    cov = cov,
    var = var,

    # residual (y | x)
    res.cov = res_cov,
    res.var = res_var,
    res.th = res_th,
    res.th.nox = res_th_nox,
    res.slopes = res_slopes,
    res.int = res_int,
    mean.x = mean_x,
    cov.x = cov_x,
    bifreq = bifreq,
    group.w = group_w,

    # convenience
    nobs = nobs,
    ntotal = sum(unlist(nobs)),
    ngroups = ngroups,
    x.idx = x_idx,

    # extra sample statistics
    icov = icov,
    cov.log.det = cov_log_det,
    res.icov = res_icov,
    res.cov.log.det = res_cov_log_det,
    ridge = ridge_eps,
    WLS.obs = wls_obs,
    WLS.V = wls_v,
    WLS.VD = wls_vd,
    NACOV = nacov,
    NACOV.user = nacov_user,

    # cluster/levels
    YLp = ylp,

    # missingness
    missing.flag = missing_flag,
    missing = missing_1,
    missing.h1 = missing_h1,
    zero.cell.tables = zero_cell_tables
  )

  # just a SINGLE warning if we have empty cells
  if ((!is.null(lavoptions$samplestats) && lavoptions$samplestats) &&
      categorical && zero_cell_warn &&
      any(sapply(zero_cell_tables, nrow) > 0L)) {
    nempty <- sum(sapply(zero_cell_tables, nrow))
    lav_msg_warn(gettextf(
      "%s bivariate tables have empty cells; to see them, use:
      lavInspect(fit, \"zero.cell.tables\")", nempty)
    )
  }

  lav_sample_stats
}


lav_samp_from_moments <- function(sample_cov = NULL,
                                         sample_mean = NULL,
                                         sample_th = NULL,
                                         sample_nobs = NULL,
                                         ov_names = NULL, # including x
                                         ov_names_x = NULL,
                                         wls_v = NULL,
                                         nacov = NULL,
                                         lavoptions = NULL) {
  # extract options
  estimator <- lavoptions$estimator
  # mimic <- lavoptions$mimic
  meanstructure <- lavoptions$meanstructure
  correlation <- lavoptions$correlation
  group_w_free <- lavoptions$group.w.free
  ridge <- lavoptions$ridge
  rescale <- lavoptions$sample.cov.rescale

  # no multilevel yet
  nlevels <- 1L

  # ridge default
  if (ridge) {
    if (is.numeric(lavoptions$ridge.constant)) {
      ridge_eps <- lavoptions$ridge.constant
    } else {
      ridge_eps <- 1e-5
    }
  } else {
    ridge_eps <- 0.0
  }

  # new in 0.6-3:
  # check if sample.cov has attributes if conditional.x = TRUE
  sample_res_slopes <- attr(sample_cov, "res.slopes")
  sample_cov_x <- attr(sample_cov, "cov.x")
  sample_mean_x <- attr(sample_cov, "mean.x")
  if (!is.null(sample_res_slopes)) {
    conditional_x <- TRUE
    # strip attributes
    attr(sample_cov, "res.slopes") <- NULL
    attr(sample_cov, "cov.x") <- NULL
    attr(sample_cov, "mean.x") <- NULL
    # make list
    if (!is.list(sample_res_slopes)) {
      sample_res_slopes <- list(sample_res_slopes)
    }
    if (!is.list(sample_cov_x)) {
      sample_cov_x <- list(sample_cov_x)
    }
    if (!is.list(sample_mean_x)) {
      sample_mean_x <- list(sample_mean_x)
    }
  } else if (!is.null(sample_cov_x)) {
    conditional_x <- FALSE
    fixed_x <- TRUE

    # strip attributes
    attr(sample_cov, "cov.x") <- NULL
    attr(sample_cov, "mean.x") <- NULL
    # make list
    if (!is.list(sample_cov_x)) {
      sample_cov_x <- list(sample_cov_x)
    }
    if (!is.list(sample_mean_x)) {
      sample_mean_x <- list(sample_mean_x)
    }
  } else if (is.null(sample_cov_x) && length(unlist(ov_names_x)) > 0L) {
    # fixed.x = TRUE, but only joint sample.cov is provided
    conditional_x <- FALSE
    fixed_x <- TRUE

    # create sample.cov.x and sample.mean.x later...
  } else {
    conditional_x <- FALSE
    fixed_x <- FALSE
  }

  # matrix -> list
  if (!is.list(sample_cov)) {
    sample_cov <- list(sample_cov)
  }

  # number of groups
  ngroups <- length(sample_cov)

  # ov.names
  if (!is.list(ov_names)) {
    ov_names <- rep(list(ov_names), ngroups)
  }
  if (!is.list(ov_names_x)) {
    ov_names_x <- rep(list(ov_names_x), ngroups)
  }

  if (!is.null(sample_mean)) {
    meanstructure <- TRUE
    if (!is.list(sample_mean)) {
      # check if sample.mean is string (between single quotes)
      if (is.character(sample_mean)) {
        sample_mean <- lav_char2num(sample_mean)
      }
      sample_mean <- list(unname(sample_mean))
    } else {
      sample_mean <- lapply(lapply(sample_mean, unname), unclass)
    }
  }

  if (!is.null(sample_th)) {
    th_idx <- attr(sample_th, "th.idx")
    attr(sample_th, "th.idx") <- NULL
    if (is.null(th_idx)) {
      lav_msg_stop(gettext("sample.th should have a th.idx attribute"))
    } else {
      if (is.list(th_idx)) {
        th_names <- lapply(th_idx, names)
        th_idx <- lapply(lapply(th_idx, unname), unclass)
      } else {
        th_names <- list(names(th_idx))
        th_idx <- list(unclass(unname(th_idx)))
      }
    }
    if (is.list(sample_th)) {
      # strip names and lavaan.vector class
      sample_th <- lapply(lapply(sample_th, unname), unclass)
    } else {
      # strip names and lavaan.vector class, make list
      sample_th <- list(unclass(unname(sample_th)))
    }
  } else {
    th_idx <- vector("list", length = ngroups)
    th_names <- vector("list", length = ngroups)
  }

  # sample statistics per group
  cov <- vector("list", length = ngroups)
  var <- vector("list", length = ngroups)
  mean <- vector("list", length = ngroups)
  th <- vector("list", length = ngroups)
  # th.idx      <- vector("list", length = ngroups)
  # th.names    <- vector("list", length = ngroups)

  # residual (y | x)
  res_cov <- vector("list", length = ngroups)
  res_var <- vector("list", length = ngroups)
  res_slopes <- vector("list", length = ngroups)
  res_int <- vector("list", length = ngroups)
  res_th <- vector("list", length = ngroups)
  res_th_nox <- vector("list", length = ngroups)

  # fixed.x / conditional.x
  mean_x <- vector("list", length = ngroups)
  cov_x <- vector("list", length = ngroups)

  bifreq <- vector("list", length = ngroups)

  # extra sample statistics per group
  icov <- vector("list", length = ngroups)
  cov_log_det <- vector("list", length = ngroups)
  res_icov <- vector("list", length = ngroups)
  res_cov_log_det <- vector("list", length = ngroups)
  wls_obs <- vector("list", length = ngroups)
  missing_1 <- vector("list", length = ngroups)
  missing_h1 <- vector("list", length = ngroups)
  missing_flag <- FALSE
  zero_cell_tables <- vector("list", length = ngroups)
  ylp <- vector("list", length = ngroups)

  # group weights
  group_w <- vector("list", length = ngroups)
  x_idx <- vector("list", length = ngroups)

  categorical <- FALSE
  if (!is.null(sample_th)) {
    categorical <- TRUE
  }

  wls_vd <- vector("list", length = ngroups)
  if (is.null(wls_v)) {
    wls_v <- vector("list", length = ngroups)
    wls_v_user <- FALSE
  } else {
    if (!is.list(wls_v)) {
      if (ngroups == 1L) {
        wls_v <- list(unclass(wls_v))
      } else {
        lav_msg_stop(gettextf("wls_v argument should be a list of length %s",
          ngroups)
        )
      }
    } else {
      if (length(wls_v) != ngroups) {
        lav_msg_stop(gettextf(
          "wls_v assumes %1$s groups; data contains %2$s groups",
          length(wls_v), ngroups))
      }
      wls_v <- lapply(wls_v, unclass)
    }

    # is wls_v full? check first
    if (is.null(dim(wls_v[[1]]))) {
      # we will assume it is the diagonal only
      wls_vd <- wls_v
      wls_v <- lapply(wls_vd, diag)
    } else {
      # create WLS.VD
      wls_vd <- lapply(wls_v, diag)
      # we could remove wls_v to save space...
    }

    wls_v_user <- TRUE
    # FIXME: check dimension of wls_v!!
  }

  if (is.null(nacov)) {
    nacov <- vector("list", length = ngroups)
    nacov_user <- FALSE
  } else {
    if (!is.list(nacov)) {
      if (ngroups == 1L) {
        nacov <- list(unclass(nacov))
      } else {
        lav_msg_stop(gettextf(
          "NACOV argument should be a list of length %s", ngroups))
      }
    } else {
      if (length(nacov) != ngroups) {
        lav_msg_stop(gettextf(
          "NACOV assumes %1$s groups; data contains %2$s groups",
          length(nacov), ngroups))
      }
      nacov <- lapply(nacov, unclass)
    }
    nacov_user <- TRUE
    # FIXME: check dimension of NACOV!!
  }

  nobs <- as.list(as.integer(sample_nobs))


  for (g in 1:ngroups) {
    # exogenous x?
    nexo <- length(ov_names_x[[g]])
    if (nexo) {
      # two cases: ov.names contains 'x' variables, or not
      if (conditional_x) {
        # ov.names.x are NOT in ov.names
        x_idx[[g]] <- which(ov_names[[g]] %in% ov_names_x[[g]])
      } else {
        if (fixed_x) {
          # ov.names.x are a subset of ov.names
          x_idx[[g]] <- match(ov_names_x[[g]], ov_names[[g]])
          stopifnot(!anyNA(x_idx[[g]]))
        } else {
          x_idx[[g]] <- integer(0L)
        }
      }
    } else {
      x_idx[[g]] <- integer(0L)
      conditional_x <- FALSE
      fixed_x <- FALSE
    }


    # group weight
    group_w[[g]] <- nobs[[g]] / sum(unlist(nobs))

    tmp_cov <- sample_cov[[g]]

    # make sure that the matrix is fully symmetric (NEEDED?)
    t_1 <- t(tmp_cov)
    tmp_cov[upper.tri(tmp_cov)] <- t_1[upper.tri(t_1)]

    # check dimnames
    if (!is.null(rownames(tmp_cov))) {
      cov_names <- rownames(tmp_cov)
    } else if (!is.null(colnames(tmp_cov))) {
      cov_names <- colnames(tmp_cov)
    } else {
      lav_msg_stop(gettext(
        "please provide row/col names for the covariance matrix!"))
    }

    # extract only the part we need (using ov.names)
    if (conditional_x) {
      idx <- match(ov_names[[g]][-x_idx[[g]]], cov_names)
    } else {
      idx <- match(ov_names[[g]], cov_names)
    }
    if (any(is.na(idx))) {
      cat("found: ", cov_names, "\n")
      cat("expected: ", ov_names[[g]], "\n")
      lav_msg_stop(gettextf(
        "rownames of covariance matrix do not match the model!
        found: %1$s expected: %2$s", lav_msg_view(cov_names),
        lav_msg_view(ov_names[[g]])))
    } else {
      tmp_cov <- tmp_cov[idx, idx, drop = FALSE]
    }

    # strip dimnames
    dimnames(tmp_cov) <- NULL

    if (is.null(sample_mean)) {
      # assume zero mean vector
      tmp_mean <- numeric(ncol(tmp_cov))
    } else {
      # extract only the part we need
      tmp_mean <- unclass(sample_mean[[g]][idx])
    }



    if (categorical) {
      # categorical + conditional.x = TRUE
      if (conditional_x) {
        th_g <- numeric(length(th_idx[[g]]))
        ord_idx <- which(th_idx[[g]] > 0)
        num_idx <- which(th_idx[[g]] == 0)
        if (length(ord_idx) > 0L) {
          th_g[ord_idx] <- sample_th[[g]][ord_idx]
        }
        if (length(num_idx) > 0L) {
          ord_var_idx <- unique(th_idx[[g]][th_idx[[g]] > 0])
          th_g[num_idx] <- -1 * sample_mean[[g]][-ord_var_idx]
        }
        res_th[[g]] <- th_g
        res_th_nox[[g]] <- sample_th[[g]]

        res_cov[[g]] <- tmp_cov
        if (ridge) {
          diag(res_cov[[g]]) <- diag(res_cov[[g]]) + ridge_eps
        }
        res_var[[g]] <- diag(tmp_cov)
        res_int[[g]] <- tmp_mean

        res_slopes[[g]] <- unclass(unname(sample_res_slopes[[g]]))
        cov_x[[g]] <- unclass(unname(sample_cov_x[[g]]))
        mean_x[[g]] <- unclass(unname(sample_mean_x[[g]]))

        # th.idx and th.names are already ok

        # categorical + conditional.x = FALSE
      } else {
        th_g <- numeric(length(th_idx[[g]]))
        ord_idx <- which(th_idx[[g]] > 0)
        num_idx <- which(th_idx[[g]] == 0)
        if (length(ord_idx) > 0L) {
          th_g[ord_idx] <- sample_th[[g]][ord_idx]
        }
        if (length(num_idx) > 0L) {
          ord_var_idx <- unique(th_idx[[g]][th_idx[[g]] > 0])
          th_g[num_idx] <- -1 * sample_mean[[g]][-ord_var_idx]
        }
        th[[g]] <- th_g

        cov[[g]] <- tmp_cov
        if (ridge) {
          diag(cov[[g]]) <- diag(cov[[g]]) + ridge_eps
        }
        var[[g]] <- diag(tmp_cov)
        mean[[g]] <- tmp_mean

        # fixed.x? (needed?)
        if (fixed_x) {
          cov_x[[g]] <- unclass(unname(sample_cov_x[[g]]))
          mean_x[[g]] <- unclass(unname(sample_mean_x[[g]]))
        }

        # th, th.idx and th.names are already ok
      }

      # multilevel
    } else if (nlevels > 1L) {
      lav_msg_stop(gettext("multilevel + sample stats not ready yet"))


      # single level
    } else {
      # single-level + continuous + conditional.x = TRUE
      if (conditional_x) {
        res_cov[[g]] <- tmp_cov
        if (ridge) {
          diag(res_cov[[g]]) <- diag(res_cov[[g]]) + ridge_eps
        }
        res_var[[g]] <- diag(tmp_cov)
        res_int[[g]] <- tmp_mean
        res_slopes[[g]] <- unclass(unname(sample_res_slopes[[g]]))
        cov_x[[g]] <- unclass(unname(sample_cov_x[[g]]))
        mean_x[[g]] <- unclass(unname(sample_mean_x[[g]]))

        # no rescale!

        # icov and cov.log.det
        # if(lavoptions$sample.icov) {
        out <- lav_samp_icov(
          cov_1 = res_cov[[g]],
          ridge = 1e-05,
          x_idx = x_idx[[g]],
          ngroups = ngroups, g = g
        )
        res_icov[[g]] <- out$icov
        res_cov_log_det[[g]] <- out$cov.log.det
        # }

        # continuous + conditional.x = FALSE
      } else {
        cov[[g]] <- tmp_cov
        mean[[g]] <- tmp_mean

        # rescale cov by (N-1)/N?
        if (rescale) {
          # we 'transform' the sample cov (divided by n-1)
          # to a sample cov divided by 'n'
          cov[[g]] <- ((nobs[[g]] - 1) / nobs[[g]]) * cov[[g]]
        }
        if (ridge) {
          diag(cov[[g]]) <- diag(cov[[g]]) + ridge_eps
        }
        var[[g]] <- diag(cov[[g]])

        # icov and cov.log.det
        # if(lavoptions$sample.icov) {
        out <- lav_samp_icov(
          cov_1 = cov[[g]],
          ridge = 1e-05,
          x_idx = x_idx[[g]],
          ngroups = ngroups,
          g = g
        )
        icov[[g]] <- out$icov
        cov_log_det[[g]] <- out$cov.log.det
        # }

        # fixed.x?
        if (fixed_x) {
          if (is.null(sample_cov_x)) {
            cov_x[[g]] <- cov[[g]][x_idx[[g]], x_idx[[g]],
              drop = FALSE
            ]
          } else {
            cov_x[[g]] <- unclass(unname(sample_cov_x[[g]]))
          }
          if (is.null(sample_mean_x)) {
            mean_x[[g]] <- mean[[g]][x_idx[[g]]]
          } else {
            mean_x[[g]] <- unclass(unname(sample_mean_x[[g]]))
          }
        }
      }

      # correlation structure?
      if (correlation) {
        cov[[g]] <- cov2cor(cov[[g]])
        var[[g]] <- rep(1, length(var[[g]]))
        if (conditional_x) {
          res_cov[[g]] <- cov2cor(res_cov[[g]])
          res_var[[g]] <- rep(1, length(res_var[[g]]))
          cov_x[[g]] <- cov2cor(cov_x[[g]])
          # FIXME: slopes? more?
        }
      }
    }

    # WLS.obs
    wls_obs[[g]] <- lav_samp_wls_obs(
      mean_g = mean[[g]],
      cov_g = cov[[g]], var_g = var[[g]], th_g = th[[g]],
      th_idx_g = th_idx[[g]], res_int_g = res_int[[g]],
      res_cov_g = res_cov[[g]], res_var_g = res_var[[g]],
      res_th_g = res_th[[g]], res_slopes_g = res_slopes[[g]],
      group_w_g = log(nobs[[g]]),
      categorical = categorical, conditional_x = conditional_x,
      meanstructure = meanstructure, correlation = correlation,
      slopestructure = conditional_x,
      group_w_free = group_w_free
    )

    # wls_v
    if (!wls_v_user) {
      if (estimator == "GLS") {
        # FIXME: in <0.5-21, we had
        # V11 <- icov[[g]]
        #    if(mimic == "Mplus") { # is this a bug in Mplus?
        #        V11 <- V11 * nobs[[g]]/(nobs[[g]]-1)
        #    }
        if (correlation) {
          gamma_nt <- lav_samp_cor_gamma_nt(
            m_cov          = cov[[g]],
            m_mean         = mean[[g]],
            #rescale        = FALSE,
            x_idx          = x_idx[[g]],     # not used yet
            fixed_x        = fixed_x,        # not used yet
            conditional_x  = conditional_x,  # not used yet
            meanstructure  = meanstructure,  # not used yet
            slopestructure = conditional_x   # not used yet
          )
          wls_v[[g]] <- lav_mat_sym_inverse(gamma_nt)
        } else {
          wls_v[[g]] <- lav_samp_gamma_inverse_nt(
            m_icov = icov[[g]],
            m_cov = cov[[g]],
            m_mean = mean[[g]],
            rescale = FALSE,
            x_idx = x_idx[[g]],
            fixed_x = fixed_x,
            conditional_x = conditional_x,
            meanstructure = meanstructure,
            slopestructure = conditional_x
          )
        }
      } else if (estimator == "ULS") {
        wls_v[[g]] <- diag(length(wls_obs[[g]]))
        wls_vd[[g]] <- rep(1, length(wls_obs[[g]]))
      } else if (estimator == "WLS" || estimator == "DWLS") {
        if (is.null(wls_v[[g]])) {
          lav_msg_stop(gettext(
            "the (D)WLS estimator is only available with full data
            or with a user-provided wls_v"))
        }
      }

      # group.w.free
      if (!is.null(wls_v[[g]]) && group_w_free) {
        # FIXME!!!
        wls_v[[g]] <- lav_mat_bdiag(matrix(1, 1, 1), wls_v[[g]])
      }
    }
  } # ngroups

  # construct SampleStats object
  lav_sample_stats <- new("lavSampleStats",
    # sample moments
    th = th,
    th.idx = th_idx,
    th.names = th_names,
    mean = mean,
    cov = cov,
    var = var,

    # residual (y | x)
    res.cov = res_cov,
    res.var = res_var,
    res.th = res_th,
    res.th.nox = res_th_nox,
    res.slopes = res_slopes,
    res.int = res_int,

    # fixed.x
    mean.x = mean_x,
    cov.x = cov_x,

    # other
    bifreq = bifreq,
    group.w = group_w,

    # convenience
    nobs = nobs,
    ntotal = sum(unlist(nobs)),
    ngroups = ngroups,
    x.idx = x_idx,

    # extra sample statistics
    icov = icov,
    cov.log.det = cov_log_det,
    res.icov = res_icov,
    res.cov.log.det = res_cov_log_det,
    ridge = ridge_eps,
    WLS.obs = wls_obs,
    WLS.V = wls_v,
    WLS.VD = wls_vd,
    NACOV = nacov,
    NACOV.user = nacov_user,

    # cluster/level
    YLp = ylp,

    # missingness
    missing.flag = missing_flag,
    missing = missing_1,
    missing.h1 = missing_h1,
    zero.cell.tables = zero_cell_tables
  )

  lav_sample_stats
}

# compute sample statistics, per missing pattern
lav_samp_mi_patterns <- function(
  y = NULL, mp = NULL, wt = NULL,lp = NULL) {
  # coerce Y to matrix
  y <- as.matrix(y)

  # handle two-level data
  if (!is.null(lp)) {
    y_orig <- y
    # z <- NULL
    if (length(lp$between.idx[[2]]) > 0L) {
      y <- y[, -lp$between.idx[[2]], drop = FALSE]
      # z_idx <- which(!duplicated(lp$cluster.idx[[2]]))
      # z <- y_orig[z_idx, lp$between.idx[[2]], drop = FALSE]
    }
  }

  if (is.null(mp)) {
    mp <- lav_data_mi_patterns(y,
      sort_freq = FALSE, coverage = FALSE,
      lp = lp
    )
  }

  yp <- vector("list", length = mp$npatterns)

  # fill in pattern statistics
  for (p in seq_len(mp$npatterns)) {
    # extract raw data for these cases
    raw_1 <- y[mp$case.idx[[p]], mp$pat[p, ], drop = FALSE]

    # more than one case
    if (mp$freq[p] > 1L) {
      if (!is.null(wt)) {
        out <- stats::cov.wt(raw_1,
          wt = wt[mp$case.idx[[p]]],
          method = "ML"
        )
        sy <- out$cov
        my <- out$center
      } else {
        my <- base::.colMeans(raw_1, m = NROW(raw_1), n = NCOL(raw_1))
        # SY <- crossprod(RAW)/Mp$freq[p] - tcrossprod(MY)
        # bad practice, better like this:
        sy <- lav_mat_cov(raw_1)
      }
    # only a single observation (no need to weight!)
    } else {
      sy <- 0
      my <- as.numeric(raw_1)
    }

    if (!is.null(wt)) {
      freq <- sum(wt[mp$case.idx[[p]]])
    } else {
      freq <- mp$freq[p]
    }

    # store sample statistics, var.idx and freq
    yp[[p]] <- list(
      SY = sy, MY = my, var.idx = mp$pat[p, ],
      freq = freq
    )

    # if clustered data, add rowsum over all cases per cluster
    if (!is.null(lp)) {
      tmp <- rowsum.default(raw_1, group = mp$j.idx[[p]], reorder = FALSE)
      yp[[p]]$ROWSUM <- tmp
    }
  }

  # add Zp as an attribute
  # if(!is.null(Lp)) {
  #    Zp <- lav_samp_mi_patterns(Y = Z, Mp = Mp$Zp)
  #    for(p in Mp$Zp$npatterns) {
  #        this.z <- Z[Mp$Zp$case.idx[[p]], drop = FALSE]
  #        Zp[[p]]$ROWSUM <- t(this.z)
  #
  #    }
  #    attr(Yp, "Zp") <- Zp
  # }

  yp
}

# compute sample statistics, per cluster
lav_samp_cl_patterns <- function(y = NULL, lp = NULL,
                                             conditional_x = FALSE) {
  # coerce Y to matrix
  y1 <- as.matrix(y)
  n <- NROW(y1)
  p <- NCOL(y1)

  if (is.null(lp)) {
    lav_msg_stop(gettext("Lp is NULL"))
  }

  # how many levels?
  nlevels <- length(lp$cluster) + 1L

  # compute some sample statistics per level
  ylp <- vector("list", length = nlevels)
  for (l in 2:nlevels) {
    ncluster_sizes <- lp$ncluster.sizes[[l]]
    cluster_size <- lp$cluster.size[[l]]
    cluster_sizes <- lp$cluster.sizes[[l]]
    nclusters <- lp$nclusters[[l]]
    both_idx <- lp$both.idx[[l]]
    within_idx <- lp$within.idx[[l]]
    between_idx <- lp$between.idx[[l]]
    cluster_idx <- lp$cluster.idx[[l]]
    # cluster_size_ns <- lp$cluster.size.ns[[l]]

    # s <- (N^2 - sum(cluster.size^2)) / (N*(nclusters - 1L))
    # same as
    s <- (n - sum(cluster_size^2) / n) / (nclusters - 1)
    # NOTE: must be (nclusters - 1), otherwise, s is not average cluster
    # size even in the balanced case

    y1_means <- colMeans(y1, na.rm = TRUE)
    y1y1 <- lav_mat_crossprod(y1)
    both_idx <- all_idx <- seq_len(p)

    if (length(within_idx) > 0L ||
      length(between_idx) > 0L) {
      both_idx <- all_idx[-c(within_idx, between_idx)]
      # hm, this assumes the 'order' is the
      # same at both levels...
    }

    # cluster-means
    y2 <- rowsum.default(y1,
      group = cluster_idx, reorder = FALSE,
      na.rm = FALSE, # must be FALSE!
    ) / cluster_size
    y2c <- t(t(y2) - y1_means)

    # compute S.w
    # center within variables by grand mean instead of group mean?
    # (YR: apparently not for S.PW)

    y2a <- y2
    # if(length(within.idx) > 0L) {
    #    for(i in 1:length(within.idx)) {
    #        Y2a[, within.idx[i]] <- Y1.means[within.idx[i]]
    #    }
    # }
    y1a <- y1 - y2a[cluster_idx, , drop = FALSE]
    s_w <- lav_mat_crossprod(y1a) / (n - nclusters)

    # S.b
    # three parts: within/within, between/between, between/within
    # standard definition of the between variance matrix
    # divides by (nclusters - 1)
    s_b <- lav_mat_crossprod(y2c * cluster_size, y2c) / (nclusters - 1)

    # check for zero variances
    if (length(both_idx) > 0L) {
      zero_idx <- which(diag(s_b)[both_idx] < 0.0001)
      if (length(zero_idx) > 0L && !anyNA(y2)) {
        lav_msg_warn(gettext(
          "(near) zero variance at between level for split variable:"),
          paste(lp$both.names[[l]][zero_idx], collapse = " ")
        )
      }
    }

    s_1 <- cov(y1, use = "pairwise.complete.obs") * (n - 1L) / n
  # missing by design?
  s_1[is.na(s_1)] <- as.numeric(NA)

    # loglik.x
    # extract 'fixed' level-1 loglik from here
    wx_idx <- lp$ov.x.idx[[1]]
    if (length(wx_idx) > 0L) {
      loglik_x_w <- lav_mvn_h1_loglik_samp(
        sample_nobs = lp$nclusters[[1]],
        sample_cov  = s_1[wx_idx, wx_idx, drop = FALSE]
      )
    } else {
      loglik_x_w <- 0
    }
    # extract 'fixed' level-2 loglik
    bx_idx <- lp$ov.x.idx[[2]]
    if (length(bx_idx) > 0L) {
      covb <- cov(y2[, bx_idx, drop = FALSE]) * (nclusters - 1) / nclusters
      loglik_x_b <- lav_mvn_h1_loglik_samp(
        sample_nobs = lp$nclusters[[2]],
        sample_cov  = covb
      )
    } else {
      loglik_x_b <- 0
    }
    loglik_x <- loglik_x_w + loglik_x_b


    s_pw_start <- s_w
    if (length(within_idx) > 0L) {
      s_pw_start[within_idx, within_idx] <-
        s_1[within_idx, within_idx, drop = FALSE]
    }

    if (length(between_idx) > 0L) {
      s_w[between_idx, ] <- 0
      s_w[, between_idx] <- 0
      s_pw_start[between_idx, ] <- 0
      s_pw_start[, between_idx] <- 0
    }

    if (length(between_idx) > 0L) {
      # this is what is needed for MUML:
      s_b[, between_idx] <-
        (s * nclusters / n) * s_b[, between_idx, drop = FALSE]
      s_b[between_idx, ] <-
        (s * nclusters / n) * s_b[between_idx, , drop = FALSE]
      s_b[between_idx, between_idx] <-
        (s * lav_mat_crossprod(
          y2c[, between_idx, drop = FALSE],
          y2c[, between_idx, drop = FALSE]
        ) / nclusters)
    }

    sigma_b <- (s_b - s_w) / s
    sigma_b[within_idx, ] <- 0
    sigma_b[, within_idx] <- 0

    # what if we have negative variances in Sigma.B?
    # this may happen if 'split' a variable that has no between variance
    zero_idx <- which(diag(sigma_b) < 1e-10)
    if (length(zero_idx) > 0L) {
      sigma_b[zero_idx, ] <- 0
      sigma_b[, zero_idx] <- 0
    }


    mu_w <- numeric(p)
    mu_w[within_idx] <- y1_means[within_idx]

    mu_b <- y1_means
    mu_b[within_idx] <- 0
    if (length(between_idx) > 0L) {
      # replace between.idx by cov(Y2)[,] elements...
      mu_b[between_idx] <- colMeans(y2[, between_idx, drop = FALSE],
        na.rm = TRUE
      )

      s2 <- (cov(y2, use = "pairwise.complete.obs") *
        (nclusters - 1L) / nclusters)

      sigma_b[between_idx, between_idx] <-
        s2[between_idx, between_idx, drop = FALSE]
    }

    # FIXME: Mu.B not quite ok for (fixed.x) x variables if they
    # occur both at level 1 AND level 2
    mu_b_start <- mu_b
    # Mu.B.start[both.idx] <- Mu.B.start[both.idx] - colMeans(Y2c[,both.idx])

    # sample statistics PER CLUSTER-SIZE

    # summary statistics for complete data, conditional.x = FALSE
    # also needed for h1 (even if conditional.x = TRUE)
    cov_d <- vector("list", length = ncluster_sizes)
    mean_d <- vector("list", length = ncluster_sizes)
    for (clz in seq_len(ncluster_sizes)) {
      nj <- cluster_sizes[clz]
      # select clusters with this size
      d_idx <- which(cluster_size == nj)
      ns <- length(d_idx)
      # NOTE:!!!!
      # reorder columns
      # to match A.inv and m.k later on in objective!!!
      tmp2 <- y2[d_idx,
        c(between_idx, sort.int(c(both_idx, within_idx))),
        drop = FALSE
      ]
      mean_d[[clz]] <- colMeans(tmp2, na.rm = TRUE)
      bad_idx <- which(!is.finite(mean_d[[clz]])) # if nrow = 1 + NA
      if (length(bad_idx) > 0L) {
        mean_d[[clz]][bad_idx] <- 0 # ugly, only for starting values
      }
      if (length(d_idx) > 1L) {
        if (any(is.na(tmp2))) {
          # if full column has NA, this will fail...
          # not needed anyway
          # out <- lav_mvn_mi_h1_est_moments(Y = tmp2,
          #          max.iter = 10L)
          # cov.d[[clz]] <- out$Sigma
          cov_d[[clz]] <- 0
        } else {
          cov_d[[clz]] <- (cov(tmp2, use = "complete.obs") *
            (ns - 1) / ns)
        }
      } else {
        cov_d[[clz]] <- 0
      }
    } # clz

    # new in 0.6-12:
    # summary statistics for complete data, conditional.x = TRUE
    # ONLY for twolevel
    if (conditional_x) {
      within_x_idx <- lp$within.x.idx[[1]]
      between_y_idx <- lp$between.y.idx[[2]]
      between_x_idx <- lp$between.x.idx[[2]]
      y1_idx <- lp$ov.y.idx[[1]]
      x1_idx <- c(within_x_idx, between_x_idx) # in that order

      # data
      y1_wb <- y1[, y1_idx, drop = FALSE]
      y2_wb <- y2[, y1_idx, drop = FALSE]
      # if (length(between_y_idx) > 0L) {
      #   y2_z <- y2[, between_y_idx, drop = FALSE]
      # }
      if (length(x1_idx) > 0L) {
        exo_wb1 <- cbind(1, y1[, x1_idx, drop = FALSE])
        exo_wb2 <- cbind(1, y2[, x1_idx, drop = FALSE])
      } else {
        exo_wb1 <- matrix(1, nrow(y1), 1L)
        exo_wb2 <- matrix(1, nrow(y2), 1L)
      }

      # sample beta.wb (level 1)
      sample_wb <- solve(crossprod(exo_wb1), crossprod(exo_wb1, y1_wb))
      sample_yhat_wb1 <- exo_wb1 %*% sample_wb
      sample_yres_wb1 <- y1_wb - sample_yhat_wb1
      sample_yyres_wb1 <- crossprod(sample_yres_wb1)
      sample_xx_wb1 <- crossprod(exo_wb1)

      # sample beta.wb (level 2)
      xx_wb2 <- crossprod(exo_wb2)
      sample_wb2 <- try(solve(xx_wb2, crossprod(exo_wb2, y2_wb)),
        silent = TRUE
      )
      if (inherits(sample_wb2, "try-error")) {
        # this may happen if the covariate is cluster-centered
        # using the observed cluster means; then the 'means' will
        # be all (near) zero, and there is no variance
        sample_wb2 <- MASS::ginv(xx_wb2) %*% crossprod(exo_wb2, y2_wb)
      }
      sample_yhat_wb2 <- exo_wb2 %*% sample_wb2
      sample_yres_wb2 <- y2_wb - sample_yhat_wb2

      # weighted by cluster.size
      sample_yyres_wb2 <- crossprod(
        sample_yres_wb2,
        sample_yres_wb2 * cluster_size
      )
      sample_yres_x_wb2 <- crossprod(
        sample_yres_wb2,
        exo_wb2 * cluster_size
      )
      sample_xx_wb2 <- crossprod(
        exo_wb2,
        exo_wb2 * cluster_size
      )


      sample_clz_y2_res <- vector("list", ncluster_sizes)
      sample_clz_y2_xx <- vector("list", ncluster_sizes)
      sample_clz_y2_b <- vector("list", ncluster_sizes)

      if (length(between_y_idx) > 0L) {
        sample_clz_zz_res <- vector("list", ncluster_sizes)
        sample_clz_zz_xx <- vector("list", ncluster_sizes)
        sample_clz_zz_b <- vector("list", ncluster_sizes)

        sample_clz_yz_res <- vector("list", ncluster_sizes)
        sample_clz_yz_xx <- vector("list", ncluster_sizes)
        sample_clz_yres_xz <- vector("list", ncluster_sizes)
        sample_clz_xwzres <- vector("list", ncluster_sizes)
      }
      for (clz in seq_len(ncluster_sizes)) {
        # cluster size
        nj <- cluster_sizes[clz]
        nj_idx <- which(cluster_size == nj)

        # Y2
        y2_clz <- y2[nj_idx, y1_idx, drop = FALSE]
        if (length(x1_idx) > 0L) {
          exo2_clz <- cbind(1, y2[nj_idx, x1_idx, drop = FALSE])
        } else {
          exo2_clz <- matrix(1, nrow(y2_clz), 1L)
        }
        xx_clz <- crossprod(exo2_clz)
        clz_y2_b <- try(solve(xx_clz, crossprod(exo2_clz, y2_clz)),
          silent = TRUE
        )
        if (inherits(clz_y2_b, "try-error")) {
          clz_y2_b <-
            MASS::ginv(xx_clz) %*% crossprod(exo2_clz, y2_clz)
        }
        clz_y2_hat <- exo2_clz %*% clz_y2_b
        clz_y2_res <- y2_clz - clz_y2_hat
        sample_clz_y2_b[[clz]] <- clz_y2_b
        sample_clz_y2_res[[clz]] <- crossprod(clz_y2_res)
        sample_clz_y2_xx[[clz]] <- crossprod(exo2_clz)

        # Z
        if (length(between_y_idx) > 0L) {
          z_clz_z <- y2[nj_idx, between_y_idx, drop = FALSE]
          if (length(between_x_idx) > 0L) {
            exo_clz_z <- cbind(
              1,
              y2[nj_idx, between_x_idx, drop = FALSE]
            )
          } else {
            exo_clz_z <- matrix(1, nrow(z_clz_z), 1L)
          }
          zz_clz <- crossprod(exo_clz_z)
          clz_zz_b <- try(
            solve(
              zz_clz,
              crossprod(exo_clz_z, z_clz_z)
            ),
            silent = TRUE
          )
          if (inherits(clz_zz_b, "try-error")) {
            clz_zz_b <-
              MASS::ginv(zz_clz) %*% crossprod(exo_clz_z, z_clz_z)
          }
          clz_z_hat <- exo_clz_z %*% clz_zz_b
          clz_z_res <- z_clz_z - clz_z_hat
          sample_clz_zz_b[[clz]] <- clz_zz_b
          sample_clz_zz_res[[clz]] <- crossprod(clz_z_res)
          sample_clz_zz_xx[[clz]] <- crossprod(exo_clz_z)

          sample_clz_yz_res[[clz]] <- crossprod(clz_y2_res, clz_z_res)
          sample_clz_yz_xx[[clz]] <- crossprod(exo2_clz, exo_clz_z)
          sample_clz_yres_xz[[clz]] <- crossprod(clz_y2_res, exo_clz_z)
          sample_clz_xwzres[[clz]] <- crossprod(exo2_clz, clz_z_res)
        }
      } # clz
    } # conditional.x

    ylp[[l]] <- list(
      Y1Y1 = y1y1,
      Y2 = y2, s = s, S.b = s_b, S.PW.start = s_pw_start,
      Sigma.W = s_w, Mu.W = mu_w,
      Sigma.B = sigma_b, Mu.B = mu_b,
      Mu.B.start = mu_b_start, loglik.x = loglik_x,
      mean.d = mean_d, cov.d = cov_d
    )

    # if conditional, add more stuff
    if (conditional_x) {
      if (length(between_y_idx) > 0L) {
        extra <- list(
          sample.wb = sample_wb,
          sample.YYres.wb1 = sample_yyres_wb1,
          sample.XX.wb1 = sample_xx_wb1,
          sample.wb2 = sample_wb2,
          sample.YYres.wb2 = sample_yyres_wb2,
          sample.YresX.wb2 = sample_yres_x_wb2,
          sample.XX.wb2 = sample_xx_wb2,
          sample.clz.Y2.res = sample_clz_y2_res,
          sample.clz.Y2.XX = sample_clz_y2_xx,
          sample.clz.Y2.B = sample_clz_y2_b,
          sample.clz.ZZ.res = sample_clz_zz_res,
          sample.clz.ZZ.XX = sample_clz_zz_xx,
          sample.clz.ZZ.B = sample_clz_zz_b,
          sample.clz.YZ.res = sample_clz_yz_res,
          sample.clz.YZ.XX = sample_clz_yz_xx,
          sample.clz.YresXZ = sample_clz_yres_xz, # zero?
          sample.clz.XWZres = sample_clz_xwzres
        )
      } else {
        extra <- list(
          sample.wb = sample_wb,
          sample.YYres.wb1 = sample_yyres_wb1,
          sample.XX.wb1 = sample_xx_wb1,
          sample.wb2 = sample_wb2,
          sample.YYres.wb2 = sample_yyres_wb2,
          sample.YresX.wb2 = sample_yres_x_wb2,
          sample.XX.wb2 = sample_xx_wb2,
          sample.clz.Y2.res = sample_clz_y2_res,
          sample.clz.Y2.XX = sample_clz_y2_xx,
          sample.clz.Y2.B = sample_clz_y2_b
        )
      }
      ylp[[l]] <- c(ylp[[l]], extra)
    }
  } # l

  ylp
}
