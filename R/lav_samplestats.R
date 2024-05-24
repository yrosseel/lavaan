# constructor for the 'lavSampleStats' class
#
# initial version: YR 25/03/2009
# major revision: YR 5/11/2011: separate data.obs and sample statistics
# YR 5/01/2016: add rescov, resvar, ... if conditional.x = TRUE

# YR 18 Jan 2021: use lavoptions

lav_samplestats_from_data <- function(lavdata = NULL,
                                      lavoptions = NULL,
                                      WLS.V = NULL,
                                      NACOV = NULL) {
  # extra info from lavoptions
  stopifnot(!is.null(lavoptions))
  missing <- lavoptions$missing
  rescale <- lavoptions$sample.cov.rescale
  estimator <- lavoptions$estimator
  mimic <- lavoptions$mimic
  meanstructure <- lavoptions$meanstructure
  correlation <- lavoptions$correlation
  conditional.x <- lavoptions$conditional.x
  fixed.x <- lavoptions$fixed.x
  group.w.free <- lavoptions$group.w.free
  se <- lavoptions$se
  test <- lavoptions$test
  ridge <- lavoptions$ridge
  zero.add <- lavoptions$zero.add
  zero.keep.margins <- lavoptions$zero.keep.margins
  zero.cell.warn <- lavoptions$zero.cell.warn
  dls.a <- lavoptions$estimator.args$dls.a
  dls.GammaNT <- lavoptions$estimator.args$dls.GammaNT
  debug <- lavoptions$debug
  verbose <- lavoptions$verbose

  # sample.icov (new in 0.6-9; ensure it exists, for older objects)
  sample.icov <- TRUE
  if (!is.null(lavoptions$sample.icov)) {
    sample.icov <- lavoptions$sample.icov
  }

  # ridge default
  if (ridge) {
    if (is.numeric(lavoptions$ridge.constant)) {
      ridge.eps <- lavoptions$ridge.constant
    } else {
      ridge.eps <- 1e-5
    }
  } else {
    ridge.eps <- 0.0
  }

  # check lavdata
  stopifnot(!is.null(lavdata))

  # lavdata slots (FIXME: keep lavdata@ names)
  X <- lavdata@X
  Mp <- lavdata@Mp
  ngroups <- lavdata@ngroups
  nlevels <- lavdata@nlevels
  nobs <- lavdata@nobs
  ov.names <- lavdata@ov.names
  ov.names.x <- lavdata@ov.names.x
  DataOv <- lavdata@ov
  eXo <- lavdata@eXo
  WT <- lavdata@weights

  # new in 0.6-6
  # if sampling weights have been used, redefine nobs:
  # per group, we define nobs == sum(wt)
  for (g in seq_len(ngroups)) {
    if (!is.null(WT[[g]])) {
      nobs[[g]] <- sum(WT[[g]])
    }
  }

  # sample.cov.robust cannot be used if sampling weights are used
  if (lavoptions$sample.cov.robust) {
    if (!is.null(WT[[1]])) {
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
  th.idx <- vector("list", length = ngroups)
  th.names <- vector("list", length = ngroups)

  # residual (y | x)
  res.cov <- vector("list", length = ngroups)
  res.var <- vector("list", length = ngroups)
  res.th <- vector("list", length = ngroups)
  res.th.nox <- vector("list", length = ngroups)
  res.slopes <- vector("list", length = ngroups)
  res.int <- vector("list", length = ngroups)

  # fixed.x
  mean.x <- vector("list", length = ngroups)
  cov.x <- vector("list", length = ngroups)

  # binary/ordinal
  bifreq <- vector("list", length = ngroups)

  # extra sample statistics per group
  icov <- vector("list", length = ngroups)
  cov.log.det <- vector("list", length = ngroups)
  res.icov <- vector("list", length = ngroups)
  res.cov.log.det <- vector("list", length = ngroups)
  WLS.obs <- vector("list", length = ngroups)
  missing. <- vector("list", length = ngroups)
  missing.h1. <- vector("list", length = ngroups)
  missing.flag. <- FALSE
  zero.cell.tables <- vector("list", length = ngroups)
  YLp <- vector("list", length = ngroups)

  # group weights
  group.w <- vector("list", length = ngroups)

  # convenience? # FIXME!
  x.idx <- vector("list", length = ngroups)


  WLS.VD <- vector("list", length = ngroups)
  if (is.null(WLS.V)) {
    WLS.V <- vector("list", length = ngroups)
    WLS.V.user <- FALSE
  } else {
    if (!is.list(WLS.V)) {
      if (ngroups == 1L) {
        WLS.V <- list(WLS.V)
      } else {
        lav_msg_stop(gettextf(
          "WLS.V argument should be a list of length %s", ngroups)
        )
      }
    } else {
      if (length(WLS.V) != ngroups) {
        lav_msg_stop(gettextf(
          "WLS.V assumes %1$s groups; data contains %2$s groups",
          length(WLS.V), ngroups))
      }
    }

    # is WLS.V full? check first
    if (is.null(dim(WLS.V[[1]]))) {
      # we will assume it is the diagonal only
      WLS.VD <- WLS.V
      WLS.V <- lapply(WLS.VD, diag)
    } else {
      # create WLS.VD
      WLS.VD <- lapply(WLS.V, diag)
    }

    WLS.V.user <- TRUE
    # FIXME: check dimension of WLS.V!!
  }

  NACOV.compute <- FALSE # since 0.6-6
  if (is.null(NACOV)) {
    NACOV <- vector("list", length = ngroups)
    NACOV.user <- FALSE
    if (se == "robust.sem" && missing == "listwise") {
      NACOV.compute <- TRUE
    }
    # note: test can be a vector...
    if (missing == "listwise" && any(test %in% c(
      "satorra.bentler",
      "mean.var.adjusted",
      "scaled.shifted"
    ))) {
      NACOV.compute <- TRUE
    }
  } else if (is.logical(NACOV)) {
    if (!NACOV) {
      NACOV.compute <- FALSE
    } else {
      NACOV.compute <- TRUE
    }
    NACOV.user <- FALSE
    NACOV <- vector("list", length = ngroups)
  } else {
    if (!is.list(NACOV)) {
      if (ngroups == 1L) {
        NACOV <- list(NACOV)
      } else {
        lav_msg_stop(gettextf(
          "NACOV argument should be a list of length ", ngroups))
      }
    } else {
      if (length(NACOV) != ngroups) {
        lav_msg_stop(gettextf(
          "NACOV assumes %1$s groups; data contains %2$s groups",
          length(NACOV), ngroups))
      }
    }
    NACOV.user <- TRUE
    # FIXME: check dimension of NACOV!!
  }



  # compute some sample statistics per group
  for (g in 1:ngroups) {
    # switch off computing all sample statistics? (housekeeping only)
    if (!is.null(lavoptions$samplestats) && !lavoptions$samplestats) {
	  next
	}
    # check nobs
    if (is.null(WT[[g]])) {
      if (nobs[[g]] < 2L) {
        if (nobs[[g]] == 0L) {
          lav_msg_stop(gettext("data contains no observations"),
            if (ngroups > 1L) gettextf("in group %s", g) else "")
        } else {
          lav_msg_stop(gettext("data contains only a single observation"),
           if (ngroups > 1L) gettextf("in group %s", g) else "")
        }
      }
    }

    # exogenous x?
    nexo <- length(ov.names.x[[g]])
    if (nexo) {
      stopifnot(nexo == NCOL(eXo[[g]]))

      # two cases: ov.names contains 'x' variables, or not
      if (conditional.x) {
        # ov.names.x are NOT in ov.names
        x.idx[[g]] <- length(ov.names[[g]]) + seq_len(nexo)
      } else {
        if (fixed.x) {
          # ov.names.x are a subset of ov.names
          x.idx[[g]] <- match(ov.names.x[[g]], ov.names[[g]])
          stopifnot(!anyNA(x.idx[[g]]))
        } else {
          x.idx[[g]] <- integer(0L)
        }
      }
    } else {
      x.idx[[g]] <- integer(0L)
      conditional.x <- FALSE
      fixed.x <- FALSE
    }

    # group weight
    group.w[[g]] <- nobs[[g]] / sum(unlist(nobs))

    # check if we have categorical data in this group
    categorical <- FALSE
    ov.types <- DataOv$type[match(ov.names[[g]], DataOv$name)]
    ov.levels <- DataOv$nlev[match(ov.names[[g]], DataOv$name)]
    CAT <- list()
    if ("ordered" %in% ov.types) {
      categorical <- TRUE
      if (nlevels > 1L) {
        lav_msg_warn(gettext("multilevel + categorical not supported yet."))
      }
    }

    if (categorical) {
      # compute CAT

      if (estimator %in% c("ML", "REML", "PML", "FML", "MML", "none", "ULS")) {
        WLS.W <- FALSE
        if (estimator == "ULS" && se == "robust.sem") { #||
          # any(test %in% c("satorra.bentler", "scaled.shifted",
          #            "mean.var.adjusted")))) {
          WLS.W <- TRUE
        }
      } else {
        WLS.W <- TRUE
      }
	  # check cat.wls.w option (new in 0.6-18)
	  if (!is.null(lavoptions$cat.wls.w) && !lavoptions$cat.wls.w) {
	    WLS.W <- FALSE # perhaps do.fit = FALSE? (eg sam())
	  }
      if (verbose) {
        cat("Estimating sample thresholds and correlations ... ")
      }

      if (conditional.x) {
        CAT <- muthen1984(
          Data = X[[g]],
          wt = WT[[g]],
          ov.names = ov.names[[g]],
          ov.types = ov.types,
          ov.levels = ov.levels,
          ov.names.x = ov.names.x[[g]],
          eXo = eXo[[g]],
          group = g, # for error messages only
          WLS.W = WLS.W,
          zero.add = zero.add,
          zero.keep.margins = zero.keep.margins,
          zero.cell.warn = FALSE,
          zero.cell.tables = TRUE,
          verbose = debug
        )
      } else {
        CAT <- muthen1984(
          Data = X[[g]],
          wt = WT[[g]],
          ov.names = ov.names[[g]],
          ov.types = ov.types,
          ov.levels = ov.levels,
          ov.names.x = NULL,
          eXo = NULL,
          group = g, # for error messages only
          WLS.W = WLS.W,
          zero.add = zero.add,
          zero.keep.margins = zero.keep.margins,
          zero.cell.warn = FALSE,
          zero.cell.tables = TRUE,
          verbose = debug
        )
      }
      # empty cell tables
      zero.cell.tables[[g]] <- CAT$zero.cell.tables
      if (verbose) cat("done\n")
    }

    if (categorical) {
      # convenience
      th.idx[[g]] <- unlist(CAT$TH.IDX)
      th.names[[g]] <- unlist(CAT$TH.NAMES)

      if (conditional.x) {
        # residual var/cov
        res.var[[g]] <- unlist(CAT$VAR)
        res.cov[[g]] <- unname(CAT$COV)
        if (ridge) {
          diag(res.cov[[g]]) <- diag(res.cov[[g]]) + ridge.eps
          res.var[[g]] <- diag(res.cov[[g]])
        }

        # th also contains the means of numeric variables
        res.th[[g]] <- unlist(CAT$TH)
        res.th.nox[[g]] <- unlist(CAT$TH.NOX)

        # for convenience, we store the intercept of numeric
        # variables in res.int
        NVAR <- NCOL(res.cov[[g]])
        mean[[g]] <- res.int[[g]] <- numeric(NVAR)
        num.idx <- which(!seq_len(NVAR) %in% th.idx[[g]])
        if (length(num.idx) > 0L) {
          NUM.idx <- which(th.idx[[g]] == 0L)
          mean[[g]][num.idx] <- res.th.nox[[g]][NUM.idx]
          res.int[[g]][num.idx] <- res.th[[g]][NUM.idx]
        }

        # slopes
        res.slopes[[g]] <- CAT$SLOPES
      } else {
        # var/cov
        var[[g]] <- unlist(CAT$VAR)
        cov[[g]] <- unname(CAT$COV)
        if (ridge) {
          diag(cov[[g]]) <- diag(cov[[g]]) + ridge.eps
          var[[g]] <- diag(cov[[g]])
        }

        # th also contains the means of numeric variables
        th[[g]] <- unlist(CAT$TH)

        # mean (numeric only)
        NVAR <- NCOL(cov[[g]])
        mean[[g]] <- numeric(NVAR)
        num.idx <- which(!seq_len(NVAR) %in% th.idx[[g]])
        if (length(num.idx) > 0L) {
          NUM.idx <- which(th.idx[[g]] == 0L)
          mean[[g]][num.idx] <- th[[g]][NUM.idx]
        }
      }

      # only for catML
      if (estimator == "catML") {
        COV <- cov2cor(lav_matrix_symmetric_force_pd(cov[[g]],
          tol = 1e-04
        ))
        # overwrite
        cov[[g]] <- COV
        out <- lav_samplestats_icov(
          COV = COV,
          x.idx = x.idx[[g]],
          ngroups = ngroups, g = g, warn = TRUE
        )
        icov[[g]] <- out$icov
        cov.log.det[[g]] <- out$cov.log.det

        # the same for res.cov if conditional.x = TRUE
        if (conditional.x) {
          RES.COV <-
            cov2cor(lav_matrix_symmetric_force_pd(res.cov[[g]],
              tol = 1e-04
            ))
          # overwrite
          res.cov[[g]] <- RES.COV
          out <- lav_samplestats_icov(
            COV = RES.COV,
            ridge = 1e-05,
            x.idx = x.idx[[g]],
            ngroups = ngroups, g = g, warn = TRUE
          )
          res.icov[[g]] <- out$icov
          res.cov.log.det[[g]] <- out$cov.log.det
        }
      }
    } # categorical

    # continuous -- multilevel
    else if (nlevels > 1L) {
      # level-based sample statistics
      YLp[[g]] <- lav_samplestats_cluster_patterns(
        Y = X[[g]],
        Lp = lavdata@Lp[[g]],
        conditional.x = lavoptions$conditional.x
      )


      if (conditional.x) {
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
        Y <- X[[g]] # contains eXo
        COV <- unname(stats::cov(Y, use = "pairwise.complete.obs"))
		# if we have missing values (missing by design?), replace them by 0
		COV[is.na(COV)] <- 0
        MEAN <- unname(colMeans(Y, na.rm = TRUE))
        var[[g]] <- diag(COV)
        # rescale cov by (N-1)/N? (only COV!)
        if (rescale) {
          # we 'transform' the sample cov (divided by n-1)
          # to a sample cov divided by 'n'
          COV <- ((nobs[[g]] - 1) / nobs[[g]]) * COV
        }
        cov[[g]] <- COV
        if (ridge) {
          diag(cov[[g]]) <- diag(cov[[g]]) + ridge.eps
          var[[g]] <- diag(cov[[g]])
        }
        mean[[g]] <- MEAN

        A <- COV[-x.idx[[g]], -x.idx[[g]], drop = FALSE]
        B <- COV[-x.idx[[g]], x.idx[[g]], drop = FALSE]
        C <- COV[x.idx[[g]], x.idx[[g]], drop = FALSE]
        # FIXME: make robust against singular C!!!
        res.cov[[g]] <- A - B %*% solve(C) %*% t(B)
        res.var[[g]] <- diag(cov[[g]])

        MY <- MEAN[-x.idx[[g]]]
        MX <- MEAN[x.idx[[g]]]
        C3 <- rbind(
          c(1, MX),
          cbind(MX, C + tcrossprod(MX))
        )
        B3 <- cbind(MY, B + tcrossprod(MY, MX))
        COEF <- unname(solve(C3, t(B3)))

        res.int[[g]] <- COEF[1, ] # intercepts
        res.slopes[[g]] <- t(COEF[-1, , drop = FALSE]) # slopes
      } else {
        # FIXME: needed?
		COV <- unname(stats::cov(X[[g]], use = "pairwise.complete.obs"))
        # if we have missing values (missing by design?), replace them by 0
        COV[is.na(COV)] <- 0
        cov[[g]] <- COV
        mean[[g]] <- unname(colMeans(X[[g]], na.rm = TRUE))
        var[[g]] <- diag(cov[[g]])

        # missing patterns
        if (missing %in% c("ml", "ml.x")) {
          missing.flag. <- TRUE
          missing.[[g]] <-
            lav_samplestats_missing_patterns(
              Y = X[[g]],
              Mp = Mp[[g]],
              wt = WT[[g]],
              Lp = lavdata@Lp[[g]]
            )
        }
      }
    } # multilevel

    # continuous -- single-level
    else {
      if (conditional.x) {
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

        Y <- cbind(X[[g]], eXo[[g]])
		COV <- unname(stats::cov(Y, use = "pairwise.complete.obs"))
        # if we have missing values (missing by design?), replace them by 0
        COV[is.na(COV)] <- 0
        MEAN <- unname(colMeans(Y, na.rm = TRUE))
        # rescale cov by (N-1)/N? (only COV!)
        if (rescale) {
          # we 'transform' the sample cov (divided by n-1)
          # to a sample cov divided by 'n'
          COV <- ((nobs[[g]] - 1) / nobs[[g]]) * COV
        }
        cov[[g]] <- COV
        var[[g]] <- diag(COV)
        if (ridge) {
          diag(cov[[g]]) <- diag(cov[[g]]) + ridge.eps
          var[[g]] <- diag(cov[[g]])
        }
        mean[[g]] <- MEAN

        A <- COV[-x.idx[[g]], -x.idx[[g]], drop = FALSE]
        B <- COV[-x.idx[[g]], x.idx[[g]], drop = FALSE]
        C <- COV[x.idx[[g]], x.idx[[g]], drop = FALSE]
        # FIXME: make robust against singular C!!!
        res.cov[[g]] <- A - B %*% solve(C) %*% t(B)
        res.var[[g]] <- diag(cov[[g]])


        MY <- MEAN[-x.idx[[g]]]
        MX <- MEAN[x.idx[[g]]]
        C3 <- rbind(
          c(1, MX),
          cbind(MX, C + tcrossprod(MX))
        )
        B3 <- cbind(MY, B + tcrossprod(MY, MX))
        COEF <- unname(solve(C3, t(B3)))

        res.int[[g]] <- COEF[1, ] # intercepts
        res.slopes[[g]] <- t(COEF[-1, , drop = FALSE]) # slopes
      } else if (missing == "two.stage" ||
        missing == "robust.two.stage") {
        missing.flag. <- FALSE # !!! just use sample statistics
        missing.[[g]] <-
          lav_samplestats_missing_patterns(
            Y = X[[g]],
            Mp = Mp[[g]],
            wt = WT[[g]]
          )
        out <- lav_mvnorm_missing_h1_estimate_moments(
          Y = X[[g]],
          wt = WT[[g]],
          Mp = Mp[[g]], Yp = missing.[[g]], verbose = verbose,
          max.iter = lavoptions$em.h1.iter.max,
          tol = lavoptions$em.h1.tol,
          warn = lavoptions$em.h1.warn
        )
        missing.h1.[[g]]$sigma <- out$Sigma
        missing.h1.[[g]]$mu <- out$Mu
        missing.h1.[[g]]$h1 <- out$fx

        # here, sample statistics == EM estimates
        cov[[g]] <- missing.h1.[[g]]$sigma
        if (ridge) {
          diag(cov[[g]]) <- diag(cov[[g]]) + ridge.eps
        }
        var[[g]] <- diag(cov[[g]])
        mean[[g]] <- missing.h1.[[g]]$mu
      } else if (missing %in% c("ml", "ml.x")) {
        missing.flag. <- TRUE
        missing.[[g]] <-
          lav_samplestats_missing_patterns(
            Y = X[[g]],
            Mp = Mp[[g]],
            wt = WT[[g]]
          )

        if (nlevels == 1L) {
          # estimate moments unrestricted model
          out <- lav_mvnorm_missing_h1_estimate_moments(
            Y = X[[g]],
            wt = WT[[g]], verbose = verbose,
            Mp = Mp[[g]], Yp = missing.[[g]],
            max.iter = lavoptions$em.h1.iter.max,
            tol = lavoptions$em.h1.tol,
            warn = lavoptions$em.h1.warn
          )
          missing.h1.[[g]]$sigma <- out$Sigma
          missing.h1.[[g]]$mu <- out$Mu
          missing.h1.[[g]]$h1 <- out$fx
        }

        if (!is.null(WT[[g]])) {
          # here, sample statistics == EM estimates
          cov[[g]] <- missing.h1.[[g]]$sigma
          if (ridge) {
            diag(cov[[g]]) <- diag(cov[[g]]) + ridge.eps
          }
          var[[g]] <- diag(cov[[g]])
          mean[[g]] <- missing.h1.[[g]]$mu
        } else {
          # NEEDED? why not just EM-based?
          COV <- unname(stats::cov(X[[g]], use = "pairwise.complete.obs"))
          # if we have missing values (missing by design?), replace them by 0
		  COV[is.na(COV)] <- 0
          cov[[g]] <- COV
          # rescale cov by (N-1)/N? (only COV!)
          if (rescale) {
            # we 'transform' the sample cov (divided by n-1)
            # to a sample cov divided by 'n'
            cov[[g]] <- ((nobs[[g]] - 1) / nobs[[g]]) * cov[[g]]
          }
          if (ridge) {
            diag(cov[[g]]) <- diag(cov[[g]]) + ridge.eps
          }
          var[[g]] <- diag(cov[[g]])
          mean[[g]] <- colMeans(X[[g]], na.rm = TRUE)
        }
      } else {
        # LISTWISE
        if (!is.null(WT[[g]])) {
          out <- stats::cov.wt(X[[g]],
            wt = WT[[g]],
            method = "ML"
          )
		  COV <- out$cov
          # if we have missing values (missing by design?), replace them by 0
          COV[is.na(COV)] <- 0
          cov[[g]] <- COV
          if (ridge) {
            diag(cov[[g]]) <- diag(cov[[g]]) + ridge.eps
          }
          var[[g]] <- diag(cov[[g]])
          mean[[g]] <- out$center
        } else if (lavoptions$sample.cov.robust) {
          # fixme: allow prob/max.it to be options
          out <- lav_cov_huber(
            Y = X[[g]], prob = 0.95,
            max.it = 200L, tol = 1e-07
          )
          cov[[g]] <- out$Sigma
          var[[g]] <- diag(cov[[g]])
          mean[[g]] <- out$Mu
        } else {
		  COV <- unname(stats::cov(X[[g]], use = "pairwise.complete.obs"))
          # if we have missing values (missing by design?), replace them by 0
          COV[is.na(COV)] <- 0
          cov[[g]] <- COV
          # rescale cov by (N-1)/N? (only COV!)
          if (rescale) {
            # we 'transform' the sample cov (divided by n-1)
            # to a sample cov divided by 'n'
            cov[[g]] <- ((nobs[[g]] - 1) / nobs[[g]]) * cov[[g]]
          }
          if (ridge) {
            diag(cov[[g]]) <- diag(cov[[g]]) + ridge.eps
          }
          var[[g]] <- diag(cov[[g]])
          mean[[g]] <- colMeans(X[[g]], na.rm = TRUE)
        }
      }

      # correlation structure?
      if (correlation) {
        cov[[g]] <- cov2cor(cov[[g]])
        var[[g]] <- rep(1, length(var[[g]]))
        if (conditional.x) {
          res.cov[[g]] <- cov2cor(res.cov[[g]])
          res.var[[g]] <- rep(1, length(res.var[[g]]))
          cov.x[[g]] <- cov2cor(cov.x[[g]])
          # FIXME: slopes? more?
        }
      }

      # icov and cov.log.det (but not if missing)
      if (sample.icov && !missing %in% c("ml", "ml.x")) {
        out <- lav_samplestats_icov(
          COV = cov[[g]], ridge = 1e-05,
          x.idx = x.idx[[g]],
          ngroups = ngroups, g = g, warn = TRUE
        )
        icov[[g]] <- out$icov
        cov.log.det[[g]] <- out$cov.log.det

        # the same for res.cov if conditional.x = TRUE
        if (conditional.x) {
          out <- lav_samplestats_icov(
            COV = res.cov[[g]],
            ridge = 1e-05,
            x.idx = x.idx[[g]],
            ngroups = ngroups, g = g, warn = TRUE
          )
          res.icov[[g]] <- out$icov
          res.cov.log.det[[g]] <- out$cov.log.det
        }
      }
    } # continuous - single level


    # WLS.obs
    if (nlevels == 1L) {
      if (estimator == "catML") {
        # correlations only (for now)
        tmp.categorical <- FALSE
        tmp.meanstructure <- FALSE
      } else {
        tmp.categorical <- categorical
        tmp.meanstructure <- meanstructure
      }
      WLS.obs[[g]] <- lav_samplestats_wls_obs(
        mean.g = mean[[g]],
        cov.g = cov[[g]], var.g = var[[g]], th.g = th[[g]],
        th.idx.g = th.idx[[g]], res.int.g = res.int[[g]],
        res.cov.g = res.cov[[g]], res.var.g = res.var[[g]],
        res.th.g = res.th[[g]], res.slopes.g = res.slopes[[g]],
        group.w.g = log(nobs[[g]]),
        categorical = tmp.categorical, conditional.x = conditional.x,
        meanstructure = tmp.meanstructure, correlation = correlation,
        slopestructure = conditional.x,
        group.w.free = group.w.free
      )
    }

    # fill in the other slots
    if (!is.null(eXo[[g]])) {
      if (!is.null(WT[[g]])) {
        if (missing != "listwise") {
          cov.x[[g]] <- missing.h1.[[g]]$sigma[x.idx[[g]],
            x.idx[[g]],
            drop = FALSE
          ]
          mean.x[[g]] <- missing.h1.[[g]]$mu[x.idx[[g]]]
        } else {
          out <- stats::cov.wt(eXo[[g]],
            wt = WT[[g]],
            method = "ML"
          )
          cov.x[[g]] <- out$cov
          mean.x[[g]] <- out$center
        }
      } else {
        cov.x[[g]] <- cov(eXo[[g]], use = "pairwise")
        if (rescale) {
          # we 'transform' the sample cov (divided by n-1)
          # to a sample cov divided by 'n'
          cov.x[[g]] <- ((nobs[[g]] - 1) / nobs[[g]]) * cov.x[[g]]
        }
        mean.x[[g]] <- colMeans(eXo[[g]])
      }
    }

    # NACOV (=GAMMA)
    if (!NACOV.user && nlevels == 1L) {
      if (estimator == "ML" && !missing.flag. && NACOV.compute) {
        if (conditional.x) {
          Y <- Y
        } else {
          Y <- X[[g]]
        }

        if (length(lavdata@cluster) > 0L) {
          cluster.idx <- lavdata@Lp[[g]]$cluster.idx[[2]]
        } else {
          cluster.idx <- NULL
        }

        NACOV[[g]] <-
          lav_samplestats_Gamma(
            Y = Y,
            x.idx = x.idx[[g]],
            cluster.idx = cluster.idx,
            fixed.x = fixed.x,
            conditional.x = conditional.x,
            meanstructure = meanstructure,
            slopestructure = conditional.x,
            gamma.n.minus.one =
              lavoptions$gamma.n.minus.one,
            unbiased = lavoptions$gamma.unbiased,
            Mplus.WLS = FALSE
          )
      } else if (estimator %in% c("WLS", "DWLS", "ULS", "DLS", "catML")) {
        if (!categorical) {
          # sample size large enough?
          nvar <- ncol(X[[g]])
          # if(conditional.x && nexo > 0L) {
          #    nvar <- nvar - nexo
          # }
          pstar <- nvar * (nvar + 1) / 2
          if (meanstructure) pstar <- pstar + nvar
          if (conditional.x && nexo > 0L) {
            pstar <- pstar + (nvar * nexo)
          }
          if (nrow(X[[g]]) < pstar) {
            lav_msg_warn(gettextf(
              "number of observations (%s) too small to compute Gamma",
              nrow(X[[g]])),
              if (ngroups > 1L) gettextf("in group %s", g) else ""
            )
          }
          if (conditional.x) {
            Y <- Y
          } else {
            Y <- X[[g]]
          }

          if (length(lavdata@cluster) > 0L) {
            cluster.idx <- lavdata@Lp[[g]]$cluster.idx[[2]]
          } else {
            cluster.idx <- NULL
          }
          NACOV[[g]] <-
            lav_samplestats_Gamma(
              Y = Y,
              x.idx = x.idx[[g]],
              cluster.idx = cluster.idx,
              fixed.x = fixed.x,
              conditional.x = conditional.x,
              meanstructure = meanstructure,
              slopestructure = conditional.x,
              gamma.n.minus.one =
                lavoptions$gamma.n.minus.one,
              unbiased =
                lavoptions$gamma.unbiased,
              Mplus.WLS = (mimic == "Mplus")
            )
        } else { # categorical case
          NACOV[[g]] <- CAT$WLS.W * nobs[[g]]
          if (lavoptions$gamma.n.minus.one) {
            NACOV[[g]] <- NACOV[[g]] * (nobs[[g]] / (nobs[[g]] - 1L))
          }
          if (estimator == "catML") {
            # remove all but the correlation part
            ntotal <- nrow(NACOV[[g]])
            pstar <- nrow(CAT$A22)
            nocor <- ntotal - pstar
            if (length(nocor) > 0L) {
              NACOV[[g]] <- NACOV[[g]][
                -seq_len(nocor),
                -seq_len(nocor)
              ]
            }
          }
        }
      } else if (estimator == "PML") {
        # no NACOV ... for now
      }

      # group.w.free
      if (!is.null(NACOV[[g]]) && group.w.free) {
        # unweight!!
        a <- group.w[[g]] * sum(unlist(nobs)) / nobs[[g]]
        # always 1!!!
        NACOV[[g]] <- lav_matrix_bdiag(matrix(a, 1, 1), NACOV[[g]])
      }
    }

    # WLS.V
    if (!WLS.V.user && nlevels == 1L) {
      if (estimator == "DLS" && dls.GammaNT == "sample" && dls.a < 1.0) {
        # compute GammaNT here
        GammaNT <- lav_samplestats_Gamma_NT(
          COV            = cov[[g]],
          MEAN           = mean[[g]],
          rescale        = FALSE,
          x.idx          = x.idx[[g]],
          fixed.x        = fixed.x,
          conditional.x  = conditional.x,
          meanstructure  = meanstructure,
          slopestructure = conditional.x
        )
      }

      if (estimator == "GLS" ||
        (estimator == "DLS" && dls.GammaNT == "sample" &&
          dls.a == 1.0)) {
        # Note: we need the 'original' COV/MEAN/ICOV
        #        sample statistics; not the 'residual' version
        WLS.V[[g]] <- lav_samplestats_Gamma_inverse_NT(
          ICOV           = icov[[g]],
          COV            = cov[[g]],
          MEAN           = mean[[g]],
          rescale        = FALSE,
          x.idx          = x.idx[[g]],
          fixed.x        = fixed.x,
          conditional.x  = conditional.x,
          meanstructure  = meanstructure,
          slopestructure = conditional.x
        )
        if (mimic == "Mplus" && !conditional.x && meanstructure) {
          # bug in Mplus? V11 rescaled by nobs[[g]]/(nobs[[g]]-1)
          nvar <- NCOL(cov[[g]])
          WLS.V[[g]][1:nvar, 1:nvar] <-
            WLS.V[[g]][1:nvar, 1:nvar,
              drop = FALSE
            ] * (nobs[[g]] / (nobs[[g]] - 1))
        }
      } else if (estimator == "ML") {
        # no WLS.V here, since function of model-implied moments
      } else if (estimator %in% c("WLS", "DWLS", "ULS", "DLS")) {
        if (!categorical) {
          if (estimator == "WLS" || estimator == "DLS") {
            if (!fixed.x) {
              if (estimator != "DLS") {
                # Gamma should be po before we invert
                ev <- eigen(NACOV[[g]], # symmetric=FALSE,
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
              if (estimator == "DLS" && dls.GammaNT == "sample") {
                if (dls.a == 1.0) {
                  # nothing to do, use GLS version
                } else {
                  W.DLS <-
                    (1 - dls.a) * NACOV[[g]] + dls.a * GammaNT
                  WLS.V[[g]] <-
                    lav_matrix_symmetric_inverse(W.DLS)
                }
              } else { # WLS
                WLS.V[[g]] <-
                  lav_matrix_symmetric_inverse(NACOV[[g]])
              }
            } else {
              # fixed.x: we have zero cols/rows
              # ginv does the trick, but perhaps this is overkill
              # just removing the zero rows/cols, invert, and
              # fill back in the zero rows/cols would do it
              # WLS.V[[g]] <- MASS::ginv(NACOV[[g]])
              if (estimator == "DLS" && dls.GammaNT == "sample") {
                W.DLS <- (1 - dls.a) * NACOV[[g]] + dls.a * GammaNT
                WLS.V[[g]] <-
                  lav_matrix_symmetric_inverse(W.DLS)
              } else { # WLS
                WLS.V[[g]] <-
                  lav_matrix_symmetric_inverse(NACOV[[g]])
              }
            }
          } else if (estimator == "DWLS") {
            dacov <- diag(NACOV[[g]])
            if (!all(is.finite(dacov))) {
              lav_msg_stop(gettext(
                "diagonal of Gamma (NACOV) contains non finite values"))
            }
            if (fixed.x) {
              # structural zeroes!
              zero.idx <- which(dacov == 0.0)
              idacov <- 1 / dacov
              idacov[zero.idx] <- 0.0
            } else {
              idacov <- 1 / dacov
            }
            WLS.V[[g]] <- diag(idacov,
              nrow = NROW(NACOV[[g]]),
              ncol = NCOL(NACOV[[g]])
            )
            WLS.VD[[g]] <- idacov
          } else if (estimator == "ULS") {
            # WLS.V[[g]] <- diag(length(WLS.obs[[g]]))
            WLS.VD[[g]] <- rep(1, length(WLS.obs[[g]]))
          }
        } else {
          if (estimator == "WLS") {
            WLS.V[[g]] <- inv.chol(CAT$WLS.W * nobs[[g]])
          } else if (estimator == "DWLS") {
            dacov <- diag(CAT$WLS.W * nobs[[g]])
            # WLS.V[[g]] <- diag(1/dacov, nrow=NROW(CAT$WLS.W),
            #                            ncol=NCOL(CAT$WLS.W))
            WLS.VD[[g]] <- 1 / dacov
          } else if (estimator == "ULS") {
            # WLS.V[[g]] <- diag(length(WLS.obs[[g]]))
            WLS.VD[[g]] <- rep(1, length(WLS.obs[[g]]))
          }
        }
      } else if (estimator == "PML" || estimator == "FML") {
        # no WLS.V here
      }

      # group.w.free (only if categorical)
      if (group.w.free && categorical) {
        if (!is.null(WLS.V[[g]])) {
          # unweight!!
          a <- group.w[[g]] * sum(unlist(nobs)) / nobs[[g]]
          # always 1!!!
          # invert
          a <- 1 / a
          WLS.V[[g]] <- lav_matrix_bdiag(matrix(a, 1, 1), WLS.V[[g]])
        }
        if (!is.null(WLS.VD[[g]])) {
          # unweight!!
          a <- group.w[[g]] * sum(unlist(nobs)) / nobs[[g]]
          # always 1!!!
          # invert
          a <- 1 / a
          WLS.VD[[g]] <- c(a, WLS.VD[[g]])
        }
      }
    }
  } # ngroups

  # remove 'CAT', unless debug -- this is to save memory
  if (!debug) {
    CAT <- list()
  }

  # construct SampleStats object
  lavSampleStats <- new("lavSampleStats",
    # sample moments
    th = th,
    th.idx = th.idx,
    th.names = th.names,
    mean = mean,
    cov = cov,
    var = var,

    # residual (y | x)
    res.cov = res.cov,
    res.var = res.var,
    res.th = res.th,
    res.th.nox = res.th.nox,
    res.slopes = res.slopes,
    res.int = res.int,
    mean.x = mean.x,
    cov.x = cov.x,
    bifreq = bifreq,
    group.w = group.w,

    # convenience
    nobs = nobs,
    ntotal = sum(unlist(nobs)),
    ngroups = ngroups,
    x.idx = x.idx,

    # extra sample statistics
    icov = icov,
    cov.log.det = cov.log.det,
    res.icov = res.icov,
    res.cov.log.det = res.cov.log.det,
    ridge = ridge.eps,
    WLS.obs = WLS.obs,
    WLS.V = WLS.V,
    WLS.VD = WLS.VD,
    NACOV = NACOV,
    NACOV.user = NACOV.user,

    # cluster/levels
    YLp = YLp,

    # missingness
    missing.flag = missing.flag.,
    missing = missing.,
    missing.h1 = missing.h1.,
    zero.cell.tables = zero.cell.tables
  )

  # just a SINGLE warning if we have empty cells
  if ((!is.null(lavoptions$samplestats) && lavoptions$samplestats) &&
      categorical && zero.cell.warn &&
      any(sapply(zero.cell.tables, nrow) > 0L)) {
    nempty <- sum(sapply(zero.cell.tables, nrow))
    lav_msg_warn(gettextf(
      "%s bivariate tables have empty cells; to see them, use:
      lavInspect(fit, \"zero.cell.tables\")", nempty)
    )
  }

  lavSampleStats
}


lav_samplestats_from_moments <- function(sample.cov = NULL,
                                         sample.mean = NULL,
                                         sample.th = NULL,
                                         sample.nobs = NULL,
                                         ov.names = NULL, # including x
                                         ov.names.x = NULL,
                                         WLS.V = NULL,
                                         NACOV = NULL,
                                         lavoptions = NULL) {
  # extract options
  estimator <- lavoptions$estimator
  mimic <- lavoptions$mimic
  meanstructure <- lavoptions$meanstructure
  group.w.free <- lavoptions$group.w.free
  ridge <- lavoptions$ridge
  rescale <- lavoptions$sample.cov.rescale

  # no multilevel yet
  nlevels <- 1L

  # ridge default
  if (ridge) {
    if (is.numeric(lavoptions$ridge.constant)) {
      ridge.eps <- lavoptions$ridge.constant
    } else {
      ridge.eps <- 1e-5
    }
  } else {
    ridge.eps <- 0.0
  }

  # new in 0.6-3:
  # check if sample.cov has attributes if conditional.x = TRUE
  sample.res.slopes <- attr(sample.cov, "res.slopes")
  sample.cov.x <- attr(sample.cov, "cov.x")
  sample.mean.x <- attr(sample.cov, "mean.x")
  if (!is.null(sample.res.slopes)) {
    conditional.x <- TRUE
    # strip attributes
    attr(sample.cov, "res.slopes") <- NULL
    attr(sample.cov, "cov.x") <- NULL
    attr(sample.cov, "mean.x") <- NULL
    # make list
    if (!is.list(sample.res.slopes)) {
      sample.res.slopes <- list(sample.res.slopes)
    }
    if (!is.list(sample.cov.x)) {
      sample.cov.x <- list(sample.cov.x)
    }
    if (!is.list(sample.mean.x)) {
      sample.mean.x <- list(sample.mean.x)
    }
  } else if (!is.null(sample.cov.x)) {
    conditional.x <- FALSE
    fixed.x <- TRUE

    # strip attributes
    attr(sample.cov, "cov.x") <- NULL
    attr(sample.cov, "mean.x") <- NULL
    # make list
    if (!is.list(sample.cov.x)) {
      sample.cov.x <- list(sample.cov.x)
    }
    if (!is.list(sample.mean.x)) {
      sample.mean.x <- list(sample.mean.x)
    }
  } else if (is.null(sample.cov.x) && length(unlist(ov.names.x)) > 0L) {
    # fixed.x = TRUE, but only joint sample.cov is provided
    conditional.x <- FALSE
    fixed.x <- TRUE

    # create sample.cov.x and sample.mean.x later...
  } else {
    conditional.x <- FALSE
    fixed.x <- FALSE
  }

  # matrix -> list
  if (!is.list(sample.cov)) {
    sample.cov <- list(sample.cov)
  }

  # number of groups
  ngroups <- length(sample.cov)

  # ov.names
  if (!is.list(ov.names)) {
    ov.names <- rep(list(ov.names), ngroups)
  }
  if (!is.list(ov.names.x)) {
    ov.names.x <- rep(list(ov.names.x), ngroups)
  }

  if (!is.null(sample.mean)) {
    meanstructure <- TRUE
    if (!is.list(sample.mean)) {
      # check if sample.mean is string (between single quotes)
      if (is.character(sample.mean)) {
        sample.mean <- char2num(sample.mean)
      }
      sample.mean <- list(unname(sample.mean))
    } else {
      sample.mean <- lapply(lapply(sample.mean, unname), unclass)
    }
  }

  if (!is.null(sample.th)) {
    th.idx <- attr(sample.th, "th.idx")
    attr(sample.th, "th.idx") <- NULL
    if (is.null(th.idx)) {
      lav_msg_stop(gettext("sample.th should have a th.idx attribute"))
    } else {
      if (is.list(th.idx)) {
        th.names <- lapply(th.idx, names)
        th.idx <- lapply(lapply(th.idx, unname), unclass)
      } else {
        th.names <- list(names(th.idx))
        th.idx <- list(unclass(unname(th.idx)))
      }
    }
    if (is.list(sample.th)) {
      # strip names and lavaan.vector class
      sample.th <- lapply(lapply(sample.th, unname), unclass)
    } else {
      # strip names and lavaan.vector class, make list
      sample.th <- list(unclass(unname(sample.th)))
    }
  } else {
    th.idx <- vector("list", length = ngroups)
    th.names <- vector("list", length = ngroups)
  }

  # sample statistics per group
  cov <- vector("list", length = ngroups)
  var <- vector("list", length = ngroups)
  mean <- vector("list", length = ngroups)
  th <- vector("list", length = ngroups)
  # th.idx      <- vector("list", length = ngroups)
  # th.names    <- vector("list", length = ngroups)

  # residual (y | x)
  res.cov <- vector("list", length = ngroups)
  res.var <- vector("list", length = ngroups)
  res.slopes <- vector("list", length = ngroups)
  res.int <- vector("list", length = ngroups)
  res.th <- vector("list", length = ngroups)
  res.th.nox <- vector("list", length = ngroups)

  # fixed.x / conditional.x
  mean.x <- vector("list", length = ngroups)
  cov.x <- vector("list", length = ngroups)

  bifreq <- vector("list", length = ngroups)

  # extra sample statistics per group
  icov <- vector("list", length = ngroups)
  cov.log.det <- vector("list", length = ngroups)
  res.icov <- vector("list", length = ngroups)
  res.cov.log.det <- vector("list", length = ngroups)
  WLS.obs <- vector("list", length = ngroups)
  missing. <- vector("list", length = ngroups)
  missing.h1. <- vector("list", length = ngroups)
  missing.flag. <- FALSE
  zero.cell.tables <- vector("list", length = ngroups)
  YLp <- vector("list", length = ngroups)

  # group weights
  group.w <- vector("list", length = ngroups)
  x.idx <- vector("list", length = ngroups)

  categorical <- FALSE
  if (!is.null(sample.th)) {
    categorical <- TRUE
  }

  WLS.VD <- vector("list", length = ngroups)
  if (is.null(WLS.V)) {
    WLS.V <- vector("list", length = ngroups)
    WLS.V.user <- FALSE
  } else {
    if (!is.list(WLS.V)) {
      if (ngroups == 1L) {
        WLS.V <- list(unclass(WLS.V))
      } else {
        lav_msg_stop(gettextf("WLS.V argument should be a list of length %s",
          ngroups)
        )
      }
    } else {
      if (length(WLS.V) != ngroups) {
        lav_msg_stop(gettextf(
          "WLS.V assumes %1$s groups; data contains %2$s groups",
          length(WLS.V), ngroups))
      }
      WLS.V <- lapply(WLS.V, unclass)
    }

    # is WLS.V full? check first
    if (is.null(dim(WLS.V[[1]]))) {
      # we will assume it is the diagonal only
      WLS.VD <- WLS.V
      WLS.V <- lapply(WLS.VD, diag)
    } else {
      # create WLS.VD
      WLS.VD <- lapply(WLS.V, diag)
      # we could remove WLS.V to save space...
    }

    WLS.V.user <- TRUE
    # FIXME: check dimension of WLS.V!!
  }

  if (is.null(NACOV)) {
    NACOV <- vector("list", length = ngroups)
    NACOV.user <- FALSE
  } else {
    if (!is.list(NACOV)) {
      if (ngroups == 1L) {
        NACOV <- list(unclass(NACOV))
      } else {
        lav_msg_stop(gettextf(
          "NACOV argument should be a list of length %s", ngroups))
      }
    } else {
      if (length(NACOV) != ngroups) {
        lav_msg_stop(gettextf(
          "NACOV assumes %1$s groups; data contains %2$s groups",
          length(NACOV), ngroups))
      }
      NACOV <- lapply(NACOV, unclass)
    }
    NACOV.user <- TRUE
    # FIXME: check dimension of NACOV!!
  }

  nobs <- as.list(as.integer(sample.nobs))


  for (g in 1:ngroups) {
    # exogenous x?
    nexo <- length(ov.names.x[[g]])
    if (nexo) {
      # two cases: ov.names contains 'x' variables, or not
      if (conditional.x) {
        # ov.names.x are NOT in ov.names
        x.idx[[g]] <- which(ov.names[[g]] %in% ov.names.x[[g]])
      } else {
        if (fixed.x) {
          # ov.names.x are a subset of ov.names
          x.idx[[g]] <- match(ov.names.x[[g]], ov.names[[g]])
          stopifnot(!anyNA(x.idx[[g]]))
        } else {
          x.idx[[g]] <- integer(0L)
        }
      }
    } else {
      x.idx[[g]] <- integer(0L)
      conditional.x <- FALSE
      fixed.x <- FALSE
    }


    # group weight
    group.w[[g]] <- nobs[[g]] / sum(unlist(nobs))

    tmp.cov <- sample.cov[[g]]

    # make sure that the matrix is fully symmetric (NEEDED?)
    T <- t(tmp.cov)
    tmp.cov[upper.tri(tmp.cov)] <- T[upper.tri(T)]

    # check dimnames
    if (!is.null(rownames(tmp.cov))) {
      cov.names <- rownames(tmp.cov)
    } else if (!is.null(colnames(tmp.cov))) {
      cov.names <- colnames(tmp.cov)
    } else {
      lav_msg_stop(gettext(
        "please provide row/col names for the covariance matrix!"))
    }

    # extract only the part we need (using ov.names)
    if (conditional.x) {
      idx <- match(ov.names[[g]][-x.idx[[g]]], cov.names)
    } else {
      idx <- match(ov.names[[g]], cov.names)
    }
    if (any(is.na(idx))) {
      cat("found: ", cov.names, "\n")
      cat("expected: ", ov.names[[g]], "\n")
      lav_msg_stop(gettextf(
        "rownames of covariance matrix do not match the model!
        found: %1$s expected: %2$s", lav_msg_view(cov.names),
        lav_msg_view(ov.names[[g]])))
    } else {
      tmp.cov <- tmp.cov[idx, idx, drop = FALSE]
    }

    # strip dimnames
    dimnames(tmp.cov) <- NULL

    if (is.null(sample.mean)) {
      # assume zero mean vector
      tmp.mean <- numeric(ncol(tmp.cov))
    } else {
      # extract only the part we need
      tmp.mean <- unclass(sample.mean[[g]][idx])
    }



    if (categorical) {
      # categorical + conditional.x = TRUE
      if (conditional.x) {
        th.g <- numeric(length(th.idx[[g]]))
        ord.idx <- which(th.idx[[g]] > 0)
        num.idx <- which(th.idx[[g]] == 0)
        if (length(ord.idx) > 0L) {
          th.g[ord.idx] <- sample.th[[g]]
        }
        if (length(num.idx) > 0L) {
          ord.var.idx <- unique(th.idx[[g]][th.idx[[g]] > 0])
          th.g[num.idx] <- -1 * sample.mean[[g]][-ord.var.idx]
        }
        res.th[[g]] <- th.g
        res.th.nox[[g]] <- sample.th[[g]]

        res.cov[[g]] <- tmp.cov
        if (ridge) {
          diag(res.cov[[g]]) <- diag(res.cov[[g]]) + ridge.eps
        }
        res.var[[g]] <- diag(tmp.cov)
        res.int[[g]] <- tmp.mean

        res.slopes[[g]] <- unclass(unname(sample.res.slopes[[g]]))
        cov.x[[g]] <- unclass(unname(sample.cov.x[[g]]))
        mean.x[[g]] <- unclass(unname(sample.mean.x[[g]]))

        # th.idx and th.names are already ok

        # categorical + conditional.x = FALSE
      } else {
        th.g <- numeric(length(th.idx[[g]]))
        ord.idx <- which(th.idx[[g]] > 0)
        num.idx <- which(th.idx[[g]] == 0)
        if (length(ord.idx) > 0L) {
          th.g[ord.idx] <- sample.th[[g]]
        }
        if (length(num.idx) > 0L) {
          ord.var.idx <- unique(th.idx[[g]][th.idx[[g]] > 0])
          th.g[num.idx] <- -1 * sample.mean[[g]][-ord.var.idx]
        }
        th[[g]] <- th.g

        cov[[g]] <- tmp.cov
        if (ridge) {
          diag(cov[[g]]) <- diag(cov[[g]]) + ridge.eps
        }
        var[[g]] <- diag(tmp.cov)
        mean[[g]] <- tmp.mean

        # fixed.x? (needed?)
        if (fixed.x) {
          cov.x[[g]] <- unclass(unname(sample.cov.x[[g]]))
          mean.x[[g]] <- unclass(unname(sample.mean.x[[g]]))
        }

        # th, th.idx and th.names are already ok
      }

      # multilevel
    } else if (nlevels > 1L) {
      lav_msg_stop(gettext("multilevel + sample stats not ready yet"))


      # single level
    } else {
      # single-level + continuous + conditional.x = TRUE
      if (conditional.x) {
        res.cov[[g]] <- tmp.cov
        if (ridge) {
          diag(res.cov[[g]]) <- diag(res.cov[[g]]) + ridge.eps
        }
        res.var[[g]] <- diag(tmp.cov)
        res.int[[g]] <- tmp.mean
        res.slopes[[g]] <- unclass(unname(sample.res.slopes[[g]]))
        cov.x[[g]] <- unclass(unname(sample.cov.x[[g]]))
        mean.x[[g]] <- unclass(unname(sample.mean.x[[g]]))

        # no rescale!

        # icov and cov.log.det
        # if(lavoptions$sample.icov) {
        out <- lav_samplestats_icov(
          COV = res.cov[[g]],
          ridge = 1e-05,
          x.idx = x.idx[[g]],
          ngroups = ngroups, g = g,
          warn = TRUE
        )
        res.icov[[g]] <- out$icov
        res.cov.log.det[[g]] <- out$cov.log.det
        # }

        # continuous + conditional.x = FALSE
      } else {
        cov[[g]] <- tmp.cov
        mean[[g]] <- tmp.mean

        # rescale cov by (N-1)/N?
        if (rescale) {
          # we 'transform' the sample cov (divided by n-1)
          # to a sample cov divided by 'n'
          cov[[g]] <- ((nobs[[g]] - 1) / nobs[[g]]) * cov[[g]]
        }
        if (ridge) {
          diag(cov[[g]]) <- diag(cov[[g]]) + ridge.eps
        }
        var[[g]] <- diag(cov[[g]])

        # icov and cov.log.det
        # if(lavoptions$sample.icov) {
        out <- lav_samplestats_icov(
          COV = cov[[g]],
          ridge = 1e-05,
          x.idx = x.idx[[g]],
          ngroups = ngroups,
          g = g,
          warn = TRUE
        )
        icov[[g]] <- out$icov
        cov.log.det[[g]] <- out$cov.log.det
        # }

        # fixed.x?
        if (fixed.x) {
          if (is.null(sample.cov.x)) {
            cov.x[[g]] <- cov[[g]][x.idx[[g]], x.idx[[g]],
              drop = FALSE
            ]
          } else {
            cov.x[[g]] <- unclass(unname(sample.cov.x[[g]]))
          }
          if (is.null(sample.mean.x)) {
            mean.x[[g]] <- mean[[g]][x.idx[[g]]]
          } else {
            mean.x[[g]] <- unclass(unname(sample.mean.x[[g]]))
          }
        }
      }
    }

    # WLS.obs
    WLS.obs[[g]] <- lav_samplestats_wls_obs(
      mean.g = mean[[g]],
      cov.g = cov[[g]], var.g = var[[g]], th.g = th[[g]],
      th.idx.g = th.idx[[g]], res.int.g = res.int[[g]],
      res.cov.g = res.cov[[g]], res.var.g = res.var[[g]],
      res.th.g = res.th[[g]], res.slopes.g = res.slopes[[g]],
      group.w.g = log(nobs[[g]]),
      categorical = categorical, conditional.x = conditional.x,
      meanstructure = meanstructure, slopestructure = conditional.x,
      group.w.free = group.w.free
    )

    # WLS.V
    if (!WLS.V.user) {
      if (estimator == "GLS") {
        # FIXME: in <0.5-21, we had
        # V11 <- icov[[g]]
        #    if(mimic == "Mplus") { # is this a bug in Mplus?
        #        V11 <- V11 * nobs[[g]]/(nobs[[g]]-1)
        #    }
        WLS.V[[g]] <- lav_samplestats_Gamma_inverse_NT(
          ICOV = icov[[g]],
          COV = cov[[g]],
          MEAN = mean[[g]],
          rescale = FALSE,
          x.idx = x.idx[[g]],
          fixed.x = fixed.x,
          conditional.x = conditional.x,
          meanstructure = meanstructure,
          slopestructure = conditional.x
        )
      } else if (estimator == "ULS") {
        WLS.V[[g]] <- diag(length(WLS.obs[[g]]))
        WLS.VD[[g]] <- rep(1, length(WLS.obs[[g]]))
      } else if (estimator == "WLS" || estimator == "DWLS") {
        if (is.null(WLS.V[[g]])) {
          lav_msg_stop(gettext(
            "the (D)WLS estimator is only available with full data
            or with a user-provided WLS.V"))
        }
      }

      # group.w.free
      if (!is.null(WLS.V[[g]]) && group.w.free) {
        # FIXME!!!
        WLS.V[[g]] <- lav_matrix_bdiag(matrix(1, 1, 1), WLS.V[[g]])
      }
    }
  } # ngroups

  # construct SampleStats object
  lavSampleStats <- new("lavSampleStats",
    # sample moments
    th = th,
    th.idx = th.idx,
    th.names = th.names,
    mean = mean,
    cov = cov,
    var = var,

    # residual (y | x)
    res.cov = res.cov,
    res.var = res.var,
    res.th = res.th,
    res.th.nox = res.th.nox,
    res.slopes = res.slopes,
    res.int = res.int,

    # fixed.x
    mean.x = mean.x,
    cov.x = cov.x,

    # other
    bifreq = bifreq,
    group.w = group.w,

    # convenience
    nobs = nobs,
    ntotal = sum(unlist(nobs)),
    ngroups = ngroups,
    x.idx = x.idx,

    # extra sample statistics
    icov = icov,
    cov.log.det = cov.log.det,
    res.icov = res.icov,
    res.cov.log.det = res.cov.log.det,
    ridge = ridge.eps,
    WLS.obs = WLS.obs,
    WLS.V = WLS.V,
    WLS.VD = WLS.VD,
    NACOV = NACOV,
    NACOV.user = NACOV.user,

    # cluster/level
    YLp = YLp,

    # missingness
    missing.flag = missing.flag.,
    missing = missing.,
    missing.h1 = missing.h1.,
    zero.cell.tables = zero.cell.tables
  )

  lavSampleStats
}

# compute sample statistics, per missing pattern
lav_samplestats_missing_patterns <- function(Y = NULL, Mp = NULL, wt = NULL,
                                             Lp = NULL) {
  # coerce Y to matrix
  Y <- as.matrix(Y)

  # handle two-level data
  if (!is.null(Lp)) {
    Y.orig <- Y
    Z <- NULL
    if (length(Lp$between.idx[[2]]) > 0L) {
      Y <- Y[, -Lp$between.idx[[2]], drop = FALSE]
      z.idx <- which(!duplicated(Lp$cluster.idx[[2]]))
      Z <- Y.orig[z.idx, Lp$between.idx[[2]], drop = FALSE]
    }
  }

  if (is.null(Mp)) {
    Mp <- lav_data_missing_patterns(Y,
      sort.freq = FALSE, coverage = FALSE,
      Lp = Lp
    )
  }

  Yp <- vector("list", length = Mp$npatterns)

  # fill in pattern statistics
  for (p in seq_len(Mp$npatterns)) {
    # extract raw data for these cases
    RAW <- Y[Mp$case.idx[[p]], Mp$pat[p, ], drop = FALSE]

    # more than one case
    if (Mp$freq[p] > 1L) {
      if (!is.null(wt)) {
        out <- stats::cov.wt(RAW,
          wt = wt[Mp$case.idx[[p]]],
          method = "ML"
        )
        SY <- out$cov
        MY <- out$center
      } else {
        MY <- base::.colMeans(RAW, m = NROW(RAW), n = NCOL(RAW))
        # SY <- crossprod(RAW)/Mp$freq[p] - tcrossprod(MY)
        # bad practice, better like this:
        SY <- lav_matrix_cov(RAW)
      }
    }
    # only a single observation (no need to weight!)
    else {
      SY <- 0
      MY <- as.numeric(RAW)
    }

    if (!is.null(wt)) {
      FREQ <- sum(wt[Mp$case.idx[[p]]])
    } else {
      FREQ <- Mp$freq[p]
    }

    # store sample statistics, var.idx and freq
    Yp[[p]] <- list(
      SY = SY, MY = MY, var.idx = Mp$pat[p, ],
      freq = FREQ
    )

    # if clustered data, add rowsum over all cases per cluster
    if (!is.null(Lp)) {
      tmp <- rowsum.default(RAW, group = Mp$j.idx[[p]], reorder = FALSE)
      Yp[[p]]$ROWSUM <- tmp
    }
  }

  # add Zp as an attribute
  # if(!is.null(Lp)) {
  #    Zp <- lav_samplestats_missing_patterns(Y = Z, Mp = Mp$Zp)
  #    for(p in Mp$Zp$npatterns) {
  #        this.z <- Z[Mp$Zp$case.idx[[p]], drop = FALSE]
  #        Zp[[p]]$ROWSUM <- t(this.z)
  #
  #    }
  #    attr(Yp, "Zp") <- Zp
  # }

  Yp
}

# compute sample statistics, per cluster
lav_samplestats_cluster_patterns <- function(Y = NULL, Lp = NULL,
                                             conditional.x = FALSE) {
  # coerce Y to matrix
  Y1 <- as.matrix(Y)
  N <- NROW(Y1)
  P <- NCOL(Y1)

  if (is.null(Lp)) {
    lav_msg_stop(gettext("Lp is NULL"))
  }

  # how many levels?
  nlevels <- length(Lp$cluster) + 1L

  # compute some sample statistics per level
  YLp <- vector("list", length = nlevels)
  for (l in 2:nlevels) {
    ncluster.sizes <- Lp$ncluster.sizes[[l]]
    cluster.size <- Lp$cluster.size[[l]]
    cluster.sizes <- Lp$cluster.sizes[[l]]
    nclusters <- Lp$nclusters[[l]]
    both.idx <- Lp$both.idx[[l]]
    within.idx <- Lp$within.idx[[l]]
    between.idx <- Lp$between.idx[[l]]
    cluster.idx <- Lp$cluster.idx[[l]]
    cluster.size.ns <- Lp$cluster.size.ns[[l]]

    # s <- (N^2 - sum(cluster.size^2)) / (N*(nclusters - 1L))
    # same as
    s <- (N - sum(cluster.size^2) / N) / (nclusters - 1)
    # NOTE: must be (nclusters - 1), otherwise, s is not average cluster
    # size even in the balanced case

    Y1.means <- colMeans(Y1, na.rm = TRUE)
    Y1Y1 <- lav_matrix_crossprod(Y1)
    both.idx <- all.idx <- seq_len(P)

    if (length(within.idx) > 0L ||
      length(between.idx) > 0L) {
      both.idx <- all.idx[-c(within.idx, between.idx)]
      # hm, this assumes the 'order' is the
      # same at both levels...
    }

    # cluster-means
    Y2 <- rowsum.default(Y1,
      group = cluster.idx, reorder = FALSE,
      na.rm = FALSE, # must be FALSE!
    ) / cluster.size
    Y2c <- t(t(Y2) - Y1.means)

    # compute S.w
    # center within variables by grand mean instead of group mean?
    # (YR: apparently not for S.PW)

    Y2a <- Y2
    # if(length(within.idx) > 0L) {
    #    for(i in 1:length(within.idx)) {
    #        Y2a[, within.idx[i]] <- Y1.means[within.idx[i]]
    #    }
    # }
    Y1a <- Y1 - Y2a[cluster.idx, , drop = FALSE]
    S.w <- lav_matrix_crossprod(Y1a) / (N - nclusters)

    # S.b
    # three parts: within/within, between/between, between/within
    # standard definition of the between variance matrix
    # divides by (nclusters - 1)
    S.b <- lav_matrix_crossprod(Y2c * cluster.size, Y2c) / (nclusters - 1)

    # check for zero variances
    if (length(both.idx) > 0L) {
      zero.idx <- which(diag(S.b)[both.idx] < 0.0001)
      if (length(zero.idx) > 0L && !anyNA(Y2)) {
        lav_msg_warn(gettext(
          "(near) zero variance at between level for splitted variable:"),
          paste(Lp$both.names[[l]][zero.idx], collapse = " ")
        )
      }
    }

    S <- cov(Y1, use = "pairwise.complete.obs") * (N - 1L) / N
	# missing by design?
	S[is.na(S)] <- as.numeric(NA)

    # loglik.x
    # extract 'fixed' level-1 loglik from here
    wx.idx <- Lp$ov.x.idx[[1]]
    if (length(wx.idx) > 0L) {
      loglik.x.w <- lav_mvnorm_h1_loglik_samplestats(
        sample.nobs = Lp$nclusters[[1]],
        sample.cov  = S[wx.idx, wx.idx, drop = FALSE]
      )
    } else {
      loglik.x.w <- 0
    }
    # extract 'fixed' level-2 loglik
    bx.idx <- Lp$ov.x.idx[[2]]
    if (length(bx.idx) > 0L) {
      COVB <- cov(Y2[, bx.idx, drop = FALSE]) * (nclusters - 1) / nclusters
      loglik.x.b <- lav_mvnorm_h1_loglik_samplestats(
        sample.nobs = Lp$nclusters[[2]],
        sample.cov  = COVB
      )
    } else {
      loglik.x.b <- 0
    }
    loglik.x <- loglik.x.w + loglik.x.b


    S.PW.start <- S.w
    if (length(within.idx) > 0L) {
      S.PW.start[within.idx, within.idx] <-
        S[within.idx, within.idx, drop = FALSE]
    }

    if (length(between.idx) > 0L) {
      S.w[between.idx, ] <- 0
      S.w[, between.idx] <- 0
      S.PW.start[between.idx, ] <- 0
      S.PW.start[, between.idx] <- 0
    }

    if (length(between.idx) > 0L) {
      # this is what is needed for MUML:
      S.b[, between.idx] <-
        (s * nclusters / N) * S.b[, between.idx, drop = FALSE]
      S.b[between.idx, ] <-
        (s * nclusters / N) * S.b[between.idx, , drop = FALSE]
      S.b[between.idx, between.idx] <-
        (s * lav_matrix_crossprod(
          Y2c[, between.idx, drop = FALSE],
          Y2c[, between.idx, drop = FALSE]
        ) / nclusters)
    }

    Sigma.B <- (S.b - S.w) / s
    Sigma.B[within.idx, ] <- 0
    Sigma.B[, within.idx] <- 0

    # what if we have negative variances in Sigma.B?
    # this may happen if 'split' a variable that has no between variance
    zero.idx <- which(diag(Sigma.B) < 1e-10)
    if (length(zero.idx) > 0L) {
      Sigma.B[zero.idx, ] <- 0
      Sigma.B[, zero.idx] <- 0
    }


    Mu.W <- numeric(P)
    Mu.W[within.idx] <- Y1.means[within.idx]

    Mu.B <- Y1.means
    Mu.B[within.idx] <- 0
    if (length(between.idx) > 0L) {
      # replace between.idx by cov(Y2)[,] elements...
      Mu.B[between.idx] <- colMeans(Y2[, between.idx, drop = FALSE],
        na.rm = TRUE
      )

      S2 <- (cov(Y2, use = "pairwise.complete.obs") *
        (nclusters - 1L) / nclusters)

      Sigma.B[between.idx, between.idx] <-
        S2[between.idx, between.idx, drop = FALSE]
    }

    # FIXME: Mu.B not quite ok for (fixed.x) x variables if they
    # occur both at level 1 AND level 2
    Mu.B.start <- Mu.B
    # Mu.B.start[both.idx] <- Mu.B.start[both.idx] - colMeans(Y2c[,both.idx])

    # sample statistics PER CLUSTER-SIZE

    # summary statistics for complete data, conditional.x = FALSE
    # also needed for h1 (even if conditional.x = TRUE)
    cov.d <- vector("list", length = ncluster.sizes)
    mean.d <- vector("list", length = ncluster.sizes)
    for (clz in seq_len(ncluster.sizes)) {
      nj <- cluster.sizes[clz]
      # select clusters with this size
      d.idx <- which(cluster.size == nj)
      ns <- length(d.idx)
      # NOTE:!!!!
      # reorder columns
      # to match A.inv and m.k later on in objective!!!
      tmp2 <- Y2[d.idx,
        c(between.idx, sort.int(c(both.idx, within.idx))),
        drop = FALSE
      ]
      mean.d[[clz]] <- colMeans(tmp2, na.rm = TRUE)
      bad.idx <- which(!is.finite(mean.d[[clz]])) # if nrow = 1 + NA
      if (length(bad.idx) > 0L) {
        mean.d[[clz]][bad.idx] <- 0 # ugly, only for starting values
      }
      if (length(d.idx) > 1L) {
        if (any(is.na(tmp2))) {
          # if full column has NA, this will fail...
          # not needed anyway
          # out <- lav_mvnorm_missing_h1_estimate_moments(Y = tmp2,
          #          max.iter = 10L)
          # cov.d[[clz]] <- out$Sigma
          cov.d[[clz]] <- 0
        } else {
          cov.d[[clz]] <- (cov(tmp2, use = "complete.obs") *
            (ns - 1) / ns)
        }
      } else {
        cov.d[[clz]] <- 0
      }
    } # clz

    # new in 0.6-12:
    # summary statistics for complete data, conditional.x = TRUE
    # ONLY for twolevel
    if (conditional.x) {
      within.x.idx <- Lp$within.x.idx[[1]]
      between.y.idx <- Lp$between.y.idx[[2]]
      between.x.idx <- Lp$between.x.idx[[2]]
      y1.idx <- Lp$ov.y.idx[[1]]
      x1.idx <- c(within.x.idx, between.x.idx) # in that order

      # data
      Y1.wb <- Y1[, y1.idx, drop = FALSE]
      Y2.wb <- Y2[, y1.idx, drop = FALSE]
      if (length(between.y.idx) > 0L) {
        Y2.z <- Y2[, between.y.idx, drop = FALSE]
      }
      if (length(x1.idx) > 0L) {
        EXO.wb1 <- cbind(1, Y1[, x1.idx, drop = FALSE])
        EXO.wb2 <- cbind(1, Y2[, x1.idx, drop = FALSE])
      } else {
        EXO.wb1 <- matrix(1, nrow(Y1), 1L)
        EXO.wb2 <- matrix(1, nrow(Y2), 1L)
      }

      # sample beta.wb (level 1)
      sample.wb <- solve(crossprod(EXO.wb1), crossprod(EXO.wb1, Y1.wb))
      sample.yhat.wb1 <- EXO.wb1 %*% sample.wb
      sample.yres.wb1 <- Y1.wb - sample.yhat.wb1
      sample.YYres.wb1 <- crossprod(sample.yres.wb1)
      sample.XX.wb1 <- crossprod(EXO.wb1)

      # sample beta.wb (level 2)
      XX.wb2 <- crossprod(EXO.wb2)
      sample.wb2 <- try(solve(XX.wb2, crossprod(EXO.wb2, Y2.wb)),
        silent = TRUE
      )
      if (inherits(sample.wb2, "try-error")) {
        # this may happen if the covariate is cluster-centered
        # using the observed cluster means; then the 'means' will
        # be all (near) zero, and there is no variance
        sample.wb2 <- MASS::ginv(XX.wb2) %*% crossprod(EXO.wb2, Y2.wb)
      }
      sample.yhat.wb2 <- EXO.wb2 %*% sample.wb2
      sample.yres.wb2 <- Y2.wb - sample.yhat.wb2

      # weighted by cluster.size
      sample.YYres.wb2 <- crossprod(
        sample.yres.wb2,
        sample.yres.wb2 * cluster.size
      )
      sample.YresX.wb2 <- crossprod(
        sample.yres.wb2,
        EXO.wb2 * cluster.size
      )
      sample.XX.wb2 <- crossprod(
        EXO.wb2,
        EXO.wb2 * cluster.size
      )


      sample.clz.Y2.res <- vector("list", ncluster.sizes)
      sample.clz.Y2.XX <- vector("list", ncluster.sizes)
      sample.clz.Y2.B <- vector("list", ncluster.sizes)

      if (length(between.y.idx) > 0L) {
        sample.clz.ZZ.res <- vector("list", ncluster.sizes)
        sample.clz.ZZ.XX <- vector("list", ncluster.sizes)
        sample.clz.ZZ.B <- vector("list", ncluster.sizes)

        sample.clz.YZ.res <- vector("list", ncluster.sizes)
        sample.clz.YZ.XX <- vector("list", ncluster.sizes)
        sample.clz.YresXZ <- vector("list", ncluster.sizes)
        sample.clz.XWZres <- vector("list", ncluster.sizes)
      }
      for (clz in seq_len(ncluster.sizes)) {
        # cluster size
        nj <- cluster.sizes[clz]
        nj.idx <- which(cluster.size == nj)

        # Y2
        Y2.clz <- Y2[nj.idx, y1.idx, drop = FALSE]
        if (length(x1.idx) > 0L) {
          EXO2.clz <- cbind(1, Y2[nj.idx, x1.idx, drop = FALSE])
        } else {
          EXO2.clz <- matrix(1, nrow(Y2.clz), 1L)
        }
        XX.clz <- crossprod(EXO2.clz)
        clz.Y2.B <- try(solve(XX.clz, crossprod(EXO2.clz, Y2.clz)),
          silent = TRUE
        )
        if (inherits(clz.Y2.B, "try-error")) {
          clz.Y2.B <-
            MASS::ginv(XX.clz) %*% crossprod(EXO2.clz, Y2.clz)
        }
        clz.Y2.hat <- EXO2.clz %*% clz.Y2.B
        clz.Y2.res <- Y2.clz - clz.Y2.hat
        sample.clz.Y2.B[[clz]] <- clz.Y2.B
        sample.clz.Y2.res[[clz]] <- crossprod(clz.Y2.res)
        sample.clz.Y2.XX[[clz]] <- crossprod(EXO2.clz)

        # Z
        if (length(between.y.idx) > 0L) {
          Z.clz.z <- Y2[nj.idx, between.y.idx, drop = FALSE]
          if (length(between.x.idx) > 0L) {
            EXO.clz.z <- cbind(
              1,
              Y2[nj.idx, between.x.idx, drop = FALSE]
            )
          } else {
            EXO.clz.z <- matrix(1, nrow(Z.clz.z), 1L)
          }
          ZZ.clz <- crossprod(EXO.clz.z)
          clz.ZZ.B <- try(
            solve(
              ZZ.clz,
              crossprod(EXO.clz.z, Z.clz.z)
            ),
            silent = TRUE
          )
          if (inherits(clz.ZZ.B, "try-error")) {
            clz.ZZ.B <-
              MASS::ginv(ZZ.clz) %*% crossprod(EXO.clz.z, Z.clz.z)
          }
          clz.Z.hat <- EXO.clz.z %*% clz.ZZ.B
          clz.Z.res <- Z.clz.z - clz.Z.hat
          sample.clz.ZZ.B[[clz]] <- clz.ZZ.B
          sample.clz.ZZ.res[[clz]] <- crossprod(clz.Z.res)
          sample.clz.ZZ.XX[[clz]] <- crossprod(EXO.clz.z)

          sample.clz.YZ.res[[clz]] <- crossprod(clz.Y2.res, clz.Z.res)
          sample.clz.YZ.XX[[clz]] <- crossprod(EXO2.clz, EXO.clz.z)
          sample.clz.YresXZ[[clz]] <- crossprod(clz.Y2.res, EXO.clz.z)
          sample.clz.XWZres[[clz]] <- crossprod(EXO2.clz, clz.Z.res)
        }
      } # clz
    } # conditional.x

    YLp[[l]] <- list(
      Y1Y1 = Y1Y1,
      Y2 = Y2, s = s, S.b = S.b, S.PW.start = S.PW.start,
      Sigma.W = S.w, Mu.W = Mu.W,
      Sigma.B = Sigma.B, Mu.B = Mu.B,
      Mu.B.start = Mu.B.start, loglik.x = loglik.x,
      mean.d = mean.d, cov.d = cov.d
    )

    # if conditional, add more stuff
    if (conditional.x) {
      if (length(between.y.idx) > 0L) {
        extra <- list(
          sample.wb = sample.wb,
          sample.YYres.wb1 = sample.YYres.wb1,
          sample.XX.wb1 = sample.XX.wb1,
          sample.wb2 = sample.wb2,
          sample.YYres.wb2 = sample.YYres.wb2,
          sample.YresX.wb2 = sample.YresX.wb2,
          sample.XX.wb2 = sample.XX.wb2,
          sample.clz.Y2.res = sample.clz.Y2.res,
          sample.clz.Y2.XX = sample.clz.Y2.XX,
          sample.clz.Y2.B = sample.clz.Y2.B,
          sample.clz.ZZ.res = sample.clz.ZZ.res,
          sample.clz.ZZ.XX = sample.clz.ZZ.XX,
          sample.clz.ZZ.B = sample.clz.ZZ.B,
          sample.clz.YZ.res = sample.clz.YZ.res,
          sample.clz.YZ.XX = sample.clz.YZ.XX,
          sample.clz.YresXZ = sample.clz.YresXZ, # zero?
          sample.clz.XWZres = sample.clz.XWZres
        )
      } else {
        extra <- list(
          sample.wb = sample.wb,
          sample.YYres.wb1 = sample.YYres.wb1,
          sample.XX.wb1 = sample.XX.wb1,
          sample.wb2 = sample.wb2,
          sample.YYres.wb2 = sample.YYres.wb2,
          sample.YresX.wb2 = sample.YresX.wb2,
          sample.XX.wb2 = sample.XX.wb2,
          sample.clz.Y2.res = sample.clz.Y2.res,
          sample.clz.Y2.XX = sample.clz.Y2.XX,
          sample.clz.Y2.B = sample.clz.Y2.B
        )
      }
      YLp[[l]] <- c(YLp[[l]], extra)
    }
  } # l

  YLp
}
