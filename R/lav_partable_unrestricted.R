# YR - 26 Nov 2013: generate partable for the unrestricted model
# YR - 19 Mar 2017: handle twolevel model
# YR - 27 May 2021: added lav_pt_unrestricted_chol so we can use
#                   a cholesky parameterization: S = LAMBDA %*% t(LAMBDA)

lav_pt_unrestricted <- function(lavobject = NULL,          # nolint start
                                      # if no object is available,
                                      lavdata = NULL,
                                      lavpta = NULL,
                                      lavoptions = NULL,
                                      lavsamplestats = NULL,
                                      lavh1 = NULL,
                                      # optional user-provided sample stats
                                      sample_cov = NULL,
                                      sample_mean = NULL,
                                      sample_slopes = NULL,
                                      sample_th = NULL,
                                      sample_th_idx = NULL,
                                      sample_cov_x = NULL,
                                      sample_mean_x = NULL) {   # nolint end
  lav_pt_indep_or_unrestricted(
    lavobject = lavobject,
    lavdata = lavdata, lavpta = lavpta, lavoptions = lavoptions,
    lavsamplestats = lavsamplestats, lavh1 = lavh1,
    sample_cov = sample_cov, sample_mean = sample_mean,
    sample_slopes = sample_slopes,
    sample_th = sample_th, sample_th_idx = sample_th_idx,
    independent = FALSE
  )
}

# generate parameter table for an independence model
# YR - 12 Sep 2017: special case of lav_pt_unrestricted()
lav_pt_independence <- function(lavobject = NULL,         # nolint start
                                      # if no object is available,
                                      lavdata = NULL,
                                      lavpta = NULL,
                                      lavoptions = NULL,
                                      lavsamplestats = NULL,
                                      lavh1 = NULL,
                                      # optional user-provided sample stats
                                      sample_cov = NULL,
                                      sample_mean = NULL,
                                      sample_slopes = NULL,
                                      sample_th = NULL,
                                      sample_th_idx = NULL,
                                      sample_cov_x = NULL,
                                      sample_mean_x = NULL) {   # nolint end
  lav_pt_indep_or_unrestricted(
    lavobject = lavobject,
    lavdata = lavdata, lavpta = lavpta, lavoptions = lavoptions,
    lavsamplestats = lavsamplestats, lavh1 = lavh1,
    sample_cov = sample_cov, sample_mean = sample_mean,
    sample_slopes = sample_slopes,
    sample_th = sample_th, sample_th_idx = sample_th_idx,
    independent = TRUE
  )
}

lav_pt_indep_or_unrestricted <- function(lavobject = NULL,
                                               # if no object is available,
                                               lavdata = NULL,
                                               lavpta = NULL,
                                               lavoptions = NULL,
                                               lavsamplestats = NULL,
                                               lavh1 = NULL,
                                           # optional user-provided sample stats
                                               sample_cov = NULL,
                                               sample_mean = NULL,
                                               sample_slopes = NULL,
                                               sample_th = NULL,
                                               sample_th_idx = NULL,
                                               sample_cov_x = NULL,
                                               sample_mean_x = NULL,
                                               independent = FALSE) {
  # grab everything from lavaan lavobject
  if (!is.null(lavobject)) {
    stopifnot(inherits(lavobject, "lavaan"))

    lavdata <- lavobject@Data
    lavoptions <- lavobject@Options
    lavsamplestats <- lavobject@SampleStats
    lavpta <- lavobject@pta
    lavh1 <- lavobject@h1
  }

  if (lavdata@data.type == "none") {
    lavsamplestats <- NULL
  }

  # conditional.x ? check res.cov[[1]] slot
  conditional_x <- FALSE
  if (!is.null(lavsamplestats) && !is.null(lavsamplestats@res.cov[[1]])) {
    conditional_x <- TRUE
  } else if (!is.null(lavoptions) && lavoptions$conditional.x) {
    conditional_x <- TRUE
  }

  # group.w.free?
  group_w_free <- FALSE
  if (!is.null(lavoptions) && lavoptions$group.w.free) {
    group_w_free <- TRUE
  }
  # we use CAPS below for the list version, so we can use 'small caps'
  # within the for() loop

  # get sample statistics, all groups
  sample_cov_1 <- sample_cov
  if (is.null(sample_cov_1) && !is.null(lavsamplestats)) {
    if (conditional_x) {
      sample_cov_1 <- lavsamplestats@res.cov
    } else {
      sample_cov_1 <- lavsamplestats@cov
    }
  }

  sample_mean_1 <- sample_mean
  if (is.null(sample_mean_1) && !is.null(lavsamplestats)) {
    if (conditional_x) {
      sample_mean_1 <- lavsamplestats@res.int
    } else {
      sample_mean_1 <- lavsamplestats@mean
    }
  }

  sample_slopes_1 <- sample_slopes
  if (conditional_x && is.null(sample_slopes_1) && !is.null(lavsamplestats)) {
    sample_slopes_1 <- lavsamplestats@res.slopes
  }

  sample_th_1 <- sample_th
  if (is.null(sample_th_1) && !is.null(lavsamplestats)) {
    if (conditional_x) {
      sample_th_1 <- lavsamplestats@res.th
    } else {
      sample_th_1 <- lavsamplestats@th
    }
  }

  sample_th_idx_1 <- sample_th_idx
  if (is.null(sample_th_idx_1) && !is.null(lavsamplestats)) {
    sample_th_idx_1 <- lavsamplestats@th.idx
  }

  sample_cov_x_1 <- sample_cov_x
  if (is.null(sample_cov_x_1) && !is.null(lavsamplestats)) {
    sample_cov_x_1 <- lavsamplestats@cov.x
  }

  sample_mean_x_1 <- sample_mean_x
  if (is.null(sample_mean_x_1) && !is.null(lavsamplestats)) {
    sample_mean_x_1 <- lavsamplestats@mean.x
  }



  ov <- lavdata@ov
  meanstructure <- lavoptions$meanstructure
  categorical <- any(ov$type == "ordered")
  ngroups <- lavdata@ngroups
  nlevels <- lavdata@nlevels
  if (lavoptions$estimator == "catML") {
    categorical <- FALSE
  }
  correlation <- FALSE
  if (!is.null(lavoptions$correlation)) {
    correlation <- lavoptions$correlation
  }

  # what with fixed.x?
  # - does not really matter; fit will be saturated anyway
  # - fixed.x = TRUE may avoid convergence issues with non-numeric
  #             x-covariates
  fixed_x <- lavoptions$fixed.x

  # if multilevel
  if (nlevels > 1L) {
    # fixed.x       <- FALSE # for now
    if (!lavoptions$estimator %in% c("WLS", "DWLS", "ULS")) {
      categorical <- FALSE # for now
      conditional_x <- FALSE # for now
    }
    # else: two-level (D)WLS supports conditional.x (slopes per block)
  }

  lhs <- rhs <- op <- character(0)
  group <- block <- level <- free <- exo <- integer(0)
  ustart <- numeric(0)

  # block number
  b <- 0L
  for (g in 1:ngroups) {
    # only for multilevel
    if (nlevels > 1L) {
      # ylp <- lavsamplestats@YLp[[g]]
      lp <- lavdata@Lp[[g]]
    }

    # local copy
    sample_cov <- sample_cov_1[[g]]
    sample_mean <- sample_mean_1[[g]]
    sample_slopes <- sample_slopes_1[[g]]
    sample_th <- sample_th_1[[g]]
    sample_th_idx <- sample_th_idx_1[[g]]
    sample_cov_x <- sample_cov_x_1[[g]]
    sample_mean_x <- sample_mean_x_1[[g]]

    # force local sample.cov to be pd -- just for starting values anyway
    if (!is.null(sample_cov) && !anyNA(sample_cov)) {
      sample_cov <- lav_mat_sym_force_pd(sample_cov)
    }

    for (l in 1:nlevels) {
      # block
      b <- b + 1L

      # ov.names for this block
      if (is.null(lavpta)) { # only data was used
      if (nlevels > 1L) {
      ov_names <- lavdata@ov.names.l[[g]][[(ngroups - 1) * g + b]]
    } else {
          ov_names <- lavdata@ov.names[[g]]
    }
        ov_names_x <- lavdata@ov.names.x[[g]]
        ov_names_nox <- ov_names[!ov_names %in% ov_names_x]
      } else {
        if (conditional_x) {
          ov_names <- lavpta$vnames$ov.nox[[b]]
        } else {
          ov_names <- lavpta$vnames$ov[[b]]
        }
        ov_names_x <- lavpta$vnames$ov.x[[b]]
        ov_names_nox <- lavpta$vnames$ov.nox[[b]]
      }

      # only for multilevel, overwrite sample.cov and sample_mean
      if (nlevels > 1L) {
        if (independent) {
          ov_names_x <- lavpta$vnames$ov.x[[b]]
          ov_names_nox <- lavpta$vnames$ov.nox[[b]]
          if (!is.null(lavh1$implied$cov.x) &&
              !is.null(lavh1$implied$cov.x[[b]]) &&
              NROW(lavh1$implied$cov.x[[b]]) > 0L) {
            # conditional.x with y-only h1 cov blocks and separately
            # stored covariate moments (two-level (D)WLS)
            sample_cov_x <- lavh1$implied$cov.x[[b]]
            sample_mean_x <- lavh1$implied$mean.x[[b]]
          } else if (length(ov_names_x) > 0L) {
            # better use lavdata@Lp[[g]]$ov.x.idx??
            # in case we have x/y mismatch across levels?
            ov_x_idx <- lavpta$vidx$ov.x[[b]]
            sample_cov_x <- lavh1$implied$cov[[b]][ov_x_idx,
              ov_x_idx,
              drop = FALSE
            ]
            sample_mean_x <- lavh1$implied$mean[[b]][ov_x_idx]
          } else {
            sample_cov_x <- NULL
            sample_mean_x <- NULL
          }
        } else {
          ov_names_x <- character(0L)
          ov_names_nox <- ov_names
        }

        if (length(lavh1) > 0L) {
          sample_cov <- lavh1$implied$cov[[b]]
          sample_mean <- lavh1$implied$mean[[b]]
          if (conditional_x && !is.null(lavh1$implied$res.slopes) &&
              !is.null(lavh1$implied$res.slopes[[b]])) {
            sample_slopes <- lavh1$implied$res.slopes[[b]]
          }
        } else {
          sample_cov <- diag(length(ov_names))
          sample_mean <- numeric(length(ov_names))
        }

        # if(l == 1L) {
        #    sample.cov  <- YLp[[2]]$Sigma.W[block.idx, block.idx,
        #                                    drop = FALSE]
        #    sample_mean <- YLp[[2]]$Mu.W[block.idx]
        # } else {
        #    sample.cov  <- YLp[[2]]$Sigma.B[block.idx, block.idx,
        #                                    drop = FALSE]
        #    sample_mean <- YLp[[2]]$Mu.B[block.idx]
        # }

        # force local sample.cov to be strictly pd (and exaggerate)
        # just for starting values anyway, but at least the first
        # evaluation will be feasible
        sample_cov <- lav_mat_sym_force_pd(sample_cov,
          tol = 1e-03
        )
      }


      # a) VARIANCES (all ov's, if !conditional.x, also exo's)
      nvar <- length(ov_names)

      lhs <- c(lhs, ov_names)
      op <- c(op, rep("~~", nvar))
      rhs <- c(rhs, ov_names)
      block <- c(block, rep(b, nvar))
      group <- c(group, rep(g, nvar))
      level <- c(level, rep(l, nvar))
      if (correlation) {
        free <- c(free, rep(0L, nvar))
      } else {
        free <- c(free, rep(1L, nvar))
      }
      exo <- c(exo, rep(0L, nvar))

      # starting values -- variances
      if (correlation) {
        ustart <- c(ustart, rep(1, nvar))
      } else if (!is.null(sample_cov)) {
        ustart <- c(ustart, diag(sample_cov))
      } else {
        ustart <- c(ustart, rep(as.numeric(NA), nvar))
      }

      # COVARIANCES!
      if (!independent) {
        pstar <- nvar * (nvar - 1) / 2
        if (pstar > 0L) { # only if more than 1 variable
          tmp <- utils::combn(ov_names, 2)
          lhs <- c(lhs, tmp[1, ]) # to fill upper.tri
          op <- c(op, rep("~~", pstar))
          rhs <- c(rhs, tmp[2, ])
          block <- c(block, rep(b, pstar))
          group <- c(group, rep(g, pstar))
          level <- c(level, rep(l, pstar))
          free <- c(free, rep(1L, pstar))
          exo <- c(exo, rep(0L, pstar))
        }

        # starting values -- covariances
        if (!is.null(sample_cov)) {
          sample_cov_vech <- lav_mat_vech(sample_cov, diagonal = FALSE)
          ustart <- c(ustart, sample_cov_vech)
          # check for 'missing by design' cells: here, the sample.cov
          # element is *exactly* zero (new in 0.6-18)
          zero_cov <- which(sample_cov_vech == 0)
          if (length(zero_cov) > 0L && !is.null(lavh1)) {
            n_tmp <- length(free)
            ones_and_zeroes <- rep(1L, pstar)
            ones_and_zeroes[zero_cov] <- 0L
            free[(n_tmp - pstar + 1):n_tmp] <- ones_and_zeroes
          }
        } else {
          ustart <- c(ustart, rep(as.numeric(NA), pstar))
        }
      }

      # ordered? fix variances, add thresholds
      ord_names <- character(0L)
      if (categorical) {
        ord_names <- ov$name[ov$type == "ordered"]
        # only for this group
        ord_names <- ov_names[which(ov_names %in% ord_names)]

        if (length(ord_names) > 0L) {
          # fix variances to 1.0
          # (two-level, theta parameterization: at the WITHIN level only;
          #  the between-level variances of the ordinal variables are free)
          if (nlevels > 1L) {
            if (l == 1L) {
              idx <- which(lhs %in% ord_names & op == "~~" & lhs == rhs &
                           block == b)
              ustart[idx] <- 1.0
              free[idx] <- 0L
            }
          } else {
            idx <- which(lhs %in% ord_names & op == "~~" & lhs == rhs)
            ustart[idx] <- 1.0
            free[idx] <- 0L
          }

          # add thresholds
          # (two-level: FREE at the between level; fixed-to-zero
          #  placeholder rows at the within level -- the same convention
          #  as lavaanify() for the user model)
          lhs_th <- character(0)
          rhs_th <- character(0)
          for (o in ord_names) {
            nth <- ov$nlev[ov$name == o] - 1L
            if (nth < 1L) next
            lhs_th <- c(lhs_th, rep(o, nth))
            rhs_th <- c(rhs_th, paste("t", seq_len(nth), sep = ""))
          }
          nel <- length(lhs_th)
          lhs <- c(lhs, lhs_th)
          rhs <- c(rhs, rhs_th)
          op <- c(op, rep("|", nel))
          block <- c(block, rep(b, nel))
          group <- c(group, rep(g, nel))
          level <- c(level, rep(l, nel))
          exo <- c(exo, rep(0L, nel))
          if (nlevels == 1L || l == nlevels) {
            free <- c(free, rep(1L, nel))
            # starting values
            if (!is.null(sample_th) && !is.null(sample_th_idx)) {
              th_start <- sample_th[sample_th_idx > 0L]
              ustart <- c(ustart, th_start)
            } else {
              ustart <- c(ustart, rep(as.numeric(NA), nel))
            }
          } else {
            # within level: fixed-to-zero placeholders
            free <- c(free, rep(0L, nel))
            ustart <- c(ustart, rep(0, nel))
          } # thresholds

          # fixed-to-zero intercepts (since 0.5.17)
          ov_int <- ord_names
          nel <- length(ov_int)
          lhs <- c(lhs, ov_int)
          op <- c(op, rep("~1", nel))
          rhs <- c(rhs, rep("", nel))
          block <- c(block, rep(b, nel))
          group <- c(group, rep(g, nel))
          level <- c(level, rep(l, nel))
          free <- c(free, rep(0L, nel))
          exo <- c(exo, rep(0L, nel))
          ustart <- c(ustart, rep(0, nel))

          # ~*~ (since 0.6-1)
          nel <- length(ov_int)
          lhs <- c(lhs, ov_int)
          op <- c(op, rep("~*~", nel))
          rhs <- c(rhs, ov_int)
          block <- c(block, rep(b, nel))
          group <- c(group, rep(g, nel))
          level <- c(level, rep(l, nel))
          free <- c(free, rep(0L, nel))
          exo <- c(exo, rep(0L, nel))
          ustart <- c(ustart, rep(1, nel))
        }
      } # categorical

      # correlation structure?
      if (!categorical && correlation) {
        nel <- nvar
        lhs <- c(lhs, ov_names)
        op <- c(op, rep("~*~", nel))
        rhs <- c(rhs, ov_names)
        block <- c(block, rep(b, nel))
        group <- c(group, rep(g, nel))
        level <- c(level, rep(l, nel))
        free <- c(free, rep(0L, nel))
        exo <- c(exo, rep(0L, nel))
        ustart <- c(ustart, rep(1, nel))
      }

      # meanstructure?
      if (meanstructure) {
        # auto-remove ordinal variables
        ov_int <- ov_names
        idx <- which(ov_int %in% ord_names)
        if (length(idx)) ov_int <- ov_int[-idx]

        nel <- length(ov_int)
        lhs <- c(lhs, ov_int)
        op <- c(op, rep("~1", nel))
        rhs <- c(rhs, rep("", nel))
        block <- c(block, rep(b, nel))
        group <- c(group, rep(g, nel))
        level <- c(level, rep(l, nel))
        # if multilevel, level=1 has fixed zeroes
        if (nlevels > 1L && l == 1L) {
          within_1 <- rep(0L, nel)
          # within-only variables have a free mean at level 1; match by
          # NAME: ov_int may exclude ordinal (and exogenous) variables,
          # while within.idx refers to the full ('tilde') variable set
          # FIXME: assuming 1 group
          within_names <- lp$ov.names[lp$within.idx[[2]]]
          within_1[ov_int %in% within_names] <- 1L
          free <- c(free, within_1)
        } else {
          free <- c(free, rep(1L, nel))
        }
        exo <- c(exo, rep(0L, nel))

        # starting values
        if (!is.null(sample_mean)) {
          sample_int_idx <- match(ov_int, ov_names)
          ustart <- c(ustart, sample_mean[sample_int_idx])
        } else {
          ustart <- c(ustart, rep(as.numeric(NA), length(ov_int)))
        }
      }


      # fixed.x exogenous variables?
      if (!conditional_x && (nx <- length(ov_names_x)) > 0L) {
        if (independent && lavoptions$baseline.fixed.x.free.cov) {
          # add covariances for eXo
          pstar <- nx * (nx - 1) / 2
          if (pstar > 0L) { # only if more than 1 variable
            tmp <- utils::combn(ov_names_x, 2)
            lhs <- c(lhs, tmp[1, ]) # to fill upper.tri
            op <- c(op, rep("~~", pstar))
            rhs <- c(rhs, tmp[2, ])
            block <- c(block, rep(b, pstar))
            group <- c(group, rep(g, pstar))
            level <- c(level, rep(l, pstar))
            free <- c(free, rep(1L, pstar))
            exo <- c(exo, rep(0L, pstar))

            # starting values
            if (!is.null(sample_cov_x)) {
              rhs_idx <- match(tmp[1, ], ov_names_x)
              lhs_idx <- match(tmp[2, ], ov_names_x)
              ustart <- c(
                ustart,
                sample_cov_x[cbind(rhs_idx, lhs_idx)]
              )
            } else {
              ustart <- c(ustart, rep(as.numeric(NA), pstar))
            }
          }
        }

        if (fixed_x) {
          # fix variances/covariances
          exo_idx <- which(rhs %in% ov_names_x &
            lhs %in% ov_names_x &
            op == "~~" & group == g) # ok
          exo[exo_idx] <- 1L
          free[exo_idx] <- 0L

          # fix means
          exo_idx <- which(lhs %in% ov_names_x &
            op == "~1" & group == g) # ok
          exo[exo_idx] <- 1L
          free[exo_idx] <- 0L
        }
      }

      # conditional.x?
      if (conditional_x && (nx <- length(ov_names_x)) > 0L) {
        # eXo variances
        nel <- length(ov_names_x)
        lhs <- c(lhs, ov_names_x)
        op <- c(op, rep("~~", nel))
        rhs <- c(rhs, ov_names_x)
        block <- c(block, rep(b, nel))
        group <- c(group, rep(g, nel))
        level <- c(level, rep(l, nel))
        if (fixed_x) {
          free <- c(free, rep(0L, nel))
          exo <- c(exo, rep(1L, nel))
        } else {
          free <- c(free, rep(1L, nel))
          exo <- c(exo, rep(0L, nel))
        }

        # starting values
        if (!is.null(sample_cov_x)) {
          ustart <- c(ustart, diag(sample_cov_x))
        } else {
          ustart <- c(ustart, rep(as.numeric(NA), nel))
        }


        # eXo covariances
        pstar <- nx * (nx - 1) / 2
        if (pstar > 0L) { # only if more than 1 variable
          tmp <- utils::combn(ov_names_x, 2)
          lhs <- c(lhs, tmp[1, ]) # to fill upper.tri
          op <- c(op, rep("~~", pstar))
          rhs <- c(rhs, tmp[2, ])
          block <- c(block, rep(b, pstar))
          group <- c(group, rep(g, pstar))
          level <- c(level, rep(l, pstar))
          if (fixed_x) {
            free <- c(free, rep(0L, pstar))
            exo <- c(exo, rep(1L, pstar))
          } else {
            free <- c(free, rep(1L, pstar))
            exo <- c(exo, rep(0L, pstar))
          }

          # starting values
          if (!is.null(sample_cov_x)) {
            rhs_idx <- match(tmp[1, ], ov_names_x)
            lhs_idx <- match(tmp[2, ], ov_names_x)
            ustart <- c(
              ustart,
              sample_cov_x[cbind(rhs_idx, lhs_idx)]
            )
          } else {
            ustart <- c(ustart, rep(as.numeric(NA), pstar))
          }
        }

        # eXo means
        if (meanstructure) {
          ov_int <- ov_names_x

          nel <- length(ov_int)
          lhs <- c(lhs, ov_int)
          op <- c(op, rep("~1", nel))
          rhs <- c(rhs, rep("", nel))
          group <- c(group, rep(g, nel))
          block <- c(block, rep(b, nel))
          level <- c(level, rep(l, nel))
          if (fixed_x) {
            free <- c(free, rep(0L, nel))
            exo <- c(exo, rep(1L, nel))
          } else {
            free <- c(free, rep(1L, nel))
            exo <- c(exo, rep(0L, nel))
          }

          # starting values
          if (!is.null(sample_mean_x)) {
            sample_int_idx <- match(ov_int, ov_names_x)
            ustart <- c(ustart, sample_mean_x[sample_int_idx])
          } else {
            ustart <- c(ustart, rep(as.numeric(NA), length(ov_int)))
          }
        }

        # slopes
        nnox <- length(ov_names_nox)
        nel <- nnox * nx

        lhs <- c(lhs, rep(ov_names_nox, times = nx))
        op <- c(op, rep("~", nel))
        rhs <- c(rhs, rep(ov_names_x, each = nnox))
        block <- c(block, rep(b, nel))
        group <- c(group, rep(g, nel))
        level <- c(level, rep(l, nel))
        if (independent) {
          if (lavoptions$baseline.conditional.x.free.slopes) {
            free <- c(free, rep(1L, nel))
          } else {
            free <- c(free, rep(0L, nel))
          }
        } else {
          free <- c(free, rep(1L, nel))
        }
        exo <- c(exo, rep(1L, nel))

        # starting values -- slopes
        if (independent) {
          # FIXME: zero slope-structure provides a fit that
          # is equal to the conditional.x = FALSE version;
          # in principle, we could just fix the slope-structure
          # to the sample-based slopes

          # to get the old behaviour:
          if (!lavoptions$baseline.conditional.x.free.slopes) {
            ustart <- c(ustart, rep(0, nel))
          } else {
            # but we probably should do:
            ustart <- c(ustart, lav_mat_vec(sample_slopes))
          }
        } else if (!is.null(sample_slopes)) {
          ustart <- c(ustart, lav_mat_vec(sample_slopes))
        } else {
          ustart <- c(ustart, rep(as.numeric(NA), nel))
        }
      } # conditional.x

      # group.w.free (new in 0.6-8)
      if (group_w_free) {
        lhs <- c(lhs, "group")
        op <- c(op, "%")
        rhs <- c(rhs, "w")
        block <- c(block, b)
        group <- c(group, g)
        level <- c(level, l)
        free <- c(free, 1L)
        exo <- c(exo, 0L)
        ustart <- c(ustart, lavsamplestats@WLS.obs[[g]][1])
      }
    } # levels
  } # ngroups

  # free counter
  idx_free <- which(free > 0)
  free[idx_free] <- seq_along(idx_free)

  list_1 <- list(
    id = seq_along(lhs),
    lhs = lhs,
    op = op,
    rhs = rhs,
    user = rep(1L, length(lhs)),
    block = block,
    group = group,
    level = level,
    free = free,
    ustart = ustart,
    exo = exo # ,
    # label       = rep("",  length(lhs))
    # eq.id       = rep(0L,  length(lhs)),
    # unco        = free
  )


  # keep level column if no levels? (no for now)
  if (nlevels < 2L) {
    list_1$level <- NULL
  }

  list_1
}

# - currently only used for continuous twolevel data
# - conditional.x not supported (yet)
lav_pt_unrestricted_chol <- function(lavobject = NULL,
                                           # if no object is available,
                                           lavdata = NULL,
                                           lavpta = NULL, # optional
                                           lavoptions = NULL,
                                           group = NULL) {
  # grab everything from lavaan lavobject
  if (!is.null(lavobject)) {
    stopifnot(inherits(lavobject, "lavaan"))

    lavdata <- lavobject@Data
    lavoptions <- lavobject@Options
    # lavsamplestats <- lavobject@SampleStats
    lavpta <- lavobject@pta
    # lavh1 <- lavobject@h1
  }

  ov <- lavdata@ov
  meanstructure <- lavoptions$meanstructure
  categorical <- any(ov$type == "ordered")
  if (categorical) {
    lav_msg_stop(gettext("categorical data not supported in this function"))
  }
  ngroups <- lavdata@ngroups
  nlevels <- lavdata@nlevels

  # select groups
  if (is.null(group)) {
    group_select <- seq_len(ngroups)
  } else {
    stopifnot(is.numeric(group), all(group <= ngroups))
    group_select <- group
    if (length(group_select) == 0L) {
      lav_msg_stop(gettext("no groups selected"))
    }
  }

  # what with fixed.x?
  # - does not really matter; fit will be saturated anyway
  # - fixed.x = TRUE may avoid convergence issues with non-numeric
  #             x-covariates
  # fixed_x <- lavoptions$fixed.x

  # if multilevel
  if (nlevels > 1L) {
    # fixed.x       <- FALSE # for now
    # conditional_x <- FALSE # for now
    categorical <- FALSE # for now
  }

  lhs <- rhs <- op <- character(0)
  group <- block <- level <- free <- exo <- integer(0)
  ustart <- lower <- numeric(0)

  # block number
  b <- 0L
  for (g in 1:ngroups) {

    # select group?
    if (! g %in% group_select) {
      next
    }

    # only for multilevel
    if (nlevels > 1L) {
      lp <- lavdata@Lp[[g]]
    }

    for (l in 1:nlevels) {
      # block
      b <- b + 1L

      if (is.null(lavpta)) {
        ov_names <- lavdata@ov.names[[b]]
        # ov_names_x <- lavdata@ov.names.x[[b]]
        # ov_names_nox <- ov_names[!ov_names %in% ov_names_x]
      } else {
        ov_names <- lavpta$vnames$ov[[b]]
        # ov_names_x <- lavpta$vnames$ov.x[[b]]
        # ov_names_nox <- lavpta$vnames$ov.nox[[b]]
      }

      # only for multilevel, overwrite sample.cov and sample_mean
      if (nlevels > 1L) {
        # ov_names_x <- character(0L)
        # ov_names_nox <- ov_names
      }

      # create lv.names == ov.names
      lv_names <- paste("f", ov_names, sep = "")

      # a) OV VARIANCES -> fixed to zero
      nvar <- length(ov_names)
      lhs <- c(lhs, ov_names)
      op <- c(op, rep("~~", nvar))
      rhs <- c(rhs, ov_names)
      block <- c(block, rep(b, nvar))
      group <- c(group, rep(g, nvar))
      level <- c(level, rep(l, nvar))
      ustart <- c(ustart, rep(0.0001, nvar)) ### Force PD!! (option?)
      free <- c(free, rep(0L, nvar))
      exo <- c(exo, rep(0L, nvar))
      lower <- c(lower, rep(0.0, nvar))

      # b) LV VARIANCES -> fixed to 1.0
      nvar <- length(lv_names)
      lhs <- c(lhs, lv_names)
      op <- c(op, rep("~~", nvar))
      rhs <- c(rhs, lv_names)
      block <- c(block, rep(b, nvar))
      group <- c(group, rep(g, nvar))
      level <- c(level, rep(l, nvar))
      ustart <- c(ustart, rep(1.0, nvar))
      free <- c(free, rep(0L, nvar))
      exo <- c(exo, rep(0L, nvar))
      lower <- c(lower, rep(1.0, nvar))

      # c) LOADINGS self
      nvar <- length(ov_names)
      lhs <- c(lhs, lv_names)
      op <- c(op, rep("=~", nvar))
      rhs <- c(rhs, ov_names)
      block <- c(block, rep(b, nvar))
      group <- c(group, rep(g, nvar))
      level <- c(level, rep(l, nvar))
      ustart <- c(ustart, rep(as.numeric(NA), nvar))
      free <- c(free, rep(1L, nvar))
      exo <- c(exo, rep(0L, nvar))
      lower <- c(lower, rep(0.0, nvar)) # lower bound!

      # d) LOADINGS other
      if (length(ov_names) > 1L) {
        tmp <- utils::combn(ov_names, 2)
        pstar <- ncol(tmp)
        lhs <- c(lhs, paste("f", tmp[1, ], sep = ""))
        op <- c(op, rep("=~", pstar))
        rhs <- c(rhs, tmp[2, ])
        block <- c(block, rep(b, pstar))
        group <- c(group, rep(g, pstar))
        level <- c(level, rep(l, pstar))
        free <- c(free, rep(1L, pstar))
        exo <- c(exo, rep(0L, pstar))
        lower <- c(lower, rep(-Inf, pstar))
        ustart <- c(ustart, rep(as.numeric(NA), pstar))
      }

      # check for zero coverage at level 1 (new in 0.6-18)
      if (lavdata@missing %in% c("ml", "ml.x") && l == 1 &&
          !is.null(lavdata@Mp[[g]])) {
        coverage <- lavdata@Mp[[g]]$coverage
        sample_cov_vech <- lav_mat_vech(coverage, diagonal = FALSE)
        zero_cov <- which(sample_cov_vech == 0)
        if (length(zero_cov) > 0L) {
          n_tmp <- length(free)
          ones_and_zeroes <- rep(1L, pstar)
          ones_and_zeroes[zero_cov] <- 0L
      inf_and_zeroes <- rep(-Inf, pstar)
      inf_and_zeroes[zero_cov] <- 0
      na_and_zeroes <- rep(as.numeric(NA), pstar)
      na_and_zeroes[zero_cov] <- 0
          free[(n_tmp - pstar + 1):n_tmp] <- ones_and_zeroes
      ustart[(n_tmp - pstar + 1):n_tmp] <- na_and_zeroes
      lower[(n_tmp - pstar + 1):n_tmp] <- inf_and_zeroes
        }
      }

      # meanstructure?
      if (meanstructure) {
        # OV
        ov_int <- ov_names

        nel <- length(ov_int)
        lhs <- c(lhs, ov_int)
        op <- c(op, rep("~1", nel))
        rhs <- c(rhs, rep("", nel))
        block <- c(block, rep(b, nel))
        group <- c(group, rep(g, nel))
        level <- c(level, rep(l, nel))
        # if multilevel, level=1 has fixed zeroes
        if (nlevels > 1L && l == 1L) {
          within_1 <- rep(0L, nel)
          within_idx <- match(lp$within.idx[[2]], lp$ov.idx[[1]])
          within_1[within_idx] <- 1L
          free <- c(free, within_1)
        } else {
          free <- c(free, rep(1L, nel))
        }
        exo <- c(exo, rep(0L, nel))
        lower <- c(lower, rep(-Inf, nel))
        ustart <- c(ustart, rep(as.numeric(NA), nel))

        # LV
        ov_int <- lv_names

        nel <- length(ov_int)
        lhs <- c(lhs, ov_int)
        op <- c(op, rep("~1", nel))
        rhs <- c(rhs, rep("", nel))
        block <- c(block, rep(b, nel))
        group <- c(group, rep(g, nel))
        level <- c(level, rep(l, nel))
        free <- c(free, rep(0L, nel))
        exo <- c(exo, rep(0L, nel))
        ustart <- c(ustart, rep(0.0, nel))
        lower <- c(lower, rep(-Inf, nel))
      }
    } # levels
  } # ngroups

  # free counter
  idx_free <- which(free > 0)
  free[idx_free] <- seq_along(idx_free)

  list_1 <- list(
    id = seq_along(lhs),
    lhs = lhs,
    op = op,
    rhs = rhs,
    user = rep(1L, length(lhs)),
    block = block,
    group = group,
    level = level,
    free = free,
    ustart = ustart,
    exo = exo,
    lower = lower # ,
    # label       = rep("",  length(lhs))
    # eq.id       = rep(0L,  length(lhs)),
    # unco        = free
  )


  # keep level column if no levels? (no for now)
  if (nlevels < 2L) {
    list_1$level <- NULL
  }

  list_1
}

# create a 'baseline' model from an existing object/partable
# - fix all (free) directed effects (factor loadings, regressions) to zero
# - fix all (free) covariances to zero
# - neutralize the (measured) latent variables by fixing their variances
#   to zero
# - keep all other (relevant) constraints
# - this *should* result in a baseline model that is always nested within
#   the original model
lav_pt_baseline <- function(lavobject = NULL,
                                  # if no object is available,
                                  lavpartable = NULL,
                                  lavh1 = NULL) {

  # grab everything from lavaan lavobject
  if (!is.null(lavobject)) {
    stopifnot(inherits(lavobject, "lavaan"))
    lavpartable <- lavobject@ParTable
    lavh1 <- lavobject@h1
  }

  # number of blocks
  nblocks <- lav_pt_nblocks(lavpartable)

  # conditional.x ? check res.cov[[1]] slot
  conditional_x <- FALSE
  if (!is.null(lavh1) && !is.null(lavh1$implied$res.cov[[1]])) {
    conditional_x <- TRUE
  }

  # get sample statistics, all groups
  if (conditional_x) {
    sample_cov <- lavh1$implied$res.cov
    # sample_mean <- lavh1$implied$res.int
  } else {
    sample_cov <- lavh1$implied$cov
    # sample_mean <- lavh1$implied$mean
  }

  # shortcut of lavpartable
  pt_1 <- lavpartable

  # keep exo=1 starting values
  exo_idx <- which(pt_1$exo == 1L)
  if (!is.null(pt_1$est)) {
    pt_1$ustart[exo_idx] <- pt_1$est[exo_idx]
  } else if (!is.null(pt_1$start)) {
    pt_1$iustart[exo_idx] <- pt_1$start[exo_idx]
  }

  # remove est/se columns, if present, but keep start column
  pt_1$est <- pt_1$se <- NULL

  # lv.names
  lv_names <- lav_pt_vnames(pt_1, "lv")

  # zero-out all directed effects and covariances
  directed_idx <- which(pt_1$op %in% c("=~", "~") & pt_1$free > 0L &
                                                            pt_1$exo != 1L)
  cov_idx <- which(pt_1$op == "~~" & pt_1$lhs != pt_1$rhs & pt_1$free > 0L)

  # neutralize latent variables
  lv_var_idx <- which(pt_1$op == "~~" & pt_1$lhs %in% lv_names &
                                        pt_1$lhs == pt_1$rhs & pt_1$free > 0L)
  #lv.int.idx <- which(PT$op == "~1" & PT$lhs %in% lv.names & PT$free > 0L)
  #zero.idx <- c(directed.idx, cov.idx, lv.var.idx, lv.int.idx)
  zero_idx <- c(directed_idx, cov_idx, lv_var_idx)

  pt_1$free[zero_idx] <- 0L
  pt_1$start[zero_idx] <- 0
  pt_1$ustart[zero_idx] <- 0

  # Question: fill in more elements in PT$start? (only ov variances for now)
  for (b in seq_len(nblocks)) {
    ov_names_num <- lav_pt_vnames(lavpartable, "ov.num", block = b)
    if (conditional_x) {
      ov_names_x <- lav_pt_vnames(lavpartable, "ov.x", block = b)
      ov_names_num <- ov_names_num[!ov_names_num %in% ov_names_x]
    }
    ovar_idx <- which(lavpartable$block == b &
      lavpartable$op == "~~" &
      lavpartable$lhs %in% ov_names_num &
      lavpartable$lhs == lavpartable$rhs)
    sample_var_idx <- match(lavpartable$lhs[ovar_idx], ov_names_num)
    pt_1$ustart[ovar_idx] <- diag(sample_cov[[b]])[sample_var_idx]
  }

  pt_1 <- lav_pt_complete(pt_1)
  pt_1
}
