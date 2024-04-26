# YR - 26 Nov 2013: generate partable for the unrestricted model
# YR - 19 Mar 2017: handle twolevel model
# YR - 27 May 2021: added lav_partable_unrestricted_chol so we can use
#                   a cholesky parameterization: S = LAMBDA %*% t(LAMBDA)

lav_partable_unrestricted <- function(lavobject = NULL,
                                      # if no object is available,
                                      lavdata = NULL,
                                      lavpta = NULL,
                                      lavoptions = NULL,
                                      lavsamplestats = NULL,
                                      lavh1 = NULL,
                                      # optional user-provided sample stats
                                      sample.cov = NULL,
                                      sample.mean = NULL,
                                      sample.slopes = NULL,
                                      sample.th = NULL,
                                      sample.th.idx = NULL,
                                      sample.cov.x = NULL,
                                      sample.mean.x = NULL) {
  lav_partable_indep_or_unrestricted(
    lavobject = lavobject,
    lavdata = lavdata, lavpta = lavpta, lavoptions = lavoptions,
    lavsamplestats = lavsamplestats, lavh1 = lavh1,
    sample.cov = sample.cov, sample.mean = sample.mean,
    sample.slopes = sample.slopes,
    sample.th = sample.th, sample.th.idx = sample.th.idx,
    independent = FALSE
  )
}

# generate parameter table for an independence model
# YR - 12 Sep 2017: special case of lav_partable_unrestricted()
lav_partable_independence <- function(lavobject = NULL,
                                      # if no object is available,
                                      lavdata = NULL,
                                      lavpta = NULL,
                                      lavoptions = NULL,
                                      lavsamplestats = NULL,
                                      lavh1 = NULL,
                                      # optional user-provided sample stats
                                      sample.cov = NULL,
                                      sample.mean = NULL,
                                      sample.slopes = NULL,
                                      sample.th = NULL,
                                      sample.th.idx = NULL,
                                      sample.cov.x = NULL,
                                      sample.mean.x = NULL) {
  lav_partable_indep_or_unrestricted(
    lavobject = lavobject,
    lavdata = lavdata, lavpta = lavpta, lavoptions = lavoptions,
    lavsamplestats = lavsamplestats, lavh1 = lavh1,
    sample.cov = sample.cov, sample.mean = sample.mean,
    sample.slopes = sample.slopes,
    sample.th = sample.th, sample.th.idx = sample.th.idx,
    independent = TRUE
  )
}

lav_partable_indep_or_unrestricted <- function(lavobject = NULL,
                                               # if no object is available,
                                               lavdata = NULL,
                                               lavpta = NULL,
                                               lavoptions = NULL,
                                               lavsamplestats = NULL,
                                               lavh1 = NULL,
                                               # optional user-provided sample stats
                                               sample.cov = NULL,
                                               sample.mean = NULL,
                                               sample.slopes = NULL,
                                               sample.th = NULL,
                                               sample.th.idx = NULL,
                                               sample.cov.x = NULL,
                                               sample.mean.x = NULL,
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
  conditional.x <- FALSE
  if (!is.null(lavsamplestats) && !is.null(lavsamplestats@res.cov[[1]])) {
    conditional.x <- TRUE
  } else if (!is.null(lavoptions) && lavoptions$conditional.x) {
    conditional.x <- TRUE
  }

  # group.w.free?
  group.w.free <- FALSE
  if (!is.null(lavoptions) && lavoptions$group.w.free) {
    group.w.free <- TRUE
  }
  # we use CAPS below for the list version, so we can use 'small caps'
  # within the for() loop

  # get sample statistics, all groups
  SAMPLE.cov <- sample.cov
  if (is.null(SAMPLE.cov) && !is.null(lavsamplestats)) {
    if (conditional.x) {
      SAMPLE.cov <- lavsamplestats@res.cov
    } else {
      SAMPLE.cov <- lavsamplestats@cov
    }
  }

  SAMPLE.mean <- sample.mean
  if (is.null(SAMPLE.mean) && !is.null(lavsamplestats)) {
    if (conditional.x) {
      SAMPLE.mean <- lavsamplestats@res.int
    } else {
      SAMPLE.mean <- lavsamplestats@mean
    }
  }

  SAMPLE.slopes <- sample.slopes
  if (conditional.x && is.null(SAMPLE.slopes) && !is.null(lavsamplestats)) {
    SAMPLE.slopes <- lavsamplestats@res.slopes
  }

  SAMPLE.th <- sample.th
  if (is.null(SAMPLE.th) && !is.null(lavsamplestats)) {
    if (conditional.x) {
      SAMPLE.th <- lavsamplestats@res.th
    } else {
      SAMPLE.th <- lavsamplestats@th
    }
  }

  SAMPLE.th.idx <- sample.th.idx
  if (is.null(SAMPLE.th.idx) && !is.null(lavsamplestats)) {
    SAMPLE.th.idx <- lavsamplestats@th.idx
  }

  SAMPLE.cov.x <- sample.cov.x
  if (is.null(SAMPLE.cov.x) && !is.null(lavsamplestats)) {
    SAMPLE.cov.x <- lavsamplestats@cov.x
  }

  SAMPLE.mean.x <- sample.mean.x
  if (is.null(SAMPLE.mean.x) && !is.null(lavsamplestats)) {
    SAMPLE.mean.x <- lavsamplestats@mean.x
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
  fixed.x <- lavoptions$fixed.x

  # if multilevel
  if (nlevels > 1L) {
    # fixed.x       <- FALSE # for now
    conditional.x <- FALSE # for now
    categorical <- FALSE # for now
  }

  lhs <- rhs <- op <- character(0)
  group <- block <- level <- free <- exo <- integer(0)
  ustart <- numeric(0)

  # block number
  b <- 0L
  for (g in 1:ngroups) {
    # only for multilevel
    if (nlevels > 1L) {
      YLp <- lavsamplestats@YLp[[g]]
      Lp <- lavdata@Lp[[g]]
    }

    # local copy
    sample.cov <- SAMPLE.cov[[g]]
    sample.mean <- SAMPLE.mean[[g]]
    sample.slopes <- SAMPLE.slopes[[g]]
    sample.th <- SAMPLE.th[[g]]
    sample.th.idx <- SAMPLE.th.idx[[g]]
    sample.cov.x <- SAMPLE.cov.x[[g]]
    sample.mean.x <- SAMPLE.mean.x[[g]]

    # force local sample.cov to be pd -- just for starting values anyway
    if (!is.null(sample.cov) && !anyNA(sample.cov)) {
      sample.cov <- lav_matrix_symmetric_force_pd(sample.cov)
    }

    for (l in 1:nlevels) {
      # block
      b <- b + 1L

      # ov.names for this block
      if (is.null(lavpta)) { # only data was used
        ov.names <- lavdata@ov.names[[g]]
        ov.names.x <- lavdata@ov.names.x[[g]]
        ov.names.nox <- ov.names[!ov.names %in% ov.names.x]
      } else {
        if (conditional.x) {
          ov.names <- lavpta$vnames$ov.nox[[b]]
        } else {
          ov.names <- lavpta$vnames$ov[[b]]
        }
        ov.names.x <- lavpta$vnames$ov.x[[b]]
        ov.names.nox <- lavpta$vnames$ov.nox[[b]]
      }

      # only for multilevel, overwrite sample.cov and sample.mean
      if (nlevels > 1L) {
        if (independent) {
          # beter use lavdata@Lp[[g]]$ov.x.idx??
          # in case we have x/y mismatch across levels?
          ov.x.idx <- lavpta$vidx$ov.x[[b]]
          ov.names.x <- lavpta$vnames$ov.x[[b]]
          ov.names.nox <- lavpta$vnames$ov.nox[[b]]
          sample.cov.x <- lavh1$implied$cov[[b]][ov.x.idx,
            ov.x.idx,
            drop = FALSE
          ]
          sample.mean.x <- lavh1$implied$mean[[b]][ov.x.idx]
        } else {
          ov.names.x <- character(0L)
          ov.names.nox <- ov.names
        }

        if (length(lavh1) > 0L) {
          sample.cov <- lavh1$implied$cov[[b]]
          sample.mean <- lavh1$implied$mean[[b]]
        } else {
          sample.cov <- diag(length(ov.names))
          sample.mean <- numeric(length(ov.names))
        }

        # if(l == 1L) {
        #    sample.cov  <- YLp[[2]]$Sigma.W[block.idx, block.idx,
        #                                    drop = FALSE]
        #    sample.mean <- YLp[[2]]$Mu.W[block.idx]
        # } else {
        #    sample.cov  <- YLp[[2]]$Sigma.B[block.idx, block.idx,
        #                                    drop = FALSE]
        #    sample.mean <- YLp[[2]]$Mu.B[block.idx]
        # }

        # force local sample.cov to be strictly pd (and exaggerate)
        # just for starting values anyway, but at least the first
        # evaluation will be feasible
        sample.cov <- lav_matrix_symmetric_force_pd(sample.cov,
          tol = 1e-03
        )
      }


      # a) VARIANCES (all ov's, if !conditional.x, also exo's)
      nvar <- length(ov.names)

      lhs <- c(lhs, ov.names)
      op <- c(op, rep("~~", nvar))
      rhs <- c(rhs, ov.names)
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
      } else if (!is.null(sample.cov)) {
        ustart <- c(ustart, diag(sample.cov))
      } else {
        ustart <- c(ustart, rep(as.numeric(NA), nvar))
      }

      # COVARIANCES!
      if (!independent) {
        pstar <- nvar * (nvar - 1) / 2
        if (pstar > 0L) { # only if more than 1 variable
          tmp <- utils::combn(ov.names, 2)
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
        if (!is.null(sample.cov)) {
          ustart <- c(ustart, lav_matrix_vech(sample.cov,
            diagonal = FALSE
          ))
        } else {
          ustart <- c(ustart, rep(as.numeric(NA), pstar))
        }
      }

      # ordered? fix variances, add thresholds
      ord.names <- character(0L)
      if (categorical) {
        ord.names <- ov$name[ov$type == "ordered"]
        # only for this group
        ord.names <- ov.names[which(ov.names %in% ord.names)]

        if (length(ord.names) > 0L) {
          # fix variances to 1.0
          idx <- which(lhs %in% ord.names & op == "~~" & lhs == rhs)
          ustart[idx] <- 1.0
          free[idx] <- 0L

          # add thresholds
          lhs.th <- character(0)
          rhs.th <- character(0)
          for (o in ord.names) {
            nth <- ov$nlev[ov$name == o] - 1L
            if (nth < 1L) next
            lhs.th <- c(lhs.th, rep(o, nth))
            rhs.th <- c(rhs.th, paste("t", seq_len(nth), sep = ""))
          }
          nel <- length(lhs.th)
          lhs <- c(lhs, lhs.th)
          rhs <- c(rhs, rhs.th)
          op <- c(op, rep("|", nel))
          block <- c(block, rep(b, nel))
          group <- c(group, rep(g, nel))
          level <- c(level, rep(l, nel))
          free <- c(free, rep(1L, nel))
          exo <- c(exo, rep(0L, nel))

          # starting values
          if (!is.null(sample.th) && !is.null(sample.th.idx)) {
            th.start <- sample.th[sample.th.idx > 0L]
            ustart <- c(ustart, th.start)
          } else {
            ustart <- c(ustart, rep(as.numeric(NA), nel))
          }

          # fixed-to-zero intercepts (since 0.5.17)
          ov.int <- ord.names
          nel <- length(ov.int)
          lhs <- c(lhs, ov.int)
          op <- c(op, rep("~1", nel))
          rhs <- c(rhs, rep("", nel))
          block <- c(block, rep(b, nel))
          group <- c(group, rep(g, nel))
          level <- c(level, rep(l, nel))
          free <- c(free, rep(0L, nel))
          exo <- c(exo, rep(0L, nel))
          ustart <- c(ustart, rep(0, nel))

          # ~*~ (since 0.6-1)
          nel <- length(ov.int)
          lhs <- c(lhs, ov.int)
          op <- c(op, rep("~*~", nel))
          rhs <- c(rhs, ov.int)
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
        lhs <- c(lhs, ov.names)
        op <- c(op, rep("~*~", nel))
        rhs <- c(rhs, ov.names)
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
        ov.int <- ov.names
        idx <- which(ov.int %in% ord.names)
        if (length(idx)) ov.int <- ov.int[-idx]

        nel <- length(ov.int)
        lhs <- c(lhs, ov.int)
        op <- c(op, rep("~1", nel))
        rhs <- c(rhs, rep("", nel))
        block <- c(block, rep(b, nel))
        group <- c(group, rep(g, nel))
        level <- c(level, rep(l, nel))
        # if multilevel, level=1 has fixed zeroes
        if (nlevels > 1L && l == 1L) {
          WITHIN <- rep(0L, nel)
          # FIXME: assuming 1 group
          within.idx <- match(Lp$within.idx[[2]], Lp$ov.idx[[1]])
          WITHIN[within.idx] <- 1L
          free <- c(free, WITHIN)
        } else {
          free <- c(free, rep(1L, nel))
        }
        exo <- c(exo, rep(0L, nel))

        # starting values
        if (!is.null(sample.mean)) {
          sample.int.idx <- match(ov.int, ov.names)
          ustart <- c(ustart, sample.mean[sample.int.idx])
        } else {
          ustart <- c(ustart, rep(as.numeric(NA), length(ov.int)))
        }
      }


      # fixed.x exogenous variables?
      if (!conditional.x && (nx <- length(ov.names.x)) > 0L) {
        if (independent && lavoptions$mimic %in% c("Mplus", "lavaan")) {
          # add covariances for eXo
          pstar <- nx * (nx - 1) / 2
          if (pstar > 0L) { # only if more than 1 variable
            tmp <- utils::combn(ov.names.x, 2)
            lhs <- c(lhs, tmp[1, ]) # to fill upper.tri
            op <- c(op, rep("~~", pstar))
            rhs <- c(rhs, tmp[2, ])
            block <- c(block, rep(b, pstar))
            group <- c(group, rep(g, pstar))
            level <- c(level, rep(l, pstar))
            free <- c(free, rep(1L, pstar))
            exo <- c(exo, rep(0L, pstar))

            # starting values
            if (!is.null(sample.cov.x)) {
              rhs.idx <- match(tmp[1, ], ov.names.x)
              lhs.idx <- match(tmp[2, ], ov.names.x)
              ustart <- c(
                ustart,
                sample.cov.x[cbind(rhs.idx, lhs.idx)]
              )
            } else {
              ustart <- c(ustart, rep(as.numeric(NA), pstar))
            }
          }
        }

        if (fixed.x) {
          # fix variances/covariances
          exo.idx <- which(rhs %in% ov.names.x &
            lhs %in% ov.names.x &
            op == "~~" & group == g) # ok
          exo[exo.idx] <- 1L
          free[exo.idx] <- 0L

          # fix means
          exo.idx <- which(lhs %in% ov.names.x &
            op == "~1" & group == g) # ok
          exo[exo.idx] <- 1L
          free[exo.idx] <- 0L
        }
      }

      # conditional.x?
      if (conditional.x && (nx <- length(ov.names.x)) > 0L) {
        # eXo variances
        nel <- length(ov.names.x)
        lhs <- c(lhs, ov.names.x)
        op <- c(op, rep("~~", nel))
        rhs <- c(rhs, ov.names.x)
        block <- c(block, rep(b, nel))
        group <- c(group, rep(g, nel))
        level <- c(level, rep(l, nel))
        if (fixed.x) {
          free <- c(free, rep(0L, nel))
          exo <- c(exo, rep(1L, nel))
        } else {
          free <- c(free, rep(1L, nel))
          exo <- c(exo, rep(0L, nel))
        }

        # starting values
        if (!is.null(sample.cov.x)) {
          ustart <- c(ustart, diag(sample.cov.x))
        } else {
          ustart <- c(ustart, rep(as.numeric(NA), nel))
        }


        # eXo covariances
        pstar <- nx * (nx - 1) / 2
        if (pstar > 0L) { # only if more than 1 variable
          tmp <- utils::combn(ov.names.x, 2)
          lhs <- c(lhs, tmp[1, ]) # to fill upper.tri
          op <- c(op, rep("~~", pstar))
          rhs <- c(rhs, tmp[2, ])
          block <- c(block, rep(b, pstar))
          group <- c(group, rep(g, pstar))
          level <- c(level, rep(l, pstar))
          if (fixed.x) {
            free <- c(free, rep(0L, pstar))
            exo <- c(exo, rep(1L, pstar))
          } else {
            free <- c(free, rep(1L, pstar))
            exo <- c(exo, rep(0L, pstar))
          }

          # starting values
          if (!is.null(sample.cov.x)) {
            rhs.idx <- match(tmp[1, ], ov.names.x)
            lhs.idx <- match(tmp[2, ], ov.names.x)
            ustart <- c(
              ustart,
              sample.cov.x[cbind(rhs.idx, lhs.idx)]
            )
          } else {
            ustart <- c(ustart, rep(as.numeric(NA), pstar))
          }
        }

        # eXo means
        if (meanstructure) {
          ov.int <- ov.names.x

          nel <- length(ov.int)
          lhs <- c(lhs, ov.int)
          op <- c(op, rep("~1", nel))
          rhs <- c(rhs, rep("", nel))
          group <- c(group, rep(g, nel))
          block <- c(block, rep(b, nel))
          level <- c(level, rep(l, nel))
          if (fixed.x) {
            free <- c(free, rep(0L, nel))
            exo <- c(exo, rep(1L, nel))
          } else {
            free <- c(free, rep(1L, nel))
            exo <- c(exo, rep(0L, nel))
          }

          # starting values
          if (!is.null(sample.mean.x)) {
            sample.int.idx <- match(ov.int, ov.names.x)
            ustart <- c(ustart, sample.mean.x[sample.int.idx])
          } else {
            ustart <- c(ustart, rep(as.numeric(NA), length(ov.int)))
          }
        }

        # slopes
        nnox <- length(ov.names.nox)
        nel <- nnox * nx

        lhs <- c(lhs, rep(ov.names.nox, times = nx))
        op <- c(op, rep("~", nel))
        rhs <- c(rhs, rep(ov.names.x, each = nnox))
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
            ustart <- c(ustart, lav_matrix_vec(sample.slopes))
          }
        } else if (!is.null(sample.slopes)) {
          ustart <- c(ustart, lav_matrix_vec(sample.slopes))
        } else {
          ustart <- c(ustart, rep(as.numeric(NA), nel))
        }
      } # conditional.x

      # group.w.free (new in 0.6-8)
      if (group.w.free) {
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
  idx.free <- which(free > 0)
  free[idx.free] <- 1:length(idx.free)

  LIST <- list(
    id = 1:length(lhs),
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
    LIST$level <- NULL
  }

  LIST
}

# - currently only used for continuous twolevel data
# - conditional.x not supported (yet)
lav_partable_unrestricted_chol <- function(lavobject = NULL,
                                           # if no object is available,
                                           lavdata = NULL,
                                           lavpta = NULL,
                                           lavoptions = NULL) {
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

  # what with fixed.x?
  # - does not really matter; fit will be saturated anyway
  # - fixed.x = TRUE may avoid convergence issues with non-numeric
  #             x-covariates
  fixed.x <- lavoptions$fixed.x

  # if multilevel
  if (nlevels > 1L) {
    # fixed.x       <- FALSE # for now
    conditional.x <- FALSE # for now
    categorical <- FALSE # for now
  }

  lhs <- rhs <- op <- character(0)
  group <- block <- level <- free <- exo <- integer(0)
  ustart <- lower <- numeric(0)

  # block number
  b <- 0L
  for (g in 1:ngroups) {
    # only for multilevel
    if (nlevels > 1L) {
      Lp <- lavdata@Lp[[g]]
    }

    for (l in 1:nlevels) {
      # block
      b <- b + 1L

      if (is.null(lavpta)) {
        ov.names <- lavdata@ov.names[[b]]
        ov.names.x <- lavdata@ov.names.x[[b]]
        ov.names.nox <- ov.names[!ov.names %in% ov.names.x]
      } else {
        ov.names <- lavpta$vnames$ov[[b]]
        ov.names.x <- lavpta$vnames$ov.x[[b]]
        ov.names.nox <- lavpta$vnames$ov.nox[[b]]
      }

      # only for multilevel, overwrite sample.cov and sample.mean
      if (nlevels > 1L) {
        ov.names.x <- character(0L)
        ov.names.nox <- ov.names
      }

      # create lv.names == ov.names
      lv.names <- paste("f", ov.names, sep = "")

      # a) OV VARIANCES -> fixed to zero
      nvar <- length(ov.names)
      lhs <- c(lhs, ov.names)
      op <- c(op, rep("~~", nvar))
      rhs <- c(rhs, ov.names)
      block <- c(block, rep(b, nvar))
      group <- c(group, rep(g, nvar))
      level <- c(level, rep(l, nvar))
      ustart <- c(ustart, rep(0.0001, nvar)) ### Force PD!! (option?)
      free <- c(free, rep(0L, nvar))
      exo <- c(exo, rep(0L, nvar))
      lower <- c(lower, rep(0.0, nvar))

      # b) LV VARIANCES -> fixed to 1.0
      nvar <- length(lv.names)
      lhs <- c(lhs, lv.names)
      op <- c(op, rep("~~", nvar))
      rhs <- c(rhs, lv.names)
      block <- c(block, rep(b, nvar))
      group <- c(group, rep(g, nvar))
      level <- c(level, rep(l, nvar))
      ustart <- c(ustart, rep(1.0, nvar))
      free <- c(free, rep(0L, nvar))
      exo <- c(exo, rep(0L, nvar))
      lower <- c(lower, rep(1.0, nvar))

      # c) LOADINGS self
      nvar <- length(ov.names)
      lhs <- c(lhs, lv.names)
      op <- c(op, rep("=~", nvar))
      rhs <- c(rhs, ov.names)
      block <- c(block, rep(b, nvar))
      group <- c(group, rep(g, nvar))
      level <- c(level, rep(l, nvar))
      ustart <- c(ustart, rep(as.numeric(NA), nvar))
      free <- c(free, rep(1L, nvar))
      exo <- c(exo, rep(0L, nvar))
      lower <- c(lower, rep(0.0, nvar)) # lower bound!

      # d) LOADINGS other
      if (length(ov.names) > 1L) {
        tmp <- utils::combn(ov.names, 2)
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

      # meanstructure?
      if (meanstructure) {
        # OV
        ov.int <- ov.names

        nel <- length(ov.int)
        lhs <- c(lhs, ov.int)
        op <- c(op, rep("~1", nel))
        rhs <- c(rhs, rep("", nel))
        block <- c(block, rep(b, nel))
        group <- c(group, rep(g, nel))
        level <- c(level, rep(l, nel))
        # if multilevel, level=1 has fixed zeroes
        if (nlevels > 1L && l == 1L) {
          WITHIN <- rep(0L, nel)
          within.idx <- match(Lp$within.idx[[2]], Lp$ov.idx[[1]])
          WITHIN[within.idx] <- 1L
          free <- c(free, WITHIN)
        } else {
          free <- c(free, rep(1L, nel))
        }
        exo <- c(exo, rep(0L, nel))
        lower <- c(lower, rep(-Inf, nel))
        ustart <- c(ustart, rep(as.numeric(NA), nel))

        # LV
        ov.int <- lv.names

        nel <- length(ov.int)
        lhs <- c(lhs, ov.int)
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
  idx.free <- which(free > 0)
  free[idx.free] <- 1:length(idx.free)

  LIST <- list(
    id = 1:length(lhs),
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
    LIST$level <- NULL
  }

  LIST
}
