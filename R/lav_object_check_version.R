# the entries of the lavaan and lavaanList objects have changed over time
# this function will check if the lavaan/lavaanList object is up to date, and
# adapt, if not.
#
# this may be useful if an older (say 0.5) lavaan object was saved, and
# passed to a function like lavPredict() in, say, lavaan 0.6-21.

# notes:
# - pre<0.5 lavaan objects are no longer supported
# - @Fit slot is ignored (as not used anymore)

# YR 16 Oct 2025 + LDW 22 Oct 2025

lav_object_check_version <- function(object = NULL) {

  is.lavaan.object <- inherits(object, "lavaan")
  if (!is.lavaan.object) {
    # check if lavaanList object, if not return input object
    if (!inherits(object, "lavaanList")) return(object)
  }

  # flag: check or not?
  check_not_needed_flag <- TRUE

  # do we have a version slot?
  if(.hasSlot(object, "version")) {
    has_version_flag <- TRUE
    lavobject_version <- object@version
    lavaanpkg_version <- version <- read.dcf(
      file = system.file("DESCRIPTION", package = "lavaan"),
      fields = "Version"
    )[1]
    if (lavobject_version != lavaanpkg_version) {
      check_not_needed_flag <- FALSE
    }
  } else {
    # <0.6
    has_version_flag <- FALSE
    check_not_needed_flag <- FALSE
  }

  # check needed?
  if (check_not_needed_flag) {
    return(object)
  }

  # ok, we have potentially an older (saved) lavaan or lavaanList object
  # check needed slots, and if missing, add them
  suppressWarnings(lavobject <- object)
  ngroups <- lav_partable_ngroups(lavobject@ParTable)
  nblocks <- lav_partable_nblocks(lavobject@ParTable)

  if (!has_version_flag) { # pre 0.6 object!
    # 0.5-10 (25 Oct 2012)
    if (!.hasSlot(lavobject@Data, "group")) {
      lavobject@Data@group <- character(0L)
    }

    # 0.5-11 (19 dec 2012)
    if (is.lavaan.object) {
      if (!.hasSlot(lavobject@SampleStats, "bifreq")) {
        lavobject@SampleStats@bifreq <- vector("list", length = ngroups)
      }
      if (!.hasSlot(lavobject, "Cache")) {
        lavobject@Cache <- list()
      }
    }

    # 0.5-12 (8 March 2013)
    if (is.lavaan.object) {
      if (!.hasSlot(lavobject@SampleStats, "ridge")) {
        lavobject@SampleStats@ridge <- 0
      }
    }

    # 0.5-14 (21 July 2013)
    if (is.lavaan.object) {
      if (!.hasSlot(lavobject@Model, "ov.x.dummy.ov.idx")) {
        lavobject@Model@ov.x.dummy.ov.idx <- vector("list", length = nblocks)
        lavobject@Model@ov.x.dummy.lv.idx <- vector("list", length = nblocks)
        lavobject@Model@ov.y.dummy.ov.idx <- vector("list", length = nblocks)
        lavobject@Model@ov.y.dummy.lv.idx <- vector("list", length = nblocks)
      }
      if (!.hasSlot(lavobject@SampleStats, "mean.x")) {
        lavobject@SampleStats@mean.x <- vector("list", length = ngroups)
        for (g in seq_len(ngroups)) {
          if (!is.null(lavobject@SampleStats@x.idx[[g]])) {
            lavobject@SampleStats@mean.x[[g]] <-
              lavobject@SampleStats@mean[[g]][lavobject@SampleStats@x.idx[[g]]]
          }
        }
      }
      if (!.hasSlot(lavobject, "pta")) {
        lavobject@pta <- list()
      }
    }

    # 0.5-15 (15 Nov 2013)
    if (is.lavaan.object) {
      if (!.hasSlot(lavobject@Data, "Rp")) {
        lavobject@Data@Rp <- vector("list", length = ngroups)
      }
    }

    # 0.5-16 (7 March 2014)
    if (is.lavaan.object) {
      if (!.hasSlot(lavobject@SampleStats, "group.w")) {
        lavobject@SampleStats@group.w <-  vector("list", length = ngroups)
        for (g in seq_len(ngroups)) {
          lavobject@SampleStats@group.w[[g]] <-
            lavobject@SampleStats@nobs[[g]] / lavobject@SampleStats@ntotal
        }
      }
      if (!.hasSlot(lavobject@Model, "group.w.free")) {
        lavobject@Model@group.w.free <- FALSE
      }
      if (!.hasSlot(lavobject@Model, "parameterization")) {
        lavobject@Model@parameterization <- "delta"
      }
      if (!.hasSlot(lavobject@Model, "link")) {
        lavobject@Model@link <- "default"
      }
    }

    # 0.5-17 (30 Sept 2014)
    if (is.lavaan.object) {
      if (!.hasSlot(lavobject@SampleStats, "WLS.VD")) {
        lavobject@SampleStats@WLS.VD <- vector("list", length = ngroups)
      }
    }

    # 0.5-18 (18 Nov 2014)
    if (is.lavaan.object) {
      if (!.hasSlot(lavobject@Model, "eq.constraints.k0")) {
        lavobject@Model@eq.constraints.k0 <- numeric(0L)
      }
      if (!.hasSlot(lavobject@Model, "ceq.linear.idx")) {
        lavobject@Model@ceq.linear.idx <- integer(0L)
        lavobject@Model@ceq.nonlinear.idx <- integer(0L)
        lavobject@Model@cin.linear.idx <- integer(0L)
        lavobject@Model@cin.nonlinear.idx <- integer(0L)
      }
    }

    # 0.5-18 (13 Jan 2015)
    if (is.lavaan.object) {
      if (!.hasSlot(lavobject@Model, "ceq.JAC")) {
        lavobject@Model@ceq.JAC <- matrix(0, nrow = 0L,
                                            ncol = lavobject@Model@nx.free)
        lavobject@Model@ceq.rhs <- numeric(0L)
        lavobject@Model@cin.JAC <- matrix(0, nrow = 0L,
                                            ncol = lavobject@Model@nx.free)
        lavobject@Model@cin.rhs <- numeric(0L)
      }
    }

    # 0.5-19 (30 Jul 2015)
    if (is.lavaan.object) {
      if (!.hasSlot(lavobject, "boot")) {
        # construct partial optim list
        optim_list <- list(x = lavobject@Fit@x,
                           dx = numeric(0L),
                           npar = lavobject@Fit@npar,
                           iterations = lavobject@Fit@iterations,
                           converged = lavobject@Fit@converged,
                           warn.txt = "",
                           parscale = rep(1, lavobject@Fit@npar),
                           fx = lavobject@Fit@fx,
                           fx.group = lavobject@Fit@fx.group,
                           logl.group = lavobject@Fit@logl.group,
                           control = lavobject@Fit@control)
        lavobject@boot <- vector("list", 0L)
        lavobject@optim <- optim_list
        lavobject@implied <- lav_model_implied(lavobject@Model)
        lavobject@vcov <- list(se = lavobject@Options$se[1],
                               information = lavobject@Options$information[1],
                               vcov = matrix(0, 0, 0)) # for now
        lavobject@test <- lavTest(lavobject)
        lavobject@external <- vector("list", 0L)
      }
    }

    # 0.5-19: est/se move to @ParTable
    if (!is.null(lavobject@Fit@est) && is.null(lavobject@ParTable$est)) {
      lavobject@ParTable$est <- lavobject@Fit@est
    }
    if (!is.null(lavobject@Fit@se) && is.null(lavobject@ParTable$se)) {
      lavobject@ParTable$se <- lavobject@Fit@se
    }


    # 0.5-21 (16 Dec 2015)
    if (is.lavaan.object) {
      if (!.hasSlot(lavobject@Model, "conditional.x")) {
        lavobject@Model@conditional.x <- FALSE
      }
    }

    # 0.5-21 (5 Jan 2016)
    if (is.lavaan.object) {
      if (!.hasSlot(lavobject@SampleStats, "x.idx")) {
        lavobject@SampleStats@x.idx <- rep(list(integer(0L)), ngroups)
      }
    }

    # 0.5-21 (8 Jan 2016)
    if (is.lavaan.object) {
      if (!.hasSlot(lavobject@SampleStats, "res.cov")) {
        lavobject@SampleStats@res.cov <- vector("list", ngroups)
        lavobject@SampleStats@res.var <- vector("list", ngroups)
        lavobject@SampleStats@res.th <- vector("list", ngroups)
        lavobject@SampleStats@res.th.nox <- vector("list", ngroups)
        lavobject@SampleStats@res.slopes <- vector("list", ngroups)
        lavobject@SampleStats@res.int <- vector("list", ngroups)
        lavobject@SampleStats@res.icov <- vector("list", ngroups)
        lavobject@SampleStats@res.icov.log.det <- vector("list", ngroups)
      }
    }

    # 0.5-21 (28 Mar 2016)
    if (is.lavaan.object) {
      if (!.hasSlot(lavobject@SampleStats, "NACOV.user")) {
        lavobject@SampleStats@NACOV.user <- FALSE
      }
    }

    #### 0.5-21, 3 Jul 2016, Class lavaanList is added ####

    # 0.5-23 (25 Jan 2017)
    if (!.hasSlot(lavobject@Model, "estimator")) {
      lavobject@Model@estimator <- lavobject@Options$estimator
    }

    # 0.5-23 (30 Jan 2017)
    if (is.lavaan.object) {
      if (!.hasSlot(lavobject@Data, "cluster")) {
        lavobject@Data@cluster <- character(0L)
        lavobject@Data@ordered <- character(0L)
      }
    } else {
      for (j in seq_along(lavobject@DataList)) {
        if (!.hasSlot(lavobject@DataList[[j]], "cluster")) {
          lavobject@DataList[[j]]@cluster <- character(0L)
          lavobject@DataList[[j]]@ordered <- character(0L)
        }
      }
    }

    # 0.5-23 (7 Feb 2017)
    if (is.lavaan.object) {
      if (!.hasSlot(lavobject@SampleStats, "zero.cell.tables")) {
        lavobject@SampleStats@zero.cell.tables <- vector("list", ngroups)
      }
    } else {
      for (j in seq_along(lavobject@SampleStatsList)) {
        if (!.hasSlot(lavobject@SampleStatsList[[j]], "zero.cell.tables")) {
          lavobject@SampleStatsList[[j]]@zero.cell.tables <- vector("list", ngroups)
        }
      }
    }

    # 0.5-23 (21 Feb 2017)
    if (!.hasSlot(lavobject@Model, "nblocks")) {
      lavobject@Model@nblocks <- nblocks
    }
    if (is.lavaan.object) {
      if (!.hasSlot(lavobject@Data, "block.label")) {
        lavobject@Data@block.label <- character(0L) # or length nblocks ?
        lavobject@Data@level.label <- character(0L)
      }
    } else {
      for (j in seq_along(lavobject@DataList)) {
        if (!.hasSlot(lavobject@DataList[[j]], "block.label")) {
          lavobject@DataList[[j]]@block.label <- character(0L) # or length nblocks ?
          lavobject@DataList[[j]]@level.label <- character(0L)
        }
      }
    }

    # 0.5-23 (24 Feb 2017)
    if (is.lavaan.object) {
      if (!.hasSlot(lavobject@Data, "nlevels")) {
        lavobject@Data@nlevels <- 1L
        lavobject@Data@Lp <- vector("list", ngroups)
      }
    } else {
      for (j in seq_along(lavobject@DataList)) {
        if (!.hasSlot(lavobject@DataList[[j]], "nlevels")) {
          lavobject@DataList[[j]]@nlevels <- 1L
          lavobject@DataList[[j]]@Lp <- vector("list", ngroups)
        }
      }
    }
    if (is.lavaan.object) {
      if (!.hasSlot(lavobject@SampleStats, "YLp")) {
        lavobject@SampleStats@YLp <- vector("list", ngroups)
      }
    } else {
      for (j in seq_along(lavobject@SampleStatsList)) {
        if (!.hasSlot(lavobject@SampleStatsList[[j]], "YLp")) {
          lavobject@SampleStatsList[[j]]@YLp <- vector("list", ngroups)
        }
      }
    }

    # 0.6-1 (8 Mar 2017)
    if (is.lavaan.object) {
      if (!.hasSlot(lavobject, "h1")) {
        lavobject@h1 <- lav_h1_implied_logl(lavdata = lavobject@Data,
          lavsamplestats = lavobject@SampleStats,
          lavpartable = lavobject@ParTable,
          lavoptions = lavobject@Options)
        lavobject@baseline <- lav_lavaan_step15_baseline(
          lavoptions = lavobject@Options,
          lavsamplestats = lavobject@SampleStats,
          lavdata = lavobject@Data,
          lavcache = lavobject@Cache,
          lavpartable = lavobject@ParTable
        )
      }
      if (!.hasSlot(lavobject@Data, "ov.names.l")) {
        lavobject@Data@ov.names.l <- vector("list", 0L)
      }
    } else {
      for (j in seq_along(lavobject@DataList)) {
        if (!.hasSlot(lavobject@DataList[[j]], "ov.names.l")) {
          lavobject@DataList[[j]]@ov.names.l <- vector("list", 0L)
        }
      }
    }

    # 0.6-1 (10 Mar 2017)
    if (!.hasSlot(lavobject@Model, "multilevel")) {
      lavobject@Model@multilevel <- FALSE
    }

    # 0.6-1 (19 Mar 2017)
    if (is.lavaan.object) {
      if (!.hasSlot(lavobject, "loglik")) {
        lavobject@Model@ceq.simple.only <- FALSE
        lavobject@Model@cin.simple.only <- FALSE
        lavobject@loglik <- lav_model_loglik(
          lavdata = lavobject@Data,
          lavsamplestats = lavobject@SampleStats,
          lavimplied = lavobject@implied,
          lavmodel = lavobject@Model,
          lavoptions = lavobject@Options
        )
      }
    }

    # 0.6-1 (1 Oct 2017)
    if (is.lavaan.object) {
      if (!.hasSlot(lavobject@Data, "weights")) {
        lavobject@Data@weights <- vector("list", ngroups)
      }
    } else {
      for (j in seq_along(lavobject@DataList)) {
        if (!.hasSlot(lavobject@DataList[[j]], "weights")) {
          lavobject@DataList[[j]]@weights <- vector("list", ngroups)
        }
      }
    }

    # 0.6-1 (3 Oct 2017)
    if (is.lavaan.object) {
      if (!.hasSlot(lavobject@Data, "sampling.weights")) {
        lavobject@Data@sampling.weights <- character(0L)
      }
    } else {
      for (j in seq_along(lavobject@DataList)) {
        if (!.hasSlot(lavobject@DataList[[j]], "sampling.weights")) {
          lavobject@DataList[[j]]@sampling.weights <- character(0L)
        }
      }
    }

    # 0.6-1 (2 May 2018)
    if (is.lavaan.object) {
      if (!.hasSlot(lavobject, "version")) lavobject@version <- "PRE 0.6"
    }
  } # no-version-flag (pre 0.6)

  ### from here on, we assume that the object is generated by lavaan 0.6-1 or
  ### higher

  # check if @test list is named
  if (is.null(names(lavobject@test))) {
    names(lavobject@test) <- sapply(lavobject@test, "[[", "test")
  }

  # 0.6-2 (12 Jun 2018)
  if (!.hasSlot(lavobject@Model, "x.free.var.idx")) {
    lavobject@Model@x.free.var.idx <- integer(0L)
  }

  # 0.6-3 (17 Sep 2018)
  if (!is.lavaan.object) {
    if (!.hasSlot(lavobject, "h1List")) {
      lavobject@h1List <- vector("list", 0L)
      lavobject@loglikList <- vector("list", 0L)
    }
  }

  # 0.6-4 (30 Mar 2019)
  if (!.hasSlot(lavobject@Model, "ov.efa.idx")) {
    lavobject@Model@ov.efa.idx <- vector("list", nblocks)
    lavobject@Model@lv.efa.idx <- vector("list", nblocks)
  }

  # 0.6-4 (11 Apr 2019)
  if (!.hasSlot(lavobject@Model, "nefa")) {
    lavobject@Model@nefa <- 0L
  }

  # 0.6-4 (24 Apr 2019)
  if (!.hasSlot(lavobject@Model, "H")) {
    lavobject@Model@H <- vector("list", 0L)
    lavobject@Model@lv.order <- vector("list", 0L)
  }

  # 0.6-4 (26 Apr 2019)
  if (!.hasSlot(lavobject@Model, "ceq.efa.JAC")) {
    lavobject@Model@ceq.efa.JAC <- matrix(0, nrow <- 0L, ncol <- 0L)
  }

  # 0.6-5 (7 Jul 2019)
  if (!is.lavaan.object) {
    if (!.hasSlot(lavobject, "baselineList")) {
      lavobject@baselineList <- vector("list", 0L)
    }
  }

  # 0.6-8 (29 Sep 2020)
  if (!.hasSlot(lavobject@Model, "rv.ov")) {
    lavobject@Model@rv.ov <- vector("list", 0L)
    lavobject@Model@rv.lv <- vector("list", 0L)
  }

  # 0.6-8 (18 Dec 2020)
  if (!.hasSlot(lavobject@Model, "estimator.args")) {
    lavobject@Model@estimator.args <- vector("list", 0L)
  }

  # 0.6-9 (15 Mar 2021)
  if (!.hasSlot(lavobject@Model, "modprop")) {
    lavobject@Model@modprop = lav_model_properties(
      GLIST = lavobject@Model@GLIST,
      lavpartable = lavobject@ParTable,
      nmat = lavobject@Model@nmat,
      m.free.idx = lavobject@Model@m.free.idx
    )
  }

  # 0.6-9 (22 Jun 2021)
  if (is.lavaan.object) {
    if (!.hasSlot(lavobject, "internal")) {
      lavobject@internal <- vector("list", 0L)
    }
  } else {
    if (!.hasSlot(lavobject, "internalList")) {
      lavobject@internalList <- vector("list", 0L)
    }
  }

  # 0.6-11 (28 Feb 2022)
  if (!.hasSlot(lavobject@Model, "nx.unco")) {
    # is not available, unco == free
    lavobject@Model@nx.unco <- lavobject@Model@nx.free
    lavobject@Model@x.unco.idx <- lavobject@Model@x.free.idx
    lavobject@Model@ceq.simple.only <- FALSE
    lavobject@Model@ceq.simple.K <- matrix(0, nrow <- 0L, ncol <- 0L)
  }

  # 0.6-13 (25 Jul 2022)
  if (!.hasSlot(lavobject@Model, "correlation")) {
    lavobject@Model@correlation <- FALSE
  }

  # 0.6-18 (25 Apr 2024)
  if (!is.lavaan.object) {
    if (!.hasSlot(lavobject, "version")) {
      lavobject@version <- "PRE 0.6.18"
    }
  }

  # 0.6-19 (27 Sep 2024)
  if (!.hasSlot(lavobject@Model, "cin.simple.only")) {
    lavobject@Model@cin.simple.only <- FALSE
  }

  # 0.6-20 (24 Jan 2025)
  if (!.hasSlot(lavobject@Model, "composites")) {
    lavobject@Model@composites <- any(lavobject@ParTable$op == "<~")
  }

  # check missing options
  object_options <- lavobject@Options
  all_options <- lavOptions()
  missing.idx <- which(!names(all_options) %in% names(object_options))
  new_options <- c(object_options, all_options[missing.idx])
  # fill in some "default" values
  if (new_options$estimator.orig == "default") {
    new_options$estimator.orig <- new_options$estimator
  }
  new_options$gamma.vcov.mplus <- new_options$mimic == "Mplus"
  new_options$gamma.wls.mplus <- new_options$mimic == "Mplus"
  new_options$gls.v11.mplus <- new_options$mimic == "Mplus"
  new_options$cinformation.expected.mplus <- new_options$mimic == "Mplus"
  new_options$h1.information.meat <- "structured"
  new_options$mega.h1.information <- "unstructured"
  lavobject@Options <- new_options

  lavobject
}
