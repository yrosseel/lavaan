# data simulation engine behind simulateData() / lavSimulateData()
#
# This file contains the single, unified engine and its workers:
#  - lav_data_simulate()     : the public router (single-level vs multilevel)
#  - lav_data_simulate_sl()  : single-level worker (the historical engine);
#                              continuous, single-level output is byte-identical
#                              to the long-standing simulateData() behaviour, and
#                              it also handles ov.var / skewness / kurtosis /
#                              standardized / mass / categorical data
#  - lav_data_simulate_ml()  : multilevel worker; calls lavaan() directly to
#                              obtain the model-implied moments per block, and
#                              supports two levels (optionally multigroup)
#  - lav_data_simulate_old() : deprecated thin wrapper kept for backward
#                              compatibility
#  - helper functions (cluster.idx handling, Vale-Maurelli nonnormal data, ...)

# ---------------------------------------------------------------------------
# lav_data_simulate(): the single, unified data-simulation engine behind
# simulateData() / lavSimulateData().
#
# It routes to one of two internal workers, but presents a single interface:
#  - multilevel requests  -> lav_data_simulate_ml() (model-implied moments per
#    block, via a lavaan fit; supports two levels, optionally multigroup)
#  - single-level requests -> lav_data_simulate_sl() (the historical engine;
#    supports ov.var, skewness/kurtosis, standardized, mass, ...). The
#    continuous, single-level output is byte-identical to the historical
#    simulateData() output.
#
# The signature mirrors the historical simulateData() signature (so existing
# code keeps working), plus 'cluster.idx' (multilevel) and 'ordered.center'.
# ---------------------------------------------------------------------------
lav_data_simulate <- function(model = NULL,                       # nolint start
                              model.type = "sem",
                              # model modifiers
                              meanstructure = FALSE,
                              int.ov.free = TRUE,
                              int.lv.free = FALSE,
                              marker.int.zero = FALSE,
                              conditional.x = FALSE,
                              composites = TRUE,
                              fixed.x = FALSE,
                              orthogonal = FALSE,
                              std.lv = TRUE,
                              auto.fix.first = FALSE,
                              auto.fix.single = FALSE,
                              auto.var = TRUE,
                              auto.cov.lv.x = TRUE,
                              auto.cov.y = TRUE,
                              ...,
                              # data properties
                              sample.nobs = 500L,
                              ov.var = NULL,
                              group.label = NULL,
                              skewness = NULL,
                              kurtosis = NULL,
                              cluster.idx = NULL,
                              # control
                              seed = NULL,
                              empirical = FALSE,
                              mass = FALSE,
                              ordered.center = TRUE,
                              return.type = "data.frame",
                              return.fit = FALSE,
                              debug = FALSE,
                              standardized = FALSE) {           # nolint end
  multilevel <- !is.null(cluster.idx) || lav_simulate_is_multilevel(model)

  if (multilevel) {
    # ------- multilevel worker -------
    # the multilevel worker does not (yet) support these (single-level) features
    if (!is.null(ov.var) || !is.null(skewness) || !is.null(kurtosis) ||
        isTRUE(standardized) || isTRUE(mass)) {
      lav_msg_warn(gettext(
        "arguments 'ov.var', 'skewness', 'kurtosis', 'standardized' and 'mass'
        are not supported for multilevel data and will be ignored"))
    }
    return(lav_data_simulate_ml(
      model = model, cmd.pop = model.type,
      # forward the model modifiers to the population-model fit
      int.ov.free = int.ov.free, int.lv.free = int.lv.free,
      marker.int.zero = marker.int.zero, conditional.x = conditional.x,
      composites = composites, fixed.x = fixed.x, orthogonal = orthogonal,
      std.lv = std.lv, auto.fix.first = auto.fix.first,
      auto.fix.single = auto.fix.single, auto.var = auto.var,
      auto.cov.lv.x = auto.cov.lv.x, auto.cov.y = auto.cov.y,
      ...,
      sample.nobs = sample.nobs, cluster.idx = cluster.idx,
      seed = seed, empirical = empirical, ordered.center = ordered.center,
      return.fit = return.fit, output = return.type))
  }

  # ------- single-level worker -------
  lav_data_simulate_sl(
    model = model, model.type = model.type,
    meanstructure = meanstructure, int.ov.free = int.ov.free,
    int.lv.free = int.lv.free, marker.int.zero = marker.int.zero,
    conditional.x = conditional.x, composites = composites,
    fixed.x = fixed.x, orthogonal = orthogonal, std.lv = std.lv,
    auto.fix.first = auto.fix.first, auto.fix.single = auto.fix.single,
    auto.var = auto.var, auto.cov.lv.x = auto.cov.lv.x,
    auto.cov.y = auto.cov.y, ...,
    sample.nobs = sample.nobs, ov.var = ov.var, group.label = group.label,
    skewness = skewness, kurtosis = kurtosis, seed = seed,
    empirical = empirical, mass = mass, return.type = return.type,
    return.fit = return.fit, debug = debug, standardized = standardized,
    ordered.center = ordered.center)
}

# exported synonym
lavSimulateData <- lav_data_simulate

# detect whether a model (syntax string or parameter table) is multilevel
lav_simulate_is_multilevel <- function(model = NULL) {
  if (is.null(model)) {
    return(FALSE)
  }
  if (is.character(model)) {
    txt <- paste(model, collapse = "\n")
    return(grepl("(^|\n)[ \t]*level[ \t]*:", txt))
  }
  if (is.list(model)) {
    if (!is.null(model$level)) {
      lev <- model$level[!is.na(model$level)]
      return(length(unique(lev[lev > 0L])) > 1L)
    }
  }
  FALSE
}

# multilevel data-simulation worker behind lav_data_simulate()
lav_data_simulate_ml <- function(model = NULL,
                              cmd.pop = "sem",
                              ...,
                              # data properties
                              sample.nobs = 1000L,
                              cluster.idx = NULL,
                              # control
                              seed = NULL,
                              empirical = FALSE,
                              ordered.center = TRUE,
                              # output
                              add.labels = TRUE,
                              return.fit = FALSE,
                              output = "data.frame") {
  if (!is.null(seed)) set.seed(seed)

  # dotdotdot
  dotdotdot <- list(...)
  dotdotdot.orig <- dotdotdot

  # remove/override some options
  dotdotdot$verbose <- FALSE
  dotdotdot$debug <- FALSE
  dotdotdot$data <- NULL
  dotdotdot$sample.cov <- NULL
  dotdotdot$do.fit <- NULL

  # always use meanstructure = TRUE
  dotdotdot$meanstructure <- TRUE

  # helper to 'fit' the (no-data) population model; for multigroup models that
  # are only defined via c() modifiers, lavaan needs sample.nobs to know the
  # number of groups, so we retry with sample.nobs if the first attempt fails
  fit_pop_model <- function(dd) {
    out <- tryCatch(do.call(cmd.pop, args = c(list(model = model), dd)),
                    error = function(e) e)
    if (inherits(out, "error")) {
      out <- do.call(cmd.pop, args = c(list(model = model,
                                            sample.nobs = sample.nobs), dd))
    }
    out
  }

  # 'fit' population model: first pretend we generate continuous data only
  # (so the model-implied moments are the latent-response moments)
  dotdotdot.con <- dotdotdot
  dotdotdot.con$ordered <- NULL
  fit.con <- fit_pop_model(dotdotdot.con)

  # categorical? refit keeping the 'ordered' argument, to obtain thresholds
  if (!is.null(dotdotdot.orig$ordered) || fit.con@Model@categorical) {
    fit.pop <- fit_pop_model(dotdotdot)
  } else {
    fit.pop <- fit.con
  }

  # extract model implied statistics and data slot
  lavimplied  <- fit.con@implied # take continuous mean/cov
  lavdata     <- fit.pop@Data
  lavmodel    <- fit.pop@Model
  lavpartable <- fit.pop@ParTable
  lavoptions  <- fit.pop@Options

  # number of groups/levels/blocks
  ngroups <- lav_pt_ngroups(lavpartable)
  nblocks <- lav_pt_nblocks(lavpartable)
  nlevels <- lavdata@nlevels

  # check/expand the sample.nobs and cluster.idx arguments
  tmp <- lav_data_simulate_nobs(sample.nobs = sample.nobs,
                                cluster.idx = cluster.idx,
                                ngroups = ngroups, nblocks = nblocks,
                                nlevels = nlevels)
  sample.nobs <- tmp$sample.nobs
  cluster.idx <- tmp$cluster.idx

  # check if ov.names are the same for each group
  if (ngroups > 1L) {
    N1 <- lavdata@ov.names[[1]]
    if (!all(sapply(lavdata@ov.names, function(x) all(x %in% N1)))) {
      if (output == "data.frame") {
        output <- "matrix"
        lav_msg_warn(gettext(
          "groups do not contain the same set of variables; changing output=
          argument to \"matrix\""))
      }
    }
  }

  # prepare data containers
  X <- vector("list", length = nblocks)

  # generate data per BLOCK
  for (b in seq_len(nblocks)) {
    if (lavoptions$conditional.x) {
      lav_msg_stop(gettext("conditional.x = TRUE is not supported (yet) by the
                            multilevel data simulation"))
    } else {
      COV <- lavimplied$cov[[b]]
      MU  <- lavimplied$mean[[b]]
    }

    # if empirical = TRUE, rescale by N/(N-1), so that estimator = ML returns
    # exact results
    if (empirical) {
      if (sample.nobs[b] < NCOL(COV)) {
        lav_msg_stop(gettextf(
          "empirical = TRUE requires sample.nobs (= %1$s) to be larger than the
          number of variables (= %2$s) in block = %3$s",
          sample.nobs[b], NCOL(COV), b))
      }
      if (nlevels > 1L && (b %% nlevels == 1L)) {
        COV <- COV * sample.nobs[b] / (sample.nobs[b] - sample.nobs[b + 1])
      } else {
        COV <- COV * sample.nobs[b] / (sample.nobs[b] - 1)
      }
    }

    # generate normal data (using sign-invariant method for reproducibility)
    tmp <- try(lav_mvrnorm(
      n = sample.nobs[b],
      mu = MU, sigma_1 = COV, empirical = empirical),
      silent = TRUE)

    if (inherits(tmp, "try-error")) {
      # something went wrong; most likely: non-positive COV?
      ev <- eigen(COV, symmetric = TRUE, only.values = TRUE)$values
      if (any(ev < 0)) {
        lav_msg_stop(gettextf(
          "model-implied covariance matrix is not positive-definite in block =
          %1$s; smallest eigen value = %2$s; change the model parameters.",
          b, round(min(ev), 5)))
      } else {
        lav_msg_stop(gettextf("data generation failed for block = %s", b))
      }
    } else {
      X[[b]] <- unname(tmp)
    }
  } # block

  if (output == "block") {
    return(X)
  }

  # if multilevel, make a copy, and create X[[g]] per group
  if (nlevels > 1L) {
    X.block <- X
    X <- vector("list", length = ngroups)
  }

  # assemble data per group
  group.values <- lav_pt_group_values(lavpartable)
  for (g in 1:ngroups) {
    # multilevel?
    if (nlevels > 1L) {
      # which block?
      bb <- (g - 1) * nlevels + 1L

      Lp <- lavdata@Lp[[g]]
      p.tilde <- length(lavdata@ov.names[[g]])
      tmp1 <- matrix(0, nrow(X.block[[bb]]), p.tilde + 1L) # one extra column
      tmp2 <- matrix(0, nrow(X.block[[bb]]), p.tilde + 1L) # for the clus id

      # level 1
      tmp1[, Lp$ov.idx[[1]]] <- X.block[[bb]]

      # level 2 (expand cluster-level values to the level-1 units)
      tmp2[, Lp$ov.idx[[2]]] <- X.block[[bb + 1L]][cluster.idx[[g]], ,
                                                   drop = FALSE]
      # final
      X[[g]] <- tmp1 + tmp2

      # cluster id
      X[[g]][, p.tilde + 1L] <- cluster.idx[[g]]
    }

    # add variable names?
    if (add.labels) {
      if (nlevels > 1L) {
        colnames(X[[g]]) <- c(lavdata@ov.names[[g]], "cluster")
      } else {
        colnames(X[[g]]) <- lavdata@ov.names[[g]]
      }
    }

    # any categorical variables?
    ov.ord <- lav_object_vnames(fit.pop, "ov.ord", group = group.values[g])
    if (is.list(ov.ord)) {
      # multilevel -> use within level only
      ov.ord <- ov.ord[[1L]]
    }
    if (length(ov.ord) > 0L) {
      ov.names <- lavdata@ov.names[[g]]

      # which block?
      bb <- (g - 1) * nlevels + 1L

      # th/names
      TH.VAL <- as.numeric(fit.pop@implied$th[[bb]])
      if (length(lavmodel@num.idx[[bb]]) > 0L) {
        NUM.idx <- which(lavmodel@th.idx[[bb]] == 0)
        TH.VAL <- TH.VAL[-NUM.idx]
      }
      th.names <- fit.pop@pta$vnames$th[[bb]]
      TH.NAMES <- sapply(strsplit(th.names, split = "|", fixed = TRUE),
                         "[[", 1L)

      # use thresholds to cut
      for (o in ov.ord) {
        o.idx  <- which(o == ov.names)
        th.idx <- which(o == TH.NAMES)
        th.val <- c(-Inf, sort(TH.VAL[th.idx]), +Inf)
        tmp <- X[[g]][, o.idx]
        if (ordered.center) {
          # center first (so the cut also works when the model-implied 'mean'
          # is nonzero); set ordered.center = FALSE for the 'old' behaviour
          tmp <- tmp - mean(tmp, na.rm = TRUE)
        }
        X[[g]][, o.idx] <- cut(tmp, th.val, labels = FALSE)
      }
    }
  }

  # output
  if (output == "matrix") {
    if (ngroups == 1L) {
      out <- X[[1L]]
    } else {
      out <- X
    }
  } else if (output == "data.frame") {
    if (ngroups == 1L) {
      out <- as.data.frame(X[[1L]], stringsAsFactors = FALSE)
    } else {
      # rbind + add group column
      out <- as.data.frame(do.call("rbind", X), stringsAsFactors = FALSE)
      out$group <- rep.int(1:ngroups, times = sapply(X, NROW))
    }
  } else if (output == "cov") {
    if (ngroups == 1L) {
      out <- cov(X[[1L]])
    } else {
      out <- lapply(X, cov)
    }
  } else {
    lav_msg_stop(gettextf("unknown option for argument output: %s", output))
  }

  if (return.fit) {
    attr(out, "fit") <- fit.pop
  }

  out
}

# check/expand the sample.nobs and cluster.idx arguments; for the multilevel
# case, return a per-block sample.nobs vector and a per-group cluster.idx list
lav_data_simulate_nobs <- function(sample.nobs = 1000L, cluster.idx = NULL,
                                   ngroups = 1L, nblocks = 1L, nlevels = 1L) {
  if (nlevels > 1L) {
    if (is.null(cluster.idx)) {
      # no cluster.idx given: build balanced clusters from sample.nobs
      # sample.nobs is interpreted (per group) as:
      #  - a single number: total number of level-1 units (clusters default
      #    to 1/10 of that, with a minimum of 2)
      #  - a vector c(n.level1, n.level2): #level-1 units and #clusters
      if (length(sample.nobs) == 1L) {
        n1 <- as.integer(sample.nobs)
        n2 <- max(2L, as.integer(round(n1 / 10)))
        nobs.l <- c(n1, n2)
      } else if (length(sample.nobs) == nlevels) {
        nobs.l <- as.integer(sample.nobs)
      } else {
        lav_msg_stop(gettext("for the multilevel case, sample.nobs should be a
                             single number, or a vector with one number per
                             level (eg c(1000, 100))"))
      }
      cluster.idx <- vector("list", ngroups)
      for (g in seq_len(ngroups)) {
        cluster.idx[[g]] <- lav_data_simulate_clusidx(nobs.l[1L], nobs.l[2L])
      }
    } else {
      # cluster.idx given
      if (!is.list(cluster.idx)) {
        cluster.idx <- rep(list(cluster.idx), ngroups)
      }
    }
    # cluster.idx values within each group should be 1, 2, ... nclus
    sample.nobs <- numeric(nblocks)
    for (g in seq_len(ngroups)) {
      cluster.idx[[g]] <- as.integer(as.factor(cluster.idx[[g]]))
      gg <- (g - 1L) * nlevels + 1L
      sample.nobs[gg]     <- length(cluster.idx[[g]])
      sample.nobs[gg + 1] <- length(unique(cluster.idx[[g]]))
    }
  } else {
    # single level
    if (length(sample.nobs) == ngroups) {
      # nothing to do
    } else if (ngroups > 1L && length(sample.nobs) == 1L) {
      sample.nobs <- rep.int(sample.nobs, ngroups)
    } else {
      lav_msg_stop(gettextf(
        "ngroups = %1$s but sample.nobs has length = %2$s",
        ngroups, length(sample.nobs)))
    }
  }

  list(sample.nobs = as.integer(sample.nobs), cluster.idx = cluster.idx)
}

# create a balanced cluster.idx vector with n1 level-1 units in n2 clusters
lav_data_simulate_clusidx <- function(n1, n2) {
  if (n2 < 1L) n2 <- 1L
  if (n2 > n1) n2 <- n1
  base <- n1 %/% n2
  rem  <- n1 %% n2
  sizes <- rep.int(base, n2)
  if (rem > 0L) sizes[seq_len(rem)] <- sizes[seq_len(rem)] + 1L
  rep.int(seq_len(n2), times = sizes)
}


# single-level data-simulation worker (the historical engine)
#
# initial version: YR 24 jan 2011
# revision for 0.4-11: YR 21 okt 2011
#
# This is the single-level worker behind lav_data_simulate(). It is kept as a
# separate function because the continuous, single-level output must remain
# byte-identical to the historical simulateData() output. The multilevel case is
# handled by lav_data_simulate_ml().
#
lav_data_simulate_sl <- function( # user-specified model    # nolint start
                         model = NULL,
                         model.type = "sem",
                         # model modifiers
                         meanstructure = FALSE,
                         int.ov.free = TRUE,
                         int.lv.free = FALSE,
                         marker.int.zero = FALSE,
                         conditional.x = FALSE,
                         composites = TRUE,
                         fixed.x = FALSE,
                         orthogonal = FALSE,
                         std.lv = TRUE,
                         auto.fix.first = FALSE,
                         auto.fix.single = FALSE,
                         auto.var = TRUE,
                         auto.cov.lv.x = TRUE,
                         auto.cov.y = TRUE,
                         ...,
                         # data properties
                         sample.nobs = 500L,
                         ov.var = NULL,
                         group.label = paste("G", 1:ngroups, sep = ""),
                         skewness = NULL,
                         kurtosis = NULL,
                         # control
                         seed = NULL,
                         empirical = FALSE,
                         mass = FALSE,
                         return.type = "data.frame",
                         return.fit = FALSE,
                         debug = FALSE,
                         standardized = FALSE,
                         ordered.center = FALSE) {         # nolint end
  if (!missing(debug)) {
    current_debug <- lav_debug()
    if (lav_debug(debug))
      on.exit(lav_debug(current_debug), TRUE)
  }
  if (!is.null(seed)) set.seed(seed)
  # if(!exists(".Random.seed", envir = .GlobalEnv))
  #    runif(1)               # initialize the RNG if necessary
  # RNGstate <- .Random.seed

  # lav_model_pt
  if (is.list(model)) {
    # two possibilities: either model is already lavaanified
    # or it is something else...
    if (!is.null(model$lhs) && !is.null(model$op) &&
      !is.null(model$rhs) && !is.null(model$free)) {
      lav <- model

      # until 0.6-5, we only used the 'ustart' column
      # but what if 'lav' is a fitted lavaan object -> use 'est'
      if (!is.null(lav$est)) {
        lav$ustart <- lav$est
        lav$se <- NULL
        lav$est <- NULL
        lav$start <- NULL
      }
    } else if (is.character(model[[1]])) {
      lav_msg_stop(gettext("model is a list, but not a parameterTable?"))
    }
  } else {
    lav <- lav_model_pt(
      model = model,
      meanstructure = meanstructure,
      int_ov_free = int.ov.free,
      int_lv_free = int.lv.free,
      marker_int_zero = marker.int.zero,
      composites = composites,
      conditional_x = conditional.x,
      fixed_x = fixed.x,
      orthogonal = orthogonal,
      std_lv = std.lv,
      auto_fix_first = auto.fix.first,
      auto_fix_single = auto.fix.single,
      auto_var = auto.var,
      auto_cov_lv_x = auto.cov.lv.x,
      auto_cov_y = auto.cov.y,
      ngroups = length(sample.nobs)
    )
  }

  group_values <- lav_pt_group_values(lav)
  if (lav_debug()) {
    cat("initial lav\n")
    print(as.data.frame(lav))
  }

  # fill in any remaining NA values (needed for unstandardize)
  # 1 for variances and (unstandardized) factor loadings, 0 otherwise
  idx <- which(lav$op == "=~" & is.na(lav$ustart))
  if (length(idx) > 0L) {
    if (standardized) {
      lav$ustart[idx] <- 0.7
    } else {
      lav$ustart[idx] <- 1.0
    }
  }

  idx <- which(lav$op == "~~" & is.na(lav$ustart) & lav$lhs == lav$rhs)
  if (length(idx) > 0L) lav$ustart[idx] <- 1.0

  idx <- which(lav$op == "~" & is.na(lav$ustart))
  if (length(idx) > 0L) {
    lav_msg_warn(gettext(
      "some regression coefficients are unspecified and will be set to zero"))
  }

  idx <- which(is.na(lav$ustart))
  if (length(idx) > 0L) lav$ustart[idx] <- 0.0

  if (lav_debug()) {
    cat("lav + default values\n")
    print(as.data.frame(lav))
  }

  # set residual variances to enforce a standardized solution
  # but only if no *residual* variances have been specified in the syntax

  if (standardized) {
    # check if factor loadings are smaller than 1.0
    lambda_idx <- which(lav$op == "=~")
    if (any(lav$ustart[lambda_idx] >= 1.0)) {
      lav_msg_warn(gettext("standardized=TRUE but factor loadings are >= 1.0"))
    }

    # check if regression coefficients are smaller than 1.0
    reg_idx <- which(lav$op == "~")
    if (any(lav$ustart[reg_idx] >= 1.0)) {
      lav_msg_warn(gettext(
        "standardized=TRUE but regression coefficients are >= 1.0"))
    }

    # for ordered observed variables, we will get '0.0', but that is ok
    # so there is no need to make a distinction between numeric/ordered
    # here??
    ngroups <- lav_pt_ngroups(lav)
    ov_names <- lav_pt_vnames(lav, "ov")
    ov_nox <- lav_pt_vnames(lav, "ov.nox")
    # lv_names <- lav_pt_vnames(lav, "lv")
    # lv_y <- lav_pt_vnames(lav, "lv.y")
    lv_nox <- lav_pt_vnames(lav, "lv.nox")
    ov_var_idx <- which(lav$op == "~~" & lav$lhs %in% ov_nox &
      lav$rhs == lav$lhs)
    lv_var_idx <- which(lav$op == "~~" & lav$lhs %in% lv_nox &
      lav$rhs == lav$lhs)
    if (any(lav$user[c(ov_var_idx, lv_var_idx)] > 0L)) {
      lav_msg_warn(gettext(
        "if residual variances are specified, please use standardized=FALSE"))
    }

    # new in 0.6-20: - use lav_lisrel_residual_variances
    #                - use lav_lisrel_comp_set_intresvar
    dotdotdot <- list(...)
    dotdotdot$sample.nobs <- sample.nobs
    dotdotdot$fixed.x <- FALSE # for now
    dotdotdot$representation <- "LISREL"
    dotdotdot$composites <- composites
    dotdotdot$correlation <- TRUE # this is the trick
    tmp_fit <- do.call("lavaan", args = c(list(model = lav), dotdotdot))
    # set/get parameters to invoke lav_lisrel_residual_variances
    tmp_lav <- tmp_fit@ParTable
    tmp_x <- lav_model_get_parameters(tmp_fit@Model)
    tmp_model <- lav_model_set_parameters(tmp_fit@Model, x = tmp_x)
    tmp_lav$ustart <- lav_model_get_parameters(tmp_model, type = "user")

    # copy residual values to lav (without assuming parameter tables look the
    # same)
    res_idx <- c(ov_var_idx, lv_var_idx)
    for (i in seq_along(res_idx)) {
      # lookup this parameter in tmp.lav
      idx_in_lav <- res_idx[i]
      this_lhs <- lav$lhs[idx_in_lav]
      idx_in_tmp <- which(tmp_lav$op == "~~" & tmp_lav$lhs == this_lhs &
                          tmp_lav$rhs == tmp_lav$lhs)
      if (length(idx_in_tmp) == 0L) {
        # hm, not found? Give a warning?
      } else {
        vals <- tmp_lav$ustart[idx_in_tmp]
        # check if we have unfortunate values
        bad_idx <- which(!is.finite(vals) | vals < 0)
        vals[bad_idx] <- 1.0 # not pretty, but safe
        lav$ustart[idx_in_lav] <- vals
      }
    }

    # this is what we did <0.6-20
    # fit <- lavaan(model = lav, sample.nobs = sample.nobs, ...)
    # Sigma.hat <- lav_model_sigma(lavmodel = fit@Model)
    # ETA <- lav_model_veta(lavmodel = fit@Model)

    # if (lav_debug()) {
    #   cat("Sigma.hat:\n")
    #   print(Sigma.hat)
    #   cat("Eta:\n")
    #   print(ETA)
    # }

    # # stage 1: standardize LV
    # if (length(lv.nox) > 0L) {
    #   for (g in 1:ngroups) {
    #     var.group <- which(lav$op == "~~" & lav$lhs %in% lv.nox &
    #       lav$rhs == lav$lhs &
    #       lav$group == group.values[g])
    #     eta.idx <- match(lv.nox, lv.names)
    #     lav$ustart[var.group] <- 1 - diag(ETA[[g]])[eta.idx]
    #   }
    # }
    # # refit
    # fit <- lavaan(model = lav, sample.nobs = sample.nobs, ...)
    # Sigma.hat <- lav_model_sigma(lavmodel = fit@Model)

    # if (lav_debug()) {
    #   cat("after stage 1:\n")
    #   cat("Sigma.hat:\n")
    #   print(Sigma.hat)
    # }

    # # stage 2: standardize OV
    # for (g in 1:ngroups) {
    #   var.group <- which(lav$op == "~~" & lav$lhs %in% ov.nox &
    #     lav$rhs == lav$lhs &
    #     lav$group == group.values[g])
    #   ov.idx <- match(ov.nox, ov.names)
    #   lav$ustart[var.group] <- 1 - diag(Sigma.hat[[g]])[ov.idx]
    # }

    # if (lav_debug()) {
    #   cat("after standardisation lav\n")
    #   print(as.data.frame(lav))
    # }
  }


  # unstandardize
  if (!is.null(ov.var)) {
    # FIXME: if ov.var is named, check the order of the elements

    # 1. unstandardize observed variables
    lav$ustart <- lav_unstandardize_ov(partable = lav, ov_var = ov.var)

    # 2. unstandardized latent variables

    if (lav_debug()) {
      cat("after unstandardisation lav\n")
      print(as.data.frame(lav))
    }
  }

  # fit the model without data
  fit <- lavaan(model = lav, sample.nobs = sample.nobs, ...)

  # the model-implied moments for the population
  sigma_hat <- lav_model_sigma(lavmodel = fit@Model)
  mu_hat <- lav_model_mu(lavmodel = fit@Model)
  if (fit@Model@categorical) {
    th <- lav_model_th(lavmodel = fit@Model)
  }

  if (lav_debug()) {
    cat("\nModel-implied moments (before Vale-Maurelli):\n")
    print(sigma_hat)
    print(mu_hat)
    if (exists("TH")) print(th)
  }

  # ngroups
  ngroups <- length(sample.nobs)

  # prepare
  x <- vector("list", length = ngroups)
  # out <- vector("list", length = ngroups)

  for (g in 1:ngroups) {
    cov_1 <- sigma_hat[[g]]

    # if empirical = TRUE, rescale by N/(N-1), so that estimator=ML
    # returns exact results
    if (empirical) {
      cov_1 <- cov_1 * sample.nobs[g] / (sample.nobs[g] - 1)
    }

    # Using sign-invariant method for cross-machine reproducibility
    if (is.null(skewness) && is.null(kurtosis)) {
      if (mass) {
        x[[g]] <- MASS::mvrnorm(
          n = sample.nobs[g],
          mu = mu_hat[[g]],
          Sigma = cov_1,
          empirical = empirical
        )
      } else {
        x[[g]] <- lav_mvrnorm(
          n = sample.nobs[g],
          mu = mu_hat[[g]],
          sigma_1 = cov_1,
          empirical = empirical
        )
      }
    } else {
      # first generate Z
      z <- lav_data_valemaurelli1983(
        n = sample.nobs[g],
        cor_1 = cov2cor(cov_1),
        skewness = skewness, # FIXME: per group?
        kurtosis = kurtosis, mass = mass
      )
      # rescale
      # Note: 'scale()' will first center, and then scale
      # but we need to first scale, and then center...
      # this was reported by Jordan Brace (9 may 2014)
      # X[[g]] <- scale(Z, center = -Mu.hat[[g]],
      #                   scale  = 1/sqrt(diag(COV)))

      # first, we scale
      tmp <- scale(z,
        center = FALSE,
        scale = 1 / sqrt(diag(cov_1))
      )[, , drop = FALSE]

      # then, we center
      x[[g]] <- sweep(tmp, MARGIN = 2, STATS = mu_hat[[g]], FUN = "+")
    }

    # any categorical variables?
    ov_ord <- lav_pt_vnames(lav, type = "ov.ord", group = group_values[g])
    if (length(ov_ord) > 0L) {
      ov_names <- lav_pt_vnames(lav, type = "ov", group = group_values[g])
      # use thresholds to cut
      for (o in ov_ord) {
        o_idx <- which(o == ov_names)
        th_idx <- which(lav$op == "|" & lav$lhs == o &
          lav$group == group_values[g])
        th_val <- c(-Inf, sort(lav$ustart[th_idx]), +Inf)
        tmp_o <- x[[g]][, o_idx]
        if (ordered.center) {
          # center first (so the cut also works when the model-implied 'mean'
          # is nonzero); ordered.center = FALSE reproduces the historical cut
          tmp_o <- tmp_o - mean(tmp_o, na.rm = TRUE)
        }
        x[[g]][, o_idx] <- as.integer(cut(tmp_o, th_val))
      }
    }

    if (return.type == "data.frame") x[[g]] <- as.data.frame(x[[g]])
  }

  if (return.type == "matrix") {
    if (ngroups == 1L) {
      x[[1L]]
    } else {
      x
    }
  } else if (return.type == "data.frame") {
    data_1 <- x[[1L]]

    # if multiple groups, add group column
    if (ngroups > 1L) {
      for (g in 2:ngroups) {
        data_1 <- rbind(data_1, x[[g]])
      }
      data_1$group <- rep(1:ngroups, times = sample.nobs)
    }
    var_names <- lav_pt_vnames(fit@ParTable, type = "ov", group = 1L)
    if (ngroups > 1L) var_names <- c(var_names, "group")
    names(data_1) <- var_names
    if (return.fit) {
      attr(data_1, "fit") <- fit
    }
    data_1
  } else if (return.type == "cov") {
    if (ngroups == 1L) {
      cov(x[[1L]])
    } else {
      cov_list <- lapply(x, cov)
      cov_list
    }
  }
}
# deprecated: kept for backward compatibility (eg as the default data_function
# of lavSimulate()); routes through the unified engine lav_data_simulate(), which
# dispatches to lav_data_simulate_sl() (single level) or lav_data_simulate_ml()
# (multilevel). New code should call lavSimulateData() / lav_data_simulate().
lav_data_simulate_old <- function(..., ordered.center = FALSE) {
  lav_data_simulate(..., ordered.center = ordered.center)
}



lav_data_valemaurelli1983 <- function(n = 100L, cor_1, skewness, kurtosis,
                                      mass = FALSE) {
  fleishman1978_abcd <- function(skewness, kurtosis) {
    system_function <- function(x, skewness, kurtosis) {
      b <- x[1L]
      c_1 <- x[2L]
      d <- x[3L]
      eq1 <- b * b + 6 * b * d + 2 * c_1 * c_1 + 15 * d * d - 1
      eq2 <- 2 * c_1 * (b * b + 24 * b * d + 105 * d * d + 2) - skewness
      eq3 <- 24 * (b * d + c_1 * c_1 * (1 + b * b + 28 * b * d) +
        d * d * (12 + 48 * b * d + 141 * c_1 * c_1 + 225 * d * d)) - kurtosis
      eq <- c(eq1, eq2, eq3)
      sum(eq * eq) ## SS
    }

    out <- nlminb(
      start = c(1, 0, 0), objective = system_function,
      scale = 10,
      control = list(trace = 0),
      skewness = skewness, kurtosis = kurtosis
    )
    if (out$convergence != 0 || out$objective > 1e-5) {
      lav_msg_warn(gettext("lav_data_valemaurelli1983 method did not converge,
                   or it did not find the roots"))
    }
    b <- out$par[1L]
    c_1 <- out$par[2L]
    d <- out$par[3L]
    a <- -c_1
    c(a, b, c_1, d)
  }

  get_icov <- function(b1, c1, d1, b2, c2, d2, r) {
    objective_function <- function(x, b1, c1, d1, b2, c2, d2, r) {
      rho <- x[1L]
      eq <- rho * (b1 * b2 + 3 * b1 * d2 + 3 * d1 * b2 + 9 * d1 * d2) +
        rho * rho * (2 * c1 * c2) + rho * rho * rho * (6 * d1 * d2) - r
      eq * eq
    }

    # gradientFunction <- function(x, bcd1, bcd2, R) {
    #
    # }

    out <- nlminb(
      start = r, objective = objective_function,
      scale = 10, control = list(trace = 0),
      b1 = b1, c1 = c1, d1 = d1, b2 = b2, c2 = c2, d2 = d2, r = r
    )
    if (out$convergence != 0 || out$objective > 1e-5)
      lav_msg_warn(gettext("no convergence"))
    rho <- out$par[1L]
    rho
  }

  # number of variables
  nvar <- ncol(cor_1)
  # check skewness
  if (is.null(skewness)) {
    sk <- rep(0, nvar)
  } else if (length(skewness) == nvar) {
    sk <- skewness
  } else if (length(skewness) == 1L) {
    sk <- rep(skewness, nvar)
  } else {
    lav_msg_stop(gettext("skewness has wrong length"))
  }

  if (is.null(kurtosis)) {
    ku <- rep(0, nvar)
  } else if (length(kurtosis) == nvar) {
    ku <- kurtosis
  } else if (length(kurtosis) == 1L) {
    ku <- rep(kurtosis, nvar)
  } else {
    lav_msg_stop(gettext("kurtosis has wrong length"))
  }

  # create Fleishman table
  ftable_1 <- matrix(0, nvar, 4L)
  for (i in 1:nvar) {
    ftable_1[i, ] <- fleishman1978_abcd(skewness = sk[i], kurtosis = ku[i])
  }

  # compute intermediate correlations between all pairs
  icor <- diag(nvar)
  for (j in 1:(nvar - 1L)) {
    for (i in (j + 1):nvar) {
      if (cor_1[i, j] == 0) next
      icor[i, j] <- icor[j, i] <-
        get_icov(ftable_1[i, 2], ftable_1[i, 3], ftable_1[i, 4],
          ftable_1[j, 2], ftable_1[j, 3], ftable_1[j, 4],
          r = cor_1[i, j]
        )
    }
  }

  if (lav_debug()) {
    cat("\nOriginal correlations (for Vale-Maurelli):\n")
    print(cor_1)
    cat("\nIntermediate correlations (for Vale-Maurelli):\n")
    print(icor)
    cat("\nEigen values ICOR:\n")
    print(eigen(icor)$values)
  }

  # generate Z (using sign-invariant method for cross-machine reproducibility)
  if (mass) {
    x <- z <- MASS::mvrnorm(n = n, mu = rep(0, nvar), Sigma = icor)
  } else {
    x <- z <- lav_mvrnorm(n = n, mu = rep(0, nvar), sigma_1 = icor)
  }

  # transform Z using Fleishman constants
  for (i in 1:nvar) {
    x[, i] <- ftable_1[i, 1L] + ftable_1[i, 2L] * z[, i] +
      ftable_1[i, 3L] * z[, i] * z[, i] +
      ftable_1[i, 4L] * z[, i] * z[, i] * z[, i]
  }

  x
}
