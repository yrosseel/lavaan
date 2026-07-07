# lavEffects(): compute (and provide standard errors for) various 'effects'
# that are functions of the estimated model parameters.
#
# In this first version, we focus on the classic LISREL-style 'total' and
# 'indirect' effects among the variables that appear in the structural part
# of the model (the BETA matrix in the LISREL representation). Both observed
# and latent variables that are involved in a regression are included.
#
# All effects are derived from the reduced-form matrix (I - BETA)^{-1}:
#
#   - direct   effect of j on i:   BETA[i, j]
#   - total    effect of j on i:   ((I - BETA)^{-1} - I)[i, j]
#   - indirect effect of j on i:   total[i, j] - direct[i, j]
#
# Standard errors are obtained either by the delta method (se_def = "delta")
# or by the Monte Carlo method (se_def = "monte.carlo", the default), both of
# which only need the point estimates and the (co)variance matrix of the free
# parameters.
#
# Future versions may add standardized effects, average treatment effects, and
# support for conditional.x = TRUE models (effects of exogenous covariates that
# are stored in the GAMMA matrix).

# YR 2026 -- first version

lavEffects <- function(object,
                       effects = c("total", "indirect"),
                       se_def = NULL,
                       level = 0.95,
                       monte_carlo = NULL,
                       boot_ci_type = "perc",
                       zstat = TRUE,
                       pvalue = TRUE,
                       ci = TRUE,
                       standardized = FALSE,
                       cov_std = TRUE,
                       add_class = TRUE,
                       output = "data.frame",
                      ...) {
  dotdotdot <- list(...)
  lav_adapt_func(environment(), dotdotdot, NULL)

  # check object
  stopifnot(inherits(object, "lavaan"))

  # are bootstrap draws available?
  boot.available <- !is.null(object@boot$coef) &&
    NROW(object@boot$coef) > 0L

  # resolve the 'se_def' argument: if the user did not specify it, honor the
  # choices that were made when the model was fitted
  if (is.null(se_def)) {
    if (identical(object@Options$se, "bootstrap") && boot.available) {
      # se = "bootstrap" -> use the available bootstrap draws
      se_def <- "bootstrap"
    } else if (identical(object@Options$se.def, "monte.carlo")) {
      # the user explicitly requested the Monte Carlo method for the
      # defined parameters of the original model
      se_def <- "monte.carlo"
    } else {
      # lavEffects() default
      se_def <- "monte.carlo"
    }
  }

  # check/expand 'effects' argument
  effects.choices <- c("total", "indirect", "direct")
  if (!is.character(effects) || length(effects) == 0L) {
    lav_msg_stop(gettext(
      "The \"effects\" argument must be a character vector with one or more of:
       \"total\", \"indirect\", \"direct\"."))
  }
  bad.idx <- which(!effects %in% effects.choices)
  if (length(bad.idx) > 0L) {
    lav_msg_stop(gettextf(
      "Unknown effect(s) in the \"effects\" argument: %s; valid choices are:
       \"total\", \"indirect\", \"direct\".",
      paste(dQuote(effects[bad.idx]), collapse = ", ")))
  }
  # keep canonical order, drop duplicates
  effects <- effects.choices[effects.choices %in% effects]

  # check 'se_def' argument
  se_def <- tolower(as.character(se_def)[1L])
  # allow a few aliases
  se_def <- switch(se_def,
    "mc" = , "montecarlo" = , "monte_carlo" =, "monte-carlo" = "monte.carlo",
    "boot" = "bootstrap",
    se_def)
  if (!se_def %in% c("monte.carlo", "delta", "bootstrap", "none")) {
    lav_msg_stop(gettextf(
      "Unknown \"se_def\" argument: %s; valid choices are \"monte.carlo\",
       \"delta\", \"bootstrap\" or \"none\".", dQuote(se_def)))
  }
  if (identical(se_def, "bootstrap") && !boot.available) {
    lav_msg_warn(gettext(
      "se_def = \"bootstrap\" was requested, but no bootstrap draws are
       available in the fitted object; standard errors are not available.
       Refit the model with se = \"bootstrap\" (or use bootstrapLavaan())."))
    se_def <- "none"
  }

  # check 'boot_ci_type' argument (only relevant when se_def = "bootstrap")
  boot_ci_type <- as.character(boot_ci_type)[1L]
  if (!boot_ci_type %in% c("norm", "basic", "perc", "bca.simple", "bca")) {
    lav_msg_stop(gettextf(
      "Unknown \"boot_ci_type\" argument: %s; valid choices are \"norm\",
       \"basic\", \"perc\", \"bca.simple\" or \"bca\".", dQuote(boot_ci_type)))
  }

  # resolve the Monte Carlo settings (R, seed): if not specified, inherit the
  # 'monte.carlo' list that was used when the model was fitted
  if (is.null(monte_carlo)) {
    monte_carlo <- object@Options$monte.carlo
  }
  mc_r <- 20000L
  mc_seed <- NULL
  if (!is.null(monte_carlo)) {
    if (!is.null(monte_carlo$R)) {
      mc_r <- as.integer(monte_carlo$R)
    }
    if (!is.null(monte_carlo$seed)) {
      mc_seed <- monte_carlo$seed
    }
  }

  # check 'output' argument
  output <- tolower(as.character(output)[1L])
  if (!output %in% c("data.frame", "list")) {
    lav_msg_stop(gettextf(
      "Unknown \"output\" argument: %s; valid choices are \"data.frame\"
       or \"list\".", dQuote(output)))
  }

  lavmodel <- object@Model
  lavdata  <- object@Data

  # only the LISREL representation is supported (the default)
  if (lavmodel@representation != "LISREL") {
    lav_msg_stop(gettextf(
      "lavEffects() is only available for the LISREL representation (not %s).",
      dQuote(lavmodel@representation)))
  }

  nblocks <- lavmodel@nblocks
  nmat    <- lavmodel@nmat
  glist_names <- names(lavmodel@GLIST)

  # build the per-block 'specs' that describe, for each block, which effects
  # are structurally non-zero and must be reported.
  #
  # The structural model is  eta = B eta + Gamma x + zeta, so we work with an
  # augmented effect matrix whose rows are the (endogenous) BETA variables
  # 'eta' and whose columns are the predictors: the eta variables themselves
  # plus, when conditional.x = TRUE, the exogenous covariates 'x' (stored in
  # the GAMMA matrix). When conditional.x = FALSE there is no GAMMA matrix and
  # the columns are simply the eta variables (square case).
  specs <- vector("list", nblocks)
  has_beta <- logical(nblocks)
  for (b in seq_len(nblocks)) {
    # which GLIST matrices belong to block b?
    mm_in_block <- seq_len(nmat[b]) + cumsum(c(0L, nmat))[b]
    beta_idx <- mm_in_block[glist_names[mm_in_block] == "beta"]
    if (length(beta_idx) == 0L) {
      has_beta[b] <- FALSE
      next
    }
    has_beta[b] <- TRUE
    beta_idx <- beta_idx[1L]

    b0 <- lavmodel@GLIST[[beta_idx]]
    nr <- nrow(b0)
    row.names.b <- lavmodel@dimNames[[beta_idx]][[1L]]
    if (is.null(row.names.b)) {
      row.names.b <- as.character(seq_len(nr))
    }
    beta_m_idx <- lavmodel@m.free.idx[[beta_idx]]
    beta_x_idx <- lavmodel@x.free.idx[[beta_idx]]

    # structural adjacency among eta: a[i, j] = TRUE if there is a (possible)
    # direct edge j -> i (the element is free, or fixed to a non-zero value)
    a <- (b0 != 0)
    if (length(beta_m_idx) > 0L) {
      a[beta_m_idx] <- TRUE
    }
    diag(a) <- FALSE

    # reachability among eta (transitive closure): t1[i, j] = TRUE if j
    # influences i through one or more directed paths
    t1 <- lav_effects_reachable(a)

    # GAMMA (exogenous covariates), if conditional.x = TRUE
    gamma.idx <- mm_in_block[glist_names[mm_in_block] == "gamma"]
    g0 <- NULL
    gamma_m_idx <- gamma_x_idx <- integer(0L)
    x.names.b <- character(0L)
    nx <- 0L
    if (length(gamma.idx) > 0L) {
      gamma.idx <- gamma.idx[1L]
      g0 <- lavmodel@GLIST[[gamma.idx]]
      nx <- ncol(g0)
      x.names.b <- lavmodel@dimNames[[gamma.idx]][[2L]]
      if (is.null(x.names.b)) {
        x.names.b <- as.character(seq_len(nx))
      }
      gamma_m_idx <- lavmodel@m.free.idx[[gamma.idx]]
      gamma_x_idx <- lavmodel@x.free.idx[[gamma.idx]]
    }

    nr.col <- nr + nx
    col.names.b <- c(row.names.b, x.names.b)

    # augmented structural patterns (rows = eta, cols = eta then x)
    #   direct  : [ a | ax ]
    #   total   : [ t1 | ((I|t1) %*% ax) ]   (reflexive reach * Gamma)
    #   indirect: [ (a %*% t1) | (t1 %*% ax) ]
    ind_ee <- ((a %*% t1) > 0L)
    diag(ind_ee) <- FALSE
    if (nx > 0L) {
      ax <- (g0 != 0)
      if (length(gamma_m_idx) > 0L) {
        ax[gamma_m_idx] <- TRUE
      }
      refl <- t1
      diag(refl) <- TRUE # reflexive reachability (I | t1)
      tot <- cbind(t1, (refl %*% ax) > 0L)
      dir <- cbind(a,  ax)
      ind <- cbind(ind_ee, (t1 %*% ax) > 0L)
    } else {
      tot <- t1
      dir <- a
      ind <- ind_ee
    }

    # per effect type, the (i, j) cells (and linear indices) to report
    pat <- list(total = tot, indirect = ind, direct = dir)
    cells <- list()
    for (e in effects) {
      P <- pat[[e]]
      ij <- which(P, arr.ind = TRUE) # col 1 = i (row/outcome), 2 = j (col/pred)
      if (nrow(ij) == 0L) {
        cells[[e]] <- NULL
        next
      }
      # order by outcome (i), then predictor (j), following variable order
      o <- order(ij[, 1L], ij[, 2L])
      ij <- ij[o, , drop = FALSE]
      cells[[e]] <- list(
        i   = ij[, 1L],
        j   = ij[, 2L],
        lin = ij[, 1L] + (ij[, 2L] - 1L) * nr # linear idx into nr x nr.col mat
      )
    }

    specs[[b]] <- list(
      nr          = nr,         # number of rows (eta / outcomes)
      nr.col      = nr.col,     # number of columns (eta + x / predictors)
      row.names   = row.names.b,
      col.names   = col.names.b,
      b0          = b0,
      beta_m_idx  = beta_m_idx,
      beta_x_idx  = beta_x_idx,
      g0          = g0,
      gamma_m_idx = gamma_m_idx,
      gamma_x_idx = gamma_x_idx,
      cells       = cells
    )
  }

  if (!any(has_beta)) {
    lav_msg_warn(gettext(
      "The model contains no regressions (empty BETA matrix); there are no
       indirect or total effects to compute."))
    empty <- lav_effects_empty_df(effects, zstat, pvalue)
    return(empty)
  }

  # build the (ordered) skeleton data.frame and remember, for each row, the
  # block and the linear cell index, so that the effect-evaluation function
  # returns its values in exactly the same order
  skeleton <- lav_effects_skeleton(specs, effects, lavdata, lavmodel)
  if (nrow(skeleton$df) == 0L) {
    empty <- lav_effects_empty_df(effects, zstat, pvalue)
    return(empty)
  }
  eff_df <- skeleton$df
  row_block <- skeleton$block
  row_eff   <- skeleton$effect
  row_lin   <- skeleton$lin

  # the master effect-evaluation function: free parameter vector -> vector of
  # all reported effects (in the order of eff_df). Written to also work with a
  # complex 'x' (needed for the complex-step jacobian).
  eff_func <- lav_effects_func(specs, effects, nblocks)

  # point estimates
  x_hat <- lav_model_get_parameters(lavmodel = lavmodel, type = "free")
  est <- eff_func(x_hat)
  eff_df$est <- as.numeric(est)

  # standard errors / confidence intervals
  se <- rep(as.numeric(NA), length(est))
  ci_lower <- ci_upper <- rep(as.numeric(NA), length(est))

  a <- (1 - level) / 2
  nboot <- NULL # for the print header (bootstrap only)

  if (se_def == "bootstrap") {
    # bootstrap-based standard errors and percentile confidence intervals,
    # using the bootstrap draws stored in the fitted object
    boot <- try(lav_inspect_boot(object, add_labels = FALSE,
      add_class = FALSE), silent = TRUE)
    error.idx <- attr(boot, "error.idx")
    if (length(error.idx) > 0L) {
      boot <- boot[-error.idx, , drop = FALSE] # drops attributes
    }
    if (inherits(boot, "try-error") || is.null(boot) || NROW(boot) == 0L) {
      lav_msg_warn(gettext(
        "Could not extract the bootstrap draws from the fitted object;
         standard errors are not available."))
    } else {
      nboot <- NROW(boot)
      boot_eff <- apply(boot, 1L, eff_func)
      if (length(est) == 1L) {
        boot_eff <- matrix(boot_eff, ncol = 1L)
      } else {
        boot_eff <- t(boot_eff)
      }
      boot_eff[!is.finite(boot_eff)] <- as.numeric(NA)
      # bootstrap standard error, using the same MLE-normalised covariance
      # (cov * (nboot - 1)/nboot) that lavaan uses elsewhere for bootstrap SEs
      se <- apply(boot_eff, 2L, sd, na.rm = TRUE) * sqrt((nboot - 1) / nboot)
      # bootstrap confidence interval; the per-quantity computation is shared
      # with parameterEstimates() via lav_bootstrap_ci() (see lav_bootstrap.R)
      if (identical(boot_ci_type, "bca")) {
        # adjusted bootstrap percentile method: the acceleration constants are
        # obtained from the empirical influence values of the effects
        design <- lav_bootstrap_bca_design(object, error.idx)
        stopifnot(nrow(design) == nrow(boot_eff))
        acc <- lav_bootstrap_acceleration(boot_eff, design)
        eff_ci <- lav_bootstrap_ci(boot_t = boot_eff, t0 = est,
          boot_ci_type = "bca", level = level, acc = acc)
      } else {
        # For the "norm" method, the bias is evaluated at the mean of the
        # bootstrapped (free) parameters, exactly as in parameterEstimates().
        bias <- if (identical(boot_ci_type, "norm")) {
          as.numeric(eff_func(colMeans(boot, na.rm = TRUE))) - est
        } else {
          NULL
        }
        eff_ci <- lav_bootstrap_ci(boot_t = boot_eff, t0 = est,
          boot_ci_type = boot_ci_type, level = level, se = se, bias = bias)
      }
      ci_lower <- eff_ci[, 1L]
      ci_upper <- eff_ci[, 2L]
    }
  } else if (se_def %in% c("delta", "monte.carlo")) {
    # (co)variance matrix of the free parameters
    vcov <- try(unclass(lavInspect(object, "vcov")), silent = TRUE)
    if (inherits(vcov, "try-error") || is.null(vcov)) {
      lav_msg_warn(gettext(
        "Could not extract the variance-covariance matrix of the parameters;
         standard errors are not available."))
    } else {
      # reduce to the 'free' (compact) space if needed (ceq.simple.only)
      vcov <- lav_model_vcov_unco_to_free(lavmodel, vcov)

      if (se_def == "delta") {
        jac <- try(lav_func_jacobian_complex(func = eff_func, x = x_hat),
          silent = TRUE)
        if (inherits(jac, "try-error")) {
          jac <- lav_func_jacobian_simple(func = eff_func, x = x_hat)
        }
        eff.cov <- jac %*% vcov %*% t(jac)
        diag.cov <- diag(eff.cov)
        diag.cov[diag.cov < 0] <- as.numeric(NA)
        se <- sqrt(diag.cov)
        # symmetric normal-theory confidence interval
        ci_lower <- est + qnorm(a) * se
        ci_upper <- est + qnorm(1 - a) * se
      } else { # monte.carlo
        mc <- lav_effects_mc_draws(x_hat, vcov, r = mc_r, seed = mc_seed)
        mc.eff <- apply(mc, 1L, eff_func)
        if (length(est) == 1L) {
          mc.eff <- matrix(mc.eff, ncol = 1L)
        } else {
          mc.eff <- t(mc.eff)
        }
        mc.eff[!is.finite(mc.eff)] <- as.numeric(NA)
        se <- apply(mc.eff, 2L, sd, na.rm = TRUE)
        # percentile confidence interval
        mc.ci <- apply(mc.eff, 2L, quantile, probs = c(a, 1 - a),
          na.rm = TRUE, names = FALSE)
        ci_lower <- mc.ci[1L, ]
        ci_upper <- mc.ci[2L, ]
      }
    }
  }
  eff_df$se <- as.numeric(se)

  if (zstat || pvalue) {
    z <- eff_df$est / eff_df$se
    if (zstat) {
      eff_df$z <- z
    }
    if (pvalue) {
      eff_df$pvalue <- 2 * (1 - pnorm(abs(z)))
    }
  }

  if (ci && se_def != "none") {
    eff_df$ci_lower <- as.numeric(ci_lower)
    eff_df$ci_upper <- as.numeric(ci_upper)
  }

  if (output == "list") {
    return(lav_effects_as_list(eff_df, specs, effects, lavmodel,
      row_block = row_block, row_eff = row_eff, row_lin = row_lin))
  }

  # standardized effects? Add std.lv/std.all/std.nox/std.user columns (point
  # estimates), just like the standardized= argument of parameterEstimates().
  std_types <- lav_effects_std_types(object, standardized)
  if (length(std_types$types) > 0L) {
    eff_df <- lav_effects_add_standardized(eff_df, object, lavmodel, specs,
      row_block = row_block, row_lin = row_lin,
      std_types = std_types$types, ov_std_user = std_types$ov_std_user,
      cov_std = cov_std)
  }

  # attributes for printing
  attr(eff_df, "se_def") <- se_def
  attr(eff_df, "level")  <- if (ci) level else NULL
  attr(eff_df, "R")      <- if (se_def == "monte.carlo") mc_r else NULL
  attr(eff_df, "nboot")  <- nboot
  attr(eff_df, "boot_ci_type") <- if (se_def == "bootstrap") boot_ci_type
    else NULL
  attr(eff_df, "ngroups") <- lavdata@ngroups

  if (add_class) {
    class(eff_df) <- c("lavaan.effects", "lavaan.data.frame", "data.frame")
  }

  eff_df
}

# resolve the 'standardized' argument (as in parameterEstimates()) into a
# vector of standardization types and an optional set of user-specified
# observed variables to standardize (std.user)
lav_effects_std_types <- function(object, standardized) {
  ov_std_user <- NULL
  if (is.logical(standardized)) {
    if (isTRUE(standardized)) {
      types <- c("std.lv", "std.all")
      if (length(lavNames(object, "ov.x")) && object@Options$fixed.x) {
        types <- c(types, "std.nox")
      }
    } else {
      types <- character(0L)
    }
  } else {
    standardized <- as.character(standardized)
    std.set <- c("std.lv", "std.all", "std.nox")
    if (length(standardized) > 0L &&
        !any(tolower(standardized) %in% std.set)) {
      # interpret 'standardized' as a vector of observed variable names
      ov.all <- lavNames(object, "ov")
      bad <- standardized[!(standardized %in% ov.all)]
      if (length(bad) > 0L) {
        lav_msg_warn(gettextf(
          "standardized= argument contains unknown observed variable(s): %s",
          paste(bad, collapse = ", ")))
      }
      ov_std_user <- standardized[standardized %in% ov.all]
      if (length(ov_std_user) == 0L) {
        types <- character(0L)
      } else {
        types <- c("std.lv", "std.user")
      }
    } else {
      types <- tolower(standardized)
      if ("std.nox" %in% types) {
        if (length(lavNames(object, "ov.x")) == 0L) {
          lav_msg_note(gettext(
            "`std.nox' unavailable without fixed exogenous predictors"))
          types <- setdiff(types, "std.nox")
        } else if (!object@Options$fixed.x) {
          lav_msg_note(gettext("`std.nox' unavailable when fixed.x = FALSE"))
          types <- setdiff(types, "std.nox")
        }
      }
    }
  }
  list(types = types, ov_std_user = ov_std_user)
}

# the (model-implied) standard-deviation vector used to standardize the BETA
# variables of a block, for a given standardization 'type'. Variables that are
# not standardized under this type get a SD of 1.
lav_effects_std_sd_vec <- function(var_names, sd_all, type, lv_names,
                                   ovx_names, ov_std_user) {
  sdv <- sd_all
  if (type == "std.lv") {
    # only the (true) latent variables are standardized
    sdv[!(var_names %in% lv_names)] <- 1
  } else if (type == "std.nox") {
    # all variables except the exogenous observed 'x'
    sdv[var_names %in% ovx_names] <- 1
  } else if (type == "std.user") {
    # latents + user-specified observed variables
    keep <- (var_names %in% lv_names) | (var_names %in% ov_std_user)
    sdv[!keep] <- 1
  }
  # type == "std.all": standardize everything (sdv = sd_all)
  sdv
}

# add standardized effect columns (point estimates) to the effects data.frame.
# The standardized effect of j on i equals effect[i,j] * sd_j / sd_i, where the
# SDs are the model-implied standard deviations of the BETA variables under the
# requested standardization type.
lav_effects_add_standardized <- function(eff_df, object, lavmodel, specs,
                                         row_block, row_lin, std_types,
                                         ov_std_user = NULL, cov_std = TRUE) {
  nblocks <- lavmodel@nblocks
  nmat <- lavmodel@nmat
  glist_names <- names(lavmodel@GLIST)
  partable <- object@ParTable

  # model-implied (co)variances of the BETA (eta) variables, per block
  veta <- lav_model_veta(lavmodel = lavmodel, glist = lavmodel@GLIST)

  # per block: 'all' SDs for the row (eta) and column (eta + exogenous x)
  # variables, plus the variable classifications needed by each type
  sd.row.all <- vector("list", nblocks)
  sd.col.all <- vector("list", nblocks)
  lv_names <- vector("list", nblocks)
  ovx_names <- vector("list", nblocks)
  for (b in seq_len(nblocks)) {
    sp <- specs[[b]]
    if (is.null(sp)) {
      next
    }
    v2 <- diag(veta[[b]])
    v2[v2 < 0] <- as.numeric(NA)
    sd.eta <- sqrt(v2)
    sd.row.all[[b]] <- sd.eta
    if (sp$nr.col > sp$nr) {
      # exogenous covariate SDs (from cov.x), in the GAMMA column order
      mm_in_block <- seq_len(nmat[b]) + cumsum(c(0L, nmat))[b]
      covx.idx <- mm_in_block[glist_names[mm_in_block] == "cov.x"]
      if (length(covx.idx) > 0L) {
        sd.x <- sqrt(diag(lavmodel@GLIST[[covx.idx[1L]]]))
      } else {
        sd.x <- rep(1, sp$nr.col - sp$nr)
      }
      sd.col.all[[b]] <- c(sd.eta, sd.x)
    } else {
      sd.col.all[[b]] <- sd.eta
    }
    lv_names[[b]] <- lav_pt_vnames(partable, "lv", block = b)
    ovx_names[[b]] <- lav_pt_vnames(partable, "ov.x", block = b)
  }

  nr.rows <- nrow(eff_df)
  for (ty in std_types) {
    # per-block SD vectors (row and column) for this standardization type
    sd.row <- vector("list", nblocks)
    sd.col <- vector("list", nblocks)
    for (b in seq_len(nblocks)) {
      if (is.null(specs[[b]])) {
        next
      }
      sd.row[[b]] <- lav_effects_std_sd_vec(specs[[b]]$row.names,
        sd.row.all[[b]], ty, lv_names[[b]], ovx_names[[b]], ov_std_user)
      sd.col[[b]] <- lav_effects_std_sd_vec(specs[[b]]$col.names,
        sd.col.all[[b]], ty, lv_names[[b]], ovx_names[[b]], ov_std_user)
    }
    std.col <- rep(as.numeric(NA), nr.rows)
    for (r in seq_len(nr.rows)) {
      b <- row_block[r]
      nr <- specs[[b]]$nr
      lin <- row_lin[r]
      i <- ((lin - 1L) %% nr) + 1L
      j <- ((lin - 1L) %/% nr) + 1L
      std.col[r] <- eff_df$est[r] * sd.col[[b]][j] / sd.row[[b]][i]
    }
    eff_df[[ty]] <- std.col
  }

  eff_df
}

# print method for the object returned by lavEffects()
lav_effects_print <- function(x, ..., nd = 3L) {
  # if the object has been subsetted (and lost the columns/attributes that
  # this method relies on), fall back to a plain data.frame print
  if (is.null(x$effect) || is.null(x$lhs) || is.null(x$est)) {
    y <- x
    class(y) <- "data.frame"
    print(y, ...)
    return(invisible(x))
  }

  # number formats (mirrors lav_parameterestimates_print)
  num.format  <- paste("%", max(8L, nd + 5L), ".", nd, "f", sep = "")
  char.format <- paste("%", max(8L, nd + 5L), "s", sep = "")

  se_def  <- attr(x, "se_def")
  level   <- attr(x, "level")
  ngroups <- attr(x, "ngroups")
  if (is.null(ngroups)) {
    ngroups <- 1L
  }

  # header
  cat("\nEffects:\n\n")
  if (!is.null(se_def)) {
    c1 <- c("Standard errors")
    c2 <- c(switch(se_def, monte.carlo = "Monte Carlo", delta = "Delta",
      bootstrap = "Bootstrap", none = "None", se_def))
    if (identical(se_def, "monte.carlo") && !is.null(attr(x, "R"))) {
      c1 <- c(c1, "Monte Carlo draws")
      c2 <- c(c2, format(attr(x, "R")))
    }
    if (identical(se_def, "bootstrap") && !is.null(attr(x, "nboot"))) {
      c1 <- c(c1, "Bootstrap draws")
      c2 <- c(c2, format(attr(x, "nboot")))
    }
    if (identical(se_def, "bootstrap") && !is.null(attr(x, "boot_ci_type"))) {
      c1 <- c(c1, "Bootstrap CI type")
      c2 <- c(c2, attr(x, "boot_ci_type"))
    }
    if (!is.null(level) && !identical(se_def, "none")) {
      c1 <- c(c1, "Confidence level")
      c2 <- c(c2, format(level))
    }
    c1 <- format(c1, width = 33L)
    c2 <- format(c2, width = 14L + max(0, (nd - 3L)) * 4L, justify = "right")
    m1 <- cbind(c1, c2, deparse.level = 0)
    colnames(m1) <- rep("", ncol(m1))
    rownames(m1) <- rep(" ", nrow(m1))
    write.table(m1, row.names = TRUE, col.names = FALSE, quote = FALSE)
  }

  effect.title <- c(total = "Total Effects", indirect = "Indirect Effects",
    direct = "Direct Effects")

  # group/level handling
  group.vec <- if (!is.null(x$group)) x$group else rep(1L, nrow(x))
  ugroups <- unique(group.vec)

  for (g in ugroups) {
    if (ngroups > 1L) {
      cat("\n\nGroup ", which(ugroups == g), " [", as.character(g), "]:\n",
        sep = "")
    }
    g.idx <- which(group.vec == g)

    for (e in c("total", "indirect", "direct")) {
      e.idx <- g.idx[x$effect[g.idx] == e]
      if (length(e.idx) == 0L) {
        next
      }
      cat("\n", effect.title[e], ":\n", sep = "")

      sub <- x[e.idx, , drop = FALSE]
      # left label: "outcome ~ predictor"
      labels <- paste(format(sub$lhs), sub$op, format(sub$rhs))

      # numeric columns to show
      cols <- list()
      cols[["Estimate"]] <- sprintf(num.format, sub$est)
      if (!is.null(sub$se)) {
        cols[["Std.Err"]] <- sprintf(num.format, sub$se)
      }
      if (!is.null(sub$z)) {
        cols[["z-value"]] <- sprintf(num.format, sub$z)
      }
      if (!is.null(sub$pvalue)) {
        cols[["P(>|z|)"]] <- sprintf(num.format, sub$pvalue)
      }
      if (!is.null(sub$ci_lower)) {
        cols[["ci_lower"]] <- sprintf(num.format, sub$ci_lower)
      }
      if (!is.null(sub$ci_upper)) {
        cols[["ci_upper"]] <- sprintf(num.format, sub$ci_upper)
      }
      # standardized columns (if any)
      for (sc in c("std.lv", "std.all", "std.nox", "std.user")) {
        if (!is.null(sub[[sc]])) {
          cols[[sc]] <- sprintf(num.format, sub[[sc]])
        }
      }

      m <- do.call(cbind, cols)
      if (!is.null(sub$se)) {
        # blank out se/z/p/ci when se is NA (no SE available)
        na.se <- which(is.na(sub$se))
        if (length(na.se) > 0L) {
          for (cn in intersect(c("Std.Err", "z-value", "P(>|z|)",
            "ci_lower", "ci_upper"), colnames(m))) {
            m[na.se, cn] <- sprintf(char.format, "")
          }
        }
        # blank out se/z/p for fixed effects (se == 0); keep the ci columns
        # (which equal the point estimate)
        zero.se <- which(!is.na(sub$se) & sub$se == 0)
        if (length(zero.se) > 0L) {
          for (cn in intersect(c("Std.Err", "z-value", "P(>|z|)"),
            colnames(m))) {
            m[zero.se, cn] <- sprintf(char.format, "")
          }
        }
      }
      rownames(m) <- paste("   ", format(labels, justify = "left"))
      print.default(m, quote = FALSE, right = TRUE)
    }
  }
  cat("\n")

  invisible(x)
}

# transitive closure (reachability) of a boolean adjacency matrix, where
# a[i, j] = TRUE means there is a direct edge j -> i. Returns t1 with
# t1[i, j] = TRUE if j reaches i via one or more directed edges. The diagonal
# is forced to FALSE (we are not interested in self-loops).
lav_effects_reachable <- function(a) {
  n <- nrow(a)
  r <- (a != 0)
  if (n > 1L) {
    # Warshall: r[i, j] |= r[i, k] & r[k, j]
    for (k in seq_len(n)) {
      r <- r | (r[, k, drop = FALSE] %*% r[k, , drop = FALSE] > 0L)
    }
  }
  diag(r) <- FALSE
  r
}

# build the master effect-evaluation function:
# free parameter vector 'x' -> numeric (or complex) vector of all reported
# effects, in the order: block, then effect type, then (outcome, predictor)
lav_effects_func <- function(specs, effects, nblocks) {
  function(x) {
    out <- vector(mode = if (is.complex(x)) "complex" else "numeric",
      length = 0L)
    for (b in seq_len(nblocks)) {
      sp <- specs[[b]]
      if (is.null(sp)) {
        next
      }
      b <- sp$b0
      if (is.complex(x)) {
        storage.mode(b) <- "complex"
      }
      if (length(sp$beta_m_idx) > 0L) {
        b[sp$beta_m_idx] <- x[sp$beta_x_idx]
      }
      diag(b) <- 0
      nr <- sp$nr
      I.B.inv <- solve(diag(nr) - b)

      # eta -> eta effects (square block)
      total.ee <- I.B.inv - diag(nr)
      if (is.null(sp$g0)) {
        total  <- total.ee
        direct <- b
      } else {
        # eta <- Gamma x effects: total = (I-B)^-1 Gamma, direct = Gamma,
        # indirect = total - direct
        g <- sp$g0
        if (is.complex(x)) {
          storage.mode(g) <- "complex"
        }
        if (length(sp$gamma_m_idx) > 0L) {
          g[sp$gamma_m_idx] <- x[sp$gamma_x_idx]
        }
        total  <- cbind(total.ee, I.B.inv %*% g)
        direct <- cbind(b, g)
      }
      indirect <- total - direct

      for (e in effects) {
        cell <- sp$cells[[e]]
        if (is.null(cell)) {
          next
        }
        MAT <- switch(e,
          total    = total,
          indirect = indirect,
          direct   = direct)
        out <- c(out, MAT[cell$lin])
      }
    }
    out
  }
}

# draw r samples from MVN(x_hat, vcov), saving/restoring the RNG state
lav_effects_mc_draws <- function(x_hat, vcov, r = 20000L, seed = NULL) {
  r <- as.integer(r)
  stopifnot(r > 0L)

  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    saved.seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit(assign(".Random.seed", saved.seed, envir = .GlobalEnv), add = TRUE)
  } else {
    on.exit({
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }

  vcov.sym <- 0.5 * (vcov + t(vcov))
  mc <- lav_mvrnorm(n = r, mu = x_hat, sigma_1 = vcov.sym,
    check_symmetry = FALSE)
  if (!is.matrix(mc)) {
    mc <- matrix(mc, nrow = r)
  }
  mc
}

# build the (ordered) skeleton data.frame: one row per reported effect, in the
# same order as lav_effects_func() returns its values
lav_effects_skeleton <- function(specs, effects, lavdata, lavmodel) {
  nblocks <- lavmodel@nblocks
  ngroups <- lavdata@ngroups
  nlevels <- lavdata@nlevels

  eff.col <- character(0L)
  lhs.col <- character(0L)
  rhs.col <- character(0L)
  blk.col <- integer(0L)
  lin.col <- integer(0L)

  for (b in seq_len(nblocks)) {
    sp <- specs[[b]]
    if (is.null(sp)) {
      next
    }
    for (e in effects) {
      cell <- sp$cells[[e]]
      if (is.null(cell)) {
        next
      }
      n.e <- length(cell$lin)
      eff.col <- c(eff.col, rep(e, n.e))
      lhs.col <- c(lhs.col, sp$row.names[cell$i]) # outcome (eta)
      rhs.col <- c(rhs.col, sp$col.names[cell$j]) # predictor (eta or x)
      blk.col <- c(blk.col, rep(b, n.e))
      lin.col <- c(lin.col, cell$lin)
    }
  }

  df <- data.frame(
    effect = eff.col,
    lhs    = lhs.col,
    op     = rep("~", length(lhs.col)),
    rhs    = rhs.col,
    stringsAsFactors = FALSE
  )

  # group/level/block columns (only when relevant)
  if (nblocks > 1L) {
    if (ngroups > 1L) {
      group.idx <- ((blk.col - 1L) %/% max(nlevels, 1L)) + 1L
      glabels <- lavdata@group.label
      if (length(glabels) >= max(group.idx)) {
        df$group <- glabels[group.idx]
      } else {
        df$group <- group.idx
      }
    }
    if (nlevels > 1L) {
      level.idx <- ((blk.col - 1L) %% nlevels) + 1L
      llabels <- lavdata@level.label
      if (length(llabels) >= max(level.idx)) {
        df$level <- llabels[level.idx]
      } else {
        df$level <- level.idx
      }
    }
    if (ngroups <= 1L && nlevels <= 1L) {
      df$block <- blk.col
    }
  }

  list(df = df, effect = eff.col, block = blk.col, lin = lin.col)
}

# an empty result (no effects) with the expected columns
lav_effects_empty_df <- function(effects, zstat, pvalue) {
  df <- data.frame(
    effect = character(0L),
    lhs    = character(0L),
    op     = character(0L),
    rhs    = character(0L),
    est    = numeric(0L),
    se     = numeric(0L),
    stringsAsFactors = FALSE
  )
  if (zstat) {
    df$z <- numeric(0L)
  }
  if (pvalue) {
    df$pvalue <- numeric(0L)
  }
  df$ci_lower <- numeric(0L)
  df$ci_upper <- numeric(0L)
  class(df) <- c("lavaan.effects", "lavaan.data.frame", "data.frame")
  df
}

# convert the result to a list of matrices (one per block and per effect type),
# using the per-row block/effect/linear-cell bookkeeping from the skeleton
lav_effects_as_list <- function(eff_df, specs, effects, lavmodel,
                                row_block, row_eff, row_lin) {
  nblocks <- lavmodel@nblocks
  out <- list()
  for (b in seq_len(nblocks)) {
    sp <- specs[[b]]
    if (is.null(sp)) {
      next
    }
    nr <- sp$nr
    nr.col <- sp$nr.col
    dn <- list(sp$row.names, sp$col.names)
    block.list <- list()
    for (e in effects) {
      rows <- which(row_block == b & row_eff == e)
      if (length(rows) == 0L) {
        next
      }
      M.est <- matrix(0, nr, nr.col, dimnames = dn)
      M.se  <- matrix(as.numeric(NA), nr, nr.col, dimnames = dn)
      M.est[row_lin[rows]] <- eff_df$est[rows]
      if (!is.null(eff_df$se)) {
        M.se[row_lin[rows]] <- eff_df$se[rows]
      }
      block.list[[e]] <- list(est = M.est, se = M.se)
    }
    if (nblocks > 1L) {
      out[[paste0("block.", b)]] <- block.list
    } else {
      out <- block.list
    }
  }
  out
}
