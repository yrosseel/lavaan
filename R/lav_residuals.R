# residual diagnostics

# two types:
# 1) residuals for summary statistics
# 2) case-wise residuals

# this (new) version written around Aug/Sept 2018 for 0.6-3
# - based on obsList (inspect_sampstat) and estList (inspect_implied)
# - pre-scaling for type = "cor.bollen" and type = "cor.bentler"
# - summary statistics: rmr, srmr, crmr, urmr, usrmr, ucrmr; standard errors,
#       confidence intervals (for u(cs)rmr),
#       z-statistics (exact test, close test), p-values
# - type = "normalized" is based on lav_model_h1_acov(), and should now work
#   for all estimators
# - type = "standardized" now uses the correct formula, and should work for
#   for all estimators
# - type = "standardized.mplus" uses the simplified Mplus/LISREL version,
#   often resulting in NAs due to negative var(resid) estimates
#   (this was "standardized" in lavaan < 0.6.3

# WARNING: only partial support for conditional.x = TRUE!!
# - in categorical case: we only compute summary statistics, using cor + th
#   (no var, slopes, ...)
# - twolevel not supported here; see lav_fit_srmr.R, where we convert to
#   the unconditional setting

# - change 0.6-6: we enforce observed.information = "h1" to ensure 'Q' is a
#                 projection matrix (see lav_residuals_acov)

# - change 0.6-13: fixed.x = TRUE is ignored (to conform with 'tradition')

setMethod(
  "residuals", "lavaan",
  function(object, type = "raw", labels = TRUE, ...) {
    dotdotdot <- list(...)
    if (length(dotdotdot) > 0L) {
      for (j in seq_along(dotdotdot)) {
        lav_msg_warn(gettextf(
          "Unknown argument %s for %s", sQuote(names(dotdotdot)[j]),
          sQuote("residuals"))
        )
      }
    }
    # lowercase type
    type <- tolower(type)

    # type = "casewise"
    if (type %in% c("casewise", "case", "obs", "observations", "ov")) {
      return(lav_residuals_casewise(object, labels = labels))
    } else {
      out <- lav_residuals(
        object = object, type = type, h1 = TRUE,
        add_type = TRUE,
        rename_cov_cor = FALSE, # should become FALSE!
        # after packages (eg jmv)
        # have adapted 0.6-3 style
        add_labels = labels, add_class = TRUE,
        drop_list_single_group = TRUE
      )
    }

    out
  }
)

setMethod(
  "resid", "lavaan",
  function(object, type = "raw", ...) {
    dotdotdot <- list(...)
    if (length(dotdotdot) > 0L) {
      for (j in seq_along(dotdotdot)) {
        lav_msg_warn(gettextf(
          "Unknown argument %s for %s", sQuote(names(dotdotdot)[j]),
          sQuote("resid"))
        )
      }
    }
    residuals(object, type = type, labels = TRUE)
  }
)


# user-visible function
lavResiduals <- function(object, type = "cor.bentler", custom.rmr = NULL,  # nolint start
                         se = FALSE, zstat = TRUE, summary = TRUE,
                         h1.acov = "unstructured",
                         add.type = TRUE, add.labels = TRUE, add.class = TRUE,
                         drop.list.single.group = TRUE,
                         maximum.number = length(res_vech),
                         output = "list") {                                # nolint end
  out <- lav_residuals(
    object = object, type = type, h1 = TRUE,
    custom_rmr = custom.rmr, se = se, zstat = zstat,
    summary = summary,
    summary_options = list(
      se = TRUE, zstat = TRUE, pvalue = TRUE,
      unbiased = TRUE, unbiased_se = TRUE, unbiased_ci = TRUE,
      unbiased_ci_level = 0.90, unbiased_zstat = TRUE,
      unbiased_test_val = 0.05, unbiased_pvalue = TRUE
    ),
    h1_acov = h1.acov, add_type = add.type,
    add_labels = add.labels, add_class = add.class,
    drop_list_single_group = drop.list.single.group
  )

  # no pretty printing yet...
  if (output == "table") {
    # new in 0.6-18: handle multiple blocks
    nblocks <- lav_pt_nblocks(object@ParTable)
    out_list <- vector("list", length = nblocks)
    for (block in seq_len(nblocks)) {
      if (nblocks == 1L) {
        res <- out$cov
      } else {
        res <- out[[block]]$cov
      }
      # extract only below-diagonal elements
      res_vech <- lav_mat_vech(res, diagonal = FALSE)

      # get names
      p <- nrow(res)
      names_1 <- colnames(res)

      nam <- expand.grid(
        names_1,
        names_1
      )[lav_mat_vech_idx(p, diagonal = FALSE), ]
      names_vech <- paste(nam[, 1], "~~", nam[, 2], sep = "")

      # create table
      tab <- data.frame(
        name = names_vech, res = round(res_vech, 3),
        stringsAsFactors = FALSE
      )

      # sort table
      idx <- sort.int(abs(tab$res),
        decreasing = TRUE,
        index.return = TRUE
      )$ix
      out_sorted <- tab[idx, ]

      # show first rows only
      maximum_number <- maximum.number
      if (maximum_number == 0L || maximum_number > length(res_vech)) {
        maximum_number <- length(res_vech)
      }
      out_list[[block]] <- out_sorted[seq_len(maximum_number), ]
    }
    if (nblocks == 1L) {
      out <- out_list[[1]]
    } else {
      out <- out_list
      names(out) <- object@Data@block.label
    }
  } else {
    # list -> nothing to do
  }

  out
}

# main function
lav_residuals <- function(object, type = "raw", h1 = TRUE, custom_rmr = NULL,
                          se = FALSE, zstat = FALSE, summary = FALSE,
                          summary_options = list(
                            se = TRUE, zstat = TRUE,
                            pvalue = TRUE, unbiased = TRUE, unbiased_se = TRUE,
                            unbiased_ci = TRUE, unbiased_ci_level = 0.90,
                            unbiased_zstat = FALSE, unbiased_test_val = 0.05,
                            unbiased_pvalue = FALSE
                          ),
                          h1_acov = "unstructured", add_type = FALSE,
                          rename_cov_cor = FALSE,
                          add_labels = FALSE, add_class = FALSE,
                          drop_list_single_group = FALSE) {
  # check object
  object <- lav_object_check_version(object)

  # working version summary options
  summary_options_1 <- summary_options

  # type
  type <- tolower(type)[1]

  # check type
  if (!type %in% c(
    "raw", "cor", "cor.bollen", "cor.bentler", "cor.eqs",
    "rmr", "srmr", "crmr",
    "normalized", "standardized", "standardized.mplus"
  )) {
    lav_msg_stop(gettext("unknown argument for type:"), dQuote(type))
  }

  # if cor, choose 'default'
  if (type == "cor") {
    if (object@Options$mimic == "EQS") {
      type <- "cor.bentler"
    } else {
      type <- "cor.bollen"
    }
  }
  if (type == "cor.eqs") {
    type <- "cor.bentler"
  }
  if (type == "rmr") {
    type <- "raw"
  }
  if (type == "srmr") {
    type <- "cor.bentler"
  }
  if (type == "crmr") {
    type <- "cor.bollen"
  }

  # slots
  lavdata <- object@Data
  lavmodel <- object@Model

  # change options if multilevel (for now)
  if (lavdata@nlevels > 1L) {
    zstat <- se <- FALSE
    summary <- FALSE
  }

  # change options if categorical (for now)
  if (lavmodel@categorical) {
    # only if conditional.x = FALSE AND no continuous endogenous variables
    # -> only the simple setting where we only have thresholds and
    #    correlations

    # As soon as we add continuous variables, we get means/variances too,
    # and we need to decide how WLS.obs/WLS.est/WLS.V will then map to
    # the output of lavInspect(fit, "implied") and
    # lavInspect(fit, "sampstat")

    if (lavmodel@conditional.x || length(unlist(lavmodel@num.idx)) > 0L) {
      zstat <- se <- FALSE
      summary <- FALSE
      summary_options_1 <- list(
        se = FALSE, zstat = FALSE,
        pvalue = FALSE, unbiased = FALSE,
        unbiased_se = FALSE,
        unbiased_ci = FALSE, unbiased_ci_level = 0.90,
        unbiased_zstat = FALSE, unbiased_test_val = 0.05,
        unbiased_pvalue = FALSE
      )
    }
  }

  # change options if conditional.x (for now)
  if (!lavmodel@categorical && lavmodel@conditional.x) {
    zstat <- se <- FALSE
    summary <- FALSE
    summary_options_1 <- list(
      se = FALSE, zstat = FALSE,
      pvalue = FALSE, unbiased = FALSE,
      unbiased_se = FALSE,
      unbiased_ci = FALSE, unbiased_ci_level = 0.90,
      unbiased_zstat = FALSE, unbiased_test_val = 0.05,
      unbiased_pvalue = FALSE
    )
  }

  # change options if estimator = "IV" (for now)
  if (lavmodel@estimator == "IV") {
    zstat <- se <- FALSE
    summary <- FALSE
    summary_options_1 <- list(
      se = FALSE, zstat = FALSE,
      pvalue = FALSE, unbiased = FALSE,
      unbiased_se = FALSE,
      unbiased_ci = FALSE, unbiased_ci_level = 0.90,
      unbiased_zstat = FALSE, unbiased_test_val = 0.05,
      unbiased_pvalue = FALSE
    )
  }

  # observed and fitted sample statistics
  obs_list <- lav_inspect_sampstat(object,
    h1 = h1,
    add_labels = add_labels, add_class = add_class,
    drop_list_single_group = FALSE
  )
  est_list <- lav_inspect_implied(object,
    add_labels = add_labels, add_class = add_class,
    drop_list_single_group = FALSE
  )
  # blocks
  nblocks <- length(obs_list)

  # pre-scale?
  if (type %in% c("cor.bentler", "cor.bollen")) {
    for (b in seq_len(nblocks)) {
      var_obs <- if (lavmodel@conditional.x) {
        diag(obs_list[[b]][["res.cov"]])
      } else {
        diag(obs_list[[b]][["cov"]])
      }
      var_est <- if (lavmodel@conditional.x) {
        diag(est_list[[b]][["res.cov"]])
      } else {
        diag(est_list[[b]][["cov"]])
      }

      # rescale obsList
      obs_list[[b]] <-
        lav_residuals_rescale(x = obs_list[[b]], diag_cov = var_obs)
      # rescale estList
      if (type == "cor.bentler") { # use obsList
        est_list[[b]] <-
          lav_residuals_rescale(x = est_list[[b]], diag_cov = var_obs)
      } else if (type == "cor.bollen") { # use estList for COV only
        est_list[[b]] <- lav_residuals_rescale(
          x = est_list[[b]],
          diag_cov = var_est, diag_cov2 = var_obs
        )
      }
    }
  }

  # compute residuals: (observed - implied)
  res_list <- vector("list", length = nblocks)
  for (b in seq_len(nblocks)) {
    res_list[[b]] <- lapply(seq_along(obs_list[[b]]),
      FUN = function(el) {
        obs_list[[b]][[el]] - est_list[[b]][[el]]
      }
    )
    # always name the elements, even if add.labels = FALSE
    names_1 <- names(obs_list[[b]])
    names(res_list[[b]]) <- names_1
  }

  # do we need seList?
  if (se || zstat) {
    se_list <- lav_residuals_se(object,
      type = type, z_type = "standardized",
      h1_acov = h1_acov,
      add_class = add_class, add_labels = add_labels
    )
  } else if (type %in% c("normalized", "standardized", "standardized.mplus")) {
    se_list <- lav_residuals_se(object,
      type = "raw", z_type = type,
      h1_acov = h1_acov,
      add_class = add_class, add_labels = add_labels
    )
  } else {
    se_list <- NULL
  }

  # normalize/standardize?
  if (type %in% c("normalized", "standardized", "standardized.mplus")) {
    for (b in seq_len(nblocks)) {
      if (add_labels) {
        names_1 <- names(res_list[[b]])
      }
      res_list[[b]] <- lapply(seq_along(res_list[[b]]),
        FUN = function(el) {
          a_1 <- res_list[[b]][[el]]
          b_1 <- se_list[[b]][[el]]
          near_zero_idx <- which(abs(a_1) < 1e-05)
          if (length(near_zero_idx) > 0L) {
            b_1[near_zero_idx] <- 1
          }
          a_1 / b_1
        }
      )
      if (add_labels) {
        names(res_list[[b]]) <- names_1
      }
    }
  }

  # add se
  res_list_orig <- res_list
  if (se) {
    for (b in seq_len(nblocks)) {
      names_res <- names(res_list[[b]])
      names_se <- paste0(names_res, ".se")
      res_list[[b]] <- c(res_list[[b]], se_list[[b]])
      names(res_list[[b]]) <- c(names_res, names_se)
    }
  }


  # add zstat
  if (zstat) {
    for (b in seq_len(nblocks)) {
      names_res <- names(res_list[[b]])
      names_z <- paste0(names(res_list_orig[[b]]), ".z")
      tmp <- lapply(seq_along(res_list_orig[[b]]),
        FUN = function(el) {
          a_1 <- res_list_orig[[b]][[el]]
          b_1 <- se_list[[b]][[el]]
          # NOTE: which threshold should we use?
          # used to be 1e-05
          # changed to 1e-04 in 0.6-4
          near_zero_idx <- which(abs(a_1) < 1e-04)
          if (length(near_zero_idx) > 0L) {
            # B[near_zero_idx] <- as.numeric(NA)
            b_1[near_zero_idx] <- 1.0
          }
          a_1 / b_1
        }
      )
      res_list[[b]] <- c(res_list[[b]], tmp)
      names(res_list[[b]]) <- c(names_res, names_z)
    }
  }

  # add summary statistics (rms, mabs)
  if (summary) {
    summary_options_1 <- lav_snake_case(summary_options_1)
    args <- c(
      list(
        object = object, type = type, h1_acov = h1_acov,
        add_class = add_class, custom_rmr = custom_rmr
      ),
      summary_options_1
    )
    sum_stat <- do.call("lav_residuals_summary", args)
    for (b in seq_len(nblocks)) {
      names_1 <- names(res_list[[b]])
      res_list[[b]] <- c(res_list[[b]], list(sum_stat[[b]][[1]])) # only 1
      names_1 <- c(names_1, "summary")
      names(res_list[[b]]) <- names_1
    }
  }

  # last: add type
  if (add_type) {
    for (b in seq_len(nblocks)) {
      names_1 <- names(res_list[[b]])
      res_list[[b]] <- c(type, res_list[[b]])
      names_1 <- c("type", names_1)
      names(res_list[[b]]) <- names_1
    }
  }

  # optional: rename 'cov' to 'cor' (if type = "cor")
  if (rename_cov_cor && type %in% c("cor.bentler", "cor.bollen")) {
    for (b in seq_len(nblocks)) {
      names_1 <- names(res_list[[b]])
      names_1 <- gsub("cov", "cor", names_1)
      names(res_list[[b]]) <- names_1
    }
  }


  # output
  out <- res_list
  if (nblocks == 1L && drop_list_single_group) {
    out <- out[[1]]
  } else {
    if (lavdata@nlevels == 1L &&
      length(lavdata@group.label) > 0L) {
      names(out) <- unlist(lavdata@group.label)
    } else if (lavdata@nlevels > 1L &&
      length(lavdata@group.label) == 0L) {
      names(out) <- lavdata@level.label
    }
  }

  out
}

# return ACOV as list per group
lav_residuals_acov <- function(object, type = "raw", z_type = "standardized",
                               h1_acov = "unstructured") {
  # check type
  if (z_type %in% c("normalized", "standardized.mplus") && type != "raw") {
    lav_msg_stop(gettextf(
      "z.type = %1$s can only be used with type = %2$s",
      dQuote(z_type), dQuote("raw")))
  }

  # slots
  lavdata <- object@Data
  lavmodel <- object@Model
  lavsamplestats <- object@SampleStats

  # return list per group
  acov_res <- vector("list", length = lavdata@ngroups)

  # compute ACOV for observed h1 sample statistics (ACOV == Gamma/N)
  if (!is.null(lavsamplestats@NACOV[[1]])) {
    nacov_obs <- lavsamplestats@NACOV
                # if this changes, tag @TDJorgensen in commit message
    acov_obs <- lapply(nacov_obs, function(x) x / lavsamplestats@ntotal)
  } else {
    acov_obs <- lav_model_h1_acov(
      lavobject = object,
      h1_information = h1_acov
    )
  }

  # shortcut for normalized
  if (z_type == "normalized") {
    acov_res <- acov_obs
    return(acov_res)
  } else {
    if (z_type == "standardized") {
      a1 <- lav_model_h1_info(object)
      if (lavmodel@estimator == "DWLS" || lavmodel@estimator == "ULS") {
        # A1 is diagonal matrix
        a1 <- lapply(a1, diag)
      }
      if (type %in% c("cor.bentler", "cor.bollen")) {
        sampstat <- lavTech(object, "sampstat")
      }
    } else if (z_type == "standardized.mplus") {
      vcov_1 <- lavTech(object, "vcov")
    }
    mm_delta <- lavTech(object, "delta")
  }

  # for each group, compute ACOV
  for (g in seq_len(lavdata@ngroups)) {
    # group weight
    gw <- object@SampleStats@nobs[[g]] / object@SampleStats@ntotal
                 # if this changes, tag @TDJorgensen in commit message

    if (z_type == "standardized.mplus") { # simplified formula
      # also used by LISREL?
      # see https://www.statmodel.com/download/StandardizedResiduals.pdf

      acov_est_g <- mm_delta[[g]] %*% vcov_1 %*% t(mm_delta[[g]])
      acov_res[[g]] <- acov_obs[[g]] - acov_est_g
    } else if (z_type == "standardized") {
      # see Ogasawara (2001) using Bentler & Dijkstra (1985) eq 1.7.4

      # NVarCov, but always 'not' robust
      #
      # new in 0.6-6: to ensure Q is a projection matrix, we
      #               force observed.information = "h1"
      #               (only needed if information is observed)
      this_options <- object@Options
      this_options$observed.information[1] <- "h1"
      a0_g_inv <- lav_model_info(
        lavmodel = lavmodel,
        lavsamplestats = object@SampleStats,
        lavdata = lavdata,
        lavcache = object@Cache,
        lavimplied = object@implied,
        lavh1 = object@h1,
        lavoptions = this_options,
        extra = FALSE,
        augmented = TRUE,
        inverted = TRUE,
        use_ginv = TRUE
      )

      acov_est_g <- gw * (mm_delta[[g]] %*% a0_g_inv %*% t(mm_delta[[g]]))
      q_1 <- diag(nrow = nrow(acov_est_g)) - acov_est_g %*% a1[[g]]
      acov_res[[g]] <- q_1 %*% acov_obs[[g]] %*% t(q_1)

      # correct ACOV.res for type = "cor.bentler" or type = "cor.bollen"
      if (type == "cor.bentler") {
        if (lavmodel@categorical) {
          if (lavmodel@conditional.x ||
            length(unlist(lavmodel@num.idx)) > 0L) {
            lav_msg_stop(gettext(
      "SE for cor.bentler not available (yet) if categorical = TRUE, and
      conditional.x = TRUE OR some endogenous variables are continuous"))
          } else {
            # nothing to do, as we already are in correlation metric
          }
        } else {
          # Ogasawara (2001), eq (13), or
          # Maydeu-Olivares (2017), eq (16)
          cov_1 <- if (lavmodel@conditional.x) {
            sampstat[[g]][["res.cov"]]
          } else {
            sampstat[[g]][["cov"]]
          }
          ss <- 1 / sqrt(diag(cov_1))
          tmp <- lav_mat_vech(tcrossprod(ss))
          g_inv_sqrt <- diag(tmp, nrow = length(tmp))
          if (lavmodel@meanstructure) {
            gg <- lav_mat_bdiag(
              diag(ss, nrow = length(ss)),
              g_inv_sqrt
            )
          } else {
            gg <- g_inv_sqrt
          }
          acov_res[[g]] <- gg %*% acov_res[[g]] %*% gg
        } # continuous
      } else if (type == "cor.bollen") {
        if (lavmodel@categorical) {
          if (lavmodel@conditional.x ||
            length(unlist(lavmodel@num.idx)) > 0L) {
            lav_msg_stop(gettext(
     "SE for cor.bentler not available (yet) if categorical = TRUE, and
     conditional.x = TRUE OR some endogenous variables are continuous"))
          } else {
            # nothing to do, as we already are in correlation metric
          }
        } else {
          # here we use the Maydeu-Olivares (2017) approach, see eq 17
          cov_1 <- if (lavmodel@conditional.x) {
            sampstat[[g]][["res.cov"]]
          } else {
            sampstat[[g]][["cov"]]
          }
          f1 <- lav_deriv_cov2cor_b(cov_1)
          if (lavmodel@meanstructure) {
            ss <- 1 / sqrt(diag(cov_1))
            ff <- lav_mat_bdiag(diag(ss, nrow = length(ss)), f1)
          } else {
            ff <- f1
          }
          acov_res[[g]] <- ff %*% acov_res[[g]] %*% t(ff)
        } # continuous
      } # cor.bollen
    } # z.type = "standardized"
  } # g

  acov_res
}

# return resList with 'se' values for each residual
lav_residuals_se <- function(object, type = "raw", z_type = "standardized",
                             h1_acov = "unstructured",
                             add_class = FALSE, add_labels = FALSE) {
  # slots
  lavdata <- object@Data
  lavmodel <- object@Model
  lavpta <- object@pta

  # return list per group
  se_list <- vector("list", length = lavdata@ngroups)

  # get ACOV per group
  acov_res <- lav_residuals_acov(
    object = object, type = type,
    z_type = z_type, h1_acov = h1_acov
  )

  # labels
  if (add_labels) {
    ov_names <- object@pta$vnames$ov
    # ov_names_res <- object@pta$vnames$ov.nox
    # ov_names_x <- object@pta$vnames$ov.x
  }

  # for each group, compute 'se' values, and fill list
  for (g in seq_len(lavdata@ngroups)) {
    nvar <- object@pta$nvar[[g]] # block or group-based?
    diag_acov <- diag(acov_res[[g]])

    # take care of negative, or non-finite diag.ACOV elements
    diag_acov[!is.finite(diag_acov)] <- NA
    diag_acov[diag_acov < 0] <- NA

    # categorical
    if (lavmodel@categorical) {
      if (lavmodel@conditional.x ||
        length(unlist(lavmodel@num.idx)) > 0L) {
        lav_msg_stop(gettext("not ready yet!"))
      }

      # COR
      nth <- length(lavmodel@th.idx[[g]])
      tmp <- sqrt(diag_acov[-(1:nth)])
      cov_se <- lav_mat_vech_rev(tmp, diagonal = FALSE)

      # MEAN
      mean_se <- rep(as.numeric(NA), nvar)

      # TH
      th_se <- sqrt(diag_acov[1:nth])

      if (add_class) {
        class(cov_se) <- c("lavaan.matrix.symmetric", "matrix")
        class(mean_se) <- c("lavaan.vector", "numeric")
        class(th_se) <- c("lavaan.vector", "numeric")
      }
      if (add_labels) {
        rownames(cov_se) <- colnames(cov_se) <- ov_names[[g]]
        names(mean_se) <- ov_names[[g]]
        names(th_se) <- lavpta$vnames$th.mean[[g]]
      }
      se_list[[g]] <- list(
        cov.se = cov_se, mean.se = mean_se,
        th.se = th_se
      )

      # continuous -- single level
    } else if (lavdata@nlevels == 1L) {
      if (lavmodel@conditional.x) {
        lav_msg_stop(gettext("not ready yet"))
      } else {
        if (lavmodel@meanstructure) {
          tmp <- sqrt(diag_acov[-(1:nvar)])
          cov_se <- lav_mat_vech_rev(tmp)
          mean_se <- sqrt(diag_acov[1:nvar])
          if (add_class) {
            class(cov_se) <- c("lavaan.matrix.symmetric", "matrix")
            class(mean_se) <- c("lavaan.vector", "numeric")
          }
          if (add_labels) {
            rownames(cov_se) <- colnames(cov_se) <- ov_names[[g]]
            names(mean_se) <- ov_names[[g]]
          }
          se_list[[g]] <- list(cov.se = cov_se, mean.se = mean_se)
        } else {
          cov_se <- lav_mat_vech_rev(sqrt(diag_acov))
          if (add_class) {
            class(cov_se) <- c("lavaan.matrix.symmetric", "matrix")
          }
          if (add_labels) {
            rownames(cov_se) <- colnames(cov_se) <- ov_names[[g]]
          }
          se_list[[g]] <- list(cov.se = cov_se)
        }
      }

      # continuous -- multilevel
    } else if (lavdata@nlevels > 1L) {
      lav_msg_stop(gettext("not ready yet"))
    }
  } # g

  se_list
}

# return summary statistics as list per group
lav_residuals_summary <- function(object, type = c("rmr", "srmr", "crmr"),
                                  h1_acov = "unstructured", custom_rmr = NULL,
                                  se = FALSE, zstat = FALSE, pvalue = FALSE,
                                  unbiased = FALSE, unbiased_se = FALSE,
                                  unbiased_ci = FALSE, unbiased_ci_level = 0.90,
                                  unbiased_zstat = FALSE,
                                  unbiased_test_val = 0.05,
                                  unbiased_pvalue = FALSE,
                                  add_class = FALSE) {
  # arguments
  if (length(custom_rmr)) {
    if (!is.list(custom_rmr)) lav_msg_stop(gettext("custom.rmr must be a list"))
    ## Each custom (S/C)RMR must have a unique name
    custom_names <- names(custom_rmr)
    if (is.null(custom_names)) lav_msg_stop(gettext(
      "custom.rmr list must have names"))
    if (length(unique(custom_names)) < length(custom_rmr)) {
      lav_msg_stop(gettext(
        "custom.rmr must have a unique name for each summary"))
    }
    ## Each list must contain a list consisting of $cov and/or
    #                                                  $mean (no $th yet)
    for (i in seq_along(custom_rmr)) {
      if (!is.list(custom_rmr[[i]])) {
        lav_msg_stop(gettext("Each element in custom.rmr must be a list"))
      }
      if (is.null(names(custom_rmr[[i]]))) {
        lav_msg_stop(gettext("The list in custom.rmr must have names"))
      }
      if (!all(names(custom_rmr[[i]]) %in% c("cov", "mean"))) {
        lav_msg_stop(gettext(
          'Elements in custom.rmr must be named "cov" and/or "mean"'))
      }
      ## below, verify dimensions match rmsList.g
    }
    # FIXME: blocks can have unique models, need another layer of lists
    #       between custom summaries and moments
  } else {
    custom_names <- NULL
  }

  if (pvalue) {
    zstat <- TRUE
  }
  if (zstat) {
    se <- TRUE
  }
  if (unbiased_pvalue) {
    unbiased_zstat <- TRUE
  }
  if (unbiased_zstat) {
    unbiased_se <- TRUE
  }

  if (!all(type %in% c(
    "rmr", "srmr", "crmr",
    "raw", "cor.bentler", "cor.bollen"
  ))) {
    lav_msg_stop(gettext("unknown type:"), dQuote(type))
  }

  # change type name
  idx <- which(type == "raw")
  if (length(idx) > 0L) {
    type[idx] <- "rmr"
  }
  idx <- which(type == "cor.bentler")
  if (length(idx) > 0L) {
    type[idx] <- "srmr"
  }
  idx <- which(type == "cor.bollen")
  if (length(idx) > 0L) {
    type[idx] <- "crmr"
  }

  # slots
  lavdata <- object@Data
  lavmodel <- object@Model

  # fixed.x/conditional.x
  # fixed_x <- lavmodel@fixed.x
  conditional_x <- lavmodel@conditional.x

  # rmr_flag <- srmr_flag <- crmr_flag <- FALSE
  if ("rmr" %in% type || "raw" %in% type) {
    # FIXME: recursive call to lav_residuals() is summary = TRUE!!
    rmr_list <- lav_residuals(object = object, type = "raw")
    if (se || unbiased) {
      rmr_list_se <- lav_residuals_acov(
        object = object, type = "raw",
        z_type = "standardized",
        h1_acov = "unstructured"
      )
    }
  }
  if ("srmr" %in% type || "cor.bentler" %in% type || "cor" %in% type) {
    srmr_list <- lav_residuals(object = object, type = "cor.bentler")
    if (se || unbiased) {
      srmr_list_se <- lav_residuals_acov(
        object = object,
        type = "cor.bentler",
        z_type = "standardized",
        h1_acov = "unstructured"
      )
    }
  }
  if ("crmr" %in% type || "cor.bollen" %in% type) {
    crmr_list <- lav_residuals(object = object, type = "cor.bollen")
    if (se || unbiased) {
      crmr_list_se <- lav_residuals_acov(
        object = object,
        type = "cor.bollen",
        z_type = "standardized",
        h1_acov = "unstructured"
      )
    }
  }

  # return list per group
  sum_stat <- vector("list", length = lavdata@ngroups)

  # for each group, compute ACOV
  for (g in seq_len(lavdata@ngroups)) {
    nvar <- object@pta$nvar[[g]] # block or group-based?

    # categorical single level
    if (lavdata@nlevels == 1L && lavmodel@categorical) {
      if ((se || unbiased) && (conditional_x ||
        length(unlist(lavmodel@num.idx)) > 0L)) {
        lav_msg_stop(gettext("not ready yet"))
      } else {
        # remove fixed.x elements:
        # seems like a good idea, but nobody likes it
        # nvar.x <- pstar.x <- 0L
        # if(lavmodel@fixed.x) {
        #     nvar.x <- lavmodel@nexo[g]
        #     pstar.x <- nvar.x * (nvar.x - 1) / 2 # note '-'
        # }

        out <- vector("list", length(type))
        names(out) <- type

        for (typ in seq_along(type)) {
          if (type[typ] == "rmr") {
            rms_list_g <- rmr_list[[g]]
            if (se || unbiased) {
              rms_list_se_g <- rmr_list_se[[g]]
            }
          } else if (type[typ] == "srmr") {
            rms_list_g <- srmr_list[[g]]
            if (se || unbiased) {
              rms_list_se_g <- srmr_list_se[[g]]
            }
          } else if (type[typ] == "crmr") {
            rms_list_g <- crmr_list[[g]]
            if (se || unbiased) {
              rms_list_se_g <- crmr_list_se[[g]]
            }
          }

          # COR
          nth <- length(lavmodel@th.idx[[g]])
          if (conditional_x) {
            stats <- lav_mat_vech(rms_list_g[["res.cov"]],
              diagonal = FALSE
            )
          } else {
            stats <- lav_mat_vech(rms_list_g[["cov"]],
              diagonal = FALSE
            )
          }

          # should pstar be p*(p+1)/2 or p*(p-1)/2
          # we use the first for SRMR and the latter for CRMR
          if (type[typ] == "crmr") {
            pstar <- length(stats)
          } else {
            pstar <- length(stats) + nvar
          }
          acov <- NULL
          if (se || unbiased) {
            acov <- rms_list_se_g[-seq_len(nth),
              -seq_len(nth),
              drop = FALSE
            ]
          }
          rms_cor <- lav_residuals_summary_rms(
            stats = stats,
            acov = acov, se = se, zstat = zstat, pvalue = pvalue,
            unbiased = unbiased, unbiased_se = unbiased_se,
            unbiased_ci = unbiased_ci,
            unbiased_ci_level = unbiased_ci_level,
            unbiased_zstat = unbiased_zstat,
            unbiased_test_val = unbiased_test_val,
            unbiased_pvalue = unbiased_pvalue,
            pstar = pstar, type = type[typ]
          )


          # THRESHOLDS
          if (conditional_x) {
            stats <- rms_list_g[["res.th"]]
          } else {
            stats <- rms_list_g[["th"]]
          }
          pstar <- length(stats)
          acov <- NULL
          if (se || unbiased) {
            acov <- rms_list_se_g[seq_len(nth),
              seq_len(nth),
              drop = FALSE
            ]
          }
          rms_th <- lav_residuals_summary_rms(
            stats = stats,
            acov = acov, se = se, zstat = zstat, pvalue = pvalue,
            unbiased = unbiased, unbiased_se = unbiased_se,
            unbiased_ci = unbiased_ci,
            unbiased_ci_level = unbiased_ci_level,
            unbiased_zstat = unbiased_zstat,
            unbiased_test_val = unbiased_test_val,
            unbiased_pvalue = unbiased_pvalue,
            pstar = pstar, type = type[typ]
          )

          # MEAN
          # STATS <- rmsList.g[["mean"]]
          stats <- numeric(0L)
          pstar <- length(stats)
          acov <- NULL
          if (se || unbiased) {
            # TODO: extract from rmsList.se.g
          }
          rms_mean <- lav_residuals_summary_rms(
            stats = stats,
            acov = acov, se = se, zstat = zstat, pvalue = pvalue,
            unbiased = unbiased, unbiased_se = unbiased_se,
            unbiased_ci = unbiased_ci,
            unbiased_ci_level = unbiased_ci_level,
            unbiased_zstat = unbiased_zstat,
            unbiased_test_val = unbiased_test_val,
            unbiased_pvalue = unbiased_pvalue,
            pstar = pstar, type = type[typ]
          )

          # VAR (not ready yet)
          # STATS <- diag(rmsList.g[["cov"]])[lavmodel@num.idx[[g]]]
          stats <- numeric(0L)
          pstar <- length(stats)
          acov <- NULL
          if (se || unbiased) {
            # TODO: extract from rmsList.se.g
          }
          # rms_var <- lav_residuals_summary_rms(
          #   STATS = stats,
          #   ACOV = acov, se = se, zstat = zstat, pvalue = pvalue,
          #   unbiased = unbiased, unbiased_se = unbiased_se,
          #   unbiased_ci = unbiased_ci,
          #   unbiased_ci_level = unbiased_ci_level,
          #   unbiased_zstat = unbiased_zstat,
          #   unbiased_test_val = unbiased_test_val,
          #   unbiased_pvalue = unbiased_pvalue,
          #   pstar = pstar, type = type[typ]
          # )

          # TOTAL -- FIXME: for conditional.x ....
          if (conditional_x) {
            stats <- c(
              lav_mat_vech(rms_list_g[["res.cov"]],
                diagonal = FALSE
              ),
              rms_list_g[["res.th"]]
            )
          } else {
            stats <- c(
              lav_mat_vech(rms_list_g[["cov"]],
                diagonal = FALSE
              ),
              rms_list_g[["th"]]
            )
            # rmsList.g[["mean"]],
            # diag(rmsList.g[["cov"]])[lavmodel@num.idx[[g]]])
          }

          # should pstar be p*(p+1)/2 or p*(p-1)/2 for COV/COR?
          # we use the first for SRMR and the latter for CRMR
          if (type[typ] == "crmr") {
            pstar <- length(stats)
          } else {
            pstar <- length(stats) + nvar
          }

          # if(lavmodel@fixed.x) {
          #    pstar <- pstar - pstar.x
          # }

          acov <- NULL
          if (se || unbiased) {
            acov <- rms_list_se_g
          }
          rms_total <- lav_residuals_summary_rms(
            stats = stats,
            acov = acov, se = se, zstat = zstat, pvalue = pvalue,
            unbiased = unbiased, unbiased_se = unbiased_se,
            unbiased_ci = unbiased_ci,
            unbiased_ci_level = unbiased_ci_level,
            unbiased_zstat = unbiased_zstat,
            unbiased_test_val = unbiased_test_val,
            unbiased_pvalue = unbiased_pvalue,
            pstar = pstar, type = type[typ]
          )

          table_1 <- as.data.frame(cbind(
            rms_cor,
            rms_th,
            # RMS.MEAN,
            # RMS.VAR,
            rms_total
          ))
          # colnames(TABLE) <- c("cor", "thresholds", "mean",
          #                     "var", "total")
          colnames(table_1) <- c("cor", "thresholds", "total")
          if (add_class) {
            class(table_1) <- c("lavaan.data.frame", "data.frame")
          }
          out[[typ]] <- table_1
        } # type
      } # not conditional.x or mixed cat/con

      # continuous -- single level
    } else if (lavdata@nlevels == 1L) {
      if ((se || unbiased) && conditional_x) {
        lav_msg_stop(gettext("not ready yet"))
      } else {
        # nvar.x <- pstar.x <- 0L
        # if(lavmodel@fixed.x) {
        #    nvar.x <- lavmodel@nexo[g]
        #    pstar.x <- nvar.x * (nvar.x + 1) / 2
        # }

        out <- vector("list", length(type))
        names(out) <- type

        for (typ in seq_along(type)) {
          if (type[typ] == "rmr") {
            rms_list_g <- rmr_list[[g]]
            if (se || unbiased) {
              rms_list_se_g <- rmr_list_se[[g]]
            }
          } else if (type[typ] == "srmr") {
            rms_list_g <- srmr_list[[g]]
            if (se || unbiased) {
              rms_list_se_g <- srmr_list_se[[g]]
            }
          } else if (type[typ] == "crmr") {
            rms_list_g <- crmr_list[[g]]
            if (se || unbiased) {
              rms_list_se_g <- crmr_list_se[[g]]
            }
          }

          # COV
          if (conditional_x) {
            stats <- lav_mat_vech(rms_list_g[["res.cov"]])
          } else {
            stats <- lav_mat_vech(rms_list_g[["cov"]])
          }
          # pstar <- ( length(STATS) - pstar.x )
          pstar <- length(stats)
          if (type[typ] == "crmr") {
            # pstar <- pstar - ( nvar - nvar.x )
            if (conditional_x) {
              pstar <- pstar - nrow(rms_list_g[["res.cov"]])
            } else {
              pstar <- pstar - nvar
            }
          }

          acov <- NULL
          if (se || unbiased) {
            acov <- if (lavmodel@meanstructure) {
              rms_list_se_g[-seq_len(nvar),
                -seq_len(nvar),
                drop = FALSE
              ]
            } else {
              rms_list_se_g
            }
          }
          rms_cov <- lav_residuals_summary_rms(
            stats = stats,
            acov = acov, se = se, zstat = zstat, pvalue = pvalue,
            unbiased = unbiased, unbiased_se = unbiased_se,
            unbiased_ci = unbiased_ci,
            unbiased_ci_level = unbiased_ci_level,
            unbiased_zstat = unbiased_zstat,
            unbiased_test_val = unbiased_test_val,
            unbiased_pvalue = unbiased_pvalue,
            pstar = pstar, type = type[typ]
          )

          # MEAN
          if (lavmodel@meanstructure) {
            if (conditional_x) {
              stats <- rms_list_g[["res.int"]]
            } else {
              stats <- rms_list_g[["mean"]]
            }
            # pstar <- ( length(STATS) - nvar.x )
            pstar <- length(stats)
            acov <- NULL
            if (se || unbiased) {
              acov <- rms_list_se_g[seq_len(nvar),
                seq_len(nvar),
                drop = FALSE
              ]
            }
            rms_mean <- lav_residuals_summary_rms(
              stats = stats,
              acov = acov,
              se = se, zstat = zstat, pvalue = pvalue,
              unbiased = unbiased, unbiased_se = unbiased_se,
              unbiased_ci = unbiased_ci,
              unbiased_ci_level = unbiased_ci_level,
              unbiased_zstat = unbiased_zstat,
              unbiased_test_val = unbiased_test_val,
              unbiased_pvalue = unbiased_pvalue,
              pstar = pstar, type = type[typ]
            )
          }

          # TOTAL
          if (lavmodel@meanstructure) {
            if (conditional_x) {
              stats <- c(
                rms_list_g[["res.int"]],
                lav_mat_vech(rms_list_g[["res.cov"]])
              )
            } else {
              stats <- c(
                rms_list_g[["mean"]],
                lav_mat_vech(rms_list_g[["cov"]])
              )
            }
            # pstar <- ( length(STATS) - ( pstar.x + nvar.x) )
            pstar <- length(stats)
            if (type[typ] == "crmr") {
              # pstar <- pstar - ( nvar - nvar.x )
              pstar <- pstar - nvar
            }
            acov <- NULL
            if (se || unbiased) {
              acov <- rms_list_se_g
            }
            rms_total <- lav_residuals_summary_rms(
              stats = stats,
              acov = acov,
              se = se, zstat = zstat, pvalue = pvalue,
              unbiased = unbiased, unbiased_se = unbiased_se,
              unbiased_ci = unbiased_ci,
              unbiased_ci_level = unbiased_ci_level,
              unbiased_zstat = unbiased_zstat,
              unbiased_test_val = unbiased_test_val,
              unbiased_pvalue = unbiased_pvalue,
              pstar = pstar, type = type[typ]
            )
          }

          # CUSTOM
          if (length(custom_rmr)) {
            if (lavmodel@fixed.x && !lavmodel@conditional.x) {
              ## save exogenous-variable indices, use to remove or set
              ## FALSE any moments that cannot have nonzero residuals
              x_idx <- which(rownames(rms_list_g$cov) %in%
                                          object@Data@ov.names.x[[g]])
            }

            rms_custom_list <- vector("list", length(custom_names))

            for (cus in custom_names) {
              ## in case there is no meanstructure
              stats <- NULL
              acov_idx <- NULL

              # MEANS?
              if (lavmodel@meanstructure) {
                if ("mean" %in% names(custom_rmr[[cus]])) {
                  ## if logical, save numeric indices
                  if (is.logical(custom_rmr[[cus]]$mean)) {
                    ## check length
                    if (length(custom_rmr[[cus]]$mean) !=
                                          length(rms_list_g[["mean"]])) {
                      lav_msg_stop(gettextf("length(custom.rmr$%s$mean) must
                            match length(lavResiduals(fit)$mean)", cus))
                    }
                    acov_idx <- which(custom_rmr[[cus]]$mean)
                    if (lavmodel@fixed.x && !lavmodel@conditional.x) {
                      acov_idx[x_idx] <- FALSE
                    }
                  } else if (!is.numeric(custom_rmr[[cus]]$mean)) {
                    lav_msg_stop(gettextf("custom.rmr$%s$mean must contain
                                  logical or numeric indices.", cus))
                  } else {
                    acov_idx <- custom_rmr[[cus]]$mean
                    if (lavmodel@fixed.x && !lavmodel@conditional.x) {
                      acov_idx <- setdiff(acov_idx, x_idx)
                    }
                    acov_idx <- acov_idx[!is.na(acov_idx)] # necessary?
                    if (max(acov_idx) > length(rms_list_g[["mean"]])) {
                      lav_msg_stop(gettextf(
                        "custom.rmr$%1$s$mean[%2$s] is an out-of-bounds index",
                        cus, which.max(acov_idx))
                      )
                    }
                  }
                  stats <- rms_list_g[["mean"]][acov_idx]
                }
              }
              # (CO)VARIANCES?
              if ("cov" %in% names(custom_rmr[[cus]])) {
                ## if numeric, create a logical matrix to obtain
                ## ACOV.idx and check for x.idx
                if (is.numeric(custom_rmr[[cus]]$cov)) {
                  cus_cov <- rms_list_g[["cov"]] == "start with all FALSE"
                  ## matrix of row/column indices?
                  if (length(dim(custom_rmr[[cus]]$cov))) {
                    if (max(custom_rmr[[cus]]$cov[, 1:2] >
                                              nrow(rms_list_g[["cov"]]))) {
                      lav_msg_stop(gettextf(
                        "numeric indices in custom.rmr$%1$s$cov cannot
                        exceed %2$s", cus, nrow(rms_list_g[["cov"]])))
                    }
                    for (rr in seq_len(nrow(custom_rmr[[cus]]$cov))) {
                      cus_cov[
                        custom_rmr[[cus]]$cov[rr, 1],
                        custom_rmr[[cus]]$cov[rr, 2]
                      ] <- TRUE
                    }
                  } else {
                    ## numeric-vector indices
                    if (max(custom_rmr[[cus]]$cov >
                                          length(rms_list_g[["cov"]]))) {
                      lav_msg_stop(gettextf(
                        "numeric indices in custom.rmr$%1$s$cov cannot
                        exceed %2$s", cus, length(rms_list_g[["cov"]])))
                    }
                    cus_cov[custom_rmr[[cus]]$cov] <- TRUE
                  }

                  ## numeric indices no longer needed, use logical
                  custom_rmr[[cus]]$cov <- cus_cov
                } else if (!is.logical(custom_rmr[[cus]]$cov)) {
                  lav_msg_stop(gettextf(
                    "custom.rmr$%s$cov must be a logical square matrix or a
                    numeric matrix of (row/column) indices.", cus))
                }

                ## check dimensions
                if (!all(dim(custom_rmr[[cus]]$cov) ==
                                                  dim(rms_list_g[["cov"]]))) {
                  lav_msg_stop(gettextf(
                    "dim(custom.rmr$%s$cov) must match
                    dim(lavResiduals(fit)$cov)", cus))
                }
                ## users can specify upper.tri or lower.tri indices
                custom_rmr[[cus]]$cov <-
                                custom_rmr[[cus]]$cov | t(custom_rmr[[cus]]$cov)
                ## but ACOV refers to lower.tri indices
                custom_rmr[[cus]]$cov[upper.tri(custom_rmr[[cus]]$cov)] <- FALSE
                ## diagonal relevant?
                if (type[typ] == "crmr") diag(custom_rmr[[cus]]$cov) <- FALSE
                ## extract lower.tri indices
                vech_idx <- which(lav_mat_vech(custom_rmr[[cus]]$cov))

                ## add residuals to STATS, indices to ACOV.idx
                stats <- c(stats,
                         lav_mat_vech(rms_list_g[["cov"]])[vech_idx])
                acov_idx <- c(acov_idx, vech_idx)
              }


              ## count residuals in summary (x.idx already removed)
              pstar <- length(stats)

              acov <- NULL
              if (se || unbiased) {
                acov <- rms_list_se_g[acov_idx, acov_idx, drop = FALSE]
              }
              rms_custom_list[[cus]] <- lav_residuals_summary_rms(
                stats = stats,
                acov = acov,
                se = se, zstat = zstat, pvalue = pvalue,
                unbiased = unbiased, unbiased_se = unbiased_se,
                unbiased_ci = unbiased_ci,
                unbiased_ci_level = unbiased_ci_level,
                unbiased_zstat = unbiased_zstat,
                unbiased_test_val = unbiased_test_val,
                unbiased_pvalue = unbiased_pvalue,
                pstar = pstar, type = type[typ]
              )

              # FIXME: update for categorical
            } # cus
            rms_custom <- do.call(rbind, rms_custom_list)
          } else {
            rms_custom <- NULL
          }


          if (lavmodel@meanstructure) {
            table_1 <- as.data.frame(cbind(
              rms_cov,
              rms_mean,
              rms_total,
              rms_custom
            ))
            colnames(table_1) <- c(
              "cov", "mean", "total",
              custom_names
            )
          } else {
            table_1 <- as.data.frame(cbind(rms_cov, rms_custom))
            colnames(table_1) <- c("cov", custom_names)
          }
          if (add_class) {
            class(table_1) <- c("lavaan.data.frame", "data.frame")
          }
          out[[typ]] <- table_1
        } # type
      } # continuous, single-level, unconditional

      # continuous -- multilevel
    } else if (lavdata@nlevels > 1L) {
      lav_msg_stop(gettext("not ready yet"))
    }

    sum_stat[[g]] <- out
  } # g

  sum_stat
}


lav_residuals_summary_rms <- function(stats = NULL, acov = NULL,
                                      se = FALSE,
                                      level = 0.90,
                                      zstat = FALSE, pvalue = FALSE,
                                      unbiased = FALSE,
                                      unbiased_se = FALSE,
                                      unbiased_ci = FALSE,
                                      unbiased_ci_level = 0.90,
                                      unbiased_zstat = FALSE,
                                      unbiased_test_val = 0.05,
                                      unbiased_pvalue = FALSE,
                                      pstar = 0, type = "rms") {
  out <- vector("list", length = 0L)


  # covariance matrix
  if (length(stats) > 0L) {
    rms <- sqrt(sum(stats * stats) / pstar)
  } else {
    rms <- 0
    se <- unbiased <- zstat <- FALSE
  }

  # default is NULL
  rms_se <- rms_z <- rms_pvalue <- NULL
  urms <- urms_se <- urms_z <- urms_pvalue <- NULL
  urms_ci_lower <- urms_ci_upper <- NULL
  if (!unbiased_zstat) {
    unbiased_test_val <- NULL
  }

  if (se || unbiased) {
    tr2 <- sum(diag(acov %*% acov))
    tr1 <- sum(diag(acov))
    if (se) {
      rms_avar <- tr2 / (tr1 * 2 * pstar)
      if (!is.finite(rms_avar) || rms_avar < .Machine$double.eps) {
        rms_se <- as.numeric(NA)
      } else {
        rms_se <- sqrt(rms_avar)
      }
    }
  }

  if (zstat) {
    e_rms <- (sqrt(tr1 / pstar) * (4 * tr1 * tr1 - tr2) / (4 * tr1 * tr1))
    rms_z <- max((rms - e_rms), 0) / rms_se
    if (pvalue) {
      rms_pvalue <- 1 - pnorm(rms_z)
    }
  }

  if (unbiased) {
    t_cov <- as.numeric(crossprod(stats))
    e_ve <- as.numeric(t(stats) %*% acov %*% stats)
    k_cov <- 1 - (tr2 + 2 * e_ve) / (4 * t_cov * t_cov)
    urms <- (1 / k_cov * sqrt(max((t_cov - tr1), 0) / pstar))
    if (unbiased_se) {
      urms_avar <- (1 / (k_cov * k_cov) * (tr2 + 2 * e_ve) /
                                               (2 * pstar * t_cov))
      if (!is.finite(urms_avar) || urms_avar < .Machine$double.eps) {
        urms_se <- as.numeric(NA)
      } else {
        urms_se <- sqrt(urms_avar)
      }
      if (unbiased_ci) {
        a <- (1 - unbiased_ci_level) / 2
        a <- c(a, 1 - a)
        fac <- stats::qnorm(a)
        urms_ci_lower <- urms + urms_se * fac[1]
        urms_ci_upper <- urms + urms_se * fac[2]
      }
      if (unbiased_zstat) {
        urms_z <- (urms - unbiased_test_val) / urms_se
        if (unbiased_pvalue) {
          urms_pvalue <- 1 - pnorm(urms_z)
        }
      }
    }
  }

  # labels
  if (type == "rmr") {
    out <- list(
      rmr = rms, rmr.se = rms_se,
      rmr.exactfit.z = rms_z, rmr.exactfit.pvalue = rms_pvalue,
      urmr = urms, urmr.se = urms_se,
      urmr.ci.lower = urms_ci_lower,
      urmr.ci.upper = urms_ci_upper,
      urmr.closefit.h0.value = unbiased_test_val,
      urmr.closefit.z = urms_z,
      urmr.closefit.pvalue = urms_pvalue
    )
  } else if (type == "srmr") {
    out <- list(
      srmr = rms, srmr.se = rms_se,
      srmr.exactfit.z = rms_z, srmr.exactfit.pvalue = rms_pvalue,
      usrmr = urms, usrmr.se = urms_se,
      usrmr.ci.lower = urms_ci_lower,
      usrmr.ci.upper = urms_ci_upper,
      usrmr.closefit.h0.value = unbiased_test_val,
      usrmr.closefit.z = urms_z,
      usrmr.closefit.pvalue = urms_pvalue
    )
  } else if (type == "crmr") {
    out <- list(
      crmr = rms, crmr.se = rms_se,
      crmr.exactfit.z = rms_z, crmr.exactfit.pvalue = rms_pvalue,
      ucrmr = urms, ucrmr.se = urms_se,
      ucrmr.ci.lower = urms_ci_lower,
      ucrmr.ci.upper = urms_ci_upper,
      ucrmr.closefit.h0.value = unbiased_test_val,
      ucrmr.closefit.z = urms_z,
      ucrmr.closefit.pvalue = urms_pvalue
    )
  }

  unlist(out)
}

# generate summary statistics for the residuals
lav_residuals_summary_old <- function(res_list = NULL,
                                      add_class = FALSE, add_labels = FALSE) {
  # per block
  nblocks <- length(res_list)

  for (b in seq_len(nblocks)) {
    # create new list, including with summary statistics interleaved
    x <- vector("list", length = 0L)
    nel <- length(res_list[[b]])
    names_1 <- names(res_list[[b]])

    for (el in seq_len(nel)) {
      el_1 <- res_list[[b]][[el]]
      if (!is.null(names_1)) {
        name <- names_1[el]
      }

      if (is.character(el_1)) {
        new_x <- list(el_1)
        if (add_labels) {
          names(new_x) <- "type"
        }
        x <- c(x, new_x)
      } else if (is.matrix(el_1) && isSymmetric(el_1)) {
        tmp <- na.omit(lav_mat_vech(el_1))
        rms <- sqrt(sum(tmp * tmp) / length(tmp))
        mabs <- mean(abs(tmp))
        tmp2 <- na.omit(lav_mat_vech(el_1, diagonal = FALSE))
        rms_nodiag <- sqrt(sum(tmp2 * tmp2) / length(tmp2))
        mabs_nodiag <- mean(abs(tmp2))
        cov_summary <- c(rms, rms_nodiag, mabs, mabs_nodiag)
        if (add_labels) {
          names(cov_summary) <-
            c("rms", "rms.nodiag", "mabs", "mabs.nodiag")
        }
        if (add_class) {
          class(cov_summary) <- c("lavaan.vector", "numeric")
        }
        new_x <- list(el_1, cov_summary)
        if (add_labels && !is.null(names_1)) {
          names(new_x) <- c(name, paste0(name, ".summary"))
        }
        x <- c(x, new_x)
      } else {
        tmp <- na.omit(el_1)
        rms <- sqrt(sum(tmp * tmp) / length(tmp))
        mabs <- mean(abs(tmp))
        mean_summary <- c(rms, mabs)
        if (add_labels) {
          names(mean_summary) <- c("rms", "mabs")
        }
        if (add_class) {
          class(mean_summary) <- c("lavaan.vector", "numeric")
        }
        new_x <- list(el_1, mean_summary)
        if (add_labels && !is.null(names_1)) {
          names(new_x) <- c(name, paste0(name, ".summary"))
        }
        x <- c(x, new_x)
      }
    } # nel

    # fill in block including summary statistics
    res_list[[b]] <- x
  } # nblocks

  res_list
}


# x is a list with sample statistics (eg output of lavInspect(fit, "sampstat")
# y is another (possibly the same) list with sample statistics
#
# to avoid many 'NAs', we set the scale-factor to 1
# if the to-be-scaled value is < 1e-05 (in absolute value)
lav_residuals_rescale <- function(x, diag_cov = NULL, diag_cov2 = NULL) {
  if (is.null(diag_cov2)) {
    diag_cov2 <- diag_cov
  }

  # make sure we can take the sqrt and invert
  diag_cov[!is.finite(diag_cov)] <- NA
  diag_cov[diag_cov < .Machine$double.eps] <- NA
  scale_cov <- tcrossprod(1 / sqrt(diag_cov))

  # for the mean, we use diag.cov2
  diag_cov2[!is.finite(diag_cov2)] <- NA
  diag_cov2[diag_cov2 < .Machine$double.eps] <- NA
  scale_mean <- 1 / sqrt(diag_cov2)

  # rescale cov
  if (!is.null(x[["cov"]])) {
    # catch (near) zero elements in x$cov
    near_zero_idx <- which(abs(x[["cov"]]) < 1e-05)
    scale_cov[near_zero_idx] <- 1
    x[["cov"]][] <- x[["cov"]] * scale_cov
  }
  if (!is.null(x[["res.cov"]])) {
    # catch (near) zero elements in x$res.cov
    near_zero_idx <- which(abs(x[["res.cov"]]) < 1e-05)
    scale_cov[near_zero_idx] <- 1
    x[["res.cov"]][] <- x[["res.cov"]] * scale_cov
  }

  # rescale int/mean
  if (!is.null(x[["res.int"]])) {
    # catch (near) zero elements in x$res.int
    near_zero_idx <- which(abs(x[["res.int"]]) < 1e-05)
    scale_mean[near_zero_idx] <- 1
    x[["res.int"]] <- x[["res.int"]] * scale_mean
  }

  if (!is.null(x[["mean"]])) {
    # catch (near) zero elements in x$mean
    near_zero_idx <- which(abs(x[["mean"]]) < 1e-05)
    scale_mean[near_zero_idx] <- 1
    x[["mean"]] <- x[["mean"]] * scale_mean
  }

  # FIXME: do something sensible for th, slopes, ...

  x
}
