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

# - change 0.7-1: refactor (issue #291): ACOV.obs (= the expensive ACOV of
#                  the observed h1 sample statistics) is computed at most once
#                  per call and reused for both the residual SEs and the
#                  summary statistics SEs (see lav_residuals_acov_obs())
# - change 0.7-1: the (undocumented, never-functional) custom.rmr argument
#                  has been removed
# - change 0.7-1: raw (rmr) residual SEs/z-statistics/usrmr are now available
#                  for the mixed (continuous + ordinal) case, as long as
#                  conditional.x = FALSE; the standardized (cor.bentler/
#                  cor.bollen) SEs for mixed data are still not ready

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
lavResiduals <- function(object, type = "cor.bentler",                    # nolint start
                         se = FALSE, zstat = TRUE, summary = TRUE,
                         h1.acov = "unstructured",
                         add.type = TRUE, add.labels = TRUE, add.class = TRUE,
                         drop.list.single.group = TRUE,
                         maximum.number = length(res_vech),
                         output = "list") {                                # nolint end
  out <- lav_residuals(
    object = object, type = type, h1 = TRUE,
    se = se, zstat = zstat,
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

# summary_options with everything switched off; used for the settings where
# we only report point estimates (no SE, no z, no unbiased statistics)
lav_residuals_summary_options_off <- function() {
  list(
    se = FALSE, zstat = FALSE,
    pvalue = FALSE, unbiased = FALSE,
    unbiased_se = FALSE,
    unbiased_ci = FALSE, unbiased_ci_level = 0.90,
    unbiased_zstat = FALSE, unbiased_test_val = 0.05,
    unbiased_pvalue = FALSE
  )
}

# main function
lav_residuals <- function(object, type = "raw", h1 = TRUE,
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
    # The standardized (cor.bentler/cor.bollen) SEs require a cov->cor
    # jacobian in the mixed metric that is not ready yet, and conditional.x
    # adds regression slopes to the moment vector (also not ready). The raw
    # residual SEs, however, are metric-agnostic (Q-projection ACOV), so we
    # support them for the unconditional mixed (continuous + ordinal) case.
    mixed <- length(unlist(lavmodel@num.idx)) > 0L
    raw_ok <- (type == "raw") && !lavmodel@conditional.x
    if ((lavmodel@conditional.x || mixed) && !raw_ok) {
      zstat <- se <- FALSE
      summary <- FALSE
      summary_options_1 <- lav_residuals_summary_options_off()
    }
  }

  # change options if conditional.x (for now)
  if (!lavmodel@categorical && lavmodel@conditional.x) {
    zstat <- se <- FALSE
    summary <- FALSE
    summary_options_1 <- lav_residuals_summary_options_off()
  }

  # change options if estimator = "IV" (for now)
  if (lavmodel@estimator == "IV") {
    zstat <- se <- FALSE
    summary <- FALSE
    summary_options_1 <- lav_residuals_summary_options_off()
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
  #
  # ACOV.obs (the ACOV of the observed h1 sample statistics) is the expensive
  # building block underlying both the residual SEs and the summary statistics
  # SEs. We compute it (at most) once here and pass it down to avoid computing
  # it twice (see issue #291).
  acov_obs <- NULL
  if (se || zstat) {
    acov_obs <- lav_residuals_acov_obs(object, h1_acov = h1_acov)
    se_list <- lav_residuals_se(object,
      type = type, z_type = "standardized",
      h1_acov = h1_acov,
      add_class = add_class, add_labels = add_labels,
      acov_obs = acov_obs
    )
  } else if (type %in% c("normalized", "standardized", "standardized.mplus")) {
    acov_obs <- lav_residuals_acov_obs(object, h1_acov = h1_acov)
    se_list <- lav_residuals_se(object,
      type = "raw", z_type = type,
      h1_acov = h1_acov,
      add_class = add_class, add_labels = add_labels,
      acov_obs = acov_obs
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
    # lav_residuals_summary() computes its ACOV.obs using h1_acov =
    # "unstructured"; reuse the one we already have if it matches (issue #291)
    summary_acov_obs <- if (identical(h1_acov, "unstructured")) {
      acov_obs
    } else {
      NULL
    }
    args <- c(
      list(
        object = object, type = type, h1_acov = h1_acov,
        add_class = add_class,
        acov_obs = summary_acov_obs
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

# compute ACOV for the observed h1 sample statistics (ACOV == Gamma/N)
#
# This is the (often expensive) building block underlying all residual
# standard errors. It only depends on 'object' and 'h1_acov' (NOT on 'type'
# or 'z_type'), so it can be computed once and reused across the various
# lav_residuals_acov() calls (see issue #291).
lav_residuals_acov_obs <- function(object, h1_acov = "unstructured") {
  lavsamplestats <- object@SampleStats

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

  acov_obs
}

# return ACOV as list per group
#
# 'acov_obs' (the ACOV of the observed h1 sample statistics) may be passed in
# by the caller to avoid recomputing it; if NULL, it is computed here.
lav_residuals_acov <- function(object, type = "raw", z_type = "standardized",
                               h1_acov = "unstructured", acov_obs = NULL) {
  # check type
  if (z_type %in% c("normalized", "standardized.mplus") && type != "raw") {
    lav_msg_stop(gettextf(
      "z.type = %1$s can only be used with type = %2$s",
      dQuote(z_type), dQuote("raw")))
  }

  # slots
  lavdata <- object@Data
  lavmodel <- object@Model

  # return list per group
  acov_res <- vector("list", length = lavdata@ngroups)

  # compute ACOV for observed h1 sample statistics (reuse if provided)
  if (is.null(acov_obs)) {
    acov_obs <- lav_residuals_acov_obs(object, h1_acov = h1_acov)
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

# index map of the (categorical) h1 sample-statistics vector, in the order
# produced by lav_samp_wls_obs() for the unconditional categorical case:
#   [ th-block (thresholds + continuous means, interleaved by variable),
#     continuous variances (num.idx),
#     covariances (vech, no diagonal) ]
#
# Returns, as positions into that vector: $th (thresholds), $mean (continuous
# means), $var (continuous variances), $cov (off-diagonal (co)variances), plus
# $num_idx (continuous variable indices) and $nvar. Only supported for the
# unconditional (conditional.x = FALSE) case.
lav_residuals_cat_idx <- function(lavmodel, g = 1L, nvar = NULL) {
  th_idx <- lavmodel@th.idx[[g]]
  nth_block <- length(th_idx)
  num_idx <- lavmodel@num.idx[[g]]
  n_num <- length(num_idx)
  if (is.null(nvar)) {
    nvar <- length(unique(th_idx[th_idx != 0])) + n_num
  }
  n_offdiag <- nvar * (nvar - 1) / 2

  list(
    th = which(th_idx != 0L), # threshold positions in the th-block
    mean = which(th_idx == 0L), # continuous-mean positions in the th-block
    var = nth_block + seq_len(n_num), # continuous variance positions
    cov = nth_block + n_num + seq_len(n_offdiag), # off-diagonal cov positions
    num_idx = num_idx,
    nvar = nvar
  )
}

# return resList with 'se' values for each residual
lav_residuals_se <- function(object, type = "raw", z_type = "standardized",
                             h1_acov = "unstructured",
                             add_class = FALSE, add_labels = FALSE,
                             acov_obs = NULL) {
  # slots
  lavdata <- object@Data
  lavmodel <- object@Model
  lavpta <- object@pta

  # return list per group
  se_list <- vector("list", length = lavdata@ngroups)

  # get ACOV per group (reuse acov_obs if provided)
  acov_res <- lav_residuals_acov(
    object = object, type = type,
    z_type = z_type, h1_acov = h1_acov, acov_obs = acov_obs
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
      if (lavmodel@conditional.x) {
        # the moment vector also contains regression slopes; not ready yet
        lav_msg_stop(gettext("not ready yet!"))
      }

      # index map of the categorical moment vector (handles the mixed
      # continuous + ordinal case; reduces to [th, cov] when no num.idx)
      idx <- lav_residuals_cat_idx(lavmodel, g, nvar = nvar)

      # COV (off-diagonal); continuous variances go on the diagonal, the
      # ordinal-variable diagonal stays zero (variance is fixed at 1)
      cov_se <- lav_mat_vech_rev(sqrt(diag_acov[idx$cov]), diagonal = FALSE)
      if (length(idx$num_idx) > 0L) {
        diag(cov_se)[idx$num_idx] <- sqrt(diag_acov[idx$var])
      }

      # MEAN (continuous variables only; ordinal means are NA)
      mean_se <- rep(as.numeric(NA), nvar)
      if (length(idx$mean) > 0L) {
        mean_se[idx$num_idx] <- sqrt(diag_acov[idx$mean])
      }

      # TH
      th_se <- sqrt(diag_acov[idx$th])

      if (add_class) {
        class(cov_se) <- c("lavaan.matrix.symmetric", "matrix")
        class(mean_se) <- c("lavaan.vector", "numeric")
        class(th_se) <- c("lavaan.vector", "numeric")
      }
      if (add_labels) {
        rownames(cov_se) <- colnames(cov_se) <- ov_names[[g]]
        names(mean_se) <- ov_names[[g]]
        # thresholds only (no continuous means); equals th.mean for pure
        # categorical models, but excludes the continuous means in mixed ones
        names(th_se) <- lavpta$vnames$th[[g]]
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

# return summary statistics (rmr/srmr/crmr) as a list per group
#
# For each requested 'type' and each group, we summarize the residuals into
# one or more RMS columns (e.g. cov/mean/total or cor/thresholds/total). The
# heavy lifting (point estimate, SE, z-statistic, unbiased version, CI, close-
# fit test) is done by lav_residuals_summary_rms(); this function only decides
# *which* residuals enter *which* column, and slices the corresponding ACOV.
lav_residuals_summary <- function(object, type = c("rmr", "srmr", "crmr"),
                                  h1_acov = "unstructured",
                                  se = FALSE, zstat = FALSE, pvalue = FALSE,
                                  unbiased = FALSE, unbiased_se = FALSE,
                                  unbiased_ci = FALSE, unbiased_ci_level = 0.90,
                                  unbiased_zstat = FALSE,
                                  unbiased_test_val = 0.05,
                                  unbiased_pvalue = FALSE,
                                  add_class = FALSE, acov_obs = NULL) {
  # imply flags
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

  # change type name to the canonical rmr/srmr/crmr
  type[type == "raw"] <- "rmr"
  type[type == "cor.bentler"] <- "srmr"
  type[type == "cor.bollen"] <- "crmr"

  # slots
  lavdata <- object@Data
  lavmodel <- object@Model
  conditional_x <- lavmodel@conditional.x

  # residual point estimates (and, if needed, their ACOV) per requested type;
  # acov_obs is reused across types to avoid recomputing it (issue #291)
  if ("rmr" %in% type) {
    rmr_list <- lav_residuals(object = object, type = "raw")
    if (se || unbiased) {
      rmr_list_se <- lav_residuals_acov(
        object = object, type = "raw", z_type = "standardized",
        h1_acov = "unstructured", acov_obs = acov_obs
      )
    }
  }
  if ("srmr" %in% type) {
    srmr_list <- lav_residuals(object = object, type = "cor.bentler")
    if (se || unbiased) {
      srmr_list_se <- lav_residuals_acov(
        object = object, type = "cor.bentler", z_type = "standardized",
        h1_acov = "unstructured", acov_obs = acov_obs
      )
    }
  }
  if ("crmr" %in% type) {
    crmr_list <- lav_residuals(object = object, type = "cor.bollen")
    if (se || unbiased) {
      crmr_list_se <- lav_residuals_acov(
        object = object, type = "cor.bollen", z_type = "standardized",
        h1_acov = "unstructured", acov_obs = acov_obs
      )
    }
  }

  # select the residual (and ACOV) lists for a given type + group
  get_rms_lists <- function(type_name, g) {
    est <- switch(type_name,
      rmr = rmr_list[[g]], srmr = srmr_list[[g]], crmr = crmr_list[[g]]
    )
    se_g <- NULL
    if (se || unbiased) {
      se_g <- switch(type_name,
        rmr = rmr_list_se[[g]], srmr = srmr_list_se[[g]],
        crmr = crmr_list_se[[g]]
      )
    }
    list(est = est, se = se_g)
  }

  # compute one RMS summary column (captures the various reporting flags)
  do_rms <- function(stats, acov, pstar, type_name) {
    lav_residuals_summary_rms(
      stats = stats, acov = acov, se = se, zstat = zstat, pvalue = pvalue,
      unbiased = unbiased, unbiased_se = unbiased_se,
      unbiased_ci = unbiased_ci, unbiased_ci_level = unbiased_ci_level,
      unbiased_zstat = unbiased_zstat, unbiased_test_val = unbiased_test_val,
      unbiased_pvalue = unbiased_pvalue, pstar = pstar, type = type_name
    )
  }

  # return list per group
  sum_stat <- vector("list", length = lavdata@ngroups)

  for (g in seq_len(lavdata@ngroups)) {
    nvar <- object@pta$nvar[[g]] # block or group-based?

    out <- vector("list", length(type))
    names(out) <- type

    # categorical -- single level
    if (lavdata@nlevels == 1L && lavmodel@categorical) {
      mixed <- length(unlist(lavmodel@num.idx)) > 0L
      # raw (rmr) SEs are metric-agnostic and supported for the unconditional
      # mixed (continuous + ordinal) case; the standardized srmr/crmr SEs and
      # the conditional.x case need extra machinery and are not ready yet
      if ((se || unbiased) &&
        (conditional_x || (mixed && !all(type == "rmr")))) {
        lav_msg_stop(gettext("not ready yet"))
      }

      # moment-vector index map; for pure-categorical models idx$cov and idx$th
      # reproduce the simple [thresholds, correlations] split, but it also
      # handles the continuous variances/means in the mixed case
      nth <- length(lavmodel@th.idx[[g]])
      cov_nm <- if (conditional_x) "res.cov" else "cov"
      th_nm <- if (conditional_x) "res.th" else "th"
      idx <- if ((se || unbiased) && !conditional_x) {
        lav_residuals_cat_idx(lavmodel, g, nvar = nvar)
      } else {
        NULL
      }

      for (typ in seq_along(type)) {
        ty <- type[typ]
        sel <- get_rms_lists(ty, g)
        rms_list_g <- sel$est
        rms_list_se_g <- sel$se

        # COR (off-diagonal (co)variances)
        # pstar is p*(p+1)/2 for SRMR but p*(p-1)/2 for CRMR
        stats <- lav_mat_vech(rms_list_g[[cov_nm]], diagonal = FALSE)
        pstar <- if (ty == "crmr") length(stats) else length(stats) + nvar
        acov <- if (se || unbiased) {
          rms_list_se_g[idx$cov, idx$cov, drop = FALSE]
        } else {
          NULL
        }
        rms_cor <- do_rms(stats, acov, pstar, ty)

        # THRESHOLDS
        stats <- rms_list_g[[th_nm]]
        pstar <- length(stats)
        acov <- if (se || unbiased) {
          rms_list_se_g[idx$th, idx$th, drop = FALSE]
        } else {
          NULL
        }
        rms_th <- do_rms(stats, acov, pstar, ty)

        # TOTAL
        stats <- c(
          lav_mat_vech(rms_list_g[[cov_nm]], diagonal = FALSE),
          rms_list_g[[th_nm]]
        )
        pstar <- if (ty == "crmr") length(stats) else length(stats) + nvar
        acov <- NULL
        if (se || unbiased) {
          # match the [cor, thresholds] ordering of 'stats' in the mixed case;
          # pure-categorical keeps the historical full-matrix ordering
          acov <- if (mixed) {
            rms_list_se_g[c(idx$cov, idx$th), c(idx$cov, idx$th), drop = FALSE]
          } else {
            rms_list_se_g
          }
        }
        rms_total <- do_rms(stats, acov, pstar, ty)

        table_1 <- as.data.frame(cbind(rms_cor, rms_th, rms_total))
        colnames(table_1) <- c("cor", "thresholds", "total")
        if (add_class) {
          class(table_1) <- c("lavaan.data.frame", "data.frame")
        }
        out[[typ]] <- table_1
      } # type

      # continuous -- single level
    } else if (lavdata@nlevels == 1L) {
      if ((se || unbiased) && conditional_x) {
        lav_msg_stop(gettext("not ready yet"))
      }

      # ACOV order: [means(nvar), covariances(vech, with diagonal)]
      cov_nm <- if (conditional_x) "res.cov" else "cov"
      mean_nm <- if (conditional_x) "res.int" else "mean"

      for (typ in seq_along(type)) {
        ty <- type[typ]
        sel <- get_rms_lists(ty, g)
        rms_list_g <- sel$est
        rms_list_se_g <- sel$se

        # COV
        stats <- lav_mat_vech(rms_list_g[[cov_nm]])
        pstar <- length(stats)
        if (ty == "crmr") {
          # CRMR excludes the (co)variance diagonal
          pstar <- pstar -
            if (conditional_x) nrow(rms_list_g[[cov_nm]]) else nvar
        }
        acov <- NULL
        if (se || unbiased) {
          acov <- if (lavmodel@meanstructure) {
            rms_list_se_g[-seq_len(nvar), -seq_len(nvar), drop = FALSE]
          } else {
            rms_list_se_g
          }
        }
        rms_cov <- do_rms(stats, acov, pstar, ty)

        if (lavmodel@meanstructure) {
          # MEAN
          stats <- rms_list_g[[mean_nm]]
          pstar <- length(stats)
          acov <- if (se || unbiased) {
            rms_list_se_g[seq_len(nvar), seq_len(nvar), drop = FALSE]
          } else {
            NULL
          }
          rms_mean <- do_rms(stats, acov, pstar, ty)

          # TOTAL
          stats <- c(
            rms_list_g[[mean_nm]],
            lav_mat_vech(rms_list_g[[cov_nm]])
          )
          pstar <- length(stats)
          if (ty == "crmr") {
            pstar <- pstar - nvar
          }
          acov <- if (se || unbiased) rms_list_se_g else NULL
          rms_total <- do_rms(stats, acov, pstar, ty)

          table_1 <- as.data.frame(cbind(rms_cov, rms_mean, rms_total))
          colnames(table_1) <- c("cov", "mean", "total")
        } else {
          table_1 <- as.data.frame(cbind(rms_cov))
          colnames(table_1) <- "cov"
        }
        if (add_class) {
          class(table_1) <- c("lavaan.data.frame", "data.frame")
        }
        out[[typ]] <- table_1
      } # type

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
