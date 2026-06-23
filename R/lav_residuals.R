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

# conditional.x = TRUE is supported (both categorical and continuous): the
# moment vector then includes a regression-slopes block. Slopes are not
# rescaled in the cor metric (they are not rescaled in the residual point
# estimate either).
#
# Multilevel (twolevel): per-element raw/standardized/normalized residual SEs
# are supported for a single group; the within-level means are fixed at zero
# (degenerate moments) and get a zero SE. The cor.bentler/cor.bollen
# standardization, the summary statistics, and multigroup + multilevel are
# not ready yet.

# - change 0.6-6: we enforce observed.information = "h1" to ensure 'Q' is a
#                 projection matrix (see lav_residuals_acov)

# - change 0.6-13: fixed.x = TRUE is ignored (to conform with 'tradition')

# - change 0.7-1: refactor (issue #291): ACOV.obs (= the expensive ACOV of
#                  the observed h1 sample statistics) is computed at most once
#                  per call and reused for both the residual SEs and the
#                  summary statistics SEs (see lav_residuals_acov_obs())
# - change 0.7-1: the (undocumented, never-functional) custom.rmr argument
#                  has been removed
# - change 0.7-1: residual SEs/z-statistics and the (u)rmr/(u)srmr/(u)crmr
#                  summaries are now available for the categorical case,
#                  including mixed (continuous + ordinal), and for both
#                  categorical and continuous conditional.x (with a
#                  regression-slopes block). The cor.bentler/cor.bollen
#                  standardization uses a cov->cor jacobian built for the
#                  categorical moment ordering (see
#                  lav_residuals_cor_jacobian_cat()).
# - change 0.7-1: per-element raw/standardized/normalized residual SEs are now
#                  available for (single-group) multilevel models; the within-
#                  level means are degenerate (fixed at zero) and lav_model_h1
#                  _acov() now drops such zero-information moments before
#                  inverting (reinserting zeros) instead of returning NaN.
# - change 0.7-1: new h1= argument (lavResiduals) to supply a user-provided
#                  saturated model (a fitted lavaan object, a lavh1 list, or a
#                  moments list); it is swapped into object@h1 and used as the
#                  observed summary statistics (see lav_residuals_resolve_h1())

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
lavResiduals <- function(object, type = "cor.bentler", h1 = NULL,         # nolint start
                         se = FALSE, zstat = TRUE, summary = TRUE,
                         h1.acov = "unstructured",
                         add.type = TRUE, add.labels = TRUE, add.class = TRUE,
                         drop.list.single.group = TRUE,
                         maximum.number = length(res_vech),
                         output = "list") {                                # nolint end
  out <- lav_residuals(
    object = object, type = type, h1 = h1,
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

# resolve a user-provided h1= argument to a 'lavh1' list (list with an
# $implied element), so it can be swapped into object@h1 and used as the
# observed (saturated) sample statistics. Accepts:
#  - a fitted lavaan object       -> its @h1 slot
#  - a lavh1 list (has $implied)  -> as-is
#  - an implied-moments list      -> wrapped as list(implied = .)
# Returns NULL when no user override is requested (h1 is NULL or logical).
lav_residuals_resolve_h1 <- function(object, h1 = NULL) {
  if (is.null(h1) || is.logical(h1)) {
    return(NULL)
  }

  if (inherits(h1, "lavaan")) {
    lavh1 <- h1@h1
    if (length(lavh1) == 0L || is.null(lavh1$implied)) {
      lav_msg_stop(gettext(
        "the h1= lavaan object has no (non-empty) h1 slot; refit it with the
         default option h1 = TRUE"))
    }
  } else if (is.list(h1)) {
    if (!is.null(h1$implied)) {
      lavh1 <- h1 # already a lavh1 list
    } else if (!is.null(h1$cov) || !is.null(h1$res.cov)) {
      lavh1 <- list(implied = h1) # an 'implied' moments list
    } else {
      lav_msg_stop(gettext(
        "the h1= list must be a lavh1 list (with an 'implied' element) or an
         implied-moments list (with a 'cov' or 'res.cov' element)"))
    }
  } else {
    lav_msg_stop(gettext(
      "h1= must be a fitted lavaan object, a lavh1 list, or a moments list"))
  }

  # check that the saturated moments are compatible with the fitted object
  cov_nm <- if (object@Model@conditional.x) "res.cov" else "cov"
  usr_cov <- lavh1$implied[[cov_nm]]
  if (is.null(usr_cov) || !is.list(usr_cov)) {
    lav_msg_stop(gettextf(
      "the h1= implied moments must contain a (per-block) list element %s",
      dQuote(cov_nm)))
  }
  ref_cov <- object@h1$implied[[cov_nm]]
  if (!is.null(ref_cov)) {
    if (length(usr_cov) != length(ref_cov)) {
      lav_msg_stop(gettext(
        "the h1= saturated model has a different number of blocks/groups than
         the fitted object"))
    }
    for (b in seq_along(ref_cov)) {
      if (!all(dim(as.matrix(usr_cov[[b]])) == dim(as.matrix(ref_cov[[b]])))) {
        lav_msg_stop(gettextf(
          "the h1= saturated model and the fitted object have incompatible
           dimensions in block %d", b))
      }
    }
  }

  lavh1
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

  # user-provided saturated model (h1=)? if so, swap it into object@h1 so the
  # observed (saturated) sample statistics -- and their ACOV -- are taken from
  # it; the model-implied moments still come from 'object'
  if (is.null(h1)) {
    h1 <- TRUE
  } else if (is.logical(h1)) {
    # logical TRUE/FALSE: use object's own saturated model (keep as-is)
  } else if (inherits(h1, "lavaan") || is.list(h1)) {
    object@h1 <- lav_residuals_resolve_h1(object, h1)
    h1 <- TRUE
  } else {
    lav_msg_stop(gettext(
      "h1= must be NULL, logical, a fitted lavaan object, a lavh1 list, or a
       moments list"))
  }

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
    # per-element raw/standardized/normalized residual SEs are supported (for
    # a single group); the cor.bentler/cor.bollen standardization and the
    # summary statistics are not ready yet for the multilevel case
    summary <- FALSE
    if (type %in% c("cor.bentler", "cor.bollen") || lavdata@ngroups > 1L) {
      zstat <- se <- FALSE
      summary_options_1 <- lav_residuals_summary_options_off()
    }
  }

  # residual SEs/summaries are supported for the unconditional case (pure
  # categorical, continuous, and mixed continuous + ordinal) as well as the
  # conditional.x case (the moment vector then includes a regression-slopes
  # block, both for categorical and continuous models); nothing to disable.

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
          mixed_or_condx <- lavmodel@conditional.x ||
            length(unlist(lavmodel@num.idx)) > 0L
          if (mixed_or_condx) {
            # mixed continuous + ordinal and/or conditional.x: transform to
            # correlation metric using the (residual) covariance matrix
            cov_1 <- if (lavmodel@conditional.x) {
              sampstat[[g]][["res.cov"]]
            } else {
              sampstat[[g]][["cov"]]
            }
            nexo <- if (lavmodel@conditional.x) lavmodel@nexo[g] else 0L
            jac <- lav_residuals_cor_jacobian_cat(
              lavmodel, g, cov_1, "cor.bentler", nexo = nexo
            )
            acov_res[[g]] <- jac %*% acov_res[[g]] %*% t(jac)
          } else {
            # pure categorical: already in correlation metric, nothing to do
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
          if (lavmodel@conditional.x) {
            # moment order [vecr(cbind(int, slopes)), vech(cov)]: intercepts
            # are rescaled by 1/SD, regression slopes are left unchanged
            d_intsl <- lav_residuals_condx_intsl_scale(
              ss, lavmodel@nexo[g], lavmodel@meanstructure
            )
            gg <- lav_mat_bdiag(
              diag(d_intsl, nrow = length(d_intsl)), g_inv_sqrt
            )
          } else if (lavmodel@meanstructure) {
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
          mixed_or_condx <- lavmodel@conditional.x ||
            length(unlist(lavmodel@num.idx)) > 0L
          if (mixed_or_condx) {
            # mixed continuous + ordinal and/or conditional.x: transform to
            # correlation metric using the (residual) covariance matrix
            cov_1 <- if (lavmodel@conditional.x) {
              sampstat[[g]][["res.cov"]]
            } else {
              sampstat[[g]][["cov"]]
            }
            nexo <- if (lavmodel@conditional.x) lavmodel@nexo[g] else 0L
            jac <- lav_residuals_cor_jacobian_cat(
              lavmodel, g, cov_1, "cor.bollen", nexo = nexo
            )
            acov_res[[g]] <- jac %*% acov_res[[g]] %*% t(jac)
          } else {
            # pure categorical: already in correlation metric, nothing to do
          }
        } else {
          # here we use the Maydeu-Olivares (2017) approach, see eq 17
          cov_1 <- if (lavmodel@conditional.x) {
            sampstat[[g]][["res.cov"]]
          } else {
            sampstat[[g]][["cov"]]
          }
          f1 <- lav_deriv_cov2cor_b(cov_1)
          if (lavmodel@conditional.x) {
            ss <- 1 / sqrt(diag(cov_1))
            d_intsl <- lav_residuals_condx_intsl_scale(
              ss, lavmodel@nexo[g], lavmodel@meanstructure
            )
            ff <- lav_mat_bdiag(diag(d_intsl, nrow = length(d_intsl)), f1)
          } else if (lavmodel@meanstructure) {
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
# produced by lav_samp_wls_obs():
#   [ th-block (thresholds + continuous means/intercepts, by variable),
#     regression slopes (vec, only if conditional.x = TRUE; nvar * nexo),
#     continuous variances (num.idx),
#     covariances (vech, no diagonal) ]
#
# Returns, as positions into that vector: $th (thresholds), $mean (continuous
# means/intercepts), $slopes (regression slopes), $var (continuous variances),
# $cov (off-diagonal (co)variances), plus $num_idx (continuous variable
# indices), $nvar and $nexo. For conditional.x = FALSE, $slopes is empty and
# nexo = 0.
lav_residuals_cat_idx <- function(lavmodel, g = 1L, nvar = NULL, nexo = 0L) {
  th_idx <- lavmodel@th.idx[[g]]
  nth_block <- length(th_idx)
  num_idx <- lavmodel@num.idx[[g]]
  n_num <- length(num_idx)
  if (is.null(nvar)) {
    nvar <- length(unique(th_idx[th_idx != 0])) + n_num
  }
  n_slopes <- nvar * nexo
  n_offdiag <- nvar * (nvar - 1) / 2

  list(
    th = which(th_idx != 0L), # threshold positions in the th-block
    mean = which(th_idx == 0L), # continuous mean/intercept positions
    slopes = if (n_slopes > 0L) nth_block + seq_len(n_slopes) else integer(0),
    var = nth_block + n_slopes + seq_len(n_num), # continuous variances
    cov = nth_block + n_slopes + n_num + seq_len(n_offdiag), # off-diag covs
    num_idx = num_idx,
    nvar = nvar,
    nexo = nexo
  )
}

# Build the cov->cor standardization jacobian for the *categorical* moment
# vector ordering [th-block, (slopes), continuous variances, vech(cov, no
# diagonal)], used to transform the raw-residual ACOV into the correlation
# metric for mixed (continuous + ordinal) models. 'cov_1' is the observed (h1)
# (residual) covariance matrix (ordinal variables have a unit diagonal, so
# 1/sqrt(diag) is 1 for them and 1/SD for the continuous ones). 'type' is
# either "cor.bentler" (observed SDs treated as fixed -> diagonal scaling) or
# "cor.bollen" (full cov2cor jacobian). 'nexo' > 0 for conditional.x (a slopes
# block is present). Returns J such that ACOV.cor = J %*% ACOV.raw %*% t(J).
#
# Thresholds and regression slopes are left unchanged (slopes are not rescaled
# in the residual point estimate either, cf. lav_residuals_rescale()). For a
# pure-categorical model (no continuous variables) this reduces to the
# identity; we therefore only use it in the mixed case.
lav_residuals_cor_jacobian_cat <- function(lavmodel, g, cov_1, type,
                                           nexo = 0L) {
  idx <- lav_residuals_cat_idx(lavmodel, g, nvar = nrow(cov_1), nexo = nexo)
  nvar <- idx$nvar
  num_idx <- idx$num_idx
  p <- length(idx$th) + length(idx$mean) + length(idx$slopes) +
    length(idx$var) + length(idx$cov)

  ss <- 1 / sqrt(diag(cov_1)) # 1 for ordinal (unit diagonal), 1/SD otherwise

  jac <- matrix(0, nrow = p, ncol = p)
  # thresholds and regression slopes: unchanged
  if (length(idx$th) > 0L) {
    jac[cbind(idx$th, idx$th)] <- 1
  }
  if (length(idx$slopes) > 0L) {
    jac[cbind(idx$slopes, idx$slopes)] <- 1
  }
  # continuous means/intercepts: scaled by 1/SD (observed SD treated as fixed)
  if (length(idx$mean) > 0L) {
    jac[cbind(idx$mean, idx$mean)] <- ss[num_idx]
  }

  if (type == "cor.bentler") {
    # diagonal scaling; observed SDs treated as fixed constants
    if (length(idx$var) > 0L) {
      jac[cbind(idx$var, idx$var)] <- ss[num_idx]^2
    }
    if (length(idx$cov) > 0L) {
      jac[cbind(idx$cov, idx$cov)] <-
        lav_mat_vech(tcrossprod(ss), diagonal = FALSE)
    }
  } else {
    # cor.bollen: full cov2cor jacobian, restricted to the moments that are
    # actually free (continuous variances + off-diagonal covariances); the
    # ordinal variances are fixed at 1 and are dropped
    f1 <- lav_deriv_cov2cor_b(cov_1)
    diagh <- lav_mat_diagh_idx(nvar)
    offdiag <- setdiff(seq_len(nvar * (nvar + 1) / 2), diagh)
    full_idx <- c(diagh[num_idx], offdiag)
    catcov <- c(idx$var, idx$cov) # variances first, then off-diagonals
    jac[catcov, catcov] <- f1[full_idx, full_idx]
  }

  jac
}

# scaling factors for the [vecr(cbind(intercepts, slopes))] block of the
# *continuous* conditional.x moment vector: intercepts are rescaled by 1/SD
# (ss), regression slopes are left unchanged (slopes are not rescaled in the
# residual point estimate either). Returns the diagonal in vecr (row-major)
# order. If meanstructure = FALSE the block contains only the slopes.
lav_residuals_condx_intsl_scale <- function(ss, nexo, meanstructure) {
  nvar <- length(ss)
  b <- if (meanstructure) {
    cbind(ss, matrix(1, nvar, nexo))
  } else {
    matrix(1, nvar, nexo)
  }
  as.vector(t(b))
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

  # return list per block (= per group for single-level data; for multilevel
  # data there is one block per level within each group)
  se_list <- vector("list", length = lavmodel@nblocks)

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
      conditional_x <- lavmodel@conditional.x
      nexo <- if (conditional_x) lavmodel@nexo[g] else 0L
      # endogenous variables only (res.cov is over the endogenous ov);
      # for conditional.x the exogenous covariates are NOT in the moment vector
      nvar_endo <- if (conditional_x) {
        length(lavpta$vnames$ov.nox[[g]])
      } else {
        nvar
      }

      # index map of the categorical moment vector (handles the mixed
      # continuous + ordinal case, and the conditional.x slopes block)
      idx <- lav_residuals_cat_idx(lavmodel, g, nvar = nvar_endo, nexo = nexo)
      num_idx <- idx$num_idx

      # COV (off-diagonal); continuous variances go on the diagonal, the
      # ordinal-variable diagonal stays zero (variance is fixed at 1)
      cov_se <- lav_mat_vech_rev(sqrt(diag_acov[idx$cov]), diagonal = FALSE)
      if (length(num_idx) > 0L) {
        diag(cov_se)[num_idx] <- sqrt(diag_acov[idx$var])
      }

      # MEAN/INTERCEPT (continuous variables only; ordinal entries are NA)
      mean_se <- rep(as.numeric(NA), nvar_endo)
      if (length(idx$mean) > 0L) {
        mean_se[num_idx] <- sqrt(diag_acov[idx$mean])
      }

      # TH (thresholds only)
      th_se <- sqrt(diag_acov[idx$th])

      if (conditional_x) {
        # SLOPES (nvar_endo x nexo, vec ordering)
        slopes_se <- matrix(sqrt(diag_acov[idx$slopes]),
          nrow = nvar_endo, ncol = nexo
        )
        # exogenous moments are conditioned on (not in the vector): no SE
        covx_se <- matrix(as.numeric(NA), nexo, nexo)
        meanx_se <- rep(as.numeric(NA), nexo)
        if (add_class) {
          class(cov_se) <- c("lavaan.matrix.symmetric", "matrix")
          class(covx_se) <- c("lavaan.matrix.symmetric", "matrix")
          class(slopes_se) <- c("lavaan.matrix", "matrix")
          class(mean_se) <- class(meanx_se) <- class(th_se) <-
            c("lavaan.vector", "numeric")
        }
        if (add_labels) {
          rownames(cov_se) <- colnames(cov_se) <- lavpta$vnames$ov.nox[[g]]
          names(mean_se) <- lavpta$vnames$ov.nox[[g]]
          names(th_se) <- lavpta$vnames$th[[g]]
          rownames(slopes_se) <- lavpta$vnames$ov.nox[[g]]
          colnames(slopes_se) <- lavpta$vnames$ov.x[[g]]
          rownames(covx_se) <- colnames(covx_se) <- lavpta$vnames$ov.x[[g]]
          names(meanx_se) <- lavpta$vnames$ov.x[[g]]
        }
        # order must match lav_inspect_sampstat(): res.cov, res.int, res.th,
        # res.slopes, cov.x, mean.x
        se_list[[g]] <- list(
          res.cov.se = cov_se, res.int.se = mean_se, res.th.se = th_se,
          res.slopes.se = slopes_se, cov.x.se = covx_se, mean.x.se = meanx_se
        )
      } else {
        if (add_class) {
          class(cov_se) <- c("lavaan.matrix.symmetric", "matrix")
          class(mean_se) <- c("lavaan.vector", "numeric")
          class(th_se) <- c("lavaan.vector", "numeric")
        }
        if (add_labels) {
          rownames(cov_se) <- colnames(cov_se) <- ov_names[[g]]
          names(mean_se) <- ov_names[[g]]
          # thresholds only (no continuous means); equals th.mean for pure
          # categorical models, excludes continuous means in mixed ones
          names(th_se) <- lavpta$vnames$th[[g]]
        }
        se_list[[g]] <- list(
          cov.se = cov_se, mean.se = mean_se,
          th.se = th_se
        )
      }

      # continuous -- single level
    } else if (lavdata@nlevels == 1L) {
      if (lavmodel@conditional.x) {
        # moment order: [vecr(cbind(res.int, res.slopes)), vech(res.cov)]
        # (intercepts only if meanstructure); res.slopes is nvar.endo x nexo
        nexo <- lavmodel@nexo[g]
        nvar_endo <- length(lavpta$vnames$ov.nox[[g]])
        ncol_b <- if (lavmodel@meanstructure) 1L + nexo else nexo
        n_intsl <- nvar_endo * ncol_b
        # row-major (vecr) reshape of the intercept/slopes block
        b_se <- matrix(sqrt(diag_acov[seq_len(n_intsl)]),
          nrow = nvar_endo, ncol = ncol_b, byrow = TRUE
        )
        cov_se <- lav_mat_vech_rev(sqrt(diag_acov[-seq_len(n_intsl)]),
          diagonal = TRUE
        )
        covx_se <- matrix(as.numeric(NA), nexo, nexo)
        meanx_se <- rep(as.numeric(NA), nexo)
        if (lavmodel@meanstructure) {
          int_se <- b_se[, 1]
          slopes_se <- b_se[, -1, drop = FALSE]
        } else {
          int_se <- NULL
          slopes_se <- b_se
        }
        if (add_class) {
          class(cov_se) <- c("lavaan.matrix.symmetric", "matrix")
          class(covx_se) <- c("lavaan.matrix.symmetric", "matrix")
          class(slopes_se) <- c("lavaan.matrix", "matrix")
          class(meanx_se) <- c("lavaan.vector", "numeric")
          if (!is.null(int_se)) class(int_se) <- c("lavaan.vector", "numeric")
        }
        if (add_labels) {
          rownames(cov_se) <- colnames(cov_se) <- lavpta$vnames$ov.nox[[g]]
          rownames(slopes_se) <- lavpta$vnames$ov.nox[[g]]
          colnames(slopes_se) <- lavpta$vnames$ov.x[[g]]
          rownames(covx_se) <- colnames(covx_se) <- lavpta$vnames$ov.x[[g]]
          names(meanx_se) <- lavpta$vnames$ov.x[[g]]
          if (!is.null(int_se)) names(int_se) <- lavpta$vnames$ov.nox[[g]]
        }
        # order must match lav_inspect_sampstat(): res.cov, res.int (if
        # meanstructure), res.slopes, cov.x, mean.x
        if (lavmodel@meanstructure) {
          se_list[[g]] <- list(
            res.cov.se = cov_se, res.int.se = int_se,
            res.slopes.se = slopes_se, cov.x.se = covx_se, mean.x.se = meanx_se
          )
        } else {
          se_list[[g]] <- list(
            res.cov.se = cov_se, res.slopes.se = slopes_se,
            cov.x.se = covx_se, mean.x.se = meanx_se
          )
        }
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
      if (lavdata@ngroups > 1L) {
        # multigroup + multilevel: not ready yet
        lav_msg_stop(gettext("not ready yet"))
      }
      # the moment vector stacks the level-specific blocks, each ordered
      # [mean, vech(cov)]; the within-level means are fixed at zero and thus
      # have a zero standard error
      offset <- 0L
      for (b in seq_len(lavmodel@nblocks)) {
        nvb <- object@pta$nvar[[b]]
        pstar_b <- nvb * (nvb + 1) / 2
        mean_se <- sqrt(diag_acov[offset + seq_len(nvb)])
        cov_se <- lav_mat_vech_rev(
          sqrt(diag_acov[offset + nvb + seq_len(pstar_b)])
        )
        offset <- offset + nvb + pstar_b
        if (add_class) {
          class(cov_se) <- c("lavaan.matrix.symmetric", "matrix")
          class(mean_se) <- c("lavaan.vector", "numeric")
        }
        if (add_labels) {
          rownames(cov_se) <- colnames(cov_se) <- ov_names[[b]]
          names(mean_se) <- ov_names[[b]]
        }
        se_list[[b]] <- list(cov.se = cov_se, mean.se = mean_se)
      }
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
      # rmr/srmr/crmr SEs are supported for the unconditional case (pure
      # categorical as well as mixed continuous + ordinal) and for the
      # conditional.x case (the moment vector then includes a slopes block,
      # which does not enter the cor/threshold summaries)

      # endogenous variables only; for conditional.x the exogenous covariates
      # are not part of the moment vector
      nexo <- if (conditional_x) lavmodel@nexo[g] else 0L
      nvar_endo <- if (conditional_x) {
        length(object@pta$vnames$ov.nox[[g]])
      } else {
        nvar
      }

      # moment-vector index map; for pure-categorical models idx$cov and idx$th
      # reproduce the simple [thresholds, correlations] split, but it also
      # handles the continuous variances/means and the conditional.x slopes
      cov_nm <- if (conditional_x) "res.cov" else "cov"
      th_nm <- if (conditional_x) "res.th" else "th"
      idx <- if (se || unbiased) {
        lav_residuals_cat_idx(lavmodel, g, nvar = nvar_endo, nexo = nexo)
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
        pstar <- if (ty == "crmr") length(stats) else length(stats) + nvar_endo
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
        # 'stats' and 'acov' must use the same ordering, otherwise the
        # unbiased e_ve = t(stats) %*% acov %*% stats quadratic form is wrong.
        # We use the canonical lavaan WLS ordering [thresholds, covariances]
        # (matching lav_samp_wls_obs() and hence the ACOV); for pure
        # categorical models this leaves the ACOV in its original order.
        stats <- c(
          rms_list_g[[th_nm]],
          lav_mat_vech(rms_list_g[[cov_nm]], diagonal = FALSE)
        )
        pstar <- if (ty == "crmr") length(stats) else length(stats) + nvar_endo
        acov <- NULL
        if (se || unbiased) {
          acov <- rms_list_se_g[c(idx$th, idx$cov), c(idx$th, idx$cov),
            drop = FALSE
          ]
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
      cov_nm <- if (conditional_x) "res.cov" else "cov"
      mean_nm <- if (conditional_x) "res.int" else "mean"
      # endogenous variables only (the cov/mean summaries are over the
      # endogenous ov; the conditional.x slopes block does not enter them)
      nvar_endo <- if (conditional_x) {
        length(object@pta$vnames$ov.nox[[g]])
      } else {
        nvar
      }

      # ACOV positions of the mean/intercept block and the covariance block.
      # Non-conditional.x order: [means(nvar), vech(cov)]. Conditional.x order:
      # [vecr(cbind(int, slopes)), vech(cov)], so the intercepts are row-major
      # interleaved with the slopes (which are skipped here).
      if (conditional_x) {
        nexo <- lavmodel@nexo[g]
        ncol_b <- if (lavmodel@meanstructure) 1L + nexo else nexo
        n_intsl <- nvar_endo * ncol_b
        mean_acov_idx <- if (lavmodel@meanstructure) {
          seq.int(1L, by = ncol_b, length.out = nvar_endo)
        } else {
          integer(0)
        }
        cov_acov_idx <- n_intsl + seq_len(nvar_endo * (nvar_endo + 1) / 2)
      } else {
        mean_acov_idx <- if (lavmodel@meanstructure) {
          seq_len(nvar_endo)
        } else {
          integer(0)
        }
        cov_acov_idx <- length(mean_acov_idx) +
          seq_len(nvar_endo * (nvar_endo + 1) / 2)
      }

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
          pstar <- pstar - nvar_endo
        }
        acov <- if (se || unbiased) {
          rms_list_se_g[cov_acov_idx, cov_acov_idx, drop = FALSE]
        } else {
          NULL
        }
        rms_cov <- do_rms(stats, acov, pstar, ty)

        if (lavmodel@meanstructure) {
          # MEAN / INTERCEPT
          stats <- rms_list_g[[mean_nm]]
          pstar <- length(stats)
          acov <- if (se || unbiased) {
            rms_list_se_g[mean_acov_idx, mean_acov_idx, drop = FALSE]
          } else {
            NULL
          }
          rms_mean <- do_rms(stats, acov, pstar, ty)

          # TOTAL (means/intercepts + covariances; slopes excluded)
          stats <- c(
            rms_list_g[[mean_nm]],
            lav_mat_vech(rms_list_g[[cov_nm]])
          )
          pstar <- length(stats)
          if (ty == "crmr") {
            pstar <- pstar - nvar_endo
          }
          acov <- if (se || unbiased) {
            rms_list_se_g[
              c(mean_acov_idx, cov_acov_idx),
              c(mean_acov_idx, cov_acov_idx),
              drop = FALSE
            ]
          } else {
            NULL
          }
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
