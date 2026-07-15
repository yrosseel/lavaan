# step 1 in SAM: fitting the measurement blocks

lav_sam_step1 <- function(cmd = "sem", mm_list = NULL, mm_args = list(),
                          fit = NULL, sam_method = "local") {
  lavoptions <- fit@Options
  lavpta <- fit@pta
  pt_1 <- fit@ParTable
  nblocks <- lavpta$nblocks
  ngroups <- lavpta$ngroups

  if (lav_verbose()) {
    cat("Fitting the measurement part:\n")
  }

  # total number of free parameters
  if (fit@Model@ceq.simple.only) {
    npar <- fit@Model@nx.unco
    pt_free <- pt_1$free
    pt_free[pt_free > 0] <- seq_len(npar)
  } else {
    npar <- fit@Model@nx.free
    pt_free <- pt_1$free
  }
  if (npar < 1L) {
    lav_msg_stop(gettext("model does not contain any free parameters"))
  }

  # do we have at least 1 'regular' (measured) latent variable?
  lv_names <- unique(unlist(fit@pta$vnames$lv.regular))
  if (length(lv_names) == 0L) {
    lav_msg_warn(gettext("structural model does not contain any (measured)
                          latent variables; consider using sem() instead"))
  }

  # check for higher-order factors (supported since 0.6-20):
  # remove the lower-order factors from lv_names
  lv_ind_names <- unique(unlist(fit@pta$vnames$lv.ind))
  if (length(lv_ind_names) > 0L) {
    lv_names <- lv_names[!lv_names %in% lv_ind_names]
  }

  # how many measurement models?
  if (!is.null(mm_list)) {
    n_mmblocks <- length(mm_list)
    # check each measurement block
    for (b in seq_len(n_mmblocks)) {
      # check if we can find all lv names in LV.names
      if (!all(unlist(mm_list[[b]]) %in% lv_names)) {
        tmp <- unlist(mm_list[[b]])
        lav_msg_stop(gettext("mm.list contains unknown latent variable(s):"),
          lav_msg_view(tmp[!tmp %in% lv_names], "none"))
      }
      # make list per block
      if (!is.list(mm_list[[b]])) {
        mm_list[[b]] <- rep(list(mm_list[[b]]), nblocks)
      } else {
        if (length(mm_list[[b]]) != nblocks) {
          lav_msg_stop(gettextf(
            "mm.list block %1$s has length %2$s but nblocks = %3$s",
            b, length(mm_list[[b]]), nblocks))
        }
      }
    }
  } else {
    mm_list_per_block <- lav_sam_get_mmlist(fit) # we get a list PER BLOCK
    if (fit@Data@nlevels > 1L) {
      # multilevel: join the within and between clusters that share
      # observed indicators into a single (two-level) measurement block
      mm_list <- lav_sam_mmlist_merge_blocks(fit, mm_list_per_block)
    } else {
      # for now, we only keep the first block
      mm_list <- mm_list_per_block[[1]]
    }
    n_mmblocks <- length(mm_list)
  }

  # adjust options for measurement models
  lavoptions_mm <- lavoptions
  # "yuan.chan" is a SAM-global test for the JOINT model only; the measurement
  # blocks use the ordinary (standard) test
  if (any(lavoptions_mm$test == "yuan.chan")) {
    lavoptions_mm$test <- "standard"
  }
  lavoptions_mm$optim.bounds <- NULL
  if (lavoptions$se %in% c("none", "bootstrap")) {
    lavoptions_mm$se <- "none"
  } else {
    # single-level clustered data: use the cluster-robust variants, so that
    # the step 1 standard errors (and Sigma.11) reflect the clustering;
    # for the Gamma (NACOV) based sandwich (se = "robust.sem"), no switch
    # is needed: the NACOV is computed with the cluster correction whenever
    # the data contains a cluster variable
    clustered_flag <- fit@Data@nlevels == 1L && length(fit@Data@cluster) > 0L
    # categorical?
    # (PML first: it is a categorical estimator, but se = "robust.sem" is
    # not valid for PML -- it needs the huber-white sandwich)
    if (lavoptions$estimator.orig == "PML") {
      lavoptions_mm$se <- "robust.huber.white"
    } else if (fit@Model@categorical) {
      lavoptions_mm$se <- "robust.sem"
    } else if (lavoptions$estimator.orig == "MLM") {
      lavoptions_mm$se <- "robust.sem"
    } else if (lavoptions$estimator.orig == "MLR") {
      lavoptions_mm$se <- if (clustered_flag) "robust.cluster"
                          else "robust.huber.white"
    } else if (clustered_flag) {
      lavoptions_mm$se <- "robust.cluster"
    } else {
      lavoptions_mm$se <- "standard" # may be overridden later
    }
  }
  # note: we keep the tests, as we need them for the summary info about MM
  lavoptions_mm$check.post <- FALSE # neg lv variances may be overridden
  lavoptions_mm$check.gradient <- FALSE # too sensitive in large model (global)
  lavoptions_mm$baseline <- FALSE
  lavoptions_mm$fit.by.level <- FALSE
  lavoptions_mm$bounds <- "wide.zerovar"

  # conditional.x: for CONTINUOUS data the measurement models are always
  # fitted unconditionally -- this does not affect the measurement parameters
  # (if the covariates act on the structural part only). For CATEGORICAL
  # data, however, the latent-response scale differs between the marginal
  # (Var(y*) = 1) and the conditional (Var(y*|x) = 1) parameterization, and
  # the marginal polychorics are misspecified when the covariates are
  # non-normal (eg binary). The measurement blocks are then fitted on the
  # same conditional scale as the joint model: conditional.x = TRUE, with
  # saturated eta ~ x regressions (added via lav_pt_subset_mm() below).
  if (fit@Model@categorical && lavoptions$conditional.x) {
    lavoptions_mm$conditional.x <- TRUE
  } else {
    lavoptions_mm$conditional.x <- FALSE
  }

  # override with user-specified mm.args
  lavoptions_mm <- modifyList(lavoptions_mm, mm_args)

  # create MM slot_options
  slot_options_mm <- lav_options_set(lavoptions_mm)

  # we assume the same number/names of lv's per group!!!
  mm_fit <- vector("list", n_mmblocks) # fitted object

  # for joint model later
  if (!lavoptions$se %in% c("none", "bootstrap")) {
    sigma_11_1 <- matrix(0, npar, npar)
    colnames(sigma_11_1) <- rownames(sigma_11_1) <-
      lav_pt_labels(fit@ParTable, type = "free")
  }
  step1_free_idx <- integer(0L)
  block_mm_idx  <- vector("list", length = n_mmblocks)
  block_ptm_idx <- vector("list", length = n_mmblocks)

  # NOTE: we should explicitly add zero-constrained LV covariances
  # to PT, and keep them zero in PTM
  add_lv_cov <- cmd != "lavaan"

  # if data.type == "moment", we need the h1 sample statistics (per group)
  if (fit@Data@data.type == "moment") {
    h1_all <- lavTech(fit, "h1", add.labels = TRUE)
  }

  # fit mm model for each measurement block
  for (mm in seq_len(n_mmblocks)) {
    if (lav_verbose()) {
      cat(
        "  block ", mm, "[",
        paste(mm_list[[mm]], collapse = " "), "]\n"
      )
    }

    # create parameter table for this measurement block only
    ptm <- lav_pt_subset_mm(
      pt_1 = pt_1,
      add_lv_cov = add_lv_cov,
      add_idx = TRUE,
      lv_names = mm_list[[mm]],
      add_exo_reg = slot_options_mm$conditional.x,
      ov_names_x = fit@pta$vnames$ov.x
    )
    mm_idx <- attr(ptm, "idx")
    attr(ptm, "idx") <- NULL
    ptm$est <- NULL
    ptm$se <- NULL
    block_mm_idx[[mm]] <- mm_idx

    # check for categorical in PTM in this mm-block
    # (a fully continuous block is also fitted with conditional.x = TRUE
    # when the joint model is categorical + conditional.x: its statistics
    # must live in the same conditional space as the joint model)
    if (!any(ptm$op == "|")) {
      slot_options_mm$categorical <- FALSE
      slot_options_mm$.categorical <- FALSE
    }

    # update slot_data for this measurement block
    # (under conditional.x, lavData keeps the exogenous covariates OUT of
    # ov.names/X -- they live in ov.names.x/eXo -- so the block ov.names
    # must be the indicators only)
    ov_names_block <- lapply(1:ngroups, function(g) {
      if (slot_options_mm$conditional.x) {
        unique(unlist(lav_pt_vnames(ptm, type = "ov.nox", group = g)))
      } else {
        unique(unlist(lav_pt_vnames(ptm, type = "ov", group = g)))
      }
    })
    slot_data_block <- lav_data_update_subset(fit@Data,
      ov_names = ov_names_block
    )
    if (!slot_options_mm$conditional.x) {
      # get rid of ov.names.x
      slot_data_block@ov.names.x <-
        lapply(seq_len(nblocks), function(x) character(0L))
      slot_data_block@eXo <-
        lapply(seq_len(nblocks), function(x) NULL)
    } else {
      # keep ALL exogenous covariates for the conditional block
      # (lav_data_update_subset() dropped them, as they are not part of
      # ov_names_block); also restore their entries in the ov table
      slot_data_block@ov.names.x <- fit@Data@ov.names.x
      slot_data_block@eXo <- fit@Data@eXo
      ov_keep_idx <- which(fit@Data@ov$name %in%
        unique(c(unlist(ov_names_block), unlist(fit@Data@ov.names.x))))
      slot_data_block@ov <- lapply(fit@Data@ov, "[", ov_keep_idx)
    }

    # if data.type == "moment", (re)create sample.cov and sample.nobs
    if (fit@Data@data.type == "moment") {
      mm_sample_mean <- NULL
      if (ngroups == 1L) {
        mm_sample_cov <- h1_all[[1L]]$cov
        if (fit@Model@meanstructure) {
          mm_sample_mean <- h1_all[[1L]]$mean
        }
        mm_sample_nobs <- fit@SampleStats@nobs[[1L]]
      } else {
        mm_sample_cov <- lapply(seq_len(ngroups), function(x) {
          h1_all[[x]]$cov[ov_names_block[[x]], ov_names_block[[x]]]
        })
        if (fit@Model@meanstructure) {
          mm_sample_mean <- lapply(seq_len(ngroups), function(x) {
            h1_all[[x]]$mean[ov_names_block[[x]]]
          })
        }
        mm_sample_nobs <- fit@SampleStats@nobs
      }
    }

    # handle single block 1-factor CFA with (only) two indicators
    # (under conditional.x the block variables include the exogenous
    # covariates -- count the indicators only)
    n_ind_block <- length(unlist(ov_names_block))
    if (slot_options_mm$conditional.x) {
      n_ind_block <- length(unlist(lapply(1:ngroups, function(g) {
        unique(unlist(lav_pt_vnames(ptm, type = "ov.nox", group = g)))
      })))
    }
    if (n_ind_block == 2L && ngroups == 1L) {
      lambda_idx <- which(ptm$op == "=~")
      # check if both factor loadings are fixed
      # (note: this assumes std.lv = FALSE)
      if (any(ptm$free[lambda_idx] != 0)) {
        ptm$free[lambda_idx] <- 0L
        ptm$ustart[lambda_idx] <- 1
        ptm$start[lambda_idx] <- 1
        # adjust free counter
        free_idx <- which(as.logical(ptm$free))
        if (length(free_idx) > 0L) {
          ptm$free[free_idx] <- seq_along(free_idx)
        }
        lav_msg_warn(gettextf(
          "measurement block [%1$s] (%2$s) contains only two indicators;
          -> fixing both factor loadings to unity",
          mm, lav_msg_view(mm_list[[mm]], "none")))
      }
    }

    # fit this measurement model only
    # (question: can we re-use even more slots?)
    if (fit@Data@data.type == "full") {
      fit_mm_block <- lavaan(
        model = ptm, slot_data = slot_data_block,
        slot_options = slot_options_mm, debug = FALSE, verbose = FALSE
      )
    } else if (fit@Data@data.type == "moment") {
      slot_options_mm$sample.cov.rescale <- FALSE
      fit_mm_block <- lavaan(
        model = ptm, slot_data = slot_data_block,
        sample_cov = mm_sample_cov, sample_mean = mm_sample_mean,
        sample_nobs = mm_sample_nobs,
        slot_options = slot_options_mm, debug = FALSE, verbose = FALSE
      )
    }

    # check convergence
    if (!lavInspect(fit_mm_block, "converged")) {
      # warning for now, but this is not good!
      lav_msg_warn(gettextf(
        "measurement model for %s did not converge!",
        lav_msg_view(mm_list[[mm]], "none")))
    }

    # check that this measurement block is identified on its own: in the
    # SAM approach every block is estimated in isolation, so a block with
    # more free parameters than sample statistics (negative degrees of
    # freedom) cannot be estimated, even if the full (joint) model is
    # identified (eg a factor with three indicators plus a residual
    # covariance between two of them). Without this check the block fit
    # produces garbage and downstream computations fail cryptically.
    blk_df <- fit_mm_block@test[[1]]$df
    if (!is.null(blk_df) && !is.na(blk_df) && blk_df < 0L) {
      lav_msg_stop(gettextf(
        "measurement block %1$s is not identified on its own: it has more
         free parameters than sample statistics (df = %2$d). The SAM
         approach estimates each measurement block in isolation, so every
         block must be identified by itself, even if the full model is
         identified. Consider simplifying the measurement block (eg
         removing residual covariances), combining blocks via mm.list, or
         using sem() instead.",
        lav_msg_view(mm_list[[mm]], "none"), blk_df))
    }

    # store fitted measurement model
    mm_fit[[mm]] <- fit_mm_block

    # fill in point estimates measurement block (including slack values)
    ptm <- mm_fit[[mm]]@ParTable
    # note: the row-numbers in PT that correspond to the rows in PTM
    # are given by mm_idx
    ptm_idx <- which((ptm$free > 0L | ptm$op %in% c(":=", "<", ">")) &
      ptm$user != 3L)
    block_ptm_idx[[mm]] <- ptm_idx
    pt_1$est[mm_idx[ptm_idx]] <- ptm$est[ptm_idx]

    # if categorical, add non-free residual variances
    if (fit_mm_block@Model@categorical || fit_mm_block@Model@correlation) {
      extra_idx <- which(ptm$op %in% c("~~", "~*~") &
        ptm$lhs == ptm$rhs &
        ptm$user == 0L &
        ptm$free == 0L &
        ptm$ustart == 1)
      if (length(extra_idx) > 0L) {
        pt_1$est[mm_idx[extra_idx]] <- ptm$est[extra_idx]
      }
    }
    # if EFA, add user=7 values (but do not add to ptm.idx)
    user7_idx <- which(ptm$user == 7L)
    if (length(user7_idx)) {
      pt_1$est[mm_idx[user7_idx]] <- ptm$est[user7_idx]
    }

    # add step1.free.idx
    par_idx <- pt_free[mm_idx[ptm_idx]]
    # store (ordered) indices in step1.free.idx
    step1_free_idx <- c(step1_free_idx, sort.int(par_idx)) # all combined

    # fill in standard errors measurement block
    if (!lavoptions$se %in% c("none", "bootstrap")) {
      ptm_free <- ptm$free
      if (fit_mm_block@Model@ceq.simple.only) {
        ptm_free[ptm_free > 0] <- seq_len(fit_mm_block@Model@nx.unco)
      }

      ptm_se_idx <- which((ptm$free > 0L) & ptm$user != 3L) # no :=, <, >
      pt_1$se[mm_idx[ptm_se_idx]] <- ptm$se[ptm_se_idx]

      # fill in variance matrix for this measurement block
      sigma_11 <- mm_fit[[mm]]@vcov$vcov
      if (is.null(sigma_11)) {
        # the block vcov could not be computed (eg the block information
        # matrix could not be inverted, typically an identification issue
        # that the df >= 0 check above cannot catch)
        lav_msg_stop(gettextf(
          "the variance matrix of measurement block %s could not be
           computed (its information matrix could not be inverted); the
           block may not be identified on its own (empirically). Consider
           simplifying the measurement block, combining blocks via
           mm.list, or using sem() instead.",
          lav_msg_view(mm_list[[mm]], "none")))
      }
      keep_idx <- ptm_free[ptm_idx]
      sigma_11_1[par_idx, par_idx] <-
        sigma_11[keep_idx, keep_idx, drop = FALSE]
    }
  } # measurement block

  # only keep 'measurement part' parameters in Sigma.11
  if (!lavoptions$se %in% c("none", "bootstrap")) {
    sigma_11_1 <- sigma_11_1[step1_free_idx, step1_free_idx, drop = FALSE]
  } else {
    sigma_11_1 <- NULL
  }

  # create STEP1 list
  step1 <- list(
    MM.FIT = mm_fit, Sigma.11 = sigma_11_1,
    step1.free.idx = step1_free_idx,
    block.mm.idx = block_mm_idx,
    block.ptm.idx = block_ptm_idx,
    PT.free = pt_free,
    mm.list = mm_list, PT = pt_1
  )

  step1
}
