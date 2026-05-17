# factor score regression

# four methods:
#  - naive (regression or Bartlett)
#  - Skrondal & Laake (2001) (regression models only)
#  - Croon (2002) (general + robust SE)
#  - simple: always use Bartlett, replace var(f) by psi estimate
#
# TODO:
#  - Hishino & Bentler: this is simple + WLS

# changes 09 dec 2018: add analytic SE ('standard')
#                      make this the new default

fsr <- function(model = NULL,
                data = NULL,
                cmd = "sem",
                fsr_method = "Croon",
                fs_method = "Bartlett",
                fs_scores = FALSE,
                mm_options = list(se = "standard", test = "standard"),
                gamma_nt = TRUE,
                lvinfo = FALSE,
                mm_list = NULL,
                ...,
                output = "lavaan") {
  # we need full data
  if (is.null(data)) {
    lav_msg_stop(gettext("full data is required for factor score regression"))
  }
  # dot dot dot
  dotdotdot <- list(...)
  # ------------- handling of warn/debug/verbose switches ----------
  if (!is.null(dotdotdot$debug)) {
    current_debug <- lav_debug()
    if (lav_debug(dotdotdot$debug))
      on.exit(lav_debug(current_debug), TRUE)
    dotdotdot$debug <- NULL
    if (lav_debug()) {
      dotdotdot$warn <- TRUE       # force warnings if debug
      dotdotdot$verbose <- TRUE    # force verbose if debug
    }
  }
  if (!is.null(dotdotdot$warn)) {
    current_warn <- lav_warn()
    if (lav_warn(dotdotdot$warn))
      on.exit(lav_warn(current_warn), TRUE)
    dotdotdot$warn <- NULL
  }
  if (!is.null(dotdotdot$verbose)) {
    current_verbose <- lav_verbose()
    if (lav_verbose(dotdotdot$verbose))
      on.exit(lav_verbose(current_verbose), TRUE)
    dotdotdot$verbose <- NULL
  }

  # check fsr.method argument
  fsr_method <- tolower(fsr_method)
  if (fsr_method == "naive") {
    # nothing to do
  } else if (fsr_method %in% c(
    "skrondal", "laake", "skrondallaake",
    "skrondal.laake", "skrondal-laake"
  )) {
    fsr_method <- "skrondal.laake"
  } else if (fsr_method == "croon") {
    # nothing to do
  } else if (fsr_method == "simple") {
    # force fs.method to Bartlett!
    fs_method <- "Bartlett"
  } else {
    lav_msg_stop(gettext("invalid option for argument fsr.method:"),
      fsr_method)
  }

  # check fs.method argument
  fs_method <- tolower(fs_method)
  if (fs_method %in% c("bartlett", "barttlett", "bartlet")) {
    fs_method <- "Bartlett"
  } else if (fs_method == "regression") {
    # nothing to do
  } else {
    lav_msg_stop(gettext("invalid option for argument fs.method:"),
      fs_method
    )
  }

  if (output %in% c("scores", "fs.scores", "fsr.scores")) {
    fs_scores <- TRUE
  }

  # change 'default' values for fsr
  if (is.null(dotdotdot$se)) {
    dotdotdot$se <- "standard"
  }
  if (is.null(dotdotdot$test)) {
    dotdotdot$test <- "standard"
  }
  if (is.null(dotdotdot$missing)) {
    dotdotdot$missing <- "ml"
  }
  if (is.null(dotdotdot$meanstructure)) {
    dotdotdot$meanstructure <- TRUE
  }


  # STEP 0: process full model, without fitting
  dotdotdot0 <- dotdotdot
  dotdotdot0$do.fit <- NULL
  dotdotdot0$se <- "none" # to avoid warning about missing="listwise"
  dotdotdot0$test <- "none" # to avoid warning about missing="listwise"

  # check for arguments that we do not want (eg sample.cov)?
  # TODO

  # initial processing of the model, no fitting
  fit_1 <- suppressWarnings(do.call(cmd,
    args = c(list(
      model = model,
      data = data,
      # meanstructure = TRUE,
      do.fit = FALSE
    ), dotdotdot0)
  ))
  lavoptions <- lavInspect(fit_1, "options")
  # restore
  lavoptions$se <- dotdotdot$se
  lavoptions$test <- dotdotdot$test
  ngroups <- lavInspect(fit_1, "ngroups")
  lavpta <- fit_1@pta
  lavpartable <- lav_partable_set_cache(fit_1@ParTable, lavpta)

  # FIXME: not ready for multiple groups yet
  if (ngroups > 1L) {
    lav_msg_stop(gettext("fsr code not ready for multiple groups (yet)"))
  }

  # if missing = "listwise", make data complete
  if (lavoptions$missing == "listwise") {
    # FIXME: make this work for multiple groups!!
    ov <- unique(unlist(lavpta$vnames$ov))
    data <- na.omit(data[, ov])
  }

  # any `regular' latent variables?
  lv_names <- unique(unlist(fit_1@pta$vnames$lv.regular))
  ov_names <- unique(unlist(fit_1@pta$vnames$ov))

  # check for higher-order factors
  good_idx <- logical(length(lv_names))
  for (f in seq_along(lv_names)) {
    # check the indicators
    fac <- lv_names[f]
    ind <- lavpartable$rhs[lavpartable$lhs == fac &
      lavpartable$op == "=~"]
    if (all(ind %in% ov_names)) {
      good_idx[f] <- TRUE
    }
    # FIXME: check for mixed lv/ov indicators
  }
  lv_names <- lv_names[good_idx]

  if (length(lv_names) == 0L) {
    lav_msg_stop(gettext(
      "model does not contain any (measured) latent variables"))
  }
  # nfac <- length(lv_names)

  # check parameter table
  pt_1 <- lav_partable_set_cache(parTable(fit_1))
  pt_1$est <- pt_1$se <- NULL

  # extract structural part
  pt_pa <- lav_partable_subset_structural_model(pt_1)

  # check if we can use skrondal & laake (no mediational terms?)
  if (fsr_method == "skrondal.laake") {
    # determine eqs.y and eqs.x names
    eqs_x_names <- unlist(fit_1@pta$vnames$eqs.x)
    eqs_y_names <- unlist(fit_1@pta$vnames$eqs.y)
    # eqs_names <- unique(c(eqs_x_names, eqs_y_names))
    if (any(eqs_x_names %in% eqs_y_names)) {
      lav_msg_stop(
        gettextf("mediational relationships are not allowed for the
                 Skrondal.Laake method; use %s instead.",
                 dQuote("Croon")))
    }
  }


  # STEP 1a: compute factor scores for each measurement model (block)

  # how many measurement models?
  if (!is.null(mm_list)) {
    if (fsr_method != "simple") {
      lav_msg_stop(gettext("mm.list only available if fsr.method = \"simple\""))
    }

    nblocks <- length(mm_list)
    # check each measurement block
    for (b in seq_len(nblocks)) {
      if (!all(mm_list[[b]] %in% lv_names)) {
        lav_msg_stop(
          gettextf("mm.list contains unknown latent variable(s): %s",
            lav_msg_view(mm_list[[b]][mm_list[[b]] %in% lv_names],
                         log_sep = "none")))
      }
    }
  } else {
    # TODO: here comes the automatic 'detection' of linked
    #       measurement models
    #
    # for now we take a single latent variable per measurement model block
    mm_list <- as.list(lv_names)
    nblocks <- length(mm_list)
  }

  # compute factor scores, per latent variable
  fs_scores_1 <- vector("list", length = ngroups)
  lvinfo_1 <- vector("list", length = ngroups)
  if (ngroups > 1L) {
    names(fs_scores_1) <- names(lvinfo_1) <- lavInspect(fit_1, "group.label")
  }
  for (g in 1:ngroups) {
    fs_scores_1[[g]] <- vector("list", length = nblocks)
    # names(FS.SCORES[[g]]) <-  lv.names
    lvinfo_1[[g]] <- vector("list", length = nblocks)
    # names(LVINFO[[g]]) <- lv.names
  }

  # adjust options
  dotdotdot2 <- dotdotdot
  dotdotdot2$se <- "none"
  dotdotdot2$test <- "none"
  dotdotdot2$debug <- FALSE    # only transmitted to lavaan call = ok
  dotdotdot2$verbose <- FALSE  # only transmitted to lavaan call = ok
  dotdotdot2$auto.cov.lv.x <- TRUE # allow correlated exogenous factors

  # override with mm.options
  dotdotdot2 <- modifyList(dotdotdot2, mm_options)

  # we assume the same number/names of lv's per group!!!
  mm_fit <- vector("list", nblocks)
  sigma2_block <- vector("list", nblocks)
  for (b in 1:nblocks) {
    # create parameter table for this measurement block only
    pt_block <-
      lav_partable_subset_measurement_model(
        pt_1 = pt_1,
        add_lv_cov = TRUE,
        lv_names = mm_list[[b]]
      )
    # fit 1-factor model
    fit_block <- do.call("lavaan",
      args = c(list(
        model = pt_block,
        data = data
      ), dotdotdot2)
    )
    # check convergence
    if (!lavInspect(fit_block, "converged")) {
      lav_msg_stop(
        gettextf("measurement model for %s did not converge.",
        lav_msg_view(mm_list[[b]]))
      )
    }
    # store fitted measurement model
    mm_fit[[b]] <- fit_block

    # fs.method?
    if (fsr_method == "skrondal.laake") {
      # dependent -> Bartlett
      if (lv_names[b] %in% eqs_y_names) {
        fs_method <- "Bartlett"
      } else {
        fs_method <- "regression"
      }
    }

    # compute factor scores
    sc <- lavPredict(fit_block, method = fs_method, fsm = TRUE)
    fsm <- attr(sc, "fsm")
    attr(sc, "fsm") <- NULL

    # warning, FSM may be a list per pattern!
    # if(fit.block@Options$missing == "ml") {
    #    # do something...
    #    ngroups <- fit.block@Data@ngroups
    #    FSM.missing <- FSM
    #    FSM <- vector("list", length = "ngroups")
    #    for(g in seq_len(ngroups)) {
    #
    #    }
    # }

    mm_lambda <- lav_model_lambda(fit_block@Model) # FIXME: remove dummy lv's?
    mm_theta <- lav_model_theta(fit_block@Model) # FIXME: remove not used ov?
    mm_psi <- lav_model_veta(fit_block@Model)

    # if ngroups = 1, make list again
    if (ngroups == 1L) {
      # because lavPredict() drops the list
      sc <- list(sc)
    }

    # store results
    for (g in 1:ngroups) {
      fs_scores_1[[g]][[b]] <- sc[[g]]
      if (fsr_method %in% c("croon", "simple")) {
        offset <- fsm[[g]] %*% mm_theta[[g]] %*% t(fsm[[g]])
        scale <- fsm[[g]] %*% mm_lambda[[g]]
        scale_inv <- solve(scale)
        scoffset <- scale_inv %*% offset %*% scale_inv

        lvinfo_1[[g]][[b]] <- list(
          lv.names = mm_list[[b]],
          fsm = fsm[[g]],
          lambda = mm_lambda[[g]],
          psi = mm_psi[[g]],
          theta = mm_theta[[g]],
          offset = offset,
          scale = scale,
          scale.inv = scale_inv,
          scoffset = scoffset
        )
      }
    } # g

    # Delta.21: list per group
    delta_21 <- lav_fsr_delta21(fit_block, fsm)

    # vcov
    sigma1_block <- vcov(fit_block)
    tmp <- matrix(0, nrow(delta_21[[1]]), nrow(delta_21[[1]]))
    lavsamplestats <- fit_block@SampleStats
    for (g in 1:ngroups) {
      fg <- lavsamplestats@nobs[[g]] / lavsamplestats@ntotal
      tmp <-
        tmp + fg * (delta_21[[g]] %*% sigma1_block %*% t(delta_21[[g]]))
    }
    sigma2_block[[b]] <- tmp
  } # measurement block

  # Sigma.2 = Delta.21 %*% Sigma.1  %*% t(Delta.21)
  sigma_2 <- lav_matrix_bdiag(sigma2_block)


  # compute empirical covariance matrix factor scores + observed variables
  # in structural part
  group_values <- lav_partable_group_values(pt_pa)
  fs_cov <- vector("list", length = ngroups)
  fsr_cov <- vector("list", length = ngroups)
  fsr_cov2 <- vector("list", length = ngroups)
  y <- vector("list", length = ngroups)
  if (lavoptions$meanstructure) {
    fs_mean <- vector("list", length = ngroups)
  } else {
    fs_mean <- NULL
  }
  for (g in seq_len(ngroups)) {
    # full data for structural model
    struc_names <- lav_object_vnames(pt_pa, "ov", group = group_values[g])
    # reorder struc.names, so that order is the same as in MM (new in 0.6-9)
    lv_idx <- which(struc_names %in% lv_names)
    struc_names[lv_idx] <- lv_names

    struc_ov_idx <- which(!struc_names %in% lv_names)
    struc_lv_idx <- which(struc_names %in% lv_names)
    lv_order <- match(lv_names, struc_names)
    if (length(struc_ov_idx) > 0L) {
      ov_idx <- which(fit_1@Data@ov.names[[g]] %in%
        struc_names[struc_ov_idx])
      y_g <- matrix(0,
        nrow = nrow(fs_scores_1[[g]][[1]]),
        ncol = length(struc_names)
      )
      y_g[, struc_lv_idx] <- do.call(
        "cbind",
        fs_scores_1[[g]]
      )[, lv_order, drop = FALSE]
      y_g[, struc_ov_idx] <- fit_1@Data@X[[g]][, ov_idx, drop = FALSE]
    } else {
      y_g <- do.call("cbind", fs_scores_1[[g]])[, lv_order, drop = FALSE]
    }
    y[[g]] <- y_g

    # sample statistics for structural model
    cov_1 <- cov(y_g) # divided by N-1
    if (lavoptions$likelihood == "normal") {
      ng <- lavInspect(fit_1, "nobs")[g]
      cov_1 <- cov_1 * (ng - 1) / ng
    }
    fs_cov[[g]] <- cov_1

    if (lavoptions$meanstructure) {
      fs_mean[[g]] <- colMeans(y_g)
    }

    # STEP 1b: if using `Croon' method: correct COV matrix:
    if (fsr_method %in% c("croon")) {
      scoffset <- lav_matrix_bdiag(lapply(lvinfo_1[[g]], "[[", "scoffset"))
      scale_inv <- lav_matrix_bdiag(lapply(lvinfo_1[[g]], "[[", "scale.inv"))

      scoffset_1 <- matrix(0,
        nrow = length(struc_names),
        ncol = length(struc_names)
      )
      scoffset_1[struc_lv_idx, struc_lv_idx] <- scoffset

      scale_inv_1 <- diag(length(struc_names))
      scale_inv_1[struc_lv_idx, struc_lv_idx] <- scale_inv

      fsr_cov[[g]] <- scale_inv_1 %*% fs_cov[[g]] %*% scale_inv_1 - scoffset_1
    } else if (fsr_method == "simple") {
      psi <- lav_matrix_bdiag(lapply(lvinfo_1[[g]], "[[", "psi"))

      fsr_cov[[g]] <- fs_cov[[g]]
      # scalar version only (for now)
      diag(fsr_cov[[g]])[struc_lv_idx] <- psi
    } else {
      fsr_cov[[g]] <- fs_cov[[g]]
    }

    # copy with different labels
    fsr_cov2[[g]] <- fsr_cov[[g]]

    # add row/col names
    rownames(fs_cov[[g]]) <- colnames(fs_cov[[g]]) <- struc_names
    rownames(fsr_cov[[g]]) <- colnames(fsr_cov[[g]]) <- struc_names
    rownames(fsr_cov2[[g]]) <- colnames(fsr_cov2[[g]]) <- struc_names
    rownames(fsr_cov2[[g]])[struc_lv_idx] <-
      colnames(fsr_cov2[[g]])[struc_lv_idx] <-
      paste(lv_names, ".si", sep = "")

    # check if FSR.COV is positive definite for all groups
    eigvals <- eigen(fsr_cov[[g]], symmetric = TRUE, only.values = TRUE)$values
    if (any(eigvals < .Machine$double.eps^(3 / 4))) {
      if (ngroups > 1L) {
        lav_msg_stop(gettextf("corrected covariance matrix of factor scores
                               is not positive definite in group %s", g))

      } else {
        lav_msg_stop(gettext("corrected covariance matrix of factor scores
                              is not positive definite"))
      }
    }
  } # g


  # STEP 1c: do we need full set of factor scores?
  if (fs_scores) {
    # transform?
    if (fsr_method %in% c("croon", "simple")) {
      for (g in 1:ngroups) {
        old_inv <- solve(fs_cov[[g]])
        old_inv_sqrt <- lav_matrix_symmetric_sqrt(old_inv)
        fsr_cov_sqrt <- lav_matrix_symmetric_sqrt(fsr_cov[[g]])
        sc <- as.matrix(y[[g]])
        sc <- sc %*% old_inv_sqrt %*% fsr_cov_sqrt
        sc <- as.data.frame(sc)
        names(sc) <- lv_names
        y[[g]] <- sc
      }
    }

    # unlist if multiple groups, add group column
    if (ngroups == 1L) {
      y <- as.data.frame(y[[1]])
    } else {
      lav_msg_fixme("fix this!")
    }
  }




  # STEP 2: fit structural model using (corrected?) factor scores

  # free all means/intercepts (of observed variables only)
  lv_names_pa <- lav_object_vnames(pt_pa, "lv")
  int_idx <- which(pt_pa$op == "~1" & !pt_pa$lhs %in% lv_names_pa)
  pt_pa$free[int_idx] <- 1L
  pt_pa$free[pt_pa$free > 0L] <- seq_len(sum(pt_pa$free > 0L))
  pt_pa$ustart[int_idx] <- NA



  # adjust lavoptions
  if (is.null(dotdotdot$do.fit)) {
    lavoptions$do.fit <- TRUE
  } else {
    lavoptions$do.fit <- dotdotdot$do.fit
  }
  if (is.null(dotdotdot$se)) {
    lavoptions$se <- "standard"
  } else {
    lavoptions$se <- dotdotdot$se
  }
  if (is.null(dotdotdot$test)) {
    lavoptions$test <- "standard"
  } else {
    lavoptions$test <- dotdotdot$test
  }
  if (is.null(dotdotdot$sample.cov.rescale)) {
    lavoptions$sample.cov.rescale <- FALSE
  } else {
    lavoptions$sample.cov.rescale <- dotdotdot$sample.cov.rescale
  }

  # fit structural model -- point estimation ONLY
  lavoptions2 <- lavoptions
  # if(lavoptions$se == "standard") {
  #    lavoptions2$se <- "external"
  # }
  # lavoptions2$test <- "none"
  lavoptions2$se <- "none"
  lavoptions2$test <- "none"
  lavoptions2$missing <- "listwise" # always complete data anyway...
  fit <- lavaan(pt_pa,
    sample.cov = fsr_cov,
    sample.mean = fs_mean,
    sample.nobs = fit_1@SampleStats@nobs,
    slotOptions = lavoptions2
  )

  # only to correct the SE, we create another model, augmented with
  # the croon parameters
  pt_pa2 <- parTable(fit)
  pt_si <- lav_fsr_pa2si(pt_pa2, lvinfo = lvinfo_1)
  idx1 <- pt_si$free[pt_si$user == 10L & pt_si$free > 0L]
  idx2 <- pt_si$free[pt_si$user != 10L & pt_si$free > 0L]

  lavoptions3 <- lavoptions2
  lavoptions3$optim.method <- "none"
  lavoptions3$test <- "standard"
  lavoptions3$se <- "none"
  lavoptions3$check.gradient <- FALSE
  lavoptions3$information <- "expected" ## FIXME: lav_model_gradient + delta
  fit_si2 <- lavaan(pt_si,
    sample.cov  = fsr_cov2,
    sample.mean = fs_mean,
    sample.nobs = fit_1@SampleStats@nobs,
    slotOptions = lavoptions3
  )
  info_all <- lavTech(fit_si2, "information") * nobs(fit)
  i33 <- info_all[idx2, idx2]
  i32 <- info_all[idx2, idx1]
  # i23 <- info_all[idx1, idx2]
  # i22 <- info_all[idx1, idx1]

  i33_inv <- lav_matrix_symmetric_inverse(i33)

  v1 <- i33_inv
  v2 <- i33_inv %*% i32 %*% sigma_2 %*% t(i32) %*% i33_inv
  vcov_1 <- v1 + v2

  # fill in standard errors step 2
  pt_pa2$se[pt_pa2$free > 0L] <- sqrt(diag(vcov_1))

  if (output == "lavaan" || output == "fsr") {
    lavoptions3$se <- "twostep"
    fit <- lavaan::lavaan(pt_pa2,
      sample.cov = fsr_cov,
      sample.mean = fs_mean,
      sample.nobs = fit_1@SampleStats@nobs,
      slotOptions = lavoptions3
    )
    fit@vcov$vcov <- vcov_1
  }

  # extra info
  extra <- list(
    FS.COV = fs_cov, FS.SCORES = y,
    FSR.COV = fsr_cov,
    LVINFO = lvinfo_1, Sigma.2 = sigma_2
  )

  # standard errors
  # lavsamplestats <- fit@SampleStats
  # lavsamplestats@NACOV <- Omega.f
  # VCOV <- lav_model_vcov(fit@Model, lavsamplestats = lavsamplestats,
  #                       lavoptions = lavoptions)
  # SE <- lav_model_vcov_se(fit@Model, fit@ParTable, VCOV = VCOV)
  # PE$se <- SE
  # tmp.se <- ifelse(PE$se == 0.0, NA, PE$se)
  # zstat <- pvalue <- TRUE
  # if(zstat) {
  #    PE$z <- PE$est / tmp.se
  #    if(pvalue) {
  #        PE$pvalue <- 2 * (1 - pnorm( abs(PE$z) ))
  #    }
  # }

  if (output == "fsr") {
    header <- paste("This is fsr (0.2) -- factor score regression using ",
      "fsr.method = ", fsr_method,
      sep = ""
    )
    out <- list(header = header, MM.FIT = mm_fit, STRUC.FIT = fit)
    if (lvinfo) {
      out$lvinfo <- extra
    }

    class(out) <- c("lavaan.fsr", "list")
  } else if (output %in% c("lavaan", "fit")) {
    out <- fit
  } else if (output == "extra") {
    out <- extra
  } else if (output == "lvinfo") {
    out <- lvinfo_1
  } else if (output %in% c("scores", "f.scores", "fs.scores")) {
    out <- y
  } else if (output %in% c(
    "FSR.COV", "fsr.cov", "croon", "cov.croon",
    "croon.cov", "COV", "cov"
  )) {
    out <- fsr_cov
  } else if (output %in% c("FS.COV", "fs.cov")) {
    out <- fs_cov
  } else {
    lav_msg_stop(gettext("unknown output= argument:"), output)
  }

  out
}
