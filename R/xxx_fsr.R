# factor score regression

# four methods:
#  - naive (regression or Bartlett)
#  - Skrondal & Laake (2001) (regression models only)
#  - Croon (2002) (general + robust SE)
#  - simple: always use Bartlett, replace var(f) by psi estimate
#
# TODO
#  - Hishino & Bentler: this is simple + WLS

# changes 09 dec 2018: add analytic SE ('standard')
#                      make this the new default

fsr <- function(model = NULL,
                data = NULL,
                cmd = "sem",
                fsr.method = "Croon",
                fs.method = "Bartlett",
                fs.scores = FALSE,
                mm.options = list(se = "standard", test = "standard"),
                Gamma.NT = TRUE,
                lvinfo = FALSE,
                mm.list = NULL,
                ...,
                output = "lavaan") {
  # we need full data
  if (is.null(data)) {
    lav_msg_stop(gettext("full data is required for factor score regression"))
  }

  # check fsr.method argument
  fsr.method <- tolower(fsr.method)
  if (fsr.method == "naive") {
    # nothing to do
  } else if (fsr.method %in% c(
    "skrondal", "laake", "skrondallaake",
    "skrondal.laake", "skrondal-laake"
  )) {
    fsr.method <- "skrondal.laake"
  } else if (fsr.method == "croon") {
    # nothing to do
  } else if (fsr.method == "simple") {
    # force fs.method to Bartlett!
    fs.method <- "Bartlett"
  } else {
    lav_msg_stop(gettext("invalid option for argument fsr.method:"),
      fsr.method)
  }

  # check fs.method argument
  fs.method <- tolower(fs.method)
  if (fs.method %in% c("bartlett", "barttlett", "bartlet")) {
    fs.method <- "Bartlett"
  } else if (fs.method == "regression") {
    # nothing to do
  } else {
    lav_msg_stop(gettext("invalid option for argument fs.method:"),
      fs.method
    )
  }

  if (output %in% c("scores", "fs.scores", "fsr.scores")) {
    fs.scores <- TRUE
  }

  # dot dot dot
  dotdotdot <- list(...)

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
  FIT <- suppressWarnings(do.call(cmd,
    args = c(list(
      model = model,
      data = data,
      # meanstructure = TRUE,
      do.fit = FALSE
    ), dotdotdot0)
  ))
  lavoptions <- lavInspect(FIT, "options")
  # restore
  lavoptions$se <- dotdotdot$se
  lavoptions$test <- dotdotdot$test
  ngroups <- lavInspect(FIT, "ngroups")
  lavpta <- FIT@pta
  lavpartable <- lav_partable_set_cache(FIT@ParTable, lavpta)

  # FIXME: not ready for multiple groups yet
  if (ngroups > 1L) {
    lav_msg_stop(gettext("fsr code not ready for multiple groups (yet)"))
  }

  # if missing = "listwise", make data complete
  if (lavoptions$missing == "listwise") {
    # FIXME: make this work for multiple groups!!
    OV <- unique(unlist(lavpta$vnames$ov))
    data <- na.omit(data[, OV])
  }

  # any `regular' latent variables?
  lv.names <- unique(unlist(FIT@pta$vnames$lv.regular))
  ov.names <- unique(unlist(FIT@pta$vnames$ov))

  # check for higher-order factors
  good.idx <- logical(length(lv.names))
  for (f in seq_len(length(lv.names))) {
    # check the indicators
    FAC <- lv.names[f]
    IND <- lavpartable$rhs[lavpartable$lhs == FAC &
      lavpartable$op == "=~"]
    if (all(IND %in% ov.names)) {
      good.idx[f] <- TRUE
    }
    # FIXME: check for mixed lv/ov indicators
  }
  lv.names <- lv.names[good.idx]

  if (length(lv.names) == 0L) {
    lav_msg_stop(gettext("model does not contain any (measured) latent variables"))
  }
  nfac <- length(lv.names)

  # check parameter table
  PT <- lav_partable_set_cache(parTable(FIT))
  PT$est <- PT$se <- NULL

  # extract structural part
  PT.PA <- lav_partable_subset_structural_model(PT)

  # check if we can use skrondal & laake (no mediational terms?)
  if (fsr.method == "skrondal.laake") {
    # determine eqs.y and eqs.x names
    eqs.x.names <- unlist(FIT@pta$vnames$eqs.x)
    eqs.y.names <- unlist(FIT@pta$vnames$eqs.y)
    eqs.names <- unique(c(eqs.x.names, eqs.y.names))
    if (any(eqs.x.names %in% eqs.y.names)) {
      lav_msg_stop(
        gettextf("mediational relationships are not allowed for the
                 Skrondal.Laake method; use %s instead.",
                 dQuote("Croon")))
    }
  }


  # STEP 1a: compute factor scores for each measurement model (block)

  # how many measurement models?
  if (!is.null(mm.list)) {
    if (fsr.method != "simple") {
      lav_msg_stop(gettext("mm.list only available if fsr.method = \"simple\""))
    }

    nblocks <- length(mm.list)
    # check each measurement block
    for (b in seq_len(nblocks)) {
      if (!all(mm.list[[b]] %in% lv.names)) {
        lav_msg_stop(
          gettextf("mm.list contains unknown latent variable(s): %s",
            lav_msg_view(mm.list[[b]][mm.list[[b]] %in% lv.names],
                         log.sep = "none")))
      }
    }
  } else {
    # TODO: here comes the automatic 'detection' of linked
    #       measurement models
    #
    # for now we take a single latent variable per measurement model block
    mm.list <- as.list(lv.names)
    nblocks <- length(mm.list)
  }

  # compute factor scores, per latent variable
  FS.SCORES <- vector("list", length = ngroups)
  LVINFO <- vector("list", length = ngroups)
  if (ngroups > 1L) {
    names(FS.SCORES) <- names(LVINFO) <- lavInspect(FIT, "group.label")
  }
  for (g in 1:ngroups) {
    FS.SCORES[[g]] <- vector("list", length = nblocks)
    # names(FS.SCORES[[g]]) <-  lv.names
    LVINFO[[g]] <- vector("list", length = nblocks)
    # names(LVINFO[[g]]) <- lv.names
  }

  # adjust options
  dotdotdot2 <- dotdotdot
  dotdotdot2$se <- "none"
  dotdotdot2$test <- "none"
  dotdotdot2$debug <- FALSE
  dotdotdot2$verbose <- FALSE
  dotdotdot2$auto.cov.lv.x <- TRUE # allow correlated exogenous factors

  # override with mm.options
  dotdotdot2 <- modifyList(dotdotdot2, mm.options)

  # we assume the same number/names of lv's per group!!!
  MM.FIT <- vector("list", nblocks)
  Sigma2.block <- vector("list", nblocks)
  for (b in 1:nblocks) {
    # create parameter table for this measurement block only
    PT.block <-
      lav_partable_subset_measurement_model(
        PT = PT,
        add.lv.cov = TRUE,
        lv.names = mm.list[[b]]
      )
    # fit 1-factor model
    fit.block <- do.call("lavaan",
      args = c(list(
        model = PT.block,
        data = data
      ), dotdotdot2)
    )
    # check convergence
    if (!lavInspect(fit.block, "converged")) {
      lav_msg_stop(
        gettextf("measurement model for %s did not converge.",
        lav_msg_view(mm.list[[b]]))
      )
    }
    # store fitted measurement model
    MM.FIT[[b]] <- fit.block

    # fs.method?
    if (fsr.method == "skrondal.laake") {
      # dependent -> Bartlett
      if (lv.names[b] %in% eqs.y.names) {
        fs.method <- "Bartlett"
      } else {
        fs.method <- "regression"
      }
    }

    # compute factor scores
    SC <- lavPredict(fit.block, method = fs.method, fsm = TRUE)
    FSM <- attr(SC, "fsm")
    attr(SC, "fsm") <- NULL

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

    LAMBDA <- computeLAMBDA(fit.block@Model) # FIXME: remove dummy lv's?
    THETA <- computeTHETA(fit.block@Model) # FIXME: remove not used ov?
    PSI <- computeVETA(fit.block@Model)

    # if ngroups = 1, make list again
    if (ngroups == 1L) {
      # because lavPredict() drops the list
      SC <- list(SC)
    }

    # store results
    for (g in 1:ngroups) {
      FS.SCORES[[g]][[b]] <- SC[[g]]
      if (fsr.method %in% c("croon", "simple")) {
        offset <- FSM[[g]] %*% THETA[[g]] %*% t(FSM[[g]])
        scale <- FSM[[g]] %*% LAMBDA[[g]]
        scale.inv <- solve(scale)
        scoffset <- scale.inv %*% offset %*% scale.inv

        LVINFO[[g]][[b]] <- list(
          lv.names = mm.list[[b]],
          fsm = FSM[[g]],
          lambda = LAMBDA[[g]],
          psi = PSI[[g]],
          theta = THETA[[g]],
          offset = offset,
          scale = scale,
          scale.inv = scale.inv,
          scoffset = scoffset
        )
      }
    } # g

    # Delta.21: list per group
    Delta.21 <- lav_fsr_delta21(fit.block, FSM)

    # vcov
    Sigma1.block <- vcov(fit.block)
    tmp <- matrix(0, nrow(Delta.21[[1]]), nrow(Delta.21[[1]]))
    lavsamplestats <- fit.block@SampleStats
    for (g in 1:ngroups) {
      fg <- lavsamplestats@nobs[[g]] / lavsamplestats@ntotal
      tmp <-
        tmp + fg * (Delta.21[[g]] %*% Sigma1.block %*% t(Delta.21[[g]]))
    }
    Sigma2.block[[b]] <- tmp
  } # measurement block

  # Sigma.2 = Delta.21 %*% Sigma.1  %*% t(Delta.21)
  Sigma.2 <- lav_matrix_bdiag(Sigma2.block)


  # compute empirical covariance matrix factor scores + observed variables
  # in structural part
  group.values <- lav_partable_group_values(PT.PA)
  FS.COV <- vector("list", length = ngroups)
  FSR.COV <- vector("list", length = ngroups)
  FSR.COV2 <- vector("list", length = ngroups)
  Y <- vector("list", length = ngroups)
  if (lavoptions$meanstructure) {
    FS.MEAN <- vector("list", length = ngroups)
  } else {
    FS.MEAN <- NULL
  }
  for (g in seq_len(ngroups)) {
    # full data for structural model
    struc.names <- lavNames(PT.PA, "ov", group = group.values[g])
    # reorder struc.names, so that order is the same as in MM (new in 0.6-9)
    lv.idx <- which(struc.names %in% lv.names)
    struc.names[lv.idx] <- lv.names

    struc.ov.idx <- which(!struc.names %in% lv.names)
    struc.lv.idx <- which(struc.names %in% lv.names)
    lv.order <- match(lv.names, struc.names)
    if (length(struc.ov.idx) > 0L) {
      ov.idx <- which(FIT@Data@ov.names[[g]] %in%
        struc.names[struc.ov.idx])
      Y.g <- matrix(0,
        nrow = nrow(FS.SCORES[[g]][[1]]),
        ncol = length(struc.names)
      )
      Y.g[, struc.lv.idx] <- do.call(
        "cbind",
        FS.SCORES[[g]]
      )[, lv.order, drop = FALSE]
      Y.g[, struc.ov.idx] <- FIT@Data@X[[g]][, ov.idx, drop = FALSE]
    } else {
      Y.g <- do.call("cbind", FS.SCORES[[g]])[, lv.order, drop = FALSE]
    }
    Y[[g]] <- Y.g

    # sample statistics for structural model
    COV <- cov(Y.g) # divided by N-1
    if (lavoptions$likelihood == "normal") {
      Ng <- lavInspect(FIT, "nobs")[g]
      COV <- COV * (Ng - 1) / Ng
    }
    FS.COV[[g]] <- COV

    if (lavoptions$meanstructure) {
      FS.MEAN[[g]] <- colMeans(Y.g)
    }

    # STEP 1b: if using `Croon' method: correct COV matrix:
    if (fsr.method %in% c("croon")) {
      scoffset <- lav_matrix_bdiag(lapply(LVINFO[[g]], "[[", "scoffset"))
      scale.inv <- lav_matrix_bdiag(lapply(LVINFO[[g]], "[[", "scale.inv"))

      SCOFFSET <- matrix(0,
        nrow = length(struc.names),
        ncol = length(struc.names)
      )
      SCOFFSET[struc.lv.idx, struc.lv.idx] <- scoffset

      SCALE.INV <- diag(length(struc.names))
      SCALE.INV[struc.lv.idx, struc.lv.idx] <- scale.inv

      FSR.COV[[g]] <- SCALE.INV %*% FS.COV[[g]] %*% SCALE.INV - SCOFFSET
    } else if (fsr.method == "simple") {
      psi <- lav_matrix_bdiag(lapply(LVINFO[[g]], "[[", "psi"))

      FSR.COV[[g]] <- FS.COV[[g]]
      # scalar version only (for now)
      diag(FSR.COV[[g]])[struc.lv.idx] <- psi
    } else {
      FSR.COV[[g]] <- FS.COV[[g]]
    }

    # copy with different labels
    FSR.COV2[[g]] <- FSR.COV[[g]]

    # add row/col names
    rownames(FS.COV[[g]]) <- colnames(FS.COV[[g]]) <- struc.names
    rownames(FSR.COV[[g]]) <- colnames(FSR.COV[[g]]) <- struc.names
    rownames(FSR.COV2[[g]]) <- colnames(FSR.COV2[[g]]) <- struc.names
    rownames(FSR.COV2[[g]])[struc.lv.idx] <-
      colnames(FSR.COV2[[g]])[struc.lv.idx] <-
      paste(lv.names, ".si", sep = "")

    # check if FSR.COV is positive definite for all groups
    eigvals <- eigen(FSR.COV[[g]], symmetric = TRUE, only.values = TRUE)$values
    if (any(eigvals < .Machine$double.eps^(3 / 4))) {
      lav_msg_stop(gettext(
        "corrected covariance matrix of factor scores is not positive definite"),
        if (ngroups > 1L) gettextf("in group %s", g) else ""
      )
    }
  } # g


  # STEP 1c: do we need full set of factor scores?
  if (fs.scores) {
    # transform?
    if (fsr.method %in% c("croon", "simple")) {
      for (g in 1:ngroups) {
        OLD.inv <- solve(FS.COV[[g]])
        OLD.inv.sqrt <- lav_matrix_symmetric_sqrt(OLD.inv)
        FSR.COV.sqrt <- lav_matrix_symmetric_sqrt(FSR.COV[[g]])
        SC <- as.matrix(Y[[g]])
        SC <- SC %*% OLD.inv.sqrt %*% FSR.COV.sqrt
        SC <- as.data.frame(SC)
        names(SC) <- lv.names
        Y[[g]] <- SC
      }
    }

    # unlist if multiple groups, add group column
    if (ngroups == 1L) {
      Y <- as.data.frame(Y[[1]])
    } else {
      lav_msg_fixme("fix this!")
    }
  }




  # STEP 2: fit structural model using (corrected?) factor scores

  # free all means/intercepts (of observed variables only)
  lv.names.pa <- lavNames(PT.PA, "lv")
  int.idx <- which(PT.PA$op == "~1" & !PT.PA$lhs %in% lv.names.pa)
  PT.PA$free[int.idx] <- 1L
  PT.PA$free[PT.PA$free > 0L] <- seq_len(sum(PT.PA$free > 0L))
  PT.PA$ustart[int.idx] <- NA



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
  fit <- lavaan(PT.PA,
    sample.cov = FSR.COV,
    sample.mean = FS.MEAN,
    sample.nobs = FIT@SampleStats@nobs,
    slotOptions = lavoptions2
  )

  # only to correct the SE, we create another model, augmented with
  # the croon parameters
  PT.PA2 <- parTable(fit)
  PT.si <- lav_fsr_pa2si(PT.PA2, LVINFO = LVINFO)
  idx1 <- PT.si$free[PT.si$user == 10L & PT.si$free > 0L]
  idx2 <- PT.si$free[PT.si$user != 10L & PT.si$free > 0L]

  lavoptions3 <- lavoptions2
  lavoptions3$optim.method <- "none"
  lavoptions3$test <- "standard"
  lavoptions3$se <- "none"
  lavoptions3$check.gradient <- FALSE
  lavoptions3$information <- "expected" ## FIXME: lav_model_gradient + delta
  fit.si2 <- lavaan(PT.si,
    sample.cov  = FSR.COV2,
    sample.mean = FS.MEAN,
    sample.nobs = FIT@SampleStats@nobs,
    slotOptions = lavoptions3
  )
  Info.all <- lavTech(fit.si2, "information") * nobs(fit)
  I33 <- Info.all[idx2, idx2]
  I32 <- Info.all[idx2, idx1]
  I23 <- Info.all[idx1, idx2]
  I22 <- Info.all[idx1, idx1]

  I33.inv <- lav_matrix_symmetric_inverse(I33)

  V1 <- I33.inv
  V2 <- I33.inv %*% I32 %*% Sigma.2 %*% t(I32) %*% I33.inv
  VCOV <- V1 + V2

  # fill in standard errors step 2
  PT.PA2$se[PT.PA2$free > 0L] <- sqrt(diag(VCOV))

  if (output == "lavaan" || output == "fsr") {
    lavoptions3$se <- "twostep"
    fit <- lavaan::lavaan(PT.PA2,
      sample.cov = FSR.COV,
      sample.mean = FS.MEAN,
      sample.nobs = FIT@SampleStats@nobs,
      slotOptions = lavoptions3
    )
    fit@vcov$vcov <- VCOV
  }

  # extra info
  extra <- list(
    FS.COV = FS.COV, FS.SCORES = Y,
    FSR.COV = FSR.COV,
    LVINFO = LVINFO, Sigma.2 = Sigma.2
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
    HEADER <- paste("This is fsr (0.2) -- factor score regression using ",
      "fsr.method = ", fsr.method,
      sep = ""
    )
    out <- list(header = HEADER, MM.FIT = MM.FIT, STRUC.FIT = fit)
    if (lvinfo) {
      out$lvinfo <- extra
    }

    class(out) <- c("lavaan.fsr", "list")
  } else if (output %in% c("lavaan", "fit")) {
    out <- fit
  } else if (output == "extra") {
    out <- extra
  } else if (output == "lvinfo") {
    out <- LVINFO
  } else if (output %in% c("scores", "f.scores", "fs.scores")) {
    out <- Y
  } else if (output %in% c(
    "FSR.COV", "fsr.cov", "croon", "cov.croon",
    "croon.cov", "COV", "cov"
  )) {
    out <- FSR.COV
  } else if (output %in% c("FS.COV", "fs.cov")) {
    out <- FS.COV
  } else {
    lav_msg_stop(gettext("unknown output= argument:"), output)
  }

  out
}
