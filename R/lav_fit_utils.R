# utility functions needed to compute various (robust) fit measures:
#
# - lav_fit_catml_dwls (for 'robust' RMSEA/CFI if data is cateogrical)
# - lav_fit_fiml_corrected (correct RMSEA/CFI if data is incomplete)

# compute scaling-factor (c.hat3) for fit.dwls, using fit.catml ingredients
# see:
#     Savalei, V. (2021) Improving Fit Indices In SEM with categorical data.
#     Multivariate Behavioral Research, 56(3), 390-407.
#

# YR Dec 2022: first version
# YR Jan 2023: catml_dwls should check if the input 'correlation' matrix
#              is positive-definite (or not)

lav_fit_catml_dwls <- function(lavobject, check.pd = TRUE) {
  # empty list
  empty.list <- list(
    XX3 = as.numeric(NA), df3 = as.numeric(NA),
    c.hat3 = as.numeric(NA), XX3.scaled = as.numeric(NA),
    XX3.null = as.numeric(NA), df3.null = as.numeric(NA),
    c.hat3.null = as.numeric(NA)
  )

  # limitations
  if (!lavobject@Model@categorical ||
    lavobject@Options$conditional.x ||
    length(unlist(lavobject@pta$vnames$ov.num)) > 0L) {
    return(empty.list)
  } else {
    lavdata <- lavobject@Data
    lavsamplestats <- lavobject@SampleStats
  }

  # check if input matrix (or matrices) are all positive definite
  # (perhaps later, we can rely on 'smoothing', but not for now
  pd.flag <- TRUE
  if (check.pd) {
    for (g in seq_len(lavdata@ngroups)) {
      COR <- lavsamplestats@cov[[g]]
      ev <- eigen(COR, symmetric = TRUE, only.values = TRUE)$values
      if (any(ev < .Machine$double.eps^(1 / 2))) {
        # non-pd!
        pd.flag <- FALSE
        # should we give a warning here? (not for now)
        # warning("lavaan WARNING: robust RMSEA/CFI could not be computed because the input correlation matrix is not positive-definite")
        # what should we do? return NA (for now)
        return(empty.list)
      }
    }
  }

  # 'refit' using estimator = "catML"
  fit.catml <- try(lav_object_catml(lavobject), silent = TRUE)
  if (inherits(fit.catml, "try-error")) {
    return(empty.list)
  }

  XX3 <- fit.catml@test[[1]]$stat
  df3 <- fit.catml@test[[1]]$df


  # compute 'k'
  V <- lavTech(fit.catml, "wls.v") # NT-ML weight matrix

  W.dwls <- lavTech(lavobject, "wls.v") # DWLS weight matrix
  Gamma <- lavTech(lavobject, "gamma") # acov of polychorics
  Delta <- lavTech(lavobject, "delta")
  E.inv <- lavTech(lavobject, "inverted.information")

  fg <- unlist(lavsamplestats@nobs) / lavsamplestats@ntotal

  # Fixme: as we only need the trace, perhaps we could do this
  # group-specific? (see lav_test_satorra_bentler_trace_original)
  V.g <- V
  W.dwls.g <- W.dwls
  Gamma.f <- Gamma
  Delta.g <- Delta
  for (g in seq_len(lavdata@ngroups)) {
    ntotal <- nrow(Gamma[[g]])
    nvar <- lavobject@Model@nvar[[g]]
    pstar <- nvar * (nvar - 1) / 2
    rm.idx <- seq_len(ntotal - pstar)

    # reduce
    Delta.g[[g]] <- Delta[[g]][-rm.idx, , drop = FALSE]
    # reduce and weight
    W.dwls.g[[g]] <- fg[g] * W.dwls[[g]][-rm.idx, -rm.idx]
    V.g[[g]] <- fg[g] * V[[g]] # should already have the right dims
    Gamma.f[[g]] <- 1 / fg[g] * Gamma[[g]][-rm.idx, -rm.idx]
  }
  # create 'big' matrices
  W.dwls.all <- lav_matrix_bdiag(W.dwls.g)
  V.all <- lav_matrix_bdiag(V.g)
  Gamma.all <- lav_matrix_bdiag(Gamma.f)
  Delta.all <- do.call("rbind", Delta.g)

  # compute trace
  WiU.all <- diag(nrow(W.dwls.all)) - Delta.all %*% E.inv %*% t(Delta.all) %*% W.dwls.all
  ks <- sum(diag(t(WiU.all) %*% V.all %*% WiU.all %*% Gamma.all))

  # convert to lavaan 'scaling.factor'
  c.hat3 <- ks / df3
  XX3.scaled <- XX3 / c.hat3

  # baseline model
  XX3.null <- fit.catml@baseline$test[[1]]$stat
  if (is.null(XX3.null)) {
    XX3.null <- as.numeric(NA)
    df3.null <- as.numeric(NA)
    kbs <- as.numeric(NA)
    c.hat3.null <- as.numeric(NA)
  } else {
    df3.null <- fit.catml@baseline$test[[1]]$df
    kbs <- sum(diag(Gamma.all))
    c.hat3.null <- kbs / df3.null
  }

  # return values
  list(
    XX3 = XX3, df3 = df3, c.hat3 = c.hat3, XX3.scaled = XX3.scaled,
    XX3.null = XX3.null, df3.null = df3.null, c.hat3.null = c.hat3.null
  )
}


# compute ingredients to compute FIML-Corrected RMSEA/CFI
# see:
#     Zhang X, Savalei V. (2022). New computations for RMSEA and CFI
#     following FIML and TS estimation with missing data. Psychological Methods.

lav_fit_fiml_corrected <- function(lavobject, version = "V3") {
  version <- toupper(version)
  if (!version %in% c("V3", "V6")) {
    lav_msg_stop(gettext("only FIML-C(V3) and FIML-C(V6) are available."))
  }

  # empty list
  empty.list <- list(
    XX3 = as.numeric(NA), df3 = as.numeric(NA),
    c.hat3 = as.numeric(NA), XX3.scaled = as.numeric(NA),
    XX3.null = as.numeric(NA), df3.null = as.numeric(NA),
    c.hat3.null = as.numeric(NA)
  )

  # limitations
  if (lavobject@Options$conditional.x ||
    lavobject@Data@nlevels > 1L ||
    !.hasSlot(lavobject, "h1") ||
    is.null(lavobject@h1$implied$cov[[1]])) {
    return(empty.list)
  } else {
    lavdata <- lavobject@Data
    lavsamplestats <- lavobject@SampleStats

    h1 <- lavTech(lavobject, "h1", add.labels = TRUE)
    COV.tilde <- lapply(h1, "[[", "cov")
    MEAN.tilde <- lapply(h1, "[[", "mean")
    sample.nobs <- unlist(lavsamplestats@nobs)
  }

  # 'refit' using 'tilde' (=EM/saturated) sample statistics
  fit.tilde <- try(lavaan(
    model = parTable(lavobject),
    sample.cov = COV.tilde,
    sample.mean = MEAN.tilde,
    sample.nobs = sample.nobs,
    sample.cov.rescale = FALSE,
    information = "observed",
    optim.method = "none",
    se = "none",
    test = "standard",
    baseline = FALSE,
    check.post = FALSE
  ), silent = TRUE)
  if (inherits(fit.tilde, "try-error")) {
    return(empty.list)
  }

  XX3 <- fit.tilde@test[[1]]$stat
  df3 <- fit.tilde@test[[1]]$df

  # compute 'k'

  # V3/V6: always use h1.information = "unstructured"!!
  lavobject@Options$h1.information <- c("unstructured", "unstructured")
  lavobject@Options$observed.information <- c("h1", "h1")
  fit.tilde@Options$h1.information <- c("unstructured", "unstructured")

  Wm <- Wm.g <- lav_model_h1_information_observed(lavobject)
  Jm <- lav_model_h1_information_firstorder(lavobject)
  Wc <- Wc.g <- lav_model_h1_information_observed(fit.tilde)
  if (version == "V3") {
    Gamma.f <- vector("list", length = lavdata@ngroups)
  }
  Delta <- lavTech(lavobject, "delta")
  E.inv <- lavTech(lavobject, "inverted.information")
  # Wmi <- Wmi.g <- lapply(Wm, solve) ## <- how wrote this? (I did)
  Wmi <- Wmi.g <- try(lapply(Wm, lav_matrix_symmetric_inverse),
    silent = TRUE
  )
  if (inherits(Wmi, "try-error")) {
    return(empty.list)
  }

  fg <- unlist(lavsamplestats@nobs) / lavsamplestats@ntotal
  # Fixme: as we only need the trace, perhaps we could do this
  # group-specific? (see lav_test_satorra_bentler_trace_original)
  for (g in seq_len(lavdata@ngroups)) {
    # group weight
    Wc.g[[g]] <- fg[g] * Wc[[g]]
    Wm.g[[g]] <- fg[g] * Wm[[g]]
    Wmi.g[[g]] <- 1 / fg[g] * Wmi[[g]]

    # Gamma
    if (version == "V3") {
      Gamma.g <- Wmi[[g]] %*% Jm[[g]] %*% Wmi[[g]]
      Gamma.f[[g]] <- 1 / fg[g] * Gamma.g
    }
  }
  # create 'big' matrices
  Wc.all <- lav_matrix_bdiag(Wc.g)
  Wm.all <- lav_matrix_bdiag(Wm.g)
  Wmi.all <- lav_matrix_bdiag(Wmi.g)
  Delta.all <- do.call("rbind", Delta)

  # compute trace
  U <- Wm.all - Wm.all %*% Delta.all %*% E.inv %*% t(Delta.all) %*% Wm.all

  # V3 or V6?
  if (version == "V3") {
    Gamma.all <- lav_matrix_bdiag(Gamma.f)
    k.fimlc <-
      sum(diag(U %*% Wmi.all %*% Wc.all %*% Wmi.all %*% U %*% Gamma.all))
  } else {
    # V6
    k.fimlc <- sum(diag(Wc.all %*% Wmi.all %*% U %*% Wmi.all))
  }

  # convert to lavaan 'scaling.factor'
  c.hat3 <- k.fimlc / df3
  XX3.scaled <- XX3 / c.hat3

  # collect temp results
  out <- list(
    XX3 = XX3, df3 = df3,
    c.hat3 = c.hat3, XX3.scaled = XX3.scaled,
    XX3.null = as.numeric(NA), df3.null = as.numeric(NA),
    c.hat3.null = as.numeric(NA)
  )

  # baseline model
  fitB <- try(lav_object_independence(lavobject), silent = TRUE)
  if (inherits(fitB, "try-error")) {
    return(out)
  }
  # 'refit' using 'tilde' (=EM/saturated) sample statistics
  fitB.tilde <- try(lavaan(
    model = parTable(fitB),
    sample.cov = COV.tilde,
    sample.mean = MEAN.tilde,
    sample.nobs = sample.nobs,
    sample.cov.rescale = FALSE,
    information = "observed",
    optim.method = "none",
    se = "none",
    test = "standard",
    baseline = FALSE,
    check.post = FALSE
  ), silent = TRUE)
  if (inherits(fitB.tilde, "try-error")) {
    return(out)
  }

  XX3.null <- fitB.tilde@test[[1]]$stat
  df3.null <- fitB.tilde@test[[1]]$df

  fitB@Options$h1.information <- c("unstructured", "unstructured")
  fitB@Options$observed.information <- c("h1", "h1")
  E.invB <- lavTech(fitB, "inverted.information")
  DeltaB <- lavTech(fitB, "Delta")
  DeltaB.all <- do.call("rbind", DeltaB)

  # trace baseline model
  UB <- Wm.all - Wm.all %*% DeltaB.all %*% E.invB %*% t(DeltaB.all) %*% Wm.all
  # V3 or V6?
  if (version == "V3") {
    kb.fimlc <-
      sum(diag(UB %*% Wmi.all %*% Wc.all %*% Wmi.all %*% UB %*% Gamma.all))
  } else {
    # V6
    kb.fimlc <- sum(diag(Wc.all %*% Wmi.all %*% UB %*% Wmi.all))
  }

  # convert to lavaan 'scaling.factor'
  c.hat3.null <- kb.fimlc / df3.null

  # return values
  list(
    XX3 = XX3, df3 = df3, c.hat3 = c.hat3, XX3.scaled = XX3.scaled,
    XX3.null = XX3.null, df3.null = df3.null, c.hat3.null = c.hat3.null
  )
}
