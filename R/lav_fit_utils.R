# utility functions needed to compute various (robust) fit measures:
#
# - lav_fit_catml_dwls (for 'robust' RMSEA/CFI if data is categorical)
# - lav_fit_fiml_corrected (correct RMSEA/CFI if data is incomplete)

# compute scaling-factor (c.hat3) for fit.dwls, using fit.catml ingredients
# see:
#     Savalei, V. (2021) Improving Fit Indices In SEM with categorical data.
#     Multivariate Behavioral Research, 56(3), 390-407.
#

# YR Dec 2022: first version
# YR Jan 2023: catml_dwls should check if the input 'correlation' matrix
#              is positive-definite (or not)

lav_fit_catml_dwls <- function(lavobject, check_pd = TRUE) {
  # empty list
  empty_list <- list(
    XX3 = as.numeric(NA), df3 = as.numeric(NA),
    c.hat3 = as.numeric(NA), XX3.scaled = as.numeric(NA),
    XX3.null = as.numeric(NA), df3.null = as.numeric(NA),
    c.hat3.null = as.numeric(NA)
  )

  # limitations
  if (!lavobject@Model@categorical ||
    lavobject@Options$conditional.x ||
    length(unlist(lavobject@pta$vnames$ov.num)) > 0L) {
    return(empty_list)
  } else {
    lavdata <- lavobject@Data
    lavsamplestats <- lavobject@SampleStats
  }

  # check if input matrix (or matrices) are all positive definite
  # (perhaps later, we can rely on 'smoothing', but not for now
  # pd_flag <- TRUE
  if (check_pd) {
    for (g in seq_len(lavdata@ngroups)) {
      cor_1 <- lavsamplestats@cov[[g]]
      ev <- eigen(cor_1, symmetric = TRUE, only.values = TRUE)$values
      if (any(ev < .Machine$double.eps^(1 / 2))) {
        # non-pd!
        # pd_flag <- FALSE
        # should we give a warning here? (not for now)
        # warning("lavaan WARNING: robust RMSEA/CFI could not be computed
        #   because the input correlation matrix is not positive-definite")
        # what should we do? return NA (for now)
        return(empty_list)
      }
    }
  }

  # 'refit' using estimator = "catML"
  fit_catml <- try(lav_object_catml(lavobject), silent = TRUE)
  if (inherits(fit_catml, "try-error")) {
    return(empty_list)
  }

  xx3 <- fit_catml@test[[1]]$stat
  df3 <- fit_catml@test[[1]]$df


  # compute 'k'
  v <- lavTech(fit_catml, "wls.v") # NT-ML weight matrix

  w_dwls <- lavTech(lavobject, "wls.v") # DWLS weight matrix
  gamma <- lavTech(lavobject, "gamma") # acov of polychorics
  delta <- lavTech(lavobject, "delta")
  e_inv <- lavTech(lavobject, "inverted.information")

  fg <- unlist(lavsamplestats@nobs) / lavsamplestats@ntotal

  # Fixme: as we only need the trace, perhaps we could do this
  # group-specific? (see lav_test_satorra_bentler_trace_original)
  v_g <- v
  w_dwls_g <- w_dwls
  gamma_f <- gamma
  delta_g <- delta
  for (g in seq_len(lavdata@ngroups)) {
    ntotal <- nrow(gamma[[g]])
    nvar <- lavobject@Model@nvar[[g]]
    pstar <- nvar * (nvar - 1) / 2
    rm_idx <- seq_len(ntotal - pstar)

    # reduce
    delta_g[[g]] <- delta[[g]][-rm_idx, , drop = FALSE]
    # reduce and weight
    w_dwls_g[[g]] <- fg[g] * w_dwls[[g]][-rm_idx, -rm_idx]
    v_g[[g]] <- fg[g] * v[[g]] # should already have the right dims
    gamma_f[[g]] <- 1 / fg[g] * gamma[[g]][-rm_idx, -rm_idx]
  }
  # create 'big' matrices
  w_dwls_all <- lav_matrix_bdiag(w_dwls_g)
  v_all <- lav_matrix_bdiag(v_g)
  gamma_all <- lav_matrix_bdiag(gamma_f)
  delta_all <- do.call("rbind", delta_g)

  # compute trace
  wi_u_all <- diag(nrow(w_dwls_all)) -
              delta_all %*% e_inv %*% t(delta_all) %*% w_dwls_all
  ks <- sum(diag(t(wi_u_all) %*% v_all %*% wi_u_all %*% gamma_all))

  # convert to lavaan 'scaling.factor'
  c_hat3 <- ks / df3
  xx3_scaled <- xx3 / c_hat3

  # baseline model
  xx3_null <- fit_catml@baseline$test[[1]]$stat
  if (is.null(xx3_null)) {
    xx3_null <- as.numeric(NA)
    df3_null <- as.numeric(NA)
    kbs <- as.numeric(NA)
    c_hat3_null <- as.numeric(NA)
  } else {
    df3_null <- fit_catml@baseline$test[[1]]$df
    kbs <- sum(diag(gamma_all))
    c_hat3_null <- kbs / df3_null
  }

  # return values
  list(
    XX3 = xx3, df3 = df3, c.hat3 = c_hat3, XX3.scaled = xx3_scaled,
    XX3.null = xx3_null, df3.null = df3_null, c.hat3.null = c_hat3_null
  )
}


# compute ingredients to compute FIML-Corrected RMSEA/CFI
# see:
#     Zhang X, Savalei V. (2022). New computations for RMSEA and CFI
#     following FIML and TS estimation with missing data. Psychological Methods.

lav_fit_fiml_corrected <- function(lavobject, baseline_model,
                                   version = "V3") {
  version <- toupper(version)
  if (!version %in% c("V3", "V6")) {
    lav_msg_stop(gettext("only FIML-C(V3) and FIML-C(V6) are available."))
  }

  # empty list
  empty_list <- list(
    XX3 = as.numeric(NA), df3 = as.numeric(NA),
    c.hat3 = as.numeric(NA), XX3.scaled = as.numeric(NA),
    XX3.null = as.numeric(NA), df3.null = as.numeric(NA),
    c.hat3.null = as.numeric(NA)
  )

  # limitations
  if (lavobject@Options$conditional.x ||
    lavobject@Data@nlevels > 1L ||
    is.null(lavobject@h1$implied$cov[[1]])) {
    return(empty_list)
  } else {
    lavdata <- lavobject@Data
    lavsamplestats <- lavobject@SampleStats

    h1 <- lavTech(lavobject, "h1", add.labels = TRUE)
    cov_tilde <- lapply(h1, "[[", "cov")
    mean_tilde <- lapply(h1, "[[", "mean")
    sample_nobs <- unlist(lavsamplestats@nobs)
  }

  # 'refit' using 'tilde' (=EM/saturated) sample statistics
  fit_tilde <- try(lavaan(
    model = parTable(lavobject),
    sample.cov = cov_tilde,
    sample.mean = mean_tilde,
    sample.nobs = sample_nobs,
    sample.cov.rescale = FALSE,
    information = "observed",
    optim.method = "none",
    se = "none",
    test = "standard",
    baseline = FALSE,
    check.post = FALSE
  ), silent = TRUE)
  if (inherits(fit_tilde, "try-error")) {
    return(empty_list)
  }

  xx3 <- fit_tilde@test[[1]]$stat
  df3 <- fit_tilde@test[[1]]$df

  # compute 'k'

  # V3/V6: always use h1.information = "unstructured"!!
  lavobject@Options$h1.information <- c("unstructured", "unstructured")
  lavobject@Options$observed.information <- c("h1", "h1")
  fit_tilde@Options$h1.information <- c("unstructured", "unstructured")
  fit_tilde@Options$observed.information <- c("h1", "h1")

  wm <- wm_g <- lav_model_h1_information_observed(lavobject)
  wc <- wc_g <- lav_model_h1_information_observed(fit_tilde)

  if (version == "V3") {
    jm <- jm_g <- lav_model_h1_information_firstorder(lavobject)
    gamma_f <- vector("list", length = lavdata@ngroups)
  }
  delta <- lavTech(lavobject, "delta")
  e_inv <- lavTech(lavobject, "inverted.information")
  wmi <- wmi_g <- try(lapply(wm, lav_matrix_symmetric_inverse),
    silent = TRUE
  )
  if (inherits(wmi, "try-error")) {
    return(empty_list)
  }

  fg <- unlist(lavsamplestats@nobs) / lavsamplestats@ntotal
  # Fixme: as we only need the trace, perhaps we could do this
  # group-specific? (see lav_test_satorra_bentler_trace_original)
  for (g in seq_len(lavdata@ngroups)) {
    # group weight
    wc_g[[g]] <- fg[g] * wc[[g]]
    wm_g[[g]] <- fg[g] * wm[[g]]
    wmi_g[[g]] <- 1 / fg[g] * wmi[[g]]

    # gamma
    if (version == "V3") {
      jm_g[[g]] <- fg[g] * jm[[g]]
      gamma_g <- wmi[[g]] %*% jm[[g]] %*% wmi[[g]]
      gamma_f[[g]] <- 1 / fg[g] * gamma_g
    }
  }
  # create 'big' matrices
  wc_all <- lav_matrix_bdiag(wc_g)
  wm_all <- lav_matrix_bdiag(wm_g)
  wmi_all <- lav_matrix_bdiag(wmi_g)
  delta_all <- do.call("rbind", delta)

  e_comp <- t(delta_all) %*% wc_all %*% delta_all
               # VS: or grab from fit.tilde, with observed.info="h1"


  # compute trace
  if (version == "V3") {
    gamma_all <- lav_matrix_bdiag(gamma_f)
    # VS: Simplification of k.fimlc to minimize matrix multiplication
    #                                                 of big matrices
    jm_all <- lav_matrix_bdiag(jm_g)

    # VS: tr11 is also used for baseline
    # VS: tr(AB) = sum(A*t(B)) is more efficient

    tr11 <- sum(wc_all * gamma_all)
    tr12 <- sum((t(delta_all) %*% jm_all %*% wmi_all %*%
                 wc_all %*% delta_all) * e_inv)
    tr22 <- sum((t(delta_all) %*% jm_all %*% delta_all %*% e_inv)
             * t(t(delta_all) %*% wc_all %*% delta_all %*% e_inv))

    k_fimlc <- tr11 - 2 * tr12 + tr22
  } else {
    # V6
    tr1 <- sum(wc_all * wmi_all)
    k_fimlc <- tr1 - sum(e_comp * e_inv)
  }

  # convert to lavaan 'scaling.factor'
  c_hat3 <- k_fimlc / df3
  xx3_scaled <- xx3 / c_hat3

  # collect temp results
  out <- list(
    XX3 = xx3, df3 = df3,
    c.hat3 = c_hat3, XX3.scaled = xx3_scaled,
    XX3.null = as.numeric(NA), df3.null = as.numeric(NA),
    c.hat3.null = as.numeric(NA)
  )

  # baseline model
  if (!is.null(baseline_model)) {
    fit_b <- baseline_model
  } else {
    fit_b <- try(lav_object_independence(lavobject), silent = TRUE)
  }

  if (inherits(fit_b, "try-error")) {
    return(out)
  }

  # 'refit' using 'tilde' (=EM/saturated) sample statistics
  fit_b_tilde <- try(lavaan(
    model = parTable(fit_b),
    sample.cov = cov_tilde,
    sample.mean = mean_tilde,
    sample.nobs = sample_nobs,
    sample.cov.rescale = FALSE,
    information = "observed",
    optim.method = "none",
    se = "none",
    test = "standard",
    baseline = FALSE,
    check.post = FALSE
  ), silent = TRUE)
  if (inherits(fit_b_tilde, "try-error")) {
    return(out)
  }

  xx3_null <- fit_b_tilde@test[[1]]$stat
  df3_null <- fit_b_tilde@test[[1]]$df

  fit_b@Options$h1.information <- c("unstructured", "unstructured")
  fit_b@Options$observed.information <- c("h1", "h1")
  fit_b_tilde@Options$h1.information <- c("unstructured", "unstructured")
  fit_b_tilde@Options$observed.information <- c("h1", "h1")

  e_inv_b <- lavTech(fit_b, "inverted.information")
  delta_b <- lavTech(fit_b, "Delta")
  delta_b_all <- do.call("rbind", delta_b)

  e_comp_b <- t(delta_b_all) %*% wc_all %*% delta_b_all #or grab from fitB.tilde

  # V3 or V6?
  if (version == "V3") {
    tr12b <-
      sum((t(delta_b_all) %*% jm_all %*% wmi_all %*% wc_all %*% delta_b_all) *
        e_inv_b)
    tr22b <-
      sum((t(delta_b_all) %*% jm_all %*% delta_b_all %*% e_inv_b) *
        t(t(delta_b_all) %*% wc_all %*% delta_b_all %*% e_inv_b))
    kb_fimlc <- tr11 - 2 * tr12b + tr22b
  } else {
    # V6
    kb_fimlc <- tr1 - sum(e_comp_b * e_inv_b)
  }

  # convert to lavaan 'scaling.factor'
  c_hat3_null <- kb_fimlc / df3_null

  # return values
  list(
    XX3 = xx3, df3 = df3, c.hat3 = c_hat3, XX3.scaled = xx3_scaled,
    XX3.null = xx3_null, df3.null = df3_null, c.hat3.null = c_hat3_null
  )
}
