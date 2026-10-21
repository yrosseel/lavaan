# importing results from an external estimation backend into a lavaan
# object (step 3 of the backend-integration plan for data-defined, ":~",
# parameters; see also lav_partable_dv.R)
#
# the object passed in is the (unfitted) object lavaan prepared for the
# backend (do.fit = FALSE): it holds the parameter table, the data, the
# sample statistics and the model matrices. The backend returns (at least)
# the estimates of the free parameters, keyed by plabel (or label), and a
# loglikelihood. We rebuild the object through lavaan's own pipeline
# (optim.method = "none", as sam() does), and then fill the slots the
# pipeline can not compute itself: @vcov, @loglik, @optim, @external.
#
# YR/Claude 21 Oct 2026: initial version

# public entry point
lavImportResults <- function(object,                         # nolint
                             est,
                             vcov = NULL,
                             se = NULL,
                             loglik = NULL,
                             converged = TRUE,
                             iterations = NA_integer_,
                             backend = "external",
                             status = NULL,
                             nobs = NULL,
                             loglik_casewise = NULL,
                             scores = NULL,
                             backend_object = NULL,
                             ...) {
  dotdotdot <- list(...)
  lav_adapt_func(environment(), dotdotdot, NULL)
  results <- list(
    est = est, vcov = vcov, se = se, loglik = loglik,
    converged = converged, iterations = iterations,
    backend = backend, status = status, nobs = nobs,
    loglik_casewise = loglik_casewise, scores = scores,
    object = backend_object
  )
  lav_object_backend_results(object = object, results = results)
}

# map a (named) estimate vector onto lavaan's free-parameter vector;
# names are plabels (preferred) or labels; an unnamed vector must have
# length nx.free and is taken to be in free-parameter order
lav_backend_map_x <- function(partable, nx_free, values, what = "est") {
  if (!is.numeric(values)) {
    lav_msg_stop(gettextf("`%s' must be a numeric vector.", what))
  }
  if (any(!is.finite(values))) {
    lav_msg_stop(gettextf("`%s' contains non-finite values.", what))
  }
  if (is.null(names(values))) {
    if (length(values) != nx_free) {
      lav_msg_stop(gettextf(
        "`%1$s' has length %2$d, but the model has %3$d free parameters;
         supply a named vector (names = plabels or labels), or a vector of
         length %3$d in free-parameter order.", what, length(values),
        nx_free))
    }
    return(as.numeric(values))
  }
  free_idx <- which(partable$free > 0L)
  plabel <- partable$plabel[free_idx]
  label <- partable$label[free_idx]
  # third key: the names coef() produces (lhs-op-rhs, with group suffixes)
  coef_names <- lav_pt_labels(partable, type = "user")[free_idx]
  x <- rep(NA_real_, nx_free)
  seen <- character(0L)
  for (nm in names(values)) {
    row_idx <- free_idx[plabel == nm]
    if (length(row_idx) == 0L) {
      row_idx <- free_idx[nzchar(label) & label == nm]
    }
    if (length(row_idx) == 0L) {
      row_idx <- free_idx[coef_names == nm]
    }
    if (length(row_idx) == 0L) {
      lav_msg_stop(gettextf(
        "`%1$s' contains the name %2$s, which is not the plabel, the label,
         or the coef() name of a free parameter.", what,
        dQuote(nm, q = FALSE)))
    }
    for (ri in row_idx) {
      k <- partable$free[ri]
      if (!is.na(x[k]) && abs(x[k] - values[[nm]]) >
            1e-8 * max(1, abs(x[k]))) {
        lav_msg_warn(gettextf(
          "`%1$s' supplies different values for the same free parameter
           (%2$s: %3$s and %4$s); using the first one.",
          what, dQuote(nm, q = FALSE), x[k], values[[nm]]))
      } else {
        x[k] <- values[[nm]]
      }
    }
    seen <- c(seen, nm)
  }
  if (anyNA(x)) {
    missing_k <- which(is.na(x))
    missing_names <- vapply(missing_k, function(k) {
      ri <- free_idx[partable$free[free_idx] == k][1L]
      if (nzchar(partable$label[ri])) partable$label[ri] else
        partable$plabel[ri]
    }, "")
    lav_msg_stop(gettextf(
      "`%1$s' does not provide a value for free parameter(s): %2$s.",
      what, lav_msg_view(missing_names, "none")))
  }
  x
}

# the number of free parameters, taking equality constraints into account
# (as in lav_model_loglik)
lav_backend_npar <- function(lavmodel) {
  npar <- lavmodel@nx.free
  if (lavmodel@ceq.simple.only) {
    # nothing to do
  } else if (!lavmodel@cin.simple.only && nrow(lavmodel@con.jac) > 0L) {
    ceq_idx <- attr(lavmodel@con.jac, "ceq.idx")
    if (length(ceq_idx) > 0L) {
      npar <- npar - qr(lavmodel@con.jac[ceq_idx, , drop = FALSE])$rank
    }
  }
  npar
}

# internal importer: 'results' is the list a backend hook returns
lav_object_backend_results <- function(object, results) {
  stopifnot(inherits(object, "lavaan"))
  object <- lav_object_check_version(object)
  lavmodel <- object@Model
  lavpartable <- object@ParTable
  nx_free <- lavmodel@nx.free
  dv_flag <- lav_model_dv_flag(lavmodel)
  backend_name <- if (is.null(results$backend)) "external" else
    as.character(results$backend)[1L]

  if (nx_free == 0L) {
    lav_msg_stop(gettext("the model has no free parameters."))
  }

  # 1. estimates
  x <- lav_backend_map_x(lavpartable, nx_free, results$est, what = "est")
  if (any(!is.finite(x))) {
    lav_msg_stop(gettext("`est' contains non-finite values."))
  }

  # 2. constraints: the backend is responsible for honouring them
  con_tol <- 1e-5
  if (!is.null(body(lavmodel@ceq.function))) {
    ceq <- try(lavmodel@ceq.function(x), silent = TRUE)
    if (!inherits(ceq, "try-error") && length(ceq) > 0L &&
        any(abs(ceq) > con_tol)) {
      lav_msg_warn(gettextf(
        "the imported estimates do not satisfy the equality constraints of
         the model (max |violation| = %s); the backend (%s) may have
         ignored them.", format(max(abs(ceq)), digits = 3), backend_name))
    }
  }
  if (!is.null(body(lavmodel@cin.function))) {
    cin <- try(lavmodel@cin.function(x), silent = TRUE)
    if (!inherits(cin, "try-error") && length(cin) > 0L &&
        any(cin < -con_tol)) {
      lav_msg_warn(gettextf(
        "the imported estimates do not satisfy the inequality constraints
         of the model (max violation = %s); the backend (%s) may have
         ignored them.", format(max(-cin), digits = 3), backend_name))
    }
  }

  # 3. vcov / se
  vcov_mat <- NULL
  se_free <- NULL
  if (!is.null(results$vcov)) {
    vcov_mat <- results$vcov
    if (!is.matrix(vcov_mat) || nrow(vcov_mat) != nx_free ||
        ncol(vcov_mat) != nx_free) {
      lav_msg_stop(gettextf(
        "`vcov' must be a %1$d x %1$d matrix (one row/column per free
         parameter, in free-parameter order).", nx_free))
    }
    if (!is.null(rownames(vcov_mat)) || !is.null(colnames(vcov_mat))) {
      # reorder by name, if the names are plabels/labels
      nms <- if (!is.null(rownames(vcov_mat))) rownames(vcov_mat) else
        colnames(vcov_mat)
      perm <- lav_backend_map_x(lavpartable, nx_free,
                                setNames(seq_len(nx_free), nms),
                                what = "vcov dimnames")
      vcov_mat <- vcov_mat[perm, perm, drop = FALSE]
    }
    dimnames(vcov_mat) <- NULL
    storage.mode(vcov_mat) <- "double"
    asym <- max(abs(vcov_mat - t(vcov_mat)))
    if (asym > 1e-6 * max(1, max(abs(vcov_mat)))) {
      lav_msg_stop(gettext("`vcov' is not symmetric."))
    }
    vcov_mat <- (vcov_mat + t(vcov_mat)) / 2
    ev <- eigen(vcov_mat, symmetric = TRUE, only.values = TRUE)$values
    if (any(ev < -1e-8 * max(1, max(abs(ev))))) {
      lav_msg_warn(gettextf(
        "`vcov' is not positive semi-definite (smallest eigenvalue %s).",
        format(min(ev), digits = 3)))
    }
    se_free <- sqrt(pmax(diag(vcov_mat), 0))
  } else if (!is.null(results$se)) {
    se_free <- lav_backend_map_x(lavpartable, nx_free, results$se,
                                 what = "se")
  }

  # 4. number of observations: the backend must have used lavaan's cases
  ntotal <- object@SampleStats@ntotal
  if (!is.null(results$nobs) && !isTRUE(all.equal(as.numeric(results$nobs),
                                                   as.numeric(ntotal)))) {
    lav_msg_stop(gettextf(
      "the backend used %1$s observations, but lavaan retained %2$s (see
       lavInspect(object, \"case.idx\")); the backend must use the same
       cases.", results$nobs, ntotal))
  }

  # 5. loglikelihood
  loglik <- results$loglik
  if (!is.null(loglik)) {
    if (!is.numeric(loglik) || length(loglik) != 1L || !is.finite(loglik)) {
      lav_msg_stop(gettext("`loglik' must be a single finite number."))
    }
  } else if (dv_flag) {
    lav_msg_warn(gettext(
      "no loglikelihood supplied: logLik(), AIC() and BIC() will not be
       available."))
  }

  # 6. rebuild the object through the pipeline, with the estimates as
  #    starting values and no optimization
  lavoptions <- object@Options
  lavoptions$optim.method <- "none"
  lavoptions$optim.parscale <- "none"
  lavoptions$optim.force.converged <- TRUE
  lavoptions$start <- "default"
  lavoptions$check.gradient <- FALSE
  lavoptions$check.start <- FALSE
  lavoptions$se <- "external"
  lavoptions$do.fit <- TRUE
  lavoptions$rotation <- "none"
  lavoptions$backend <- "lavaan" # the dry run must not call the backend
  if (dv_flag) {
    lavoptions$test <- "none"
    lavoptions$baseline <- FALSE
    lavoptions$fit.by.level <- FALSE
  }
  free_idx <- which(lavpartable$free > 0L)
  lavpartable$ustart[free_idx] <- x[lavpartable$free[free_idx]]
  lavpartable$se <- rep(NA_real_, length(lavpartable$lhs))
  lavpartable$est <- NULL
  lavpartable$start <- NULL

  current_warn <- lav_warn()
  lav_warn(FALSE) # the dry run would warn about se = "external"
  on.exit(lav_warn(current_warn), add = TRUE)
  fit <- lavaan(lavpartable,
    slot_options = lavoptions,
    slot_data = object@Data,
    slot_sample_stats = object@SampleStats,
    verbose = FALSE
  )
  lav_warn(current_warn)

  # 7. fill the slots the pipeline can not compute
  # 7a. parameter table: standard errors (incl. := rows, delta method)
  if (!is.null(vcov_mat)) {
    fit@ParTable$se <- lav_model_vcov_se(
      lavmodel = fit@Model, lavpartable = fit@ParTable,
      vcov = vcov_mat, lavoptions = fit@Options
    )
  } else if (!is.null(se_free)) {
    se_user <- rep(NA_real_, length(fit@ParTable$lhs))
    free_idx <- which(fit@ParTable$free > 0L)
    se_user[free_idx] <- se_free[fit@ParTable$free[free_idx]]
    # fixed parameters: 0; defined (:=) parameters: unknown without vcov
    se_user[fit@ParTable$free == 0L & fit@ParTable$op != ":="] <- 0
    fit@ParTable$se <- se_user
  }
  # host rows of ":~" definitions: no scalar estimate exists
  if (dv_flag && !is.null(fit@ParTable$dv)) {
    host_idx <- which(nzchar(fit@ParTable$dv))
    fit@ParTable$est[host_idx] <- NA_real_
    fit@ParTable$se[host_idx] <- NA_real_
  }

  # 7b. vcov slot
  fit@vcov <- list(
    se = "external",
    information = "external",
    vcov = vcov_mat
  )

  # 7c. loglikelihood slot
  if (!is.null(loglik)) {
    npar <- lav_backend_npar(fit@Model)
    fit@loglik <- list(
      loglik = loglik,
      loglik.group = if (object@Data@ngroups == 1L) loglik else NULL,
      npar = npar,
      ntotal = ntotal,
      AIC = lav_fit_aic(logl = loglik, npar = npar),
      BIC = lav_fit_bic(logl = loglik, npar = npar, n = ntotal),
      BIC2 = lav_fit_sabic(logl = loglik, npar = npar, n = ntotal),
      estimator = fit@Options$estimator,
      conditional.x = fit@Options$conditional.x,
      fixed.x = fit@Options$fixed.x,
      backend = backend_name
    )
  } else if (dv_flag) {
    fit@loglik <- list()
  }

  # 7d. optimization slot
  fit@optim$iterations <- if (is.null(results$iterations) ||
                              is.na(results$iterations)) NA_integer_ else
    as.integer(results$iterations)
  fit@optim$converged <- isTRUE(results$converged)
  if (dv_flag) {
    fit@optim$fx <- if (is.null(loglik)) NA_real_ else -2 * loglik
    fit@optim$fx.group <- NULL
  }
  fit@optim$backend <- list(
    name = backend_name,
    status = results$status,
    converged = isTRUE(results$converged)
  )

  # 7e. everything else the backend supplied
  fit@external$backend <- list(
    name = backend_name,
    status = results$status,
    loglik_casewise = results$loglik_casewise,
    scores = results$scores,
    object = results$object
  )

  # 7f. the call of the object lavaan prepared for the backend
  fit@call <- object@call

  if (!isTRUE(results$converged)) {
    lav_msg_warn(gettextf(
      "the backend (%s) reports that the estimation did not converge;
       the estimates may be unreliable.", backend_name))
  }

  fit
}

# header rows for summary(): replace the information block by the
# backend information (cfr. lav_noniter_se_rows)
lav_backend_se_rows <- function(object = NULL) {
  backend <- object@optim$backend
  if (is.null(backend)) {
    return(NULL)
  }
  if (identical(object@Options$se, "bootstrap")) {
    n_req <- object@Options$bootstrap$R
    if (is.null(n_req)) {
      n_req <- NA
    }
    n_ok <- if (is.null(object@boot$coef)) NA else {
      nrow(object@boot$coef) - length(attr(object@boot$coef, "error.idx"))
    }
    return(c("Standard errors" = "Bootstrap",
             "Number of requested bootstrap draws" = as.character(n_req),
             "Number of successful bootstrap draws" = as.character(n_ok),
             "Estimation backend" = backend$name))
  }
  out <- c("Standard errors" = "External",
           "Estimation backend" = backend$name)
  if (is.null(object@vcov$vcov)) {
    out["Standard errors"] <- "None"
  }
  out
}

# se = "bootstrap" for a model estimated by a backend: the generic
# nonparametric bootstrap (lav_bootstrap_internal) refits every resampled
# data set through lavaan(), which hands each replicate to the backend
# again (the stored options carry the backend name and do.fit = TRUE)
lav_backend_bootstrap <- function(fit, lavoptions = NULL) {
  if (is.null(lavoptions)) {
    lavoptions <- fit@Options
  }
  r <- 1000L
  if (!is.null(lavoptions$bootstrap) && !is.null(lavoptions$bootstrap$R)) {
    r <- as.integer(lavoptions$bootstrap$R)
  }
  if ("bollen.stine" %in% lavoptions$test || "yuan" %in% lavoptions$test) {
    lav_msg_stop(gettext(
      "model-based (Bollen-Stine, Yuan) bootstrap is not available for
       models estimated by an external backend; only the ordinary
       nonparametric bootstrap is."))
  }
  lavoptions_b <- fit@Options
  lavoptions_b$do.fit <- TRUE
  lavoptions_b$se <- "none"
  lavoptions_b$parallel <- lavoptions$parallel
  lavoptions_b$ncpus <- lavoptions$ncpus
  lavoptions_b$cl <- lavoptions$cl
  lavoptions_b$iseed <- lavoptions$iseed
  lavoptions_b$bootstrap <- lavoptions$bootstrap

  if (lav_verbose()) {
    cat("computing VCOV for      se = bootstrap (backend ",
        fit@optim$backend$name, ") ...", sep = "")
  }
  coef_1 <- lav_bootstrap_internal(
    object = NULL,
    lavmodel = fit@Model,
    lavsamplestats = fit@SampleStats,
    lavpartable = fit@ParTable,
    lavoptions = lavoptions_b,
    lavdata = fit@Data,
    r = r,
    check_post = lavoptions$check.post,
    type = "ordinary",
    fun = "coef"
  )
  if (lav_verbose()) {
    cat(" done.\n")
  }
  coef_orig <- coef_1
  error_idx <- attr(coef_1, "error.idx")
  if (length(error_idx) > 0L) {
    lav_msg_warn(gettextf(
      "%s bootstrap runs failed or did not converge.", length(error_idx)))
    coef_1 <- coef_1[-error_idx, , drop = FALSE]
  }
  notok <- length(attr(coef_1, "nonadmissible"))
  if (notok > 0L) {
    lav_msg_warn(gettextf(
      "%s bootstrap runs resulted in nonadmissible solutions.", notok))
  }
  nboot <- nrow(coef_1)
  if (nboot < 2L) {
    lav_msg_warn(gettext(
      "fewer than 2 successful bootstrap draws; no standard errors."))
    return(fit)
  }
  vcov_mat <- cov(coef_1) * (nboot - 1) / nboot
  dimnames(vcov_mat) <- NULL

  fit@boot <- list(coef = coef_orig)
  fit@vcov <- list(se = "bootstrap", information = "external",
                   vcov = vcov_mat)
  fit@Options$se <- "bootstrap"
  fit@ParTable$se <- lav_model_vcov_se(
    lavmodel = fit@Model, lavpartable = fit@ParTable,
    vcov = vcov_mat, boot = coef_1, lavoptions = fit@Options
  )
  if (lav_model_dv_flag(fit@Model) && !is.null(fit@ParTable$dv)) {
    host_idx <- which(nzchar(fit@ParTable$dv))
    fit@ParTable$se[host_idx] <- NA_real_
  }
  fit
}

# ---------------------------------------------------------------------
# the backend registry (step 4 of the plan): named estimation backends
#
# a backend is a descriptor: a name, a function
#
#     fun(object, args)
#
# that receives the lavaan object prepared for it (do.fit = FALSE) and the
# backend.args list, and returns the results list lavImportResults() takes
# (est, vcov/se, loglik, converged, iterations, status, nobs,
# loglik_casewise, scores, object), plus metadata: the package that
# provides it (loaded on demand, also on parallel workers), a version,
# the model classes it supports ("all", "dv" = data-defined parameters,
# "standard" = ordinary models), and free-form capabilities.
#
# the `backend' option is a string resolved through this registry (the
# form that survives the lavaan 2 transition: no R closure involved); a
# bare function is accepted as development sugar and registered under an
# anonymous name for the session.
lav_backend_registry <- new.env(parent = emptyenv())

lavBackendRegister <- function(name,                           # nolint
                               fun,
                               package = NULL,
                               version = NULL,
                               supports = "all",
                               capabilities = list(),
                               description = "",
                               overwrite = TRUE) {
  if (!is.character(name) || length(name) != 1L || !nzchar(name)) {
    lav_msg_stop(gettext("`name' must be a single (non-empty) string."))
  }
  if (tolower(name) == "lavaan") {
    lav_msg_stop(gettext("the name \"lavaan\" is reserved."))
  }
  if (!is.function(fun)) {
    lav_msg_stop(gettext("`fun' must be a function(object, args)."))
  }
  supports <- lav_keyword_canonical(supports)
  bad <- supports[!supports %in% c("all", "dv", "standard")]
  if (length(bad) > 0L) {
    lav_msg_stop(gettextf(
      "`supports' must be a subset of %s.",
      lav_msg_view(c("all", "dv", "standard"), "none")))
  }
  key <- tolower(name)
  if (!overwrite && exists(key, envir = lav_backend_registry,
                           inherits = FALSE)) {
    lav_msg_stop(gettextf("a backend named %s is already registered.",
                          dQuote(name, q = FALSE)))
  }
  descriptor <- list(
    name = name,
    fun = fun,
    package = if (is.null(package)) NA_character_ else as.character(package)[1L],
    version = if (is.null(version)) NA_character_ else as.character(version)[1L],
    supports = supports,
    capabilities = capabilities,
    description = as.character(description)[1L]
  )
  assign(key, descriptor, envir = lav_backend_registry)
  invisible(descriptor)
}

lavBackends <- function() {                                      # nolint
  keys <- sort(ls(lav_backend_registry))
  out <- data.frame(
    name = character(0L), package = character(0L),
    version = character(0L), supports = character(0L),
    description = character(0L), stringsAsFactors = FALSE
  )
  for (key in keys) {
    d <- get(key, envir = lav_backend_registry)
    out[nrow(out) + 1L, ] <- list(
      d$name, d$package, d$version, paste(d$supports, collapse = ","),
      d$description
    )
  }
  out
}

# look up a backend by name; if unknown, try to load a package of that
# name (its .onLoad() is expected to register the backend)
lav_backend_get <- function(name) {
  key <- tolower(name)
  if (!exists(key, envir = lav_backend_registry, inherits = FALSE)) {
    if (requireNamespace(name, quietly = TRUE) &&
        exists(key, envir = lav_backend_registry, inherits = FALSE)) {
      # registered by the package's .onLoad()
    } else {
      registered <- gettext("none")
      if (length(ls(lav_backend_registry)) > 0L) {
        registered <- lav_msg_view(lavBackends()$name, "none")
      }
      lav_msg_stop(gettextf(
        "unknown estimation backend %1$s; registered backends: %2$s (see
         lavBackendRegister()).", dQuote(name, q = FALSE), registered))
    }
  }
  get(key, envir = lav_backend_registry)
}

# resolve the `backend' option (a string, or a function) to a descriptor
lav_backend_resolve <- function(backend) {
  if (is.function(backend)) {
    # already registered under an anonymous name?
    for (key in ls(lav_backend_registry)) {
      d <- get(key, envir = lav_backend_registry)
      if (startsWith(key, "function#") && identical(d$fun, backend)) {
        return(d)
      }
    }
    k <- sum(startsWith(ls(lav_backend_registry), "function#")) + 1L
    return(lavBackendRegister(
      name = paste0("function#", k), fun = backend,
      description = gettext("anonymous backend function (this session)")))
  }
  if (!is.character(backend) || length(backend) != 1L) {
    lav_msg_stop(gettext(
      "the `backend' option must be a single string (the name of a
       registered backend) or a function(object, args)."))
  }
  lav_backend_get(backend)
}

# run a backend on the prepared object, and import its results
lav_backend_run <- function(object, descriptor, args = list()) {
  dv_flag <- lav_model_dv_flag(object@Model)
  needed <- if (dv_flag) "dv" else "standard"
  if (!any(descriptor$supports %in% c("all", needed))) {
    lav_msg_stop(gettextf(
      "backend %1$s does not support this model class (%2$s); it supports:
       %3$s.", dQuote(descriptor$name, q = FALSE), needed,
      lav_msg_view(descriptor$supports, "none")))
  }
  if (!is.na(descriptor$package) &&
      !requireNamespace(descriptor$package, quietly = TRUE)) {
    lav_msg_stop(gettextf(
      "backend %1$s needs package %2$s, which is not available.",
      dQuote(descriptor$name, q = FALSE), dQuote(descriptor$package, q = FALSE)))
  }
  if (lav_verbose()) {
    cat(sprintf("estimation by backend %s ...", descriptor$name))
  }
  results <- descriptor$fun(object, args)
  if (lav_verbose()) {
    cat(" done.\n")
  }
  if (!is.list(results) || is.null(results$est)) {
    lav_msg_stop(gettextf(
      "backend %s did not return a list with (at least) an `est' element.",
      dQuote(descriptor$name, q = FALSE)))
  }
  if (is.null(results$backend)) {
    results$backend <- descriptor$name
  }
  fit <- lav_object_backend_results(object = object, results = results)
  fit@Options$backend <- descriptor$name
  fit
}
