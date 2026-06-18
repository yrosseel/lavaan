# initial version: YR 03/05/2017

# major change: YR 14/06/2022 for 0.6-12
# - summary() is now silent if not printed
# - here, we only collect the necessary ingredients, and store them in a
#   a list
# - the result is a S3 class lavaan.summary
# - the actual printing is done by lav_summary_print (see lav_print.R)

# YR 26 July 2022: add fm.args= argument to change the way (some) fit measures
#                  are computed
# YR 24 Sept 2022: add efa= argument
# YR 19 Nov  2023: add remove.unused= argument
# TDJ 28 March 2024: deprecate std.nox= argument ("std.nox" can be %in%
#                                                          standardized=)

# Is at least one defined (:=) parameter a NONLINEAR function of the free
# parameters? The (first-order) delta method is exact for linear definitions
# (sums/differences), so the se.def note in lav_summary_print() is only
# relevant when a nonlinear definition is present. Detected numerically: a
# function is linear iff its Jacobian does not depend on the parameter values.
lav_object_def_nonlinear <- function(object) {
  if (sum(object@ParTable$op == ":=") == 0L) {
    return(FALSE)
  }
  lavmodel <- object@Model
  if (!is.function(lavmodel@def.function)) {
    return(FALSE)
  }
  x <- try(lav_model_get_parameters(lavmodel, type = "free"), silent = TRUE)
  if (inherits(x, "try-error") || length(x) < 1L) {
    return(FALSE)
  }
  jac <- function(p) {
    out <- try(lav_func_jacobian_complex(lavmodel@def.function, p),
               silent = TRUE)
    if (inherits(out, "try-error")) {
      out <- try(lav_func_jacobian_simple(lavmodel@def.function, p),
                 silent = TRUE)
    }
    out
  }
  j1 <- jac(x)
  j2 <- jac(x + 1e-3 * (abs(x) + 1)) # deterministic perturbation
  if (inherits(j1, "try-error") || inherits(j2, "try-error")) {
    return(FALSE)
  }
  ok <- is.finite(j1) & is.finite(j2)
  isTRUE(any(abs(j1[ok] - j2[ok]) > 1e-8))
}

# create summary of a lavaan object
lav_object_summary <- function(object, header = TRUE,
                               fit_measures = FALSE,
                               baseline_model = NULL,
                               h1_model = NULL,
                               fm_args =
                                 list(
                                   standard.test = "default",
                                   scaled.test = "default",
                                   rmsea.ci.level = 0.90,
                                   rmsea.h0.closefit = 0.05,
                                   rmsea.h0.notclosefit = 0.08
                                 ),
                               estimates = TRUE,
                               ci = FALSE,
                               fmi = FALSE,
                               standardized = FALSE,
                               std = standardized,
                               remove_system_eq = TRUE,
                               remove_eq = TRUE,
                               remove_ineq = TRUE,
                               remove_def = FALSE,
                               remove_nonfree = FALSE,
                               remove_step1 = TRUE,
                               remove_unused = TRUE,
                               plabel = FALSE,
                               cov_std = TRUE,
                               rsquare = FALSE,
                               efa = FALSE,
                               efa_args =
                                 list(
                                   lambda = TRUE,
                                   theta = TRUE,
                                   psi = TRUE,
                                   eigenvalues = TRUE,
                                   sumsq.table = TRUE,
                                   lambda.structure = FALSE,
                                   fs.determinacy = FALSE,
                                   se = FALSE,
                                   zstat = FALSE,
                                   pvalue = FALSE
                                 ),
                               modindices = FALSE,
                               srmr_close_h0 = NULL) {

  # check object
  object <- lav_object_check_version(object)

  # default fm.args
  default_fm_args <- list(
    standard.test = "default",
    scaled.test = "default",
    rmsea.ci.level = 0.90,
    rmsea.close.h0 = 0.05,
    rmsea.notclose.h0 = 0.08,
    robust = TRUE,
    cat.check.pd = TRUE
  )
  if (is.logical(fit_measures)) {
    if (fit_measures) {
      fit_measures <- "default"
    } else {
      fit_measures <- "none"
    }
  }
  if (!missing(fm_args)) {
    lav_deprecated_args("fit.measures", "fm.args")
    fm_args <- modifyList(default_fm_args, fm_args)
  } else {
    fm_args <- default_fm_args
  }
  if (is.list(fit_measures)) {
    if (is.null(names(fit_measures)) ||
        is.null(fit_measures$fit.measures)) {
      lav_msg_stop(gettextf(
        "If %s is a list, it must contain a named element %s.",
        "fit.measures"
      ))
    }
    temp <- fit_measures$fit.measures
    fit_measures$fit.measures <- NULL
    fm_args <- modifyList(default_fm_args, fit_measures)
    fit_measures <- temp
  }


  # return a list with the main ingredients
  res <- list()

  # this is to avoid partial matching of 'std' with std.nox
  if (is.logical(std) && is.logical(standardized)) {
    standardized <- std || standardized
  } else {
    # At least 1 is not logical. Retain only valid standardization options.
    standardized <- intersect(union(tolower(std), tolower(standardized)),
                              c("std.lv", "std.all", "std.nox"))
  }

  # create the 'short' summary
  if (header) {
    # 1. collect header information
    version <- object@version

    res$header <- list(
      lavaan.version = version,
      sam.approach = !is.null(object@internal$sam.method),
      optim.method = object@Options$optim.method,
      optim.iterations = object@optim$iterations,
      optim.converged = object@optim$converged
    )

    # sam or sem?
    if (!is.null(object@internal$sam.method)) {
      # SAM version

      # 2. sam header
      res$sam.header <-
        list(
          sam.method = object@internal$sam.method,
          sam.local.options = object@internal$sam.local.options,
          sam.mm.list = object@internal$sam.mm.list,
          sam.mm.estimator = object@internal$sam.mm.estimator,
          sam.struc.estimator = object@internal$sam.struc.estimator
        )

      # 3. no EFA (for now)?

      # 4. summarize lavdata
      res$data <- lav_data_summary_short(object@Data)

      # 5a. sam local test statistics
      res$sam <- list(
        sam.method = object@internal$sam.method,
        sam.mm.table = object@internal$sam.mm.table,
        sam.mm.rel = object@internal$sam.mm.rel,
        sam.struc.fit = object@internal$sam.struc.fit,
        ngroups = object@Data@ngroups,
        group.label = object@Data@group.label,
        nlevels = object@Data@nlevels,
        level.label = object@Data@level.label,
        block.label = object@Data@block.label
      )

      # 5b. global test statistics (for global only)
      if (object@internal$sam.method == "global") {
        test <- object@test
        # attach the meta info that lav_test_print() needs (the SAM branch
        # does not go through the SEM path that normally sets this)
        if (is.null(attr(test, "info"))) {
          lavdata <- object@Data
          lavoptions <- object@Options
          attr(test, "info") <-
            list(
              ngroups = lavdata@ngroups,
              group.label = lavdata@group.label,
              information = lavoptions$information,
              h1.information = lavoptions$h1.information,
              observed.information = lavoptions$observed.information
            )
        }
        res$test <- test
      }
    } else {
      # SEM version

      # 2. summarize optim info (including estimator)
      nrow_ceq_jac <- nrow(object@Model@ceq.JAC)
      #if (object@Model@ceq.simple.only) {
        # not needed, as nrow.ceq.jac is already zero
      #  nrow.ceq.jac <- 0L
      #}
      cin_simple_only <- FALSE
      ceq_simple_only <- FALSE
      cin_simple_only <- object@Model@cin.simple.only
      ceq_simple_only <- object@Model@ceq.simple.only
      nrow_cin_jac <- nrow(object@Model@cin.JAC)
      if (cin_simple_only) {
        nrow_cin_jac <- 0L
      }
      if (ceq_simple_only && cin_simple_only) {
        nrow_con_jac <- 0L
        con_jac_rank <- 0L
      } else {
        nrow_con_jac <- nrow(object@Model@con.jac)
        con_jac_rank <- qr(object@Model@con.jac)$rank
      }

      res$optim <- list(
        estimator = object@Options$estimator,
        estimator.args = object@Options$estimator.args,
        optim.method = object@Options$optim.method,
        npar = object@Model@nx.free,
        eq.constraints = object@Model@eq.constraints,
        nrow.ceq.jac = nrow_ceq_jac,
        nrow.cin.jac = nrow_cin_jac,
        nrow.con.jac = nrow_con_jac,
        con.jac.rank = con_jac_rank
      )


      # 3. if EFA/ESEM, summarize rotation info
      if (object@Model@nefa > 0L) {
        res$rotation <-
          list(
            rotation = object@Options$rotation,
            rotation.args = object@Options$rotation.args
          )
      }

      # 4. summarize lavdata
      res$data <- lav_data_summary_short(object@Data)

      # 5. test statistics
      test <- object@test
      # TDJ: check for user-supplied h1 model
      if (!is.null(object@external$h1.model)) {
        stopifnot(inherits(object@external$h1.model, "lavaan"))
        ## update @test slot
        test <- lav_update_test_custom_h1(
          lav_obj_h0 = object, lav_obj_h1 = object@external$h1.model)@test
      }
      # double check if we have attr(TEST, "info") (perhaps old object?)
      if (is.null(attr(test, "info"))) {
        lavdata <- object@Data
        lavoptions <- object@Options
        attr(test, "info") <-
          list(
            ngroups = lavdata@ngroups,
            group.label = lavdata@group.label,
            information = lavoptions$information,
            h1.information = lavoptions$h1.information,
            observed.information = lavoptions$observed.information
          )
      }
      res$test <- test
    } # regular sem
  } # header

  # efa-related info
  if (efa) {
    res$efa <- lav_efa_summary(object, efa_args = efa_args)
  } # efa

  # only if requested, add the additional fit measures
  if (fit_measures != "none") {
    # some early warnings (to avoid a hard stop)
    if (object@Data@data.type == "none") {
      lav_msg_warn(gettext(
        "fit measures not available if there is no data"))
    } else if (length(object@Options$test) == 1L &&
      object@Options$test == "none") {
      lav_msg_warn(gettext(
        "fit measures not available if test = \"none\""))
    } else if (object@optim$npar > 0L && !object@optim$converged) {
      lav_msg_warn(gettext(
        "fit measures not available if model did not converge"))
    } else {
      fit <- lav_fit(object,
       fit_measures = c(list(fit.measures = fit_measures), fm_args),
       baseline_model = baseline_model,
       h1_model = h1_model)
      res$fit <- fit
    }
  }

  # main ingredient: the parameter table
  if (estimates) {
    pe <- lavParameterEstimates(object,
      ci = ci, standardized = standardized,
      rsquare = rsquare, fmi = fmi,
      cov.std = cov_std,
      remove.eq = remove_eq, remove.system.eq = remove_system_eq,
      remove.ineq = remove_ineq, remove.def = remove_def,
      remove.nonfree = remove_nonfree,
      remove.step1 = remove_step1,
      remove.unused = remove_unused,
      plabel = plabel,
      output = "text",
      header = TRUE
    )
    res$pe <- as.data.frame(pe)

    # flag: nonlinear defined (:=) parameters whose SE/CI rely on the
    # (first-order) delta method -> suggest se.def="mc" or bootstrap
    # (printed by lav_summary_print(); silent otherwise). Stored as an
    # attribute of the parameter table (res$pe) instead of a separate
    # list element, to avoid disturbing code that relies on partial
    # matching of names(summary(object)) (e.g. $p -> $pe).
    attr(res$pe, "delta.note") <-
      identical(object@Options$se.def, "default") &&
      !identical(object@Options$se, "bootstrap") &&
      !identical(object@Options$se, "none") &&
      lav_object_def_nonlinear(object)
  }

  # modification indices?
  if (is.logical(modindices)) {
    if (modindices) modindices <- list(standardized = TRUE, cov.std = cov_std)
  }
  if (is.list(modindices)) {
    mi <- do.call("modificationIndices", c(list(object = object), modindices))
    res$mi <- mi
  }

  # create lavaan.summary S3 class
  class(res) <- c("lavaan.summary", "list")

  res
}
