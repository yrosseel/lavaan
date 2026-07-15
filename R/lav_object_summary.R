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
# YR 19 Nov  2023: add remove_unused= argument
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
                               residuals = FALSE,
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
                                   sumsq_table = TRUE,
                                   lambda_structure = FALSE,
                                   fs_determinacy = FALSE,
                                   se = FALSE,
                                   zstat = FALSE,
                                   pvalue = FALSE
                                 ),
                               modindices = FALSE,
                               srmr_close_h0 = NULL,
                              ...) {
   dotdotdot <- list(...)
   lav_adapt_func(environment(), dotdotdot, NULL)
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
    cat.nonpd = "na"
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
    fit_measures <- lav_fit_measures_list_alias(fit_measures)
    if (is.null(names(fit_measures)) ||
        is.null(fit_measures$fit.measures)) {
      lav_msg_stop(gettextf(
        "If %s is a list, it must contain a named element %s.",
        "fit.measures", "fit.measures (or fit_measures)"
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
    # At least 1 is not logical.
    tmp <- union(
      if (is.logical(std)) character(0L) else as.character(std),
      if (is.logical(standardized)) character(0L) else as.character(standardized)
    )
    std_types <- c("std.lv", "std.all", "std.nox")
    if (length(tmp) > 0L && !any(tolower(tmp) %in% std_types)) {
      # a vector of observed variable names: standardize only those (a
      # generalization of std.nox); pass through unchanged
      standardized <- tmp
    } else {
      # Retain only valid standardization options.
      standardized <- intersect(tolower(tmp), std_types)
    }
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
      } else if (!identical(fit_measures, "none") &&
                 !is.null(lav_sam_struc_object(object))) {
        # 5c. local family + fit.measures: print the (structural) chi-square
        # test statistic(s) as usual, instead of the one-line
        # 'Summary Information Structural part:' (issue #517)
        fit_pa <- lav_sam_struc_object(object)
        test <- fit_pa@test
        if (is.null(attr(test, "info"))) {
          attr(test, "info") <-
            list(
              ngroups = fit_pa@Data@ngroups,
              group.label = fit_pa@Data@group.label,
              information = fit_pa@Options$information,
              h1.information = fit_pa@Options$h1.information,
              observed.information = fit_pa@Options$observed.information
            )
        }
        # the header is printed before the test statistics, making clear
        # that the model test (and the fit measures below it) concern the
        # structural part only
        attr(test, "header") <- gettext(
          "Fit measures of the structural part (given the measurement model):")
        res$test <- test
        res$sam$sam.struc.fit <- NULL
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
  # (fit_measures may be a vector of measure names; only the single
  #  "none" turns this section off)
  if (!identical(fit_measures, "none")) {
    # some early warnings (to avoid a hard stop)
    if (object@Data@data.type == "none") {
      lav_msg_warn(gettext(
        "fit measures not available if there is no data"))
    } else if (lav_sam_local_flag(object) &&
      is.null(lav_sam_struc_object(object))) {
      lav_msg_warn(gettext(
        "fit measures not available: this sam object was created by an
         older version of lavaan; please rerun sam()"))
    } else if (length(object@Options$test) == 1L &&
      object@Options$test == "none" &&
      length(object@Model@rv.ov) == 0L &&
      length(object@Model@rv.lv) == 0L) {
      # (random-slope models do provide the loglikelihood-based
      #  measures, despite test = "none")
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
      # local sam objects: the 'structural part' header is already printed
      # once, before the test statistics (see res$test above)
      if (!is.null(attr(res$test, "header"))) {
        attr(fit, "header") <- NULL
      }
      res$fit <- fit

      # multilevel: add level-specific fit measures (based on the stored
      # partially saturated fits; see fit.by.level option)
      if (object@Data@nlevels > 1L) {
        res$fit.by.level <- lav_fit_by_level_fm(object)
      }
    }
  }

  # only if requested, add the residuals (and their summary). 'residuals' may be
  # TRUE/FALSE, or a (named) list of arguments passed on to lavResiduals().
  if ((is.list(residuals) || isTRUE(residuals)) &&
      object@optim$npar > 0L && !object@optim$converged) {
    lav_msg_warn(gettext("residuals not available if model did not converge"))
    residuals <- FALSE
  }
  # sam objects: warn-and-skip where lavResiduals() would stop (see the
  # sam branch in lavResiduals())
  if ((is.list(residuals) || isTRUE(residuals)) &&
      !is.null(object@internal$sam.method)) {
    if (lav_sam_local_flag(object) &&
        is.null(lav_sam_struc_object(object))) {
      lav_msg_warn(gettext(
        "residuals not available: this sam object was created by an older
         version of lavaan; please rerun sam()"))
      residuals <- FALSE
    } else if (!lav_sam_local_flag(object)) {
      lav_msg_warn(gettext(
        "residual summary statistics are not available if sam.method =
        \"global\"; use residuals() to inspect the raw residuals"))
      residuals <- FALSE
    }
  }
  if (is.list(residuals) || isTRUE(residuals)) {
    res_args <- if (is.list(residuals)) residuals else list()
    # defaults (overridable via the list form), matching lavResiduals()
    res_defaults <- list(
      type = "cor.bentler", se = FALSE, zstat = TRUE, n.largest = 5L,
      combine = FALSE, usrmr.ci.level = 0.90, usrmr.close.h0 = 0.05
    )
    res_args <- modifyList(res_defaults, res_args)
    # the summary needs both the element-wise residuals and the summary tables
    res_args <- c(
      res_args,
      list(
        object = object, summary = TRUE, elementwise = TRUE, output = "list"
      )
    )
    res_list <- do.call(lavResiduals, res_args)
    res$residuals <- res_list

    # for local sam objects, lavResiduals() delegated to the structural part
    res_pt <- object@ParTable
    if (lav_sam_local_flag(object)) {
      res_pt <- lav_sam_struc_object(object)@ParTable
    }

    # pre-build the text representation used by lav_summary_print()
    res_text <- lav_residuals_list_to_text(
      res_list, type = res_args$type,
      nblocks = lav_pt_nblocks(res_pt), drop_single = TRUE,
      n_largest = res_args$n.largest, show_se = isTRUE(res_args$se),
      show_z = isTRUE(res_args$zstat), ci_level = res_args$usrmr.ci.level,
      combine = isTRUE(res_args$combine), do_summary = TRUE, do_largest = TRUE
    )
    attr(res$residuals, "text") <- res_text

    # if the fit measures are also shown AND the residual summary is the SRMR,
    # drop the SRMR from the fit measures to avoid printing it twice
    if (!is.null(res$fit) && !is.null(res_text$summary_tables)) {
      r_metric <- rownames(res_text$summary_tables[[1]])[1]
      drop_srmr <- c("srmr", "srmr_nomean", "srmr_within", "srmr_between")
      if (identical(r_metric, "srmr") &&
          any(names(res$fit) %in% drop_srmr)) {
        fit_attr <- attributes(res$fit)
        res$fit <- res$fit[setdiff(names(res$fit), drop_srmr)]
        # restore the class (and any other attributes) lost by subsetting
        for (a in setdiff(names(fit_attr), "names")) {
          attr(res$fit, a) <- fit_attr[[a]]
        }
      }
    }
  }

  # main ingredient: the parameter table
  if (estimates) {
    pe <- lavParameterEstimates(object,
      ci = ci, standardized = standardized,
      rsquare = rsquare, fmi = fmi,
      cov_std = cov_std,
      remove_eq = remove_eq, remove_system_eq = remove_system_eq,
      remove_ineq = remove_ineq, remove_def = remove_def,
      remove_nonfree = remove_nonfree,
      remove_step1 = remove_step1,
      remove_unused = remove_unused,
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
    if (modindices) modindices <- list(standardized = TRUE, cov_std = cov_std)
  }
  if (is.list(modindices)) {
    mi <- do.call("modificationIndices", c(list(object = object), modindices))
    res$mi <- mi
  }

  # create lavaan.summary S3 class
  class(res) <- c("lavaan.summary", "list")

  res
}
