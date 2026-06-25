# `methods' for fitted lavaan objects
#
# standard (S4) methods:
# - show()
# - summary()
# - coef()
# - fitted.values() + fitted()
# - vcov()
# - logLik()
# - nobs()
# - update()
# - anova()

# lavaan-specific methods:
#
# - lavParameterEstimates()
# - standardizedSolution()
# - parameterTable()
# - varTable()


setMethod(
  "show", "lavaan",
  function(object) {
    # efa?
    efa_flag <- object@Options$model.type == "efa"

    # show only basic information
    res <- lav_object_summary(object,
      fit_measures = FALSE,
      estimates = FALSE,
      modindices = FALSE,
      efa = efa_flag
    )
    if (efa_flag) {
      # print (standardized) loadings only
      class(res) <- c("lavaan.efa", "list")
      print(res)
    } else {
      # print lavaan header
      print(res)
    }
    invisible(res)
  }
)

setMethod(
  "summary", "lavaan",
  function(object, header = TRUE,                                # nolint start
           fit.measures = FALSE,
           residuals = FALSE,
           estimates = TRUE,
           ci = FALSE,
           fmi = FALSE,
           standardized = FALSE,
           std = standardized,
           std.nox = FALSE, # TODO: remove deprecated argument in early 2025
           remove.system.eq = TRUE,
           remove.eq = TRUE,
           remove.ineq = TRUE,
           remove.def = FALSE,
           remove.nonfree = FALSE,
           remove.step1 = TRUE,
           remove.unused = TRUE,
           plabel = FALSE,
           cov.std = TRUE,
           rsquare = FALSE,
           baseline.model = NULL,
           h1.model = NULL,
           fm.args = list(
             standard.test = "default",
             scaled.test = "default",
             rmsea.ci.level = 0.90,
             rmsea.h0.closefit = 0.05,
             rmsea.h0.notclosefit = 0.08,
             robust = TRUE,
             cat.check.pd = TRUE
           ),
           modindices = FALSE,
           srmr.close.h0 = NULL,
           nd = 3L, cutoff = 0.3, dot.cutoff = 0.1, ...) {      # nolint end
    dotdotdot <- list(...)
    if (length(dotdotdot) > 0L) {
      for (j in seq_along(dotdotdot)) {
        lav_msg_warn(gettextf(
          "Unknown argument %s for %s", sQuote(names(dotdotdot)[j]),
          sQuote("summary"))
        )
      }
    }
    # rename modified argument(s)
    fit_measures <- fit.measures

    # efa?
    efa_flag <- object@Options$model.type == "efa"

    if (is.logical(fit_measures)) {
      fit_measures <- if (fit_measures) "default" else "none"
    }
    if (!is.list(fit_measures))
        fit_measures <- list(fit.measures = fit_measures)
    if (is.logical(fit_measures$fit.measures)) {
      fit_measures$fit.measures <- if (fit_measures$fit.measures) "default" else "none"
    }
    if (!missing(fm.args)) {
      lav_deprecated_args("fit.measures", "fm.args")
      fit_measures <- c(fit_measures, fm.args)
    }
    res <- lav_object_summary(
      object = object, header = header,
      fit_measures = fit_measures, residuals = residuals,
      estimates = estimates,
      baseline_model = baseline.model,
      h1_model = h1.model,
      ci = ci, fmi = fmi, std = std, standardized = standardized,
      remove_system_eq = remove.system.eq,
      remove_eq = remove.eq, remove_ineq = remove.ineq,
      remove_def = remove.def, remove_nonfree = remove.nonfree,
      remove_step1 = remove.step1, remove_unused = remove.unused,
      plabel = plabel, cov_std = cov.std,
      rsquare = rsquare, efa = efa_flag,
      modindices = modindices,
      srmr_close_h0 = srmr.close.h0
    )
    # res has class c("lavaan.summary", "list")

    # what about nd? only used if we actually print; save as attribute
    attr(res, "nd") <- nd

    # if efa, add cutoff and dot.cutoff, and change class
    if (efa_flag) {
      # class(res) <- c("lavaan.summary.efa", "list")
      attr(res, "cutoff") <- cutoff
      attr(res, "dot.cutoff") <- dot.cutoff
    }

    res
  }
)


setMethod(
  "coef", "lavaan",
  function(object, type = "free", labels = TRUE, ...) {
    dotdotdot <- list(...)
    if (length(dotdotdot) > 0L) {
      for (j in seq_along(dotdotdot)) {
        lav_msg_warn(gettextf(
          "Unknown argument %s for %s", sQuote(names(dotdotdot)[j]),
            sQuote("coef"))
        )
      }
    }
    # check object
    object <- lav_object_check_version(object)

    lav_inspect_coef(
      object = object, type = type,
      add_labels = labels, add_class = TRUE
    )
  }
)

# Warn when computing a *standardized* defined (:=) parameter is ambiguous.
#
# If a (user) label is attached to two or more model parameters -- ie they are
# constrained to be equal -- and that same label is referenced inside a defined
# parameter (:=), then the standardized value of the defined parameter is not
# well-defined: equal *unstandardized* coefficients need not have equal
# *standardized* coefficients (the latter also depend on the standard deviations
# of the variables involved, which typically differ across groups/equations).
# lavaan must still report a single number, and silently uses the value from the
# first occurrence of the label. This check issues an informative warning, as
# suggested by Keith Markus on the lavaan discussion list:
#   https://groups.google.com/g/lavaan/c/AWGORQnmKBg
#
# 'est.std' (optional) is the vector of standardized estimates aligned with the
# rows of 'partable'. When supplied, a label is only flagged if the standardized
# values of its (equal) parameters actually differ -- so the warning never fires
# when standardization happens to be invariant.
lav_partable_check_std_def <- function(partable, est.std = NULL,
                                       tol = 1e-08) {
  # we need labels, and at least one := definition
  if (is.null(partable$label)) {
    return(invisible(NULL))
  }
  def_idx <- which(partable$op == ":=")
  if (length(def_idx) == 0L) {
    return(invisible(NULL))
  }

  # labels attached to 2+ *model* parameters (these are the equality
  # constraints); exclude ==/</>/:= rows (the latter carries the def name)
  model_idx <- which(!partable$op %in% c("==", "<", ">", ":=") &
                       nzchar(partable$label))
  if (length(model_idx) == 0L) {
    return(invisible(NULL))
  }
  lab <- partable$label[model_idx]
  shared <- unique(lab[duplicated(lab)])
  if (length(shared) == 0L) {
    return(invisible(NULL))
  }

  def_rhs <- partable$rhs[def_idx]
  def_lhs <- partable$lhs[def_idx]
  flagged_lab <- character(0L)
  flagged_def <- character(0L)
  for (lb in shared) {
    # if standardized values are available, only flag when they truly differ
    if (!is.null(est.std)) {
      vals <- est.std[model_idx[lab == lb]]
      vals <- vals[is.finite(vals)]
      if (length(vals) < 2L ||
          diff(range(vals)) <= tol * max(1, max(abs(vals)))) {
        next
      }
    }
    # is the label referenced (as a whole token) inside any := definition?
    pat <- paste0("\\b", gsub("([][{}().|^$*+?\\\\])", "\\\\\\1", lb), "\\b")
    hit <- grepl(pat, def_rhs)
    if (any(hit)) {
      flagged_lab <- c(flagged_lab, lb)
      flagged_def <- unique(c(flagged_def, def_lhs[hit]))
    }
  }
  if (length(flagged_lab) == 0L) {
    return(invisible(NULL))
  }

  lav_msg_warn(gettextf(
    "the standardized value of the defined parameter(s) %s is ambiguous: it
     uses the label(s) %s, which are constrained equal across groups (or
     equations). Equal unstandardized coefficients can have different
     standardized coefficients, so lavaan reports the standardized defined
     parameter based on the first occurrence only. Interpret with caution.",
    paste(sQuote(flagged_def), collapse = ", "),
    paste(sQuote(flagged_lab), collapse = ", ")))

  invisible(list(def = flagged_def, label = flagged_lab))
}

standardizedSolution <-                                      # nolint start
  standardizedsolution <- function(object,
                                   type = "std.all",
                                   se = TRUE,
                                   zstat = TRUE,
                                   pvalue = TRUE,
                                   ci = TRUE,
                                   level = 0.95,
                                   boot.ci.type = "perc",
                                   cov.std = TRUE,
                                   remove.eq = TRUE,
                                   remove.ineq = TRUE,
                                   remove.def = FALSE,
                                   remove.aux = TRUE,
                                   partable = NULL,
                                   GLIST = NULL,
                                   est = NULL,
                                   output = "data.frame") { # nolint end

    # check object
    object <- lav_object_check_version(object)

    # check type: either one of the std.* keywords, or a vector of (observed)
    # variable names, in which case only the parameters involving these
    # variables are standardized (a generalization of "std.nox")
    ov_std_user <- NULL
    std_types <- c("std.all", "std.lv", "std.nox")
    if (length(type) > 1L || !any(tolower(type) %in% std_types)) {
      ov_all <- lav_object_vnames(object, "ov")
      bad <- type[!(type %in% ov_all)]
      if (length(bad) > 0L) {
        lav_msg_warn(gettextf(
          "type= argument contains unknown observed variable(s): %s",
          paste(bad, collapse = ", ")))
      }
      ov_std_user <- type[type %in% ov_all]
      if (length(ov_std_user) == 0L) {
        lav_msg_stop(gettext(
          "type= must be one of std.lv, std.all, std.nox, or a vector of
           observed variable names."))
      }
      type <- "std.user"
    } else {
      type <- tolower(type)
      stopifnot(type %in% std_types)
    }

    # check output= argument
    output <- tolower(output)
    if (output %in% c("data.frame", "table")) {
      output <- "data.frame"
    } else if (output %in% c("text", "pretty")) {
      output <- "text"
    } else {
      lav_msg_stop(gettextf(
        "output must be %s or %s", sQuote("data.frame"), sQuote("text"))
      )
    }

    # no zstat + pvalue if estimator is Bayes
    if (object@Options$estimator == "Bayes") {
      zstat <- pvalue <- FALSE
    }

    # no se if class is not lavaan
    # using class() -- can't use inherits(), as this includes blavaan
    if (class(object)[1L] != "lavaan") {
      if (missing(se) || !se) {
        se <- FALSE
        zstat <- FALSE
        pvalue <- FALSE
      }
    }

    if (is.null(partable)) {
      tmp_partable <- lavInspect(object, "list")
    } else {
      tmp_partable <- partable
    }
    tmp_list <- tmp_partable[, c("lhs", "op", "rhs", "exo")]
    if (!is.null(tmp_partable$group)) {
      tmp_list$group <- tmp_partable$group
    }
    if (!is.null(tmp_partable$block)) {
      tmp_list$block <- tmp_partable$block
    }
    if (sum(nchar(tmp_partable$label)) != 0L) {
      tmp_list$label <- tmp_partable$label
    }

    # add std and std.all columns
    if (type == "std.lv") {
      tmp_list$est.std <- lav_standardize_lv(object,
        est = est, glist = GLIST,
        partable = partable, cov_std = cov.std
      )
    } else if (type == "std.all") {
      tmp_list$est.std <- lav_standardize_all(object,
        est = est, glist = GLIST,
        partable = partable, cov_std = cov.std
      )
    } else if (type == "std.nox") {
      tmp_list$est.std <- lav_standardize_all_nox(object,
        est = est, glist = GLIST,
        partable = partable, cov_std = cov.std
      )
    } else if (type == "std.user") {
      tmp_list$est.std <- lav_standardize_all(object,
        est = est, glist = GLIST,
        partable = partable, cov_std = cov.std, ov_std = ov_std_user
      )
    }

    # warn if a defined (:=) parameter standardizes a label that is constrained
    # equal across groups/equations (its standardized value is then ambiguous)
    lav_partable_check_std_def(tmp_partable, est.std = tmp_list$est.std)

    # bootstrap? If the model was fitted with se = "bootstrap" (and bootstrap
    # draws are available), we provide bootstrap standard errors and bootstrap
    # confidence intervals for the standardized parameters (just like
    # parameterEstimates() does for the unstandardized parameters). Every row
    # of the standardized solution is a function of the free parameters, so we
    # simply re-standardize each bootstrap draw.
    boot_std <- NULL
    boot_error_idx <- integer(0L)
    tmp_boot_free <- NULL
    tmp_fun_x <- NULL
    if (object@Options$se == "bootstrap" && se &&
        !is.null(object@boot$coef) && NROW(object@boot$coef) > 0L) {
      stopifnot(boot.ci.type %in%
        c("norm", "basic", "perc", "bca.simple", "bca"))
      # standardization as a function of the free parameter vector
      if (type == "std.lv") {
        tmp_fun_x <- lav_standardize_lv_x
      } else if (type == "std.all") {
        tmp_fun_x <- lav_standardize_all_x
      } else if (type == "std.nox") {
        tmp_fun_x <- lav_standardize_all_nox_x
      } else if (type == "std.user") {
        tmp_fun_x <- function(x, lavobject) {
          lav_standardize_all_x(x, lavobject = lavobject, ov_std = ov_std_user)
        }
      }
      tmp_boot_free <- lav_inspect_boot(object)
      boot_error_idx <- attr(tmp_boot_free, "error.idx")
      if (length(boot_error_idx) > 0L) {
        tmp_boot_free <- tmp_boot_free[-boot_error_idx, , drop = FALSE]
      }
      boot_std <- try(apply(tmp_boot_free, 1L, tmp_fun_x, lavobject = object),
                      silent = TRUE)
      if (inherits(boot_std, "try-error") || is.null(boot_std)) {
        boot_std <- NULL
      } else {
        if (is.matrix(boot_std)) {
          boot_std <- t(boot_std)
        } else {
          boot_std <- as.matrix(boot_std)
        }
        boot_std[!is.finite(boot_std)] <- as.numeric(NA)
      }
    }

    if (object@Options$se != "none" && se && !is.null(boot_std)) {
      # bootstrap standard errors (MLE-normalised covariance, as elsewhere)
      nboot <- nrow(boot_std)
      tmp <- apply(boot_std, 2L, sd, na.rm = TRUE) * sqrt((nboot - 1) / nboot)
      # catch near-zero SEs
      zero_idx <- which(tmp < .Machine$double.eps^(1 / 4))
      if (length(zero_idx) > 0L) {
        tmp[zero_idx] <- 0.0
      }
      tmp_list$se <- tmp
      if (zstat) {
        tmp_se <- ifelse(tmp_list$se == 0.0, NA, tmp_list$se)
        tmp_list$z <- tmp_list$est.std / tmp_se
      }
      if (zstat && pvalue) {
        tmp_list$pvalue <- 2 * (1 - pnorm(abs(tmp_list$z)))
      }
    } else if (object@Options$se != "none" && se) {
      # add 'se' for standardized parameters
      tmp_vcov <- try(lav_inspect_vcov(object,
        standardized = TRUE,
        type = type, ov_std = ov_std_user, free_only = FALSE,
        add_labels = FALSE,
        add_class = FALSE
      ))
      if (inherits(tmp_vcov, "try-error") || is.null(tmp_vcov)) {
        tmp_list$se <- rep(NA, length(tmp_list$lhs))
        if (zstat) {
          tmp_list$z <- rep(NA, length(tmp_list$lhs))
        }
        if (pvalue) {
          tmp_list$pvalue <- rep(NA, length(tmp_list$lhs))
        }
      } else {
        tmp <- diag(tmp_vcov)
        # catch negative values
        min_idx <- which(tmp < 0)
        if (length(min_idx) > 0L) {
          tmp[min_idx] <- as.numeric(NA)
        }
        # now, we can safely take the square root
        tmp <- sqrt(tmp)

        # catch near-zero SEs
        zero_idx <- which(tmp < .Machine$double.eps^(1 / 4)) # was 1/2 < 0.6
        # was 1/3 < 0.6-9
        if (length(zero_idx) > 0L) {
          tmp[zero_idx] <- 0.0
        }
        tmp_list$se <- tmp

        # add 'z' column
        if (zstat) {
          tmp_se <- ifelse(tmp_list$se == 0.0, NA, tmp_list$se)
          tmp_list$z <- tmp_list$est.std / tmp_se
        }
        if (zstat && pvalue) {
          tmp_list$pvalue <- 2 * (1 - pnorm(abs(tmp_list$z)))
        }
      }
    }

    # simple symmetric confidence interval
    if (se && object@Options$se != "none" && ci) {
      # next three lines based on confint.lm
      a <- (1 - level) / 2
      a <- c(a, 1 - a)
      fac <- qnorm(a)
      ci <- tmp_list$est.std + tmp_list$se %o% fac

      if (!is.null(boot_std)) {
        # bootstrap confidence intervals (shared with parameterEstimates() and
        # lavEffects() via lav_bootstrap_ci(); see lav_bootstrap.R). Every row
        # is a derived quantity, so the CI is computed uniformly for all rows.
        t0 <- tmp_list$est.std
        if (boot.ci.type == "bca") {
          design <- lav_bootstrap_bca_design(object, boot_error_idx)
          stopifnot(nrow(design) == nrow(boot_std))
          acc <- lav_bootstrap_acceleration(boot_std, design)
          ci <- lav_bootstrap_ci(boot_t = boot_std, t0 = t0,
            boot.ci.type = "bca", level = level, acc = acc)
        } else {
          # For "norm", the bias is evaluated at the mean of the bootstrapped
          # (free) parameters, exactly as in parameterEstimates().
          bias <- if (boot.ci.type == "norm") {
            as.numeric(tmp_fun_x(colMeans(tmp_boot_free, na.rm = TRUE),
                                 lavobject = object)) - t0
          } else {
            NULL
          }
          ci <- lav_bootstrap_ci(boot_t = boot_std, t0 = t0,
            boot.ci.type = boot.ci.type, level = level,
            se = tmp_list$se, bias = bias)
        }
      } else {
        # Monte Carlo: asymmetric percentile-based CI for standardized
        # defined parameters (Preacher & Selig 2012)
        mc_coef <- lav_inspect_mc(object)
        def_idx <- which(tmp_list$op == ":=")
        if (!is.null(mc_coef) && length(def_idx) > 0L) {
          if (type == "std.lv") {
            tmp_fun_mc <- lav_standardize_lv_x
          } else if (type == "std.all") {
            tmp_fun_mc <- lav_standardize_all_x
          } else if (type == "std.nox") {
            tmp_fun_mc <- lav_standardize_all_nox_x
          } else if (type == "std.user") {
            tmp_fun_mc <- function(x, lavobject) {
              lav_standardize_all_x(x, lavobject = lavobject,
                                    ov_std = ov_std_user)
            }
          }
          mc_std <- try(apply(mc_coef, 1L, tmp_fun_mc,
                              lavobject = object), silent = TRUE)
          if (!inherits(mc_std, "try-error")) {
            if (is.matrix(mc_std)) {
              mc_std <- t(mc_std)
            } else {
              mc_std <- as.matrix(mc_std)
            }
            alpha <- c((1 - level) / 2, 1 - (1 - level) / 2)
            mc_qq <- apply(mc_std[, def_idx, drop = FALSE], 2L,
                           quantile, probs = alpha, na.rm = TRUE, type = 7)
            ci[def_idx, ] <- t(mc_qq)
          }
        }
      }

      tmp_list$ci.lower <- ci[, 1]
      tmp_list$ci.upper <- ci[, 2]
    }


    # if single group, remove group column
    if (object@Data@ngroups == 1L) tmp_list$group <- NULL

    # remove == rows?
    if (remove.eq) {
      eq_idx <- which(tmp_list$op == "==")
      if (length(eq_idx) > 0L) {
        tmp_list <- tmp_list[-eq_idx, ]
      }
    }
    # remove <> rows?
    if (remove.ineq) {
      ineq_idx <- which(tmp_list$op %in% c("<", ">"))
      if (length(ineq_idx) > 0L) {
        tmp_list <- tmp_list[-ineq_idx, ]
      }
    }
    # remove := rows?
    if (remove.def) {
      def_idx <- which(tmp_list$op == ":=")
      if (length(def_idx) > 0L) {
        tmp_list <- tmp_list[-def_idx, ]
      }
    }
    # remove auxiliary (saturated-correlates) rows? (FIML aux= variables)
    if (remove.aux && length(object@Options$aux) > 0L) {
      aux_idx <- which(tmp_list$lhs %in% object@Options$aux |
        tmp_list$rhs %in% object@Options$aux)
      if (length(aux_idx) > 0L) {
        tmp_list <- tmp_list[-aux_idx, ]
      }
    }

    # remove attribute for data order
    attr(tmp_list, "ovda") <- NULL

    if (output == "text") {
      class(tmp_list) <- c(
        "lavaan.parameterEstimates", "lavaan.data.frame",
        "data.frame"
      )
      # tmp_list$exo is needed for printing, don't remove it
      attr(tmp_list, "group.label") <- object@Data@group.label
      attr(tmp_list, "level.label") <- object@Data@level.label
      # attr(tmp_list, "header") <- FALSE
    } else {
      tmp_list$exo <- NULL
      tmp_list$block <- NULL
      class(tmp_list) <- c("lavaan.data.frame", "data.frame")
    }

    tmp_list
  }

lavParameterEstimates <- function(object,                      # nolint start
                                 # select columns
                                 se = TRUE,
                                 zstat = TRUE,
                                 pvalue = TRUE,
                                 ci = TRUE,
                                 standardized = FALSE,
                                 fmi = FALSE,
                                 plabel = FALSE,
                                 # control
                                 level = 0.95,
                                 boot.ci.type = "perc",
                                 cov.std = TRUE,
                                 fmi.options = list(),
                                 # add rows
                                 rsquare = FALSE,
                                 # remove rows
                                 remove.system.eq = TRUE,
                                 remove.eq = TRUE,
                                 remove.ineq = TRUE,
                                 remove.def = FALSE,
                                 remove.nonfree = FALSE,
                                 remove.step1 = TRUE,
                                 remove.unused = FALSE,
                                 remove.aux = TRUE,
                                 # output
                                 add.attributes = FALSE,
                                 output = "data.frame",
                                 header = FALSE) {           # nolint end

    # lavaan.fsr?
    if (inherits(object, "lavaan.fsr")) {
      return(object$PE)
    }

    # check object
    object <- lav_object_check_version(object)

    # deprecated add.attributes (for psycho/blavaan)
    if (add.attributes) {
      output <- "text"
    }


    # no se if class is not lavaan
    # can't use inherits(), as this would return TRUE if object is from blavaan
    if (class(object)[1L] != "lavaan") {
      if (missing(se) || !se) {
        se <- FALSE
        zstat <- FALSE
        pvalue <- FALSE
      }
    }

    # check output= argument
    output <- tolower(output)
    if (output %in% c("data.frame", "table")) {
      output <- "data.frame"
      header <- FALSE
    } else if (output %in% c("text", "pretty")) {
      output <- "text"
    } else {
      lav_msg_stop(gettextf(
        "output must be %s or %s", sQuote("data.frame"), sQuote("text"))
      )
    }

    # check fmi
    if (fmi) {
      if (inherits(object, "lavaanList")) {
        lav_msg_warn(gettext(
          "fmi not available for object of class \"lavaanList\""))
        fmi <- FALSE
      }
      if (object@Options$se != "standard") {
        lav_msg_warn(gettext(
          "fmi only available if se = \"standard\""))
        fmi <- FALSE
      }
      if (object@Options$estimator != "ML") {
        lav_msg_warn(gettext(
          "fmi only available if estimator = \"ML\""))
        fmi <- FALSE
      }
      if (!object@SampleStats@missing.flag) {
        lav_msg_warn(gettext(
          "fmi only available if missing = \"(fi)ml\""))
        fmi <- FALSE
      }
      if (!object@optim$converged) {
        lav_msg_warn(gettext(
          "fmi not available; model did not converge"))
        fmi <- FALSE
      }
    }

    # no zstat + pvalue if estimator is Bayes
    if (object@Options$estimator == "Bayes") {
      zstat <- pvalue <- FALSE
    }

    tmp_partable <- as.data.frame(object@ParTable, stringsAsFactors = FALSE)
    tmp_list <- tmp_partable[, c("lhs", "op", "rhs", "free")]
    if (!is.null(tmp_partable$user)) {
      tmp_list$user <- tmp_partable$user
    }
    if (!is.null(tmp_partable$block)) {
      tmp_list$block <- tmp_partable$block
    } else {
      tmp_list$block <- rep(1L, length(tmp_list$lhs))
    }
    if (!is.null(tmp_partable$level)) {
      tmp_list$level <- tmp_partable$level
    } else {
      tmp_list$level <- rep(1L, length(tmp_list$lhs))
    }
    if (!is.null(tmp_partable$group)) {
      tmp_list$group <- tmp_partable$group
    } else {
      tmp_list$group <- rep(1L, length(tmp_list$lhs))
    }
    if (!is.null(tmp_partable$step)) {
      tmp_list$step <- tmp_partable$step
    }
    if (!is.null(tmp_partable$efa)) {
      tmp_list$efa <- tmp_partable$efa
    }
    if (!is.null(tmp_partable$label)) {
      tmp_list$label <- tmp_partable$label
    } else {
      tmp_list$label <- rep("", length(tmp_list$lhs))
    }
    if (!is.null(tmp_partable$exo)) {
      tmp_list$exo <- tmp_partable$exo
    } else {
      tmp_list$exo <- rep(0L, length(tmp_list$lhs))
    }
    if (inherits(object, "lavaanList")) {
      # per default: nothing!
      # if("partable" %in% object@meta$store.slots) {
      #    COF <- sapply(object@ParTableList, "[[", "est")
      #    tmp_list$est <- rowMeans(COF)
      # }
      tmp_list$est <- NULL
    } else if (!is.null(tmp_partable$est)) {
      tmp_list$est <- tmp_partable$est
    } else {
      tmp_list$est <- lav_model_get_parameters(object@Model,
        type = "user",
        extra = TRUE
      )
    }
    if (!is.null(tmp_partable$lower)) {
      tmp_list$lower <- tmp_partable$lower
      tmp_list$lower[tmp_list$free > 0L] <-
        tmp_list$lower[tmp_list$free > 0L] + .Machine$double.eps
    }
    if (!is.null(tmp_partable$upper)) {
      tmp_list$upper <- tmp_partable$upper
      tmp_list$upper[tmp_list$free > 0L] <-
        tmp_list$upper[tmp_list$free > 0L] - .Machine$double.eps
    }


    # add se, zstat, pvalue
    if (se && object@Options$se != "none") {
      tmp_list$se <- lav_inspect_se(object)
      # handle tiny SEs
      tmp_list$se <- ifelse(tmp_list$se < sqrt(.Machine$double.eps),
        0, tmp_list$se
      )
      tmp_se <- ifelse(tmp_list$se < sqrt(.Machine$double.eps), NA, tmp_list$se)
      if (zstat) {
        tmp_list$z <- tmp_list$est / tmp_se
        if (pvalue) {
          tmp_list$pvalue <- 2 * (1 - pnorm(abs(tmp_list$z)))
          # remove p-value if bounds have been used
          if (!is.null(tmp_partable$lower)) {
            b_idx <- which(abs(tmp_partable$lower - tmp_partable$est) <
              sqrt(.Machine$double.eps) &
              tmp_partable$free > 0L)
            if (length(b_idx) > 0L) {
              tmp_list$pvalue[b_idx] <- as.numeric(NA)
            }
          }
          if (!is.null(tmp_partable$upper)) {
            b_idx <- which(abs(tmp_partable$upper - tmp_partable$est) <
              sqrt(.Machine$double.eps) &
              tmp_partable$free > 0L)
            if (length(b_idx) > 0L) {
              tmp_list$pvalue[b_idx] <- as.numeric(NA)
            }
          }
        }
      }
    }

    # extract bootstrap data (if any)
    if (object@Options$se == "bootstrap" ||
      "bootstrap" %in% object@Options$test ||
      "bollen.stine" %in% object@Options$test) {
      tmp_boot <- lav_inspect_boot(object)
      bootstrap_seed <- attr(tmp_boot, "seed") # for bca
      error_idx <- attr(tmp_boot, "error.idx")
      if (length(error_idx) > 0L) {
        tmp_boot <- tmp_boot[-error_idx, , drop = FALSE] # drops attributes
      }
    } else {
      tmp_boot <- NULL
    }

    bootstrap_successful <- NROW(tmp_boot) # should be zero if NULL

    # confidence interval
    if (se && object@Options$se != "none" && ci) {
      # next three lines based on confint.lm
      a <- (1 - level) / 2
      a <- c(a, 1 - a)
      if (object@Options$se != "bootstrap") {
        fac <- qnorm(a)
        ci <- tmp_list$est + tmp_list$se %o% fac
        # Monte Carlo: replace CI for defined parameters with
        # asymmetric percentile-based bounds (Preacher & Selig 2012)
        mc_coef <- lav_inspect_mc(object)
        def_idx <- which(object@ParTable$op == ":=")
        if (!is.null(mc_coef) && length(def_idx) > 0L) {
          mc_def <- apply(mc_coef, 1L, object@Model@def.function)
          if (length(def_idx) == 1L) {
            mc_def <- as.matrix(mc_def)
          } else {
            mc_def <- t(mc_def)
          }
          mc_def[!is.finite(mc_def)] <- as.numeric(NA)
          alpha <- c((1 - level) / 2, 1 - (1 - level) / 2)
          mc_qq <- apply(mc_def, 2L, quantile, probs = alpha,
                         na.rm = TRUE, type = 7)
          ci[def_idx, ] <- t(mc_qq)
        }
      } else if (object@Options$se == "bootstrap") {
        # bootstrap confidence interval; the order-statistic interpolation
        # (boot package 'norm.inter') is shared with lavEffects() via
        # lav_bootstrap_norm_inter() (see lav_bootstrap.R)
        stopifnot(!is.null(tmp_boot))
        stopifnot(boot.ci.type %in% c(
          "norm", "basic", "perc",
          "bca.simple", "bca"
        ))
        if (boot.ci.type == "norm") {
          fac <- qnorm(a)
          boot_x <- colMeans(tmp_boot, na.rm = TRUE)
          boot_est <-
            lav_model_get_parameters(object@Model,
              glist = lav_model_x2glist(object@Model, boot_x),
              type = "user", extra = TRUE
            )
          bias_est <- (boot_est - tmp_list$est)
          ci <- (tmp_list$est - bias_est) + tmp_list$se %o% fac
        } else if (boot.ci.type %in% c("basic", "perc", "bca.simple")) {
          # the per-quantity bound computation is shared with lavEffects()
          # via lav_bootstrap_ci() (see lav_bootstrap.R)
          ci <- cbind(tmp_list$est, tmp_list$est)

          # free parameters
          free_idx <- which(object@ParTable$free &
            !duplicated(object@ParTable$free))
          ci[free_idx, ] <- lav_bootstrap_ci(
            boot_t = tmp_boot, t0 = tmp_list$est[free_idx],
            boot.ci.type = boot.ci.type, level = level
          )

          # defined parameters
          def_idx <- which(object@ParTable$op == ":=")
          if (length(def_idx) > 0L) {
            boot_def <- apply(tmp_boot, 1L, object@Model@def.function)
            if (length(def_idx) == 1L) {
              boot_def <- as.matrix(boot_def)
            } else {
              boot_def <- t(boot_def)
            }
            ci[def_idx, ] <- lav_bootstrap_ci(
              boot_t = boot_def, t0 = tmp_list$est[def_idx],
              boot.ci.type = boot.ci.type, level = level
            )
          }

          # TODO: add cin/ceq?
        } else if (boot.ci.type == "bca") { # new in 0.6-12
          # we assume that the 'ordinary' (nonparametric) bootstrap was used.
          # The empirical-influence design matrix and the acceleration/CI
          # computations are shared with lavEffects() (see lav_bootstrap.R).
          design <- lav_bootstrap_bca_design(object, error_idx)
          stopifnot(nrow(design) == nrow(tmp_boot))

          ci <- cbind(tmp_list$est, tmp_list$est)

          # free parameters
          free_idx <- which(object@ParTable$free &
            !duplicated(object@ParTable$free))
          stopifnot(length(free_idx) == ncol(tmp_boot))
          acc <- lav_bootstrap_acceleration(tmp_boot, design)
          ci[free_idx, ] <- lav_bootstrap_ci(
            boot_t = tmp_boot, t0 = tmp_list$est[free_idx],
            boot.ci.type = "bca", level = level, acc = acc
          )

          # defined parameters
          def_idx <- which(object@ParTable$op == ":=")
          if (length(def_idx) > 0L) {
            boot_def <- apply(tmp_boot, 1L, object@Model@def.function)
            if (length(def_idx) == 1L) {
              boot_def <- as.matrix(boot_def)
            } else {
              boot_def <- t(boot_def)
            }
            acc <- lav_bootstrap_acceleration(boot_def, design)
            ci[def_idx, ] <- lav_bootstrap_ci(
              boot_t = boot_def, t0 = tmp_list$est[def_idx],
              boot.ci.type = "bca", level = level, acc = acc
            )
          }

          # TODO: add cin/ceq?
        }
      }

      tmp_list$ci.lower <- ci[, 1]
      tmp_list$ci.upper <- ci[, 2]
    }

    # standardized estimates?
    # 28 March 2024: TDJ adds option to select specific types
    # observed variables to standardize when standardized= is a vector of
    # (observed) variable names (a generalization of std.nox); see below
    ov_std_user <- NULL
    if (is.logical(standardized)) {
      if (standardized) {
        standardized <- c("std.lv", "std.all")
        if (length(lav_object_vnames(object, "ov.x")) &&
                                 object@Options$fixed.x) {
          standardized <- c(standardized, "std.nox")
        }
      } else {
        standardized <- character(0)
      } # corresponds to standardized=FALSE
    } else {
      # !is.logical(standardized)
      standardized <- as.character(standardized)
      std_types <- c("std.lv", "std.all", "std.nox")
      if (length(standardized) > 0L &&
          !any(tolower(standardized) %in% std_types)) {
        # interpret 'standardized' as a vector of observed variable names:
        # standardize ONLY the parameters involving these variables (a
        # generalization of "std.nox", where the exogenous 'x' are the ones
        # left unstandardized)
        ov_all <- lav_object_vnames(object, "ov")
        bad <- standardized[!(standardized %in% ov_all)]
        if (length(bad) > 0L) {
          lav_msg_warn(gettextf(
            "standardized= argument contains unknown observed variable(s): %s",
            paste(bad, collapse = ", ")))
        }
        ov_std_user <- standardized[standardized %in% ov_all]
        if (length(ov_std_user) == 0L) {
          standardized <- character(0)
        } else {
          standardized <- c("std.lv", "std.user")
        }
      } else {
        standardized <- tolower(standardized)
        if ("std.nox" %in% standardized) {
          # sanity checks
          if (length(lav_object_vnames(object, "ov.x")) == 0) {
            lav_msg_note(gettext(
              "`std.nox' unavailable without fixed exogenous predictors"))
            standardized <- setdiff(standardized, "std.nox")
          } else if (!object@Options$fixed.x) {
            lav_msg_note(gettext("`std.nox' unavailable when fixed.x = FALSE"))
            standardized <- setdiff(standardized, "std.nox")
          }
        }
      }
    }
    # Then add each requested type
    # (original source code, but now independently conditional)
    if ("std.lv" %in% standardized) {
      tmp_list$std.lv <- lav_standardize_lv(object, cov_std = cov.std)
    }
    if ("std.all" %in% standardized) {
      tmp_list$std.all <- lav_standardize_all(object,
        est_std = tmp_list$est.std,
        cov_std = cov.std
      )
    }
    if ("std.nox" %in% standardized) {
      tmp_list$std.nox <- lav_standardize_all_nox(object,
        est_std = tmp_list$est.std,
        cov_std = cov.std
      )
    }
    if ("std.user" %in% standardized) {
      tmp_list$std.user <- lav_standardize_all(object,
        est_std = tmp_list$est.std,
        cov_std = cov.std, ov_std = ov_std_user
      )
    }

    # warn if a defined (:=) parameter standardizes a label that is constrained
    # equal across groups/equations (its standardized value is then ambiguous);
    # base the numeric comparison on whichever standardized column is available
    if (length(standardized) > 0L) {
      tmp_std <- tmp_list$std.all
      if (is.null(tmp_std)) tmp_std <- tmp_list$std.nox
      if (is.null(tmp_std)) tmp_std <- tmp_list$std.user
      if (is.null(tmp_std)) tmp_std <- tmp_list$std.lv
      lav_partable_check_std_def(tmp_list, est.std = tmp_std)
    }

    # rsquare?
    if (rsquare) {
      r2 <- lavTech(object, "rsquare", add.labels = TRUE)
      tmp_names <- unlist(lapply(r2, names))
      nel <- length(tmp_names)
      if (nel == 0L) {
        lav_msg_warn(
          gettext("rsquare = TRUE, but there are no dependent variables"))
      } else {
        if (lav_pt_nlevels(tmp_list) == 1L) {
          block <- rep(seq_along(r2), sapply(r2, length))
          first_block_idx <- which(!duplicated(tmp_list$block) &
            tmp_list$block > 0L)
          gval <- tmp_list$group[first_block_idx]
          if (length(gval) > 0L) {
            group <- rep(gval, sapply(r2, length))
          } else {
            # single block, single group
            group <- rep(1L, length(block))
          }
          r2 <- data.frame(
            lhs = tmp_names, op = rep("r2", nel),
            rhs = tmp_names, block = block, group = group,
            est = unlist(r2), stringsAsFactors = FALSE
          )
        } else {
          # add level column
          block <- rep(seq_along(r2), sapply(r2, length))
          first_block_idx <- which(!duplicated(tmp_list$block) &
            tmp_list$block > 0L)
          # always at least two blocks
          gval <- tmp_list$group[first_block_idx]
          group <- rep(gval, sapply(r2, length))
          lval <- tmp_list$level[first_block_idx]
          level <- rep(lval, sapply(r2, length))
          r2 <- data.frame(
            lhs = tmp_names, op = rep("r2", nel),
            rhs = tmp_names,
            block = block, group = group,
            level = level,
            est = unlist(r2), stringsAsFactors = FALSE
          )
        }
        # add step column if needed
        if (!is.null(tmp_list$step)) {
          r2$step <- 2L # per default
          # simplification: we assume that only the
          # observed indicators of latent variables are step 1
          ov_ind <- unlist(object@pta$vnames$ov.ind)
          step1_idx <- which(r2$lhs %in% ov_ind)
          r2$step[step1_idx] <- 1L
        }
        tmp_list <- lav_pt_merge(pt1 = tmp_list, pt2 = r2, warn = FALSE)
      }
    }

    # fractional missing information (if estimator="fiml")
    if (fmi) {
      se_orig <- tmp_list$se

      # new in 0.6-6, use 'EM' based (unstructured) sample statistics
      # otherwise, it would be as if we use expected info, while the
      # original use observed, producing crazy results
      if (object@Data@ngroups > 1L) {
        em_cov <- lapply(lavInspect(object, "sampstat.h1"), "[[", "cov")
        em_mean <- lapply(lavInspect(object, "sampstat.h1"), "[[", "mean")
      } else {
        em_cov <- lavInspect(object, "sampstat.h1")$cov
        em_mean <- lavInspect(object, "sampstat.h1")$mean
      }

      tmp_pt <- parTable(object)
      tmp_pt$ustart <- tmp_pt$est
      tmp_pt$start <- tmp_pt$est <- NULL

      this_options <- object@Options
      if (!is.null(fmi.options) && is.list(fmi.options)) {
        # modify original options
        this_options <- modifyList(this_options, fmi.options)
      }
      # override
      this_options$optim.method <- "none"
      this_options$sample.cov.rescale <- FALSE
      this_options$check.gradient <- FALSE
      this_options$baseline <- FALSE
      this_options$h1 <- FALSE
      this_options$test <- FALSE

      fit_complete <- lavaan(
        model = tmp_pt,
        sample_cov = em_cov,
        sample_mean = em_mean,
        sample_nobs = lavInspect(object, "nobs"),
        slot_options = this_options
      )

      se_comp <- lavParameterEstimates(fit_complete,
        ci = FALSE, fmi = FALSE,
        zstat = FALSE, pvalue = FALSE, remove.system.eq = FALSE,
        remove.eq = FALSE, remove.ineq = FALSE,
        remove.def = FALSE, remove.nonfree = FALSE, remove.unused = FALSE,
        rsquare = rsquare, add.attributes = FALSE
      )$se

      se_comp <- ifelse(se_comp == 0.0, as.numeric(NA), se_comp)
      tmp_list$fmi <- 1 - (se_comp * se_comp) / (se_orig * se_orig)
    }

    # if single level, remove level column
    if (object@Data@nlevels == 1L) tmp_list$level <- NULL

    # if single group, remove group column
    if (object@Data@ngroups == 1L) tmp_list$group <- NULL

    # if single everything, remove block column
    if (object@Data@nlevels == 1L &&
      object@Data@ngroups == 1L) {
      tmp_list$block <- NULL
    }

    # if no user-defined labels, remove label column
    if (sum(nchar(object@ParTable$label)) == 0L) {
      tmp_list$label <- NULL
    }

    # if plabel = TRUE, add it (new in 0.6-20)
    if (plabel) {
      if (!is.null(tmp_partable$plabel)) {
        tmp_list$plabel <- tmp_partable$plabel
      } else {
        if (!is.null(tmp_partable$id)) {
          tmp_list$plabel <- paste(".p", tmp_list$id, ".", sep = "")
        } else {
          tmp_list$plabel <- paste(".p", seq_along(tmp_list$plabel),
                                   ".", sep = "")
        }
      }
    }

    # remove non-free parameters? (but keep ==, >, < and :=)
    if (remove.nonfree) {
      nonfree_idx <- which(tmp_list$free == 0L &
        !tmp_list$op %in% c("==", ">", "<", ":="))
      if (length(nonfree_idx) > 0L) {
        tmp_list <- tmp_list[-nonfree_idx, ]
      }
    }

    # remove 'unused' parameters
    # these are parameters that are automatically added (user == 0),
    # but with their final (est) values fixed to their default values
    # (typically 1 or 0).
    # currently only intercepts and scaling-factors (for now)
    # should we also remove fixed-to-1 variances? (parameterization = theta)?
    if (remove.unused) {
      # intercepts
      int_idx <- which(tmp_list$op == "~1" &
        tmp_list$user == 0L &
        tmp_list$free == 0L &
        tmp_list$est == 0)
      if (length(int_idx) > 0L) {
        tmp_list <- tmp_list[-int_idx, ]
      }

      # scaling factors
      scaling_idx <- which(tmp_list$op == "~*~" &
        tmp_list$user == 0L &
        tmp_list$free == 0L &
        tmp_list$est == 1)
      if (length(scaling_idx) > 0L) {
        tmp_list <- tmp_list[-scaling_idx, ]
      }
    }


    # remove 'free' column
    tmp_list$free <- NULL

    # remove == rows?
    if (remove.eq) {
      eq_idx <- which(tmp_list$op == "==" & tmp_list$user == 1L)
      if (length(eq_idx) > 0L) {
        tmp_list <- tmp_list[-eq_idx, ]
      }
    }
    if (remove.system.eq) {
      eq_idx <- which(tmp_list$op == "==" & tmp_list$user != 1L)
      if (length(eq_idx) > 0L) {
        tmp_list <- tmp_list[-eq_idx, ]
      }
    }
    # remove <> rows?
    if (remove.ineq) {
      ineq_idx <- which(tmp_list$op %in% c("<", ">"))
      if (length(ineq_idx) > 0L) {
        tmp_list <- tmp_list[-ineq_idx, ]
      }
    }
    # remove := rows?
    if (remove.def) {
      def_idx <- which(tmp_list$op == ":=")
      if (length(def_idx) > 0L) {
        tmp_list <- tmp_list[-def_idx, ]
      }
    }

    # remove step 1 rows?
    if (remove.step1 && !is.null(tmp_list$step)) {
      step1_idx <- which(tmp_list$step == 1L)
      if (length(step1_idx) > 0L) {
        tmp_list <- tmp_list[-step1_idx, ]
      }
      # remove step column
      tmp_list$step <- NULL
    }

    # remove auxiliary (saturated-correlates) rows? (FIML aux= variables)
    if (remove.aux && length(object@Options$aux) > 0L) {
      aux_idx <- which(tmp_list$lhs %in% object@Options$aux |
        tmp_list$rhs %in% object@Options$aux)
      if (length(aux_idx) > 0L) {
        tmp_list <- tmp_list[-aux_idx, ]
      }
    }

    # remove attribute for data order
    attr(tmp_list, "ovda") <- NULL

    # remove tmp_list$user
    tmp_list$user <- NULL

    if (output == "text") {
      class(tmp_list) <- c(
        "lavaan.parameterEstimates", "lavaan.data.frame",
        "data.frame"
      )
      if (header) {
        attr(tmp_list, "categorical") <- object@Model@categorical
        attr(tmp_list, "parameterization") <- object@Model@parameterization
        attr(tmp_list, "information") <- object@Options$information[1]
        attr(tmp_list, "information.meat") <- object@Options$information.meat
        attr(tmp_list, "se") <- object@Options$se
        attr(tmp_list, "group.label") <- object@Data@group.label
        attr(tmp_list, "level.label") <- object@Data@level.label
        attr(tmp_list, "bootstrap") <- object@Options$bootstrap
        attr(tmp_list, "bootstrap.successful") <- bootstrap_successful
        attr(tmp_list, "missing") <- object@Options$missing
        attr(tmp_list, "observed.information") <-
          object@Options$observed.information[1]
        attr(tmp_list, "h1.information") <- object@Options$h1.information[1]
        attr(tmp_list, "h1.information.meat") <-
          object@Options$h1.information.meat
        attr(tmp_list, "header") <- header
        # FIXME: add more!!
      }
    } else {
      tmp_list$exo <- NULL
      tmp_list$lower <- tmp_list$upper <- NULL
      class(tmp_list) <- c("lavaan.data.frame", "data.frame")
    }

    tmp_list
  }
parameterEstimates <- lavParameterEstimates     # synonym   # nolint

parameterTable <- parametertable <- parTable <- partable <- # nolint
  function(object) {
    # check object
    object <- lav_object_check_version(object)

    # convert to data.frame
    out <- as.data.frame(object@ParTable, stringsAsFactors = FALSE)

    class(out) <- c("lavaan.data.frame", "data.frame")
    out
  }

varTable <- vartable <- function(object, ov.names = names(object),  # nolint start
                                 ov.names.x = NULL,
                                 ordered = NULL, factor = NULL,
                                 as.data.frame. = TRUE) {           # nolint end

  if (inherits(object, "lavaan")) {
    # check object
    object <- lav_object_check_version(object)
    tmp_var <- object@Data@ov
  } else if (inherits(object, "lavData")) {
    tmp_var <- object@ov
  } else if (inherits(object, "data.frame")) {
    tmp_var <- lav_dataframe_vartable(
      frame = object, ov_names = ov.names,
      ov_names_x = ov.names.x,
      ordered = ordered, factor = factor,
      as_data_frame = FALSE
    )
  } else {
    lav_msg_stop(gettext("object must of class lavaan or a data.frame"))
  }

  if (as.data.frame.) {
    tmp_var <- as.data.frame(tmp_var,
      stringsAsFactors = FALSE,
      row.names = seq_along(tmp_var$name)
    )
    class(tmp_var) <- c("lavaan.data.frame", "data.frame")
  }

  tmp_var
}


setMethod(
  "fitted.values", "lavaan",
  function(object, type = "moments", labels = TRUE, ...) {
    dotdotdot <- list(...)
    if (length(dotdotdot) > 0L) {
      for (j in seq_along(dotdotdot)) {
        lav_msg_warn(gettextf(
          "Unknown argument %s for %s", sQuote(names(dotdotdot)[j]),
          sQuote("fitted.values"))
        )
      }
    }
    # lowercase type
    type <- tolower(type)

    # check object
    object <- lav_object_check_version(object)

    # catch type="casewise"
    if (type %in% c("casewise", "case", "obs", "observations", "ov")) {
      return(lavPredict(object, type = "ov", label = labels))
    }

    lav_inspect_implied(object,
      add_labels = labels, add_class = TRUE,
      drop_list_single_group = TRUE
    )
  }
)


setMethod(
  "fitted", "lavaan",
  function(object, type = "moments", labels = TRUE, ...) {
    dotdotdot <- list(...)
    if (length(dotdotdot) > 0L) {
      for (j in seq_along(dotdotdot)) {
        lav_msg_warn(gettextf(
          "Unknown argument %s for %s", sQuote(names(dotdotdot)[j]),
          sQuote("fitted"))
        )
      }
    }
    fitted.values(object, type = type, labels = labels)
  }
)

setMethod(
  "vcov", "lavaan",
  function(object, type = "free", labels = TRUE, remove.duplicated = FALSE, # nolint start
           standardized = NULL, free.only = TRUE, ...) {                    # nolint end
    dotdotdot <- list(...)
    if (length(dotdotdot) > 0L) {
      for (j in seq_along(dotdotdot)) {
        lav_msg_warn(gettextf(
          "Unknown argument %s for %s", sQuote(names(dotdotdot)[j]),
          sQuote("vcov"))
        )
      }
    }
    # check object
    object <- lav_object_check_version(object)

    # check for convergence first!
    if (object@optim$npar > 0L && !object@optim$converged) {
      lav_msg_stop(gettext("model did not converge"))
    }

    if (object@Options$se == "none") {
      lav_msg_stop(gettext("vcov not available if se=\"none\""))
    }

    ## verify there are any user-defined parameters
    #FIXME? smarter to check @ParTable for $op == ":="?
    if (is.null(formals(object@Model@def.function))) {
      type <- "free" # avoids error in lav_inspect_vcov_def()
    }

    if (!is.null(standardized)) {
      standardized <- tolower(standardized[1])
      stopifnot(standardized %in% c("std.lv", "std.all", "std.nox"))
    }

    if (type == "user" || type == "joint" || type == "all" || type == "full" ||
      type == "complete") {
      if (remove.duplicated) {
        lav_msg_stop(gettext(
          "argument \"remove.duplicated\" not supported if type = \"user\""
        ))
      }
      tmp_varcov <- lav_inspect_vcov_def(object,
        joint = TRUE,
        standardized = !is.null(standardized),
        type = ifelse(is.null(standardized), "std.all", standardized),
        add_labels = labels,
        add_class = TRUE
      )
    } else if (type == "free") {
      tmp_varcov <- lav_inspect_vcov(object,
        standardized = !is.null(standardized),
        type = ifelse(is.null(standardized), "std.all", standardized),
        free_only = free.only,
        add_labels = labels,
        add_class = TRUE,
        remove_duplicated = remove.duplicated
      )
    } else {
      lav_msg_stop(gettext("type argument should be \"user\" or \"free\""))
    }

    tmp_varcov
  }
)


# logLik (so that we can use the default AIC/BIC functions from stats4(
setMethod(
  "logLik", "lavaan",
  function(object, ...) {
    dotdotdot <- list(...)
    if (length(dotdotdot) > 0L) {
      for (j in seq_along(dotdotdot)) {
        lav_msg_warn(gettextf(
          "Unknown argument %s for %s", sQuote(names(dotdotdot)[j]),
          sQuote("logLik"))
        )
      }
    }
    # check object
    object <- lav_object_check_version(object)

    if (object@Options$estimator != "ML") {
      lav_msg_warn(gettext("logLik only available if estimator is ML"))
    }
    if (object@optim$npar > 0L && !object@optim$converged) {
      lav_msg_warn(gettext("model did not converge"))
    }

    # new in 0.6-1: we use the @loglik slot (instead of fitMeasures)
    tmp_logl <- object@loglik
    logl <- tmp_logl$loglik
    attr(logl, "df") <- tmp_logl$npar ### note: must be npar, not df!!
    attr(logl, "nobs") <- tmp_logl$ntotal
    class(logl) <- "logLik"
    logl
  }
)

# nobs
if (!exists("nobs", envir = asNamespace("stats4"))) {
  setGeneric("nobs", function(object, ...) standardGeneric("nobs"))
}
setMethod(
  "nobs", signature(object = "lavaan"),
  function(object, ...) {
    dotdotdot <- list(...)
    if (length(dotdotdot) > 0L) {
      for (j in seq_along(dotdotdot)) {
        lav_msg_warn(gettextf(
          "Unknown argument %s for %s", sQuote(names(dotdotdot)[j]),
          sQuote("nobs"))
        )
      }
    }
    object@SampleStats@ntotal
  }
)

# see: src/library/stats/R/update.R
setMethod(
  "update", signature(object = "lavaan"),
  function(object, model, add, ..., evaluate = TRUE) {
    # check object
    object <- lav_object_check_version(object)

    call <- object@call
    if (is.null(call)) {
      lav_msg_stop(gettext("need an object with call slot"))
    }
    current_call_names <- names(call)
    lavaan_formals <- names(formals(lavaan::lavaan))
    current_formals <- lav_snake_case(current_call_names) %in% lavaan_formals
    current_call_names[current_formals] <-
              lav_snake_case(current_call_names[current_formals])
    names(call) <- current_call_names

    extras <- match.call(expand.dots = FALSE)$...

    if (!missing(model)) {
      # call$formula <- update.formula(formula(object), formula.)
      call$model <- model
    } else if (exists(as.character(call$model))) {
      call$model <- eval(call$model, parent.frame())
    } else if (is.character(call$model)) {
      ## do nothing
      ## call$model <- call$model
    } else {
      call$model <- parTable(object)
      call$model$est <- NULL
      call$model$se <- NULL
    }
    if (!is.null(call$slot_par_table) && is.list(call$model)) {
      call$slot_par_table <- call$model
    }

    if (length(extras) > 0) {
      ## check for call$slot_options conflicts
      if (!is.null(call$slot_options)) {
        same_names <- intersect(names(lavOptions()), names(extras))
        for (i in same_names) {
          call$slot_options[[i]] <- extras[[i]]
          extras[i] <- NULL # not needed if they are in slotOptions
        }
      }
      existing <- !is.na(match(names(extras), names(call)))
      for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
      if (any(!existing)) {
        call <- c(as.list(call), extras[!existing])
        call <- as.call(call)
      }
    }

    if (missing(add) && !evaluate) {
      return(call)
    }
    ## for any of the other 3 scenarios, we need the updated fit

    ## Check if "add" and "model" are both strings; combine them
    if (missing(add)) {
      add_allready_in_partable <- TRUE # because nothing to add
    } else {
      if (is.character(add) && is.character(call$model)) {
        call$model <- c(call$model, add)
        add_allready_in_partable <- TRUE
      } else {
        add_allready_in_partable <- FALSE
      }
    }
    newfit <- eval(call, parent.frame())
    if (add_allready_in_partable && evaluate) {
      return(newfit)
    }

    ## only remaining situations: "add" exists, but either "add" or "model"
    ## is a parameter table, so update the parameter table in the call
    if (!(mode(add) %in% c("list", "character"))) {
      lav_msg_stop(
        gettext("'add' argument must be model syntax or parameter table.
                See ?lav_model_pt help page.")
      )
    }
    tmp_pt <- lav_object_extended(newfit, add = add)@ParTable
    tmp_pt$user <- NULL # get rid of "10" category used in lavTestScore()
    ## group == 0L in new rows
    tmp_pt$group[tmp_pt$group == 0L] <- tmp_pt$block[tmp_pt$group == 0L]
    # tmp_pt$plabel == "" in new rows.  Consequences?
    tmp_pt$est <- NULL
    tmp_pt$se <- NULL
    call$model <- tmp_pt

    if (evaluate) {
      eval(call, parent.frame())
    } else {
      call
    }
  }
)

setMethod(
  "anova", signature(object = "lavaan"),
  function(object, ...) {
    # NOTE: if we add additional arguments, it is not the same generic
    # anova() function anymore, and match.call will be screwed up

    # NOTE: we need to extract the names of the models from match.call here,
    #       otherwise, we lose them in the call stack

    mcall <- match.call(expand.dots = TRUE)
    dots <- list(...)

    modp <- if (length(dots)) {
      sapply(dots, inherits, "lavaan")
    } else {
      logical(0)
    }
    tmp_names <- sapply(as.list(mcall)[c(FALSE, TRUE, modp)], deparse)

    lavTestLRT(object = object, ..., model.names = tmp_names)
  }
)
