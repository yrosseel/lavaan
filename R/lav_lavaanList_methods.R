# methods
setMethod(
  "show", "lavaanList",
  function(object) {
    # check object
    object <- lav_object_check_version(object)
    # show only basic information
    lav_lavaanlist_short_summary(object, print = TRUE)
  }
)

lav_lavaanlist_short_summary <- function(object, print = TRUE) {
  txt <- sprintf(
    "lavaanList (%s) -- based on %d datasets (%d converged)\n",
    object@version,
    object@meta$ndat,
    sum(object@meta$ok)
  )

  if (print) {
    cat(txt)
  }

  invisible(txt)
}

setMethod(
  "summary", "lavaanList",
  function(object, header = TRUE,
           estimates = TRUE,
           print = TRUE,
           nd = 3L,
           simulate.args = list(est.bias = TRUE,                     # nolint
                                se.bias  = TRUE,
                                prop.sig = TRUE,
                                coverage = TRUE,
                                level    = 0.95,
                                trim     = 0),
           ...) {
    dotdotdot <- list(...)
    if (length(dotdotdot) > 0L) {
      for (j in seq_along(dotdotdot)) {
        lav_msg_warn(gettextf(
          "Unknown argument %s for %s", sQuote(names(dotdotdot)[j]),
          sQuote("summary"))
        )
      }
    }
    lav_lavaanlist_summary(object,
      header = header, estimates = estimates,
      print = print, nd = nd, simulate_args = simulate.args
    )
  }
)

lav_lavaanlist_summary <- function(object,
                                   header = TRUE,
                                   estimates = TRUE,
                                   simulate_args = list(est.bias = TRUE,
                                                        se.bias  = TRUE,
                                                        prop.sig = TRUE,
                                                        coverage = TRUE,
                                                        level    = 0.95,
                                                        trim     = 0),
                                   zstat = TRUE,
                                   pvalue = TRUE,
                                   print = TRUE,
                                   nd = 3L) {
  # check object
  object <- lav_object_check_version(object)

  out <- list()

  if (header) {
    out$header <- lav_lavaanlist_short_summary(object, print = print)

    # if(print) {
    #    # show only basic information
    #    lav_lavaanlist_short_summary(object)
    # }
  }
  if (print) {
    output <- "text"
  } else {
    output <- "data.frame"
  }

  if (estimates && "partable" %in% object@meta$store.slots) {
    pe <- lavParameterEstimates(object,
      se = FALSE,
      remove.system.eq = FALSE, remove.eq = FALSE,
      remove.ineq = FALSE, remove.def = FALSE,
      remove.nonfree = FALSE, remove.unused = FALSE,
      remove.step1 = FALSE, # in case we used sam()
      # zstat = FALSE, pvalue = FALSE, ci = FALSE,
      standardized = FALSE,
      output = output
    )

    # scenario 1: simulation
    if (!is.null(object@meta$lavSimulate)) {

      # default behavior
      sim_args <- list(est.bias = TRUE,
                       se.bias  = TRUE,
                       prop.sig = TRUE,
                       coverage = TRUE,
                       level    = 0.95,
                       trim     = 0)
      sim_args <- modifyList(sim_args, simulate_args)
      #if (!sim.args$est.bias) {
      #  sim.args$se.bias <- FALSE
      #}
      if (!sim_args$se.bias) {
        sim_args$prop.sig <- FALSE
      }

      pe$est.true <- object@meta$est.true
      nel <- length(pe$est.true)

      # always compute EST
      est <- lav_lavaanlist_partable(object, what = "est", type = "all")

      # sometimes compute SE
      if (sim_args$se.bias || sim_args$prop.sig || sim_args$coverage) {
        se <- lav_lavaanlist_partable(object, what = "se", type = "all")
      }

      # est.bias?
      if (sim_args$est.bias) {
        ave_1 <- apply(est, 1L, mean, na.rm = TRUE, trim = sim_args$trim)

        # remove things like equality constraints
        if (length(ave_1) > nel) {
          ave_1 <- ave_1[seq_len(nel)]
        }
        ave_1[!is.finite(ave_1)] <- as.numeric(NA)
        pe$est.ave <- ave_1
        pe$est.bias <- pe$est.ave - pe$est.true

        # FIXME: should we also add bias^2 and MSE?
        # and what about relative bias?
      }


      # SE?
      if (sim_args$se.bias) {
        se_obs_1 <- apply(est, 1L, lav_sample_trimmed_sd, na.rm = TRUE,
                        trim = sim_args$trim)
        if (length(se_obs_1) > nel) {
          se_obs_1 <- se_obs_1[seq_len(nel)]
        }
        se_obs_1[!is.finite(se_obs_1)] <- as.numeric(NA)
        pe$se.obs <- se_obs_1
        se_ave <- apply(se, 1L, mean, na.rm = TRUE, trim = sim_args$trim)
        if (length(se_ave) > nel) {
          se_ave <- se_ave[seq_len(nel)]
        }
        se_ave[!is.finite(se_ave)] <- as.numeric(NA)
        pe$se.ave <- se_ave
        se_obs <- se_obs_1
        se_obs[se_obs < .Machine$double.eps ^ (1 / 3)] <- as.numeric(NA)
        pe$se.bias <- pe$se.ave / se_obs # use ratio!
        pe$se.bias[!is.finite(pe$se.bias)] <- as.numeric(NA)
      }

      if (sim_args$prop.sig) {
        se[se < sqrt(.Machine$double.eps)] <- as.numeric(NA)
        wald_1 <- est / se
        wald <- apply(wald_1, 1L, mean, na.rm = TRUE,
                      trim = sim_args$trim)
        wald[!is.finite(wald)] <- as.numeric(NA)
        pe$wald <- wald
        pval <- 2 * (1 - pnorm(abs(wald_1)))
        propsig <- apply(pval, 1L, function(x) {
                        x_ok <- x[is.finite(x)]
                        nx <- length(x_ok)
                        sum(x_ok < (1 - sim_args$level)) / nx
                       })
        propsig[!is.finite(propsig)] <- as.numeric(NA)
        pe$prop.sig <- propsig
      }

      if (sim_args$coverage) {
        # next three lines based on confint.lm
        a <- (1 - sim_args$level) / 2
        a <- c(a, 1 - a)
        fac <- qnorm(a)
        ci_lower_1 <- est + fac[1] * se
        ci_upper_1 <- est + fac[2] * se
        ci_lower <- apply(ci_lower_1, 1L, mean, na.rm = TRUE,
                          trim = sim_args$trim)
        ci_upper <- apply(ci_upper_1, 1L, mean, na.rm = TRUE,
                          trim = sim_args$trim)
        ci_lower[!is.finite(ci_lower)] <- as.numeric(NA)
        ci_upper[!is.finite(ci_upper)] <- as.numeric(NA)
        pe$ci.lower <- ci_lower
        pe$ci.upper <- ci_upper
        # columnwise comparison
        inside_flag <- (ci_lower_1 <= pe$est.true) & (ci_upper_1 >= pe$est.true)
        coverage <- apply(inside_flag, 1L, mean, na.rm = TRUE)
        coverage[!is.finite(coverage)] <- as.numeric(NA)
        pe$coverage <- coverage
      }

      # if sam(), should we keep or remove the step1 values?
      # keep them for now

      # scenario 2: bootstrap
    } else if (!is.null(object@meta$lavBootstrap)) {
      # print the average value for est
      est <- lav_lavaanlist_partable(object, what = "est", type = "all")
      pe$est.ave <- rowMeans(est, na.rm = TRUE)

      # scenario 3: multiple imputation
    } else if (!is.null(object@meta$lavMultipleImputation)) {
      # pool est: take the mean
      est <- lav_lavaanlist_partable(object, what = "est", type = "all")
      m <- NCOL(est)
      pe$est <- rowMeans(est, na.rm = TRUE)

      # pool se

      # between-imputation variance
      # B.var <- apply(EST, 1L, var)
      est1 <- rowMeans(est, na.rm = TRUE)
      est2 <- rowMeans(est^2, na.rm = TRUE)
      b_var <- (est2 - est1 * est1) * m / (m - 1)

      # within-imputation variance
      se <- lav_lavaanlist_partable(object, what = "se", type = "all")
      w_var <- rowMeans(se^2, na.rm = TRUE)

      # total variance: T.var = W.var + B.var + B.var/m
      pe$se <- sqrt(w_var + b_var + (b_var / m))

      tmp_se <- ifelse(pe$se == 0.0, NA, pe$se)
      if (zstat) {
        pe$z <- pe$est / tmp_se
        if (pvalue) {
          pe$pvalue <- 2 * (1 - pnorm(abs(pe$z)))
        }
      }

      # scenario 4: multiple groups/sets
    } else if (!is.null(object@meta$lavMultipleGroups)) {
      # show individual estimates, for each group
      est <- lav_lavaanlist_partable(object, what = "est", type = "all")
      est <- as.list(as.data.frame(est))
      names(est) <- object@meta$group.label
      attr_1 <- attributes(pe)
      names_1 <- c(names(pe), names(est))
      pe <- c(pe, est)
      attributes(pe) <- attr_1
      names(pe) <- names_1
    } else {
      # scenario 5: just a bunch of fits, using different datasets
      # print the average value for est
      est <- lav_lavaanlist_partable(object, what = "est", type = "all")
      pe$est.ave <- rowMeans(est, na.rm = TRUE)

      # more?
    }

    # remove ==,<,>
    rm_idx <- which(pe$op %in% c("==", "<", ">"))
    if (length(rm_idx) > 0L) {
      pe <- pe[-rm_idx, ]
    }

    out$pe <- pe

    if (print) {
      # print pe?
      print(pe, nd = nd)
    }
  } else {
    cat("available slots (per dataset) are:\n")
    print(object@meta$store.slots)
  }

  invisible(out)
}

setMethod(
  "coef", "lavaanList",
  function(object, type = "free", labels = TRUE, ...) {
    # check object
    object <- lav_object_check_version(object)
    dotdotdot <- list(...)
    if (length(dotdotdot) > 0L) {
      for (j in seq_along(dotdotdot)) {
        lav_msg_warn(gettextf(
          "Unknown argument %s for %s", sQuote(names(dotdotdot)[j]),
          sQuote("coef"))
        )
      }
    }
    lav_lavaanlist_partable(
      object = object, what = "est", type = type,
      labels = labels
    )
  }
)

lav_lavaanlist_partable <- function(object, what = "est",
                                    type = "free", labels = TRUE) {
  # check object
  object <- lav_object_check_version(object)
  if ("partable" %in% object@meta$store.slots) {
    if (what %in% names(object@ParTableList[[1]])) {
      out <- sapply(object@ParTableList, "[[", what)
    } else {
      lav_msg_stop(gettextf(
        "column `%s' not found in the first element of the ParTableList slot.",
        what))
    }
  } else {
    lav_msg_stop(gettext("no ParTable slot stored in lavaanList object"))
  }

  if (type == "user" || type == "all") {
    type <- "user"
    idx <- seq_along(object@ParTable$lhs)
  } else if (type == "free") {
    idx <- which(object@ParTable$free > 0L & !duplicated(object@ParTable$free))
  } else {
    lav_msg_stop(gettext("argument `type' must be one of free or user"))
  }

  out <- out[idx, , drop = FALSE]

  if (labels) {
    rownames(out) <- lav_partable_labels(object@ParTable, type = type)
  }

  out
}
