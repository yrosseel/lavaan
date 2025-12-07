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
           simulate.args = list(est.bias = TRUE,
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
      print = print, nd = nd, simulate.args = simulate.args
    )
  }
)

lav_lavaanlist_summary <- function(object,
                                   header = TRUE,
                                   estimates = TRUE,
                                   simulate.args = list(est.bias = TRUE,
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
    pe <- parameterEstimates(object,
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
      sim.args <- list(est.bias = TRUE,
                       se.bias  = TRUE,
                       prop.sig = TRUE,
                       coverage = TRUE,
                       level    = 0.95,
                       trim     = 0)
      sim.args <- modifyList(sim.args, simulate.args)
      #if (!sim.args$est.bias) {
      #  sim.args$se.bias <- FALSE
      #}
      if (!sim.args$se.bias) {
        sim.args$prop.sig <- FALSE
      }

      pe$est.true <- object@meta$est.true
      nel <- length(pe$est.true)

      # always compute EST
      EST <- lav_lavaanlist_partable(object, what = "est", type = "all")

      # sometimes compute SE
      if (sim.args$se.bias || sim.args$prop.sig || sim.args$coverage) {
        SE <- lav_lavaanlist_partable(object, what = "se", type = "all")
      }

      # est.bias?
      if (sim.args$est.bias) {
        AVE <- apply(EST, 1L, mean, na.rm = TRUE, trim = sim.args$trim)

        # remove things like equality constraints
        if (length(AVE) > nel) {
          AVE <- AVE[seq_len(nel)]
        }
        AVE[!is.finite(AVE)] <- as.numeric(NA)
        pe$est.ave <- AVE
        pe$est.bias <- pe$est.ave - pe$est.true
 
        # FIXME: should we also add bias^2 and MSE?
        # and what about relative bias?
      }

      
      # SE?
      if (sim.args$se.bias) {
        SE.OBS <- apply(EST, 1L, lav_sample_trimmed_sd, na.rm = TRUE,
                        trim = sim.args$trim)
        if (length(SE.OBS) > nel) {
          SE.OBS <- SE.OBS[seq_len(nel)]
        }
        SE.OBS[!is.finite(SE.OBS)] <- as.numeric(NA)
        pe$se.obs <- SE.OBS
        SE.AVE <- apply(SE, 1L, mean, na.rm = TRUE, trim = sim.args$trim)
        if (length(SE.AVE) > nel) {
          SE.AVE <- SE.AVE[seq_len(nel)]
        }
        SE.AVE[!is.finite(SE.AVE)] <- as.numeric(NA)
        pe$se.ave <- SE.AVE
        se.obs <- SE.OBS
        se.obs[se.obs < .Machine$double.eps^(1/3)] <- as.numeric(NA)
        pe$se.bias <- pe$se.ave / se.obs # use ratio!
        pe$se.bias[!is.finite(pe$se.bias)] <- as.numeric(NA)
      }

      if (sim.args$prop.sig) {
        SE[SE < sqrt(.Machine$double.eps)] <- as.numeric(NA)
        WALD <- EST/SE
        wald <- apply(WALD, 1L, mean, na.rm = TRUE,
                      trim = sim.args$trim)
        wald[!is.finite(wald)] <- as.numeric(NA)
        pe$wald <- wald
        PVAL <- 2 * (1 - pnorm(abs(WALD)))
        propsig <- apply(PVAL, 1L, function(x) {
                        x.ok <- x[is.finite(x)]
                        nx <- length(x.ok)
                        sum(x.ok < (1 - sim.args$level))/nx
                       })
        propsig[!is.finite(propsig)] <- as.numeric(NA)
        pe$prop.sig <- propsig
      }

      if (sim.args$coverage) {
        # next three lines based on confint.lm
        a <- (1 - sim.args$level) / 2
        a <- c(a, 1 - a)
        fac <- qnorm(a)
        CI.LOWER <- EST + fac[1]*SE
        CI.UPPER <- EST + fac[2]*SE
        ci.lower <- apply(CI.LOWER, 1L, mean, na.rm = TRUE,
                          trim = sim.args$trim)
        ci.upper <- apply(CI.UPPER, 1L, mean, na.rm = TRUE,
                          trim = sim.args$trim)
        ci.lower[!is.finite(ci.lower)] <- as.numeric(NA)
        ci.upper[!is.finite(ci.upper)] <- as.numeric(NA)
        pe$ci.lower <- ci.lower
        pe$ci.upper <- ci.upper
        # columnwise comparison
        inside.flag <- (CI.LOWER <= pe$est.true) & (CI.UPPER >= pe$est.true)
        coverage <- apply(inside.flag, 1L, mean, na.rm = TRUE)
        coverage[!is.finite(coverage)] <- as.numeric(NA)
        pe$coverage <- coverage
      }

      # if sam(), should we keep or remove the step1 values?
      # keep them for now

      # scenario 2: bootstrap
    } else if (!is.null(object@meta$lavBootstrap)) {
      # print the average value for est
      EST <- lav_lavaanlist_partable(object, what = "est", type = "all")
      pe$est.ave <- rowMeans(EST, na.rm = TRUE)

      # scenario 3: multiple imputation
    } else if (!is.null(object@meta$lavMultipleImputation)) {
      # pool est: take the mean
      EST <- lav_lavaanlist_partable(object, what = "est", type = "all")
      m <- NCOL(EST)
      pe$est <- rowMeans(EST, na.rm = TRUE)

      # pool se

      # between-imputation variance
      # B.var <- apply(EST, 1L, var)
      est1 <- rowMeans(EST, na.rm = TRUE)
      est2 <- rowMeans(EST^2, na.rm = TRUE)
      B.var <- (est2 - est1 * est1) * m / (m - 1)

      # within-imputation variance
      SE <- lav_lavaanlist_partable(object, what = "se", type = "all")
      W.var <- rowMeans(SE^2, na.rm = TRUE)

      # total variance: T.var = W.var + B.var + B.var/m
      pe$se <- sqrt(W.var + B.var + (B.var / m))

      tmp.se <- ifelse(pe$se == 0.0, NA, pe$se)
      if (zstat) {
        pe$z <- pe$est / tmp.se
        if (pvalue) {
          pe$pvalue <- 2 * (1 - pnorm(abs(pe$z)))
        }
      }

      # scenario 4: multiple groups/sets
    } else if (!is.null(object@meta$lavMultipleGroups)) {
      # show individual estimates, for each group
      EST <- lav_lavaanlist_partable(object, what = "est", type = "all")
      EST <- as.list(as.data.frame(EST))
      ngroups <- length(EST)
      names(EST) <- object@meta$group.label
      ATTR <- attributes(pe)
      NAMES <- c(names(pe), names(EST))
      pe <- c(pe, EST)
      attributes(pe) <- ATTR
      names(pe) <- NAMES
    }

    # scenario 5: just a bunch of fits, using different datasets
    else {
      # print the average value for est
      EST <- lav_lavaanlist_partable(object, what = "est", type = "all")
      pe$est.ave <- rowMeans(EST, na.rm = TRUE)

      # more?
    }

    # remove ==,<,>
    rm.idx <- which(pe$op %in% c("==", "<", ">"))
    if (length(rm.idx) > 0L) {
      pe <- pe[-rm.idx, ]
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
      OUT <- sapply(object@ParTableList, "[[", what)
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
    idx <- 1:length(object@ParTable$lhs)
  } else if (type == "free") {
    idx <- which(object@ParTable$free > 0L & !duplicated(object@ParTable$free))
  } else {
    lav_msg_stop(gettext("argument `type' must be one of free or user"))
  }

  OUT <- OUT[idx, , drop = FALSE]

  if (labels) {
    rownames(OUT) <- lav_partable_labels(object@ParTable, type = type)
  }

  OUT
}
