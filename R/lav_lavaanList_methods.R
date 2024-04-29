# methods
setMethod(
  "show", "lavaanList",
  function(object) {
    # show only basic information
    lav_lavaanList_short_summary(object, print = TRUE)
  }
)

lav_lavaanList_short_summary <- function(object, print = TRUE) {
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
           nd = 3L) {
    lav_lavaanList_summary(object,
      header = header, estimates = estimates,
      print = print, nd = nd
    )
  }
)

lav_lavaanList_summary <- function(object,
                                   header = TRUE,
                                   estimates = TRUE,
                                   est.bias = TRUE,
                                   se.bias = TRUE,
                                   zstat = TRUE,
                                   pvalue = TRUE,
                                   print = TRUE,
                                   nd = 3L) {
  out <- list()

  if (header) {
    out$header <- lav_lavaanList_short_summary(object, print = print)

    # if(print) {
    #    # show only basic information
    #    lav_lavaanList_short_summary(object)
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
      # zstat = FALSE, pvalue = FALSE, ci = FALSE,
      standardized = FALSE,
      output = output
    )

    # scenario 1: simulation
    if (!is.null(object@meta$lavSimulate)) {
      pe$est.true <- object@meta$est.true
      nel <- length(pe$est.true)

      # EST
      EST <- lav_lavaanList_partable(object, what = "est", type = "all")
      AVE <- rowMeans(EST, na.rm = TRUE)

      # remove things like equality constraints
      if (length(AVE) > nel) {
        AVE <- AVE[seq_len(nel)]
      }
      pe$est.ave <- AVE
      if (est.bias) {
        pe$est.bias <- pe$est.ave - pe$est.true
      }

      # SE?
      if (se.bias) {
        SE.OBS <- apply(EST, 1L, sd, na.rm = TRUE)
        if (length(SE.OBS) > nel) {
          SE.OBS <- SE.OBS[seq_len(nel)]
        }
        pe$se.obs <- SE.OBS
        SE <- lav_lavaanList_partable(object, what = "se", type = "all")
        SE.AVE <- rowMeans(SE, na.rm = TRUE)
        if (length(SE.AVE) > nel) {
          SE.AVE <- SE.AVE[seq_len(nel)]
        }
        pe$se.ave <- SE.AVE
        pe$se.bias <- pe$se.ave - pe$se.obs
      }

      # scenario 2: bootstrap
    } else if (!is.null(object@meta$lavBootstrap)) {
      # print the average value for est
      EST <- lav_lavaanList_partable(object, what = "est", type = "all")
      pe$est.ave <- rowMeans(EST, na.rm = TRUE)

      # scenario 3: multiple imputation
    } else if (!is.null(object@meta$lavMultipleImputation)) {
      # pool est: take the mean
      EST <- lav_lavaanList_partable(object, what = "est", type = "all")
      m <- NCOL(EST)
      pe$est <- rowMeans(EST, na.rm = TRUE)

      # pool se

      # between-imputation variance
      # B.var <- apply(EST, 1L, var)
      est1 <- rowMeans(EST, na.rm = TRUE)
      est2 <- rowMeans(EST^2, na.rm = TRUE)
      B.var <- (est2 - est1 * est1) * m / (m - 1)

      # within-imputation variance
      SE <- lav_lavaanList_partable(object, what = "se", type = "all")
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
      EST <- lav_lavaanList_partable(object, what = "est", type = "all")
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
      EST <- lav_lavaanList_partable(object, what = "est", type = "all")
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
  function(object, type = "free", labels = TRUE) {
    lav_lavaanList_partable(
      object = object, what = "est", type = type,
      labels = labels
    )
  }
)

lav_lavaanList_partable <- function(object, what = "est",
                                    type = "free", labels = TRUE) {
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
