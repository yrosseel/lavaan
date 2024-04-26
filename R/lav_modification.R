# univariate modification indices
#

modindices <- function(object,
                       standardized = TRUE,
                       cov.std = TRUE,
                       information = "expected",
                       # power statistics?
                       power = FALSE,
                       delta = 0.1,
                       alpha = 0.05,
                       high.power = 0.75,
                       # customize output
                       sort. = FALSE,
                       minimum.value = 0.0,
                       maximum.number = nrow(LIST),
                       free.remove = TRUE,
                       na.remove = TRUE,
                       op = NULL) {
  # check if model has converged
  if (object@optim$npar > 0L && !object@optim$converged) {
    lav_msg_warn(gettext("model did not converge"))
  }

  # not ready for estimator = "PML"
  if (object@Options$estimator == "PML") {
    lav_msg_stop(gettext(
      "modification indices for estimator PML are not implemented yet."))
  }

  # new in 0.6-17: check if the model contains equality constraints
  if (object@Model@eq.constraints) {
    lav_msg_warn(
      gettext("the modindices() function ignores equality constraints;
              use lavTestScore() to assess the impact of releasing one or
              multiple constraints.")
    )
  }

  # sanity check
  if (power) {
    standardized <- TRUE
  }

  # extended list (fixed-to-zero parameters)
  strict.exo <- FALSE
  if (object@Model@conditional.x) {
    strict.exo <- TRUE
  }
  FULL <- lav_partable_full(
    partable = lav_partable_set_cache(object@ParTable, object@pta),
    free = TRUE, start = TRUE,
    strict.exo = strict.exo
  )
  FULL$free <- rep(1L, nrow(FULL))
  FULL$user <- rep(10L, nrow(FULL))

  FIT <- lav_object_extended(object, add = FULL, all.free = TRUE)
  LIST <- FIT@ParTable


  # compute information matrix 'extended model'
  # ALWAYS use *expected* information (for now)
  Information <- lavTech(FIT, paste("information", information, sep = "."))

  # compute gradient 'extended model'
  score <- lavTech(FIT, "gradient.logl")

  # Saris, Satorra & Sorbom 1987
  # partition Q into Q_11, Q_22 and Q_12/Q_21
  # which elements of Q correspond with 'free' and 'nonfree' parameters?
  model.idx <- LIST$free[LIST$free > 0L & LIST$user != 10L]
  extra.idx <- LIST$free[LIST$free > 0L & LIST$user == 10L]

  # catch empty extra.idx (no modification indices!)
  if (length(extra.idx) == 0L) {
    # 2 possibilities: either model is saturated, or we have constraints
    if (object@test[[1]]$df == 0) {
      lav_msg_warn(gettext(
        "list with extra parameters is empty; model is saturated"))
    } else {
      lav_msg_warn(gettext(
        "list with extra parameters is empty; to release equality
        constraints, use lavTestScore()"))
    }
    LIST <- data.frame(
      lhs = character(0), op = character(0),
      rhs = character(0), group = integer(0),
      mi = numeric(0), epc = numeric(0),
      sepc.lv = numeric(0), sepc.all = numeric(0),
      sepc.nox = numeric(0)
    )
    return(LIST)
  }

  # partition
  I11 <- Information[extra.idx, extra.idx, drop = FALSE]
  I12 <- Information[extra.idx, model.idx, drop = FALSE]
  I21 <- Information[model.idx, extra.idx, drop = FALSE]
  I22 <- Information[model.idx, model.idx, drop = FALSE]

  # ALWAYS use *expected* information (for now)
  I22.inv <- try(
    lavTech(object, paste("inverted.information",
      information,
      sep = "."
    )),
    silent = TRUE
  )
  # just in case...
  if (inherits(I22.inv, "try-error")) {
    lav_msg_stop(gettext(
      "could not compute modification indices; information matrix is singular"))
  }

  V <- I11 - I12 %*% I22.inv %*% I21
  V.diag <- diag(V)
  # dirty hack: catch very small or negative values in diag(V)
  # this is needed eg when parameters are not identified if freed-up;
  idx <- which(V.diag < .Machine$double.eps^(1 / 3)) # was 1/2 <0.6-14
  if (length(idx) > 0L) {
    V.diag[idx] <- as.numeric(NA)
  }

  # create and fill in mi
  if (object@Data@nlevels == 1L) {
    N <- object@SampleStats@ntotal
    if (object@Model@estimator %in% ("ML")) {
      score <- -1 * score # due to gradient.logl
    }
  } else {
    # total number of clusters (over groups)
    N <- 0
    for (g in 1:object@SampleStats@ngroups) {
      N <- N + object@Data@Lp[[g]]$nclusters[[2]]
    }
    # score <- score * (2 * object@SampleStats@ntotal) / N
    score <- score / 2 # -2 * LRT
  }
  mi <- numeric(length(score))
  mi[extra.idx] <- N * (score[extra.idx] * score[extra.idx]) / V.diag
  if (length(model.idx) > 0L) {
    mi[model.idx] <- N * (score[model.idx] * score[model.idx]) / diag(I22)
  }

  LIST$mi <- rep(as.numeric(NA), length(LIST$lhs))
  LIST$mi[LIST$free > 0] <- mi

  # handle equality constraints (if any)
  # eq.idx <- which(LIST$op == "==")
  # if(length(eq.idx) > 0L) {
  #    OUT <- lavTestScore(object, warn = FALSE)
  #    LIST$mi[ eq.idx ] <- OUT$uni$X2
  # }

  # scaled?
  # if(length(object@test) > 1L) {
  #    LIST$mi.scaled <- LIST$mi / object@test[[2]]$scaling.factor
  # }

  # EPC
  d <- (-1 * N) * score
  # needed? probably not; just in case
  d[which(abs(d) < 1e-15)] <- 1.0
  LIST$epc[LIST$free > 0] <- mi / d


  # standardize?
  if (standardized) {
    EPC <- LIST$epc

    if (cov.std) {
      # replace epc values for variances by est values
      var.idx <- which(LIST$op == "~~" & LIST$lhs == LIST$rhs &
        LIST$exo == 0L)
      EPC[var.idx] <- LIST$est[var.idx]
    }

    # two problems:
    #   - EPC of variances can be negative, and that is
    #     perfectly legal
    #   - EPC (of variances) can be tiny (near-zero), and we should
    #     not divide by tiny variables
    small.idx <- which(LIST$op == "~~" &
      LIST$lhs == LIST$rhs &
      abs(EPC) < sqrt(.Machine$double.eps))
    if (length(small.idx) > 0L) {
      EPC[small.idx] <- as.numeric(NA)
    }

    # get the sign
    EPC.sign <- sign(LIST$epc)

    LIST$sepc.lv <- EPC.sign * lav_standardize_lv(object,
      partable = LIST,
      est = abs(EPC),
      cov.std = cov.std
    )
    if (length(small.idx) > 0L) {
      LIST$sepc.lv[small.idx] <- 0
    }
    LIST$sepc.all <- EPC.sign * lav_standardize_all(object,
      partable = LIST,
      est = abs(EPC),
      cov.std = cov.std
    )
    if (length(small.idx) > 0L) {
      LIST$sepc.all[small.idx] <- 0
    }
    LIST$sepc.nox <- EPC.sign * lav_standardize_all_nox(object,
      partable = LIST,
      est = abs(EPC),
      cov.std = cov.std
    )
    if (length(small.idx) > 0L) {
      LIST$sepc.nox[small.idx] <- 0
    }
  }

  # power?
  if (power) {
    LIST$delta <- delta
    # FIXME: this is using epc in unstandardized metric
    #        this would be much more useful in standardized metric
    #        we need a lav_standardize_all.reverse function...
    LIST$ncp <- (LIST$mi / (LIST$epc * LIST$epc)) * (delta * delta)
    LIST$power <- 1 - pchisq(qchisq((1.0 - alpha), df = 1),
      df = 1, ncp = LIST$ncp
    )
    LIST$decision <- character(length(LIST$power))

    # five possibilities (Table 6 in Saris, Satorra, van der Veld, 2009)
    mi.significant <- ifelse(1 - pchisq(LIST$mi, df = 1) < alpha,
      TRUE, FALSE
    )
    high.power <- LIST$power > high.power
    # FIXME: sepc.all or epc??
    # epc.high <- abs(LIST$sepc.all) > LIST$delta
    epc.high <- abs(LIST$epc) > LIST$delta

    LIST$decision[which(!mi.significant & !high.power)] <- "(i)"
    LIST$decision[which(mi.significant & !high.power)] <- "**(m)**"
    LIST$decision[which(!mi.significant & high.power)] <- "(nm)"
    LIST$decision[which(mi.significant & high.power &
      !epc.high)] <- "epc:nm"
    LIST$decision[which(mi.significant & high.power &
      epc.high)] <- "*epc:m*"

    # LIST$decision[ which(mi.significant &  high.power) ] <- "epc"
    # LIST$decision[ which(mi.significant & !high.power) ] <- "***"
    # LIST$decision[ which(!mi.significant & !high.power) ] <- "(i)"
  }

  # remove rows corresponding to 'fixed.x' exogenous parameters
  # exo.idx <- which(LIST$exo == 1L & nchar(LIST$plabel) > 0L)
  # if(length(exo.idx) > 0L) {
  #    LIST <- LIST[-exo.idx,]
  # }

  # remove some columns
  LIST$id <- LIST$ustart <- LIST$exo <- LIST$label <- LIST$plabel <- NULL
  LIST$start <- LIST$free <- LIST$est <- LIST$se <- LIST$prior <- NULL
  LIST$upper <- LIST$lower <- NULL

  if (power) {
    LIST$sepc.lv <- LIST$sepc.nox <- NULL
  }

  # create data.frame
  LIST <- as.data.frame(LIST, stringsAsFactors = FALSE)
  class(LIST) <- c("lavaan.data.frame", "data.frame")

  # remove rows corresponding to 'old' free parameters
  if (free.remove) {
    old.idx <- which(LIST$user != 10L)
    if (length(old.idx) > 0L) {
      LIST <- LIST[-old.idx, ]
    }
  }

  # remove rows corresponding to 'equality' constraints
  eq.idx <- which(LIST$op == "==")
  if (length(eq.idx) > 0L) {
    LIST <- LIST[-eq.idx, ]
  }

  # remove even more columns
  LIST$user <- NULL

  # remove block/group/level is only single block
  if (lav_partable_nblocks(LIST) == 1L) {
    LIST$block <- NULL
    LIST$group <- NULL
    LIST$level <- NULL
  }

  # sort?
  if (sort.) {
    LIST <- LIST[order(LIST$mi, decreasing = TRUE), ]
  }
  if (minimum.value > 0.0) {
    LIST <- LIST[!is.na(LIST$mi) & LIST$mi > minimum.value, ]
  }
  if (maximum.number < nrow(LIST)) {
    LIST <- LIST[seq_len(maximum.number), ]
  }
  if (na.remove) {
    idx <- which(is.na(LIST$mi))
    if (length(idx) > 0) {
      LIST <- LIST[-idx, ]
    }
  }
  if (!is.null(op)) {
    idx <- LIST$op %in% op
    if (length(idx) > 0) {
      LIST <- LIST[idx, ]
    }
  }

  # add header
  # TODO: small explanation of the columns in the header?
  #    attr(LIST, "header") <-
  # c("modification indices for newly added parameters only; to\n",
  #   "see the effects of releasing equality constraints, use the\n",
  #   "lavTestScore() function")

  LIST
}

# aliases
modificationIndices <- modificationindices <- modindices
