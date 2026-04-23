# univariate modification indices
#

modindices <- function(object,                         # nolint start
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
                       maximum.number = nrow(list_1),
                       free.remove = TRUE,
                       na.remove = TRUE,
                       op = NULL) {                    # nolint end
  # rename modified argument
  high_power <- high.power

  # check object
  object <- lav_object_check_version(object)

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
  strict_exo <- FALSE
  if (object@Model@conditional.x) {
    strict_exo <- TRUE
  }
  full <- lav_partable_full(
    partable = lav_partable_set_cache(object@ParTable, object@pta),
    free = TRUE, start = TRUE,
    strict.exo = strict_exo
  )
  full$free <- rep(1L, nrow(full))
  full$user <- rep(10L, nrow(full))

  fit <- lav_object_extended(object, add = full, all.free = TRUE)
  list_1 <- fit@ParTable


  # compute information matrix 'extended model'
  # ALWAYS use *expected* information (for now)
  information_1 <- lavTech(fit, paste("information", information, sep = "."))

  # compute gradient 'extended model'
  score <- lavTech(fit, "gradient.logl")

  # Saris, Satorra & Sorbom 1987
  # partition Q into Q_11, Q_22 and Q_12/Q_21
  # which elements of Q correspond with 'free' and 'nonfree' parameters?
  model_idx <- list_1$free[list_1$free > 0L & list_1$user != 10L]
  extra_idx <- list_1$free[list_1$free > 0L & list_1$user == 10L]

  # catch empty extra.idx (no modification indices!)
  if (length(extra_idx) == 0L) {
    # 2 possibilities: either model is saturated, or we have constraints
    if (object@test[[1]]$df == 0) {
      lav_msg_warn(gettext(
        "list with extra parameters is empty; model is saturated"))
    } else {
      lav_msg_warn(gettext(
        "list with extra parameters is empty; to release equality
        constraints, use lavTestScore()"))
    }
    list_1 <- data.frame(
      lhs = character(0), op = character(0),
      rhs = character(0), group = integer(0),
      mi = numeric(0), epc = numeric(0),
      sepc.lv = numeric(0), sepc.all = numeric(0),
      sepc.nox = numeric(0)
    )
    return(list_1)
  }

  # partition
  i11 <- information_1[extra_idx, extra_idx, drop = FALSE]
  i12 <- information_1[extra_idx, model_idx, drop = FALSE]
  i21 <- information_1[model_idx, extra_idx, drop = FALSE]
  i22 <- information_1[model_idx, model_idx, drop = FALSE]

  # ALWAYS use *expected* information (for now)
  i22_inv <- try(
    lavTech(object, paste("inverted.information",
      information,
      sep = "."
    )),
    silent = TRUE
  )
  # just in case...
  if (inherits(i22_inv, "try-error")) {
    lav_msg_stop(gettext(
      "could not compute modification indices; information matrix is singular"))
  }

  v <- i11 - i12 %*% i22_inv %*% i21
  v_diag <- diag(v)
  # dirty hack: catch very small or negative values in diag(V)
  # this is needed eg when parameters are not identified if freed-up;
  idx <- which(v_diag < .Machine$double.eps^(1 / 3)) # was 1/2 <0.6-14
  if (length(idx) > 0L) {
    v_diag[idx] <- as.numeric(NA)
  }

  # create and fill in mi
  if (object@Data@nlevels == 1L) {
    n <- object@SampleStats@ntotal
    #if (object@Model@estimator %in% ("ML")) {
      score <- -1 * score # due to gradient.logl
    #}
  } else {
    # total number of clusters (over groups)
    n <- 0
    for (g in 1:object@SampleStats@ngroups) {
      n <- n + object@Data@Lp[[g]]$nclusters[[2]]
    }
    # score <- score * (2 * object@SampleStats@ntotal) / N
    score <- score / 2 # -2 * LRT
  }
  mi <- numeric(length(score))
  mi[extra_idx] <- n * (score[extra_idx] * score[extra_idx]) / v_diag
  if (length(model_idx) > 0L) {
    mi[model_idx] <- n * (score[model_idx] * score[model_idx]) / diag(i22)
  }

  list_1$mi <- rep(as.numeric(NA), length(list_1$lhs))
  list_1$mi[list_1$free > 0] <- mi

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
  d <- (-1 * n) * score
  # needed? probably not; just in case
  d[which(abs(d) < 1e-15)] <- 1.0
  list_1$epc[list_1$free > 0] <- mi / d


  # standardize?
  if (standardized) {
    epc <- list_1$epc

    if (cov.std) {
      # replace epc values for variances by est values
      var_idx <- which(list_1$op == "~~" & list_1$lhs == list_1$rhs &
        list_1$exo == 0L)
      epc[var_idx] <- list_1$est[var_idx]
    }

    # two problems:
    #   - EPC of variances can be negative, and that is
    #     perfectly legal
    #   - EPC (of variances) can be tiny (near-zero), and we should
    #     not divide by tiny variables
    small_idx <- which(list_1$op == "~~" &
      list_1$lhs == list_1$rhs &
      abs(epc) < sqrt(.Machine$double.eps))
    if (length(small_idx) > 0L) {
      epc[small_idx] <- as.numeric(NA)
    }

    # get the sign
    epc_sign <- sign(list_1$epc)

    list_1$sepc.lv <- epc_sign * lav_standardize_lv(object,
      partable = list_1,
      est = abs(epc),
      cov.std = cov.std
    )
    if (length(small_idx) > 0L) {
      list_1$sepc.lv[small_idx] <- 0
    }
    list_1$sepc.all <- epc_sign * lav_standardize_all(object,
      partable = list_1,
      est = abs(epc),
      cov.std = cov.std
    )
    if (length(small_idx) > 0L) {
      list_1$sepc.all[small_idx] <- 0
    }
    list_1$sepc.nox <- epc_sign * lav_standardize_all_nox(object,
      partable = list_1,
      est = abs(epc),
      cov.std = cov.std
    )
    if (length(small_idx) > 0L) {
      list_1$sepc.nox[small_idx] <- 0
    }
  }

  # power?
  if (power) {
    list_1$delta <- delta
    # FIXME: this is using epc in unstandardized metric
    #        this would be much more useful in standardized metric
    #        we need a lav_standardize_all.reverse function...
    list_1$ncp <- (list_1$mi / (list_1$epc * list_1$epc)) * (delta * delta)
    list_1$power <- 1 - pchisq(qchisq((1.0 - alpha), df = 1),
      df = 1, ncp = list_1$ncp
    )
    list_1$decision <- character(length(list_1$power))

    # five possibilities (Table 6 in Saris, Satorra, van der Veld, 2009)
    mi_significant <- ifelse(1 - pchisq(list_1$mi, df = 1) < alpha,
      TRUE, FALSE
    )
    high_power <- list_1$power > high_power
    # FIXME: sepc.all or epc??
    # epc.high <- abs(LIST$sepc.all) > LIST$delta
    epc_high <- abs(list_1$epc) > list_1$delta

    list_1$decision[which(!mi_significant & !high_power)] <- "(i)"
    list_1$decision[which(mi_significant & !high_power)] <- "**(m)**"
    list_1$decision[which(!mi_significant & high_power)] <- "(nm)"
    list_1$decision[which(mi_significant & high_power &
      !epc_high)] <- "epc:nm"
    list_1$decision[which(mi_significant & high_power &
      epc_high)] <- "*epc:m*"

    # LIST$decision[ which(mi.significant &  high_power) ] <- "epc"
    # LIST$decision[ which(mi.significant & !high_power) ] <- "***"
    # LIST$decision[ which(!mi.significant & !high_power) ] <- "(i)"
  }

  # remove rows corresponding to 'fixed.x' exogenous parameters
  # exo.idx <- which(LIST$exo == 1L & nchar(LIST$plabel) > 0L)
  # if(length(exo.idx) > 0L) {
  #    LIST <- LIST[-exo.idx,]
  # }

  # remove some columns
  list_1$id <- list_1$ustart <- list_1$exo <- list_1$label <-
                                              list_1$plabel <- NULL
  list_1$start <- list_1$free <- list_1$est <- list_1$se <-
                                               list_1$prior <- NULL
  list_1$upper <- list_1$lower <- NULL

  if (power) {
    list_1$sepc.lv <- list_1$sepc.nox <- NULL
  }

  # create data.frame
  list_1 <- as.data.frame(list_1, stringsAsFactors = FALSE)
  class(list_1) <- c("lavaan.data.frame", "data.frame")

  # remove rows corresponding to 'old' free parameters
  if (free.remove) {
    old_idx <- which(list_1$user != 10L)
    if (length(old_idx) > 0L) {
      list_1 <- list_1[-old_idx, ]
    }
  }

  # remove rows corresponding to 'equality' constraints
  eq_idx <- which(list_1$op == "==")
  if (length(eq_idx) > 0L) {
    list_1 <- list_1[-eq_idx, ]
  }

  # remove even more columns
  list_1$user <- NULL

  # remove block/group/level is only single block
  if (lav_partable_nblocks(list_1) == 1L) {
    list_1$block <- NULL
    list_1$group <- NULL
    list_1$level <- NULL
  }

  # sort?
  if (sort.) {
    list_1 <- list_1[order(list_1$mi, decreasing = TRUE), ]
  }
  if (minimum.value > 0.0) {
    list_1 <- list_1[!is.na(list_1$mi) & list_1$mi > minimum.value, ]
  }
  if (maximum.number < nrow(list_1)) {
    list_1 <- list_1[seq_len(maximum.number), ]
  }
  if (na.remove) {
    idx <- which(is.na(list_1$mi))
    if (length(idx) > 0) {
      list_1 <- list_1[-idx, ]
    }
  }
  if (!is.null(op)) {
    idx <- list_1$op %in% op
    if (length(idx) > 0) {
      list_1 <- list_1[idx, ]
    }
  }

  # add header
  # TODO: small explanation of the columns in the header?
  #    attr(LIST, "header") <-
  # c("modification indices for newly added parameters only; to\n",
  #   "see the effects of releasing equality constraints, use the\n",
  #   "lavTestScore() function")

  list_1
}

# aliases
modificationIndices <- modificationindices <- modindices     # nolint
