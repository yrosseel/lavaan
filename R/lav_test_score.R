# classic score test (= Lagrange Multiplier test)
#
# this function can run in two modes:
#
# MODE 1: 'add'
#   add new parameters that are currently not included in de model
#   (aka fixed to zero), but should be released
#
# MODE 2: 'release' (the default)
#   release existing "==" constraints
#
lavTestScore <- function(object, add = NULL, release = NULL,
                         univariate = TRUE, cumulative = FALSE,
                         epc = FALSE, standardized = epc, cov.std = epc,
                         verbose = FALSE, warn = TRUE,
                         information = "expected") {
  # check object
  stopifnot(inherits(object, "lavaan"))
  lavoptions <- object@Options

  if (object@optim$npar > 0L && !object@optim$converged) {
    lav_msg_stop(gettext("model did not converge"))
  }

  # check for inequality constraints
  PT <- object@ParTable
  if (any(PT$op == ">" | PT$op == "<")) {
    lav_msg_stop(gettext(
      "lavTestScore() does not handle inequality constraints (yet)"))
  }

  # check arguments
  if (cumulative) {
    univariate <- TRUE
  }


  # Mode 1: ADDING new parameters
  if (!is.null(add) && all(nchar(add) > 0L)) {
    # check release argument
    if (!is.null(release)) {
      lav_msg_stop(gettext(
        "`add' and `release' arguments cannot be used together."))
    }

    # extend model with extra set of parameters
    FIT <- lav_object_extended(object, add = add)

    score <- lavTech(FIT, "gradient.logl")
    Information <- lavTech(FIT, paste("information", information, sep = "."))

    npar <- object@Model@nx.free
    nadd <- FIT@Model@nx.free - npar

    # R
    R.model <- object@Model@con.jac[, , drop = FALSE]
    if (nrow(R.model) > 0L) {
      R.model <- cbind(R.model, matrix(0, nrow(R.model), ncol = nadd))
      R.add <- cbind(matrix(0, nrow = nadd, ncol = npar), diag(nadd))
      R <- rbind(R.model, R.add)

      Z <- cbind(
        rbind(Information, R.model),
        rbind(t(R.model), matrix(0, nrow(R.model), nrow(R.model)))
      )
      Z.plus <- MASS::ginv(Z)
      J.inv <- Z.plus[1:nrow(Information), 1:nrow(Information)]

      r.idx <- seq_len(nadd) + nrow(R.model)
    } else {
      R <- cbind(matrix(0, nrow = nadd, ncol = npar), diag(nadd))
      J.inv <- MASS::ginv(Information)

      r.idx <- seq_len(nadd)
    }

    # lhs/rhs
    lhs <- lav_partable_labels(FIT@ParTable)[FIT@ParTable$user == 10L]
    op <- rep("==", nadd)
    rhs <- rep("0", nadd)
    Table <- data.frame(
      lhs = lhs, op = op, rhs = rhs,
      stringsAsFactors = FALSE
    )
    class(Table) <- c("lavaan.data.frame", "data.frame")
  } else {
    # MODE 2: releasing constraints

    R <- object@Model@con.jac[, , drop = FALSE]
    if (nrow(R) == 0L) {
      lav_msg_stop(gettext("no equality constraints found in model."))
    }

    score <- lavTech(object, "gradient.logl")
    Information <- lavTech(
      object,
      paste("information", information, sep = ".")
    )
    J.inv <- MASS::ginv(Information) # FIXME: move into if(is.null(release))?
    #                 else written over with Z1.plus if(is.numeric(release))
    # R <- object@Model@con.jac[,]

    if (is.null(release)) {
      # ALL constraints
      r.idx <- seq_len(nrow(R))
    } else if (is.numeric(release)) {
      r.idx <- release
      if (max(r.idx) > nrow(R)) {
        lav_msg_stop(gettextf(
          "maximum constraint number (%1$s) is larger than number of
          constraints (%2$s)", max(r.idx), nrow(R)))
      }

      # neutralize the non-needed constraints
      R1 <- R[-r.idx, , drop = FALSE]
      Z1 <- cbind(
        rbind(Information, R1),
        rbind(t(R1), matrix(0, nrow(R1), nrow(R1)))
      )
      Z1.plus <- MASS::ginv(Z1)
      J.inv <- Z1.plus[1:nrow(Information), 1:nrow(Information)]
    } else if (is.character(release)) {
      lav_msg_stop(gettext("not implemented yet"))
    }

    # lhs/rhs
    eq.idx <- which(object@ParTable$op == "==")
    if (length(eq.idx) > 0L) {
      lhs <- object@ParTable$lhs[eq.idx][r.idx]
      op <- rep("==", length(r.idx))
      rhs <- object@ParTable$rhs[eq.idx][r.idx]
    }
    Table <- data.frame(
      lhs = lhs, op = op, rhs = rhs,
      stringsAsFactors = FALSE
    )
    class(Table) <- c("lavaan.data.frame", "data.frame")
  }

  if (object@Data@nlevels == 1L) {
    N <- object@SampleStats@ntotal
    if (lavoptions$mimic == "EQS") {
      N <- N - 1
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

  if (lavoptions$se == "standard") {
    stat <- as.numeric(N * score %*% J.inv %*% score)
  } else {
    # generalized score test
    if (warn) {
      lav_msg_warn(gettext("se is not `standard'; not implemented yet;
                           falling back to ordinary score test"))
    }

    # NOTE!!!
    # we can NOT use VCOV here, because it reflects the constraints,
    # and the whole point is to test for these constraints...

    stat <- as.numeric(N * score %*% J.inv %*% score)
  }

  # compute df, taking into account that some of the constraints may
  # be needed to identify the model (and hence Information is singular)
  # Information.plus <- Information + crossprod(R)
  # df <- qr(R[r.idx,,drop = FALSE])$rank +
  #          ( qr(Information)$rank - qr(Information.plus)$rank )
  df <- nrow(R[r.idx, , drop = FALSE])
  pvalue <- 1 - pchisq(stat, df = df)

  # total score test
  TEST <- data.frame(
    test = "score", X2 = stat, df = df, p.value = pvalue,
    stringsAsFactors = FALSE
  )
  class(TEST) <- c("lavaan.data.frame", "data.frame")
  attr(TEST, "header") <- "total score test:"

  OUT <- list(test = TEST)

  if (univariate) {
    TS <- numeric(nrow(R))
    EPC.uni <- numeric(nrow(R)) # ignored in release= mode
    for (r in r.idx) {
      R1 <- R[-r, , drop = FALSE]
      Z1 <- cbind(
        rbind(Information, R1),
        rbind(t(R1), matrix(0, nrow(R1), nrow(R1)))
      )
      Z1.plus <- MASS::ginv(Z1)
      Z1.plus1 <- Z1.plus[1:nrow(Information), 1:nrow(Information)]
      TS[r] <- as.numeric(N * t(score) %*% Z1.plus1 %*% score)
      if (epc && !is.null(add)) {
        # EPC.uni[r] <- -1 * utils::tail(as.numeric(score %*%  Z1.plus1),
        #                               n = nrow(R))[r]
        # to keep the 'sign' consistent with modindices(), which
        # uses epc = 'new - old'
        EPC.uni[r] <- +1 * utils::tail(as.numeric(score %*% Z1.plus1),
          n = nrow(R)
        )[r]
      }
    }

    Table2 <- Table
    Table2$X2 <- TS[r.idx]
    Table2$df <- rep(1, length(r.idx))
    Table2$p.value <- 1 - pchisq(Table2$X2, df = Table2$df)
    if (epc && !is.null(add)) {
      Table2$epc <- EPC.uni[r.idx]
    }
    attr(Table2, "header") <- "univariate score tests:"
    OUT$uni <- Table2
  }

  if (cumulative) {
    TS.order <- sort.int(TS, index.return = TRUE, decreasing = TRUE)$ix
    ROW.order <- sort.int(TS[r.idx], index.return = TRUE, decreasing = TRUE)$ix
    TS <- numeric(length(r.idx))
    for (r in 1:length(r.idx)) {
      rcumul.idx <- TS.order[1:r]

      R1 <- R[-rcumul.idx, , drop = FALSE]
      Z1 <- cbind(
        rbind(Information, R1),
        rbind(t(R1), matrix(0, nrow(R1), nrow(R1)))
      )
      Z1.plus <- MASS::ginv(Z1)
      Z1.plus1 <- Z1.plus[1:nrow(Information), 1:nrow(Information)]
      TS[r] <- as.numeric(N * t(score) %*% Z1.plus1 %*% score)
    }

    Table3 <- Table[ROW.order, ]
    Table3$X2 <- TS
    Table3$df <- seq_len(length(TS))
    Table3$p.value <- 1 - pchisq(Table3$X2, df = Table3$df)
    attr(Table3, "header") <- "cumulative score tests:"
    OUT$cumulative <- Table3
  }

  if (epc) {
    # EPC <- vector("list", length = length(r.idx))
    # for(i in 1:length(r.idx)) {
    #    r <- r.idx[i]
    #    R1 <- R[-r,,drop = FALSE]
    #    Z1 <- cbind( rbind(Information, R1),
    #                 rbind(t(R1), matrix(0,nrow(R1),nrow(R1))) )
    #    Z1.plus <- MASS::ginv(Z1)
    #    Z1.plus1 <- Z1.plus[ 1:nrow(Information), 1:nrow(Information) ]
    #    EPC[[i]] <- -1 * as.numeric(score %*%  Z1.plus1)
    # }
    #
    # OUT$EPC <- EPC

    # alltogether
    R1 <- R[-r.idx, , drop = FALSE]
    Z1 <- cbind(
      rbind(Information, R1),
      rbind(t(R1), matrix(0, nrow(R1), nrow(R1)))
    )
    Z1.plus <- MASS::ginv(Z1)
    Z1.plus1 <- Z1.plus[1:nrow(Information), 1:nrow(Information)]
    # EPC.all <- -1 * as.numeric(score %*%  Z1.plus1)
    # to keep the 'sign' consistent with modindices(), which
    # uses epc = 'new - old'
    EPC.all <- +1 * as.numeric(score %*% Z1.plus1)

    # create epc table for the 'free' parameters
    if (!is.null(add) && all(nchar(add) > 0L)) {
      LIST <- parTable(FIT)
    } else {
      ## release mode
      LIST <- parTable(object)
    }
    if (lav_partable_ngroups(LIST) == 1L) {
      LIST$group <- NULL
    }
    nonpar.idx <- which(LIST$op %in% c("==", ":=", "<", ">"))
    if (length(nonpar.idx) > 0L) {
      LIST <- LIST[-nonpar.idx, ]
    }

    LIST$est[LIST$free > 0 & LIST$user != 10] <- lav_object_inspect_coef(object, type = "free")
    LIST$est[LIST$user == 10L] <- 0
    LIST$epc <- rep(as.numeric(NA), length(LIST$lhs))
    LIST$epc[LIST$free > 0] <- EPC.all
    LIST$epv <- LIST$est + LIST$epc

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

    LIST$free[LIST$user == 10L] <- 0
    LIST$user <- NULL
    # remove some more columns
    LIST$id <- LIST$ustart <- LIST$exo <- LIST$start <- LIST$se <- LIST$prior <- NULL
    if (lav_partable_nblocks(LIST) == 1L) {
      LIST$block <- NULL
      LIST$group <- NULL
      LIST$level <- NULL
    }

    attr(LIST, "header") <- "expected parameter changes (epc) and expected parameter values (epv):"

    OUT$epc <- LIST
  }

  OUT
}
