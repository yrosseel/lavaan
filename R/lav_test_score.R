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
# - YR 8 Jan 2026: use ceq.JAC instead of con.jac (as we do not support
#                  inequality constraints anyway, and we do not want to
#                  release the EFA constraints)

lavTestScore <- function(object, add = NULL, release = NULL,       # nolint start
                         univariate = TRUE, cumulative = FALSE,
                         epc = FALSE, standardized = epc, cov.std = epc,
                         verbose = FALSE, warn = TRUE,
                         information = "expected") {               # nolint end
  # check object
  object <- lav_object_check_version(object)

  if (!missing(warn)) {
    current_warn <- lav_warn()
    if (lav_warn(warn))
      on.exit(lav_warn(current_warn))
  }
  if (!missing(verbose)) {
    current_verbose <- lav_verbose()
    if (lav_verbose(verbose))
      on.exit(lav_verbose(current_verbose), TRUE)
  }
  # check object
  stopifnot(inherits(object, "lavaan"))
  lavoptions <- object@Options

  if (object@optim$npar > 0L && !object@optim$converged) {
    lav_msg_stop(gettext("model did not converge"))
  }

  # check for inequality constraints
  pt_1 <- object@ParTable
  if (any(pt_1$op == ">" | pt_1$op == "<")) {
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
    fit <- lav_object_extended(object, add = add)

    score <- lavTech(fit, "gradient.logl")
    information_1 <- lavTech(fit, paste("information", information, sep = "."))

    npar <- object@Model@nx.free
    nadd <- fit@Model@nx.free - npar

    # R
    r_model <- object@Model@ceq.JAC[, , drop = FALSE]
    if (nrow(r_model) > 0L) {
      r_model <- cbind(r_model, matrix(0, nrow(r_model), ncol = nadd))
      r_add <- cbind(matrix(0, nrow = nadd, ncol = npar), diag(nadd))
      r_1 <- rbind(r_model, r_add)

      z <- cbind(
        rbind(information_1, r_model),
        rbind(t(r_model), matrix(0, nrow(r_model), nrow(r_model)))
      )
      z_plus <- MASS::ginv(z)
      j_inv <- z_plus[seq_len(nrow(information_1)),
                      seq_len(nrow(information_1))]

      r_idx <- seq_len(nadd) + nrow(r_model)
    } else {
      r_1 <- cbind(matrix(0, nrow = nadd, ncol = npar), diag(nadd))
      j_inv <- MASS::ginv(information_1)

      r_idx <- seq_len(nadd)
    }

    # lhs/rhs
    lhs <- lav_partable_labels(fit@ParTable)[fit@ParTable$user == 10L]
    op <- rep("==", nadd)
    rhs <- rep("0", nadd)
    table_1 <- data.frame(
      lhs = lhs, op = op, rhs = rhs,
      stringsAsFactors = FALSE
    )
    class(table_1) <- c("lavaan.data.frame", "data.frame")
  } else {
    # MODE 2: releasing constraints

    r_1 <- object@Model@ceq.JAC[, , drop = FALSE]
    if (nrow(r_1) == 0L) {
      lav_msg_stop(gettext("no equality constraints found in model."))
    }

    score <- lavTech(object, "gradient.logl")
    information_1 <- lavTech(
      object,
      paste("information", information, sep = ".")
    )
    j_inv <- MASS::ginv(information_1) # FIXME: move into if(is.null(release))?
    #                 else written over with Z1.plus if(is.numeric(release))
    # R <- object@Model@con.jac[,]

    if (is.null(release)) {
      # ALL constraints
      r_idx <- seq_len(nrow(r_1))
    } else if (is.numeric(release)) {
      r_idx <- release
      if (max(r_idx) > nrow(r_1)) {
        lav_msg_stop(gettextf(
          "maximum constraint number (%1$s) is larger than number of
          constraints (%2$s)", max(r_idx), nrow(r_1)))
      }

      # neutralize the non-needed constraints
      r1 <- r_1[-r_idx, , drop = FALSE]
      z1 <- cbind(
        rbind(information_1, r1),
        rbind(t(r1), matrix(0, nrow(r1), nrow(r1)))
      )
      z1_plus <- MASS::ginv(z1)
      j_inv <- z1_plus[seq_len(nrow(information_1)),
                       seq_len(nrow(information_1))]
    } else if (is.character(release)) {
      lav_msg_stop(gettext("not implemented yet"))
    }

    # lhs/rhs
    eq_idx <- which(object@ParTable$op == "==")
    if (length(eq_idx) > 0L) {
      lhs <- object@ParTable$lhs[eq_idx][r_idx]
      op <- rep("==", length(r_idx))
      rhs <- object@ParTable$rhs[eq_idx][r_idx]
    }
    table_1 <- data.frame(
      lhs = lhs, op = op, rhs = rhs,
      stringsAsFactors = FALSE
    )
    class(table_1) <- c("lavaan.data.frame", "data.frame")
  }

  if (object@Data@nlevels == 1L) {
    n <- object@SampleStats@ntotal
    if (lavoptions$mimic == "EQS") {
      n <- n - 1
    }
  } else {
    # total number of clusters (over groups)
    n <- 0
    for (g in 1:object@SampleStats@ngroups) {
      n <- n + object@Data@Lp[[g]]$nclusters[[2]]
    }
    # score <- score * (2 * object@SampleStats@ntotal) / N
    score <- score / 2 # -2 * LRT
  }

  if (lavoptions$se == "standard") {
    stat <- as.numeric(n * score %*% j_inv %*% score)
  } else {
    # generalized score test
    lav_msg_warn(gettext("se is not `standard'; not implemented yet;
                         falling back to ordinary score test"))

    # NOTE!!!
    # we can NOT use VCOV here, because it reflects the constraints,
    # and the whole point is to test for these constraints...

    stat <- as.numeric(n * score %*% j_inv %*% score)
  }

  # compute df, taking into account that some of the constraints may
  # be needed to identify the model (and hence Information is singular)
  # Information.plus <- Information + crossprod(R)
  # df <- qr(R[r_idx,,drop = FALSE])$rank +
  #          ( qr(Information)$rank - qr(Information.plus)$rank )
  df <- nrow(r_1[r_idx, , drop = FALSE])
  pvalue <- 1 - pchisq(stat, df = df)

  # total score test
  test <- data.frame(
    test = "score", X2 = stat, df = df, p.value = pvalue,
    stringsAsFactors = FALSE
  )
  class(test) <- c("lavaan.data.frame", "data.frame")
  attr(test, "header") <- "total score test:"

  out <- list(test = test)

  if (univariate) {
    ts_1 <- numeric(nrow(r_1))
    epc_uni <- numeric(nrow(r_1)) # ignored in release= mode
    for (r in r_idx) {
      r1 <- r_1[-r, , drop = FALSE]
      z1 <- cbind(
        rbind(information_1, r1),
        rbind(t(r1), matrix(0, nrow(r1), nrow(r1)))
      )
      z1_plus <- MASS::ginv(z1)
      z1_plus1 <- z1_plus[seq_len(nrow(information_1)),
                          seq_len(nrow(information_1))]
      ts_1[r] <- as.numeric(n * t(score) %*% z1_plus1 %*% score)
      if (epc && !is.null(add)) {
        # EPC.uni[r] <- -1 * utils::tail(as.numeric(score %*%  Z1.plus1),
        #                               n = nrow(R))[r]
        # to keep the 'sign' consistent with modindices(), which
        # uses epc = 'new - old'
        epc_uni[r] <- +1 * utils::tail(as.numeric(score %*% z1_plus1),
          n = nrow(r_1)
        )[r]
      }
    }

    table2 <- table_1
    table2$X2 <- ts_1[r_idx]
    table2$df <- rep(1, length(r_idx))
    table2$p.value <- 1 - pchisq(table2$X2, df = table2$df)
    if (epc && !is.null(add)) {
      table2$epc <- epc_uni[r_idx]
    }
    attr(table2, "header") <- "univariate score tests:"
    out$uni <- table2
  }

  if (cumulative) {
    ts_order <- sort.int(ts_1, index.return = TRUE, decreasing = TRUE)$ix
    row_order <-
         sort.int(ts_1[r_idx], index.return = TRUE, decreasing = TRUE)$ix
    ts_1 <- numeric(length(r_idx))
    for (r in seq_along(r_idx)) {
      rcumul_idx <- ts_order[1:r]

      r1 <- r_1[-rcumul_idx, , drop = FALSE]
      z1 <- cbind(
        rbind(information_1, r1),
        rbind(t(r1), matrix(0, nrow(r1), nrow(r1)))
      )
      z1_plus <- MASS::ginv(z1)
      z1_plus1 <- z1_plus[seq_len(nrow(information_1)),
                          seq_len(nrow(information_1))]
      ts_1[r] <- as.numeric(n * t(score) %*% z1_plus1 %*% score)
    }

    table3 <- table_1[row_order, ]
    table3$X2 <- ts_1
    table3$df <- seq_along(ts_1)
    table3$p.value <- 1 - pchisq(table3$X2, df = table3$df)
    attr(table3, "header") <- "cumulative score tests:"
    out$cumulative <- table3
  }

  if (epc) {
    # EPC <- vector("list", length = length(r_idx))
    # for(i in seq_along(r_idx)) {
    #    r <- r_idx[i]
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
    r1 <- r_1[-r_idx, , drop = FALSE]
    z1 <- cbind(
      rbind(information_1, r1),
      rbind(t(r1), matrix(0, nrow(r1), nrow(r1)))
    )
    z1_plus <- MASS::ginv(z1)
    z1_plus1 <- z1_plus[seq_len(nrow(information_1)),
                           seq_len(nrow(information_1))]
    # EPC.all <- -1 * as.numeric(score %*%  Z1.plus1)
    # to keep the 'sign' consistent with modindices(), which
    # uses epc = 'new - old'
    epc_all <- +1 * as.numeric(score %*% z1_plus1)

    # create epc table for the 'free' parameters
    if (!is.null(add) && all(nchar(add) > 0L)) {
      list_1 <- parTable(fit)
    } else {
      ## release mode
      list_1 <- parTable(object)
    }
    if (lav_partable_ngroups(list_1) == 1L) {
      list_1$group <- NULL
    }
    nonpar_idx <- which(list_1$op %in% c("==", ":=", "<", ">"))
    if (length(nonpar_idx) > 0L) {
      list_1 <- list_1[-nonpar_idx, ]
    }

    list_1$est[list_1$free > 0 & list_1$user != 10] <-
            lav_object_inspect_coef(object, type = "free")
    list_1$est[list_1$user == 10L] <- 0
    list_1$epc <- rep(as.numeric(NA), length(list_1$lhs))
    list_1$epc[list_1$free > 0] <- epc_all
    list_1$epv <- list_1$est + list_1$epc

    if (standardized) {
      epc_1 <- list_1$epc

      if (cov.std) {
        # replace epc values for variances by est values
        var_idx <- which(list_1$op == "~~" & list_1$lhs == list_1$rhs &
          list_1$exo == 0L)
        epc_1[var_idx] <- list_1$est[var_idx]
      }

      # two problems:
      #   - EPC of variances can be negative, and that is
      #     perfectly legal
      #   - EPC (of variances) can be tiny (near-zero), and we should
      #     not divide by tiny variables
      small_idx <- which(list_1$op == "~~" &
        list_1$lhs == list_1$rhs &
        abs(epc_1) < sqrt(.Machine$double.eps))
      if (length(small_idx) > 0L) {
        epc_1[small_idx] <- as.numeric(NA)
      }

      # get the sign
      epc_sign <- sign(list_1$epc)

      list_1$sepc.lv <- epc_sign * lav_standardize_lv(object,
        partable = list_1,
        est = abs(epc_1),
        cov.std = cov.std
      )
      if (length(small_idx) > 0L) {
        list_1$sepc.lv[small_idx] <- 0
      }
      list_1$sepc.all <- epc_sign * lav_standardize_all(object,
        partable = list_1,
        est = abs(epc_1),
        cov.std = cov.std
      )
      if (length(small_idx) > 0L) {
        list_1$sepc.all[small_idx] <- 0
      }
      list_1$sepc.nox <- epc_sign * lav_standardize_all_nox(object,
        partable = list_1,
        est = abs(epc_1),
        cov.std = cov.std
      )
      if (length(small_idx) > 0L) {
        list_1$sepc.nox[small_idx] <- 0
      }
    }

    list_1$free[list_1$user == 10L] <- 0
    list_1$user <- NULL
    # remove some more columns
    list_1$id <- list_1$ustart <- list_1$exo <-
           list_1$start <- list_1$se <- list_1$prior <- NULL
    if (lav_partable_nblocks(list_1) == 1L) {
      list_1$block <- NULL
      list_1$group <- NULL
      list_1$level <- NULL
    }

    attr(list_1, "header") <-
      "expected parameter changes (epc) and expected parameter values (epv):"

    out$epc <- list_1
  }

  out
}
