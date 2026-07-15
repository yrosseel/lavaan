# classic score test (= Lagrange Multiplier test)
#
# this function can run in two modes:
#
# MODE 1: 'add'
#   add new parameters that are currently not included in the model
#   (aka fixed to zero), but should be released
#
# MODE 2: 'release' (the default)
#   release existing "==" constraints
#
# - YR 8 Jan 2026: use ceq.JAC instead of con.jac (as we do not support
#                  inequality constraints anyway, and we do not want to
#                  release the EFA constraints)
# - YR 11 Jul 2026: add the scaled, adjusted and generalized ('robust')
#                   versions of the score test (Satorra, 2000); see
#                   lav_test_satorra2000.R

# For ceq.simple models, the simple equality constraints are absorbed (no '=='
# rows; ceq.JAC is empty). Reconstruct an explicit constraint jacobian in the
# 'unco' parameter space (one +1/-1 row per equality, ref == other), together
# with matching lhs/rhs labels, so that lavTestScore() can release them. The
# 'unco' parameter order matches which(partable$free > 0), i.e. the space used
# by lavTech(object, "gradient.logl") and "information.*".
lav_test_score_ceq_simple <- function(object) {
  pt <- object@ParTable
  free_rows <- which(pt$free > 0L) # one entry per 'unco' parameter
  nunco <- length(free_rows)
  comp_idx <- pt$free[free_rows]   # compact (nx.free) index per unco parameter
  # use plabels (unique per parameter) so the two sides of each equality are
  # distinguishable, exactly as for explicit '==' constraints (ceq.simple=FALSE)
  labels <- pt$plabel[free_rows]

  r_list <- list()
  lhs <- character(0L)
  rhs <- character(0L)
  dup_comp <- unique(comp_idx[duplicated(comp_idx)])
  for (k in dup_comp) {
    members <- which(comp_idx == k) # unco positions sharing free index k
    ref <- members[1L]
    for (o in members[-1L]) {
      row <- numeric(nunco)
      row[ref] <- 1
      row[o] <- -1
      r_list[[length(r_list) + 1L]] <- row
      lhs <- c(lhs, labels[ref])
      rhs <- c(rhs, labels[o])
    }
  }
  r_mat <- if (length(r_list) > 0L) {
    do.call(rbind, r_list)
  } else {
    matrix(0, 0L, nunco)
  }
  list(R = r_mat, lhs = lhs, op = rep("==", length(lhs)), rhs = rhs)
}

# the [1:npar, 1:npar] block of the inverse of the information matrix
# bordered with the constraint rows in r1 (the constraints that are KEPT,
# i.e. not released): the 'constrained' inverse that neutralizes the
# released constraints only.
#
# computed via the null-space identity
#   J.inv = Z solve(Z' I Z) Z',  Z = orthonormal basis of null(r1),
# which equals the [1:npar, 1:npar] block of
#   ginv(rbind(cbind(I, t(r1)), cbind(r1, 0)))
# (see also lav_model_info_augment_invert); one solve of order
# npar - nrow(r1) instead of the SVD of the (npar + nrow(r1))^2 bordered
# matrix -- this matters when the constraints are released one at a time.
# if Z'IZ cannot be solved (the model is not identified even with the
# kept constraints), we fall back to the explicit bordered Moore-Penrose
# route
lav_test_score_iinv <- function(information = NULL, r1 = NULL) {
  if (is.null(r1) || nrow(r1) == 0L) {
    return(MASS::ginv(information))
  }
  z <- lav_mat_ortho_complement(t(r1))
  ziz <- crossprod(z, information %*% z)
  out <- try(z %*% solve(ziz, t(z)), silent = TRUE)
  if (inherits(out, "try-error")) {
    z1 <- cbind(
      rbind(information, r1),
      rbind(t(r1), matrix(0, nrow(r1), nrow(r1)))
    )
    z1_plus <- MASS::ginv(z1)
    out <- z1_plus[seq_len(nrow(information)), seq_len(nrow(information))]
  } else {
    out <- (out + t(out)) / 2
  }
  out
}

lavTestScore <- function(object, add = NULL, release = NULL,       # nolint
                         univariate = TRUE, cumulative = FALSE,
                         epc = FALSE, standardized = epc, cov_std = epc,
                         verbose = FALSE, warn = TRUE,
                         information = "expected", ...) {
  dotdotdot <- list(...)
  lav_adapt_func(environment(), dotdotdot, NULL)
  # check object
  object <- lav_object_check_version(object)

  # sam objects: the score test needs the (joint) information/scores, which
  # a two-step (sam) procedure does not provide (issue #517)
  if (!is.null(object@internal$sam.method)) {
    lav_msg_stop(gettext(
      "lavTestScore() is not available for models fitted with sam();
       consider modindices() to inspect the structural part."))
  }

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

  # PML: the expected (h1) information is not available; use the
  # observed (pairwise) information instead
  if (object@Model@estimator == "PML" && information == "expected") {
    information <- "observed"
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

    # the object providing the ingredients for the robust versions
    # (Satorra, 2000): score/information/meat all live in the parameter
    # space of the extended model
    robust_object <- fit

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

      # keep the model constraints
      j_inv <- lav_test_score_iinv(information_1, r_model)

      r_idx <- seq_len(nadd) + nrow(r_model)
    } else {
      r_1 <- cbind(matrix(0, nrow = nadd, ncol = npar), diag(nadd))
      j_inv <- lav_test_score_iinv(information_1)

      r_idx <- seq_len(nadd)
    }

    # lhs/rhs
    lhs <- lav_pt_labels(fit@ParTable)[fit@ParTable$user == 10L]
    op <- rep("==", nadd)
    rhs <- rep("0", nadd)
    table_1 <- data.frame(
      lhs = lhs, op = op, rhs = rhs,
      stringsAsFactors = FALSE
    )
    class(table_1) <- c("lavaan.data.frame", "data.frame")
  } else {
    # MODE 2: releasing constraints

    robust_object <- object

    r_1 <- object@Model@ceq.JAC[, , drop = FALSE]
    ceq_simple_con <- NULL
    if (nrow(r_1) == 0L && object@Model@ceq.simple.only) {
      # simple equality constraints (ceq.simple) are absorbed: ceq.JAC is empty.
      # reconstruct an explicit constraint jacobian in the 'unco' parameter
      # space (one +1/-1 row per equality), matching the space of the score and
      # information below, so the constraints can be released individually.
      ceq_simple_con <- lav_test_score_ceq_simple(object)
      r_1 <- ceq_simple_con$R
    }
    if (nrow(r_1) == 0L) {
      lav_msg_stop(gettext("no equality constraints found in model."))
    }

    score <- lavTech(object, "gradient.logl")
    information_1 <- lavTech(
      object,
      paste("information", information, sep = ".")
    )
    # R <- object@Model@con.jac[,]

    if (is.null(release)) {
      # ALL constraints
      r_idx <- seq_len(nrow(r_1))
      j_inv <- lav_test_score_iinv(information_1)
    } else if (is.numeric(release)) {
      r_idx <- release
      if (max(r_idx) > nrow(r_1)) {
        lav_msg_stop(gettextf(
          "maximum constraint number (%1$s) is larger than number of
          constraints (%2$s)", max(r_idx), nrow(r_1)))
      }

      # neutralize the non-needed constraints
      j_inv <- lav_test_score_iinv(
        information_1, r_1[-r_idx, , drop = FALSE]
      )
    } else if (is.character(release)) {
      lav_msg_stop(gettext("not implemented yet"))
    }

    # lhs/rhs
    if (!is.null(ceq_simple_con)) {
      lhs <- ceq_simple_con$lhs[r_idx]
      op <- rep("==", length(r_idx))
      rhs <- ceq_simple_con$rhs[r_idx]
    } else {
      eq_idx <- which(object@ParTable$op == "==")
      if (length(eq_idx) > 0L) {
        lhs <- object@ParTable$lhs[eq_idx][r_idx]
        op <- rep("==", length(r_idx))
        rhs <- object@ParTable$rhs[eq_idx][r_idx]
      }
    }
    table_1 <- data.frame(
      lhs = lhs, op = op, rhs = rhs,
      stringsAsFactors = FALSE
    )
    class(table_1) <- c("lavaan.data.frame", "data.frame")
  }

  # PML: lavTech(object, "gradient.logl") returns the TOTAL (summed over
  # cases) pairwise gradient, while the information matrices are on the
  # unit scale (see lav_model_info_firstorder()); rescale the score so
  # that both are on the same (unit) scale
  if (object@Model@estimator == "PML") {
    if (length(object@Data@sampling.weights) == 0L) {
      score <- score / object@SampleStats@ntotal
    } else {
      score <- score / sum(unlist(object@Data@weights))
    }
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

  # the standard (normal-theory) score test statistic
  #
  # NOTE: we can NOT use VCOV here, because it reflects the constraints,
  # and the whole point is to test for these constraints...
  stat <- as.numeric(n * score %*% j_inv %*% score)

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

  # scaled, adjusted and generalized ('robust') versions (Satorra, 2000),
  # whenever the object provides the ingredients for the sandwich: the
  # (unit-scale) meat B, in the same parameter space as information_1
  b_meat <- lav_test_meat(robust_object)
  if (!is.null(b_meat)) {
    a_mat <- r_1[r_idx, , drop = FALSE]
    aj <- a_mat %*% j_inv
    m1 <- aj %*% t(a_mat)              # A J A'
    m2 <- aj %*% b_meat %*% t(aj)      # A J B J A'
    s2000 <- lav_test_satorra2000(
      stat = stat, df = df, m1 = m1, m2 = m2,
      v = as.numeric(aj %*% score), n = n
    )
    test <- rbind(test, data.frame(
      test = c("score.scaled", "score.adjusted", "score.robust"),
      X2 = c(
        s2000$stat.scaled, s2000$stat.adjusted,
        s2000$stat.robust
      ),
      df = c(s2000$df.scaled, s2000$df.adjusted, s2000$df.robust),
      p.value = c(
        s2000$p.value.scaled, s2000$p.value.adjusted,
        s2000$p.value.robust
      ),
      stringsAsFactors = FALSE
    ))
  } else if (!lavoptions$se %in% c("standard", "none")) {
    lav_msg_warn(gettextf(
      "the scaled, adjusted and robust versions of the score test are not
      available when se = %s; only the standard score test is computed.",
      dQuote(lavoptions$se)))
  }

  class(test) <- c("lavaan.data.frame", "data.frame")
  attr(test, "header") <- "total score test:"

  out <- list(test = test)

  if (univariate) {
    ts_1 <- numeric(nrow(r_1))
    ts_scaled <- rep(as.numeric(NA), nrow(r_1))
    epc_uni <- numeric(nrow(r_1)) # ignored in release= mode
    for (r in r_idx) {
      z1_plus1 <- lav_test_score_iinv(
        information_1, r_1[-r, , drop = FALSE]
      )
      ts_1[r] <- as.numeric(n * t(score) %*% z1_plus1 %*% score)
      if (!is.null(b_meat)) {
        # single restriction: the scaled, adjusted and generalized
        # versions all coincide (Satorra, 2000, eq. 25-28)
        a_r <- r_1[r, , drop = FALSE]
        aj_r <- a_r %*% z1_plus1
        m1_r <- as.numeric(aj_r %*% t(a_r))
        m2_r <- as.numeric(aj_r %*% b_meat %*% t(aj_r))
        if (is.finite(m1_r) && is.finite(m2_r) &&
          m1_r > 0 && m2_r > 0) {
          ts_scaled[r] <- ts_1[r] * m1_r / m2_r
        }
      }
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
    if (!is.null(b_meat)) {
      table2$X2.scaled <- ts_scaled[r_idx]
      table2$p.value.scaled <- 1 - pchisq(table2$X2.scaled, df = 1)
    }
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
    ts_scaled <- rep(as.numeric(NA), length(r_idx))
    for (r in seq_along(r_idx)) {
      rcumul_idx <- ts_order[1:r]

      z1_plus1 <- lav_test_score_iinv(
        information_1, r_1[-rcumul_idx, , drop = FALSE]
      )
      ts_1[r] <- as.numeric(n * t(score) %*% z1_plus1 %*% score)
      if (!is.null(b_meat)) {
        a_c <- r_1[rcumul_idx, , drop = FALSE]
        aj_c <- a_c %*% z1_plus1
        m1_c <- aj_c %*% t(a_c)
        m2_c <- aj_c %*% b_meat %*% t(aj_c)
        s2000 <- lav_test_satorra2000(
          stat = ts_1[r], df = r,
          m1 = m1_c, m2 = m2_c
        )
        ts_scaled[r] <- s2000$stat.scaled
      }
    }

    table3 <- table_1[row_order, ]
    table3$X2 <- ts_1
    table3$df <- seq_along(ts_1)
    table3$p.value <- 1 - pchisq(table3$X2, df = table3$df)
    if (!is.null(b_meat)) {
      table3$X2.scaled <- ts_scaled
      table3$p.value.scaled <- 1 - pchisq(table3$X2.scaled, df = table3$df)
    }
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
    z1_plus1 <- lav_test_score_iinv(
      information_1, r_1[-r_idx, , drop = FALSE]
    )
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
    if (lav_pt_ngroups(list_1) == 1L) {
      list_1$group <- NULL
    }
    nonpar_idx <- which(list_1$op %in% c("==", ":=", "<", ">"))
    if (length(nonpar_idx) > 0L) {
      list_1 <- list_1[-nonpar_idx, ]
    }

    list_1$est[list_1$free > 0 & list_1$user != 10] <-
            lav_inspect_coef(object, type = "free")
    list_1$est[list_1$user == 10L] <- 0
    list_1$epc <- rep(as.numeric(NA), length(list_1$lhs))
    list_1$epc[list_1$free > 0] <- epc_all
    list_1$epv <- list_1$est + list_1$epc

    if (standardized) {
      epc_1 <- list_1$epc

      if (cov_std) {
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
        cov_std = cov_std
      )
      if (length(small_idx) > 0L) {
        list_1$sepc.lv[small_idx] <- 0
      }
      list_1$sepc.all <- epc_sign * lav_standardize_all(object,
        partable = list_1,
        est = abs(epc_1),
        cov_std = cov_std
      )
      if (length(small_idx) > 0L) {
        list_1$sepc.all[small_idx] <- 0
      }
      list_1$sepc.nox <- epc_sign * lav_standardize_all_nox(object,
        partable = list_1,
        est = abs(epc_1),
        cov_std = cov_std
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
    if (lav_pt_nblocks(list_1) == 1L) {
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
