# classic Wald test
#
# NOTE: does not handle redundant constraints yet!
#
# - YR 11 Jul 2026: add the scaled and adjusted versions of the Wald test
#                   (Satorra, 2000); see lav_test_satorra2000.R. Note that
#                   the main statistic is already the generalized ('robust')
#                   Wald test whenever the object's se is sandwich-based,
#                   as it is based on the object's vcov.
#

lavTestWald <- function(object, constraints = NULL, verbose = FALSE) { # nolint
  # check object
  object <- lav_object_check_version(object)

  if (!missing(verbose)) {
    current_verbose <- lav_verbose()
    if (lav_verbose(verbose))
      on.exit(lav_verbose(current_verbose), TRUE)
  }
  if (object@optim$npar > 0L && !object@optim$converged) {
    lav_msg_stop(gettext("model did not converge"))
  }

  if (is.null(constraints) || all(nchar(constraints) == 0L)) {
    lav_msg_stop(gettext("constraints are empty"))
  }

  # extract slots
  lavoptions <- object@Options
  lavmodel <- object@Model
  lavpartable <- data.frame(object@ParTable)

  # remove == constraints from parTable
  eq_idx <- which(lavpartable$op == "==")
  if (length(eq_idx) > 0L) {
    lavpartable <- lavpartable[-eq_idx, ]
  }
  partable <- as.list(lavpartable)

  # parse constraints
  flat <- lavParseModelString(constraints, parser = lavoptions$parser)
  con <- attr(flat, "constraints")
  list_1 <- list()
  if (length(con) > 0L) {
    lhs <- unlist(lapply(con, "[[", "lhs"))
    op <- unlist(lapply(con, "[[", "op"))
    rhs <- unlist(lapply(con, "[[", "rhs"))
    list_1$lhs <- c(list_1$lhs, lhs)
    list_1$op <- c(list_1$op, op)
    list_1$rhs <- c(list_1$rhs, rhs)
  } else {
    lav_msg_stop(gettext(
      "no equality constraints found in constraints argument"))
  }

  # theta = free parameters only
  theta <- lav_model_get_parameters(lavmodel)

  # build constraint function
  ceq_function <- lav_pt_con_ceq(
    partable = partable,
    con = list_1, debug = FALSE
  )
  # compute jacobian restrictions
  jac <- try(lav_func_jacobian_complex(func = ceq_function, x = theta),
    silent = TRUE
  )
  if (inherits(jac, "try-error")) { # eg. pnorm()
    jac <- lav_func_jacobian_simple(func = ceq_function, x = theta)
  }

  # check for linear redundant rows in JAC
  out <- lav_mat_rref(t(jac))
  ran_k <- length(out$pivot)
  if (ran_k < nrow(jac)) {
    lav_msg_warn(gettext(
      "Jacobian of constraints is rank deficient. Some constraints",
      " may be redundant, and have been removed."))
  }
  jac <- jac[out$pivot, , drop = FALSE]

  if (lav_verbose()) {
    cat("Restriction matrix (jacobian):\n")
    print(jac)
    cat("\n")
  }

  # linear restriction
  theta_r <- ceq_function(theta)
  theta_r <- theta_r[out$pivot]

  if (lav_verbose()) {
    cat("Restricted theta values:\n")
    print(theta_r)
    cat("\n")
  }

  # get VCOV
  # VCOV <- vcov(object, labels = FALSE)
  # avoid S4 dispatch
  vcov_1 <- lav_inspect_vcov(object,
    standardized = FALSE,
    free_only = TRUE,
    add_labels = FALSE,
    add_class = FALSE,
    remove_duplicated = FALSE
  )

  # the constraint jacobian is in the compact (nx.free) space; when
  # ceq.simple.only, vcov_1 is in the larger 'unco' space, so reduce it
  vcov_1 <- lav_model_vcov_unco_to_free(lavmodel, vcov_1)

  # restricted vcov
  vcov_r <- jac %*% vcov_1 %*% t(jac)

  # fixme: what if VCOV.r is singular?
  # Wald test statistic
  # (this is the generalized/robust Wald test of Satorra (2000) whenever
  # the object's se is sandwich-based, as vcov_1 is then a sandwich)
  wald <- as.numeric(t(theta_r) %*% solve(vcov_r) %*% theta_r)

  # df
  wald_df <- ran_k

  # p-value based on chisq
  wald_pvalue <- 1 - pchisq(wald, df = wald_df)

  # prepare output
  out <- list(
    stat = wald, df = wald_df, p.value = wald_pvalue,
    se = lavoptions$se
  )

  # scaled, adjusted and generalized ('robust') versions (Satorra, 2000);
  # these need both the standard (normal-theory) and a sandwich vcov;
  # not for rotated (EFA/ESEM) solutions, as the stored vcov is
  # post-processed in that case (see step 16)
  se_flavor <- lav_test_robust_se_flavor(object)
  if (!is.null(se_flavor) &&
    !(lavmodel@nefa > 0L && lavoptions$rotation != "none")) {
    # helper: (re)compute the vcov for a given se setting, in the same
    # (nx.free) parameter space as the constraint jacobian
    wald_vcov <- function(se) {
      lavoptions2 <- lavoptions
      lavoptions2$se <- se
      lavoptions2$check.vcov <- FALSE
      lavh1 <- object@h1
      if (length(lavh1) == 0L) {
        lavh1 <- NULL
      }
      vv <- try(
        lav_model_vcov(
          lavmodel = lavmodel,
          lavsamplestats = object@SampleStats,
          lavoptions = lavoptions2,
          lavdata = object@Data,
          lavpartable = object@ParTable,
          lavcache = object@Cache,
          lavimplied = object@implied,
          lavh1 = lavh1
        ),
        silent = TRUE
      )
      if (inherits(vv, "try-error") || is.null(vv)) {
        return(NULL)
      }
      lav_model_vcov_unco_to_free(lavmodel, vv)
    }

    # sandwich vcov (reuse the object's vcov if it is already of this
    # flavor) and normal-theory vcov
    if (lavoptions$se == se_flavor) {
      vcov_rob <- vcov_1
    } else {
      vcov_rob <- wald_vcov(se_flavor)
    }
    if (lavoptions$se == "standard") {
      vcov_nt <- vcov_1
    } else {
      vcov_nt <- wald_vcov("standard")
    }

    if (!is.null(vcov_rob) && !is.null(vcov_nt)) {
      m1 <- jac %*% vcov_nt %*% t(jac)   # A J A' / n
      m2 <- jac %*% vcov_rob %*% t(jac)  # A J B J A' / n
      # the standard (normal-theory) Wald statistic
      wald_nt <- as.numeric(t(theta_r) %*% MASS::ginv(m1) %*% theta_r)
      s2000 <- lav_test_satorra2000(
        stat = wald_nt, df = wald_df, m1 = m1, m2 = m2,
        v = theta_r, n = 1
      )
      out$stat.standard <- wald_nt
      out$stat.scaled <- s2000$stat.scaled
      out$p.value.scaled <- s2000$p.value.scaled
      out$scaling.factor <- s2000$scaling.factor
      out$stat.adjusted <- s2000$stat.adjusted
      out$df.adjusted <- s2000$df.adjusted
      out$p.value.adjusted <- s2000$p.value.adjusted
      out$stat.robust <- s2000$stat.robust
      out$p.value.robust <- s2000$p.value.robust
    }
  }

  out
}
