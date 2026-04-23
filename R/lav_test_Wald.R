# classic Wald test
#
# NOTE: does not handle redundant constraints yet!
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
  ceq_function <- lav_partable_constraints_ceq(
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
  out <- lav_matrix_rref(t(jac))
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
  vcov_1 <- lav_object_inspect_vcov(object,
    standardized = FALSE,
    free.only = TRUE,
    add.labels = FALSE,
    add.class = FALSE,
    remove.duplicated = FALSE
  )

  # restricted vcov
  vcov_r <- jac %*% vcov_1 %*% t(jac)

  # fixme: what if VCOV.r is singular?
  # Wald test statistic
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

  out
}
