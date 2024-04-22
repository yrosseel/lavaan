# classic Wald test
#
# NOTE: does not handle redundant constraints yet!
#

lavTestWald <- function(object, constraints = NULL, verbose = FALSE) {
  if (object@optim$npar > 0L && !object@optim$converged) {
    lav_msg_stop(gettext("model did not converge"))
  }

  if (is.null(constraints) || all(nchar(constraints) == 0L)) {
    lav_msg_stop(gettext("constraints are empty"))
  }

  # extract slots
  lavoptions <- object@Options
  lavmodel <- object@Model
  lavpartable <- object@ParTable

  # remove == constraints from parTable
  eq.idx <- which(lavpartable$op == "==")
  if (length(eq.idx) > 0L) {
    lavpartable <- lavpartable[-eq.idx, ]
  }
  partable <- as.list(lavpartable)

  # parse constraints
  FLAT <- lavParseModelString(constraints, parser = lavoptions$parser)
  CON <- attr(FLAT, "constraints")
  LIST <- list()
  if (length(CON) > 0L) {
    lhs <- unlist(lapply(CON, "[[", "lhs"))
    op <- unlist(lapply(CON, "[[", "op"))
    rhs <- unlist(lapply(CON, "[[", "rhs"))
    LIST$lhs <- c(LIST$lhs, lhs)
    LIST$op <- c(LIST$op, op)
    LIST$rhs <- c(LIST$rhs, rhs)
  } else {
    lav_msg_stop(gettext(
      "no equality constraints found in constraints argument"))
  }

  # theta = free parameters only
  theta <- lav_model_get_parameters(lavmodel)

  # build constraint function
  ceq.function <- lav_partable_constraints_ceq(
    partable = partable,
    con = LIST, debug = FALSE
  )
  # compute jacobian restrictions
  JAC <- try(lav_func_jacobian_complex(func = ceq.function, x = theta),
    silent = TRUE
  )
  if (inherits(JAC, "try-error")) { # eg. pnorm()
    JAC <- lav_func_jacobian_simple(func = ceq.function, x = theta)
  }

  if (verbose) {
    cat("Restriction matrix (jacobian):\n")
    print(JAC)
    cat("\n")
  }

  # linear restriction
  theta.r <- ceq.function(theta)

  if (verbose) {
    cat("Restricted theta values:\n")
    print(theta.r)
    cat("\n")
  }

  # get VCOV
  # VCOV <- vcov(object, labels = FALSE)
  # avoid S4 dispatch
  VCOV <- lav_object_inspect_vcov(object,
    standardized = FALSE,
    free.only = TRUE,
    add.labels = FALSE,
    add.class = FALSE,
    remove.duplicated = FALSE
  )

  # restricted vcov
  VCOV.r <- JAC %*% VCOV %*% t(JAC)

  # fixme: what if VCOV.r is singular?
  # Wald test statistic
  Wald <- as.numeric(t(theta.r) %*% solve(VCOV.r) %*% theta.r)

  # df
  Wald.df <- nrow(JAC)

  # p-value based on chisq
  Wald.pvalue <- 1 - pchisq(Wald, df = Wald.df)

  # prepare output
  out <- list(
    stat = Wald, df = Wald.df, p.value = Wald.pvalue,
    se = lavoptions$se
  )

  out
}
