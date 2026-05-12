# deprecated: only kept in order to avoid some older packages
lav_model_fit <- function(lavpartable = NULL,
                          lavmodel = NULL,
                          lavimplied = NULL,
                          x = NULL,
                          vcov_1 = NULL,
                          test = NULL) {
  stopifnot(is.list(lavpartable), inherits(lavmodel, "lavModel"))

  # extract information from 'x'
  iterations <- attr(x, "iterations")
  converged <- attr(x, "converged")
  fx <- attr(x, "fx")
  fx_group <- attr(fx, "fx.group")
  if (!is.null(attr(fx, "logl.group"))) {
    logl_group <- attr(fx, "logl.group")
    logl <- sum(logl_group)
  } else {
    logl_group <- as.numeric(NA)
    logl <- as.numeric(NA)
  }
  # print(fx.group)
  control <- attr(x, "control")
  attributes(fx) <- NULL
  x_copy <- x # we are going to change it (remove attributes)
  attributes(x_copy) <- NULL
  est <- lav_model_get_parameters(lavmodel = lavmodel, type = "user")

  # did we compute standard errors?
  if (is.null(lavpartable$se)) {
    if (is.null(vcov_1)) {
      se <- rep(as.numeric(NA), lavmodel@nx.user)
      se[lavpartable$free == 0L] <- 0
    } else {
      se <- lav_model_vcov_se(
        lavmodel = lavmodel,
        lavpartable = lavpartable,
        VCOV = vcov_1,
        BOOT = attr(vcov_1, "BOOT.COEF")
      )
    }
  } else {
    se <- as.numeric(lavpartable$se) # could be logical NA
  }

  # did we compute test statistics
  if (is.null(test)) {
    test_1 <- list()
  } else {
    test_1 <- test
  }

  # for convenience: compute lavmodel-implied Sigma and Mu
  if (is.null(lavimplied) || length(lavimplied) == 0L) {
    implied <- lav_model_implied(lavmodel)
  } else {
    implied <- lavimplied
  }

  # if bootstrapped parameters, add attr to 'est'
  if (!is.null(attr(vcov_1, "BOOT.COEF"))) {
    attr(est, "BOOT.COEF") <- attr(vcov_1, "BOOT.COEF")
  }

  # partrace?
  if (!is.null(attr(x, "partrace"))) {
    partrace <- attr(x, "partrace")
  } else {
    partrace <- matrix(0, 0L, 0L)
  }

  new("Fit",
    npar       = max(lavpartable$free),
    x          = x_copy,
    partrace   = partrace,
    start      = lavpartable$start, # needed? (package stremo!)
    est        = est, # at least 5 packages!!
    se         = se,
    fx         = fx,
    fx.group   = fx_group,
    logl       = logl,
    logl.group = logl_group,
    iterations = iterations,
    converged  = converged,
    control    = control,
    Sigma.hat  = if (lavmodel@conditional.x) implied$res.cov else implied$cov,
    Mu.hat     = if (lavmodel@conditional.x) implied$res.int else implied$mean,
    TH         = if (lavmodel@conditional.x) implied$res.th else implied$th,
    test       = test_1
  )
}
