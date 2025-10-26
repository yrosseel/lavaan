lav_samplestats_icov <- function(COV = NULL, ridge = 0.0, x.idx = integer(0L),
                                 ngroups = 1L, g = 1L) {

  cS <- tryCatch(chol(COV), error = function(e) NULL)

  # what if this fails...
  # ridge exogenous part only (if any); this may help for GLS (but not ML)
  if (is.null(cS)) {
    if (length(x.idx) > 0L && ridge > 0) {
      # maybe, we can fix it by gently ridging the exo variances
      ridge.eps <- ridge
      diag(COV)[x.idx] <- diag(COV)[x.idx] + ridge.eps

      # try again
      cS <- tryCatch(chol(COV), error = function(e) NULL)

      if (is.null(cS)) {
        # fatal stop after all
        lav_msg_stop(gettext(
          "sample covariance matrix is not positive-definite"))
      } else {
        cov <- chol2inv(cS)
        d <- diag(cS)
        cov.log.det <- 2 * sum(log(d))
        # give a warning
        if (ngroups > 1) {
          lav_msg_warn(gettextf(
            "sample covariance matrix in group: %s is not
            positive-definite", g))
        } else {
          lav_msg_warn(gettext(
            "sample covariance matrix is not positive-definite"))
        }
      }
    } else {
      # fatal stop
      lav_msg_stop(gettext("sample covariance matrix is not positive-definite"))
    }
  } else {
    icov <- chol2inv(cS)
    d <- diag(cS)
    cov.log.det <- 2 * sum(log(d))
  }

  list(icov = icov, cov.log.det = cov.log.det)
}
