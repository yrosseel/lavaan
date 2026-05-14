lav_samplestats_icov <- function(cov_1 = NULL, ridge = 0.0, x_idx = integer(0L),
                                 ngroups = 1L, g = 1L) {

  c_s <- tryCatch(chol(cov_1), error = function(e) NULL)

  # what if this fails...
  # ridge exogenous part only (if any); this may help for GLS (but not ML)
  if (is.null(c_s)) {
    if (length(x_idx) > 0L && ridge > 0) {
      # maybe, we can fix it by gently ridging the exo variances
      ridge_eps <- ridge
      diag(cov_1)[x_idx] <- diag(cov_1)[x_idx] + ridge_eps

      # try again
      c_s <- tryCatch(chol(cov_1), error = function(e) NULL)

      if (is.null(c_s)) {
        # fatal stop after all
        lav_msg_stop(gettext(
          "sample covariance matrix is not positive-definite"))
      } else {
        icov <- chol2inv(c_s)
        d <- diag(c_s)
        cov_log_det <- 2 * sum(log(d))
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
    icov <- chol2inv(c_s)
    d <- diag(c_s)
    cov_log_det <- 2 * sum(log(d))
  }

  list(icov = icov, cov.log.det = cov_log_det)
}
