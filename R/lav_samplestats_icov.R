lav_samplestats_icov <- function(COV = NULL, ridge = 0.0, x.idx = integer(0L),
                                 ngroups = 1L, g = 1L, warn = TRUE) {
  tmp <- try(inv.chol(COV, logdet = TRUE), silent = TRUE)

  # what if this fails...
  # ridge exogenous part only (if any); this may help for GLS (but not ML)
  if (inherits(tmp, "try-error")) {
    if (length(x.idx) > 0L && ridge > 0) {
      # maybe, we can fix it by gently ridging the exo variances
      ridge.eps <- ridge
      diag(COV)[x.idx] <- diag(COV)[x.idx] + ridge.eps

      # try again
      tmp <- try(inv.chol(COV, logdet = TRUE), silent = TRUE)

      if (inherits(tmp, "try-error")) {
        # fatal stop after all
        lav_msg_stop(gettext(
          "sample covariance matrix is not positive-definite"))
      } else {
        cov.log.det <- attr(tmp, "logdet")
        attr(tmp, "logdet") <- NULL
        icov <- tmp

        # give a warning
        if (warn) {
          if (ngroups > 1) {
            lav_msg_warn(gettextf(
              "sample covariance matrix in group: %s is not
              positive-definite", g))
          } else {
            lav_msg_warn(gettext(
              "sample covariance matrix is not positive-definite"))
          }
        }
      }
    } else {
      # fatal stop
      lav_msg_stop(gettext("sample covariance matrix is not positive-definite"))
    }
  } else {
    cov.log.det <- attr(tmp, "logdet")
    attr(tmp, "logdet") <- NULL
    icov <- tmp
  }

  list(icov = icov, cov.log.det = cov.log.det)
}
