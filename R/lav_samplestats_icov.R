lav_samplestats_icov <- function(COV = NULL, ridge = 0.0, x.idx = integer(0L),
                                 ngroups = 1L, g = 1L, warn = TRUE) {

    tmp <- try(inv.chol(COV, logdet = TRUE), silent = TRUE)

    # what if this fails...
    if(inherits(tmp, "try-error")) {

        if(length(x.idx) > 0L) {
            # maybe, we can fix it by gently ridging the exo variances
            ridge.eps <- ridge
            diag(COV)[x.idx] <- diag(COV)[x.idx] + ridge.eps

            # try again
            tmp <- try(inv.chol(COV, logdet = TRUE), silent = TRUE)

            if(inherits(tmp, "try-error")) {
                # fatal stop after all
                stop("lavaan ERROR: sample covariance matrix is not positive-definite")
            } else {
                cov.log.det <- attr(tmp, "logdet")
                attr(tmp, "logdet") <- NULL
                icov <- tmp

                # give a warning
                if(warn) {
                    if(ngroups > 1) {
                        warning("lavaan WARNING sample covariance matrix in group: ",
                                g, " is not positive-definite")
                    } else {
                        warning("lavaan WARNING: sample covariance matrix is not positive-definite")
                    }
                }
            }
        } else {
            # fatal stop
            stop("lavaan ERROR: sample covariance matrix is not positive-definite")
        }
    } else {
        cov.log.det <- attr(tmp, "logdet")
        attr(tmp, "logdet") <- NULL
        icov <- tmp
    }

    list(icov = icov, cov.log.det = cov.log.det)
}
