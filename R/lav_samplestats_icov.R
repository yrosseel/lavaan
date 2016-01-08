lav_samplestats_icov <- function(COV = NULL, ridge = 0.0,
                                 ngroups = 1L, g = 1L, warn = TRUE) {

    tmp <- try(inv.chol(COV, logdet = TRUE), silent = TRUE)

    # what if this fails...
    if(inherits(tmp, "try-error")) {
        if(warn) {
            if(ngroups > 1) {
                warning("lavaan WARNING sample covariance can not be inverted in group: ", g)
            } else {
                warning("lavaan WARNING: sample covariance can not be inverted")
            }
        }

        # ok, try ridging for exogenous x only
        ## FIXME -- only x (but all for now)
        ridge.eps <- ridge
        diag(COV) <- diag(COV) + ridge.eps

        # try again
        tmp <- try(inv.chol(COV, logdet = TRUE), silent = TRUE)
        if(inherits(tmp, "try-error")) {
            # emergency values
            icov <- MASS::ginv(COV)
            cov.log.det <- log(.Machine$double.eps)
        } else {
            cov.log.det <- attr(tmp, "logdet")
            attr(tmp, "logdet") <- NULL
            icov        <- tmp
        }
    } else {
        cov.log.det <- attr(tmp, "logdet")
        attr(tmp, "logdet") <- NULL
        icov        <- tmp
    }

    list(icov = icov, cov.log.det = cov.log.det)
}
