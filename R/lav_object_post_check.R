# check if a fitted model is admissible
lav_object_post_check <- function(object, verbose = FALSE) {

    stopifnot(inherits(object, "lavaan"))
    lavpartable    <- object@ParTable
    lavmodel       <- object@Model
    lavdata        <- object@Data
    lavfit         <- object@Fit

    result.ok <- TRUE

    # 1. check for Heywood cases, negative (residual) variances, ...
    var.idx <- which(lavpartable$op == "~~" &
                     #!lavpartable$lhs %in% unlist(lavpta$vnames$ov.ord) &
                     lavpartable$lhs == lavpartable$rhs)
    if(length(var.idx) > 0L && any(lavfit@est[var.idx] < 0.0)) {
        result.ok <- FALSE
        warning("lavaan WARNING: some estimated variances are negative")
    }

    # 2. is cov.lv (PSI) positive definite?
    if(length(lavNames(lavpartable, type="lv.regular")) > 0L) {
        ETA <- lavTech(object, "cov.lv")
        for(g in 1:lavdata@ngroups) {
            txt.group <- ifelse(lavdata@ngroups > 1L,
                                paste("in group", g, ".", sep=""), "")
            eigvals <- eigen(ETA[[g]], symmetric=TRUE, only.values=TRUE)$values
            if(any(eigvals < -1 * .Machine$double.eps^(3/4))) {
                warning("lavaan WARNING: covariance matrix of latent variables is not positive definite;", txt.group, " use inspect(fit,\"cov.lv\") to investigate.")
                result.ok <- FALSE
            }
        }
    }

    # 3. is THETA positive definite (but only for numeric variables)
    THETA <- lavTech(object, "theta")
    for(g in 1:lavdata@ngroups) {
        num.idx <- lavmodel@num.idx[[g]]
        if(length(num.idx) > 0L) {
            txt.group <- ifelse(lavdata@ngroups > 1L,
                                paste("in group", g, ".", sep=""), "")
            eigvals <- eigen(THETA[[g]][num.idx, num.idx, drop=FALSE],
                             symmetric = TRUE,
                             only.values = TRUE)$values
            if(any(eigvals < -1 * .Machine$double.eps^(3/4))) {
                warning("lavaan WARNING: observed variable error term matrix (theta) is not positive definite;", txt.group, " use inspect(fit,\"theta\") to investigate.")
                result.ok <- FALSE
            }
        }
    }

    result.ok
}

