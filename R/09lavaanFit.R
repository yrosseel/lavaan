Fit <- function(user=NULL, start, model, x=NULL, VCOV=NULL, TEST=NULL) {

    stopifnot(is.list(user), length(user$lhs) == length(start),
              class(model) == "Model")

    # extract information from 'x'
    iterations = attr(x, "iterations")
    converged  = attr(x, "converged")
    fx         = attr(x, "fx")
    fx.group   = attr(fx, "fx.group")
    attributes(fx) <- NULL
    attributes(x) <- NULL
    est <- getModelParameters(model, type="user")

    # impute computed values for 'variable definitions'
    def.idx <- which(user$op == ":=")
    if(length(def.idx) > 0L) {
        def.est <- model@def.function(x)
        est[def.idx] <- def.est
    }

    # set est values for "==" and "<>" elements to NA
    con.idx <- which(user$op == "==" | 
                     user$op == "<"  | 
                     user$op == ">")
    est[con.idx] <- NA

    # did we compute standard errors?
    se <- numeric( length(est) )
    if(!is.null(VCOV)) { 
        x.se <- sqrt( diag(VCOV) )
        GLIST <- x2GLIST(model, x=x.se, type="free")
        se <- getModelParameters(model, GLIST=GLIST, type="user")
        # fixed parameters -> se = 0.0
        se[ which(user$free.uncon == 0L) ] <- 0.0
        #se[ user$free & !duplicated(user$free) ] <- x.se
    }
    se[con.idx] <- NA

    # did we compute test statistics
    if(is.null(TEST)) {
        test <- list()
    } else {
        test <- TEST
    }

    # for convenience: compute model-implied Sigma and Mu
    Sigma.hat <- computeSigmaHat(model)
       Mu.hat <-    computeMuHat(model)

    # if bootstrapped parameters, add attr to 'est'
    if(!is.null(attr(VCOV, "BOOT"))) {
        attr(est, "BOOT") <- attr(VCOV, "BOOT")
    }

    new("Fit",
        npar       = max(user$free),
        x          = x,
        start      = start,
        est        = est,
        se         = se,
        fx         = fx,
        fx.group   = fx.group,
        iterations = iterations,
        converged  = converged,
        Sigma.hat  = Sigma.hat,
        Mu.hat     = Mu.hat,
        test       = test
       )
}
