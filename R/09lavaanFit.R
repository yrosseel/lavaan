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

    # did we compute standard errors?
    se <- numeric( length(est) )
    if(!is.null(VCOV)) { 
        x.se <- sqrt( diag(VCOV) )
        GLIST <- x2GLIST(model, x=x.se, type="free")
        se <- getModelParameters(model, GLIST=GLIST, type="user", 
                                 extra=FALSE) # no def/cin/ceq entries!
        # fixed parameters -> se = 0.0
        se[ which(user$free.uncon == 0L) ] <- 0.0

        # defined parameters: 
        def.idx <- which(user$op == ":=")
        if(length(def.idx) > 0L) {
            if(!is.null(attr(VCOV, "BOOT"))) {
                BOOT <- attr(VCOV, "BOOT")
                BOOT.def <- apply(BOOT, 1, model@def.function)
                if(length(def.idx) == 1L) {
                    BOOT.def <- as.matrix(BOOT.def)
                } else {
                    BOOT.def <- t(BOOT.def)
                }
                def.cov <- cov(BOOT.def )
            } else {
                # regular delta method
                nvar <- length(x)
                JAC <- jacobian(func = model@def.function, x = x, 
                                method = "Richardson")
                def.cov <- JAC %*% VCOV %*% t(JAC)
            }
            se[def.idx] <- sqrt(diag(def.cov))
        }
    }

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
