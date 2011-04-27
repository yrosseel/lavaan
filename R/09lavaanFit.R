Fit <- function(user=NULL, start, model, x=NULL, VCOV=NULL, TEST=NULL) {

    stopifnot(is.list(user), length(user$lhs) == length(start),
              class(model) == "Model")

    est <- start
    if(is.null(x)) {
        x <- numeric(0L)
        se <- numeric(0L)
        converged  <- FALSE
        iterations <- 0L
        fx <- as.numeric(NA)
        fx.group <- as.numeric(NA)
    } else {
        iterations = attr(x, "iterations")
        converged  = attr(x, "converged")
        fx         = attr(x, "fx")
        fx.group   = attr(fx, "fx.group")
        attributes(fx) <- NULL
        attributes(x) <- NULL

        est <- getModelParameters(model, type="user")
        se <- numeric( length(est) )
        if(!is.null(VCOV) && nrow(VCOV) == length(x)) { 
            x.se <- sqrt( diag(VCOV) )
            GLIST <- x2GLIST(model, x=x.se, type="free")
            se <- getModelParameters(model, GLIST=GLIST, type="user")
            # fixed parameters -> se = 0.0
            se[ which(user$free.uncon == 0L) ] <- 0.0
            #se[ user$free & !duplicated(user$free) ] <- x.se
        }
    }

    # compute model-implied Sigma and Mu
    Sigma.hat <- computeSigmaHat(model)
       Mu.hat <-    computeMuHat(model)

    # test statistics
    if(is.null(TEST)) {
        test <- list()
    } else {
        test <- TEST
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
