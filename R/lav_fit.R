lav_model_fit <- function(lavpartable = NULL,
                          lavmodel    = NULL, 
                          x           = NULL, 
                          VCOV        = NULL, 
                          TEST        = NULL) {

    stopifnot(is.list(lavpartable), class(lavmodel) == "Model")

    # extract information from 'x'
    iterations = attr(x, "iterations")
    converged  = attr(x, "converged")
    fx         = attr(x, "fx")
    fx.group   = attr(fx, "fx.group")
    if(!is.null(attr(fx, "logl.group"))) {
        logl.group = attr(fx, "logl.group")
        logl       = sum(logl.group)
    } else {
        logl.group = as.numeric(NA)
        logl = as.numeric(NA)
    }
    #print(fx.group)
    control    = attr(x, "control")
    attributes(fx) <- NULL
    x.copy <- x # we are going to change it (remove attributes)
    attributes(x.copy) <- NULL
    est <- lav_model_get_parameters(lavmodel = lavmodel, type = "user")

    # did we compute standard errors?
    se <- numeric( length(est) )
    if(!is.null(VCOV)) { 
        x.var <- diag(VCOV)
        # check for negative values (what to do: NA or 0.0?)
        x.var[x.var < 0] <- as.numeric(NA)
        x.se <- sqrt( x.var )
        GLIST <- lav_model_x2GLIST(lavmodel = lavmodel, x = x.se, type = "free")
        se <- lav_model_get_parameters(lavmodel = lavmodel, GLIST = GLIST, 
                                       type = "user", extra = FALSE) # no def/cin/ceq entries!
        # fixed parameters -> se = 0.0
        se[ which(lavpartable$free == 0L) ] <- 0.0

        # defined parameters: 
        def.idx <- which(lavpartable$op == ":=")
        if(length(def.idx) > 0L) {
            if(!is.null(attr(VCOV, "BOOT.COEF"))) {
                BOOT <- attr(VCOV, "BOOT.COEF")
                BOOT.def <- apply(BOOT, 1L, lavmodel@def.function)
                if(length(def.idx) == 1L) {
                    BOOT.def <- as.matrix(BOOT.def)
                } else {
                    BOOT.def <- t(BOOT.def)
                }
                def.cov <- cov(BOOT.def )
            } else {
                # regular delta method
                JAC <- try(lav_func_jacobian_complex(func = lavmodel@def.function, x = x),
                           silent=TRUE)
                if(inherits(JAC, "try-error")) { # eg. pnorm()
                    JAC <- lav_func_jacobian_simple(func = lavmodel@def.function, x = x)
                }
                def.cov <- JAC %*% VCOV %*% t(JAC)
            }
            # check for negative se's
            diag.def.cov <- diag(def.cov)
            diag.def.cov[ diag.def.cov < 0 ] <- as.numeric(NA)
            se[def.idx] <- sqrt(diag.def.cov)
        }
    }

    # did we compute test statistics
    if(is.null(TEST)) {
        test <- list()
    } else {
        test <- TEST
    }

    # for convenience: compute lavmodel-implied Sigma and Mu
    Sigma.hat <- computeSigmaHat(lavmodel = lavmodel)
       Mu.hat <-    computeMuHat(lavmodel = lavmodel)
    if(lavmodel@categorical) {
        TH <- computeTH(lavmodel = lavmodel)
    } else {
        TH <- vector("list", length = lavmodel@ngroups)
    }


    # if bootstrapped parameters, add attr to 'est'
    if(!is.null(attr(VCOV, "BOOT.COEF"))) {
        attr(est, "BOOT.COEF") <- attr(VCOV, "BOOT.COEF")
    }

    # partrace?
    if(!is.null(attr(x, "partrace"))) {
        PARTRACE <- attr(x, "partrace")
    } else {
        PARTRACE <- matrix(0, 0L, 0L)
    }

    new("Fit",
        npar       = max(lavpartable$free),
        x          = x.copy,
        partrace   = PARTRACE,
        start      = lavpartable$start, # needed?
        est        = est,
        se         = se,
        fx         = fx,
        fx.group   = fx.group,
        logl       = logl,
        logl.group = logl.group,
        iterations = iterations,
        converged  = converged,
        control    = control,
        Sigma.hat  = Sigma.hat,
        Mu.hat     = Mu.hat,
        TH         = TH,
        test       = test
       )
}
