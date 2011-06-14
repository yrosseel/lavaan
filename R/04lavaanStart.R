# StartingValues.R
#
# YR 30/11/2010: initial version
# YR 08/06/2011: add fabin3 start values for factor loadings

# fill in the 'ustart' column in a User data.frame with reasonable
# starting values, using the sample data

StartingValues <- function(user       = NULL, 
                           sample     = NULL,
                           model.type = "sem",
                           debug      = FALSE) {

    # check arguments
    stopifnot(is.list(user), class(sample) == "Sample")

    # info from sample
    ngroups <- sample@ngroups

    # info from user model
    ov.names    <- vnames(user, "ov")
    lv.names    <- vnames(user, "lv"); nfac <- length(lv.names)
    ov.names.x  <- vnames(user, "ov.x")
    
    # 0. everyting is zero
    start <- numeric( length(user$ustart) )

    # 1. =~ factor loadings: 1.0 
    start[ which(user$op == "=~") ] <- 1.0
    # or slightly better using 2sls
    for(g in 1:ngroups) {
        if(!sum( user$ustart[ user$op == "=~" & user$group == g], 
                na.rm=TRUE) == length(lv.names)) {
            next
        }
        # only if all latent variables have a reference item,
        # we use the fabin3 estimator (2sls) of Hagglund (1982)
        # per factor
        for(f in lv.names) {
            free.idx <- which( user$lhs == f & user$op == "=~"
                                             & user$group == g
                                             & user$free > 0L)
            if(length(free.idx) < 2L) next
            user.idx <- which( user$lhs == f & user$op == "=~" 
                                             & user$group == g )
            # no second order
            if(any(user$rhs[user.idx] %in% lv.names)) next

            # get observed indicators for this latent variable
            ov.idx <- match(user$rhs[user.idx], ov.names)
            if(length(ov.idx) > 2L && !any(is.na(ov.idx))) {
                if(sample@missing.flag[g]) {
                    COV <- sample@missing[[g]]$sigma[ov.idx,ov.idx]
                } else {
                    COV <- sample@cov[[g]][ov.idx,ov.idx]
                }
                start[user.idx] <- fabin3.uni(COV)
            }
        }
    }

    # 2. residual lv variances for latent variables
    lv.var.idx <- which(user$op == "~~"        & 
                        user$lhs %in% lv.names & 
                        user$lhs == user$rhs)
    start[lv.var.idx] <- 0.05

    # group/sample dependent
    for(g in 1:ngroups) {
        # 1. residual ov variances (including exo, to be overriden)
        ov.var.idx <- which(user$group == g         & 
                            user$op    == "~~"      & 
                            user$lhs %in% ov.names  & 
                            user$lhs == user$rhs)
        start[ov.var.idx] <- (1.0 - 0.50) * sample@var[[g]]

        # 2. intercepts
        ov.int.idx <- which(user$group == g         &
                            user$op == "~1"         & 
                            user$lhs %in% ov.names)
        start[ov.int.idx] <- sample@mean[[g]]

        # 3. exogenous `fixed.x' covariates
        if(length(ov.names.x) > 0) {
            exo.idx <- which(user$group == g        &
                             user$op == "~~"        & 
                             user$rhs %in% ov.names.x)
            row.idx <- match(user$lhs[exo.idx], ov.names)
            col.idx <- match(user$rhs[exo.idx], ov.names)
            start[exo.idx] <- sample@cov[[g]][ cbind(row.idx, col.idx) ]
        }
    }

    # growth models:
    # - compute starting values for mean latent variables
    # - compute starting values for variance latent variables
    if(model.type == "growth") {
        ### DEBUG ONLY
        #lv.var.idx <- which(user$op == "~~"                &
        #                user$lhs %in% lv.names &
        #                user$lhs == user$rhs)
        #start[lv.var.idx] <- c(2.369511, 0.7026852)
        
        ### DEBUG ONLY
        #lv.int.idx <- which(user$op == "~1"         &
        #                    user$lhs %in% lv.names)
        #start[lv.int.idx] <- c(0.617156788, 1.005192793)
    }
  
    # override if user-specified starting values are provided
    user.idx <- which(!is.na(user$ustart))
    start[user.idx] <- user$ustart[user.idx]

    if(debug) {
        cat("lavaan DEBUG: lavaanStart\n")
        print( start )
    }

    start
}
