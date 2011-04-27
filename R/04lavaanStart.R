# StartingValues.R
#
# initial version: YR 30/11/2010

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
    lv.names    <- vnames(user, "lv")
    ov.names.x  <- vnames(user, "ov.x")
    
    # 0. everyting is zero
    start <- numeric( length(user$ustart) )

    # 1. =~ factor loadings: 1.0
    start[ which(user$op == "=~") ] <- 1.0

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
