# StartingValues.R
#
# YR 30/11/2010: initial version
# YR 08/06/2011: add fabin3 start values for factor loadings

# fill in the 'ustart' column in a User data.frame with reasonable
# starting values, using the sample data

StartingValues <- function(start.method = "default",
                           user         = NULL, 
                           sample       = NULL,
                           model.type   = "sem",
                           mimic        = "lavaan",
                           debug        = FALSE) {

    # check arguments
    stopifnot(is.list(user), class(sample) == "SampleStats")

    # shortcut for 'simple'
    if(start.method == "simple") {
        start <- numeric( length(user$ustart) )
        start[ which(user$op == "=~") ] <- 1.0    
        var.idx <- which(user$op == "~~" & user$lhs == user$rhs)
        start[var.idx] <- 1.0
        user.idx <- which(!is.na(user$ustart))
        start[user.idx] <- user$ustart[user.idx]
        return(start)
    }

    # check start.method
    if(mimic == "lavaan") {
        start.initial <- "lavaan"
    } else if(mimic == "Mplus") {
        start.initial <- "mplus"
    } else {
        # FIXME: use LISREL/EQS/AMOS/.... schems
        start.initial <- "lavaan"
    }
    start.user    <- NULL
    if(is.character(start.method)) {
        start.method <- tolower(start.method)
        if(start.method == "default") {
            # nothing to do
        } else if(start.method %in% c("simple", "lavaan", "mplus")) { 
            start.initial <- start.method
        } else {
            stop("lavaan ERROR: unknown value for start argument")
        }
    } else if(is.list(start.method)) {
        start.user <- start.method
    } else if(class(start.method) == "lavaan") {
        start.user <- parameterEstimates(start.method)
    }
    # check model list elements, if provided
    if(!is.null(start.user)) {
        if(is.null(start.user$lhs) ||
           is.null(start.user$op)  ||
           is.null(start.user$rhs)) {
            stop("lavaan ERROR: problem with start argument: model list does not contain all elements: lhs/op/rhs")
        }
        if(!is.null(start.user$est)) {
            # excellent, we got an est column; nothing to do
        } else if(!is.null(start.user$start)) {
            # no est column, but we use the start column
            start.user$est <- start.user$start
        } else if(!is.null(start.user$ustart)) {
            # no ideal, but better than nothing
            start.user$est <- start.user$ustart
        } else {
            stop("lavaan ERROR: problem with start argument: could not find est/start column in model list")
        }
    }   

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

    if(start.initial %in% c("lavaan", "mplus") && 
       model.type %in% c("sem", "cfa")) {
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
        sample.var.idx <- match(user$lhs[ov.var.idx], ov.names)
        if(start.initial == "mplus") {
            start[ov.var.idx] <- (1.0 - 0.50)*sample@var[[1L]][sample.var.idx]
        } else {
            start[ov.var.idx] <- (1.0 - 0.50)*sample@var[[g]][sample.var.idx]
        }

        # 2. intercepts
        ov.int.idx <- which(user$group == g         &
                            user$op == "~1"         & 
                            user$lhs %in% ov.names)
        sample.var.idx <- match(user$lhs[ov.int.idx], ov.names)
        if(start.initial == "mplus") {
            start[ov.int.idx] <- sample@mean[[1L]][sample.var.idx]
        } else {
            start[ov.int.idx] <- sample@mean[[g]][sample.var.idx]
        }

        # 3. exogenous `fixed.x' covariates
        if(length(ov.names.x) > 0) {
            exo.idx <- which(user$group == g          &
                             user$op == "~~"          & 
                             user$lhs %in% ov.names.x &
                             user$rhs %in% ov.names.x)
            row.idx <- match(user$lhs[exo.idx], ov.names)
            col.idx <- match(user$rhs[exo.idx], ov.names)
            start[exo.idx] <- sample@cov[[g]][ cbind(row.idx, col.idx) ]
        }
    }

    # growth models:
    # - compute starting values for mean latent variables
    # - compute starting values for variance latent variables
    if(start.initial %in% c("lavaan", "mplus") && 
       model.type == "growth") {
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

    # override if a user list with starting values is provided 
    # we only look at the 'est' column for now
    if(!is.null(start.user)) {
        # FIXME: avoid for loop!!!
        for(i in 1:length(user$lhs)) {
            # find corresponding parameters
            lhs <- user$lhs[i]; op <- user$op[i]; rhs <- user$rhs[i]
            start.user.idx <- which(start.user$lhs == lhs &
                                    start.user$op  ==  op &
                                    start.user$rhs == rhs)
            if(length(start.user.idx) == 1L && 
               is.finite(start.user$est[start.user.idx])) {
                start[i] <- start.user$est[start.user.idx]
            }
        }
    }
  
    # override if the model syntax contains explicit starting values
    user.idx <- which(!is.na(user$ustart))
    start[user.idx] <- user$ustart[user.idx]

    

    if(debug) {
        cat("lavaan DEBUG: lavaanStart\n")
        print( start )
    }

    start
}
