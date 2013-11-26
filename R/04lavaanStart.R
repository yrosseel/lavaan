# StartingValues.R
#
# YR 30/11/2010: initial version
# YR 08/06/2011: add fabin3 start values for factor loadings

# fill in the 'ustart' column in a User data.frame with reasonable
# starting values, using the sample data

StartingValues <- function(start.method = "default",
                           partable     = NULL, 
                           samplestats  = NULL,
                           model.type   = "sem",
                           mimic        = "lavaan",
                           debug        = FALSE) {

    # check arguments
    stopifnot(is.list(partable), class(samplestats) == "lavSampleStats")

    # categorical?
    categorical <- any(partable$op == "|")
    #ord.names <- unique(partable$lhs[ partable$op == "|" ])

    # shortcut for 'simple'
    if(identical(start.method, "simple")) {
        start <- numeric( length(partable$ustart) )
        start[ which(partable$op == "=~") ] <- 1.0
        start[ which(partable$op == "~*~") ] <- 1.0
        ov.names.ord <- vnames(partable, "ov.ord")
        var.idx <- which(partable$op == "~~" & partable$lhs == partable$rhs &
                         !(partable$lhs %in% ov.names.ord))
        start[var.idx] <- 1.0
        user.idx <- which(!is.na(partable$ustart))
        start[user.idx] <- partable$ustart[user.idx]
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
        start.method. <- tolower(start.method)
        if(start.method. == "default") {
            # nothing to do
        } else if(start.method. %in% c("simple", "lavaan", "mplus")) { 
            start.initial <- start.method.
        } else {
            stop("lavaan ERROR: unknown value for start argument")
        }
    } else if(is.list(start.method)) {
        start.user <- start.method
    } else if(inherits(start.method, "lavaan")) {
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


    # global settings
    # 0. everyting is zero
    start <- numeric( length(partable$ustart) )

    # 1. =~ factor loadings: 
    if(categorical) {
        # if std.lv=TRUE, more likely initial Sigma.hat is positive definite
        # 0.8 is too large
        start[ which(partable$op == "=~") ] <- 0.7
    } else {
        start[ which(partable$op == "=~") ] <- 1.0
    }

    # 2. residual lv variances for latent variables
    lv.names    <- vnames(partable, "lv") # all groups
    lv.var.idx <- which(partable$op == "~~"        &
                        partable$lhs %in% lv.names &
                        partable$lhs == partable$rhs)
    start[lv.var.idx] <- 0.05

    # 3. latent response scales (if any)
    delta.idx <- which(partable$op == "~*~")
    start[delta.idx] <- 1.0


    # group-specific settings
    ngroups <- samplestats@ngroups

    for(g in 1:ngroups) {

        # info from user model for this group
        if(categorical) {
            ov.names     <- vnames(partable, "ov.nox", group=g)
            ov.names.num <- vnames(partable, "ov.num", group=g)
        } else {
            ov.names.num <- ov.names <- vnames(partable, "ov", group=g)
        }
        lv.names    <- vnames(partable, "lv",   group=g)
        ov.names.x  <- vnames(partable, "ov.x", group=g)

        # g1) factor loadings
        if(start.initial %in% c("lavaan", "mplus") && 
           model.type %in% c("sem", "cfa") &&
           #!categorical &&
           sum( partable$ustart[ partable$op == "=~" & partable$group == g],
                                   na.rm=TRUE) == length(lv.names) ) {
            # only if all latent variables have a reference item,
            # we use the fabin3 estimator (2sls) of Hagglund (1982)
            # per factor
            # 9 Okt 2013: if only 2 indicators, we use the regression
            # coefficient (y=marker, x=2nd indicator)
            for(f in lv.names) {
                free.idx <- which( partable$lhs == f & partable$op == "=~"
                                                 & partable$group == g
                                                 & partable$free > 0L)
                 
                user.idx <- which( partable$lhs == f & partable$op == "=~" 
                                                 & partable$group == g )
                # no second order
                if(any(partable$rhs[user.idx] %in% lv.names)) next

                # get observed indicators for this latent variable
                ov.idx <- match(partable$rhs[user.idx], ov.names)
                if(length(ov.idx) > 2L && !any(is.na(ov.idx))) {
                    if(samplestats@missing.flag) {
                        COV <- samplestats@missing.h1[[g]]$sigma[ov.idx,ov.idx]
                    } else {
                        COV <- samplestats@cov[[g]][ov.idx,ov.idx]
                    }
                    start[user.idx] <- fabin3.uni(COV)
                } else if(length(free.idx) == 1L && length(ov.idx) == 2L) {
                    REG2 <- ( samplestats@cov[[g]][ov.idx[1],ov.idx[2]] /
                              samplestats@cov[[g]][ov.idx[1],ov.idx[1]] )
                    start[free.idx] <- REG2
                }

                # standardized?
                var.f.idx <- which(partable$lhs == f & partable$op == "~~" &
                                   partable$rhs == f)
                if(partable$free[var.f.idx] == 0 &&
                   partable$ustart[var.f.idx] == 1) {
                   # make sure factor loadings are between -0.7 and 0.7
                    x <- start[user.idx]
                    start[user.idx] <- (x / max(abs(x))) * 0.7
                }
            }
        }

        if(model.type == "unrestricted") {
           # fill in 'covariances' from samplestats
            cov.idx <- which(partable$group == g             &
                             partable$op    == "~~"          &
                             partable$lhs != partable$rhs)
            lhs.idx <- match(partable$lhs[cov.idx], ov.names)
            rhs.idx <- match(partable$rhs[cov.idx], ov.names)
            start[cov.idx] <- samplestats@cov[[g]][ cbind(lhs.idx, rhs.idx) ]
        }

        # 2g) residual ov variances (including exo, to be overriden)
        ov.var.idx <- which(partable$group == g             & 
                            partable$op    == "~~"          & 
                            partable$lhs %in% ov.names.num  & 
                            partable$lhs == partable$rhs)
        sample.var.idx <- match(partable$lhs[ov.var.idx], ov.names.num)
        if(model.type == "unrestricted") {
            start[ov.var.idx] <- diag(samplestats@cov[[g]])[sample.var.idx]
        } else {
            if(start.initial == "mplus") {
                start[ov.var.idx] <- 
                    (1.0 - 0.50)*samplestats@var[[1L]][sample.var.idx]
            } else {
                # start[ov.var.idx] <- 
                #     (1.0 - 0.50)*samplestats@var[[g]][sample.var.idx]
                start[ov.var.idx] <- 
                    (1.0 - 0.50)*diag(samplestats@cov[[g]])[sample.var.idx]
            }
        }

        # 3g) intercepts
        ov.int.idx <- which(partable$group == g         &
                            partable$op == "~1"         & 
                            partable$lhs %in% ov.names)
        sample.int.idx <- match(partable$lhs[ov.int.idx], ov.names)
        if(samplestats@missing.flag) {
            start[ov.int.idx] <- samplestats@missing.h1[[g]]$mu[sample.int.idx]
        } else {
            start[ov.int.idx] <- samplestats@mean[[g]][sample.int.idx]
        }
        
        # thresholds
        th.idx <- which(partable$group == g & partable$op == "|")
        if(length(th.idx) > 0L) {
            th.names.partable <- paste(partable$lhs[th.idx], "|",
                                       partable$rhs[th.idx], sep="")
            th.names.sample   <- 
                samplestats@th.names[[g]][ samplestats@th.idx[[g]] > 0L ]
            # th.names.sample should identical to
           # vnames(partable, "th", group = g)
           th.values <- samplestats@th.nox[[g]][ samplestats@th.idx[[g]] > 0L ]
            start[th.idx] <- th.values[match(th.names.partable,
                                             th.names.sample)]
        }
        

        # 4g) exogenous `fixed.x' covariates
        if(!categorical && length(ov.names.x) > 0) {
            exo.idx <- which(partable$group == g          &
                             partable$op == "~~"          & 
                             partable$lhs %in% ov.names.x &
                             partable$rhs %in% ov.names.x)
            row.idx <- match(partable$lhs[exo.idx], ov.names)
            col.idx <- match(partable$rhs[exo.idx], ov.names)
            if(samplestats@missing.flag) {
                start[exo.idx] <- 
                    samplestats@missing.h1[[g]]$sigma[cbind(row.idx,col.idx)]
            } else {
                start[exo.idx] <- samplestats@cov[[g]][cbind(row.idx,col.idx)]
            }
        }
    }

    # group weights
    group.idx <- which(partable$lhs == "group" &
                       partable$op  == "%")
    if(length(group.idx) > 0L) {
        prop <- rep(1/ngroups, ngroups)
        # use last group as reference
        start[group.idx] <- log(prop/prop[ngroups])
    }

    # growth models:
    # - compute starting values for mean latent variables
    # - compute starting values for variance latent variables
    if(start.initial %in% c("lavaan", "mplus") && 
       model.type == "growth") {
        ### DEBUG ONLY
        #lv.var.idx <- which(partable$op == "~~"                &
        #                partable$lhs %in% lv.names &
        #                partable$lhs == partable$rhs)
        #start[lv.var.idx] <- c(2.369511, 0.7026852)
        
        ### DEBUG ONLY
        #lv.int.idx <- which(partable$op == "~1"         &
        #                    partable$lhs %in% lv.names)
        #start[lv.int.idx] <- c(0.617156788, 1.005192793)
    }

    # override if a user list with starting values is provided 
    # we only look at the 'est' column for now
    if(!is.null(start.user)) {
        # FIXME: avoid for loop!!!
        for(i in 1:length(partable$lhs)) {
            # find corresponding parameters
            lhs <- partable$lhs[i]; op <- partable$op[i]; rhs <- partable$rhs[i]
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
    user.idx <- which(!is.na(partable$ustart))
    start[user.idx] <- partable$ustart[user.idx]

    if(debug) {
        cat("lavaan DEBUG: lavaanStart\n")
        print( start )
    }

    start
}
