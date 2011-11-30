# public version
lavaanNames <- function(object, type="ov", group=NULL) {

    if(class(object) == "lavaan") {
         user <- object@User
    } else if(class(object) == "list" ||
              class(object) == "data.frame") {
        user <- object
    }

    stopifnot(is.list(user), 
              type %in% c("lv",   "ov",
                          "lv.x", "ov.x",
                          "lv.y", "ov.y",
                          "ov.nox"))
    vnames(user, type=type, group=group)
}

# internal version
vnames <- function(user, type=NULL, group=NULL) {

    stopifnot(is.list(user), !missing(type),
              type %in% c("lv",   "ov", 
                          "lv.x", "ov.x",
                          "lv.y", "ov.y",
                          "ov.nox"))

    # select single group only?
    if(!is.null(group)) {
        group <- as.integer(group)
        idx <- which(user$group %in% group)
        user$lhs <- user$lhs[idx]
        user$op  <- user$op[idx]
        user$rhs <- user$rhs[idx]
    }

    # regular latent variables: lhs =~
    if(type == "lv") {
        out <- unique( user$lhs[ user$op == "=~" ] )
    } else 

    # observed variables 
    if(type == "ov") {
        lv.names <- unique( user$lhs[ user$op == "=~" ] )
        
        # order is important!
        # 1. indicators, which are not latent variables themselves
        v.ind <- unique( user$rhs[ user$op == "=~" ] )
        ov.ind <- v.ind[ !v.ind %in% lv.names ]

        # 2. dependent ov's
        eqs.y <- unique( user$lhs[ user$op == "~" ] )
        ov.y <- eqs.y[ !eqs.y %in% c(lv.names, ov.ind) ]

        # 2. independent ov's
        eqs.x <- unique( user$rhs[ user$op == "~" ] )
        ov.x <- eqs.x[ !eqs.x %in% c(lv.names, ov.ind, ov.y) ]
        
        out <- c(ov.ind, ov.y, ov.x)
    } else

    # exogenous `x' covariates
    if(type == "ov.x") {
        lv.names <- unique( user$lhs[ user$op == "=~" ] )
        
        # 1. indicators, which are not latent variables themselves
        v.ind <- unique( user$rhs[ user$op == "=~" ] )
        ov.ind <- v.ind[ !v.ind %in% lv.names ]

        # 2. dependent ov's
        eqs.y <- unique( user$lhs[ user$op == "~" ] )
        ov.y <- eqs.y[ !eqs.y %in% c(lv.names, ov.ind) ]

        # 2. independent ov's
        eqs.x <- unique( user$rhs[ user$op == "~" ] )
        ov.x <- eqs.x[ !eqs.x %in% c(lv.names, ov.ind, ov.y) ]

        # correction: is any of these ov.names.x mentioned as a variance,
        #             covariance, or intercept? If so, omit from ov.x!
        # FIXME: better use fixed.x=FALSE flag
        #if(is.null(user$user)) { # FLAT!
        #    vars <- c(user$lhs[user$op == "~1"],
        #              user$lhs[user$op == "~~"],
        #              user$rhs[user$op == "~~"])
        #} else {
        #    vars <- c(user$lhs[user$op == "~1" & user$user == 1],
        #              user$lhs[user$op == "~~" & user$user == 1],
        #              user$rhs[user$op == "~~" & user$user == 1])
        #}
        #idx.no.x <- which(ov.x %in% vars)
        #if(length(idx.no.x)) ov.x <- ov.x[-idx.no.x]
        
        out <- ov.x
    } else

    # ov's withouth ov.x
    if(type == "ov.nox") {
        out <- vnames(user, "ov", group=group)
        ov.names.x <- vnames(user, "ov.x", group=group)
        idx <- which(out %in% ov.names.x)
        if(length(idx)) out <- out[-idx]
    } else


    # exogenous lv's
    if(type == "lv.x") {
        lv.names <- unique( user$lhs[ user$op == "=~" ] )
        v.ind    <- unique( user$rhs[ user$op == "=~" ] )
        eqs.y    <- unique( user$lhs[ user$op == "~"  ] )

        tmp <- lv.names[ !lv.names %in% c(v.ind, eqs.y) ]
        # make sure order is the same as lv.names
        out <- lv.names[ match(tmp, lv.names) ]
    } else
 
    # dependent ov (but not also indicator or x)
    if(type == "ov.y") {
        ov.names <- vnames(user, "ov", group=group)
        lv.names <- unique( user$lhs[ user$op == "=~" ] )
        v.ind    <- unique( user$rhs[ user$op == "=~" ] )
        eqs.x <- unique( user$rhs[ user$op == "~" ] )
        eqs.y    <- unique( user$lhs[ user$op == "~"  ] )

        tmp <- eqs.y[ !eqs.y %in% c(v.ind, eqs.x, lv.names) ]
        # make sure order is the same as ov.names
        out <- ov.names[ match(tmp, ov.names) ]
    } else

    # dependent lv (but not also indicator or x)
    if(type == "lv.y") {
        lv.names <- unique( user$lhs[ user$op == "=~" ] )
        v.ind    <- unique( user$rhs[ user$op == "=~" ] )
        eqs.x <- unique( user$rhs[ user$op == "~" ] )
        eqs.y    <- unique( user$lhs[ user$op == "~"  ] )

        tmp <- eqs.y[ !eqs.y %in% c(v.ind, eqs.x) &
                       eqs.y %in% lv.names ]
        # make sure order is the same as lv.names
        out <- lv.names[ match(tmp, lv.names) ]
    }

    out
}

getNDAT <- function(user, group=NULL) {

    ov.names <- vnames(user, "ov", group=group)
    nvar <- length(ov.names)

    ngroups <- max(user$group)
    meanstructure <- any(user$op == "~1")
    fixed.x <- any(user$fixed.x > 0 & user$free == 0)

    pstar <- nvar*(nvar+1)/2; if(meanstructure) pstar <- pstar + nvar
    ndat  <- ngroups*pstar

    # correction for fixed.x?
    if(fixed.x) {
        ov.names.x <- vnames(user, "ov.x", group=group)
        nvar.x <- length(ov.names.x)
        pstar.x <- nvar.x * (nvar.x + 1) / 2
        if(meanstructure) pstar.x <- pstar.x + nvar.x
        ndat <- ndat - (ngroups * pstar.x)
    }

    ndat
}

getNPAR <- function(user) {
    npar <- max(user$free)
    npar
}

getDF <- function(user, group=NULL) {

    ov.names <- vnames(user, "ov", group=group)
    nvar <- length(ov.names)
    npar <- max(user$free)

    ngroups <- max(user$group)
    meanstructure <- any(user$op == "~1")
    fixed.x <- any(user$fixed.x > 0 & user$free == 0)

    pstar <- nvar*(nvar+1)/2; if(meanstructure) pstar <- pstar + nvar
    ndat  <- ngroups*pstar

    # correction for fixed.x?
    if(fixed.x) {
        ov.names.x <- vnames(user, "ov.x", group=group)
        nvar.x <- length(ov.names.x)
        pstar.x <- nvar.x * (nvar.x + 1) / 2
        if(meanstructure) pstar.x <- pstar.x + nvar.x
        ndat <- ndat - (ngroups * pstar.x)
    }

    # degrees of freedom
    df <- ndat - npar

    as.integer(df)
}

getParameterLabels <- function(user, type="user") {
    # default labels
    label <- paste(user$lhs, user$op, user$rhs, sep="")
    
    # handle multiple groups
    ngroups <- max(user$group)
    if(ngroups > 1L) {
        for(g in 2:ngroups) {
            label[user$group == g] <- 
                paste(label[user$group == 1], ".g", g, sep="")
        }
    }

    # user-specified labels
    user.idx <- which(nchar(user$label) > 0L)
    label[user.idx] <- user$label[user.idx]

    # which labels do we need?
    if(type == "user") {
        idx <- 1:length(label)
    } else if(type == "free") {
        idx <- which(user$free > 0L & !duplicated(user$free))
    } else if(type == "unco") {
        idx <- which(user$free.uncon > 0L & !duplicated(user$free.uncon))
    } else {
        stop("argument `type' must be one of free, unco, or user")
    }
   
    label[idx]
}


getUserListFull <- function(user=NULL, group=NULL) {

    # meanstructure + number of groups
    meanstructure <- any(user$op == "~1")
    ngroups <- max(user$group)

    # extract `names' of various types of variables:
    lv.names     <- vnames(user, type="lv",  group=group)   # latent variables
    ov.names     <- vnames(user, type="ov",  group=group)   # observed variables
    ov.names.x   <- vnames(user, type="ov.x",group=group)   # exogenous x covariates
    ov.names.nox <- vnames(user, type="ov.nox",group=group) # ov's without exo's
    lv.names.x   <- vnames(user, type="lv.x",group=group)   # exogenous lv
    ov.names.y   <- vnames(user, type="ov.y",group=group)   # dependent ov
    lv.names.y   <- vnames(user, type="lv.y",group=group)   # dependent lv
    lvov.names.y <- c(ov.names.y, lv.names.y)


    # 1 "=~"
    l.lhs <- r.rhs <- op <- character(0)
    l.lhs <- rep(lv.names, each=length(ov.names.nox))
    l.rhs <- rep(ov.names.nox, times=length(lv.names))
    l.op  <- rep("=~", length(l.lhs))

    # 2a. "~~" ov ## FIXME: ov.names.nox or ov.names??
    ov.lhs <- ov.rhs <- ov.op <- character(0)
    nx <- length(ov.names)
    idx <- lower.tri(matrix(0, nx, nx), diag=TRUE)
    ov.lhs <- rep(ov.names,  each=nx)[idx] # fill upper.tri
    ov.rhs <- rep(ov.names, times=nx)[idx]
    ov.op  <- rep("~~", length(ov.lhs))

    # 2b. "~~" lv
    lv.lhs <- lv.rhs <- lv.op <- character(0)
    nx <- length(lv.names)
    idx <- lower.tri(matrix(0, nx, nx), diag=TRUE)
    lv.lhs <- rep(lv.names,  each=nx)[idx] # fill upper.tri
    lv.rhs <- rep(lv.names, times=nx)[idx]
    lv.op  <- rep("~~", length(lv.lhs))

    # 3 regressions?
    r.lhs <- r.rhs <- r.op <- character(0)
    if(any(user$op == "~")) {
        # FIXME: better choices?
        eqs.names <- unique( c(user$lhs[ user$op == "~" ],
                               user$rhs[ user$op == "~" ]) )
        r.lhs <- rep(eqs.names, each=length(eqs.names))
        r.rhs <- rep(eqs.names, times=length(eqs.names))
        # remove self-arrows
        idx <- which(r.lhs == r.rhs)
        r.lhs <- r.lhs[-idx]
        r.rhs <- r.rhs[-idx]
        r.op <- rep("~", length(r.rhs))
    }

    # 4. intercetps
    int.lhs <- int.rhs <- int.op <- character(0)
    if(meanstructure) {
        int.lhs <- c(ov.names, lv.names)
        int.rhs <- rep("",   length(int.lhs))
        int.op  <- rep("~1", length(int.lhs))
    }

    # combine
    lhs <- c(l.lhs, ov.lhs, lv.lhs, r.lhs, int.lhs)
    rhs <- c(l.rhs, ov.rhs, lv.rhs, r.rhs, int.rhs)
     op <- c(l.op,  ov.op,  lv.op,  r.op,  int.op)


    # multiple groups!
    group <- 1L
    if(ngroups > 1) {
        group   <- rep(1:ngroups, each=length(lhs))
        lhs     <- rep(lhs,     times=ngroups)
        op      <- rep(op,      times=ngroups)
        rhs     <- rep(rhs,     times=ngroups)
    }

    LIST <- data.frame(lhs=lhs, op=op, rhs=rhs, group=group,
                       stringsAsFactors=FALSE)
        
}

