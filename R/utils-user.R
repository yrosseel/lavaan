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
        group.idx <- which(user$group %in% group)
        if(is.data.frame(user)) {
            user <- user[group.idx,]
        } else {
            user.old <- user; user <- list()
            user$lhs   <- user.old$lhs[group.idx]
            user$op    <- user.old$op[group.idx]
            user$rhs   <- user.old$rhs[group.idx]
            user$group <- user.old$group[group.idx]
            user$exo   <- user.old$exo[group.idx]
        }
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

        # 3. independent ov's
        eqs.x <- unique( user$rhs[ user$op == "~" ] )
        ov.x <- eqs.x[ !eqs.x %in% c(lv.names, ov.ind, ov.y) ]

        out <- c(ov.ind, ov.y, ov.x)

        # 4. orphaned covariances
        ov.cov <- c(user$lhs[ user$op == "~~" & !user$lhs %in% lv.names ], 
                    user$rhs[ user$op == "~~" & !user$rhs %in% lv.names ])
        ov.int <- user$lhs[ user$op == "~1" & !user$lhs %in% lv.names ]
        extra <- unique(c(ov.cov, ov.int))
        extra.idx <- which(!extra %in% out)
        out <- c(out, extra[extra.idx])
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

        # extra
        if(!is.null(user$exo)) {
            ov.cov <- c(user$lhs[ user$op == "~~" & user$exo == 1L],
                        user$rhs[ user$op == "~~" & user$exo == 1L])
            ov.int <- user$lhs[ user$op == "~1" & user$exo == 1L ]
            extra <- unique(c(ov.cov, ov.int))
            extra.idx <- which(!extra %in% out)
            out <- c(out, extra[extra.idx])
        }
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
    fixed.x <- any(user$exo > 0L & user$free == 0L)

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
    fixed.x <- any(user$exo > 0L & user$free == 0L)

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

getParameterLabels <- function(user, group.equal="", group.partial="", 
                               type="user") {
    # default labels
    label <- paste(user$lhs, user$op, user$rhs, sep="")
    
    # handle multiple groups
    ngroups <- max(user$group)
    if(ngroups > 1L) {
        for(g in 2:ngroups) {
            label[user$group == g] <- 
                paste(label[user$group == g], ".g", g, sep="")
        }
    }
 
    #cat("DEBUG: label start:\n"); print(label); cat("\n")
    #cat("group.equal = ", group.equal, "\n")
    #cat("group.partial = ", group.partial, "\n")

    # use group.equal so that equal sets of parameters get the same label
    if(ngroups > 1L && length(group.equal) > 0L) {

        if("intercepts" %in% group.equal ||
           "residuals"  %in%  group.equal ||
           "residual.covariances" %in%  group.equal) {
            ov.names.nox <- vector("list", length=ngroups)
            for(g in 1:ngroups)
                ov.names.nox[[g]] <- vnames(user, "ov.nox", group=g)
        }
        if("means" %in% group.equal ||
           "lv.variances" %in% group.equal ||
           "lv.covariances" %in% group.equal) {
            lv.names <- vector("list", length=ngroups)
            for(g in 1:ngroups)
                lv.names[[g]] <- vnames(user, "lv", group=g)
        }

        # g1.flag: TRUE if included, FALSE if not
        g1.flag <- logical(length(which(user$group == 1L)))

        # LOADINGS
        if("loadings" %in% group.equal)
            g1.flag[ user$op == "=~" & user$group == 1L  ] <- TRUE
        # INTERCEPTS (OV)
        if("intercepts" %in% group.equal)
            g1.flag[ user$op == "~1"  & user$group == 1L  &
                     user$lhs %in% ov.names.nox[[1L]] ] <- TRUE
        # MEANS (LV)
        if("means" %in% group.equal)
            g1.flag[ user$op == "~1" & user$group == 1L &
                     user$lhs %in% lv.names[[1L]] ] <- TRUE
        # REGRESSIONS
        if("regressions" %in% group.equal)
            g1.flag[ user$op == "~" & user$group == 1L ] <- TRUE
        # RESIDUAL variances (FIXME: OV ONLY!)
        if("residuals" %in% group.equal)
            g1.flag[ user$op == "~~" & user$group == 1L &
                     user$lhs %in% ov.names.nox[[1L]] &
                     user$lhs == user$rhs ] <- TRUE
        # RESIDUAL covariances (FIXME: OV ONLY!)
        if("residual.covariances" %in% group.equal)
            g1.flag[ user$op == "~~" & user$group == 1L &
                     user$lhs %in% ov.names.nox[[1L]] &
                     user$lhs != user$rhs ] <- TRUE
        # LV VARIANCES
        if("lv.variances" %in% group.equal)
            g1.flag[ user$op == "~~" & user$group == 1L &
                     user$lhs %in% lv.names[[1L]] &
                     user$lhs == user$rhs ] <- TRUE
        # LV COVARIANCES
        if("lv.covariances" %in% group.equal)
            g1.flag[ user$op == "~~" & user$group == 1L &
                     user$lhs %in% lv.names[[1L]] &
                     user$lhs != user$rhs ] <- TRUE

        # if group.partial, set corresponding flag to FALSE
        if(length(group.partial) > 0L) {
            g1.flag[ label %in% group.partial &
                     user$group == 1L ] <- FALSE
        }

        # for each (constrained) parameter in 'group 1', find a similar one
        # in the other groups (we assume here that the models need
        # NOT be the same across groups!
        g1.idx <- which(g1.flag)
        for(i in 1:length(g1.idx)) {
            ref.idx <- g1.idx[i]
            idx <- which(user$lhs == user$lhs[ref.idx] &
                         user$op  == user$op[ ref.idx] &
                         user$rhs == user$rhs[ref.idx] &
                         user$group > 1L)
            label[idx] <- label[ref.idx]
        }
    }

    #cat("DEBUG: g1.idx = ", g1.idx, "\n")
    #cat("DEBUG: label after group.equal:\n"); print(label); cat("\n")

    # user-specified labels -- override everything!!
    user.idx <- which(nchar(user$label) > 0L)
    label[user.idx] <- user$label[user.idx]

    #cat("DEBUG: user.idx = ", user.idx, "\n")
    #cat("DEBUG: label after user.idx:\n"); print(label); cat("\n")

    # which labels do we need?
    if(type == "user") {
        idx <- 1:length(label)
    } else if(type == "free") {
        idx <- which(user$free > 0L & !duplicated(user$free))
    } else if(type == "unco") {
        idx <- which(user$unco > 0L & !duplicated(user$unco))
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

getLIST <- function(FLAT=NULL,
                    meanstructure   = FALSE,
                    int.ov.free     = FALSE,
                    int.lv.free     = FALSE,
                    orthogonal      = FALSE,
                    std.lv          = FALSE,
                    fixed.x         = TRUE,
                    auto.fix.first  = FALSE,
                    auto.fix.single = FALSE,
                    auto.var        = FALSE,
                    auto.cov.lv.x   = FALSE,
                    auto.cov.y      = FALSE,
                    group.equal     = NULL,
                    ngroups         = 1L) {

    ### DEFAULT elements: parameters that are typically not specified by
    ###                   users, but should typically be considered, 
    ###                   either free or fixed

    # extract `names' of various types of variables:
    lv.names     <- vnames(FLAT, type="lv")     # latent variables
    ov.names     <- vnames(FLAT, type="ov")     # observed variables
    ov.names.x   <- vnames(FLAT, type="ov.x")   # exogenous x covariates 
    ov.names.nox <- vnames(FLAT, type="ov.nox")
    lv.names.x   <- vnames(FLAT, type="lv.x")   # exogenous lv
    ov.names.y   <- vnames(FLAT, type="ov.y")   # dependent ov
    lv.names.y   <- vnames(FLAT, type="lv.y")   # dependent lv
    #lvov.names.y <- c(ov.names.y, lv.names.y)
    lvov.names.y <- c(lv.names.y, ov.names.y)

    lhs <- rhs <- character(0)

    # 1. default (residual) variances and covariances

    # a) (residual) VARIANCES (all ov's, except exo and regular lv's)
    if(auto.var) {
        lhs <- c(lhs, ov.names.nox, lv.names)
        rhs <- c(rhs, ov.names.nox, lv.names)
    }

    # b) `independent` latent variable COVARIANCES (lv.names.x)
    if(auto.cov.lv.x && length(lv.names.x) > 1L) {
        tmp <- combn(lv.names.x, 2)
        lhs <- c(lhs, tmp[1,]) # to fill upper.tri
        rhs <- c(rhs, tmp[2,])
    }

    # c) `dependent` latent variables COVARIANCES (lv.y.idx + ov.y.lv.idx)
    if(auto.cov.y && length(lvov.names.y) > 1L) {
        tmp <- combn(lvov.names.y, 2L)
        lhs <- c(lhs, tmp[1,]) # to fill upper.tri
        rhs <- c(rhs, tmp[2,])
    }

    # d) exogenous x covariates: VARIANCES + COVARIANCES
    if((nx <- length(ov.names.x)) > 0L) {
        idx <- lower.tri(matrix(0, nx, nx), diag=TRUE)
        lhs <- c(lhs, rep(ov.names.x,  each=nx)[idx]) # fill upper.tri
        rhs <- c(rhs, rep(ov.names.x, times=nx)[idx])
    }
    op <- rep("~~", length(lhs))

    # 2. INTERCEPTS
    if(meanstructure) {
        int.lhs <- c(ov.names, lv.names)
        lhs <- c(lhs, int.lhs)
        rhs <- c(rhs, rep("",   length(int.lhs)))
        op  <- c(op,  rep("~1", length(int.lhs)))
    }

    DEFAULT <- data.frame(lhs=lhs, op=op, rhs=rhs,
                          mod.idx=rep(0L, length(lhs)),
                          stringsAsFactors=FALSE)


    # 4. USER: user-specified elements
    lhs     <- FLAT$lhs
     op     <- FLAT$op
    rhs     <- FLAT$rhs
    mod.idx <- FLAT$mod.idx

    lv.names     <- vnames(FLAT, type="lv")     # latent variables
    ov.names     <- vnames(FLAT, type="ov")     # observed variables

    # check order of covariances: we only fill the upper.tri!
    cov.idx <- which(op == "~~" & lhs != rhs)
    for(i in cov.idx) {
        lv.ov.names <- c(lv.names, ov.names) ### FIXME!!! OK??
        lv.idx <- match(c(lhs[i], rhs[i]), lv.ov.names)
        if(lv.idx[1] > lv.idx[2]) { # swap!
            tmp <- lhs[i]; lhs[i] <- rhs[i]; rhs[i] <- tmp
        }
        if(lhs[i] %in% lv.names && rhs[i] %in% lv.names) {
            lv.idx <- match(c(lhs[i], rhs[i]), lv.names)
            if(lv.idx[1] > lv.idx[2]) { # swap!
                tmp <- lhs[i]; lhs[i] <- rhs[i]; rhs[i] <- tmp
            }
        } else if(lhs[i] %in% ov.names && rhs[i] %in% ov.names) {
            ov.idx <- match(c(lhs[i], rhs[i]), ov.names)
            if(ov.idx[1] > ov.idx[2]) { # swap!
                tmp <- lhs[i]; lhs[i] <- rhs[i]; rhs[i] <- tmp
            }
        } else { # mixed!! # we allow this since 0.4-10
            lv.ov.names <- c(lv.names, ov.names) ### FIXME!!! OK??
            lv.idx <- match(c(lhs[i], rhs[i]), lv.ov.names)
            if(lv.idx[1] > lv.idx[2]) { # swap!
                tmp <- lhs[i]; lhs[i] <- rhs[i]; rhs[i] <- tmp
            }
        }
    }

    USER <- data.frame(lhs=lhs, op=op, rhs=rhs, mod.idx=mod.idx,
                       stringsAsFactors=FALSE)

    # check for duplicated elements in DEFAULT
    # - FIXME: can we not avoid this somehow??
    # - for example, if the user model includes 'x1 ~~ x1'
    #   or 'x1 ~ 1' 
    # - remove them from DEFAULT
    TMP <- rbind(DEFAULT[,1:3], USER[,1:3])
    idx <- which(duplicated(TMP, fromLast=TRUE)) # idx should be in DEFAULT
    if(length(idx)) {
        for(i in idx) {
            flat.idx <- which(USER$lhs   == DEFAULT$lhs[i] &
                              USER$op    == DEFAULT$op[i]  &
                              USER$rhs   == DEFAULT$rhs[i])
            if(length(flat.idx) != 1L) {
                cat("[lavaan DEBUG] idx in TMP: i = ", i, "\n"); print(TMP[i,])
                cat("[lavaan DEBUG] idx in DEFAULT: i = ", i, "\n"); print(DEFAULT[i,])
               cat("[lavaan DEBUG] flat.idx:"); print(flat.idx)
            }
    }
        DEFAULT <- DEFAULT[-idx,]
    }

    # now that we have removed all duplicated elements, we can construct
    # the LIST for a single group
    lhs     <- c(USER$lhs, DEFAULT$lhs)
    op      <- c(USER$op,  DEFAULT$op)
    rhs     <- c(USER$rhs, DEFAULT$rhs)
    user    <- c(rep(1L, length(USER$lhs)),
                 rep(0L, length(DEFAULT$lhs)))
    mod.idx <- c(USER$mod.idx, DEFAULT$mod.idx)
    free    <- rep(1L,  length(lhs))
    ustart  <- rep(as.numeric(NA), length(lhs))
    #label   <- paste(lhs, op, rhs, sep="")
    label   <- rep(character(1), length(lhs))
    exo     <- rep(0L, length(lhs))

    # 1. fix metric of regular latent variables
    if(std.lv) {
        # fix metric by fixing the variance of the latent variable
        lv.var.idx <- which(op == "~~" &
                            lhs %in% lv.names & lhs == rhs)
        ustart[lv.var.idx] <- 1.0
          free[lv.var.idx] <- 0L
    }
    if(auto.fix.first) {
        # fix metric by fixing the loading of the first indicator
        mm.idx <- which(op == "=~")
        first.idx <- mm.idx[which(!duplicated(lhs[mm.idx]))]
        ustart[first.idx] <- 1.0
          free[first.idx] <- 0L
    }

    # 2. fix residual variance of single indicators to zero
    if(auto.var && auto.fix.single) {
        mm.idx <- which(op == "=~")
        T <- table(lhs[mm.idx])
        if(any(T == 1L)) {
            # ok, we have a LV with only a single indicator
            lv.names.single <- names(T)[T == 1L]
            # get corresponding indicator if unique
            lhs.mm <- lhs[mm.idx]; rhs.mm <- rhs[mm.idx]
            single.ind <- rhs.mm[which(lhs.mm %in% lv.names.single &
                                       !(duplicated(rhs.mm) |
                                         duplicated(rhs.mm, fromLast=TRUE)))]
            # is the indicator unique?
            if(length(single.ind)) {
                var.idx <- which(op == "~~" & lhs %in% single.ind
                                            & rhs %in% single.ind
                                            & lhs == rhs
                                            & user == 0L)
                ustart[var.idx] <- 0.0
                  free[var.idx] <- 0L
            }
        }
    }

    # 3. orthogonal=TRUE?
    if(orthogonal) {
        # FIXME: only lv.x.idx for now
        lv.cov.idx <- which(op == "~~" &
                            lhs %in% lv.names &
                            lhs != rhs &
                            user == 0L)
        ustart[lv.cov.idx] <- 0.0
          free[lv.cov.idx] <- 0L
    }

    # 4. intercepts
    if(meanstructure) {
        if(int.ov.free == FALSE) {
            # zero intercepts/means observed variables
                   ov.int.idx <- which(op == "~1" &
                                       lhs %in% ov.names &
                                       user == 0L)
            ustart[ov.int.idx] <- 0.0
              free[ov.int.idx] <- 0L
        }
        if(int.lv.free == FALSE) {
            # zero intercepts/means latent variables
                   lv.int.idx <- which(op == "~1" &
                                       lhs %in% lv.names &
                                       user == 0L)
            ustart[lv.int.idx] <- 0.0
              free[lv.int.idx] <- 0L
        }
    }

    # 5. handle exogenous `fixed.x' covariates
    if(length(ov.names.x) > 0 && fixed.x) {
        # 1. variances/covariances
               exo.idx  <- which(op == "~~" &
                                 rhs %in% ov.names.x &
                                 user == 0L)
        ustart[exo.idx] <- as.numeric(NA) # should be overriden later!
          free[exo.idx] <- 0L
           exo[exo.idx] <- 1L

        # 2. intercepts
               exo.int.idx  <- which(op == "~1" &
                                     lhs %in% ov.names.x &
                                     user == 0L)
        ustart[exo.int.idx] <- as.numeric(NA) # should be overriden later!
          free[exo.int.idx] <- 0L
           exo[exo.int.idx] <- 1L
    }

    # 6. multiple groups?
    group <- rep(1L, length(lhs))
    if(ngroups > 1) {
        group   <- rep(1:ngroups, each=length(lhs))
        user    <- rep(user,    times=ngroups)
        lhs     <- rep(lhs,     times=ngroups)
        op      <- rep(op,      times=ngroups)
        rhs     <- rep(rhs,     times=ngroups)
        free    <- rep(free,    times=ngroups)
        ustart  <- rep(ustart,  times=ngroups)
        mod.idx <- rep(mod.idx, times=ngroups)
        label   <- rep(label,   times=ngroups)
        exo     <- rep(exo,     times=ngroups)

        # specific changes per group
        for(g in 2:ngroups) {
            # label
            # label[group == g] <- paste(label[group == 1], ".g", g, sep="")

            # free/fix intercepts
            if(meanstructure) {
                int.idx  <- which(op == "~1" &
                                  lhs %in% lv.names &
                                  user == 0L &
                                  group == g)
                if(int.lv.free == FALSE && g > 1 &&
                   "intercepts" %in% group.equal &&
                   !("means" %in% group.equal) ) {
                      free[ int.idx ] <- 1L
                    ustart[ int.idx ] <- as.numeric(NA)
                }
            }
        } # g
    } # ngroups

    # construct LIST
    #LIST  <- data.frame(
    LIST   <- list(     id          = 1:length(lhs),
                        lhs         = lhs,
                        op          = op,
                        rhs         = rhs,
                        user        = user,
                        group       = group,
                        mod.idx     = mod.idx,
                        free        = free,
                        ustart      = ustart,
                        exo         = exo,
                        label       = label,
                        eq.id       = rep(0L,  length(lhs)),
                        unco        = rep(0L,  length(lhs))
                   )
    #                   stringsAsFactors=FALSE)

    LIST
}

independenceModel <- function(ov.names=NULL, ov.names.x=NULL, sample.cov=NULL,
                              meanstructure=FALSE, sample.mean=NULL) {

    ngroups <- length(ov.names)
    ov.names.nox <- lapply(as.list(1:ngroups), function(g) 
                    ov.names[[g]][ !ov.names[[g]] %in% ov.names.x[[g]] ])


    lhs <- rhs <- op <- character(0)
    group <- free <- exo <- integer(0)
    ustart <- numeric(0)

    for(g in 1:ngroups) {

        # a) VARIANCES (all ov's, except exo's)
        nvar  <- length(ov.names.nox[[g]])
        lhs   <- c(lhs, ov.names.nox[[g]])
         op   <- c(op, rep("~~", nvar))
        rhs   <- c(rhs, ov.names.nox[[g]])
        group <- c(group, rep(g,  nvar))
        free  <- c(free,  rep(1L, nvar))
        exo   <- c(exo,   rep(0L, nvar))

        # starting values
        if(!is.null(sample.cov)) {
            sample.var.idx <- match(ov.names.nox[[g]], ov.names[[g]])
            ustart <- c(ustart, diag(sample.cov[[g]])[sample.var.idx])
        } else {
            ustart <- c(ustart, rep(as.numeric(NA), nvar))
        }

        # meanstructure?
        if(meanstructure) {
            lhs   <- c(lhs, ov.names.nox[[g]])
             op   <- c(op, rep("~1", nvar))
            rhs   <- c(rhs, rep("", nvar))
            group <- c(group, rep(g,  nvar))
            free  <- c(free,  rep(1L, nvar))
            exo   <- c(exo,   rep(0L, nvar))
            # starting values
            if(!is.null(sample.mean)) {
                sample.int.idx <- match(ov.names.nox[[g]], ov.names[[g]])
                ustart <- c(ustart, sample.mean[[g]][sample.int.idx])
            } else {
                ustart <- c(ustart, rep(as.numeric(NA), nvar))
            }
        }


        # fixed.x exogenous variables?
        if((nx <- length(ov.names.x[[g]])) > 0L) {
            idx <- lower.tri(matrix(0, nx, nx), diag=TRUE)
            nel <- sum(idx)
            lhs    <- c(lhs, rep(ov.names.x[[g]],  each=nx)[idx]) # upper.tri
             op    <- c(op, rep("~~", nel))
            rhs    <- c(rhs, rep(ov.names.x[[g]], times=nx)[idx])
            free   <- c(free,  rep(0L, nel))
            group  <- c(group, rep(g,  nel))
            exo    <- c(exo,   rep(1L, nel))
            ustart <- c(ustart, rep(as.numeric(NA), nel))

            # meanstructure?
            if(meanstructure) {
                lhs    <- c(lhs,    ov.names.x[[g]])
                 op    <- c(op,     rep("~1", nx))
                rhs    <- c(rhs,    rep("", nx))
                group  <- c(group,  rep(g,  nx))
                free   <- c(free,   rep(0L, nx))
                exo    <- c(exo,    rep(1L, nx))
                ustart <- c(ustart, rep(as.numeric(NA), nx))
            }
        }

    } # ngroups

    # free counter
    idx.free <- which(free > 0)
    free[idx.free] <- 1:length(idx.free)

    LIST   <- list(     id          = 1:length(lhs),
                        lhs         = lhs,
                        op          = op,
                        rhs         = rhs,
                        user        = rep(1L,  length(lhs)),
                        group       = group,
                        mod.idx     = rep(0L,  length(lhs)),
                        free        = free,
                        ustart      = ustart,
                        exo         = exo,
                        label       = rep("",  length(lhs)),
                        eq.id       = rep(0L,  length(lhs)),
                        unco        = free
                   )
    LIST
}


