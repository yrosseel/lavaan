# lav_partable (old name: utils-user.R)
#
# functions to generate/compute/extract information from the lavaan
# `parameter table'
#
# YR. 29 june 2013 (as lav_partable)

# return 'attributes' of a lavaan partable -- generate a new set if necessary
lav_partable_attributes <- function(partable, pta=NULL) {

    if(is.null(pta)) {
        # attached to partable?
        pta <- attributes(partable)
        if(!is.null(pta$vnames)) {
            return(pta)
        } else {
            pta <- list()
        }
    }

    if(is.null(pta$vnames)) {
        pta$vnames <- lav_partable_vnames(partable, type="all", group="list")
    }

    pta
}


getIDX <- function(partable, type="th", group=NULL) {

    stopifnot(is.list(partable), !missing(type),
              type %in% c("th"))

    if(type == "th") {
        ovn <- lav_partable_vnames(partable, type="ov.nox", group=group)
        ov.num <- lav_partable_vnames(partable, type="ov.num", group=group)
        th  <- lav_partable_vnames(partable, type="th.mean", group=group)
        th[th %in% ov.num] <- "__NUM__"
        th1 <- gsub("\\|t[0-9]*","",th)
        out <- match(th1, ovn)
        out[is.na(out)] <- 0
    }

    out
}

getNDAT <- function(partable, group=NULL) {

    ngroups <- max(partable$group)
    meanstructure <- any(partable$op == "~1")
    fixed.x <- any(partable$exo > 0L & partable$free == 0L)
    categorical <- any(partable$op == "|")
    if(categorical) meanstructure <- TRUE

    if(categorical) {
        ov.names <- lav_partable_vnames(partable, "ov.nox", group=group)
    } else {
        ov.names <- lav_partable_vnames(partable, "ov", group=group)
    }
    nvar <- length(ov.names)

    pstar <- nvar*(nvar+1)/2; if(meanstructure) pstar <- pstar + nvar
    ndat  <- ngroups*pstar

    # correction for fixed.x?
    if(!categorical && fixed.x) {
        ov.names.x <- lav_partable_vnames(partable, "ov.x", group=group)
        nvar.x <- length(ov.names.x)
        pstar.x <- nvar.x * (nvar.x + 1) / 2
        if(meanstructure) pstar.x <- pstar.x + nvar.x
        ndat <- ndat - (ngroups * pstar.x)
    }

    # correction for ordinal data?
    if(categorical) {
        ov.names.x <- lav_partable_vnames(partable, "ov.x", group=group)
        nexo     <- length(ov.names.x)
        ov.ord   <- lav_partable_vnames(partable, "ov.ord", group=group)
        nvar.ord <- length(ov.ord)
        th       <- lav_partable_vnames(partable, "th", group=group)
        nth      <- length(th)
        # no variances
        ndat <- ndat - (ngroups * nvar.ord)
        # no means
        ndat <- ndat - (ngroups * nvar.ord)
        # but additional thresholds
        ndat <- ndat + (ngroups * nth)
        # add slopes
        ndat <- ndat + (ngroups * nvar * nexo)
    }

    ndat
}

getNPAR <- function(partable) {
    npar <- max(partable$free)
    npar
}

getDF <- function(partable, group=NULL) {

    npar <- getNPAR(partable)
    ndat <- getNDAT(partable, group=group)

    # degrees of freedom
    df <- ndat - npar

    as.integer(df)
}

getParameterLabels <- function(partable, group.equal="", group.partial="", 
                               type="user") {
    # default labels
    label <- paste(partable$lhs, partable$op, partable$rhs, sep="")
    
    # handle multiple groups
    ngroups <- max(partable$group)
    if(ngroups > 1L) {
        for(g in 2:ngroups) {
            label[partable$group == g] <- 
                paste(label[partable$group == g], ".g", g, sep="")
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
                ov.names.nox[[g]] <- lav_partable_vnames(partable, "ov.nox", group=g)
        }
        if("thresholds" %in% group.equal) {
            ov.names.ord <- vector("list", length=ngroups)
            for(g in 1:ngroups)
                ov.names.ord[[g]] <- lav_partable_vnames(partable, "ov.ord", group=g)
        }
        if("means" %in% group.equal ||
           "lv.variances" %in% group.equal ||
           "lv.covariances" %in% group.equal) {
            lv.names <- vector("list", length=ngroups)
            for(g in 1:ngroups)
                lv.names[[g]] <- lav_partable_vnames(partable, "lv", group=g)
        }

        # g1.flag: TRUE if included, FALSE if not
        g1.flag <- logical(length(which(partable$group == 1L)))

        # LOADINGS
        if("loadings" %in% group.equal)
            g1.flag[ partable$op == "=~" & partable$group == 1L  ] <- TRUE
        # INTERCEPTS (OV)
        if("intercepts" %in% group.equal)
            g1.flag[ partable$op == "~1"  & partable$group == 1L  &
                     partable$lhs %in% ov.names.nox[[1L]] ] <- TRUE
        # THRESHOLDS (OV-ORD)
        if("thresholds" %in% group.equal)
            g1.flag[ partable$op == "|"  & partable$group == 1L  &
                     partable$lhs %in% ov.names.ord[[1L]] ] <- TRUE
        # MEANS (LV)
        if("means" %in% group.equal)
            g1.flag[ partable$op == "~1" & partable$group == 1L &
                     partable$lhs %in% lv.names[[1L]] ] <- TRUE
        # REGRESSIONS
        if("regressions" %in% group.equal)
            g1.flag[ partable$op == "~" & partable$group == 1L ] <- TRUE
        # RESIDUAL variances (FIXME: OV ONLY!)
        if("residuals" %in% group.equal)
            g1.flag[ partable$op == "~~" & partable$group == 1L &
                     partable$lhs %in% ov.names.nox[[1L]] &
                     partable$lhs == partable$rhs ] <- TRUE
        # RESIDUAL covariances (FIXME: OV ONLY!)
        if("residual.covariances" %in% group.equal)
            g1.flag[ partable$op == "~~" & partable$group == 1L &
                     partable$lhs %in% ov.names.nox[[1L]] &
                     partable$lhs != partable$rhs ] <- TRUE
        # LV VARIANCES
        if("lv.variances" %in% group.equal)
            g1.flag[ partable$op == "~~" & partable$group == 1L &
                     partable$lhs %in% lv.names[[1L]] &
                     partable$lhs == partable$rhs ] <- TRUE
        # LV COVARIANCES
        if("lv.covariances" %in% group.equal)
            g1.flag[ partable$op == "~~" & partable$group == 1L &
                     partable$lhs %in% lv.names[[1L]] &
                     partable$lhs != partable$rhs ] <- TRUE

        # if group.partial, set corresponding flag to FALSE
        if(length(group.partial) > 0L) {
            g1.flag[ label %in% group.partial &
                     partable$group == 1L ] <- FALSE
        }

        # for each (constrained) parameter in 'group 1', find a similar one
        # in the other groups (we assume here that the models need
        # NOT be the same across groups!
        g1.idx <- which(g1.flag)
        for(i in 1:length(g1.idx)) {
            ref.idx <- g1.idx[i]
            idx <- which(partable$lhs == partable$lhs[ref.idx] &
                         partable$op  == partable$op[ ref.idx] &
                         partable$rhs == partable$rhs[ref.idx] &
                         partable$group > 1L)
            label[idx] <- label[ref.idx]
        }
    }

    #cat("DEBUG: g1.idx = ", g1.idx, "\n")
    #cat("DEBUG: label after group.equal:\n"); print(label); cat("\n")

    # user-specified labels -- override everything!!
    user.idx <- which(nchar(partable$label) > 0L)
    label[user.idx] <- partable$label[user.idx]

    #cat("DEBUG: user.idx = ", user.idx, "\n")
    #cat("DEBUG: label after user.idx:\n"); print(label); cat("\n")

    # which labels do we need?
    if(type == "user") {
        idx <- 1:length(label)
    } else if(type == "free") {
        idx <- which(partable$free > 0L & !duplicated(partable$free))
    } else if(type == "unco") {
        idx <- which(partable$unco > 0L & !duplicated(partable$unco))
    } else {
        stop("argument `type' must be one of free, unco, or user")
    }

    label[idx]
}


getUserListFull <- function(partable=NULL, group=NULL) {

    # meanstructure + number of groups
    meanstructure <- any(partable$op == "~1")
    ngroups <- max(partable$group)

    # extract `names' of various types of variables:
    lv.names     <- lav_partable_vnames(partable, type="lv",  group=group)   # latent variables
    ov.names     <- lav_partable_vnames(partable, type="ov",  group=group)   # observed variables
    ov.names.x   <- lav_partable_vnames(partable, type="ov.x",group=group)   # exogenous x covariates
    ov.names.nox <- lav_partable_vnames(partable, type="ov.nox",group=group) # ov's without exo's
    lv.names.x   <- lav_partable_vnames(partable, type="lv.x",group=group)   # exogenous lv
    ov.names.y   <- lav_partable_vnames(partable, type="ov.y",group=group)   # dependent ov
    lv.names.y   <- lav_partable_vnames(partable, type="lv.y",group=group)   # dependent lv
    lvov.names.y <- c(ov.names.y, lv.names.y)
    ov.names.ord <- lav_partable_vnames(partable, type="ov.ord", group=group)


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
    if(any(partable$op == "~")) {
        # FIXME: better choices?
        eqs.names <- unique( c(partable$lhs[ partable$op == "~" ],
                               partable$rhs[ partable$op == "~" ]) )
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

    # 5. thresholds
    th.lhs <- th.rhs <- th.op <- character(0)
    if(length(ov.names.ord) > 0L) {
        tmp <- strsplit(lav_partable_vnames(partable, "th", group=group), "\\|")
        th.lhs <- sapply(tmp, function(x) x[1])
        th.rhs <- sapply(tmp, function(x) x[2])
        th.op  <- rep("|", length(th.lhs))
    }

    # 6. scaling parameters
    delta.lhs <- delta.rhs <- delta.op <- character(0)
    if(ngroups > 1L && length(ov.names.ord) > 0L) {
        delta.lhs <- ov.names.ord
        delta.rhs <- ov.names.ord
        delta.op  <- rep("~*~", length(delta.lhs))
    }

    # combine
    lhs <- c(l.lhs, ov.lhs, lv.lhs, r.lhs, int.lhs, th.lhs, delta.lhs)
    rhs <- c(l.rhs, ov.rhs, lv.rhs, r.rhs, int.rhs, th.rhs, delta.rhs)
     op <- c(l.op,  ov.op,  lv.op,  r.op,  int.op,  th.op,  delta.op)


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
                    auto.th         = FALSE,
                    auto.delta      = FALSE,
                    varTable        = NULL,
                    group.equal     = NULL,
                    ngroups         = 1L) {

    categorical <- FALSE

    ### DEFAULT elements: parameters that are typically not specified by
    ###                   users, but should typically be considered, 
    ###                   either free or fixed

    # extract `names' of various types of variables:
    lv.names     <- lav_partable_vnames(FLAT, type="lv")     # latent variables
    lv.names.r   <- lav_partable_vnames(FLAT, type="lv.regular") # regular latent variables
    ov.names     <- lav_partable_vnames(FLAT, type="ov")     # observed variables
    ov.names.x   <- lav_partable_vnames(FLAT, type="ov.x")   # exogenous x covariates 
    ov.names.nox <- lav_partable_vnames(FLAT, type="ov.nox")
    lv.names.x   <- lav_partable_vnames(FLAT, type="lv.x")   # exogenous lv
    ov.names.y   <- lav_partable_vnames(FLAT, type="ov.y")   # dependent ov
    lv.names.y   <- lav_partable_vnames(FLAT, type="lv.y")   # dependent lv
    #lvov.names.y <- c(ov.names.y, lv.names.y)
    lvov.names.y <- c(lv.names.y, ov.names.y)


    # get 'ordered' variables, either from FLAT or varTable
    ov.names.ord1 <- lav_partable_vnames(FLAT, type="ov.ord")
    # check if we have "|" for exogenous variables
    if(length(ov.names.ord1) > 0L) {
        idx <- which(ov.names.ord1 %in% ov.names.x)
        if(length(idx) > 0L) {
            warning("lavaan WARNING: thresholds are defined for exogenous variables: ", paste(ov.names.ord1[idx], collapse=" "))
        }
    }
 
    if(!is.null(varTable)) {
        ov.names.ord2 <- as.character(varTable$name[ varTable$type == "ordered" ])
        # remove fixed.x variables
        idx <- which(ov.names.ord2 %in% ov.names.x)
        if(length(idx) > 0L)
            ov.names.ord2 <- ov.names.ord2[-idx]
    } else {
        ov.names.ord2 <- character(0)
    }
    #### FIXME!!!!! ORDER!
    ov.names.ord <- unique(c(ov.names.ord1, ov.names.ord2))

    if(length(ov.names.ord) > 0L)
        categorical <- TRUE

    lhs <- rhs <- character(0)

    # 1. THRESHOLDS (based on varTable)
    nth <- 0L
    if(auto.th && length(ov.names.ord2) > 0L) {
        for(o in ov.names.ord2) {
            nth  <- varTable$nlev[ varTable$name == o ] - 1L
            if(nth < 1L) next
            lhs <- c(lhs, rep(o, nth))
            rhs <- c(rhs, paste("t", seq_len(nth), sep=""))
        }
        nth <- length(lhs)
    }

    # 2. default (residual) variances and covariances

    # a) (residual) VARIANCES (all ov's except exo, and regular lv's)
    if(auto.var) {
        ov.var <- ov.names.nox
        # auto-remove ordinal variables
        #idx <- match(ov.names.ord, ov.var)
        #if(length(idx)) ov.var <- ov.var[-idx]
        lhs <- c(lhs, ov.var, lv.names.r)
        rhs <- c(rhs, ov.var, lv.names.r)
    }

    # b) `independent` latent variable COVARIANCES (lv.names.x)
    if(auto.cov.lv.x && length(lv.names.x) > 1L) {
        tmp <- utils::combn(lv.names.x, 2)
        lhs <- c(lhs, tmp[1,]) # to fill upper.tri
        rhs <- c(rhs, tmp[2,])
    }

    # c) `dependent` latent variables COVARIANCES (lv.y.idx + ov.y.lv.idx)
    if(auto.cov.y && length(lvov.names.y) > 1L) {
        tmp <- utils::combn(lvov.names.y, 2L)
        lhs <- c(lhs, tmp[1,]) # to fill upper.tri
        rhs <- c(rhs, tmp[2,])
    }

    # d) exogenous x covariates: VARIANCES + COVARIANCES
    if(!categorical && (nx <- length(ov.names.x)) > 0L) {
        idx <- lower.tri(matrix(0, nx, nx), diag=TRUE)
        lhs <- c(lhs, rep(ov.names.x,  each=nx)[idx]) # fill upper.tri
        rhs <- c(rhs, rep(ov.names.x, times=nx)[idx])
    }
 
    # create 'op' (thresholds come first, then variances)
    op <- rep("~~", length(lhs)); op[seq_len(nth)] <- "|"

    # LATENT RESPONSE SCALES (DELTA)
    if(auto.delta && auto.th && length(ov.names.ord) > 0L && 
       length(lv.names) > 0L &&
       (ngroups > 1L || any(FLAT$op == "~*~"))) {
        lhs <- c(lhs, ov.names.ord)
        rhs <- c(rhs, ov.names.ord)
         op <- c(op,  rep("~*~", length(ov.names.ord)))
    }

    # 3. INTERCEPTS
    if(meanstructure) {
        if(categorical) {
            ov.int <- ov.names.nox
        } else {
            ov.int <- ov.names
        }
        # auto-remove ordinal variables
        idx <- which(ov.int %in% ov.names.ord)
        if(length(idx)) ov.int <- ov.int[-idx]

        int.lhs <- c(ov.int, lv.names)
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

    lv.names     <- lav_partable_vnames(FLAT, type="lv")     # latent variables
    ov.names     <- lav_partable_vnames(FLAT, type="ov")     # observed variables

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
        if(categorical) {
            # zero intercepts/means ordinal variables
                   ov.int.idx <- which(op == "~1" &
                                       lhs %in% ov.names.ord &
                                       user == 0L)
            ustart[ov.int.idx] <- 0.0
              free[ov.int.idx] <- 0L
        }
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

    # 5b. residual variances of ordinal variables?
    if(length(ov.names.ord) > 0L) {
        ord.idx <- which(lhs %in% ov.names.ord &
                         op == "~~" &
                         lhs == rhs)
        ustart[ord.idx] <- 1L ## FIXME!! or 0?? (0 breaks ex3.12)
          free[ord.idx] <- 0L
    }

    # 5c latent response scales of ordinal variables?
    if(length(ov.names.ord) > 0L) {
        delta.idx <- which(op == "~*~")
        ustart[delta.idx] <- 1.0
          free[delta.idx] <- 0L
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
                   ("intercepts" %in% group.equal ||
                    "thresholds" %in% group.equal) &&
                   !("means" %in% group.equal) ) {
                      free[ int.idx ] <- 1L
                    ustart[ int.idx ] <- as.numeric(NA)
                }
            }

            # latent response scaling
            if(any(op == "~*~" & group == g) &&
               ("thresholds" %in% group.equal)) {
                delta.idx <- which(op == "~*~" & group == g)
                  free[ delta.idx ] <- 1L
                ustart[ delta.idx ] <- as.numeric(NA)
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

independenceModel <- function(ov.names=NULL, ov=NULL, 
                              ov.names.x=NULL, sample.cov=NULL,
                              meanstructure=FALSE, sample.mean=NULL,
                              sample.th=NULL,
                              fixed.x=TRUE) {

    ngroups <- length(ov.names)
    ov.names.nox <- lapply(as.list(1:ngroups), function(g) 
                    ov.names[[g]][ !ov.names[[g]] %in% ov.names.x[[g]] ])

    lhs <- rhs <- op <- character(0)
    group <- free <- exo <- integer(0)
    ustart <- numeric(0)

    categorical <- any(ov$type == "ordered")

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

        # ordered? fix variances, add thresholds
        ord.names <- character(0L)
        if(categorical) {
            ord.names <- ov$name[ ov$type == "ordered" ]
            # only for this group
            ord.names <- ov.names[[g]][ which(ov.names[[g]] %in% ord.names) ]
            
            if(length(ord.names) > 0L) {
                # fix variances to 1.0
                idx <- which(lhs %in% ord.names & op == "~~" & lhs == rhs)
                ustart[idx] <- 1.0
                free[idx] <- 0L

                # add thresholds
                lhs.th <- character(0); rhs.th <- character(0)
                for(o in ord.names) {
                    nth  <- ov$nlev[ ov$name == o ] - 1L
                    if(nth < 1L) next
                    lhs.th <- c(lhs.th, rep(o, nth))
                    rhs.th <- c(rhs.th, paste("t", seq_len(nth), sep=""))
                }
                nel   <- length(lhs.th)
                lhs   <- c(lhs, lhs.th)
                rhs   <- c(rhs, rhs.th)
                 op   <- c(op, rep("|", nel))
                group <- c(group, rep(g, nel))
                 free <- c(free, rep(1L, nel))
                  exo <- c(exo, rep(0L, nel))
               th.start <- rep(0, nel)
               #th.start <- sample.th[[g]] ### FIXME::: ORDER??? ONLY ORD!!!
               ustart <- c(ustart, th.start)
            }
        }
        # meanstructure?
        if(meanstructure) {
            # auto-remove ordinal variables
            ov.int <- ov.names.nox[[g]]
            idx <- which(ov.int %in% ord.names)
            if(length(idx)) ov.int <- ov.int[-idx]

            nel <- length(ov.int)
            lhs   <- c(lhs, ov.int)
             op   <- c(op, rep("~1", nel))
            rhs   <- c(rhs, rep("", nel))
            group <- c(group, rep(g,  nel))
            free  <- c(free,  rep(1L, nel))
            exo   <- c(exo,   rep(0L, nel))
            # starting values
            if(!is.null(sample.mean)) {
                sample.int.idx <- match(ov.int, ov.names[[g]])
                ustart <- c(ustart, sample.mean[[g]][sample.int.idx])
            } else {
                ustart <- c(ustart, rep(as.numeric(NA), length(ov.int)))
            }
        }

        # fixed.x exogenous variables?
        if(!categorical && (nx <- length(ov.names.x[[g]])) > 0L) {
            idx <- lower.tri(matrix(0, nx, nx), diag=TRUE)
            nel <- sum(idx)
            lhs    <- c(lhs, rep(ov.names.x[[g]],  each=nx)[idx]) # upper.tri
             op    <- c(op, rep("~~", nel))
            rhs    <- c(rhs, rep(ov.names.x[[g]], times=nx)[idx])
            if(fixed.x) {
                free   <- c(free,  rep(0L, nel))
            } else {
                free   <- c(free,  rep(1L, nel))
            }
            exo    <- c(exo,   rep(1L, nel))
            group  <- c(group, rep(g,  nel))
            ustart <- c(ustart, rep(as.numeric(NA), nel))

            # meanstructure?
            if(meanstructure) {
                lhs    <- c(lhs,    ov.names.x[[g]])
                 op    <- c(op,     rep("~1", nx))
                rhs    <- c(rhs,    rep("", nx))
                group  <- c(group,  rep(g,  nx))
                if(fixed.x) {
                    free   <- c(free,   rep(0L, nx))
                } else {
                    free   <- c(free,   rep(1L, nx))
                }
                exo    <- c(exo,    rep(1L, nx))
                ustart <- c(ustart, rep(as.numeric(NA), nx))
            }
        }

        if(categorical && (nx <- length(ov.names.x[[g]])) > 0L) {
            # add regressions
            lhs <- c(lhs, rep("dummy", nx))
             op <- c( op, rep("~", nx))
            rhs <- c(rhs, ov.names.x[[g]])
            # add 3 dummy lines
            lhs <- c(lhs, "dummy"); op <- c(op, "=~"); rhs <- c(rhs, "dummy")
            lhs <- c(lhs, "dummy"); op <- c(op, "~~"); rhs <- c(rhs, "dummy")
            lhs <- c(lhs, "dummy"); op <- c(op, "~1"); rhs <- c(rhs, "")

            exo    <- c(exo,   rep(0L, nx + 3L))
            group  <- c(group, rep(g,  nx + 3L))
            free   <- c(free,  rep(0L, nx + 3L))
            ustart <- c(ustart, rep(0, nx)); ustart <- c(ustart, c(0,1,0))
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

# given two models, M1 and M0, where M0 is nested in M1,
# create a function 'af(x)' where 'x' is the full parameter vector of M1
# and af(x) returns the evaluated restrictions under M0).
# The 'jacobian' of this function 'A' will be used in the anova
# anova() function, and elsewhere
getConstraintsFunction <- function(p1, p0) {

    # check for inequality constraints
    if(any(c(p0$op,p1$op) %in% c(">","<"))) 
        stop("lavaan ERROR: anova() can not handle inequality constraints; use InformativeTesting() instead")

    npar.p1 <- max(p1$free)
    npar.p0 <- max(p0$free)
    npar.diff <- npar.p1 - npar.p0

    con.function <- function() NULL
    formals(con.function) <- alist(x=, ...=)
    BODY.txt <- paste("{\nout <- rep(NA, ", npar.diff, ")\n", sep="")

    # for each free parameter in p1, we 'check' is it is somehow 
    # restricted in p0
    ncon <- 0L; EQ.ID <- integer(); EQ.P1 <- integer()
    for(i in seq_len(npar.p1)) {
        idx <- which(p1$free == i)[1L]
        lhs <- p1$lhs[idx]; op <- p1$op[idx]; rhs <- p1$rhs[idx]
        group <- p1$group[idx]
        p0.idx <- which(p0$lhs == lhs & p0$op == op & p0$rhs == rhs &
                        p0$group == group)
        if(length(p0.idx) == 0L) {
            # this parameter is not to be found in M0, we will
            # assume that is fixed to zero
            ncon <- ncon + 1L 
            BODY.txt <- paste(BODY.txt,
                "out[", ncon, "] = x[", i, "] - 0\n", sep="")
            next
        }
        if(p0$free[p0.idx] == 0L) {
            ncon <- ncon + 1L
            # simple fixed value
            BODY.txt <- paste(BODY.txt,
                "out[", ncon, "] = x[", i, "] - ",
                                   p0$ustart[p0.idx], "\n", sep="")
        } 
        # how to deal with *new* equality constraints?
        # check if p0 has an eq.id not yet in EQ.ID (while p1 eq.id is empty)
        if(p0$eq.id[p0.idx] != 0 && p1$eq.id[idx] == 0) {
            # if not in EQ.ID, put it there and continue
            EQ <- p0$eq.id[p0.idx]
            if(EQ %in% EQ.ID) {
                # add constraint
                ncon <- ncon + 1L
                BODY.txt <- paste(BODY.txt,
                                  "out[", ncon, "] = x[", EQ.P1[EQ.ID == EQ], 
                                                "] - x[", i, "]\n", sep="")
            } else {
                EQ.ID <- c(EQ.ID, EQ)
                EQ.P1 <- c(EQ.P1,  i)
            }
        }
    }

    # add all NEW constraints using "==" (not already in p1)
    p0.con.idx <- which(p0$op == "==")
    if(length(p0.con.idx) > 0L) {
        # first, remove those also present in p1
        del.idx <- integer(0L)
        for(con in 1:length(p0.con.idx)) {
            p0.idx <- p0.con.idx[con]
            p1.idx <- which(p1$op == "==" &
                            p1$lhs == p0$lhs[p0.idx] &
                            p1$rhs == p0$rhs[p0.idx])
            if(length(p1.idx) > 0L)
                del.idx <- c(del.idx, con)
        }
        if(length(del.idx) > 0L)
            p0.con.idx <- p0.con.idx[-del.idx]   
    }
    eq.idx <- p0.con.idx
    if(length(eq.idx) > 0L) {
        def.idx <- which(p0$op == ":=")
        # first come the variable definitions
        if(length(def.idx) > 0L) {
            for(i in 1:length(def.idx)) {
                lhs <- p0$lhs[ def.idx[i] ]
                rhs <- p0$rhs[ def.idx[i] ]
                def.string <- rhs
                # coerce to expression to extract variable names
                def.labels <- all.vars( parse(file="", text=def.string) )
                # get corresponding 'x' indices
                def.x.idx  <- p0$free[match(def.labels, p0$label)]
                def.x.lab  <- paste("x[", def.x.idx, "]",sep="")
                # put both the labels and the expression in the function BODY
                BODY.txt <- paste(BODY.txt,
                    paste(def.labels, "=",def.x.lab, collapse=";"),"\n",
                    lhs, " = ", def.string, "\n", sep="")
            }
        }

        for(i in 1:length(eq.idx)) {
            lhs <- p0$lhs[ eq.idx[i] ]
            rhs <- p0$rhs[ eq.idx[i] ]
            if(rhs == "0") {
                eq.string <- lhs
            } else {
                eq.string <- paste(lhs, "- (", rhs, ")", sep="")
            }
            # coerce to expression to extract variable names
            eq.labels <- all.vars( parse(file="", text=eq.string) )
            # get corresponding 'x' indices
            if(length(def.idx) > 0L) {
                # remove def.names from ineq.labels
                def.names <- as.character(p0$lhs[def.idx])
                d.idx <- which(eq.labels %in% def.names)
                if(length(d.idx) > 0) eq.labels <- eq.labels[-d.idx]
            }
            if(length(eq.labels) > 0L) {
                eq.x.idx  <- p0$free[match(eq.labels, p0$label)]
                if(any(is.na(eq.x.idx))) {
                    stop("lavaan ERROR: unknown label(s) in equality constraint: ",
                         paste(eq.labels[which(is.na(eq.x.idx))], collapse=" "))
                }
                if(any(eq.x.idx == 0)) {
                    stop("lavaan ERROR: non-free parameter(s) in inequality constraint: ",
                        paste(eq.labels[which(eq.x.idx == 0)], collapse=" "))
                }
                eq.x.lab  <- paste("x[", eq.x.idx, "]",sep="")
                # put both the labels and the expression in the function BODY
                BODY.txt <- paste(BODY.txt,
                    paste(eq.labels, "=", eq.x.lab, collapse=";"),"\n",
                    "out[", i, "] = ", eq.string, "\n", sep="")
            } else {
                BODY.txt <- paste(BODY.txt,
                    "out[", i, "] = ", eq.string, "\n", sep="")
            }
        }
               
    }

    # wrap function
    BODY.txt <- paste(BODY.txt, "return(out)\n}\n", sep="")
    body(con.function) <- parse(file="", text=BODY.txt)

    con.function
}

# user visible function to add 'matrix' entries in the parameter table
lavMatrixRepresentation <- function(partable, representation = "LISREL",
                                    as.data.frame. = TRUE) {

    # get model matrices
    if(representation == "LISREL") {
        REP <- representation.LISREL(partable, target=NULL, extra=FALSE)
    } else {
        stop("lavaan ERROR: only representation \"LISREL\" has been implemented.")
    }
    
    partable$mat <- REP$mat
    partable$row <- REP$row
    partable$col <- REP$col

    if(as.data.frame.) {
        partable <- as.data.frame(partable, stringsAsFactors=FALSE)
        class(partable) <- c("lavaan.data.frame", "data.frame")
    }

    partable
}

