# create `full' parameter table, containing (almost) all parameters
# that could be free
#
# main motivation: univariate scores tests (modification indices)
#
lav_partable_full <- function(partable = NULL,
                              strict.exo = FALSE,
                              free = FALSE, start = FALSE) {

    # check minimum requirements: lhs, op, rhs
    stopifnot( !is.null(partable$lhs), 
               !is.null(partable$op), 
               !is.null(partable$rhs) )

    # meanstructure
    meanstructure <- any(partable$op == "~1")

    # number of blocks
    nblocks <- lav_partable_nblocks(partable)

    # extract `names' of various types of variables:
    lv.names     <- lav_partable_vnames(partable, type="lv")   
    ov.names     <- lav_partable_vnames(partable, type="ov")
    ov.names.x   <- lav_partable_vnames(partable, type="ov.x")
    ov.names.nox <- lav_partable_vnames(partable, type="ov.nox")
    lv.names.x   <- lav_partable_vnames(partable, type="lv.x")
    ov.names.y   <- lav_partable_vnames(partable, type="ov.y") 
    lv.names.y   <- lav_partable_vnames(partable, type="lv.y")
    lvov.names.y <- c(ov.names.y, lv.names.y)
    ov.names.ord <- lav_partable_vnames(partable, type="ov.ord")
    ov.names.ind <- lav_partable_vnames(partable, type="ov.ind")

    # eqs.y, eqs.x
    if(any(partable$op == "~")) {
        eqs.names <- unique( c(partable$lhs[partable$op == "~"],
                               partable$rhs[partable$op == "~"]) )
        eqs.y <- eqs.names
        if(strict.exo) {
            x.idx <- which(eqs.names %in% ov.names.x)
            if(length(x.idx) > 0L) {
                eqs.y <- eqs.names[-x.idx]
            }
        }
        eqs.x <- eqs.names
    } else {
        eqs.y <- character(0L)
        eqs.x <- character(0L)
    }


    # 1 "=~"
    l.lhs <- r.rhs <- op <- character(0)
    l.lhs <- rep(lv.names, each=length(ov.names.nox))
    l.rhs <- rep(ov.names.nox, times=length(lv.names))

    # remove factor ~ eqs.y combinations, if any
    # because they also appear as a regression
    bad.idx <- which( l.lhs %in% lv.names &
                      l.rhs %in% eqs.y)
    if(length(bad.idx) > 0L) {
        l.lhs <- l.lhs[-bad.idx]
        l.rhs <- l.rhs[-bad.idx]
    }

    l.op  <- rep("=~", length(l.lhs))

    # 2a. "~~" ov ## FIXME: ov.names.nox or ov.names??
    ov.lhs <- ov.rhs <- ov.op <- character(0)
    #if(strict.exo) {
        OV <- ov.names.nox
    #} else {
    #    OV <- ov.names
    #}
    nx <- length(OV)
    idx <- lower.tri(matrix(0, nx, nx), diag=TRUE)
    ov.lhs <- rep(OV,  each=nx)[idx] # fill upper.tri
    ov.rhs <- rep(OV, times=nx)[idx]
    ov.op  <- rep("~~", length(ov.lhs))

    # exo ~~
    if(!strict.exo && length(ov.names.x) > 0L) {
        OV <- ov.names.x
        nx <- length(OV)
        idx <- lower.tri(matrix(0, nx, nx), diag=TRUE)
        more.lhs <- rep(OV,  each=nx)[idx] # fill upper.tri
        more.rhs <- rep(OV, times=nx)[idx]
        ov.lhs <- c(ov.lhs, more.lhs)
        ov.rhs <- c(ov.rhs, more.rhs)
        ov.op  <- c(ov.op,  rep("~~", length(more.lhs)))
    }

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

        r.lhs <- rep(eqs.y, each  = length(eqs.x))
        r.rhs <- rep(eqs.x, times = length(eqs.y))

        # remove self-arrows
        idx <- which(r.lhs == r.rhs)
        if(length(idx) > 0L) {
            r.lhs <- r.lhs[-idx]
            r.rhs <- r.rhs[-idx]
        }

        # remove indicator ~ factor if they exist
        bad.idx <- which(r.lhs %in% ov.names.ind &
                         r.rhs %in% lv.names)
        if(length(bad.idx) > 0L) {
            r.lhs <- r.lhs[-bad.idx]
            r.rhs <- r.rhs[-bad.idx]
        }

        r.op <- rep("~", length(r.rhs))
    }

    # 4. intercepts
    int.lhs <- int.rhs <- int.op <- character(0)
    if(meanstructure) {
        if(strict.exo) {
            int.lhs <- c(ov.names.nox, lv.names)
        } else {
            int.lhs <- c(ov.names, lv.names)
        }
        int.rhs <- rep("",   length(int.lhs))
        int.op  <- rep("~1", length(int.lhs))
    }

    # 5. thresholds
    th.lhs <- th.rhs <- th.op <- character(0)
    if(length(ov.names.ord) > 0L) {
        tmp <- strsplit(lav_partable_vnames(partable, "th"), "\\|")
        th.lhs <- sapply(tmp, function(x) x[1])
        th.rhs <- sapply(tmp, function(x) x[2])
        th.op  <- rep("|", length(th.lhs))
    }

    # 6. scaling parameters
    delta.lhs <- delta.rhs <- delta.op <- character(0)
    if(nblocks > 1L && length(ov.names.ord) > 0L) {
        delta.lhs <- ov.names.ord
        delta.rhs <- ov.names.ord
        delta.op  <- rep("~*~", length(delta.lhs))
    }

    # combine
    lhs <- c(l.lhs, ov.lhs, lv.lhs, r.lhs, int.lhs, th.lhs, delta.lhs)
    rhs <- c(l.rhs, ov.rhs, lv.rhs, r.rhs, int.rhs, th.rhs, delta.rhs)
     op <- c(l.op,  ov.op,  lv.op,  r.op,  int.op,  th.op,  delta.op)


    # multiple blocks!
    block <- 1L
    if(nblocks > 1) {
        block   <- rep(1:nblocks, each = length(lhs))
        group   <- rep(1:nblocks, each = length(lhs)) # group == blocks for now
        lhs     <- rep(lhs,     times = nblocks)
        op      <- rep(op,      times = nblocks)
        rhs     <- rep(rhs,     times = nblocks)
    } else {
        group   <- block
    }

    LIST <- data.frame(lhs = lhs, op = op, rhs = rhs, block = block,
                       group = group, stringsAsFactors = FALSE)

    if(free) {
        LIST$free <- rep(0L, nrow(LIST))
    }

    if(start) {
        LIST$start <- rep(0, nrow(LIST))
    }

    LIST
}

