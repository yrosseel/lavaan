# YR - 26 Nov 2013: generate partable for the unrestricted model
# YR - 19 Mar 2017: handle twolevel model
lav_partable_unrestricted <- function(lavobject      = NULL,
                                      # if no object is available,
                                      lavdata        = NULL,
                                      lavpta         = NULL,
                                      lavoptions     = NULL,
                                      lavsamplestats = NULL,
                                      # optional user-provided sample stats
                                      sample.cov     = NULL,
                                      sample.mean    = NULL,
                                      sample.slopes  = NULL,
                                      sample.th      = NULL,
                                      sample.th.idx  = NULL) {

    # grab everything from lavaan lavobject
    if(!is.null(lavobject)) {
        stopifnot(inherits(lavobject, "lavaan"))

        lavdata <- lavobject@Data
        lavoptions <- lavobject@Options
        lavsamplestats <- lavobject@SampleStats
        lavpta <- lavobject@pta
    }

    # conditional.x ? check res.cov[[1]] slot
    conditional.x <- FALSE
    if(!is.null(lavsamplestats) && !is.null(lavsamplestats@res.cov[[1]])) {
        conditional.x <- TRUE
    } else if(!is.null(lavoptions) && lavoptions$conditional.x) {
        conditional.x <- TRUE
    }

    # get sample statistics, all groups
    SAMPLE.cov <- sample.cov
    if(is.null(SAMPLE.cov) && !is.null(lavsamplestats)) {
        if(conditional.x) {
            SAMPLE.cov <- lavsamplestats@res.cov
        } else {
            SAMPLE.cov <- lavsamplestats@cov
        }
    }

    SAMPLE.mean <- sample.mean
    if(is.null(SAMPLE.mean) && !is.null(lavsamplestats)) {
        if(conditional.x) {
            SAMPLE.mean <- lavsamplestats@res.int
        } else {
            SAMPLE.mean <- lavsamplestats@mean
        }
    }

    SAMPLE.slopes <- sample.slopes
    if(conditional.x && is.null(SAMPLE.slopes) && !is.null(lavsamplestats)) {
        SAMPLE.slopes <- lavsamplestats@res.slopes
    }

    SAMPLE.th <- sample.th
    if(is.null(SAMPLE.th) && !is.null(lavsamplestats)) {
        if(conditional.x) {
             SAMPLE.th <- lavsamplestats@res.th
         } else {
             SAMPLE.th <- lavsamplestats@th
         }
    }

    SAMPLE.th.idx <- sample.th.idx
    if(is.null(SAMPLE.th.idx) && !is.null(lavsamplestats)) {
         SAMPLE.th.idx <- lavsamplestats@th.idx
    }

    #ov.names      <- lavdata@ov.names
    ov            <- lavdata@ov
    #ov.names.x    <- lavdata@ov.names.x
    meanstructure <- lavoptions$meanstructure
    categorical   <- any(ov$type == "ordered")
    ngroups       <- lavdata@ngroups
    nlevels       <- lavdata@nlevels

    # what with fixed.x? 
    # - does not really matter; fit will be saturated any way
    # - fixed.x = TRUE may avoid convergence issues with non-numeric 
    #             x-covariates
    #if(lavoptions$mimic %in% c("lavaan", "Mplus")) {
        fixed.x = lavoptions$fixed.x
    #} else if(lavoptions$mimic == "EQS") {
        # always ignore fixed.x
    #    ov.names.x = NULL
    #    fixed.x = FALSE
    #} else if(lavoptions$mimic == "LISREL") {
    #    # always ignore fixed.x??? CHECKME!!
    #    ov.names.x = NULL
    #    fixed.x = FALSE
    #}

    # if multilevel, ALWAYS fixed.x = FALSE
    if(nlevels > 1L) {
        fixed.x       <- FALSE
        conditional.x <- FALSE
        categorical   <- FALSE # for now
    }

    lhs <- rhs <- op <- character(0)
    group <- block <- level <- free <- exo <- integer(0)
    ustart <- numeric(0)

    # block number
    b <- 0L
    for(g in 1:ngroups) {

        # only for multilevel
        if(nlevels > 1L) {
            OV.NAMES <- lavdata@ov.names[[g]]
            YLp <-  lavsamplestats@YLp[[g]]
        }

        sample.cov    <- SAMPLE.cov[[g]]
        sample.mean   <- SAMPLE.mean[[g]]
        sample.slopes <- SAMPLE.slopes[[g]]
        sample.th     <- SAMPLE.th[[g]]
        sample.th.idx <- SAMPLE.th.idx[[g]]

        # force sample.cov to be pd
        if(!is.null(sample.cov)) {
            sample.cov <- lav_matrix_symmetric_force_pd(sample.cov)
        }

        for(l in 1:nlevels) {

            # block
            b <- b + 1L

            # ov.names for this block
            if(is.null(lavpta)) { # only data was used
                ov.names     <- lavdata@ov.names[[g]]
                ov.names.x   <- lavdata@ov.names.x[[g]]
                ov.names.nox <- ov.names[!ov.names %in% ov.names.x]
            } else {
                ov.names     <- lavpta$vnames$ov[[b]]
                ov.names.x   <- lavpta$vnames$ov.x[[b]]
                ov.names.nox <- lavpta$vnames$ov.nox[[b]]
            }

            # only for multilevel, overwrite sample.cov and sample.mean
            if(nlevels > 1L) {

                # NO FIXED.X!! (because we can not compute samplestats for
                # them in the unbalanced case)
                ov.names.x <- character(0L)
                ov.names.nox <- ov.names

                block.idx <- match(ov.names, OV.NAMES)

                if(l == 1L) {
                    sample.cov  <- YLp[[2]]$Sigma.W[block.idx, block.idx, 
                                                    drop = FALSE]
                    sample.mean <- NULL
                } else {
                    sample.cov  <- YLp[[2]]$Sigma.B[block.idx, block.idx, 
                                                    drop = FALSE]
                    sample.mean <- YLp[[2]]$Mu.B[block.idx]
                }
            } 


            # a) VARIANCES (all ov's, if !conditional.x, also exo's)
            nvar  <- length(ov.names)
    
            lhs   <- c(lhs, ov.names)
             op   <- c(op, rep("~~", nvar))
            rhs   <- c(rhs, ov.names)
            block <- c(block, rep(b,  nvar))
            group <- c(group, rep(g,  nvar))
            level <- c(level, rep(l,  nvar))
            free  <- c(free,  rep(1L, nvar))
            exo   <- c(exo,   rep(0L, nvar))

            # starting values -- variances
            if(!is.null(sample.cov)) {
                ustart <- c(ustart, diag(sample.cov))
            } else {
                ustart <- c(ustart, rep(as.numeric(NA), nvar))
            }

            # COVARIANCES!
            pstar <- nvar*(nvar-1)/2
            if(pstar > 0L) { # only if more than 1 variable
                tmp <- utils::combn(ov.names, 2)
                lhs <- c(lhs, tmp[1,]) # to fill upper.tri
                 op <- c(op,   rep("~~", pstar))
                rhs <- c(rhs, tmp[2,])
                block <- c(block, rep(b,  pstar))
                group <- c(group, rep(g,  pstar))
                level <- c(level, rep(l,  pstar))
                free  <- c(free,  rep(1L, pstar))
                exo   <- c(exo,   rep(0L, pstar))
            }
    
            # starting values -- covariances
            if(!is.null(sample.cov)) {
                ustart <- c(ustart, lav_matrix_vech(sample.cov, 
                                                    diagonal = FALSE))
            } else {
                ustart <- c(ustart, rep(as.numeric(NA), pstar))
            }

            # ordered? fix variances, add thresholds
            ord.names <- character(0L)
            if(categorical) {
                ord.names <- ov$name[ ov$type == "ordered" ]
                # only for this group
                ord.names <- ov.names[ which(ov.names %in% ord.names) ]
                
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
                    block <- c(block, rep(b,  nel))
                    group <- c(group, rep(g,  nel))
                    level <- c(level, rep(l,  nel))
                     free <- c(free, rep(1L, nel))
                      exo <- c(exo, rep(0L, nel))
    
                    # starting values
                    if(!is.null(sample.th) && !is.null(sample.th.idx)) {
                        th.start <- sample.th[ sample.th.idx > 0L ]
                        ustart <- c(ustart, th.start)
                    } else {
                        ustart <- c(ustart, rep(as.numeric(NA), nel))
                    }
                }
            }

            # meanstructure?
            if(meanstructure) {
                # auto-remove ordinal variables
                ov.int <- ov.names
                idx <- which(ov.int %in% ord.names)
                if(length(idx)) ov.int <- ov.int[-idx]

                nel <- length(ov.int)
                lhs   <- c(lhs, ov.int)
                 op   <- c(op, rep("~1", nel))
                rhs   <- c(rhs, rep("", nel))
                block <- c(block, rep(b,  nel))
                group <- c(group, rep(g,  nel))
                level <- c(level, rep(l,  nel))
                # if multilevel, level=1 has fixed zeroes
                if(nlevels > 1L && l == 1L) {
                    free  <- c(free,  rep(0L, nel))
                } else {
                    free  <- c(free,  rep(1L, nel))
                }
                exo   <- c(exo,   rep(0L, nel))
                # starting values
                if(nlevels > 1L && l == 1L) {
                    ustart <- c(ustart, rep(0, length(ov.int)))
                } else if(!is.null(sample.mean)) {
                    sample.int.idx <- match(ov.int, ov.names)
                    ustart <- c(ustart, sample.mean[sample.int.idx])
                } else {
                    ustart <- c(ustart, rep(as.numeric(NA), length(ov.int)))
                }
            }

            # categorical? insert means as fixed-to-zero parameters
            # since 0.5-17
            # always add ~*~ parameters (since 0.6-1)
            if(categorical) {
                ov.int <- ov.names
                idx <- which(ov.int %in% ord.names)
                ov.int <- ov.int[idx]

                nel <- length(ov.int)
                lhs   <- c(lhs, ov.int)
                 op   <- c(op, rep("~1", nel))
                rhs   <- c(rhs, rep("", nel))
                block <- c(block, rep(b,  nel))
                group <- c(group, rep(g,  nel))
                level <- c(level, rep(l,  nel))
                free  <- c(free,  rep(0L, nel))
                exo   <- c(exo,   rep(0L, nel))
               ustart <- c(ustart, rep(0, nel))

                # ~*~
                nel <- length(ov.int)
                lhs   <- c(lhs, ov.int)
                 op   <- c(op, rep("~*~", nel))
                rhs   <- c(rhs, ov.int)
                block <- c(block, rep(b,  nel))
                group <- c(group, rep(g,  nel))
                level <- c(level, rep(l,  nel))
                free  <- c(free,  rep(0L, nel))
                exo   <- c(exo,   rep(0L, nel))
               ustart <- c(ustart, rep(1, nel))
            }
 

            # fixed.x exogenous variables?
            if(!conditional.x && 
               fixed.x && (nx <- length(ov.names.x)) > 0L) {
                # fix variances/covariances 
                exo.idx <- which(rhs %in% ov.names.x &
                                 lhs %in% ov.names.x &
                                 op == "~~" & group == g)
                exo[exo.idx] <- 1L
                free[exo.idx] <- 0L
    
                # fix means
                exo.idx <- which(rhs %in% ov.names.x &
                                 op == "~1" & group == g)
                exo[exo.idx] <- 1L
                free[exo.idx] <- 0L
            }

            # conditional.x?
            if(conditional.x && (nx <- length(ov.names.x)) > 0L) {
                nnox <- length(ov.names.nox)
                nel  <- nnox * nx
    
                lhs <- c(lhs, rep(ov.names.nox, times = nx))
                 op <- c(op,  rep("~", nel))
                rhs <- c(rhs, rep(ov.names.x, each = nnox))
                block <- c(block, rep(b,  nel))
                group <- c(group, rep(g,  nel))
                level <- c(level, rep(l,  nel))
                free <- c(free,  rep(1L, nel))
                exo <- c(exo,   rep(1L, nel))

                # starting values -- slopes
                if(!is.null(sample.slopes)) {
                    ustart <- c(ustart, lav_matrix_vec(sample.slopes))
                } else {
                    ustart <- c(ustart, rep(as.numeric(NA), nel))
                }
            }
        } # levels
    } # ngroups

    # free counter
    idx.free <- which(free > 0)
    free[idx.free] <- 1:length(idx.free)

    LIST   <- list(     id          = 1:length(lhs),
                        lhs         = lhs,
                        op          = op,
                        rhs         = rhs,
                        user        = rep(1L,  length(lhs)),
                        block       = block,
                        group       = group,
                        level       = level,
                        free        = free,
                        ustart      = ustart,
                        exo         = exo #,
                        #label       = rep("",  length(lhs))
                        #eq.id       = rep(0L,  length(lhs)),
                        #unco        = free
                   )

    
    # keep level column if no levels? (no for now)
    if(nlevels < 2L) {
        LIST$level <- NULL
    }

    LIST

}
