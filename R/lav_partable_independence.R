# generate parameter table for an independence model
lav_partable_independence <- function(lavobject      = NULL, 
                                      lavdata        = NULL,
                                      lavoptions     = NULL,
                                      lavsamplestats = NULL,
                                      # optional user-provided sample stats
                                      sample.cov       = NULL,
                                      sample.mean      = NULL,
                                      sample.th        = NULL, 
                                      sample.th.idx    = NULL) {

    # grab everything from lavaan lavobject
    if(!is.null(lavobject)) {
        stopifnot(inherits(lavobject, "lavaan"))

        lavdata <- lavobject@Data
        lavoptions <- lavobject@Options
        lavsamplestats <- lavobject@SampleStats
    }

    # conditional.x ? check res.cov[[1]] slot
    conditional.x <- FALSE
    if(!is.null(lavsamplestats) && !is.null(lavsamplestats@res.cov[[1]])) {
        conditional.x <- TRUE
    } else if(!is.null(lavoptions) && lavoptions$conditional.x) {
        conditional.x <- TRUE
    }


    # if user-based moments are given, use these
    if(is.null(sample.cov) && !is.null(lavsamplestats)) {
        if(conditional.x) {
            sample.cov <- lavsamplestats@res.cov
        } else {
            sample.cov <- lavsamplestats@cov
        }
    }
    if(is.null(sample.mean) && !is.null(lavsamplestats)) {
        if(conditional.x) {
            sample.mean <- lavsamplestats@res.int
        } else {
            sample.mean <- lavsamplestats@mean
        }
    }
    if(is.null(sample.th) && !is.null(lavsamplestats)) {
         if(conditional.x) {
             sample.th <- lavsamplestats@res.th
         } else {
             sample.th <- lavsamplestats@th
         }
    }

    if(is.null(sample.th.idx) && !is.null(lavsamplestats)) {
         sample.th.idx <- lavsamplestats@th.idx
    }

    ov.names         = lavdata@ov.names
    ov               = lavdata@ov
    meanstructure    = lavoptions$meanstructure
    parameterization = lavoptions$parameterization

    # what with fixed.x?
    if(lavoptions$mimic %in% c("lavaan", "Mplus")) {
        fixed.x = lavoptions$fixed.x
        ov.names.x = lavdata@ov.names.x
    } else if(lavoptions$mimic == "EQS") {
        # always ignore fixed.x
        fixed.x = FALSE
        ov.names.x = NULL
    } else if(lavoptions$mimic == "LISREL") {
        # always ignore fixed.x??? CHECKME!!
        fixed.x = FALSE
        ov.names.x = NULL
    }

    ngroups <- length(ov.names)


    lhs <- rhs <- op <- character(0)
    group <- free <- exo <- integer(0)
    ustart <- numeric(0)

    categorical <- any(ov$type == "ordered")

    for(g in 1:ngroups) {

        # a) VARIANCES (all ov's, if !conditional.x, also exo's)
        nvar  <- length(ov.names[[g]])
        lhs   <- c(lhs, ov.names[[g]])
         op   <- c(op, rep("~~", nvar))
            rhs   <- c(rhs, ov.names[[g]])
        group <- c(group, rep(g,  nvar))
        free  <- c(free,  rep(1L, nvar))
        exo   <- c(exo,   rep(0L, nvar))

        # starting values
        if(!is.null(sample.cov)) {
            ustart <- c(ustart, diag(sample.cov[[g]]))
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

                # starting values
                if(!is.null(sample.th) && !is.null(sample.th.idx)) {
                    th.start <- sample.th[[g]][ sample.th.idx[[g]] > 0L ]
                    ustart <- c(ustart, th.start)
                } else {
                    ustart <- c(ustart, rep(as.numeric(NA), nel))
                }

                # add delta
                if(parameterization == "theta") {
                   lhs.delta <- character(0); rhs.delta <- character(0)
                   lhs.delta <- ov.names[[g]]
                   nel   <- length(lhs.delta)
                   lhs   <- c(lhs, lhs.delta)
                   rhs   <- c(rhs, lhs.delta)
                    op   <- c(op, rep("~*~", nel))
                   group <- c(group, rep(g, nel))
                    free <- c(free, rep(0L, nel))
                     exo <- c(exo, rep(0L, nel))
                  delta.start <- rep(1, nel)
                  ustart <- c(ustart, delta.start)
                }

                # add mean/intercept, but fix to zero
                lhs.int <- ord.names
                nel   <- length(lhs.int)
                rhs.int <- rep("", nel)
                lhs   <- c(lhs, lhs.int)
                rhs   <- c(rhs, rhs.int)
                 op   <- c(op, rep("~1", nel))
                group <- c(group, rep(g, nel))
                 free <- c(free, rep(0L, nel))
                  exo <- c(exo, rep(0L, nel))
               int.start <- rep(0, nel)
               ustart <- c(ustart, int.start)
            }
        }

        # meanstructure?
        if(meanstructure) {
            # auto-remove ordinal variables
            ov.int <- ov.names[[g]]
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
        if(!conditional.x && (nx <- length(ov.names.x[[g]])) > 0L) {
            # fix variances
            exo.idx <- which(rhs %in% ov.names.x[[g]] &
                             lhs %in% ov.names.x[[g]] &
                             op == "~~" & group == g)
            if(fixed.x) {
                exo[exo.idx] <- 1L
                free[exo.idx] <- 0L
            }

            # fix means
            exo.idx <- which(lhs %in% ov.names.x[[g]] &
                             op == "~1" & group == g)
            if(fixed.x) {
                exo[exo.idx] <- 1L
                free[exo.idx] <- 0L
            }

            # add covariances
            pstar <- nx*(nx-1)/2
            if(pstar > 0L) { # only if more than 1 variable
                tmp <- utils::combn(ov.names.x[[g]], 2)
                lhs <- c(lhs, tmp[1,]) # to fill upper.tri
                 op <- c(op,   rep("~~", pstar))
                rhs <- c(rhs, tmp[2,])
                group <- c(group, rep(g,  pstar))
                if(fixed.x) {
                    free  <- c(free,  rep(0L, pstar))
                    exo   <- c(exo,   rep(1L, pstar))
                } else {
                    free  <- c(free,  rep(1L, pstar))
                    exo   <- c(exo,   rep(0L, pstar))
                }

                # starting values
                if(!is.null(sample.cov)) {
                    rhs.idx <- match(tmp[1,], ov.names[[g]])
                    lhs.idx <- match(tmp[2,], ov.names[[g]])
                    ustart <- c(ustart, 
                                sample.cov[[g]][ cbind(rhs.idx, lhs.idx) ])
                } else {
                    ustart <- c(ustart, rep(as.numeric(NA), pstar))
                }
            }
        }

        if(conditional.x && (nx <- length(ov.names.x[[g]])) > 0L) {
            # add regressions
            lhs <- c(lhs, rep("dummy", nx))
             op <- c( op, rep("~", nx))
            rhs <- c(rhs, ov.names.x[[g]])

            # add dummy latent
            lhs <- c(lhs,"dummy"); op <- c(op, "=~"); rhs <- c(rhs, "dummy")
            lhs <- c(lhs,"dummy"); op <- c(op, "~~"); rhs <- c(rhs, "dummy") 

            exo <- c(exo,    rep(1L, nx)); exo <- c(exo,   c(0L,0L))
          group <- c(group,  rep(g,  nx + 2L))
           free <- c(free,   rep(0L, nx + 2L))
         ustart <- c(ustart, rep(0,  nx + 2L))

            if(meanstructure) {
                lhs <- c(lhs,"dummy"); op <- c(op, "~1"); rhs <- c(rhs, "")

                exo    <- c(exo,    0L)
                group  <- c(group,  g)
                free   <- c(free,   0L)
                ustart <- c(ustart, 0)
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
                        label       = rep("",  length(lhs))
                        #eq.id       = rep(0L,  length(lhs)),
                        #unco        = free
                   )
    LIST

}
