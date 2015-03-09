# YR - 26 Nov 2013: generate partable for the unrestricted model
lav_partable_unrestricted <- function(lavobject      = NULL,
                                      lavdata        = NULL,
                                      lavoptions     = NULL,
                                      lavsamplestats = NULL,
                                      # optional user-provided sample stats
                                      sample.cov     = NULL,
                                      sample.mean    = NULL,
                                      sample.th      = NULL,
                                      sample.th.idx  = NULL) {

    # grab everything from lavaan lavobject
    if(!is.null(lavobject)) {
        stopifnot(inherits(lavobject, "lavaan"))

        lavdata <- lavobject@Data
        lavoptions <- lavobject@Options
        lavsamplestats <- lavobject@SampleStats
    }

    # if user-based moments are given, use these
    if(is.null(sample.cov) && !is.null(lavsamplestats)) {
        sample.cov <- lavsamplestats@cov
    }
    if(is.null(sample.mean) && !is.null(lavsamplestats)) {
        sample.mean <- lavsamplestats@mean
    }
    if(is.null(sample.th) && !is.null(lavsamplestats)) {
         sample.th <- lavsamplestats@th
    }
    if(is.null(sample.th.idx) && !is.null(lavsamplestats)) {
         sample.th.idx <- lavsamplestats@th.idx
    }

    ov.names      = lavdata@ov.names
    ov            = lavdata@ov
    ov.names.x    = lavdata@ov.names.x
    meanstructure = lavoptions$meanstructure

    # what with fixed.x?
    if(lavoptions$mimic %in% c("lavaan", "Mplus")) {
        fixed.x = lavoptions$fixed.x
    } else if(lavoptions$mimic == "EQS") {
        # always ignore fixed.x
        ov.names.x = NULL
        fixed.x = FALSE
    } else if(lavoptions$mimic == "LISREL") {
        # always ignore fixed.x??? CHECKME!!
        ov.names.x = NULL
        fixed.x = FALSE
    }

    ngroups <- length(ov.names)
    ov.names.nox <- lapply(as.list(1:ngroups), function(g) 
                    ov.names[[g]][ !ov.names[[g]] %in% ov.names.x[[g]] ])



    lhs <- rhs <- op <- character(0)
    group <- free <- exo <- integer(0)
    ustart <- numeric(0)

    categorical <- any(ov$type == "ordered")
    if(categorical) {
        OV <- ov.names.nox
    } else {
        OV <- ov.names
    }
    

    for(g in 1:ngroups) {

        # a) VARIANCES (all ov's, if !categorical, also exo's)
        nvar  <- length(OV[[g]])
        lhs   <- c(lhs, OV[[g]])
         op   <- c(op, rep("~~", nvar))
        rhs   <- c(rhs, OV[[g]])
        group <- c(group, rep(g,  nvar))
        free  <- c(free,  rep(1L, nvar))
        exo   <- c(exo,   rep(0L, nvar))

        # starting values -- variances
        if(!is.null(sample.cov)) {
            sample.var.idx <- match(OV[[g]], ov.names[[g]])
            ustart <- c(ustart, diag(sample.cov[[g]])[sample.var.idx])
        } else {
            ustart <- c(ustart, rep(as.numeric(NA), nvar))
        }

        # COVARIANCES!
        pstar <- nvar*(nvar-1)/2
        if(pstar > 0L) { # only if more than 1 variable
            tmp <- utils::combn(OV[[g]], 2)
            lhs <- c(lhs, tmp[1,]) # to fill upper.tri
             op <- c(op,   rep("~~", pstar))
            rhs <- c(rhs, tmp[2,])
            group <- c(group, rep(g,  pstar))
            free  <- c(free,  rep(1L, pstar))
            exo   <- c(exo,   rep(0L, pstar))
        }

        # starting values -- variances
        if(!is.null(sample.cov)) {
            sample.var.idx <- match(OV[[g]], ov.names[[g]])
            COV <-  sample.cov[[g]][sample.var.idx, sample.var.idx]
            ### CHECK ME!!! upper tri??
            ustart <- c(ustart, COV[lower.tri(COV, diag=FALSE)])
        } else {
            ustart <- c(ustart, rep(as.numeric(NA), pstar))
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
            }
        }

        # meanstructure?
        if(meanstructure) {
            # auto-remove ordinal variables
            ov.int <- OV[[g]]
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

        # categorical? insert means as fixed-to-zero parameters
        # since 0.5-17
        if(categorical) {
            ov.int <- OV[[g]]
            idx <- which(ov.int %in% ord.names)
            ov.int <- ov.int[idx]

            nel <- length(ov.int)
            lhs   <- c(lhs, ov.int)
             op   <- c(op, rep("~1", nel))
            rhs   <- c(rhs, rep("", nel))
            group <- c(group, rep(g,  nel))
            free  <- c(free,  rep(0L, nel))
            exo   <- c(exo,   rep(0L, nel))
           ustart <- c(ustart, rep(0L, nel))
        }
 

        # fixed.x exogenous variables?
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
                        label       = rep("",  length(lhs))
                        #eq.id       = rep(0L,  length(lhs)),
                        #unco        = free
                   )
    LIST

}
