# YR - 26 Nov 2013: generate partable for the unrestricted model
lav_partable_unrestricted <- function(lavobject      = NULL,
                                      lavdata        = NULL,
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
    }

    # if user-based moments are given, use these
    if(is.null(sample.cov) && !is.null(lavsamplestats)) {
        sample.cov <- lavsamplestats@cov
    }
    if(is.null(sample.mean) && !is.null(lavsamplestats)) {
        sample.mean <- lavsamplestats@mean
    }
    if(is.null(sample.slopes) && !is.null(lavsamplestats)) {
        sample.slopes <- lavsamplestats@res.slopes
    }
    if(is.null(sample.th) && !is.null(lavsamplestats)) {
         sample.th <- lavsamplestats@th
    }
    if(is.null(sample.th.idx) && !is.null(lavsamplestats)) {
         sample.th.idx <- lavsamplestats@th.idx
    }

    ov.names      <- lavdata@ov.names
    ov            <- lavdata@ov
    ov.names.x    <- lavdata@ov.names.x
    meanstructure <- lavoptions$meanstructure
    categorical   <- any(ov$type == "ordered")
    conditional.x <- lavoptions$conditional.x
    ngroups <- length(ov.names)

    # what with fixed.x? 
    # - does not really matter; fit will be saturated any way
    # - fixed.x = TRUE may avoid convergence issues with non-numeric 
    #             x-covariates
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

    if(conditional.x) {
        ov.names.nox <- lapply(seq_len(ngroups), function(g)
                ov.names[[g]][ !ov.names[[g]] %in% ov.names.x[[g]] ])
    }

    lhs <- rhs <- op <- character(0)
    group <- free <- exo <- integer(0)
    ustart <- numeric(0)

    for(g in 1:ngroups) {

        # a) VARIANCES (all ov's, if !conditional.x, also exo's)
        nvar  <- length(ov.names[[g]])

        lhs   <- c(lhs, ov.names[[g]])
         op   <- c(op, rep("~~", nvar))
        rhs   <- c(rhs, ov.names[[g]])
        group <- c(group, rep(g,  nvar))
        free  <- c(free,  rep(1L, nvar))
        exo   <- c(exo,   rep(0L, nvar))

        # starting values -- variances
        if(!is.null(sample.cov)) {
            ustart <- c(ustart, diag(sample.cov[[g]]))
        } else {
            ustart <- c(ustart, rep(as.numeric(NA), nvar))
        }

        # COVARIANCES!
        pstar <- nvar*(nvar-1)/2
        if(pstar > 0L) { # only if more than 1 variable
            tmp <- utils::combn(ov.names[[g]], 2)
            lhs <- c(lhs, tmp[1,]) # to fill upper.tri
             op <- c(op,   rep("~~", pstar))
            rhs <- c(rhs, tmp[2,])
            group <- c(group, rep(g,  pstar))
            free  <- c(free,  rep(1L, pstar))
            exo   <- c(exo,   rep(0L, pstar))
        }

        # starting values -- covariances
        if(!is.null(sample.cov)) {
            ustart <- c(ustart, lav_matrix_vech(sample.cov[[g]], 
                                                diagonal = FALSE))
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

        # categorical? insert means as fixed-to-zero parameters
        # since 0.5-17
        if(categorical) {
            ov.int <- ov.names[[g]]
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
        if(!conditional.x && 
           fixed.x && (nx <- length(ov.names.x[[g]])) > 0L) {
            # fix variances/covariances 
            exo.idx <- which(rhs %in% ov.names.x[[g]] &
                             lhs %in% ov.names.x[[g]] &
                             op == "~~" & group == g)
            exo[exo.idx] <- 1L
            free[exo.idx] <- 0L

            # fix means
            exo.idx <- which(rhs %in% ov.names.x[[g]] &
                             op == "~1" & group == g)
            exo[exo.idx] <- 1L
            free[exo.idx] <- 0L
        }

        # conditional.x?
        if(conditional.x && (nx <- length(ov.names.x[[g]])) > 0L) {
            nnox <- length(ov.names.nox[[g]])
            nel  <- nnox * nx

            lhs <- c(lhs, rep(ov.names.nox[[g]], times = nx))
             op <- c(op,  rep("~", nel))
            rhs <- c(rhs, rep(ov.names.x[[g]], each = nnox))
          group <- c(group, rep(g,  nel))
           free <- c(free,  rep(1L, nel))
            exo <- c(exo,   rep(1L, nel))

            # starting values -- slopes
            if(!is.null(sample.slopes)) {
                ustart <- c(ustart, lav_matrix_vec(sample.slopes[[g]]))
            } else {
                ustart <- c(ustart, rep(as.numeric(NA), nel))
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
                        #mod.idx     = rep(0L,  length(lhs)),
                        free        = free,
                        ustart      = ustart,
                        exo         = exo #,
                        #label       = rep("",  length(lhs))
                        #eq.id       = rep(0L,  length(lhs)),
                        #unco        = free
                   )


    LIST

}
