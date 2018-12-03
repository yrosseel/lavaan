# twostep SEM
# experimental version (for simulation purposes only)
# YR -- 9 July 2017
#
# step 1: fit all measurement blocks (here: one per latent variable)
# step 2: fit structural part, keeping all measurement model parameters fixed
#
# standard errors structural part are computed as in PML (Gong & Samaniego)
# see also Bakk, Oberski & Vermunt (2014)


twostep <- function(model      = NULL,
                    cmd        = "sem",
                    ...,
                    add.class = TRUE,
                    output    = "lavaan") { # lavaan, PT or list

    # dot dot dot
    dotdotdot <- list(...)

    # STEP 0: process full model, without fitting
    dotdotdot0 <- dotdotdot
    dotdotdot0$do.fit <- NULL
    dotdotdot0$se     <- "none"
    dotdotdot0$test   <- "none"

    # check for arguments that we do not want?
    # TODO

    # initial processing of the model, no fitting
    FIT <- do.call(cmd,
                   args =  c(list(model  = model,
                                  do.fit = FALSE), dotdotdot0) )

    # restore options
    lavoptions <- lavInspect(FIT, "options")
    lavoptions$do.fit <- TRUE

    if(!is.null(dotdotdot$se)) {
        lavoptions$se   <- dotdotdot$se
    } else {
        lavoptions$se   <- "standard"
    }

    if(!is.null(dotdotdot$test)) {
        lavoptions$test <- dotdotdot$test
    } else {
        lavoptions$test <- "standard"
    }


    # extract other info
    ngroups    <- lavInspect(FIT, "ngroups")
    lavpta     <- FIT@pta



    # first, check if we have lv's
    lv.names <- unique(unlist(FIT@pta$vnames$lv.regular))
    if(length(lv.names) == 0L) {
        stop("lavaan ERROR: model does not contain any latent variables")
    }
    nfac     <- length(lv.names)

    # extract parameter table
    PT <- FIT@ParTable

    # total number of free parameters
    npar <- lav_partable_npar(PT)
    if(npar < 1L) {
        stop("lavaan ERROR: model does not contain any free parameters")
    }


    # make a copy
    PT.orig <- PT

    # est equals ustart by default
    PT$est <- PT$ustart

    # clear se values
    PT$se <- rep(as.numeric(NA), length(PT$lhs))
    PT$se[ PT$free == 0L & !is.na(PT$ustart) ] <- 0.0


    # STEP 1: fit measurent model for each meaesurement block
    #         FIXME: measurement block == single lv for now!!
    MM <- vector("list", length = nfac)
    MM.INFO <- matrix(0, npar, npar)
    step1.idx <- integer(0L)
    for(f in seq_len(nfac)) {

        # which parameters are related to this measurement block?
        mm.idx <- lav_partable_subset_measurement_model(PT = PT,
                      lavpta = lavpta, lv.names = lv.names[f], idx.only = TRUE)

        if(length(mm.idx) == 0L) {
            # empty measurement block (single-indicator lv?)
            next
        }

        # adapt parameter table:
        # - only parameters related to this measurement block are 'free'
        # - everything else is set to zero (except variances, which are
        #   are set to 1)
        PTM <- PT
        PTM$free[ !seq_len(length(PT$lhs)) %in% mm.idx &
                   PT$free > 0L ] <- 0L
        PTM$free[ PTM$free > 0L ] <- seq_len( sum(PTM$free > 0L) )

        # set all other non-fixed ustart values to 0
        PTM$ustart[ PTM$free == 0L &
                    is.na(PTM$ustart) ] <- 0
        # set all other non-fixed ustart variance values to 1
        PTM$ustart[ PTM$free == 0L & PTM$op == "~~" & PTM$lhs == PTM$rhs &
                    PTM$ustart == 0 ] <- 1

        # fit this measurement block, store the fitted object in the MM list
        MM[[f]] <- lavaan::lavaan(model = PTM, ...) ## FIXME: reuse slots!

        # fill in point estimates measurement block
        PT$est[ seq_len(length(PT$lhs)) %in% mm.idx &
                PT$free > 0L ] <- MM[[f]]@ParTable$est[ PTM$free > 0L ]

        # fill in standard errors measurement block
        PT$se[ seq_len(length(PT$lhs)) %in% mm.idx &
               PT$free > 0L ] <- MM[[f]]@ParTable$se[ PTM$free > 0L ]

        # compute `total' information for this measurement block
        mm.info <- lavTech(MM[[f]], "information") * lavTech(FIT, "nobs")

        # fill in `total' information matrix
        par.idx <- PT.orig$free[ seq_len(length(PT$lhs)) %in% mm.idx &
                                 PT$free > 0L ]
        MM.INFO[par.idx, par.idx] <- mm.info

        # store indices in step1.idx
        step1.idx <- c(step1.idx, par.idx)
    }

    # the measurement model parameters now become fixed ustart values
    PT$ustart <- PT$est




    # STEP 2: fit structural part only, holding the measurement model fixed
    reg.idx <- lav_partable_subset_structural_model(PT = PT,
                      lavpta = lavpta, idx.only = TRUE)

    # remove 'exogenous' factor variances (if any) from reg.idx
    lv.names.x <- lv.names[ lv.names %in% unlist(lavpta$vnames$eqs.x) ]
    if(lavoptions$fixed.x && length(lv.names.x) > 0L) {
        var.idx <- which(PT$lhs %in% lv.names.x &
                         PT$op == "~~" &
                         PT$lhs == PT$rhs)
        rm.idx <- which(reg.idx %in% var.idx)
        if(length(rm.idx) > 0L) {
            reg.idx <- reg.idx[ -rm.idx ]
        }
    }

    # adapt parameter table for structural part
    PTS <- PT
    PTS$free[ !seq_len(length(PT$lhs)) %in% reg.idx &
                   PT$free > 0L ] <- 0L
    PTS$free[ PTS$free > 0L ] <- seq_len( sum(PTS$free > 0L) )

    # set 'ustart' values for free STRUC parameter to NA
    PTS$ustart[ PTS$free > 0L ] <- as.numeric(NA)

    # estimate structural part
    STRUC <- lavaan::lavaan(model = PTS, ...)  ### FIXME: reuse slots

    # fill in point estimates structural part
    PT$est[ seq_len(length(PT$lhs)) %in% reg.idx &
            PT$free > 0L ] <- STRUC@ParTable$est[ PTS$free > 0L ]

    # construct JOINT model
    JOINT <- lavaan::lavaan(PT, ..., optim.method = "none",
                            check.gradient = FALSE, # often not ok
                            se = "external")

    # TOTAL information
    INFO <- lavInspect(JOINT, "information") * nobs(FIT)

    #step1.idx <- PT$free[ !seq_len(length(PT$lhs)) %in% reg.idx &
    #                       PT$free > 0L ]
    step2.idx <- PT$free[  seq_len(length(PT$lhs)) %in% reg.idx &
                           PT$free > 0L ]

    I.12 <- INFO[step1.idx, step2.idx]
    I.22 <- INFO[step2.idx, step2.idx]
    I.21 <- INFO[step2.idx, step1.idx]

    # compute Sigma.11
    # note: step1.idx and step2.idx have some overlap, but removing
    #       the step2.idx parameters from step1.idx BEFORE we take
    #       the inverse, results in std.errors that are lower than ML,
    #       in addition, fixed.x = TRUE/FALSE give different results
    #       for SE of beta
    #
    # the solution seems to be to keep ALL step1.idx elements in MM.INFO
    # for the inversion, and then set rows/cols of Sigma.11 to zero if
    # the overlap with reg.idx (YR, 28 nov 2018; lavaan 0.6-4)

    MM.INFO <- MM.INFO[step1.idx, step1.idx, drop = FALSE]
    Sigma.11 <- solve(MM.INFO)

    # overlap? set corresponding rows/cols of Sigma.11 to zero
    both.idx <- which(step1.idx %in% step2.idx)
    if(length(both.idx) > 0L) {
        Sigma.11[both.idx,] <- 0
        Sigma.11[,both.idx] <- 0
    }

    # V2
    I.22.inv <- solve(I.22)
    V2 <- I.22.inv

    # V1
    V1 <- I.22.inv %*% I.21 %*% Sigma.11 %*% I.12 %*% I.22.inv

    # V for second step
    V <- V2 + V1

    # fill in standard errors step 2
    PT$se[ seq_len(length(PT$lhs)) %in% reg.idx &
               PT$free > 0L ] <- sqrt( diag(V) )

    if(output == "lavaan") {
        FINAL <- lavaan::lavaan(PT, ..., optim.method = "none",
                            se = "external", check.gradient = FALSE)
        return(FINAL)
    } else if(output == "PT") {
        # for pretty printing only
        if(add.class) {
            class(PT) <- c("lavaan.data.frame", "data.frame")
        }
        return(PT)
    } else {
        return( list(MM = MM, STRUC = STRUC, JOINT = JOINT,
                     V = V, V1 = V1, V2 = V2, Sigma.11 = Sigma.11,
                     MM.INFO = MM.INFO, PT = PT) )
    }
}

