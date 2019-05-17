# twostep SEM
# experimental version (for simulation purposes only)
# YR -- 9 July 2017
#
# step 1: fit all measurement blocks (here: one per latent variable)
# step 2: fit structural part, keeping all measurement model parameters fixed
#
# standard errors structural part are computed as in PML (Gong & Samaniego)
# see also Bakk, Oberski & Vermunt (2014)

# YR -- 16 May 2019: - add mm.list
#                    - more overlap with xxx_sam()


twostep <- function(model      = NULL,
                    data       = NULL,
                    cmd        = "sem",
                    mm.list    = NULL,
                    mm.args    = list(),
                    struc.args = list(),
                    ...,         # global options
                    output     = "global") { # list, structural, global

    # check arguments

    # output
    output <- tolower(output)
    if(output == "list") {
        # nothing to do
    } else if(output %in% c("struc", "structural", "structural.part")) {
        output <- "structural"
        stop("lavaan ERROR: output=structural not ready yet.")
    } else if(output %in% c("joint", "full", "global", "total")) {
        output <- "global"
    } else {
        stop("lavaan ERROR: output should be one of list, structural or global.")
    }

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
                                  data   = data,
                                  do.fit = FALSE), dotdotdot0) )
    lavoptions <- lavInspect(FIT, "options")

    # restore options
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


    # what have we learned?
    ngroups    <- lavInspect(FIT, "ngroups")
    lavpta     <- FIT@pta
    PT         <- FIT@ParTable


    # any `regular' latent variables? (across groups!)
    LV.names <- unique(unlist(FIT@pta$vnames$lv.regular))
    OV.names <- unique(unlist(FIT@pta$vnames$ov))

    # check for higher-order factors
    LV.IND.names <- unique(unlist(FIT@pta$vnames$lv.ind))
    if(length(LV.IND.names) > 0L) {
        ind.idx <- match(LV.IND.names, LV.names)
        LV.names <- LV.names[-ind.idx]
    }

    # do we have at least 1 'regular' (measured) latent variable?
    if(length(LV.names) == 0L) {
        stop("lavaan ERROR: model does not contain any (measured) latent variables; use sem() instead")
    }
    nfac <- length(LV.names)


    # total number of free parameters
    npar <- lav_partable_npar(PT)
    if(npar < 1L) {
        stop("lavaan ERROR: model does not contain any free parameters")
    }

    # check parameter table
    PT$est <- PT$se <- NULL
    # est equals start by default (except exo values)
    PT$est <- PT$start
    if(any(PT$exo > 0L)) {
        PT$est[PT$exo > 0L] <- PT$start[PT$exo > 0L]
    }

    # clear se values (needed here?)
    PT$se <- rep(as.numeric(NA), length(PT$lhs))
    PT$se[ PT$free == 0L & !is.na(PT$ustart) ] <- 0.0


    # STEP 1: fit each measurement model (block)
    # how many measurement models?
    if(!is.null(mm.list)) {
        nMMblocks <- length(mm.list)
        # check each measurement block
        for(b in seq_len(nMMblocks)) {
            if(!all(mm.list[[b]] %in% LV.names)) {
              stop("lavaan ERROR: mm.list contains unknown latent variable(s):",
                paste( mm.list[[b]][ !mm.list[[b]] %in% LV.names ], sep = " "),
                "\n")
            }
        }
    } else {
        # TODO: here comes the automatic 'detection' of linked
        #       measurement models
        #
        # for now we take a single latent variable per measurement model block
        mm.list <- as.list(LV.names)
        nMMblocks <- length(mm.list)
    }

    # adjust options for measurement models
    dotdotdot.mm <- dotdotdot
    #dotdotdot.mm$se <- "none"
    dotdotdot.mm$test <- "none" # global view!
    dotdotdot.mm$debug <- FALSE
    dotdotdot.mm$verbose <- FALSE
    dotdotdot.mm$check.post <- FALSE # neg lv variances may be overriden

    # override with mm.args
    dotdotdot.mm <- modifyList(dotdotdot.mm, mm.args)



    # STEP 1: fit measurent model for each measurement block
    #         FIXME: measurement block == single lv for now!!
    MM.FIT <- vector("list", nMMblocks)

    #MM.INFO <- matrix(0, npar, npar)
    Sigma.11 <- matrix(0, npar, npar)
    step1.idx <- integer(0L)

    for(mm in seq_len(nMMblocks)) {

        # which parameters are related to this measurement block?
        mm.idx <- lav_partable_subset_measurement_model(PT = PT,
                      lavpta = lavpta, lv.names = mm.list[[mm]],
                      idx.only = TRUE)

        if(length(mm.idx) == 0L) {
            # empty measurement block (single-indicator lv?)
            next
        }

        # adapt parameter table:
        # - only parameters related to this measurement block are 'free'
        # - everything else is set to zero (except variances, which are
        #   are set to 1)
        # - if multiple latent variables in this measurement block, we
        #   must ADD and free the covariances among these latent variables
        # - remove non-needed constraints
        PTM <- PT
        if(length(mm.list[[mm]]) > 1L) {
            PTM <- lav_partable_add_lv_cov(PT = PT,
                          lavpta = lavpta, lv.names = mm.list[[mm]])
        }
        con.idx <- which(PTM$op %in% c("==","<",">"))
        if(length(con.idx) > 0L) {
            needed.idx <- which(con.idx %in% mm.idx)
            if(length(needed.idx) > 0L) {
                con.idx <- con.idx[-needed.idx]
            }
            if(length(con.idx) > 0L) {
                PTM <- as.data.frame(PTM, stringsAsFactors = FALSE)
                PTM <- PTM[-con.idx, ]
            }
        }
        PTM$est <- NULL
        PTM$se <- NULL

        PTM$free[ !seq_len(length(PTM$lhs)) %in% mm.idx &
                   PTM$free > 0L ] <- 0L
        PTM$free[ PTM$user == 3L ] <- 1L
        PTM$free[ PTM$free > 0L ] <- seq_len( sum(PTM$free > 0L) )

        # set all other non-fixed ustart values to 0
        PTM$ustart[ PTM$free == 0L & is.na(PTM$ustart) ] <- 0

        # set all other non-fixed ustart variance values to 1
        PTM$ustart[ PTM$free == 0L & PTM$op == "~~" & PTM$lhs == PTM$rhs &
                    PTM$ustart == 0 ] <- 1

        # fit this measurement block, store the fitted object in the MM list
        fit.mm.block <- do.call("lavaan",
                                args =  c(list(model  = PTM,
                                               data   = data), dotdotdot.mm) )
        # check convergence
        if(!lavInspect(fit.mm.block, "converged")) {
            # warning for now
            warning("lavaan WARNING: measurement model for ",
                 paste(mm.list[[mm]], collapse = " "), " did not converge.")
        }

        # store fitted measurement model
        MM.FIT[[mm]] <- fit.mm.block

        # fill in point estimates measurement block
        PTM <- MM.FIT[[mm]]@ParTable
        PT$est[ seq_len(length(PT$lhs)) %in% mm.idx & PT$free > 0L ] <-
            PTM$est[ PTM$free > 0L & PTM$user != 3L]

        # fill in standard errors measurement block
        PT$se[ seq_len(length(PT$lhs)) %in% mm.idx & PT$free > 0L ] <-
            PTM$se[ PTM$free > 0L & PTM$user != 3L]

        # compute `total' information for this measurement block
        #mm.info <- lavTech(MM.FIT[[mm]], "information") * lavTech(FIT, "nobs")
        sigma.11 <- MM.FIT[[mm]]@vcov$vcov

        # fill in `total' information matrix
        par.idx <- PT$free[ seq_len(length(PT$lhs)) %in% mm.idx & PT$free > 0L ]
        keep.idx <- PTM$free[ PTM$free > 0 & PTM$user != 3L ]
        #MM.INFO[par.idx, par.idx] <- mm.info[keep.idx, keep.idx, drop = FALSE]
        Sigma.11[par.idx, par.idx] <- sigma.11[keep.idx, keep.idx, drop = FALSE]

        # store indices in step1.idx
        step1.idx <- c(step1.idx, par.idx)

    } # measurement block
    Sigma.11 <- Sigma.11[step1.idx, step1.idx, drop = FALSE]

    # do we have any parameters left?
    if(length(step1.idx) >= npar) {
        warning("lavaan WARNING: ",
                "no free parameters left for structural part.\n",
                "        Returning measurement part only.")
        if(output == "list") {
            return(list(MM.FIT = MM.FIT, Sigma.11 = Sigma.11, PT = PT))
        } else {
            if(nMMblocks == 1L) {
                return(MM.FIT[[1]])
            } else {
                return(MM.FIT)
            }
        }
    }




    # STEP 2: fit structural part only, holding the measurement model fixed

    # the measurement model parameters now become fixed ustart values
    PT$ustart[PT$free > 0] <- PT$est[PT$free > 0]

    reg.idx <- lav_partable_subset_structural_model(PT = PT,
                      lavpta = lavpta, idx.only = TRUE)

    # remove 'exogenous' factor variances (if any) from reg.idx
    lv.names.x <- LV.names[ LV.names %in% unlist(lavpta$vnames$eqs.x)  &
                           !LV.names %in% unlist(lavpta$vnames$eqs.y) ]
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

    # remove constraints we don't need
    con.idx <- which(PTS$op %in% c("==","<",">"))
    if(length(con.idx) > 0L) {
        needed.idx <- which(con.idx %in% reg.idx)
        if(length(needed.idx) > 0L) {
            con.idx <- con.idx[-needed.idx]
        }
        if(length(con.idx) > 0L) {
            PTS <- as.data.frame(PTS, stringsAsFactors = FALSE)
            PTS <- PTS[-con.idx, ]
        }
    }
    PTS$est <- NULL
    PTS$se <- NULL

    PTS$free[ !seq_len(length(PTS$lhs)) %in% reg.idx & PTS$free > 0L ] <- 0L
    PTS$free[ PTS$free > 0L ] <- seq_len( sum(PTS$free > 0L) )

    # set 'ustart' values for free FIT.PA parameter to NA
    PTS$ustart[ PTS$free > 0L ] <- as.numeric(NA)

    # adjust options
    lavoptions.PA <- lavoptions
    lavoptions.PA <- modifyList(lavoptions.PA, struc.args)

    # override, not matter what
    lavoptions.PA$do.fit <- TRUE

    # estimate structural part
    FIT.PA <- lavaan::lavaan(model = PTS,
                             slotData = FIT@Data,
                             slotSampleStats = FIT@SampleStats,
                             slotOptions = lavoptions.PA)

    # fill in point estimates structural part
    PTS <- FIT.PA@ParTable
    PT$est[ seq_len(length(PT$lhs)) %in% reg.idx & PT$free > 0L ] <-
        PTS$est[ PTS$free > 0L ]
    step2.idx <- PT$free[  seq_len(length(PT$lhs)) %in% reg.idx &
                           PT$free > 0L ]


    # compute standard errors
    # current approach:
    # - create 'global' model, to get the 'joint' information matrix
    # - partition information matrix (step 1, step 2)
    # - apply two-step correction for second step
    # - 'insert' these corrected SEs (and vcov) for the parameters of FIT.PA

    lavoptions.joint <- lavoptions
    lavoptions.joint$optim.method = "none"
    lavoptions.joint$check.gradient = FALSE
    #lavoptions.joint$se = "none"
    #lavoptions.joint$test = "none"
    JOINT <- lavaan::lavaan(PT, slotOptions = lavoptions.joint,
                            slotSampleStats = FIT@SampleStats,
                            slotData = FIT@Data)

    # TOTAL information
    INFO <- lavInspect(JOINT, "information") * nobs(FIT)
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

    #MM.INFO <- MM.INFO[step1.idx, step1.idx, drop = FALSE]
    #Sigma.11 <- solve(MM.INFO)

    # overlap? set corresponding rows/cols of Sigma.11 to zero
    both.idx <- which(step1.idx %in% step2.idx)
    if(length(both.idx) > 0L) {
        Sigma.11[both.idx,] <- 0
        Sigma.11[,both.idx] <- 0
    }

    # V2
    I.22.inv <- solve(I.22)

    # method below has the advantage that we can use a 'robust' vcov
    # but does not work if we have equality constraints in the MM
    # -> D will be singular
    #A <- JOINT@vcov$vcov[ step2.idx,  step2.idx]
    #B <- JOINT@vcov$vcov[ step2.idx, -step2.idx]
    #C <- JOINT@vcov$vcov[-step2.idx,  step2.idx]
    #D <- JOINT@vcov$vcov[-step2.idx, -step2.idx]
    #I.22.inv <- A - B %*% solve(D) %*% C
    V2 <- I.22.inv

    # V1
    V1 <- I.22.inv %*% I.21 %*% Sigma.11 %*% I.12 %*% I.22.inv

    # V for second step
    VCOV <- V2 + V1

    # fill in standard errors step 2
    PT$se[ seq_len(length(PT$lhs)) %in% reg.idx &
               PT$free > 0L ] <- sqrt( diag(VCOV) )

    if(output == "global") {
        FINAL <- JOINT
        FINAL@ParTable <- PT
        if(FINAL@Options$se == "none") {
        } else {
            FINAL@Options$se <- "twostep"
            FINAL@vcov$vcov[step2.idx, step2.idx] <- VCOV
        }
        res <- FINAL
    } else {
        res <- list(MM.FIT = MM.FIT, FIT.PA = FIT.PA, JOINT = JOINT,
                    VCOV = VCOV, V1 = V1, V2 = V2, Sigma.11 = Sigma.11,
                    MM.INFO = MM.INFO, PT = PT)
    }

    res
}

