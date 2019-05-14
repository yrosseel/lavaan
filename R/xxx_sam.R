# SAM: a Structural After Measurement approach
#
# Yves Rosseel & Wen-Wei Loh, Feb-May 2019

# alternative for FSR+Croon
#  - but no need to compute factor scores or corrections

# restrictions
#  - only if LAMBDA is of full column rank (eg no SRM, no bi-factor, no MTMM)
#  - if multiple groups: each group has the same set of latent variables!

# YR 12 May 2019 - first version
#

sam <- function(model      = NULL,
                data       = NULL,
                cmd        = "sem",
                mm.list    = NULL,
                mm.args    = list(se = "standard", test = "standard"),
                M.method   = "GLS",
                struc.args = list(se = "twostep", test = "standard"),
                ...,         # global options
                output     = "lavaan") {

    # check input
    if(!M.method %in% c("GLS", "ML", "ULS")) {
        stop("lavaan ERROR: M.method should be one of GLS, ML or ULS.")
    }

    # dot dot dot
    dotdotdot <- list(...)

    # STEP 0: process full model, without fitting
    dotdotdot0 <- dotdotdot
    dotdotdot0$do.fit <- NULL
    dotdotdot0$se     <- "none"      # to avoid warning about missing="listwise"
    dotdotdot0$test   <- "none"      # to avoid warning about missing="listwise"
    if(!is.null(dotdotdot0$auto.cov.lv.x) &&
       !dotdotdot0$auto.cov.lv.x) {
        warning("lavaan WARNING: sam forces auto.cov.lv.x = TRUE")
    }
    dotdotdot0$auto.cov.lv.x <- TRUE # force allow correlated exogenous factors

    # initial processing of the model, no fitting
    FIT <- do.call(cmd,
                   args =  c(list(model  = model,
                                  data   = data,
                                  do.fit = FALSE), dotdotdot0) )
    lavoptions <- lavInspect(FIT, "options")

    # what have we learned?
    ngroups     <- lavInspect(FIT, "ngroups")
    lavpta      <- FIT@pta
    lavpartable <- FIT@ParTable

    # if missing = "listwise", make data complete, to avoid different
    # datasets per measurement block (NEEDED?)
    if(lavoptions$missing == "listwise") {
        # FIXME: make this work for multiple groups!!
        OV <- unique(unlist(lavpta$vnames$ov))
        data <- na.omit(data[,OV])
    }

    # any `regular' latent variables? (across groups!)
    LV.names <- unique(unlist(FIT@pta$vnames$lv.regular))
    OV.names <- unique(unlist(FIT@pta$vnames$ov))

    # check for higher-order factors
    good.idx <- logical( length(LV.names) )
    for(f in seq_len(length(LV.names))) {
        # check the indicators
        FAC <- LV.names[f]
        IND <- lavpartable$rhs[ lavpartable$lhs == FAC &
                                lavpartable$op  == "=~" ]
        if(all(IND %in% OV.names)) {
            good.idx[f] <- TRUE
        }
        # FIXME: check for mixed lv/ov indicators
    }
    LV.names <- LV.names[ good.idx ]

    # do we have at least 1 'regular' (measured) latent variable?
    if(length(LV.names) == 0L) {
        stop("lavaan ERROR: model does not contain any (measured) latent variables; use sem() instead")
    }
    nfac <- length(LV.names)

    # check parameter table
    PT <- parTable(FIT)
    PT.orig <- PT
    PT$est <- PT$se <- NULL
    # est equals ustart by default (except exo values)
    PT$est <- PT$ustart
    if(any(PT$exo > 0L)) {
        PT$est[PT$exo > 0L] <- PT$start[PT$exo > 0L]
    }

    # clear se values (needed here?)
    PT$se <- rep(as.numeric(NA), length(PT$lhs))
    PT$se[ PT$free == 0L & !is.na(PT$ustart) ] <- 0.0
    # total number of free parameters
    npar <- lav_partable_npar(PT)
    if(npar < 1L) {
        stop("lavaan ERROR: model does not contain any free parameters")
    }


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
    dotdotdot.mm$se <- "none"
    dotdotdot.mm$test <- "none"
    dotdotdot.mm$debug <- FALSE
    dotdotdot.mm$verbose <- FALSE

    # override with mm.args
    dotdotdot.mm <- modifyList(dotdotdot.mm, mm.args)

    # we assume the same number/names of lv's per group!!!
    MM.FIT <- vector("list", nMMblocks)         # fitted object
    LAMBDA.list <- vector("list", nMMblocks)
    THETA.list  <- vector("list", nMMblocks)
    NU.list     <- vector("list", nMMblocks)
    LV.idx.list <- vector("list", nMMblocks)
    OV.idx.list <- vector("list", nMMblocks)

    # for joint model later
    MM.INFO <- matrix(0, npar, npar)
    step1.idx <- integer(0L)

    for(mm in seq_len(nMMblocks)) {

        # create parameter table for this measurement block only
        PT.block <-
            lav_partable_subset_measurement_model(PT = PT,
                                                  lavpta = lavpta,
                                                  add.lv.cov = TRUE,
                                                  add.idx = TRUE,
                                                  lv.names = mm.list[[mm]])
        mm.idx <- attr(PT.block, "idx"); attr(PT.block, "idx") <- NULL

        ind.names <- lav_partable_vnames(PT.block, "ov.ind")
        LV.idx.list[[mm]] <- match(mm.list[[mm]], FIT@Model@dimNames[[1]][[2]])
        OV.idx.list[[mm]] <- match(ind.names, FIT@Model@dimNames[[1]][[1]])

        # fit this measurement model only
        fit.mm.block <- do.call("lavaan",
                                args =  c(list(model  = PT.block,
                                               data   = data), dotdotdot.mm) )
        # check convergence
        if(!lavInspect(fit.mm.block, "converged")) {
            # fatal for now
            stop("lavaan ERROR: measurement model for ",
                 paste(mm.list[[mm]], collapse = " "), " did not converge.")
        }

        # store fitted measurement model
        MM.FIT[[mm]] <- fit.mm.block

        # store LAMBDA/THETA
        LAMBDA.list[[mm]] <- computeLAMBDA(fit.mm.block@Model )
         THETA.list[[mm]] <- computeTHETA( fit.mm.block@Model )
        if(lavoptions$meanstructure) {
            NU.list[[mm]] <- computeNU( fit.mm.block@Model,
                                        lavsamplestats = FIT@SampleStats )
        }

        # fill in point estimates measurement block
        PTM <- MM.FIT[[mm]]@ParTable
        PT$est[ seq_len(length(PT$lhs)) %in% mm.idx & PT$free > 0L ] <-
            PTM$est[ PTM$free > 0L ]

        # compute `total' information for this measurement block
        mm.info <- lavTech(MM.FIT[[mm]], "information") * lavTech(FIT, "nobs")

        # fill in `total' information matrix
        par.idx <- PT.orig$free[ seq_len(length(PT$lhs)) %in% mm.idx &
                                 PT$free > 0L ]
        MM.INFO[par.idx, par.idx] <- mm.info

        # store indices in step1.idx
        step1.idx <- c(step1.idx, par.idx)

    } # measurement block

    # store MM fits (for now) in output
    out <- list()
    out$MM.FIT <- MM.FIT


    ## STEP 2: compute Var(eta) and E(eta) per group
    VETA <- vector("list", ngroups)
    if(lavoptions$meanstructure) {
        EETA <- vector("list", ngroups)
    } else {
        EETA <- NULL
    }

    # global MM matrices
    LAMBDA <- computeLAMBDA(FIT@Model)
    THETA  <- computeTHETA(FIT@Model, fix = FALSE) # keep zeroes for dummy lv
    if(lavoptions$meanstructure) {
        NU <- computeNU(FIT@Model, lavsamplestats = FIT@SampleStats)
    }
    M <- vector("list", ngroups)


    for(g in seq_len(ngroups)) {
        # assemble global LAMBDA/THETA (per group)
        for(mm in seq_len(nMMblocks)) {
            ov.idx <- OV.idx.list[[mm]]
            lv.idx <- LV.idx.list[[mm]]
            LAMBDA[[g]][ov.idx, lv.idx] <- LAMBDA.list[[mm]][[g]]
             THETA[[g]][ov.idx, ov.idx] <-  THETA.list[[mm]][[g]]
            if(lavoptions$meanstructure) {
                NU[[g]][ov.idx, 1] <- NU.list[[mm]][[g]]
            }
        }

        # check if LAMBDA has full column rank
        if(qr(LAMBDA[[g]])$rank < ncol(LAMBDA[[g]])) {
            print(LAMBDA[[g]])
            stop("lavaan ERROR: LAMBDA has no full column rank.")
        }

        # get sample statistics for this group
        ybar <- FIT@SampleStats@mean[[g]]
        cov  <- FIT@SampleStats@cov[[g]]

        # compute 'M'
        if(M.method == "GLS") {
            icov <- FIT@SampleStats@icov[[g]]
            Mg <- ( solve(t(LAMBDA[[g]]) %*% icov %*% LAMBDA[[g]]) %*%
                        t(LAMBDA[[g]]) %*% icov )
        } else if(M.method == "ML") {
            zero.theta.idx <- which(diag(THETA[[g]]) == 0)
            if(length(zero.theta.idx) > 0L) {
                tmp <- THETA[[g]][-zero.theta.idx, -zero.theta.idx,
                                  drop = FALSE]
                tmp.inv <- solve(tmp)
                THETA.inv <- THETA[[g]]
                THETA.inv[-zero.theta.idx, -zero.theta.idx] <- tmp.inv
                diag(THETA.inv)[zero.theta.idx] <- 1
            } else {
                THETA.inv <- solve(THETA[[g]])
            }
            Mg <- ( solve(t(LAMBDA[[g]]) %*% THETA.inv %*% LAMBDA[[g]]) %*%
                        t(LAMBDA[[g]]) %*% THETA.inv )
        } else if(M.method == "ULS") {
            Mg <- solve(t(LAMBDA[[g]]) %*%  LAMBDA[[g]]) %*% t(LAMBDA[[g]])
        }

        # compute VETA
        VETA[[g]] <- Mg %*% (cov - THETA[[g]]) %*% t(Mg)

        # names
        psi.idx <- which(names(FIT@Model@GLIST) == "psi")
        dimnames(VETA[[g]]) <- FIT@Model@dimNames[[psi.idx]]

        # compute EETA
        if(lavoptions$meanstructure) {
            EETA[[g]] <- M %*% (ybar - NU[[g]])
        }

        # store M
        M[[g]] <- Mg
    }

    # label groups
    if(ngroups > 1L) {
        names(VETA) <- FIT@Data@group.label
    }

    # store LAMBDA/THETA/NU per group
    out$LAMBDA <- LAMBDA
    out$THETA  <- THETA
    if(lavoptions$meanstructure) {
        out$NU     <- NU
    }

    # store EETA/VETA
    out$VETA <- VETA
    out$EETA <- EETA

    # store M
    out$M <- M


    # STEP 3: compute sample.cov and sample.mean for structural part

    # extract structural part
    PT.PA <- lav_partable_subset_structural_model(PT, lavpta = lavpta,
                                                  add.idx = TRUE)
    reg.idx <- attr(PT.PA, "idx"); attr(PT.PA, "idx") <- NULL
    group.values <- lav_partable_group_values(PT.PA)

    # adjust options
    lavoptions.PA <- lavoptions
    lavoptions.PA <- modifyList(lavoptions.PA, mm.args)

    # override, not matter what
    lavoptions.PA$do.fit <- TRUE
    lavoptions.PA$missing <- "listwise"
    lavoptions.PA$sample.cov.rescale <- FALSE

    # fit structural model
    FIT.PA <- lavaan::lavaan(PT.PA,
                             sample.cov = VETA,
                             sample.mean = EETA, # NULL if no meanstructure
                             sample.nobs = FIT@Data@nobs,
                             slotOptions = lavoptions.PA)

    # fill in point estimates structural part
    PTS <- FIT.PA@ParTable
    PT$est[ seq_len(length(PT$lhs)) %in% reg.idx & PT$free > 0L ] <-
        PTS$est[ PTS$free > 0L ]
    step2.idx <- PT$free[  seq_len(length(PT$lhs)) %in% reg.idx &
                           PT$free > 0L ]


    # compute standard errors
    # current approach:
    # - create 'global' model, only to get the 'joint' information matrix
    # - partition information matrix (step 1, step 2)
    # - apply two-step correction for second step
    # - 'insert' these corrected SEs (and vcov) in FIT.PA

    lavoptions.joint <- lavoptions
    lavoptions.joint$optim.method = "none"
    lavoptions.joint$check.gradient = FALSE
    lavoptions.joint$se = "none"
    lavoptions.joint$test = "none"
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
    VCOV <- V2 + V1

    # fill in standard errors step 2
    PTS$se[ PTS$free > 0L ] <- sqrt( diag(VCOV) )

    # overwrite slots in FIT.PA (better way?)
    FIT.PA@Options$se <- "twostep"
    FIT.PA@ParTable <- PTS
    FIT.PA@vcov$vcov <- VCOV


    # store FIT.PA
    out$FIT.PA <- FIT.PA


    # prepare output
    if(output == "list") {
        res <- out
    } else if(output == "lavaan") {
        res <- FIT.PA
    } else {
        res <- out
    }

    res
}



