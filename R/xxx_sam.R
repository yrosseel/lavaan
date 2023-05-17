# SAM: a Structural After Measurement approach
#
# Yves Rosseel & Wen-Wei Loh, Feb-May 2019

# local vs global sam
# local sam = alternative for FSR+Croon
# - but no need to compute factor scores or corrections
# gloal sam = (old) twostep
# - but we can also take a 'local' perspective

# restrictions:
#
# local and global:
#  - all (measured) latent variables must have indicators that are observed
# local:
#  - only if LAMBDA is of full column rank (eg no SRM, no bi-factor, no MTMM)
#  - if multiple groups: each group has the same set of latent variables!
#  - global approach is used to compute corrected two-step standard errors

# YR 12 May 2019 - first version
# YR 22 May 2019 - merge sam/twostep (call it 'local' vs 'global' sam)

# YR 27 June 2021 - prepare for `public' release
#                 - add Fuller (1987) 'lambda' correction if (MSM - MTM) is not
#                   positive definite
#                 - se = "none" now works
#                 - store 'local' information in @internal slot (for printing)

# YR 16 Oct 2021 - if an indicator is also a predictor/outcome in the
#                  structural part, treat it as an observed predictor
#                  without measurement error in the second step
#                  (ie, set THETA element to zero)

# YR 03 Dec 2022 - allow for sam.method = "fsr" and se = "naive"
#                - add alpha.correction= argument (for small sample correction)


# twostep = wrapper for global sam
twostep <- function(model = NULL, data = NULL, cmd = "sem",
                    mm.list = NULL, mm.args = list(), struc.args = list(),
                    ...,         # global options
                    output = "lavaan") {

    sam(model = model, data = data, cmd = cmd, mm.list = mm.list,
        mm.args = mm.args, struc.args = struc.args,
        sam.method = "global", # or global
                ...,         # global options
        output = output)
}

# fsr = wrapper for local sam
# TODO


sam <- function(model            = NULL,
                data             = NULL,
                cmd              = "sem",
                se               = "twostep",
                mm.list          = NULL,
                mm.args          = list(bounds = "wide.zerovar",
                                        se = "standard"),
                struc.args       = list(estimator = "ML",
                                        #fixed.x = TRUE,
                                        se = "standard"),
                sam.method       = "local", # or "global", or "fsr"
                alpha.correction = 0L,
                ...,           # global options
                local.options    = list(M.method = "ML",
                                        veta.force.pd = TRUE,
                                        twolevel.method = "h1"),
                                        # h1, anova, mean
                global.options   = list(), # not used for now
                output           = "lavaan") {

    # default local.options
    local.opt <- list(M.method = "ML",
                      veta.force.pd = TRUE,
                      twolevel.method = "h1")
    local.options <- modifyList(local.opt, local.options, keep.null = FALSE)

    # check arguments
    if(!is.null(local.opt[["M.method"]])) {
        local.M.method <- toupper(local.options[["M.method"]])
    }
    if(!local.M.method %in% c("GLS", "ML", "ULS")) {
        stop("lavaan ERROR: local option M.method should be one of GLS, ML or ULS.")
    }
    local.twolevel.method <- tolower(local.options[["twolevel.method"]])
    if(!local.twolevel.method %in% c("h1", "anova", "mean")) {
        stop("lavaan ERROR: local option twolevel.method should be one of h1, anova or mean.")
    }

    # output
    output <- tolower(output)
    if(output == "list" || output == "lavaan") {
        # nothing to do
    } else {
        stop("lavaan ERROR: output should be \"list\" or \"lavaan.\"")
    }

    # handle dot dot dot
    dotdotdot <- list(...)

    # STEP 0: process full model, without fitting
    dotdotdot0 <- dotdotdot
    dotdotdot0$do.fit <- NULL
    if(sam.method %in% c("local", "fsr")) {
        dotdotdot0$sample.icov <- FALSE # if N < nvar
    }
    dotdotdot0$se     <- "none"
    dotdotdot0$test   <- "none"
    dotdotdot0$verbose <- FALSE # no output for this 'dummy' FIT
    dotdotdot0$ceq.simple <- TRUE # if not the default yet

    # initial processing of the model, no fitting
    FIT <- do.call(cmd,
                   args =  c(list(model  = model,
                                  data   = data,
                                  do.fit = FALSE), dotdotdot0) )
    lavoptions <- lavInspect(FIT, "options")

    # restore options
    lavoptions$do.fit <- TRUE
    if(sam.method %in% c("local", "fsr")) {
        lavoptions$sample.icov <- TRUE
    }
    # se
    lavoptions$se   <- se
    # test
    if(!is.null(dotdotdot$test)) {
        lavoptions$test <- dotdotdot$test
    } else {
        lavoptions$test <- "standard"
    }
    # verbose
    if(!is.null(dotdotdot$verbose)) {
        lavoptions$verbose <- dotdotdot$verbose
    }

    if(lavoptions$verbose) {
        cat("This is sam using sam.method = ", sam.method, ".\n", sep = "")
    }


    # what have we learned?
    lavpta  <- FIT@pta
    ngroups <- lavpta$ngroups
    nlevels <- lavpta$nlevels
    nblocks <- lavpta$nblocks
    PT      <- FIT@ParTable


    # local only
    if(sam.method %in% c("local", "fsr")) {
        # if missing = "listwise", make data complete, to avoid different
        # datasets per measurement block
        if(lavoptions$missing == "listwise") {
            # FIXME: make this work for multiple groups!!
            OV <- unique(unlist(lavpta$vnames$ov))
            # add group/cluster/sample.weights variables (if any)
            OV <- c(OV, FIT@Data@group, FIT@Data@cluster,
                    FIT@Data@sampling.weights)
            data <- na.omit(data[,OV])
        }
    }

    # any `regular' latent variables? (across groups!)
    LV.names <- unique(unlist(FIT@pta$vnames$lv.regular))
    OV.names <- unique(unlist(FIT@pta$vnames$ov))

    # check for higher-order factors
    # 0.6-11: hard stop for now, as we do not support them (yet)!
    LV.IND.names <- unique(unlist(FIT@pta$vnames$lv.ind))
    if(length(LV.IND.names) > 0L) {
        stop("lavaan ERROR: model contains indicators that are also latent variables:\n\t", paste(LV.IND.names, collapse = " "))
        #ind.idx <- match(LV.IND.names, LV.names)
        #LV.names <- LV.names[-ind.idx]
    }

    # do we have at least 1 'regular' (measured) latent variable?
    if(length(LV.names) == 0L) {
        stop("lavaan ERROR: model does not contain any (measured) latent variables; use sem() instead")
    }
    nfac <- length(LV.names)

    # total number of free parameters
    if(FIT@Model@ceq.simple.only) {
        npar <- FIT@Model@nx.unco
        PT.free <- PT$free
        PT.free[ PT.free > 0 ] <- seq_len(npar)
    } else {
        npar <- FIT@Model@nx.free
        PT.free <- PT$free
    }
    if(npar < 1L) {
        stop("lavaan ERROR: model does not contain any free parameters")
    }

    # check parameter table
    PT$est <- PT$se <- NULL
    # est equals ustart by default (except exo values)
    PT$est <- PT$ustart
    if(any(PT$exo > 0L)) {
        PT$est[PT$exo > 0L] <- PT$start[PT$exo > 0L]
    }

    # clear se values (needed here?) only for global approach to compute SE
    PT$se <- rep(as.numeric(NA), length(PT$lhs))
    PT$se[ PT$free == 0L & !is.na(PT$ustart) ] <- 0.0

    # how many measurement models?
    if(!is.null(mm.list)) {
        nMMblocks <- length(mm.list)
        # check each measurement block
        for(b in seq_len(nMMblocks)) {
            # check if we can find all lv names in LV.names
            if(!all(unlist(mm.list[[b]]) %in% LV.names)) {
              tmp <- unlist(mm.list[[b]])
              stop("lavaan ERROR: mm.list contains unknown latent variable(s): ",
                paste( tmp[ !tmp %in% LV.names ], collapse = " "),
                "\n")
            }
            # make list per block
            if(!is.list(mm.list[[b]])) {
                mm.list[[b]] <- rep(list(mm.list[[b]]), nblocks)
            } else {
                if(length(mm.list[[b]]) != nblocks) {
                    stop("lavaan ERROR: mm.list block ", b, " has length ",
                         length(mm.list[[b]]), " but nblocks = ", nblocks)
                }
            }
        }
    } else {
        # TODO: here comes the automatic 'detection' of linked
        #       measurement models
        #
        # for now we take a single latent variable per measurement model block
        mm.list <- as.list(LV.names)
        nMMblocks <- length(mm.list)
        for(b in seq_len(nMMblocks)) {
            # make list per block
            mm.list[[b]] <- rep(list(mm.list[[b]]), nblocks)
        }
    }




    # STEP 1: fit each measurement model (block)

    if(lavoptions$verbose) {
        cat("Fitting the measurement part:\n")
    }

    # adjust options for measurement models
    dotdotdot.mm <- dotdotdot
    #dotdotdot.mm$se <- "none"
    #if(sam.method == "global") {
    #    dotdotdot.mm$test <- "none"
    #}
    # we need the tests to create summary info about MM
    dotdotdot.mm$debug <- FALSE
    dotdotdot.mm$verbose <- FALSE
    dotdotdot.mm$check.post <- FALSE # neg lv variances may be overriden
    dotdotdot.mm$check.gradient <- FALSE # too sensitive in large model (global)
    dotdotdot.mm$baseline <- FALSE

    # override with mm.args
    dotdotdot.mm <- modifyList(dotdotdot.mm, mm.args)

    # we assume the same number/names of lv's per group!!!
    MM.FIT <- vector("list", nMMblocks)         # fitted object

    # local only
    if(sam.method %in% c("local", "fsr")) {
        LAMBDA.list <- vector("list", nMMblocks)
        THETA.list  <- vector("list", nMMblocks)
        NU.list     <- vector("list", nMMblocks)
        LV.idx.list <- vector("list", nMMblocks)
        OV.idx.list <- vector("list", nMMblocks)
    }

    # for joint model later
    if(lavoptions$se != "none") {
        Sigma.11 <- matrix(0, npar, npar)
    }
    step1.idx <- integer(0L)


    # NOTE: we should explicitly add zero-constrained LV covariances
    # to PT, and keep them zero in PTM
    if(cmd == "lavaan") {
        add.lv.cov <- FALSE
    } else {
        add.lv.cov <- TRUE
    }


    for(mm in seq_len(nMMblocks)) {

        if(lavoptions$verbose) {
            cat("  block ", mm, "[",
                paste(mm.list[[mm]], collapse = " "), "]\n")
        }

        if(sam.method %in% c("local", "fsr")) {
            #s LV.idx.list/OV.idx.list: list per block
            LV.idx.list[[mm]] <- vector("list", nblocks)
            OV.idx.list[[mm]] <- vector("list", nblocks)
        }

        # create parameter table for this measurement block only
        PTM <- lav_partable_subset_measurement_model(PT = PT,
                                                     lavpta = lavpta,
                                                     add.lv.cov = add.lv.cov,
                                                     add.idx = TRUE,
                                                     lv.names = mm.list[[mm]])
        mm.idx <- attr(PTM, "idx"); attr(PTM, "idx") <- NULL
        PTM$est <- NULL
        PTM$se <- NULL

        # fit this measurement model only
        fit.mm.block <- do.call("lavaan",
                                args =  c(list(model  = PTM,
                                               data   = data), dotdotdot.mm) )
        # check convergence
        if(!lavInspect(fit.mm.block, "converged")) {
            # fatal for now
            stop("lavaan ERROR: measurement model for ",
                 paste(mm.list[[mm]], collapse = " "), " did not converge.")
        }

        # store fitted measurement model
        MM.FIT[[mm]] <- fit.mm.block

        if(sam.method %in% c("local", "fsr")) {
            # store LAMBDA/THETA
            LAMBDA.list[[mm]] <- computeLAMBDA(fit.mm.block@Model )
             THETA.list[[mm]] <- computeTHETA( fit.mm.block@Model )
            if(lavoptions$meanstructure) {
                NU.list[[mm]] <- computeNU( fit.mm.block@Model,
                                            lavsamplestats = FIT@SampleStats )
            }

            # store indices
            for(bb in seq_len(nblocks)) {
                lambda.idx <- which(names(FIT@Model@GLIST) == "lambda")[bb]
                ind.names <- lav_partable_vnames(PTM, "ov.ind", block = bb)
                LV.idx.list[[mm]][[bb]] <- match(mm.list[[mm]][[bb]],
                    FIT@Model@dimNames[[lambda.idx]][[2]])
                OV.idx.list[[mm]][[bb]] <- match(ind.names,
                    FIT@Model@dimNames[[lambda.idx]][[1]])
            }
        }

        # fill in point estimates measurement block (including slack values)
        PTM <- MM.FIT[[mm]]@ParTable
        PT$est[ seq_len(length(PT$lhs)) %in% mm.idx &
               (PT$free > 0L | PT$op %in% c(":=", ">", "<")) ] <-
            PTM$est[ (PTM$free > 0L | PTM$op %in% c(":=", "<", ">")) &
                     PTM$user != 3L ]

        # fill in standard errors measurement block
        if(lavoptions$se != "none") {

            if(fit.mm.block@Model@ceq.simple.only) {
                PTM.free <- PTM$free
                PTM.free[ PTM.free > 0 ] <- seq_len(fit.mm.block@Model@nx.unco)
            } else {
                PTM.free <- PTM$free
            }


            PT$se[ seq_len(length(PT$lhs)) %in% mm.idx & PT$free > 0L ] <-
                PTM$se[ PTM$free > 0L & PTM$user != 3L]

            # compute variance matrix for this measurement block
            sigma.11 <- MM.FIT[[mm]]@vcov$vcov

            # fill in variance matrix
            par.idx <- PT.free[ seq_len(length(PT$lhs)) %in% mm.idx &
                                PT$free > 0L ]
            keep.idx <- PTM.free[ PTM$free > 0 & PTM$user != 3L ]
            Sigma.11[par.idx, par.idx] <-
                sigma.11[keep.idx, keep.idx, drop = FALSE]

            # store indices in step1.idx
            step1.idx <- c(step1.idx, par.idx)
        }

    } # measurement block

    # only keep 'measurement part' parameters in Sigma.11
    if(lavoptions$se != "none") {
        Sigma.11 <- Sigma.11[step1.idx, step1.idx, drop = FALSE]
    }

    # store MM fits (for now) in output
    out <- list()
    out$MM.FIT <- MM.FIT

    if(sam.method %in% c("local", "fsr")) {
        if(lavoptions$verbose) {
            cat("Constructing the mapping matrix using the ",
                local.M.method, " method ... ", sep = "")
        }

        # assemble global LAMBDA/THETA (per block)
        LAMBDA <- computeLAMBDA(FIT@Model, handle.dummy.lv = FALSE)
        THETA  <- computeTHETA(FIT@Model, fix = FALSE) # keep dummy lv
        if(lavoptions$meanstructure) {
            NU <- computeNU(FIT@Model, lavsamplestats = FIT@SampleStats)
        }
        for(b in seq_len(nblocks)) {
            for(mm in seq_len(nMMblocks)) {
                ov.idx <- OV.idx.list[[mm]][[b]]
                lv.idx <- LV.idx.list[[mm]][[b]]
                LAMBDA[[b]][ov.idx, lv.idx] <- LAMBDA.list[[mm]][[b]]
                 THETA[[b]][ov.idx, ov.idx] <-  THETA.list[[mm]][[b]]
                # new in 0.6-10: check if any indicators are also involved
                # in the structural part; if so, set THETA row/col to zero
                # and make sure LAMBDA element is correctly set
                # (we also need to adjust M)
                dummy.ov.idx <- FIT@Model@ov.y.dummy.ov.idx[[b]]
                dummy.lv.idx <- FIT@Model@ov.y.dummy.lv.idx[[b]]
                if(length(dummy.ov.idx)) {
                    THETA[[b]][dummy.ov.idx,] <- 0
                    THETA[[b]][,dummy.ov.idx] <- 0
                    LAMBDA[[b]][dummy.ov.idx,] <- 0
                    LAMBDA[[b]][cbind(dummy.ov.idx, dummy.lv.idx)] <- 1
                }
                if(lavoptions$meanstructure) {
                    NU[[b]][ov.idx, 1] <- NU.list[[mm]][[b]]
                    if(length(dummy.ov.idx)) {
                        NU[[b]][dummy.ov.idx, 1] <- 0
                    }
                }
            }

            # check if LAMBDA has full column rank
            if(qr(LAMBDA[[b]])$rank < ncol(LAMBDA[[b]])) {
                print(LAMBDA[[b]])
                stop("lavaan ERROR: LAMBDA has no full column rank. Please use sam.method = global")
            }
        } # b

        # store LAMBDA/THETA/NU per block
        out$LAMBDA <- LAMBDA
        out$THETA  <- THETA
        if(lavoptions$meanstructure) {
            out$NU     <- NU
        }
    }



    ## STEP 1b: compute Var(eta) and E(eta) per block
    ##          only needed for local approach!
    if(sam.method %in% c("local", "fsr")) {
        VETA   <- vector("list", nblocks)
        REL    <- vector("list", nblocks)
        alpha  <- vector("list", nblocks)
        lambda <- vector("list", nblocks)
        if(lavoptions$meanstructure) {
            EETA <- vector("list", nblocks)
        } else {
            EETA <- NULL
        }
        M <- vector("list", nblocks)

        # compute VETA/EETA per block
        if(nlevels > 1L && local.twolevel.method == "h1") {
            H1 <- lav_h1_implied_logl(lavdata = FIT@Data,
                                      lavsamplestats = FIT@SampleStats,
                                      lavoptions     = FIT@Options)
        }

        for(b in seq_len(nblocks)) {

            # get sample statistics for this block
            if(nlevels > 1L) {
                if(ngroups > 1L) {
                    this.level <- (b - 1L) %% ngroups + 1L
                } else {
                    this.level <- b
                }
                this.group <- floor(b/nlevels + 0.5)

                if(this.level == 1L) {

                    if(local.twolevel.method == "h1") {
                        COV  <- H1$implied$cov[[1]]
                        YBAR <- H1$implied$mean[[1]]
                    } else if(local.twolevel.method == "anova" ||
                              local.twolevel.method == "mean") {
                        COV  <- FIT@SampleStats@YLp[[this.group]][[2]]$Sigma.W
                        YBAR <- FIT@SampleStats@YLp[[this.group]][[2]]$Mu.W
                    }

                    # reduce
                    ov.idx <- FIT@Data@Lp[[this.group]]$ov.idx[[this.level]]
                    COV <- COV[ov.idx, ov.idx, drop = FALSE]
                    YBAR <- YBAR[ov.idx]
                } else if(this.level == 2L) {
                    if(local.twolevel.method == "h1") {
                        COV  <- H1$implied$cov[[2]]
                        YBAR <- H1$implied$mean[[2]]
                    } else if(local.twolevel.method == "anova") {
                        COV  <- FIT@SampleStats@YLp[[this.group]][[2]]$Sigma.B
                        YBAR <- FIT@SampleStats@YLp[[this.group]][[2]]$Mu.B
                    } else if(local.twolevel.method == "mean") {
                        S.PW <- FIT@SampleStats@YLp[[this.group]][[2]]$Sigma.W
                        NJ   <- FIT@SampleStats@YLp[[this.group]][[2]]$s
                        Y2   <- FIT@SampleStats@YLp[[this.group]][[2]]$Y2
                        # grand mean
                        MU.Y <- ( FIT@SampleStats@YLp[[this.group]][[2]]$Mu.W +                                   FIT@SampleStats@YLp[[this.group]][[2]]$Mu.B )
                        Y2c <- t( t(Y2) - MU.Y ) # MUST be centered
                        YB <- crossprod(Y2c)/nrow(Y2c)
                        COV  <- YB - 1/NJ * S.PW
                        YBAR <- FIT@SampleStats@YLp[[this.group]][[2]]$Mu.B
                    }

                    # reduce
                    ov.idx <- FIT@Data@Lp[[this.group]]$ov.idx[[this.level]]
                    COV <- COV[ov.idx, ov.idx, drop = FALSE]
                    YBAR <- YBAR[ov.idx]
                } else {
                    stop("lavaan ERROR: level 3 not supported (yet).")
                }
            } else {
                YBAR <- FIT@h1$implied$mean[[b]] # EM version if missing="ml"
                COV  <- FIT@h1$implied$cov[[b]]
                if(local.M.method == "GLS") {
                    if(FIT@Options$sample.cov.rescale) {
                       # get unbiased S
                       N <- FIT@SampleStats@nobs[[b]]
                       COV.unbiased <- COV * N/(N-1)
                       ICOV <- solve(COV.unbiased)
                    } else {
                        ICOV <- solve(COV)
                    }
                }
            }

            # compute mapping matrix 'M'
            Mb <- lav_sam_mapping_matrix(LAMBDA = LAMBDA[[b]],
                                         THETA = THETA[[b]],
                                         S = COV, S.inv = ICOV,
                                         method = local.M.method,
                                         warn = lavoptions$warn)

            # handle observed-only variables
            dummy.ov.idx <- FIT@Model@ov.y.dummy.ov.idx[[b]]
            dummy.lv.idx <- FIT@Model@ov.y.dummy.lv.idx[[b]]
            if(length(dummy.ov.idx)) {
                Mb[dummy.lv.idx,] <- 0
                Mb[cbind(dummy.lv.idx, dummy.ov.idx)] <- 1
            }

            # compute VETA
            if(sam.method == "local") {
                tmp <- lav_sam_veta(M = Mb, S = COV, THETA = THETA[[b]],
                           alpha.correction = alpha.correction,
                           lambda.correction = local.options[["veta.force.pd"]],
                           N <- nobs(FIT), extra = TRUE)
                VETA[[b]] <- tmp[,,drop = FALSE] # drop attributes
                alpha[[b]]  <- attr(tmp, "alpha")
                lambda[[b]] <- attr(tmp, "lambda")
            } else {
                # FSR -- no correction
                VETA[[b]] <- Mb %*% COV %*% t(Mb)
            }

            # standardize? not really needed, but we may have 1.0000001
            # as variances, and this may lead to false convergence
            if(FIT@Options$std.lv) {
                VETA[[b]] <- stats::cov2cor(VETA[[b]])
            }

            # names
            psi.idx <- which(names(FIT@Model@GLIST) == "psi")[b]
            dimnames(VETA[[b]]) <- FIT@Model@dimNames[[psi.idx]]

            # compute EETA
            if(lavoptions$meanstructure) {
                EETA[[b]] <- lav_sam_eeta(M = Mb, YBAR = YBAR, NU = NU[[b]])
            }

            # compute model-based reliability
            MSM <- Mb %*% COV %*% t(Mb)
            #REL[[b]] <- diag(VETA[[b]] %*% solve(MSM)) # CHECKme!
            REL[[b]] <- diag(VETA[[b]]) / diag(MSM)

            # store M
            M[[b]] <- Mb

        } # blocks

        # label blocks
        if(nblocks > 1L) {
            names(VETA) <- FIT@Data@block.label
            names(REL)  <- FIT@Data@block.label
        }

        # store EETA/VETA/M/alpha/lambda
        out$VETA <- VETA; out$EETA <- EETA; out$REL  <- REL
        out$M <- M; out$lambda <- lambda; out$alpha  <- alpha

        if(lavoptions$verbose) {
            cat("done.\n")
        }

    } # local


    # do we have any parameters left?
    if(length(unique(step1.idx)) >= npar) {
        warning("lavaan WARNING: ",
                "no free parameters left for structural part.\n",
                "        Returning measurement part only.")
        if(output == "list") {
            return(out)
        } else {
            if(nMMblocks == 1L) {
                return(MM.FIT[[1]])
            } else {
                return(MM.FIT)
            }
        }
    }



    ####################################
    # STEP 2: estimate structural part #
    ####################################

    # adjust options
    lavoptions.PA <- lavoptions
    #lavoptions.PA$fixed.x <- TRUE # may be false if indicator is predictor
    lavoptions.PA$fixed.x <- FALSE # until we fix this...
    lavoptions.PA$verbose <- FALSE # must be in struc.args
    lavoptions.PA$categorical <- FALSE
    lavoptions.PA$.caegorical <- FALSE
    lavoptions.PA <- modifyList(lavoptions.PA, struc.args)

    # override, no matter what
    lavoptions.PA$do.fit <- TRUE

    if(sam.method %in% c("local", "fsr")) {
        lavoptions.PA$missing <- "listwise"
        #if(lavoptions$se != "naive") {
        #    lavoptions.PA$se <- "none" # sample statistics input
        #}
        lavoptions.PA$sample.cov.rescale <- FALSE
        #lavoptions.PA$baseline <- FALSE
        lavoptions.PA$h1 <- FALSE
        #lavoptions.PA$implied <- FALSE
        lavoptions.PA$loglik <- FALSE
    } else {
        #lavoptions.PA$baseline <- FALSE
        lavoptions.PA$h1 <- FALSE
        #lavoptions.PA$implied <- FALSE
        lavoptions.PA$loglik <- FALSE
    }

    # construct PTS
    if(sam.method %in% c("local", "fsr")) {
        # extract structural part
        PTS <- lav_partable_subset_structural_model(PT, lavpta = lavpta,
                   add.idx = TRUE,
                   fixed.x = lavoptions.PA$fixed.x,
                   add.exo.cov = TRUE)
        # remove est/se/start columns
        PTS$est   <- NULL
        PTS$se    <- NULL
        PTS$start <- NULL

        if(nlevels > 1L) {
            PTS$level <- NULL
            PTS$group <- NULL
            PTS$group <- PTS$block
            NOBS <- FIT@Data@Lp[[1]]$nclusters
        } else {
            NOBS <- FIT@Data@nobs
        }
        # if meanstructure, 'free' user=0 intercepts?
        if(lavoptions.PA$meanstructure) {
            extra.int.idx <- which(PTS$op == "~1" & PTS$user == 0L &
                                   PTS$exo == 0L)
            if(length(extra.int.idx) > 0L) {
                PTS$free[  extra.int.idx ] <- 1L
                PTS$ustart[extra.int.idx ] <- as.numeric(NA)
                PTS$free[ PTS$free > 0L ] <-
                    seq_len( length(PTS$free[ PTS$free > 0L ]) )
            }
        } else {
            extra.int.idx <- integer(0L)
        }
        reg.idx <- attr(PTS, "idx"); attr(PTS, "idx") <- NULL
    } else {
        # the measurement model parameters now become fixed ustart values
        PT$ustart[PT$free > 0] <- PT$est[PT$free > 0]

        reg.idx <- lav_partable_subset_structural_model(PT = PT,
                          lavpta = lavpta, idx.only = TRUE)

        # remove 'exogenous' factor variances (if any) from reg.idx
        lv.names.x <- LV.names[ LV.names %in% unlist(lavpta$vnames$eqs.x)  &
                               !LV.names %in% unlist(lavpta$vnames$eqs.y) ]
        if(lavoptions.PA$fixed.x && length(lv.names.x) > 0L) {
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
        con.idx <- which(PTS$op %in% c("==","<",">",":="))
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

        extra.int.idx <- integer(0L)
    } # global

    # fit structural model
    if(lavoptions$verbose) {
        cat("Fitting the structural part ... \n")
    }
    if(sam.method %in% c("local", "fsr")) {
        FIT.PA <- lavaan::lavaan(PTS,
                                 sample.cov = VETA,
                                 sample.mean = EETA, # NULL if no meanstructure
                                 sample.nobs = NOBS,
                                 slotOptions = lavoptions.PA)

    } else {
        FIT.PA <- lavaan::lavaan(model = PTS,
                                 slotData = FIT@Data,
                                 slotSampleStats = FIT@SampleStats,
                                 slotOptions = lavoptions.PA)
    }
    if(lavoptions$verbose) {
        cat("Fitting the structural part ... done.\n")
    }
    # store FIT.PA
    out$FIT.PA <- FIT.PA

    # fill in point estimates structural part
    PTS <- FIT.PA@ParTable
    p2def.idx <- seq_len(length(PT$lhs)) %in% reg.idx &
                 (PT$free > 0 | PT$op %in% c(":=", "<", ">"))

    # which parameters from PTS do we wish to fill in:
    # - all 'free' parameters
    # - :=, <, > (if any)
    # - but NOT elements in extra.int.idx
    # - and NOT element with user=3 (add.exo.cov = TRUE)
    pts.idx <- which( (PTS$free > 0L | (PTS$op %in% c(":=", "<", ">"))) &
                      !PTS$user == 3L &
                      !seq_len(length(PTS$lhs)) %in% extra.int.idx )

    # find corresponding rows in PT
    PTS2 <- as.data.frame(PTS, stringsAsFactors = FALSE)
    pt.idx <- lav_partable_map_id_p1_in_p2(PTS2[pts.idx,], PT,
                                           exclude.nonpar = FALSE)
    # fill in
    PT$est[ pt.idx ] <- PTS$est[ pts.idx ]


    # create step2.idx
    p2.idx <- seq_len(length(PT$lhs)) %in% reg.idx & PT$free > 0 # no def!
    step2.idx <- PT.free[ p2.idx ]

    # add 'step' column in PT
    PT$step <- rep(1L, length(PT$lhs))
    PT$step[seq_len(length(PT$lhs)) %in% reg.idx] <- 2L



    ################################################################
    # Step 3: assemble results in a 'dummy' JOINT model for output #
    ################################################################
    if(lavoptions$verbose) {
        cat("Assembling results for output ... ")
    }
    lavoptions.joint <- lavoptions
    lavoptions.joint$optim.method <- "none"
    lavoptions.joint$optim.force.converged <- TRUE
    PT$ustart <- PT$est # as this is used if optim.method == "none"
    lavoptions.joint$check.gradient <- FALSE
    lavoptions.joint$check.start <- FALSE
    lavoptions.joint$check.post <- FALSE
    if(sam.method %in% c("local", "fsr")) {
        lavoptions.joint$baseline <- FALSE
        lavoptions.joint$sample.icov <- FALSE
        lavoptions.joint$h1 <- FALSE
        lavoptions.joint$test <- "none"
        lavoptions.joint$estimator <- "none"
    } else {
        lavoptions.joint$test <- lavoptions$test
        lavoptions.joint$estimator <- lavoptions$estimator
    }
    lavoptions.joint$se   <- "none"
    lavoptions.joint$store.vcov <- FALSE # we do this manually
    lavoptions.joint$verbose <- FALSE

    JOINT <- lavaan::lavaan(PT, slotOptions = lavoptions.joint,
                            slotSampleStats = FIT@SampleStats,
                            slotData = FIT@Data)
    if(lavoptions$verbose) {
        cat("done.\n")
    }


    ###################################
    # Step 4: compute standard errors #
    ###################################

    # current approach:
    # - create 'global' model, only to get the 'joint' information matrix
    # - partition information matrix (step 1, step 2)
    # - apply two-step correction for second step
    # - 'insert' these corrected SEs (and vcov) in FIT.PA
    # compute information matrix

    if(lavoptions$se == "none") {
        # nothing to do...
    } else {
        if(lavoptions$verbose) {
            cat("Computing ", lavoptions$se, " standard errors ... ", sep = "")
        }
        JOINT@Model@estimator <- FIT@Options$estimator # could be DWLS!
        JOINT@Options$se <- lavoptions$se # always set to standard?
        if(JOINT@Model@ceq.simple.only) {
            VCOV.ALL <-  matrix(0, JOINT@Model@nx.unco,
                                   JOINT@Model@nx.unco)
        } else {
            VCOV.ALL <-  matrix(0, JOINT@Model@nx.free,
                                   JOINT@Model@nx.free)
        }
        VCOV.ALL[step1.idx, step1.idx] <- Sigma.11
        JOINT@vcov <- list(se = lavoptions$se,
                           information = lavoptions$information,
                           vcov = VCOV.ALL)

        INFO <- lavInspect(JOINT, "information")
        I.12 <- INFO[step1.idx, step2.idx]
        I.22 <- INFO[step2.idx, step2.idx]
        I.21 <- INFO[step2.idx, step1.idx]

        # compute Sigma.11
        # overlap? set corresponding rows/cols of Sigma.11 to zero
        both.idx <- which(step1.idx %in% step2.idx)
        if(length(both.idx) > 0L) {
            Sigma.11[both.idx,] <- 0
            Sigma.11[,both.idx] <- 0
        }

        # V2
        if(nlevels > 1L) {
            # FIXME: not ok for multigroup multilevel
            N <- FIT@Data@Lp[[1]]$nclusters[[2]] # first group only
        } else {
            N <- nobs(FIT)
        }

        # invert augmented information, for I.22 block only
        # new in 0.6-16 (otherwise, eq constraints in struc part are ignored)
        I.22.inv <-
            lav_model_information_augment_invert(lavmodel = FIT.PA@Model,
                     information = I.22, inverted = TRUE)

        # method below has the advantage that we can use a 'robust' vcov
        # for the joint model;
        # but does not work if we have equality constraints in the MM!
        # -> D will be singular
        #A <- JOINT@vcov$vcov[ step2.idx,  step2.idx]
        #B <- JOINT@vcov$vcov[ step2.idx, -step2.idx]
        #C <- JOINT@vcov$vcov[-step2.idx,  step2.idx]
        #D <- JOINT@vcov$vcov[-step2.idx, -step2.idx]
        #I.22.inv <- A - B %*% solve(D) %*% C

        if(lavoptions$se == "standard") {
            VCOV <- 1/N * I.22.inv
            out$VCOV <- VCOV
        } else if(lavoptions$se == "naive") {
            VCOV <- FIT.PA@vcov$vcov
            out$VCOV <- VCOV
        } else {

            # FIXME:
            V2 <- 1/N * I.22.inv # not the same as FIT.PA@vcov$vcov!!
            #V2 <- JOINT@vcov$vcov[ step2.idx,  step2.idx]
            #V2 <- FIT.PA@vcov$vcov

            # V1
            #V1 <- I.22.inv %*% I.21 %*% Sigma.11 %*% I.12 %*% I.22.inv


#            theta1totheta2 <- function(x) {
#
#
#                lavoptions.PA$se <- "none"
#                FIT.PA <- lavaan::lavaan(PTS,
#                                 sample.cov = VETA,
#                                 sample.mean = EETA, # NULL if no meanstructure
#                                 sample.nobs = NOBS,
#                                 slotOptions = lavoptions.PA)
#                out <- coef(FIT.PA)
#                out
#            }
#

            V1 <- I.22.inv %*% I.21 %*% Sigma.11 %*% I.12 %*% I.22.inv



            # V for second step
            if(alpha.correction > 0) {
                alpha.N1 <- alpha.correction / (N - 1)
                if(alpha.N1 > 1.0) {
                    alpha.N1 <- 1.0
                } else if (alpha.N1 < 0.0) {
                    alpha.N1 <- 0.0
                }
                VCOV.naive <- FIT.PA@vcov$vcov
                VCOV.corrected <- V2 + V1
                VCOV <- alpha.N1 * VCOV.naive + (1 - alpha.N1) * VCOV.corrected
            } else {
                VCOV <- V2 + V1
            }

            # store in out
            out$V2 <- V2
            out$V1 <- V1
            out$VCOV <- VCOV
        }
        if(lavoptions$verbose) {
            cat("done.\n")
        }
    } # se != "none



    ##################
    # Step 5: Output #
    ##################

    # assemble final lavaan objects
    if(output == "lavaan") {
        sam.mm.table <- data.frame(
            Block  = seq_len(length(mm.list)),
            Latent = sapply(MM.FIT, function(x) {
                      paste(unique(unlist(x@pta$vnames$lv)), collapse=",")}),
            Nind = sapply(MM.FIT, function(x) {
                       length(unique(unlist(x@pta$vnames$ov)))}),
            #Estimator = sapply(MM.FIT, function(x) { x@Model@estimator} ),
            Chisq  = sapply(MM.FIT, function(x) {x@test[[1]]$stat}),
            Df     = sapply(MM.FIT, function(x) {x@test[[1]]$df}) )
            #pvalue = sapply(MM.FIT, function(x) {x@test[[1]]$pvalue}) )
        class(sam.mm.table) <- c("lavaan.data.frame", "data.frame")


        # extra info for @internal slot
        if(sam.method %in% c("local", "fsr")) {
            sam.struc.fit <- try(fitMeasures(FIT.PA,
                                               c("chisq", "df", # "pvalue",
                                                 "cfi", "rmsea", "srmr")),
                                 silent = TRUE)
            if(inherits(sam.struc.fit, "try-error")) {
                sam.struc.fit <- "(unable to obtain fit measures)"
                names(sam.struc.fit) <- "warning"
            }
            sam.mm.rel <- REL
        } else {
            sam.struc.fit <- "no local fit measures available for structural part if sam.method is global"
            names(sam.struc.fit) <- "warning"
            sam.mm.rel <- numeric(0L)
        }


        SAM <- list(sam.method          = sam.method,
                    sam.local.options   = local.options,
                    sam.global.options  = global.options,
                    sam.mm.list         = mm.list,
                    sam.mm.estimator    = MM.FIT[[1]]@Model@estimator,
                    sam.mm.args         = mm.args,
                    sam.mm.ov.names     = lapply(MM.FIT, function(x) {
                                                 x@pta$vnames$ov }),
                    sam.mm.table        = sam.mm.table,
                    sam.mm.rel          = sam.mm.rel,
                    sam.struc.estimator = FIT.PA@Model@estimator,
                    sam.struc.args      = struc.args,
                    sam.struc.fit       = sam.struc.fit
                   )
        JOINT@internal <- SAM

        # fill in twostep standard errors
        if(JOINT@Options$se != "none") {
            JOINT@Options$se <- lavoptions$se
            JOINT@vcov$se    <- lavoptions$se
            JOINT@vcov$vcov[step2.idx, step2.idx] <- VCOV
            PT$se <- lav_model_vcov_se(lavmodel = JOINT@Model,
                                       lavpartable = PT,
                                       VCOV = JOINT@vcov$vcov)
            JOINT@ParTable <- PT
        }

        # fill information from FIT.PA
        JOINT@Options$optim.method <- FIT.PA@Options$optim.method
        if(sam.method %in% c("local", "fsr")) {
            JOINT@optim <- FIT.PA@optim
            JOINT@test  <- FIT.PA@test
        }

    } # output = "lavaan"


    # prepare output
    if(output == "lavaan") {
        res <- JOINT
    } else {
        res <- out
    }

    if(lavoptions$verbose) {
        cat("End of sam.\n")
    }

    res
}

