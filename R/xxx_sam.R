# SAM: a Structural After Measurement approach
#
# Yves Rosseel & Wen-Wei Loh, Feb-May 2019

# local vs global sam
# local sam = alternative for FSR+Croon
# - but no need to compute factor scores or corrections
# gloal sam = (old) twostep
# - but we can also take a 'local' perspective

# restrictions
# local:
#  - only if LAMBDA is of full column rank (eg no SRM, no bi-factor, no MTMM)
#  - if multiple groups: each group has the same set of latent variables!
#  - global approach is used to compute corrected two-step standard errors
# global:
#  - (none)?

# YR 12 May 2019 - first version
# YR 22 May 2019 - merge sam/twostep (call it 'local' vs 'global' sam)


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


sam <- function(model          = NULL,
                data           = NULL,
                cmd            = "sem",
                mm.list        = NULL,
                mm.args        = list(),
                struc.args     = list(),
                sam.method     = "local", # or global
                ...,           # global options
                local.options  = list(M.method = "ML",
                                      veta.force.pd = FALSE,
                                      veta.force.pd.tol = 1e-03,
                                      twolevel.method = "h1"), # h1, anova, mean
                global.options = list(), # not used for now
                output         = "lavaan") {

    # default local.options
    local.opt <- list(M.method = "ML",
                      veta.force.pd = FALSE,
                      veta.force.pd.tol = 1e-03,
                      twolevel.method = "h1")
    local.options <- modifyList(local.opt, local.options, keep.null = FALSE)

    # check arguments
    if(!is.null(local.opt[["M.method"]])) {
        local.M.method <- toupper(local.options[["M.method"]])
    }
    if(!local.M.method %in% c("GLS", "ML", "ULS")) {
        stop("lavaan ERROR: local.M.method should be one of GLS, ML or ULS.")
    }
    local.twolevel.method <- tolower(local.options[["twolevel.method"]])
    if(!local.twolevel.method %in% c("h1", "anova", "mean")) {
        stop("lavaan ERROR: local.twolevel.method should be one of h1, anova or mean.")
    }

    # output
    output <- tolower(output)
    if(output == "list" || output == "lavaan") {
        # nothing to do
    } else {
        stop("lavaan ERROR: output should be one list or lavaan.")
    }

    # handle dot dot dot
    dotdotdot <- list(...)

    # STEP 0: process full model, without fitting
    dotdotdot0 <- dotdotdot
    dotdotdot0$do.fit <- NULL
    dotdotdot0$se     <- "none"
    dotdotdot0$test   <- "none"

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
    lavpta  <- FIT@pta
    ngroups <- lavpta$ngroups
    nlevels <- lavpta$nlevels
    nblocks <- lavpta$nblocks
    PT      <- FIT@ParTable


    # local only
    if(sam.method == "local") {
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
              stop("lavaan ERROR: mm.list contains unknown latent variable(s):",
                paste( tmp[ !tmp %in% LV.names ], sep = " "),
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

    # adjust options for measurement models
    dotdotdot.mm <- dotdotdot
    #dotdotdot.mm$se <- "none"
    if(sam.method == "global") {
        dotdotdot.mm$test <- "none"
    }
    dotdotdot.mm$debug <- FALSE
    dotdotdot.mm$verbose <- FALSE
    dotdotdot.mm$check.post <- FALSE # neg lv variances may be overriden

    # override with mm.args
    dotdotdot.mm <- modifyList(dotdotdot.mm, mm.args)

    # we assume the same number/names of lv's per group!!!
    MM.FIT <- vector("list", nMMblocks)         # fitted object

    # local only
    if(sam.method == "local") {
        LAMBDA.list <- vector("list", nMMblocks)
        THETA.list  <- vector("list", nMMblocks)
        NU.list     <- vector("list", nMMblocks)
        LV.idx.list <- vector("list", nMMblocks)
        OV.idx.list <- vector("list", nMMblocks)
    }

    # for joint model later
    #MM.INFO <- matrix(0, npar, npar)
    Sigma.11 <- matrix(0, npar, npar)
    step1.idx <- integer(0L)


    # local only
    if(sam.method == "local") {
        # NOTE: we should explicitly add zero-constrained LV covariances
        # to PT, and keep them zero in PTM
        if(cmd == "lavaan") {
            add.lv.cov <- FALSE
        } else {
            add.lv.cov <- TRUE
        }
    }


    for(mm in seq_len(nMMblocks)) {

        if(sam.method == "local") {
            # LV.idx.list/OV.idx.list: list per block
            LV.idx.list[[mm]] <- vector("list", nblocks)
            OV.idx.list[[mm]] <- vector("list", nblocks)

            # create parameter table for this measurement block only
            PTM <-
                lav_partable_subset_measurement_model(PT = PT,
                                                      lavpta = lavpta,
                                                      add.lv.cov = add.lv.cov,
                                                      add.idx = TRUE,
                                                      lv.names = mm.list[[mm]])
            mm.idx <- attr(PTM, "idx"); attr(PTM, "idx") <- NULL
            PTM$est <- NULL
            PTM$se <- NULL
        } else {
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
            # - remove non-needed defs
            PTM <- PT
            if(any(sapply(mm.list[[mm]], length) > 1L)) {
                PTM <- lav_partable_add_lv_cov(PT = PT,
                              lavpta = lavpta, lv.names = mm.list[[mm]])
            }
            con.idx <- which(PTM$op %in% c("==","<",">",":="))
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
        }



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

        if(sam.method == "local") {
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

    # store MM fits (for now) in output
    out <- list()
    out$MM.FIT <- MM.FIT

    # do we have any parameters left?
    if(length(step1.idx) >= npar) {
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

    if(sam.method == "local") {
        # assemble global LAMBDA/THETA (per block)
        LAMBDA <- computeLAMBDA(FIT@Model)
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
                if(lavoptions$meanstructure) {
                    NU[[b]][ov.idx, 1] <- NU.list[[mm]][[b]]
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
    if(sam.method == "local") {
        VETA <- vector("list", nblocks)
        if(lavoptions$meanstructure) {
            EETA <- vector("list", nblocks)
        } else {
            EETA <- NULL
        }
        M <- vector("list", nblocks)

        # compute VETA/EETA per block
        if(nlevels > 1L && local.twolevel.method == "h1") {
            out <- lav_h1_implied_logl(lavdata = FIT@Data,
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
                        COV  <- out$implied$cov[[1]]
                        YBAR <- out$implied$mean[[1]]
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
                        COV  <- out$implied$cov[[2]]
                        YBAR <- out$implied$mean[[2]]
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
                YBAR <- FIT@SampleStats@mean[[b]]
                COV  <- FIT@SampleStats@cov[[b]]
                if(local.M.method == "GLS") {
                    ICOV <- FIT@SampleStats@icov[[b]]
                }
            }

            # compute 'M'
            if(local.M.method == "GLS") {
                Mg <- ( solve(t(LAMBDA[[b]]) %*% ICOV %*% LAMBDA[[b]]) %*%
                            t(LAMBDA[[b]]) %*% ICOV )
            } else if(local.M.method == "ML") {
                zero.theta.idx <- which(diag(THETA[[b]]) == 0)
                if(length(zero.theta.idx) > 0L) {
                    tmp <- THETA[[b]][-zero.theta.idx, -zero.theta.idx,
                                      drop = FALSE]
                    tmp.inv <- solve(tmp)
                    THETA.inv <- THETA[[b]]
                    THETA.inv[-zero.theta.idx, -zero.theta.idx] <- tmp.inv
                    diag(THETA.inv)[zero.theta.idx] <- 1
                } else {
                    THETA.inv <- solve(THETA[[b]])
                }
                Mg <- ( solve(t(LAMBDA[[b]]) %*% THETA.inv %*% LAMBDA[[b]]) %*%
                            t(LAMBDA[[b]]) %*% THETA.inv )
            } else if(local.M.method == "ULS") {
                Mg <- solve(t(LAMBDA[[b]]) %*%  LAMBDA[[b]]) %*% t(LAMBDA[[b]])
            }

            # compute VETA
            VETA[[b]] <- Mg %*% (COV - THETA[[b]]) %*% t(Mg)

            # force pd?
            if(local.options[["veta.force.pd"]]) {
                VETA[[b]] <- lav_matrix_symmetric_force_pd( VETA[[b]],
                                 tol = local.options[["veta.force.pd.tol"]] )
            }

            # names
            psi.idx <- which(names(FIT@Model@GLIST) == "psi")[b]
            dimnames(VETA[[b]]) <- FIT@Model@dimNames[[psi.idx]]

            # compute EETA
            if(lavoptions$meanstructure) {
                EETA[[b]] <- Mg %*% (YBAR - NU[[b]])
            }

            # store M
            M[[b]] <- Mg

        } # blocks

        # label groups (if not multilevel)
        # FIXME: we need block.names after all...
        if(ngroups > 1L && nblocks == ngroups)  {
            names(VETA) <- FIT@Data@group.label
        }

        # store EETA/VETA
        out$VETA <- VETA
        out$EETA <- EETA

        # store M
        out$M <- M

    } # local




    ####################################
    # STEP 2: estimate structural part #
    ####################################
    if(sam.method == "local") {
        # extract structural part
        PTS <- lav_partable_subset_structural_model(PT, lavpta = lavpta,
                                                    add.idx = TRUE)
        if(nlevels > 1L) {
            PTS$level <- NULL
            PTS$group <- NULL
            PTS$group <- PTS$block
            NOBS <- FIT@Data@Lp[[1]]$nclusters
        } else {
            NOBS <- FIT@Data@nobs
        }
        # if meanstructure, 'free' user=0 intercepts?
        if(lavoptions$meanstructure) {
            extra.int.idx <- which(PTS$op == "~1" & PTS$user == 0L)
            if(length(extra.int.idx) > 0L) {
                PTS$free[ extra.int.idx ] <- 1L
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

    # adjust options
    lavoptions.PA <- lavoptions
    lavoptions.PA <- modifyList(lavoptions.PA, struc.args)

    # override, not matter what
    lavoptions.PA$do.fit <- TRUE

    if(sam.method == "local") {
        lavoptions.PA$missing <- "listwise"
        lavoptions.PA$se <- "none" # sample statistics input
        lavoptions.PA$sample.cov.rescale <- FALSE
    } else {
        # global (so PA will not be used for for final model)  
        lavoptions.PA$baseline <- FALSE
        lavoptions.PA$h1 <- FALSE
        lavoptions.PA$implied <- FALSE
    }

    # fit structural model
    if(lavoptions.PA$verbose) {
        cat("Fitting Structural Part:\n")
    }
    if(sam.method == "local") {
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
    if(lavoptions.PA$verbose) {
        cat("Done.\n")
    }
    # store FIT.PA
    out$FIT.PA <- FIT.PA

    # fill in point estimates structural part
    PTS <- FIT.PA@ParTable
    p2def.idx <- seq_len(length(PT$lhs)) %in% reg.idx &
                 (PT$free > 0 | PT$op == ":=")

    # which parameters from PTS do we wish to fill in:
    # - all 'free' parameters
    # - := (if any)
    # - but NOT elements in extra.int.idx
    pts.idx <- which( (PTS$free > 0L | PTS$op == ":=") &
                      !seq_len(length(PTS$lhs)) %in% extra.int.idx )

    # find corresponding rows in PT
    PTS2 <- as.data.frame(PTS, stringsAsFactors = FALSE)
    pt.idx <- lav_partable_map_id_p1_in_p2(PTS2[pts.idx,], PT,
                                           exclude.nonpar = FALSE)
    # fill in
    PT$est[ pt.idx ] <- PTS$est[ pts.idx ]


    # fill in any def parameters
    #def.idx <- which(PT$op == ":=")
    #if(length(def.idx) > 0L) {
    #    # call def.function
    #    est.def <- FIT@Model@def.function( PT$est[ PT$free > 0L ] )
    #    if(length(est.def) == length(def.idx)) {
    #        PT$est[ def.idx ] <- est.def
    #    } # else warning?
    #}

    # create step2.idx
    p2.idx <- seq_len(length(PT$lhs)) %in% reg.idx & PT$free > 0 # no def!
    step2.idx <- PT$free[ p2.idx ]



    ###########################
    # compute standard errors #
    ###########################

    # current approach:
    # - create 'global' model, only to get the 'joint' information matrix
    # - partition information matrix (step 1, step 2)
    # - apply two-step correction for second step
    # - 'insert' these corrected SEs (and vcov) in FIT.PA

    lavoptions.joint <- lavoptions
    lavoptions.joint$optim.method = "none"
    PT$ustart <- PT$est # as this is used if optim.method == "none"
    lavoptions.joint$check.gradient = FALSE
    lavoptions.joint$check.start = FALSE
    lavoptions.joint$check.post = FALSE
    #lavoptions.joint$baseline = FALSE
    #lavoptions.joint$h1 = FALSE
    #lavoptions.joint$test <- FIT.PA@Options$test
    #lavoptions.joint$se   <- "none"
    if(sam.method == "local") {
        lavoptions.joint$test <- "none" # local view
    }
    JOINT <- lavaan::lavaan(PT, slotOptions = lavoptions.joint,
                            slotSampleStats = FIT@SampleStats,
                            slotData = FIT@Data)
    # TOTAL information
    INFO <- lavInspect(JOINT, "information")
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
    if(nlevels > 1L) {
        # FIXME: not ok for multigroup multilevel
        N <- FIT@Data@Lp[[1]]$nclusters[[2]] # first group only
    } else {
        N <- nobs(FIT)
    }
    I.22.inv <- 1/N * solve(I.22)

    # method below has the advantage that we can use a 'robust' vcov
    # for the joint model;
    # but does not work if we have equality constraints in the MM!
    # -> D will be singular
    #A <- JOINT@vcov$vcov[ step2.idx,  step2.idx]
    #B <- JOINT@vcov$vcov[ step2.idx, -step2.idx]
    #C <- JOINT@vcov$vcov[-step2.idx,  step2.idx]
    #D <- JOINT@vcov$vcov[-step2.idx, -step2.idx]
    #I.22.inv <- A - B %*% solve(D) %*% C

    # FIXME:
    V2 <- I.22.inv
    #V2 <- JOINT@vcov$vcov[ step2.idx,  step2.idx]

    # V1
    V1 <- (I.22.inv %*% I.21 %*% Sigma.11 %*% I.12 %*% I.22.inv) * N * N

    # V for second step
    VCOV <- V2 + V1

    # store in out
    out$V2 <- V2
    out$V1 <- V1
    out$VCOV <- VCOV



    ##########
    # Output #
    ##########

    # assemble final lavaan objects
    if(sam.method == "local" && output == "lavaan") {
        # overwrite slots in FIT.PA (better way?)
        #PTS$se[ PTS$free > 0L &
        #        !seq_len(length(PTS$lhs)) %in% extra.int.idx ] <-
        #    sqrt( diag(VCOV) )
        PTS$se <- lav_model_vcov_se(lavmodel = FIT.PA@Model,
                                    lavpartable = PTS,
                                    VCOV = VCOV)
        FIT.PA@Options$se <- "twostep"
        FIT.PA@ParTable <- PTS
        FIT.PA@vcov$se <- "twostep"
        FIT.PA@vcov$vcov <- VCOV
        FINAL <- FIT.PA
    } else if(sam.method == "global" && output == "lavaan") {
        #PT$se[ seq_len(length(PT$lhs)) %in% reg.idx &
        #       PT$free > 0L ] <- sqrt( diag(VCOV) )
        FINAL <- JOINT

        if(FINAL@Options$se == "none") {
        } else {
            FINAL@Options$se <- "twostep"
            FINAL@vcov$se <- "twostep"
            FINAL@vcov$vcov[step2.idx, step2.idx] <- VCOV
            PT$se <- lav_model_vcov_se(lavmodel = JOINT@Model,
                                   lavpartable = PT,
                                   VCOV = FINAL@vcov$vcov)
            FINAL@ParTable <- PT
        }
        res <- FINAL
    }


    # prepare output
    if(output == "lavaan") {
        res <- FINAL
    } else {
        res <- out
    }

    res
}



