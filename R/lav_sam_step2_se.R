# compute two-step standard errors for SAM models

lav_sam_step2_se <- function(FIT = NULL, JOINT = NULL,
                             STEP1 = NULL, STEP2 = NULL,
                             local.options = list()) {

    # current approach for se = "twostep":
    # - create 'global' model, only to get the 'joint' information matrix
    # - partition information matrix (step 1, step 2)
    # - apply two-step correction for second step
    # - 'insert' these corrected SEs (and vcov) in FIT.PA

    out <- list()
    Sigma.11 <- STEP1$Sigma.11
    step1.free.idx <- STEP1$step1.free.idx
    step2.free.idx <- STEP2$step2.free.idx
    lavoptions <- FIT@Options
    nlevels <- FIT@pta$nlevels
    FIT.PA <- STEP2$FIT.PA
    extra.int.idx <- STEP2$extra.int.idx

    # catch empty step2.free.idx
    if(length(step2.free.idx) == 0L) {
        # no (free) structural parameters at all!
        out <- list(V1 = matrix(0, 0, 0), V2 = matrix(0, 0, 0),
                    VCOV = matrix(0, 0, 0))
        return(out)
    }

    if(lavoptions$se == "none") {
        # nothing to do...
    } else {
        if(lavoptions$verbose) {
            cat("Computing ", lavoptions$se, " standard errors ... ", sep = "")
        }

        INFO <- lavInspect(JOINT, "information")
        I.12 <- INFO[step1.free.idx, step2.free.idx]
        I.22 <- INFO[step2.free.idx, step2.free.idx]
        I.21 <- INFO[step2.free.idx, step1.free.idx]

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
        if(inherits(I.22.inv, "try-error")) {
            # hm, not good
            if(lavoptions$se != "naive") {
                warning("lavaan WARNING: problem inverting information matrix (I.22);\n\t\t  -> switching to naive standard errors!")
                lavoptions$se <- "naive"
            }
        }

        # method below has the advantage that we can use a 'robust' vcov
        # for the joint model;
        # but does not work if we have equality constraints in the MM!
        # -> D will be singular
        #A <- JOINT@vcov$vcov[ step2.free.idx,  step2.free.idx]
        #B <- JOINT@vcov$vcov[ step2.free.idx, -step2.free.idx]
        #C <- JOINT@vcov$vcov[-step2.free.idx,  step2.free.idx]
        #D <- JOINT@vcov$vcov[-step2.free.idx, -step2.free.idx]
        #I.22.inv <- A - B %*% solve(D) %*% C

        if(lavoptions$se == "standard") {
            VCOV <- 1/N * I.22.inv
            out$VCOV <- VCOV
        } else if(lavoptions$se == "naive") {
            if(is.null(FIT.PA@vcov$vcov)) {
                FIT.PA@Options$se <- "standard"
                VCOV.naive <- lavTech(FIT.PA, "vcov")
            } else {
                VCOV.naive <- FIT.PA@vcov$vcov
            }
            if(length(extra.int.idx) > 0L) {
                rm.idx <- FIT.PA@ParTable$free[extra.int.idx]
                VCOV.naive <- VCOV.naive[-rm.idx, -rm.idx]
            }
            out$VCOV <- VCOV.naive
        } else {
           # twostep

            # FIXME:
            V2 <- 1/N * I.22.inv # not the same as FIT.PA@vcov$vcov!!
            #V2 <- JOINT@vcov$vcov[ step2.free.idx,  step2.free.idx]
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
            if(!is.null(local.options$alpha.correction) &&
               local.options$alpha.correction > 0) {
                alpha.N1 <- local.options$alpha.correction / (N - 1)
                if(alpha.N1 > 1.0) {
                    alpha.N1 <- 1.0
                } else if (alpha.N1 < 0.0) {
                    alpha.N1 <- 0.0
                }
                if(is.null(FIT.PA@vcov$vcov)) {
                    FIT.PA@Options$se <- "standard"
                    VCOV.naive <- lavTech(FIT.PA, "vcov")
                } else {
                    VCOV.naive <- FIT.PA@vcov$vcov
                }
                if(length(extra.int.idx) > 0L) {
                     rm.idx <- FIT.PA@ParTable$free[extra.int.idx]
                    VCOV.naive <- VCOV.naive[-rm.idx, -rm.idx]
                }
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

    out
}
