# factor score regression

# four methods:
#  - naive (regression or Bartlett)
#  - Skrondal & Laake (2001) (regression models only)
#  - Croon (2002) (general + robust SE)
#  - simple: always use Bartlett, replace var(f) by psi estimate
#
# TODO
#  - Hishino & Bentler: this is simple + WLS

# changes 28 nov 2018: add analytic SE ('standad') if fs.method = "Bartlett"
#                      make this the new default

fsr <- function(model      = NULL,
                data       = NULL,
                cmd        = "sem",
                fsr.method = "Croon",
                fs.method  = "Bartlett",
                fs.scores  = FALSE,
                mm.options = list(se = "standard", test = "standard"),
                Gamma.NT   = TRUE,
                lvinfo     = FALSE,
                mm.list    = NULL,
                ...,
                output     = "fsr") {

    # we need full data
    if(is.null(data)) {
        stop("lavaan ERROR: full data is required for factor score regression")
    }

    # check fsr.method argument
    fsr.method <- tolower(fsr.method)
    if(fsr.method == "naive") {
        # nothing to do
    } else if(fsr.method %in% c("skrondal", "laake", "skrondallaake",
                            "skrondal.laake", "skrondal-laake")) {
        fsr.method <- "skrondal.laake"
    } else if(fsr.method == "croon") {
        # nothing to do
    } else if(fsr.method == "simple") {
        # force fs.method to Bartlett!
        fs.method <- "Bartlett"
    } else {
        stop("lavaan ERROR: invalid option for argument fsr.method: ",
             fsr.method)
    }

    # check fs.method argument
    fs.method <- tolower(fs.method)
    if(fs.method %in% c("bartlett", "barttlett", "bartlet")) {
        fs.method <- "Bartlett"
    } else if(fs.method == "regression") {
        # nothing to do
    } else {
        stop("lavaan ERROR: invalid option for argument fs.method: ",
             fs.method)
    }

    if(output %in% c("scores", "fs.scores", "fsr.scores")) {
        fs.scores <- TRUE
    }

    # dot dot dot
    dotdotdot <- list(...)

    # change 'default' values for fsr
    if(is.null(dotdotdot$se)) {
        dotdotdot$se <- "standard"
    }
    if(is.null(dotdotdot$test)) {
        dotdotdot$test <- "satorra.bentler"
    }
    if(is.null(dotdotdot$missing)) {
        dotdotdot$missing <- "ml"
    }
    if(is.null(dotdotdot$meanstructure)) {
        dotdotdot$meanstructure <- TRUE
    }


    # STEP 0: process full model, without fitting
    dotdotdot0 <- dotdotdot
    dotdotdot0$do.fit <- NULL
    dotdotdot0$se <- "none"   # to avoid warning about missing="listwise"
    dotdotdot0$test <- "none" # to avoid warning about missing="listwise"

    # check for arguments that we do not want (eg sample.cov)?
    # TODO

    # initial processing of the model, no fitting
    FIT <- do.call(cmd,
                   args =  c(list(model  = model,
                                  data   = data,
                                  #meanstructure = TRUE,
                                  do.fit = FALSE), dotdotdot0) )
    lavoptions <- lavInspect(FIT, "options")
    # restore
    lavoptions$se   <- dotdotdot$se
    lavoptions$test <- dotdotdot$test
    ngroups     <- lavInspect(FIT, "ngroups")
    lavpta      <- FIT@pta
    lavpartable <- FIT@ParTable

    # FIXME: not ready for multiple groups yet
    if(ngroups > 1L) {
        stop("lavaan ERROR: fsr code not ready for multiple groups (yet)")
    }

    # if missing = "listwise", make data complete
    if(lavoptions$missing == "listwise") {
        # FIXME: make this work for multiple groups!!
        OV <- unique(unlist(lavpta$vnames$ov))
        data <- na.omit(data[,OV])
    }


    # any `regular' latent variables?
    lv.names <- unique(unlist(FIT@pta$vnames$lv.regular))
    ov.names <- unique(unlist(FIT@pta$vnames$ov))

    # check for higher-order factors
    good.idx <- logical( length(lv.names) )
    for(f in seq_len(length(lv.names))) {
        # check the indicators
        FAC <- lv.names[f]
        IND <- lavpartable$rhs[ lavpartable$lhs == FAC &
                                lavpartable$op  == "=~" ]
        if(all(IND %in% ov.names)) {
            good.idx[f] <- TRUE
        }
        # FIXME: check for mixed lv/ov indicators
    }
    lv.names <- lv.names[ good.idx ]

    if(length(lv.names) == 0L) {
        stop("lavaan ERROR: model does not contain any (measured) latent variables")
    }
    nfac     <- length(lv.names)

    # check parameter table
    PT <- parTable(FIT)
    PT$est <- PT$se <- NULL

    # find the structural regressions in the parameter table
    #eqs.idx <- which(PT$op == "~" & (PT$lhs %in% lv.names |
    #                                 PT$rhs %in% lv.names))
    # FIXME: we should allow for just correlations too?
    #if(length(eqs.idx) == 0L) {
    #    stop("lavaan ERROR: regressions do not involve any latent variables")
    #}

    # determine eqs.y and eqs.x names
    eqs.x.names <- unlist(FIT@pta$vnames$eqs.x)
    eqs.y.names <- unlist(FIT@pta$vnames$eqs.y)
    eqs.names <- unique( c(eqs.x.names, eqs.y.names) )

    # check if we can use skrondal & laake (no mediational terms?)
    if(fsr.method == "skrondal.laake") {
        if(any(eqs.x.names %in% eqs.y.names)) {
            stop("lavaan ERROR: mediational relationships are not allowed for the Skrondal.Laake method; use ", sQuote("Croon"), " instead.")
        }
    }


    # STEP 1a: compute factor scores for each measurement model (block)

    # how many measurement models?
    if(!is.null(mm.list)) {

        if(fsr.method != "simple") {
            stop("lavaan ERROR: mm.list only available if fsr.method = \"simple\"")
        }

        nblocks <- length(mm.list)
        # check each measurement block
        for(b in seq_len(nblocks)) {
            if(!all(mm.list[[b]] %in% lv.names)) {
              stop("lavaan ERROR: mm.list contains unknown latent variable(s):",
                paste( mm.list[[b]][ mm.list[[b]] %in% lv.names ], sep = " "),
                "\n")
            }
        }
    } else {
        # TODO: here comes the automatic 'detection' of linked
        #       measurement models
        #
        # for now we take a single latent variable per measurement model block
        mm.list <- as.list(lv.names)
        nblocks <- length(mm.list)
    }

    # compute factor scores, per latent variable
    FS.SCORES  <- vector("list", length = ngroups)
    LVINFO     <- vector("list", length = ngroups)
    if(ngroups > 1L) {
        names(FS.SCORES) <- names(LVINFO) <- lavInspect(FIT, "group.label")
    }
    for(g in 1:ngroups) {
        FS.SCORES[[g]]  <- vector("list", length = nblocks)
        #names(FS.SCORES[[g]]) <-  lv.names
        LVINFO[[g]] <- vector("list", length = nblocks)
        #names(LVINFO[[g]]) <- lv.names
    }

    # adjust options
    dotdotdot2 <- dotdotdot
    dotdotdot2$se <- "none"
    dotdotdot2$test <- "none"
    dotdotdot2$debug <- FALSE
    dotdotdot2$verbose <- FALSE
    dotdotdot2$auto.cov.lv.x <- TRUE # allow correlated exogenous factors

    # override with mm.options
    dotdotdot2 <- modifyList(dotdotdot2, mm.options)

    # we assume the same number/names of lv's per group!!!
    MM.FIT <- vector("list", nblocks)
    for(b in 1:nblocks) {

        # create parameter table for this measurement block only
        PT.block <-
            lav_partable_subset_measurement_model(PT = PT,
                                                  lavpta = lavpta,
                                                  add.lv.cov = TRUE,
                                                  lv.names = mm.list[[b]])
        # fit 1-factor model
        fit.block <- do.call("lavaan",
                            args =  c(list(model  = PT.block,
                                           data   = data), dotdotdot2) )
        # check convergence
        if(!lavInspect(fit.block, "converged")) {
            stop("lavaan ERROR: measurement model for ",
                 paste(mm.list[[b]], collapse = " "), " did not converge.")
        }
        # check fit
        MM.FIT[[b]] <- fit.block

        # fs.method?
        if(fsr.method == "skrondal.laake") {
            # dependent -> Bartlett
            if(lv.names[b] %in% eqs.y.names) {
                fs.method <- "Bartlett"
            } else {
                fs.method <- "regression"
            }
        }

        # compute factor scores
        if(fsr.method %in% c("croon", "simple") ||
           lavoptions$se == "robust.sem") {
            # we use lavPredict() here to remove unwanted dummy lv's, if any
            SC <- lavPredict(fit.block, method = fs.method, fsm = TRUE)
            FSM <- attr(SC, "fsm"); attr(SC, "fsm") <- NULL
            LAMBDA <- computeLAMBDA(fit.block@Model) # FIXME: remove dummy lv's?
            THETA  <- computeTHETA(fit.block@Model)  # FIXME: remove not used ov?
            PSI <- computeVETA(fit.block@Model)
        } else {
            SC <- lavPredict(fit.block, method = fs.method, fsm = FALSE)
        }
        # if ngroups = 1, make list again
        if(ngroups == 1L) {
            # because lavPredict() drops the list
            SC <- list(SC)
        }


        # store results
        for(g in 1:ngroups) {
            FS.SCORES[[g]][[b]] <- SC[[g]]
            if(fsr.method %in% c("croon", "simple")) {
                LVINFO[[g]][[b]] <- list(fsm = FSM[[g]],
                                         lambda = LAMBDA[[g]],
                                         psi    = PSI[[g]],
                                         theta  = THETA[[g]])
            }
        } # g

    } # measurement block

    # cbind factor scores
    FS.SCORES <- lapply(1:ngroups, function(g) {
        SC <- as.data.frame(FS.SCORES[[g]])
        SC
    })

    # compute empirical covariance matrix factor scores
    FS.COV <- lapply(1:ngroups, function(g) {
        COV <- cov(FS.SCORES[[g]]) ## divided by N-1!!!
        if(lavoptions$likelihood == "normal") {
            Ng <- lavInspect(FIT, "nobs")[g]
            COV <- COV * (Ng - 1) / Ng
        }
        COV
    })
    if(lavoptions$meanstructure) {
        FS.MEAN <- lapply(1:ngroups, function(g) { colMeans(FS.SCORES[[g]]) })
    } else {
        FS.MEAN <- NULL
    }

    # STEP 1b: if using `Croon' method: correct COV matrix:
    if(fsr.method %in% c("croon")) {
         FSR.COV <- lav_fsr_croon_correction(FS.COV    = FS.COV,
                                             LVINFO    = LVINFO,
                                             fs.method = fs.method)
    } else if(fsr.method == "simple") {
        FSR.COV <- lav_fsr_simple_correction(FS.COV    = FS.COV,
                                             LVINFO    = LVINFO,
                                             mm.list   = mm.list,
                                             force.pd  = FALSE)
    } else {
        FSR.COV <- FS.COV
    }

    # check if FSR.COV is positive definite for all groups
    for(g in 1:ngroups) {
        txt.group <- ifelse(ngroups > 1L, paste(" in group ", g, sep=""), "")
        eigvals <- eigen(FSR.COV[[g]], symmetric=TRUE, only.values=TRUE)$values
            if(any(eigvals < .Machine$double.eps^(3/4))) {
                stop(
"lavaan ERROR: corrected covariance matrix of factor scores\n",
"                is not positive definite", txt.group, ";\n")
            }
    }


    # STEP 1c: do we need full set of factor scores?
    if(fs.scores) {
        # transform?
        if(fsr.method %in% c("croon", "simple")) {
            for(g in 1:ngroups) {
                OLD.inv <- solve(FS.COV[[g]])
                OLD.inv.sqrt <- lav_matrix_symmetric_sqrt(OLD.inv)
                FSR.COV.sqrt <- lav_matrix_symmetric_sqrt(FSR.COV[[g]])
                SC <- as.matrix(FS.SCORES[[g]])
                SC <- SC %*% OLD.inv.sqrt %*% FSR.COV.sqrt
                SC <- as.data.frame(SC)
                names(SC) <- lv.names
                FS.SCORES[[g]] <- SC
            }
        }

        # unlist if multiple groups, add group column
        if(ngroups == 1L) {
            FS.SCORES <- as.data.frame(FS.SCORES[[1]])
        } else {
            stop("fix this!")
        }
    }




    # STEP 2: fit structural model using (corrected?) factor scores

    PT.PA <- lav_partable_subset_structural_model(PT, lavpta = lavpta)

    # free all means/intercepts (of observed variables only)
    lv.names.pa <- lavNames(PT.PA, "lv")
    int.idx <- which(PT.PA$op == "~1" & !PT.PA$lhs %in% lv.names.pa)
    PT.PA$free[int.idx] <- 1L
    PT.PA$ustart[int.idx] <- NA
    # adjust lavoptions

    #if(is.null(dotdotdot$missing)) {
    #    lavoptions$missing <- "listwise" # factor scores are always complete
    #} else {
    #    lavoptions$missing <- dotdotdot$missing
    #}
    #if(is.null(dotdotdot$information)) {
    #    lavoptions$information <- "expected"
    #} else {
    #    lavoptions$information <- dotdotdot$information
    #}
    if(is.null(dotdotdot$do.fit)) {
        lavoptions$do.fit <- TRUE
    } else {
        lavoptions$do.fit <- dotdotdot$do.fit
    }
    if(is.null(dotdotdot$se)) {
        lavoptions$se <- "robust.sem"
    } else {
        lavoptions$se <- dotdotdot$se
    }
    if(is.null(dotdotdot$test)) {
        lavoptions$test <- "satorra.bentler"
    } else {
        lavoptions$test <- dotdotdot$test
    }
    if(is.null(dotdotdot$sample.cov.rescale)) {
        lavoptions$sample.cov.rescale <- FALSE
    } else {
        lavoptions$sample.cov.rescale <- dotdotdot$sample.cov.rescale
    }

    # take care of NACOV, in case we want correct standard errors
    if(lavoptions$se == "robust.sem") {
        Omega.f <- vector("list", length = ngroups)
        for(g in 1:ngroups) {
            DATA <- FIT@Data@X[[g]]
            if(Gamma.NT) {
                if(lavoptions$missing == "listwise") {
                    Omega.y <- lav_samplestats_Gamma_NT(Y = DATA,
                                   meanstructure = lavoptions$meanstructure,
                                   rescale = TRUE, fixed.x = FALSE)
                } else if(lavoptions$missing %in% c("ml", "ml.x")) {
                    # we assume UNSTRUCTURED Mu and Sigma!!
                    MU    <- FIT@SampleStats@missing.h1[[g]]$mu
                    SIGMA <- FIT@SampleStats@missing.h1[[g]]$sigma
                    if(lavoptions$information == "expected") {
                        Info <- lav_mvnorm_missing_information_expected(
                                    Y = DATA, Mp = FIT@Data@Mp[[g]],
                                    wt = FIT@Data@weights[[g]],
                                    Mu = MU, Sigma = SIGMA,
                                    x.idx = FIT@SampleStats@x.idx[[g]])
                    } else {
                        Info <-
                            lav_mvnorm_missing_information_observed_samplestats(
                                Yp = FIT@SampleStats@missing[[g]],
                                # wt not needed
                                Mu = MU, Sigma = SIGMA,
                                x.idx = FIT@SampleStats@x.idx[[g]])
                    }
                    Omega.y <- lav_matrix_symmetric_inverse(Info)
                } else {
                    stop("lavaan ERROR: can not handle missing = ",
                         lavoptions$missing)
                }

            } else {
                if(lavoptions$missing == "listwise") {
                    Omega.y <- lav_samplestats_Gamma(Y = DATA,
                                   meanstructure = lavoptions$meanstructure,
                                   fixed.x = FALSE)
                } else if(lavoptions$missing %in% c("ml", "ml.x")) {
                    # we assume UNSTRUCTURED Mu and Sigma!!
                    MU    <- FIT@SampleStats@missing.h1[[g]]$mu
                    SIGMA <- FIT@SampleStats@missing.h1[[g]]$sigma
                    Omega.y <- lav_mvnorm_missing_h1_omega_sw(Y =
                             DATA, Mp = FIT@Data@Mp[[g]],
                             Yp = FIT@SampleStats@missing[[g]],
                             Mu = MU, Sigma = SIGMA,
                             x.idx = FIT@SampleStats@x.idx[[g]],
                             information = lavoptions$information)
                } else {
                    stop("lavaan ERROR: can not handle missing = ",
                         lavoptions$missing)
                }
            }

            # factor score matrices
            A <- lav_matrix_bdiag(lapply(LVINFO[[g]], "[[", "fsm"))

            # compensate for Croon correction
            if(fs.method == "regression") {
                if(!exists("OLD.inv.sqrt")) {
                    OLD.inv <- solve(FS.COV[[g]])
                    OLD.inv.sqrt <- lav_matrix_symmetric_sqrt(OLD.inv)
                }
                if(!exists("FSR.COV.sqrt")) {
                    FSR.COV.sqrt <- lav_matrix_symmetric_sqrt(FSR.COV[[g]])
                }
                A <- OLD.inv.sqrt %*% FSR.COV.sqrt %*% A
            }

            # mean + vech(sigma)
            A22 <- lav_matrix_duplication_post(
                   lav_matrix_duplication_ginv_pre(A %x% A))
            if(lavoptions$meanstructure) {
                A11 <- A
                A.tilde <- lav_matrix_bdiag(A11, A22)
            } else {
                A.tilde <- A22
            }
            Omega.f[[g]] <- A.tilde %*% Omega.y %*% t(A.tilde)
        } # g
    } else {
        Omega.f <- NULL
    }


    # fit structural model
    lavoptions2 <- lavoptions
    #if(lavoptions$se == "standard") {
    #    lavoptions2$se <- "external"
    #}
    #lavoptions2$test <- "none"
    lavoptions2$missing <- "listwise" # always complete data anyway...
    fit <- lavaan(PT.PA,
                  sample.cov = FSR.COV,
                  sample.mean = FS.MEAN,
                  sample.nobs = FIT@SampleStats@nobs,
                  NACOV       = Omega.f,
                  slotOptions = lavoptions2)

    # extra info
    extra <- list( FS.COV =  FS.COV,  FS.SCORES =  FS.SCORES,
                   FSR.COV = FSR.COV,
                   LVINFO = LVINFO)

    # standard errors
    #lavsamplestats <- fit@SampleStats
    #lavsamplestats@NACOV <- Omega.f
    #VCOV <- lav_model_vcov(fit@Model, lavsamplestats = lavsamplestats,
    #                       lavoptions = lavoptions)
    #SE <- lav_model_vcov_se(fit@Model, fit@ParTable, VCOV = VCOV)
    #PE$se <- SE
    #tmp.se <- ifelse(PE$se == 0.0, NA, PE$se)
    #zstat <- pvalue <- TRUE
    #if(zstat) {
    #    PE$z <- PE$est / tmp.se
    #    if(pvalue) {
    #        PE$pvalue <- 2 * (1 - pnorm( abs(PE$z) ))
    #    }
    #}

    if(output == "fsr") {
        #PE <- parameterEstimates(fit, add.attributes = TRUE, ci = FALSE)
        HEADER <- paste("This is fsr (0.2) -- factor score regression using ",
                        "fsr.method = ", fsr.method, sep = "")
        out <- list(header = HEADER, MM.FIT = MM.FIT, STRUC.FIT = fit)
        if(lvinfo) {
            out$lvinfo <- extra
        }

        class(out) <- c("lavaan.fsr", "list")

    } else if(output %in% c("lavaan", "fit")) {
        out <- fit
    } else if(output == "extra") {
        out <- extra
    } else if(output == "lvinfo") {
        out <- LVINFO
    } else if(output %in% c("scores", "f.scores", "fs.scores")) {
        out <- FS.SCORES
    } else if(output %in% c("FSR.COV", "fsr.cov", "croon", "cov.croon",
                            "croon.cov", "COV", "cov")) {
        out <- FSR.COV
    } else if(output %in% c("FS.COV", "fs.cov")) {
        out <- FS.COV
    } else {
        stop("lavaan ERROR: unknown output= argument: ", output)
    }

    out
}



