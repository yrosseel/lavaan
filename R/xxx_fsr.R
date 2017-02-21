# factor score regression

# three methods:
#  - naive (regression or Bartlett)
#  - Skrondal & Laake (2001) (regression models only)
#  - Croon (2002) (general + robust SE)

fsr <- function(model      = NULL, 
                data       = NULL, 
                cmd        = "sem",
                fsr.method = "Croon", 
                fs.method  = "Bartlett", 
                fs.scores  = FALSE,
                Gamma.NT   = TRUE,
                lvinfo     = FALSE,
                ...) {
   
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
    
    # dot dot dot
    dotdotdot <- list(...)

    # change 'default' values for fsr
    if(is.null(dotdotdot$se)) {
        dotdotdot$se <- "none"
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
    ngroups    <- lavInspect(FIT, "ngroups")
    lavpta     <- FIT@pta

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
    if(length(lv.names) == 0L) {
        stop("lavaan ERROR: model does not contain any latent variables")
    }
    nfac     <- length(lv.names)

    # check parameter table
    PT <- parTable(FIT)
    PT$est <- PT$se <- NULL

    # find the structural regressions in the parameter table
    eqs.idx <- which(PT$op == "~" & (PT$lhs %in% lv.names |
                                     PT$rhs %in% lv.names))
    if(length(eqs.idx) == 0L) {
        stop("lavaan ERROR: regressions do not involve any latent variables")
    }

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


    # STEP 1a: compute factor scores for each latent variable

    # compute factor scores, per latent variable
    FS.SCORES  <- vector("list", length = ngroups)
    LVINFO     <- vector("list", length = ngroups)
    if(ngroups > 1L) {
        names(FS.SCORES) <- names(LVINFO) <- lavInspect(FIT, "group.label")
    }
    for(g in 1:ngroups) {
        FS.SCORES[[g]]  <- vector("list", length = nfac)
        names(FS.SCORES[[g]]) <-  lv.names
        LVINFO[[g]] <- vector("list", length = nfac)
        names(LVINFO[[g]]) <- lv.names
    }

    # adjust options
    dotdotdot2 <- dotdotdot
    dotdotdot2$se <- "none"
    dotdotdot2$test <- "none"
    dotdotdot2$debug <- FALSE
    dotdotdot2$verbose <- FALSE
    dotdotdot2$auto.cov.lv.x <- TRUE # allow correlated exogenous factors
 
    # we assume the same number/names of lv's per group!!!
    for(f in 1:nfac) {

        # create parameter table for this factor only
        PT.1fac <- lav_partable_subset_measurement_model(PT = PT,
                                                         lavpta = lavpta,
                                                         lv.names = lv.names[f])
        # fit 1-factor model
        fit.1fac <- do.call("lavaan",
                            args =  c(list(model  = PT.1fac,
                                           data   = data), dotdotdot2) )

        # fs.method?
        if(fsr.method == "skrondal.laake") {
            # dependent -> Bartlett
            if(lv.names[f] %in% eqs.y.names) {
                fs.method <- "Bartlett"
            } else {
                fs.method <- "regression"
            }
        }

        # compute factor scores
        if(fsr.method %in% c("croon") || 
           lavoptions$se == "robust.sem") {
            SC <- lav_predict_eta(fit.1fac, method = fs.method, fsm = TRUE)
            FSM <- attr(SC, "fsm"); attr(SC, "fsm") <- NULL
            LAMBDA <- computeLAMBDA(fit.1fac@Model)
            THETA  <- computeTHETA(fit.1fac@Model)
        } else {
            SC <- lav_predict_eta(fit.1fac, method = fs.method, fsm = FALSE)
        }

        # store results
        for(g in 1:ngroups) {
            FS.SCORES[[g]][[f]] <- SC[[g]]
            if(fsr.method %in% c("croon") ||
               lavoptions$se == "robust.sem") {
                LVINFO[[g]][[f]] <- list(fsm = FSM[[g]], lambda = LAMBDA[[g]],
                                                         theta  = THETA[[g]])
            }
        } # g

    } # nfac


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
    } else {
        FSR.COV <- FS.COV 
    }

    # STEP 1c: do we need full set of factor scores?
    if(fs.scores) {
        # transform?
        if(fsr.method == "croon") {
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

    # free all means/intercepts
    int.idx <- which(PT.PA$op == "~1")
    PT.PA$free[int.idx] <- 1L
    PT.PA$ustart[int.idx] <- NA

    # adjust lavoptions
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
                } else if(lavoptions$missing == "ml") {
                    # we assume UNSTRUCTURED Mu and Sigma!!
                    MU    <- FIT@SampleStats@missing.h1[[g]]$mu
                    SIGMA <- FIT@SampleStats@missing.h1[[g]]$sigma
                    if(lavoptions$information == "expected") {
                        Info <- lav_mvnorm_missing_information_expected(
                                    Y = DATA, Mp = FIT@Data@Mp[[g]],
                                    Mu = MU, Sigma = SIGMA)
                    } else {
                        Info <- lav_mvnorm_missing_information_observed_samplestats(
                                Yp = FIT@SampleStats@missing[[g]],
                                Mu = MU, Sigma = SIGMA)
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
                } else if(lavoptions$missing == "ml") {
                    # we assume UNSTRUCTURED Mu and Sigma!!
                    MU    <- FIT@SampleStats@missing.h1[[g]]$mu
                    SIGMA <- FIT@SampleStats@missing.h1[[g]]$sigma
                    Omega.y <- lav_mvnorm_missing_h1_omega_sw(Y =
                             DATA, Mp = FIT@Data@Mp[[g]],
                             Yp = FIT@SampleStats@missing[[g]],
                             Mu = MU, Sigma = SIGMA,
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
    #lavoptions2$se <- "none"
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

    PE <- parameterEstimates(fit, add.attributes = TRUE)

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
   
    out <- list(header = "This is fsr (0.1) -- factor score regression.",
                PE = PE)

    if(lvinfo) {
        out$lvinfo <- extra
    }

    class(out) <- c("lavaan.fsr", "list")

    out
}



