# bootstrap based NVCOV
lav_model_nvcov_bootstrap <- function(lavmodel = NULL, lavsamplestats = NULL, 
                                      lavoptions = NULL, lavdata = NULL,
                                      lavcache = NULL, lavpartable = NULL, 
                                      control=list()) {

    # number of bootstrap draws
    if(!is.null(lavoptions$bootstrap)) {
        R <- lavoptions$bootstrap
    } else {
        R <- 1000L
    }
  
    boot.type <- "ordinary"
    if(lavoptions$test == "bollen.stine") boot.type <- "bollen.stine"

    TEST <- NULL
    COEF <- bootstrap.internal(object          = NULL,
                               lavmodel.       = lavmodel, 
                               lavsamplestats. = lavsamplestats, 
                               lavpartable.    = lavpartable, 
                               lavoptions.     = lavoptions, 
                               lavdata.        = lavdata,
                               R               = R, 
                               verbose         = lavoptions$verbose,
                               type            = boot.type,
                               FUN  = ifelse(boot.type == "bollen.stine",
                                          "coeftest", "coef"),
                               warn            = -1L,
                               parallel        = control$parallel,
                               ncpus           = control$ncpus,
                               cl              = control$cl)
    if(boot.type == "bollen.stine") {
        nc <- ncol(COEF)
        TEST <- COEF[,nc]
        COEF <- COEF[,-nc]
    }

    # FIXME: cov rescale? Yes for now
    nboot <- nrow(COEF)
    NVarCov <- lavsamplestats@ntotal * (cov(COEF) * (nboot-1)/nboot )

    # save COEF and TEST (if any)
    attr(NVarCov, "BOOT.COEF") <- COEF
    attr(NVarCov, "BOOT.TEST") <- TEST
 
    NVarCov
}


# robust `sem' NVCOV (see Browne, 1984,  bentler & dijkstra 1985)
lav_model_nvcov_robust_sem <- function(lavmodel = NULL, lavsamplestats = NULL,
                                       lavdata = NULL, lavcache = NULL,
                                       estimator = "ML", mimic = "lavaan") {

    # compute inverse of the expected(!) information matrix
    if(estimator == "ML" && mimic == "Mplus") {
        # YR - 11 aug 2010 - what Mplus seems to do is (see Muthen apx 4 eq102)
        # - WLS.V is not based on Sigma.hat and Mu.hat (as it
        #   should be?), but on lavsamplestats@cov and lavsamplestats@mean...
        # - Gamma is not identical to what is used for WLS; closer to EQS
        # - N/N-1 bug in G11 for NVarCov (but not test statistic)
        # - we divide by N-1! (just like EQS)
        E.inv <- lav_model_information_expected_MLM(lavmodel = lavmodel, 
                                           lavsamplestats = lavsamplestats,
                                           extra          = TRUE,
                                           augmented      = TRUE,
                                           inverted       = TRUE)
    } else {
        E.inv <- lav_model_information_expected(lavmodel = lavmodel, 
                                        lavsamplestats = lavsamplestats, 
                                        lavdata        = lavdata, 
                                        estimator      = estimator, 
                                        extra          = TRUE,
                                        augmented      = TRUE,
                                        inverted       = TRUE)
    }

    # check if E.inv is ok
    if(inherits(E.inv, "try-error")) { 
        return(E.inv)
    }

    Delta <- attr(E.inv, "Delta")
    WLS.V <- attr(E.inv, "WLS.V")

    # Gamma
    Gamma <- lavsamplestats@NACOV
    if(estimator == "ML" && mimic == "Mplus") {
        # 'fix' G11 part of Gamma (NOTE: this is NOT needed for SB test 
        # statistic
        for(g in 1:lavsamplestats@ngroups) {
            gg1 <- (lavsamplestats@nobs[[g]]-1)/lavsamplestats@nobs[[g]]
            G11 <- gg1 * lavsamplestats@cov[[g]]
            Gamma[[g]][1:nrow(G11), 1:nrow(G11)] <- G11
        } # g
    }
   

    tDVGVD <- matrix(0, ncol=ncol(E.inv), nrow=nrow(E.inv))
    for(g in 1:lavsamplestats@ngroups) {
        fg  <-  lavsamplestats@nobs[[g]]   /lavsamplestats@ntotal
        fg1 <- (lavsamplestats@nobs[[g]]-1)/lavsamplestats@ntotal
        # fg twice for WLS.V, 1/fg1 once for GaMMA
        # if fg==fg1, there would be only one fg, as in Satorra 1999 p.8
        # t(Delta) * WLS.V %*% Gamma %*% WLS.V %*% Delta
        if(estimator == "DWLS" || estimator == "ULS") {
            # diagonal weight matrix
            WD <- WLS.V[[g]] * Delta[[g]]
        } else {
            # full weight matrix
            WD <- WLS.V[[g]] %*% Delta[[g]]
        }
        tDVGVD <- tDVGVD + fg*fg/fg1 * crossprod(WD, Gamma[[g]] %*% WD)
    } # g
    NVarCov <- (E.inv %*% tDVGVD %*% E.inv)

    # to be reused by lavaanTest
    attr(NVarCov, "E.inv") <- E.inv
    attr(NVarCov, "Delta") <- Delta
    attr(NVarCov, "WLS.V") <- WLS.V

    NVarCov
}

lav_model_nvcov_robust_sandwich <- function(lavmodel      = lavmodel, 
                                           lavsamplestats = NULL, 
                                           lavdata        = NULL,
                                           information    = "observed", 
                                           lavcache       = NULL, 
                                           estimator      = "ML") {

    # sandwich estimator: A.inv %*% B %*% t(A.inv)
    # where A.inv == E.inv
    #       B == outer product of case-wise scores

    # inverse observed/expected information matrix
    E.inv <- lav_model_information(lavmodel       = lavmodel,
                                   lavsamplestats = lavsamplestats,
                                   lavdata        = lavdata,
                                   estimator      = estimator,
                                   lavcache       = lavcache,
                                   information    = information,
                                   extra          = FALSE,
                                   augmented      = TRUE,
                                   inverted       = TRUE)

    # check if E.inv is ok
    if(inherits(E.inv, "try-error")) {
        return(E.inv)
    }

    # outer product of case-wise scores
    B0 <- 
        lav_model_information_firstorder(lavmodel       = lavmodel,
                                         lavsamplestats = lavsamplestats,
                                         lavdata        = lavdata,
                                         estimator      = estimator,
                                         lavcache       = lavcache,
                                         extra          = TRUE,
                                         check.pd       = FALSE,
                                         augmented      = FALSE,
                                         inverted       = FALSE)

    # compute sandwich estimator
    NVarCov <- E.inv %*% B0 %*% E.inv

    attr(NVarCov, "B0.group") <- attr(B0, "B0.group")
    attr(NVarCov, "E.inv") <- E.inv

    NVarCov
}


lav_model_vcov <- function(lavmodel       = NULL, 
                           lavsamplestats = NULL, 
                           lavoptions     = NULL, 
                           lavdata        = NULL, 
                           lavpartable    = NULL, 
                           lavcache       = NULL, 
                           control=list()) {

    estimator   <- lavoptions$estimator
    likelihood  <- lavoptions$likelihood
    information <- lavoptions$information
    se          <- lavoptions$se
    verbose     <- lavoptions$verbose
    mimic       <- lavoptions$mimic
  
    if(se == "none") return(matrix(0,0,0))

    # some require meanstructure (for now)
    if(se %in% c("first.order", "robust.sem", "robust.huber.white") && 
       !lavoptions$meanstructure) {
        stop("se (", se, ") requires meanstructure (for now)")
    }

    if(se == "standard") {
        NVarCov <- lav_model_information(lavmodel       = lavmodel,
                                         lavsamplestats = lavsamplestats,
                                         lavdata        = lavdata,
                                         estimator      = estimator,
                                         lavcache       = lavcache,
                                         information    = information,
                                         extra          = FALSE,
                                         augmented      = TRUE,
                                         inverted       = TRUE)

    } else if(se == "first.order") {
        NVarCov <- 
            lav_model_information_firstorder(lavmodel = lavmodel,
                                             lavsamplestats = lavsamplestats,
                                             lavdata        = lavdata,
                                             estimator      = estimator,
                                             lavcache       = lavcache,
                                             extra          = TRUE,
                                             check.pd       = FALSE,
                                             augmented      = TRUE,
                                             inverted       = TRUE)
    
    } else if(se == "robust.sem") {
        NVarCov <-
            lav_model_nvcov_robust_sem(lavmodel       = lavmodel,
                                       lavsamplestats = lavsamplestats,
                                       estimator      = estimator,
                                       mimic          = mimic,
                                       lavcache       = lavcache,
                                       lavdata        = lavdata)

    } else if(se == "robust.huber.white") {
        NVarCov <-
            lav_model_nvcov_robust_sandwich(lavmodel = lavmodel,
                                            lavsamplestats = lavsamplestats,
                                            lavdata        = lavdata,
                                            information    = information,
                                            lavcache       = lavcache,
                                            estimator      = estimator)

    } else if(se == "bootstrap") {
        NVarCov <- try( lav_model_nvcov_bootstrap(lavmodel       = lavmodel,
                                        lavsamplestats = lavsamplestats,
                                        lavoptions     = lavoptions,
                                        lavdata        = lavdata,
                                        lavcache       = lavcache,
                                        lavpartable    = lavpartable,
                                        control        = control),
                        silent=TRUE )
    }

    if(! inherits(NVarCov, "try-error") ) {

        # denominator!
        if(estimator %in% c("ML","PML","FML") && likelihood == "normal") {
            N <- lavsamplestats@ntotal
        } else {
            N <- lavsamplestats@ntotal - lavsamplestats@ngroups
        }

        #if(estimator %in% c("PML", "MML")) {
        #    VarCov <- NVarCov
        #} else {
            VarCov <- 1/N * NVarCov
        #}

    } else {
        warning("lavaan WARNING: could not compute standard errors!\n  lavaan NOTE: this may be a symptom that the model is not identified.\n")
        VarCov <- NULL
    }

    VarCov
}


