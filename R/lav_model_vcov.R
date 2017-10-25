# bootstrap based NVCOV
lav_model_nvcov_bootstrap <- function(lavmodel       = NULL, 
                                      lavsamplestats = NULL, 
                                      lavoptions     = NULL, 
                                      lavimplied     = NULL,
                                      lavdata        = NULL,
                                      lavcache       = NULL, 
                                      lavpartable    = NULL) {

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
                               warn            = -1L)
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
lav_model_nvcov_robust_sem <- function(lavmodel       = NULL, 
                                       lavsamplestats = NULL,
                                       lavdata        = NULL, 
                                       lavcache       = NULL,
                                       lavimplied     = NULL,
                                       lavoptions     = NULL,
                                       use.ginv       = FALSE) {

    # compute inverse of the expected(!) information matrix
    if(lavmodel@estimator == "ML" && lavoptions$mimic == "Mplus") {
        # YR - 11 aug 2010 - what Mplus seems to do is (see Muthen apx 4 eq102)
        # - A1 is not based on Sigma.hat and Mu.hat, 
        # but on lavsamplestats@cov and lavsamplestats@mean... ('unstructured')
        # - Gamma is not identical to what is used for WLS; closer to EQS
        # - N/N-1 bug in G11 for NVarCov (but not test statistic)
        # - we divide by N-1! (just like EQS)
        E.inv <- lav_model_information_expected_MLM(lavmodel = lavmodel, 
                                           lavsamplestats = lavsamplestats,
                                           extra          = TRUE,
                                           augmented      = TRUE,
                                           inverted       = TRUE,
                                           use.ginv       = use.ginv)
    } else {
        E.inv <- lav_model_information(lavmodel       = lavmodel,
                                       lavsamplestats = lavsamplestats,
                                       lavdata        = lavdata,
                                       lavimplied     = lavimplied,
                                       lavoptions     = lavoptions,
                                       extra          = TRUE,
                                       augmented      = TRUE,
                                       inverted       = TRUE,
                                       use.ginv       = use.ginv)
    }

    # check if E.inv is ok
    if(inherits(E.inv, "try-error")) { 
        return(E.inv)
    }

    Delta <- attr(E.inv, "Delta")
    WLS.V <- attr(E.inv, "WLS.V")

    # Gamma
    Gamma <- lavsamplestats@NACOV
    if(lavmodel@estimator == "ML" && 
       lavoptions$mimic == "Mplus" && !lavsamplestats@NACOV.user) {
        # 'fix' G11 part of Gamma (NOTE: this is NOT needed for SB test 
        # statistic
        for(g in 1:lavsamplestats@ngroups) {
            gg1 <- (lavsamplestats@nobs[[g]]-1)/lavsamplestats@nobs[[g]]
            if(lavmodel@conditional.x) {
                nvar <- NCOL(lavsamplestats@res.cov[[g]])
            } else {
                nvar <- NCOL(lavsamplestats@cov[[g]])
            }
            G11 <- Gamma[[g]][1:nvar, 1:nvar, drop = FALSE]
            Gamma[[g]][1:nvar, 1:nvar] <- G11 * gg1
        } # g
    }
   

    tDVGVD <- matrix(0, ncol=ncol(E.inv), nrow=nrow(E.inv))
    for(g in 1:lavsamplestats@ngroups) {
        fg  <-  lavsamplestats@nobs[[g]]   /lavsamplestats@ntotal
        if(lavoptions$mimic == "Mplus") {
            fg1 <- (lavsamplestats@nobs[[g]]-1)/lavsamplestats@ntotal
        } else {
            # from 0.6 onwards, we use fg1 == fg, to be more consistent with
            # lav_test()
            fg1 <- fg
        }
        # fg twice for WLS.V, 1/fg1 once for GaMMA
        # if fg==fg1, there would be only one fg, as in Satorra 1999 p.8
        # t(Delta) * WLS.V %*% Gamma %*% WLS.V %*% Delta
        if(lavmodel@estimator == "DWLS" || lavmodel@estimator == "ULS") {
            # diagonal weight matrix
            WD <- WLS.V[[g]] * Delta[[g]]
        } else {
            # full weight matrix
            WD <- WLS.V[[g]] %*% Delta[[g]]
        }
        tDVGVD <- tDVGVD + fg*fg/fg1 * crossprod(WD, Gamma[[g]] %*% WD)
    } # g
    NVarCov <- (E.inv %*% tDVGVD %*% E.inv)

    # to be reused by lav_test()
    attr(NVarCov, "E.inv") <- E.inv
    attr(NVarCov, "Delta") <- Delta
    attr(NVarCov, "WLS.V") <- WLS.V

    NVarCov
}

lav_model_nvcov_robust_sandwich <- function(lavmodel       = NULL,
                                            lavsamplestats = NULL, 
                                            lavdata        = NULL,
                                            lavoptions     = NULL,
                                            lavimplied     = NULL,
                                            lavcache       = NULL, 
                                            use.ginv       = FALSE) {

    # sandwich estimator: A.inv %*% B %*% t(A.inv)
    # where A.inv == E.inv
    #       B == outer product of case-wise scores

    # inverse observed/expected information matrix
    E.inv <- lav_model_information(lavmodel       = lavmodel,
                                   lavsamplestats = lavsamplestats,
                                   lavdata        = lavdata,
                                   lavcache       = lavcache,
                                   lavimplied     = lavimplied,
                                   lavoptions     = lavoptions,
                                   extra          = FALSE,
                                   augmented      = TRUE,
                                   inverted       = TRUE,
                                   use.ginv       = use.ginv)

    # check if E.inv is ok
    if(inherits(E.inv, "try-error")) {
        return(E.inv)
    }

    # outer product of case-wise scores
    B0 <- 
        lav_model_information_firstorder(lavmodel       = lavmodel,
                                         lavsamplestats = lavsamplestats,
                                         lavdata        = lavdata,
                                         lavcache       = lavcache,
                                         lavimplied     = lavimplied,
                                         lavoptions     = lavoptions,
                                         extra          = TRUE,
                                         check.pd       = FALSE,
                                         augmented      = FALSE,
                                         inverted       = FALSE,
                                         use.ginv       = use.ginv)

    # compute sandwich estimator
    NVarCov <- E.inv %*% B0 %*% E.inv

    attr(NVarCov, "B0.group") <- attr(B0, "B0.group")
    attr(NVarCov, "E.inv") <- E.inv

    NVarCov
}

# two stage 
# - two.stage: Gamma = I_1^{-1}
# - robust.two.stage: Gamma = incomplete Gamma (I_1^{-1} J_1 I_1^{-1})
# where I_1 and J_1 are based on the (saturated) model h1 
# (either unstructured, or structured)
#
# references:
#
# - Savalei \& Bentler (2009) eq (6) for se = "two.stage"
# - Savalei \& Falk (2014) eq  (3)   for se = "robust.two.stage"
# - Yuan \& Bentler (2000)
lav_model_nvcov_two_stage <- function(lavmodel       = NULL, 
                                      lavsamplestats = NULL,
                                      lavoptions     = NULL,
                                      lavimplied     = NULL,
                                      lavdata        = NULL,
                                      use.ginv       = FALSE) {

    # expected OR observed, depending on lavoptions$information
    if(is.null(lavoptions) && is.null(lavoptions$information)) {
        lavoptions <- list(information = "observed",
                           observed.information = "h1",
                           h1.information = "structured")
    }

    # restrictions:
    # only works if:
    # - information is expected,
    # - or information is observed but with observed.information == "h1"
    if(lavoptions$information == "observed" && 
       lavoptions$observed.information != "h1") {
            stop("lavaan ERROR: two.stage + observed information currently only works with observed.information = ", dQuote("h1"))
    }
    # no weights (yet)
    if(!is.null(lavdata@weights[[1]])) {
        stop("lavaan ERROR: two.stage + sampling.weights is not supported yet")
    }
    # no fixed.x (yet)
    if(!is.null(lavsamplestats@x.idx) && 
       length(lavsamplestats@x.idx[[1]]) > 0L) {
        stop("lavaan ERROR: two.stage + fixed.x = TRUE is not supported yet")
    }


    # information matrix
    E.inv <- lav_model_information(lavmodel       = lavmodel,
                                   lavsamplestats = lavsamplestats,
                                   lavdata        = lavdata,
                                   lavoptions     = lavoptions,
                                   lavimplied     = lavimplied,
                                   extra          = TRUE,
                                   augmented      = TRUE,
                                   inverted       = TRUE,
                                   use.ginv       = use.ginv)
    Delta <- attr(E.inv, "Delta")
    WLS.V <- attr(E.inv, "WLS.V") # this is 'H' or 'A1' in the literature
    attr(E.inv, "Delta") <- NULL
    attr(E.inv, "WLS.V") <- NULL

    # check if E.inv is ok
    if(inherits(E.inv, "try-error")) {
        return(E.inv)
    }

    # check WLS.V = A1
    if(is.null(WLS.V)) {
        stop("lavaan ERROR: WLS.V/H/A1 is NULL, observed.information = hessian?")
    }

    # Gamma
    Gamma <- vector("list", length = lavsamplestats@ngroups)

    # handle multiple groups
    tDVGVD <- matrix(0, ncol=ncol(E.inv), nrow=nrow(E.inv))
    for(g in 1:lavsamplestats@ngroups) {
        fg  <-  lavsamplestats@nobs[[g]]   /lavsamplestats@ntotal
        #fg1 <- (lavsamplestats@nobs[[g]]-1)/lavsamplestats@ntotal
        fg1 <- fg
        # fg twice for WLS.V, 1/fg1 once for GaMMA
        # if fg==fg1, there would be only one fg, as in Satorra 1999 p.8
        # t(Delta) * WLS.V %*% Gamma %*% WLS.V %*% Delta
        WD <- WLS.V[[g]] %*% Delta[[g]]

        # to compute (incomplete) GAMMA, should we use 
        # structured or unstructured mean/sigma?
        #
        # we use the same setting as to compute 'H' (the h1 information matrix)
        # so that at Omega = H if data is complete
        if(lavoptions$h1.information == "unstructured") {
            MU    <- lavsamplestats@missing.h1[[g]]$mu
            SIGMA <- lavsamplestats@missing.h1[[g]]$sigma
        } else {
            MU    <- lavimplied$mean[[g]]
            SIGMA <- lavimplied$cov[[g]]   
        }

        # compute 'Gamma' (or Omega.beta)
        if(lavoptions$se == "two.stage") {
            # this is Savalei & Bentler (2009)
            if(lavoptions$information == "expected") {
                Info <- lav_mvnorm_missing_information_expected(
                            Y = lavdata@X[[g]], Mp = lavdata@Mp[[g]], 
                            Mu = MU, Sigma = SIGMA)
            } else {
                Info <- lav_mvnorm_missing_information_observed_samplestats(
                            Yp = lavsamplestats@missing[[g]], 
                            Mu = MU, Sigma = SIGMA)
            }
            Gamma[[g]] <- lav_matrix_symmetric_inverse(Info)
        } else { # we assume "robust.two.stage"
                 # NACOV is here incomplete Gamma
                 # Savalei & Falk (2014)
                 #
            Gamma[[g]] <- lav_mvnorm_missing_h1_omega_sw(Y = 
                             lavdata@X[[g]], Mp = lavdata@Mp[[g]], 
                             Yp = lavsamplestats@missing[[g]], 
                             Mu = MU, Sigma = SIGMA,
                             information = lavoptions$information)
        }
 
        # compute
        tDVGVD <- tDVGVD + fg*fg/fg1 * crossprod(WD, Gamma[[g]] %*% WD)
    } # g

    NVarCov <- (E.inv %*% tDVGVD %*% E.inv)

    # to be reused by lavaanTest
    attr(NVarCov, "Delta") <- Delta
    attr(NVarCov, "Gamma") <- Gamma
    #if(lavoptions$h1.information.se == lavoptions$h1.information.test) {
        attr(NVarCov, "E.inv") <- E.inv 
        attr(NVarCov, "WLS.V") <- WLS.V
    #}

    NVarCov
}



lav_model_vcov <- function(lavmodel       = NULL, 
                           lavsamplestats = NULL, 
                           lavoptions     = NULL, 
                           lavdata        = NULL, 
                           lavpartable    = NULL, 
                           lavcache       = NULL,
                           lavimplied     = NULL,
                           use.ginv       = FALSE) {

    likelihood  <- lavoptions$likelihood
    information <- lavoptions$information
    se          <- lavoptions$se
    verbose     <- lavoptions$verbose
    mimic       <- lavoptions$mimic
  
    # special cases
    if(se == "none" || se == "external") return(matrix(0,0,0))

    # some require meanstructure (for now)
    #if(se %in% c("first.order", "robust.sem", "robust.huber.white") && 
    #   !lavoptions$meanstructure) {
    #    stop("se (", se, ") requires meanstructure (for now)")
    #}

    if(se == "standard") {
        NVarCov <- lav_model_information(lavmodel       = lavmodel,
                                         lavsamplestats = lavsamplestats,
                                         lavdata        = lavdata,
                                         lavcache       = lavcache,
                                         lavimplied     = lavimplied,
                                         lavoptions     = lavoptions,
                                         extra          = FALSE,
                                         augmented      = TRUE,
                                         inverted       = TRUE,
                                         use.ginv       = use.ginv)

    } else if(se == "first.order") {
        NVarCov <- 
            lav_model_information_firstorder(lavmodel = lavmodel,
                                             lavsamplestats = lavsamplestats,
                                             lavdata        = lavdata,
                                             lavcache       = lavcache,
                                             lavimplied     = lavimplied,
                                             lavoptions     = lavoptions,
                                             extra          = TRUE,
                                             check.pd       = FALSE,
                                             augmented      = TRUE,
                                             inverted       = TRUE,
                                             use.ginv       = use.ginv)
    
    } else if(se == "robust.sem") {
        NVarCov <-
            lav_model_nvcov_robust_sem(lavmodel       = lavmodel,
                                       lavsamplestats = lavsamplestats,
                                       lavcache       = lavcache,
                                       lavdata        = lavdata,
                                       lavimplied     = lavimplied,
                                       lavoptions     = lavoptions,
                                       use.ginv       = use.ginv)

    } else if(se == "robust.huber.white") {
        NVarCov <-
            lav_model_nvcov_robust_sandwich(lavmodel = lavmodel,
                                            lavsamplestats = lavsamplestats,
                                            lavdata        = lavdata,
                                            lavcache       = lavcache,
                                            lavimplied     = lavimplied,
                                            lavoptions     = lavoptions,
                                            use.ginv       = use.ginv)

    } else if(se %in% c("two.stage", "robust.two.stage")) {
        NVarCov <-
            lav_model_nvcov_two_stage(lavmodel       = lavmodel,
                                      lavsamplestats = lavsamplestats,
                                      lavoptions     = lavoptions,
                                      lavdata        = lavdata,
                                      lavimplied     = lavimplied,
                                      use.ginv       = use.ginv)

    } else if(se == "bootstrap") {
        NVarCov <- try( lav_model_nvcov_bootstrap(lavmodel       = lavmodel,
                                        lavsamplestats = lavsamplestats,
                                        lavoptions     = lavoptions,
                                        lavdata        = lavdata,
                                        lavimplied     = lavimplied,
                                        lavcache       = lavcache,
                                        lavpartable    = lavpartable),
                        silent=TRUE )
    } else {
        warning("lavaan WARNING: unknown se type: ", se)
    }

    if(! inherits(NVarCov, "try-error") ) {

        # denominator!
        if(lavmodel@estimator %in% c("ML","PML","FML") && 
           likelihood == "normal") {
            N <- lavsamplestats@ntotal
        } else {
            N <- lavsamplestats@ntotal - lavsamplestats@ngroups
        }

        VarCov <- 1/N * NVarCov

    } else {
        warning("lavaan WARNING: could not compute standard errors!\n  lavaan NOTE: this may be a symptom that the model is not identified.\n")
        VarCov <- NULL
    }

    VarCov
}

lav_model_vcov_se <- function(lavmodel, lavpartable, VCOV = NULL,
                              BOOT = NULL) {

    # 0. special case
    if(is.null(VCOV)) {
        se <- rep(as.numeric(NA), lavmodel@nx.user)
        se[ lavpartable$free == 0L ] <- 0.0
        return(se)
    }

    # 1. free parameters only
        x.var <- diag(VCOV)
        # check for negative values (what to do: NA or 0.0?)
        x.var[x.var < 0] <- as.numeric(NA)
        x.se <- sqrt( x.var )
        GLIST <- lav_model_x2GLIST(lavmodel = lavmodel, x = x.se, type = "free")

        # se for full parameter table, but with 0.0 entries for def/ceq/cin
        # elements
        se <- lav_model_get_parameters(lavmodel = lavmodel, GLIST = GLIST,
                                       type = "user", extra = FALSE)


    # 2. fixed parameters -> se = 0.0
        se[ which(lavpartable$free == 0L) ] <- 0.0


    # 3. defined parameters: 
        def.idx <- which(lavpartable$op == ":=")
        if(length(def.idx) > 0L) {
            if(!is.null(BOOT)) {
                BOOT.def <- apply(BOOT, 1L, lavmodel@def.function)
                if(length(def.idx) == 1L) {
                    BOOT.def <- as.matrix(BOOT.def)
                } else {
                    BOOT.def <- t(BOOT.def)
                }
                def.cov <- cov(BOOT.def )
            } else {
                # regular delta method
                x <- lav_model_get_parameters(lavmodel = lavmodel, type = "free")
                 JAC <- try(lav_func_jacobian_complex(func = lavmodel@def.function, x = x),
                           silent=TRUE)
                if(inherits(JAC, "try-error")) { # eg. pnorm()
                    JAC <- lav_func_jacobian_simple(func = lavmodel@def.function, x = x)
                }
                def.cov <- JAC %*% VCOV %*% t(JAC)
            }
            # check for negative se's
            diag.def.cov <- diag(def.cov)
            diag.def.cov[ diag.def.cov < 0 ] <- as.numeric(NA)
            se[def.idx] <- sqrt(diag.def.cov)
        }

    se
}
