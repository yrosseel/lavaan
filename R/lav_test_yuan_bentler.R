lav_test_yuan_bentler <- function(lavobject      = NULL,
                                  lavsamplestats = NULL,
                                  lavmodel       = NULL,
                                  lavimplied     = NULL,
                                  lavh1          = NULL,
                                  lavoptions     = NULL,
                                  lavdata        = NULL,
                                  TEST.unscaled  = NULL,
                                  E.inv          = NULL,
                                  B0.group       = NULL,
                                  test           = "yuan.bentler",
                                  mimic          = "lavaan",
                                  #method         = "default",
                                  return.ugamma  = FALSE) {

    TEST <- list()

    if(!is.null(lavobject)) {
        lavsamplestats <- lavobject@SampleStats
        lavmodel       <- lavobject@Model
        lavoptions     <- lavobject@Options
        lavpartable    <- lavobject@ParTable
        lavimplied     <- lavobject@implied
        lavh1          <- lavobject@h1
        lavdata        <- lavobject@Data
        TEST$standard  <- lavobject@test[[1]]
    } else {
        TEST$standard  <- TEST.unscaled
    }

    # check test
    if(!all(test %in% c("yuan.bentler",
                        "yuan.bentler.mplus"))) {
        warning("lavaan WARNING: test must be one of `yuan.bentler', or `yuan.bentler.mplus'; will use `yuan.bentler' only")
        test <- "yuan.bentler"
    }

    # information
    information <- lavoptions$information
    # x.idx
    if(lavoptions$conditional.x) {
        x.idx <- NULL
    } else {
        x.idx <- lavsamplestats@x.idx
    }
    # ndat
    ndat <- numeric(lavsamplestats@ngroups)


    if(is.null(E.inv)) {
        E.inv <- try(lav_model_information(lavmodel       = lavmodel,
                                           lavsamplestats = lavsamplestats,
                                           lavdata        = lavdata,
                                           lavimplied     = lavimplied,
                                           lavoptions     = lavoptions,
                                           extra          = FALSE,
                                           augmented      = TRUE,
                                           inverted       = TRUE),
                      silent = TRUE)
        if(inherits(E.inv, "try-error")) {
            TEST <- list(test = test, stat = as.numeric(NA),
                stat.group = rep(as.numeric(NA), lavsamplestats@ngroups),
                df = TEST$standard$df, refdistr = TEST$standard$refdistr,
                pvalue = as.numeric(NA), scaling.factor = as.numeric(NA))
            warning("lavaan WARNING: could not invert information matrix\n")
            return(TEST)
        }
    }

    # mean and variance adjusted?
    Satterthwaite <- FALSE # for now
    #if(any(test %in% c("mean.var.adjusted", "scaled.shifted"))) {
    #    Satterthwaite <- TRUE
    #}

    # FIXME: should we not always use 'unstructured' here?
    # if the model is, say, the independence model, the
    # 'structured' information (A1) will be so far away from B1
    # that we will end with 'NA'
    h1.options <- lavoptions
    if(test == "yuan.bentler.mplus") {
        # always 'unstructured' H1 information
        h1.options$h1.information <- "unstructured"
    }

    # A1 is usually expected or observed
    A1.group <- lav_model_h1_information(lavmodel       = lavmodel,
                                         lavsamplestats = lavsamplestats,
                                         lavdata        = lavdata,
                                         lavimplied     = lavimplied,
                                         lavh1          = lavh1,
                                         lavoptions     = h1.options)
    # B1 is always first.order
    B1.group <- lav_model_h1_information_firstorder(lavmodel = lavmodel,
                                         lavsamplestats = lavsamplestats,
                                         lavdata        = lavdata,
                                         lavimplied     = lavimplied,
                                         lavh1          = lavh1,
                                         lavoptions     = h1.options)

    if(test == "yuan.bentler.mplus") {
        if(is.null(B0.group)) {
            B0 <- lav_model_information_firstorder(lavmodel = lavmodel,
                                             lavsamplestats = lavsamplestats,
                                             lavdata        = lavdata,
                                             lavh1          = lavh1,
                                             lavoptions     = lavoptions,
                                             extra          = TRUE,
                                             check.pd       = FALSE,
                                             augmented      = FALSE,
                                             inverted       = FALSE)
           B0.group <- attr(B0, "B0.group")
        }
        trace.UGamma <-
            lav_test_yuan_bentler_mplus_trace(lavsamplestats = lavsamplestats,
                                              A1.group       = A1.group,
                                              B1.group       = B1.group,
                                              B0.group       = B0.group,
                                              E.inv          = E.inv,
                                              x.idx          = x.idx,
                                              meanstructure  = lavmodel@meanstructure)
    } else if(test == "yuan.bentler") {

        Delta <- computeDelta(lavmodel = lavmodel)
        trace.UGamma <- lav_test_yuan_bentler_trace(
            lavsamplestats = lavsamplestats,
            meanstructure  = lavmodel@meanstructure,
            A1.group       = A1.group,
            B1.group       = B1.group,
            Delta          = Delta,
            E.inv          = E.inv,
            x.idx          = x.idx,
            Satterthwaite  = TRUE) # for now
    }

    # unscaled test
    df <- TEST$standard$df
    chisq.group <- TEST$standard$stat.group

    scaling.factor       <- sum(trace.UGamma) / df
    if(scaling.factor < 0) scaling.factor <- as.numeric(NA)
    chisq.scaled         <- sum(chisq.group / scaling.factor)
    pvalue.scaled        <- 1 - pchisq(chisq.scaled, df)

    ndat <- sum(attr(trace.UGamma, "h1.ndat"))
    npar <- lavmodel@nx.free

    scaling.factor.h1    <- sum( attr(trace.UGamma, "h1") ) / ndat
    scaling.factor.h0    <- sum( attr(trace.UGamma, "h0") ) / npar
    trace.UGamma2        <- attr(trace.UGamma, "trace.UGamma2")
    attributes(trace.UGamma) <- NULL

    if("yuan.bentler" %in% test) {
        TEST$yuan.bentler <-
            list(test              = test,
                 stat              = chisq.scaled,
                 stat.group        = (chisq.group / scaling.factor),
                 df                = df,
                 pvalue            = pvalue.scaled,
                 scaling.factor    = scaling.factor,
                 scaling.factor.h1 = scaling.factor.h1,
                 scaling.factor.h0 = scaling.factor.h0,
                 trace.UGamma      = trace.UGamma,
                 trace.UGamma2     = trace.UGamma2)

    } else if("yuan.bentler.mplus" %in% test) {
        TEST$yuan.bentler.mplus <-
            list(test              = test,
                 stat              = chisq.scaled,
                 stat.group        = (chisq.group / scaling.factor),
                 df                = df,
                 pvalue            = pvalue.scaled,
                 scaling.factor    = scaling.factor,
                 scaling.factor.h1 = scaling.factor.h1,
                 scaling.factor.h0 = scaling.factor.h0,
                 trace.UGamma      = trace.UGamma,
                 trace.UGamma2     = as.numeric(NA))
    }

    TEST
}


lav_test_yuan_bentler_trace <- function(lavsamplestats =lavsamplestats,
                                        meanstructure  = TRUE,
                                        A1.group       = NULL,
                                        B1.group       = NULL,
                                        Delta          = NULL,
                                        E.inv          = NULL,
                                        x.idx          = list(integer(0)),
                                        Satterthwaite  = FALSE) {

    # we always assume a meanstructure (nope, not any longer, since 0.6)
    #meanstructure <- TRUE

    trace.UGamma  <- numeric( lavsamplestats@ngroups )
    trace.UGamma2 <- numeric( lavsamplestats@ngroups )
    trace.h1      <- numeric( lavsamplestats@ngroups )
    trace.h0      <- numeric( lavsamplestats@ngroups )
    h1.ndat       <- numeric( lavsamplestats@ngroups )

    for(g in 1:lavsamplestats@ngroups) {

        # group weight
        fg <- lavsamplestats@nobs[[g]]/lavsamplestats@ntotal

        A1 <- A1.group[[g]]
        B1 <- B1.group[[g]]
        DELTA <- Delta[[g]]

        # mask independent 'fixed-x' variables
        zero.idx <- which(diag(A1) == 0)
        if(length(zero.idx) > 0L) {
            A1.inv <- matrix(0, nrow(A1), ncol(A1))
            a1 <- A1[-zero.idx, -zero.idx]
            a1.inv <- solve(a1)
            A1.inv[-zero.idx, -zero.idx] <- a1.inv
        } else {
            A1.inv <- solve(A1)
        }
        h1.ndat[g] <- ncol(A1) - length(zero.idx)

        D.Einv.tD <- DELTA %*% tcrossprod(E.inv, DELTA)

        trace.h1[g] <- sum( B1 * t( A1.inv ) )
        # fg cancels out: trace.h1[g] <- sum( fg*B1 * t( 1/fg*A1.inv ) )
        trace.h0[g] <- fg * sum( B1 * D.Einv.tD )
        trace.UGamma[g] <- trace.h1[g] - trace.h0[g]

        if(Satterthwaite) {
            A1invB1 <- A1.inv %*% B1
            UG <- A1invB1 - fg * (A1invB1 %*% D.Einv.tD %*% A1)
            trace.UGamma2[g] <- sum(UG * t(UG))
        }
    }

    # traces
    trace.UGamma <- sum(trace.UGamma)
    attr(trace.UGamma, "h1") <- trace.h1
    attr(trace.UGamma, "h0") <- trace.h0
    attr(trace.UGamma, "h1.ndat") <- h1.ndat

    if(Satterthwaite) {
        attr(trace.UGamma, "trace.UGamma2") <- sum(trace.UGamma2)
    }

    trace.UGamma
}

lav_test_yuan_bentler_mplus_trace <- function(lavsamplestats=NULL,
                                              A1.group = NULL,
                                              B1.group = NULL,
                                              B0.group=NULL,
                                              E.inv=NULL,
                                              x.idx=list(integer(0)),
                                              meanstructure = TRUE) {
    # typical for Mplus:
    # - do NOT use the YB formula, but use an approximation
    #   relying  on A0 ~= Delta' A1 Delta and the same for B0
    #
    # NOTE: if A0 is based on the hessian, then A0 only approximates
    #       Delta' A1 Delta
    #
    # - always use h1.information = "unstructured"!!!

    ngroups <- lavsamplestats@ngroups

    trace.UGamma <- numeric( lavsamplestats@ngroups )
    trace.h1     <- numeric( lavsamplestats@ngroups )
    trace.h0     <- numeric( lavsamplestats@ngroups )
    h1.ndat      <- numeric( lavsamplestats@ngroups )

    for(g in 1:lavsamplestats@ngroups) {

        # group weight
        fg <- lavsamplestats@nobs[[g]]/lavsamplestats@ntotal

        A1 <- A1.group[[g]]
        B1 <- B1.group[[g]]

        # mask independent 'fixed-x' variables
        zero.idx <- which(diag(A1) == 0)
        if(length(zero.idx) > 0L) {
            A1.inv <- matrix(0, nrow(A1), ncol(A1))
            a1 <- A1[-zero.idx, -zero.idx]
            a1.inv <- solve(a1)
            A1.inv[-zero.idx, -zero.idx] <- a1.inv
        } else {
            A1.inv <- solve(A1)
        }
        h1.ndat[g] <- ncol(A1) - length(zero.idx)

        # if data is complete, why not just A1 %*% Gamma?
        trace.h1[g]     <- sum( B1 * t( A1.inv ) )
        trace.h0[g]     <- fg * sum( B0.group[[g]] * t(E.inv) )
        trace.UGamma[g] <- (trace.h1[g] - trace.h0[g])
    }

    # we take the sum here
    trace.UGamma <- sum(trace.UGamma)

    attr(trace.UGamma, "h1") <- trace.h1
    attr(trace.UGamma, "h0") <- trace.h0
    attr(trace.UGamma, "h1.ndat") <- h1.ndat

    trace.UGamma
}


