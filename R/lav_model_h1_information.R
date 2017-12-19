# the information matrix of the unrestricted (H1) model
# taking into account:
#   - the estimator (ML or (D)WLS/ULS)
#   - missing or not
#   - fixed.x = TRUE or FALSE
#   - conditional.x = TRUE or FALSE
#   - h1.information is "structured" or "unstructured"
#
# Note: this replaces the (old) lav_model_wls_v() function
#
# - YR 22 Okt 2017
# - YR 03 Dec 2017: add lavh1, implied is either lavimplied or lavh1
#                   add support for clustered data

lav_model_h1_information <- function(lavobject      = NULL,
                                     lavmodel       = NULL,
                                     lavsamplestats = NULL,
                                     lavdata        = NULL,
                                     lavimplied     = NULL,
                                     lavh1          = NULL,
                                     lavcache       = NULL,
                                     lavoptions     = NULL) {

    if(!is.null(lavobject) && inherits(lavobject, "lavaan")) {
        lavmodel       <- lavobject@Model
        lavsamplestats <- lavobject@SampleStats
        lavdata        <- lavobject@Data
        lavimplied     <- lavobject@implied
        lavh1          <- lavobject@h1
        lavcache       <- lavobject@Cache
        lavoptions     <- lavobject@Options
    }


    estimator   <- lavmodel@estimator
    information <- lavoptions$information

    # compute information matrix
    if(information == "observed") {
        I1 <- lav_model_h1_information_observed(lavmodel = lavmodel,
            lavsamplestats = lavsamplestats, lavdata = lavdata,
            lavimplied = lavimplied, lavh1 = lavh1,
            lavcache = lavcache, lavoptions = lavoptions)
    } else if(information == "expected") {
        I1 <- lav_model_h1_information_expected(lavmodel = lavmodel,
            lavsamplestats = lavsamplestats, lavdata = lavdata,
            lavimplied = lavimplied, lavh1 = lavh1,
            lavcache = lavcache, lavoptions = lavoptions)
    } else if(information == "first.order") {
        I1 <- lav_model_h1_information_firstorder(lavmodel = lavmodel,
            lavsamplestats = lavsamplestats, lavdata = lavdata,
            lavimplied = lavimplied, lavh1 = lavh1,
            lavcache = lavcache, lavoptions = lavoptions)
    }

    # I1 information, as a list per group
    I1
}

# fisher/expected information of H1
lav_model_h1_information_expected <- function(lavobject      = NULL,
                                              lavmodel       = NULL,
                                              lavsamplestats = NULL,
                                              lavdata        = NULL,
                                              lavoptions     = NULL,
                                              lavimplied     = NULL,
                                              lavh1          = NULL,
                                              lavcache       = NULL) {

    if(!is.null(lavobject) && inherits(lavobject, "lavaan")) {
        lavmodel       <- lavobject@Model
        lavsamplestats <- lavobject@SampleStats
        lavdata        <- lavobject@Data
        lavimplied     <- lavobject@implied
        lavh1          <- lavobject@h1
        lavcache       <- lavobject@Cache
        lavoptions     <- lavobject@Options
    }


    estimator <- lavmodel@estimator

    # structured of unstructured? (since 0.5-23)
    if(!is.null(lavoptions) &&
       !is.null(lavoptions$h1.information) &&
       lavoptions$h1.information == "unstructured") {
        structured <- FALSE
    } else {
        structured <- TRUE
    }

    # 1. WLS.V (=A1) for GLS/WLS
    if(lavmodel@estimator == "GLS"  || lavmodel@estimator == "WLS") {
        A1 <- lavsamplestats@WLS.V
    }

    # 2. DWLS/ULS diagonal @WLS.VD slot
    else if(lavmodel@estimator == "DWLS"  || lavmodel@estimator == "ULS") {
        # diagonal only!!
        A1 <- lavsamplestats@WLS.VD
    }

    # 3. ML
    else if(lavmodel@estimator == "ML" ||
            lavmodel@estimator == "NTRLS") {
        A1 <- vector("list", length=lavsamplestats@ngroups)

        # structured? compute model-implied statistics
        if(structured && is.null(lavimplied)) {
            lavimplied <- lav_model_implied(lavmodel)
        }

        for(g in 1:lavsamplestats@ngroups) {

            if(lavsamplestats@missing.flag) {
                # mvnorm
                # FIXME: allow for meanstructure = FALSE
                # FIXME: allow for conditional.x = TRUE
                # FIXME: allow for wt
                # FIXME: allow for x.idx?
                if(lavmodel@meanstructure && structured) {
                    MEAN <- lavimplied$mean[[g]]
                } else {
                    MEAN <- lavsamplestats@missing.h1[[g]]$mu
                }

                if(structured) {
                    A1[[g]] <- 
                      lav_mvnorm_missing_information_expected(
                          Y = lavdata@X[[g]],
                          Mp = lavdata@Mp[[g]],
                          wt = lavdata@weights[[g]],
                          Mu = MEAN,
                          # meanstructure = lavmodel@meanstructure,
                          Sigma = lavimplied$cov[[g]])
                } else {
                    A1[[g]] <-
                      lav_mvnorm_missing_information_expected(
                          Y = lavdata@X[[g]],
                          Mp = lavdata@Mp[[g]],
                          Mu = MEAN,
                          # meanstructure = lavmodel@meanstructure,
                          Sigma = lavsamplestats@missing.h1[[g]]$sigma)
                }
            } else {
                if(lavmodel@conditional.x) {
                    # mvreg
                    if(lavmodel@meanstructure && structured) {
                        RES.INT    <- lavimplied$res.int[[g]]
                        RES.SLOPES <- lavimplied$res.slopes[[g]]
                    } else {
                        RES.INT    <- lavsamplestats@res.int[[g]]
                        RES.SLOPES <- lavsamplestats@res.slopes[[g]]
                    }

                    if(structured) {
                        A1[[g]] <- lav_mvreg_information_expected(
                            sample.mean.x     = lavsamplestats@mean.x[[g]],
                            sample.cov.x      = lavsamplestats@cov.x[[g]],
                            sample.nobs       = lavsamplestats@nobs[[g]],
                            res.int           = RES.INT,
                            res.slopes        = RES.SLOPES,
                            #wt               = lavdata@weights[[g]],
                            #meanstructure    = lavmodel@meanstructure,
                            res.cov           = lavimplied$res.cov[[g]])
                    } else {
                        A1[[g]] <- lav_mvreg_information_expected(
                            sample.mean.x     = lavsamplestats@mean.x[[g]],
                            sample.cov.x      = lavsamplestats@cov.x[[g]],
                            sample.nobs       = lavsamplestats@nobs[[g]],
                            res.int           = lavsamplestats@res.int[[g]],
                            res.slopes        = lavsamplestats@res.slopes[[g]],
                            #wt               = lavdata@weights[[g]],
                            #meanstructure    = lavmodel@meanstructure,
                            res.cov           = lavsamplestats@res.cov[[g]])
                    }

                } else {
                    # conditional.x = FALSE
                    # mvnorm
                    if(lavmodel@meanstructure && structured) {
                        MEAN <- lavimplied$mean[[g]]
                    } else {
                        MEAN <- lavsamplestats@mean[[g]]
                    }

                    if(structured) {
                        A1[[g]] <- lav_mvnorm_information_expected(
                                  Sigma         = lavimplied$cov[[g]],
                                  #wt = lavdata@weights[[g]], # not needed
                                  x.idx         = lavsamplestats@x.idx[[g]],
                                  meanstructure = lavmodel@meanstructure)
                    } else {
                        A1[[g]] <- lav_mvnorm_h1_information_expected(
                                  sample.cov.inv = lavsamplestats@icov[[g]],
                                  #wt = lavdata@weights[[g]], not needed
                                  x.idx          = lavsamplestats@x.idx[[g]],
                                  meanstructure  = lavmodel@meanstructure)
                    }
                } # conditional.x
            } # missing

            # stochastic group weight
            if(lavmodel@group.w.free) {
                # unweight!!
                a <- exp(lavimplied$group.w[[g]]) / lavsamplestats@nobs[[g]]
                A1[[g]] <- lav_matrix_bdiag( matrix(a, 1L, 1L), A1[[g]])
            }

        } # g
    } # ML

    A1
}

lav_model_h1_information_observed <- function(lavobject      = NULL,
                                              lavmodel       = NULL,
                                              lavsamplestats = NULL,
                                              lavdata        = NULL,
                                              lavimplied     = NULL,
                                              lavh1          = NULL,
                                              lavcache       = NULL,
                                              lavoptions     = NULL) {

    if(!is.null(lavobject) && inherits(lavobject, "lavaan")) {
        lavmodel       <- lavobject@Model
        lavsamplestats <- lavobject@SampleStats
        lavdata        <- lavobject@Data
        lavimplied     <- lavobject@implied
        lavh1          <- lavobject@h1
        lavcache       <- lavobject@Cache
        lavoptions     <- lavobject@Options
    }

    estimator <- lavmodel@estimator

    # structured?
    if(!is.null(lavoptions) &&
       !is.null(lavoptions$h1.information) &&
       lavoptions$h1.information == "unstructured") {
        structured <- FALSE
    } else {
        structured <- TRUE
    }
 
    # 1. WLS.V (=A1) for GLS/WLS
    if(lavmodel@estimator == "GLS"  || lavmodel@estimator == "WLS") {
        A1 <- lavsamplestats@WLS.V
    }

    # 2. DWLS/ULS diagonal @WLS.VD slot
    else if(lavmodel@estimator == "DWLS"  || lavmodel@estimator == "ULS") {
        # diagonal only!!
        A1 <- lavsamplestats@WLS.VD
    }

    # 3. ML
    else if(lavmodel@estimator == "ML") {
        A1 <- vector("list", length=lavsamplestats@ngroups)
  
        # structured? compute model-implied statistics
        if(structured && is.null(lavimplied)) {
            lavimplied <- lav_model_implied(lavmodel)
        }

        for(g in 1:lavsamplestats@ngroups) {

            if(lavsamplestats@missing.flag) {
                # mvnorm
                # FIXME: allow for meanstructure = FALSE
                # FIXME: allow for conditional.x = TRUE
                # FIXME: allow for wt
                # FIXME: allow for x.idx?
                if(lavmodel@meanstructure && structured) {
                    MEAN <- lavimplied$mean[[g]]
                } else {
                    MEAN <- lavsamplestats@missing.h1[[g]]$mu
                }

                if(structured) {
                    A1[[g]] <- 
                      lav_mvnorm_missing_information_observed_samplestats(
                          Yp = lavsamplestats@missing[[g]],
                          #wt = lavdata@weights[[g]], ?
                          Mu = MEAN,
                          # meanstructure = lavmodel@meanstructure,
                          Sigma = lavimplied$cov[[g]])
                } else {
                    A1[[g]] <-
                      lav_mvnorm_missing_information_observed_samplestats(
                          Yp = lavsamplestats@missing[[g]],
                          #wt = lavdata@weights[[g]], ?
                          Mu = MEAN,
                          # meanstructure = lavmodel@meanstructure,
                          Sigma = lavsamplestats@missing.h1[[g]]$sigma)
                }
            } else {
                if(lavmodel@conditional.x) {
                    # mvreg
                    if(lavmodel@meanstructure && structured) {
                        RES.INT    <- lavimplied$res.int[[g]]
                        RES.SLOPES <- lavimplied$res.slopes[[g]]
                    } else {
                        RES.INT    <- lavsamplestats@res.int[[g]]
                        RES.SLOPES <- lavsamplestats@res.slopes[[g]]
                    }

                    if(structured) {
                        A1[[g]] <- lav_mvreg_information_observed_samplestats(
                            sample.res.int    = lavsamplestats@res.int[[g]],
                            sample.res.slopes = lavsamplestats@res.slopes[[g]],
                            sample.res.cov    = lavsamplestats@res.cov[[g]],
                            sample.mean.x     = lavsamplestats@mean.x[[g]],
                            sample.cov.x      = lavsamplestats@cov.x[[g]],
                            res.int           = RES.INT,
                            res.slopes        = RES.SLOPES,
                            #wt               = lavdata@weights[[g]],
                            #meanstructure    = lavmodel@meanstructure,
                            res.cov           = lavimplied$res.cov[[g]])
                    } else {
                        A1[[g]] <- lav_mvreg_information_observed_samplestats(
                            sample.res.int    = lavsamplestats@res.int[[g]],
                            sample.res.slopes = lavsamplestats@res.slopes[[g]],
                            sample.res.cov    = lavsamplestats@res.cov[[g]],
                            sample.mean.x     = lavsamplestats@mean.x[[g]],
                            sample.cov.x      = lavsamplestats@cov.x[[g]],
                            res.int           = lavsamplestats@res.int[[g]],
                            res.slopes        = lavsamplestats@res.slopes[[g]],
                            #wt               = lavdata@weights[[g]],
                            #meanstructure    = lavmodel@meanstructure,
                            res.cov           = lavsamplestats@res.cov[[g]])
                    }

                } else {
                    # conditional.x = FALSE
                    # mvnorm
                    if(lavmodel@meanstructure && structured) {
                        MEAN <- lavimplied$mean[[g]]
                    } else {
                        MEAN <- lavsamplestats@mean[[g]]
                    }

                    if(structured) {
                        A1[[g]] <- lav_mvnorm_information_observed_samplestats(
                                  sample.mean   = lavsamplestats@mean[[g]],
                                  sample.cov    = lavsamplestats@cov[[g]],
                                  Mu            = MEAN,
                                  Sigma         = lavimplied$cov[[g]],
                                  #wt = lavdata@weights[[g]], # not needed
                                  x.idx         = lavsamplestats@x.idx[[g]],
                                  meanstructure = lavmodel@meanstructure)
                    } else {
                        A1[[g]] <- lav_mvnorm_h1_information_observed_samplestats(
                                  sample.mean    = lavsamplestats@mean[[g]],
                                  sample.cov     = lavsamplestats@cov[[g]],
                                  sample.cov.inv = lavsamplestats@icov[[g]],
                                  #wt = lavdata@weights[[g]], not needed
                                  x.idx          = lavsamplestats@x.idx[[g]],
                                  meanstructure  = lavmodel@meanstructure)
                    }
                } # conditional.x
            } # missing

            # stochastic group weight
            if(lavmodel@group.w.free) {
                # unweight!!
                a <- exp(lavimplied$group.w[[g]]) / lavsamplestats@nobs[[g]]
                A1[[g]] <- lav_matrix_bdiag( matrix(a,1,1), A1[[g]])
            }

        } # g
    } # ML

    A1
}

# outer product of the case-wise scores (gradients)
lav_model_h1_information_firstorder <- function(lavobject      = NULL,
                                                lavmodel       = NULL,
                                                lavsamplestats = NULL,
                                                lavdata        = NULL,
                                                lavimplied     = NULL,
                                                lavh1          = NULL,
                                                lavcache       = NULL,
                                                lavoptions     = NULL) {

    if(!is.null(lavobject) && inherits(lavobject, "lavaan")) {
        lavmodel       <- lavobject@Model
        lavsamplestats <- lavobject@SampleStats
        lavdata        <- lavobject@Data
        lavimplied     <- lavobject@implied
        lavh1          <- lavobject@h1
        lavcache       <- lavobject@Cache
        lavoptions     <- lavobject@Options
    }

    estimator <- lavmodel@estimator
    if(!estimator %in% c("ML", "PML")) {
        stop("lavaan ERROR: information = \"first.order\" not available for estimator ", sQuote(estimator))
    }

    if(!is.null(lavoptions) &&
       !is.null(lavoptions$h1.information) &&
       lavoptions$h1.information == "unstructured") {
        structured <- FALSE
    } else {
        structured <- TRUE
    }

    # structured? compute model-implied statistics
    if(estimator == "PML" || structured) {
        if(is.null(lavimplied)) {
            lavimplied <- lav_model_implied(lavmodel)
        }
    }

    # structured? lavimplied vs lavh1
    if(structured) {
        implied <- lavimplied
    } else {
        implied <- lavh1$implied
    }

    B1 <- vector("list", length=lavsamplestats@ngroups)
    for(g in 1:lavdata@ngroups) {
        if(estimator == "PML") {
            # slow approach: compute outer product of case-wise scores

            if(lavmodel@conditional.x) {
                SIGMA <- implied$res.cov[[g]]
                TH    <- implied$res.th[[g]]
                PI    <- implied$res.slopes[[g]]
            } else {
                SIGMA <- implied$cov[[g]]
                TH    <- implied$th[[g]]
                PI    <- implied$slopes[[g]]
            }
            SC <- pml_deriv1(Sigma.hat  = SIGMA,
                             TH         = TH,
                             th.idx     = lavmodel@th.idx[[g]],
                             num.idx    = lavmodel@num.idx[[g]],
                             X          = lavdata@X[[g]],
                             eXo        = lavdata@eXo[[g]],
                             PI         = PI,
                             lavcache   = lavcache[[g]],
                             missing    = lavdata@missing,
                             scores     = TRUE,
                             negative   = FALSE)
            # information H1
            B1[[g]] <- crossprod(SC)

        } else if(estimator == "ML" && lavdata@nlevels > 1L) {

            MU.W    <- implied$mean[[ (g-1)*lavdata@nlevels + 1L ]]
            MU.B    <- implied$mean[[ (g-1)*lavdata@nlevels + 2L ]]
            SIGMA.W <- implied$cov[[  (g-1)*lavdata@nlevels + 1L ]]
            SIGMA.B <- implied$cov[[  (g-1)*lavdata@nlevels + 2L ]]

            # clustered data
            B1[[g]] <- lav_mvnorm_cluster_information_firstorder(
                           Y1           = lavdata@X[[g]],
                           YLp          = lavsamplestats@YLp[[g]],
                           Lp           = lavdata@Lp[[g]],
                           Mu.W         = MU.W,
                           Sigma.W      = SIGMA.W,
                           Mu.B         = MU.B,
                           Sigma.B      = SIGMA.B,
                           divide.by.two = TRUE)

        } else if(estimator == "ML") {
            if(lavsamplestats@missing.flag) {
                # mvnorm
                # FIXME: allow for meanstructure = FALSE
                # FIXME: allow for conditional.x = TRUE
                # FIXME: allow for wt
                # FIXME: allow for x.idx?
                if(lavmodel@meanstructure && structured) {
                    MEAN <- lavimplied$mean[[g]]
                } else {
                    MEAN <- lavsamplestats@missing.h1[[g]]$mu
                }
 
                B1[[g]] <- lav_mvnorm_missing_information_firstorder(
                               Y = lavdata@X[[g]],
                              Mp = lavdata@Mp[[g]], wt = lavdata@weights[[g]],
                              Mu = MEAN,
                              # meanstructure = lavmodel@meanstructure,
                              Sigma = implied$cov[[g]])

            } else {
                if(lavmodel@conditional.x) {
                    # mvreg
                    if(lavmodel@meanstructure && structured) {
                        RES.INT    <- lavimplied$res.int[[g]]
                        RES.SLOPES <- lavimplied$res.slopes[[g]]
                    } else {
                        RES.INT    <- lavsamplestats@res.int[[g]]
                        RES.SLOPES <- lavsamplestats@res.slopes[[g]]
                    }

                    B1[[g]] <- lav_mvreg_information_firstorder(
                                  Y              = lavdata@X[[g]],
                                  eXo            = lavdata@eXo[[g]],
                                  res.int        = RES.INT,
                                  res.slopes     = RES.SLOPES,
                                  #wt            = lavdata@weights[[g]],
                                  #meanstructure = lavmodel@meanstructure,
                                  res.cov        = implied$res.cov[[g]])
                } else {
                    # conditional.x = FALSE
                    # mvnorm
                    if(lavmodel@meanstructure && structured) {
                        MEAN <- lavimplied$mean[[g]]
                    } else {
                        # NOTE: the information matrix will be the same (minus
                        # the meanstructure block), but once INVERTED, the
                        # standard errors will be (slightly) smaller!!!
                        # This is only visibile when estimator = "MLF"
                        # (or information = "first.order")
                        MEAN <- lavsamplestats@mean[[g]] # saturated
                    }

                    if(structured) {
                        B1[[g]] <- lav_mvnorm_information_firstorder(
                                  Y = lavdata@X[[g]],
                                  Mu = MEAN, Sigma = lavimplied$cov[[g]],
                                  wt = lavdata@weights[[g]],
                                  x.idx = lavsamplestats@x.idx[[g]],
                                  meanstructure = lavmodel@meanstructure)
                    } else {
                        B1[[g]] <- lav_mvnorm_h1_information_firstorder(
                                  Y = lavdata@X[[g]],
                                  sample.cov.inv = lavsamplestats@icov[[g]],
                                  Gamma = lavsamplestats@NACOV[[g]],
                                  wt = lavdata@weights[[g]],
                                  x.idx = lavsamplestats@x.idx[[g]],
                                  meanstructure = lavmodel@meanstructure)
                    }
                } # mvnorm
            } # missing
        } # ML

        # stochastic group weight
        if(lavmodel@group.w.free) {
            # unweight!!
            a <- exp(lavimplied$group.w[[g]]) / lavsamplestats@nobs[[g]]
            B1[[g]] <- lav_matrix_bdiag( matrix(a,1,1), B1[[g]])
        }

    } # g

    B1
}

