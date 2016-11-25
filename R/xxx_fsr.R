# factor score regression

# three methods:
#  - naive (regression or Bartlett)
#  - Skrondal & Laake (2001) (regression models only)
#  - Croon (2002) (general + robust SE)

fsr <- function(model = NULL, data = NULL, cmd = "sem",
                fsr.method = "Croon", fs.method  = "Bartlett", ...) {
   
    # we need full data
    if(is.null(data)) {
        stop("lavaan ERROR: full data is required for factor score regression")
    }

    # check arguments
    fsr.method <- tolower(fsr.method)
    if(fsr.method == "naive") {
        # nothing to do
    } else if(fsr.method %in% c("skrondal", "laake", "skrondallaake",
                            "skrondal.laake", "skrondal-laake")) {
        fsr.method <- "skrondal.laake"
    } else if(fsr.method == "croon") {
        # nothing to do
    }

    fs.method <- tolower(fs.method)

    # dot dot dot
    dotdotdot <- list(...)

    # check dotdotdot
    if(!is.null(dotdotdot$meanstructure)) {
        dotdotdot$meanstructure <- NULL
    }

    if(!is.null(dotdotdot$do.fit)) {
        dotdotdot$do.fit <- NULL
    }
  

    # check for arguments that we do not want (eg sample.cov)?
    # TODO

    # process model, no fitting
    FIT <- do.call(cmd, 
                   args =  c(list(model  = model, 
                                  data   = data,
                                  meanstructure = TRUE,
                                  do.fit = FALSE), dotdotdot) )
    ngroups <- lavInspect(FIT, "ngroups")

    # any `regular' latent variables?
    lv.names <- unique(unlist(FIT@pta$vnames$lv.regular))
    nfac     <- length(lv.names)
    if(length(lv.names) == 0L) {
        stop("lavaan ERROR: model does not contain any latent variables")
    }

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
            stop("lavaan ERROR: mediational relationship are not allowed for the Skrondal.Laake method; use ", sQuote("Croon"), " instead.")
        }
    }


    # STEP 1:
    # compute factor scores, per latent variable
    FS.SCORES  <- vector("list", length = ngroups)
    FSR.SCORES <- vector("list", length = ngroups)
    FS.COV     <- vector("list", length = ngroups)
    FSR.COV    <- vector("list", length = ngroups)
    LVINFO     <- vector("list", length = ngroups)
    if(ngroups > 1L) {
        names(FS.SCORES) <- names(LVINFO) <- names(FS.COV) <- 
        names(FSR.SCORES) <- names(FSR.COV) <-
            lavInspect(FIT, "group.label")
    }

    for(g in 1:ngroups) {
        FS.SCORES[[g]]  <- vector("list", length = nfac)
        FSR.SCORES[[g]] <- vector("list", length = nfac)
        names(FS.SCORES[[g]]) <- names(FSR.SCORES[[g]]) <- lv.names

        LVINFO[[g]] <- vector("list", length = nfac)
        names(LVINFO[[g]]) <- lv.names
    }

    # adjust options
    dotdotdot2 <- dotdotdot
    dotdotdot2$se <- "none"
    dotdotdot2$test <- "none"
    dotdotdot2$debug <- FALSE
    dotdotdot2$verbose <- FALSE
 

    # we assume the same number/names of lv's per group!!!
    for(f in 1:nfac) {

        FAC <- lv.names[f]
        IND <- PT$rhs[ PT$op == "=~" & PT$lhs == FAC ]

        # check number of indicators
        if(length(IND) < 3L) {
            stop("lavaan ERROR: fsr currently needs at least 3 indicators per factor")
        }

        # the latent variable definitions
        op.idx <- which(PT$op == "=~" & PT$lhs == FAC)

        # the residuals + factor variance
        var.idx <- which(PT$op == "~~" & PT$lhs %in% c(IND,FAC) &
                                         PT$rhs %in% c(IND,FAC) &
                                         PT$lhs == PT$rhs)

        # any residual covariances among the indicators
        cov.idx <- which(PT$op == "~~" & PT$lhs %in% IND &
                                         PT$rhs %in% IND &
                                         PT$lhs != PT$rhs)

        # means/intercepts
        # ov.int <- which(PT$op == "~1" & PT$lhs %in% IND)
        # lv.int <- which(PT$op == "~1" & PT$lhs == FAC)


        # any regression where lhs is an indicator
        reg.idx <- which(PT$op == "~" & PT$lhs %in% IND & 
                                        !PT$rhs %in% lv.names)

        # and their variances...
        reg.names <- PT$rhs[ reg.idx ]
        var2.idx <- which(PT$op == "~~" & PT$lhs %in% reg.names &
                                          PT$rhs %in% reg.names &
                                          PT$lhs == PT$rhs)

        # eq constraints?
        # TODO!!
        keep.idx <- c(op.idx, var.idx, cov.idx, reg.idx, var2.idx) #,
                      #ov.int, lv.int)
        PT.1fac <- PT[keep.idx, , drop = FALSE]

        # clean up
        PT.1fac <- lav_partable_complete(PT.1fac)

        # fit 1-factor model
        #fit.1fac <- lavaan(PT.1fac, data = data, ...)
        fit.1fac <- do.call("lavaan",
                            args =  c(list(model  = PT.1fac,
                                           data   = data), dotdotdot2) )

        # fs.method?
        if(fsr.method == "skrondal.laake") {
            # dependent -> Bartlett
            if(FAC %in% eqs.y.names) {
                fs.method <- "Bartlett"
            } else {
                fs.method <- "regression"
            }
        }

        if(fsr.method %in% c("croon") || FIT@Options$se == "robust.sem") {
            fsm <- TRUE
        } else {
            fsm <- FALSE
        }

        # compute factor scores
        SC <- lav_predict_eta(fit.1fac, method = fs.method, fsm = fsm)

        for(g in 1:ngroups) {

            if(fsr.method %in% c("croon") || FIT@Options$se == "robust.sem") {
                FSM <- attr(SC, "fsm")
                attr(SC, "fsm") <- NULL
                FS.SCORES[[g]][[f]] <- SC[[g]]
  
                lambda.idx <- which(names(fit.1fac@Model@GLIST) == "lambda")
                theta.idx  <- which(names(fit.1fac@Model@GLIST) == "theta")
                LVINFO[[g]][[f]] <- 
                    list(fsm = FSM[[g]], 
                         lambda = fit.1fac@Model@GLIST[[lambda.idx[g]]],
                         theta  = fit.1fac@Model@GLIST[[theta.idx[g]]])
            } else {
                FS.SCORES[[g]][[f]] <- SC[[g]]
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
        if(FIT@Options$likelihood == "normal") {
            Ng <- lavInspect(FIT, "nobs")[g]
            COV <- COV * (Ng - 1) / Ng
        }
        COV
    })
    FSR.COV <- FS.COV


    # STEP 2:
    # construct path analysis model (removing all measurement elements)
    PT.PA <- PT
    PT.PA$est <- PT.PA$se <- NULL

    # extract all regressions
    reg.idx <- which(PT$op == "~" & PT$lhs %in% eqs.names &
                                    PT$rhs %in% eqs.names)

    # the variances
    var.idx <- which(PT$op == "~~" & PT$lhs %in% eqs.names &
                                     PT$rhs %in% eqs.names &
                                     PT$lhs == PT$rhs)

    # optionally covariances (exo!)
    cov.idx <- which(PT$op == "~~" & PT$lhs %in% eqs.names &
                                     PT$rhs %in% eqs.names &
                                     PT$lhs != PT$rhs)

    # means/intercepts
    int.idx <- which(PT$op == "~1" & PT$lhs %in% eqs.names)

    keep.idx <- c(reg.idx, var.idx, cov.idx, int.idx)
    PT.PA <- PT.PA[keep.idx, , drop = FALSE]

    # free all means/intercepts
    int.idx <- which(PT.PA$op == "~1")
    PT.PA$free[int.idx] <- 1L
    PT.PA$ustart[int.idx] <- NA

    # what about equality constraints?
    # TODO
    # what about inequality constraints?
    if(any(PT.PA$op %in% c(">", "<"))) {
        stop("lavaan ERROR: fsr does not support inequality constraints")
    }

    PT.PA <- lav_partable_complete(PT.PA)

    # adjust using Croon method
    if(fsr.method %in% c("croon")) {
        # compute 'corrected' COV for the factor scores 
        # using the Croon method
        for(g in 1:ngroups) {

            # correct covariances only
            if(fs.method != "bartlett") {
                for(i in 1:(nfac-1)) {
                    LHS <- lv.names[i]

                    A.y <- LVINFO[[g]][[LHS]]$fsm
                    lambda.y <- LVINFO[[g]][[LHS]]$lambda

                    for(j in (i+1):nfac) {
                        RHS <- lv.names[j]

                        A.x <- LVINFO[[g]][[RHS]]$fsm
                        lambda.x <- LVINFO[[g]][[RHS]]$lambda
                
                        # always 1 if Bartlett
                        A.xy <- as.numeric(crossprod(A.x %*% lambda.x,
                                                     A.y %*% lambda.y))

                        # corrected covariance
                        FSR.COV[[g]][i,j] <- FSR.COV[[g]][j,i] <- 
                            FS.COV[[g]][LHS,RHS] / A.xy
                    }
                }
            }

            # correct variances
            for(i in 1:nfac) {
                RHS <- lv.names[i]
                
                A.x <- LVINFO[[g]][[RHS]]$fsm
                lambda.x <- LVINFO[[g]][[RHS]]$lambda
                theta.x <- LVINFO[[g]][[RHS]]$theta

                if(fs.method == "bartlett") {
                    A.xx <- 1.0
                } else {
                    A.xx <- as.numeric(crossprod(A.x %*% lambda.x))
                }

                offset.x <- as.numeric(A.x %*% theta.x %*% t(A.x))

                FSR.COV[[g]][i,i] <- (FS.COV[[g]][RHS, RHS] - offset.x)/A.xx
            }
        } # g

    } # croon


    # Step 3: sem using factor scores
    if(fsr.method == "naive") {

        # FIXME!!! rbind FS.SCORES in the multiple group case!!!!
        if(ngroups > 1L) {
            stop("lavaan ERROR: fsr code not ready for multiple groups (yet)")
        }
        FS.SCORES <- as.data.frame(FS.SCORES[[1]])

        if(FIT@Options$se == "robust.sem") {
            # compute Omega.y (using NT for now)
            DATA <- lavInspect(FIT, "data")
            Omega.y <- lav_samplestats_Gamma_NT(Y             = DATA,
                                                meanstructure = TRUE,
                                                rescale       = TRUE,
                                                fixed.x       = FALSE)

            ## FIXME: this should be the Jacobian!!
            A <- lav_matrix_bdiag(lapply(LVINFO[[1]], "[[", "fsm"))
            A11 <- A
            A22 <- lav_matrix_duplication_post(
                   lav_matrix_duplication_ginv_pre(A %x% A))
            A.tilde <- lav_matrix_bdiag(A11, A22)
            Omega.f <- A.tilde %*% Omega.y %*% t(A.tilde)

            # add factor scores to data.frame
            fit <- lavaan(PT.PA, data = cbind(data, FS.SCORES),
                          meanstructure = TRUE,
                          NACOV = Omega.f,
                          se = "robust",
                          fixed.x = FALSE)
        } else {
            # add factor scores to data.frame
            fit <- lavaan(PT.PA, data = cbind(data, FS.SCORES), ...)
        }

    } else if(fsr.method == "skrondal.laake") {

        # FIXME!!! rbind FS.SCORES in the multiple group case!!!!
        if(ngroups > 1L) {
            stop("lavaan ERROR: fsr code not ready for multiple groups (yet)")
        }
        FS.SCORES <- as.data.frame(FS.SCORES[[1]])

        # apply bias-avoiding method
        fit <- lavaan(PT.PA, data = cbind(data, FS.SCORES), ...)

    } else if(fsr.method == "croon") {

        # apply bias-correcting method

        # transform FS.SCORES
        for(g in 1:ngroups) {
            OLD.inv <- solve(FS.COV[[g]]) 
            OLD.inv.sqrt <- lav_matrix_symmetric_sqrt(OLD.inv)
            FSR.COV.sqrt <- lav_matrix_symmetric_sqrt(FSR.COV[[g]])
            SC <- as.matrix(FS.SCORES[[g]])
            SC <- SC %*% OLD.inv.sqrt %*% FSR.COV.sqrt
            SC <- as.data.frame(SC)
            names(SC) <- lv.names
            FSR.SCORES[[g]] <- SC
        }

        # FIXME!!! rbind FS.SCORES in the multiple group case!!!!
        if(ngroups > 1L) {
            stop("lavaan ERROR: fsr code not ready for multiple groups (yet)")
        }
        FSR.SCORES <- as.data.frame(FSR.SCORES[[1]])

        # compute Omega.y (using NT for now)
        DATA <- lavInspect(FIT, "data")
        Omega.y <- lav_samplestats_Gamma_NT(Y             = DATA, 
                                            meanstructure = TRUE, 
                                            rescale       = TRUE,
                                            fixed.x       = FALSE)
 
        # factor score matrices
        A <- lav_matrix_bdiag(lapply(LVINFO[[1]], "[[", "fsm"))
        # compensate for Croon correction
        if(fs.method == "regression") {
            A <- OLD.inv.sqrt %*% FSR.COV.sqrt %*% A
        }

        # mean + vech(sigma)
        A11 <- A
        A22 <- lav_matrix_duplication_post(
               lav_matrix_duplication_ginv_pre(A %x% A))
        A.tilde <- lav_matrix_bdiag(A11, A22)
        Omega.f <- A.tilde %*% Omega.y %*% t(A.tilde)

        # add factor scores to data.frame
        fit <- lavaan(PT.PA, data = cbind(data, FSR.SCORES),
                      meanstructure = TRUE,
                      NACOV = Omega.f,
                      se = "robust",
                      fixed.x = FALSE)

    } else {
        stop("lavaan ERROR: fsr.method [", fsr.method, "] unknown", sep="")
    }

    # use 'external' slot to stores info
    fit@external <- list( FS.COV =  FS.COV,  FS.SCORES =  FS.SCORES, 
                         FSR.COV = FSR.COV, FSR.SCORES = FSR.SCORES,
                         LVINFO = LVINFO)
            
    fit
}

