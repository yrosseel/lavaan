# factor score regression

# four methods:
#  - regression fsr (biased)
#  - Bartlett fsr (biased)
#  - Skrondal & Laake (2001)
#  - Croon (2002)

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
    } else if(fsr.method == "croonb") {
        # nothing to do
    }

    fs.method <- tolower(fs.method)

    # dot dot dot
    dotdotdot <- list(...)

    # check for arguments that we do not want (eg sample.cov)?
    # TODO

    # process model, no fitting
    FIT <- do.call(cmd, 
                   args =  c(list(model  = model, 
                                  data   = data,
                                  do.fit = FALSE), dotdotdot) )

    # any `regular' latent variables?
    lv.names <- unique(unlist(FIT@pta$vnames$lv.regular))
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
            stop("lavaan ERROR: mediational relationship are (currently) not allowed for the skrondal.laake method")
        }
    }


    # STEP 1:
    # compute factor scores, per latent variable
    SCORES <- vector("list", length = length(lv.names))
    if(fsr.method %in% c("croon", "croonb")) { 
        CROON <- vector("list", length = length(lv.names))
    }
    for(f in 1:length(lv.names)) {
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
        keep.idx <- c(op.idx, var.idx, cov.idx, reg.idx, var2.idx)
        PT.1fac <- PT[keep.idx, , drop = FALSE]

        # clean up
        PT.1fac <- lav_partable_complete(PT.1fac)

        # fit 1-factor model
        fit.1fac <- lavaan(PT.1fac, data = data, ...)

        # fs.method?
        if(fsr.method == "skrondal.laake") {
            # dependent -> Bartlett
            if(FAC %in% eqs.y.names) {
                fs.method <- "Bartlett"
            } else {
                fs.method <- "regression"
            }
        }

        if(fsr.method == "croonb") {
            fs.method <- "Bartlett"
        }

        if(fsr.method %in% c("croon" ,"croonb")) {
            fsm <- TRUE
        } else {
            fsm <- FALSE
        }

        # compute factor scores
        SC <- lavPredict(fit.1fac, type = "lv", method = fs.method, fsm = fsm)

        if(fsr.method %in% c("croon", "croonb")) {
            FSM <- attr(SC, "fsm")
            attr(SC, "fsm") <- NULL
            SCORES[[f]] <- SC

            lambda.idx <- which(names(fit.1fac@Model@GLIST) == "lambda")
            theta.idx  <- which(names(fit.1fac@Model@GLIST) == "theta")
            CROON[[f]] <- list(fsm = FSM, 
                               lambda = fit.1fac@Model@GLIST[lambda.idx],
                               theta  = fit.1fac@Model@GLIST[theta.idx])
        } else {
            SCORES[[f]] <- SC
        }
        
    }

    # FIXME!!! rbind SCORES in the multiple group case!!!!
    if(FIT@Data@ngroups > 1L) {
        stop("lavaan ERROR: fsr code not ready for multiple groups (yet)")
    }

    names(SCORES) <- lv.names
    SCORES <- as.data.frame(SCORES)
    if(fsr.method %in% c("croon", "croonb")) {
        names(CROON) <- lv.names
    }


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

    keep.idx <- c(reg.idx, var.idx, cov.idx)
    PT.PA <- PT.PA[keep.idx, , drop = FALSE]

    # what about equality constraints?
    # TODO
    # what about inequality constraints?
    if(any(PT.PA$op %in% c(">", "<"))) {
        stop("lavaan ERROR: fsr does not support inequality constraints")
    }

    PT.PA <- lav_partable_complete(PT.PA)

    if(fsr.method == "naive") {
        # add factor scores to data.frame
        fit <- lavaan(PT.PA, data = cbind(data, SCORES), ...)
    } else if(fsr.method == "skrondal.laake") {
        # apply bias-avoiding method
        fit <- lavaan(PT.PA, data = cbind(data, SCORES), ...)
    } else if(fsr.method == "croon") {
        # apply bias-correcting method
        #fit <- lavaan(PT.PA, data = cbind(data, SCORES), ...)

        # here, we first do this manually; later, we will create a  
        # lavaan object

        PE <- vector("list", length = FIT@Data@ngroups)
        for(g in FIT@Data@ngroups) {

            # FIXME: we assume only 1 group
            COV.SC <- cov(SCORES) ## divided by N-1!!!
            N <- nobs(FIT)
            COV.SC <- COV.SC * (N-1) / N
            all.names <- rownames(COV.SC)

            NEW <- COV.SC

            # correct covariances only
            if(fs.method != "bartlett") {
                for(i in 1:ncol(COV.SC)) {
                    LHS <- all.names[i]

                    A.y <- CROON[[LHS]]$fsm[[g]]
                    lambda.y <- CROON[[LHS]]$lambda[[g]]

                    for(j in (i+1):nrow(COV.SC)) {
                        RHS <- all.names[j]

                        A.x <- CROON[[RHS]]$fsm[[g]]
                        lambda.x <- CROON[[RHS]]$lambda[[g]]
                
                        # always 1 if Bartlett
                        A.xy <- as.numeric(crossprod(A.x %*% lambda.x,
                                                     A.y %*% lambda.y))

                        # corrected covariance
                        NEW[i,j] <- NEW[j,i] <- COV.SC[LHS,RHS] / A.xy
                    }
                }
            }

            # correct variances
            for(i in 1:ncol(COV.SC)) {
                RHS <- all.names[i]
                if(!RHS %in% eqs.x.names) next
                
                A.x <- CROON[[RHS]]$fsm[[g]]
                lambda.x <- CROON[[RHS]]$lambda[[g]]
                theta.x <- CROON[[RHS]]$theta[[g]]

                if(fs.method == "bartlett") {
                    A.xx <- 1.0
                } else {
                    A.xx <- as.numeric(crossprod(A.x %*% lambda.x))
                }

                offset.x <- as.numeric(A.x %*% theta.x %*% t(A.x))

                NEW[i,i] <- (COV.SC[RHS, RHS] - offset.x)/A.xx
            }



            lhs <- op <- rhs <- character(0L)
            est <- se <- numeric(0L)
            for(y in eqs.y.names) {
                reg.idx <- which(PT.PA$op == "~" & PT.PA$group == g 
                                                 & PT.PA$lhs == y)
                nx <- length(reg.idx)
                x.names <- PT.PA$rhs[reg.idx]
                EST <- solve(NEW[x.names, x.names, drop = FALSE], 
                             NEW[x.names, y,       drop = FALSE])

                lhs <- c(lhs, rep(y, nx))
                 op <- c(op, rep("~", nx))
                rhs <- c(rhs, x.names)
                est <- c(est, as.numeric(EST))
                 se <- c(se, rep(NA, nx))
            }
            tmp.se <- ifelse( se == 0.0, NA, se)
            z <- est/tmp.se
            pvalue <- 2 * (1 - pnorm( abs(z) ))
            PE[[g]] <- data.frame(lhs = lhs, op = op, rhs = rhs,
                                  est = est, se = se, z = z, pvalue = pvalue)
            class(PE[[g]]) <- c("lavaan.data.frame", "data.frame")
            attr(PE[[g]], "header") <- 
                paste("lavaan factor score regression\n",
                      "  fsr.method = ", fsr.method, "\n",
                      "  fs.method  = ", fs.method, "\n",  sep = "")
        }
        if(FIT@Data@ngroups == 1L) {
            fit <- PE[[1]]
        } else {
            fit <- PE
        }
    } else if(fsr.method == "croonb") {
 
        N <- nobs(FIT)

        # rescale x-scores
        for(g in FIT@Data@ngroups) {
            for(RHS in eqs.x.names) {
                A.x <- CROON[[RHS]]$fsm[[g]]
                theta.x <- CROON[[RHS]]$theta[[g]]
                offset.x <- as.numeric(A.x %*% theta.x %*% t(A.x))
                OLD.varx <- var(SCORES[,RHS]) # divided by N-1!!!
                OLD.varx <- OLD.varx * (N - 1) / N
                NEW.varx <- OLD.varx - offset.x
                scale.factor <- NEW.varx / OLD.varx

                SCORES[,RHS] <- scale.factor * SCORES[,RHS]
            }
        }
        # add factor scores to data.frame
        fit <- lavaan(PT.PA, data = cbind(data, SCORES), ...)
    } else {
        stop("lavaan ERROR: fsr.method [", fsr.method, "] unknown", sep="")
    }
            
    fit
}

