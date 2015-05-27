# classic score test (= Lagrange Multiplier test)
#
# 'add' contains new parameters that are currently not included in de model
# (aka fixed to zero), but should be released
#
# if empty, we test the internal == constraints
#
lavTestScore <- function(object, add = NULL, verbose = FALSE, warn = TRUE,
                         univariate = TRUE, cumulative = FALSE) {

    # check object
    stopifnot(inherits(object, "lavaan"))
    lavoptions <- object@Options

    if(object@Fit@npar > 0L && !object@Fit@converged) {
        stop("lavaan ERROR: model did not converge")
    }

    # check arguments
    if(cumulative) {
        univariate <- TRUE
    }

    # add.flag?
    if(is.null(add) || nchar(add) == 0L) {
        add.flag <- FALSE
    } else {
        add.flag <- TRUE
    }

    if(add.flag) {
        # extend model with extra set of parameters
        FIT <- lav_object_extended(object, add = add)
        score <- lavTech(FIT, "gradient")
        information <- lavTech(FIT, "information.expected")

        npar <- object@Model@nx.free
        nadd <- FIT@Model@nx.free - npar

        # R
        R <- cbind(matrix(0, nrow = nadd, ncol = npar), diag(nadd))
        J.inv <- MASS:::ginv(information)

        # lhs/rhs
        lhs <- paste(".c", seq_len(nadd), sep = "")
        rhs <- rep("0", nadd)
        LABEL <- FIT@ParTable$label[ FIT@ParTable$user == 10L ]
        label.idx <- which(nchar(LABEL) > 0L)
        if(length(label.idx) > 0L) {
            lhs[ label.idx ] <- LABEL[ label.idx ]
        }
    } else {
        score <- lavTech(object, "gradient")
        information <- lavTech(object, "information.expected")
        J.inv <- MASS::ginv(information)
        R <- object@Model@con.jac[,]

        # lhs/rhs
        eq.idx <- which(object@ParTable$op == "==")
        if(length(eq.idx) > 0L) {
            lhs <- object@ParTable$lhs[eq.idx]
            rhs <- object@ParTable$rhs[eq.idx]
        }
    }

    if(nrow(R) == 0L) {
        stop("lavaan ERROR: now equality constraints found in model")
    }

    N <- nobs(object)
    if(lavoptions$mimic == "EQS") {
        N <- N - 1
    }
    
    if(lavoptions$se == "standard") {
        stat <- as.numeric(N * score %*% J.inv %*% score)
    } else {
        # generalized score test
        if(warn) {
            warning("lavaan WARNING: se is not `standard'; not implemented yet; falling back to ordinary score test")
        }
 
        # NOTE!!!
        # we can NOT use VCOV here, because it reflects the constraints,
        # and the whole point is to test for these constraints...
        
        stat <- as.numeric(N * score %*% J.inv %*% score)
    }

    # compute df, taking into account that some of the constraints may
    # be needed to identify the model (and hence information is singular)
    information.plus <- information + crossprod(R)
    df <- qr(R)$rank + qr(information)$rank - qr(information.plus)$rank
    pvalue <- 1 - pchisq(stat, df=df)

    OUT <- list(stat = stat, df = df, p.value = pvalue, se = lavoptions$se,
                lhs = lhs, rhs = rhs)

    if(univariate) {
        TS <- numeric( nrow(R) )
        for(r in 1:nrow(R)) {
            R1 <- R[-r,,drop = FALSE]
            Z1 <- cbind( rbind(information, R1),
                         rbind(t(R1), matrix(0,nrow(R1),nrow(R1))) )
            Z1.plus <- MASS::ginv(Z1)
            Z1.plus1 <- Z1.plus[ 1:nrow(information), 1:nrow(information) ]
            TS[r] <- as.numeric(N * t(score) %*%  Z1.plus1 %*% score)
        }

        OUT$TS.univariate <- TS
    }

    if(cumulative) {
        TS.order <- sort.int(TS, index.return = TRUE, decreasing = TRUE)$ix
        TS <- numeric( nrow(R) )
        for(r in 1:nrow(R)) {
            r.idx <- TS.order[1:r]

            R1 <- R[-r.idx,,drop = FALSE]
            Z1 <- cbind( rbind(information, R1),
                         rbind(t(R1), matrix(0,nrow(R1),nrow(R1))) )
            Z1.plus <- MASS::ginv(Z1)
            Z1.plus1 <- Z1.plus[ 1:nrow(information), 1:nrow(information) ]
            TS[r] <- as.numeric(N * t(score) %*%  Z1.plus1 %*% score)
        }

        OUT$TS.order <- TS.order
        OUT$TS.cumulative <- TS
    }

    OUT
}
