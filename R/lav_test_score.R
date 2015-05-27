# classic score test (= Lagrange Multiplier test)
#
# this function can run in two modes:
# 
# MODE 1: 'add'
#   add new parameters that are currently not included in de model
#   (aka fixed to zero), but should be released
#
# MODE 2: 'release' (the default)
#   release existing "==" constraints
#
lavTestScore <- function(object, add = NULL, release = NULL,
                         univariate = TRUE, cumulative = FALSE,
                         epc = FALSE, verbose = FALSE, warn = TRUE) {

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


    # Mode 1: ADDING new parameters
    if(!is.null(add) && nchar(add) > 0L) {
        # check release argument
        if(!is.null(release)) {
            stop("lavaan ERROR: `add' and `release' arguments can be used together.")
        }

        # extend model with extra set of parameters
        FIT <- lav_object_extended(object, add = add)
        score <- lavTech(FIT, "gradient")
        information <- lavTech(FIT, "information.expected")

        npar <- object@Model@nx.free
        nadd <- FIT@Model@nx.free - npar

        # R
        R.model <- object@Model@con.jac[,]
        if(nrow(R.model) > 0L) {
            R.model <- cbind(R.model, matrix(0, nrow(R.model), ncol = nadd))
            R.add   <- cbind(matrix(0, nrow = nadd, ncol = npar), diag(nadd))
            R       <- rbind(R.model, R.add)

            Z <- cbind(rbind(information, R.model),
                       rbind(t(R.model),matrix(0,nrow(R.model),nrow(R.model))))
            Z.plus <- MASS::ginv(Z)
            J.inv  <- Z.plus[ 1:nrow(information), 1:nrow(information) ]

            r.idx <- seq_len(nadd) + nrow(R.model)
        } else {
            R <- cbind(matrix(0, nrow = nadd, ncol = npar), diag(nadd))
            J.inv <- MASS::ginv(information)

            r.idx <- seq_len(nadd)
        }

        # lhs/rhs
        lhs <- lav_partable_labels(FIT@ParTable)[ FIT@ParTable$user == 10L ]
        rhs <- rep("0", nadd)
    } else {
    # MODE 2: releasing constraints

        R <- object@Model@con.jac[,]
        if(nrow(R) == 0L) {
            stop("lavaan ERROR: no equality constraints found in model.")
        }

        score <- lavTech(object, "gradient")
        information <- lavTech(object, "information.expected")
        J.inv <- MASS::ginv(information)
        R <- object@Model@con.jac[,]

        if(is.null(release)) {
            # ALL constraints
            r.idx <- seq_len( nrow(R) )
        } else if(is.numeric(release)) {
            if(length(release) == 1L) {
                r.idx <- seq_len( release )
            } else {
                r.idx <- release
            }

            if(max(r.idx) > nrow(R)) {
                stop("lavaan ERROR: maximum constraint number (", max(r.idx),
                     ") is larger than number of constraints (", nrow(R), ")")
            }
        } else if(is.character(release)) {
            stop("not implemented yet")
        }

        # lhs/rhs
        eq.idx <- which(object@ParTable$op == "==")
        if(length(eq.idx) > 0L) {
            lhs <- object@ParTable$lhs[eq.idx][r.idx]
            rhs <- object@ParTable$rhs[eq.idx][r.idx]
        }
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
    df <- qr(R[r.idx,,drop = FALSE])$rank + 
              ( qr(information)$rank - qr(information.plus)$rank )
    pvalue <- 1 - pchisq(stat, df=df)

    OUT <- list(stat = stat, df = df, p.value = pvalue, se = lavoptions$se,
                lhs = lhs, rhs = rhs)

    if(univariate) {
        TS <- numeric( nrow(R) )
        for(r in r.idx) {
            R1 <- R[-r,,drop = FALSE]
            Z1 <- cbind( rbind(information, R1),
                         rbind(t(R1), matrix(0,nrow(R1),nrow(R1))) )
            Z1.plus <- MASS::ginv(Z1)
            Z1.plus1 <- Z1.plus[ 1:nrow(information), 1:nrow(information) ]
            TS[r] <- as.numeric(N * t(score) %*%  Z1.plus1 %*% score)
        }

        OUT$TS.univariate <- TS[r.idx]
    }

    if(cumulative) {
        TS.order <- sort.int(TS, index.return = TRUE, decreasing = TRUE)$ix
        TS <- numeric( length(r.idx) )
        for(r in 1:length(r.idx)) {
            rcumul.idx <- TS.order[1:r]

            R1 <- R[-rcumul.idx,,drop = FALSE]
            Z1 <- cbind( rbind(information, R1),
                         rbind(t(R1), matrix(0,nrow(R1),nrow(R1))) )
            Z1.plus <- MASS::ginv(Z1)
            Z1.plus1 <- Z1.plus[ 1:nrow(information), 1:nrow(information) ]
            TS[r] <- as.numeric(N * t(score) %*%  Z1.plus1 %*% score)
        }

        OUT$TS.order <- TS.order
        OUT$TS.cumulative <- TS
    }

    if(epc) {
        EPC <- vector("list", length = length(r.idx))
        for(i in 1:length(r.idx)) {
            r <- r.idx[i]
            R1 <- R[-r,,drop = FALSE]
            Z1 <- cbind( rbind(information, R1),
                         rbind(t(R1), matrix(0,nrow(R1),nrow(R1))) )
            Z1.plus <- MASS::ginv(Z1)
            Z1.plus1 <- Z1.plus[ 1:nrow(information), 1:nrow(information) ]
            EPC[[i]] <- -1 * as.numeric(score %*%  Z1.plus1)
        }

        OUT$EPC <- EPC
    }

    OUT
}
