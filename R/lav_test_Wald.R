# classic Wald test
#
lavTestWald <- function(object, constraints = NULL, verbose = FALSE) {

    if(object@Fit@npar > 0L && !object@Fit@converged)
        stop("lavaan ERROR: model did not converge")

    # remove == constraints from parTable
    PT <- as.data.frame(object@ParTable, stringsAsFactors = FALSE)
    eq.idx <- which(PT$op == "==")
    if(length(eq.idx) > 0L) {
        PT <- PT[-eq.idx,]
    }
    partable <- as.list(PT)

    # parse constraints
    FLAT <- lavParseModelString( constraints ); CON <- attr(FLAT, "constraints")
    LIST <- list()
    if(length(CON) > 0L) {
        lhs = unlist(lapply(CON, "[[", "lhs"))
         op = unlist(lapply(CON, "[[",  "op"))
        rhs = unlist(lapply(CON, "[[", "rhs"))
        LIST$lhs        <- c(LIST$lhs,        lhs)
        LIST$op         <- c(LIST$op,         op)
        LIST$rhs        <- c(LIST$rhs,        rhs)
    } else {
        stop("lavaan ERROR: no equality constraints found in constraints argument")
    }

    # theta = free parameters only
    theta <- object@Fit@x

    # build constraint function
    ceq.function <- lav_partable_constraints_ceq(partable = partable,
                                                 con = LIST, debug = FALSE)
    # compute jacobian restrictions
    JAC <- try(lavJacobianC(func = ceq.function, x = theta), silent=TRUE)
    if(inherits(JAC, "try-error")) { # eg. pnorm()
        JAC <- lavJacobianD(func = ceq.function, x = theta)
    }

    if(verbose) {
        cat("Restriction matrix (jacobian):\n"); print(JAC)
    }

    # linear restriction
    theta.r <- ceq.function( theta )

    if(verbose) {
        cat("Restricted theta values:\n"); print(theta.r)
    }

    # get VCOV
    VCOV <- vcov(object, labels = FALSE)

    # restricted vcov
    info.r  <- JAC %*% VCOV %*% t(JAC)
   
    # Wald test statistic
    Wald <- as.numeric(t(theta.r) %*% solve( info.r ) %*% theta.r)

    # df
    Wald.df <- nrow(JAC)

    # p-value based on chisq
    Wald.pvalue <- 1 - pchisq(Wald, df=Wald.df)

    list(stat=Wald, df=Wald.df, p.value=Wald.pvalue)
}
