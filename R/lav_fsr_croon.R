# use the `Croon' method to correct the covariance matrix
# of the factor scores
lav_fsr_croon_correction <- function(FS.COV, LVINFO, fs.method = "bartlett") {

    # ngroups
    ngroups <- length(FS.COV)

    # FSR.COV
    FSR.COV <- FS.COV

    for(g in 1:ngroups) {

        # number of factors - lv.names
        nfac <- nrow(FS.COV[[g]])
        lv.names <- names(LVINFO[[g]])

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

    FSR.COV
}

