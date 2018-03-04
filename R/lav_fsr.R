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

        # correct covariances only
        if(fs.method != "bartlett") {
            for(i in 1:(nfac-1)) {
                if(i == 0L) break

                A.y <- LVINFO[[g]][[i]]$fsm
                lambda.y <- LVINFO[[g]][[i]]$lambda

                if(nfac > 1L) {
                    for(j in (i+1):nfac) {
    
                        A.x <- LVINFO[[g]][[j]]$fsm
                        lambda.x <- LVINFO[[g]][[j]]$lambda

                        # always 1 if Bartlett
                        A.xy <- as.numeric(crossprod(A.x %*% lambda.x,
                                                     A.y %*% lambda.y))

                        # corrected covariance
                        FSR.COV[[g]][i,j] <- FSR.COV[[g]][j,i] <-
                        FS.COV[[g]][i,j] / A.xy
                    }
                } # nfac > 1L
            } # i
        } 

        # correct variances
        for(i in 1:nfac) {
            A.x <- LVINFO[[g]][[i]]$fsm
            lambda.x <- LVINFO[[g]][[i]]$lambda
            theta.x <- LVINFO[[g]][[i]]$theta

            if(fs.method == "bartlett") {
                A.xx <- 1.0
            } else {
                A.xx <- as.numeric(crossprod(A.x %*% lambda.x))
            }

            offset.x <- as.numeric(A.x %*% theta.x %*% t(A.x))

            FSR.COV[[g]][i,i] <- (FS.COV[[g]][i, i] - offset.x)/A.xx
        }
    } # g

    FSR.COV
}


# simple version of Croon's correction:
# - we always assume ML/Bartlett factor scores, so there is no need to
#   adjust the covariances
# - we simply replace the variances (of the observed factor scores) by
#   the estimated variances of the factors (PSI elements)
lav_fsr_simple_correction <- function(FS.COV, LVINFO, mm.list = NULL,
                                      force.pd = FALSE) {

    # ngroups
    ngroups <- length(FS.COV)

    # FSR.COV
    FSR.COV <- FS.COV

    for(g in 1:ngroups) {

        # number of measurement blocks
        nblocks <- length(mm.list)

        # correct variances only
        for(b in seq_len(nblocks)) {
            psi.x <- LVINFO[[g]][[b]]$psi
            idx <- match(mm.list[[b]], unlist(mm.list, use.names = FALSE))
            FSR.COV[[g]][idx, idx] <- psi.x
        }

        # force pd?
        if(force.pd) {
            tmp <- FSR.COV[[g]]
            tmp2 <- lav_matrix_symmetric_force_pd(tmp)
            FSR.COV[[g]][] <- tmp2
        }

    } # g

    FSR.COV
}

