# these functions need serious revision, and move to lav_samplestats_gamma.R


# compute Omega.beta: the 'incomplete Gamma' using either 
# expected or observed 'information' (for Abeta, aka the information matrix
# of the saturated model), using either structured (model-implied) or
# unstructured estimates of Beta=(mu, Sigma)
#
# using the formulas of Savalei & Bentler (2009) and Yuan & Bentler (2000)
compute.Omega.Beta <- function(Sigma.hat=NULL, Mu.hat=NULL, X=NULL, M=NULL,
                               information="observed") {

    out <- compute.Abeta.Bbeta(Sigma.hat=Sigma.hat, Mu.hat=Mu.hat, X=X, M=M, 
                               Abeta=TRUE, Bbeta=TRUE, information=information)
    Abeta <- out$Abeta
    Bbeta <- out$Bbeta

    Abeta.inv <- solve(Abeta)
    Omega.Sw <- Abeta.inv %*% Bbeta %*% Abeta.inv

    Omega.Sw
}


compute.Abeta.Bbeta <- function(Sigma.hat=NULL, Mu.hat=NULL, 
                                X=NULL, M=NULL, Abeta=TRUE, Bbeta=TRUE,
                                information="observed") {
    if(is.null(Mu.hat)) {
        stop("Mu.hat=NULL not implemented yet; use meanstructure=TRUE")
    }

    if(is.null(M)) {
        type <- "full"
        if(is.null(X) || !is.matrix(X)) {
            stop("X is null or non-matrix")
        }
        ntotal <- nrow(X)
        nvar <- ncol(X)
        npatterns <- 1L
    } else {
        type <- "missing"
        ntotal <- sum(sapply(M, "[[", "freq"))
        nvar <- length(M[[1]]$var.idx)
        npatterns <- length(M)
    }

    if(Abeta) {
        Aj22 <- matrix(0, nvar*nvar, nvar*nvar)
        Aj11 <- matrix(0, nvar, nvar)
        Aj12 <- matrix(0, nvar, nvar*nvar)
    }

    if(Bbeta) {
        nstar <- nvar + nvar * (nvar + 1)/2
        Bj <- matrix(0, nstar, nstar)
    }

    for(j in 1:npatterns) {
        if(type == "missing") {
            Xp <- M[[j]][["X"]]
            SX <- M[[j]][["SY"]]
            MX <- M[[j]][["MY"]]
            nobs <- M[[j]][["freq"]]
            var.idx <- M[[j]][["var.idx"]]
            taoj <- diag(nvar)[var.idx, , drop = FALSE]
        } else {
            nobs <- ntotal
            Xp <- X
            MX <- colMeans(X)
            SX <- crossprod(X)/nobs - tcrossprod(MX)
            var.idx <- rep(TRUE, nvar)
            taoj <- diag(nvar)
        }
        Sigma.inv <- inv.chol(Sigma.hat[var.idx, var.idx], logdet = FALSE)
        Mj <- t(taoj) %*% Sigma.inv %*% taoj
        Mu <- Mu.hat[var.idx]
        TT <- SX + tcrossprod(MX - Mu)
        if(Abeta) {
            if(information == "observed") {
                #cat("Abeta: observed\n")
                hjbar <- t(taoj) %*% Sigma.inv %*% (MX - Mu)
                Hjbar <- t(taoj) %*% (Sigma.inv %*% TT %*% Sigma.inv) %*% taoj
                Aj22 <- Aj22 + nobs * (Mj %x% (Hjbar - 0.5 * Mj))
                Aj11 <- Aj11 + nobs * Mj
                Aj12 <- Aj12 + nobs * (Mj %x% t(hjbar))
            } else {
                #cat("Abeta: expected\n")
                Aj22 <- Aj22 + nobs * 0.5 * ( Mj %x% Mj )
                Aj11 <- Aj11 + nobs * Mj
            }
        }
        if(Bbeta) {
            for(i in 1:nobs) {
                dx.Mu    <- matrix(0, nvar, 1)
                dx.Sigma <- matrix(0, nvar, nvar)

                diff.i <- t(Xp[i,] - Mu)
                dx.Mu[var.idx, 1] <- -1 * t(diff.i %*% Sigma.inv)
                dx.Sigma[var.idx, var.idx] <- (-1 * (Sigma.inv %*% 
                  (crossprod(diff.i) - Sigma.hat[var.idx, var.idx]) %*% 
                  Sigma.inv))
                diag(dx.Sigma) <- diag(dx.Sigma)/2
                # in lavaan: first the means, then the covariances
                dx <- c(dx.Mu, lav_matrix_vech(dx.Sigma))
                Bj <- Bj + tcrossprod(dx)
            }
        }
    }

    A.beta <- NULL
    B.beta <- NULL
    if (Abeta) {
        Abeta22 <- lav_matrix_duplication_pre_post(Aj22)
        Abeta12 <- lav_matrix_duplication_post(Aj12); Abeta21 <- t(Abeta12)
        Abeta11 <- Aj11
        A.beta  <- 1/ntotal * rbind( cbind(Abeta11, Abeta12), 
                                     cbind(Abeta21, Abeta22)  )
    }
    if (Bbeta) {
        B.beta <- 1/ntotal * Bj
    }

    list(Abeta = A.beta, Bbeta = B.beta)
}


# shortcut to only get Abeta
compute.Abeta <- function(Sigma.hat=NULL, Mu.hat=NULL, lavsamplestats=NULL,
                          lavdata=NULL, group=1L, information="observed") { 
    if(lavsamplestats@missing.flag) {
        X <- NULL
        M <- lavsamplestats@missing[[group]]
    } else {
        X <- lavdata@X[[group]]
        M <- NULL
    }
    out <- compute.Abeta.Bbeta(Sigma.hat=Sigma.hat, Mu.hat=Mu.hat, 
                               X=X, M=M, Abeta=TRUE, Bbeta=FALSE, 
                               information=information)
    Abeta <- out$Abeta
    Abeta
}

compute.Abeta.complete <- function(Sigma.hat=NULL, meanstructure=TRUE) {

    inv <- attr(Sigma.hat, "inv")
    if(is.matrix(inv)) {
        Sigma.hat.inv <- inv
    } else {
        Sigma.hat.inv <- inv.chol(Sigma.hat)
    }

    A <- 0.5 * lav_matrix_duplication_pre_post(Sigma.hat.inv %x% Sigma.hat.inv)

    if(meanstructure) {
        A11 <- Sigma.hat.inv
        A <- lav_matrix_bdiag(A11, A)
    }

    A
}

# shortcut if Sigma and Mu are simply the sample counterparts
compute.A1.sample <- function(lavsamplestats, group=1L, meanstructure=TRUE,
                              idx=NULL, information=NULL) {

    # note: for complete data, the type of information does not matter
    # but for incomplete data, it makes a difference!

    if(lavsamplestats@missing.flag) {
        # incomplete data
        A1 <- compute.Abeta(Sigma.hat=lavsamplestats@missing.h1[[group]]$sigma,
                            Mu.hat=lavsamplestats@missing.h1[[group]]$mu,
                            lavsamplestats=lavsamplestats, group=group,
                            information=information)
    } else {
        # complete data
        sample.icov <- lavsamplestats@icov[[group]]
        sample.mean <- lavsamplestats@mean[[group]]

        # do only subset
        if(length(idx) > 0) {
            sample.icov <- sample.icov[idx,idx]
            sample.mean <- sample.mean[idx]
        }

        A1 <- 0.5 * lav_matrix_duplication_pre_post(sample.icov %x% sample.icov)

        if(meanstructure) {
            A11 <- sample.icov
            A1 <- lav_matrix_bdiag(A11, A1)
        }
    }

    A1
}

# shortcut to only get Bbeta
compute.Bbeta <- function(Sigma.hat=NULL, Mu.hat=NULL, lavsamplestats=NULL, 
                          lavdata=NULL, group=1L) {
    if(lavsamplestats@missing.flag) {
        X <- NULL
        M <- lavsamplestats@missing[[group]]
    } else {
        X <- lavdata@X[[group]]
        M <- NULL
    }
    out <- compute.Abeta.Bbeta(Sigma.hat=Sigma.hat, Mu.hat=Mu.hat, 
                               X=X, M=M, Abeta=FALSE, Bbeta=TRUE)
    Bbeta <- out$Bbeta
    Bbeta
}

compute.B1.sample <- function(lavsamplestats=NULL, lavdata=NULL,
                              group=1L, meanstructure=TRUE) {

    if(lavsamplestats@missing.flag) {
        stopifnot(meanstructure == TRUE)
        B1 <- compute.Bbeta(Sigma.hat=lavsamplestats@missing.h1[[group]]$sigma,
                            Mu.hat=lavsamplestats@missing.h1[[group]]$mu,
                            lavsamplestats=lavsamplestats, group=group)
    } else {
        # complete lavdata and sample values only: B1 = A1 %*% Gamma %*% A1
        A1 <- compute.A1.sample(lavsamplestats=lavsamplestats, group=group, 
                                meanstructure=meanstructure)
        Gamma <- lav_samplestats_Gamma(Y = lavdata@X[[group]],
                                       meanstructure = meanstructure)
        B1 <- A1 %*% Gamma %*% A1
    }

    B1
}


### NOT USED, NOT (ENTIRELY CORRECT) !!!
### ONLY CORRECT FOR SATURATED MU!!!
# compute.Bbeta.complete <- function(Sigma.hat=NULL, Mu.hat=NULL, X=NULL,
#                                    sample.cov=NULL, sample.mean=NULL,
#                                    meanstructure=TRUE) {
# 
#     # alternative `analytic' version (5 times faster?)

#     # see Yuan & Hayashi 1996 page 406 (Beta hat)
# 
#     # ONLY CORRECT FOR SATURATED MU!!!
# 
#     diff.mean <- as.matrix(sample.mean - Mu.hat)
#     TT <- sample.cov + tcrossprod(diff.mean)
#     Sigma.hat.inv <- attr(Sigma.hat, "inv")
#     W <- 0.5 * lav_matrix_duplication_pre_post(Sigma.hat.inv %x% Sigma.hat.inv)
# 
#     if(meanstructure) {
#         G11 <- TT
#     #    G12 <- compute.third.moment(X) 
#         G22 <- compute.Gamma1(X)
# 
#         B11 <- Sigma.hat.inv %*% G11 %*% Sigma.hat.inv
#         B12 <- Sigma.hat.inv %*% G12 %*% W
#         B21 <- t(B12)
# 
#         diff <- tcrossprod(lav_matrix_vech(TT) - lav_matrix_vech(Sigma.hat))
#         G0 <- G22 + diff
# 
#         B22 <- W %*% G0 %*% W
# 
#         B1 <- rbind( cbind(B11, B12),
#                      cbind(B21, B22) )
#     } else { 
#         # this gives slightly different results, compared
#         # to meanstructure=TRUE and is probably not correct?
#         # not used for now
#         G <- compute.Gamma1(X)
#         diff <- tcrossprod(lav_matrix_vech(sample.cov) - lav_matrix_vech(Sigma.hat))
#         G0 <- G + diff
#         B1 <- W %*% G0 %*% W
#     }
#  
#     B1
# }
