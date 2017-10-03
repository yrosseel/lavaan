# the multivariate normal distribution

# 1) loglikelihood (from raw data, or sample statistics)
# 2) derivatives with respect to mu, Sigma, vech(Sigma)
# 3) casewise scores with respect to mu, vech(Sigma), mu + vech(Sigma)
# 4) hessian mu + vech(Sigma)
# 5) information h0 mu + vech(Sigma)
#    5a: (unit)    expected information
#    5b: (unit)    observed information
#    5c: (unit) first.order information
# 6) inverted information h0 mu + vech(Sigma)
#    6a: (unit) inverted expected information
#    6b: /
#    6c: /
# 7) ACOV h0 mu + vech(Sigma)
#    7a: 1/N * inverted expected    information
#    7b: 1/N * inverted observed    information
#    7c: 1/N * inverted first-order information
#    7d: sandwich acov

# YR 07 Feb 2016: first version
# YR 24 Mar 2016: added firstorder information, hessian logl
# YR 19 Jan 2017: added lav_mvnorm_inverted_information_expected

# 0. densities 
lav_mvnorm_dmvnorm <- function(Y             = NULL,
                               wt            = NULL,
                               Mu            = NULL,
                               Sigma         = NULL,
                               Sigma.inv     = NULL,
                               Sinv.method   = "eigen",
                               log           = TRUE) {

    if(is.matrix(Y)) {
        if(is.null(Mu) && is.null(Sigma) && is.null(Sigma.inv)) {
            out <- lav_mvnorm_loglik_data_z(Y = Y, casewise = TRUE)
        } else {
            out <- lav_mvnorm_loglik_data(Y = Y, Mu = Mu, Sigma = Sigma,
                                          casewise = TRUE,
                                          Sinv.method = Sinv.method)
        }
    } else {
        # just one
        P <- length(Y); LOG.2PI <- log(2 * pi)

        if(is.null(Mu) && is.null(Sigma) && is.null(Sigma.inv)) {
            # mahalanobis distance
            DIST <- sum(Y * Y)
            out <- -(P * LOG.2PI + DIST)/2
        } else {
            if(is.null(Sigma.inv)) {
                Sigma.inv <- lav_matrix_symmetric_inverse(S = Sigma, 
                    logdet = TRUE, Sinv.method = Sinv.method)
                logdet <- attr(Sigma.inv, "logdet")
            } else {
                logdet <- attr(Sigma.inv, "logdet")
                if(is.null(logdet)) {
                    # compute - ln|Sigma.inv|
                    ev <- eigen(Sigma.inv, symmetric = TRUE, only.values = TRUE)
                    logdet <- -1 * sum(log(ev$values))
                }
            }

            # mahalanobis distance
            Yc <- Y - Mu
            DIST <- sum(Yc %*% Sigma.inv * Yc)
            out <- -(P * LOG.2PI + logdet + DIST)/2
        }
    }

    if(!is.null(wt)) {
        out <- out * wt
    }

    if(!log) {
        out <- exp(out)
    }

    out
}

# 1. likelihood

# 1a: input is raw data
# (note casewise = TRUE same as: dmvnorm(Y, mean, sigma, log = TRUE))
lav_mvnorm_loglik_data <- function(Y             = NULL,
                                   wt            = NULL,
                                   Mu            = NULL,
                                   Sigma         = NULL,
                                   casewise      = FALSE,
                                   Sinv.method   = "eigen") {

    if(!is.null(wt)) {
        N <- sum(wt)
    } else {
        N <- NROW(Y)
    }

    P <- NCOL(Y); Mu <- as.numeric(Mu)

    if(casewise) {
        LOG.2PI <- log(2 * pi)

        # invert Sigma
        if(Sinv.method == "chol") {
            cS <- chol(Sigma); icS <- backsolve(cS, diag(P))
            Yc <- t( t(Y) - Mu )
            DIST <- rowSums((Yc %*% icS)^2)
            logdet <- -2 * sum(log(diag(icS)))
        } else {
            Sigma.inv <- lav_matrix_symmetric_inverse(S = Sigma, logdet = TRUE,
                                                      Sinv.method = Sinv.method)
            logdet <- attr(Sigma.inv, "logdet")
            # mahalanobis distance
            Yc <- t( t(Y) - Mu )
            DIST <- rowSums(Yc %*% Sigma.inv * Yc)
        }

        loglik <- -(P * LOG.2PI + logdet + DIST)/2

        # weights
        if(!is.null(wt)) {
            loglik <- loglik * wt
        }

    } else {
        # invert Sigma
        Sigma.inv <- lav_matrix_symmetric_inverse(S = Sigma, logdet = TRUE,
                                                  Sinv.method = Sinv.method)
        if(!is.null(wt)) {
            out <- stats::cov.wt(Y, wt = wt, method = "ML")
            sample.mean <- out$center
            sample.cov  <- out$cov
        } else {
            sample.mean <- colMeans(Y)
            sample.cov <- 1/N*crossprod(Y) - tcrossprod(sample.mean)
        }
        loglik <- lav_mvnorm_loglik_samplestats(sample.mean = sample.mean,
                                                sample.cov  = sample.cov,
                                                sample.nobs = N,
                                                Mu          = Mu,
                                                Sigma.inv   = Sigma.inv)
    }

    loglik
}



# 1b: input are sample statistics (mean, cov, N) only
lav_mvnorm_loglik_samplestats <- function(sample.mean  = NULL,
                                          sample.cov   = NULL,
                                          sample.nobs  = NULL,
                                          Mu           = NULL,
                                          Sigma        = NULL,
                                          Sinv.method  = "eigen",
                                          Sigma.inv    = NULL) {

    P <- length(sample.mean); N <- sample.nobs
    Mu <- as.numeric(Mu); sample.mean <- as.numeric(sample.mean)
    LOG.2PI <- log(2 * pi)

    if(is.null(Sigma.inv)) {
        Sigma.inv <- lav_matrix_symmetric_inverse(S = Sigma, logdet = TRUE,
                                                  Sinv.method = Sinv.method)
        logdet <- attr(Sigma.inv, "logdet")
    } else {
        logdet <- attr(Sigma.inv, "logdet")
        if(is.null(logdet)) {
            # compute - ln|Sigma.inv|
            ev <- eigen(Sigma.inv, symmetric = TRUE, only.values = TRUE)
            logdet <- -1 * sum(log(ev$values))
        }
    }

    # tr(Sigma^{-1} %*% S)
    DIST1 <- sum(Sigma.inv * sample.cov)
    # (ybar - mu)^T %*% Sigma.inv %*% (ybar - mu)
    Diff <- as.numeric(sample.mean - Mu)
    DIST2 <- sum(as.numeric(crossprod(Diff, Sigma.inv)) * Diff)

    loglik <- -N/2 * (P * LOG.2PI + logdet + DIST1 + DIST2)

    loglik
}

# 1c special case: Mu = 0, Sigma = I
lav_mvnorm_loglik_data_z <- function(Y             = NULL,
                                     wt            = NULL,
                                     casewise      = FALSE) {
    if(!is.null(wt)) {
        N <- sum(wt)
    } else {
        N <- NROW(Y)
    }

    P <- NCOL(Y); LOG.2PI <- log(2 * pi)
   
    if(casewise) {
        DIST <- rowSums(Y * Y)
        loglik <- -(P * LOG.2PI + DIST)/2
        if(!is.null(wt)) {
            loglik <- loglik * wt
        }
    } else {
        if(!is.null(wt)) {
            out <- stats::cov.wt(Y, wt = wt, method = "ML")
            sample.mean <- out$center
            sample.cov  <- out$cov
        } else {
            sample.mean <- colMeans(Y)
            sample.cov <- 1/N*crossprod(Y) - tcrossprod(sample.mean)
        }

        DIST1 <- sum(diag(sample.cov))
        DIST2 <- sum(sample.mean * sample.mean)

        loglik <- -N/2 * (P * LOG.2PI + DIST1 + DIST2)
    }

    loglik
}





# 2. Derivatives

# 2a: derivative logl with respect to mu
lav_mvnorm_dlogl_dmu <- function(Y           = NULL,
                                 wt          = NULL,
                                 Mu          = NULL,
                                 Sigma       = NULL,
                                 Sinv.method = "eigen",
                                 Sigma.inv   = NULL) {
    Mu <- as.numeric(Mu)
  
    if(is.null(Sigma.inv)) {
        # invert Sigma
        Sigma.inv <- lav_matrix_symmetric_inverse(S = Sigma, logdet = FALSE,
                                                  Sinv.method = Sinv.method)
    }

    # substract 'Mu' from Y
    Yc <- t( t(Y) - Mu )

    # weights
    if(!is.null(wt)) {
        Yc <- Yc * wt
    }

    # derivative
    dmu <- as.numeric(Sigma.inv %*% colSums(Yc))

    dmu
}

# 2b: derivative logl with respect to Sigma (full matrix, ignoring symmetry)
lav_mvnorm_dlogl_dSigma <- function(Y           = NULL,
                                    wt          = NULL,
                                    Mu          = NULL,
                                    Sigma       = NULL,
                                    Sinv.method = "eigen",
                                    Sigma.inv   = NULL) {

    if(!is.null(wt)) {
        N <- sum(wt)
    } else {
        N <- NROW(Y)
    }

    Mu <- as.numeric(Mu)

    if(is.null(Sigma.inv)) {
        # invert Sigma
        Sigma.inv <- lav_matrix_symmetric_inverse(S = Sigma, logdet = FALSE,
                                                  Sinv.method = Sinv.method)
    }

    # substract 'Mu' from Y
    Yc <- t( t(Y) - Mu )

    # W.tilde
    if(!is.null(wt)) {
        out <- stats::cov.wt(Y, wt = wt, method = "ML")
        SY <- out$cov
        MY <- out$center
        W.tilde <- SY + tcrossprod(MY - Mu)
    } else {
        W.tilde <- crossprod(Yc) / N
    }

    # derivative
    dSigma <- -(N/2)* (Sigma.inv - (Sigma.inv %*% W.tilde %*% Sigma.inv))

    dSigma
}

# 2c: derivative logl with respect to vech(Sigma)
lav_mvnorm_dlogl_dvechSigma <- function(Y           = NULL,
                                        wt          = NULL,
                                        Mu          = NULL,
                                        Sigma       = NULL,
                                        Sinv.method = "eigen",
                                        Sigma.inv   = NULL) {

    if(!is.null(wt)) {
        N <- sum(wt)
    } else {
        N <- NROW(Y)
    }
    
    Mu <- as.numeric(Mu)

    if(is.null(Sigma.inv)) {
        # invert Sigma
        Sigma.inv <- lav_matrix_symmetric_inverse(S = Sigma, logdet = FALSE,
                                                  Sinv.method = Sinv.method)
    }

    # substract 'Mu' from Y
    Yc <- t( t(Y) - Mu )

    # W.tilde
    if(!is.null(wt)) {
        out <- stats::cov.wt(Y, wt = wt, method = "ML")
        SY <- out$cov
        MY <- out$center
        W.tilde <- SY + tcrossprod(MY - Mu)
    } else {
        W.tilde <- crossprod(Yc) / N
    }

    # derivative (avoiding kronecker product)
    dSigma <- -(N/2)* (Sigma.inv - (Sigma.inv %*% W.tilde %*% Sigma.inv))
    dvechSigma <- as.numeric( lav_matrix_duplication_pre( 
                                  as.matrix(lav_matrix_vec(dSigma)) ) )

    dvechSigma
}

# 3. Casewise scores

# 3a: casewise scores with respect to mu
lav_mvnorm_scores_mu <- function(Y           = NULL,
                                 wt          = NULL,
                                 Mu          = NULL,
                                 Sigma       = NULL,
                                 Sinv.method = "eigen",
                                 Sigma.inv   = NULL) {
    Mu <- as.numeric(Mu)

    if(is.null(Sigma.inv)) {
        # invert Sigma
        Sigma.inv <- lav_matrix_symmetric_inverse(S = Sigma, logdet = FALSE,
                                                  Sinv.method = Sinv.method)
    }

    # substract Mu
    Yc <- t( t(Y) - Mu )

    # postmultiply with Sigma.inv
    SC <- Yc %*% Sigma.inv

    # weights
    if(!is.null(wt)) {
        SC <- SC * wt
    }

    SC
}

# 3b: casewise scores with respect to vech(Sigma)
lav_mvnorm_scores_vech_sigma <- function(Y           = NULL,
                                         wt          = NULL,
                                         Mu          = NULL,
                                         Sigma       = NULL,
                                         Sinv.method = "eigen",
                                         Sigma.inv   = NULL) {
    P <- NCOL(Y); Mu <- as.numeric(Mu)

    if(is.null(Sigma.inv)) {
        # invert Sigma
        Sigma.inv <- lav_matrix_symmetric_inverse(S = Sigma, logdet = FALSE,
                                                  Sinv.method = Sinv.method)
    }

    # vech(Sigma.inv)
    isigma <- lav_matrix_vech(Sigma.inv)

    # substract Mu
    Yc <- t( t(Y) - Mu )

    # postmultiply with Sigma.inv
    Yc <- Yc %*% Sigma.inv

    # tcrossprod
    idx1 <- lav_matrix_vech_col_idx(P); idx2 <- lav_matrix_vech_row_idx(P)
    Z <- Yc[,idx1] * Yc[,idx2]

    # substract isigma from each row
    SC <- t( t(Z) - isigma )

    # adjust for vech
    SC[,lav_matrix_diagh_idx(P)] <- SC[,lav_matrix_diagh_idx(P)] / 2

    # weights
    if(!is.null(wt)) {
        SC <- SC * wt
    }

    SC
}

# 3c: casewise scores with respect to mu + vech(Sigma)
lav_mvnorm_scores_mu_vech_sigma <- function(Y           = NULL,
                                            wt          = NULL,
                                            Mu          = NULL,
                                            Sigma       = NULL,
                                            Sinv.method = "eigen",
                                            Sigma.inv   = NULL) {
    P <- NCOL(Y); Mu <- as.numeric(Mu)

    if(is.null(Sigma.inv)) {
        # invert Sigma
        Sigma.inv <- lav_matrix_symmetric_inverse(S = Sigma, logdet = FALSE,
                                                  Sinv.method = Sinv.method)
    }

    # vech(Sigma.inv)
    isigma <- lav_matrix_vech(Sigma.inv)
    
    # substract Mu
    Yc <- t( t(Y) - Mu )

    # postmultiply with Sigma.inv
    Yc <- Yc %*% Sigma.inv
    
    # tcrossprod
    idx1 <- lav_matrix_vech_col_idx(P); idx2 <- lav_matrix_vech_row_idx(P)
    Z <- Yc[,idx1] * Yc[,idx2]
    
    # substract isigma from each row
    SC <- t( t(Z) - isigma )

    # adjust for lav_matrix_duplication_pre (not vech!)
    SC[,lav_matrix_diagh_idx(P)] <- SC[,lav_matrix_diagh_idx(P)] / 2
    
    out <- cbind(Yc, SC)

    # weights
    if(!is.null(wt)) {
        out <- out * wt
    }

    out
}


# 4. hessian of logl

# 4a: hessian logl Mu and vech(Sigma) from raw data
lav_mvnorm_logl_hessian_data <- function(Y           = NULL,
                                         wt          = NULL,
                                         Mu          = NULL,
                                         Sigma       = NULL,
                                         Sinv.method = "eigen",
                                         Sigma.inv   = NULL) {
    if(!is.null(wt)) {
        N <- sum(wt)
    } else {
        N <- NROW(Y)
    }

    # observed information
    observed <- lav_mvnorm_information_observed_data(Y = Y, wt = wt, Mu = Mu,
        Sigma = Sigma, Sinv.method = Sinv.method, Sigma.inv = Sigma.inv)

    -N*observed
}

# 4b: hessian Mu and vech(Sigma) from samplestats
lav_mvnorm_logl_hessian_samplestats <-
    function(sample.mean  = NULL,
             sample.cov   = NULL,
             sample.nobs  = NULL,
             Mu           = NULL,
             Sigma        = NULL,
             Sinv.method  = "eigen",
             Sigma.inv    = NULL) {

    N <- sample.nobs

    # observed information
    observed <- lav_mvnorm_information_observed_samplestats(sample.mean = 
        sample.mean, sample.cov = sample.cov, Mu = Mu, 
        Sigma = Sigma, Sinv.method = Sinv.method, Sigma.inv = Sigma.inv)
    
    -N*observed
}

# 5) Information h0

# 5a: unit expected information h0 Mu and vech(Sigma)
lav_mvnorm_information_expected <- function(Y             = NULL, # unused!
                                            wt            = NULL, # unused!
                                            Mu            = NULL, # unused!
                                            Sigma         = NULL,
                                            Sinv.method   = "eigen",
                                            Sigma.inv     = NULL,
                                            meanstructure = TRUE) {

    if(is.null(Sigma.inv)) {
        # invert Sigma
        Sigma.inv <- lav_matrix_symmetric_inverse(S = Sigma, logdet = FALSE,
                                                  Sinv.method = Sinv.method)
    }

    I22 <- 0.5 * lav_matrix_duplication_pre_post(Sigma.inv %x% Sigma.inv)

    if(meanstructure) {
        I11 <- Sigma.inv
        out <- lav_matrix_bdiag(I11, I22)
    } else {
        out <- I22
    }

    out
}

# 5b: unit observed information h0
lav_mvnorm_information_observed_data <- function(Y           = NULL,
                                                 wt          = NULL,
                                                 Mu          = NULL,
                                                 Sigma       = NULL,
                                                 Sinv.method = "eigen",
                                                 Sigma.inv   = NULL) {

    if(!is.null(wt)) {
        N <- sum(wt)
        out <- stats::cov.wt(Y, wt = wt, method = "ML")
        sample.cov  <- out$cov
        sample.mean <- out$center
    } else {
        N <- NROW(Y)
        # sample statistics
        sample.mean <- colMeans(Y)
        sample.cov <- 1/N*crossprod(Y) - tcrossprod(sample.mean)
    }

    lav_mvnorm_information_observed_samplestats(sample.mean = sample.mean,
        sample.cov = sample.cov, Mu = Mu, Sigma = Sigma,
        Sinv.method = Sinv.method, Sigma.inv = Sigma.inv)
}

# 5b-bis: observed information h0 from sample statistics
lav_mvnorm_information_observed_samplestats <-
    function(sample.mean  = NULL,
             sample.cov   = NULL,
             Mu           = NULL,
             Sigma        = NULL,
             Sinv.method  = "eigen",
             Sigma.inv    = NULL) {

    sample.mean <- as.numeric(sample.mean); Mu <- as.numeric(Mu)

    if(is.null(Sigma.inv)) {
        # invert Sigma
        Sigma.inv <- lav_matrix_symmetric_inverse(S = Sigma, logdet = FALSE,
                                                  Sinv.method = Sinv.method)
    }

    W.tilde <- sample.cov + tcrossprod(sample.mean - Mu)
    
    I11 <- Sigma.inv
    I21 <- lav_matrix_duplication_pre( (Sigma.inv %*% (sample.mean - Mu)) %x%
                                        Sigma.inv )
    I12 <- t(I21)
    
    AAA <- Sigma.inv %*% (2*W.tilde - Sigma) %*% Sigma.inv
    I22 <- (1/2) * lav_matrix_duplication_pre_post(Sigma.inv %x% AAA)

    rbind( cbind(I11, I12),
           cbind(I21, I22) )
}

# 5c: unit first-order information h0
lav_mvnorm_information_firstorder <- function(Y             = NULL,
                                              wt            = NULL,
                                              Mu            = NULL,
                                              Sigma         = NULL,
                                              Sinv.method   = "eigen",
                                              Sigma.inv     = NULL,
                                              meanstructure = TRUE) {

    if(!is.null(wt)) {
        N <- sum(wt)
    } else {
        N <- NROW(Y)
    }

    if(meanstructure) {
        SC <- lav_mvnorm_scores_mu_vech_sigma(Y = Y, wt = wt, 
                  Mu = Mu, Sigma = Sigma,
                  Sinv.method = Sinv.method, Sigma.inv = Sigma.inv)
    } else {
        SC <- lav_mvnorm_scores_vech_sigma(Y = Y, wt = wt,
                  Mu = Mu, Sigma = Sigma,
                  Sinv.method = Sinv.method, Sigma.inv = Sigma.inv)
    }

    crossprod(SC)/N
}


# 6: inverted information h0

# 6a: inverted unit expected information h0 Mu and vech(Sigma)
lav_mvnorm_inverted_information_expected <- function(Y     = NULL, # unused!
                                                     wt    = NULL, # unused!
                                                     Mu    = NULL, # unused!
                                                     Sigma         = NULL,
                                                     meanstructure = TRUE) {

    I22 <- 2 * lav_matrix_duplication_ginv_pre_post(Sigma %x% Sigma)

    if(meanstructure) {
        I11 <- Sigma
        out <- lav_matrix_bdiag(I11, I22)
    } else {
        out <- I22
    }

    out
}

# 6b:  inverted unit observed information h0

# one could use the inverse of a partitioned matrix, but that does not
# seem to help much... unless we can find an expression for solve(I22)

# 6c: inverted unit first-order information h0
# /


# 7) ACOV h0 mu + vech(Sigma)
# not implemented, as too trivial

#    7a: 1/N * inverted expected    information

#    7b: 1/N * inverted observed    information

#    7c: 1/N * inverted first-order information

#    7d: sandwich acov

