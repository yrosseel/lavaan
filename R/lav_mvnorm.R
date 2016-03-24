# the multivariate normal distribution

# 1) loglikelihood (from raw data, or sample statistics)
# 2) derivatives with respect to mu, Sigma, vech(Sigma)
# 3) casewise scores with respect to mu, vech(Sigma), mu + vech(Sigma)
# 4) 4a: (unit) information of mu + vech(Sigma)
#    4b: (unit) observed information

# YR 07 Feb 2016: first version

# 1) likelihood

# 1a: input is raw data
# (note casewise = TRUE same as: dmvnorm(Y, mean, sigma, log = TRUE))
lav_mvnorm_loglik_data <- function(Y             = NULL,
                                   Mu            = NULL,
                                   Sigma         = NULL,
                                   casewise      = FALSE,
                                   Sinv.method   = "eigen") {
    P <- NCOL(Y); N <- NROW(Y); Mu <- as.numeric(Mu)

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

    } else {
        # invert Sigma
        Sigma.inv <- lav_matrix_symmetric_inverse(S = Sigma, logdet = TRUE,
                                                  Sinv.method = Sinv.method)
        sample.mean <- colMeans(Y)
        sample.cov <- 1/N*crossprod(Y) - tcrossprod(sample.mean)
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



# 2. Derivatives

# 2a: derivative logl with respect to mu
lav_mvnorm_dlogl_dmu <- function(Y           = NULL,
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

    # derivative
    dmu <- as.numeric(Sigma.inv %*% colSums(Yc))

    dmu
}

# 2b: derivative logl with respect to Sigma (full matrix, ignoring symmetry)
lav_mvnorm_dlogl_dSigma <- function(Y           = NULL,
                                    Mu          = NULL,
                                    Sigma       = NULL,
                                    Sinv.method = "eigen",
                                    Sigma.inv   = NULL) {
    N <- NROW(Y); Mu <- as.numeric(Mu)

    if(is.null(Sigma.inv)) {
        # invert Sigma
        Sigma.inv <- lav_matrix_symmetric_inverse(S = Sigma, logdet = FALSE,
                                                  Sinv.method = Sinv.method)
    }

    # substract 'Mu' from Y
    Yc <- t( t(Y) - Mu )

    # W.tilde
    W.tilde <- crossprod(Yc) / N

    # derivative
    dSigma <- -(N/2)* (Sigma.inv - (Sigma.inv %*% W.tilde %*% Sigma.inv))

    dSigma
}

# 2c: derivative logl with respect to vech(Sigma)
lav_mvnorm_dlogl_dvechSigma <- function(Y           = NULL,
                                        Mu          = NULL,
                                        Sigma       = NULL,
                                        Sinv.method = "eigen",
                                        Sigma.inv   = NULL) {
    N <- NROW(Y); Mu <- as.numeric(Mu)

    if(is.null(Sigma.inv)) {
        # invert Sigma
        Sigma.inv <- lav_matrix_symmetric_inverse(S = Sigma, logdet = FALSE,
                                                  Sinv.method = Sinv.method)
    }

    # substract 'Mu' from Y
    Yc <- t( t(Y) - Mu )

    # W.tilde
    W.tilde <- crossprod(Yc) / N

    # derivative (avoiding kronecker product)
    dSigma <- -(N/2)* (Sigma.inv - (Sigma.inv %*% W.tilde %*% Sigma.inv))
    dvechSigma <- as.numeric( lav_matrix_duplication_pre( 
                                  as.matrix(lav_matrix_vec(dSigma)) ) )

    dvechSigma
}

# 3. Casewise scores

# 3a: casewise scores with respect to mu
lav_mvnorm_scores_mu <- function(Y           = NULL,
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

    SC
}

# 3b: casewise scores with respect to vech(Sigma)
lav_mvnorm_scores_vech_sigma <- function(Y           = NULL,
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

    SC
}

# 3c: casewise scores with respect to mu + vech(Sigma)
lav_mvnorm_scores_mu_vech_sigma <- function(Y           = NULL,
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
    
    cbind(Yc, SC)
}

# 4) Information 

# 4a: unit expected information Mu and vech(Sigma)
lav_mvnorm_information_expected <- function(Sigma.inv = NULL) {

    I11 <- Sigma.inv
    I22 <- 0.5 * lav_matrix_duplication_pre_post(Sigma.inv %x% Sigma.inv)

    lav_matrix_bdiag(I11, I22)
}

# 4b: unit observed inforrmation Mu and vech(Sigma) from raw data
lav_mvnorm_information_observed_data <- function(Y           = NULL,
                                                 Mu          = NULL,
                                                 Sigma       = NULL,
                                                 Sinv.method = "eigen",
                                                 Sigma.inv   = NULL) {
    N <- NROW(Y)

    # sample statistics
    sample.mean <- colMeans(Y)
    sample.cov <- cov(Y) * (N-1)/N

    lav_mvnorm_information_observed_samplestats(sample.mean = sample.mean,
        sample.cov = sample.cov, sample.nobs = N, Mu = Mu, Sigma = Sigma,
        Sinv.method = Sinv.method, Sigma.inv = Sigma.inv)
}

# 4b-bis: hessian Mu and vech(Sigma) from samplestats
lav_mvnorm_information_observed_samplestats <-
    function(sample.mean  = NULL,
             sample.cov   = NULL,
             sample.nobs  = NULL,
             Mu           = NULL,
             Sigma        = NULL,
             Sinv.method  = "eigen",
             Sigma.inv    = NULL) {

    sample.mean <- as.numeric(sample.mean); Mu <- as.numeric(Mu)
    N <- sample.nobs

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

