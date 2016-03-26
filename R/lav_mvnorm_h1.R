# the multivariate normal distribution, unrestricted (h1)
# - everything is evalued under the MLEs: Mu = ybar, Sigma = S 

# 1) loglikelihood h1 (from raw data, or sample statistics)
# 4) hessian h1 around MLEs
# 5) information h1 (restricted Sigma/mu)
#    5a: (unit)    expected information h1 (A1)
#    5b: (unit)    observed information h1 (A1)
#    5c: (unit) first.order information h1 (B1 = A1 %*% Gamma %*% A1)


# YR 25 March 2016: first version

# 1. likelihood h1

# 1a: input is raw data
lav_mvnorm_h1_loglik_data <- function(Y             = NULL,
                                      casewise      = FALSE,
                                      Sinv.method   = "eigen") {
    P <- NCOL(Y); N <- NROW(Y)

    # sample statistics
    sample.mean <- colMeans(Y)
    sample.cov <- 1/N*crossprod(Y) - tcrossprod(sample.mean)

    if(casewise) {
        LOG.2PI <- log(2 * pi)

        # invert sample.cov
        if(Sinv.method == "chol") {
            cS <- chol(sample.cov); icS <- backsolve(cS, diag(P))
            Yc <- t( t(Y) - sample.mean )
            DIST <- rowSums((Yc %*% icS)^2)
            logdet <- -2 * sum(log(diag(icS)))
        } else {
            sample.cov.inv <- lav_matrix_symmetric_inverse(S = sample.cov, 
                                  logdet = TRUE, Sinv.method = Sinv.method)
            logdet <- attr(sample.cov.inv, "logdet")
            # mahalanobis distance
            Yc <- t( t(Y) - sample.mean )
            DIST <- rowSums(Yc %*% sample.cov.inv * Yc)
        }

        loglik <- -(P * LOG.2PI + logdet + DIST)/2

    } else {
        # invert sample.cov
        sample.cov.inv <- lav_matrix_symmetric_inverse(S = sample.cov, 
                              logdet = TRUE, Sinv.method = Sinv.method)
        loglik <- 
            lav_mvnorm_h1_loglik_samplestats(sample.mean    = sample.mean,
                                             sample.cov     = sample.cov,
                                             sample.nobs    = N,
                                             sample.cov.inv = sample.cov.inv)
    }

    loglik
}



# 1b: input are sample statistics (mean, cov, N) only
lav_mvnorm_h1_loglik_samplestats <- function(sample.mean    = NULL,
                                             sample.cov     = NULL,
                                             sample.nobs    = NULL,
                                             Sinv.method    = "eigen",
                                             sample.cov.inv = NULL) {

    P <- length(sample.mean); N <- sample.nobs
    sample.mean <- as.numeric(sample.mean)
    LOG.2PI <- log(2 * pi)

    if(is.null(sample.cov.inv)) {
        sample.cov.inv <- lav_matrix_symmetric_inverse(S = sample.cov, 
                              logdet = TRUE, Sinv.method = Sinv.method)
        logdet <- attr(sample.cov.inv, "logdet")
    } else {
        logdet <- attr(sample.cov.inv, "logdet")
        if(is.null(logdet)) {
            # compute - ln|S.inv|
            ev <- eigen(sample.cov.inv, symmetric = TRUE, only.values = TRUE)
            logdet <- -1 * sum(log(ev$values))
        }
    }

    loglik <- -N/2 * (P * LOG.2PI + logdet + P)

    loglik
}


# 4. hessian of logl (around MLEs of Mu and Sigma)

# 4a: hessian logl Mu and vech(Sigma) from raw data
lav_mvnorm_h1_logl_hessian_data <- function(Y              = NULL,
                                            Sinv.method    = "eigen",
                                            sample.cov.inv = NULL) {
    N <- NROW(Y)

    # observed information
    observed <- lav_mvnorm_h1_information_observed_data(Y = Y, 
                    Sinv.method = Sinv.method, sample.cov.inv = sample.cov.inv)

    -N*observed
}

# 4b: hessian Mu and vech(Sigma) from samplestats
lav_mvnorm_h1_logl_hessian_samplestats <-
    function(sample.mean    = NULL, # unused!
             sample.cov     = NULL,
             sample.nobs    = NULL,
             Sinv.method    = "eigen",
             sample.cov.inv = NULL) {

    N <- sample.nobs

    # observed information
    observed <- lav_mvnorm_h1_information_observed_samplestats(sample.mean = 
        sample.mean, sample.cov = sample.cov, Sinv.method = Sinv.method, 
        sample.cov.inv = sample.cov.inv)
    
    -N*observed
}



# 5) Information h1 (not expected == observed if data is complete!)

# 5a: unit expected information h1
lav_mvnorm_h1_information_expected <- function(Y              = NULL,
                                               sample.cov     = NULL,
                                               Sinv.method    = "eigen",
                                               sample.cov.inv = NULL) {
    if(is.null(sample.cov.inv)) {

        if(is.null(sample.cov)) {
            # sample statistics
            sample.mean <- colMeans(Y); N <- NCOL(Y)
            sample.cov <- 1/N*crossprod(Y) - tcrossprod(sample.mean)
        }

        # invert sample.cov
        sample.cov.inv <- lav_matrix_symmetric_inverse(S = sample.cov, 
                              logdet = FALSE, Sinv.method = Sinv.method)
    }

    I11 <- sample.cov.inv
    I22 <- 0.5 * lav_matrix_duplication_pre_post(sample.cov.inv %x% 
                                                 sample.cov.inv)

    lav_matrix_bdiag(I11, I22)
}

# 5b: unit observed information h1
lav_mvnorm_h1_information_observed_data <- function(Y              = NULL,
                                                    Sinv.method    = "eigen",
                                                    sample.cov.inv = NULL) {

    lav_mvnorm_h1_information_expected(Y = Y, Sinv.method = Sinv.method,
                                       sample.cov.inv = sample.cov.inv)
}

# 5b-bis: observed information h1 from sample statistics
lav_mvnorm_h1_information_observed_samplestats <-
    function(sample.mean    = NULL, # unused!
             sample.cov     = NULL, 
             Sinv.method    = "eigen",
             sample.cov.inv = NULL) {

    if(is.null(sample.cov.inv)) {
        # invert sample.cov
        sample.cov.inv <- lav_matrix_symmetric_inverse(S = sample.cov, 
                              logdet = FALSE, Sinv.method = Sinv.method)
    }

    I11 <- sample.cov.inv
    I22 <- 0.5 * lav_matrix_duplication_pre_post(sample.cov.inv %x% 
                                                 sample.cov.inv)

    lav_matrix_bdiag(I11, I22)
}

# 5c: unit first-order information h1
#     note: first order information h1 == A1 %*% Gamma %*% A1 
#           (where A1 = obs/exp information h1)
lav_mvnorm_h1_information_firstorder <- function(Y              = NULL,
                                                 sample.cov     = NULL,
                                                 Sinv.method    = "eigen",
                                                 sample.cov.inv = NULL,
                                                 Gamma          = NULL) {
    # Gamma
    if(is.null(Gamma)) {
        Gamma <- lav_samplestats_Gamma(Y, meanstructure = TRUE)
    }

    # sample.cov.in
    if(is.null(sample.cov.inv)) {
        # invert sample.cov
        sample.cov.inv <- lav_matrix_symmetric_inverse(S = sample.cov,
                              logdet = FALSE, Sinv.method = Sinv.method)
    }

    # A1
    A1 <- lav_mvnorm_h1_information_expected(Y = Y, Sinv.method = Sinv.method,
                                             sample.cov.inv = sample.cov.inv)

    A1 %*% Gamma %*% A1
}


