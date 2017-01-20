# the multivariate linear model using maximum likelihood
# - loglikelihood (from raw data, or sample statitics)
# - derivatives with respect to Beta, Sigma, vech(Sigma)
# - casewise scores with respect to Beta, vech(Sigma), Beta + vech(Sigma)
# - (unit) information of Beta + vech(Sigma)
# - hessian of Beta + vech(Sigma)

# 1. input is raw data
lav_mvreg_loglik_data <- function(Y           = NULL,
                                  X           = NULL, # includes intercept
                                  Beta        = NULL,
                                  Sigma       = NULL,
                                  casewise    = FALSE,
                                  Sinv.method = "eigen") {
    Q <- NCOL(Y); N <- NROW(Y)

    if(casewise) {
        LOG.2PI <- log(2 * pi)

        # invert Sigma
        if(Sinv.method == "chol") {
            cS <- chol(Sigma); icS <- backsolve(cS, diag(Q))
            logdet <- -2 * sum(log(diag(icS)))

            RES <- Y - X %*% Beta
            DIST <- rowSums((RES %*% icS)^2)
        } else {
            Sigma.inv <- lav_matrix_symmetric_inverse(S = Sigma, logdet = TRUE,
                                                      Sinv.method = Sinv.method)
            logdet <- attr(Sigma.inv, "logdet")

            RES <- Y - X %*% Beta
            DIST <- rowSums(RES %*% Sigma.inv * RES)
        }

        loglik <- -(Q * LOG.2PI + logdet + DIST)/2

    } else {
        # invert Sigma
        Sigma.inv <- lav_matrix_symmetric_inverse(S = Sigma, logdet = TRUE,
                                                  Sinv.method = Sinv.method)
        logdet <- attr(Sigma.inv, "logdet")

        RES <- Y - X %*% Beta
        # TOTAL <- TR( (Y - X%*%Beta) %*% Sigma.inv %*% t(Y - X%*%Beta) )
        TOTAL <- sum( rowSums(RES %*% Sigma.inv * RES) )
        loglik <- -(N*Q/2)*log(2*pi) - (N/2)*logdet - (1/2)*TOTAL
    }

    loglik
}


# 2. input are sample statistics (beta, cov, N) only
lav_mvreg_loglik_samplestats <- function(sample.res.beta = NULL,
                                         sample.res.cov  = NULL,
                                         sample.XX       = NULL,
                                         sample.nobs     = NULL,
                                         Beta            = NULL,
                                         Sigma           = NULL,
                                         Sinv.method     = "eigen",
                                         Sigma.inv       = NULL) {

    Q <- NCOL(sample.res.cov); N <- sample.nobs
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
    DIST1 <- sum(Sigma.inv * sample.res.cov)

    # tr( Sigma^{-1} (B-beta)' X'X (B-beta) 
    Diff <- sample.res.beta - Beta
    DIST2 <- sum(Sigma.inv * crossprod(Diff, (1/N)*sample.XX) %*% Diff)

    loglik <- -(N/2) * (Q*log(2*pi) + logdet + DIST1 + DIST2)

    loglik
}

# derivative logl with respect to Beta
# version 1: using Y/X
#    lav_matrix_vec( t(X) %*% RES %*% Sigma.inv )
# version 2: using B/S
#    lav_matrix_vec(XX %*% (B - Beta) %*% Sigma.inv)
lav_mvreg_dlogl_dbeta <- function(Y           = NULL,
                                  X           = NULL,
                                  Beta        = NULL,
                                  Sigma       = NULL,
                                  Sinv.method = "eigen",
                                  Sigma.inv   = NULL) {
    if(is.null(Sigma.inv)) {
        # invert Sigma
        Sigma.inv <- lav_matrix_symmetric_inverse(S = Sigma, logdet = FALSE,
                                                  Sinv.method = Sinv.method)
    }

    # substract 'X %*% Beta' from Y
    RES <- Y - X %*% Beta

    # derivative
    dbeta <- as.numeric( t(X) %*% RES %*% Sigma.inv )

    dbeta
}

# derivative logl with respect to Sigma (full matrix, ignoring symmetry)
lav_mvreg_dlogl_dSigma <- function(Y           = NULL,
                                   X           = NULL,
                                   Beta        = NULL,
                                   Sigma       = NULL,
                                   Sinv.method = "eigen",
                                   Sigma.inv   = NULL) {
    N <- NROW(Y)

    if(is.null(Sigma.inv)) {
        # invert Sigma
        Sigma.inv <- lav_matrix_symmetric_inverse(S = Sigma, logdet = FALSE,
                                                  Sinv.method = Sinv.method)
    }

    # substract 'X %*% Beta' from Y
    RES <- Y - X %*% Beta

    # W.tilde
    W.tilde <- crossprod(RES)/N

    # derivative
    dSigma <- -(N/2)* (Sigma.inv - (Sigma.inv %*% W.tilde %*% Sigma.inv))

    dSigma
}

# derivative logl with respect to vech(Sigma)
lav_mvreg_dlogl_dvechSigma <- function(Y           = NULL,
                                       X           = NULL,
                                       Beta        = NULL,
                                       Sigma       = NULL,
                                       Sinv.method = "eigen",
                                       Sigma.inv   = NULL) {
    N <- NROW(Y)

    if(is.null(Sigma.inv)) {
        # invert Sigma
        Sigma.inv <- lav_matrix_symmetric_inverse(S = Sigma, logdet = FALSE,
                                                  Sinv.method = Sinv.method)
    }

    # substract 'X %*% Beta' from Y
    RES <- Y - X %*% Beta

    # W.tilde
    W.tilde <- crossprod(RES)/N

    # derivative
    dSigma <- -(N/2)* (Sigma.inv - (Sigma.inv %*% W.tilde %*% Sigma.inv))
    dvechSigma <- as.numeric( lav_matrix_duplication_pre(
                                  as.matrix(lav_matrix_vec(dSigma)) ) )

    dvechSigma
}


# casewise scores with respect to Beta
lav_mvreg_scores_beta <- function(Y           = NULL,
                                  X           = NULL,
                                  Beta        = NULL,
                                  Sigma       = NULL,
                                  Sinv.method = "eigen",
                                  Sigma.inv   = NULL) {
    Q <- NCOL(Y); P <- NCOL(X)

    if(is.null(Sigma.inv)) {
        # invert Sigma
        Sigma.inv <- lav_matrix_symmetric_inverse(S = Sigma, logdet = FALSE,
                                                  Sinv.method = Sinv.method)
    }

    # substract Mu
    RES <- Y - X %*% Beta

    # post-multiply with Sigma.inv
    RES <- RES %*% Sigma.inv

    SC.Beta <- X[,  rep(1:P, times = Q), drop = FALSE] * 
               RES[,rep(1:Q,  each = P), drop = FALSE]

    SC.Beta
}


# casewise scores with respect to vech(Sigma)
lav_mvreg_scores_vech_sigma <- function(Y           = NULL,
                                        X           = NULL,
                                        Beta        = NULL,
                                        Sigma       = NULL,
                                        Sinv.method = "eigen",
                                        Sigma.inv   = NULL) {
    Q <- NCOL(Y)
    
    if(is.null(Sigma.inv)) {
        # invert Sigma
        Sigma.inv <- lav_matrix_symmetric_inverse(S = Sigma, logdet = FALSE,
                                                  Sinv.method = Sinv.method)
    }                       
        
    # vech(Sigma.inv)
    isigma <- lav_matrix_vech(Sigma.inv)
    
    # substract X %*% Beta
    RES <- Y - X %*% Beta

    # postmultiply with Sigma.inv
    RES <- RES %*% Sigma.inv
    
    # tcrossprod 
    idx1 <- lav_matrix_vech_col_idx(Q)
    idx2 <- lav_matrix_vech_row_idx(Q)
    Z <- RES[,idx1] * RES[,idx2]

    # substract isigma from each row
    SC <- t( t(Z) - isigma )

    # adjust for vech (and avoiding the 1/2 factor)
    SC[,lav_matrix_diagh_idx(Q)] <- SC[,lav_matrix_diagh_idx(Q)] / 2

    SC
}


# casewise scores with respect to beta + vech(Sigma)
lav_mvreg_scores_beta_vech_sigma <- function(Y           = NULL,
                                             X           = NULL,
                                             Beta        = NULL,
                                             Sigma       = NULL,
                                             Sinv.method = "eigen",
                                             Sigma.inv   = NULL) {
    Q <- NCOL(Y); P <- NCOL(X)
    
    if(is.null(Sigma.inv)) {
        # invert Sigma
        Sigma.inv <- lav_matrix_symmetric_inverse(S = Sigma, logdet = FALSE,
                                                  Sinv.method = Sinv.method)
    }                       
        
    # vech(Sigma.inv)
    isigma <- lav_matrix_vech(Sigma.inv)
    
    # substract X %*% Beta
    RES <- Y - X %*% Beta

    # postmultiply with Sigma.inv
    RES <- RES %*% Sigma.inv

    SC.Beta <- X[,  rep(1:P, times = Q), drop = FALSE] *
               RES[,rep(1:Q,  each = P), drop = FALSE]
    
    # tcrossprod 
    idx1 <- lav_matrix_vech_col_idx(Q)
    idx2 <- lav_matrix_vech_row_idx(Q)
    Z <- RES[,idx1] * RES[,idx2]

    # substract isigma from each row
    SC <- t( t(Z) - isigma )

    # adjust for vech (and avoiding the 1/2 factor)
    SC[,lav_matrix_diagh_idx(Q)] <- SC[,lav_matrix_diagh_idx(Q)] / 2

    cbind(SC.Beta, SC)
}

# information Beta and vech(Sigma)
lav_mvreg_information_beta_vech_sigma_samplestats <- 
    function(Sigma.inv = NULL,
             sample.XX = NULL,
             sample.nobs = NULL) {

    XXN <- (1/sample.nobs) * sample.XX

    I11 <- Sigma.inv %x% XXN
    I22 <- 0.5 * lav_matrix_duplication_pre_post(Sigma.inv %x% Sigma.inv)

    lav_matrix_bdiag(I11, I22)
}


# hessian Beta and vech(Sigma)
lav_mvreg_hessian_beta_vech_sigma <- function(Y           = NULL,
                                              X           = NULL,
                                              Beta        = NULL,
                                              Sigma       = NULL,
                                              Sinv.method = "eigen",
                                              Sigma.inv   = NULL) {

    N <- NROW(Y)

    if(is.null(Sigma.inv)) {
        # invert Sigma
        Sigma.inv <- lav_matrix_symmetric_inverse(S = Sigma, logdet = FALSE,
                                                  Sinv.method = Sinv.method)
    }

    RES <- Y - X %*% Beta
    W.tilde <- 1/N * crossprod(RES)

    H11 <- Sigma.inv %x% ((1/N) * crossprod(X))
    H21 <- lav_matrix_duplication_pre( Sigma.inv %x% 
       (Sigma.inv %*% ((1/N) * crossprod(RES, X))) )
    H12 <- t(H21)

    AAA <- Sigma.inv %*% (2*W.tilde - Sigma) %*% Sigma.inv
    H22 <- (1/2) * lav_matrix_duplication_pre_post(Sigma.inv %x% AAA)

    H <- -N * rbind( cbind(H11, H12),
                     cbind(H21, H22) )

    H
}



