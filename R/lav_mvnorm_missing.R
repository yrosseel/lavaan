# the multivariate normal distribution + missing values
#  (so-called 'FIML')

# 1) loglikelihood (from raw data, or sample statitics)
# 2) derivatives with respect to mu, Sigma, vech(Sigma)
# 3) casewise scores with respect to mu, vech(Sigma), mu + vech(Sigma)
# 4) hessian of mu + vech(Sigma)
# 5) (unit) information of mu + vech(Sigma)


# YR 09 Feb 2016: first version


# 1) likelihood

# 1a: input is raw data
#  - two strategies: 1) using missing patterns (pattern = TRUE)
#                    2) truly case per case    (pattern = FALSE)
#    depending on the sample size, missing patterns, etc... one can be 
#    (much) faster than the other
lav_mvnorm_missing_loglik_data <- function(Y           = NULL,
                                           Mu          = NULL,
                                           Sigma       = NULL,
                                           casewise    = FALSE,
                                           pattern     = TRUE,
                                           Sinv.method = "eigen") {

    if(pattern) {
        llik <- lav_mvnorm_missing_llik_pattern(Y = Y, Mu = Mu, 
                    Sigma = Sigma, Sinv.method = Sinv.method)
    } else {
        llik <- lav_mvnorm_missing_llik_casewise(Y = Y, Mu = Mu, 
                    Sigma = Sigma, Sinv.method = Sinv.method)
    }

    if(casewise) {
        loglik <- llik
    } else {
        loglik <- sum(llik, na.rm = TRUE)
    }
    
    loglik
}

# 1b: input are sample statistics (mean, cov, N) per pattern
lav_mvnorm_missing_loglik_samplestats <- function(Yp          = NULL,
                                                  Mu          = NULL,
                                                  Sigma       = NULL,
                                                  Sinv.method = "eigen") {

    LOG.2PI <- log(2*pi); pat.N <- length(Yp);  P <- length(Yp[[1]]$var.idx)

    # global inverse + logdet
    Sigma.inv <- lav_matrix_symmetric_inverse(S = Sigma, logdet = TRUE,
                                              Sinv.method = Sinv.method)
    Sigma.logdet <- attr(Sigma.inv, "logdet")

    # DIST/logdet per pattern
    DIST <- logdet <- P.LOG.2PI <- numeric(pat.N)

    # for each pattern, compute sigma.inv/logdet; compute DIST for all
    # observations of this pattern
    for(p in seq_len(pat.N)) {

        # observed variables for this pattern
        var.idx <- Yp[[p]]$var.idx

        # missing values for this pattern
        na.idx <- which(!var.idx)

        # constant
        P.LOG.2PI[p] <- sum(var.idx) * LOG.2PI * Yp[[p]]$freq

        # invert Sigma for this pattern
        if(length(na.idx) > 0L) {
            sigma.inv <- lav_matrix_symmetric_inverse_update(S.inv = Sigma.inv,
                rm.idx = na.idx, logdet = TRUE, S.logdet = Sigma.logdet)
            logdet[p] <- attr(sigma.inv, "logdet") * Yp[[p]]$freq
        } else {
            sigma.inv <- Sigma.inv
            logdet[p] <- Sigma.logdet * Yp[[p]]$freq
        }

        TT <- Yp[[p]]$SY + tcrossprod(Yp[[p]]$MY - Mu[var.idx])
        DIST[p] <- sum(sigma.inv * TT) * Yp[[p]]$freq
    }

    # compute loglikelihoods per pattern
    loglik <- sum(-(P.LOG.2PI + logdet + DIST)/2)

    loglik
}

## casewise loglikelihoods

# casewise Sinv.method
lav_mvnorm_missing_llik_casewise <- function(Y           = NULL,
                                             Mu          = NULL,
                                             Sigma       = NULL,
                                             Sinv.method = "eigen") {

    P <- NCOL(Y); N <- NROW(Y); LOG.2PI <- log(2*pi); Mu <- as.numeric(Mu)

    # global inverse + logdet
    Sigma.inv <- lav_matrix_symmetric_inverse(S = Sigma, logdet = TRUE,
                                              Sinv.method = Sinv.method)
    Sigma.logdet <- attr(Sigma.inv, "logdet")

    # subtract Mu
    Yc <- t( t(Y) - Mu )

    # DIST/logdet per case
    DIST <- logdet <- P.LOG.2PI <- rep(as.numeric(NA), N)

    # missing pattern per case
    OBS <- !is.na(Y); P.i <- rowSums(OBS)

    # constant
    P.LOG.2PI <- P.i * LOG.2PI

    # complete cases first (only an advantage if we have mostly complete 
    # observations)
    other.idx <- seq_len(N)
    complete.idx <- which(P.i == P)
    if(length(complete.idx) > 0L) {
        other.idx <- other.idx[-complete.idx]
        DIST[complete.idx] <-
            rowSums(Yc[complete.idx,,drop = FALSE] %*% Sigma.inv *
                    Yc[complete.idx,,drop = FALSE])
        logdet[complete.idx] <- Sigma.logdet
    }

    # non-complete cases
    for(i in other.idx) {
        na.idx <- which(!OBS[i,])

        # catch empty cases
        if(length(na.idx) == P) next

        # invert Sigma for this pattern
        sigma.inv <- lav_matrix_symmetric_inverse_update(S.inv = Sigma.inv,
                rm.idx = na.idx, logdet = TRUE, S.logdet = Sigma.logdet)
        logdet[i] <- attr(sigma.inv, "logdet")

        # distance for this case
        DIST[i] <- sum(sigma.inv * crossprod(Yc[i, OBS[i,], drop = FALSE]))
    }

    # compute casewise loglikelihoods
    llik <- -(P.LOG.2PI + logdet + DIST)/2

    llik
}

# pattern-based, but casewise loglikelihoods
lav_mvnorm_missing_llik_pattern <- function(Y           = NULL,
                                            Mp          = NULL,
                                            Mu          = NULL,
                                            Sigma       = NULL,
                                            Sinv.method = "eigen") {

    P <- NCOL(Y); N <- NROW(Y); LOG.2PI <- log(2*pi); Mu <- as.numeric(Mu)

    # global inverse + logdet
    Sigma.inv <- lav_matrix_symmetric_inverse(S = Sigma, logdet = TRUE,
                                              Sinv.method = Sinv.method)
    Sigma.logdet <- attr(Sigma.inv, "logdet")

    # subtract Mu
    Yc <- t( t(Y) - Mu )

    # DIST/logdet per case
    DIST <- logdet <- P.LOG.2PI <- rep(as.numeric(NA), N)

    # missing patterns
    if(is.null(Mp)) {
        Mp <- lav_data_missing_patterns(Y)
    }

    # for each pattern, compute sigma.inv/logdet; compute DIST for all
    # observations of this pattern
    for(p in seq_len(Mp$npatterns)) {

        # observed values for this pattern
        var.idx <- Mp$pat[p,]

        # missing values for this pattern
        na.idx <- which(!var.idx)

        # identify cases with this pattern
        case.idx <- Mp$case.idx[[p]]

        # constant
        P.LOG.2PI[case.idx] <- sum(var.idx) * LOG.2PI

        # invert Sigma for this pattern
        if(length(na.idx) > 0L) {
            sigma.inv <- lav_matrix_symmetric_inverse_update(S.inv = Sigma.inv,
                rm.idx = na.idx, logdet = TRUE, S.logdet = Sigma.logdet)
            logdet[case.idx] <- attr(sigma.inv, "logdet")
        } else {
            sigma.inv <- Sigma.inv
            logdet[case.idx] <- Sigma.logdet
        }

        if(Mp$freq[p] == 1L) {
            DIST[case.idx] <- sum(sigma.inv * 
                crossprod(Yc[case.idx, var.idx, drop = FALSE]))
        } else {
            DIST[case.idx] <-
                rowSums(Yc[case.idx, var.idx, drop = FALSE] %*% sigma.inv *
                        Yc[case.idx, var.idx, drop = FALSE])
        }  
    }

    # compute casewise loglikelihoods
    llik <- -(P.LOG.2PI + logdet + DIST)/2

    llik
}




# 2. Derivatives

# 2a: derivative logl with respect to mu
lav_mvnorm_missing_dlogl_dmu <- function(Y           = NULL,
                                         Mu          = NULL,
                                         Sigma       = NULL,
                                         Sigma.inv   = NULL,
                                         Sinv.method = "eigen") {

    SC <- lav_mvnorm_missing_scores_mu(Y = Y, Mu = Mu, Sigma = Sigma,
              Sigma.inv = Sigma.inv, Sinv.method = Sinv.method)

    colSums(SC, na.rm = TRUE)
}

# 2abis: using samplestats 
lav_mvnorm_missing_dlogl_dmu_samplestats <- function(Yp          = NULL,
                                                     Mu          = NULL,
                                                     Sigma       = NULL,
                                                     Sigma.inv   = NULL,
                                                     Sinv.method = "eigen") {
    pat.N <- length(Yp);  P <- length(Yp[[1]]$var.idx)

    if(is.null(Sigma.inv)) {
        Sigma.inv <- lav_matrix_symmetric_inverse(S = Sigma, logdet = FALSE,
                                                  Sinv.method = Sinv.method)
    }

    # dmu
    dmu <- numeric(P)

    # for each pattern, compute sigma.inv
    for(p in seq_len(pat.N)) {

        # observed variables for this pattern
        var.idx <- Yp[[p]]$var.idx

        # missing values for this pattern
        na.idx <- which(!var.idx)

        # invert Sigma for this pattern
        if(length(na.idx) > 0L) {
            sigma.inv <- lav_matrix_symmetric_inverse_update(S.inv = Sigma.inv,
                             rm.idx = na.idx, logdet = FALSE)
        } else {
            sigma.inv <- Sigma.inv
        }

        # dmu for this pattern
        dmu.pattern <- as.numeric(sigma.inv %*% (Yp[[p]]$MY - Mu[var.idx]))

        # update mu
        dmu[var.idx] <- dmu[var.idx] + (dmu.pattern * Yp[[p]]$freq)
    }

    dmu
}



# 2b: derivative logl with respect to Sigma (full matrix, ignoring symmetry)
lav_mvnorm_missing_dlogl_dSigma <- function(Y           = NULL,
                                            Mp          = NULL,
                                            Mu          = NULL,
                                            Sigma       = NULL,
                                            Sigma.inv   = NULL,
                                            Sinv.method = "eigen") {
    P <- NCOL(Y); N <- NROW(Y); Mu <- as.numeric(Mu)

    if(is.null(Sigma.inv)) {
        # invert Sigma
        Sigma.inv <- lav_matrix_symmetric_inverse(S = Sigma, logdet = FALSE,
                                                  Sinv.method = Sinv.method)
    }

    # subtract Mu
    Yc <- t( t(Y) - Mu )

    # dvechSigma
    dSigma <- matrix(0, P, P)

    # missing patterns
    if(is.null(Mp)) {
        Mp <- lav_data_missing_patterns(Y)
    }

    # for each pattern
    for(p in seq_len(Mp$npatterns)) {

        # observed values for this pattern
        var.idx <- Mp$pat[p,]

        # missing values for this pattern
        na.idx <- which(!var.idx)

        # cases with this pattern
        case.idx <- Mp$case.idx[[p]]

        # invert Sigma for this pattern
        if(length(na.idx) > 0L) {
            sigma.inv <- lav_matrix_symmetric_inverse_update(S.inv = Sigma.inv,
                             rm.idx = na.idx, logdet = FALSE)
        } else {
            sigma.inv <- Sigma.inv
        }

        if(length(case.idx) > 1L) {
            W.tilde <- crossprod(Yc[case.idx, var.idx, drop = FALSE])/Mp$freq[p]
        } else {
            W.tilde <- tcrossprod(Yc[case.idx, var.idx])
        }

        # dSigma for this pattern
        dSigma.pattern <- matrix(0, P, P)
        dSigma.pattern[var.idx, var.idx] <- -(1/2) * (sigma.inv - 
                        (sigma.inv %*% W.tilde %*% sigma.inv))

        # update dSigma
        dSigma <- dSigma + (dSigma.pattern * Mp$freq[p])
    }

    dSigma
}

# 2bbis: using samplestats
lav_mvnorm_missing_dlogl_dSigma_samplestats <- function(Yp          = NULL,
                                                        Mu          = NULL,
                                                        Sigma       = NULL,
                                                        Sigma.inv   = NULL,
                                                        Sinv.method = "eigen") {
    pat.N <- length(Yp);  P <- length(Yp[[1]]$var.idx)

    if(is.null(Sigma.inv)) {
        # invert Sigma
        Sigma.inv <- lav_matrix_symmetric_inverse(S = Sigma, logdet = FALSE,
                                                  Sinv.method = Sinv.method)
    }

    # dvechSigma
    dSigma <- matrix(0, P, P)

    # for each pattern
    for(p in seq_len(pat.N)) {

        # observed variables for this pattern
        var.idx <- Yp[[p]]$var.idx

        # missing values for this pattern
        na.idx <- which(!var.idx)

        # invert Sigma for this pattern
        if(length(na.idx) > 0L) {
            sigma.inv <- lav_matrix_symmetric_inverse_update(S.inv = Sigma.inv,
                             rm.idx = na.idx, logdet = FALSE)
        } else {
            sigma.inv <- Sigma.inv
        }

        W.tilde <- Yp[[p]]$SY + tcrossprod(Yp[[p]]$MY - Mu[var.idx])

        # dSigma for this pattern
        dSigma.pattern <- matrix(0, P, P)
        dSigma.pattern[var.idx, var.idx] <- -(1/2) * (sigma.inv - 
                        (sigma.inv %*% W.tilde %*% sigma.inv))

        # update dSigma
        dSigma <- dSigma + (dSigma.pattern * Yp[[p]]$freq)
    }

    dSigma
}


# 2c: derivative logl with respect to vech(Sigma)
lav_mvnorm_missing_dlogl_dvechSigma <- function(Y           = NULL,
                                                Mu          = NULL,
                                                Sigma       = NULL,
                                                Sigma.inv   = NULL,
                                                Sinv.method = "eigen") {

    FULL <- lav_mvnorm_missing_dlogl_dSigma(Y = Y, Mu = Mu, Sigma = Sigma,
              Sigma.inv = Sigma.inv, Sinv.method = Sinv.method)
    as.numeric( lav_matrix_duplication_pre( as.matrix(lav_matrix_vec(FULL)) ) )
}

# 2cbis: using samplestats
lav_mvnorm_missing_dlogl_dvechSigma_samplestats <- 
    function(Yp          = NULL,
             Mu          = NULL,
             Sigma       = NULL,
             Sigma.inv   = NULL,
             Sinv.method = "eigen") {
    pat.N <- length(Yp);  P <- length(Yp[[1]]$var.idx)

    if(is.null(Sigma.inv)) {
        # invert Sigma
        Sigma.inv <- lav_matrix_symmetric_inverse(S = Sigma, logdet = FALSE,
                                                  Sinv.method = Sinv.method)
    }

    # dvechSigma
    dvechSigma <- numeric(P*(P+1)/2)

    # for each pattern
    for(p in seq_len(pat.N)) {

        # observed variables for this pattern
        var.idx <- Yp[[p]]$var.idx

        # missing values for this pattern
        na.idx <- which(!var.idx)

        # invert Sigma for this pattern
        if(length(na.idx) > 0L) {
            sigma.inv <- lav_matrix_symmetric_inverse_update(S.inv = Sigma.inv,
                             rm.idx = na.idx, logdet = FALSE)
        } else {
            sigma.inv <- Sigma.inv
        }

        W.tilde <- Yp[[p]]$SY + tcrossprod(Yp[[p]]$MY - Mu[var.idx])

        # dSigma for this pattern
        dSigma.pattern <- matrix(0, P, P)
        dSigma.pattern[var.idx, var.idx] <- -(1/2) * (sigma.inv - 
                        (sigma.inv %*% W.tilde %*% sigma.inv))

        # convert to vechSigma
        dvechSigma.pattern <- as.numeric( lav_matrix_duplication_pre(
                                  as.matrix(lav_matrix_vec(dSigma.pattern)) ) )

        # update dvechSigma
        dvechSigma <- dvechSigma + (dvechSigma.pattern * Yp[[p]]$freq)
    }

    dvechSigma
}




# 3. Casewise scores

# 3a: casewise scores with respect to mu
lav_mvnorm_missing_scores_mu <- function(Y           = NULL,
                                         Mp          = NULL,
                                         Mu          = NULL,
                                         Sigma       = NULL,
                                         Sigma.inv   = NULL,
                                         Sinv.method = "eigen") {

    P <- NCOL(Y); N <- NROW(Y); Mu <- as.numeric(Mu)

    if(is.null(Sigma.inv)) {
        # invert Sigma
        Sigma.inv <- lav_matrix_symmetric_inverse(S = Sigma, logdet = FALSE,
                                                  Sinv.method = Sinv.method)
    }

    # missing patterns
    if(is.null(Mp)) {
        Mp <- lav_data_missing_patterns(Y)
    }

    # subtract Mu
    Yc <- t( t(Y) - Mu )

    # dmu per case
    dmu <- matrix(as.numeric(NA), N, P)

    # for each pattern, compute sigma.inv
    for(p in seq_len(Mp$npatterns)) {

        # observed values for this pattern
        var.idx <- Mp$pat[p,]

        # missing values for this pattern
        na.idx <- which(!var.idx)

        case.idx <- Mp$case.idx[[p]]

        # invert Sigma for this pattern
        if(length(na.idx) > 0L) {
            sigma.inv <- lav_matrix_symmetric_inverse_update(S.inv = Sigma.inv,
                             rm.idx = na.idx, logdet = FALSE)
        } else {
            sigma.inv <- Sigma.inv
        }

        # compute dMu for all observations of this pattern
        dmu[case.idx, var.idx] <-
            Yc[case.idx, var.idx, drop = FALSE] %*% sigma.inv
    }

    dmu
}

# 3b: casewise scores with respect to vech(Sigma)
lav_mvnorm_missing_scores_vech_sigma <- function(Y           = NULL,
                                                 Mp          = NULL,
                                                 Mu          = NULL,
                                                 Sigma       = NULL,
                                                 Sigma.inv   = NULL,
                                                 Sinv.method = "eigen") {

    P <- NCOL(Y); N <- NROW(Y); Mu <- as.numeric(Mu)

    if(is.null(Sigma.inv)) {
        # invert Sigma
        Sigma.inv <- lav_matrix_symmetric_inverse(S = Sigma, logdet = FALSE,
                                                  Sinv.method = Sinv.method)
    }
    
    # for the tcrossprod
    idx1 <- lav_matrix_vech_col_idx(P); idx2 <- lav_matrix_vech_row_idx(P)

    # vech(Sigma.inv)
    iSigma <- lav_matrix_vech(Sigma.inv)

    # missing patterns
    if(is.null(Mp)) {
         Mp <- lav_data_missing_patterns(Y)
    }

    # subtract Mu
    Yc <- t( t(Y) - Mu )

    # SC
    SC <- matrix(as.numeric(NA), nrow = N, ncol = length(iSigma))

    # for each pattern
    for(p in seq_len(Mp$npatterns)) {

        # observed values for this pattern
        var.idx <- Mp$pat[p,]

        # missing values for this pattern
        na.idx <- which(!var.idx)

        # cases with this pattern
        case.idx <- Mp$case.idx[[p]]

        # invert Sigma for this pattern
        if(length(na.idx) > 0L) {
            sigma.inv <- lav_matrix_symmetric_inverse_update(S.inv = Sigma.inv,
                             rm.idx = na.idx, logdet = FALSE)
            tmp <- matrix(0, P, P)
            tmp[var.idx, var.idx] <- sigma.inv
            isigma <- lav_matrix_vech(tmp)
        } else {
            sigma.inv <- Sigma.inv
            isigma <- iSigma
        }

        # postmultiply these cases with sigma.inv
        Yc[case.idx, var.idx] <- Yc[case.idx, var.idx] %*% sigma.inv

        # tcrossprod
        SC[case.idx,] <- Yc[case.idx, idx1] * Yc[case.idx, idx2]

        # substract isigma from each row
        SC[case.idx,] <- t( t(SC[case.idx,,drop = FALSE]) - isigma )
    }

    # adjust for vech
    SC[,lav_matrix_diagh_idx(P)] <- SC[,lav_matrix_diagh_idx(P)] / 2

    SC
}

# 3c: casewise scores with respect to mu + vech(Sigma)
lav_mvnorm_missing_scores_mu_vech_sigma <- function(Y           = NULL,
                                                    Mp          = NULL,
                                                    Mu          = NULL,
                                                    Sigma       = NULL,
                                                    Sigma.inv   = NULL,
                                                    Sinv.method = "eigen") {

    P <- NCOL(Y); N <- NROW(Y); Mu <- as.numeric(Mu)

    if(is.null(Sigma.inv)) {
        # invert Sigma
        Sigma.inv <- lav_matrix_symmetric_inverse(S = Sigma, logdet = FALSE,
                                                  Sinv.method = Sinv.method)
    }
    
    # for the tcrossprod
    idx1 <- lav_matrix_vech_col_idx(P); idx2 <- lav_matrix_vech_row_idx(P)

    # vech(Sigma.inv)
    iSigma <- lav_matrix_vech(Sigma.inv)

    # missing patterns
    if(is.null(Mp)) {
        Mp <- lav_data_missing_patterns(Y)
    }

    # subtract Mu
    Yc <- t( t(Y) - Mu )

    # dmu per case
    dmu <- matrix(as.numeric(NA), N, P)

    # SC
    SC <- matrix(as.numeric(NA), nrow = N, ncol = length(iSigma))

    # for each pattern, compute Yc %*% sigma.inv
    for(p in seq_len(Mp$npatterns)) {

        # observed values for this pattern
        var.idx <- Mp$pat[p,]

        # missing values for this pattern
        na.idx <- which(!var.idx)

        # cases with this pattern
        case.idx <- Mp$case.idx[[p]]

        # invert Sigma for this pattern
        if(length(na.idx) > 0L) {
            sigma.inv <- lav_matrix_symmetric_inverse_update(S.inv = Sigma.inv,
                             rm.idx = na.idx, logdet = FALSE)
            tmp <- matrix(0, P, P)
            tmp[var.idx, var.idx] <- sigma.inv
            isigma <- lav_matrix_vech(tmp)
        } else {
            sigma.inv <- Sigma.inv
            isigma <- iSigma
        }

        # compute dMu for all observations of this pattern
        dmu[case.idx, var.idx] <-
            Yc[case.idx, var.idx, drop = FALSE] %*% sigma.inv

        # postmultiply these cases with sigma.inv
        Yc[case.idx, var.idx] <- Yc[case.idx, var.idx] %*% sigma.inv

        # tcrossprod
        SC[case.idx,] <- Yc[case.idx, idx1] * Yc[case.idx, idx2]

        # substract isigma from each row
        SC[case.idx,] <- t( t(SC[case.idx,,drop = FALSE]) - isigma )
    }

    # adjust for vech
    SC[,lav_matrix_diagh_idx(P)] <- SC[,lav_matrix_diagh_idx(P)] / 2

    cbind(dmu, SC)
}


# 4) Hessian of logl
lav_mvnorm_missing_logl_hessian_data <- function(Y           = NULL,
                                                 Mp          = NULL,
                                                 Mu          = NULL,
                                                 Sigma       = NULL,
                                                 Sinv.method = "eigen",
                                                 Sigma.inv   = NULL) {
    # missing patterns
    Yp <- lav_samplestats_missing_patterns(Y = Y, Mp = Mp)

    lav_mvnorm_missing_logl_hessian_samplestats(Yp = Yp, Mu = Mu,
        Sigma = Sigma, Sinv.method = Sinv.method, Sigma.inv = Sigma.inv)
}

lav_mvnorm_missing_logl_hessian_samplestats <-
    function(Yp          = NULL,
             Mu          = NULL,
             Sigma       = NULL,
             Sinv.method = "eigen",
             Sigma.inv   = NULL) {

    pat.N <- length(Yp);  P <- length(Yp[[1]]$var.idx)

    if(is.null(Sigma.inv)) {
        Sigma.inv <- lav_matrix_symmetric_inverse(S = Sigma, logdet = FALSE,
                                                  Sinv.method = Sinv.method)
    }

    H11 <- matrix(0, P, P)
    H21 <- matrix(0, P*(P+1)/2, P)
    H22 <- matrix(0, P*(P+1)/2, P*(P+1)/2)

    # for each pattern, compute sigma.inv
    for(p in seq_len(pat.N)) {

        # observed variables
        var.idx  <- Yp[[p]]$var.idx
        pat.freq <- Yp[[p]]$freq

        # missing values for this pattern
        na.idx <- which(!var.idx)

        # invert Sigma for this pattern
        if(length(na.idx) > 0L) {
            sigma.inv <- lav_matrix_symmetric_inverse_update(S.inv = Sigma.inv,
                             rm.idx = na.idx, logdet = FALSE)
        } else {
            sigma.inv <- Sigma.inv
        }

        S.inv <- matrix(0, P, P)
        S.inv[var.idx, var.idx] <- sigma.inv

        tmp21 <- matrix(0,P,1)
        tmp21[var.idx,1] <- sigma.inv %*% (Yp[[p]]$MY - Mu[var.idx])

        W.tilde <- Yp[[p]]$SY + tcrossprod(Yp[[p]]$MY - Mu[var.idx])
        AAA <- ( sigma.inv %*% 
                 (2*W.tilde - Sigma[var.idx,var.idx,drop = FALSE]) %*% 
                 sigma.inv )
        tmp22 <- matrix(0, P, P)
        tmp22[var.idx, var.idx] <- AAA

        i11 <- S.inv
        i21 <- lav_matrix_duplication_pre( tmp21 %x% S.inv )
        i22 <- (1/2) * lav_matrix_duplication_pre_post(S.inv %x% tmp22)

        H11 <- H11 + pat.freq * i11
        H21 <- H21 + pat.freq * i21
        H22 <- H22 + pat.freq * i22
    }

    H12 <- t(H21)

    -1 * rbind( cbind(H11, H12),
                cbind(H21, H22) )
}




# 5) Information 

# 5a: expected unit information Mu and vech(Sigma)
#     (only useful under MCAR)
# (old term: Abeta, expected)
lav_mvnorm_missing_information_expected <- function(Y           = NULL,
                                                    Mp          = NULL,
                                                    Mu          = NULL,# unused
                                                    Sigma       = NULL,
                                                    Sigma.inv   = NULL,
                                                    Sinv.method = "eigen") {

    P <- NCOL(Y)

    if(is.null(Sigma.inv)) {
        Sigma.inv <- lav_matrix_symmetric_inverse(S = Sigma, logdet = FALSE,
                                                  Sinv.method = Sinv.method)
    }

    # missing patterns
    if(is.null(Mp)) {
        Mp <- lav_data_missing_patterns(Y)
    }

    # N
    N <- sum(Mp$freq) # removed empty cases!

    I11 <- matrix(0, P, P)
    I22 <- matrix(0, P*(P+1)/2, P*(P+1)/2)

    # for each pattern, compute sigma.inv
    for(p in seq_len(Mp$npatterns)) {

        # observed variables
        var.idx <- Mp$pat[p,]

        # missing values for this pattern
        na.idx <- which(!var.idx)

        # invert Sigma for this pattern
        if(length(na.idx) > 0L) {
            sigma.inv <- lav_matrix_symmetric_inverse_update(S.inv = Sigma.inv,
                             rm.idx = na.idx, logdet = FALSE)
        } else {
            sigma.inv <- Sigma.inv
        }

        S.inv <- matrix(0, P, P)
        S.inv[var.idx, var.idx] <- sigma.inv

        S2.inv <- 0.5 * lav_matrix_duplication_pre_post(S.inv %x% S.inv)

        I11 <- I11 + Mp$freq[p] * S.inv
        I22 <- I22 + Mp$freq[p] * S2.inv
    }

    lav_matrix_bdiag(I11, I22)/N
}

# 5b: unit observed information Mu and vech(Sigma) from raw data
# (old term: Abeta, observed)
lav_mvnorm_missing_information_observed_data <- function(Y           = NULL,
                                                         Mp          = NULL,
                                                         Mu          = NULL,
                                                         Sigma       = NULL,
                                                         Sinv.method = "eigen",
                                                         Sigma.inv   = NULL) {
    # missing patterns
    if(is.null(Mp)) {
        Mp <- lav_data_missing_patterns(Y)
    }

    # N
    N <- sum(Mp$freq)

    # observed information
    observed <- lav_mvnorm_missing_logl_hessian_data(Y = Y, Mp = Mp, Mu = Mu,
                    Sigma = Sigma, Sinv.method = Sinv.method, 
                    Sigma.inv = Sigma.inv)

    -observed/N
}

# 5b-bis: unit observed information Mu and vech(Sigma) from samplestats
lav_mvnorm_missing_information_observed_samplestats <-
    function(Yp          = NULL, 
             Mu          = NULL,
             Sigma       = NULL,
             Sinv.method = "eigen",
             Sigma.inv   = NULL) {

    N <- sum(sapply(Yp, "[[", "freq")) # implicitly: removed empty cases!
    
    # observed information
    observed <- lav_mvnorm_missing_logl_hessian_samplestats(Yp = Yp, Mu = Mu,
                    Sigma = Sigma, Sinv.method = Sinv.method,
                    Sigma.inv = Sigma.inv)
    
    -observed/N
}

# 5c: unit first-order information Mu and vech(Sigma) from raw data
# (old term: Bbeta)
lav_mvnorm_missing_information_firstorder <- function(Y           = NULL,
                                                      Mp          = NULL,
                                                      Mu          = NULL,
                                                      Sigma       = NULL,
                                                      Sinv.method = "eigen",
                                                      Sigma.inv   = NULL) {
    # missing patterns
    if(is.null(Mp)) {
        Mp <- lav_data_missing_patterns(Y)
    }

    # N
    N <- sum(Mp$freq)

    SC <- lav_mvnorm_missing_scores_mu_vech_sigma(Y = Y, Mp = Mp, Mu = Mu, 
              Sigma = Sigma, Sinv.method = Sinv.method, Sigma.inv = Sigma.inv)

    lav_matrix_crossprod(SC)/N
}

# 5d: both unit first-order information and expected/observed information
#     from raw data, in one go for efficiency
lav_mvnorm_missing_information_both <- function(Y           = NULL,
                                                Mp          = NULL,
                                                Mu          = NULL,
                                                Sigma       = NULL,
                                                Sinv.method = "eigen",
                                                Sigma.inv   = NULL,
                                                information = "observed") {
    P <- NCOL(Y); Mu <- as.numeric(Mu)

    if(is.null(Sigma.inv)) {
        # invert Sigma
        Sigma.inv <- lav_matrix_symmetric_inverse(S = Sigma, logdet = FALSE,
                                                  Sinv.method = Sinv.method)
    }
    
    # for the tcrossprod
    idx1 <- lav_matrix_vech_col_idx(P); idx2 <- lav_matrix_vech_row_idx(P)

    # vech(Sigma.inv)
    iSigma <- lav_matrix_vech(Sigma.inv)

    # missing patterns
    if(is.null(Mp)) {
        Mp <- lav_data_missing_patterns(Y)
    }

    if(information == "observed") {
        Yp <- lav_samplestats_missing_patterns(Y = Y, Mp = Mp)
    }

    # N
    N <- sum(Mp$freq)

    # subtract Mu
    Yc <- t( t(Y) - Mu )

    # dmu per case
    dmu <- matrix(as.numeric(NA), nrow = NROW(Y), ncol = P)

    # SC
    SC <- matrix(as.numeric(NA), nrow = NROW(Y), ncol = length(iSigma))

    # expected/observed information
    I11 <- matrix(0, P, P)
    I22 <- matrix(0, P*(P+1)/2, P*(P+1)/2)
    if(information == "observed") {
        I21 <- matrix(0, P*(P+1)/2, P)
    }

    # for each pattern, compute Yc %*% sigma.inv
    for(p in seq_len(Mp$npatterns)) {

        # observed values for this pattern
        var.idx <- Mp$pat[p,]

        # missing values for this pattern
        na.idx <- which(!var.idx)

        # cases with this pattern
        case.idx <- Mp$case.idx[[p]]

        # invert Sigma for this pattern
        if(length(na.idx) > 0L) {
            sigma.inv <- lav_matrix_symmetric_inverse_update(S.inv = Sigma.inv,
                             rm.idx = na.idx, logdet = FALSE)
            tmp <- matrix(0, P, P)
            tmp[var.idx, var.idx] <- sigma.inv
            isigma <- lav_matrix_vech(tmp)
        } else {
            sigma.inv <- Sigma.inv
            isigma <- iSigma
        }

        # information
        S.inv <- matrix(0, P, P)
        S.inv[var.idx, var.idx] <- sigma.inv

        if(information == "expected") {
            S2.inv <- 0.5 * lav_matrix_duplication_pre_post(S.inv %x% S.inv)

            I11 <- I11 + Mp$freq[p] * S.inv
            I22 <- I22 + Mp$freq[p] * S2.inv
        } else {
            pat.freq <- Yp[[p]]$freq
        
            tmp21 <- matrix(0,P,1)
            tmp21[var.idx,1] <- sigma.inv %*% (Yp[[p]]$MY - Mu[var.idx])
        
            W.tilde <- Yp[[p]]$SY + tcrossprod(Yp[[p]]$MY - Mu[var.idx])
            AAA <- ( sigma.inv %*%  
                     (2*W.tilde - Sigma[var.idx,var.idx,drop = FALSE]) %*%
                     sigma.inv )
            tmp22 <- matrix(0, P, P)
            tmp22[var.idx, var.idx] <- AAA
        
            i11 <- S.inv
            i21 <- lav_matrix_duplication_pre( tmp21 %x% S.inv )
            i22 <- (1/2) * lav_matrix_duplication_pre_post(S.inv %x% tmp22)
        
            I11 <- I11 + pat.freq * i11
            I21 <- I21 + pat.freq * i21
            I22 <- I22 + pat.freq * i22
        }

        # compute dMu for all observations of this pattern
        dmu[case.idx, var.idx] <-
            Yc[case.idx, var.idx, drop = FALSE] %*% sigma.inv

        # postmultiply these cases with sigma.inv
        Yc[case.idx, var.idx] <- Yc[case.idx, var.idx] %*% sigma.inv

        # tcrossprod
        SC[case.idx,] <- Yc[case.idx, idx1] * Yc[case.idx, idx2]

        # substract isigma from each row
        SC[case.idx,] <- t( t(SC[case.idx,,drop = FALSE]) - isigma )
    }

    # adjust for vech
    SC[,lav_matrix_diagh_idx(P)] <- SC[,lav_matrix_diagh_idx(P)] / 2

    # add dmu
    SC <- cbind(dmu, SC)

    # first order information
    Bbeta <- lav_matrix_crossprod(SC)/N

    # expected/observed information
    if(information == "expected") {
        Abeta <- lav_matrix_bdiag(I11, I22)/N
    } else {
        Abeta <-  rbind( cbind(I11, t(I21) ),
                         cbind(I21,   I22) )/N
    }

    list(Abeta = Abeta, Bbeta = Bbeta)
}



# impute missing cells, under the normal model, pattern-based
##
## FIXME: 
## 1) adopt coding style (see above) for this function
## 2) adopt estimate.moments.EM
lav_mvnorm_missing_impute_pattern <- function(Y           = NULL,
                                              Mp          = NULL,
                                              Mu          = NULL,
                                              Sigma       = NULL,
                                              Sinv.method = "eigen") {

    # complete data
    Y.complete <- Y

    for(p in seq_len(Mp$npatterns)) {

        var.idx <- Mp$pat[p,]
        if(all(var.idx)) {
            next
        }

        # extract raw data for these cases
        X <- Y[Mp$case.idx[[p]], Mp$pat[p, ], drop = FALSE]

        # partition Mu (1=missing, 2=complete)
        Mu_1 <- Mu[!var.idx]
        Mu_2 <- Mu[ var.idx]

        # partition Sigma (1=missing, 2=complete)
        Sigma_11 <- Sigma[!var.idx, !var.idx, drop=FALSE]
        Sigma_12 <- Sigma[!var.idx,  var.idx, drop=FALSE]
        Sigma_21 <- Sigma[ var.idx, !var.idx, drop=FALSE]
        Sigma_22 <- Sigma[ var.idx,  var.idx, drop=FALSE]
        Sigma_22.inv <- try(inv.chol(Sigma_22, logdet=FALSE), silent = TRUE)
        if(inherits(Sigma_22.inv, "try-error")) {
            stop("lavaan ERROR: Sigma_22.inv cannot be inverted")
        }
        #Sigma_22.inv <- solve(Sigma_22)

        # estimate missing values in this pattern
        Diff <- apply(X, 1, '-', Mu_2)
        X_missing2 <- t(Sigma_12 %*% Sigma_22.inv %*% Diff)
        X_missing <- t(apply(X_missing2, 1, '+', Mu_1))

        # complete data for this pattern
        Y.complete[Mp$case.idx[[p]], !var.idx] <- X_missing
    }

    Y.complete
}

