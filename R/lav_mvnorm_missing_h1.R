# the Multivariate normal distribution, unrestricted (h1), missing values

# 1) loglikelihood --> same as h0 but where Mu and Sigma are unrestricted
# 2) 3) 4) 5)      --> (idem)

# YR 26 Mar 2016: first version
# YR 20 Jan 2017: added _h1_omega_sw()

# here, we estimate Mu and Sigma from Y with missing values, assuming normality
# this is a rewrite of the 'estimate.moments.EM' function in <= 0.5-22
lav_mvnorm_missing_h1_estimate_moments <- function(Y           = NULL,
                                                   Mp          = NULL,
                                                   Yp          = NULL,
                                                   wt          = NULL,
                                                   Sinv.method = "eigen",
                                                   verbose     = FALSE,
                                                   max.iter    = 500L,
                                                   tol         = 1e-05) {

    # check input
    Y <- as.matrix(Y); P <- NCOL(Y); N <- NROW(Y)

    # missing patterns
    if(is.null(Mp)) {
        Mp <- lav_data_missing_patterns(Y)
    }
    if(is.null(Yp)) {
        Yp <- lav_samplestats_missing_patterns(Y = Y, Mp = Mp, wt = wt)
    }

    # remove empty cases
    if(length(Mp$empty.idx) > 0L) {
        N <- N - length(Mp$empty.idx)
    }

    # verbose?
    if(verbose) {
        cat("\n")
        cat("lav_mvnorm_missing_h1_estimate_moments: start EM steps\n")
    }

    # starting values; zero covariances to guarantee a pd matrix
    Mu0 <- base::.colMeans(Y, m = N, n = P, na.rm = TRUE)
    var0 <- base::.colMeans(Y*Y, m = N, n = P, na.rm = TRUE) - Mu0*Mu0
    Sigma0 <- diag(x = var0, nrow = P)
    Mu <- Mu0; Sigma <- Sigma0

    # report
    if(verbose) {
        #fx0 <- estimator.FIML(Sigma.hat=Sigma, Mu.hat=Mu, M=Yp)
        fx0 <- lav_mvnorm_missing_loglik_samplestats(Yp = Yp,
                                                     Mu = Mu, Sigma = Sigma,
                                                     log2pi = FALSE,
                                                     minus.two = TRUE)/N
        cat("  EM iteration:", sprintf("%4d", 0),
            " fx = ", sprintf("%15.10f", fx0),
            "\n")
    }

    # EM steps
    for(i in 1:max.iter) {

        # E-step
        Estep <- lav_mvnorm_missing_estep(Y = Y, Mp = Mp, wt = wt,
                                          Mu = Mu, Sigma = Sigma,
                                          Sinv.method = Sinv.method)
        T1 <- Estep$T1
        T2 <- Estep$T2

        # M-step
        Mu    <- T1/N
        Sigma <- T2/N - tcrossprod(Mu)

        # check if Sigma is near-pd (+ poor fix)
        ev <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)
        tol <- 1e-6 # FIXME!
        if(any(ev$values < tol)) {
            #too.small <- which( ev$values < tol )
            #ev$values[too.small] <- tol
            #ev$values <- ev$values + tol
            #Sigma <- ev$vectors %*% diag(ev$values) %*% t(ev$vectors)

            # ridge
            diag(Sigma) <- diag(Sigma) + max(diag(Sigma))*1e-08
        }

        # max absolute difference in parameter values
        DELTA <- max(abs(c(Mu,  lav_matrix_vech(Sigma)) -
                         c(Mu0, lav_matrix_vech(Sigma0))))

        # report fx
        if(verbose) {
            #fx <- estimator.FIML(Sigma.hat=Sigma, Mu.hat=Mu, M=Yp)
            fx <- lav_mvnorm_missing_loglik_samplestats(Yp = Yp,
                                                        Mu = Mu, Sigma = Sigma,
                                                        log2pi = FALSE,
                                                        minus.two = TRUE)/N
            cat("  EM iteration:", sprintf("%4d", i),
                " fx = ", sprintf("%15.10f", fx),
                " delta par = ", sprintf("%9.8f", DELTA),
                "\n")
        }

        # convergence check: using parameter values:
        if(DELTA < tol)
            break

        # again
        Mu0 <- Mu; Sigma0 <- Sigma

    } # EM iterations

    if(verbose) {
        cat("\nSigma:\n"); print(Sigma)
        cat("\nMu:\n"); print(Mu)
        cat("\n")
    }

    # compute fx if we haven't already
    if(!verbose) {
        #fx <- estimator.FIML(Sigma.hat = Sigma, Mu.hat = Mu, M = Yp)
        fx <- lav_mvnorm_missing_loglik_samplestats(Yp = Yp,
                                                    Mu = Mu, Sigma = Sigma,
                                                    log2pi = FALSE,
                                                    minus.two = TRUE)/N
    }

    list(Sigma = Sigma, Mu = Mu, fx = fx)
}

# compute N times ACOV(Mu, vech(Sigma))
# in the literature: - `Omega_{SW}'
#                    - `Gamma for incomplete data'
#                    - (N times the) sandwich estimator for acov(mu,vech(Sigma))
lav_mvnorm_missing_h1_omega_sw <- function(Y           = NULL,
                                           Mp          = NULL,
                                           wt          = NULL,
                                           cluster.idx = NULL,
                                           Yp          = NULL,
                                           Sinv.method = "eigen",
                                           Mu          = NULL,
                                           Sigma       = NULL,
                                           x.idx       = NULL,
                                           Sigma.inv   = NULL,
                                           information = "observed") {

    # missing patterns
    if(is.null(Mp)) {
        Mp <- lav_data_missing_patterns(Y)
    }

    # sample stats per pattern
    if(is.null(Yp) && (information == "observed" || is.null(Sigma))) {
        Yp <- lav_samplestats_missing_patterns(Y = Y, Mp = Mp, wt = wt)
    }

    # Sigma and Mu
    if(is.null(Sigma) || is.null(Mu)) {
        out <- lav_mvnorm_missing_h1_estimate_moments(Y = Y, Mp = Mp, Yp = Yp)
        Mu <- out$Mu
        Sigma <- out$Sigma
    }

    # information matrices
    info <- lav_mvnorm_missing_information_both(Y = Y, Mp = Mp, Mu = Mu,
                wt = wt, cluster.idx = cluster.idx,
                Sigma = Sigma, x.idx = x.idx, Sinv.method = Sinv.method,
                Sigma.inv = Sigma.inv, information = information)

    A <- info$Abeta
    A.inv <- lav_matrix_symmetric_inverse(S = A, logdet = FALSE,
                                          Sinv.method = Sinv.method)
    B <- info$Bbeta

    # sandwich
    SW <- A.inv %*% B %*% A.inv

    SW
}


