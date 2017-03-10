# loglikelihood clustered data

lav_mvnorm_cluster_loglik_samplestats_2l <- function(YLp          = NULL,
                                                     Lp           = NULL,
                                                     Sigma.W      = NULL,
                                                     Mu.B         = NULL,
                                                     Sigma.B      = NULL,
                                                     Sinv.method  = "eigen",
                                                     method       = "size",
                                                     log2pi       = FALSE,
                                                     minus.two    = TRUE) {

    between.idx <- Lp$between.idx[[2]]

    # map to matrices needed for loglik
    if(length(between.idx) > 0L) {
        mu.x <- Mu.B[ between.idx ]
        mu.y <- Mu.B[-between.idx ]

        sigma.xx <- Sigma.B[ between.idx, between.idx, drop = FALSE]
        sigma.yx <- Sigma.B[-between.idx, between.idx, drop = FALSE]
        sigma.2  <- Sigma.B[-between.idx,-between.idx, drop = FALSE]
    } else {
        mu.x <- numeric(0L)
        mu.y <- Mu.B
        sigma.xx <- matrix(0, 0L, 0L)
        sigma.yx <- matrix(0, nrow(Sigma.B), 0L)
        sigma.2 <- Sigma.B
    }
    sigma.1 <- Sigma.W

    if(method == "size") {
        loglik <- .loglik_cluster_2l_size(YLp = YLp, Lp = Lp,
                      mu.y = mu.y, mu.x = mu.x, sigma.1 = sigma.1,
                      sigma.2 = sigma.2, sigma.xx = sigma.xx,
                      sigma.yx = sigma.yx, Sinv.method = Sinv.method)
    } else if(method == "cluster") {
        loglik <- .loglik_cluster_2l_cluster(YLp = YLp, Lp = Lp,
                      mu.y = mu.y, mu.x = mu.x, sigma.1 = sigma.1,
                      sigma.2 = sigma.2, sigma.xx = sigma.xx,
                      sigma.yx = sigma.yx, Sinv.method = Sinv.method)
    } else {
        stop("lavaan ERROR: method `", method, " not supported")
    }

    loglik
}


.loglik_cluster_2l_cluster <- function(YLp         = NULL,
                                       Lp          = NULL,
                                       mu.y        = NULL,
                                       mu.x        = NULL,
                                       sigma.1     = NULL,
                                       sigma.2     = NULL,
                                       sigma.xx    = NULL,
                                       sigma.yx    = NULL,
                                       Sinv.method = "eigen") {

}

.loglik_cluster_2l_size <- function(YLp         = NULL,
                                    Lp          = NULL,
                                    mu.y        = NULL,
                                    mu.x        = NULL,
                                    sigma.1     = NULL,
                                    sigma.2     = NULL,
                                    sigma.xx    = NULL,
                                    sigma.yx    = NULL,
                                    Sinv.method = "eigen") {

    between.idx <- Lp$between.idx[[2]]

    # inverses
    sigma.1.inv <- lav_matrix_symmetric_inverse(S = sigma.1, logdet = TRUE,
                                                Sinv.method = Sinv.method)
    sigma.1.logdet <- attr(sigma.1.inv, "logdet")
    attr(sigma.1.inv, "logdet") <- NULL

    cluster.sizes  <- Lp$cluster.sizes[[2L]]
    ncluster.sizes <- Lp$ncluster.sizes[[2L]]
    B  <- numeric(ncluster.sizes)
    LB <- numeric(ncluster.sizes)

    GD     <- YLp[[2]]$GD
    cov.d  <- YLp[[2]]$cov.d
    mean.d <- YLp[[2]]$mean.d

    if(length(between.idx) > 0L) {
        sigma.xx.inv <- 
            lav_matrix_symmetric_inverse(S = sigma.xx, logdet = FALSE,
                                         Sinv.method = Sinv.method)

        # components A matrix
        A11 <- sigma.xx
        A12 <- t(sigma.yx)
        A21 <- sigma.yx
        m.k <- c(mu.x, mu.y)

        S.PW <- YLp[[2]]$Sigma.W[-between.idx,-between.idx, drop = FALSE]
    } else {
        S.PW <- YLp[[2]]$Sigma.W
        m.k <- mu.y
    }


    # min 2* logliklihood
    loglik <- 0
    for(clz in seq_len(ncluster.sizes)) {
        nd <- cluster.sizes[clz]

        # construct A matrix (Vk.mini)
        if(length(between.idx) > 0L) {
            A22 <- (1/nd * sigma.1) + sigma.2
            A <- rbind( cbind(A11, A12), cbind(A21, A22) )

            A.11 <- solve(A11 - A12 %*% solve(A22) %*% A21)
            A.22 <- solve(A22 - A21 %*% sigma.xx.inv %*% A12)
            A.12 <- - sigma.xx.inv %*% A12 %*% A.22
            A.21 <- t(A.12)
            A.inv <- rbind( cbind(A.11, A.12), cbind(A.21, A.22) )
        } else {
            A <- (1/nd * sigma.1) + sigma.2
            A.inv <- solve(A)
        }

        COV  <- cov.d[[clz]]
        MEAN <- mean.d[[clz]]
        TT <- COV + tcrossprod(MEAN - m.k)

        B[clz] <- GD[[clz]] * sum( A.inv * TT ) # eq 33, second 

        # logdet (between only)
        LB[clz] <- GD[[clz]]*log(det(A)) + GD[[clz]]*nrow(sigma.1)*log(nd)
    }

    LW <- ( (Lp$nclusters[[1]] - Lp$nclusters[[2]]) *
            (sigma.1.logdet + sum(sigma.1.inv * S.PW)) )
    loglik <- sum(LB) + sum(B) + LW
    
    loglik
}
