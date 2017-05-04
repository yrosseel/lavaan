# loglikelihood clustered data

# YR: first version around Feb 2017


# take model-implied mean+variance matrices, and reorder/augment them
# to facilitate computing of (log)likelihood in the two-level case
#
# - sigma.w and sigma.b: same dimensions, level-1 variables only
# - sigma.zz: level-2 variables only
# - sigma.yz: cov(level-1, level-2)
# - mu.y: level-1 variables only
# - my.z: level-2 variables only
lav_mvnorm_cluster_implied22l <- function(Lp           = NULL,
                                          implied      = NULL,
                                          Mu.W         = NULL,
                                          Mu.B         = NULL,
                                          Sigma.W      = NULL,
                                          Sigma.B      = NULL) {

    if(!is.null(implied)) {
        Sigma.W <- implied$cov[[1]]
        Mu.W    <- implied$mean[[1]]

        Sigma.B <- implied$cov[[2]]
        Mu.B    <- implied$mean[[2]]
    }

    # between.idx 
    between.idx <- Lp$between.idx[[2]]

    # ov.idx per level
    ov.idx      <- Lp$ov.idx

    # 'tilde' matrices: ALL variables within and between
    p.tilde <- length( unique(c(ov.idx[[1]], ov.idx[[2]])) )

    # Sigma.W.tilde
    Sigma.W.tilde <- matrix(0, p.tilde, p.tilde)
    Sigma.W.tilde[ ov.idx[[1]], ov.idx[[1]] ] <- Sigma.W

    # Sigma.B.tilde
    Sigma.B.tilde <- matrix(0, p.tilde, p.tilde)
    Sigma.B.tilde[ ov.idx[[2]], ov.idx[[2]] ] <- Sigma.B

    # Mu.W.tilde
    Mu.W.tilde <- numeric( p.tilde )
    Mu.W.tilde[ ov.idx[[1]] ] <- Mu.W

    # Mu.B.tilde
    Mu.B.tilde <- numeric( p.tilde )
    Mu.B.tilde[ ov.idx[[2]] ] <- Mu.B

    # add Mu.W[within.idx] to Mu.B
    # Mu.B.tilde[ Lp$within.idx[[2]] ] <- Mu.W.tilde[ Lp$within.idx[[2]] ]
    Mu.B.tilde[ ov.idx[[1]] ] <- ( Mu.B.tilde[ ov.idx[[1]] ] +
                                   Mu.W.tilde[ ov.idx[[1]] ] )

    # map to matrices needed for loglik
    if(length(between.idx) > 0L) {
        mu.z <- Mu.B.tilde[ between.idx ]
        mu.y <- Mu.B.tilde[-between.idx ]

        sigma.zz <- Sigma.B.tilde[ between.idx, between.idx, drop = FALSE]
        sigma.yz <- Sigma.B.tilde[-between.idx, between.idx, drop = FALSE]
        sigma.b  <- Sigma.B.tilde[-between.idx,-between.idx, drop = FALSE]
        sigma.w  <- Sigma.W.tilde[-between.idx,-between.idx, drop = FALSE]
    } else {
        mu.z <- numeric(0L)
        mu.y <- Mu.B.tilde
        sigma.zz <- matrix(0, 0L, 0L)
        sigma.yz <- matrix(0, nrow(Sigma.B.tilde), 0L)
        sigma.b <- Sigma.B.tilde
        sigma.w <- Sigma.W.tilde
    }

    list(sigma.w = sigma.w, sigma.b = sigma.b, sigma.zz = sigma.zz,
         sigma.yz = sigma.yz, mu.z = mu.z, mu.y = mu.y)
}

# Mu.W, Mu.B, Sigma.W, Sigma.B are the model-implied statistics
# (not yet reordered)
lav_mvnorm_cluster_loglik_samplestats_2l <- function(YLp          = NULL,
                                                     Lp           = NULL,
                                                     Mu.W         = NULL,
                                                     Sigma.W      = NULL,
                                                     Mu.B         = NULL,
                                                     Sigma.B      = NULL,
                                                     Sinv.method  = "eigen",
                                                     log2pi       = FALSE,
                                                     minus.two    = TRUE) {

    # map implied to 2l matrices
    out <- lav_mvnorm_cluster_implied22l(Lp = Lp, Mu.W = Mu.W, Mu.B = Mu.B,
               Sigma.W = Sigma.W, Sigma.B = Sigma.B)
    mu.y <- out$mu.y; mu.z <- out$mu.z
    sigma.w <- out$sigma.w; sigma.b <- out$sigma.b 
    sigma.zz <- out$sigma.zz; sigma.yz <- out$sigma.yz

    # Lp
    nclusters       <- Lp$nclusters[[2]]
    cluster.size    <- Lp$cluster.size[[2]]
    between.idx     <- Lp$between.idx[[2]]
    cluster.sizes   <- Lp$cluster.sizes[[2]]
    ncluster.sizes  <- Lp$ncluster.sizes[[2]]
    cluster.size.ns <- Lp$cluster.size.ns[[2]]

    # Y1
    if(length(between.idx) > 0L) {
        S.PW <- YLp[[2]]$Sigma.W[-between.idx, -between.idx, drop = FALSE]
    } else {
        S.PW <- YLp[[2]]$Sigma.W
    }

    # Y2
    #Y2 <- YLp[[2]]$Y2
    #Y2.cm <- t( t(Y2) - mu.b )
    cov.d <- YLp[[2]]$cov.d
    mean.d <- YLp[[2]]$mean.d

    # common parts:
    sigma.w.inv <- lav_matrix_symmetric_inverse(S = sigma.w,
                       logdet = TRUE, Sinv.method = Sinv.method)
    sigma.w.logdet <- attr(sigma.w.inv, "logdet")
    attr(sigma.w.inv, "logdet") <- NULL

    if(length(between.idx) > 0L) {
        sigma.zz.inv <- lav_matrix_symmetric_inverse(S = sigma.zz,
                           logdet = TRUE, Sinv.method = Sinv.method)
        sigma.zz.logdet <- attr(sigma.zz.inv, "logdet")
        attr(sigma.zz.inv, "logdet") <- NULL
        sigma.yz.zi <- sigma.yz %*% sigma.zz.inv
        sigma.zi.zy <- t(sigma.yz.zi)
        sigma.b.z <- sigma.b - sigma.yz %*% sigma.zi.zy
    } else {
        sigma.zz.logdet <- 0
        sigma.b.z <- sigma.b
    }

    # min 2* logliklihood
    L <- numeric(ncluster.sizes) # logdet
    B <- numeric(ncluster.sizes) # between qf
    for(clz in seq_len(ncluster.sizes)) {

        # cluster size
        nj <- cluster.sizes[clz]

        # data between
        Y2Yc <- (cov.d[[clz]] + tcrossprod(mean.d[[clz]] - c(mu.z, mu.y)))

        # FIXME: avoid reorder/b.idx, so can use between.idx
        if(length(between.idx) > 0L) {
            b.idx <- seq_len(length(Lp$between.idx[[2]]))
            Y2Yc.zz <- Y2Yc[ b.idx, b.idx, drop = FALSE]
            Y2Yc.yz <- Y2Yc[-b.idx, b.idx, drop = FALSE]
            Y2Yc.yy <- Y2Yc[-b.idx,-b.idx, drop = FALSE]
        } else {
            Y2Yc.yy <- Y2Yc
        }

        # construct sigma.j
        sigma.j <- (nj * sigma.b.z) + sigma.w
        sigma.j.inv <- lav_matrix_symmetric_inverse(S = sigma.j,
                       logdet = TRUE, Sinv.method = Sinv.method)
        sigma.j.logdet <- attr(sigma.j.inv, "logdet")
        attr(sigma.j.inv, "logdet") <- NULL
        
        # logdet -- between only
        L[clz] <- (sigma.zz.logdet + sigma.j.logdet)

        if(length(between.idx) > 0L) {
            # part 1 -- zz
            sigma.ji.yz.zi <- sigma.j.inv %*% sigma.yz.zi
            Vinv.11 <- sigma.zz.inv + nj*(sigma.zi.zy %*% sigma.ji.yz.zi)
            q.zz <- sum(Vinv.11 * Y2Yc.zz)

            # part 2 -- yz
            q.yz <- - nj * sum(sigma.ji.yz.zi * Y2Yc.yz)
        } else {
            q.zz <- q.yz <- 0
        }

        # part 5 -- yyc
        q.yyc <- -nj * sum(sigma.j.inv * Y2Yc.yy )

        # qf -- between only
        B[clz] <- q.zz + 2*q.yz - q.yyc
    }
    # q.yya + q.yyb
    q.W <- sum(cluster.size - 1) * sum(sigma.w.inv * S.PW)
    # logdet within part
    L.W  <- sum(cluster.size - 1) * sigma.w.logdet

    # -2*times logl (without the constant)
    loglik <- sum(L * cluster.size.ns) + sum(B * cluster.size.ns) + q.W + L.W

    # functions below compute -2 * logl
    if(!minus.two) {
        loglik <- loglik / (-2)
    }

    # constant
    # Note: total 'N' = (nobs * #within vars) + (nclusters * #between vars)
    if(log2pi) {
        LOG.2PI <- log(2 * pi)
        nWithin <- length(c(Lp$both.idx[[2]], Lp$within.idx[[2]]))
        nBetween <- length(Lp$between.idx[[2]])
        P <- Lp$nclusters[[1]]*nWithin + Lp$nclusters[[2]]*nBetween
        constant <- -(P * LOG.2PI)/2
        loglik <- loglik + constant
    }

    loglik
}


