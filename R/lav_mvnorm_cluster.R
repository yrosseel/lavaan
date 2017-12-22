# loglikelihood clustered data

# YR: first version around Feb 2017


# take model-implied mean+variance matrices, and reorder/augment them
# to facilitate computing of (log)likelihood in the two-level case
#
# - sigma.w and sigma.b: same dimensions, level-1 variables only
# - sigma.zz: level-2 variables only
# - sigma.yz: cov(level-1, level-2)
# - mu.y: level-1 variables only (mu.w + mu.b)
# - mu.w: y within  part
# - mu.b: y between part
# - my.z: level-2 variables only
lav_mvnorm_cluster_implied22l <- function(Lp           = NULL,
                                          implied      = NULL,
                                          Mu.W         = NULL,
                                          Mu.B         = NULL,
                                          Sigma.W      = NULL,
                                          Sigma.B      = NULL) {

    if(!is.null(implied)) {
        # FIXME: only for single-group analysis!
        Sigma.W <- implied$cov[[1]]
        Mu.W    <- implied$mean[[1]]

        Sigma.B <- implied$cov[[2]]
        Mu.B    <- implied$mean[[2]]
    }

    # within/between.idx 
    between.idx <- Lp$between.idx[[2]]
    within.idx  <- Lp$within.idx[[2]]
    both.idx    <- Lp$both.idx[[2]]

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
    Mu.WB.tilde <- numeric( p.tilde) 
    Mu.WB.tilde[ within.idx ] <- Mu.W.tilde[ within.idx ]
    Mu.WB.tilde[ both.idx ] <- ( Mu.B.tilde[ both.idx ] +
                                 Mu.W.tilde[ both.idx ] )

    # map to matrices needed for loglik
    if(length(within.idx) > 0L) {
        Mu.B.tilde[within.idx] <- 0
    }
    if(length(between.idx) > 0L) {
        mu.z <- Mu.B.tilde[  between.idx ]
        mu.y <- Mu.WB.tilde[-between.idx ]
        mu.w <- Mu.W.tilde[ -between.idx ]
        mu.b <- Mu.B.tilde[ -between.idx ]
        sigma.zz <- Sigma.B.tilde[ between.idx, between.idx, drop = FALSE]
        sigma.yz <- Sigma.B.tilde[-between.idx, between.idx, drop = FALSE]
        sigma.b  <- Sigma.B.tilde[-between.idx,-between.idx, drop = FALSE]
        sigma.w  <- Sigma.W.tilde[-between.idx,-between.idx, drop = FALSE]
    } else {
        mu.z <- numeric(0L)
        mu.y <- Mu.WB.tilde
        mu.w <- Mu.W.tilde
        mu.b <- Mu.B.tilde
        sigma.zz <- matrix(0, 0L, 0L)
        sigma.yz <- matrix(0, nrow(Sigma.B.tilde), 0L)
        sigma.b <- Sigma.B.tilde
        sigma.w <- Sigma.W.tilde
    }

    list(sigma.w = sigma.w, sigma.b = sigma.b, sigma.zz = sigma.zz,
         sigma.yz = sigma.yz, mu.z = mu.z, mu.y = mu.y, mu.w = mu.w,
         mu.b = mu.b)
}

lav_mvnorm_cluster_2l2implied <- function(Lp, 
                                          sigma.w = NULL,
                                          sigma.b = NULL,
                                          sigma.zz = NULL,
                                          sigma.yz = NULL,
                                          mu.z     = NULL,
                                          mu.y     = NULL,
                                          mu.w     = NULL,
                                          mu.b     = NULL) {

    # between.idx 
    between.idx <- Lp$between.idx[[2]]
    within.idx  <- Lp$within.idx[[2]]

    # ov.idx per level
    ov.idx      <- Lp$ov.idx

    # 'tilde' matrices: ALL variables within and between
    p.tilde <- length( unique(c(ov.idx[[1]], ov.idx[[2]])) )

    # if we have mu.y, convert to mu.w and mu.b
    if(!is.null(mu.y)) {
        mu.b <- mu.y
        mu.w.tilde <- numeric( p.tilde )
        mu.w.tilde[ ov.idx[[1]] ] <- mu.y
        if(length(within.idx) > 0L) {
            mu.w.tilde[  -within.idx ] <- 0
        } else {
            mu.w.tilde[] <- 0
        }
        mu.w <- mu.w.tilde[ ov.idx[[1]] ]
    }

    Mu.W.tilde <- numeric( p.tilde )
########## DEBUG ##############
    #if(length(within.idx) > 0) {
        Mu.W.tilde[ ov.idx[[1]] ] <- mu.w
    #}
###############################
    Mu.W <- Mu.W.tilde[ ov.idx[[1]] ]

    # Mu.B
    Mu.B.tilde <- numeric(p.tilde)
    Mu.B.tilde[ ov.idx[[1]] ] <- mu.b
    Mu.B.tilde[ between.idx ] <- mu.z
    if(length(within.idx) > 0) {
        Mu.B.tilde[within.idx] <- 0
    }
    Mu.B <- Mu.B.tilde[ ov.idx[[2]] ]

    # Sigma.W
    Sigma.W <- sigma.w

    # Sigma.B
    Sigma.B.tilde <- matrix(0, p.tilde, p.tilde)
    Sigma.B.tilde[ ov.idx[[1]], ov.idx[[1]] ] <- sigma.b
    Sigma.B.tilde[ ov.idx[[1]], between.idx ] <- sigma.yz
    Sigma.B.tilde[ between.idx, ov.idx[[1]] ] <- t(sigma.yz)
    Sigma.B.tilde[ between.idx, between.idx ] <- sigma.zz
    Sigma.B <- Sigma.B.tilde[ ov.idx[[2]], ov.idx[[2]], drop = FALSE ]

    list(Mu.W = Mu.W, Mu.B = Mu.B, Sigma.W = Sigma.W, Sigma.B = Sigma.B)
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

        # check: what if sigma.j is non-pd?
        if(is.na(sigma.j.logdet)) {
            # stop, and return NA right away
            #return(as.numeric(NA))
            # FORCE?
            #sigma.j <- lav_matrix_symmetric_force_pd(sigma.j)
            #sigma.j.inv <- lav_matrix_symmetric_inverse(S = sigma.j,
            #           logdet = TRUE, Sinv.method = Sinv.method)
            #sigma.j.logdet <- attr(sigma.j.inv, "logdet")
            #attr(sigma.j.inv, "logdet") <- NULL
        }
        
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


# first derivate -2*logl wrt Mu.W, Mu.B, Sigma.W, Sigma.B
lav_mvnorm_cluster_dlogl_2l_samplestats <- function(YLp          = NULL,
                                                    Lp           = NULL,
                                                    Mu.W         = NULL,
                                                    Sigma.W      = NULL,
                                                    Mu.B         = NULL,
                                                    Sigma.B      = NULL,
                                                    Sinv.method  = "eigen") {

    # map implied to 2l matrices
    out <- lav_mvnorm_cluster_implied22l(Lp = Lp, Mu.W = Mu.W, Mu.B = Mu.B,
               Sigma.W = Sigma.W, Sigma.B = Sigma.B)
    mu.y <- out$mu.y; mu.z <- out$mu.z
    sigma.w <- out$sigma.w; sigma.b <- out$sigma.b 
    sigma.zz <- out$sigma.zz; sigma.yz <- out$sigma.yz

    # Lp
    nclusters       <- Lp$nclusters[[2]]
    cluster.size    <- Lp$cluster.size[[2]]
    cluster.sizes   <- Lp$cluster.sizes[[2]]
    cluster.idx     <- Lp$cluster.idx[[2]]
    between.idx     <- Lp$between.idx[[2]]
    ncluster.sizes  <- Lp$ncluster.sizes[[2]]
    cluster.size.ns <- Lp$cluster.size.ns[[2]]

    # Y1
    if(length(between.idx) > 0L) {
        S.PW <- YLp[[2]]$Sigma.W[-between.idx, -between.idx, drop = FALSE]
    } else {
        S.PW <- YLp[[2]]$Sigma.W
    }

    # Y2
    cov.d <- YLp[[2]]$cov.d
    mean.d <- YLp[[2]]$mean.d

    # common parts:
    sigma.w.inv <- lav_matrix_symmetric_inverse(S = sigma.w,
                       logdet = FALSE, Sinv.method = Sinv.method)

    # both level-1 and level-2
    G.muy      <- matrix(0, ncluster.sizes, length(mu.y))
    G.Sigma.w  <- matrix(0, ncluster.sizes, length(lav_matrix_vech(sigma.w)))
    G.Sigma.b  <- matrix(0, ncluster.sizes, length(lav_matrix_vech(sigma.b)))

    if(length(between.idx) > 0L) {

        G.muz      <- matrix(0, ncluster.sizes, length(mu.z))
        G.Sigma.zz <- matrix(0, ncluster.sizes, 
                                length(lav_matrix_vech(sigma.zz)))
        G.Sigma.yz <- matrix(0, ncluster.sizes, length(lav_matrix_vec(sigma.yz)))

        sigma.zz.inv <- lav_matrix_symmetric_inverse(S = sigma.zz,
                           logdet = FALSE, Sinv.method = Sinv.method)
        sigma.yz.zi <- sigma.yz %*% sigma.zz.inv
        sigma.zi.zy <- t(sigma.yz.zi)
        sigma.b.z <- sigma.b - sigma.yz %*% sigma.zi.zy


        for(clz in seq_len(ncluster.sizes)) {
            # cluster size
            nj <- cluster.sizes[clz]

            # level-2 vectors
            b.idx <- seq_len(length(Lp$between.idx[[2]]))
            zyc <- mean.d[[clz]] - c(mu.z, mu.y)
            yc <- zyc[-b.idx]
            zc <- zyc[ b.idx]

            # level-2 crossproducts
            Y2Yc <- (cov.d[[clz]] + tcrossprod(mean.d[[clz]] - c(mu.z, mu.y)))
            b.idx <- seq_len(length(Lp$between.idx[[2]]))
            Y2Yc.zz <- Y2Yc[ b.idx, b.idx, drop = FALSE]
            Y2Yc.yz <- Y2Yc[-b.idx, b.idx, drop = FALSE]
            Y2Yc.yy <- Y2Yc[-b.idx,-b.idx, drop = FALSE]

            # construct sigma.j
            sigma.j <- (nj * sigma.b.z) + sigma.w
            sigma.j.inv <- lav_matrix_symmetric_inverse(S = sigma.j,
                               logdet = FALSE, Sinv.method = Sinv.method)
            sigma.ji.yz.zi <- sigma.j.inv %*% sigma.yz.zi
            sigma.zi.zy.ji <- t(sigma.ji.yz.zi)

            # common parts
            jYZj <- nj * (   sigma.j.inv %*%
                            (sigma.yz.zi %*% Y2Yc.zz %*% t(sigma.yz.zi)
                               - Y2Yc.yz %*% t(sigma.yz.zi)
                             - t(Y2Yc.yz %*% t(sigma.yz.zi)) + Y2Yc.yy)
                         %*% sigma.j.inv )

            Z1 <- Y2Yc.zz %*% t(sigma.ji.yz.zi) %*% sigma.yz
            YZ1 <- t(Y2Yc.yz) %*% sigma.j.inv %*% sigma.yz


            # Mu.Z
            G.muz[clz,] <- -2 * as.numeric(
                (sigma.zz.inv + nj*(sigma.zi.zy.ji %*% sigma.yz.zi)) %*% zc
                 -nj*sigma.zi.zy.ji %*% yc)

            # MU.Y
            G.muy[clz,] <- 2*nj * as.numeric(zc %*% sigma.zi.zy.ji -
                                            yc %*% sigma.j.inv)

            # SIGMA.W (between part)
            g.sigma.w <- sigma.j.inv - jYZj
            tmp <- g.sigma.w*2; diag(tmp) <- diag(g.sigma.w)
            G.Sigma.w[clz,] <- lav_matrix_vech(tmp)

            # SIGMA.B
            g.sigma.b <- nj * (sigma.j.inv - jYZj)
            tmp <- g.sigma.b*2; diag(tmp) <- diag(g.sigma.b)
            G.Sigma.b[clz,] <- lav_matrix_vech(tmp)

            # SIGMA.ZZ
            g.sigma.zz <- ( sigma.zz.inv + nj * sigma.zz.inv %*% (
                              t(sigma.yz) %*% (sigma.j.inv - jYZj) %*% sigma.yz
                            -(1/nj * Y2Yc.zz + t(Z1) + Z1 - t(YZ1) - YZ1) ) %*%
                                                sigma.zz.inv )

            tmp <- g.sigma.zz*2; diag(tmp) <- diag(g.sigma.zz)
            G.Sigma.zz[clz,] <- lav_matrix_vech(tmp)

            # SIGMA.ZY
            g.sigma.yz <- 2 * nj * (
                          (sigma.j.inv %*%
                              (sigma.yz.zi %*% Y2Yc.zz - sigma.yz - Y2Yc.yz)
                               + jYZj %*% sigma.yz) %*% sigma.zz.inv )

            G.Sigma.yz[clz,] <- lav_matrix_vec(g.sigma.yz)
        }

        # level-1
        d.mu.y     <- colSums(G.muy * cluster.size.ns)
        d.sigma.w1 <- lav_matrix_vech_reverse(colSums(G.Sigma.w * 
                                                      cluster.size.ns))
        d.sigma.b  <- lav_matrix_vech_reverse(colSums(G.Sigma.b * 
                                                      cluster.size.ns))

        # level-2
        d.mu.z     <- colSums(G.muz * cluster.size.ns)
        d.sigma.zz <- lav_matrix_vech_reverse(colSums(G.Sigma.zz * 
                                                      cluster.size.ns))
        d.sigma.yz <- matrix(colSums(G.Sigma.yz * cluster.size.ns), 
                             nrow(sigma.yz), ncol(sigma.yz))
    } # between.idx 

    else { # no level-2 variables

        for(clz in seq_len(ncluster.sizes)) {
            # cluster size
            nj <- cluster.sizes[clz]

            # level-2 vectors
            yc <- mean.d[[clz]] - mu.y

            # level-2 crossproducts
            Y2Yc.yy <- (cov.d[[clz]] + tcrossprod(mean.d[[clz]] - mu.y))

            # construct sigma.j
            sigma.j <- (nj * sigma.b) + sigma.w
            sigma.j.inv <- lav_matrix_symmetric_inverse(S = sigma.j,
                               logdet = FALSE, Sinv.method = Sinv.method)
            # common part 
            jYYj <- nj * sigma.j.inv %*% Y2Yc.yy %*% sigma.j.inv

            # MU.Y
            G.muy[clz,] <- -2*nj * as.numeric(yc %*% sigma.j.inv)

            # SIGMA.W (between part)
            g.sigma.w <- sigma.j.inv - jYYj
            tmp <- g.sigma.w*2; diag(tmp) <- diag(g.sigma.w)
            G.Sigma.w[clz,] <- lav_matrix_vech(tmp)

            # SIGMA.B
            g.sigma.b <- nj * (sigma.j.inv - jYYj)
            tmp <- g.sigma.b*2; diag(tmp) <- diag(g.sigma.b)
            G.Sigma.b[clz,] <- lav_matrix_vech(tmp)
        }

        # level-1
        d.mu.y     <- colSums(G.muy * cluster.size.ns)
        d.sigma.w1 <- lav_matrix_vech_reverse(colSums(G.Sigma.w * 
                                                      cluster.size.ns))
        d.sigma.b  <- lav_matrix_vech_reverse(colSums(G.Sigma.b * 
                                                      cluster.size.ns))
        # level-2
        d.mu.z     <- numeric(0L)
        d.sigma.zz <- matrix(0, 0L, 0L)
        d.sigma.yz <- matrix(0, 0L, 0L)
    }

    # Sigma.W (bis)
    d.sigma.w2 <- (Lp$nclusters[[1]] - nclusters) * ( sigma.w.inv
                   - sigma.w.inv %*% S.PW %*% sigma.w.inv )
    tmp <- d.sigma.w2*2; diag(tmp) <- diag(d.sigma.w2)
    d.sigma.w2 <- tmp

    d.sigma.w <- d.sigma.w1 + d.sigma.w2

    # rearrange 
    dout <- lav_mvnorm_cluster_2l2implied(Lp = Lp,
                sigma.w = d.sigma.w, sigma.b = d.sigma.b,
                sigma.yz = d.sigma.yz, sigma.zz = d.sigma.zz,
                mu.y = d.mu.y, mu.z = d.mu.z)

    dout
}

# cluster-wise scores -2*logl wrt Mu.W, Mu.B, Sigma.W, Sigma.B
lav_mvnorm_cluster_scores_2l <- function(Y1           = NULL,
                                         YLp          = NULL,
                                         Lp           = NULL,
                                         Mu.W         = NULL,
                                         Sigma.W      = NULL,
                                         Mu.B         = NULL,
                                         Sigma.B      = NULL,
                                         Sinv.method  = "eigen") {

    # map implied to 2l matrices
    out <- lav_mvnorm_cluster_implied22l(Lp = Lp, Mu.W = Mu.W, Mu.B = Mu.B,
               Sigma.W = Sigma.W, Sigma.B = Sigma.B)
    mu.y <- out$mu.y; mu.z <- out$mu.z
    sigma.w <- out$sigma.w; sigma.b <- out$sigma.b 
    sigma.zz <- out$sigma.zz; sigma.yz <- out$sigma.yz

    # Lp
    nclusters       <- Lp$nclusters[[2]]
    cluster.size    <- Lp$cluster.size[[2]]
    cluster.idx     <- Lp$cluster.idx[[2]]
    between.idx     <- Lp$between.idx[[2]]

    # Y1
    if(length(between.idx) > 0L) {
        Y1w <- Y1[,-Lp$between.idx[[2]], drop = FALSE]
    } else {
        Y1w <- Y1
    }
    Y1w.cm <- t( t(Y1w) - mu.y )

    # Y2
    Y2 <- YLp[[2]]$Y2
    # NOTE: ORDER mu.b must match Y2
    mu.b <- numeric(ncol(Y2))
    if(length(between.idx) > 0L) {
        mu.b[-Lp$between.idx[[2]]] <- mu.y
        mu.b[ Lp$between.idx[[2]]] <- mu.z
    } else {
        mu.b <- mu.y
    }
    Y2.cm <- t( t(Y2) - mu.b )

    # common parts:
    sigma.w.inv <- lav_matrix_symmetric_inverse(S = sigma.w,
                       logdet = FALSE, Sinv.method = Sinv.method)

    # both level-1 and level-2
    G.muy      <- matrix(0, nclusters, length(mu.y))
    G.Sigma.w  <- matrix(0, nclusters, length(lav_matrix_vech(sigma.w)))
    G.Sigma.b  <- matrix(0, nclusters, length(lav_matrix_vech(sigma.b)))
    G.muz      <- matrix(0, nclusters, length(mu.z))
    G.Sigma.zz <- matrix(0, nclusters, length(lav_matrix_vech(sigma.zz)))
    G.Sigma.yz <- matrix(0, nclusters, length(lav_matrix_vec(sigma.yz)))

    if(length(between.idx) > 0L) {

        sigma.zz.inv <- lav_matrix_symmetric_inverse(S = sigma.zz,
                           logdet = FALSE, Sinv.method = Sinv.method)
        sigma.yz.zi <- sigma.yz %*% sigma.zz.inv
        sigma.zi.zy <- t(sigma.yz.zi)
        sigma.b.z <- sigma.b - sigma.yz %*% sigma.zi.zy


        for(cl in seq_len(nclusters)) {
            # cluster size
            nj <- cluster.size[cl]

            # data within for the cluster (centered by mu.y)
            Y1m <- Y1w.cm[cluster.idx == cl,, drop = FALSE]
            yc <- Y2.cm[cl,-Lp$between.idx[[2]]]
            zc <- Y2.cm[cl, Lp$between.idx[[2]]]

            # data between
            Y2Yc <- tcrossprod(Y2.cm[cl,])
            Y2Yc.zz <- Y2Yc[Lp$between.idx[[2]],
                        Lp$between.idx[[2]], drop = FALSE]
            Y2Yc.yz <- Y2Yc[-Lp$between.idx[[2]],
                             Lp$between.idx[[2]], drop = FALSE]
            Y2Yc.yy <- Y2Yc[-Lp$between.idx[[2]],
                            -Lp$between.idx[[2]], drop = FALSE]

            # construct sigma.j
            sigma.j <- (nj * sigma.b.z) + sigma.w
            sigma.j.inv <- lav_matrix_symmetric_inverse(S = sigma.j,
                               logdet = FALSE, Sinv.method = Sinv.method)
            sigma.ji.yz.zi <- sigma.j.inv %*% sigma.yz.zi
            sigma.zi.zy.ji <- t(sigma.ji.yz.zi)

            # common parts
            jYZj <- nj * (   sigma.j.inv %*%
                            (sigma.yz.zi %*% Y2Yc.zz %*% t(sigma.yz.zi)
                               - Y2Yc.yz %*% t(sigma.yz.zi)
                             - t(Y2Yc.yz %*% t(sigma.yz.zi)) + Y2Yc.yy)
                         %*% sigma.j.inv )

            Z1 <- Y2Yc.zz %*% t(sigma.ji.yz.zi) %*% sigma.yz
            YZ1 <- t(Y2Yc.yz) %*% sigma.j.inv %*% sigma.yz


            # Mu.Z
            G.muz[cl,] <- -2 * as.numeric(
                (sigma.zz.inv + nj*(sigma.zi.zy.ji %*% sigma.yz.zi)) %*% zc
                 -nj*sigma.zi.zy.ji %*% yc)

            # MU.Y
            G.muy[cl,] <- 2*nj * as.numeric(zc %*% sigma.zi.zy.ji -
                                            yc %*% sigma.j.inv)

            # SIGMA.W
            g.sigma.w <- ( (nj-1) * sigma.w.inv
                - sigma.w.inv %*% (crossprod(Y1m) - nj*Y2Yc.yy) %*% sigma.w.inv
                + sigma.j.inv - jYZj )

            tmp <- g.sigma.w*2; diag(tmp) <- diag(g.sigma.w)
            G.Sigma.w[cl,] <- lav_matrix_vech(tmp)

            # SIGMA.B
            g.sigma.b <- nj * (sigma.j.inv - jYZj)

            tmp <- g.sigma.b*2; diag(tmp) <- diag(g.sigma.b)
            G.Sigma.b[cl,] <- lav_matrix_vech(tmp)


            # SIGMA.ZZ
            g.sigma.zz <- ( sigma.zz.inv + nj * sigma.zz.inv %*% (
                              t(sigma.yz) %*% (sigma.j.inv - jYZj) %*% sigma.yz
                            -(1/nj * Y2Yc.zz + t(Z1) + Z1 - t(YZ1) - YZ1) ) %*%
                                                sigma.zz.inv )

            tmp <- g.sigma.zz*2; diag(tmp) <- diag(g.sigma.zz)
            G.Sigma.zz[cl,] <- lav_matrix_vech(tmp)

            # SIGMA.ZY
            g.sigma.yz <- 2 * nj * (
                          (sigma.j.inv %*%
                              (sigma.yz.zi %*% Y2Yc.zz - sigma.yz - Y2Yc.yz)
                               + jYZj %*% sigma.yz) %*% sigma.zz.inv )

            G.Sigma.yz[cl,] <- lav_matrix_vec(g.sigma.yz)
        }
    } # between.idx 

    else { # no level-2 variables

        for(cl in seq_len(nclusters)) {
            # cluster size
            nj <- cluster.size[cl]

            # data within for the cluster (centered by mu.y)
            Y1m <- Y1w.cm[cluster.idx == cl,, drop = FALSE]
            yc <- Y2.cm[cl,]

            # data between
            Y2Yc.yy <- tcrossprod(Y2.cm[cl,])

            # construct sigma.j
            sigma.j <- (nj * sigma.b) + sigma.w
            sigma.j.inv <- lav_matrix_symmetric_inverse(S = sigma.j,
                               logdet = FALSE, Sinv.method = Sinv.method)
            # common part 
            jYYj <- nj * sigma.j.inv %*% Y2Yc.yy %*% sigma.j.inv

            # MU.Y
            G.muy[cl,] <- -2*nj * as.numeric(yc %*% sigma.j.inv)

            # SIGMA.W
            g.sigma.w <- ( (nj-1) * sigma.w.inv
                - sigma.w.inv %*% (crossprod(Y1m) - nj*Y2Yc.yy) %*% sigma.w.inv
                + sigma.j.inv - jYYj )
            tmp <- g.sigma.w*2; diag(tmp) <- diag(g.sigma.w)
            G.Sigma.w[cl,] <- lav_matrix_vech(tmp)

            # SIGMA.B
            g.sigma.b <- nj * (sigma.j.inv - jYYj)
            tmp <- g.sigma.b*2; diag(tmp) <- diag(g.sigma.b)
            G.Sigma.b[cl,] <- lav_matrix_vech(tmp)
        }
    }

    # rearrange columns to Mu.W, Mu.B, Sigma.W, Sigma.B
    ov.idx  <- Lp$ov.idx
    p.tilde <- length( unique(c(ov.idx[[1]], ov.idx[[2]])) )

    # Mu.W
    Mu.W <- G.muy

    # Mu.B
    Mu.B.tilde <- matrix(0, nclusters, p.tilde)
    Mu.B.tilde[, ov.idx[[1]] ] <- G.muy
    if(length(between.idx) > 0L) {
        Mu.B.tilde[, between.idx ] <- G.muz
    }
    Mu.B <- Mu.B.tilde[, ov.idx[[2]], drop = FALSE]

    # Sigma.W
    Sigma.W <- G.Sigma.w

    # Sigma.B
    if(length(between.idx) > 0L) {
        p.tilde.star <- p.tilde*(p.tilde+1)/2
        B.tilde <- lav_matrix_vech_reverse(seq_len(p.tilde.star))

        Sigma.B.tilde <- matrix(0, nclusters, p.tilde.star)

        col.idx <- lav_matrix_vech( B.tilde[ ov.idx[[1]], ov.idx[[1]],
                                             drop = FALSE ] )
        Sigma.B.tilde[ , col.idx ] <- G.Sigma.b

        col.idx <- lav_matrix_vec( B.tilde[ ov.idx[[1]], between.idx,
                                            drop = FALSE ] )
        Sigma.B.tilde[ , col.idx ] <- G.Sigma.yz

        col.idx <- lav_matrix_vech( B.tilde[ between.idx, between.idx,
                                             drop = FALSE ] )
        Sigma.B.tilde[ , col.idx ] <- G.Sigma.zz

        col.idx <- lav_matrix_vech( B.tilde[ ov.idx[[2]], ov.idx[[2]],
                                             drop = FALSE ] )
        Sigma.B <- Sigma.B.tilde[ , col.idx, drop = FALSE ]
    } else {
        Sigma.B <- G.Sigma.b
    }

    SCORES <- cbind(Mu.W, Sigma.W, Mu.B, Sigma.B)

    SCORES
}

lav_mvnorm_cluster_information_firstorder <- function(Y1           = NULL,
                                                      YLp          = NULL,
                                                      Lp           = NULL,
                                                      Mu.W         = NULL,
                                                      Sigma.W      = NULL,
                                                      Mu.B         = NULL,
                                                      Sigma.B      = NULL,
                                                      divide.by.two = FALSE,
                                                      Sinv.method  = "eigen") {
    N <- NROW(Y1)

    SCORES <- lav_mvnorm_cluster_scores_2l(Y1           = Y1,
                                           YLp          = YLp,
                                           Lp           = Lp,
                                           Mu.W         = Mu.W,
                                           Sigma.W      = Sigma.W,
                                           Mu.B         = Mu.B,
                                           Sigma.B      = Sigma.B,
                                           Sinv.method  = Sinv.method)

    # divide by 2 (if we want scores wrt objective function)
    if(divide.by.two) {
        SCORES <- SCORES / 2
    }

    # unit information
    information <- crossprod(SCORES)/N

    information
}

# estimate ML estimates of Mu.W, Mu.B, Sigma.W, Sigma.B
# using the EM algorithm
#
# per cluster-SIZE
#
lav_mvnorm_cluster_em_sat <- function(YLp            = NULL,
                                      Lp             = NULL,
                                      verbose        = TRUE,
                                      tol            = 1e-04,
                                      max.iter       = 5000) {


    # sample stats    
    nclusters      <- Lp$nclusters[[2]]
    cluster.size   <- Lp$cluster.size[[2]]
    cluster.idx    <- Lp$cluster.idx[[2]]
    within.idx     <- Lp$within.idx[[2]]
    between.idx    <- Lp$between.idx[[2]]
    cluster.sizes  <- Lp$cluster.sizes[[2]]
    ncluster.sizes <- Lp$ncluster.sizes[[2]]
    n.s            <- Lp$cluster.size.ns[[2]]

    Y2 <- YLp[[2]]$Y2
    Y1Y1 <- YLp[[2]]$Y1Y1

    # take care of Mu.W (within only) and Mu.B (between only) estimates
    Mu.W.tilde <- YLp[[2]]$Mu.W
    Mu.B.tilde <- YLp[[2]]$Mu.B
    if(length(between.idx) > 0) {
        Mu.W <- Mu.W.tilde[-between.idx]
    } else {
        Mu.W <- Mu.W.tilde
    }
    if(length(within.idx) > 0) {
        Mu.B <- Mu.B.tilde[-within.idx]
    } else {
        Mu.B <- Mu.B.tilde
    }

    # starting values for Sigma
    ov.idx <- Lp$ov.idx
    Sigma.W <- diag( length(ov.idx[[1]]) )
    #if(length(between.idx) > 0) {
    #    Sigma.W <- YLp[[2]]$Sigma.W[-between.idx, -between.idx]
    #} else {
    #    Sigma.W <- YLp[[2]]$Sigma.W
    #}
    Sigma.B <- diag( length(ov.idx[[2]]) )

    # initial logl
    fx <- lav_mvnorm_cluster_loglik_samplestats_2l(YLp = YLp,
              Lp = Lp, Mu.W = Mu.W, Sigma.W = Sigma.W,
              Mu.B = Mu.B, Sigma.B = Sigma.B,
              Sinv.method = "eigen", log2pi = TRUE, minus.two = FALSE)

    # if verbose, report
    if(verbose) {
        cat("EM iter:", sprintf("%6d", 0),
            " fx = ", sprintf("%15.10f", fx),
            "\n")
    }

    # translate to internal matrices
    out <- lav_mvnorm_cluster_implied22l(Lp = Lp,
              Mu.W = Mu.W, Mu.B = Mu.B,
               Sigma.W = Sigma.W, Sigma.B = Sigma.B)
    mu.y <- out$mu.y; mu.z <- out$mu.z
    sigma.w <- out$sigma.w; sigma.b <- out$sigma.b
    sigma.zz <- out$sigma.zz; sigma.yz <- out$sigma.yz

    nvar.y <- ncol(sigma.w)
    nvar.z <- ncol(sigma.zz)

    # sigma.zz can be compute beforehand
    if(length(between.idx) > 0L) {
        Z <- Y2[,between.idx,drop=FALSE]
        zbar <- colMeans(Y2)[between.idx]
        sigma.zz <- 1/nclusters * crossprod(Z) - tcrossprod(zbar)

        Y1Y1 <- Y1Y1[-between.idx, -between.idx, drop=FALSE]
    }

    # EM iterations
    fx.old <- fx
    for(i in 1:max.iter) {

        # E-step

        if(length(between.idx) > 0L) {
            sigma.1 <- cbind(sigma.yz, sigma.b)
            mu <- c(mu.z, mu.y)
        } else {
            sigma.1 <- sigma.b
            mu <- mu.y
        }

        C2.W <- matrix(0, nvar.y, nvar.y)
        C4.W <- matrix(0, nvar.y, nvar.y)
        V.J  <- matrix(0, nvar.y, nvar.y)
        A.J  <- numeric(nvar.y)
        ZY.J <- matrix(0, nvar.z, nvar.y)
        for(clz in seq_len(ncluster.sizes)) {

            # cluster size
            nj <- cluster.sizes[clz]

            # tj
            if(length(between.idx) > 0L) {
                ybar.j <- Y2[cluster.size == nj, -between.idx, drop = FALSE]
                b.j <- cbind(Y2[cluster.size == nj, between.idx, drop = FALSE],
                             Y2[cluster.size == nj,-between.idx, drop = FALSE])
            } else {
                ybar.j <- b.j <- Y2[cluster.size == nj, , drop = FALSE]
           }

            # compute Sigma.j
            sigma.j <- sigma.w + nj * sigma.b
            if(length(between.idx) > 0L) {
                omega.j <- rbind( cbind(sigma.zz, t(sigma.yz)),
                                  cbind(sigma.yz, 1/nj * sigma.j) )
            } else {
                omega.j <- 1/nj * sigma.j
            }
            omega.j.inv <- solve(omega.j)
            sigma.1.j.inv <- sigma.1 %*% omega.j.inv

            # conditional expectation of 'v.j'
            b.jc <- t( t(b.j) - mu )
            tmp <- b.jc %*% t(sigma.1.j.inv)
            a.j <- t(t(tmp) + mu.y)

            # conditional expection of 'v.j %*% t(v.j)'
            V.j <- n.s[clz]*(sigma.b - (sigma.1.j.inv %*% t(sigma.1))) +
                           crossprod(a.j)

            # conditional expection of '(x - v.j) %*% t(x - v.j)'
            # contains 4 parts: xx, xv, vx, vv
            #C.j  <- ( (1/nj) * crossprod(Y1.j)   # xx
            #          - crossprod(ybar.j, a.j)  # xv
            #          - crossprod(a.j, ybar.j)  # vx
            #          + V.j )                   # vv

            #C1.j <-  (1/nj) * crossprod(Y1.j)   # xx
            C2.j <-  - crossprod(ybar.j, a.j)  # xv
            #C3.j <-  - crossprod(a.j, ybar.j)  # vx
            #C4.j <-  V.j                    # vv

            # sum over clusters
            V.J <- V.J + V.j
            A.J <- A.J + colSums(a.j)

            # between only
            if(length(between.idx) > 0L) {
                ZY.J <- ZY.J + crossprod(Y2[cluster.size == nj,
                                         between.idx,drop = FALSE], a.j)
            }

            # if Sigma.j is the same for all j
            C2.W <- C2.W + nj*C2.j
            C4.W <- C4.W + nj*V.j
        }

        # force symmetry for V.j 
        # this only seems needed because of t(sigma.1) above, if between.idx > 0
        V.J <- (V.J + t(V.J))/2

        abar <- 1/nclusters * A.J
        C.B  <- 1/nclusters * V.J

        C1.W <- Y1Y1
        C2.W <- C2.W
        C3.W <- t(C2.W)
        C4.W <- C4.W
        C.W  <- 1/sum(cluster.size) * (C1.W + C2.W + C3.W + C4.W)

        # between only
        if(length(between.idx) > 0L) {
            A <- 1/nclusters * ZY.J - tcrossprod(zbar, abar)
        }

        # M-step
        sigma.w <- C.W
        sigma.b <- C.B - tcrossprod(abar)
        mu.y    <- abar

        if(length(between.idx) > 0L) {
            sigma.yz <- t(A)
        }

        # back to model-implied dimensions, for fx
        implied <- lav_mvnorm_cluster_2l2implied(Lp = Lp,
                       sigma.w = sigma.w, sigma.b = sigma.b,
                       sigma.zz = sigma.zz, sigma.yz = sigma.yz,
                       mu.z = mu.z, mu.y = mu.y)
        fx <- lav_mvnorm_cluster_loglik_samplestats_2l(YLp = YLp,
              Lp = Lp, Mu.W = implied$Mu.W, Sigma.W = implied$Sigma.W,
              Mu.B = implied$Mu.B, Sigma.B = implied$Sigma.B,
              Sinv.method = "eigen", log2pi = TRUE, minus.two = FALSE)

        # diff
        fx.delta <- fx - fx.old

        if(fx.delta < 0) {
            warning("lavaan WARNING: last EM iteration failed to improve fx")
        }

        if(verbose) {
            cat("EM iter:", sprintf("%6d", i),
                " fx = ", sprintf("%15.10f", fx),
                #" fx2 = ", sprintf("%15.10f", fx2),
                " fx.delta = ", sprintf("%15.10f", fx.delta),
                "\n")
        }

        # check for convergence
        if(fx.delta < tol) {
            break
        } else {
            fx.old <- fx
        }

    } # EM iterations

    list(Sigma.W = implied$Sigma.W, Sigma.B = implied$Sigma.B,
         Mu.W = implied$Mu.W, Mu.B = implied$Mu.B, logl = fx)
}


# this is a temporary solution - 21 Dec 2017
# based on EM_LeePoon1998_v5.R
# final version will be different
lav_mvnorm_cluster_em_h0 <- function(YLp            = NULL,
                                 Lp             = NULL,
                                 Y1             = NULL,
                                 lavpartable    = NULL,
                                 h1             = NULL,
                                 ov.names.l     = NULL,
                                 verbose        = TRUE,
                                 verbose.x      = FALSE,
                                 tol            = 1e-04,
                                 max.iter       = 5000,
                                 mstep.iter.max = 10000L,
                                 mstep.rel.tol  = 1e-10) {
    # sample stats
    nclusters      <- Lp$nclusters[[2]]
    cluster.size   <- Lp$cluster.size[[2]]
    cluster.idx    <- Lp$cluster.idx[[2]]
    within.idx     <- Lp$within.idx[[2]]
    between.idx    <- Lp$between.idx[[2]]
    both.idx       <- Lp$both.idx[[2]]
    within.names   <- Lp$within.names[[2]]
    between.names  <- Lp$between.names[[2]]
#    cluster.sizes  <- Lp$cluster.sizes[[2]]
#    ncluster.sizes <- Lp$ncluster.sizes[[2]]
#    n.s            <- Lp$cluster.size.ns[[2]]

    nvar <- NCOL(Y1)
    nobs <- NROW(Y1)


    # data
    Y1Y1 <- YLp[[2]]$Y1Y1
    Y2   <- YLp[[2]]$Y2
    Y1.mean <- colMeans(Y1)

    #############################################################
    if(length(within.idx) > 0L && !any(lavpartable$exo == 1L)) {
        pt.within.int.idx <- which(lavpartable$op == "~1" &
                                   lavpartable$lhs %in% within.names)
        x.within.int.idx <- lavpartable$free[ pt.within.int.idx ]
    }
    if(length(between.idx) > 0L && !any(lavpartable$exo == 1L)) {
        pt.between.int.idx <- which(lavpartable$op == "~1" &
                                    lavpartable$lhs %in% between.names)
        x.between.int.idx <- lavpartable$free[ pt.between.int.idx ]
        pt.between.var.idx <- which(lavpartable$op == "~~" &
                                    lavpartable$lhs %in% between.names,
                                    lavpartable$rhs %in% between.names)
        # FIXME:
        # make sure that intercepts/(co)variances of exogenous (but free)
        # between-only variables are already OK in lavpartable$start!!
    }
    ############################################################

    # use h1 to start
    if(!is.null(h1)) {
        Mu.W <- h1$implied$mean[[1]]
        Mu.B <- h1$implied$mean[[2]]
        Sigma.W <- h1$implied$cov[[1]]
        Sigma.B <- h1$implied$cov[[2]]
    } else {
        # initial values
        ov.idx <- Lp$ov.idx
        Mu.W <- numeric( length(ov.idx[[1]]) )
        Sigma.W <- diag( length(ov.idx[[1]]) )
        Mu.B <- numeric( length(ov.idx[[2]]) )
        Sigma.B <- diag( length(ov.idx[[2]]) )
    }
    #x.current <- fit@optim$x
    x.current <- NULL
    # report initial fx
    fx <- lav_mvnorm_cluster_loglik_samplestats_2l(YLp = YLp,
              Lp = Lp, Mu.W = Mu.W, Sigma.W = Sigma.W,
              Mu.B = Mu.B, Sigma.B = Sigma.B,
              Sinv.method = "eigen", log2pi = TRUE, minus.two = FALSE)

    # if verbose, report
    if(verbose) {
        cat("EM iter:", sprintf("%4d", 0),
            " fx =", sprintf("%14.10f", fx),
            "\n")
    }

    ###############################################
    # override Mu.W if within-only x are included #
    # --> we fill them in at the very end         #
    ###############################################
    if(length(within.idx) > 0L) {
        Mu.W <- numeric( length(Mu.W) )
    }

    # translate to internal matrices
    out <- lav_mvnorm_cluster_implied22l(Lp = Lp,
              Mu.W = Mu.W, Mu.B = Mu.B,
               Sigma.W = Sigma.W, Sigma.B = Sigma.B)
    mu.y <- out$mu.y; mu.z <- out$mu.z
    sigma.w <- out$sigma.w; sigma.b <- out$sigma.b
    sigma.zz <- out$sigma.zz; sigma.yz <- out$sigma.yz

    nvar.y <- ncol(sigma.w)
    nvar.z <- ncol(sigma.zz)

    # mu.z and sigma.zz can be computed beforehand
    if(length(between.idx) > 0L) {
        Z <- Y2[, between.idx, drop = FALSE]
        mu.z <- colMeans(Y2)[between.idx]
        sigma.zz <- 1/nclusters * crossprod(Z) - tcrossprod(mu.z)

        Y1Y1 <- Y1Y1[-between.idx, -between.idx, drop=FALSE]
    }

    C.j <- vector("list", length = nclusters)
    D.j <- vector("list", length = nclusters)
    M.j <- matrix(0, nrow = nclusters, ncol = nvar.y)
    ZY.j <- vector("list", length = nclusters)

    # EM iterations
    fx.old <- fx
    for(i in 1:max.iter) {

        if(length(between.idx) > 0L) {
            sigma.1 <- cbind(sigma.yz, sigma.b)
            mu <- c(mu.z, mu.y)
        } else {
            sigma.1 <- sigma.b
            mu <- mu.y
        }

        # E-step
        for(cl in seq_len(nclusters)) {
            nj <- cluster.size[cl]

            # data
            if(length(between.idx) > 0L) {
                # z comes first!
                b.j    <- c(Y2[cl, between.idx],
                            Y2[cl,-between.idx])
                ybar.j <- Y2[cl,-between.idx]
                y1j <- Y1[cluster.idx == cl, -between.idx, drop = FALSE]
            } else {
                ybar.j <- b.j <- Y2[cl,]
                y1j <- Y1[cluster.idx == cl,,drop = FALSE]
            }
            sigma.j <- sigma.w + nj*sigma.b
            if(length(between.idx) > 0L) {
                omega.j <- rbind( cbind(sigma.zz, t(sigma.yz)),
                                  cbind(sigma.yz, 1/nj * sigma.j) )
            } else {
                omega.j <- 1/nj * sigma.j
            }
            omega.j.inv <- solve(omega.j)

            # E(v|y)
            Ev <- as.numeric(mu.y + (sigma.1 %*% omega.j.inv %*% (b.j - mu)))

            # Cov(v|y)
            Covv <- sigma.b - (sigma.1 %*% omega.j.inv %*% t(sigma.1))

            # force symmetry
            Covv <- (Covv + t(Covv))/2

            # E(vv|y) = Cov(v|y) + E(v|y)E(v|y)^T
            Evv <- Covv + tcrossprod(Ev)

            # C.j for this cluster
            tmp1 <- 1/nj * crossprod(y1j)  # yy
            tmp2 <- tcrossprod(ybar.j, Ev) # yv
            tmp3 <- t(tmp2)                # vy
            tmp4 <- Evv                    # vv
            c.j <- (tmp1 - tmp2 - tmp3 + tmp4)

            # D.g
            tmp1 <- Evv                   # vv
            tmp2 <- tcrossprod(Ev, mu.y)  # vy
            tmp3 <- t(tmp2)               # yv
            tmp4 <- tcrossprod(mu.y)      # yy
            d.j <- (tmp1 - tmp2 - tmp3 + tmp4)

            C.j[[cl]] <- nj * c.j # multiplied by nj!
          D.j[[cl]] <- d.j
            M.j[cl,] <- Ev
            # between only
            if(length(between.idx) > 0L) {
                ZY.j[[cl]] <- tcrossprod(Y2[cl,between.idx], Ev)
            }
        }
        M.b <- 1/nclusters * colSums(M.j)
        C.b <- 1/nclusters * Reduce("+", D.j)
        C.w <- 1/nobs * Reduce("+", C.j)

        # between only
        if(length(between.idx) > 0L) {
            A <- 1/nclusters * Reduce("+", ZY.j) - tcrossprod(mu.z, M.b)
        }

        # end of E-step

        # make symmetric (not needed here)
        C.b <- (C.b + t(C.b))/2
        C.w <- (C.w + t(C.w))/2

        # M-step
        sigma.w <- C.w
        sigma.b <- C.b
        mu.y    <- M.b # note mu.y[within.idx] = 0
        if(length(between.idx) > 0L) {
            sigma.yz <- t(A)
        }

        # back to model-implied dimensions
        implied <- lav_mvnorm_cluster_2l2implied(Lp = Lp,
                       sigma.w = sigma.w, sigma.b = sigma.b,
                       sigma.zz = sigma.zz, sigma.yz = sigma.yz,
                       mu.z = mu.z, mu.y = mu.y)
        rownames(implied$Sigma.W) <- colnames(implied$Sigma.W) <- ov.names.l[[1]]
        rownames(implied$Sigma.B) <- colnames(implied$Sigma.B) <- ov.names.l[[2]]

        # fit two-group model
        local.partable <- lavpartable
        level.idx <- which(names(local.partable) == "level")
        names(local.partable)[level.idx] <- "group"
        local.partable$est <- NULL
        local.partable$se  <- NULL
        #local.partable$start <- NULL

        # ensure 'fixed' values of exo covariates are in ustart
        exo.idx <- which(lavpartable$exo == 1L)
        if(length(exo.idx) > 0L) {
            local.partable$ustart[ exo.idx ] <- lavpartable$start[ exo.idx ]
        }

        # give current values as starting values
        free.idx <- which(lavpartable$free > 0L)
        if(!is.null(x.current)) {
            local.partable$ustart[ free.idx ] <- x.current
        } else {
            local.partable$ustart[ free.idx ] <- lavpartable$start[ free.idx ]

            # REMOVE within intercepts (fix to zero)!
            if(length(within.idx) > 0L && !any(lavpartable$exo == 1L)) {
                local.partable$ustart[ pt.within.int.idx ] <- 0
            }

            # between-only intercepts (mu.z) (needed?)
            if(length(between.idx) > 0L && !any(lavpartable$exo == 1L)) {
                local.partable$ustart[ pt.between.int.idx ] <- mu.z
            }

            # between-only (co)variances (sigma.zz) (needed?)
           # between-only (co)variances (sigma.zz) (needed?)

        }

        local.fit <- lavaan(local.partable,
                            sample.cov   = list(within  = implied$Sigma.W,
                                                between = implied$Sigma.B),
                            sample.mean  = list(within  = implied$Mu.W,
                                                between = implied$Mu.B),
                            sample.nobs  = Lp$nclusters,
                            sample.cov.rescale = FALSE,
                            #verbose = TRUE, 
                            control = list(iter.max = mstep.iter.max,
                                           rel.tol  = mstep.rel.tol),
                            fixed.x = any(lavpartable$exo == 1L),
                            estimator = "ML",
                            se = "none",
                            warn = FALSE,
                            #debug = TRUE,
                            test = "none")
        # end of M-step

        implied2 <- local.fit@implied
        fx <- lav_mvnorm_cluster_loglik_samplestats_2l(YLp = YLp,
          Lp = Lp, Mu.W = implied2$mean[[1]], Sigma.W = implied2$cov[[1]],
          Mu.B = implied2$mean[[2]], Sigma.B = implied2$cov[[2]],
          Sinv.method = "eigen", log2pi = TRUE, minus.two = FALSE)

        # diff
        fx.delta <- fx - fx.old

        if(verbose) {
            cat("EM iter:", sprintf("%4d", i),
                " fx =", sprintf("%17.10f", fx),
                " fx.delta =", sprintf("%9.8f", fx.delta),
                " mstep.iter =", sprintf("%3d",
                                  lavInspect(local.fit, "iterations")),
                "\n")
        }

        # convergence check: using parameter values only (for now)
        if(fx.delta < tol && i > 1L) {
            break
        } else {
            fx.old <- fx
            x.current <- local.fit@optim$x
            if(verbose.x) {
                print(round(x.current, 3))
            }
        }

        # translate to internal matrices
        out <- lav_mvnorm_cluster_implied22l(Lp = Lp,
                   Mu.W = implied2$mean[[1]], Sigma.W = implied2$cov[[1]],
                   Mu.B = implied2$mean[[2]], Sigma.B = implied2$cov[[2]])
        mu.y <- out$mu.y; sigma.w <- out$sigma.w; sigma.b <- out$sigma.b
        sigma.yz <- out$sigma.yz
        #mu.z <- out$mu.z
        #sigma.zz <- out$sigma.zz

        ####################################
        # remove within.idx part from mu.y #
        ####################################
        if(length(within.idx) > 0L) {
            w.idx <- match(within.names, ov.names.l[[1]])
            mu.y[w.idx] <- 0
        }

    } # EM iterations

    final.x <- local.fit@optim$x

    ##############################
    # add within.only intercepts #
    ##############################
    if(length(within.idx) > 0L && !any(lavpartable$exo == 1L)) {
        final.x[ x.within.int.idx ] <- Y1.mean[ within.idx ]
    }

    # add attributes
    attr(final.x, "iterations") <- i
    attr(final.x, "converged") <- TRUE # always true for EM...
    attr(final.x, "control") <- list(em.iter.max = max.iter,
                               em.tol = tol)
    attr(fx, "fx.group") <- fx # single group for now
    attr(final.x, "fx") <- fx

    final.x
}

