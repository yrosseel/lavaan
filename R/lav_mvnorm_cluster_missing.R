# loglikelihood clustered/twolevel data in the presence of missing data

# YR: first version around March 2020



# Mu.W, Mu.B, Sigma.W, Sigma.B are the model-implied statistics
lav_mvnorm_cluster_missing_loglik_samplestats_2l <- function(Yp  = NULL,
                                                     Y2           = NULL,
                                                     Lp           = NULL,
                                                     Mp           = NULL,
                                                     Mu.W         = NULL,
                                                     Sigma.W      = NULL,
                                                     Mu.B         = NULL,
                                                     Sigma.B      = NULL,
                                                     Sinv.method  = "eigen",
                                                     log2pi       = FALSE,
                                                     loglik.x     = 0,
                                                     minus.two    = TRUE) {
    # map implied to 2l matrices
    out <- lav_mvnorm_cluster_implied22l(Lp = Lp, Mu.W = Mu.W, Mu.B = Mu.B,
               Sigma.W = Sigma.W, Sigma.B = Sigma.B)
    mu.y <- out$mu.y; mu.z <- out$mu.z
    sigma.w <- out$sigma.w; sigma.b <- out$sigma.b
    sigma.zz <- out$sigma.zz; sigma.yz <- out$sigma.yz

    # Lp
    nclusters       <- Lp$nclusters[[2]]
    between.idx     <- Lp$between.idx[[2]]

    # global
    sigma.w.inv <- solve.default(sigma.w)
    sigma.w.logdet <- log(det(sigma.w))

    # y
    ny <- ncol(sigma.w)
    ny.diag.idx <- lav_matrix_diag_idx(ny)

    # z
    if(length(between.idx) > 0L) {
        Z <- Y2[, between.idx, drop = FALSE]
        Z.c <- t( t(Z) - mu.z )

        sigma.zy <- t(sigma.yz)
        sigma.zz.inv <- solve.default(sigma.zz)
        sigma.zz.logdet <- log(det(sigma.zz))
        sigma.zi.zy <- sigma.zz.inv %*% sigma.zy
    }

    # Z per missing pattern
    if(length(between.idx) > 0L) {
        Zp  <- Mp$Zp
        SIGMA.B.Z <- vector("list", length = Zp$npatterns + 1L) # +1 for empty
        SIGMA.B.Z.cl <- integer(nclusters) # which sigma.b.z per cluster

        V.zz.logdet <- q.zz.a <- 0
        GJ <- matrix(0, nrow(Y2), ny)
        for(p in seq_len(Zp$npatterns)) {
            freq <- Zp$freq[p]; na.idx <- which(!Zp$pat[p,])
            j.idx <- Zp$case.idx[[p]] # cluster indices with this pattern
            SIGMA.B.Z.cl[j.idx] <- p

            if(length(na.idx) > 0L) {
                zp <- sigma.zz[-na.idx, -na.idx, drop = FALSE]
                zp.inv <- lav_matrix_symmetric_inverse_update(
                              S.inv = sigma.zz.inv, rm.idx = na.idx,
                              logdet = TRUE, S.logdet = sigma.zz.logdet)
                zp.logdet <- attr(zp.inv, "logdet")
                V.zz.logdet <- V.zz.logdet + (zp.logdet*freq)

                this.z <- Z.c[Zp$case.idx[[p]], -na.idx, drop = FALSE]
                zi.delta.p <- this.z %*% zp.inv
                q.zz.a <- q.zz.a + sum(zi.delta.p * this.z) # TR

                Z.G.ZY <- zp.inv %*% sigma.zy[-na.idx, , drop = FALSE]
                SIGMA.B.Z[[p]] <-
                    (sigma.b - sigma.yz[, -na.idx, drop = FALSE] %*% Z.G.ZY)
                GJ[j.idx, ] <- this.z %*% Z.G.ZY
            } else {
                # complete case
                V.zz.logdet <- V.zz.logdet + (sigma.zz.logdet * freq)

                this.z <- Z.c[Zp$case.idx[[p]], , drop = FALSE]
                zi.delta.p <- this.z %*% sigma.zz.inv
                q.zz.a <- q.zz.a + sum(zi.delta.p * this.z) # TR

                SIGMA.B.Z[[p]] <- ( sigma.b -
                                    sigma.yz %*% sigma.zz.inv %*% sigma.zy )
                GJ[j.idx, ] <- this.z %*% sigma.zi.zy
            }
        } # p

        # add empty patterns (if any)
        if(length(Zp$empty.idx) > 0L) {
            j.idx <- Zp$empty.idx
            SIGMA.B.Z.cl[j.idx] <- p + 1L
            SIGMA.B.Z[[ p + 1L ]] <- sigma.b
            GJ[j.idx,] <- numeric(ny) # must be zero vector
        }

    }


    # Y per missing pattern
    q.yy.a <- W.logdet <- 0
    PJ <- matrix(0, nrow = nclusters, ncol = ny)
    .A.j <- rep(list(matrix(0, ny, ny)), nclusters)
    for(p in seq_len(Mp$npatterns)) {
        freq <- Mp$freq[p]; na.idx <- which(!Mp$pat[p,])
        j.idx <- Mp$j.idx[[p]]; j1.idx <- Mp$j1.idx[[p]]
        TAB <- integer(nclusters); TAB[j1.idx] <- Mp$j.freq[[p]]

        if(length(na.idx) > 0L) {
            wp <- sigma.w[-na.idx, -na.idx, drop = FALSE]
            wp.inv <- lav_matrix_symmetric_inverse_update(
                          S.inv = sigma.w.inv, rm.idx = na.idx,
                          logdet = TRUE, S.logdet = sigma.w.logdet)
            wp.logdet <- attr(wp.inv, "logdet")
            W.logdet <- W.logdet + (wp.logdet * freq)

            TT <- Yp[[p]]$SY + tcrossprod(Yp[[p]]$MY - mu.y[-na.idx])
            q.yy.a <- q.yy.a + sum(wp.inv * TT) * Yp[[p]]$freq

            RS <- Yp[[p]]$ROWSUM
            MU <- matrix(mu.y[-na.idx], nrow(RS),
                         ncol = length(mu.y[-na.idx]), byrow = TRUE)
            tmp <- (RS - MU*Mp$j.freq[[p]]) %*% wp.inv
            PJ[j1.idx, -na.idx] <- ( PJ[j1.idx, -na.idx, drop = FALSE] + tmp )

            A.j <- matrix(0, ny, ny); A.j[-na.idx, -na.idx] <- wp.inv
            for(j in j1.idx) {
                .A.j[[j]] <- .A.j[[j]] + (A.j * TAB[j])
            }
        } else {
            # complete case
            W.logdet <- W.logdet + (sigma.w.logdet * freq)

            TT <- Yp[[p]]$SY + tcrossprod(Yp[[p]]$MY - mu.y)
            q.yy.a <- q.yy.a + sum(sigma.w.inv * TT) * Yp[[p]]$freq

            RS <- Yp[[p]]$ROWSUM
            MU <- matrix(mu.y, nrow(RS), ncol = length(mu.y), byrow = TRUE)
            tmp <- (RS - MU*Mp$j.freq[[p]]) %*% sigma.w.inv
            PJ[j1.idx, ] <- ( PJ[j1.idx, ,drop = FALSE] + tmp )

            for(j in j1.idx) {
                .A.j[[j]] <- .A.j[[j]] + (sigma.w.inv * TAB[j])
            }
        }
    } # p

    # per cluster
    q.zz.b <- q.zy <- q.yy.b <- IBZA.logdet <- 0
    for(cl in seq_len(nclusters)) {
        if(length(between.idx) > 0L) {
            g.j <- GJ[cl,]; p.j <- PJ[cl,]

            sigma.b.z <- SIGMA.B.Z[[ SIGMA.B.Z.cl[cl] ]]
            IBZA <- sigma.b.z %*% .A.j[[cl]]
            IBZA[ny.diag.idx] <- IBZA[ny.diag.idx] + 1

            IBZA.inv.g  <- solve.default(IBZA, g.j) # IBZA is NOT symmetric!
            IBZA.inv.BZ <- solve.default(IBZA, sigma.b.z)
            tmp <- determinant.matrix(IBZA, logarithm = TRUE)
            IBZA.logdet <- IBZA.logdet + tmp$modulus[1] * tmp$sign
            A.IBZA.inv.g  <- .A.j[[cl]] %*% IBZA.inv.g

            q.zz.b <- q.zz.b + sum(g.j * A.IBZA.inv.g)
            q.zy   <- q.zy   - sum(p.j * IBZA.inv.g)
            q.yy.b <- q.yy.b + sum(p.j * (IBZA.inv.BZ %*% p.j))
        } else {
           p.j <- PJ[cl,]
           IBZA <- sigma.b %*% .A.j[[cl]]
           IBZA[ny.diag.idx] <- IBZA[ny.diag.idx] + 1
           IBZA.inv.BZ <- solve.default(IBZA, sigma.b)
           tmp <- determinant.matrix(IBZA, logarithm = TRUE)
           IBZA.logdet <- IBZA.logdet + tmp$modulus[1] * tmp$sign

           q.yy.b <- q.yy.b + sum(p.j * (IBZA.inv.BZ %*% p.j))
        }
    }

    if(length(between.idx) > 0L) {
        P <- Mp$nel + Mp$Zp$nel
        DIST <- (q.yy.a - q.yy.b) + 2*q.zy + q.zz.a + q.zz.b
        LOGDET <- W.logdet + IBZA.logdet + V.zz.logdet
    } else {
        P <- Mp$nel
        DIST <- q.yy.a - q.yy.b
        LOGDET <- W.logdet + sum(IBZA.logdet)
    }

    if(log2pi && !minus.two) {
        LOG.2PI <- log(2 * pi)
        loglik <- -(P * LOG.2PI + LOGDET + DIST)/2
    } else {
        loglik <- DIST + LOGDET
    }

    # loglik.x (only if loglik is requested)
    if(length(unlist(Lp$ov.x.idx)) > 0L && log2pi && !minus.two) {
        loglik <- loglik - loglik.x
    }

    loglik
}


