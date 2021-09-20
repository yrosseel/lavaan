# loglikelihood clustered/twolevel data -- conditional.x = TRUE

# YR: first version around Sept 2021

# take model-implied mean+variance matrices, and reorder/augment them
# to facilitate computing of (log)likelihood in the two-level case

# when conditional.x = TRUE:
# - sigma.w and sigma.b: same dimensions, level-1 'Y' variables only
# - sigma.zz: level-2 variables only
# - sigma.yz: cov(level-1, level-2)
# - int.w: intercepts y within part
# - int.b: intercepts y between part
# - int.z: intercepts z (between-only)
# - slopes.w: slopes y within part
# - slopes.b: slopes y between part
# - slopes.z: slopes z (between-only)
# - beta.w: beta y within part
# - beta.b: beta y between part
# - beta.z: beta z (between-only)
lav_mvreg_cluster_implied22l <- function(Lp           = NULL,
                                         implied      = NULL,
                                         Res.Int.W    = NULL,
                                         Res.Int.B    = NULL,
                                         Res.Pi.W     = NULL,
                                         Res.Pi.B     = NULL,
                                         Res.Sigma.W  = NULL,
                                         Res.Sigma.B  = NULL) {
    if(!is.null(implied)) {
        # FIXME: only for single-group analysis!
        Res.Sigma.W <- implied$res.cov[[1]]
        Res.Int.W   <- implied$res.int[[1]]
        Res.Pi.W    <- implied$res.slopes[[1]]

        Res.Sigma.B <- implied$res.cov[[2]]
        Res.Int.B   <- implied$res.int[[2]]
        Res.Pi.B    <- implied$res.slopes[[2]]
    }

    # within/between idx
    within.y.idx  <- Lp$within.y.idx[[2]]
    within.x.idx  <- Lp$within.x.idx[[2]]
    between.y.idx <- Lp$between.y.idx[[2]]
    between.x.idx <- Lp$between.x.idx[[2]]
    both.idx      <- Lp$both.idx[[2]]
    # NOTE: these indices assume ALL ovs are present at both levels...

    # indices for level-2 ov's only!
    Z.IDX <- which(Lp$ov.y.idx[[2]] %in% Lp$between.idx[[2]])


    # within y-part
    sigma.w <- Res.Sigma.W
    int.w   <- Res.Int.W
    pi.w    <- Res.Pi.W
    beta.w  <- rbind(matrix(int.w, nrow = 1), t(pi.w))

    # between y-part (no-z)
    if(length(Z.IDX) > 0L) {
        SIGMA.B <- implied$res.cov[[2]][-Z.IDX, -Z.IDX, drop = FALSE]
        INT.B   <- implied$res.int[[2]][-Z.IDX, , drop = FALSE]
        PI.B    <- implied$res.slopes[[2]][-Z.IDX, , drop = FALSE]
    } else {
        SIGMA.B <- implied$res.cov[[2]]
        INT.B   <- implied$res.int[[2]]
        PI.B    <- implied$res.slopes[[2]]
    }
    sigma.b <- matrix(0, nrow(sigma.w), ncol(sigma.w))
    sigma.b[both.idx, both.idx] <- SIGMA.B

    int.b   <- matrix(0, nrow(sigma.w), 1L)
    int.b[both.idx, 1L] <- INT.B

    pi.b    <- matrix(0, nrow(sigma.w), ncol(PI.B))
    pi.b[both.idx, ] <- PI.B

    beta.b  <- rbind(matrix(int.b, nrow = 1L), t(pi.b))

    if(length(Z.IDX) > 0L) {
        sigma.zz <- implied$res.cov[[2]][Z.IDX, Z.IDX, drop = FALSE]
        int.z    <- implied$res.int[[2]][Z.IDX, , drop = FALSE]
        pi.z     <- implied$res.slopes[[2]][Z.IDX, , drop = FALSE]
        beta.z   <- rbind(matrix(int.z, nrow = 1), t(pi.z))

        SIGMA.YZ <- implied$res.cov[[2]][-Z.IDX, Z.IDX, drop = FALSE]
        sigma.yz <- matrix(0, nrow(sigma.w), ncol(SIGMA.YZ))
        sigma.yz[both.idx,] <- SIGMA.YZ
    } else {
        sigma.zz <- matrix(0, 0L, 0L)
        sigma.yz <- matrix(0, 0L, 0L)
        int.z <- numeric(0L)
        pi.z <- matrix(0, 0L, 0L)
        beta.z <- matrix(0, 0L, 0L)

        sigma.yz <- matrix(0, nrow(sigma.b), 0L)
    }

    list(sigma.w = sigma.w, sigma.b = sigma.b, sigma.zz = sigma.zz,
         sigma.yz = sigma.yz,
         int.w = int.w, pi.w = pi.w, beta.w = beta.w,
         int.b = int.b, pi.b = pi.b, beta.b = beta.b,
         int.z = int.z, pi.z = pi.z, beta.z = beta.z)
}




