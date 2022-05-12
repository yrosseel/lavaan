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
    within.idx  <- Lp$within.idx[[2]]
    within.x.idx  <- Lp$within.x.idx[[2]]
    between.idx <- Lp$between.idx[[2]]
    between.y.idx <- Lp$between.y.idx[[2]]
    between.x.idx <- Lp$between.x.idx[[2]]

#    # ov.idx per level
#    ov.idx      <- Lp$ov.idx
#
#    # remove 'x' variables from ov.idx[[1]]
#    if(length(within.x.idx) > 0L) {
#        x.idx <- which(ov.idx[[1]] %in% within.x.idx)
#        ov.idx1.nox <- ov.idx[[1]][-x.idx]
#    } else {
#        ov.idx1.nox <- ov.idx[[1]]
#    }
#
#    # remove 'x' variables from ov.idx[[2]]
#    if(length(between.x.idx) > 0L) {
#        x.idx <- which(ov.idx[[2]] %in% between.x.idx)
#        ov.idx2.nox <- ov.idx[[2]][-x.idx]
#    } else {
#        ov.idx2.nox <- ov.idx[[2]]
#    }

    # only 'y'
    ov.idx <- Lp$ov.y.idx

    # two levels only (for now)
    ov.idx1.nox <- ov.idx[[1]]
    ov.idx2.nox <- ov.idx[[2]]

    # 'tilde' matrices: ALL variables within and between
    p.tilde <- length( unique(c(ov.idx[[1]], ov.idx[[2]])) )

    # Sigma.W.tilde
    Sigma.W.tilde <- matrix(0, p.tilde, p.tilde)
    Sigma.W.tilde[ ov.idx1.nox, ov.idx1.nox ] <- Res.Sigma.W

    # INT.W.tilde
    INT.W.tilde <- matrix(0, p.tilde, 1L)
    INT.W.tilde[ ov.idx1.nox, 1L ] <- Res.Int.W

    # PI.W.tilde
    PI.W.tilde <- matrix(0, p.tilde, ncol(Res.Pi.W))
    PI.W.tilde[ ov.idx1.nox, ] <- Res.Pi.W

    BETA.W.tilde <- rbind(t(INT.W.tilde), t(PI.W.tilde))



    # Sigma.B.tilde
    Sigma.B.tilde <- matrix(0, p.tilde, p.tilde)
    Sigma.B.tilde[ ov.idx2.nox, ov.idx2.nox ] <- Res.Sigma.B

    # INT.B.tilde
    INT.B.tilde <- matrix(0, p.tilde, 1L)
    INT.B.tilde[ ov.idx2.nox, 1L ] <- Res.Int.B

    # PI.B.tilde
    PI.B.tilde <- matrix(0, p.tilde, ncol(Res.Pi.B))
    PI.B.tilde[ ov.idx2.nox, ] <- Res.Pi.B

    BETA.B.tilde <- rbind(t(INT.B.tilde), t(PI.B.tilde))


#    if(length(between.idx) > 0L) {
#        beta.z <- BETA.B.tilde[, between.y.idx, drop = FALSE]
#        beta.b <- BETA.B.tilde[, -c(between.idx, within.x.idx), drop = FALSE]
#        beta.w <- BETA.W.tilde[, -c(between.idx, within.x.idx), drop = FALSE]
#        sigma.zz <- Sigma.B.tilde[ between.y.idx, between.y.idx, drop = FALSE]
#        sigma.yz <- Sigma.B.tilde[ -c(between.idx, within.x.idx),
#                                   between.y.idx, drop = FALSE]
#        sigma.b  <- Sigma.B.tilde[ -c(between.idx, within.x.idx),
#                                   -c(between.idx, within.x.idx), drop = FALSE]
#        sigma.w  <- Sigma.W.tilde[ -c(between.idx, within.x.idx),
#                                   -c(between.idx, within.x.idx), drop = FALSE]
#    } else {
#        beta.z <- matrix(0, 0L, 0L)
#        sigma.zz <- matrix(0, 0L, 0L)
#        if(length(within.x.idx) > 0L) {
#            beta.b <- BETA.B.tilde[, -within.x.idx, drop = FALSE]
#            beta.w <- BETA.W.tilde[, -within.x.idx, drop = FALSE]
#            sigma.b <- Sigma.B.tilde[-within.x.idx, -within.x.idx, drop = FALSE]
#            sigma.w <- Sigma.W.tilde[-within.x.idx, -within.x.idx, drop = FALSE]
#        } else { # should never happen; there is always either at least
#                 # 1 element in within.x.idx or between.idx (if conditional.x)
#            beta.b <- BETA.B.tilde
#            beta.w <- BETA.W.tilde
#            sigma.b <- Sigma.B.tilde
#            sigma.w <- Sigma.W.tilde
#        }
#        sigma.yz <- matrix(0, nrow(sigma.w), 0L)
#    }

    if(length(between.y.idx) > 0L) {
        beta.z   <- BETA.B.tilde[,  between.y.idx, drop = FALSE]
        beta.b   <- BETA.B.tilde[, -between.y.idx, drop = FALSE]
        beta.w   <- BETA.W.tilde[, -between.y.idx, drop = FALSE]
        sigma.zz <- Sigma.B.tilde[  between.y.idx, between.y.idx, drop = FALSE]
        sigma.yz <- Sigma.B.tilde[ -between.y.idx, between.y.idx, drop = FALSE]
        sigma.b  <- Sigma.B.tilde[ -between.y.idx,
                                   -between.y.idx, drop = FALSE]
        sigma.w  <- Sigma.W.tilde[ -between.y.idx,
                                   -between.y.idx, drop = FALSE]
    } else {
        beta.z   <- matrix(0, 0L, 0L)
        sigma.zz <- matrix(0, 0L, 0L)
        beta.b   <- BETA.B.tilde
        beta.w   <- BETA.W.tilde
        sigma.b  <- Sigma.B.tilde
        sigma.w  <- Sigma.W.tilde
        sigma.yz <- matrix(0, nrow(sigma.w), 0L)
    }


    # beta.wb # FIXme: not correct if some 'x' are splitted (overlap)
    # but because we ALWAYS treat splitted-x as 'y', this is not a problem
    beta.wb <- rbind(beta.w, beta.b[-1,,drop = FALSE])
    beta.wb[1,] <- beta.wb[1,,drop = FALSE] + beta.b[1,,drop = FALSE]

    list(sigma.w = sigma.w, sigma.b = sigma.b, sigma.zz = sigma.zz,
         sigma.yz = sigma.yz, beta.w = beta.w, beta.b = beta.b, beta.z = beta.z,
         beta.wb = beta.wb)
}




