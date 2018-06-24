# special functions for the one-factor model 
# YR 24 June 2018

# 1-factor model with (only) three indicators:
# no iterations needed; can be solved analytically

# denote s11, s22, s33 the diagonal elements, and  
# s21, s31, s32 the off-diagonal elements

# under the 1-factor model; typically, either psi == 1, or l1 == 1
# - s11 == l1^2*psi + theta1
# - s22 == l2^2*psi + theta2
# - s33 == l3^2*psi + theta3
# - s21 == l2*l1*psi
# - s31 == l3*l1*psi
# - s32 == l3*l2*psi
# 6 unknowns, 6 knowns

# note: if the triad of covariances is negative, there is no
# `valid' solution, for example:
#
# > S
#      [,1] [,2] [,3]
# [1,]  1.0  0.6  0.3
# [2,]  0.6  1.0 -0.1
# [3,]  0.3 -0.1  1.0
#
# (note: all eigenvalues are positive)

lav_cfa_1fac_3ind <- function(sample.cov, std.lv = FALSE, 
                              warn.neg.triad = TRUE) {
    
    # check sample cov 
    stopifnot(is.matrix(sample.cov))
    nRow <- NROW(sample.cov); nCol <- NCOL(sample.cov)
    stopifnot(nRow == 3L, nCol == 3L)

    s11 <- sample.cov[1,1]; s22 <- sample.cov[2,2]; s33 <- sample.cov[3,3]
    stopifnot(s11 > 0, s22 > 0, s33 > 0)

    s21 <- sample.cov[2,1]; s31 <- sample.cov[3,1]; s32 <- sample.cov[3,2]
    # note: s21*s31*s32 should be positive!
    if(s21 * s31 * s32 < 0 && warn.neg.triad) {
        warning("lavaan WARNING: product of the three covariances is negative!")
    }

    # solution
    if(std.lv) {
        psi <- 1
        l1.square <- (s21*s31)/s32
        l2.square <- (s21*s32)/s31
        l3.square <- (s31*s32)/s21
        theta1 <- s11 - l1.square
        theta2 <- s22 - l2.square
        theta3 <- s33 - l3.square
        l1 <- sign(l1.square) * sqrt( abs(l1.square) )
        l2 <- sign(l2.square) * sqrt( abs(l2.square) )
        l3 <- sign(l3.square) * sqrt( abs(l3.square) )
    } else {
        psi <- (s21*s31)/s32
        l1 <- 1
        l2 <- s32/s31 # l2 <- s21/psi
        l3 <- s32/s21 # l3 <- s31/psi
        theta1 <- s11 - psi
        theta2 <- s22 - l2*l2*psi
        theta3 <- s33 - l3*l3*psi
    }
    
    list(lambda = c(l1,l2,l3), psi = psi, theta = c(theta1, theta2, theta3))
}

# FABIN (Hagglund, 1982)
# 1-factor only (in this case fabin2 == fabin3)
lav_cfa_1fac_fabin <- function(S) {

    nvar <- NCOL(S)
    if(nvar < 3) {
        return( rep(1, nvar) )
    }

    out <- numeric( nvar ); out[1L] <- 1.0
    for(i in 2:nvar) {
        idx3 <- (1:nvar)[-c(i, 1L)]
        s23 <- S[i, idx3]
        S31 <- S13 <- S[idx3, 1L]
        out[i] <- ( s23 / S31 )
    }

    out
}

