# various rotation criteria and their gradients
# initial version YR 5 April 2019

# references:
#
# Bernaards, C. A., & Jennrich, R. I. (2005). Gradient projection algorithms
# and software for arbitrary rotation criteria in factor analysis. Educational
# and Psychological Measurement, 65(5), 676-696.
# website: http://www.stat.ucla.edu/research/gpa/splusfunctions.net
#
# Browne, M. W. (2001). An overview of analytic rotation in exploratory factor
# analysis. Multivariate behavioral research, 36(1), 111-150.
#
# Mulaik, S. A. (2010). Foundations of factor analysis (Second Edition).
# Boca Raton: Chapman and Hall/CRC.

# Note: this is YR's implementation, not a copy of the GPArotation
#       package
#
# Why did I write my own functions (and not use the GPArotation):
# - to better understand what is going on
# - to have direct access to the gradient functions
# - to avoid yet another dependency
# - to simplify further experiments


# Orthomax family (Harman, 1960)
#
# gamma = 0   -> quartimax
# gamma = 1/2 -> biquartimax
# gamma = 1/P -> equamax
# gamma = 1   -> varimax
#
lav_matrix_rotate_orthomax <- function(LAMBDA = NULL, gamma = 1, grad = FALSE) {
    L2 <- LAMBDA * LAMBDA
    # center L2 column-wise
    cL2 <- t( t(L2) - gamma * colMeans(L2) )
    out <- -1 * sum(L2 * cL2)/4

    if(grad) {
        attr(out, "grad") <- -1 * LAMBDA * cL2
    }

    out
}

# Crawford-Ferguson (1970) family
#
# combine penalization for 1) row complexity, and 2) column complexity
# if combined with orthogonal rotation, this is equivalent to the
# orthomax family:
#
#     quartimax -> gamma = 0 (only row complexity)
#     varimax   -> gamma = 1/nrow
#     equamax   -> gamma = ncol/(2*nrow)
#     parsimax  -> gamma = (ncol - 1)/(nrow + ncol - 2)
#     factor parsimony -> gamma = 1 (only column complexity)
#
# the Crawford-Ferguson family is also equivalent to the oblimin family
# if the latter is restricted to orthogonal rotation
#
lav_matrix_rotate_cf <- function(LAMBDA = NULL, gamma = 0, grad = FALSE) {
    # check if gamma is between 0 and 1?
    nRow <- nrow(LAMBDA)
    nCol <- ncol(LAMBDA)
    ROW1 <- matrix(1.0, nCol, nCol); diag(ROW1) <- 0.0
    COL1 <- matrix(1.0, nRow, nRow); diag(COL1) <- 0.0

    L2 <- LAMBDA * LAMBDA
    LR <- L2 %*% ROW1
    LC <- COL1 %*% L2

    f1 <- sum(L2 * LR)/4
    f2 <- sum(L2 * LC)/4

    out <- (1 - gamma)*f1 + gamma*f2

    if(grad) {
        attr(out, "grad") <- ((1 - gamma) * LAMBDA * LR) + (gamma * LAMBDA * LC)
    }

    out
}

# Oblimin family (Carroll, 1960; Harman, 1976)
#
# quartimin   -> gamma = 0
# biquartimin -> gamma = 1/2
# covarimin   -> gamma = 1
#
# if combined with orthogonal rotation, this is equivalent to the
# orthomax family (they have the same optimizers):
#
# gamma = 0   -> quartimax
# gamma = 1/2 -> biquartimax
# gamma = 1   -> varimax
# gamma = P/2 -> equamax
#
lav_matrix_rotate_oblimin <- function(LAMBDA = NULL, gamma = 0,
                                           grad = FALSE) {
    nRow <- nrow(LAMBDA)
    nCol <- ncol(LAMBDA)
    ROW1 <- matrix(1.0, nCol, nCol); diag(ROW1) <- 0.0

    L2 <- LAMBDA * LAMBDA
    LR <- L2 %*% ROW1
    Jp <- matrix(1, nRow, nRow)/nRow

    # see Jennrich (2002, p. 11)
    tmp <- (diag(nRow) - gamma * Jp) %*% LR

    # same as t( t(L2) - gamma * colMeans(L2) ) %*% ROW1

    out <- sum(L2 * tmp)/4

    if(grad) {
        attr(out, "grad") <- LAMBDA * tmp
    }

    out
}


# quartimax criterion
# Carroll (1953); Saunders (1953) Neuhaus & Wrigley (1954); Ferguson (1954)
# we use here the equivalent 'Ferguson, 1954' variant
# (See Mulaik 2010, p. 303)
lav_matrix_rotate_quartimax <- function(LAMBDA = NULL, grad = FALSE) {
    L2 <- LAMBDA * LAMBDA
    out <- -1 * sum(L2 * L2)/4

    if(grad) {
        attr(out, "grad") <- -1 * LAMBDA * L2
    }

    out
}


# varimax criterion
# Kaiser (1958, 1959)
#
# special case of the Orthomax family (Harman, 1960), where gamma = 1
# see Jennrich (2001, p. 296)
lav_matrix_rotate_varimax <- function(LAMBDA = NULL, grad = FALSE) {
    L2 <- LAMBDA * LAMBDA
    # center L2 column-wise
    cL2 <- t( t(L2) - colMeans(L2) )
    out <- -1 * abs(sum(L2 * cL2))/4  # abs needed?

    if(grad) {
        attr(out, "grad") <- -1 * LAMBDA * cL2
    }

    out
}

# quartimin criterion (part of Carroll's oblimin family
lav_matrix_rotate_quartimin <- function(LAMBDA = NULL, grad = FALSE) {
    nCol <- ncol(LAMBDA)
    ROW1 <- matrix(1.0, nCol, nCol); diag(ROW1) <- 0.0

    L2 <- LAMBDA * LAMBDA
    LR <- L2 %*% ROW1

    out <- sum(L2 * LR)/4

    if(grad) {
        attr(out, "grad") <- LAMBDA * LR
    }

    out
}

# Browne's (2001) version of Yates (1984) geomin criterion
#
# we use the exp/log trick as in Bernaard & Jennrich (2005, p. 687)
lav_matrix_rotate_geomin <- function(LAMBDA = NULL, epsilon = 0.01,
                                          grad = FALSE) {
    nCol <- ncol(LAMBDA)

    L2 <- LAMBDA * LAMBDA
    L2 <- L2 + epsilon

    if(epsilon < sqrt(.Machine$double.eps)) {
        # Yates's original formula
        tmp <- apply(L2, 1, prod)^(1/nCol)
    } else {
        tmp <- exp( rowSums(log(L2)) / nCol )
    }

    out <- sum(tmp)

    if(grad) {
        attr(out, "grad") <- (2/nCol) * LAMBDA/L2 * tmp
    }

    out
}



# gradient check
ilav_matrix_rotate_grad_test <- function(crit = NULL, ...,
                                         LAMBDA = NULL,
                                         nRow = 20L, nCol = 5L) {
    # test matrix
    if(is.null(LAMBDA)) {
        LAMBDA <- matrix(rnorm(nRow*nCol), nRow, nCol)
    }

    ff <- function(x, ...) {
        Lambda <- matrix(x, nRow, nCol)
        crit(Lambda, ..., grad = FALSE)
    }

    GQ1 <- matrix(numDeriv::grad(func = ff, x = as.vector(LAMBDA), ...),
                  nRow, nCol)
    GQ2 <- attr(crit(LAMBDA, ..., grad = TRUE), "grad")

    all.equal(GQ1, GQ2)
}

ilav_matrix_rotate_grad_test_all <- function() {

    # Orthomax family with various values for gamma
    for(gamma in seq(0,1,0.2)) {
        check <- ilav_matrix_rotate_grad_test(crit = lav_matrix_rotate_orthomax,
                                              gamma = gamma)
        if(is.logical(check) && check) {
            cat("orthomax + gamma = ", sprintf("%3.1f", gamma),
                ": OK\n")
        } else {
            cat("orthomax + gamma = ", sprintf("%3.1f", gamma),
                ": FAILED\n")
        }
    }

    # Crawford-Ferguson with various values for gamma
    for(gamma in seq(0,1,0.2)) {
        check <- ilav_matrix_rotate_grad_test(crit = lav_matrix_rotate_cf,
                                              gamma = gamma)
        if(is.logical(check) && check) {
            cat("Crawford-Ferguson + gamma = ", sprintf("%3.1f", gamma),
                ": OK\n")
        } else {
            cat("Crawford-Ferguson + gamma = ", sprintf("%3.1f", gamma),
                ": FAILED\n")
        }
    }

    # Oblimin family with various values for gamma
    for(gamma in seq(0,1,0.2)) {
        check <- ilav_matrix_rotate_grad_test(crit = lav_matrix_rotate_oblimin,
                                              gamma = gamma)
        if(is.logical(check) && check) {
            cat("Oblimin + gamma = ", sprintf("%3.1f", gamma),
                ": OK\n")
        } else {
            cat("Oblimin + gamma = ", sprintf("%3.1f", gamma),
                ": FAILED\n")
        }
    }

    # quartimax
    check <- ilav_matrix_rotate_grad_test(crit = lav_matrix_rotate_quartimax)
    if(is.logical(check) && check) {
        cat("quartimax: OK\n")
    } else {
        cat("quartimax: FAILED\n")
    }

    # varimax
    check <- ilav_matrix_rotate_grad_test(crit = lav_matrix_rotate_varimax)
    if(is.logical(check) && check) {
        cat("varimax: OK\n")
    } else {
        cat("varimax: FAILED\n")
    }

    # quartimin
    check <- ilav_matrix_rotate_grad_test(crit = lav_matrix_rotate_quartimin)
    if(is.logical(check) && check) {
        cat("quartimin: OK\n")
    } else {
        cat("quartimin: FAILED\n")
    }

    # geomin
    check <- ilav_matrix_rotate_grad_test(crit = lav_matrix_rotate_geomin)
    if(is.logical(check) && check) {
        cat("geomin: OK\n")
    } else {
        cat("geomin: FAILED\n")
    }
}

