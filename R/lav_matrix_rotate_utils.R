# collection of functions that deal with rotation matrices
# YR  3 April 2019 -- initial version
# YR  6 Jan 2023: add promax
# YR 30 Jan 2024: orthogonal rotation matrix now reaches the full space

# generate random orthogonal rotation matrix
#
# reference for the orthogonal case:
#
# Stewart, G. W. (1980). The Efficient Generation of Random Orthogonal Matrices
# with an Application to Condition Estimators. SIAM Journal on Numerical
# Analysis, 17(3), 403-409. http://www.jstor.org/stable/2156882
#
lav_matrix_rotate_gen <- function(M = 10L, orthogonal = TRUE) {
  # catch M=1
  if (M == 1L) {
    return(matrix(1, 1, 1))
  }

  if (orthogonal) {
    # create random normal matrix
    tmp <- matrix(rnorm(M * M), nrow = M, ncol = M)
    # use QR decomposition
    qr.out <- qr(tmp)
    Q <- qr.Q(qr.out)
    R <- qr.R(qr.out)
    # ... "normalized so that the diagonal elements of R are positive"
    sign.diag.r <- sign(diag(R))
    out <- Q * rep(sign.diag.r, each = M)
  } else {
    # just normalize *columns* of tmp -> crossprod(out) has 1 on diagonal
    out <- t(t(tmp) / sqrt(diag(crossprod(tmp))))
  }

  out
}

# check if ROT is an orthogonal matrix if orthogonal = TRUE, or normal if
# orthogonal = FALSE
lav_matrix_rotate_check <- function(ROT = NULL, orthogonal = TRUE,
                                    tolerance = sqrt(.Machine$double.eps)) {
  # we assume ROT is a matrix
  M <- nrow(ROT)

  # crossprod
  RR <- crossprod(ROT)

  # target
  if (orthogonal) {
    # ROT^T %*% ROT = I
    target <- diag(M)
  } else {
    # diagonal should be 1
    target <- RR
    diag(target) <- 1
  }

  # compare for near-equality
  res <- all.equal(target = target, current = RR, tolerance = tolerance)

  # return TRUE or FALSE
  if (is.logical(res) && res) {
    out <- TRUE
  } else {
    out <- FALSE
  }

  out
}

# get weights vector needed to weight the rows using Kaiser normalization
lav_matrix_rotate_kaiser_weights <- function(A = NULL) {
  normalize <- 1 / sqrt(rowSums(A * A))
  idxZero <- which(normalize == 0)
  # catch rows with all zero (thanks to Coen Bernaards for suggesting this)
  normalize[idxZero] <- normalize[idxZero] + .Machine$double.eps
  normalize
}

# get weights vector needed to weight the rows using Cureton & Mulaik (1975)
# standardization
# see also Browne (2001) page 128-129
#
# Note: the 'final' weights are mutliplied by the Kaiser weights (see CEFA)
#
lav_matrix_rotate_cm_weights <- function(A = NULL) {
  P <- nrow(A)
  M <- ncol(A)

  # first principal component of AA'
  A.eigen <- eigen(tcrossprod(A), symmetric = TRUE)
  a <- A.eigen$vectors[, 1] * sqrt(A.eigen$values[1])

  Kaiser.weights <- 1 / sqrt(rowSums(A * A))
  a.star <- abs(a * Kaiser.weights) # always between 0 and 1

  m.sqrt.inv <- 1 / sqrt(M)
  acos.m.sqrt.inv <- acos(m.sqrt.inv)

  delta <- numeric(P)
  delta[a.star < m.sqrt.inv] <- pi / 2

  tmp <- (acos.m.sqrt.inv - acos(a.star)) / (acos.m.sqrt.inv - delta) * (pi / 2)

  # add constant (see Cureton & Mulaik, 1975, page 187)
  cm <- cos(tmp) * cos(tmp) + 0.001

  # final weights = weighted by Kaiser weights
  cm * Kaiser.weights
}

# taken from the stats package, but skipping varimax (already done):
lav_matrix_rotate_promax <- function(x, m = 4, varimax.ROT = NULL) {
  # this is based on promax() from factanal.R in /src/library/stats/R

  # 1. create 'ideal' pattern matrix
  Q <- x * abs(x)^(m - 1)

  # 2. regress x on Q to obtain 'rotation matrix' (same as 'procrustes')
  U <- lm.fit(x, Q)$coefficients

  # 3. rescale so that solve(crossprod(U)) has 1 on the diagonal

  d <- diag(solve(t(U) %*% U))
  U <- U %*% diag(sqrt(d))
  dimnames(U) <- NULL

  # 4. create rotated factor matrix
  z <- x %*% U

  # 5. update rotation amtrix
  U <- varimax.ROT %*% U # here we plugin the rotation matrix from varimax

  list(loadings = z, rotmat = U)
}
