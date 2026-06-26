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
lav_mat_rotate_gen <- function(m = 10L, orthogonal = TRUE) {
  # catch M=1
  if (m == 1L) {
    return(matrix(1, 1, 1))
  }

  # create random normal matrix
  tmp <- matrix(rnorm(m * m), nrow = m, ncol = m)

  if (orthogonal) {
    # use QR decomposition
    qr_out <- qr(tmp)
    q_1 <- qr.Q(qr_out)
    r <- qr.R(qr_out)
    # ... "normalized so that the diagonal elements of R are positive"
    sign_diag_r <- sign(diag(r))
    out <- q_1 * rep(sign_diag_r, each = m)
  } else {
    # just normalize *columns* of tmp -> crossprod(out) has 1 on diagonal
    out <- t(t(tmp) / sqrt(diag(crossprod(tmp))))
  }

  out
}

# check if ROT is an orthogonal matrix if orthogonal = TRUE, or normal if
# orthogonal = FALSE
lav_mat_rotate_check <- function(rot = NULL, orthogonal = TRUE,
                                    tolerance = sqrt(.Machine$double.eps)) {
  # we assume ROT is a matrix
  m <- nrow(rot)

  # crossprod
  rr <- crossprod(rot)

  # target
  if (orthogonal) {
    # ROT^T %*% ROT = I
    target <- diag(m)
  } else {
    # diagonal should be 1
    target <- rr
    diag(target) <- 1
  }

  # compare for near-equality
  res <- all.equal(target = target, current = rr, tolerance = tolerance)

  # return TRUE or FALSE
  if (is.logical(res) && res) {
    out <- TRUE
  } else {
    out <- FALSE
  }

  out
}

# get weights vector needed to weight the rows using Kaiser normalization
lav_mat_rotate_kaiser_weights <- function(a = NULL) {
  normalize <- 1 / sqrt(rowSums(a * a))
  idx_zero <- which(normalize == 0)
  # catch rows with all zero (thanks to Coen Bernaards for suggesting this)
  normalize[idx_zero] <- normalize[idx_zero] + .Machine$double.eps
  normalize
}

# get weights vector needed to weight the rows using Cureton & Mulaik (1975)
# standardization
# see also Browne (2001) page 128-129
#
# Note: the 'final' weights are multiplied by the Kaiser weights (see CEFA)
#
lav_mat_rotate_cm_weights <- function(a = NULL) {
  p <- nrow(a)
  m <- ncol(a)

  # first principal component of AA'
  a_eigen <- eigen(tcrossprod(a), symmetric = TRUE)
  a_1 <- a_eigen$vectors[, 1] * sqrt(a_eigen$values[1])

  kaiser_weights <- 1 / sqrt(rowSums(a * a))
  a_star <- abs(a_1 * kaiser_weights) # always between 0 and 1

  m_sqrt_inv <- 1 / sqrt(m)
  acos_m_sqrt_inv <- acos(m_sqrt_inv)

  delta <- numeric(p)
  delta[a_star < m_sqrt_inv] <- pi / 2

  tmp <- (acos_m_sqrt_inv - acos(a_star)) / (acos_m_sqrt_inv - delta) * (pi / 2)

  # add constant (see Cureton & Mulaik, 1975, page 187)
  cm <- cos(tmp) * cos(tmp) + 0.001

  # final weights = weighted by Kaiser weights
  cm * kaiser_weights
}

# build an orthogonal Givens (plane) rotation matrix: an m x m identity with a
# 2x2 rotation by 'theta' in the (col1, col2) plane.
lav_mat_givens_orth <- function(m = 2L, col1 = 1L, col2 = 2L, theta = 0) {
  rot <- diag(m)
  rot[col1, col1] <- base::cos(theta)
  rot[col1, col2] <- base::sin(theta)
  rot[col2, col1] <- -1 * base::sin(theta)
  rot[col2, col2] <- base::cos(theta)
  rot
}

# build an oblique (plane) rotation matrix for the (col1, col2) plane, given the
# shear 'delta' and the current factor correlation phi12. rot[col1, col1] equals
# the scaling factor 'gamma' (= sqrt(abs(gamma2))) used to update PHI.
lav_mat_givens_obliq <- function(m = 2L, col1 = 1L, col2 = 2L, delta = 0,
                                 phi12 = 0) {
  rot <- diag(m)
  gamma2 <- 1 + (2 * delta * phi12) + (delta * delta)
  rot[col1, col1] <- sqrt(abs(gamma2))
  rot[col1, col2] <- -1 * delta
  rot[col2, col1] <- 0
  rot[col2, col2] <- 1
  rot
}

# resolve the rotation criterion shared by lav_mat_rotate() and
# lav_mat_rotate_mg(): map the (already lower-cased) method name to its gradient
# function ('method_fname'), fill in the method-specific cf_gamma for the
# Crawford-Ferguson family, validate target/target_mask matrices (dimensions
# p x m), and promote a 'target.strict' with NA entries to 'pst'. 'group'
# selects the per-group matrix when target/target_mask are given as a list.
# Returns the (possibly updated) method, method_fname and method_args.
lav_mat_rotate_resolve_method <- function(method = NULL, method_args = list(),
                                          p = 0L, m = 0L, group = 1L) {
  # determine method function name
  if (method %in% c(
    "cf-quartimax", "cf-varimax", "cf-equamax",
    "cf-parsimax", "cf-facparsim"
  )) {
    method_fname <- "lav_mat_rotate_cf"
    method_args$cf_gamma <- switch(method,
      "cf-quartimax" = 0,
      "cf-varimax"   = 1 / p,
      "cf-equamax"   = m / (2 * p),
      "cf-parsimax"  = (m - 1) / (p + m - 2),
      "cf-facparsim" = 1
    )
  } else if (method %in% c("bi-quartimin", "biquartimin")) {
    method_fname <- "lav_mat_rotate_biquartimin"
  } else if (method %in% c("bi-geomin", "bigeomin")) {
    method_fname <- "lav_mat_rotate_bigeomin"
  } else if (method == "target.strict") {
    method_fname <- "lav_mat_rotate_target"
  } else {
    method_fname <- paste("lav_mat_rotate_", method, sep = "")
  }

  # check if rotation method exists
  check <- try(get(method_fname), silent = TRUE)
  if (inherits(check, "try-error")) {
    lav_msg_stop(gettext("unknown rotation method:"), method_fname)
  }

  # if target, check target matrix
  if (method == "target.strict" || method == "pst") {
    target <- method_args$target
    if (is.list(target)) {
      method_args$target <- target <- target[[group]]
    }
    # check dimension of target/A
    if (nrow(target) != p) {
      lav_msg_stop(gettext("nrow(target) != nrow(A)"))
    }
    if (ncol(target) != m) {
      lav_msg_stop(gettext("ncol(target) != ncol(A)"))
    }
  }
  if (method == "pst") {
    target_mask <- method_args$target_mask
    if (is.list(target_mask)) {
      method_args$target_mask <- target_mask <- target_mask[[group]]
    }
    # check dimension of target_mask/A
    if (nrow(target_mask) != p) {
      lav_msg_stop(gettext("nrow(target_mask) != nrow(A)"))
    }
    if (ncol(target_mask) != m) {
      lav_msg_stop(gettext("col(target_mask) != ncol(A)"))
    }
  }
  # we keep this here, so lav_mat_rotate() can be used independently
  if (method == "target.strict" && anyNA(target)) {
    method <- "pst"
    method_fname <- "lav_mat_rotate_pst"
    target_mask <- matrix(1, nrow = nrow(target), ncol = ncol(target))
    target_mask[is.na(target)] <- 0
    method_args$target_mask <- target_mask
  }

  list(method = method, method_fname = method_fname,
       method_args = method_args)
}

# taken from the stats package, but skipping varimax (already done):
lav_mat_rotate_promax <- function(x, m = 4, varimax_rot = NULL) {
  # this is based on promax() from factanal.R in /src/library/stats/R

  # 1. create 'ideal' pattern matrix
  q_1 <- x * abs(x) ^ (m - 1)

  # 2. regress x on Q to obtain 'rotation matrix' (same as 'procrustes')
  u <- lm.fit(x, q_1)$coefficients

  # 3. rescale so that solve(crossprod(U)) has 1 on the diagonal

  d <- diag(solve(t(u) %*% u))
  u <- u %*% diag(sqrt(d))
  dimnames(u) <- NULL

  # 4. create rotated factor matrix
  z <- x %*% u

  # 5. update rotation matrix
  u <- varimax_rot %*% u # here we plugin the rotation matrix from varimax

  list(loadings = z, rotmat = u)
}
