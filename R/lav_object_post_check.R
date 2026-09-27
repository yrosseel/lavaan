# check if a fitted model is admissible
lav_object_post_check <- function(object) {
  stopifnot(inherits(object, "lavaan"))
  lavpartable <- object@ParTable
  lavmodel <- object@Model
  lavdata <- object@Data

  var_ov_ok <- var_lv_ok <- result_ok <- TRUE
  var_na <- FALSE

  # block labels (only used in the warning messages; note the leading space)
  nblocks <- lavmodel@nblocks
  block_txt <- rep("", nblocks)
  if (nblocks > 1L) {
    block_label <- lavdata@block.label
    if (length(block_label) != nblocks) {
      block_label <- as.character(seq_len(nblocks))
    }
    if (lavdata@nlevels > 1L) {
      block_txt <- gettextf(" in block %s", block_label)
    } else {
      block_txt <- gettextf(" in group %s", block_label)
    }
  }

  # 1a. check for negative variances ov
  var_idx <- which(lavpartable$op == "~~" &
    lavpartable$lhs %in% lav_object_vnames(object, "ov") &
    lavpartable$lhs == lavpartable$rhs)
  if (any(is.na(lavpartable$est[var_idx]))) {
    # perhaps estimator = "IV" + stage 1 only
    var_na <- TRUE
  } else if (length(var_idx) > 0L && any(lavpartable$est[var_idx] < 0.0)) {
    result_ok <- var_ov_ok <- FALSE
    lav_msg_warn(gettext("some estimated ov variances are negative"))
  }

  # 1b. check for negative variances lv
  var_idx <- which(lavpartable$op == "~~" &
    lavpartable$lhs %in% lav_object_vnames(object, "lv") &
    lavpartable$lhs == lavpartable$rhs)
  if (any(is.na(lavpartable$est[var_idx]))) {
    # perhaps estimator = "IV" + stage 1 only
    var_na <- TRUE
  } else if (length(var_idx) > 0L && any(lavpartable$est[var_idx] < 0.0)) {
    result_ok <- var_lv_ok <- FALSE
    lav_msg_warn(gettext("some estimated lv variances are negative"))
  }

  # 2. is cov.lv (VETA, regular latent variables only) positive definite?
  # (only if we did not already warn for negative variances)
  veta_ok <- rep(TRUE, nblocks)
  if (!var_na && var_lv_ok &&
      length(lav_object_vnames(lavpartable, type = "lv.regular")) > 0L) {
    eta <- lavTech(object, "cov.lv")
    for (b in seq_len(nblocks)) {
      if (nrow(eta[[b]]) == 0L) next
      if (!lav_object_post_check_psd(eta[[b]])) {
        lav_msg_warn(gettextf(
          "covariance matrix of latent variables is not positive definite%s;
          use lavInspect(fit, \"cov.lv\") to investigate.", block_txt[b]
        ))
        result_ok <- veta_ok[b] <- FALSE
      }
    }
  }

  # 3. is THETA positive definite (but only for numeric variables)
  # and if we have not already warned for negative ov variances
  theta_ok <- rep(TRUE, nblocks)
  if (!var_na && var_ov_ok) {
    mm_theta <- lavTech(object, "theta")
    for (b in seq_len(nblocks)) {
      num_idx <- lavmodel@num.idx[[b]]
      if (length(num_idx) > 0L) {
        theta <- mm_theta[[b]][num_idx, num_idx, drop = FALSE]
        if (!lav_object_post_check_psd(theta)) {
          lav_msg_warn(gettextf(
            "the covariance matrix of the residuals of the observed variables
            (theta) is not positive definite%s; use lavInspect(fit, \"theta\")
            to investigate.", block_txt[b]))
          result_ok <- theta_ok[b] <- FALSE
        }
      }
    }
  }

  # 4. is PSI positive definite? PSI is the covariance matrix of the
  # residuals in the structural part of the model. Unlike cov.lv, it also
  # contains the dummy latent variables that represent the observed
  # endogenous variables: their residual (co)variances live in PSI, not
  # in THETA. For numeric endogenous variables, lavTech(, "theta") merges
  # this block back into THETA, so check 3 covers them; but ordered
  # endogenous variables are excluded from check 3, and cov.lv drops all
  # dummy latent variables before its eigenvalue check (check 2). A block
  # of impossible residual covariances among ordered endogenous variables
  # therefore went unnoticed (issue #633).
  # If cov.lv is not positive definite, PSI cannot be either (same
  # inertia), so we only check PSI for blocks where checks 2 and 3 passed
  # (avoiding duplicate warnings).
  if (!var_na && var_ov_ok && var_lv_ok &&
      lavmodel@representation == "LISREL") {
    psi_idx <- which(names(lavmodel@GLIST) == "psi")
    for (b in seq_len(nblocks)) {
      if (!veta_ok[b] || !theta_ok[b] || b > length(psi_idx)) next
      psi <- lavmodel@GLIST[[psi_idx[b]]]
      if (nrow(psi) == 0L || anyNA(psi)) next
      if (!lav_object_post_check_psd(psi)) {
        lav_msg_warn(gettextf(
          "the covariance matrix of the residuals in the structural part of
          the model (psi) is not positive definite%s; use
          lavInspect(fit, \"est\")$psi to investigate.", block_txt[b]))
        result_ok <- FALSE
      }
    }
  }

  result_ok
}

# is a symmetric matrix positive semi-definite (numerically)?
# we tolerate a tiny negative eigenvalue: smaller (in absolute value) than
# .Machine$double.eps^(3/4), scaled by the largest eigenvalue (or 1 if the
# matrix is small), so that the tolerance follows the scale of the matrix
lav_object_post_check_psd <- function(mat) {
  eigvals <- eigen(mat, symmetric = TRUE, only.values = TRUE)$values
  tol <- .Machine$double.eps^(3 / 4) * max(1, abs(eigvals[1L]))
  !any(eigvals < -1 * tol)
}
