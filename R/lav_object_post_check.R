# check if a fitted model is admissible
lav_object_post_check <- function(object) {
  stopifnot(inherits(object, "lavaan"))
  lavpartable <- object@ParTable
  lavmodel <- object@Model
  lavdata <- object@Data

  var_ov_ok <- var_lv_ok <- result_ok <- TRUE
  var_na <- FALSE

  # 1a. check for negative variances ov
  var_idx <- which(lavpartable$op == "~~" &
    lavpartable$lhs %in% lav_partable_vnames(lavpartable, "ov") &
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
    lavpartable$lhs %in% lav_partable_vnames(lavpartable, "lv") &
    lavpartable$lhs == lavpartable$rhs)
  if (any(is.na(lavpartable$est[var_idx]))) {
    # perhaps estimator = "IV" + stage 1 only
    var_na <- TRUE
  } else if (length(var_idx) > 0L && any(lavpartable$est[var_idx] < 0.0)) {
    result_ok <- var_lv_ok <- FALSE
    lav_msg_warn(gettext("some estimated lv variances are negative"))
  }

  # 2. is cov.lv (PSI) positive definite? (only if we did not already warn
  # for negative variances)
  if (!var_na && var_lv_ok &&
      length(lav_object_vnames(lavpartable, type = "lv.regular")) > 0L) {
    eta <- lavTech(object, "cov.lv")
    for (g in 1:lavdata@ngroups) {
      if (nrow(eta[[g]]) == 0L) next
      txt_group <- if (lavdata@ngroups > 1L) gettextf("in group %s", g) else ""
      eigvals <- eigen(eta[[g]], symmetric = TRUE, only.values = TRUE)$values
      if (any(eigvals < -1 * .Machine$double.eps^(3 / 4))) {
        lav_msg_warn(gettextf(
          "covariance matrix of latent variables is not positive definite %s;
          use lavInspect(fit, \"cov.lv\") to investigate.", txt_group
        ))
        result_ok <- FALSE
      }
    }
  }

  # 3. is THETA positive definite (but only for numeric variables)
  # and if we have not already warned for negative ov variances
  if (!var_na && var_ov_ok) {
    mm_theta <- lavTech(object, "theta")
    for (g in 1:lavdata@ngroups) {
      num_idx <- lavmodel@num.idx[[g]]
      if (length(num_idx) > 0L) {
        txt_group <- ""
        if (lavdata@ngroups > 1L) txt_group <- gettextf("in group %s", g)
        eigvals <- eigen(mm_theta[[g]][num_idx, num_idx, drop = FALSE],
          symmetric = TRUE,
          only.values = TRUE
        )$values
        if (any(eigvals < -1 * .Machine$double.eps^(3 / 4))) {
          lav_msg_warn(gettextf(
            "the covariance matrix of the residuals of the observed variables
            (theta) is not positive definite %s; use lavInspect(fit, \"theta\")
            to investigate.", txt_group))
          result_ok <- FALSE
        }
      }
    }
  }

  result_ok
}
