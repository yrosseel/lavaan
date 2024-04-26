# check if a fitted model is admissible
lav_object_post_check <- function(object, verbose = FALSE) {
  stopifnot(inherits(object, "lavaan"))
  lavpartable <- object@ParTable
  lavmodel <- object@Model
  lavdata <- object@Data

  var.ov.ok <- var.lv.ok <- result.ok <- TRUE

  # 1a. check for negative variances ov
  var.idx <- which(lavpartable$op == "~~" &
    lavpartable$lhs %in% lavNames(object, "ov") &
    lavpartable$lhs == lavpartable$rhs)
  if (length(var.idx) > 0L && any(lavpartable$est[var.idx] < 0.0)) {
    result.ok <- var.ov.ok <- FALSE
    lav_msg_warn(gettext("some estimated ov variances are negative"))
  }

  # 1b. check for negative variances lv
  var.idx <- which(lavpartable$op == "~~" &
    lavpartable$lhs %in% lavNames(object, "lv") &
    lavpartable$lhs == lavpartable$rhs)
  if (length(var.idx) > 0L && any(lavpartable$est[var.idx] < 0.0)) {
    result.ok <- var.lv.ok <- FALSE
    lav_msg_warn(gettext("some estimated lv variances are negative"))
  }

  # 2. is cov.lv (PSI) positive definite? (only if we did not already warn
  # for negative variances)
  if (var.lv.ok && length(lavNames(lavpartable, type = "lv.regular")) > 0L) {
    ETA <- lavTech(object, "cov.lv")
    for (g in 1:lavdata@ngroups) {
      if (nrow(ETA[[g]]) == 0L) next
      txt.group <- if (lavdata@ngroups > 1L) gettextf("in group %s", g) else ""
      eigvals <- eigen(ETA[[g]], symmetric = TRUE, only.values = TRUE)$values
      if (any(eigvals < -1 * .Machine$double.eps^(3 / 4))) {
        lav_msg_warn(gettextf(
          "covariance matrix of latent variables is not positive definite %s;
          use lavInspect(fit, \"cov.lv\") to investigate.", txt.group
        ))
        result.ok <- FALSE
      }
    }
  }

  # 3. is THETA positive definite (but only for numeric variables)
  # and if we not already warned for negative ov variances
  if (var.ov.ok) {
    THETA <- lavTech(object, "theta")
    for (g in 1:lavdata@ngroups) {
      num.idx <- lavmodel@num.idx[[g]]
      if (length(num.idx) > 0L) {
        txt.group <- if (lavdata@ngroups > 1L) gettextf("in group %s", g) else ""
        eigvals <- eigen(THETA[[g]][num.idx, num.idx, drop = FALSE],
          symmetric = TRUE,
          only.values = TRUE
        )$values
        if (any(eigvals < -1 * .Machine$double.eps^(3 / 4))) {
          lav_msg_warn(gettextf(
            "the covariance matrix of the residuals of the observed variables
            (theta) is not positive definite %s; use lavInspect(fit, \"theta\")
            to investigate.", txt.group))
          result.ok <- FALSE
        }
      }
    }
  }

  result.ok
}
