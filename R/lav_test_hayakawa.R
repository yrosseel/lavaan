# Corrected (mean and variance) adjusted test statistics
#
# Hayakawa, K. (2018). Corrected goodness-of-fit test in covariance
# structure analysis. Psychological Methods, 24(3), 371-389.
#
# The mean-and-variance adjusted test statistic (test = "mean.var.adjusted")
# divides the test statistic by tr(UGamma^2)/tr(UGamma), with degrees of
# freedom tr(UGamma)^2/tr(UGamma^2). But the naive (plug-in) estimator of
# tr(UGamma^2) is severely (upward) biased unless p* = p(p+1)/2 is small
# relative to N; as a result, the adjusted test underrejects (often towards
# zero) when the number of variables grows. Hayakawa (2018) replaces the
# naive estimator by the unbiased estimator of a trace of a squared
# covariance matrix of Srivastava (2005) and Himeno & Yamada (2014), which
# remains accurate in high dimensions:
#
#   a2c = [ (N-2)(N-1) tr(H^2) - N(N-1) tr(D^2) + (tr H)^2 ] /
#         [ N(N-1)(N-2)(N-3) ]
#
# where H = sum_i y_i y_i', D = diag(y_1'y_1, ..., y_N'y_N),
# y_i = U^{1/2} (s_i - sbar), and s_i is the casewise (saturated) moment
# vector of observation i.
#
# Note that U^{1/2} is never needed: with Y = [y_1, ..., y_N] and
# Zc the N x p~ matrix with rows (s_i - sbar)', the N x N Gram matrix
# M = Zc U Zc' = Y'Y satisfies tr(H) = tr(M), tr(H^2) = tr(M^2), and
# y_i'y_i = diag(M), so all three ingredients come from M directly.
#
# The corrected tests are exposed as test = "mean.var.adjusted.corrected"
# and (a natural extension, not in the paper) "scaled.shifted.corrected";
# they combine with scaled.test = "browne.residual.nt.model" to give the
# RLS-based versions recommended by Hayakawa (2018).

# eligibility: the corrected trace estimator assumes iid casewise moment
# vectors from a single sample, computed from complete continuous data
lav_test_hayakawa_check <- function(lavoptions = NULL, lavdata = NULL,
                                    lavmodel = NULL) {
  context <- gettext("corrected adjusted tests")

  if (!is.null(lavdata)) {
    if (lavdata@ngroups > 1L) {
      lav_msg_stop(gettextf(
        "%s are only available for single-group models.", context))
    }
    if (lavdata@nlevels > 1L) {
      lav_msg_stop(gettextf(
        "%s are only available for single-level models.", context))
    }
    if (length(lavdata@X) == 0L || is.null(lavdata@X[[1]])) {
      lav_msg_stop(gettextf(
        "%s require the raw data.", context))
    }
  }

  if (!is.null(lavmodel) && lavmodel@categorical) {
    lav_msg_stop(gettextf(
      "%s are only available for continuous data.", context))
  }

  if (!is.null(lavoptions)) {
    if (!isTRUE(lavoptions$estimator == "ML")) {
      lav_msg_stop(gettextf(
        "%1$s require estimator = %2$s; found %3$s.",
        context, dQuote("ML"), dQuote(lavoptions$estimator)))
    }
    if (!isTRUE(lavoptions$missing[1] == "listwise")) {
      lav_msg_stop(gettextf(
        "%1$s require missing = %2$s; found %3$s.",
        context, dQuote("listwise"), dQuote(lavoptions$missing[1])))
    }
    if (isTRUE(lavoptions$conditional.x)) {
      lav_msg_stop(gettextf(
        "%s are not available if conditional.x = TRUE.", context))
    }
  }

  invisible(TRUE)
}

# unbiased estimator of tr((U Omega)^2), Omega = Var(s_i)
# (Srivastava 2005; Himeno & Yamada 2014; Hayakawa 2018, eq. 19)
#
# u:      the U matrix (as in UfromUGamma), single group
# m_y:    raw data matrix (N x p), complete cases
# meanstructure: if TRUE, s_i = c(y_i, vech((y_i - ybar)(y_i - ybar)'));
#                if FALSE, only the vech part
#
# returns a list with both traces computed from the same casewise Gram
# matrix, exactly as in Hayakawa (2018):
#  - trace.UGamma  = tr(U Omega.hat) with the (unbiased) N-1 divisor
#  - trace.UGamma2 = the corrected (unbiased) estimator of tr((U Omega)^2)
lav_test_hayakawa_trace2 <- function(u = NULL, m_y = NULL,
                                     meanstructure = FALSE) {
  m_y <- unname(as.matrix(m_y))
  n <- nrow(m_y)
  p <- ncol(m_y)

  if (n < 4L) {
    lav_msg_stop(gettext(
      "the corrected trace estimator requires at least 4 observations."))
  }

  # casewise (saturated) moment vectors s_i
  m_yc <- t(t(m_y) - colMeans(m_y))
  idx1 <- lav_mat_vech_col_idx(p)
  idx2 <- lav_mat_vech_row_idx(p)
  z <- m_yc[, idx1, drop = FALSE] * m_yc[, idx2, drop = FALSE]
  if (meanstructure) {
    z <- cbind(m_y, z)
  }

  # sanity: U must match the moment vector
  if (ncol(z) != NROW(u)) {
    lav_msg_stop(gettextf(
      "dimension of U (%1$d) does not match the casewise moment vector
       (%2$d).", NROW(u), ncol(z)))
  }

  # center (rows become (s_i - sbar)')
  zc <- t(t(z) - colMeans(z))

  # N x N Gram matrix M = Zc U Zc' (= Y'Y with y_i = U^{1/2}(s_i - sbar);
  # no matrix square root needed)
  m_gram <- tcrossprod(zc %*% u, zc)

  tr_h <- sum(diag(m_gram))
  tr_h2 <- sum(m_gram * m_gram) # symmetric: tr(M^2)
  tr_d2 <- sum(diag(m_gram)^2)

  trace_ugamma2_c <-
    ((n - 2) * (n - 1) * tr_h2 - n * (n - 1) * tr_d2 + tr_h * tr_h) /
    (n * (n - 1) * (n - 2) * (n - 3))

  list(
    trace.UGamma = tr_h / (n - 1),
    trace.UGamma2 = trace_ugamma2_c
  )
}
