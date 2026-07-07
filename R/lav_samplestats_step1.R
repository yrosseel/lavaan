# STEP 1 of the Muthen (1984) three-stage estimator: univariate models
#
# for each (endogenous) variable: means/thresholds/intercepts, slopes
# (if conditional.x), variances (numeric variables only), and the
# corresponding casewise scores
#
# YR/LDW 2026: refactored; empty-category compensation moved to
#              lav_samp_step1_empty_cells(); unobserved *middle* categories
#              now give a clear error message (they were never supported:
#              the old code stopped with an obscure R error)

lav_samp_step1 <- function(y,
                                  wt = NULL, # new in 0.6-6
                                  ov_names = NULL,
                                  ov_types = NULL,
                                  ov_levels = NULL,
                                  ov_names_x = character(0L),
                                  exo = NULL,
                                  scores_flag = TRUE, # scores?
                                  allow_empty_cell = TRUE,
                                           # allow empty categories?
                                  group = 1L) { # for error message


  # just in case Y is a vector
  y <- as.matrix(y)

  nvar <- NCOL(y)
  n <- NROW(y)
  n_th <- ov_levels - 1L
  n_th[n_th == -1L] <- 1L # numeric variables (nlev = 0): one 'mean' slot
  nth <- sum(n_th)
  th_end_idx <- cumsum(n_th)
  th_start_idx <- th_end_idx - (n_th - 1L)

  # variable types; default = numeric
  nexo <- length(ov_names_x)
  if (nexo > 0L) stopifnot(NCOL(exo) == nexo)

  # means/thresholds/intercepts, slopes, variances
  th <- vector("list", length = nvar)
  th_nox <- vector("list", length = nvar)
  th_names <- vector("list", length = nvar)
  th_idx_1 <- vector("list", length = nvar)
  # which declared threshold slots come from the fit (FALSE = pseudo value
  # for an empty category; see lav_samp_step1_empty_cells)
  th_keep <- vector("list", length = nvar)
  slopes <- matrix(as.numeric(NA), nrow = nvar, ncol = nexo) # if conditional.x
  var_1 <- numeric(length = nvar) # continuous variables only

  # SCORES
  sc_var <- matrix(0, n, nvar)
  sc_sl <- matrix(0, n, nvar * nexo)
  sc_th <- matrix(0, n, nth)

  # fitted objects
  fit_1 <- vector("list", length = nvar)

  # stage one - TH/SLOPES/VAR only
  for (i in 1:nvar) {
    th_idx <- th_start_idx[i]:th_end_idx[i]
    sl_idx <- seq(i, by = nvar, length.out = nexo)
    if (ov_types[i] == "numeric") {
      fit <- lav_uvreg_fit(y = y[, i], x = exo, wt = wt)
      if (any(is.na(fit$theta))) {
        lav_msg_stop(gettextf(
          "linear regression failed for %1$s;
          X may not be of full rank in group %2$s", ov_names[i], group))
      }
      fit_1[[i]] <- fit
      # compute mean and variance
      th[[i]] <- th_nox[[i]] <- fit$theta[1L]
      var_1[i] <- fit$theta[fit$var_idx]
      th_names[[i]] <- ov_names[i]
      th_idx_1[[i]] <- 0L
      th_keep[[i]] <- TRUE
      if (scores_flag) {
        scores <- lav_uvreg_sc(y = y[, i], x = exo, wt = wt)
        sc_th[, th_idx] <- scores[, 1L]
        sc_var[, i] <- scores[, fit$var_idx]
      }
      if (nexo > 0L) {
        slopes[i, ] <- fit$theta[-c(1L, fit$var_idx)]
        if (scores_flag) {
          sc_sl[, sl_idx] <- scores[, -c(1L, fit$var_idx), drop = FALSE]
        }
        th_nox[[i]] <- mean(y[, i], na.rm = TRUE)
      }
    } else if (ov_types[i] == "ordered") {
      # check if we have enough categories in this group
      # FIXME: should we be more tolerant here???
      y_freq <- tabulate(y[, i], nbins = ov_levels[i])
      if (any(y_freq == 0L) && !allow_empty_cell) {
        lav_msg_stop(gettextf(
          "some categories of variable `%1$s' are empty in group %2$s;
          frequencies are [%3$s]", ov_names[i], group,
          lav_msg_view(y_freq, "none")))
      }
      fit <- lav_uvord_fit(y = y[, i], x = exo, wt = wt)
      if (any(is.na(fit$theta))) {
        lav_msg_stop(gettextf(
          "probit regression failed for %1$s; X may not be of full rank
          in group %2$s", ov_names[i], group))
      }
      fit_1[[i]] <- fit
      th[[i]] <- fit$theta[fit$th_idx]
      fit_nox <- lav_uvord_th(y = y[, i], wt = wt)
      th_nox[[i]] <- fit_nox
      if (scores_flag) {
        scores <- lav_uvord_sc(y = y[, i], x = exo, wt = wt)
      }

      th_keep[[i]] <- rep(TRUE, n_th[i])
      if (allow_empty_cell) {
        if (any(y_freq == 0L)) {
          empty <- lav_samp_step1_empty_cells(
            y_freq = y_freq, nlev = ov_levels[i], fit = fit,
            fit_nox = fit_nox, ov_name = ov_names[i], group = group
          )
          th[[i]] <- empty$th
          th_nox[[i]] <- empty$th_nox
          th_keep[[i]] <- empty$keep_idx
          if (scores_flag) sc_th[, th_idx[empty$keep_idx]] <-
                                                scores[, fit$th_idx, drop = FALSE]
        } else {
          # no empty cells for this variable: fill the scores as usual
          # (up to 0.7-1, they were left at zero, giving a rank-deficient
          #  A11 and a broken WLS.W whenever allow.empty.cell = TRUE)
          if (scores_flag) sc_th[, th_idx] <- scores[, fit$th_idx,
                                                     drop = FALSE]
        }
        # note: local 'fit' only -- fit_1[[i]] (assigned above) keeps the
        # fitted (possibly shorter) th_idx, which step 2 relies on
        fit$th_idx <- 1:n_th[i]
      } else {
        if (scores_flag) sc_th[, th_idx] <- scores[, fit$th_idx, drop = FALSE]
      }
      slopes[i, ] <- fit$theta[fit$slope_idx]
      if (scores_flag) {
        sc_sl[, sl_idx] <- scores[, fit$slope_idx, drop = FALSE]
      }
      var_1[i] <- 1.0
      th_names[[i]] <- paste(ov_names[i], "|t", seq_along(th[[i]]),
        sep = ""
      )
      th_idx_1[[i]] <- rep(i, length(th[[i]]))
    } else {
      lav_msg_stop(gettext("unknown ov.types:"), ov_types[i])
    }
  }

  list(
    FIT = fit_1, VAR = var_1, SLOPES = slopes,
    TH = th, TH.NOX = th_nox, TH.IDX = th_idx_1, TH.NAMES = th_names,
    TH.KEEP = th_keep,
    SC.TH = sc_th, SC.VAR = sc_var, SC.SL = sc_sl,
    th.start.idx = th_start_idx, th.end.idx = th_end_idx
  )
}

# compensate for empty categories of an ordinal variable (allow_empty_cell)
#
# lav_uvord_fit() (implicitly) drops thresholds when *extreme* categories
# are unobserved: leading empty categories disappear in the 1..ncat
# recoding, and trailing empty categories shorten tabulate(y). The declared
# number of thresholds is nlev - 1, so we place the fitted thresholds in
# their declared slots (keep_idx) and fill the unidentified slots with
# pseudo-values: -4 (first), +4 (last), previous + .01 (in between).
#
# Unobserved *middle* categories are NOT supported: the fit then keeps a
# (degenerate) threshold for the empty category, and no consistent mapping
# to the declared slots exists.
lav_samp_step1_empty_cells <- function(y_freq = NULL, nlev = NULL,
                                       fit = NULL, fit_nox = NULL,
                                       ov_name = NULL, group = 1L) {
  keep_idx <- rep(TRUE, (nlev - 1))
  mis_idx <- !keep_idx
  zidx <- y_freq == 0L
  dz <- diff(zidx) == 1L
  if (y_freq[nlev] == 0L) {
    # trailing empty categories
    wdz <- which(dz)
    nhi <- nlev - wdz[length(wdz)]
    keep_idx[(nlev - nhi) : (nlev - 1)] <- FALSE
    mis_idx[(nlev - nhi) : (nlev - 1)] <- TRUE
  }
  if (any(dz[-length(dz)])) {
    # middle empty categories
    wdz <- which(dz[-length(dz)])
    keep_idx[wdz + 1] <- FALSE
    mis_idx[wdz + 1] <- TRUE
  }
  if (y_freq[1] == 0L) {
    # leading empty categories
    nlow <- which(diff(zidx) == -1)[1]
    keep_idx[1:nlow] <- FALSE
    mis_idx[1:nlow] <- TRUE
  }

  if (sum(keep_idx) != length(fit$th_idx)) {
    lav_msg_stop(gettextf(
      "cannot handle the empty categories of variable %1$s in group %2$s
      (frequencies are [%3$s]): only unobserved categories at the extremes
      are supported.", ov_name, group, lav_msg_view(y_freq, "none")))
  }

  th <- th_nox <- rep(0, nlev - 1)
  th[keep_idx] <- fit$theta[fit$th_idx]
  th_nox[keep_idx] <- fit_nox[keep_idx]
  for (k in which(mis_idx)) {
    if (k == 1) {
      th[k] <- -4
      th_nox[k] <- -4
    } else if (k == (nlev - 1)) {
      th[k] <- 4
      th_nox[k] <- 4
    } else {
      th[k] <- th[(k - 1)] + .01
      th_nox[k] <- th_nox[(k - 1)] + .01
    }
  }

  list(th = th, th_nox = th_nox, keep_idx = keep_idx)
}
