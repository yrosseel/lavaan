lav_samplestats_step1 <- function(y,
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
  n_th[n_th == -1L] <- 1L
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
      if (scores_flag) {
        scores <- lav_uvreg_scores(y = y[, i], x = exo, wt = wt)
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
      if (length(y_freq) != ov_levels[i] && !allow_empty_cell) {
        lav_msg_stop(gettextf(
          "variable %1$s has fewer categories (%2$s) than
          expected (%3$s) in group %4$s", ov_names[i],
          length(y_freq), ov_levels[i], group))
      }
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
        scores <- lav_uvord_scores(y = y[, i], x = exo, wt = wt)
      }

      if (allow_empty_cell) {
        if (any(y_freq == 0L)) {
          ## lav_uvord_fit drops thresholds if extreme categories
          #                           are missing, but not otherwise
          exidx <- rep(TRUE, (ov_levels[i] - 1))
          misidx <- !exidx
          zidx <- y_freq == 0L
          dz <- diff(zidx) == 1L
          if (y_freq[ov_levels[i]] == 0L) {
            wdz <- which(dz)
            nhi <- ov_levels[i] - wdz[length(wdz)]
            exidx[(ov_levels[i] - nhi) : (ov_levels[i] - 1)] <- FALSE
            misidx[(ov_levels[i] - nhi) : (ov_levels[i] - 1)] <- TRUE
          }
          if (any(dz[-length(dz)])) {
            wdz <- which(dz[-length(dz)])
            exidx[wdz + 1] <- FALSE
            misidx[wdz + 1] <- TRUE
          }
          if (y_freq[1] == 0L) {
            nlow <- which(diff(zidx) == -1)[1]
            exidx[1:nlow] <- FALSE
            misidx[1:nlow] <- TRUE
          }
          th[[i]] <- th_nox[[i]] <- rep(0, ov_levels[i] - 1)
          th[[i]][exidx] <- fit$theta[fit$th_idx]
          th_nox[[i]][exidx] <- fit_nox[exidx]
          for (k in which(misidx)) {
            if (k == 1) {
              th[[i]][k] <- -4
              th_nox[[i]][k] <- -4
            } else if (k == (ov_levels[i] - 1)) {
              th[[i]][k] <- 4
              th_nox[[i]][k] <- 4
            } else {
              th[[i]][k] <- th[[i]][(k - 1)] + .01
              th_nox[[i]][k] <- th_nox[[i]][(k - 1)] + .01
            }
          }
          if (scores_flag) sc_th[, th_idx[!misidx]] <-
                                                scores[, fit$th_idx, drop = FALSE]
        } else if (length(y_freq) != ov_levels[i]) {
          nz <- ov_levels[i] - length(y_freq)
          th[[i]] <- c(th[[i]], th[[i]][length(y_freq)] + (1:nz) * .01)
          if (scores_flag) sc_th[, th_idx[seq_along(y_freq)]] <-
                                                scores[, fit$th_idx, drop = FALSE]
        }
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
    SC.TH = sc_th, SC.VAR = sc_var, SC.SL = sc_sl,
    th.start.idx = th_start_idx, th.end.idx = th_end_idx
  )
}
