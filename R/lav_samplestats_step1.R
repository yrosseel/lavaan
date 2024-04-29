lav_samplestats_step1 <- function(Y,
                                  wt = NULL, # new in 0.6-6
                                  ov.names = NULL,
                                  ov.types = NULL,
                                  ov.levels = NULL,
                                  ov.names.x = character(0L),
                                  eXo = NULL,
                                  scores.flag = TRUE, # scores?
                                  group = 1L) { # for error message


  # just in case Y is a vector
  Y <- as.matrix(Y)

  nvar <- NCOL(Y)
  N <- NROW(Y)
  nTH <- ov.levels - 1L
  nTH[nTH == -1L] <- 1L
  nth <- sum(nTH)
  th.end.idx <- cumsum(nTH)
  th.start.idx <- th.end.idx - (nTH - 1L)

  # variable types; default = numeric
  nexo <- length(ov.names.x)
  if (nexo > 0L) stopifnot(NCOL(eXo) == nexo)

  # means/thresholds/intercepts, slopes, variances
  TH <- vector("list", length = nvar)
  TH.NOX <- vector("list", length = nvar)
  TH.NAMES <- vector("list", length = nvar)
  TH.IDX <- vector("list", length = nvar)
  SLOPES <- matrix(as.numeric(NA), nrow = nvar, ncol = nexo) # if conditional.x
  VAR <- numeric(length = nvar) # continuous variables only

  # SCORES
  SC.VAR <- matrix(0, N, nvar)
  SC.SL <- matrix(0, N, nvar * nexo)
  SC.TH <- matrix(0, N, nth)

  # fitted objects
  FIT <- vector("list", length = nvar)

  # stage one - TH/SLOPES/VAR only
  for (i in 1:nvar) {
    th.idx <- th.start.idx[i]:th.end.idx[i]
    sl.idx <- seq(i, by = nvar, length.out = nexo)
    if (ov.types[i] == "numeric") {
      fit <- lav_uvreg_fit(y = Y[, i], X = eXo, wt = wt)
      if (any(is.na(fit$theta))) {
        lav_msg_stop(gettextf(
          "linear regression failed for %1$s;
          X may not be of full rank in group %2$s", ov.names[i], group))
      }
      FIT[[i]] <- fit
      # compute mean and variance
      TH[[i]] <- TH.NOX[[i]] <- fit$theta[1L]
      VAR[i] <- fit$theta[fit$var.idx]
      TH.NAMES[[i]] <- ov.names[i]
      TH.IDX[[i]] <- 0L
      if (scores.flag) {
        scores <- lav_uvreg_scores(y = Y[, i], X = eXo, wt = wt)
        SC.TH[, th.idx] <- scores[, 1L]
        SC.VAR[, i] <- scores[, fit$var.idx]
      }
      if (nexo > 0L) {
        SLOPES[i, ] <- fit$theta[-c(1L, fit$var.idx)]
        if (scores.flag) {
          SC.SL[, sl.idx] <- scores[, -c(1L, fit$var.idx), drop = FALSE]
        }
        TH.NOX[[i]] <- mean(Y[, i], na.rm = TRUE)
      }
    } else if (ov.types[i] == "ordered") {
      # check if we have enough categories in this group
      # FIXME: should we more tolerant here???
      y.freq <- tabulate(Y[, i], nbins = ov.levels[i])
      if (length(y.freq) != ov.levels[i]) {
        lav_msg_stop(gettextf(
          "variable %1$s has fewer categories (%2$s) than
          expected (%3$s) in group %4$s", ov.names[i],
          length(y.freq), ov.levels[i], group))
      }
      if (any(y.freq == 0L)) {
        lav_msg_stop(gettextf(
          "some categories of variable `%1$s' are empty in group %2$s;
          frequencies are [%3$s]", ov.names[i], group,
          lav_msg_view(y.freq, "none")))
      }
      fit <- lav_uvord_fit(y = Y[, i], X = eXo, wt = wt)
      if (any(is.na(fit$theta))) {
        lav_msg_stop(gettextf(
          "probit regression failed for %1$s; X may not be of full rank
          in group %2$s", ov.names[i], group))
      }
      FIT[[i]] <- fit
      TH[[i]] <- fit$theta[fit$th.idx]
      TH.NOX[[i]] <- lav_uvord_th(y = Y[, i], wt = wt)
      if (scores.flag) {
        scores <- lav_uvord_scores(y = Y[, i], X = eXo, wt = wt)
        SC.TH[, th.idx] <- scores[, fit$th.idx, drop = FALSE]
      }
      SLOPES[i, ] <- fit$theta[fit$slope.idx]
      if (scores.flag) {
        SC.SL[, sl.idx] <- scores[, fit$slope.idx, drop = FALSE]
      }
      VAR[i] <- 1.0
      TH.NAMES[[i]] <- paste(ov.names[i], "|t", 1:length(TH[[i]]),
        sep = ""
      )
      TH.IDX[[i]] <- rep(i, length(TH[[i]]))
    } else {
      lav_msg_stop(gettext("unknown ov.types:"), ov.types[i])
    }
  }

  list(
    FIT = FIT, VAR = VAR, SLOPES = SLOPES,
    TH = TH, TH.NOX = TH.NOX, TH.IDX = TH.IDX, TH.NAMES = TH.NAMES,
    SC.TH = SC.TH, SC.VAR = SC.VAR, SC.SL = SC.SL,
    th.start.idx = th.start.idx, th.end.idx = th.end.idx
  )
}
