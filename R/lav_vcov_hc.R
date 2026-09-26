# information.meat.hc: small-sample corrections of the first-order
# (casewise) sandwich standard errors
#
# YR 26 Sept 2026
#
# Every first-order sandwich in lavaan has the form V = B M B, with bread B
# and meat M = (1/n) sum_i a_i a_i', where a_i is the contribution of case
# i to the estimating function: the casewise score (se =
# "robust.huber.white"), the casewise moment contribution mapped through
# W Delta (se = "robust.sem"), or the casewise influence contribution of
# the latent moments (the local sam() SEs). Writing IF_i = B a_i for the
# casewise influence rows, V = (1/n) sum_i IF_i IF_i'. This is the analogue
# of the HC0 estimator of a regression, and it inherits HC0's small-sample
# downward bias. The corrections mirror the regression ones:
#
#  - "HC1": V * n/(n - p), with p the number of free parameters (a global
#           degrees-of-freedom factor; the only correction that applies to
#           parameters without a regression design)
#  - "HC2": the elements of IF_i that belong to the regression coefficients
#           (including the intercept) of equation j are multiplied by
#           (1 - h_ij)^(-1/2), with h_ij the leverage of case i in the
#           design Z_j = [1, casewise values of the predictors of equation j]
#           (observed predictors: the data; latent predictors: Bartlett
#           factor scores; in local sam(): the step-1 factor scores, with
#           products for latent interaction terms)
#  - "HC3": idem with (1 - h_ij)^(-1)
#
# For an observed-variable regression (fixed.x = TRUE) HC2/HC3 reproduce
# sandwich::vcovHC(lm(...), type = "HC2"/"HC3"). Parameters that are not
# regression coefficients keep their HC0 influence rows under HC2/HC3.
#
# Implementation: the correction only touches the casewise part of the
# meat, so every consumer computes
#     V_hc = V + (1/n) sum_i [ (w_i * IF_i)(w_i * IF_i)' - IF_i IF_i' ]
# where w_i is the vector of per-parameter weights of case i (1 for the
# parameters that are not adjusted). Any other, non-casewise ingredient of
# the meat (eg the finite-sample terms of the unbiased Gamma, or the
# additive part of Gamma.eta in sam()) is retained as is, and "HC0"
# leaves everything untouched.

# the leverage exponent: (1 - h)^(-delta/2)
lav_hc_delta <- function(hc = "HC0") {
  switch(toupper(hc), HC2 = 1, HC3 = 2, 0)
}

# does this hc value ask for a leverage adjustment?
lav_hc_leverage_flag <- function(hc = "HC0") {
  lav_hc_delta(hc) > 0
}

# number of free parameters for the HC1 factor: the free parameters of the
# model, minus the number of (explicit) equality constraints; for
# ceq.simple.only models nx.free is already the collapsed count
lav_hc_npar <- function(lavmodel = NULL) {
  npar <- lavmodel@nx.free
  if (!lavmodel@ceq.simple.only) {
    npar <- npar - length(lavmodel@ceq.linear.idx) -
      length(lavmodel@ceq.nonlinear.idx)
  }
  npar
}

# sam(): the number of free structural parameters for the HC1 factor --
# the regression coefficients, intercepts and residual (co)variances of
# the endogenous variables; the (saturated) moments of the exogenous
# variables are not counted (the report's N/(N - 15) for 14 paths + 1
# residual variance)
lav_hc_npar_struc <- function(lavpartable = NULL) {
  pt <- lavpartable
  grp <- pt$group
  if (is.null(grp)) {
    grp <- rep.int(1L, length(pt$lhs))
  }
  keep <- logical(length(pt$lhs))
  for (g in unique(grp)) {
    endo <- unique(pt$lhs[pt$op == "~" & grp == g])
    keep[grp == g] <- pt$op[grp == g] == "~" |
      (pt$op[grp == g] %in% c("~~", "~1") &
       (pt$lhs[grp == g] %in% endo | pt$rhs[grp == g] %in% endo))
  }
  length(unique(pt$free[keep & pt$free > 0L]))
}

# the HC1 factor n/(n - p)
lav_hc1_factor <- function(n = 0L, npar = 0L) {
  if (npar >= n) {
    lav_msg_warn(gettextf(
      "information_meat_hc = \"HC1\": the number of free parameters (%1$s)
       is not smaller than the number of observations (%2$s); no
       correction is applied.", npar, n))
    return(1)
  }
  n / (n - npar)
}

# the regression equations of a parameter table, per group: for every
# endogenous variable with (at least one) '~' row, its predictors and the
# free-parameter indices (the 'free' column) of its regression
# coefficients (the '~' rows plus the '~1' intercept row of the same
# variable)
lav_hc_equations <- function(lavpartable = NULL, ngroups = 1L) {
  pt <- lavpartable
  free_idx <- pt$free
  grp <- pt$group
  if (is.null(grp) || ngroups == 1L) {
    grp <- rep.int(1L, length(pt$lhs))
  }
  out <- vector("list", ngroups)
  for (g in seq_len(ngroups)) {
    reg_idx <- which(pt$op == "~" & grp == g)
    eqs <- list()
    for (y in unique(pt$lhs[reg_idx])) {
      idx <- reg_idx[pt$lhs[reg_idx] == y]
      int_idx <- which(pt$op == "~1" & pt$lhs == y & grp == g)
      fidx <- c(free_idx[idx], free_idx[int_idx])
      fidx <- unique(fidx[fidx > 0L])
      if (length(fidx) == 0L) {
        next
      }
      eqs[[y]] <- list(lhs = y, rhs = unique(pt$rhs[idx]), free = fidx)
    }
    out[[g]] <- eqs
  }
  out
}

# leverage of the rows of the design [1, z]; missing values in z are
# replaced by the column means (a missing predictor value then contributes
# nothing to the leverage of that case, beyond the intercept); leverages
# (numerically) equal to one would give an infinite weight: they are
# capped (with a warning)
lav_hc_leverage <- function(z = NULL, eq_name = "") {
  z <- as.matrix(z)
  if (anyNA(z)) {
    for (j in seq_len(ncol(z))) {
      na_idx <- which(is.na(z[, j]))
      if (length(na_idx) > 0L) {
        z[na_idx, j] <- mean(z[, j], na.rm = TRUE)
      }
    }
    z[is.na(z)] <- 0 # a column with no observed values at all
  }
  z <- cbind(1, z)
  qz <- qr(z)
  q <- qr.Q(qz)[, seq_len(qz$rank), drop = FALSE]
  h <- rowSums(q * q)
  cap <- 1 - 1e-8
  if (any(h > cap)) {
    lav_msg_warn(gettextf(
      "information_meat_hc: %1$s case(s) have a leverage of one in the
       regression equation of %2$s; their leverage is capped.",
      sum(h > cap), eq_name))
    h[h > cap] <- cap
  }
  h
}

# the casewise values of the predictors of an equation: observed
# variables from the data, latent variables from the factor scores, and
# products ("a:b") of the latter for (latent) interaction terms
lav_hc_design <- function(pred_names = character(0L), ov_data = NULL,
                          lv_scores = NULL, n = 0L) {
  get_col <- function(name) {
    if (!is.null(ov_data) && name %in% colnames(ov_data)) {
      return(ov_data[, name])
    }
    if (!is.null(lv_scores) && name %in% colnames(lv_scores)) {
      return(lv_scores[, name])
    }
    if (grepl(":", name, fixed = TRUE)) {
      parts <- strsplit(name, ":", fixed = TRUE)[[1]]
      return(Reduce(`*`, lapply(parts, get_col)))
    }
    lav_msg_stop(gettextf(
      "information_meat_hc: no casewise values are available for the
       predictor %s.", dQuote(name, q = FALSE)))
  }
  z <- matrix(0, nrow = n, ncol = length(pred_names))
  for (j in seq_along(pred_names)) {
    z[, j] <- get_col(pred_names[j])
  }
  z
}

# the n x npar matrix of per-case, per-parameter weights (1 = no
# adjustment) for one group; free_map (optional) maps the free-parameter
# indices of the equations to the columns of the influence rows (eg the
# step-2 parameters of a sam() model); unco_map (optional, ceq.simple.only
# models, whose vcov machinery works in the 'unco' space) gives the
# (collapsed) free index of every unco column: the weights are computed
# in the collapsed space and expanded, so that all copies of a parameter
# constrained equal across groups share the leverage of the equation the
# case contributes to. A parameter that appears in several equations of
# the same group (equality constraints across equations) takes the
# leverage of the first equation it belongs to.
lav_hc_weights <- function(hc = "HC0", eqs = list(), ov_data = NULL,
                           lv_scores = NULL, n = 0L, npar = 0L,
                           free_map = NULL, unco_map = NULL,
                           eq_classes = NULL) {
  delta <- lav_hc_delta(hc)
  if (!is.null(unco_map)) {
    w_free <- lav_hc_weights(
      hc = hc, eqs = eqs, ov_data = ov_data, lv_scores = lv_scores,
      n = n, npar = max(unco_map), free_map = free_map,
      eq_classes = eq_classes
    )
    return(w_free[, unco_map, drop = FALSE])
  }
  w <- matrix(1, nrow = n, ncol = npar)
  if (delta == 0 || length(eqs) == 0L) {
    return(w)
  }
  assigned <- logical(npar)
  for (eq in eqs) {
    cols <- eq$free
    if (!is.null(free_map)) {
      cols <- free_map[cols]
    }
    cols <- cols[!is.na(cols) & cols > 0L & cols <= npar]
    cols <- cols[!assigned[cols]]
    if (length(cols) == 0L) {
      next
    }
    z <- lav_hc_design(pred_names = eq$rhs, ov_data = ov_data,
                       lv_scores = lv_scores, n = n)
    h <- lav_hc_leverage(z, eq_name = eq$lhs)
    w[, cols] <- (1 - h)^(-delta / 2)
    assigned[cols] <- TRUE
  }
  # parameters tied by explicit equality constraints (eg the same
  # regression coefficient in several groups): every member of the class
  # represents the same parameter, so a case's weight (from the equation
  # it contributes to) applies to all of them
  for (cl in eq_classes) {
    cols <- cl
    if (!is.null(free_map)) {
      cols <- free_map[cl]
    }
    cols <- cols[!is.na(cols) & cols > 0L & cols <= npar]
    if (length(cols) < 2L) {
      next
    }
    wc <- rep(1, n)
    for (cc in cols) {
      idx <- which(wc == 1 & w[, cc] != 1)
      wc[idx] <- w[idx, cc]
    }
    w[, cols] <- wc
  }
  w
}

# free parameters tied by (explicit) equality constraints: the same
# non-empty label, or an '==' row between labels/plabels; returns a list
# of integer vectors (free indices) with at least two members. (For
# ceq.simple.only models the tied parameters already share one free
# index, and this returns an empty list.)
lav_hc_eq_classes <- function(lavpartable = NULL) {
  pt <- lavpartable
  free_rows <- which(pt$free > 0L)
  if (length(free_rows) == 0L) {
    return(list())
  }
  # union-find over the free indices
  parent <- seq_len(max(pt$free))
  find <- function(i) {
    while (parent[i] != i) {
      i <- parent[i]
    }
    i
  }
  union <- function(i, j) {
    ri <- find(i)
    rj <- find(j)
    if (ri != rj) {
      parent[max(ri, rj)] <<- min(ri, rj)
    }
  }
  # 1. shared labels
  if (!is.null(pt$label)) {
    lab <- pt$label[free_rows]
    for (l in unique(lab[nzchar(lab)])) {
      idx <- unique(pt$free[free_rows[lab == l]])
      if (length(idx) > 1L) {
        for (k in idx[-1]) union(idx[1], k)
      }
    }
  }
  # 2. explicit '==' rows (lhs/rhs are labels or plabels)
  eq_rows <- which(pt$op == "==")
  if (length(eq_rows) > 0L) {
    to_free <- function(name) {
      idx <- integer(0L)
      if (!is.null(pt$plabel)) {
        idx <- which(pt$plabel == name & pt$free > 0L)
      }
      if (length(idx) == 0L && !is.null(pt$label)) {
        idx <- which(pt$label == name & pt$free > 0L)
      }
      unique(pt$free[idx])
    }
    for (r in eq_rows) {
      f1 <- to_free(pt$lhs[r])
      f2 <- to_free(pt$rhs[r])
      if (length(f1) == 1L && length(f2) == 1L) {
        union(f1, f2)
      }
    }
  }
  roots <- vapply(seq_along(parent), find, integer(1L))
  classes <- split(seq_along(parent), roots)
  classes[vapply(classes, length, integer(1L)) > 1L]
}

# ceq.simple.only models: the (collapsed) free index of every 'unco'
# column (NULL otherwise)
lav_hc_unco_map <- function(lavmodel = NULL, lavpartable = NULL) {
  if (!lavmodel@ceq.simple.only) {
    return(NULL)
  }
  lavpartable$free[lavpartable$free > 0L]
}

# the leverage adjustment of a crossproduct: (w * if_rows)'(w * if_rows) -
# if_rows' if_rows
lav_hc_crossprod_diff <- function(if_rows = NULL, w = NULL) {
  if (is.null(w)) {
    return(matrix(0, ncol(if_rows), ncol(if_rows)))
  }
  crossprod(w * if_rows) - crossprod(if_rows)
}

# Bartlett factor scores (per group, with the latent variable names as
# column names) of the latent variables that appear as predictors in the
# regression equations; NULL if no latent predictor is needed
lav_hc_lv_scores <- function(lavmodel = NULL, lavdata = NULL,
                             lavsamplestats = NULL, lavimplied = NULL,
                             lavpartable = NULL, eqs = NULL) {
  ngroups <- lavdata@ngroups
  lv_names_all <- character(0L)
  lambda_idx <- which(names(lavmodel@GLIST) == "lambda")
  for (g in seq_len(ngroups)) {
    lv_names_all <- c(lv_names_all,
                      lavmodel@dimNames[[lambda_idx[g]]][[2L]])
  }
  needed <- FALSE
  for (g in seq_len(ngroups)) {
    for (eq in eqs[[g]]) {
      preds <- unlist(strsplit(eq$rhs, ":", fixed = TRUE))
      if (any(!preds %in% lavdata@ov.names[[g]] &
              preds %in% lv_names_all)) {
        needed <- TRUE
      }
    }
  }
  if (!needed) {
    return(NULL)
  }
  if (any(lavdata@ov$type != "numeric")) {
    lav_msg_stop(gettextf(
      "information_meat_hc = %s (factor-score leverage of latent
       predictors) is only available for continuous data.",
      dQuote("HC2/HC3", q = FALSE)))
  }
  lavpta <- lav_pt_attributes(lavpartable)
  fs <- lav_predict_eta_normal(
    lavmodel = lav_model_delta_absorb(lavmodel), lavdata = lavdata,
    lavsamplestats = lavsamplestats, lavimplied = lavimplied,
    lavpta = lavpta, method = "ml"
  )
  for (g in seq_len(ngroups)) {
    colnames(fs[[g]]) <- lavmodel@dimNames[[lambda_idx[g]]][[2L]]
  }
  fs
}

# the casewise (model-parameter) scores of a (single-level, continuous or
# incomplete) ML fit, for the se = "robust.huber.white" sandwich: one row
# per (original) case, in the ordering of lavdata@case.idx; ceq.simple.only
# models return the 'unco' space scores (matching the vcov machinery)
lav_hc_model_scores <- function(lavmodel = NULL, lavdata = NULL,
                                lavsamplestats = NULL, lavimplied = NULL,
                                lavoptions = NULL) {
  if (lavmodel@estimator != "ML") {
    lav_msg_stop(gettextf(
      "information_meat_hc = %1$s is only available for estimator = ML
       with se = \"robust.huber.white\" (not for estimator = %2$s).",
      dQuote(lavoptions$information.meat.hc, q = FALSE),
      dQuote(lavmodel@estimator, q = FALSE)))
  }
  if (length(lavdata@sampling.weights) > 0L) {
    lav_msg_stop(gettextf(
      "information_meat_hc = %s is not available with sampling weights.",
      dQuote(lavoptions$information.meat.hc, q = FALSE)))
  }
  ntab <- unlist(lavdata@norig)
  ntot <- sum(ntab)
  ngroups <- lavsamplestats@ngroups
  npar <- lavmodel@nx.free
  if (lavmodel@ceq.simple.only) {
    npar <- lavmodel@nx.unco
  }
  # the per-group implied moments, in the layout of fitted()
  moments <- vector("list", ngroups)
  for (g in seq_len(ngroups)) {
    moments[[g]] <- list(
      cov = lavimplied$cov[[g]], mean = lavimplied$mean[[g]],
      res.cov = lavimplied$res.cov[[g]], res.int = lavimplied$res.int[[g]],
      res.slopes = lavimplied$res.slopes[[g]]
    )
  }
  if (ngroups == 1L) {
    moments <- moments[[1]]
  }
  sc <- lav_sc_ml(
    ntab = ntab, ntot = ntot, npar = npar, moments = moments,
    lavdata = lavdata, lavsamplestats = lavsamplestats,
    lavmodel = lavmodel, lavoptions = lavoptions, scaling = FALSE
  )
  sc[is.na(sc)] <- 0 # empty cases
  sc
}

# the casewise (centered) moment contributions of group g -- the rows
# whose crossproduct is the (ADF) Gamma that the se = "robust.sem"
# sandwich consumes -- recomputed with the fit-time (NACOV) recipe; NULL
# when the setting has no such rows (categorical, incomplete, weighted,
# clustered or two-level data, correlation structures, no raw data)
lav_hc_gamma_rows <- function(lavdata = NULL, lavoptions = NULL,
                              lavsamplestats = NULL, g = 1L) {
  recipe <- lav_gamma_recipe(lavoptions = lavoptions, lavdata = lavdata,
                             nacov_compute = TRUE)
  if (recipe$categorical || recipe$multilevel || recipe$correlation ||
      length(lavdata@cluster) > 0L || recipe$wt ||
      !recipe$missing %in% c("listwise") ||
      lavdata@data.type != "full" || lavsamplestats@NACOV.user) {
    return(NULL)
  }
  if (recipe$conditional.x) {
    y <- cbind(lavdata@X[[g]], lavdata@eXo[[g]])
  } else {
    y <- lavdata@X[[g]]
  }
  if (!is.matrix(y) || anyNA(y)) {
    return(NULL)
  }
  n <- nrow(y)
  n_xi <- if (isTRUE(recipe$n.minus.one)) n - 1L else n
  zc <- lav_samp_gamma_zc(
    m_y = y, x_idx = lavsamplestats@x.idx[[g]],
    fixed_x = recipe$fixed.x, conditional_x = recipe$conditional.x,
    meanstructure = recipe$meanstructure,
    slopestructure = recipe$conditional.x, n_xi = n_xi
  )
  if (isTRUE(recipe$group.w.free)) {
    zc <- cbind(0, zc) # the group-weight statistic has no casewise row
  }
  zc
}

# the (nvcov-scale) leverage adjustment of the robust.sem-type sandwich
#   nvcov = E.inv [ sum_g fg^2/fg1_g WD_g' Gamma_g WD_g ] E.inv,
#   Gamma_g = crossprod(rows_g) / den_g (+ non-casewise terms)
# with WD_g = W_g Delta_g; rows and w are per-group lists (a NULL w_g
# leaves group g unchanged)
lav_hc_sem_adjust <- function(e_inv = NULL, delta = NULL, wls_v = NULL,
                              rows = NULL, w = NULL, nobs = NULL,
                              ntotal = NULL, diag_wls_v = FALSE,
                              den = NULL, fg1 = NULL) {
  npar <- ncol(e_inv)
  adj <- matrix(0, npar, npar)
  nobs <- unlist(nobs)
  for (g in seq_along(rows)) {
    if (is.null(w[[g]]) || is.null(rows[[g]])) {
      next
    }
    if (diag_wls_v) {
      wd <- wls_v[[g]] * delta[[g]]
    } else {
      wd <- wls_v[[g]] %*% delta[[g]]
    }
    if_rows <- rows[[g]] %*% (wd %*% e_inv)
    fg <- nobs[g] / ntotal
    fg1_g <- if (is.null(fg1)) fg else fg1[g]
    den_g <- if (is.null(den)) nobs[g] else den[g]
    adj <- adj + (fg * fg / fg1_g) / den_g *
      lav_hc_crossprod_diff(if_rows, w[[g]])
  }
  adj
}

# se = "robust.huber.white": the leverage-adjusted nvcov
#   nvcov = E.inv B0 E.inv, B0 = crossprod(sc) / ntotal
lav_hc_nvcov_robust_sandwich <- function(nvar_cov = NULL, e_inv = NULL,
                                         lavmodel = NULL, lavdata = NULL,
                                         lavsamplestats = NULL,
                                         lavimplied = NULL,
                                         lavoptions = NULL,
                                         lavpartable = NULL) {
  hc <- lavoptions$information.meat.hc
  sc <- lav_hc_model_scores(
    lavmodel = lavmodel, lavdata = lavdata, lavsamplestats = lavsamplestats,
    lavimplied = lavimplied, lavoptions = lavoptions
  )
  if_rows <- sc %*% e_inv
  ngroups <- lavdata@ngroups
  eqs <- lav_hc_equations(lavpartable, ngroups = ngroups)
  unco_map <- lav_hc_unco_map(lavmodel, lavpartable)
  eq_classes <- lav_hc_eq_classes(lavpartable)
  lv_scores <- lav_hc_lv_scores(
    lavmodel = lavmodel, lavdata = lavdata, lavsamplestats = lavsamplestats,
    lavimplied = lavimplied, lavpartable = lavpartable, eqs = eqs
  )
  w <- matrix(1, nrow(if_rows), ncol(if_rows))
  for (g in seq_len(ngroups)) {
    ov_data <- lavdata@X[[g]]
    colnames(ov_data) <- lavdata@ov.names[[g]]
    if (lavmodel@conditional.x && !is.null(lavdata@eXo[[g]])) {
      exo <- lavdata@eXo[[g]]
      colnames(exo) <- lavdata@ov.names.x[[g]]
      ov_data <- cbind(ov_data, exo)
    }
    w_g <- lav_hc_weights(
      hc = hc, eqs = eqs[[g]], ov_data = ov_data,
      lv_scores = lv_scores[[g]], n = nrow(ov_data), npar = ncol(if_rows),
      unco_map = unco_map, eq_classes = eq_classes
    )
    w[lavdata@case.idx[[g]], ] <- w_g
  }
  nvar_cov + lav_hc_crossprod_diff(if_rows, w) / lavsamplestats@ntotal
}

# se = "robust.sem": the leverage-adjusted nvcov (continuous, complete,
# single-level data; the rows of the ADF Gamma are recomputed)
lav_hc_nvcov_robust_sem <- function(nvar_cov = NULL, e_inv = NULL,
                                    delta = NULL, wls_v = NULL,
                                    lavmodel = NULL, lavdata = NULL,
                                    lavsamplestats = NULL,
                                    lavimplied = NULL, lavoptions = NULL,
                                    lavpartable = NULL) {
  hc <- lavoptions$information.meat.hc
  ngroups <- lavdata@ngroups
  rows <- vector("list", ngroups)
  for (g in seq_len(ngroups)) {
    rows[g] <- list(lav_hc_gamma_rows(
      lavdata = lavdata, lavoptions = lavoptions,
      lavsamplestats = lavsamplestats, g = g
    ))
    if (is.null(rows[[g]])) {
      lav_msg_stop(gettextf(
        "information_meat_hc = %s with se = \"robust.sem\" is only
         available for continuous, complete, single-level (unweighted)
         raw data; use se = \"robust.huber.white\" for incomplete data.",
        dQuote(hc, q = FALSE)))
    }
    if (ncol(rows[[g]]) != nrow(delta[[g]])) {
      lav_msg_stop(gettext(
        "information_meat_hc: the casewise moment contributions do not
         match the model statistics (robust.sem)."))
    }
  }
  eqs <- lav_hc_equations(lavpartable, ngroups = ngroups)
  unco_map <- lav_hc_unco_map(lavmodel, lavpartable)
  eq_classes <- lav_hc_eq_classes(lavpartable)
  lv_scores <- lav_hc_lv_scores(
    lavmodel = lavmodel, lavdata = lavdata, lavsamplestats = lavsamplestats,
    lavimplied = lavimplied, lavpartable = lavpartable, eqs = eqs
  )
  w <- vector("list", ngroups)
  for (g in seq_len(ngroups)) {
    ov_data <- lavdata@X[[g]]
    colnames(ov_data) <- lavdata@ov.names[[g]]
    if (lavmodel@conditional.x && !is.null(lavdata@eXo[[g]])) {
      exo <- lavdata@eXo[[g]]
      colnames(exo) <- lavdata@ov.names.x[[g]]
      ov_data <- cbind(ov_data, exo)
    }
    w[[g]] <- lav_hc_weights(
      hc = hc, eqs = eqs[[g]], ov_data = ov_data,
      lv_scores = lv_scores[[g]], n = nrow(ov_data), npar = ncol(e_inv),
      unco_map = unco_map, eq_classes = eq_classes
    )
  }
  recipe <- lav_gamma_recipe(lavoptions = lavoptions, lavdata = lavdata,
                             nacov_compute = TRUE)
  nobs <- unlist(lavsamplestats@nobs)
  den <- if (isTRUE(recipe$n.minus.one)) nobs - 1L else nobs
  fg1 <- NULL
  if (isTRUE(lavoptions$gamma.vcov.mplus)) {
    fg1 <- (nobs - 1) / lavsamplestats@ntotal
  }
  adj <- lav_hc_sem_adjust(
    e_inv = e_inv, delta = delta, wls_v = wls_v, rows = rows, w = w,
    nobs = nobs, ntotal = lavsamplestats@ntotal,
    diag_wls_v = lavmodel@estimator %in% c("DWLS", "ULS"),
    den = den, fg1 = fg1
  )
  nvar_cov + adj
}

# the ntot x npar weight matrix of a single-level (joint) model, with the
# rows in the ordering of the casewise scores (lavdata@case.idx); free_map
# (optional) maps the free-parameter indices of the parameter table to
# the columns of the influence rows
lav_hc_weights_all <- function(hc = "HC0", lavmodel = NULL, lavdata = NULL,
                               lavsamplestats = NULL, lavimplied = NULL,
                               lavpartable = NULL, npar = 0L,
                               free_map = NULL, ntot = 0L) {
  ngroups <- lavdata@ngroups
  eqs <- lav_hc_equations(lavpartable, ngroups = ngroups)
  eq_classes <- lav_hc_eq_classes(lavpartable)
  lv_scores <- lav_hc_lv_scores(
    lavmodel = lavmodel, lavdata = lavdata, lavsamplestats = lavsamplestats,
    lavimplied = lavimplied, lavpartable = lavpartable, eqs = eqs
  )
  w <- matrix(1, ntot, npar)
  for (g in seq_len(ngroups)) {
    ov_data <- lavdata@X[[g]]
    colnames(ov_data) <- lavdata@ov.names[[g]]
    if (lavmodel@conditional.x && !is.null(lavdata@eXo[[g]])) {
      exo <- lavdata@eXo[[g]]
      colnames(exo) <- lavdata@ov.names.x[[g]]
      ov_data <- cbind(ov_data, exo)
    }
    w_g <- lav_hc_weights(
      hc = hc, eqs = eqs[[g]], ov_data = ov_data,
      lv_scores = lv_scores[[g]], n = nrow(ov_data), npar = npar,
      free_map = free_map, eq_classes = eq_classes
    )
    w[lavdata@case.idx[[g]], ] <- w_g
  }
  w
}

# local sam(): the step-1 (Bartlett-type) factor scores of group g, with
# the latent variable names as column names (products for interaction
# terms are formed by lav_hc_design()); with missing = "ml" the scores are
# computed per missing pattern, and the unscoreable cases get NA rows
# (they receive the mean-imputed leverage of lav_hc_leverage())
lav_sam_hc_fs <- function(step1 = NULL, fit = NULL, g = 1L) {
  y <- fit@Data@X[[g]]
  if (!is.matrix(y)) {
    lav_msg_stop(gettext(
      "information_meat_hc: the leverage adjustment of the local standard
       errors needs the raw data."))
  }
  m <- step1$M[[g]]
  lv_names <- step1$LV.NAMES[[g]]
  if (is.null(m) || nrow(m) != length(lv_names) || ncol(m) != ncol(y)) {
    lav_msg_stop(gettext(
      "information_meat_hc: the step-1 mapping matrix does not conform with
       the data (local sam)."))
  }
  if (!anyNA(y)) {
    fs <- y %*% t(m)
  } else {
    # per-pattern factor scores at the step-1 estimates (as in
    # lav_sam_gamma_add())
    lavmodel <- fit@Model
    pt <- step1$PT
    lavmodel_0 <- lav_model_set_parameters(lavmodel,
      x = pt$est[pt$free > 0 & !duplicated(pt$free)])
    mm_idx <- lav_model_group_mm_indices(lavmodel@nmat)[[g]]
    glist_g <- lavmodel_0@GLIST[mm_idx]
    lambda_0 <- glist_g$lambda
    int_idx <- fit@pta$vidx$lv.interaction[[g]]
    if (length(int_idx) > 0L) {
      lambda_0 <- lambda_0[, -int_idx, drop = FALSE]
    }
    nu_0 <- glist_g$nu
    if (is.null(nu_0)) {
      nu_0 <- numeric(ncol(y))
    }
    out_mi <- lav_sam_fs_missing(
      y = y, mm_lambda = lambda_0, mm_theta = glist_g$theta,
      mm_nu = nu_0, s = step1$COV[[g]],
      method = step1$local.options$m_method
    )
    fs <- out_mi$fs
    fs[!out_mi$ok, ] <- NA
  }
  colnames(fs) <- lv_names
  fs
}

# local sam(): the leverage-adjusted (HC2/HC3) local vcov, in the
# parameter space of the structural fit (FIT.PA)
#   vcov = (1/N) E.inv [ sum_g fg WD_g' Gamma.eta_g WD_g ] E.inv
# with Gamma.eta_g = crossprod(rows_g) / n_g + non-casewise terms; the
# casewise rows are the ones stored in step1$Gamma.eta.rows
lav_sam_hc_local <- function(fit_pa = NULL, step1 = NULL, fit = NULL,
                             hc = "HC2", step2_rm_idx = integer(0L)) {
  rows <- step1$Gamma.eta.rows
  ngroups <- fit_pa@Data@ngroups
  if (is.null(rows) || length(rows) != ngroups ||
      any(vapply(rows, is.null, logical(1L)))) {
    lav_msg_stop(gettextf(
      "information_meat_hc = %s is not available for the local standard
       errors in this setting (the casewise contributions to Gamma.eta
       are not available).", dQuote(hc, q = FALSE)))
  }
  tmp <- lav_model_nvcov_robust_sem(
    lavmodel = fit_pa@Model, lavsamplestats = fit_pa@SampleStats,
    lavcache = fit_pa@cache, lavdata = fit_pa@Data,
    lavimplied = fit_pa@implied, lavh1 = fit_pa@h1,
    lavoptions = fit_pa@Options, use_ginv = FALSE,
    attr_delta = TRUE, attr_e_inv = TRUE, attr_wls_v = TRUE)
  e_inv <- attr(tmp, "E.inv")
  delta <- attr(tmp, "Delta")
  wls_v <- attr(tmp, "WLS.V")
  ntotal <- fit_pa@SampleStats@ntotal
  nobs <- unlist(fit_pa@SampleStats@nobs)
  eqs <- lav_hc_equations(fit_pa@ParTable, ngroups = ngroups)
  unco_map <- lav_hc_unco_map(fit_pa@Model, fit_pa@ParTable)
  eq_classes <- lav_hc_eq_classes(fit_pa@ParTable)
  w <- vector("list", ngroups)
  for (g in seq_len(ngroups)) {
    if (ncol(rows[[g]]) != nrow(delta[[g]])) {
      lav_msg_stop(gettext(
        "internal error: the casewise contributions to Gamma.eta do not
         match the structural statistics (information_meat_hc)."))
    }
    fs <- lav_sam_hc_fs(step1 = step1, fit = fit, g = g)
    if (nrow(fs) != nrow(rows[[g]])) {
      lav_msg_stop(gettext(
        "internal error: case alignment failure between the factor scores
         and the casewise contributions to Gamma.eta (information_meat_hc)."))
    }
    ov_data <- fit@Data@X[[g]]
    colnames(ov_data) <- fit@Data@ov.names[[g]]
    w[[g]] <- lav_hc_weights(
      hc = hc, eqs = eqs[[g]], ov_data = ov_data, lv_scores = fs,
      n = nrow(fs), npar = ncol(e_inv), unco_map = unco_map,
      eq_classes = eq_classes
    )
  }
  adj <- lav_hc_sem_adjust(
    e_inv = e_inv, delta = delta, wls_v = wls_v, rows = rows, w = w,
    nobs = nobs, ntotal = ntotal,
    diag_wls_v = fit_pa@Model@estimator %in% c("DWLS", "ULS")
  )
  vcov_pa <- fit_pa@vcov$vcov
  if (is.null(vcov_pa)) {
    vcov_pa <- lav_sam_step2_se_vcov_pa(fit_pa)
  }
  vcov_pa <- vcov_pa + adj / ntotal
  if (length(step2_rm_idx) > 0L) {
    vcov_pa <- vcov_pa[-step2_rm_idx, -step2_rm_idx, drop = FALSE]
  }
  vcov_pa
}
