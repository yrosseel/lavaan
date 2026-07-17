# shared building blocks for the test-statistic machinery (lav_test*.R)
#
# YR 11 July 2026: initial version, collecting code that used to be
#                  duplicated across lav_test.R, lav_test_satorra_bentler.R,
#                  lav_test_yuan_bentler.R, lav_test_browne.R and
#                  lav_test_diff.R

# degrees of freedom of the user model, adding back the rank of the
# equality-constraint jacobian (inequality constraints, active or not,
# are ignored)
lav_test_df <- function(lavpartable = NULL, lavmodel = NULL) {
  df <- lav_pt_df(lavpartable)

  if (!lavmodel@cin.simple.only && nrow(lavmodel@con.jac) > 0L) {
    ceq_idx <- attr(lavmodel@con.jac, "ceq.idx")
    if (length(ceq_idx) > 0L) {
      neq <- qr(lavmodel@con.jac[ceq_idx, , drop = FALSE])$rank
      df <- df + neq
    }
  } else if (lavmodel@ceq.simple.only) {
    ndat <- lav_pt_ndat(lavpartable)
    npar <- max(lavpartable$free)
    df <- ndat - npar
  }

  df
}

# the meta-information attribute attached to the TEST list; used by
# lav_test_print() for independent printing
lav_test_info_attr <- function(lavdata = NULL, lavoptions = NULL) {
  list(
    ngroups = lavdata@ngroups, group.label = lavdata@group.label,
    information = lavoptions$information,
    h1.information = lavoptions$h1.information,
    observed.information = lavoptions$observed.information
  )
}

# do the test-specific information options (second element) match the ones
# used for the standard errors (first element)? if not, E.inv (and Delta,
# WLS.V) must be recomputed with the test-specific options; returns the
# lavoptions with the first elements overwritten by the second, and the
# recompute flag as an attribute
lav_test_change_options <- function(lavoptions = NULL) {
  if (length(lavoptions$information) == 1L &&
    length(lavoptions$h1.information) == 1L &&
    length(lavoptions$observed.information) == 1L) {
    recompute <- FALSE
  } else if (
    (lavoptions$information[1] == lavoptions$information[2]) &&
      (lavoptions$h1.information[1] == lavoptions$h1.information[2]) &&
      (lavoptions$information[2] == "expected" ||
        (lavoptions$observed.information[1] ==
          lavoptions$observed.information[2]))) {
    recompute <- FALSE
  } else {
    recompute <- TRUE
    # change information options
    lavoptions$information[1] <- lavoptions$information[2]
    lavoptions$h1.information[1] <- lavoptions$h1.information[2]
    lavoptions$observed.information[1] <- lavoptions$observed.information[2]
  }

  attr(lavoptions, "recompute") <- recompute
  lavoptions
}

# a 'failed' scaled-test entry: the standard test augmented with NA
# scaling quantities; used when the information matrix cannot be inverted
# (na_standard = TRUE) or when df <= 0 (na_standard = FALSE, no shift)
lav_test_scaled_na <- function(standard = NULL, test1 = NULL,
                               ngroups = 1L,
                               na_standard = TRUE, shift = TRUE) {
  if (na_standard) {
    standard$stat <- as.numeric(NA)
    standard$stat.group <- rep(as.numeric(NA), ngroups)
    standard$pvalue <- as.numeric(NA)
  }
  if (shift) {
    out <- c(standard,
      scaling.factor = as.numeric(NA),
      shift.parameter = as.numeric(NA),
      label = character(0)
    )
  } else {
    out <- c(standard,
      scaling.factor = as.numeric(NA),
      label = character(0)
    )
  }
  # to prevent lavTestLRT error when a robust test is detected for some
  # but not all models
  out$test <- test1

  list(standard = standard, failed = out)
}

# label for the scaled test statistics; the Mplus-flavored labels are used
# when information.expected.mplus is TRUE
lav_test_scaled_label <- function(test = "satorra.bentler",
                                  estimator = "ML",
                                  mplus_flag = FALSE) {
  if (test == "satorra.bentler") {
    plain <- "Satorra-Bentler correction"
    mplus <- c(
      ML = "Satorra-Bentler correction (Mplus variant)",
      DWLS = "Satorra-Bentler correction (WLSM)",
      ULS = "Satorra-Bentler correction (ULSM)"
    )
  } else if (test == "mean.var.adjusted") {
    plain <- "mean and variance adjusted correction"
    mplus <- c(
      ML = "mean and variance adjusted correction (MLMV)",
      DWLS = "mean and variance adjusted correction (WLSMV)",
      ULS = "mean and variance adjusted correction (ULSMV)"
    )
  } else if (test == "scaled.shifted") {
    plain <- "simple second-order correction"
    mplus <- c(
      ML = "simple second-order correction (MLMV)",
      DWLS = "simple second-order correction (WLSMV)",
      ULS = "simple second-order correction (ULSMV)"
    )
  } else {
    return(character(0L))
  }

  if (mplus_flag && estimator %in% names(mplus)) {
    unname(mplus[estimator])
  } else {
    plain
  }
}

# stack the per-group ingredients for the 'assembled' (explicit) UGamma
# computation: fg-weighted V blocks, Gamma_g / fg blocks (the
# Cov(sqrt(ntotal) s_g) stacking convention -- see the SCALING CONVENTIONS
# note in lav_samplestats_gamma.R), and the row-stacked Delta;
# a caller-supplied gamma_full (possibly with cross-group blocks) is used
# as-is instead of the block-diagonal assembly
lav_test_ug_stack <- function(m_gamma = NULL, delta = NULL, wls_v = NULL,
                              nobs = NULL, ntotal = NULL,
                              gamma_full = NULL) {
  ngroups <- length(wls_v)
  fg <- unlist(nobs) / ntotal

  v_g <- vector("list", ngroups)
  for (g in seq_len(ngroups)) {
    if (is.matrix(wls_v[[g]])) {
      v_g[[g]] <- fg[g] * wls_v[[g]]
    } else {
      v_g[[g]] <- fg[g] * diag(wls_v[[g]], nrow = length(wls_v[[g]]))
    }
  }

  if (is.null(gamma_full)) {
    gamma_f <- lav_gamma_rescale_ntotal(m_gamma,
      nobs = nobs, ntotal = ntotal
    )
    gamma_all <- lav_mat_bdiag(gamma_f)
  } else {
    gamma_all <- gamma_full
  }

  list(
    v_all = lav_mat_bdiag(v_g),
    gamma_all = gamma_all,
    delta_all = do.call("rbind", delta)
  )
}

# streamed computation of tr(U Gamma) and tr((U Gamma)^2), with
#
#   U = V - V Delta E.inv Delta' V
#
# assembled over groups in the Cov(sqrt(ntotal) s_g) stacking convention
# (fg-weighted V blocks, Gamma_g / fg blocks), WITHOUT ever forming a
# pstar x pstar product. With WD_g = V_g Delta_g, GWD_g = Gamma_g WD_g and
# B_g = WD_g' GWD_g (= Delta' V Gamma V Delta, npar x npar):
#
#   tr(U Gamma)     = sum_g [ tr(V_g Gamma_g) - fg_g tr(E.inv B_g) ]
#   tr((U Gamma)^2) = t1 - 2 t2 + t3, with
#      t1 = sum_g tr((V_g Gamma_g)^2)
#      t2 = tr(E.inv sum_g fg_g GWD_g' V_g GWD_g)
#      t3 = tr((E.inv C)^2),  C = sum_g fg_g B_g
#
# (the cross-group blocks of U Gamma enter tr((U Gamma)^2) only through
# E.inv, which is what the C accumulation captures)
#
# the V blocks may be given as vectors (diagonal weight matrices, e.g.
# DWLS/ULS) or full matrices; this matters in the categorical setting and
# whenever the number of variables grows, where pstar is large
#
# old_approach = TRUE reproduces the < 0.6-13 behavior where also
# tr((U Gamma)^2) was accumulated per group:
#   tr((U_g Gamma_g)^2) = t1_g - 2 fg_g tr(E.inv T2_g)
#                              + fg_g^2 tr((E.inv B_g)^2)
lav_test_ug_trace_stream <- function(m_gamma = NULL, delta = NULL,
                                     wls_v = NULL, e_inv = NULL,
                                     fg = NULL,
                                     satterthwaite = FALSE,
                                     old_approach = FALSE) {
  ngroups <- length(m_gamma)
  npar <- NCOL(delta[[1]])

  trace_group <- numeric(ngroups)
  trace2_group <- numeric(ngroups) # old_approach only
  t1 <- 0
  m_t2 <- m_c <- matrix(0, npar, npar)

  for (g in seq_len(ngroups)) {
    v <- wls_v[[g]]
    gamma_g <- m_gamma[[g]]
    delta_g <- delta[[g]]

    if (is.matrix(v)) {
      wd <- v %*% delta_g
      trv <- sum(v * gamma_g) # both symmetric
    } else {
      wd <- v * delta_g
      trv <- sum(v * diag(gamma_g))
    }
    gwd <- gamma_g %*% wd
    b_g <- crossprod(wd, gwd)
    trace_group[g] <- trv - fg[g] * sum(e_inv * b_g)

    if (satterthwaite) {
      if (is.matrix(v)) {
        vg <- v %*% gamma_g
        t1_g <- sum(vg * t(vg))
        t2_g <- crossprod(gwd, v %*% gwd)
      } else {
        t1_g <- drop(crossprod(v, (gamma_g * gamma_g) %*% v))
        t2_g <- crossprod(gwd, v * gwd)
      }
      if (old_approach) {
        eb <- e_inv %*% b_g
        trace2_group[g] <- t1_g - 2 * fg[g] * sum(e_inv * t2_g) +
          fg[g] * fg[g] * sum(eb * t(eb))
      } else {
        t1 <- t1 + t1_g
        m_t2 <- m_t2 + fg[g] * t2_g
        m_c <- m_c + fg[g] * b_g
      }
    }
  } # g

  trace_ugamma <- sum(trace_group)
  trace_ugamma2 <- as.numeric(NA)
  if (satterthwaite) {
    if (old_approach) {
      trace_ugamma2 <- sum(trace2_group)
    } else {
      m_ec <- e_inv %*% m_c
      trace_ugamma2 <- t1 - 2 * sum(e_inv * m_t2) + sum(m_ec * t(m_ec))
    }
  }

  list(
    trace.UGamma = trace_ugamma,
    trace.UGamma2 = trace_ugamma2,
    trace.UGamma.group = trace_group
  )
}

# is this test statistic available? the @test slot may contain a
# placeholder entry with a NA statistic: e.g. the 'standard' test when
# estimator = "IV", where no standard chi-square is defined
lav_test_is_available <- function(test_element = NULL) {
  !is.null(test_element$stat) &&
    length(test_element$stat) > 0L &&
    !is.na(test_element$stat[1])
}

# which test statistic(s) to return when lavTest() is called without a
# 'test' argument? normally the 'standard' one, but if that one is not
# available (see above), fall back to the test statistic(s) that are
# (e.g. Browne's residual based test when estimator = "IV")
lav_test_default_names <- function(test_1 = NULL) {
  test_names <- names(test_1)

  # no names at all (e.g. test = "none") -> nothing to fall back to
  if (is.null(test_names)) {
    return("standard")
  }

  standard_idx <- which(test_names == "standard")
  if (length(standard_idx) > 0L &&
      lav_test_is_available(test_1[[standard_idx[1]]])) {
    return("standard")
  }

  ok_idx <- which(vapply(test_1, lav_test_is_available, logical(1L)))
  if (length(ok_idx) == 0L) {
    return("standard")
  }

  test_names[ok_idx]
}

# which (non-scaled) test statistic shall we scale? (the scaled.test
# option; e.g. scaled.test = "browne.residual.nt.model" gives the
# RLS-based scaled tests)
lav_test_scaled_base <- function(test_1 = NULL, lavoptions = NULL,
                                 test = NULL) {
  unscaled_test <- test_1[[1]]
  if (lavoptions$scaled.test != "standard") {
    idx <- which(names(test_1) == lavoptions$scaled.test)
    if (length(idx) > 0L) {
      unscaled_test <- test_1[[idx[1]]]
    } else {
      lav_msg_warn(gettextf(
        "scaled.test [%1$s] not found among available (non scaled) tests:
        %2$s. Using standard test instead.",
        lavoptions$scaled.test, lav_msg_view(test)))
    }
  }
  unscaled_test
}
