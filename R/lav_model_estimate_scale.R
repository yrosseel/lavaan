# parameter scaling for the optimizer (optim.parscale = "standardized")
#
# the quasi-Newton optimizers are not scale invariant: on a 'badly scaled'
# problem (a small change in one parameter changes the objective a lot,
# while a large change in another parameter hardly matters -- typically
# because the observed variables have very different variances) the
# search is slow, or fails. The remedy is a diagonal reparameterization:
# the optimizer works with u = D p, where p is the (packed) parameter
# vector and D = diag(scale) a positive diagonal matrix, chosen such that
# the elements of u are of comparable magnitude.
#
# lav_model_est_parscale() returns, for each free parameter, the factor
# that transforms it to its standardized-solution metric (as if the data
# were standardized): loadings by sd(lv)/sd(ov), regressions by
# sd(x)/sd(y), (co)variances by 1/(sd * sd), and so on. The variances of
# the latent variables are not known in advance; a latent variable
# inherits the scale of its marker indicator (the fixed-1 loading), which
# is the right order of magnitude (an upper bound); under std.lv = TRUE
# all latent variances are 1.
#
# lav_model_est_parscale_pack() maps these factors to the packed
# (optimizer) metric when linear equality constraints are handled by the
# null-space packing p = t(K) (x - k0).
#
# NOTE: the reparameterization is applied AFTER the equality-constraint
# packing, and undone BEFORE the (equality or inequality) constraint
# functions are evaluated (see lav_model_est()); the constraints are thus
# always imposed on the original parameters, whatever the scaling. (Up to
# 0.7-1 the free parameters were scaled before packing, which silently
# distorted all equality constraints except a == b.)

lav_model_est_parscale <- function(lavmodel = NULL,
                                   lavpartable = NULL,
                                   lavsamplestats = NULL,
                                   lavdata = NULL,
                                   lavh1 = NULL,
                                   lavoptions = NULL) {
  # only the parameter rows are needed (the constraint and := rows would
  # be 'standardized' too, which is pointless here, and fails when the
  # model's inequality function also holds the box-bound rows)
  con_idx <- which(lavpartable$op %in% c("==", "<", ">", ":="))
  if (length(con_idx) > 0L) {
    n_rows <- length(lavpartable$lhs)
    lavpartable <- lapply(lavpartable, function(col) {
      if (length(col) == n_rows) col[-con_idx] else col
    })
  }

  # observed variances
  if (lavdata@nlevels > 1L) {
    if (length(lavh1) > 0L) {
      ov_var <- lapply(lavh1$implied$cov, diag)
    } else {
      ov_var <- lapply(
        do.call(c, lapply(lavdata@Lp, "[[", "ov.idx")),
        function(x) rep(1, length(x))
      )
    }
  } else {
    if (lavoptions$conditional.x) {
      ov_var <- lavsamplestats@res.var
    } else {
      ov_var <- lavsamplestats@var
    }
  }

  if (lavoptions$std.lv) {
    parscale <- lav_standardize_all(
      lavobject = NULL,
      est = rep(1, length(lavpartable$lhs)),
      est_std = rep(1, length(lavpartable$lhs)),
      cov_std = FALSE, ov_var = ov_var,
      lavmodel = lavmodel, lavpartable = lavpartable,
      cov_x = lavsamplestats@cov.x
    )
  } else {
    # latent variances: under the marker (fixed-1 loading) convention, a
    # latent variable inherits the scale of its marker indicator, so the
    # marker's observed variance is the right order of magnitude;
    # higher-order factors follow their marker chain; a latent variable
    # without a marker falls back to 1.0
    lv_var <- vector("list", lavmodel@nblocks)
    block_values <- lav_pt_block_values(lavpartable)
    for (b in seq_len(lavmodel@nblocks)) {
      mm_in_block <- 1:lavmodel@nmat[b] + cumsum(c(0, lavmodel@nmat))[b]
      mm_lambda <- lavmodel@GLIST[mm_in_block]$lambda
      n_lv <- ncol(mm_lambda)
      lv_var[[b]] <- rep(1.0, n_lv)
      lv_names_b <- unique(unlist(lav_pt_vnames(lavpartable, "lv",
        block = block_values[b])))
      ov_names_b <- unique(unlist(lav_pt_vnames(lavpartable, "ov",
        block = block_values[b])))
      if (length(lv_names_b) == 0L || n_lv == 0L) {
        next
      }
      s2 <- setNames(rep(NA_real_, length(lv_names_b)), lv_names_b)
      # a few passes to resolve higher-order marker chains
      for (rep_i in 1:4) {
        for (l in lv_names_b) {
          if (!is.na(s2[[l]])) next
          m_idx <- which(lavpartable$op == "=~" &
            lavpartable$lhs == l &
            lavpartable$block == block_values[b] &
            lavpartable$free == 0L &
            !is.na(lavpartable$ustart) &
            lavpartable$ustart != 0)
          if (length(m_idx) == 0L) next
          mk <- lavpartable$rhs[m_idx[1]]
          if (mk %in% ov_names_b) {
            pos <- match(mk, ov_names_b)
            if (!is.na(pos) && pos <= length(ov_var[[b]])) {
              s2[[l]] <- ov_var[[b]][pos] / lavpartable$ustart[m_idx[1]]^2
            }
          } else if (mk %in% lv_names_b && !is.na(s2[[mk]])) {
            s2[[l]] <- s2[[mk]] / lavpartable$ustart[m_idx[1]]^2
          }
        }
        if (!anyNA(s2)) break
      }
      s2[is.na(s2)] <- 1.0
      # the columns of LAMBDA follow the lv names of the block (the model
      # matrices carry no dimnames)
      if (length(lv_names_b) == n_lv) {
        lv_var[[b]] <- as.numeric(s2)
      }
    }

    parscale <- lav_standardize_all(
      lavobject = NULL,
      est = rep(1, length(lavpartable$lhs)),
      cov_std = FALSE, ov_var = ov_var, lv_var = lv_var,
      lavmodel = lavmodel, lavpartable = lavpartable,
      cov_x = lavsamplestats@cov.x
    )
  }

  # free parameters only (one representative per set of simple equality
  # constraints)
  if (lavmodel@ceq.simple.only) {
    parscale <- parscale[lavpartable$free > 0L &
      !duplicated(lavpartable$free)]
  } else {
    parscale <- parscale[lavpartable$free > 0L]
  }

  # sanity: non-finite or non-positive factors fall back to 1.0
  bad <- !is.finite(parscale) | parscale <= 0
  if (any(bad)) {
    parscale[bad] <- 1.0
  }

  parscale
}

# map the per-parameter factors to the packed metric: if p = t(K) (x - k0)
# (null-space packing of linear equality constraints), the magnitude of
# p_i is sqrt(sum_j K_ji^2 / s_j^2) when x_j has magnitude 1/s_j
lav_model_est_parscale_pack <- function(parscale = NULL, lavmodel = NULL) {
  if (lavmodel@eq.constraints && ncol(lavmodel@eq.constraints.K) > 0L &&
      nrow(lavmodel@eq.constraints.K) == length(parscale)) {
    k_pack <- lavmodel@eq.constraints.K
    scale <- 1 / sqrt(colSums(k_pack^2 / parscale^2))
    scale[!is.finite(scale) | scale <= 0] <- 1.0
  } else {
    scale <- parscale
  }
  scale
}
