# Standardize parameter estimates of a fitted lavaan model.
#
# Three "modes" of standardization are supported, all implemented on top of a
# single observed-variable engine (lav_standardize_all):
#   - std.lv  : only the latent variables are standardized
#   - std.all : latent *and* observed variables are standardized
#   - std.nox : like std.all, but the (observed) exogenous 'x' variables are
#               left unstandardized
# std.nox is in fact a special case of a more general mechanism: standardize
# all observed variables *except* a given set (here: the exogenous 'x'). The
# 'ov_std' argument of lav_standardize_all() exposes this generality, so that a
# caller can request that only a user-specified subset of observed variables is
# standardized.

# ---------------------------------------------------------------------------
# internal helpers
# ---------------------------------------------------------------------------

# Resolve the standard set of ingredients (lavmodel, lavpartable, partable,
# est, glist, cov_x) from either a fitted lavobject or the individual pieces.
lav_standardize_setup <- function(lavobject = NULL, lavmodel = NULL,
                                  lavpartable = NULL, partable = NULL,
                                  est = NULL, glist = NULL, cov_x = NULL,
                                  need_cov_x = FALSE) {
  if (is.null(lavobject)) {
    stopifnot(!is.null(lavmodel), !is.null(lavpartable))
    if (is.null(est)) {
      if (!is.null(lavpartable$est)) {
        est <- lavpartable$est
                # if this changes, tag @TDJorgensen in commit message
      } else {
        lav_msg_stop(gettext("could not find `est' in lavpartable"))
      }
    }
  } else {
    lavmodel <- lavobject@Model
    lavpartable <- lavobject@ParTable
    if (is.null(est)) {
      est <- lav_inspect_est(lavobject)
    }
    if (need_cov_x && lavmodel@conditional.x && is.null(cov_x)) {
      if (!is.null(lavobject@implied$cov.x[[1]])) {
        cov_x <- lavobject@implied$cov.x
                # if this changes, tag @TDJorgensen in commit message
      } else {
        # perhaps a lavaanList object: extract it from GLIST per block
        cov_x <- vector("list", length = lavmodel@nblocks)
        for (b in seq_len(lavmodel@nblocks)) {
          mm_in_block <- (seq_len(lavmodel@nmat[b]) +
            cumsum(c(0, lavmodel@nmat))[b])
          mlist <- lavmodel@GLIST[mm_in_block]
          cov_x[[b]] <- mlist[["cov.x"]]
        }
      }
    }
  }

  if (is.null(partable)) {
    partable <- lavpartable
  }
  if (is.null(glist)) {
    glist <- lavmodel@GLIST
  }

  list(lavmodel = lavmodel, lavpartable = lavpartable, partable = partable,
       est = est, glist = glist, cov_x = cov_x)
}

# Evaluate the standardized values of the defined (:=), equality (==) and
# inequality (< >) rows, given the (already standardized) 'out' vector.
lav_standardize_constraints <- function(out, partable, lavmodel) {
  x <- out[partable$free & !duplicated(partable$free)]

  idx <- which(partable$op == ":=")
  if (length(idx) > 0L) {
    out[idx] <- lavmodel@def.function(x)
  }
  idx <- which(partable$op == "==")
  if (length(idx) > 0L) {
    out[idx] <- lavmodel@ceq.function(x)
  }
  idx <- which(partable$op == "<" | partable$op == ">")
  if (length(idx) > 0L) {
    out[idx] <- lavmodel@cin.function(x)
  }

  out
}

# Observed-variable standard deviations (and variances) for block 'g'.
# Returns the (possibly conditional.x-extended) ov.names together with the
# vectors of SDs ('ov') and variances ('ov2').
lav_standardize_ov_sd <- function(lavpartable, g, vy = NULL, ov_var = NULL,
                                  conditional_x = FALSE, cov_x = NULL) {
  ov_names <- lav_pt_vnames(lavpartable, "ov", block = g) # not user

  if (is.null(ov_var)) {
    ov2 <- vy[[g]]
    # replace zero values by NA (but keep negative values)
    zero_idx <- which(abs(ov2) < .Machine$double.eps)
    if (length(zero_idx) > 0L) {
      ov2[zero_idx] <- as.numeric(NA)
    }
    # replace negative values by NA (for sqrt)
    tmp_ov2 <- ov2
    neg_idx <- which(tmp_ov2 < 0)
    if (length(neg_idx) > 0L) {
      tmp_ov2[neg_idx] <- as.numeric(NA)
    }
    ov <- sqrt(tmp_ov2)
  } else {
    ov2 <- ov_var[[g]]
    ov <- sqrt(ov2)
  }

  if (conditional_x) {
    # extend OV with ov.names.x
    ov_names_x <- lav_pt_vnames(lavpartable, "ov.x", block = g)
    ov_names_nox <- lav_pt_vnames(lavpartable, "ov.nox", block = g)
    ov_names <- c(ov_names_nox, ov_names_x)
    ov2 <- c(ov2, diag(cov_x[[g]]))
    ov <- c(ov, sqrt(diag(cov_x[[g]])))
  }

  list(ov_names = ov_names, ov = ov, ov2 = ov2)
}

# Given a parameter vector 'x', set it as the model parameters and return the
# (glist, est) pair to feed into the standardization engines. Handles the
# (EFA) rotation case. Used by the bootstrap/Monte-Carlo machinery.
lav_standardize_prepare_x <- function(x, lavobject, rotation = FALSE) {
  lavmodel <- lav_model_set_parameters(lavmodel = lavobject@Model, x = x)

  if (rotation) {
    lavmodel@GLIST <- lavTech(lavobject, "est.unrotated") # unrotated!
    est_rot <- lav_model_efa_rotate_x(
      x = x,
      lavmodel = lavmodel, # unrotated!
      lavoptions = lavobject@Options,
      init_rot = lavmodel@H,
      type = "user",
      extra = TRUE
    )
    glist <- attr(est_rot, "extra")$glist
    attributes(est_rot) <- NULL
    est <- est_rot
  } else {
    glist <- lavmodel@GLIST
                    # if this changes, tag @TDJorgensen in commit message
    est <- lav_model_get_parameters(lavmodel, type = "user")
  }

  list(glist = glist, est = est)
}

# ---------------------------------------------------------------------------
# 'x' wrappers (parameter vector -> standardized values), used for bootstrap
# and Monte-Carlo standard errors.
# ---------------------------------------------------------------------------

lav_standardize_lv_x <- function(x, lavobject, partable = NULL, cov_std = TRUE,
                                 lv_var = NULL, rotation = FALSE) {
  ge <- lav_standardize_prepare_x(x, lavobject, rotation = rotation)
  lav_standardize_lv(
    lavobject = lavobject, partable = partable, est = ge$est,
    glist = ge$glist, cov_std = cov_std, lv_var = lv_var
  )
}

lav_standardize_all_x <- function(x, lavobject, partable = NULL, cov_std = TRUE,
                                  ov_std = NULL, nox = FALSE, rotation = FALSE) {
  ge <- lav_standardize_prepare_x(x, lavobject, rotation = rotation)
  lav_standardize_all(
    lavobject = lavobject, partable = partable, est = ge$est,
    est_std = NULL, glist = ge$glist, cov_std = cov_std,
    ov_std = ov_std, nox = nox
  )
}

lav_standardize_all_nox_x <- function(x, lavobject, partable = NULL,
                                      cov_std = TRUE, rotation = FALSE) {
  lav_standardize_all_x(
    x = x, lavobject = lavobject, partable = partable, cov_std = cov_std,
    nox = TRUE, rotation = rotation
  )
}

# ---------------------------------------------------------------------------
# std.lv : standardize the latent variables only
# ---------------------------------------------------------------------------

lav_standardize_lv <- function(lavobject = NULL,
                               partable = NULL, est = NULL, glist = NULL,
                               cov_std = TRUE, lv_var = NULL,
                               lavmodel = NULL, lavpartable = NULL) {
  s <- lav_standardize_setup(
    lavobject = lavobject, lavmodel = lavmodel, lavpartable = lavpartable,
    partable = partable, est = est, glist = glist
  )
  lavmodel <- s$lavmodel
  lavpartable <- s$lavpartable
  partable <- s$partable
  est <- s$est
  glist <- s$glist

  out <- est
  n <- length(est)
  stopifnot(n == length(partable$lhs))

  # working copy of est: for interaction terms (e.g., X:Z), the (co)variance
  # and mean rows are 'centered' below (the mean-driven part of the product
  # moments is removed) before any scaling takes place. The cov_std block
  # further down needs these centered -- but still unscaled -- variances to
  # compute its sqrt(variance) divisors, so we keep them in est_c.
  est_c <- est

  nmat <- lavmodel@nmat
  mm_idx <- lav_model_group_mm_indices(nmat)

  # compute ETA
  if (is.null(lv_var)) {
    lv_eta <- lav_model_veta(
      lavmodel = lavmodel,
      glist = glist
    )

    lv_eeta <- lav_model_eeta(
      lavmodel = lavmodel,
      glist = glist,
      lavsamplestats = lavobject@SampleStats
    )
  }

  for (g in 1:lavmodel@nblocks) {
    ov_names <- lav_pt_vnames(lavpartable, "ov", block = g) # not user,
    # which may be incomplete
    lv_names <- lav_pt_vnames(lavpartable, "lv", block = g)

    # shortcut: no latents in this block, nothing to do
    if (length(lv_names) == 0L) {
      next
    }

    if (is.null(lv_var)) {
      veta_g <- lv_eta[[g]] # full model-implied V(ETA); only kept as a
                            # psi fallback for the centering below (RAM)
      eta2 <- diag(veta_g)
      eeta <- lv_eeta[[g]]
    } else {
      veta_g <- NULL
      eta2 <- lv_var[[g]]
      eeta <- numeric(length(eta2))
    }
    # change negative values to NA
    eta2[eta2 < 0] <- as.numeric(NA)
    eta <- sqrt(eta2)

    # Interaction/quadratic term correction (FV)
    # (based on Kelava & Brandt, 2022; Brandt et al., 2015)
    # For interaction terms A:B, the standardized coefficient should use
    # SD(A)*SD(B) instead of SD(A:B) (Eq. 13 in Brandt et al., 2015).
    lv_int_names <- lav_pt_vnames(lavpartable, "lv.interaction",
      block = g
    )

    if (length(lv_int_names) > 0L) {

      # component indices of each interaction term (int = A:B, a = A,
      # b = B); collected for the centering of the interaction-term
      # moments below
      int_names_ok <- character(0L)
      int_comp <- list()

      for (int_name in lv_int_names) {
        components <- strsplit(int_name, ":", fixed = TRUE)[[1L]]

        # only two-way terms: a three-way (or higher) product would need
        # its own decomposition, and the product of *three* SDs below
        if (length(components) != 2L) next

        a <- components[1]
        b <- components[2]

        idx_int <- match(int_name, lv_names)
        idx_a   <- match(a, lv_names)
        idx_b   <- match(b, lv_names)

        # Need a, b and interaction term before moving on
        if (is.na(idx_a) || is.na(idx_b) || is.na(idx_int)) next

        int_names_ok <- c(int_names_ok, int_name)
        int_comp <- c(int_comp,
          list(list(int = idx_int, a = idx_a, b = idx_b)))

        exp_a <- eeta[idx_a]
        exp_b <- eeta[idx_b]

        # Do we need to shift simple main effects?
        if (exp_a != 0 || exp_b != 0) {
          # KS, 08/05/26
          #
          # This covers the general case of y ~ x + z + w + x:z + x:x
          # This should naturally extend to other conditions with more
          # interaction terms. e.g., y ~ x + z + w + x:z + x:x + w:z
          #
          # For y = b0 + b1*x + b2*z + b3*xz + b4x^2
          # b1 and b2 are the slope at the x=0 and z=0. When mean-centering
          # we have to re-define b1 and b2 to the expected slopes
          # (i.e., the partial derivatives) at x=mean(x), z=mean(z).
          #
          # Computing the partial derivatives we get
          #
          #   dy/dx = b1 + b3*z + 2*b4*x
          #   dy/dz = b2 + b3*x
          #
          # With x=mean(x) and z=mean(x), we get
          #
          #   dy/dx = b1 + b3 * mean(z) + 2 * b4 * mean(x)
          #   dy/dz = b2 + b3 * mean(x)
          #
          # Since we iteratively update out[<idx>] for b1 and b2, we carry
          # out all the necessary conditions, by the end of the loop.
          #
          # Iteration: 1, interaction term: x:z
          #   b1 <- b1 + b3 * mean(z)
          #   b2 <- b2 + b3 * mean(x)
          #
          # Iteration: 2, quadratic term: x:x
          #   b1 <- b1 + b4 * mean(x)
          #   b1 <- b1 + b4 * mean(x)
          #
          #   The above is thus equivalent to b1 + 2 * b4 * mean(x)
          #
          # Putting it all together we get
          #
          #   b1 <- b1 + b3 * mean(z) + 2 * b4 * mean(x)
          #   b2 <- b2 + b3 * mean(x)
          #

          # Find dependent variables
          lv_int_dep <- unique(partable$lhs[
            partable$op == "~" & partable$rhs == int_name &
            partable$block == g
          ])

          for (dep in lv_int_dep) {
            idx_beta_a <- which(
              partable$lhs == dep & partable$op == "~" &
              partable$rhs == a & partable$block == g
            )

            idx_beta_b <- which(
              partable$lhs == dep & partable$op == "~" &
              partable$rhs == b & partable$block == g
            )

            # read from est (not partable$est): in the *_x wrappers (delta
            # method, bootstrap) est is the perturbed parameter vector, and
            # the gradient of the shifted main effects wrt beta_ab must not
            # be lost
            beta_ab <- est[
              partable$lhs == dep & partable$op == "~" &
              partable$rhs == int_name & partable$block == g
            ]

            out[idx_beta_a] <- out[idx_beta_a] + beta_ab * exp_b
            out[idx_beta_b] <- out[idx_beta_b] + beta_ab * exp_a
          }
        }

        # Replace ETA for interaction terms with SD(A)*SD(B)
        eta[idx_int] <- eta[idx_a] * eta[idx_b]
      }

      # FV, 05/07/26
      #
      # Center the (co)variances and means of the interaction terms:
      # their psi/alpha elements describe the *uncentered* product A:B,
      # whereas the standardized solution treats the interaction as the
      # product of the *centered* (standardized) components. Without this
      # correction the standardized moments of the product terms depend
      # on how the mean structure is parameterized.
      #
      # No distributional assumptions are needed: with ac = A - E(A) and
      # bc = B - E(B), we have the exact identity
      #
      # ac*bc = A:B - E(B)*A - E(A)*B + E(A)*E(B)
      #
      # so the centered product is a linear combination of variables that
      # are already elements of eta. Collect the weights in a matrix K
      # (one row per variable; regular lvs keep their identity row), e.g.
      # for eta = (X, Z, X:Z, X:X) with mx = E(X), mz = E(Z):
      #
      # row X:Z: (-mz, -mx, 1, 0)     row X:X: (-2*mx, 0, 0, 1)
      #
      # (a quadratic term accumulates -mx twice in the column of X). Then
      # K %*% PSI %*% t(K) holds all centered (co)variances at once
      # (constants do not affect them), and the centered means follow
      # from E(ac*bc) = E(A:B) - E(A)*E(B).
      #
      # K is applied to PSI (not VETA), because the "~~" rows below *are*
      # psi elements: for exogenous pairs (the regular case) PSI equals
      # VETA, and for a residual covariance with an endogenous lv,
      # PSI[X, Y] = Cov(X, zeta_Y) is exactly what the centering needs.
      # This must happen *before* the generic scaling below, which then
      # divides the centered moments by SD(A)*SD(B) (see eta[idx_int]).
      if (length(int_comp) > 0L && !is.null(veta_g) && any(eeta != 0)) {
        psi_g <- glist[mm_idx[[g]]]$psi
        if (is.null(psi_g)) { # e.g., RAM representation
          psi_g <- veta_g
        }
        # build K and the correction matrix: shift_mat holds, for every pair
        # of latent variables, the mean-driven part of their (co)variance
        # (zero unless a product term is involved)
        kmat <- diag(nrow(psi_g))
        for (p in int_comp) {
          kmat[p$int, p$a] <- kmat[p$int, p$a] - eeta[p$b]
          kmat[p$int, p$b] <- kmat[p$int, p$b] - eeta[p$a]
        }
        shift_mat <- psi_g - kmat %*% psi_g %*% t(kmat)

        # "~~" rows involving at least one product term: this covers the
        # variance of a product term (X:Z ~~ X:Z), its covariance with a
        # regular lv (X:Z ~~ V), and product-product covariances
        # (X:Z ~~ X:X); rows whose other side is not a latent variable are
        # left alone (keep)
        idx_cov <- which(partable$op == "~~" & partable$block == g &
          (partable$lhs %in% int_names_ok |
           partable$rhs %in% int_names_ok))
        l_idx <- match(partable$lhs[idx_cov], lv_names)
        r_idx <- match(partable$rhs[idx_cov], lv_names)
        keep <- which(!is.na(l_idx) & !is.na(r_idx))
        shift <- shift_mat[cbind(l_idx[keep], r_idx[keep])]
        # est_c is also corrected: the cov_std block below derives its
        # sqrt(variance) divisors from est_c (not out, which is scaled by
        # then)
        out[idx_cov[keep]] <- out[idx_cov[keep]] - shift
        est_c[idx_cov[keep]] <- est_c[idx_cov[keep]] - shift

        # "~1" rows of the product terms: E(ac*bc) = E(A:B) - E(A)*E(B);
        # after the generic scaling this equals, e.g., cor(X, Z) for X:Z
        idx_mean <- which(partable$op == "~1" & partable$block == g &
          partable$lhs %in% int_names_ok)
        for (row_i in idx_mean) {
          p <- int_comp[[match(partable$lhs[row_i], int_names_ok)]]
          out[row_i] <- out[row_i] - eeta[p$a] * eeta[p$b]
        }
      }
    }

    # 1a. "=~" regular indicators
    idx <- which(partable$op == "=~" & !(partable$rhs %in% lv_names) &
      partable$block == g)
    out[idx] <- out[idx] * eta[match(partable$lhs[idx], lv_names)]

    # 1b. "=~" regular higher-order lv indicators
    idx <- which(partable$op == "=~" & !(partable$rhs %in% ov_names) &
      partable$block == g)
    out[idx] <- (out[idx] * eta[match(partable$lhs[idx], lv_names)]
      / eta[match(partable$rhs[idx], lv_names)])

    # 2. "~" regressions (and "<~")
    idx <- which((partable$op == "~" | partable$op == "<~") &
      partable$lhs %in% lv_names &
      partable$block == g)
    out[idx] <- out[idx] / eta[match(partable$lhs[idx], lv_names)]

    idx <- which((partable$op == "~" | partable$op == "<~") &
      partable$rhs %in% lv_names &
      partable$block == g)
    out[idx] <- out[idx] * eta[match(partable$rhs[idx], lv_names)]

    # 3b. "~~" lv variances
    # ATTENTION: in Mplus 4.1, the off-diagonal residual covariances
    #            were computed by the formula cov(i,j) / sqrt(i.var*j.var)
    #            where i.var and j.var were diagonal elements of ETA
    #
    #            in Mplus 6.1 (but also AMOS and EQS), the i.var and j.var
    #            elements are the 'PSI' diagonal elements!!
    rv_idx <- which(partable$op == "~~" & partable$rhs %in% lv_names &
      partable$lhs == partable$rhs &
      partable$block == g)
    out[rv_idx] <- (out[rv_idx] / eta[match(partable$lhs[rv_idx], lv_names)]
      / eta[match(partable$rhs[rv_idx], lv_names)])

    # covariances lv
    # three types:
    # - only lhs is LV (and fixed.x = FALSE)
    # - only rhs is LV (and fixed.x = FALSE)
    # - both lhs and rhs are LV (regular case)
    if (cov_std) {
      # est_c: est with the interaction-term variances centered (see above)
      if (!is.complex(est_c[rv_idx])) {
        rv <- sqrt(abs(est_c[rv_idx])) # abs in case of heywood cases
      } else {
        rv <- sqrt(est_c[rv_idx])
      }
      rv_names <- partable$lhs[rv_idx]
    }

    # left
    idx_lhs <- which(partable$op == "~~" &
      partable$lhs %in% lv_names &
      partable$lhs != partable$rhs &
      partable$block == g)
    if (length(idx_lhs) > 0L) {
      if (cov_std == FALSE) {
        out[idx_lhs] <-
          (out[idx_lhs] / eta[match(partable$lhs[idx_lhs], lv_names)])
      } else {
        out[idx_lhs] <-
          (out[idx_lhs] / rv[match(partable$lhs[idx_lhs], rv_names)])
      }
    }

    # right
    idx_rhs <- which(partable$op == "~~" &
      partable$rhs %in% lv_names &
      partable$lhs != partable$rhs &
      partable$block == g)
    if (length(idx_rhs) > 0L) {
      if (cov_std == FALSE) {
        out[idx_rhs] <-
          (out[idx_rhs] / eta[match(partable$rhs[idx_rhs], lv_names)])
      } else {
        out[idx_rhs] <-
          (out[idx_rhs] / rv[match(partable$rhs[idx_rhs], rv_names)])
      }
    }

    # 4b. "~1" lv
    idx <- which(partable$op == "~1" & partable$lhs %in% lv_names &
      partable$block == g)
    out[idx] <- out[idx] / eta[match(partable$lhs[idx], lv_names)]
  }

  # 5. defined/(in)equality constraints
  lav_standardize_constraints(out, partable, lavmodel)
}

# ---------------------------------------------------------------------------
# std.all / std.nox / user subset : standardize observed (and latent)
# variables.
#
# 'ov_std' and 'nox' control *which* observed variables are standardized:
#   - nox = FALSE, ov_std = NULL      -> all observed variables (std.all)
#   - nox = TRUE                      -> all but the exogenous 'x' (std.nox)
#   - nox = FALSE, ov_std = c(...)    -> only the observed variables listed in
#                                        'ov_std' (a generalization of std.nox)
# In all cases the latent variables are fully standardized (via
# lav_standardize_lv()).
# ---------------------------------------------------------------------------

lav_standardize_all <- function(lavobject = NULL,
                                partable = NULL, est = NULL, est_std = NULL,
                                glist = NULL, cov_std = TRUE, ov_var = NULL,
                                lv_var = NULL,
                                lavmodel = NULL, lavpartable = NULL,
                                cov_x = NULL, ov_std = NULL, nox = FALSE) {
  s <- lav_standardize_setup(
    lavobject = lavobject, lavmodel = lavmodel, lavpartable = lavpartable,
    partable = partable, est = est, glist = glist, cov_x = cov_x,
    need_cov_x = TRUE
  )
  lavmodel <- s$lavmodel
  lavpartable <- s$lavpartable
  partable <- s$partable
  est <- s$est
  glist <- s$glist
  cov_x <- s$cov_x

  if (is.null(est_std)) {
    est_std <- lav_standardize_lv(
      lavobject = lavobject,
      partable = partable, est = est, glist = glist,
      cov_std = cov_std, lv_var = lv_var, lavmodel = lavmodel,
      lavpartable = lavpartable
    )
  }

  out <- est_std
  n <- length(est_std)
  stopifnot(n == length(partable$lhs))

  vy <- lav_model_vy(
    lavmodel = lavmodel, glist = glist,
    diagonal_only = TRUE
  )

  for (g in 1:lavmodel@nblocks) {
    lv_names <- lav_pt_vnames(lavpartable, "lv", block = g)

    tmp <- lav_standardize_ov_sd(
      lavpartable = lavpartable, g = g, vy = vy, ov_var = ov_var,
      conditional_x = lavmodel@conditional.x, cov_x = cov_x
    )
    ov_names <- tmp$ov_names
    ov <- tmp$ov
    ov2 <- tmp$ov2

    # which observed variables should NOT be standardized in this block?
    # (std.all: none; std.nox: the exogenous 'x'; user: all but 'ov_std')
    if (nox) {
      ov_nostd <- lav_pt_vnames(lavpartable, "ov.x", block = g)
    } else if (is.null(ov_std)) {
      ov_nostd <- character(0L)
    } else {
      ov_nostd <- setdiff(ov_names, ov_std)
    }

    # 1a. "=~" regular indicators
    idx <- which(partable$op == "=~" & !(partable$rhs %in% lv_names) &
      !(partable$rhs %in% ov_nostd) &
      partable$block == g)
    out[idx] <- out[idx] / ov[match(partable$rhs[idx], ov_names)]

    # 1b. "=~" regular higher-order lv indicators -> handled in std.lv
    # 1c. "=~" indicators that are both in ov and lv -> not standardized here

    # 2. "~" regressions (and "<~")
    idx <- which((partable$op == "~" | partable$op == "<~") &
      partable$lhs %in% ov_names & !(partable$lhs %in% ov_nostd) &
      partable$block == g)
    out[idx] <- out[idx] / ov[match(partable$lhs[idx], ov_names)]

    idx <- which((partable$op == "~" | partable$op == "<~") &
      partable$rhs %in% ov_names & !(partable$rhs %in% ov_nostd) &
      partable$block == g)
    out[idx] <- out[idx] * ov[match(partable$rhs[idx], ov_names)]

    # 3a. "~~" ov
    # ATTENTION: in Mplus 4.1, the off-diagonal residual covariances
    #            were computed by the formula cov(i,j) / sqrt(i.var*j.var)
    #            where i.var and j.var were diagonal elements of OV
    #
    #            in Mplus 6.1 (but also AMOS and EQS), the i.var and j.var
    #            elements are the 'THETA' diagonal elements!!

    # variances
    rv_idx <- which(partable$op == "~~" & !(partable$lhs %in% lv_names) &
      !(partable$lhs %in% ov_nostd) &
      partable$lhs == partable$rhs &
      partable$block == g)
    out[rv_idx] <- (out[rv_idx] /
      ov2[match(partable$lhs[rv_idx], ov_names)])

    # covariances ov
    # three types:
    # - only lhs is OV (and fixed.x = FALSE)
    # - only rhs is OV (and fixed.x = FALSE)
    # - both lhs and rhs are OV (regular case)
    if (cov_std) {
      if (!is.complex(est[rv_idx])) {
        rv <- sqrt(abs(est[rv_idx]))
      } else {
        rv <- sqrt(est[rv_idx])
      }
      rv_names <- partable$lhs[rv_idx]
    }

    # left
    idx_lhs <- which(partable$op == "~~" &
      !(partable$lhs %in% lv_names) & !(partable$lhs %in% ov_nostd) &
      partable$lhs != partable$rhs &
      partable$block == g)
    if (length(idx_lhs) > 0L) {
      if (cov_std == FALSE) {
        out[idx_lhs] <-
          (out[idx_lhs] / ov[match(partable$lhs[idx_lhs], ov_names)])
      } else {
        out[idx_lhs] <-
          (out[idx_lhs] / rv[match(partable$lhs[idx_lhs], rv_names)])
      }
    }

    # right
    idx_rhs <- which(partable$op == "~~" &
      !(partable$rhs %in% lv_names) & !(partable$rhs %in% ov_nostd) &
      partable$lhs != partable$rhs &
      partable$block == g)
    if (length(idx_rhs) > 0L) {
      if (cov_std == FALSE) {
        out[idx_rhs] <-
          (out[idx_rhs] / ov[match(partable$rhs[idx_rhs], ov_names)])
      } else {
        out[idx_rhs] <-
          (out[idx_rhs] / rv[match(partable$rhs[idx_rhs], rv_names)])
      }
    }

    # 3b. "~~" lv -> handled in std.lv

    # 4a. "~1" ov
    idx <- which(partable$op == "~1" & !(partable$lhs %in% lv_names) &
      !(partable$lhs %in% ov_nostd) &
      partable$block == g)
    out[idx] <- out[idx] / ov[match(partable$lhs[idx], ov_names)]

    # 4b. "~1" lv -> handled in std.lv

    # 4c. "|" thresholds
    idx <- which(partable$op == "|" & !(partable$lhs %in% lv_names) &
      !(partable$lhs %in% ov_nostd) &
      partable$block == g)
    out[idx] <- out[idx] / ov[match(partable$lhs[idx], ov_names)]

    # 4d. "~*~" scales
    idx <- which(partable$op == "~*~" & !(partable$lhs %in% lv_names) &
      !(partable$lhs %in% ov_nostd) &
      partable$block == g)
    out[idx] <- 1.0
  }

  # 5. defined/(in)equality constraints
  lav_standardize_constraints(out, partable, lavmodel)
}

# std.nox: standardize all observed variables except the exogenous 'x'.
# Thin wrapper around lav_standardize_all().
lav_standardize_all_nox <- function(lavobject = NULL,
                                    partable = NULL, est = NULL, est_std = NULL,
                                    glist = NULL, cov_std = TRUE, ov_var = NULL,
                                    lv_var = NULL,
                                    lavmodel = NULL, lavpartable = NULL,
                                    cov_x = NULL) {
  lav_standardize_all(
    lavobject = lavobject, partable = partable, est = est, est_std = est_std,
    glist = glist, cov_std = cov_std, ov_var = ov_var, lv_var = lv_var,
    lavmodel = lavmodel, lavpartable = lavpartable, cov_x = cov_x,
    nox = TRUE
  )
}

# ---------------------------------------------------------------------------
# unstandardize observed variables (used by e.g. simulateData)
# ---------------------------------------------------------------------------

lav_unstandardize_ov <- function(partable, ov_var = NULL, cov_std = TRUE) {
  # check if ustart is missing; if so, look for est
  if (is.null(partable$ustart)) {
    partable$ustart <- partable$est
  }

  # check if block is missing
  if (is.null(partable$block)) {
    partable$block <- rep(1L, length(partable$ustart))
  }

  stopifnot(!any(is.na(partable$ustart)))
  est <- out <- partable$ustart

  # nblocks
  nblocks <- lav_pt_nblocks(partable)

  # if ov.var is NOT a list, make a list
  if (!is.list(ov_var)) {
    tmp <- ov_var
    ov_var <- vector("list", length = nblocks)
    ov_var[1:nblocks] <- list(tmp)
  }

  for (g in 1:nblocks) {
    ov_names <- lav_pt_vnames(partable, "ov", block = g) # not user
    lv_names <- lav_pt_vnames(partable, "lv", block = g)

    ov <- sqrt(ov_var[[g]])

    # 1a. "=~" regular indicators
    idx <- which(partable$op == "=~" & !(partable$rhs %in% lv_names) &
      partable$block == g)
    out[idx] <- out[idx] * ov[match(partable$rhs[idx], ov_names)]

    # 1b. "=~" regular higher-order lv indicators
    # 1c. "=~" indicators that are both in ov and lv

    # 2. "~" regressions (and "<~")
    idx <- which((partable$op == "~" | partable$op == "<~") &
      partable$lhs %in% ov_names &
      partable$block == g)
    out[idx] <- out[idx] * ov[match(partable$lhs[idx], ov_names)]

    idx <- which((partable$op == "~" | partable$op == "<~") &
      partable$rhs %in% ov_names &
      partable$block == g)
    out[idx] <- out[idx] / ov[match(partable$rhs[idx], ov_names)]

    # 3a. "~~" ov
    # ATTENTION: in Mplus 4.1, the off-diagonal residual covariances
    #            were computed by the formula cov(i,j) / sqrt(i.var*j.var)
    #            where i.var and j.var were diagonal elements of OV
    #
    #            in Mplus 6.1 (but also AMOS and EQS), the i.var and j.var
    #            elements are the 'THETA' diagonal elements!!

    # variances
    rv_idx <- which(partable$op == "~~" & !(partable$lhs %in% lv_names) &
      partable$lhs == partable$rhs &
      partable$block == g)
    out[rv_idx] <- (out[rv_idx] * ov[match(partable$lhs[rv_idx], ov_names)]
      * ov[match(partable$rhs[rv_idx], ov_names)])

    # covariances
    idx <- which(partable$op == "~~" & !(partable$lhs %in% lv_names) &
      partable$lhs != partable$rhs &
      partable$block == g)
    if (length(idx) > 0L) {
      if (cov_std == FALSE) {
        out[idx] <- (out[idx] * ov[match(partable$lhs[idx], ov_names)]
          * ov[match(partable$rhs[idx], ov_names)])
      } else {
        if (!is.complex(out[rv_idx])) {
          rv <- sqrt(abs(out[rv_idx]))
        } else {
          rv <- sqrt(out[rv_idx])
        }
        rv_names <- partable$lhs[rv_idx]
        out[idx] <- (out[idx] * rv[match(partable$lhs[idx], rv_names)]
          * rv[match(partable$rhs[idx], rv_names)])
      }
    }

    # 3b. "~~" lv

    # 4a. "~1" ov
    idx <- which(partable$op == "~1" & !(partable$lhs %in% lv_names) &
      partable$block == g)
    out[idx] <- out[idx] * ov[match(partable$lhs[idx], ov_names)]

    # 4b. "~1" lv
  }

  # 5a ":="
  # 5b "=="
  # 5c. "<" or ">"

  out
}
