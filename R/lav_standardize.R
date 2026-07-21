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

# scaling corrections for two-way latent product terms (interactions
# X:Z, quadratic terms X:X) in block 'g' of lav_standardize_lv(), based on
# Brandt et al. (2015, Eq. 13) and Kelava & Brandt (2022). The standardized
# solution treats a product term as the product of the centered,
# standardized components, whereas the model parameters describe the raw
# product A:B. Three corrections are needed, all exact (no distributional
# assumptions):

# (1) scaling SD (FV): the SD of a standardized product term is
# SD(A)*SD(B), not SD(A:B). eta[A:B] is replaced accordingly, so that the
# generic scaling in the caller applies the right factor everywhere A:B
# appears.

# (2) simple main effects (KS, 08/05/26): with uncentered components, b1 in
# y ~ b1*x + b2*z + b3*x:z + b4*x:x is the slope at x = z = 0. The
# standardized solution reports the slope at the means, i.e. the partial
# derivative dy/dx = b1 + b3*z + 2*b4*x evaluated at E(x), E(z). Looping
# over the product terms, each term adds beta_ab * E(other component) to
# the main effect of each of its components.

# (3) product moments (FV, 05/07/26): the psi/alpha rows of A:B describe the
# uncentered product. Without centering them first, the standardized
# moments would depend on how the mean structure is parameterized. With
# ac = A - E(A) and bc = B - E(B) we have the exact identity
# ac*bc = A:B - E(B)*A - E(A)*B + E(A)*E(B), so the centered product is
# a linear combination of elements of eta. Collect the weights in a
# matrix K (one row per variable, regular lvs keep their identity row),
# e.g. for eta = (X, Z, X:Z, X:X) with mx = E(X), mz = E(Z):
# row X:Z = (-mz, -mx, 1, 0), row X:X = (-2*mx, 0, 0, 1). Then
# K %*% PSI %*% t(K) holds all centered (co)variances at once (constants
# do not affect them), and the centered means follow from
# E(ac*bc) = E(A:B) - E(A)*E(B). K acts on psi_g (not VETA): the "~~"
# rows are psi elements, and for a residual covariance with an
# endogenous lv, PSI[X, Y] = Cov(X, zeta_Y) is exactly what the
# centering needs. The corrected rows are also written into est_c,
# because the cov_std block in the caller derives its sqrt(variance)
# divisors from est_c (not out, which is scaled by then).

# Returns the corrected (out, est_c, eta)
lav_standardize_lv_int <- function(out, est_c, eta, est, partable, lv_names,
                                   lv_int_names, eeta, psi_g, g) {
  # component indices of each usable product term (int = A:B, a = A, b = B)
  int_names_ok <- character(0L)
  int_comp <- list()
  for (int_name in lv_int_names) {
    components <- strsplit(int_name, ":", fixed = TRUE)[[1L]]
    # only two-way terms: a three-way (or higher) product would need its own
    # decomposition, and the product of three SDs
    if (length(components) != 2L) next
    idx_int <- match(int_name, lv_names)
    idx_a   <- match(components[1L], lv_names)
    idx_b   <- match(components[2L], lv_names)
    # need a, b and the interaction term before moving on
    if (is.na(idx_a) || is.na(idx_b) || is.na(idx_int)) next
    int_names_ok <- c(int_names_ok, int_name)
    int_comp <- c(int_comp,
      list(list(int = idx_int, a = idx_a, b = idx_b)))
  }
  if (length(int_comp) == 0L) {
    return(list(out = out, est_c = est_c, eta = eta))
  }

  # only used by the shift in (2), which needs nonzero means, but guarding
  # it with any(eeta != 0) would cost as much as building it
  reg_row <- partable$op == "~" & partable$block == g
  for (i in seq_along(int_comp)) {
    p <- int_comp[[i]]

    # (1) replace ETA for the product term with SD(A)*SD(B)
    eta[p$int] <- eta[p$a] * eta[p$b]

    # (2) shift the simple main effects to the mean-evaluated slopes
    exp_a <- eeta[p$a]
    exp_b <- eeta[p$b]
    if (exp_a != 0 || exp_b != 0) {
      int_row <- reg_row & partable$rhs == int_names_ok[i]
      for (dep in unique(partable$lhs[int_row])) {
        # beta_ab is read from est (not partable$est): in the *_x wrappers
        # (delta method, bootstrap) est is the perturbed parameter vector, and
        # the gradient of the shifted main effects wrt beta_ab must not be lost
        beta_ab <- est[int_row & partable$lhs == dep]
        dep_row <- reg_row & partable$lhs == dep
        idx_beta_a <- which(dep_row & partable$rhs == lv_names[p$a])
        idx_beta_b <- which(dep_row & partable$rhs == lv_names[p$b])
        out[idx_beta_a] <- out[idx_beta_a] + beta_ab * exp_b
        out[idx_beta_b] <- out[idx_beta_b] + beta_ab * exp_a
      }
    }
  }

  # (3) center the (co)variances and means of the product terms
  if (!is.null(psi_g) && any(eeta != 0)) {
    kmat <- diag(nrow(psi_g))
    for (p in int_comp) {
      kmat[p$int, p$a] <- kmat[p$int, p$a] - eeta[p$b]
      kmat[p$int, p$b] <- kmat[p$int, p$b] - eeta[p$a]
    }
    # shift_mat holds, for every pair of latent variables, the mean-driven
    # part of their (co)variance (zero unless a product term is involved)
    shift_mat <- psi_g - kmat %*% psi_g %*% t(kmat)

    # "~~" rows involving at least one product term: this covers the
    # variance of a product term (X:Z ~~ X:Z), its covariance with a
    # regular lv (X:Z ~~ V), and product-product covariances (X:Z ~~ X:X).
    # Rows whose other side is not a lv are left alone (keep)
    idx_cov <- which(partable$op == "~~" & partable$block == g &
      (partable$lhs %in% int_names_ok |
       partable$rhs %in% int_names_ok))
    l_idx <- match(partable$lhs[idx_cov], lv_names)
    r_idx <- match(partable$rhs[idx_cov], lv_names)
    keep <- which(!is.na(l_idx) & !is.na(r_idx))
    shift <- shift_mat[cbind(l_idx[keep], r_idx[keep])]
    out[idx_cov[keep]] <- out[idx_cov[keep]] - shift
    # of est_c only the variance entries are read back later, but centering
    # all touched rows keeps est_c fully consistent and a diagonal-only
    # filter would cost as much as it saves
    est_c[idx_cov[keep]] <- est_c[idx_cov[keep]] - shift

    # "~1" rows of the product terms: E(ac*bc) = E(A:B) - E(A)*E(B). After
    # the generic scaling this equals, e.g., cor(X, Z) for X:Z
    idx_mean <- which(partable$op == "~1" & partable$block == g &
      partable$lhs %in% int_names_ok)
    for (row_i in idx_mean) {
      p <- int_comp[[match(partable$lhs[row_i], int_names_ok)]]
      out[row_i] <- out[row_i] - eeta[p$a] * eeta[p$b]
    }
  }

  list(out = out, est_c = est_c, eta = eta)
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

  # any latent product terms (in any block)? only their corrections (see
  # lav_standardize_lv_int) need E(ETA), psi and the mm indices. Regular
  # models skip these, and E(ETA) is LISREL-only anyway (RAM would stop)
  lv_int_any <- length(lav_pt_vnames(lavpartable, "lv.interaction")) > 0L
  if (lv_int_any) {
    mm_idx <- lav_model_group_mm_indices(lavmodel@nmat)
  }

  # compute ETA
  if (is.null(lv_var)) {
    lv_eta <- lav_model_veta(
      lavmodel = lavmodel,
      glist = glist
    )

    if (lv_int_any) {
      lv_eeta <- lav_model_eeta(
        lavmodel = lavmodel,
        glist = glist,
        lavsamplestats = lavobject@SampleStats
      )
    }
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
      veta_g <- lv_eta[[g]] # full model-implied V(ETA), only kept as a
                            # psi fallback for the centering below (RAM)
      eta2 <- diag(veta_g)
    } else {
      veta_g <- NULL
      eta2 <- lv_var[[g]]
    }
    # change negative values to NA
    eta2[eta2 < 0] <- as.numeric(NA)
    eta <- sqrt(eta2)

    # two-way latent product terms (X:Z, X:X): shift the simple main
    # effects, scale by SD(A)*SD(B), and center the product moments (see
    # lav_standardize_lv_int() for the rationale)
    if (lv_int_any) {
      lv_int_names <- lav_pt_vnames(lavpartable, "lv.interaction",
        block = g
      )
    } else {
      lv_int_names <- character(0L)
    }
    if (length(lv_int_names) > 0L) {
      if (is.null(lv_var)) {
        eeta <- lv_eeta[[g]]
      } else {
        eeta <- numeric(length(eta2))
      }
      psi_g <- glist[mm_idx[[g]]]$psi
      if (is.null(psi_g)) { # e.g., RAM representation
        psi_g <- veta_g
      }
      tmp <- lav_standardize_lv_int(
        out = out, est_c = est_c, eta = eta, est = est, partable = partable,
        lv_names = lv_names, lv_int_names = lv_int_names, eeta = eeta,
        psi_g = psi_g, g = g
      )
      out <- tmp$out
      est_c <- tmp$est_c
      eta <- tmp$eta
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
      if (lavmodel@fixed.x) {
        ov_nostd <- lav_pt_vnames(lavpartable, "ov.x", block = g)
      } else {
        # if fixed.x = FALSE, the exogenous 'x' are treated as random
        # variables like any other, and std.nox reduces to std.all
        ov_nostd <- character(0L)
      }
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

    # D-augmented ML correlation mode (FREE ~*~ scales): the intercepts
    # live in the ORIGINAL metric (the implied mean is not rescaled by
    # the ~*~ scales), while VY is in the unit-diagonal P metric -- so
    # the "~1" rows are divided by the ~*~ scales as well. (Using the
    # est values keeps the delta-method standard errors consistent.)
    if (!lavmodel@categorical && length(idx) > 0L &&
        lav_model_delta_free(lavmodel)) {
      scl_idx <- which(partable$op == "~*~" & partable$block == g &
        partable$lhs == partable$rhs)
      if (length(scl_idx) > 0L) {
        m2 <- match(partable$lhs[idx], partable$lhs[scl_idx])
        ok2 <- which(!is.na(m2))
        if (length(ok2) > 0L) {
          out[idx[ok2]] <- out[idx[ok2]] / est[scl_idx][m2[ok2]]
        }
      }
    }

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
