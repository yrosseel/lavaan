# collect information about the model that we can use
# (eg. is theta diagonal or not, is the structural model recursive or not,
#  is the model just a regression model, etc)
#
# initial version: YR 15 March 2021
# - YR 05 Oct 2021: use det(I - B) to check if B is acyclic
# - YR 11 Nov 2021: if no latents, and conditional.x = TRUE, we may have no
#                   beta matrix


# note: there is no 'lavmodel' yet, because we call this in lav_model.R
lav_model_properties <- function(glist, lavpartable = NULL,
                                 nmat = NULL, m_free_idx = NULL) {
  lavpta <- lav_partable_attributes(lavpartable)
  nblocks <- lavpta$nblocks

  # is the model a univariate/multivariate linear multiple regression
  # model (per block)?
  uvreg <- logical(nblocks)
  uvord <- logical(nblocks)
  mvreg <- logical(nblocks)
  acyclic <- rep(as.logical(NA), nblocks)
  bowfree <- rep(as.logical(NA), nblocks)
  nexo <- integer(nblocks)

  for (g in seq_len(nblocks)) {
    # at least 1 regression
    if (length(lavpta$vnames$eqs.y[[g]]) == 0L) {
      next
    }

    # find beta index for this block
    mm_in_block <- 1:nmat[g] + cumsum(c(0L, nmat))[g]
    mlist <- glist[mm_in_block]
    beta_idx <- which(names(mlist) == "beta") + cumsum(c(0L, nmat))[g]
    psi_idx <- which(names(mlist) == "psi") + cumsum(c(0L, nmat))[g]

    if (length(beta_idx) > 0L) {
      # 1. acyclic?
      m_b <- glist[[beta_idx]]
      # keep fixed values (if any); fill in 1 in all 'free' positions
      m_b[m_free_idx[[beta_idx]]] <- 1
      acyclic[g] <- lav_graph_is_acyclic(m_b)

      # 2. bow-free?
      b_one <- as.integer(m_b != 0)
      psi <- glist[[psi_idx]]
      # keep fixed values (if any); fill in 1 in all 'free' positions
      psi[m_free_idx[[psi_idx]]] <- 1
      psi_one <- as.integer(psi != 0)
      both_one <- b_one + psi_one
      if (any(both_one > 1)) {
        bowfree[g] <- FALSE
      } else {
        bowfree[g] <- TRUE
      }
    } else {
      # perhaps conditional.x = TRUE?
      # if there is no BETA, then we only have Gamma, and the
      # system must be acyclic
      acyclic[g] <- TRUE
      # and also bowfree
      bowfree[g] <- TRUE
    }


    # no latent variables, at least 1 dependent variable
    if (lavpta$nfac[[g]] > 0L) {
      next
    }

    # no mediators
    if (length(lavpta$vnames$eqs.y[[g]]) !=
      length(lavpta$vnames$ov.y[[g]])) {
      next
    }

    # categorical y?
    if (length(lavpta$vnames$ov.ord[[g]]) > 0L) {
      # we only flag the univariate version
      if (length(lavpta$vnames$ov.ord[[g]]) == 1L &&
        length(lavpta$vnames$ov.y[[g]]) == 1L &&
        lavpta$vnames$ov.ord[[g]][1] == lavpta$vnames$ov.y[[g]][1]) {
        uvord[g] <- TRUE
      }

      # mvreg?
    } else {
      if (length(lavpta$vnames$ov.y[[g]]) > 1L) {
        mvreg[g] <- TRUE
      } else {
        uvreg[g] <- TRUE
      }
    }

    nexo[g] <- length(lavpta$vnames$eqs.x[[g]])
  } # g

  modprop <- list(
    uvreg = uvreg, uvord = uvord, mvreg = mvreg,
    nexo = nexo, acyclic = acyclic, bowfree = bowfree
  )

  modprop
}
