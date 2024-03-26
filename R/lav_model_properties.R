# collect information about the model that we can use
# (eg. is theta diagonal or not, is the structurual model recursive or not,
#  is the model just a regression model, etc)
#
# initial version: YR 15 March 2021
# - YR 05 Oct 2021: use det(I - B) to check if B is acyclic
# - YR 11 Nov 2021: if no latents, and conditional.x = TRUE, we may have no
#                   beta matrix


# note: there is no 'lavmodel' yet, because we call this in lav_model.R
lav_model_properties <- function(GLIST, lavpartable = NULL,
                                 nmat = NULL, m.free.idx = NULL) {
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
    mm.in.block <- 1:nmat[g] + cumsum(c(0L, nmat))[g]
    MLIST <- GLIST[mm.in.block]
    beta.idx <- which(names(MLIST) == "beta") + cumsum(c(0L, nmat))[g]
    psi.idx <- which(names(MLIST) == "psi") + cumsum(c(0L, nmat))[g]

    if (length(beta.idx) > 0L) {
      # 1. acyclic?
      B <- GLIST[[beta.idx]]
      # keep fixed values (if any); fill in 1 in all 'free' positions
      B[m.free.idx[[beta.idx]]] <- 1
      IminB <- diag(nrow(B)) - B
      # if B is acyclic, we should be able to permute the rows/cols of B
      # so that B is upper/lower triangular, and so det(I-B) = 1
      if (det(IminB) == 1) {
        acyclic[g] <- TRUE
      } else {
        acyclic[g] <- FALSE
      }

      # 2. bow-free?
      B.one <- as.integer(B != 0)
      Psi <- GLIST[[psi.idx]]
      # keep fixed values (if any); fill in 1 in all 'free' positions
      Psi[m.free.idx[[psi.idx]]] <- 1
      Psi.one <- as.integer(Psi != 0)
      Both.one <- B.one + Psi.one
      if (any(Both.one > 1)) {
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
