# utility functions for pairwise maximum likelihood

# stub for fml_deriv1
fml_deriv1 <- function(Sigma.hat = NULL, # model-based var/cov/cor
                       TH = NULL, # model-based thresholds + means
                       th.idx = NULL, # threshold idx per variable
                       num.idx = NULL, # which variables are numeric
                       X = NULL, # data
                       eXo = NULL, # external covariates
                       lavcache = NULL, # housekeeping stuff
                       scores = FALSE, # return case-wise scores
                       negative = TRUE) {
  lav_msg_stop(gettext("not implemented"))
}

# the first derivative of the pairwise logLik function with respect to the
# thresholds/slopes/var/correlations; together with DELTA, we can use the
# chain rule to get the gradient
# this is adapted from code written by Myrsini Katsikatsou
# first attempt - YR 5 okt 2012
# HJ 18/10/23: Modification for complex design and completely observed data (no
# missing) with only ordinal indicators to get the right gradient for the
# optimisation and Hessian computation.
pml_deriv1 <- function(Sigma.hat = NULL, # model-based var/cov/cor
                       Mu.hat = NULL, # model-based means
                       TH = NULL, # model-based thresholds + means
                       th.idx = NULL, # threshold idx per variable
                       num.idx = NULL, # which variables are numeric
                       X = NULL, # data
                       eXo = NULL, # external covariates
                       wt = NULL, # case weights (not used yet)
                       lavcache = NULL, # housekeeping stuff
                       PI = NULL, # slopes
                       missing = "listwise", # how to deal with missings
                       scores = FALSE, # return case-wise scores
                       negative = TRUE) { # multiply by -1

  # diagonal of Sigma.hat is not necessarily 1, even for categorical vars
  Sigma.hat2 <- Sigma.hat
  if (length(num.idx) > 0L) {
    diag(Sigma.hat2)[-num.idx] <- 1
  } else {
    diag(Sigma.hat2) <- 1
  }
  Cor.hat <- cov2cor(Sigma.hat2) # to get correlations (rho!)
  cors <- lav_matrix_vech(Cor.hat, diagonal = FALSE)

  if (any(abs(cors) > 1)) {
    # what should we do now... force cov2cor?
    # cat("FFFFOOOORRRRRCEEE PD!\n")
    # Sigma.hat <- Matrix::nearPD(Sigma.hat)
    # Sigma.hat <- as.matrix(Sigma.hat$mat)
    # Sigma.hat <- cov2cor(Sigma.hat)
    # cors <- Sigma.hat[lower.tri(Sigma.hat)]
    idx <- which(abs(cors) > 0.99)
    cors[idx] <- 0.99 # clip
    # cat("CLIPPING!\n")
  }

  nvar <- nrow(Sigma.hat)
  pstar <- nvar * (nvar - 1) / 2
  ov.types <- rep("ordered", nvar)
  if (length(num.idx) > 0L) ov.types[num.idx] <- "numeric"
  if (!is.null(eXo)) {
    nexo <- ncol(eXo)
  } else {
    nexo <- 0
  }


  if (all(ov.types == "numeric")) {
    N.TH <- nvar
  } else {
    N.TH <- length(th.idx)
  }
  N.SL <- nvar * nexo
  N.VAR <- length(num.idx)
  N.COR <- pstar

  # add num.idx to th.idx
  if (length(num.idx) > 0L) {
    th.idx[th.idx == 0] <- num.idx
  }

  # print(Sigma.hat); print(TH); print(th.idx); print(num.idx); print(str(X))

  # shortcut for ordinal-only/no-exo case
  if (!scores && all(ov.types == "ordered") && nexo == 0L) {
    # >>>>>>>> HJ/MK PML CODE >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    if (is.null(wt)) {
      n.xixj.vec <- lavcache$bifreq
    } else {
      n.xixj.vec <- lavcache$sum_obs_weights_xixj_ab_vec
    }
    gradient <- grad_tau_rho(
      no.x = nvar,
      all.thres = TH,
      index.var.of.thres = th.idx,
      rho.xixj = cors,
      n.xixj.vec = n.xixj.vec,
      out.LongVecInd = lavcache$long
    )

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    if (missing == "available.cases") {
      uniPI <- univariateExpProbVec(TH = TH, th.idx = th.idx)
      tmp <- lavcache$uniweights / uniPI

      var.idx <- split(th.idx, th.idx)
      var.idx <- unlist(lapply(var.idx, function(x) {
        c(x, x[1])
      }))

      tmp.varwise <- split(tmp, var.idx)
      tmp1 <- unlist(lapply(
        tmp.varwise,
        function(x) {
          c(x[-length(x)])
        }
      ))
      tmp2 <- unlist(lapply(tmp.varwise, function(x) {
        c(x[-1])
      }))

      uni.der.tau <- dnorm(TH) * (tmp1 - tmp2)
      nTH <- length(TH)
      gradient[1:nTH] <- gradient[1:nTH] + uni.der.tau
    }

    if (negative) {
      gradient <- -1 * gradient
    }
    return(gradient)
  }

  # in this order: TH/MEANS + SLOPES + VAR + COR
  GRAD.size <- N.TH + N.SL + N.VAR + N.COR

  # scores or gradient?
  if (scores) {
    SCORES <- matrix(0, nrow(X), GRAD.size) # we will sum up over all pairs
  } else {
    GRAD <- matrix(0, pstar, GRAD.size) # each pair is a row
  }
  PSTAR <- matrix(0, nvar, nvar) # utility matrix, to get indices
  PSTAR[lav_matrix_vech_idx(nvar, diagonal = FALSE)] <- 1:pstar
  N <- length(X[, 1])

  for (j in seq_len(nvar - 1L)) {
    for (i in (j + 1L):nvar) {
      # cat(" i = ", i, " j = ", j, "\n") # debug only
      pstar.idx <- PSTAR[i, j]
      cor.idx <- N.TH + N.SL + N.VAR + PSTAR[i, j]
      th.idx_i <- which(th.idx == i)
      th.idx_j <- which(th.idx == j)
      if (nexo > 0L) {
        sl.idx_i <- N.TH + seq(i, by = nvar, length.out = nexo)
        sl.idx_j <- N.TH + seq(j, by = nvar, length.out = nexo)

        if (length(num.idx) > 0L) {
          var.idx_i <- N.TH + N.SL + match(i, num.idx)
          var.idx_j <- N.TH + N.SL + match(j, num.idx)
        }
      } else {
        if (length(num.idx) > 0L) {
          var.idx_i <- N.TH + match(i, num.idx)
          var.idx_j <- N.TH + match(j, num.idx)
        }
      }
      if (ov.types[i] == "numeric" && ov.types[j] == "numeric") {
        if (nexo > 1L) {
          lav_msg_stop(gettext(
            "mixed + exo in PML not implemented;
            try optim.gradient = \"numerical\""))
        }

        SC <- lav_mvnorm_scores_mu_vech_sigma(
          Y = X[, c(i, j)],
          Mu = Mu.hat[c(i, j)], Sigma = Sigma.hat[c(i, j), c(i, j)]
        )

        if (scores) {
          if (all(ov.types == "numeric") && nexo == 0L) {
            # MU1 + MU2
            SCORES[, c(i, j)] <- SCORES[, c(i, j)] + SC[, c(1, 2)]
            # VAR1 + COV_12 + VAR2
            var.idx <- (nvar +
              lav_matrix_vech_match_idx(nvar, idx = c(i, j)))
            SCORES[, var.idx] <- SCORES[, var.idx] + SC[, c(3, 4, 5)]
          } else { # mixed ordered/continuous
            # MU
            mu.idx <- c(th.idx_i, th.idx_j)
            SCORES[, mu.idx] <- SCORES[, mu.idx] + (-1) * SC[, c(1, 2)]
            # VAR+COV
            var.idx <- c(var.idx_i, cor.idx, var.idx_j)
            SCORES[, var.idx] <- SCORES[, var.idx] + SC[, c(3, 4, 5)]
          }
        } else {
          if (all(ov.types == "numeric") && nexo == 0L) {
            mu.idx <- c(i, j)
            sigma.idx <- (nvar +
              lav_matrix_vech_match_idx(nvar, idx = c(i, j)))
            # MU1 + MU2
            GRAD[pstar.idx, mu.idx] <-
              colSums(SC[, c(1, 2)], na.rm = TRUE)
          } else {
            mu.idx <- c(th.idx_i, th.idx_j)
            sigma.idx <- c(var.idx_i, cor.idx, var.idx_j)
            # MU (reverse sign!)
            GRAD[pstar.idx, mu.idx] <-
              -1 * colSums(SC[, c(1, 2)], na.rm = TRUE)
          }
          # SIGMA
          GRAD[pstar.idx, sigma.idx] <-
            colSums(SC[, c(3, 4, 5)], na.rm = TRUE)
        } # gradient only
      } else if (ov.types[i] == "numeric" && ov.types[j] == "ordered") {
        # polyserial correlation
        if (nexo > 1L) {
          lav_msg_stop(gettext(
            "mixed + exo in PML not implemented;
            try optim.gradient = \"numerical\""))
        }

        SC.COR.UNI <- lav_bvmix_cor_scores(
          Y1 = X[, i], Y2 = X[, j],
          eXo = NULL, wt = wt,
          evar.y1 = Sigma.hat[i, i],
          beta.y1 = Mu.hat[i],
          th.y2 = TH[th.idx == j],
          sl.y2 = NULL,
          rho = Cor.hat[i, j],
          sigma.correction = TRUE
        )

        if (scores) {
          # MU
          SCORES[, th.idx_i] <- (SCORES[, th.idx_i] +
            -1 * SC.COR.UNI$dx.mu.y1)
          # TH
          SCORES[, th.idx_j] <- (SCORES[, th.idx_j] +
            SC.COR.UNI$dx.th.y2)
          # VAR
          SCORES[, var.idx_i] <- (SCORES[, var.idx_i] +
            SC.COR.UNI$dx.var.y1)
          # COR
          SCORES[, cor.idx] <- (SCORES[, cor.idx] +
            SC.COR.UNI$dx.rho)
        } else {
          # MU
          GRAD[pstar.idx, th.idx_i] <-
            -1 * sum(SC.COR.UNI$dx.mu.y1, na.rm = TRUE)
          # TH
          GRAD[pstar.idx, th.idx_j] <-
            colSums(SC.COR.UNI$dx.th.y2, na.rm = TRUE)
          # VAR
          GRAD[pstar.idx, var.idx_i] <-
            sum(SC.COR.UNI$dx.var.y1, na.rm = TRUE)
          # COR
          GRAD[pstar.idx, cor.idx] <-
            sum(SC.COR.UNI$dx.rho, na.rm = TRUE)
        } # grad only
      } else if (ov.types[j] == "numeric" && ov.types[i] == "ordered") {
        # polyserial correlation
        if (nexo > 1L) {
          lav_msg_stop(gettext(
            "mixed + exo in PML not implemented;
            try optim.gradient = \"numerical\""))
        }

        SC.COR.UNI <- lav_bvmix_cor_scores(
          Y1 = X[, j], Y2 = X[, i],
          eXo = NULL, wt = wt,
          evar.y1 = Sigma.hat[j, j],
          beta.y1 = Mu.hat[j],
          th.y2 = TH[th.idx == i],
          rho = Cor.hat[i, j],
          sigma.correction = TRUE
        )

        if (scores) {
          # MU
          SCORES[, th.idx_j] <- (SCORES[, th.idx_j] +
            -1 * SC.COR.UNI$dx.mu.y1)
          # TH
          SCORES[, th.idx_i] <- (SCORES[, th.idx_i] +
            SC.COR.UNI$dx.th.y2)
          # VAR
          SCORES[, var.idx_j] <- (SCORES[, var.idx_j] +
            SC.COR.UNI$dx.var.y1)
          # COR
          SCORES[, cor.idx] <- (SCORES[, cor.idx] +
            SC.COR.UNI$dx.rho)
        } else {
          # MU
          GRAD[pstar.idx, th.idx_j] <-
            -1 * sum(SC.COR.UNI$dx.mu.y1, na.rm = TRUE)
          # TH
          GRAD[pstar.idx, th.idx_i] <-
            colSums(SC.COR.UNI$dx.th.y2, na.rm = TRUE)
          # VAR
          GRAD[pstar.idx, var.idx_j] <-
            sum(SC.COR.UNI$dx.var.y1, na.rm = TRUE)
          # COR
          GRAD[pstar.idx, cor.idx] <-
            sum(SC.COR.UNI$dx.rho, na.rm = TRUE)
        } # grad only
      } else if (ov.types[i] == "ordered" && ov.types[j] == "ordered") {
        # polychoric correlation
        if (nexo == 0L) {
          SC.COR.UNI <-
            lav_bvord_cor_scores(
              Y1 = X[, i], Y2 = X[, j],
              eXo = NULL, wt = wt,
              rho = Sigma.hat[i, j],
              fit.y1 = NULL, # fixme
              fit.y2 = NULL, # fixme
              th.y1 = TH[th.idx == i],
              th.y2 = TH[th.idx == j],
              sl.y1 = NULL,
              sl.y2 = NULL,
              na.zero = TRUE
            )
        } else {
          SC.COR.UNI <-
            pc_cor_scores_PL_with_cov(
              Y1 = X[, i],
              Y2 = X[, j],
              eXo = eXo,
              Rho = Sigma.hat[i, j],
              th.y1 = TH[th.idx == i],
              th.y2 = TH[th.idx == j],
              sl.y1 = PI[i, ],
              sl.y2 = PI[j, ],
              missing.ind = missing
            )
        }

        if (scores) {
          # TH
          SCORES[, th.idx_i] <- SCORES[, th.idx_i] + SC.COR.UNI$dx.th.y1
          SCORES[, th.idx_j] <- SCORES[, th.idx_j] + SC.COR.UNI$dx.th.y2

          # SL
          if (nexo > 0L) {
            SCORES[, sl.idx_i] <- SCORES[, sl.idx_i] + SC.COR.UNI$dx.sl.y1
            SCORES[, sl.idx_j] <- SCORES[, sl.idx_j] + SC.COR.UNI$dx.sl.y2
          }
          # NO VAR
          # RHO
          SCORES[, cor.idx] <- SCORES[, cor.idx] + SC.COR.UNI$dx.rho
        } else {
          # TH
          if (length(th.idx_i) > 1L) {
            GRAD[pstar.idx, th.idx_i] <-
              colSums(SC.COR.UNI$dx.th.y1, na.rm = TRUE)
          } else {
            GRAD[pstar.idx, th.idx_i] <-
              sum(SC.COR.UNI$dx.th.y1, na.rm = TRUE)
          }
          if (length(th.idx_j) > 1L) {
            GRAD[pstar.idx, th.idx_j] <-
              colSums(SC.COR.UNI$dx.th.y2, na.rm = TRUE)
          } else {
            GRAD[pstar.idx, th.idx_j] <-
              sum(SC.COR.UNI$dx.th.y2, na.rm = TRUE)
          }

          # SL
          if (nexo > 0L) {
            if (length(sl.idx_i) > 1L) {
              GRAD[pstar.idx, sl.idx_i] <-
                colSums(SC.COR.UNI$dx.sl.y1, na.rm = TRUE)
            } else {
              GRAD[pstar.idx, sl.idx_i] <-
                sum(SC.COR.UNI$dx.sl.y1, na.rm = TRUE)
            }
            if (length(sl.idx_j) > 1L) {
              GRAD[pstar.idx, sl.idx_j] <-
                colSums(SC.COR.UNI$dx.sl.y2, na.rm = TRUE)
            } else {
              GRAD[pstar.idx, sl.idx_j] <-
                sum(SC.COR.UNI$dx.sl.y2, na.rm = TRUE)
            }
          }
          # NO VAR

          # RHO
          GRAD[pstar.idx, cor.idx] <-
            sum(SC.COR.UNI$dx.rho, na.rm = TRUE)
        }

        # GRAD2 <- numDeriv::grad(func = pc_logl_x,
        #                        x = c(Sigma.hat[i,j],
        #                              TH[ th.idx == i ],
        #                              TH[ th.idx == j]),
        #                        Y1  = X[,i],
        #                        Y2  = X[,j],
        #                        eXo = eXo,
        #                        nth.y1 = sum( th.idx == i ),
        #                        nth.y2 = sum( th.idx == j ))
      }
    }
  }

  if (missing == "available.cases" && all(ov.types == "ordered")) {
    if (nexo == 0L) {
      UNI_SCORES <- matrix(0, nrow(X), N.TH)
      for (i in seq_len(nvar)) {
        th.idx_i <- which(th.idx == i)
        derY1 <- uni_scores(
          Y1 = X[, i], th.y1 = TH[th.idx == i],
          eXo = NULL, sl.y1 = NULL,
          weights.casewise = lavcache$uniweights.casewise
        )
        UNI_SCORES[, th.idx_i] <- derY1$dx.th.y1
      }
    } else {
      UNI_SCORES <- matrix(0, nrow(X), ncol = (N.TH + N.SL))
      for (i in seq_len(nvar)) {
        th.idx_i <- which(th.idx == i)
        sl.idx_i <- N.TH + seq(i, by = nvar, length.out = nexo)
        derY1 <- uni_scores(
          Y1 = X[, i], th.y1 = TH[th.idx == i],
          eXo = eXo, sl.y1 = PI[i, ],
          weights.casewise = lavcache$uniweights.casewise
        )
        UNI_SCORES[, th.idx_i] <- derY1$dx.th.y1
        UNI_SCORES[, sl.idx_i] <- derY1$dx.sl.y1
      }
      if (scores) {
        SCORES <- SCORES[, 1:(N.TH + N.SL)] + UNI_SCORES
      } else {
        uni_gradient <- colSums(UNI_SCORES)
      }
    }
  }

  # do we need scores?
  if (scores) {
    return(SCORES)
  }

  # DEBUG
  # :print(GRAD)
  ###########



  # gradient is sum over all pairs
  gradient <- colSums(GRAD, na.rm = TRUE)

  if (missing == "available.cases" && all(ov.types == "ordered")) {
    if (nexo == 0L) {
      gradient[1:N.TH] <- gradient + uni_gradient
    } else {
      gradient[1:(N.TH + N.SL)] <- gradient + uni_gradient
    }
  }

  # we multiply by -1 because we minimize
  if (negative) {
    gradient <- -1 * gradient
  }

  gradient
}


### all code below written by Myrsini Katsikatsou


# The function grad_tau_rho
# input:
# no.x -  is scalar, the number of ordinal variables
# all.thres  - is vector containing the thresholds of all variables in the
#              following order: thres_var1, thres_var2,..., thres_var_p
#              within each variable the thresholds are in ascending order
#              Note that all.thres do NOT contain tau_0=-Inf and tau_last=Inf
#              for all variables.
# index.var.of.thres - a vector keeping track to which variable the thresholds
#                      in all.thres belongs to, it is of the form
#                      (1,1,1..., 2,2,2,...,  p,p,p,...)
# rho.xixj - is the vector of all correlations where j runs faster than i
#            i.e. the order is rho_12, rho_13, ..., rho_1p, rho_23, ..., rho_2p,
#            etc.
# n.xixj.vec -  a vector with the observed frequency for every combination
#               of categories and every pair. The frequencies are given in
#               the same order as the expected probabilities in the output of
#               pairwiseExpProbVec output
# out.LongVecInd - it is the output of function LongVecInd
# the output: it gives the elements of der.L.to.tau and der.L.to.rho in this
#             order. The elements of der.L.to.tau where the elements are
#             ordered as follows: the thresholds of each variable with respect
#             to ascending order of the variable index (i.e. thres_var1,
#             thres_var2, etc.) and within each variable the thresholds in
#             ascending order.
#             The elements of vector der.L.to.rho are der.Lxixj.to.rho.xixj
#             where j runs faster than i.

# The function depends on four other functions: LongVecTH.Rho,
# pairwiseExpProbVec, derLtoRho, and derLtoTau, all given below.

# if n.xixj.ab is either an array or a list the following should be done
# n.xixj.vec <- if(is.array(n.xixj.ab)) {
#                 c(n.xixj.ab)
#              } else if(is.list(n.xixj.ab)){
#                 unlist(n.xixj.ab)
#              }


grad_tau_rho <- function(no.x, all.thres, index.var.of.thres, rho.xixj,
                         n.xixj.vec, out.LongVecInd) {
  out.LongVecTH.Rho <- LongVecTH.Rho(
    no.x = no.x, all.thres = all.thres,
    index.var.of.thres = index.var.of.thres,
    rho.xixj = rho.xixj
  )
  pi.xixj <- pairwiseExpProbVec(
    ind.vec = out.LongVecInd,
    th.rho.vec = out.LongVecTH.Rho
  )

  out.derLtoRho <- derLtoRho(
    ind.vec = out.LongVecInd,
    th.rho.vec = out.LongVecTH.Rho,
    n.xixj = n.xixj.vec, pi.xixj = pi.xixj, no.x = no.x
  )

  out.derLtoTau <- derLtoTau(
    ind.vec = out.LongVecInd,
    th.rho.vec = out.LongVecTH.Rho,
    n.xixj = n.xixj.vec, pi.xixj = pi.xixj,
    no.x = no.x
  )

  grad <- c(out.derLtoTau, out.derLtoRho)
  attr(grad, "pi.xixj") <- pi.xixj

  grad
}
################################################################################




# The input of the function  LongVecInd:

# no.x is scalar, the number of ordinal variables

# all.thres is vector containing the thresholds of all variables in the
# following order: thres_var1, thres_var2,..., thres_var_p
# within each variable the thresholds are in ascending order
# Note that all.thres does NOT contain the first and the last threshold of the
# variables, i.e. tau_0=-Inf and tau_last=Inf

# index.var.of.thres is a vector keeping track to which variable the thresholds
# in all.thres belongs to, it is of the form (1,1,1..., 2,2,2,...,  p,p,p,...)

# The output of the function:
# it is a list of vectors keeping track of the indices
# of thresholds, of variables, and of pairs, and two T/F vectors indicating
# if the threshold index corresponds to the last threshold of a variable; all
# these for all pairs of variables. All are needed for the
# computation of expected probabilities, der.L.to.rho, and der.L.to.tau

# all duplications of indices are done as follows: within each pair of variables,
# xi-xj, if for example we want to duplicate the indices of the thresholds,
# tau^xi_a and tau^xj_b, then index a runs faster than b, i.e. for each b we
# take all different tau^xi's, and then we proceed to the next b and do the
# same. In other words if it was tabulated we fill the table columnwise.

# All pairs xi-xj are taken with index j running faster than i.

# Note that each variable may have a different number of categories, that's why
# for example we take lists below.

LongVecInd <- function(no.x, all.thres, index.var.of.thres) {
  no.thres.of.each.var <- tapply(all.thres, index.var.of.thres, length)
  index.pairs <- utils::combn(no.x, 2)
  no.pairs <- ncol(index.pairs)

  # index.thres.var1.of.pair and index.thres.var2.of.pair contain the indices of
  # of all thresholds (from tau_0 which is -Inf to tau_last which is Inf)
  # for any pair of variables appropriately duplicated so that the two vectors
  # together give all possible combinations of thresholds indices
  # Since here the threshold indices 0 and "last" are included, the vectors are
  # longer than the vectors thres.var1.of.pair and thres.var2.of.pair above.
  index.thres.var1.of.pair <- vector("list", no.pairs)
  index.thres.var2.of.pair <- vector("list", no.pairs)

  # index.var1.of.pair and index.var2.of.pair keep track the index of the
  # variable that the thresholds in index.thres.var1.of.pair and
  # index.thres.var2.of.pair belong to, respectively. So, these two variables
  # are of same length as that of index.thres.var1.of.pair and
  # index.thres.var2.of.pair
  index.var1.of.pair <- vector("list", no.pairs)
  index.var2.of.pair <- vector("list", no.pairs)

  # index.pairs.extended gives the index of the pair for each pair of variables
  # e.g. pair of variables 1-2 has index 1, variables 1-3 has index 2, etc.
  # The vector is of the same length as index.thres.var1.of.pair,
  # index.thres.var2.of.pair, index.var1.of.pair, and index.var2.of.pair
  index.pairs.extended <- vector("list", no.pairs)

  for (i in 1:no.pairs) {
    no.thres.var1.of.pair <- no.thres.of.each.var[index.pairs[1, i]]
    no.thres.var2.of.pair <- no.thres.of.each.var[index.pairs[2, i]]

    index.thres.var1.of.pair[[i]] <- rep(0:(no.thres.var1.of.pair + 1),
      times = (no.thres.var2.of.pair + 2)
    )
    index.thres.var2.of.pair[[i]] <- rep(0:(no.thres.var2.of.pair + 1),
      each = (no.thres.var1.of.pair + 2)
    )
    length.vec <- length(index.thres.var1.of.pair[[i]])
    index.var1.of.pair[[i]] <- rep(index.pairs[1, i], length.vec)
    index.var2.of.pair[[i]] <- rep(index.pairs[2, i], length.vec)
    index.pairs.extended[[i]] <- rep(i, length.vec)
  }

  index.thres.var1.of.pair <- unlist(index.thres.var1.of.pair)
  index.thres.var2.of.pair <- unlist(index.thres.var2.of.pair)
  index.var1.of.pair <- unlist(index.var1.of.pair)
  index.var2.of.pair <- unlist(index.var2.of.pair)
  index.pairs.extended <- unlist(index.pairs.extended)

  # indicator vector (T/F) showing which elements of index.thres.var1.of.pair
  # correspond to the last thresholds of variables. The length is the same as
  # that of index.thres.var1.of.pair.
  last.thres.var1.of.pair <- index.var1.of.pair == 1 &
    index.thres.var1.of.pair == (no.thres.of.each.var[1] + 1)
  # we consider up to variable (no.x-1) because in pairs xi-xj where j runs
  # faster than i, the last variable is not included in the column of xi's
  for (i in 2:(no.x - 1)) {
    new.condition <- index.var1.of.pair == i &
      index.thres.var1.of.pair == (no.thres.of.each.var[i] + 1)
    last.thres.var1.of.pair <- last.thres.var1.of.pair | new.condition
  }

  # indicator vector (T/F) showing which elements of index.thres.var2.of.pair
  # correspond to the last thresholds of variables. Notet that in pairs xi-xj
  # where j runs faster than i, the first variable is not included in the column
  # of xj's. That's why we start with variable 2. The length is the same as
  # that of index.thres.var1.of.pair.
  last.thres.var2.of.pair <- index.var2.of.pair == 2 &
    index.thres.var2.of.pair == (no.thres.of.each.var[2] + 1)
  for (i in 3:no.x) {
    new.condition <- index.var2.of.pair == i &
      index.thres.var2.of.pair == (no.thres.of.each.var[i] + 1)
    last.thres.var2.of.pair <- last.thres.var2.of.pair | new.condition
  }

  list(
    index.thres.var1.of.pair = index.thres.var1.of.pair,
    index.thres.var2.of.pair = index.thres.var2.of.pair,
    index.var1.of.pair = index.var1.of.pair,
    index.var2.of.pair = index.var2.of.pair,
    index.pairs.extended = index.pairs.extended,
    last.thres.var1.of.pair = last.thres.var1.of.pair,
    last.thres.var2.of.pair = last.thres.var2.of.pair
  )
}
################################################################################


# The input of the function  LongVecTH.Rho:

# no.x is scalar, the number of ordinal variables

# all.thres is vector containing the thresholds of all variables in the
# following order: thres_var1, thres_var2,..., thres_var_p
# within each variable the thresholds are in ascending order
# Note that all.thres does NOT contain the first and the last threshold of the
# variables, i.e. tau_0=-Inf and tau_last=Inf

# index.var.of.thres is a vector keeping track to which variable the thresholds
# in all.thres belongs to, it is of the form (1,1,1..., 2,2,2,...,  p,p,p,...)

# rho.xixj is the vector of all corrlations where j runs faster than i
# i.e. the order is rho_12, rho_13, ..., rho_1p, rho_23, ..., rho_2p, etc.

# The output of the function:
# it is a list of vectors with thresholds and rho's duplicated appropriately,
# all needed for the computation of expected probabilities,
# der.L.to.rho, and der.L.to.tau

# all duplications below are done as follows: within each pair of variables,
# xi-xj, if for example we want to duplicate their thresholds, tau^xi_a and
# tau^xj_b, then index a runs faster than b, i.e. for each b we take all
# different tau^xi's, and then we proceed to the next b and do the same.
# In other words if it was tabulated we fill the table columnwise.

# All pairs xi-xj are taken with index j running faster than i.

# Note that each variable may have a different number of categories, that's why
# for example we take lists below.

LongVecTH.Rho <- function(no.x, all.thres, index.var.of.thres, rho.xixj) {
  no.thres.of.each.var <- tapply(all.thres, index.var.of.thres, length)
  index.pairs <- utils::combn(no.x, 2)
  no.pairs <- ncol(index.pairs)

  # create the long vectors needed for the computation of expected probabilities
  # for each cell and each pair of variables. The vectors thres.var1.of.pair and
  # thres.var2.of.pair together give all the possible combinations of the
  # thresholds of any two variables. Note the combinations (-Inf, -Inf),
  # (-Inf, Inf), (Inf, -Inf), (Inf, Inf) are NOT included. Only the combinations
  # of the middle thresholds (tau_1 to tau_(last-1)).
  # thres.var1.of.pair and thres.var2.of.pair give the first and the second
  # argument, respectively, in functions pbivnorm and dbinorm
  thres.var1.of.pair <- vector("list", no.pairs)
  thres.var2.of.pair <- vector("list", no.pairs)

  # Extending the rho.vector accordingly so that it will be the the third
  # argument in pbivnorm and dbinorm functions. It is of same length as
  # thres.var1.of.pair and thres.var2.of.pair.
  rho.vector <- vector("list", no.pairs)

  # thres.var1.for.dnorm.in.der.pi.to.tau.xi and
  # thres.var2.for.dnorm.in.der.pi.to.tau.xj give the thresholds of almost
  # all variables appropriately duplicated so that the vectors can be used
  # as input in dnorm() to compute der.pi.xixj.to.tau.xi and
  # der.pi.xixj.to.tau.xj.
  # thres.var1.for.dnorm.in.der.pi.to.tau.xi does not contain the thresholds of
  # the last variable and thres.var2.for.dnorm.in.der.pi.to.tau.xj those of
  # the first variable
  thres.var1.for.dnorm.in.der.pi.to.tau.xi <- vector("list", no.pairs)
  thres.var2.for.dnorm.in.der.pi.to.tau.xj <- vector("list", no.pairs)

  for (i in 1:no.pairs) {
    single.thres.var1.of.pair <- all.thres[index.var.of.thres == index.pairs[1, i]]
    single.thres.var2.of.pair <- all.thres[index.var.of.thres == index.pairs[2, i]]
    # remember that the first (-Inf) and last (Inf) thresholds are not included
    # so no.thres.var1.of.pair is equal to number of categories of var1 minus 1
    # similarly for no.thres.var2.of.pair
    no.thres.var1.of.pair <- no.thres.of.each.var[index.pairs[1, i]]
    no.thres.var2.of.pair <- no.thres.of.each.var[index.pairs[2, i]]

    thres.var1.of.pair[[i]] <- rep(single.thres.var1.of.pair,
      times = no.thres.var2.of.pair
    )
    thres.var2.of.pair[[i]] <- rep(single.thres.var2.of.pair,
      each = no.thres.var1.of.pair
    )
    rho.vector[[i]] <- rep(rho.xixj[i], length(thres.var1.of.pair[[i]]))

    thres.var1.for.dnorm.in.der.pi.to.tau.xi[[i]] <-
      rep(single.thres.var1.of.pair, times = (no.thres.var2.of.pair + 1))
    thres.var2.for.dnorm.in.der.pi.to.tau.xj[[i]] <-
      rep(single.thres.var2.of.pair, each = (no.thres.var1.of.pair + 1))
  }

  thres.var1.of.pair <- unlist(thres.var1.of.pair)
  thres.var2.of.pair <- unlist(thres.var2.of.pair)
  rho.vector <- unlist(rho.vector)
  thres.var1.for.dnorm.in.der.pi.to.tau.xi <-
    unlist(thres.var1.for.dnorm.in.der.pi.to.tau.xi)
  thres.var2.for.dnorm.in.der.pi.to.tau.xj <-
    unlist(thres.var2.for.dnorm.in.der.pi.to.tau.xj)

  # thres.var2.for.last.cat.var1 and thres.var1.for.last.cat.var2 are needed
  # for the computation of expected probabilities. In the computation of
  # \Phi_2(tau1, tau2; rho) when either tau1 or tau2 are Inf then it is enought
  # to compute pnorm() with the non-infinite tau as an argument
  # In particular when the first variable of the pair has tau_last= Inf
  # and the second a non-infite threshold we compute
  # pnorm(thres.var2.for.last.cat.var1). Similarly, when the second variable of
  # the pair has tau_last=Inf and the first a non-infite threshold we compute
  # pnorm(thres.var1.for.last.cat.var2).
  thres.var2.for.last.cat.var1 <- vector("list", (no.x - 1))
  thres.var1.for.last.cat.var2 <- vector("list", (no.x - 1))
  for (i in 1:(no.x - 1)) {
    thres.var2.for.last.cat.var1[[i]] <-
      c(all.thres[index.var.of.thres %in% (i + 1):no.x])
    thres.var1.for.last.cat.var2[[i]] <- rep(all.thres[index.var.of.thres == i],
      times = (no.x - i)
    )
  }
  thres.var2.for.last.cat.var1 <- unlist(thres.var2.for.last.cat.var1)
  thres.var1.for.last.cat.var2 <- unlist(thres.var1.for.last.cat.var2)


  list(
    thres.var1.of.pair = thres.var1.of.pair, # these 3 of same length
    thres.var2.of.pair = thres.var2.of.pair,
    rho.vector = rho.vector,

    # the following of length dependning on the number of categories
    thres.var1.for.dnorm.in.der.pi.to.tau.xi =
      thres.var1.for.dnorm.in.der.pi.to.tau.xi,
    thres.var2.for.dnorm.in.der.pi.to.tau.xj =
      thres.var2.for.dnorm.in.der.pi.to.tau.xj,
    thres.var2.for.last.cat.var1 = thres.var2.for.last.cat.var1,
    thres.var1.for.last.cat.var2 = thres.var1.for.last.cat.var2
  )
}
#########################################################



#########################################################

# The function  pairwiseExpProbVec
# input: ind.vec - the output of function LongVecInd
#        th.rho.vec - the output of function LongVecTH.Rho
# output: it gives the elements of pairwiseTablesExpected()$pi.tables
# table-wise and column-wise within each table. In other words if
# pi^xixj_ab is the expected probability for the pair of variables xi-xj
# and categories a and b, then index a runs the fastest of all, followed by b,
# then by j, and lastly by i.

pairwiseExpProbVec <- function(ind.vec, th.rho.vec) {
  prob.vec <- rep(NA, length(ind.vec$index.thres.var1.of.pair))

  prob.vec[ind.vec$index.thres.var1.of.pair == 0 |
    ind.vec$index.thres.var2.of.pair == 0] <- 0

  prob.vec[ind.vec$last.thres.var1.of.pair &
    ind.vec$last.thres.var2.of.pair] <- 1

  prob.vec[ind.vec$last.thres.var1.of.pair &
    ind.vec$index.thres.var2.of.pair != 0 &
    !ind.vec$last.thres.var2.of.pair] <-
    pnorm(th.rho.vec$thres.var2.for.last.cat.var1)

  prob.vec[ind.vec$last.thres.var2.of.pair &
    ind.vec$index.thres.var1.of.pair != 0 &
    !ind.vec$last.thres.var1.of.pair] <-
    pnorm(th.rho.vec$thres.var1.for.last.cat.var2)

  prob.vec[is.na(prob.vec)] <- pbivnorm(
    th.rho.vec$thres.var1.of.pair,
    th.rho.vec$thres.var2.of.pair,
    th.rho.vec$rho.vector
  )

  cum.term1 <- prob.vec[ind.vec$index.thres.var1.of.pair != 0 &
    ind.vec$index.thres.var2.of.pair != 0]

  cum.term2 <- prob.vec[ind.vec$index.thres.var1.of.pair != 0 &
    !ind.vec$last.thres.var2.of.pair]

  cum.term3 <- prob.vec[ind.vec$index.thres.var2.of.pair != 0 &
    !ind.vec$last.thres.var1.of.pair]

  cum.term4 <- prob.vec[!ind.vec$last.thres.var1.of.pair &
    !ind.vec$last.thres.var2.of.pair]

  PI <- cum.term1 - cum.term2 - cum.term3 + cum.term4

  # added by YR 11 nov 2012 to avoid Nan/-Inf
  # log(.Machine$double.eps) = -36.04365
  # all elements should be strictly positive
  PI[PI < .Machine$double.eps] <- .Machine$double.eps

  PI
}


# derLtoRho
# input: ind.vec - the output of function LongVecInd
#        th.rho.vec - the output of function LongVecTH.Rho
#        n.xixj - a vector with the observed frequency for every combination
#                 of categories and every pair. The frequencies are given in
#                 the same order as the expected probabilities in the output of
#                 pairwiseExpProbVec output
#        pi.xixj - the output of pairwiseExpProbVec function
#        no.x    - the number of ordinal variables
# output: the vector of der.L.to.rho, each element corresponds to
#         der.Lxixj.to.rho.xixj where j runs faster than i

derLtoRho <- function(ind.vec, th.rho.vec, n.xixj, pi.xixj, no.x) {
  prob.vec <- rep(NA, length(ind.vec$index.thres.var1.of.pair))

  prob.vec[ind.vec$index.thres.var1.of.pair == 0 |
    ind.vec$index.thres.var2.of.pair == 0 |
    ind.vec$last.thres.var1.of.pair |
    ind.vec$last.thres.var2.of.pair] <- 0

  prob.vec[is.na(prob.vec)] <- dbinorm(th.rho.vec$thres.var1.of.pair,
    th.rho.vec$thres.var2.of.pair,
    rho = th.rho.vec$rho.vector
  )

  den.term1 <- prob.vec[ind.vec$index.thres.var1.of.pair != 0 &
    ind.vec$index.thres.var2.of.pair != 0]

  den.term2 <- prob.vec[ind.vec$index.thres.var1.of.pair != 0 &
    !ind.vec$last.thres.var2.of.pair]

  den.term3 <- prob.vec[ind.vec$index.thres.var2.of.pair != 0 &
    !ind.vec$last.thres.var1.of.pair]

  den.term4 <- prob.vec[!ind.vec$last.thres.var1.of.pair &
    !ind.vec$last.thres.var2.of.pair]

  der.pi.xixj.to.rho.xixj <- den.term1 - den.term2 - den.term3 + den.term4
  prod.terms <- (n.xixj / pi.xixj) * der.pi.xixj.to.rho.xixj

  # to get der.Lxixj.to.rho.xixj we should all the elements of
  # der.pi.xixj.to.rho.xixj which correspond to the pair xi-xj, to do so:
  xnew <- lapply(
    ind.vec[c("index.pairs.extended")],
    function(y) {
      y[ind.vec$index.thres.var1.of.pair != 0 &
        ind.vec$index.thres.var2.of.pair != 0]
    }
  )
  # der.L.to.rho is:
  tapply(prod.terms, xnew$index.pairs.extended, sum)
}
###########################################################################


# derLtoTau
# input: ind.vec - the output of function LongVecInd
#        th.rho.vec - the output of function LongVecTH.Rho
#        n.xixj - a vector with the observed frequency for every combination
#                 of categories and every pair. The frequencies are given in
#                 the same order as the expected probabilities in the output of
#                 pairwiseExpProbVec output
#        pi.xixj - the output of pairwiseExpProbVec function
# output: the vector of der.L.to.tau where the elements are ordered as follows:
#         the thresholds of each variable with respect to ascending order of
#         the variable index (i.e. thres_var1, thres_var2, etc.) and within
#         each variable the thresholds in ascending order.


derLtoTau <- function(ind.vec, th.rho.vec, n.xixj, pi.xixj, no.x = 0L) {
  # to compute der.pi.xixj.to.tau.xi
  xi <- lapply(
    ind.vec[c(
      "index.thres.var2.of.pair",
      "last.thres.var2.of.pair"
    )],
    function(y) {
      y[!(ind.vec$index.thres.var1.of.pair == 0 |
        ind.vec$last.thres.var1.of.pair)]
    }
  )

  cum.prob.vec <- rep(NA, length(xi$index.thres.var2.of.pair))
  cum.prob.vec[xi$index.thres.var2.of.pair == 0] <- 0
  cum.prob.vec[xi$last.thres.var2.of.pair] <- 1
  denom <- sqrt(1 - (th.rho.vec$rho.vector * th.rho.vec$rho.vector))
  cum.prob.vec[is.na(cum.prob.vec)] <-
    pnorm((th.rho.vec$thres.var2.of.pair -
      th.rho.vec$rho.vector * th.rho.vec$thres.var1.of.pair) /
      denom)
  den.prob.vec <- dnorm(th.rho.vec$thres.var1.for.dnorm.in.der.pi.to.tau.xi)
  der.pi.xixj.to.tau.xi <- den.prob.vec *
    (cum.prob.vec[xi$index.thres.var2.of.pair != 0] -
      cum.prob.vec[!xi$last.thres.var2.of.pair])

  # to compute der.pi.xixj.to.tau.xj
  xj <- lapply(
    ind.vec[c(
      "index.thres.var1.of.pair",
      "last.thres.var1.of.pair"
    )],
    function(y) {
      y[!(ind.vec$index.thres.var2.of.pair == 0 |
        ind.vec$last.thres.var2.of.pair)]
    }
  )

  cum.prob.vec <- rep(NA, length(xj$index.thres.var1.of.pair))
  cum.prob.vec[xj$index.thres.var1.of.pair == 0] <- 0
  cum.prob.vec[xj$last.thres.var1.of.pair] <- 1
  denom <- sqrt(1 - (th.rho.vec$rho.vector * th.rho.vec$rho.vector))
  cum.prob.vec[is.na(cum.prob.vec)] <-
    pnorm((th.rho.vec$thres.var1.of.pair -
      th.rho.vec$rho.vector * th.rho.vec$thres.var2.of.pair) /
      denom)
  den.prob.vec <- dnorm(th.rho.vec$thres.var2.for.dnorm.in.der.pi.to.tau.xj)
  der.pi.xixj.to.tau.xj <- den.prob.vec *
    (cum.prob.vec[xj$index.thres.var1.of.pair != 0] -
      cum.prob.vec[!xj$last.thres.var1.of.pair])

  # to compute der.Lxixj.tau.xi and der.Lxixj.tau.xi
  n.over.pi <- n.xixj / pi.xixj
  # get the appropriate differences of n.over.pi for der.Lxixj.to.tau.xi  and
  # der.Lxixj.to.tau.xj
  x3a <- lapply(ind.vec, function(y) {
    y[!(ind.vec$index.thres.var1.of.pair == 0 |
      ind.vec$index.thres.var2.of.pair == 0)]
  })
  diff.n.over.pi.to.xi <- n.over.pi[!x3a$last.thres.var1.of.pair] -
    n.over.pi[x3a$index.thres.var1.of.pair != 1]
  diff.n.over.pi.to.xj <- n.over.pi[!x3a$last.thres.var2.of.pair] -
    n.over.pi[x3a$index.thres.var2.of.pair != 1]
  # terms.der.Lxixj.to.tau.xi  and terms.der.Lxixj.to.tau.xj
  terms.der.Lxixj.to.tau.xi <- diff.n.over.pi.to.xi * der.pi.xixj.to.tau.xi
  terms.der.Lxixj.to.tau.xj <- diff.n.over.pi.to.xj * der.pi.xixj.to.tau.xj
  # to add appropriately elements of terms.der.Lxixj.to.tau.xi
  x3b <- lapply(
    ind.vec[c(
      "index.pairs.extended",
      "index.thres.var1.of.pair"
    )],
    function(y) {
      y[!(ind.vec$index.thres.var1.of.pair == 0 |
        ind.vec$last.thres.var1.of.pair |
        ind.vec$index.thres.var2.of.pair == 0)]
    }
  )
  # to add appropriately elements of terms.der.Lxixj.to.tau.xj
  x4b <- lapply(
    ind.vec[c(
      "index.pairs.extended",
      "index.thres.var2.of.pair"
    )],
    function(y) {
      y[!(ind.vec$index.thres.var2.of.pair == 0 |
        ind.vec$last.thres.var2.of.pair |
        ind.vec$index.thres.var1.of.pair == 0)]
    }
  )
  ind.pairs <- utils::combn(no.x, 2)
  # der.Lxixj.to.tau.xi is a matrix, nrow=no.pairs, ncol=max(no.of.free.thres)
  # thus, there are NA's, similarly for der.Lxixj.to.tau.xj
  der.Lxixj.to.tau.xi <- tapply(
    terms.der.Lxixj.to.tau.xi,
    list(
      x3b$index.pairs.extended,
      x3b$index.thres.var1.of.pair
    ),
    sum
  )
  der.Lxixj.to.tau.xj <- tapply(
    terms.der.Lxixj.to.tau.xj,
    list(
      x4b$index.pairs.extended,
      x4b$index.thres.var2.of.pair
    ),
    sum
  )

  # to add appropriately the terms of der.Lxixj.to.tau.xi and
  # der.Lxixj.to.tau.xj
  split.der.Lxixj.to.tau.xi <- split(
    as.data.frame(der.Lxixj.to.tau.xi),
    ind.pairs[1, ]
  )
  sums.der.Lxixj.to.tau.xi <- lapply(
    split.der.Lxixj.to.tau.xi,
    function(x) {
      y <- apply(x, 2, sum)
      y[!is.na(y)]
    }
  )
  # Note: NA exist in the case where the ordinal variables have different
  # number of response categories
  split.der.Lxixj.to.tau.xj <- split(
    as.data.frame(der.Lxixj.to.tau.xj),
    ind.pairs[2, ]
  )
  sums.der.Lxixj.to.tau.xj <- lapply(
    split.der.Lxixj.to.tau.xj,
    function(x) {
      y <- apply(x, 2, sum)
      y[!is.na(y)]
    }
  )
  # to get der.L.to.tau
  c(
    sums.der.Lxixj.to.tau.xi[[1]],
    c(unlist(sums.der.Lxixj.to.tau.xi[2:(no.x - 1)]) +
      unlist(sums.der.Lxixj.to.tau.xj[1:(no.x - 2)])),
    sums.der.Lxixj.to.tau.xj[[no.x - 1]]
  )
}
