# utility functions for pairwise maximum likelihood

# stub for lav_pml_fml_dploglik_dimplied
lav_pml_fml_dploglik_dimplied <- function(
                       sigma_hat = NULL, # model-based var/cov/cor
                       th = NULL, # model-based thresholds + means
                       th_idx = NULL, # threshold idx per variable
                       num_idx = NULL, # which variables are numeric
                       x = NULL, # data
                       exo = NULL, # external covariates
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
lav_pml_dploglik_dimplied <- function(
                       sigma_hat = NULL, # model-based var/cov/cor
                       mu_hat = NULL, # model-based means
                       th = NULL, # model-based thresholds + means
                       th_idx = NULL, # threshold idx per variable
                       num_idx = NULL, # which variables are numeric
                       x = NULL, # data
                       exo = NULL, # external covariates
                       wt = NULL, # case weights (not used yet)
                       lavcache = NULL, # housekeeping stuff
                       pi0 = NULL, # slopes
                       missing = "listwise", # how to deal with missings
                       scores = FALSE, # return case-wise scores
                       negative = TRUE) { # multiply by -1

  # diagonal of Sigma.hat is not necessarily 1, even for categorical vars
  sigma_hat2 <- sigma_hat
  if (length(num_idx) > 0L) {
    diag(sigma_hat2)[-num_idx] <- 1
  } else {
    diag(sigma_hat2) <- 1
  }
  cor_hat <- cov2cor(sigma_hat2) # to get correlations (rho!)
  cors <- lav_mat_vech(cor_hat, diagonal = FALSE)

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

  nvar <- nrow(sigma_hat)
  pstar <- nvar * (nvar - 1) / 2
  ov_types <- rep("ordered", nvar)
  if (length(num_idx) > 0L) ov_types[num_idx] <- "numeric"
  if (!is.null(exo)) {
    nexo <- ncol(exo)
  } else {
    nexo <- 0
  }


  if (all(ov_types == "numeric")) {
    n_th <- nvar
  } else {
    n_th <- length(th_idx)
  }
  n_sl <- nvar * nexo
  n_var <- length(num_idx)
  n_cor <- pstar

  # add num.idx to th.idx
  if (length(num_idx) > 0L) {
    th_idx[th_idx == 0] <- num_idx
  }

  # print(Sigma.hat); print(TH); print(th.idx); print(num.idx); print(str(X))

  # shortcut for ordinal-only/no-exo case
  if (!scores && all(ov_types == "ordered") && nexo == 0L) {
    # >>>>>>>> HJ/MK PML CODE >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    if (is.null(wt)) {
      n_xixj_vec <- lavcache$bifreq
    } else {
      n_xixj_vec <- lavcache$sum_obs_weights_xixj_ab_vec
    }
    gradient <- lav_pml_grad_tau_rho(
      no_x = nvar,
      all_thres = th,
      index_var_of_thres = th_idx,
      rho_xixj = cors,
      n_xixj_vec = n_xixj_vec,
      out_lav_pml_longvec_ind = lavcache$long
    )

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    if (missing == "available.cases") {
      uni_pi <- lav_pml_th_uni_prob(th = th, th_idx = th_idx)
      tmp <- lavcache$uniweights / uni_pi

      var_idx <- split(th_idx, th_idx)
      var_idx <- unlist(lapply(var_idx, function(x_1) {
        c(x_1, x_1[1])
      }))

      tmp_varwise <- split(tmp, var_idx)
      tmp1 <- unlist(lapply(
        tmp_varwise,
        function(x_1) {
          c(x_1[-length(x_1)])
        }
      ))
      tmp2 <- unlist(lapply(tmp_varwise, function(x_1) {
        c(x_1[-1])
      }))

      uni_der_tau <- dnorm(th) * (tmp1 - tmp2)
      n_th_1 <- length(th)
      gradient[1:n_th_1] <- gradient[1:n_th_1] + uni_der_tau
    }

    if (negative) {
      gradient <- -1 * gradient
    }
    return(gradient)
  }

  # in this order: TH/MEANS + SLOPES + VAR + COR
  grad_size <- n_th + n_sl + n_var + n_cor

  # scores or gradient?
  if (scores) {
    scores_1 <- matrix(0, nrow(x), grad_size) # we will sum up over all pairs
  } else {
    grad <- matrix(0, pstar, grad_size) # each pair is a row
  }
  pstar_1 <- matrix(0, nvar, nvar) # utility matrix, to get indices
  pstar_1[lav_mat_vech_idx(nvar, diagonal = FALSE)] <- 1:pstar
  # n <- length(x[, 1])

  for (j in seq_len(nvar - 1L)) {
    for (i in (j + 1L):nvar) {
      # cat(" i = ", i, " j = ", j, "\n") # debug only
      pstar_idx <- pstar_1[i, j]
      cor_idx <- n_th + n_sl + n_var + pstar_1[i, j]
      th_idx_i <- which(th_idx == i)
      th_idx_j <- which(th_idx == j)
      if (nexo > 0L) {
        sl_idx_i <- n_th + seq(i, by = nvar, length.out = nexo)
        sl_idx_j <- n_th + seq(j, by = nvar, length.out = nexo)

        if (length(num_idx) > 0L) {
          var_idx_i <- n_th + n_sl + match(i, num_idx)
          var_idx_j <- n_th + n_sl + match(j, num_idx)
        }
      } else {
        if (length(num_idx) > 0L) {
          var_idx_i <- n_th + match(i, num_idx)
          var_idx_j <- n_th + match(j, num_idx)
        }
      }
      if (ov_types[i] == "numeric" && ov_types[j] == "numeric") {
        if (nexo > 1L) {
          lav_msg_stop(gettext(
            "mixed + exo in PML not implemented;
            try optim.gradient = \"numerical\""))
        }

        sc <- lav_mvn_sc_mu_sigma(
          y = x[, c(i, j)],
          mu = mu_hat[c(i, j)], sigma_1 = sigma_hat[c(i, j), c(i, j)]
        )

        if (scores) {
          if (all(ov_types == "numeric") && nexo == 0L) {
            # MU1 + MU2
            scores_1[, c(i, j)] <- scores_1[, c(i, j)] + sc[, c(1, 2)]
            # VAR1 + COV_12 + VAR2
            var_idx <- (nvar +
              lav_mat_vech_match_idx(nvar, idx = c(i, j)))
            scores_1[, var_idx] <- scores_1[, var_idx] + sc[, c(3, 4, 5)]
          } else { # mixed ordered/continuous
            # MU
            mu_idx <- c(th_idx_i, th_idx_j)
            scores_1[, mu_idx] <- scores_1[, mu_idx] + (-1) * sc[, c(1, 2)]
            # VAR+COV
            var_idx <- c(var_idx_i, cor_idx, var_idx_j)
            scores_1[, var_idx] <- scores_1[, var_idx] + sc[, c(3, 4, 5)]
          }
        } else {
          if (all(ov_types == "numeric") && nexo == 0L) {
            mu_idx <- c(i, j)
            sigma_idx <- (nvar +
              lav_mat_vech_match_idx(nvar, idx = c(i, j)))
            # MU1 + MU2
            grad[pstar_idx, mu_idx] <-
              colSums(sc[, c(1, 2)], na.rm = TRUE)
          } else {
            mu_idx <- c(th_idx_i, th_idx_j)
            sigma_idx <- c(var_idx_i, cor_idx, var_idx_j)
            # MU (reverse sign!)
            grad[pstar_idx, mu_idx] <-
              -1 * colSums(sc[, c(1, 2)], na.rm = TRUE)
          }
          # SIGMA
          grad[pstar_idx, sigma_idx] <-
            colSums(sc[, c(3, 4, 5)], na.rm = TRUE)
        } # gradient only
      } else if (ov_types[i] == "numeric" && ov_types[j] == "ordered") {
        # polyserial correlation
        if (nexo > 1L) {
          lav_msg_stop(gettext(
            "mixed + exo in PML not implemented;
            try optim.gradient = \"numerical\""))
        }

        sc_cor_uni <- lav_bvmix_cor_sc(
          y1 = x[, i], y2 = x[, j],
          exo = NULL, wt = wt,
          evar_y1 = sigma_hat[i, i],
          beta_y1 = mu_hat[i],
          th_y2 = th[th_idx == j],
          sl_y2 = NULL,
          rho = cor_hat[i, j],
          sigma_correction = TRUE
        )

        if (scores) {
          # MU
          scores_1[, th_idx_i] <- (scores_1[, th_idx_i] +
            -1 * sc_cor_uni$dx_mu_y1)
          # TH
          scores_1[, th_idx_j] <- (scores_1[, th_idx_j] +
            sc_cor_uni$dx_th_y2)
          # VAR
          scores_1[, var_idx_i] <- (scores_1[, var_idx_i] +
            sc_cor_uni$dx_var_y1)
          # COR
          scores_1[, cor_idx] <- (scores_1[, cor_idx] +
            sc_cor_uni$dx_rho)
        } else {
          # MU
          grad[pstar_idx, th_idx_i] <-
            -1 * sum(sc_cor_uni$dx_mu_y1, na.rm = TRUE)
          # TH
          grad[pstar_idx, th_idx_j] <-
            colSums(sc_cor_uni$dx_th_y2, na.rm = TRUE)
          # VAR
          grad[pstar_idx, var_idx_i] <-
            sum(sc_cor_uni$dx_var_y1, na.rm = TRUE)
          # COR
          grad[pstar_idx, cor_idx] <-
            sum(sc_cor_uni$dx_rho, na.rm = TRUE)
        } # grad only
      } else if (ov_types[j] == "numeric" && ov_types[i] == "ordered") {
        # polyserial correlation
        if (nexo > 1L) {
          lav_msg_stop(gettext(
            "mixed + exo in PML not implemented;
            try optim.gradient = \"numerical\""))
        }

        sc_cor_uni <- lav_bvmix_cor_sc(
          y1 = x[, j], y2 = x[, i],
          exo = NULL, wt = wt,
          evar_y1 = sigma_hat[j, j],
          beta_y1 = mu_hat[j],
          th_y2 = th[th_idx == i],
          rho = cor_hat[i, j],
          sigma_correction = TRUE
        )

        if (scores) {
          # MU
          scores_1[, th_idx_j] <- (scores_1[, th_idx_j] +
            -1 * sc_cor_uni$dx_mu_y1)
          # TH
          scores_1[, th_idx_i] <- (scores_1[, th_idx_i] +
            sc_cor_uni$dx_th_y2)
          # VAR
          scores_1[, var_idx_j] <- (scores_1[, var_idx_j] +
            sc_cor_uni$dx_var_y1)
          # COR
          scores_1[, cor_idx] <- (scores_1[, cor_idx] +
            sc_cor_uni$dx_rho)
        } else {
          # MU
          grad[pstar_idx, th_idx_j] <-
            -1 * sum(sc_cor_uni$dx_mu_y1, na.rm = TRUE)
          # TH
          grad[pstar_idx, th_idx_i] <-
            colSums(sc_cor_uni$dx_th_y2, na.rm = TRUE)
          # VAR
          grad[pstar_idx, var_idx_j] <-
            sum(sc_cor_uni$dx_var_y1, na.rm = TRUE)
          # COR
          grad[pstar_idx, cor_idx] <-
            sum(sc_cor_uni$dx_rho, na.rm = TRUE)
        } # grad only
      } else if (ov_types[i] == "ordered" && ov_types[j] == "ordered") {
        # polychoric correlation
        if (nexo == 0L) {
          sc_cor_uni <-
            lav_bvord_cor_sc(
              y1 = x[, i], y2 = x[, j],
              exo = NULL, wt = wt,
              rho = sigma_hat[i, j],
              fit_y1 = NULL, # fixme
              fit_y2 = NULL, # fixme
              th_y1 = th[th_idx == i],
              th_y2 = th[th_idx == j],
              sl_y1 = NULL,
              sl_y2 = NULL,
              na_zero = TRUE
            )
        } else {
          sc_cor_uni <-
            lav_pml_dbilogl_dpar_x(
              y1 = x[, i],
              y2 = x[, j],
              exo = exo,
              rho = sigma_hat[i, j],
              th_y1 = th[th_idx == i],
              th_y2 = th[th_idx == j],
              sl_y1 = pi0[i, ],
              sl_y2 = pi0[j, ],
              missing_ind = missing
            )
        }

        if (scores) {
          # TH
          scores_1[, th_idx_i] <- scores_1[, th_idx_i] + sc_cor_uni$dx_th_y1
          scores_1[, th_idx_j] <- scores_1[, th_idx_j] + sc_cor_uni$dx_th_y2

          # SL
          if (nexo > 0L) {
            scores_1[, sl_idx_i] <- scores_1[, sl_idx_i] + sc_cor_uni$dx_sl_y1
            scores_1[, sl_idx_j] <- scores_1[, sl_idx_j] + sc_cor_uni$dx_sl_y2
          }
          # NO VAR
          # RHO
          scores_1[, cor_idx] <- scores_1[, cor_idx] + sc_cor_uni$dx_rho
        } else {
          # TH
          if (length(th_idx_i) > 1L) {
            grad[pstar_idx, th_idx_i] <-
              colSums(sc_cor_uni$dx_th_y1, na.rm = TRUE)
          } else {
            grad[pstar_idx, th_idx_i] <-
              sum(sc_cor_uni$dx_th_y1, na.rm = TRUE)
          }
          if (length(th_idx_j) > 1L) {
            grad[pstar_idx, th_idx_j] <-
              colSums(sc_cor_uni$dx_th_y2, na.rm = TRUE)
          } else {
            grad[pstar_idx, th_idx_j] <-
              sum(sc_cor_uni$dx_th_y2, na.rm = TRUE)
          }

          # SL
          if (nexo > 0L) {
            if (length(sl_idx_i) > 1L) {
              grad[pstar_idx, sl_idx_i] <-
                colSums(sc_cor_uni$dx_sl_y1, na.rm = TRUE)
            } else {
              grad[pstar_idx, sl_idx_i] <-
                sum(sc_cor_uni$dx_sl_y1, na.rm = TRUE)
            }
            if (length(sl_idx_j) > 1L) {
              grad[pstar_idx, sl_idx_j] <-
                colSums(sc_cor_uni$dx_sl_y2, na.rm = TRUE)
            } else {
              grad[pstar_idx, sl_idx_j] <-
                sum(sc_cor_uni$dx_sl_y2, na.rm = TRUE)
            }
          }
          # NO VAR

          # RHO
          grad[pstar_idx, cor_idx] <-
            sum(sc_cor_uni$dx_rho, na.rm = TRUE)
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

  if (missing == "available.cases" && all(ov_types == "ordered")) {
    # the available.cases objective adds a (weighted) UNIVARIATE likelihood
    # part to the pairwise likelihood; its scores/gradient must be added to
    # the threshold (and slope) entries. Note (fixed July 2026): the scores
    # branch used to (a) drop the univariate part entirely when nexo == 0,
    # and (b) overwrite scores_1 with only its first n_th + n_sl columns
    # when nexo > 0 (losing the correlation columns) -- which made the
    # casewise scores (and hence the first-order information) inconsistent
    # with the estimating function.
    if (nexo == 0L) {
      uni_scores <- matrix(0, nrow(x), n_th)
      for (i in seq_len(nvar)) {
        th_idx_i <- which(th_idx == i)
        der_y1 <- lav_pml_uni_sc(
          y1 = x[, i], th_y1 = th[th_idx == i],
          exo = NULL, sl_y1 = NULL,
          weights_casewise = lavcache$uniweights.casewise
        )
        uni_scores[, th_idx_i] <- der_y1$dx.th.y1
      }
      if (scores) {
        scores_1[, seq_len(n_th)] <-
          scores_1[, seq_len(n_th), drop = FALSE] + uni_scores
      } else {
        uni_gradient <- colSums(uni_scores)
      }
    } else {
      uni_scores <- matrix(0, nrow(x), ncol = (n_th + n_sl))
      for (i in seq_len(nvar)) {
        th_idx_i <- which(th_idx == i)
        sl_idx_i <- n_th + seq(i, by = nvar, length.out = nexo)
        der_y1 <- lav_pml_uni_sc(
          y1 = x[, i], th_y1 = th[th_idx == i],
          exo = exo, sl_y1 = pi0[i, ],
          weights_casewise = lavcache$uniweights.casewise
        )
        uni_scores[, th_idx_i] <- der_y1$dx.th.y1
        uni_scores[, sl_idx_i] <- der_y1$dx.sl.y1
      }
      if (scores) {
        scores_1[, seq_len(n_th + n_sl)] <-
          scores_1[, seq_len(n_th + n_sl), drop = FALSE] + uni_scores
      } else {
        uni_gradient <- colSums(uni_scores)
      }
    }
  }

  # do we need scores?
  if (scores) {
    return(scores_1)
  }

  # DEBUG
  # :print(GRAD)
  ###########



  # gradient is sum over all pairs
  gradient <- colSums(grad, na.rm = TRUE)

  if (missing == "available.cases" && all(ov_types == "ordered")) {
    # (fixed July 2026: the RHS used the FULL gradient vector instead of
    # its leading block, garbling the result through recycling; note that
    # the nexo == 0 case is normally unreachable here -- it is handled by
    # the shortcut at the top of this function)
    if (nexo == 0L) {
      gradient[seq_len(n_th)] <-
        gradient[seq_len(n_th)] + uni_gradient
    } else {
      gradient[seq_len(n_th + n_sl)] <-
        gradient[seq_len(n_th + n_sl)] + uni_gradient
    }
  }

  # we multiply by -1 because we minimize
  if (negative) {
    gradient <- -1 * gradient
  }

  gradient
}


### all code below written by Myrsini Katsikatsou


# The function lav_pml_grad_tau_rho
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
#               lav_pml_expprob_vec output
# out.lav_pml_longvec_ind - it is the output of function lav_pml_longvec_ind
# the output: it gives the elements of der.L.to.tau and der.L.to.rho in this
#             order. The elements of der.L.to.tau where the elements are
#             ordered as follows: the thresholds of each variable with respect
#             to ascending order of the variable index (i.e. thres_var1,
#             thres_var2, etc.) and within each variable the thresholds in
#             ascending order.
#             The elements of vector der.L.to.rho are der.Lxixj.to.rho.xixj
#             where j runs faster than i.

# The function depends on four other functions: lav_pml_longvec_th_rho,
# lav_pml_expprob_vec, lav_pml_dl_drho, and lav_pml_dl_dtau, all given below.

# if n.xixj.ab is either an array or a list the following should be done
# n.xixj.vec <- if(is.array(n.xixj.ab)) {
#                 c(n.xixj.ab)
#              } else if(is.list(n.xixj.ab)){
#                 unlist(n.xixj.ab)
#              }


lav_pml_grad_tau_rho <- function(no_x, all_thres, index_var_of_thres, rho_xixj,
                         n_xixj_vec, out_lav_pml_longvec_ind) {
  out_lav_pml_longvec_th_rho <- lav_pml_longvec_th_rho(
    no_x = no_x, all_thres = all_thres,
    index_var_of_thres = index_var_of_thres,
    rho_xixj = rho_xixj
  )
  pi_xixj <- lav_pml_expprob_vec(
    ind_vec = out_lav_pml_longvec_ind,
    th_rho_vec = out_lav_pml_longvec_th_rho
  )

  out_lav_pml_dl_drho <- lav_pml_dl_drho(
    ind_vec = out_lav_pml_longvec_ind,
    th_rho_vec = out_lav_pml_longvec_th_rho,
    n_xixj = n_xixj_vec, pi_xixj = pi_xixj, no_x = no_x
  )

  out_lav_pml_dl_dtau <- lav_pml_dl_dtau(
    ind_vec = out_lav_pml_longvec_ind,
    th_rho_vec = out_lav_pml_longvec_th_rho,
    n_xixj = n_xixj_vec, pi_xixj = pi_xixj,
    no_x = no_x
  )

  grad <- c(out_lav_pml_dl_dtau, out_lav_pml_dl_drho)
  attr(grad, "pi.xixj") <- pi_xixj

  grad
}
################################################################################




# The input of the function  lav_pml_longvec_ind:

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

# all duplications of indices are done as follows: within each pair of
# variables,  xi-xj, if for example we want to duplicate the indices of
# the thresholds, tau^xi_a and tau^xj_b, then index a runs faster than b,
# i.e. for each b we take all different tau^xi's, and then we proceed to
# the next b and do the same. In other words if it was tabulated we fill
# the table columnwise.

# All pairs xi-xj are taken with index j running faster than i.

# Note that each variable may have a different number of categories, that's why
# for example we take lists below.

lav_pml_longvec_ind <- function(no_x, all_thres, index_var_of_thres) {
  no_thres_of_each_var <- tapply(all_thres, index_var_of_thres, length)
  index_pairs <- utils::combn(no_x, 2)
  no_pairs <- ncol(index_pairs)

  # index.thres.var1.of.pair and index.thres.var2.of.pair contain the indices of
  # of all thresholds (from tau_0 which is -Inf to tau_last which is Inf)
  # for any pair of variables appropriately duplicated so that the two vectors
  # together give all possible combinations of thresholds indices
  # Since here the threshold indices 0 and "last" are included, the vectors are
  # longer than the vectors thres.var1.of.pair and thres.var2.of.pair above.
  index_thres_var1_of_pair <- vector("list", no_pairs)
  index_thres_var2_of_pair <- vector("list", no_pairs)

  # index.var1.of.pair and index.var2.of.pair keep track the index of the
  # variable that the thresholds in index.thres.var1.of.pair and
  # index.thres.var2.of.pair belong to, respectively. So, these two variables
  # are of same length as that of index.thres.var1.of.pair and
  # index.thres.var2.of.pair
  index_var1_of_pair <- vector("list", no_pairs)
  index_var2_of_pair <- vector("list", no_pairs)

  # index.pairs.extended gives the index of the pair for each pair of variables
  # e.g. pair of variables 1-2 has index 1, variables 1-3 has index 2, etc.
  # The vector is of the same length as index.thres.var1.of.pair,
  # index.thres.var2.of.pair, index.var1.of.pair, and index.var2.of.pair
  index_pairs_extended <- vector("list", no_pairs)

  for (i in 1:no_pairs) {
    no_thres_var1_of_pair <- no_thres_of_each_var[index_pairs[1, i]]
    no_thres_var2_of_pair <- no_thres_of_each_var[index_pairs[2, i]]

    index_thres_var1_of_pair[[i]] <- rep(0:(no_thres_var1_of_pair + 1),
      times = (no_thres_var2_of_pair + 2)
    )
    index_thres_var2_of_pair[[i]] <- rep(0:(no_thres_var2_of_pair + 1),
      each = (no_thres_var1_of_pair + 2)
    )
    length_vec <- length(index_thres_var1_of_pair[[i]])
    index_var1_of_pair[[i]] <- rep(index_pairs[1, i], length_vec)
    index_var2_of_pair[[i]] <- rep(index_pairs[2, i], length_vec)
    index_pairs_extended[[i]] <- rep(i, length_vec)
  }

  index_thres_var1_of_pair <- unlist(index_thres_var1_of_pair)
  index_thres_var2_of_pair <- unlist(index_thres_var2_of_pair)
  index_var1_of_pair <- unlist(index_var1_of_pair)
  index_var2_of_pair <- unlist(index_var2_of_pair)
  index_pairs_extended <- unlist(index_pairs_extended)

  # indicator vector (T/F) showing which elements of index.thres.var1.of.pair
  # correspond to the last thresholds of variables. The length is the same as
  # that of index.thres.var1.of.pair.
  last_thres_var1_of_pair <- index_var1_of_pair == 1 &
    index_thres_var1_of_pair == (no_thres_of_each_var[1] + 1)
  # we consider up to variable (no.x-1) because in pairs xi-xj where j runs
  # faster than i, the last variable is not included in the column of xi's
  for (i in 2:(no_x - 1)) {
    new_condition <- index_var1_of_pair == i &
      index_thres_var1_of_pair == (no_thres_of_each_var[i] + 1)
    last_thres_var1_of_pair <- last_thres_var1_of_pair | new_condition
  }

  # indicator vector (T/F) showing which elements of index.thres.var2.of.pair
  # correspond to the last thresholds of variables. Notet that in pairs xi-xj
  # where j runs faster than i, the first variable is not included in the column
  # of xj's. That's why we start with variable 2. The length is the same as
  # that of index.thres.var1.of.pair.
  last_thres_var2_of_pair <- index_var2_of_pair == 2 &
    index_thres_var2_of_pair == (no_thres_of_each_var[2] + 1)
  for (i in 3:no_x) {
    new_condition <- index_var2_of_pair == i &
      index_thres_var2_of_pair == (no_thres_of_each_var[i] + 1)
    last_thres_var2_of_pair <- last_thres_var2_of_pair | new_condition
  }

  list(
    index.thres.var1.of.pair = index_thres_var1_of_pair,
    index.thres.var2.of.pair = index_thres_var2_of_pair,
    index.var1.of.pair = index_var1_of_pair,
    index.var2.of.pair = index_var2_of_pair,
    index.pairs.extended = index_pairs_extended,
    last.thres.var1.of.pair = last_thres_var1_of_pair,
    last.thres.var2.of.pair = last_thres_var2_of_pair
  )
}
################################################################################


# The input of the function  lav_pml_longvec_th_rho:

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

lav_pml_longvec_th_rho <- function(no_x, all_thres,
                                   index_var_of_thres, rho_xixj) {
  no_thres_of_each_var <- tapply(all_thres, index_var_of_thres, length)
  index_pairs <- utils::combn(no_x, 2)
  no_pairs <- ncol(index_pairs)

  # create the long vectors needed for the computation of expected probabilities
  # for each cell and each pair of variables. The vectors thres.var1.of.pair and
  # thres.var2.of.pair together give all the possible combinations of the
  # thresholds of any two variables. Note the combinations (-Inf, -Inf),
  # (-Inf, Inf), (Inf, -Inf), (Inf, Inf) are NOT included. Only the combinations
  # of the middle thresholds (tau_1 to tau_(last-1)).
  # thres.var1.of.pair and thres.var2.of.pair give the first and the second
  # argument, respectively, in functions pbivnorm and lav_dbinorm
  thres_var1_of_pair <- vector("list", no_pairs)
  thres_var2_of_pair <- vector("list", no_pairs)

  # Extending the rho.vector accordingly so that it will be the the third
  # argument in pbivnorm and lav_dbinorm functions. It is of same length as
  # thres.var1.of.pair and thres.var2.of.pair.
  rho_vector <- vector("list", no_pairs)

  # thres.var1.for.dnorm.in.der.pi.to.tau.xi and
  # thres.var2.for.dnorm.in.der.pi.to.tau.xj give the thresholds of almost
  # all variables appropriately duplicated so that the vectors can be used
  # as input in dnorm() to compute der.pi.xixj.to.tau.xi and
  # der.pi.xixj.to.tau.xj.
  # thres.var1.for.dnorm.in.der.pi.to.tau.xi does not contain the thresholds of
  # the last variable and thres.var2.for.dnorm.in.der.pi.to.tau.xj those of
  # the first variable
  thres_var1_for_dnorm_in_der_pi_to_tau_xi <- vector("list", no_pairs) # nolint
  thres_var2_for_dnorm_in_der_pi_to_tau_xj <- vector("list", no_pairs) # nolint

  for (i in 1:no_pairs) {
    single_thres_var1_of_pair <-
       all_thres[index_var_of_thres == index_pairs[1, i]]
    single_thres_var2_of_pair <-
       all_thres[index_var_of_thres == index_pairs[2, i]]
    # remember that the first (-Inf) and last (Inf) thresholds are not included
    # so no.thres.var1.of.pair is equal to number of categories of var1 minus 1
    # similarly for no.thres.var2.of.pair
    no_thres_var1_of_pair <- no_thres_of_each_var[index_pairs[1, i]]
    no_thres_var2_of_pair <- no_thres_of_each_var[index_pairs[2, i]]

    thres_var1_of_pair[[i]] <- rep(single_thres_var1_of_pair,
      times = no_thres_var2_of_pair
    )
    thres_var2_of_pair[[i]] <- rep(single_thres_var2_of_pair,
      each = no_thres_var1_of_pair
    )
    rho_vector[[i]] <- rep(rho_xixj[i], length(thres_var1_of_pair[[i]]))

    thres_var1_for_dnorm_in_der_pi_to_tau_xi[[i]] <-   # nolint
      rep(single_thres_var1_of_pair, times = (no_thres_var2_of_pair + 1))
    thres_var2_for_dnorm_in_der_pi_to_tau_xj[[i]] <-   # nolint
      rep(single_thres_var2_of_pair, each = (no_thres_var1_of_pair + 1))
  }

  thres_var1_of_pair <- unlist(thres_var1_of_pair)
  thres_var2_of_pair <- unlist(thres_var2_of_pair)
  rho_vector <- unlist(rho_vector)
  thres_var1_for_dnorm_in_der_pi_to_tau_xi <-          # nolint
    unlist(thres_var1_for_dnorm_in_der_pi_to_tau_xi)
  thres_var2_for_dnorm_in_der_pi_to_tau_xj <-          # nolint
    unlist(thres_var2_for_dnorm_in_der_pi_to_tau_xj)

  # thres.var2.for.last.cat.var1 and thres.var1.for.last.cat.var2 are needed
  # for the computation of expected probabilities. In the computation of
  # \Phi_2(tau1, tau2; rho) when either tau1 or tau2 are Inf then it is enought
  # to compute pnorm() with the non-infinite tau as an argument
  # In particular when the first variable of the pair has tau_last= Inf
  # and the second a non-infite threshold we compute
  # pnorm(thres.var2.for.last.cat.var1). Similarly, when the second variable of
  # the pair has tau_last=Inf and the first a non-infite threshold we compute
  # pnorm(thres.var1.for.last.cat.var2).
  thres_var2_for_last_cat_var1 <- vector("list", (no_x - 1))
  thres_var1_for_last_cat_var2 <- vector("list", (no_x - 1))
  for (i in 1:(no_x - 1)) {
    thres_var2_for_last_cat_var1[[i]] <-
      c(all_thres[index_var_of_thres %in% (i + 1):no_x])
    thres_var1_for_last_cat_var2[[i]] <- rep(all_thres[index_var_of_thres == i],
      times = (no_x - i)
    )
  }
  thres_var2_for_last_cat_var1 <- unlist(thres_var2_for_last_cat_var1)
  thres_var1_for_last_cat_var2 <- unlist(thres_var1_for_last_cat_var2)


  list(
    thres.var1.of.pair = thres_var1_of_pair, # these 3 of same length
    thres.var2.of.pair = thres_var2_of_pair,
    rho.vector = rho_vector,

    # the following of length dependning on the number of categories
    thres.var1.for.dnorm.in.der.pi.to.tau.xi =
      thres_var1_for_dnorm_in_der_pi_to_tau_xi,
    thres.var2.for.dnorm.in.der.pi.to.tau.xj =
      thres_var2_for_dnorm_in_der_pi_to_tau_xj,
    thres.var2.for.last.cat.var1 = thres_var2_for_last_cat_var1,
    thres.var1.for.last.cat.var2 = thres_var1_for_last_cat_var2
  )
}
#########################################################



#########################################################

# The function  lav_pml_expprob_vec
# input: ind.vec - the output of function lav_pml_longvec_ind
#        th.rho.vec - the output of function lav_pml_longvec_th_rho
# output: it gives the elements of pairwiseTablesExpected()$pi.tables
# table-wise and column-wise within each table. In other words if
# pi^xixj_ab is the expected probability for the pair of variables xi-xj
# and categories a and b, then index a runs the fastest of all, followed by b,
# then by j, and lastly by i.

lav_pml_expprob_vec <- function(ind_vec, th_rho_vec) {
  prob_vec <- rep(NA, length(ind_vec$index.thres.var1.of.pair))

  prob_vec[ind_vec$index.thres.var1.of.pair == 0 |
    ind_vec$index.thres.var2.of.pair == 0] <- 0

  prob_vec[ind_vec$last.thres.var1.of.pair &
    ind_vec$last.thres.var2.of.pair] <- 1

  prob_vec[ind_vec$last.thres.var1.of.pair &
    ind_vec$index.thres.var2.of.pair != 0 &
    !ind_vec$last.thres.var2.of.pair] <-
    pnorm(th_rho_vec$thres.var2.for.last.cat.var1)

  prob_vec[ind_vec$last.thres.var2.of.pair &
    ind_vec$index.thres.var1.of.pair != 0 &
    !ind_vec$last.thres.var1.of.pair] <-
    pnorm(th_rho_vec$thres.var1.for.last.cat.var2)

  prob_vec[is.na(prob_vec)] <- pbivnorm(
    th_rho_vec$thres.var1.of.pair,
    th_rho_vec$thres.var2.of.pair,
    th_rho_vec$rho.vector
  )

  cum_term1 <- prob_vec[ind_vec$index.thres.var1.of.pair != 0 &
    ind_vec$index.thres.var2.of.pair != 0]

  cum_term2 <- prob_vec[ind_vec$index.thres.var1.of.pair != 0 &
    !ind_vec$last.thres.var2.of.pair]

  cum_term3 <- prob_vec[ind_vec$index.thres.var2.of.pair != 0 &
    !ind_vec$last.thres.var1.of.pair]

  cum_term4 <- prob_vec[!ind_vec$last.thres.var1.of.pair &
    !ind_vec$last.thres.var2.of.pair]

  pi0 <- cum_term1 - cum_term2 - cum_term3 + cum_term4

  # added by YR 11 nov 2012 to avoid Nan/-Inf
  # log(.Machine$double.eps) = -36.04365
  # all elements should be strictly positive
  pi0[pi0 < .Machine$double.eps] <- .Machine$double.eps

  pi0
}


# lav_pml_dl_drho
# input: ind.vec - the output of function lav_pml_longvec_ind
#        th.rho.vec - the output of function lav_pml_longvec_th_rho
#        n.xixj - a vector with the observed frequency for every combination
#                 of categories and every pair. The frequencies are given in
#                 the same order as the expected probabilities in the output of
#                 lav_pml_expprob_vec output
#        pi.xixj - the output of lav_pml_expprob_vec function
#        no.x    - the number of ordinal variables
# output: the vector of der.L.to.rho, each element corresponds to
#         der.Lxixj.to.rho.xixj where j runs faster than i

lav_pml_dl_drho <- function(ind_vec, th_rho_vec, n_xixj, pi_xixj, no_x) {
  prob_vec <- rep(NA, length(ind_vec$index.thres.var1.of.pair))

  prob_vec[ind_vec$index.thres.var1.of.pair == 0 |
    ind_vec$index.thres.var2.of.pair == 0 |
    ind_vec$last.thres.var1.of.pair |
    ind_vec$last.thres.var2.of.pair] <- 0

  prob_vec[is.na(prob_vec)] <- lav_dbinorm(th_rho_vec$thres.var1.of.pair,
    th_rho_vec$thres.var2.of.pair,
    rho = th_rho_vec$rho.vector
  )

  den_term1 <- prob_vec[ind_vec$index.thres.var1.of.pair != 0 &
    ind_vec$index.thres.var2.of.pair != 0]

  den_term2 <- prob_vec[ind_vec$index.thres.var1.of.pair != 0 &
    !ind_vec$last.thres.var2.of.pair]

  den_term3 <- prob_vec[ind_vec$index.thres.var2.of.pair != 0 &
    !ind_vec$last.thres.var1.of.pair]

  den_term4 <- prob_vec[!ind_vec$last.thres.var1.of.pair &
    !ind_vec$last.thres.var2.of.pair]

  der_pi_xixj_to_rho_xixj <- den_term1 - den_term2 - den_term3 + den_term4
  prod_terms <- (n_xixj / pi_xixj) * der_pi_xixj_to_rho_xixj

  # to get der.Lxixj.to.rho.xixj we should all the elements of
  # der.pi.xixj.to.rho.xixj which correspond to the pair xi-xj, to do so:
  xnew <- lapply(
    ind_vec[c("index.pairs.extended")],
    function(y) {
      y[ind_vec$index.thres.var1.of.pair != 0 &
        ind_vec$index.thres.var2.of.pair != 0]
    }
  )
  # der.L.to.rho is:
  tapply(prod_terms, xnew$index.pairs.extended, sum)
}
###########################################################################


# lav_pml_dl_dtau
# input: ind.vec - the output of function lav_pml_longvec_ind
#        th.rho.vec - the output of function lav_pml_longvec_th_rho
#        n.xixj - a vector with the observed frequency for every combination
#                 of categories and every pair. The frequencies are given in
#                 the same order as the expected probabilities in the output of
#                 lav_pml_expprob_vec output
#        pi.xixj - the output of lav_pml_expprob_vec function
# output: the vector of der.L.to.tau where the elements are ordered as follows:
#         the thresholds of each variable with respect to ascending order of
#         the variable index (i.e. thres_var1, thres_var2, etc.) and within
#         each variable the thresholds in ascending order.


lav_pml_dl_dtau <- function(ind_vec, th_rho_vec, n_xixj, pi_xixj, no_x = 0L) {
  # to compute der.pi.xixj.to.tau.xi
  xi <- lapply(
    ind_vec[c(
      "index.thres.var2.of.pair",
      "last.thres.var2.of.pair"
    )],
    function(y) {
      y[!(ind_vec$index.thres.var1.of.pair == 0 |
        ind_vec$last.thres.var1.of.pair)]
    }
  )

  cum_prob_vec <- rep(NA, length(xi$index.thres.var2.of.pair))
  cum_prob_vec[xi$index.thres.var2.of.pair == 0] <- 0
  cum_prob_vec[xi$last.thres.var2.of.pair] <- 1
  denom <- sqrt(1 - (th_rho_vec$rho.vector * th_rho_vec$rho.vector))
  cum_prob_vec[is.na(cum_prob_vec)] <-
    pnorm((th_rho_vec$thres.var2.of.pair -
      th_rho_vec$rho.vector * th_rho_vec$thres.var1.of.pair) /
      denom)
  den_prob_vec <- dnorm(th_rho_vec$thres.var1.for.dnorm.in.der.pi.to.tau.xi)
  der_pi_xixj_to_tau_xi <- den_prob_vec *
    (cum_prob_vec[xi$index.thres.var2.of.pair != 0] -
      cum_prob_vec[!xi$last.thres.var2.of.pair])

  # to compute der.pi.xixj.to.tau.xj
  xj <- lapply(
    ind_vec[c(
      "index.thres.var1.of.pair",
      "last.thres.var1.of.pair"
    )],
    function(y) {
      y[!(ind_vec$index.thres.var2.of.pair == 0 |
        ind_vec$last.thres.var2.of.pair)]
    }
  )

  cum_prob_vec <- rep(NA, length(xj$index.thres.var1.of.pair))
  cum_prob_vec[xj$index.thres.var1.of.pair == 0] <- 0
  cum_prob_vec[xj$last.thres.var1.of.pair] <- 1
  denom <- sqrt(1 - (th_rho_vec$rho.vector * th_rho_vec$rho.vector))
  cum_prob_vec[is.na(cum_prob_vec)] <-
    pnorm((th_rho_vec$thres.var1.of.pair -
      th_rho_vec$rho.vector * th_rho_vec$thres.var2.of.pair) /
      denom)
  den_prob_vec <- dnorm(th_rho_vec$thres.var2.for.dnorm.in.der.pi.to.tau.xj)
  der_pi_xixj_to_tau_xj <- den_prob_vec *
    (cum_prob_vec[xj$index.thres.var1.of.pair != 0] -
      cum_prob_vec[!xj$last.thres.var1.of.pair])

  # to compute der.Lxixj.tau.xi and der.Lxixj.tau.xi
  n_over_pi <- n_xixj / pi_xixj
  # get the appropriate differences of n.over.pi for der.Lxixj.to.tau.xi  and
  # der.Lxixj.to.tau.xj
  x3a <- lapply(ind_vec, function(y) {
    y[!(ind_vec$index.thres.var1.of.pair == 0 |
      ind_vec$index.thres.var2.of.pair == 0)]
  })
  diff_n_over_pi_to_xi <- n_over_pi[!x3a$last.thres.var1.of.pair] -
    n_over_pi[x3a$index.thres.var1.of.pair != 1]
  diff_n_over_pi_to_xj <- n_over_pi[!x3a$last.thres.var2.of.pair] -
    n_over_pi[x3a$index.thres.var2.of.pair != 1]
  # terms.der.Lxixj.to.tau.xi  and terms.der.Lxixj.to.tau.xj
  terms_der_lxixj_to_tau_xi <- diff_n_over_pi_to_xi * der_pi_xixj_to_tau_xi
  terms_der_lxixj_to_tau_xj <- diff_n_over_pi_to_xj * der_pi_xixj_to_tau_xj
  # to add appropriately elements of terms.der.Lxixj.to.tau.xi
  x3b <- lapply(
    ind_vec[c(
      "index.pairs.extended",
      "index.thres.var1.of.pair"
    )],
    function(y) {
      y[!(ind_vec$index.thres.var1.of.pair == 0 |
        ind_vec$last.thres.var1.of.pair |
        ind_vec$index.thres.var2.of.pair == 0)]
    }
  )
  # to add appropriately elements of terms.der.Lxixj.to.tau.xj
  x4b <- lapply(
    ind_vec[c(
      "index.pairs.extended",
      "index.thres.var2.of.pair"
    )],
    function(y) {
      y[!(ind_vec$index.thres.var2.of.pair == 0 |
        ind_vec$last.thres.var2.of.pair |
        ind_vec$index.thres.var1.of.pair == 0)]
    }
  )
  ind_pairs <- utils::combn(no_x, 2)
  # der.Lxixj.to.tau.xi is a matrix, nrow=no.pairs, ncol=max(no.of.free.thres)
  # thus, there are NA's, similarly for der.Lxixj.to.tau.xj
  der_lxixj_to_tau_xi <- tapply(
    terms_der_lxixj_to_tau_xi,
    list(
      x3b$index.pairs.extended,
      x3b$index.thres.var1.of.pair
    ),
    sum
  )
  der_lxixj_to_tau_xj <- tapply(
    terms_der_lxixj_to_tau_xj,
    list(
      x4b$index.pairs.extended,
      x4b$index.thres.var2.of.pair
    ),
    sum
  )

  # to add appropriately the terms of der.Lxixj.to.tau.xi and
  # der.Lxixj.to.tau.xj
  split_der_lxixj_to_tau_xi <- split(
    as.data.frame(der_lxixj_to_tau_xi),
    ind_pairs[1, ]
  )
  sums_der_lxixj_to_tau_xi <- lapply(
    split_der_lxixj_to_tau_xi,
    function(x) {
      y <- apply(x, 2, sum)
      y[!is.na(y)]
    }
  )
  # Note: NA exist in the case where the ordinal variables have different
  # number of response categories
  split_der_lxixj_to_tau_xj <- split(
    as.data.frame(der_lxixj_to_tau_xj),
    ind_pairs[2, ]
  )
  sums_der_lxixj_to_tau_xj <- lapply(
    split_der_lxixj_to_tau_xj,
    function(x) {
      y <- apply(x, 2, sum)
      y[!is.na(y)]
    }
  )
  # to get der.L.to.tau
  c(
    sums_der_lxixj_to_tau_xi[[1]],
    c(unlist(sums_der_lxixj_to_tau_xi[2:(no_x - 1)]) +
      unlist(sums_der_lxixj_to_tau_xj[1:(no_x - 2)])),
    sums_der_lxixj_to_tau_xj[[no_x - 1]]
  )
}
