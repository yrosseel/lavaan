# fitting function for standard ML
lav_model_objective_ml <- function(sigma_hat = NULL, mu_hat = NULL,
                         data_cov = NULL, data_mean = NULL,
                         data_cov_log_det = NULL,
                         meanstructure = FALSE) {
  # FIXME: WHAT IS THE BEST THING TO DO HERE??
  # CURRENTLY: return Inf  (at least for nlminb, this works well)
  if (!attr(sigma_hat, "po")) {
    return(Inf)
  }


  sigma_hat_inv <- attr(sigma_hat, "inv")
  sigma_hat_log_det <- attr(sigma_hat, "log.det")
  nvar <- ncol(sigma_hat)

  if (!meanstructure) {
    fx <- (sigma_hat_log_det + sum(data_cov * sigma_hat_inv) -
      data_cov_log_det - nvar)
  } else {
    w_tilde <- data_cov + tcrossprod(data_mean - mu_hat)
    fx <- (sigma_hat_log_det + sum(w_tilde * sigma_hat_inv) -
      data_cov_log_det - nvar)
  }

  # no negative values
  if (is.finite(fx) && fx < 0.0) fx <- 0.0

  fx
}

# fitting function for standard ML
lav_model_objective_ml_res <- function(sigma_hat = NULL, mu_hat = NULL,
                                       pi0 = NULL, res_cov = NULL,
                                       res_int = NULL, res_slopes = NULL,
                                       res_cov_log_det = NULL,
                                       cov_x = NULL, mean_x = NULL) {
  if (!attr(sigma_hat, "po")) {
    return(Inf)
  }

  # augmented mean.x + cov.x matrix
  c3 <- rbind(
    c(1, mean_x),
    cbind(mean_x, cov_x + tcrossprod(mean_x))
  )

  sigma_hat_inv <- attr(sigma_hat, "inv")
  sigma_hat_log_det <- attr(sigma_hat, "log.det")
  nvar <- ncol(sigma_hat)

  # sigma
  objective_sigma <- (sigma_hat_log_det + sum(res_cov * sigma_hat_inv) -
    res_cov_log_det - nvar)
  # beta
  obs <- t(cbind(res_int, res_slopes))
  est <- t(cbind(mu_hat, pi0))
  diff_1 <- obs - est
  objective_beta <- sum(sigma_hat_inv * crossprod(diff_1, c3) %*% diff_1)

  fx <- objective_sigma + objective_beta

  # no negative values
  if (is.finite(fx) && fx < 0.0) fx <- 0.0

  fx
}


# fitting function for restricted ML
lav_model_objective_reml <- function(sigma_hat = NULL, mu_hat = NULL,
                           data_cov = NULL, data_mean = NULL,
                           data_cov_log_det = NULL,
                           meanstructure = FALSE,
                           group = 1L, lavmodel = NULL,
                           lavsamplestats = NULL, lavdata = NULL) {
  if (!attr(sigma_hat, "po")) {
    return(Inf)
  }

  sigma_hat_inv <- attr(sigma_hat, "inv")
  sigma_hat_log_det <- attr(sigma_hat, "log.det")
  nvar <- ncol(sigma_hat)

  if (!meanstructure) {
    fx <- (sigma_hat_log_det + sum(data_cov * sigma_hat_inv) -
      data_cov_log_det - nvar)
  } else {
    w_tilde <- data_cov + tcrossprod(data_mean - mu_hat)
    fx <- (sigma_hat_log_det + sum(w_tilde * sigma_hat_inv) -
      data_cov_log_det - nvar)
  }

  lambda_idx <- which(names(lavmodel@GLIST) == "lambda")
  mm_lambda <- lavmodel@GLIST[[lambda_idx[group]]]
  data_cov_inv <- lavsamplestats@icov[[group]]
  reml_h0 <- log(det(t(mm_lambda) %*% sigma_hat_inv %*% mm_lambda))
  reml_h1 <- log(det(t(mm_lambda) %*% data_cov_inv %*% mm_lambda))
  nobs <- lavsamplestats@nobs[[group]]

  # fx <- (Sigma.hat.log.det + tmp - data.cov.log.det - nvar) +
  #                        1/Ng * (reml.h0  - reml.h1)
  fx <- fx + (1 / nobs * (reml_h0 - reml_h1))

  # no negative values
  if (is.finite(fx) && fx < 0.0) fx <- 0.0

  fx
}

# 'classic' fitting function for GLS
# used again since 0.6-10 (we used the much slower lav_model_objective_wls before)
lav_model_objective_gls <- function(sigma_hat = NULL,
                                    mu_hat = NULL,
                                    data_cov = NULL,
                                    data_cov_inv = NULL,
                                    data_mean = NULL,
                                    meanstructure = FALSE,
                                    correlation = FALSE) {
  tmp <- data_cov_inv %*% (data_cov - sigma_hat)
  # tmp is not perfectly symmetric, so we use t(tmp) on the next line
  # to obtain the same value as lav_model_objective_wls
  fx <- 0.5 * sum(tmp * t(tmp))

  if (correlation) {
    # Bentler & Savalei (2010) eq 1.31
    dd <- as.matrix(diag(tmp))
    tt <- diag(nrow(data_cov)) + data_cov * data_cov_inv
    fx <- fx - drop(t(dd) %*% solve(tt) %*% dd)
  }

  if (meanstructure) {
    tmp2 <- sum(data_cov_inv * tcrossprod(data_mean - mu_hat))
    fx <- fx + tmp2
  }

  # no negative values
  if (is.finite(fx) && fx < 0.0) fx <- 0.0

  fx
}

# general WLS estimator (Muthen, Appendix 4, eq 99 single group)
# full weight (WLS.V) matrix
lav_model_objective_wls <- function(wls_est = NULL,
                                 wls_obs = NULL, wls_v = NULL) {
  # diff <- as.matrix(WLS.obs - WLS.est)
  # fx <- as.numeric( t(diff) %*% WLS.V %*% diff )

  # since 0.5-17, we use crossprod twice
  diff <- wls_obs - wls_est
  fx <- as.numeric(crossprod(crossprod(wls_v, diff), diff))
  # todo alternative: using chol(WLS.V)

  # no negative values
  if (is.finite(fx) && fx < 0.0) fx <- 0.0

  fx
}

# diagonally weighted LS (DWLS)
lav_model_objective_dwls <- function(wls_est = NULL,
                            wls_obs = NULL, wls_vd = NULL) {
  diff <- wls_obs - wls_est
  fx <- sum(diff * diff * wls_vd)

  # no negative values
  if (is.finite(fx) && fx < 0.0) fx <- 0.0

  fx
}

# Full Information ML estimator (FIML) handling the missing values
lav_model_objective_fiml <- function(sigma_hat = NULL, mu_hat = NULL, yp = NULL,
                           h1 = NULL, n = NULL) {
  if (is.null(n)) {
    n <- sum(sapply(yp, "[[", "freq"))
  }

  # Note: we ignore x.idx (if any)
  fx <- lav_mvn_mi_loglik_samp(
    yp = yp,
    mu = mu_hat, sigma_1 = sigma_hat,
    log2pi = FALSE,
    minus_two = TRUE
  ) / n

  # ajust for h1
  if (!is.null(h1)) {
    fx <- fx - h1

    # no negative values
    if (is.finite(fx) && fx < 0.0) fx <- 0.0
  }

  fx
}

# pairwise maximum likelihood
# this is adapted from code written by Myrsini Katsikatsou
#
# some changes:
# - no distinction between x/y (ksi/eta)
# - 29/03/2016: adapt for exogenous covariates
# - 21/09/2016: added code for missing = doubly.robust (contributed by
#   Myrsini Katsikatsou)
# - HJ 18/10/2023: For sampling weights the lavcache$bifreq are weighted
lav_model_objective_pml <- function(sigma_hat = NULL, # model-based var/cov/cor
                          mu_hat = NULL, # model-based means
                          th = NULL, # model-based thresholds + means
                          pi0 = NULL, # slopes
                          th_idx = NULL, # threshold idx per variable
                          num_idx = NULL, # which variables are numeric
                          x = NULL, # raw data
                          exo = NULL, # eXo data
                          wt = NULL, # case weights
                          lavcache = NULL, # housekeeping stuff
                          missing = NULL) { # how to deal with missings?

  # YR 3 okt 2012
  # - the idea is to compute for each pair of variables, the model-based
  #   probability (or likelihood in mixed case) (that we observe the data
  #   for this pair under the model)
  # - if we have exogenous variables + conditional.x, do this for each case
  # - after taking logs, the sum over the cases gives the
  #   log probablity/likelihood for this pair
  # - the sum over all pairs gives the final PL based logl

  # first of all: check if all correlations are within [-1,1]
  # if not, return Inf; (at least with nlminb, this works well)

  # diagonal of Sigma.hat is not necessarily 1, even for categorical vars
  sigma_hat2 <- sigma_hat
  if (length(num_idx) > 0L) {
    diag(sigma_hat2)[-num_idx] <- 1
  } else {
    diag(sigma_hat2) <- 1
  }
  # all positive variances? (for continuous variables)
  if (any(diag(sigma_hat2) < 0)) {
    out <- +Inf
    attr(out, "logl") <- as.numeric(NA)
    return(out)
  }
  cor_hat <- cov2cor(sigma_hat2) # to get correlations (rho!)
  cors <- lav_mat_vech(cor_hat, diagonal = FALSE)

  if (length(cors) > 0L && (any(abs(cors) > 1) ||
    any(is.na(cors)))) {
    # question: what is the best approach here??
    out <- +Inf
    attr(out, "logl") <- as.numeric(NA)
    return(out)
  }

  nvar <- nrow(sigma_hat)
  if (is.null(exo)) {
    nexo <- 0L
  } else {
    nexo <- NCOL(exo)
  }
  pstar <- nvar * (nvar - 1) / 2
  ov_types <- rep("ordered", nvar)
  if (length(num_idx) > 0L) {
    ov_types[num_idx] <- "numeric"
  }


  ##### Three cases:
  ##### 1) all ordered, no exogenous (fast!)
  ##### 2) mixed ordered + continuous, no exogenous
  ##### 3) mixed ordered + continuous, exogenous (conditional.x = TRUE)




  ##### Case 1:
  #####  all ordered
  #####  no exogenous covariates
  #####
  if (all(ov_types == "ordered") && nexo == 0L) {
    # prepare for Myrsini's vectorization scheme
    long2 <- lav_pml_longvec_th_rho(
      no_x = nvar,
      all_thres = th,
      index_var_of_thres = th_idx,
      rho_xixj = cors
    )
    # get expected probability per table, per pair
    pairwise_pi <- lav_pml_expprob_vec(
      ind_vec = lavcache$long,
      th_rho_vec = long2
    )
    pairwise_pi_orig <- pairwise_pi # for doubly.robust

    # get frequency per table, per pair
    logl <- sum(lavcache$bifreq * log(pairwise_pi))

    # >>>>>>>> HJ/MK PML CODE >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    # FYI the bifreq are already weighted so this will work. Alternatively:
    if (!is.null(wt)) {
      logl <- sum(lavcache$sum_obs_weights_xixj_ab_vec * log(pairwise_pi))
    }

    # more convenient fit function
    prop <- lavcache$bifreq / lavcache$nobs
    freq <- lavcache$bifreq
    if (!is.null(wt)) {
      prop <- lavcache$sum_obs_weights_xixj_ab_vec / sum(wt)
      freq <- lavcache$sum_obs_weights_xixj_ab_vec
    }

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    # remove zero props # FIXME!!! or add 0.5???
    # zero.idx <- which(prop == 0.0)
    zero_idx <- which((prop == 0.0) | !is.finite(prop))
    if (length(zero_idx) > 0L) {
      freq <- freq[-zero_idx]
      prop <- prop[-zero_idx]
      pairwise_pi <- pairwise_pi[-zero_idx]
    }
    ## Fmin <- sum( prop*log(prop/pairwisePI) )
    fmin <- sum(freq * log(prop / pairwise_pi)) # to avoid 'N'

    if (missing == "available.cases" || missing == "doubly.robust") {
      uni_pi <- lav_pml_th_uni_prob(th = th, th_idx = th_idx)

      # shortcuts
      unifreq <- lavcache$unifreq
      uninobs <- lavcache$uninobs
      uniweights <- lavcache$uniweights

      logl <- logl + sum(uniweights * log(uni_pi))

      uniprop <- unifreq / uninobs

      # remove zero props
      # uni.zero.idx <- which(uniprop == 0.0)
      uni_zero_idx <- which((uniprop == 0.0) | !is.finite(uniprop))
      if (length(uni_zero_idx) > 0L) {
        uniprop <- uniprop[-uni_zero_idx]
        uni_pi <- uni_pi[-uni_zero_idx]
        uniweights <- uniweights[-uni_zero_idx]
      }

      fmin <- fmin + sum(uniweights * log(uniprop / uni_pi))
    }

    if (missing == "doubly.robust") {
      # COMPUTE THE SUM OF THE EXPECTED BIVARIATE CONDITIONAL LIKELIHOODS
      # SUM_{i,j} [ E_{Yi,Yj|y^o}(lnf(Yi,Yj))) ]

      # First compute the terms of the summand. Since the cells of
      # pairwiseProbGivObs are zero for the pairs of variables that at least
      # one of the variables is observed (hence not contributing to the summand)
      # there is no need to construct an index vector for summing appropriately
      # within each individual.
      log_pairwise_pi_orig <- log(pairwise_pi_orig)
      pairwise_prob_giv_obs <- lavcache$pairwiseProbGivObs
      tmp_prod <- t(t(pairwise_prob_giv_obs) * log_pairwise_pi_orig)

      sum_elnfij_casewise <- apply(tmp_prod, 1, sum)
      sum_elnfij <- sum(sum_elnfij_casewise)
      logl <- logl + sum_elnfij
      fmin <- fmin - sum_elnfij

      # COMPUTE THE THE SUM OF THE EXPECTED UNIVARIATE CONDITIONAL LIKELIHOODS
      # SUM_{i,j} [ E_{Yj|y^o}(lnf(Yj|yi))) ]

      # First compute the model-implied conditional univariate probabilities
      # p(y_i=a|y_j=b). Let ModProbY1Gy2 be the vector of these
      # probabilities. The order the probabilities
      # are listed in the vector ModProbY1Gy2 is as follows:
      # y1|y2, y1|y3, ..., y1|yp, y2|y1, y2|y3, ..., y2|yp,
      # ..., yp|y1, yp|y2, ..., yp|y(p-1). Within each pair of variables the
      # index "a" which represents the response category of variable yi runs
      # faster than "b" which represents the response category of the given
      # variable yj.
      # The computation of these probabilities are based on the model-implied
      # bivariate probabilities p(y_i=a,y_j=b). To do the appropriate summations
      # and divisions we need some index vectors to keep track of the index i,j,
      # a, and b, as well as the pair index. These index vectors should be
      # computed once and stored in lavcache. About where in the lavaan code
      # we will add the computations and how they will be done please see the
      # file "new objects in lavcache for DR-PL.r"

      idx_pairs <- lavcache$idx.pairs
      idx_cat_y2_split <- lavcache$idx.cat.y2.split
      idx_cat_y1_split <- lavcache$idx.cat.y1.split
      idx_y1 <- lavcache$idx.Y1
      idx_gy2 <- lavcache$idx.Gy2
      # idx_cat_y1 <- lavcache$idx.cat.Y1
      idx_cat_gy2 <- lavcache$idx.cat.Gy2
      id_uni_pr_giv_obs <- lavcache$id.uniPrGivObs
      # the latter keeps track which variable each column of the matrix
      # univariateProbGivObs refers to

      # For the function lav_pml_bivprob_unicondprob see the .r file
      # with the same name.
      mod_prob_y1gy2 <- lav_pml_bivprob_unicondprob(
        biv_prob = pairwise_pi_orig,
        nvar = nvar,
        idx_pairs = idx_pairs,
        idx_y1 = idx_y1,
        idx_gy2 = idx_gy2,
        idx_cat_y1_split = idx_cat_y1_split,
        idx_cat_y2_split = idx_cat_y2_split
      )

      log_mod_prob_y1gy2 <- log(mod_prob_y1gy2)

      # Let univariateProbGivObs be the matrix of the conditional univariate
      # probabilities Pr(y_i=a|y^o) that has been computed in advance and are
      # fed to the DR-PL function. The rows represent different individuals,
      # i.e. nrow=nobs, and the columns different probabilities. The columns
      # are listed as follows: a runs faster than i.

      # Note that the number of columns of univariateProbGivObs is not the
      # same with the length(log_ModProbY1Gy2), actually
      # ncol(univariateProbGivObs) < length(log_ModProbY1Gy2).
      # For this we use the following commands in order to multiply correctly.

      # Compute for each case the product  Pr(y_i=a|y^o) * log[ p(y_i=a|y_j=b) ]
      # i.e. univariateProbGivObs * log_ModProbY1Gy2
      univariate_prob_giv_obs <- lavcache$univariateProbGivObs
      nobs <- nrow(x)
      uniweights_casewise <- lavcache$uniweights.casewise
      id_cases_with_missing <- which(uniweights_casewise > 0)
      no_cases_with_missing <- length(id_cases_with_missing)
      no_obs_casewise <- nvar - uniweights_casewise
      idx_missing_var <- apply(x, 1, function(x) {
        which(is.na(x))
      })
      idx_observed_var <- lapply(idx_missing_var, function(x) {
        c(1:nvar)[-x]
      })
      idx_cat_observed_var <- sapply(1:nobs, function(i) {
        x[i, idx_observed_var[[i]]]
      })
      elnyi_givyjb_casewise <- sapply(1:no_cases_with_missing, function(i) {
        tmp_id_case <- id_cases_with_missing[i]
        tmp_no_mis <- uniweights_casewise[tmp_id_case]
        tmp_idx_mis <- idx_missing_var[[tmp_id_case]]
        tmp_idx_obs <- idx_observed_var[[tmp_id_case]]
        tmp_no_obs <- no_obs_casewise[tmp_id_case]
        tmp_idx_cat_obs <- idx_cat_observed_var[[tmp_id_case]]
        tmp_uni_prob_giv_obs_i <- univariate_prob_giv_obs[tmp_id_case, ]
        sapply(1:tmp_no_mis, function(k) {
          tmp_idx_mis_var <- tmp_idx_mis[k]
          tmp_uni_prob_giv_obs_ik <-
            tmp_uni_prob_giv_obs_i[id_uni_pr_giv_obs == tmp_idx_mis_var]
          tmp_log_mod_prob_y1gy2 <- sapply(1:tmp_no_obs, function(z) {
            log_mod_prob_y1gy2[idx_y1 == tmp_idx_mis_var &
              idx_gy2 == tmp_idx_obs[z] &
              idx_cat_gy2 == tmp_idx_cat_obs[z]]
          })
          sum(tmp_log_mod_prob_y1gy2 * tmp_uni_prob_giv_obs_ik)
        })
      })
      elnyi_givyjb <- sum(unlist(elnyi_givyjb_casewise))
      logl <- logl + elnyi_givyjb
      # for the Fmin function
      fmin <- fmin - elnyi_givyjb
    } # end of if (missing =="doubly.robust")


    ##### Case 2:
    #####  mixed ordered + numeric
    #####  no exogenous covariates
    #####
  } else if (nexo == 0L) {
    # mixed ordered/numeric variables, but no exogenous covariates
    # - no need to compute 'casewise' (log)likelihoods

    pstar_1 <- matrix(0, nvar, nvar) # utility matrix, to get indices
    pstar_1[lav_mat_vech_idx(nvar, diagonal = FALSE)] <- 1:pstar
    # n <- NROW(x)

    log_lik_pair <- numeric(pstar) # logl per pair (summed over cases)
    for (j in seq_len(nvar - 1L)) {
      for (i in (j + 1L):nvar) {
        pstar_idx <- pstar_1[i, j]
        if (ov_types[i] == "numeric" &&
          ov_types[j] == "numeric") {
          log_lik <- lav_mvn_loglik_data(
            y = x[, c(i, j)], wt = wt, mu = mu_hat[c(i, j)],
            sigma_1 = sigma_hat[c(i, j), c(i, j)], casewise = TRUE
          )
          log_lik_pair[pstar_idx] <- sum(log_lik, na.rm = TRUE)
        } else if (ov_types[i] == "numeric" &&
          ov_types[j] == "ordered") {
          # polyserial correlation
          log_lik <- lav_bvmix_lik(
            y1 = x[, i], y2 = x[, j],
            wt = wt,
            evar_y1 = sigma_hat[i, i],
            beta_y1 = mu_hat[i],
            th_y2 = th[th_idx == j],
            rho = cor_hat[i, j], take_log = TRUE
          )
          log_lik_pair[pstar_idx] <- sum(log_lik, na.rm = TRUE)
        } else if (ov_types[j] == "numeric" &&
          ov_types[i] == "ordered") {
          # polyserial correlation
          log_lik <- lav_bvmix_lik(
            y1 = x[, j], y2 = x[, i],
            wt = wt,
            evar_y1 = sigma_hat[j, j],
            beta_y1 = mu_hat[j],
            th_y2 = th[th_idx == i],
            rho = cor_hat[i, j], take_log = TRUE
          )
          log_lik_pair[pstar_idx] <- sum(log_lik, na.rm = TRUE)
        } else if (ov_types[i] == "ordered" &&
          ov_types[j] == "ordered") {
          # polychoric correlation
          pairwise_pi <- lav_bvord_noexo_pi(
            rho = cor_hat[i, j],
            th_y1 = th[th_idx == i],
            th_y2 = th[th_idx == j]
          )
          # avoid zeroes
          pairwise_pi[pairwise_pi < .Machine$double.eps] <-
            .Machine$double.eps
          # note: missing values are just not counted
          freq_1 <- lav_bvord_freq(x[, i], x[, j], wt = wt)
          log_lik_pair[pstar_idx] <- sum(freq_1 * log(pairwise_pi))
        }
      }
    } # all pairs

    na_idx <- which(is.na(log_lik_pair))
    if (length(na_idx) > 0L) {
      lav_msg_warn(gettext("some pairs produce NA values for logl:"),
        lav_msg_view(round(log_lik_pair, 3), "none")
      )
    }

    # sum over pairs
    logl <- sum(log_lik_pair)

    # Fmin
    fmin <- (-1) * logl


    ##### Case 3:
    #####  mixed ordered + numeric
    #####  exogenous covariates
    #####  (conditional.x = TRUE)
  } else {
    lik <- matrix(0, nrow(x), pstar) # likelihood per case, per pair
    pstar_1 <- matrix(0, nvar, nvar) # utility matrix, to get indices
    pstar_1[lav_mat_vech_idx(nvar, diagonal = FALSE)] <- 1:pstar
    # n <- NROW(x)

    for (j in seq_len(nvar - 1L)) {
      for (i in (j + 1L):nvar) {
        pstar_idx <- pstar_1[i, j]
        # cat("pstar.idx =", pstar.idx, "i = ", i, " j = ", j, "\n")
        if (ov_types[i] == "numeric" &&
            ov_types[j] == "numeric") {
          # ordinary pearson correlation
          lik[, pstar_idx] <-
            lav_bvreg_lik(
              y1 = x[, i], y2 = x[, j], exo = exo,
              wt = wt,
              evar_y1 = sigma_hat[i, i],
              beta_y1 = c(mu_hat[i], pi0[i, ]),
              evar_y2 = sigma_hat[j, j],
              beta_y2 = c(mu_hat[j], pi0[j, ]),
              rho = cor_hat[i, j]
            )
        } else if (ov_types[i] == "numeric" &&
                   ov_types[j] == "ordered") {
          # polyserial correlation
          ### FIXME: th.y2 should go into ps_lik!!!
          lik[, pstar_idx] <-
            lav_bvmix_lik(
              y1 = x[, i], y2 = x[, j], exo = exo,
              wt = wt,
              evar_y1 = sigma_hat[i, i],
              beta_y1 = c(mu_hat[i], pi0[i, ]),
              th_y2 = th[th_idx == j],
              sl_y2 = pi0[j, ],
              rho = cor_hat[i, j]
            )
        } else if (ov_types[j] == "numeric" &&
                   ov_types[i] == "ordered") {
          # polyserial correlation
          ### FIXME: th.y1 should go into ps_lik!!!
          lik[, pstar_idx] <-
            lav_bvmix_lik(
              y1 = x[, j], y2 = x[, i], exo = exo,
              wt = wt,
              evar_y1 = sigma_hat[j, j],
              beta_y1 = c(mu_hat[j], pi0[j, ]),
              th_y2 = th[th_idx == i],
              sl_y2 = pi0[i, ],
              rho = cor_hat[i, j]
            )
        } else if (ov_types[i] == "ordered" &&
                   ov_types[j] == "ordered") {
          lik[, pstar_idx] <-
            lav_pml_bi_lik_x(
              y1 = x[, i],
              y2 = x[, j],
              rho = sigma_hat[i, j],
              th_y1 = th[th_idx == i],
              th_y2 = th[th_idx == j],
              exo = exo,
              pi_y1 = pi0[i, ],
              pi_y2 = pi0[j, ],
              missing_ind = missing
            )
        }
      }
    } # all pairs

    # check for zero likelihoods/probabilities
    # FIXME: or should we replace them with a tiny number?
    if (any(lik == 0.0, na.rm = TRUE)) {
      out <- +Inf
      attr(out, "logl") <- as.numeric(NA)
      return(out)
    }

    # loglikelihood
    log_lik_cases <- log(lik)

    # sum over cases
    log_lik_pairs <- colSums(log_lik_cases, na.rm = TRUE)

    # sum over pairs
    logl <- logl_pairs <- sum(log_lik_pairs)

    if (missing == "available.cases" && all(ov_types == "ordered") &&
      nexo != 0L) {
      uni_lik <- matrix(0, nrow(x), ncol(x))
      for (i in seq_len(nvar)) {
        uni_lik[, i] <- lav_pml_uni_lik(
          y1 = x[, i],
          th_y1 = th[th_idx == i],
          exo = exo,
          pi_y1 = pi0[i, ]
        )
      }

      if (any(uni_lik == 0.0, na.rm = TRUE)) {
        out <- +Inf
        attr(out, "logl") <- as.numeric(NA)
        return(out)
      }

      uni_log_lik_cases <- log(uni_lik) * lavcache$uniweights.casewise

      # sum over cases
      uni_log_lik_varwise <- colSums(uni_log_lik_cases)

      # sum over variables
      uni_log_lik <- sum(uni_log_lik_varwise)

      # add with the pairwise part of LogLik
      logl <- logl_pairs + uni_log_lik
    }

    # we minimise
    fmin <- (-1) * logl
  }


  # here, we should have two quantities: logl and Fmin

  # function value as returned to the minimizer
  fx <- fmin

  # attach 'loglikelihood'
  attr(fx, "logl") <- logl

  fx
}

# full information maximum likelihood
# underlying multivariate normal approach (see Joreskog & Moustaki, 2001)
#
lav_model_objective_fml <- function(sigma_hat = NULL, # model-based var/cov/cor
                          th = NULL, # model-based thresholds + means
                          th_idx = NULL, # threshold idx per variable
                          num_idx = NULL, # which variables are numeric
                          x = NULL, # raw data
                          lavcache = NULL) { # patterns

  # YR 27 aug 2013
  # just for fun, and to compare with PML for small models

  # first of all: check if all correlations are within [-1,1]
  # if not, return Inf; (at least with nlminb, this works well)
  cors <- sigma_hat[lower.tri(sigma_hat)]

  if (any(abs(cors) > 1)) {
    return(+Inf)
  }

  nvar <- nrow(sigma_hat)
  # pstar <- nvar * (nvar - 1) / 2
  ov_types <- rep("ordered", nvar)
  if (length(num_idx) > 0L) ov_types[num_idx] <- "numeric"
  mean_1 <- rep(0, nvar)

  # shortcut for all ordered - per pattern
  if (all(ov_types == "ordered")) {
    pat <- lavcache$pat
    npatterns <- nrow(pat)
    freq <- as.numeric(rownames(pat))
    pi0 <- numeric(npatterns)
    th_var <- lapply(1:nvar, function(x) c(-Inf, th[th_idx == x], +Inf))
    # FIXME!!! ok to set diagonal to 1.0?
    diag(sigma_hat) <- 1.0
    for (r in 1:npatterns) {
      # compute probability for each pattern
      lower <- sapply(1:nvar, function(x) th_var[[x]][pat[r, x]])
      upper <- sapply(1:nvar, function(x) th_var[[x]][pat[r, x] + 1L])


      # how accurate must we be here???
      pi0[r] <- sadmvn(lower, upper,
        mean = mean_1, varcov = sigma_hat,
        maxpts = 10000 * nvar, abseps = 1e-07
      )
    }
    # sum (log)likelihood over all patterns
    # LogLik <- sum(log(PI) * freq)

    # more convenient fit function
    prop <- freq / sum(freq)
    # remove zero props # FIXME!!! or add 0.5???
    zero_idx <- which(prop == 0.0)
    if (length(zero_idx) > 0L) {
      prop <- prop[-zero_idx]
      pi0 <- pi0[-zero_idx]
    }
    fmin <- sum(prop * log(prop / pi0))
  } else { # case-wise
    pi0 <- numeric(nobs)
    for (i in 1:nobs) {
      # compute probability for each case
      pi0[i] <- lav_msg_stop(gettext("not implemented"))
    }
    # sum (log)likelihood over all observations
    # log_lik <- sum(log(pi0))
    lav_msg_stop(gettext("not implemented"))
  }

  # function value as returned to the minimizer
  # fx <- -1 * LogLik
  fx <- fmin

  fx
}

lav_model_objective_mml <- function(lavmodel = NULL,
                          mm_theta = NULL,
                          th = NULL,
                          glist = NULL,
                          group = 1L,
                          lavdata = NULL,
                          sample_mean = NULL,
                          sample_mean_x = NULL,
                          lavcache = NULL) {
  # compute case-wise likelihoods
  lik <- lav_model_lik_mml(
    lavmodel = lavmodel, mm_theta = mm_theta, th = th,
    glist = glist, group = group, lavdata = lavdata,
    sample_mean = sample_mean, sample_mean_x = sample_mean_x,
    lavcache = lavcache
  )

  # log + sum over observations
  logl <- sum(log(lik))

  # function value as returned to the minimizer
  fx <- -logl

  fx
}

lav_model_objective_2l <- function(lavmodel = NULL,
                         glist = NULL,
                         y1 = NULL, # only for missing
                         lp = NULL,
                         mp = NULL,
                         lavsamplestats = NULL,
                         lavcache_group = NULL, # only for random slopes
                         group = 1L) {
  # random slopes? (rv() modifier) -> use the random-slope kernel
  if (length(lavmodel@rv.ov) > 0L || length(lavmodel@rv.lv) > 0L) {
    rs <- lavcache_group$rs
    if (is.null(rs)) {
      lav_msg_fixme(
        "no rs element found in lavcache; this should not happen")
    }
    loglik <- lav_mvn_cl_rs_m2ll(
      lavmodel = lavmodel, glist = glist, rs = rs,
      log2pi = FALSE, minus_two = TRUE
    )
    objective <- as.numeric(loglik) / (lavsamplestats@ntotal * 2)
    return(objective)
  }

  # compute model-implied statistics for all blocks
  implied <- lav_model_implied(lavmodel, glist = glist)

  # here, we assume only 2!!! levels, at [[1]] and [[2]]
  if (lavmodel@conditional.x) {
    res_sigma_w <- implied$res.cov[[(group - 1) * 2 + 1]]
    res_int_w   <- implied$res.int[[(group - 1) * 2 + 1]]
    res_pi_w    <- implied$res.slopes[[(group - 1) * 2 + 1]]

    res_sigma_b <- implied$res.cov[[(group - 1) * 2 + 2]]
    res_int_b   <- implied$res.int[[(group - 1) * 2 + 2]]
    res_pi_b    <- implied$res.slopes[[(group - 1) * 2 + 2]]
  } else {
    sigma_w <- implied$cov[[(group - 1) * 2 + 1]]
    mu_w    <- implied$mean[[(group - 1) * 2 + 1]]
    sigma_b <- implied$cov[[(group - 1) * 2 + 2]]
    mu_b    <- implied$mean[[(group - 1) * 2 + 2]]
  }

  if (lavsamplestats@missing.flag) {
    if (lavmodel@conditional.x) {
      lav_msg_stop(gettext("multilevel + conditional.x is not ready yet for
                           fiml; rerun with conditional.x = FALSE"))
    }
    # SIGMA.B <- Sigma.B[Lp$both.idx[[2]], Lp$both.idx[[2]], drop = FALSE]
    # if(any(diag(SIGMA.B) < 0)) {
    #    return(+Inf)
    # }
    # COR.B <- cov2cor(SIGMA.B)
    # if(any(abs(lav_mat_vech(COR.B, diagonal = FALSE)) > 1)) {
    #   return(+Inf)
    # }

    y2 <- lavsamplestats@YLp[[group]][[2]]$Y2
    # yp <- lavsamplestats@missing[[group]]
    loglik <- lav_mvn_cl_mi_loglik_samp_2l(
      y1 = y1,
      y2 = y2, lp = lp, mp = mp,
      mu_w = mu_w, sigma_w = sigma_w,
      mu_b = mu_b, sigma_b = sigma_b,
      log2pi = FALSE, minus_two = TRUE
    )
  } else {
    ylp <- lavsamplestats@YLp[[group]]
    if (lavmodel@conditional.x) {
      loglik <- lav_mvreg_cl_loglik_samp_2l(
        ylp = ylp, lp = lp,
        res_sigma_w = res_sigma_w,
        res_int_w = res_int_w, res_pi_w = res_pi_w,
        res_sigma_b = res_sigma_b,
        res_int_b = res_int_b, res_pi_b = res_pi_b,
        log2pi = FALSE, minus_two = TRUE
      )
    } else {
      loglik <- lav_mvn_cl_loglik_samp_2l(
        ylp = ylp, lp = lp,
        mu_w = mu_w, sigma_w = sigma_w,
        mu_b = mu_b, sigma_b = sigma_b,
        log2pi = FALSE, minus_two = TRUE
      )
    }
  }

  # minimize (we already did -2*)
  objective <- 1 * loglik

  # divide by (N*2)
  objective <- objective / (lavsamplestats@ntotal * 2)

  # should be strictly positive
  # if(objective < 0) {
  #   objective <- +Inf
  # }

  objective
}
