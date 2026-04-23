# This code was contributed by Myrsini Katsikatsou (LSE) -- September 2016
#
# lav_pml_bivprob_unicondprob()
# lav_pml_object_uni_pairwise_prob()
# lav_pml_th_rho_generalised()
# pairwiseExpProbVec_GivenObs_UncMod()

lav_pml_bivprob_unicondprob <- function(biv_prob, nvar,
                                                 idx_pairs,
                                                 idx_y1,
                                                 idx_gy2,
                                                 idx_cat_y1_split,
                                                 idx_cat_y2_split) {
  biv_prob_split <- split(biv_prob, idx_pairs)
  lngth <- 2 * length(biv_prob)
  idx_vec_el <- 1:lngth
  prob_y1gy2 <- rep(NA, lngth)
  no_pairs <- nvar * (nvar - 1) / 2
  idx2_pairs <- combn(nvar, 2)

  for (k in 1:no_pairs) {
    y2sums <- tapply(biv_prob_split[[k]], idx_cat_y2_split[[k]], sum)
    y2sums_mult <- y2sums[idx_cat_y2_split[[k]]]
    y1gy2 <- biv_prob_split[[k]] / y2sums_mult
    tmp_idx_vec_el <- idx_vec_el[(idx_y1 == idx2_pairs[1, k]) &
      (idx_gy2 == idx2_pairs[2, k])]
    prob_y1gy2[tmp_idx_vec_el] <- y1gy2
  }

  for (k in 1:no_pairs) {
    y1sums <- tapply(biv_prob_split[[k]], idx_cat_y1_split[[k]], sum)
    y1sums_mult <- y1sums[idx_cat_y1_split[[k]]]
    y2gy1 <- biv_prob_split[[k]] / y1sums_mult
    reordered_y2gy1 <- y2gy1[order(idx_cat_y1_split[[k]])]
    tmp_idx_vec_el <- idx_vec_el[(idx_y1 == idx2_pairs[2, k]) &
      (idx_gy2 == idx2_pairs[1, k])]
    prob_y1gy2[tmp_idx_vec_el] <- reordered_y2gy1
  }

  prob_y1gy2
}

# The input of the function is a lavobject, which, in turn, is the output of the
# sem function having specified estimator="PML", missing="available.cases"

# The output of the function is a list of two lists: the pairwiseProbGivObs
# list and the univariateProbGivObs list. Each of the two lists consists of
# G matrices where G is the number of groups in a multigroup analysis.
# If G=1 each of the lists contains only one matrix that can be called as
# pairwiseProbGivObs[[1]], and univariateProbGivObs[[1]].

# Each of the matrices in the pairwiseProbGivObs list is of dimension:
# nrow=sample size, ncol=sum of the number of response categories for all
# pairs of variables (i.e. the length of the vector pxixj.ab where
# i<j=1,...,p, a=1,...,Ci, b=1,...,Cj;
# a which is the index for the response category for yi variable runs the
# fastest, then b which is the index for the response category for yj variable,
# then j, and last i.)
# The cells in a matrix of the pairwiseProbGivObs list have the value 0 except
# for those cells that correspond to the pairs of variables where both variables
# are missing. Those cells have the value of the bivariate conditional
# probability for the given pair for all their response categories.
# The bivariate probabilities are computed as follows:
# the information in the observed variables is summarised in a factor score
# for each individual and the bivariate probability given the estimated factor
# scores is computed.


# Each of the matrices in the univariateProbGivObs list is of dimension:
# nrow=sample size, ncol=sum of the number of response categories for all
# variables. The columns are indexed with i and a, where i=1,...,p, and
# a=1,...,Ci, the response categories for yi variable; a runs faster than i.
# The cells in a matrix of the univariateProbGivObs list have the value 0
# except for those cells that correspond to variables with missing values.
# Those cells have the value of the univariate conditional probability for the
# given variable for all its response categories. The univariate conditional
# probabilities are computed as follows:
# given that the bivariate conditional probabilities have been computed we sum
# over the response categories of each variable at a time (i.e. we compute the
# marginals).

# Version 3 - first compute univariate and then bivariate probabilities

lav_pml_object_uni_pairwise_prob <- function(lavobject) {  # nolint
  # compute yhat where yaht=nu + Lamda*eta + K*x where the parameter
  # estimates are used and the factor scores for eta
  # Below yhat is a list if lavobject@Data@ngroups >1, it is a list of
  # G matrices where G the number of groups and the matrices are of
  # dimension nrow=sample size and ncol=number of items.
  # If lavobject@Data@ngroups=1 then yhat is a matrix.
  yhat <- lavPredict(object = lavobject, type = "yhat")

  # compute bivariate probabilities
  ngroups <- lavobject@Data@ngroups
  univariate_prob <- vector("list", length = ngroups)
  pairwise_prob <- vector("list", length = ngroups)
  # save the indices of the Theta matrices for the groups stored in GLIST
  idx_theta_mat <- which(names(lavobject@Model@GLIST) == "theta")

  for (g in seq_len(ngroups)) { # g<-1

    if (ngroups > 1L) {
      yhat_group <- yhat[[g]]
    } else {
      yhat_group <- yhat
    }

    nsize <- lavobject@Data@nobs[[g]]
    nvar <- lavobject@Model@nvar[[g]]
    data_1 <- lavobject@Data@X[[g]]
    th <- lavobject@Fit@TH[[g]]
    th_idx <- lavobject@Model@th.idx[[g]]
    theta <- lavobject@Model@GLIST[idx_theta_mat[g]]$theta
    error_stddev <- diag(theta)^0.5

    # for the computation of the univariate probabilities
    nlev <- lavobject@Data@ov$nlev
    idx_uniy <- rep(1:nvar, times = nlev)

    # indices vectors for the computation of bivariate probabilities
    idx_pairs_yiyj <- combn(1:nvar, 2)
    no_biv_resp_cat_yiyj <- sapply(seq_len(ncol(idx_pairs_yiyj)),
     function(x) {
      prod(nlev[idx_pairs_yiyj[, x]])
    })
    idx_y1 <- unlist(
      mapply(rep, idx_pairs_yiyj[1, ], each = no_biv_resp_cat_yiyj)
    )
    idx_y2 <- unlist(
      mapply(rep, idx_pairs_yiyj[2, ], each = no_biv_resp_cat_yiyj)
    )


    univariate_prob[[g]] <- matrix(0, nrow = nsize, ncol = sum(nlev))
    pairwise_prob[[g]] <- matrix(0,
      nrow = nsize,
      ncol = length(lavobject@Cache[[g]]$bifreq)
    )

    idx_miss_var_casewise <- apply(data_1, 1, function(x) {
      which(is.na(x))
    })

    for (i in 1:nsize) {
      idx_miss_var <- idx_miss_var_casewise[[i]]
      no_miss_var <- length(idx_miss_var)

      if (no_miss_var > 0L) {
        # compute the univariate probabilities
        th_list <- split(th, th_idx)
        tmp_th <- th_list[idx_miss_var]
        tmp_lower_th <- unlist(lapply(tmp_th, function(x) {
          c(-Inf, x)
        }))
        tmp_upper_th <- unlist(lapply(tmp_th, function(x) {
          c(x, Inf)
        }))

        idx_items <- rep(c(1:no_miss_var), times = nlev[idx_miss_var])
        tmp_mean <- yhat_group[i, idx_miss_var]
        tmp_mean_extended <- tmp_mean[idx_items]
        tmp_stddev <- error_stddev[idx_miss_var]
        tmp_stddev_extended <- tmp_stddev[idx_items]

        tmp_uni_prob <- pnorm((tmp_upper_th - tmp_mean_extended) /
          tmp_stddev_extended) -
          pnorm((tmp_lower_th - tmp_mean_extended) /
            tmp_stddev_extended)
        idx_columns_uni <- which(idx_uniy %in% idx_miss_var)
        univariate_prob[[g]][i, idx_columns_uni] <- tmp_uni_prob

        # compute the bivariate probabilities
        if (no_miss_var > 1L) {
          idx_pairs_miss <- combn(idx_miss_var, 2)
          no_pairs <- ncol(idx_pairs_miss)
          idx_pairs_v2 <- combn(no_miss_var, 2)
          idx_columns <- unlist(lapply(1:no_pairs, function(x) {
            which((idx_y1 == idx_pairs_miss[1, x]) &
              (idx_y2 == idx_pairs_miss[2, x]))
          }))

          if (all(theta[t(idx_pairs_miss)] == 0)) {
            # items independence given eta
            tmp_uni_prob_list <- split(tmp_uni_prob, idx_items)
            pairwise_prob[[g]][i, idx_columns] <-
              unlist(lapply(1:no_pairs, function(x) {
                c(outer(
                  tmp_uni_prob_list[[idx_pairs_v2[1, x]]],
                  tmp_uni_prob_list[[idx_pairs_v2[2, x]]]
                ))
              }))
          } else { # when correlation between measurement errors

            tmp_th_idx <- th_idx[th_idx %in% idx_miss_var]
            # recode so that it is always 1,1,..,1, 2,...,2, etc.
            tmp_th_idx_recoded <- rep(c(1:no_miss_var),
                                  times = table(tmp_th_idx))
            tmp_th <- th[th_idx %in% idx_miss_var]

            tmp_ind_vec <- lav_pml_longvec_ind(
              no.x = no_miss_var,
              all.thres = tmp_th,
              index.var.of.thres = tmp_th_idx_recoded
            )

            tmp_th_rho_vec <- lav_pml_th_rho_generalised(
              no_x = no_miss_var,
              th = tmp_th,
              th_idx = tmp_th_idx_recoded,
              cov_xixj = theta[t(idx_pairs_miss)],
              mean_x = yhat_group[i, idx_miss_var],
              stddev_x = error_stddev[idx_miss_var]
            )

            tmp_biv_prob <- lav_pml_expprob_vec(
              ind.vec = tmp_ind_vec,
              th.rho.vec = tmp_th_rho_vec
            )

            pairwise_prob[[g]][i, idx_columns] <- tmp_biv_prob
          } # end of else of if( all( Theta[t(idx.pairsMiss)]==0 ) )
          # which checks item local independence
        } # end of if( noMissVar>1L )

        # cat(i,  "\n")
      } # end of if(noMissVar>0L)
    } # end of for(i in 1:nsize)
  } # end of for(g in seq_len(lavobject@Data@ngroups))

  list(
    univariateProbGivObs = univariate_prob,
    pairwiseProbGivObs = pairwise_prob
  )
} # end of the function lav_pml_object_uni_pairwise_prob

##################################################################



# lav_pml_th_rho_generalised function is defined as follows
lav_pml_th_rho_generalised <- function(no_x, th, th_idx,
                                      cov_xixj, mean_x, stddev_x) {
  all_std_thres <- (th - mean_x[th_idx]) / stddev_x[th_idx]
  id_pairs <- utils::combn(no_x, 2)
  cor_xixj <- cov_xixj / (stddev_x[id_pairs[1, ]] * stddev_x[id_pairs[2, ]])

  lav_pml_longvec_th_rho(
    no.x = no_x,
    all.thres = all_std_thres,
    index.var.of.thres = th_idx,
    rho.xixj = cor_xixj
  )
}

# lav_pml_th_rho_generalised is a generalisation  of the function
# lav_pml_longvec_th_rho . The latter assumes that all y* follow standard
# normal so the thresholds are automatically the standardised ones.
# lav_pml_th_rho_generalised does not assume that, each of y*'s can follow
# a normal distribution with mean mu and standard deviation sigma.
# lav_pml_th_rho_generalised has the following input arguments:
# no.x (same as in lav_pml_longvec_th_rho),
# TH (similar to the TH in lav_pml_longvec_th_rho but here they are the
#  unstandardised thresholds, i.e. of the normal distribution with
#  mean mu and standard deviation sigma)
# th.idx (same as index.var.of.thres in lav_pml_longvec_th_rho)
# cov.xixj which are the polychoric covariances of the pairs of underlying
# variables provided in a similar fashion as rho.xixj in lav_pml_longvec_th_rho)
# mean.x  is a vector including the means of y*'s provided
#   in the order mean.x1, mean.x2, ...., mean.xp
# stddev.x  is a vector including the standard deviations of y*'s provided
#   in the order stddev.x1, stddev.x2, ...., stddev.xp

# The output of the new function is similar to that of lav_pml_longvec_th_rho
# #############################################



# lavobject is the output of lavaan function where either the unconstrained
# or a hypothesized model has been fitted
lav_pml_object_uni_pairwise_unconstr <- function(lavobject) {  # nolint
  ngroups <- lavobject@Data@ngroups
  th <- lavobject@implied$th # these are the standardized thresholds
  # mean and variance of y* have been taken into account
  th_idx <- lavobject@SampleStats@th.idx
  sigma_hat <- lavobject@implied$cov

  univariate_prob <- vector("list", length = ngroups)
  pairwise_prob <- vector("list", length = ngroups)

  for (g in 1:ngroups) {
    sigma_hat_g <- sigma_hat[[g]]
    # is Sigma.hat always a correlation matrix?
    cor_hat_g <- cov2cor(sigma_hat_g)
    cors <- cor_hat_g[lower.tri(cor_hat_g)]
    if (any(abs(cors) > 1)) {
      lav_msg_warn(gettext(
      "some model-implied correlations are larger than 1.0"))
    }
    nvar <- nrow(sigma_hat_g)
    mean_1 <- rep(0, nvar)
    th_g <- th[[g]]
    th_idx_g <- th_idx[[g]]

    nlev <- lavobject@Data@ov$nlev

    # create index vector to keep track which variable each column of
    # univariateProb matrix refers to
    idx_uniy <- rep(1:nvar, times = nlev)

    # create index vector to keep track which variables each column of
    # pairwiseProb matrix refers to
    idx_pairs_yiyj <- combn(1:nvar, 2)
    no_biv_resp_cat_yiyj <- sapply(seq_len(ncol(idx_pairs_yiyj)),
     function(x) {
      prod(nlev[idx_pairs_yiyj[, x]])
    })
    idx_y1 <- unlist(
      mapply(rep, idx_pairs_yiyj[1, ], each = no_biv_resp_cat_yiyj)
    )
    idx_y2 <- unlist(
      mapply(rep, idx_pairs_yiyj[2, ], each = no_biv_resp_cat_yiyj)
    )

    data_1 <- lavobject@Data@X[[g]]
    nsize <- nrow(data_1)

    # create the lists of matrices
    univariate_prob[[g]] <- matrix(0, nrow = nsize, ncol = sum(nlev))
    pairwise_prob[[g]] <- matrix(0,
      nrow = nsize,
      ncol = length(lavobject@Cache[[g]]$bifreq)
    )

    idx_miss_var_casewise <- apply(data_1, 1, function(x) {
      which(is.na(x))
    })

    for (i in 1:nsize) {
      idx_miss_var <- idx_miss_var_casewise[[i]]
      no_miss_var <- length(idx_miss_var)

      if (no_miss_var > 0L) {
        # compute the denominator of the conditional probability
        th_var <- lapply(1:nvar, function(x) c(-Inf, th_g[th_idx_g == x], +Inf))
        lower <- sapply(1:nvar, function(x) th_var[[x]][data_1[i, x]])
        upper <- sapply(1:nvar, function(x) th_var[[x]][data_1[i, x] + 1L])
        lower_denom <- lower[-idx_miss_var]
        upper_denom <- upper[-idx_miss_var]
        mean_i <- mean_1[-idx_miss_var]
        corhat_i <- cor_hat_g[-idx_miss_var, -idx_miss_var, drop = FALSE]
        denom <- sadmvn(lower_denom, upper_denom, mean = mean_i,
           varcov = corhat_i)[1]
      } # end of if( noMissVar>0L )

      if (no_miss_var == 1L) { # only univariate probabilities for one item
        # compute the numerator
        th_miss_var <- c(-Inf, th_g[th_idx_g == idx_miss_var], +Inf)
        # for all response categories of the missing item
        no_cat <- nlev[idx_miss_var]
        numer <- sapply(1:no_cat, function(x) {
          lower[idx_miss_var] <- th_miss_var[x]
          upper[idx_miss_var] <- th_miss_var[x + 1L]
          sadmvn(lower, upper, mean = mean_1, varcov = cor_hat_g)[1]
        })
        idx_columns_uni <- which(idx_uniy %in% idx_miss_var)
        univariate_prob[[g]][i, idx_columns_uni] <- numer / denom
      } # end of if( noMissVar==1L )

      if (no_miss_var > 1L) {
        # compute the bivariate probabilities and based on them
        # calculate the univariate ones

        # form all possible pairs of items with missing values
        idx_pairs_miss <- combn(idx_miss_var, 2)
        no_pairs <- ncol(idx_pairs_miss)
        for (j in 1:no_pairs) {
          idx_missy1y2 <- idx_pairs_miss[, j]
          idx_missy1 <- idx_missy1y2[1]
          idx_missy2 <- idx_missy1y2[2]
          idx_miss_rest_items <- idx_miss_var[!(idx_miss_var %in% idx_missy1y2)]
          th_missy1 <- c(-Inf, th_g[th_idx_g == idx_missy1], +Inf)
          th_missy2 <- c(-Inf, th_g[th_idx_g == idx_missy2], +Inf)
          no_cat_missy1 <- nlev[idx_missy1]
          no_cat_missy2 <- nlev[idx_missy2]
          no_biv_resp_cat <- no_cat_missy1 * no_cat_missy2
          mat_biv_resp_cat <- matrix(1:no_biv_resp_cat,
            nrow = no_cat_missy1,
            ncol = no_cat_missy2
          )

          numer <- sapply(1:no_biv_resp_cat, function(x) {
            idx_y1_cat <- which(mat_biv_resp_cat == x, arr.ind = TRUE)[1]
            idx_y2_cat <- which(mat_biv_resp_cat == x, arr.ind = TRUE)[2]
            lower[idx_missy1y2] <-
              c(th_missy1[idx_y1_cat], th_missy2[idx_y2_cat])
            upper[idx_missy1y2] <-
              c(th_missy1[idx_y1_cat + 1L], th_missy2[idx_y2_cat + 1L])
            lower_tmp <- lower
            upper_tmp <- upper
            mean_tmp <- mean_1
            cor_hat_g_tmp <- cor_hat_g
            if (length(idx_miss_rest_items) > 0) {
              lower_tmp <- lower[-idx_miss_rest_items]
              upper_tmp <- upper[-idx_miss_rest_items]
              mean_tmp <- mean_1[-idx_miss_rest_items]
              cor_hat_g_tmp <-
                    cor_hat_g[-idx_miss_rest_items, -idx_miss_rest_items]
            }
            sadmvn(lower_tmp, upper_tmp,
              mean = mean_tmp, varcov = cor_hat_g_tmp
            )[1]
          })

          idx_columns <- which((idx_y1 == idx_missy1) &
            (idx_y2 == idx_missy2))
          tmp_biv <- numer / denom
          pairwise_prob[[g]][i, idx_columns] <- tmp_biv

          # compute the univariateProb based on the above bivariate
          # probabilities
          if (j == 1L) {
            univariate_prob[[g]][i, which(idx_uniy %in% idx_missy1)] <-
              apply(mat_biv_resp_cat, 1, function(x) {
                sum(tmp_biv[x])
              })

            univariate_prob[[g]][i, which(idx_uniy %in% idx_missy2)] <-
              apply(mat_biv_resp_cat, 2, function(x) {
                sum(tmp_biv[x])
              })
          }

          if (j > 1L && j < no_miss_var) {
            univariate_prob[[g]][i, which(idx_uniy %in% idx_missy2)] <-
              apply(mat_biv_resp_cat, 2, function(x) {
                sum(tmp_biv[x])
              })
          }
        } # end of for(j in 1:no.pairs ) #no.pairs is that of missing items
      } # end of if( noMissVar>1L )
    } # end of for(i in 1:nsize)
  } # end of for(g in 1:ngroups)

  list(
    univariateProbGivObs = univariate_prob,
    pairwiseProbGivObs = pairwise_prob
  )
} # end of function
