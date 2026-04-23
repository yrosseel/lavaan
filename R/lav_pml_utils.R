# contributed by Myrsini Katsikatsou (March 2016)

# the function lav_pml_bi_lik_x gives the value of the bivariate likelihood
# for a specific pair of ordinal variables casewise when covariates are present
# and estimator=="PML" (the bivariate likelihood is essentially the bivariate
# probability of the observed response pattern of two ordinal variables)

# Input arguments:
# Y1 is a vector, includes the observed values for the first variable for all
#  cases/units, Y1 is ordinal
# Y2 similar to Y1
# Rho is the polychoric correlation of Y1 and Y2
# th.y1 is the vector of the thresholds for Y1* excluding the first and
# the last thresholds which are -Inf and Inf
# th.y2 is similar to th.y1
# eXo is the data for the covariates in a matrix format where
# nrows= no of cases, ncols= no of covariates
# PI.y1 is a vector, includes the regression coefficients of the covariates
# for the first variable, Y1, the length of the vector is the no of covariates;
# to obtain this vector apply the function
#                                   lav_model_pi()[row_correspondin_to_Y1, ]
# PI.y2 is similar to PI.y2
# missing.ind is of "character" value, taking the values listwise, pairwise,
# available_cases; to obtain a value use lavdata@missing

# Output:
# It is a vector, length= no of cases, giving the bivariate likelihood
#  for each case.
lav_pml_bi_lik_x <- function(y1, y2, rho,
                               th_y1, th_y2,
                               exo,
                               pi_y1, pi_y2,
                               missing_ind) {
  th_y1 <- c(-100, th_y1, 100)
  th_y2 <- c(-100, th_y2, 100)
  pred_y1 <- c(exo %*% pi_y1)
  pred_y2 <- c(exo %*% pi_y2)

  th_y1_upper <- th_y1[y1 + 1L] - pred_y1
  th_y1_lower <- th_y1[y1] - pred_y1
  th_y2_upper <- th_y2[y2 + 1L] - pred_y2
  th_y2_lower <- th_y2[y2] - pred_y2

  if (missing_ind == "listwise") { # I guess this is the default which
    # also handles the case of complete data
    biv_prob <- pbivnorm(th_y1_upper, th_y2_upper, rho = rho) -
      pbivnorm(th_y1_lower, th_y2_upper, rho = rho) -
      pbivnorm(th_y1_upper, th_y2_lower, rho = rho) +
      pbivnorm(th_y1_lower, th_y2_lower, rho = rho)
    lik <- biv_prob
  } else if (missing_ind %in% c(
    "pairwise",
    "available.cases",
    "available_cases"
  )) {
    # index of cases with complete pairs
    cp_idx <- which(complete.cases(cbind(y1, y2)))

    th_y1_upper <- th_y1_upper[cp_idx]
    th_y1_lower <- th_y1_lower[cp_idx]
    th_y2_upper <- th_y2_upper[cp_idx]
    th_y2_lower <- th_y2_lower[cp_idx]

    biv_prob <- pbivnorm(th_y1_upper, th_y2_upper, rho = rho) -
      pbivnorm(th_y1_lower, th_y2_upper, rho = rho) -
      pbivnorm(th_y1_upper, th_y2_lower, rho = rho) +
      pbivnorm(th_y1_lower, th_y2_lower, rho = rho)

    # lik <- numeric( length(Y1) )
    lik <- rep(as.numeric(NA), length(y1))
    lik[cp_idx] <- biv_prob
  }
  lik
}

#################################################################


# The function  lav_pml_uni_lik gives the value of the univariate likelihood
# for a specific ordinal variable, casewise (which is essentially the
# probability for the observed response category for each case).
# The input arguments are explained before the function lav_pml_bi_lik_x above.
# Output:
# It is a vector, length= no of cases, giving the univariate likelihoods
# for each case.

lav_pml_uni_lik <- function(y1, th_y1, exo = NULL, pi_y1 = NULL) {
  th_y1 <- c(-100, th_y1, 100)
  if (!is.null(exo)) {
    pred_y1 <- c(exo %*% pi_y1)
  }

  if (is.null(exo)) {
    th_y1_upper <- th_y1[y1 + 1L]
    th_y1_lower <- th_y1[y1]
  } else {
    th_y1_upper <- th_y1[y1 + 1L] - pred_y1
    th_y1_lower <- th_y1[y1] - pred_y1
  }

  pml_uni_lik <- pnorm(th_y1_upper) - pnorm(th_y1_lower)

  pml_uni_lik[is.na(pml_uni_lik)] <- 0

  pml_uni_lik
}

#################################################################

# The function lav_pml_th_uni_prob gives the model-based univariate
# probabilities for all ordinal indicators and for all of their response
# categories, i.e. pi(xi=a), where a=1,...,ci and i=1,...,p with a index
# running faster than i index.
# Input arguments:
# TH is a vector giving the thresholds for all variables, tau_ia, with a running
#    faster than i (the first and the last thresholds which are -Inf and Inf are
#    not included). TH can be given by the lavaan function lav_model_th .
# th.idx is a vector of same length as TH which gives the value of the i index,
#      namely which variable each thresholds refers to. This can be obtained by
#      lavmodel@th.idx .
# Output:
# It is a vector, lenght= Sum_i(ci), i.e. the sum of the response categories of
# all ordinal variables. The vector contains the model-based univariate
# probabilities pi(xi=a).

lav_pml_th_uni_prob <- function(th = th, th_idx = th_idx) {
  th_split <- split(th, th_idx)
  th_lower <- unlist(lapply(th_split, function(x) {
    c(-100, x)
  }), use.names = FALSE)
  th_upper <- unlist(lapply(th_split, function(x) {
    c(x, 100)
  }), use.names = FALSE)
  prob <- pnorm(th_upper) - pnorm(th_lower)
  # to avoid Nan/-Inf
  prob[prob < .Machine$double.eps] <- .Machine$double.eps
  prob
}

#############################################################################


# The function lav_pml_dbilogl_dpar_x computes the derivatives of a bivariate
# log-likelihood of two ordinal variables casewise with respect to thresholds,
# slopes (reduced-form regression coefficients for the covariates), and
# polychoric correlation.
# The function lav_dbinorm of lavaan is used.
# The function gives the right result for both listwise and pairwise deletion,
# and the case of complete data.
# Input arguments are explained before the function lav_pml_bi_lik_x
# defined above. The only difference is that PI.y1 and PI.y2 are (accidentally)
# renamed here as sl.y1 and sl.y2
# Output:
# It is a list containing the following
# a) the derivatives w.r.t. the thresholds of the first variable casewise.
#    This is a matrix, nrows=no of cases, ncols= no of thresholds of variable 1.
# b) the derivatives w.r.t. the thresholds of the second variable casewise.
#    This is a matrix, nrows=no of cases, ncols= no of thresholds of variable 2.
# c) the derivatives w.r.t slopes for variable 1.  This is a matrix, where
#    nrows=no of cases, ncols= no of covariates.
# d) the derivatives w.r.t slopes for variable 2.  This is a matrix, where
#    nrows=no of cases, ncols= no of covariates.
# e) the derivative w.r.t the polychoric correlation of the two variables.
#    This is a vector of length= no of cases.


lav_pml_dbilogl_dpar_x <- function(y1, y2, exo, rho,
                                      th_y1, th_y2,
                                      sl_y1, sl_y2,
                                      missing_ind) {
  nth_y1 <- length(th_y1)
  nth_y2 <- length(th_y2)

  start_th_y1 <- th_y1
  start_th_y2 <- th_y2

  nobs_1 <- length(y1)

  r <- sqrt(1 - rho * rho)
  th_y1 <- c(-100, th_y1, 100)
  th_y2 <- c(-100, th_y2, 100)
  pred_y1 <- c(exo %*% sl_y1)
  pred_y2 <- c(exo %*% sl_y2)

  th_y1_z1 <- th_y1[y1 + 1L] - pred_y1
  th_y1_z2 <- th_y1[y1] - pred_y1
  th_y2_z1 <- th_y2[y2 + 1L] - pred_y2
  th_y2_z2 <- th_y2[y2] - pred_y2

  # lik, i.e. the bivariate probability case-wise
  lik <- lav_pml_bi_lik_x(
    y1 = y1, y2 = y2,
    rho = rho,
    th_y1 = start_th_y1,
    th_y2 = start_th_y2,
    exo = exo,
    pi_y1 = sl_y1,
    pi_y2 = sl_y2,
    missing_ind = missing_ind
  )


  # w.r.t. th.y1, mean tau tilde
  # derivarive bivariate prob w.r.t. tau^xi_ci, see formula in paper 2012
  y1_z1 <- dnorm(th_y1_z1) * (pnorm((th_y2_z1 - rho * th_y1_z1) / r) -
    pnorm((th_y2_z2 - rho * th_y1_z1) / r))
  # derivarive bivariate prob w.r.t. tau^xi_(ci-1),
  y1_z2 <- (-1) * (dnorm(th_y1_z2) * (pnorm((th_y2_z1 - rho * th_y1_z2) / r) -
    pnorm((th_y2_z2 - rho * th_y1_z2) / r)))


  # allocate the derivatives at the right column casewise
  idx_y1_z1 <-
    matrix(1:nth_y1, nrow = nobs_1, ncol = nth_y1, byrow = TRUE) == y1
  idx_y1_z2 <-
    matrix(1:nth_y1, nrow = nobs_1, ncol = nth_y1, byrow = TRUE) == (y1 - 1L)
  der_table_y1 <- idx_y1_z1 * y1_z1 + idx_y1_z2 * y1_z2

  # der of pl w.r.t. th.y1
  dx_th_tilde_y1 <- der_table_y1 / lik
  dx_th_tilde_y1[is.na(dx_th_tilde_y1)] <- 0

  # w.r.t. th.y2, mean tau tilde
  # derivarive bivariate prob w.r.t. tau^xi_ci, see formula in paper 2012
  y2_z1 <- dnorm(th_y2_z1) * (pnorm((th_y1_z1 - rho * th_y2_z1) / r) -
    pnorm((th_y1_z2 - rho * th_y2_z1) / r))
  # derivarive bivariate prob w.r.t. tau^xi_(ci-1),
  y2_z2 <- (-1) * (dnorm(th_y2_z2) * (pnorm((th_y1_z1 - rho * th_y2_z2) / r) -
    pnorm((th_y1_z2 - rho * th_y2_z2) / r)))
  # allocate the derivatives at the right column casewise
  idx_y2_z1 <-
    matrix(1:nth_y2, nrow = nobs_1, ncol = nth_y2, byrow = TRUE) == y2
  idx_y2_z2 <-
    matrix(1:nth_y2, nrow = nobs_1, ncol = nth_y2, byrow = TRUE) == (y2 - 1L)
  der_table_y2 <- idx_y2_z1 * y2_z1 + idx_y2_z2 * y2_z2

  # der of pl w.r.t. th.y2
  dx_th_tilde_y2 <- der_table_y2 / lik
  dx_th_tilde_y2[is.na(dx_th_tilde_y2)] <- 0



  # w.r.t. rho
  # derivarive bivariate prob w.r.t. rho, see formula in paper 2012
  dbivprob_wrt_rho <- (lav_dbinorm(th_y1_z1, th_y2_z1, rho) -
    lav_dbinorm(th_y1_z2, th_y2_z1, rho) -
    lav_dbinorm(th_y1_z1, th_y2_z2, rho) +
    lav_dbinorm(th_y1_z2, th_y2_z2, rho))
  # der of pl w.r.t. rho
  dx_rho <- dbivprob_wrt_rho / lik
  dx_rho[is.na(dx_rho)] <- 0


  # der of pl w.r.t. slopes (also referred to PI obtained by
  #     lav_model_pi function)
  row_sums_y1 <- rowSums(dx_th_tilde_y1)
  row_sums_y2 <- rowSums(dx_th_tilde_y2)
  dx_sl_y1 <- (-1) * exo * row_sums_y1
  dx_sl_y2 <- (-1) * exo * row_sums_y2


  list(
    dx.th.y1 = dx_th_tilde_y1, # note that dx.th.tilde=dx.th
    dx.th.y2 = dx_th_tilde_y2,
    dx.sl.y1 = dx_sl_y1,
    dx.sl.y2 = dx_sl_y2,
    dx.rho = dx_rho
  )
}

###############################################################


# The function lav_pml_uni_scores gives, casewise, the derivative of a
# univariate log-likelihood w.r.t. thresholds and slopes if present weighted by
# the casewise uni-weights as those defined in AC-PL (essentially the number of
# missing values per case).
# The function closely follows the "logic" of the function
# lav_pml_dbilogl_dpar_x defined above.
# Input arguments are as before plus: weights.casewise given by
# lavcavhe$uniweights.casewise .
# Output:
# A list including the following:
# a) the derivatives w.r.t. the thresholds of the variable. This is a matrix,
#    nrows=no of cases, ncols= no of thresholds of variable 1.
# b) the derivatives w.r.t slopes for the variable. If covariates are present,
#    this is a matrix, nrows=no of cases, ncols= no of covariates.
#    Otherwise it takes the value NULL.


lav_pml_uni_scores <- function(y1, th_y1, exo = NULL, sl_y1 = NULL,
                       weights_casewise) {
  nth_y1 <- length(th_y1)
  nobs_1 <- length(y1)
  th_y1 <- c(-100, th_y1, 100)

  if (is.null(exo)) {
    th_y1_z1 <- th_y1[y1 + 1L]
    th_y1_z2 <- th_y1[y1]
  } else {
    pred_y1 <- c(exo %*% sl_y1)
    th_y1_z1 <- th_y1[y1 + 1L] - pred_y1
    th_y1_z2 <- th_y1[y1] - pred_y1
  }

  # lik, i.e. the univariate probability case-wise
  lik <- lav_pml_uni_lik( # Y1 = X[,i],
    y1 = y1,
    # th_y1 = TH[th.idx==i],
    th_y1 = th_y1,
    exo = exo,
    # pi_y1 = PI[i,])
    pi_y1 = sl_y1
  )

  # w.r.t. th.y1
  # derivarive of the univariate prob w.r.t. to the upper limit threshold
  y1_z1 <- dnorm(th_y1_z1)
  # derivarive of the univariate prob w.r.t. to the lower limit threshold
  y1_z2 <- (-1) * dnorm(th_y1_z2)

  # allocate the derivatives at the right column casewise
  idx_y1_z1 <-
    matrix(1:nth_y1, nrow = nobs_1, ncol = nth_y1, byrow = TRUE) == y1
  idx_y1_z2 <-
    matrix(1:nth_y1, nrow = nobs_1, ncol = nth_y1, byrow = TRUE) == (y1 - 1L)
  der_table_y1 <- idx_y1_z1 * y1_z1 + idx_y1_z2 * y1_z2

  # der of pl w.r.t. th.y1
  dx_th_tilde_y1 <- der_table_y1 * (weights_casewise / lik)
  dx_th_tilde_y1[is.na(dx_th_tilde_y1)] <- 0

  # der of pl w.r.t. slopes (also referred to PI obtained by
  #     lav_model_pi function)
  dx_sl_y1 <- NULL
  if (!is.null(exo)) {
    row_sums_y1 <- rowSums(dx_th_tilde_y1)
    dx_sl_y1 <- (-1) * exo * row_sums_y1
  }

  list(
    dx.th.y1 = dx_th_tilde_y1, # note that dx.th.tilde=dx.th
    dx.sl.y1 = dx_sl_y1
  )
}
