# contributed by Myrsini Katsikatsou (March 2016)

# the function lav_pml_bi_lik_x gives the value of the bivariate likelihood
# for a specific pair of ordinal variables casewise when covariates are present and estimator=="PML"
# (the bivariate likelihood is essentially the bivariate probability of the
# observed response pattern of two ordinal variables)

# Input arguments:
# Y1 is a vector, includes the observed values for the first variable for all cases/units,
# Y1 is ordinal
# Y2 similar to Y1
# Rho is the polychoric correlation of Y1 and Y2
# th.y1 is the vector of the thresholds for Y1* excluding the first and
# the last thresholds which are -Inf and Inf
# th.y2 is similar to th.y1
# eXo is the data for the covariates in a matrix format where nrows= no of cases,
# ncols= no of covariates
# PI.y1 is a vector, includes the regression coefficients of the covariates
# for the first variable, Y1, the length of the vector is the no of covariates;
# to obtain this vector apply the function lavaan:::lav_model_pi()[row_correspondin_to_Y1, ]
# PI.y2 is similar to PI.y2
# missing.ind is of "character" value, taking the values listwise, pairwise, available_cases;
# to obtain a value use lavdata@missing

# Output:
# It is a vector, length= no of cases, giving the bivariate likelihood for each case.
lav_pml_bi_lik_x <- function(Y1, Y2, Rho,
                               th.y1, th.y2,
                               eXo,
                               PI.y1, PI.y2,
                               missing.ind) {
  th.y1 <- c(-100, th.y1, 100)
  th.y2 <- c(-100, th.y2, 100)
  pred.y1 <- c(eXo %*% PI.y1)
  pred.y2 <- c(eXo %*% PI.y2)

  th.y1.upper <- th.y1[Y1 + 1L] - pred.y1
  th.y1.lower <- th.y1[Y1] - pred.y1
  th.y2.upper <- th.y2[Y2 + 1L] - pred.y2
  th.y2.lower <- th.y2[Y2] - pred.y2

  if (missing.ind == "listwise") { # I guess this is the default which
    # also handles the case of complete data
    biv_prob <- pbivnorm(th.y1.upper, th.y2.upper, rho = Rho) -
      pbivnorm(th.y1.lower, th.y2.upper, rho = Rho) -
      pbivnorm(th.y1.upper, th.y2.lower, rho = Rho) +
      pbivnorm(th.y1.lower, th.y2.lower, rho = Rho)
    lik <- biv_prob
  } else if (missing.ind %in% c(
    "pairwise",
    "available.cases",
    "available_cases"
  )) {
    # index of cases with complete pairs
    CP.idx <- which(complete.cases(cbind(Y1, Y2)))

    th.y1.upper <- th.y1.upper[CP.idx]
    th.y1.lower <- th.y1.lower[CP.idx]
    th.y2.upper <- th.y2.upper[CP.idx]
    th.y2.lower <- th.y2.lower[CP.idx]

    biv_prob <- pbivnorm(th.y1.upper, th.y2.upper, rho = Rho) -
      pbivnorm(th.y1.lower, th.y2.upper, rho = Rho) -
      pbivnorm(th.y1.upper, th.y2.lower, rho = Rho) +
      pbivnorm(th.y1.lower, th.y2.lower, rho = Rho)

    # lik <- numeric( length(Y1) )
    lik <- rep(as.numeric(NA), length(Y1))
    lik[CP.idx] <- biv_prob
  }
  lik
}

#################################################################


# The function  lav_pml_uni_lik gives the value of the univariate likelihood for a
# specific ordinal variable, casewise (which is essentially the probability for
# the observed response category for each case).
# The input arguments are explained before the function lav_pml_bi_lik_x above.
# Output:
# It is a vector, length= no of cases, giving the univariate likelihoods for each case.

lav_pml_uni_lik <- function(Y1, th.y1, eXo = NULL, PI.y1 = NULL) {
  th.y1 <- c(-100, th.y1, 100)
  if (!is.null(eXo)) {
    pred.y1 <- c(eXo %*% PI.y1)
  }

  if (is.null(eXo)) {
    th.y1.upper <- th.y1[Y1 + 1L]
    th.y1.lower <- th.y1[Y1]
  } else {
    th.y1.upper <- th.y1[Y1 + 1L] - pred.y1
    th.y1.lower <- th.y1[Y1] - pred.y1
  }

  lav_pml_uni_lik <- pnorm(th.y1.upper) - pnorm(th.y1.lower)

  lav_pml_uni_lik[is.na(lav_pml_uni_lik)] <- 0

  lav_pml_uni_lik
}

#################################################################

# The function lav_pml_th_uni_prob gives the model-based univariate probabilities
# for all ordinal indicators and for all of their response categories, i.e. pi(xi=a), where
# a=1,...,ci and i=1,...,p with a index running faster than i index.
# Input arguments:
# TH is a vector giving the thresholds for all variables, tau_ia, with a running
#    faster than i (the first and the last thresholds which are -Inf and Inf are
#    not included). TH can be given by the lavaan function lav_model_th .
# th.idx is a vector of same length as TH which gives the value of the i index,
#        namely which variable each thresholds refers to. This can be obtained by
#        lavmodel@th.idx .
# Output:
# It is a vector, lenght= Sum_i(ci), i.e. the sum of the response categories of
# all ordinal variables. The vector contains the model-based univariate probabilities pi(xi=a).

lav_pml_th_uni_prob <- function(TH = TH, th.idx = th.idx) {
  TH.split <- split(TH, th.idx)
  TH.lower <- unlist(lapply(TH.split, function(x) {
    c(-100, x)
  }), use.names = FALSE)
  TH.upper <- unlist(lapply(TH.split, function(x) {
    c(x, 100)
  }), use.names = FALSE)
  prob <- pnorm(TH.upper) - pnorm(TH.lower)
  # to avoid Nan/-Inf
  prob[prob < .Machine$double.eps] <- .Machine$double.eps
  prob
}

#############################################################################


# The function lav_pml_dbilogl_dpar_x computes the derivatives of a bivariate
# log-likelihood of two ordinal variables casewise with respect to thresholds,
# slopes (reduced-form regression coefficients for the covariates), and polychoric correlation.
# The function lav_dbinorm of lavaan is used.
# The function gives the right result for both listwise and pairwise deletion,
# and the case of complete data.
# Input arguments are explained before the function lav_pml_bi_lik_x defined above.
# The only difference is that PI.y1 and PI.y2 are (accidentally) renamed here as sl.y1 and sl.y2
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


lav_pml_dbilogl_dpar_x <- function(Y1, Y2, eXo, Rho,
                                      th.y1, th.y2,
                                      sl.y1, sl.y2,
                                      missing.ind) {
  nth.y1 <- length(th.y1)
  nth.y2 <- length(th.y2)

  start.th.y1 <- th.y1
  start.th.y2 <- th.y2

  Nobs <- length(Y1)

  R <- sqrt(1 - Rho * Rho)
  th.y1 <- c(-100, th.y1, 100)
  th.y2 <- c(-100, th.y2, 100)
  pred.y1 <- c(eXo %*% sl.y1)
  pred.y2 <- c(eXo %*% sl.y2)

  th.y1.z1 <- th.y1[Y1 + 1L] - pred.y1
  th.y1.z2 <- th.y1[Y1] - pred.y1
  th.y2.z1 <- th.y2[Y2 + 1L] - pred.y2
  th.y2.z2 <- th.y2[Y2] - pred.y2

  # lik, i.e. the bivariate probability case-wise
  lik <- lav_pml_bi_lik_x(
    Y1 = Y1, Y2 = Y2,
    Rho = Rho,
    th.y1 = start.th.y1,
    th.y2 = start.th.y2,
    eXo = eXo,
    PI.y1 = sl.y1,
    PI.y2 = sl.y2,
    missing.ind = missing.ind
  )


  # w.r.t. th.y1, mean tau tilde
  # derivarive bivariate prob w.r.t. tau^xi_ci, see formula in paper 2012
  y1.Z1 <- dnorm(th.y1.z1) * (pnorm((th.y2.z1 - Rho * th.y1.z1) / R) -
    pnorm((th.y2.z2 - Rho * th.y1.z1) / R))
  # derivarive bivariate prob w.r.t. tau^xi_(ci-1),
  y1.Z2 <- (-1) * (dnorm(th.y1.z2) * (pnorm((th.y2.z1 - Rho * th.y1.z2) / R) -
    pnorm((th.y2.z2 - Rho * th.y1.z2) / R)))


  # allocate the derivatives at the right column casewise
  idx.y1.z1 <- matrix(1:nth.y1, nrow = Nobs, ncol = nth.y1, byrow = TRUE) == Y1
  idx.y1.z2 <- matrix(1:nth.y1, nrow = Nobs, ncol = nth.y1, byrow = TRUE) == (Y1 - 1L)
  der.table.y1 <- idx.y1.z1 * y1.Z1 + idx.y1.z2 * y1.Z2

  # der of pl w.r.t. th.y1
  dx.th.tilde.y1 <- der.table.y1 / lik
  dx.th.tilde.y1[is.na(dx.th.tilde.y1)] <- 0

  # w.r.t. th.y2, mean tau tilde
  # derivarive bivariate prob w.r.t. tau^xi_ci, see formula in paper 2012
  y2.Z1 <- dnorm(th.y2.z1) * (pnorm((th.y1.z1 - Rho * th.y2.z1) / R) -
    pnorm((th.y1.z2 - Rho * th.y2.z1) / R))
  # derivarive bivariate prob w.r.t. tau^xi_(ci-1),
  y2.Z2 <- (-1) * (dnorm(th.y2.z2) * (pnorm((th.y1.z1 - Rho * th.y2.z2) / R) -
    pnorm((th.y1.z2 - Rho * th.y2.z2) / R)))
  # allocate the derivatives at the right column casewise
  idx.y2.z1 <- matrix(1:nth.y2, nrow = Nobs, ncol = nth.y2, byrow = TRUE) == Y2
  idx.y2.z2 <- matrix(1:nth.y2, nrow = Nobs, ncol = nth.y2, byrow = TRUE) == (Y2 - 1L)
  der.table.y2 <- idx.y2.z1 * y2.Z1 + idx.y2.z2 * y2.Z2

  # der of pl w.r.t. th.y2
  dx.th.tilde.y2 <- der.table.y2 / lik
  dx.th.tilde.y2[is.na(dx.th.tilde.y2)] <- 0



  # w.r.t. rho
  # derivarive bivariate prob w.r.t. rho, see formula in paper 2012
  dbivprob.wrt.rho <- (lav_dbinorm(th.y1.z1, th.y2.z1, Rho) -
    lav_dbinorm(th.y1.z2, th.y2.z1, Rho) -
    lav_dbinorm(th.y1.z1, th.y2.z2, Rho) +
    lav_dbinorm(th.y1.z2, th.y2.z2, Rho))
  # der of pl w.r.t. rho
  dx.rho <- dbivprob.wrt.rho / lik
  dx.rho[is.na(dx.rho)] <- 0


  # der of pl w.r.t. slopes (also referred to PI obtained by lav_model_pi function)
  row.sums.y1 <- rowSums(dx.th.tilde.y1)
  row.sums.y2 <- rowSums(dx.th.tilde.y2)
  dx.sl.y1 <- (-1) * eXo * row.sums.y1
  dx.sl.y2 <- (-1) * eXo * row.sums.y2


  list(
    dx.th.y1 = dx.th.tilde.y1, # note that dx.th.tilde=dx.th
    dx.th.y2 = dx.th.tilde.y2,
    dx.sl.y1 = dx.sl.y1,
    dx.sl.y2 = dx.sl.y2,
    dx.rho = dx.rho
  )
}

###############################################################


# The function lav_pml_uni_scores gives, casewise, the derivative of a univariate
# log-likelihood w.r.t. thresholds and slopes if present weighted by the
# casewise uni-weights as those defined in AC-PL (essentially the number of missing values per case).
# The function closely follows the "logic" of the function lav_pml_dbilogl_dpar_x defined above.
# Input arguments are as before plus: weights.casewise given by
# lavcavhe$uniweights.casewise .
# Output:
# A list including the following:
# a) the derivatives w.r.t. the thresholds of the variable. This is a matrix,
#    nrows=no of cases, ncols= no of thresholds of variable 1.
# b) the derivatives w.r.t slopes for the variable. If covariates are present,
#    this is a matrix, nrows=no of cases, ncols= no of covariates.
#    Otherwise it takes the value NULL.


lav_pml_uni_scores <- function(Y1, th.y1, eXo = NULL, sl.y1 = NULL,
                       weights.casewise) {
  nth.y1 <- length(th.y1)
  start.th.y1 <- th.y1
  Nobs <- length(Y1)
  th.y1 <- c(-100, th.y1, 100)

  if (is.null(eXo)) {
    th.y1.z1 <- th.y1[Y1 + 1L]
    th.y1.z2 <- th.y1[Y1]
  } else {
    pred.y1 <- c(eXo %*% sl.y1)
    th.y1.z1 <- th.y1[Y1 + 1L] - pred.y1
    th.y1.z2 <- th.y1[Y1] - pred.y1
  }

  # lik, i.e. the univariate probability case-wise
  lik <- lav_pml_uni_lik( # Y1 = X[,i],
    Y1 = Y1,
    # th.y1 = TH[th.idx==i],
    th.y1 = th.y1,
    eXo = eXo,
    # PI.y1 = PI[i,])
    PI.y1 = sl.y1
  )

  # w.r.t. th.y1
  # derivarive of the univariate prob w.r.t. to the upper limit threshold
  y1.Z1 <- dnorm(th.y1.z1)
  # derivarive of the univariate prob w.r.t. to the lower limit threshold
  y1.Z2 <- (-1) * dnorm(th.y1.z2)

  # allocate the derivatives at the right column casewise
  idx.y1.z1 <- matrix(1:nth.y1, nrow = Nobs, ncol = nth.y1, byrow = TRUE) == Y1
  idx.y1.z2 <- matrix(1:nth.y1, nrow = Nobs, ncol = nth.y1, byrow = TRUE) == (Y1 - 1L)
  der.table.y1 <- idx.y1.z1 * y1.Z1 + idx.y1.z2 * y1.Z2

  # der of pl w.r.t. th.y1
  dx.th.tilde.y1 <- der.table.y1 * (weights.casewise / lik)
  dx.th.tilde.y1[is.na(dx.th.tilde.y1)] <- 0

  # der of pl w.r.t. slopes (also referred to PI obtained by lav_model_pi function)
  dx.sl.y1 <- NULL
  if (!is.null(eXo)) {
    row.sums.y1 <- rowSums(dx.th.tilde.y1)
    dx.sl.y1 <- (-1) * eXo * row.sums.y1
  }

  list(
    dx.th.y1 = dx.th.tilde.y1, # note that dx.th.tilde=dx.th
    dx.sl.y1 = dx.sl.y1
  )
}
