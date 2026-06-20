# old not-exported functions which are temporarily kept to enable
# publishing lavaan 0.7-1 on CRAN
lav_warning_nonexported <- function(funcname) {
  return(invisible(NULL))
  lav_msg_warn(gettextf(
    "The function %s is an internal function of lavaan and NOT meant to be
    called directly from outside lavaan. This function will be removed in
    the near future.", funcname
  ))
}

lav_partable_ngroups <- function(partable) {
  lav_warning_nonexported("lav_partable_ngroups")
  lav_pt_ngroups(partable)
}
lav_mvnorm_loglik_samplestats <- function(
  sample.mean = NULL,
  sample.cov = NULL,
  sample.nobs = NULL,
  Mu = NULL,
  Sigma = NULL,
  x.idx = integer(0L),
  x.mean = NULL,
  x.cov = NULL,
  Sinv.method = "eigen",
  Sigma.inv = NULL
) {
  lav_warning_nonexported("lav_mvnorm_loglik_samplestats")
  lav_mvn_loglik_samp(
    sample_mean = sample.mean,
    sample_cov = sample.cov,
    sample_nobs = sample.nobs,
    mu = Mu,
    sigma_1 = Sigma,
    x_idx = x.idx,
    x_mean = x.mean,
    x_cov = x.cov,
    sinv_method = Sinv.method,
    sigma_inv = Sigma.inv
  )
}
lav_mvnorm_missing_loglik_samplestats <- function(
  Yp = NULL,
  Mu = NULL,
  Sigma = NULL,
  x.idx = integer(0L),
  x.mean = NULL,
  x.cov = NULL,
  Sinv.method = "eigen",
  log2pi = TRUE,
  minus.two = FALSE
) {
  lav_warning_nonexported("lav_mvnorm_missing_loglik_samplestats")
  lav_mvn_mi_loglik_samp(
    yp = Yp,
    mu = Mu,
    sigma_1 = Sigma,
    x_idx = x.idx,
    x_mean = x.mean,
    x_cov = x.cov,
    sinv_method = Sinv.method,
    log2pi = log2pi,
    minus_two = minus.two
  )
}

lav_mvnorm_cluster_implied22l <- function(
  lp = NULL,
  implied = NULL,
  Mu.W = NULL,
  Mu.B = NULL,
  Sigma.W = NULL,
  Sigma.B = NULL
) {
  lav_warning_nonexported("lav_mvnorm_cluster_implied22l")
  lav_mvn_cl_implied22l(
    lp = lp,
    implied = implied,
    mu_w = Mu.W,
    mu_b = Mu.B,
    sigma_w = Sigma.W,
    sigma_b = Sigma.B
  )
}

lav_matrix_symmetric_inverse_update <- function(
  S.inv,
  rm.idx = integer(0L),
  logdet = FALSE,
  S.logdet = NULL
) {
  lav_warning_nonexported("lav_matrix_symmetric_inverse_update")
  lav_mat_sym_inverse_update(
    s_inv = S.inv,
    rm_idx = rm.idx,
    logdet = logdet,
    s_logdet = S.logdet
  )
}

lav_mvnorm_dmvnorm <- function(Y = NULL,
                               wt = NULL,
                               Mu = NULL,
                               Sigma = NULL,
                               Sigma.inv = NULL,
                               Sinv.method = "eigen",
                               x.idx = integer(0L),
                               x.mean = NULL,
                               x.cov = NULL,
                               log = TRUE) {
  lav_warning_nonexported("lav_matrix_symmetric_inverse_update")
  lav_mvn_dmvnorm(
    y = Y,
    wt = wt,
    mu = Mu,
    sigma_1 = Sigma,
    sigma_inv = Sigma.inv,
    sinv_method = Sinv.method,
    x_idx = x.idx,
    x_mean = x.mean,
    x_cov = x.cov,
    log = log
  )
}

lav_partable_subset_measurement_model <- function(
  PT = NULL,
  lv.names = NULL,
  add.lv.cov = TRUE,
  add.ind.predictors = FALSE,
  add.idx = FALSE,
  idx.only = FALSE
) {
  lav_warning_nonexported("lav_partable_subset_measurement_model")
  lav_pt_subset_mm(
    pt_1 = PT,
    lv_names = lv.names,
    add_lv_cov = add.lv.cov,
    add_ind_predictors = add.ind.predictors,
    add_idx = add.idx,
    idx_only = idx.only
  )
}

lav_partable_subset_structural_model <- function(
  PT = NULL,
  add.idx = FALSE,
  idx.only = FALSE,
  add.exo.cov = FALSE,
  fixed.x = FALSE,
  conditional.x = FALSE,
  free.fixed.var = FALSE,
  meanstructure = FALSE
) {
  lav_warning_nonexported("lav_partable_subset_structural_model")
  lav_pt_subset_sm(
    pt_1 = PT,
    add_idx = add.idx,
    idx_only = idx.only,
    add_exo_cov = add.exo.cov,
    fixed_x = fixed.x,
    conditional_x = conditional.x,
    free_fixed_var = free.fixed.var,
    meanstructure = meanstructure)
}
