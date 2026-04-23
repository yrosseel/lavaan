# All code below is written by Myrsini Katsikatsou (Feb 2015)

# The following function refers to PLRT for nested models and equality
# constraints. Namely, it is developed to test either of the following
# hypotheses:
# a) H0 states that some parameters are equal to 0
# b) H0 states that some parameters are equal to some others.
# Note that for the latter I haven't checked if it is ok when equality
# constraints are imposed on parameters that refer to different groups in a
# multi-group analysis.
# All the code below has been developed for a single-group analysis.

# Let fit_objH0 and fit_objH1 be the outputs of lavaan() function when we fit
# a model under the null hypothesis and under the alternative, respectively.
# The argument equal_constr is logical (T/F) and it is TRUE if  equality
# constraints are imposed on subsets of the parameters.

# The main idea of the code below is that we consider the parameter vector
# under the alternative H1 evaluated at the values derived under H0 and for
# these values we should evaluate the Hessian, the variability matrix (denoted
# by J) and Godambe matrix.

lav_pml_test_plrt <- function(fit_obj_h0, fit_obj_h1) {
  # sanity check, perhaps we misordered H0 and H1 in the function call??
  if (fit_obj_h1@test[[1]]$df > fit_obj_h0@test[[1]]$df) {
    tmp <- fit_obj_h0
    fit_obj_h0 <- fit_obj_h1
    fit_obj_h1 <- tmp
  }

  plrt <- 2 * (fit_obj_h1@optim$logl - fit_obj_h0@optim$logl)

  # create a new object 'objH1_h0': the object 'H1', but where
  # the parameter values are from H0
  obj_h1_h0 <- lav_test_diff_m10(m1 = fit_obj_h1, m0 = fit_obj_h0, test = FALSE)

  # EqMat # YR: from 0.6-2, use lav_test_diff_a() (again)
  #             this should allow us to test models that are
  #             nested in the covariance matrix sense, but not
  #             in the parameter (table) sense
  eq_mat <- lav_test_diff_a(m1 = fit_obj_h1, m0 = fit_obj_h0)
  if (obj_h1_h0@Model@eq.constraints) {
    eq_mat <- eq_mat %*% t(obj_h1_h0@Model@eq.constraints.K)
  }
  # if (equal_constr == TRUE) {
  #    EqMat <- fit_objH0@Model@ceq.JAC
  # } else {
  #    PT0 <- fit_objH0@ParTable
  #    PT1 <- fit_objH1@ParTable
  #    h0.par.idx <- which(PT1$free > 0 & !(PT0$free > 0))
  #    tmp.ind <- PT1$free[ h0.par.idx ]
  #
  #    no.par0 <- length(tmp.ind)
  #    tmp.ind2 <- cbind(1:no.par0, tmp.ind ) # matrix indices
  #    EqMat <- matrix(0, nrow=no.par0, ncol=fit_objH1@Model@nx.free)
  #    EqMat[tmp.ind2] <- 1
  # }

  # DEBUG YR -- eliminate the constraints also present in H1
  #          -- if we do this, there is no need to use MASS::ginv later
  # JAC0 <- fit_objH0@Model@ceq.JAC
  # JAC1 <- fit_objH1@Model@ceq.JAC
  # unique.idx <- which(apply(JAC0, 1, function(x) {
  #                    !any(apply(JAC1, 1, function(y) { all(x == y) })) }))
  # if(length(unique.idx) > 0L) {
  #    EqMat <- EqMat[unique.idx,,drop = FALSE]
  # }

  # Observed information (= for PML, this is Hessian / N)
  hes_theta0 <- lavTech(obj_h1_h0, "information.observed")

  # handle possible constraints in H1 (and therefore also in objH1_h0)
  inv_hes_theta0 <-
    lav_model_information_augment_invert(
      lavmodel = obj_h1_h0@Model,
      information = hes_theta0,
      inverted = TRUE
    )

  # the estimated variability matrix is given (=unit information first order)
  j_theta0 <- lavTech(obj_h1_h0, "first.order")

  # the Inverse of the G matrix
  inv_g <- inv_hes_theta0 %*% j_theta0 %*% inv_hes_theta0

  minv_gt_m <- eq_mat %*% inv_g %*% t(eq_mat)
  minv_ht_m <- eq_mat %*% inv_hes_theta0 %*% t(eq_mat)
  # Inv_MinvHtM <- solve(MinvHtM)
  inv_minv_ht_m <- MASS::ginv(minv_ht_m)
  tmp_prod <- minv_gt_m %*% inv_minv_ht_m
  tmp_prod2 <- tmp_prod %*% tmp_prod
  sum_eig <- sum(diag(tmp_prod))
  sum_eigsq <- sum(diag(tmp_prod2))

  fsma_plrt <- (sum_eig / sum_eigsq) * plrt
  adj_df <- (sum_eig * sum_eig) / sum_eigsq
  pvalue <- 1 - pchisq(fsma_plrt, df = adj_df)

  list(FSMA.PLRT = fsma_plrt, adj.df = adj_df, pvalue = pvalue)
}


# for testing: this is the 'original' (using m.el.idx and x.el.idx)
lav_pml_test_plrt2 <- function(fit_obj_h0, fit_obj_h1) {
  if (fit_obj_h1@test[[1]]$df > fit_obj_h0@test[[1]]$df) {
    tmp <- fit_obj_h0
    fit_obj_h0 <- fit_obj_h1
    fit_obj_h1 <- tmp
  }

  if (fit_obj_h0@Model@eq.constraints) {
    equal_constr <- TRUE
  } else {
    equal_constr <- FALSE
  }

  nsize <- fit_obj_h0@SampleStats@ntotal
  plrt <- 2 * nsize * (fit_obj_h0@optim$fx - fit_obj_h1@optim$fx)
  npar <- fit_obj_h1@optim$npar
  my_m_el_idx2 <- fit_obj_h1@Model@m.free.idx
  my_x_el_idx2 <- fit_obj_h1@Model@x.free.idx
  my_m_el_idx <- my_m_el_idx2
  my_x_el_idx <- my_x_el_idx2

  # my_m_el_idx2 <- fit_objH1@Model@m.free.idx
  # my_m_el_idx2 gives the POSITION index of the free parameters within each
  # parameter matrix under H1 model.
  # The index numbering restarts from 1 when we move to a new parameter matrix.
  # Within each matrix the index numbering "moves" columnwise.

  # my_x_el_idx2 <- fit_objH1@Model@x.free.idx
  # my_x_el_idx2 ENUMERATES the free parameters within each parameter matrix.
  # The numbering continues as we move from one parameter matrix to the
  # next one.

  # In the case of the symmetric matrices, Theta and Psi,in some functions below
  # we need to give as input my_m_el_idx2 and my_x_el_idx2 after
  # we have eliminated the information about the redundant parameters
  # (those placed above the main diagonal).
  # That's why I do the following:

  # my_m_el_idx <- my_m_el_idx2
  # my_x_el_idx <- my_x_el_idx2
  # Psi, the variance - covariance matrix of factors
  # if( length(my_x_el_idx2[[3]])!=0 & any(table(my_x_el_idx2[[3]])>1)) {
  #  nfac <- ncol(fit_objH1@Model@GLIST$lambda) #number of factors
  #  tmp  <- matrix(c(1:(nfac^2)), nrow= nfac, ncol= nfac )
  #  tmp_keep <- tmp[lower.tri(tmp, diag=TRUE)]
  #  my_m_el_idx[[3]] <- my_m_el_idx[[3]][my_m_el_idx[[3]] %in% tmp_keep]
  #  my_x_el_idx[[3]] <- unique( my_x_el_idx2[[3]] )
  # }

  # for Theta, the variance-covariance matrix of measurement errors
  # if( length(my_x_el_idx2[[2]])!=0 & any(table(my_x_el_idx2[[2]])>1)) {
  #  nvar <- fit_objH1@Model@nvar #number of indicators
  #  tmp  <- matrix(c(1:(nvar^2)), nrow= nvar, ncol= nvar )
  #  tmp_keep <- tmp[lower.tri(tmp, diag=TRUE)]
  #  my_m_el_idx[[2]] <- my_m_el_idx[[2]][my_m_el_idx[[2]] %in% tmp_keep]
  #  my_x_el_idx[[2]] <- unique( my_x_el_idx2[[2]] )
  # }

  # below the commands to find the row-column indices of the Hessian that
  # correspond to the parameters to be tested equal to 0
  # tmp.ind contains these indices
  # MY.m.el.idx2.H0 <- fit_objH0@Model@m.free.idx
  # tmp.ind <- c()
  # for(i in 1:6) {
  #   tmp.ind <- c(tmp.ind ,
  #               my_x_el_idx2[[i]] [!(my_m_el_idx2[[i]]  %in%
  #                                    MY.m.el.idx2.H0[[i]] )  ]  )
  # }
  # next line added by YR
  # tmp.ind <- unique(tmp.ind)

  # YR: use partable to find which parameters are restricted in H0
  #     (this should work in multiple groups too)
  # h0.par.idx <- which(   PT.H1.extended$free[PT.H1.extended$user < 2] > 0  &
  #                     !(PT.H0.extended$free[PT.H0.extended$user < 2] > 0)   )
  # tmp.ind <- PT.H1.extended$free[ h0.par.idx ]
  # print(tmp.ind)
  if (length(my_x_el_idx2[[3]]) != 0 && any(table(my_x_el_idx2[[3]]) > 1)) {
    nfac <- ncol(fit_obj_h1@Model@GLIST$lambda)
    tmp <- matrix(c(1:(nfac * nfac)), nrow = nfac, ncol = nfac)
    tmp_keep <- tmp[lower.tri(tmp, diag = TRUE)]
    my_m_el_idx[[3]] <- my_m_el_idx[[3]][my_m_el_idx[[3]] %in% tmp_keep]
    my_x_el_idx[[3]] <- unique(my_x_el_idx2[[3]])
  }

  if (length(my_x_el_idx2[[2]]) != 0 && any(table(my_x_el_idx2[[2]]) > 1)) {
    nvar <- fit_obj_h1@Model@nvar
    tmp <- matrix(c(1:(nvar * nvar)), nrow = nvar, ncol = nvar)
    tmp_keep <- tmp[lower.tri(tmp, diag = TRUE)]
    my_m_el_idx[[2]] <- my_m_el_idx[[2]][my_m_el_idx[[2]] %in% tmp_keep]
    my_x_el_idx[[2]] <- unique(my_x_el_idx2[[2]])
  }
  my_m_el_idx2_h0 <- fit_obj_h0@Model@m.free.idx

  tmp_ind <- c()
  for (i in 1:6) {
    tmp_ind <- c(tmp_ind, my_x_el_idx2[[i]][!(my_m_el_idx2[[i]] %in%
      my_m_el_idx2_h0[[i]])])
  }
  tmp_ind <- unique(tmp_ind)

  # if the models are nested because of equality constraints among the
  # parameters, we need to construct the matrix of derivatives of function
  # g(theta) with respect to theta where g(theta) is the function that
  # represents the equality constraints. g(theta) is
  # an rx1 vector where r are the equality constraints. In the null hypothesis
  # we test H0: g(theta)=0. The matrix of derivatives is of dimension:
  # nrows= number of free non-redundant parameters under H0, namely
  # NparH0 <- fit_objH0[[1]]@optim$npar , and ncols= number of free
  # non-redundant parameters under H1, namely
  #    NparH1 <- fit_objH0[[1]]@optim$npar.
  # The matrix of derivatives of g(theta) is composed of 0's, 1's, -1's, and in
  # the rows that refer to odd number of parameters that are equal there is
  # one -2.
  # The 1's, -1's (and possibly -2) are the contrast coefficients of the
  # parameters.
  # The sum of the rows should be equal to 0.
  # if(equal_constr==TRUE) {
  #    EqMat <- fit_objH0@Model@ceq.JAC
  # } else {
  #   no.par0 <- length(tmp.ind)
  #   tmp.ind2 <- cbind(1:no.par0, tmp.ind)
  #   EqMat <- matrix(0, nrow = no.par0, ncol = Npar)
  # EqMat[tmp.ind2] <- 1
  # }

  if (equal_constr == TRUE) {
    eq_mat <- fit_obj_h0@Model@ceq.JAC
  } else {
    no_par0 <- length(tmp_ind)
    tmp_ind2 <- cbind(1:no_par0, tmp_ind)
    eq_mat <- matrix(0, nrow = no_par0, ncol = npar)
    eq_mat[tmp_ind2] <- 1
  }

  obj <- fit_obj_h0

  # Compute the sum of the eigenvalues and the sum of the squared eigenvalues
  # so that the adjustment to PLRT can be applied.
  # Here a couple of functions (e.g. lav_pml_object_inspect_hessian) which are
  # modifications of lavaan functions (e.g. getHessian) are needed. These are
  # defined in the end of the file.

  # the quantity below follows the same logic as getHessian of lavaan 0.5-18
  # and it actually gives N*Hessian. That's why the command following the
  # command below.
  # NHes.theta0 <- lav_pml_object_inspect_hessian (object = obj@Model,
  #                           samplestats = obj@SampleStats ,
  #                           x = obj@Data@X ,
  #                           estimator = "PML",
  #                           lavcache = obj@Cache,
  #                           my_m_el_idx = my_m_el_idx,
  #                           my_x_el_idx = my_x_el_idx,
  #                           my_m_el_idx2 = my_m_el_idx2,
  #                              # input for lav_pml_object_x2glist
  #                           my_x_el_idx2 = my_x_el_idx2,
  #                              # input for lav_pml_object_x2glist
  #                           Npar = Npar,
  #                           equal_constr=equal_constr)
  nhes_theta0 <- lav_pml_object_inspect_hessian(
    object = obj@Model, samplestats = obj@SampleStats,
    x = obj@Data@X, estimator = "PML", lavcache = obj@Cache,
    my_m_el_idx = my_m_el_idx, my_x_el_idx = my_x_el_idx,
    my_m_el_idx2 = my_m_el_idx2, my_x_el_idx2 = my_x_el_idx2,
    npar = npar, equal_constr = equal_constr
  )
  hes_theta0 <- nhes_theta0 / nsize
  # Inv.Hes.theta0 <- solve(Hes.theta0)
  inv_hes_theta0 <- MASS::ginv(hes_theta0)

  nj_theta0 <- lav_pml_object_information_firstorder(
    object = obj, my_m_el_idx = my_m_el_idx,
    my_x_el_idx = my_x_el_idx, equal_constr = equal_constr
  )
  j_theta0 <- nj_theta0 / (nsize * nsize)


  inv_g <- inv_hes_theta0 %*% j_theta0 %*% inv_hes_theta0
  minv_gt_m <- eq_mat %*% inv_g %*% t(eq_mat)
  minv_ht_m <- eq_mat %*% inv_hes_theta0 %*% t(eq_mat)
  # Inv_MinvHtM <- solve(MinvHtM)    #!!! change names
  inv_minv_ht_m <- MASS::ginv(minv_ht_m)
  tmp_prod <- minv_gt_m %*% inv_minv_ht_m # !!! change names
  tmp_prod2 <- tmp_prod %*% tmp_prod
  sum_eig <- sum(diag(tmp_prod))
  sum_eigsq <- sum(diag(tmp_prod2))


  fsma_plrt <- (sum_eig / sum_eigsq) * plrt
  adj_df <- (sum_eig * sum_eig) / sum_eigsq
  pvalue <- 1 - pchisq(fsma_plrt, df = adj_df)
  list(FSMA.PLRT = fsma_plrt, adj.df = adj_df, pvalue = pvalue)
}




###############################################################################
# auxiliary functions used above, they are all copy from the corresponding
# functions of lavaan where parts no needed were deleted and some parts were
# modified.
# I mark the modifications with comments.


# library(lavaan)

# To run an example for the functions below the following input is needed.
# obj <- fit.objH0[[i]]
# object <- obj@Model
# samplestats = obj@SampleStats
# X = obj@Data@X
# estimator = "PML"
# lavcache = obj@Cache
# my_m_el_idx = my_m_el_idx
# my_x_el_idx = my_x_el_idx
# my_m_el_idx2 = my_m_el_idx2 # input for lav_pml_object_x2glist
# my_x_el_idx2 = my_x_el_idx2 # input for lav_pml_object_x2glist
# Npar = Npar
# equal_constr =TRUE

lav_pml_object_inspect_hessian <- function(object, samplestats, x,
                         estimator = "PML", lavcache,
                         my_m_el_idx, my_x_el_idx,
                         my_m_el_idx2, my_x_el_idx2,
                                 # input for lav_pml_object_x2glist
                         npar, # Npar is the number of parameters under H1
                         equal_constr) { # takes TRUE/ FALSE
  if (equal_constr) { # !!! added line
  }
  hessian <- matrix(0, npar, npar) #

  # !!!! MYfunction below
  x_1 <- lav_pml_object_inspect_parameters(
    object = object,
    glist = NULL, n = npar, # N the number of parameters to consider
    my_m_el_idx = my_m_el_idx,
    my_x_el_idx = my_x_el_idx
  )

  for (j in 1:npar) {
    h_j <- 1e-05
    x_left <- x_left2 <- x_right <- x_right2 <- x_1
    x_left[j] <- x_1[j] - h_j
    x_left2[j] <- x_1[j] - 2 * h_j
    x_right[j] <- x_1[j] + h_j
    x_right2[j] <- x_1[j] + 2 * h_j
    # !!!! MYfunction below : lav_pml_object_inspect_gradient and
    #                         lav_pml_object_x2glist
    g_left <- lav_pml_object_inspect_gradient(
      object = object,
      glist = lav_pml_object_x2glist(
        object = object, x = x_left,
        my_m_el_idx = my_m_el_idx2,
        my_x_el_idx = my_x_el_idx2
      ),
      samplestats = samplestats, x = x,
      lavcache = lavcache, estimator = "PML",
      my_m_el_idx = my_m_el_idx,
      my_x_el_idx = my_x_el_idx,
      equal_constr = equal_constr
    )

    g_left2 <- lav_pml_object_inspect_gradient(
      object = object,
      glist = lav_pml_object_x2glist(
        object = object, x = x_left2,
        my_m_el_idx = my_m_el_idx2,
        my_x_el_idx = my_x_el_idx2
      ),
      samplestats = samplestats, x = x,
      lavcache = lavcache, estimator = "PML",
      my_m_el_idx = my_m_el_idx,
      my_x_el_idx = my_x_el_idx,
      equal_constr = equal_constr
    )

    g_right <- lav_pml_object_inspect_gradient(
      object = object,
      glist = lav_pml_object_x2glist(
        object = object, x = x_right,
        my_m_el_idx = my_m_el_idx2,
        my_x_el_idx = my_x_el_idx2
      ),
      samplestats = samplestats, x = x,
      lavcache = lavcache, estimator = "PML",
      my_m_el_idx = my_m_el_idx,
      my_x_el_idx = my_x_el_idx,
      equal_constr = equal_constr
    )

    g_right2 <- lav_pml_object_inspect_gradient(
      object = object,
      glist = lav_pml_object_x2glist(
        object = object, x = x_right2,
        my_m_el_idx = my_m_el_idx2,
        my_x_el_idx = my_x_el_idx2
      ),
      samplestats = samplestats, x = x,
      lavcache = lavcache, estimator = "PML",
      my_m_el_idx = my_m_el_idx,
      my_x_el_idx = my_x_el_idx,
      equal_constr = equal_constr
    )

    hessian[, j] <- (g_left2 - 8 * g_left + 8 * g_right - g_right2) / (12 * h_j)
  }
  hessian <- (hessian + t(hessian)) / 2
  # (-1) * Hessian
  hessian
}
#############################################################################




##################################  lav_pml_object_inspect_parameters
# different input arguments: my_m_el_idx, my_x_el_idx
lav_pml_object_inspect_parameters <-         # nolint
  function(object, glist = NULL, n, # N the number of parameters to consider
           my_m_el_idx, my_x_el_idx) {
  if (is.null(glist)) {
    glist <- object@GLIST
  }

  x <- numeric(n)

  for (mm in seq_along(object@GLIST)) { # mm<-1
    m_idx <- my_m_el_idx[[mm]] # !!!!! different here and below
    x_idx <- my_x_el_idx[[mm]]
    x[x_idx] <- glist[[mm]][m_idx]
  }
  x
}
#############################################################################




#############################  lav_pml_object_inspect_gradient
# the difference are the input arguments my_m_el_idx, my_x_el_idx
# used  in  lavaan:::lav_model_delta
lav_pml_object_inspect_gradient <-                  # nolint
  function(object, glist, samplestats = NULL, x = NULL,
           lavcache = NULL, estimator = "PML",
           my_m_el_idx, my_x_el_idx, equal_constr) {
  if (equal_constr) { # added line
  }
  num_idx <- object@num.idx
  th_idx <- object@th.idx
  if (is.null(glist)) {
    glist <- object@GLIST
  }
  sigma_hat <-
    lav_model_sigma(object, GLIST = glist, extra = (estimator == "ML"))
  mu_hat <- lav_model_mu(object, GLIST = glist)
  th <- lav_model_th(object, GLIST = glist)
  g <- 1
  d1 <- lav_pml_dploglik_dimplied(
    Sigma.hat = sigma_hat[[g]], Mu.hat = mu_hat[[g]],
    TH = th[[g]], th.idx = th_idx[[g]], num.idx = num_idx[[g]],
    X = x[[g]], lavcache = lavcache[[g]]
  )

  # !?  if(equal_constr) { #delete the following three commented lines, wrong
  #     Delta <- lavaan:::lav_model_delta (lavmodel= object, GLIST. = GLIST)
  #  } else {
  delta <- lav_model_delta(
    lavmodel = object, GLIST. = glist,
    m.el.idx. = my_m_el_idx,
    x.el.idx. = my_x_el_idx
  )
  # }

  # !!! that was before: as.numeric(t(d1) %*% Delta[[g]])/samplestats@nobs[[g]]
  as.numeric(t(d1) %*% delta[[g]])
           # !!! modified to follow current computeGradient() function of lavaan
  # !!! which gives minus the gradient of PL-loglik
  }

###############################################################################


##################################  lav_pml_object_x2glist
# difference in input arguments my_m_el_idx, my_x_el_idx

lav_pml_object_x2glist <- function(object, x = NULL, my_m_el_idx, my_x_el_idx) {
  glist <- object@GLIST
  for (mm in seq_along(glist)) {
    m_el_idx <- my_m_el_idx[[mm]]
    x_el_idx <- my_x_el_idx[[mm]]
    glist[[mm]][m_el_idx] <- x[x_el_idx]
  }
  glist
}
############################################################################


##### lav_pml_object_information_firstorder function
# difference from corresponding of lavaan: I use lav_pml_model_vcov_firstorder
lav_pml_object_information_firstorder <-    # nolint
  function(object, my_m_el_idx, my_x_el_idx, equal_constr) {
  nacov <- lav_pml_model_vcov_firstorder(
    lavmodel = object@Model,
    lavsamplestats = object@SampleStats,
    lavdata = object@Data,
    estimator = "PML",
    my_m_el_idx = my_m_el_idx,
    my_x_el_idx = my_x_el_idx,
    equal_constr = equal_constr
  )
  if (equal_constr) { # added lines
  }
  b0 <- attr(nacov, "B0")
  # !!!! Note below that I don't multiply  with nsize
  # !!! so what I get is J matrix divided by n
  # if (object@Options$estimator == "PML") {
  #    B0 <- B0 * object@SampleStats@ntotal
  # }
  # !!!!!!!!!!!!!!!!!!! added the following lines so that the output of
  # !!! lav_pml_object_information_firstorder is in line with that of
  #                    lavaan 0.5-18 getVariability
  # !! what's the purpose of the following lines?
  if (object@Options$estimator == "PML") {
    b0 <- b0 * object@SampleStats@ntotal
  }

  b0
  }

##############################################################################
# example
# obj <- fit.objH0[[i]]
# object <- obj@Model
# samplestats = obj@SampleStats
# X = obj@Data@X
# estimator = "PML"
# lavcache = obj@Cache
# my_m_el_idx = my_m_el_idx
# my_x_el_idx = my_x_el_idx
# my_m_el_idx2 = my_m_el_idx2 # input for lav_pml_object_x2glist
# my_x_el_idx2 = my_x_el_idx2 # input for lav_pml_object_x2glist
# Npar = Npar
# equal_constr =TRUE


lav_pml_model_vcov_firstorder <- function(lavmodel, lavsamplestats = NULL,
                                lavdata = NULL, lavcache = NULL,
                                estimator = "PML",
                                my_m_el_idx, my_x_el_idx,
                                equal_constr) { # equal_constr takes TRUE/FALSE
  if (equal_constr) { # added lines
  }
  b0_group <- vector("list",
                     lavsamplestats@ngroups) # in my case list of length 1

  # !?   if (equal_constr) {
  #      ###the following three lines are commented because they are wrong
  #      Delta <- lavaan:::lav_model_delta(lavmodel, GLIST. = NULL)
  #   } else {
  delta <- lav_model_delta(lavmodel,
    GLIST. = NULL,
    m.el.idx. = my_m_el_idx, # !!!!! different here and below
    x.el.idx. = my_x_el_idx
  )
  #  }
  sigma_hat <- lav_model_sigma(lavmodel)
  mu_hat <- lav_model_mu(lavmodel)
  th <- lav_model_th(lavmodel)
  g <- 1

  sc <- lav_pml_dploglik_dimplied(
    Sigma.hat = sigma_hat[[g]], TH = th[[g]],
    Mu.hat = mu_hat[[g]], th.idx = lavmodel@th.idx[[g]],
    num.idx = lavmodel@num.idx[[g]],
    X = lavdata@X[[g]], lavcache = lavcache,
    scores = TRUE, negative = FALSE
  )
  group_sc <- sc %*% delta[[g]]
  b0_group[[g]] <- lav_matrix_crossprod(group_sc)
  # !!!! B0.group[[g]] <- B0.group[[g]]/lavsamplestats@ntotal
  #  !!! skip so that the result is in line with the 0.5-18 version of lavaan

  b0 <- b0_group[[1]]

  m_e <- b0

  eigvals <- eigen(m_e, symmetric = TRUE, only.values = TRUE)$values
  if (any(eigvals < -1 * .Machine$double.eps^(3 / 4))) {
    lav_msg_warn(gettext(
      "matrix based on first order outer product of the derivatives is not
      positive definite; the standard errors may not be thrustworthy"))
  }
  nvar_cov <- MASS::ginv(m_e)

  attr(nvar_cov, "B0") <- b0
  attr(nvar_cov, "B0.group") <- b0_group
  nvar_cov
}
