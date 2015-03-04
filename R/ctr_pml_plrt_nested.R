# All code below is written by Myrsini Katsikatsou (Feb 2015)

#The following function refers to PLRT for nested models and equality constraints.
# Namely, it is developed to test either of the following hypotheses:
# a) H0 states that some parameters are equal to 0
# b) H0 states that some parameters are equal to some others.
#Note that for the latter I haven't checked if it is ok when equality constraints
#are imposed on parameters that refer to different groups in a multi-group 
#analysis. All the code below has been developed for a single-group analysis.

# Let fit_objH0 and fit_objH1 be the outputs of lavaan() function when we fit
# a model under the null hypothesis and under the alternative, respectively.
# The argument equalConstr is logical (T/F) and it is TRUE if  equality constraints
# are imposed on subsets of the parameters. 

# The main idea of the code below is that we consider the parameter vector
# under the alternative H1 evaluated at the values derived under H0 and for these
# values we should evaluate the Hessian, the variability matrix (denoted by J)
# and Godambe matrix.
  

ctr_pml_plrt_nested <- function(fit_objH0, fit_objH1) {

    stopifnot(fit_objH1@Data@ngroups == 1L)


    # sanity check, perhaps we misordered H0 and H1 in the function call??
    if(fit_objH1@Fit@test[[1]]$df > fit_objH0@Fit@test[[1]]$df) {
        tmp <- fit_objH0
        fit_objH0 <- fit_objH1
        fit_objH1 <- tmp
    }

    # check if we have equality constraints
    if(fit_objH0@Model@eq.constraints) {
        equalConstr = TRUE
    } else {
        equalConstr = FALSE
    }
  
    nsize <- fit_objH0@SampleStats@ntotal
    PLRT <- 2 * (fit_objH1@Fit@logl - fit_objH0@Fit@logl)

    # create a new object 'objH1_h0': the object 'H1', but where
    # the parameter values are from H0
    #objH1_h0 <- update(fit_objH1, 
    #                   start = parameterEstimates(fit0), 
    #                   control = list(optim.method="none", 
    #                                  optim.force.converged = TRUE), 
    #                   se = "none", test = "none", verbose = FALSE)

    Options <- fit_objH1@Options
    Options$verbose <- FALSE; Options$se <- "none"; Options$test <- "none"

    PT.H0 <- fit_objH0@ParTable
    PT.H1 <- fit_objH1@ParTable

    # `extend' PT.H1 partable to include all `fixed-to-zero parameters'
    PT.H1.FULL <- lav_partable_full(PT.H1, free = TRUE, start = TRUE)
    PT.H1.extended <- lav_partable_merge(PT.H1, PT.H1.FULL, 
                                         remove = TRUE, warn = FALSE)

    # `extend' PT.H0 partable to include all `fixed-to-zero parameters'
    PT.H0.FULL <- lav_partable_full(PT.H0, free = TRUE, start = TRUE)
    PT.H0.extended <- lav_partable_merge(PT.H0, PT.H0.FULL,
                                         remove = TRUE, warn = FALSE)

    # `extend' PE of H0 to include all `fixed-to-zero parameters'
    PE.H0 <- parameterEstimates(fit_objH0)
    PE.H0.FULL <- lav_partable_full(PE.H0)
    PE.H0.extended <- lav_partable_merge(PE.H0, PE.H0.FULL,
                                         remove = TRUE, warn = FALSE)

    

    objH1_h0 <- lavaan(model = PT.H1.extended,
                       start = PE.H0.extended,
                       control=list(optim.method          = "none",
                                    optim.force.converged = TRUE) ,
                       slotOptions     = Options,
                       slotSampleStats = fit_objH1@SampleStats,
                       slotData        = fit_objH1@Data, 
                       slotCache       = fit_objH1@Cache)

  
    # number of parameters under the alternative hypothesis H1
    Npar <- fit_objH1@Fit@npar  


   #MY.m.el.idx2 <- fit_objH1@Model@m.free.idx
   # MY.m.el.idx2 gives the POSITION index of the free parameters within each 
   # parameter matrix under H1 model.
   # The index numbering restarts from 1 when we move to a new parameter matrix.
   # Within each matrix the index numbering "moves" columnwise.

   #MY.x.el.idx2 <- fit_objH1@Model@x.free.idx
   # MY.x.el.idx2 ENUMERATES the free parameters within each parameter matrix.
   # The numbering continues as we move from one parameter matrix to the next one.

   # In the case of the symmetric matrices, Theta and Psi,in some functions below 
   # we need to give as input MY.m.el.idx2 and MY.x.el.idx2 after
   # we have eliminated the information about the redundant parameters
   # (those placed above the main diagonal).
   # That's why I do the following:

   #MY.m.el.idx <- MY.m.el.idx2
   #MY.x.el.idx <- MY.x.el.idx2
   # Psi, the variance - covariance matrix of factors
   #if( length(MY.x.el.idx2[[3]])!=0 & any(table(MY.x.el.idx2[[3]])>1)) {
   #  nfac <- ncol(fit_objH1@Model@GLIST$lambda) #number of factors
   #  tmp  <- matrix(c(1:(nfac^2)), nrow= nfac, ncol= nfac )
   #  tmp_keep <- tmp[lower.tri(tmp, diag=TRUE)]
   #  MY.m.el.idx[[3]] <- MY.m.el.idx[[3]][MY.m.el.idx[[3]] %in% tmp_keep]
   #  MY.x.el.idx[[3]] <- unique( MY.x.el.idx2[[3]] )
   #}

   #for Theta, the variance-covariance matrix of measurement errors
   # if( length(MY.x.el.idx2[[2]])!=0 & any(table(MY.x.el.idx2[[2]])>1)) {
   #  nvar <- fit_objH1@Model@nvar #number of indicators
   #  tmp  <- matrix(c(1:(nvar^2)), nrow= nvar, ncol= nvar )
   #  tmp_keep <- tmp[lower.tri(tmp, diag=TRUE)]
   #  MY.m.el.idx[[2]] <- MY.m.el.idx[[2]][MY.m.el.idx[[2]] %in% tmp_keep]
   #  MY.x.el.idx[[2]] <- unique( MY.x.el.idx2[[2]] )
   # }
 
   #below the commands to find the row-column indices of the Hessian that correspond to
   #the parameters to be tested equal to 0
   #tmp.ind contains these indices
   # MY.m.el.idx2.H0 <- fit_objH0@Model@m.free.idx
   # tmp.ind <- c()
   # for(i in 1:6) {
   #   tmp.ind <- c(tmp.ind ,
   #               MY.x.el.idx2[[i]] [!(MY.m.el.idx2[[i]]  %in%
   #                                    MY.m.el.idx2.H0[[i]] )  ]  )
   # }
   # next line added by YR
   # tmp.ind <- unique(tmp.ind)

   # YR: use partable to find which parameters are restricted in H0
   #     (this should work in multiple groups too)
   #h0.par.idx <- which(   PT.H1.extended$free[PT.H1.extended$user < 2] > 0  & 
   #                     !(PT.H0.extended$free[PT.H0.extended$user < 2] > 0)   )
   #tmp.ind <- PT.H1.extended$free[ h0.par.idx ]
   #print(tmp.ind)


 # if the models are nested because of equality constraints among the parameters, we need
 # to construct the matrix of derivatives of function g(theta) with respect to theta
 # where g(theta) is the function that represents the equality constraints. g(theta) is
 # an rx1 vector where r are the equality constraints. In the null hypothesis
 # we test H0: g(theta)=0. The matrix of derivatives is of dimension:
 # nrows= number of free non-redundant parameters under H0, namely 
 # NparH0 <- fit_objH0[[1]]@Fit@npar , and ncols= number of free non-redundant
 # parameters under H1, namely NparH1 <- fit_objH0[[1]]@Fit@npar.
 # The matrix of derivatives of g(theta) is composed of 0's, 1's, -1's, and 
 # in the rows that refer to odd number of parameters that are equal there is one -2.
 # The 1's, -1's (and possibly -2) are the contrast coefficients of the parameters.
 # The sum of the rows should be equal to 0.
 #if(equalConstr==TRUE) {
 #    EqMat <- fit_objH0@Model@ceq.JAC
 #} else {
 #   no.par0 <- length(tmp.ind)
 #   tmp.ind2 <- cbind(1:no.par0, tmp.ind)
 #   EqMat <- matrix(0, nrow = no.par0, ncol = Npar)
 # EqMat[tmp.ind2] <- 1
# }

  # use a numeric approximation of `EqMat'
  Delta1 <- computeDelta(fit_objH1@Model)[[1]] # FIXME multiple groups???
  Delta0 <- computeDelta(fit_objH0@Model)[[1]]
  H <- solve(t(Delta1) %*% Delta1) %*% t(Delta1) %*% Delta0
  EqMat <- t(lav_matrix_orthogonal_complement(H))


 # Compute the sum of the eigenvalues and the sum of the squared eigenvalues
 # so that the adjustment to PLRT can be applied.
 # Here a couple of functions (e.g. MYgetHessian) which are modifications of 
 # lavaan functions (e.g. getHessian) are needed. These are defined in the end of the file.  
 # obj <- fit_objH0

 #the quantity below follows the same logic as getHessian of lavaan 0.5-18
 #and it actually gives N*Hessian. That's why the command following the command below. 
 # NHes.theta0 <- MYgetHessian (object = obj@Model,
 #                           samplestats = obj@SampleStats ,
 #                           X = obj@Data@X ,
 #                           estimator = "PML",
 #                           lavcache = obj@Cache,
 #                           MY.m.el.idx = MY.m.el.idx,
 #                           MY.x.el.idx = MY.x.el.idx,
 #                           MY.m.el.idx2 = MY.m.el.idx2, # input for MYx2GLIST
 #                           MY.x.el.idx2 = MY.x.el.idx2, # input for MYx2GLIST
 #                           Npar = Npar,
 #                           equalConstr=equalConstr)
 NHes.theta0 <- lavTech(objH1_h0, "hessian")
 Hes.theta0 <- NHes.theta0/ nsize #!!!!! to get the Hessian
 Inv.Hes.theta0 <- solve(Hes.theta0)

 #N times the estimated variability matrix is given :
 #NJ.theta0 <- MYgetVariability(object = obj,
 #                            MY.m.el.idx = MY.m.el.idx,
 #                            MY.x.el.idx = MY.x.el.idx,
 #                            equalConstr=equalConstr)
 #J.theta0.divbyN <-  NJ.theta0/(nsize^2) #!!!!!! added due to changes of lavaan functions
 NJ.theta0 <- lavTech(objH1_h0, "first.order")
 # J.theta0 <- NJ.theta0/(nsize^2)
 J.theta0 <-  NJ.theta0 / nsize
 
 # below the Inverse of the G matrix divided by N
 Inv.G <- Inv.Hes.theta0 %*% J.theta0 %*% Inv.Hes.theta0

 MInvGtM <- EqMat %*% Inv.G %*% t(EqMat)
 MinvHtM <- EqMat %*% Inv.Hes.theta0 %*% t(EqMat)
 Inv_MinvHtM <- solve(MinvHtM)    #!!! change names
 tmp.prod <- MInvGtM %*% Inv_MinvHtM  #!!! change names
 tmp.prod2 <- tmp.prod %*% tmp.prod
 sum.eig <- sum(diag(tmp.prod))
 sum.eigsq <- sum(diag(tmp.prod2))

 FSMA.PLRT <- (sum.eig/sum.eigsq) * PLRT
 adj.df <- (sum.eig^2)/sum.eigsq
 pvalue <- 1 - pchisq(FSMA.PLRT, df = adj.df)
 list(FSMA.PLRT = FSMA.PLRT, adj.df = adj.df, pvalue = pvalue)
}


# for testing: this is the 'original' (using m.el.idx and x.el.idx)
ctr_pml_plrt_nested2 <- function (fit_objH0, fit_objH1) {

    if (fit_objH1@Fit@test[[1]]$df > fit_objH0@Fit@test[[1]]$df) {
        tmp <- fit_objH0
        fit_objH0 <- fit_objH1
        fit_objH1 <- tmp
    }

    if (fit_objH0@Model@eq.constraints) {
        equalConstr = TRUE
    }   else {
        equalConstr = FALSE
    }

    nsize <- fit_objH0@Data@nobs[[1]]
    PLRT <- 2 * nsize * (fit_objH0@Fit@fx - fit_objH1@Fit@fx)
    Npar <- fit_objH1@Fit@npar
    MY.m.el.idx2 <- fit_objH1@Model@m.free.idx
    MY.x.el.idx2 <- fit_objH1@Model@x.free.idx
    MY.m.el.idx <- MY.m.el.idx2
    MY.x.el.idx <- MY.x.el.idx2

    if (length(MY.x.el.idx2[[3]]) != 0 & any(table(MY.x.el.idx2[[3]]) > 1)) {
        nfac <- ncol(fit_objH1@Model@GLIST$lambda)
        tmp <- matrix(c(1:(nfac^2)), nrow = nfac, ncol = nfac)
        tmp_keep <- tmp[lower.tri(tmp, diag = TRUE)]
        MY.m.el.idx[[3]] <- MY.m.el.idx[[3]][MY.m.el.idx[[3]] %in% tmp_keep]
        MY.x.el.idx[[3]] <- unique(MY.x.el.idx2[[3]])
    }

    if (length(MY.x.el.idx2[[2]]) != 0 & any(table(MY.x.el.idx2[[2]]) > 1)) {
        nvar <- fit_objH1@Model@nvar
        tmp <- matrix(c(1:(nvar^2)), nrow = nvar, ncol = nvar)
        tmp_keep <- tmp[lower.tri(tmp, diag = TRUE)]
        MY.m.el.idx[[2]] <- MY.m.el.idx[[2]][MY.m.el.idx[[2]] %in% tmp_keep]
        MY.x.el.idx[[2]] <- unique(MY.x.el.idx2[[2]])
    }
    MY.m.el.idx2.H0 <- fit_objH0@Model@m.free.idx

    tmp.ind <- c()
    for (i in 1:6) {
        tmp.ind <- c(tmp.ind, MY.x.el.idx2[[i]][!(MY.m.el.idx2[[i]] %in%
            MY.m.el.idx2.H0[[i]])])
    }
    tmp.ind <- unique(tmp.ind)

    if (equalConstr == TRUE) {
        EqMat <- fit_objH0@Model@ceq.JAC
    } else {
        no.par0 <- length(tmp.ind)
        tmp.ind2 <- cbind(1:no.par0, tmp.ind )
        EqMat <- matrix(0, nrow=no.par0, ncol=Npar)
        EqMat[tmp.ind2] <- 1
    }

    obj <- fit_objH0

    NHes.theta0 <- MYgetHessian(object = obj@Model, samplestats = obj@SampleStats,
        X = obj@Data@X, estimator = "PML", lavcache = obj@Cache,
        MY.m.el.idx = MY.m.el.idx, MY.x.el.idx = MY.x.el.idx,
        MY.m.el.idx2 = MY.m.el.idx2, MY.x.el.idx2 = MY.x.el.idx2,
        Npar = Npar, equalConstr = equalConstr)
    Hes.theta0 <- NHes.theta0/nsize
    Inv.Hes.theta0 <- solve(Hes.theta0)

    NJ.theta0 <- MYgetVariability(object = obj, MY.m.el.idx = MY.m.el.idx,
        MY.x.el.idx = MY.x.el.idx, equalConstr = equalConstr)
    J.theta0 <- NJ.theta0/(nsize^2)


    Inv.G <- Inv.Hes.theta0 %*% J.theta0 %*% Inv.Hes.theta0
    MInvGtM <- EqMat %*% Inv.G %*% t(EqMat)
    MinvHtM <- EqMat %*% Inv.Hes.theta0 %*% t(EqMat)
    Inv_MinvHtM <- solve(MinvHtM)    #!!! change names
    tmp.prod <- MInvGtM %*% Inv_MinvHtM  #!!! change names
    tmp.prod2 <- tmp.prod %*% tmp.prod
    sum.eig <- sum(diag(tmp.prod))
    sum.eigsq <- sum(diag(tmp.prod2))


    FSMA.PLRT <- (sum.eig/sum.eigsq) * PLRT
    adj.df <- (sum.eig^2)/sum.eigsq
    pvalue <- 1 - pchisq(FSMA.PLRT, df = adj.df)
    list(FSMA.PLRT = FSMA.PLRT, adj.df = adj.df, pvalue = pvalue)
}




###################################################################################
# auxiliary functions used above, they are all copy from the corresponding functions
# of lavaan where parts no needed were deleted and some parts were modified.
# I mark the modifications with comments.


# library(lavaan)

# To run an example for the functions below the following input is needed.
# obj <- fit.objH0[[i]] 
# object <- obj@Model 
# samplestats = obj@SampleStats 
# X = obj@Data@X 
# estimator = "PML" 
# lavcache = obj@Cache
# MY.m.el.idx = MY.m.el.idx 
# MY.x.el.idx = MY.x.el.idx
# MY.m.el.idx2 = MY.m.el.idx2 # input for MYx2GLIST
# MY.x.el.idx2 = MY.x.el.idx2 # input for MYx2GLIST
# Npar = Npar
# equalConstr =TRUE

MYgetHessian <- function (object, samplestats , X ,
                           estimator = "PML", lavcache,
                           MY.m.el.idx, MY.x.el.idx,
                           MY.m.el.idx2, MY.x.el.idx2, # input for MYx2GLIST
                           Npar,     #Npar is the number of parameters under H1
                           equalConstr  ) { # takes TRUE/ FALSE  
    if(equalConstr){   #!!! added line
    }
    Hessian <- matrix(0, Npar, Npar) #

    #!!!! MYfunction below
    x <- MYgetModelParameters(object=object,
                              GLIST = NULL, N=Npar, #N the number of parameters to consider
                              MY.m.el.idx=MY.m.el.idx,
                              MY.x.el.idx= MY.x.el.idx)

    for (j in 1:Npar) {
        h.j <- 1e-05
        x.left <- x.left2 <- x.right <- x.right2 <- x
        x.left[j] <- x[j] - h.j
        x.left2[j] <- x[j] - 2 * h.j
        x.right[j] <- x[j] + h.j
        x.right2[j] <- x[j] + 2 * h.j
        #!!!! MYfunction below : MYcomputeGradient and  MYx2GLIST
        g.left <- MYcomputeGradient(object=object,
                                    GLIST = MYx2GLIST(object=object, x = x.left,
                                                      MY.m.el.idx=MY.m.el.idx2,
                                                      MY.x.el.idx= MY.x.el.idx2),
                                    samplestats = samplestats, X = X,
                                    lavcache = lavcache, estimator = "PML",
                                    MY.m.el.idx=MY.m.el.idx,
                                    MY.x.el.idx= MY.x.el.idx,
                                    equalConstr = equalConstr   )

        g.left2 <- MYcomputeGradient(object=object,
                                    GLIST = MYx2GLIST(object=object, x = x.left2,
                                                      MY.m.el.idx=MY.m.el.idx2,
                                                      MY.x.el.idx= MY.x.el.idx2),
                                    samplestats = samplestats, X = X,
                                    lavcache = lavcache, estimator = "PML",
                                    MY.m.el.idx=MY.m.el.idx,
                                    MY.x.el.idx= MY.x.el.idx,
                                    equalConstr = equalConstr  )

        g.right <- MYcomputeGradient(object=object,
                                    GLIST = MYx2GLIST(object=object, x = x.right,
                                                      MY.m.el.idx=MY.m.el.idx2,
                                                      MY.x.el.idx= MY.x.el.idx2),
                                    samplestats = samplestats, X = X,
                                    lavcache = lavcache, estimator = "PML",
                                    MY.m.el.idx=MY.m.el.idx,
                                    MY.x.el.idx= MY.x.el.idx,
                                    equalConstr = equalConstr  )

       g.right2 <- MYcomputeGradient(object=object,
                                    GLIST = MYx2GLIST(object=object, x = x.right2,
                                                      MY.m.el.idx=MY.m.el.idx2,
                                                      MY.x.el.idx= MY.x.el.idx2),
                                    samplestats = samplestats, X = X,
                                    lavcache = lavcache, estimator = "PML",
                                    MY.m.el.idx=MY.m.el.idx,
                                    MY.x.el.idx= MY.x.el.idx,
                                    equalConstr = equalConstr )

        Hessian[, j] <- (g.left2 - 8 * g.left + 8 * g.right - g.right2)/(12 * h.j)
    }
    Hessian <- (Hessian + t(Hessian))/2
    #(-1) * Hessian
    Hessian
}
#############################################################################




##################################  MYgetModelParameters
#different input arguments: MY.m.el.idx, MY.x.el.idx
MYgetModelParameters  <- function (object, GLIST = NULL, N, #N the number of parameters to consider
                                   MY.m.el.idx, MY.x.el.idx) {
    if (is.null(GLIST)) {
        GLIST <- object@GLIST
    }

    x <- numeric(N)

    for (mm in 1:length(object@GLIST)) {      # mm<-1
        m.idx <- MY.m.el.idx[[mm]]   #!!!!! different here and below
        x.idx <- MY.x.el.idx[[mm]]
        x[x.idx] <- GLIST[[mm]][m.idx]
    }
    x
}
#############################################################################




#############################  MYcomputeGradient
#the difference are the input arguments MY.m.el.idx, MY.x.el.idx
#used  in  lavaan:::computeDelta
MYcomputeGradient <- function (object, GLIST, samplestats = NULL, X = NULL,
                               lavcache = NULL, estimator = "PML", 
                               MY.m.el.idx, MY.x.el.idx, equalConstr  ) {
   if(equalConstr){  #added line
    }
   num.idx <- object@num.idx
   th.idx <- object@th.idx
   if (is.null(GLIST)) {
        GLIST <- object@GLIST
    }
   Sigma.hat <- lavaan:::computeSigmaHat(object, GLIST = GLIST, extra = (estimator ==  "ML"))
   TH <- lavaan:::computeTH(object, GLIST = GLIST)
   g<-1
   d1 <- lavaan:::pml_deriv1(Sigma.hat = Sigma.hat[[g]], TH = TH[[g]],
                             th.idx = th.idx[[g]], num.idx = num.idx[[g]],
                             X = X[[g]], lavcache = lavcache[[g]])

 #!?  if(equalConstr) { #delete the following three commented lines, wrong
 #     Delta <- lavaan:::computeDelta (lavmodel= object, GLIST. = GLIST)
 #  } else {
      Delta <- lavaan:::computeDelta (lavmodel= object, GLIST. = GLIST,
                                      m.el.idx. = MY.m.el.idx , 
                                      x.el.idx. = MY.x.el.idx)
  # }

 #!!!!! that was before: as.numeric(t(d1) %*% Delta[[g]])/samplestats@nobs[[g]]
 as.numeric(t(d1) %*% Delta[[g]]) #!!! modified to follow current computeGradient() function of lavaan
 #!!! which gives minus the gradient of PL-loglik
}

###############################################################################


##################################  MYx2GLIST
#difference in input arguments MY.m.el.idx, MY.x.el.idx

MYx2GLIST <- function (object, x = NULL, MY.m.el.idx, MY.x.el.idx) {
    GLIST <- object@GLIST
    for (mm in 1:length(GLIST)) {
         m.el.idx <- MY.m.el.idx[[mm]]
         x.el.idx <- MY.x.el.idx[[mm]]
         GLIST[[mm]][m.el.idx] <- x[x.el.idx]
    }
    GLIST
}
############################################################################


#####MYgetVariability function
#difference from corresponding of lavaan: I use MYNvcov.first.order
MYgetVariability <- function (object, MY.m.el.idx, MY.x.el.idx, equalConstr ) {
    NACOV <- MYNvcov.first.order(lavmodel=object@Model,
                                 lavsamplestats = object@SampleStats,
                                 lavdata = object@Data,
                                 estimator = "PML",
                                 MY.m.el.idx=MY.m.el.idx,
                                 MY.x.el.idx= MY.x.el.idx, 
                                 equalConstr = equalConstr)
    if(equalConstr){  #added lines
    }
    B0 <- attr(NACOV, "B0")
    #!!!! Note below that I don't multiply  with nsize
    #!!! so what I get is J matrix divided by n
    #if (object@Options$estimator == "PML") {
    #    B0 <- B0 * object@SampleStats@ntotal
    #}
    #!!!!!!!!!!!!!!!!!!! added the following lines so that the output of
    #!!!!! MYgetVariability is in line with that of lavaan 0.5-18 getVariability
    #!! what's the purpose of the following lines?
     if (object@Options$estimator == "PML") {
        B0 <- B0 * object@SampleStats@ntotal
    }

    B0
}

##############################################################################
# example
# obj <- fit.objH0[[i]] 
# object <- obj@Model 
# samplestats = obj@SampleStats 
# X = obj@Data@X 
# estimator = "PML" 
# lavcache = obj@Cache
# MY.m.el.idx = MY.m.el.idx 
# MY.x.el.idx = MY.x.el.idx
# MY.m.el.idx2 = MY.m.el.idx2 # input for MYx2GLIST
# MY.x.el.idx2 = MY.x.el.idx2 # input for MYx2GLIST
# Npar = Npar
# equalConstr =TRUE


MYNvcov.first.order <- function (lavmodel, lavsamplestats = NULL,
                                 lavdata = NULL, lavcache = NULL,
                                 estimator = "PML",  
                                 MY.m.el.idx, MY.x.el.idx,
                                 equalConstr ) {    #equalConstr takes TRUE/FALSE
    if(equalConstr){ #added lines
    }
    B0.group <- vector("list", lavsamplestats@ngroups)  #in my case list of length 1

 #!?   if (equalConstr) {     ###the following three lines are commented because they are wrong
 #      Delta <- lavaan:::computeDelta(lavmodel, GLIST. = NULL)
 #   } else {
       Delta <- lavaan:::computeDelta(lavmodel,
                                       GLIST. = NULL,
                                       m.el.idx. = MY.m.el.idx,#!!!!! different here and below
                                       x.el.idx. = MY.x.el.idx)
  #  }
    Sigma.hat <- lavaan:::computeSigmaHat(lavmodel)
    TH <- lavaan:::computeTH(lavmodel)
    g <-1

    SC <- lavaan:::pml_deriv1(Sigma.hat = Sigma.hat[[g]], TH = TH[[g]],
                              th.idx = lavmodel@th.idx[[g]], 
                              num.idx = lavmodel@num.idx[[g]],
                              X = lavdata@X[[g]], lavcache = lavcache,
                              scores = TRUE, negative = FALSE)
    group.SC <- SC %*% Delta[[g]]
    B0.group[[g]] <- crossprod(group.SC)
    #!!!! B0.group[[g]] <- B0.group[[g]]/lavsamplestats@ntotal  !!! skip so that the result
    # is in line with the 0.5-18 version of lavaan

    B0 <- B0.group[[1]]

    E <- B0

    eigvals <- eigen(E, symmetric = TRUE, only.values = TRUE)$values
    if (any(eigvals < -1 * .Machine$double.eps^(3/4))) {
            warning("lavaan WARNING: matrix based on first order outer product of the derivatives is not positive definite; the standard errors may not be thrustworthy")
    }
    NVarCov <- MASS::ginv(E)

    attr(NVarCov, "B0") <- B0
    attr(NVarCov, "B0.group") <- B0.group
    NVarCov
}

