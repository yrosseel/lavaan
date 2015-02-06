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
  
  nsize <- fit_objH0@Data@nobs[[1]]   #[[1]] to be substituted by [[g]]?
  PLRT <- 2*nsize*(fit_objH0@Fit@fx - fit_objH1@Fit@fx)
  #!!!!!! we keep the multiplication by nsize because lavfit@fx gives the objective function divided by N

  
  # number of parameters under the alternative hypothesis H1
  Npar <- fit_objH1@Fit@npar  


   MY.m.el.idx2 <- fit_objH1@Model@m.free.idx
   # MY.m.el.idx2 gives the POSITION index of the free parameters within each 
   # parameter matrix under H1 model.
   # The index numbering restarts from 1 when we move to a new parameter matrix.
   # Within each matrix the index numbering "moves" columnwise.

   MY.x.el.idx2 <- fit_objH1@Model@x.free.idx
   # MY.x.el.idx2 ENUMERATES the free parameters within each parameter matrix.
   # The numbering continues as we move from one parameter matrix to the next one.

   # In the case of the symmetric matrices, Theta and Psi,in some functions below 
   # we need to give as input MY.m.el.idx2 and MY.x.el.idx2 after
   # we have eliminated the information about the redundant parameters
   # (those placed above the main diagonal).
   # That's why I do the following:

   MY.m.el.idx <- MY.m.el.idx2
   MY.x.el.idx <- MY.x.el.idx2
   # Psi, the variance - covariance matrix of factors
   if( length(MY.x.el.idx2[[3]])!=0 & any(table(MY.x.el.idx2[[3]])>1)) {
     nfac <- ncol(fit_objH1@Model@GLIST$lambda) #number of factors
     tmp  <- matrix(c(1:(nfac^2)), nrow= nfac, ncol= nfac )
     tmp_keep <- tmp[lower.tri(tmp, diag=TRUE)]
     MY.m.el.idx[[3]] <- MY.m.el.idx[[3]][MY.m.el.idx[[3]] %in% tmp_keep]
     MY.x.el.idx[[3]] <- unique( MY.x.el.idx2[[3]] )
   }

   #for Theta, the variance-covariance matrix of measurement errors
    if( length(MY.x.el.idx2[[2]])!=0 & any(table(MY.x.el.idx2[[2]])>1)) {
     nvar <- fit_objH1@Model@nvar #number of indicators
     tmp  <- matrix(c(1:(nvar^2)), nrow= nvar, ncol= nvar )
     tmp_keep <- tmp[lower.tri(tmp, diag=TRUE)]
     MY.m.el.idx[[2]] <- MY.m.el.idx[[2]][MY.m.el.idx[[2]] %in% tmp_keep]
     MY.x.el.idx[[2]] <- unique( MY.x.el.idx2[[2]] )
    }
 
   #below the commands to find the row-column indices of the Hessian that correspond to
   #the parameters to be tested equal to 0
   #tmp.ind contains these indices
   MY.m.el.idx2.H0 <- fit_objH0@Model@m.free.idx
   tmp.ind <- c()
   for(i in 1:6) {
     tmp.ind <- c(tmp.ind ,
                  MY.x.el.idx2[[i]] [!(MY.m.el.idx2[[i]]  %in%
                                       MY.m.el.idx2.H0[[i]] )  ]  )
   }
   # next line added by YR
   tmp.ind <- unique(tmp.ind)

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
 if(equalConstr==TRUE) {
     EqMat <- t(fit_objH0@Model@eq.constraints.K)
     #NoRowsEqMat <- nrow(EqMat)
     #NoColsEqMat <- ncol(EqMat)
     #for(i in 1:NoRowsEqMat){
     #  id_of_ones <- c(1:NoColsEqMat)[ EqMat[i,]==1] #vector with position index where 1's are
     #  tmp_values <- rep(c(1,-1), (length(id_of_ones) %/% 2)) #create pairs of the values 1, -1
     #  if((length(id_of_ones) %% 2)!=0) {  #if the number of parameters in row i is odd, 
     #      #then substitute the last -1 with -2 and add one more 1
     #      tmp_values <- tmp_values[-length(tmp_values)]
     #      tmp_values <- c(tmp_values, -2, 1) 
     #  }
     #  EqMat[i,id_of_ones] <- tmp_values 
     #}
 }


 # Compute the sum of the eigenvalues and the sum of the squared eigenvalues
 # so that the adjustment to PLRT can be applied.
 # Here a couple of functions (e.g. MYgetHessian) which are modifications of 
 # lavaan functions (e.g. getHessian) are needed. These are defined in the end of the file.  
 obj <- fit_objH0

 #the quantity below follows the same logic as getHessian of lavaan 0.5-18
 #and it actually gives N*Hessian. That's why the command following the command below. 
 NHes.theta0 <- MYgetHessian (object = obj@Model,
                            samplestats = obj@SampleStats ,
                            X = obj@Data@X ,
                            estimator = "PML",
                            lavcache = obj@Cache,
                            MY.m.el.idx = MY.m.el.idx,
                            MY.x.el.idx = MY.x.el.idx,
                            MY.m.el.idx2 = MY.m.el.idx2, # input for MYx2GLIST
                            MY.x.el.idx2 = MY.x.el.idx2, # input for MYx2GLIST
                            Npar = Npar,
                            equalConstr=equalConstr)
 Hes.theta0 <- NHes.theta0/ nsize #!!!!! to get the Hessian
 if(equalConstr==TRUE) { 
   Inv.Hes.theta0.final <- solve(Hes.theta0)
 } else {
   #Below manipulations so that the first rows-columns of the Hessian will correspond
   #to the parameters that are tested to be equal to 0.
   tmp.Hes.theta0.without <- Hes.theta0[-tmp.ind,-tmp.ind]
   #i.e. temporarily delete the rows and columns which correspond
   #to the parameters tested to be 0.

   tmp.row <- Hes.theta0[tmp.ind, -tmp.ind]
   # i.e. the corresponding row without the elements corresponding to the parameters
   
   Hes.theta0.perm <- rbind(tmp.row, tmp.Hes.theta0.without)
   if (length(tmp.ind)>1) { #if more than one parameters are tested
      tmp.col <- rbind(Hes.theta0[tmp.ind,tmp.ind], t(tmp.row))
   } else if (length(tmp.ind)==1) { #if one parameter is tested
      tmp.col <- c(Hes.theta0[tmp.ind,tmp.ind], tmp.row)
   }
   Hes.theta0.perm <- cbind(tmp.col , Hes.theta0.perm)

   #invert the Hessian and take only its first part that corresponds to
   # thetested parameters
   Inv.Hes.theta0.final <- solve( Hes.theta0.perm)
   H.to.psipsi.attheta0 <- Inv.Hes.theta0.final[1:length(tmp.ind), 1:length(tmp.ind)] #!!! correction here
 }

 #N times the estimated variability matrix is given :
 NJ.theta0 <- MYgetVariability(object = obj,
                             MY.m.el.idx = MY.m.el.idx,
                             MY.x.el.idx = MY.x.el.idx,
                             equalConstr=equalConstr)
 J.theta0.divbyN <-  NJ.theta0/(nsize^2) #!!!!!! added due to changes of lavaan functions
 if(equalConstr==TRUE) { 
   J.theta0.final <- J.theta0.divbyN
 } else {
   #similar manipulations as in the case of Hessian so that the rows and columns of J
   #referring to the parameters to be tested  appear first
   tmp.J.theta0.without <- J.theta0.divbyN[-tmp.ind,-tmp.ind]
   tmp.row <- J.theta0.divbyN[tmp.ind, -tmp.ind]
   J.theta0.perm <- rbind(tmp.row, tmp.J.theta0.without)
   if (length(tmp.ind)>1) {
      tmp.col <- rbind(J.theta0.divbyN[tmp.ind,tmp.ind], t(tmp.row))
   } else if (length(tmp.ind)==1) {
      tmp.col <- c(J.theta0.divbyN[tmp.ind,tmp.ind], tmp.row)
   }
   J.theta0.final <- cbind(tmp.col , J.theta0.perm)
 }
 
 
 # below the Inverse of the G matrix divided by N
 Inv.G <- Inv.Hes.theta0.final %*% J.theta0.final %*% Inv.Hes.theta0.final

 if(equalConstr==TRUE) { 
   MInvGtM <- EqMat %*% Inv.G %*% t(EqMat)
   MinvHtM <- EqMat %*% Inv.Hes.theta0.final %*% t(EqMat)  #!!!! corrected
   Amat <- MinvHtM #Alpha matrix in my notes
   InvAmat <- solve(Amat) 
   tmp.prod  <-  MInvGtM %*% InvAmat
   tmp.prod2 <- tmp.prod %*% tmp.prod
   sum.eig   <- sum( diag(tmp.prod) )
   sum.eigsq <- sum( diag(tmp.prod2) )
 } else {
   Inv.G.to.psipsi.attheta0 <- Inv.G[1:length(tmp.ind), 1:length(tmp.ind)]
   if (length(tmp.ind)>1) { #if more than one parameters are tested
     tmp.prod  <-  Inv.G.to.psipsi.attheta0 %*% solve( H.to.psipsi.attheta0)
     tmp.prod2 <- tmp.prod %*% tmp.prod
     sum.eig   <- sum( diag(tmp.prod) )
     sum.eigsq <- sum( diag(tmp.prod2) )
   } else if (length(tmp.ind)==1) { #if one parameter is tested
      #because here one element, the eigenvalue coincides with the following value
      sum.eig <- Inv.G.to.psipsi.attheta0 / H.to.psipsi.attheta0
      sum.eigsq <- sum.eig^2
   }
 }
 
 FSMA.PLRT <- (sum.eig    / sum.eigsq) * PLRT # adjusted PLRT
 adj.df    <- (sum.eig^2) / sum.eigsq # adjusted df
 pvalue <- 1-pchisq(FSMA.PLRT, df=adj.df )

 list(FSMA.PLRT = FSMA.PLRT, adj.df = adj.df, pvalue=pvalue)
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

