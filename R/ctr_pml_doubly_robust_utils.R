# This code was contributed by Myrsini Katsikatsou (LSE) -- September 2016

# The input of the function is a lavobject, which, in turn, is the output of the
# sem function having specified estimator="PML", missing="available.cases" 

#The output of the function is a list of two lists: the pairwiseProbGivObs list and
# the univariateProbGivObs list. Each of the two lists consists of G matrices where G
# is the number of groups in a multigroup analysis. If G=1 each of the lists 
# contains only one matrix that can be called as pairwiseProbGivObs[[1]], and
# univariateProbGivObs[[1]].

# Each of the matrices in the pairwiseProbGivObs list is of dimension: nrow=sample size, 
#ncol=sum of the number of response categories for all pairs of variables 
#(i.e. the length of the vector pxixj.ab where i<j=1,...,p, a=1,...,Ci, b=1,...,Cj; 
#a which is the index for the response category for yi variable runs the fastest, 
#then b which is the index for the response category for yj variable, 
#then j, and last i.) 
#The cells in a matrix of the pairwiseProbGivObs list have the value 0 except for
#those cells that correspond to the pairs of variables where both variables 
#are missing. Those cells have the value of the bivariate conditional probability
#for the given pair for all their response categories. The bivariate 
#probabilities are computed as follows: 
#the information in the observed variables is summarised in a factor score 
#for each individual and the bivariate probability given the estimated factor 
#scores is computed.


# Each of the matrices in the univariateProbGivObs list is of dimension: 
# nrow=sample size, ncol=sum of the number of response categories for all 
# variables. The columns are indexed with i and a, where i=1,...,p, and
# a=1,...,Ci, the response categories for yi variable; a runs faster than i. 
#The cells in a matrix of the univariateProbGivObs list have the value 0 except for
#those cells that correspond to variables with missing values. 
#Those cells have the value of the univariate conditional probability for the
#given variable for all its response categories. The univariate conditional 
#probabilities are computed as follows: 
#given that the bivariate conditional probabilities have been computed we sum over
#the response categories of each variable at a time (i.e. we compute the marginals).

#Version 3 - first compute univariate and then bivariate probabilities

pairwiseExpProbVec_GivenObs <- function(lavobject) {

 #compute yhat where yaht=nu + Lamda*eta + K*x where the parameter estimates are 
 #used and the factor scores for eta 
 #Below yhat is a list if lavobject@Data@ngroups >1, it is a list of G matrices
 #where G the number of groups and the matrices are fo dimension 
 #nrow=sample size and ncol=number of items.
 #If lavobject@Data@ngroups=1 then yhat is a matrix.
 yhat <- lavPredict(object=lavobject, type = "yhat" )

 #compute bivariate probabilities
 ngroups <- lavobject@Data@ngroups
 univariateProb <- vector("list", length=ngroups)
 pairwiseProb <- vector("list", length=ngroups)
 #save the indices of the Theta matrices for the groups stored in GLIST
 idx.ThetaMat <- which(names(lavobject@Model@GLIST)=="theta") 

 for(g in seq_len(ngroups)) { # g<-1
 
  if(ngroups>1L){
    yhat_group <- yhat[[g]]
  } else {
    yhat_group <- yhat
  }

  nsize  <- lavobject@Data@nobs[[g]]
  nvar   <- lavobject@Model@nvar[[g]]
  Data   <- lavobject@Data@X[[g]]
  TH     <- lavobject@Fit@TH[[g]]
  th.idx <- lavobject@Model@th.idx[[g]]
  Theta  <- lavobject@Model@GLIST[ idx.ThetaMat[g] ]$theta 
  error.stddev <- diag(Theta)^0.5 

  #for the computation of the univariete probabilities
  nlev <- lavobject@Data@ov$nlev
  idx.uniy <- rep(1:nvar, times=nlev)

  univariateProb[[g]] <- matrix(0, nrow = nsize, ncol = sum(nlev) ) 
  pairwiseProb[[g]]   <- matrix(0, nrow = nsize, 
                        ncol = length(lavobject@Cache[[g]]$bifreq))
  
  idx.MissVar.casewise <- apply(Data, 1, function(x) { which(is.na(x)) } )
  
  for(i in 1:nsize){
    idx.MissVar <- idx.MissVar.casewise[[i]]
    noMissVar   <- length(idx.MissVar)
    
    if( noMissVar>0L ) {
    #compute the univariate probabilities
    TH.list <- split(TH,th.idx)
    tmp.TH  <- TH.list[idx.MissVar]
    tmp.lowerTH <- unlist(lapply(tmp.TH, function(x){c(-100,x)}))
    tmp.upperTH <- unlist(lapply(tmp.TH, function(x){c(x,100) }))
    
    idx.items <- rep(c(1:noMissVar), times=nlev[idx.MissVar])
    tmp.mean    <- yhat_group[i,idx.MissVar]
    tmp.mean.extended <- tmp.mean[idx.items]
    tmp.stddev  <- error.stddev[idx.MissVar]
    tmp.stddev.extended <- tmp.stddev[idx.items]
    
    tmp.uniProb <- pnorm( (tmp.upperTH - tmp.mean.extended )/
                           tmp.stddev.extended ) -
                   pnorm( (tmp.lowerTH - tmp.mean.extended )/
                           tmp.stddev.extended )
    idx.columnsUni <- which(idx.uniy %in% idx.MissVar)
    univariateProb[[g]][i, idx.columnsUni] <- tmp.uniProb
   
  #compute the bivariate probabilities 
  if( noMissVar>1L ) {
   idx.pairsMiss <- combn(idx.MissVar ,2)
   no.pairs <- ncol(idx.pairsMiss) 
   idx.pairsV2 <- combn(noMissVar, 2)
   
   idx.columns <- c() 
   for(j in 1:no.pairs ) {
      tmp.idx <- which( (idx.y1 == idx.pairsMiss[1,j]) & 
                        (idx.y2 == idx.pairsMiss[2,j]) )
      idx.columns <- c(idx.columns, tmp.idx) 
   }
   
   if( all( Theta[t(idx.pairsMiss)]==0 ) ){ #items independence given eta 
     tmp.uniProb.list <- split(tmp.uniProb, idx.items)
     tmp.bivProb <- c()
     for(v1 in 1:(noMissVar-1) ) {
       for(v2 in (v1+1):noMissVar){
         tmp.bivProb <- c(tmp.bivProb, 
                            c( outer(tmp.uniProb.list[[v1]] ,
                                     tmp.uniProb.list[[v2]]  ) ) )
       }
     }
     pairwiseProb[[g]][i, idx.columns] <- tmp.bivProb 
   }  else { #when correlation between measurement errors
     
     tmp.th.idx <- th.idx[th.idx %in% idx.MissVar]
     #recode so that it is always 1,1,..,1, 2,...,2, etc.
     tmp.th.idx.recoded <- rep(c(1:noMissVar), times=table(tmp.th.idx))
     tmp.TH <- TH[th.idx %in% idx.MissVar]
    
     tmp.ind.vec <- LongVecInd(no.x = noMissVar, 
                                   all.thres = tmp.TH, 
                          index.var.of.thres = tmp.th.idx.recoded)
                 
     tmp.th.rho.vec <-  LongVecTH.Rho.Generalised(
                           no.x = noMissVar, 
                             TH = tmp.TH, 
                         th.idx = tmp.th.idx.recoded,
                       cov.xixj = Theta[t(pairsMiss)] ,
                         mean.x = yhat_group[i,idx.MissVar],
                       stddev.x = error.stddev[idx.MissVar] ) 
 
     tmp.bivProb <- pairwiseExpProbVec(ind.vec = tmp.ind.vec ,
                                             th.rho.vec = tmp.th.rho.vec)
    
     pairwiseProb[[g]][i, idx.columns] <- tmp.bivProb 
   } #end of else of if( all( Theta[t(idx.pairsMiss)]==0 ) ) 
     # which checks item local independence 
  } #end of if( noMissVar>1L )
   
     #cat(i,  "\n")
  } #end of if(noMissVar>0L)   
  } #end of for(i in 1:nsize)
  
 } #end of for(g in seq_len(lavobject@Data@ngroups))                           
 
 list(univariateProbGivObs = univariateProb,
        pairwiseProbGivObs = pairwiseProb)
} # end of the function pairwiseExpProbVec_GivenObs

##################################################################



# LongVecTH.Rho.Generalised function is defined as follows
LongVecTH.Rho.Generalised <- function(no.x, TH, th.idx,
                                      cov.xixj, mean.x, stddev.x ) {
 all.std.thres <- (TH - mean.x[th.idx]) / stddev.x[th.idx]
 id.pairs <- utils::combn(no.x,2)
 cor.xixj <- cov.xixj /( stddev.x[id.pairs[1,]] * stddev.x[id.pairs[2,]])
 
 LongVecTH.Rho(no.x = no.x, 
                   all.thres = all.std.thres, 
          index.var.of.thres = th.idx,
                    rho.xixj = cor.xixj)
}

# LongVecTH.Rho.Generalised is a generalisation  of the function
#  lavaan:::LongVecTH.Rho . The latter assumes that all y* follow standard
#  normal so the thresholds are automatically the standardised ones.
# LongVecTH.Rho.Generalised does not assume that, each of y*'s can follow
# a normal distribution with mean mu and standard deviation sigma.
# LongVecTH.Rho.Generalised has the following input arguments:
# no.x (same as in lavaan:::LongVecTH.Rho), 
# TH (similar to the TH in lavaan:::LongVecTH.Rho but here they are the unstandardised thresholds, i.e. of the normal distribution with mean mu and standard deviation sigma)
# th.idx (same as index.var.of.thres in lavaan:::LongVecTH.Rho)
# cov.xixj which are the polychoric covariances of the pairs of underlying variables provided in a similar fashion as rho.xixj in lavaan:::LongVecTH.Rho)
# mean.x  is a vector including the means of y*'s provided in the order mean.x1, mean.x2, ...., mean.xp
# stddev.x  is a vector including the standard deviations of y*'s provided in the order stddev.x1, stddev.x2, ...., stddev.xp

# The output of the new function is similar to that of lavaan:::LongVecTH.Rho#############################################


compute_uniCondProb_based_on_bivProb <- function(bivProb, nvar,
                                                 idx.pairs,
                                                 idx.Y1,
                                                 idx.Gy2,
                                                 idx.cat.y1.split,
                                                 idx.cat.y2.split) {
  bivProb.split <- split(bivProb, idx.pairs)
  lngth <- 2*length(bivProb)
  idx.vec.el <- 1:lngth
  ProbY1Gy2 <- rep(NA, lngth)
  no.pairs <- nvar*(nvar-1)/2
  idx2.pairs <- combn(nvar,2)

  for(k in 1:no.pairs){
      y2Sums <- tapply(bivProb.split[[k]], idx.cat.y2.split[[k]], sum)
      y2Sums.mult <- y2Sums[idx.cat.y2.split[[k]] ]
      Y1Gy2 <- bivProb.split[[k]]/ y2Sums.mult
      tmp.idx.vec.el <- idx.vec.el[(idx.Y1  == idx2.pairs[1,k]) &
                                   (idx.Gy2 == idx2.pairs[2,k])]
      ProbY1Gy2[tmp.idx.vec.el] <- Y1Gy2
  }

  for(k in 1:no.pairs){
      y1Sums <- tapply(bivProb.split[[k]], idx.cat.y1.split[[k]], sum)
      y1Sums.mult <- y1Sums[idx.cat.y1.split[[k]] ]
      Y2Gy1 <- bivProb.split[[k]]/ y1Sums.mult
      tmp.idx.vec.el <- idx.vec.el[(idx.Y1  == idx2.pairs[2,k]) &
                                   (idx.Gy2 == idx2.pairs[1,k])]
      ProbY1Gy2[tmp.idx.vec.el] <- Y2Gy1
  }

  ProbY1Gy2
}







