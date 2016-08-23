## Author: Carl F. Falk
## Version date for first lavaan commit: 2016-08-22
## First version assumed no access to internal lavaan functions.
## A re-write could make this more efficient and flexible or for
## integration with parTable
##
## General purpose likelihood-based confidence interval function using lavaan
## fitmodel - a fitted lavaan model
## pindx - index numbers or characters corresponding to parameter estimates
##  Should match order or parameters in ParTable or label in ParTable.
## alpha - confidence level (e.g., .05 -> .95 confidence interval)
## bound - "lower" for lower bound, "upper" for upper bound
## func - function to compute on model parameters
## optimizer - choice of optimizer; Rsolnp and optim currently implemented
## p.start - These values are used as starting values in optimization
##  and could be different than those in the initially fitted model
## diffmethod - method to use for difference test in lavTestLRT
lci<-function(fitmodel, pindx, alpha=.05, bound=c("lower","upper"),
              func=function(x){x}, optimizer=c("Rsolnp","optim"),
              p.start=NULL,diffmethod="default",...){
  if(optimizer[1]=="Rsolnp"){
    require(Rsolnp)
  }
  
  ## input checking
  if(!is.null(p.start)){
    if(length(pindx)!=length(p.start)){
      stop("pindx and p.start should be the same length")
    }
  }
  if(class(fitmodel)!="lavaan"){
    stop("fitmodel must be a fitted lavaan model")
  }
  if(alpha <=0 | alpha >=1){
    stop("alpha must be between 0 and 1")
  }
  
  ## extract parameter table for modification  
  ParTable<-parTable(fitmodel) #fitmodel@ParTable 

  ## Get parameter estimates and
  ## get parameter index if named by a label;
  ## Otherwise use combination of lhs, op, rhs
  p.est<-vector('numeric')
  tmp.indx<-vector("numeric")
  lhsoprhs<-paste(ParTable$lhs,ParTable$op,ParTable$rhs, sep="")
  for(i in 1:length(pindx)){
    if(is.character(pindx[i])){
      indx<-which(ParTable$label==pindx[i]|lhsoprhs==pindx[i])
      if(length(indx)>1){
        warning(paste("pindx ", pindx[i],
                      " matched more than one parameter; using the first match"))
      } else if (length(indx)==0){
        warning(paste("pindx ", pindx[i],
                      " did not match any known parameters"))
      }
      tmp.indx[i]<-indx[1]
    } else{
      tmp.indx[i]<-pindx[i]
    }
    p.est<-c(p.est, ParTable$est[tmp.indx[i]])
  }
  if(is.character(pindx)){
    pindx<-tmp.indx
  } 
  
  ## point estimate based on fitted model
  est<-do.call(func,list(p.est))
  result<-list()
  result$est<-est
    
  if(is.null(p.start)){
    p.start<-p.est
  }
  
  for(b in bound){
    if(optimizer[1]=="Rsolnp"){
      capture.output(LCI<-try(solnp(p.start,fit.diff.test.pfunc,fitmodel=fitmodel,
                     pindx=pindx, func=func, alpha=alpha, bound=b,
                     diffmethod=diffmethod,...)))
      if(class(LCI)!="try-error"){
        p.est.bound<-LCI$pars
        p.est.bound<-do.call(func,list(p.est.bound))
        conv<-LCI$convergence
      } else{
        p.est.bound<-NA
      }
    } else if (optimizer[1]=="optim"){
      LCI<-try(optim(p.start, fit.diff.test.pfunc,fitmodel=fitmodel,
                     pindx=pindx, func=func, alpha=alpha, bound=b,
                     diffmethod=diffmethod,...))
      if(class(LCI)!="try-error"){
        p.est.bound<-LCI$par
        p.est.bound<-do.call(func,list(p.est.bound))
        conv<-LCI$convergence
      } else{
        p.est.bound<-NA
      }
    }
    
    if(b=="upper"){
      result$upper<-p.est.bound
      result$convupper<-conv
    }
    if(b=="lower"){
      result$lower<-p.est.bound
      result$convlower<-conv
    }
  }
  
  if(!is.null(result$upper)&!is.null(result$lower)){
    if(result$upper < result$lower |
       result$upper < result$est | result$lower > result$est){
      warning("Optimization may have failed.
              Check boundaries vs. each other and point estimate")
      result$boundcode<-1
    } else {
      result$boundcode<-0
    }    
  } else {
    warning("One or more boundaries are null.
            Either only one was requested or optimization may have failed.")
    result$boundcode<-2
  }
  
  return(result)
}

## Optimization for function of model parameters
fit.diff.test.pfunc<-function(p,fitmodel,pindx,func=function(x){x},
                              alpha=.05,bound=c("lower","upper"),
                              diffmethod="default",...){
  
  ParTable<-parTable(fitmodel) ## extract parameter table for modification
  
  ParTable$ustart<-ParTable$est ## Starting values at MLE
  ParTable$start<-NULL
  ParTable$est<-NULL
  ParTable$se<-NULL
  #myModel1<-ParTable
  
  ## In case of multiple parameters in p and pindx
  for(i in 1:length(pindx)){
    ParTable$free[pindx[i]]<-0
    ParTable$ustart[pindx[i]]<-p[i]
  }
  
  ParTable$id<-1:length(ParTable$id)

  ## Attempt to fit model
  prevmodel<-as.list(fitmodel@call)
  prevmodel$model<-ParTable
  M1<-try(do.call(as.character(prevmodel[1]),prevmodel[-1]),silent=TRUE)

  if(class(M1)=="try-error"){
    return(NA)
  } else {
    if(!M1@Fit@converged){
      return(NA)
    } else {
      D1<-lavTestLRT(fitmodel,M1,method=diffmethod)$`Chisq diff`[2]
    }
  }
  
  ## Compute function of model parameters
  g<-do.call(func,list(p))
  
  ## Check against critical value of chi-square
  crit<-qchisq(1-alpha,1) ## hardcoded 1 df for now
  fit<-(D1-crit)^2
  
  if(bound=="lower"){
    fit<-fit+g
  } else if (bound=="upper"){
    fit<-fit-g
  } else {
    stop("No valid bound specified")
  }
  return(fit)
}






