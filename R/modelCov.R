#takes a model in lavaan syntax and the user's data and returns the covariance matrix
#of observed variables.  Useful so that the user can do things like diagnose errors
#in the cov matrix, use cov2cor to look at the correlation matrix, try and invert the
#sample covariance matrix, etc.

#last updated 5/25/2011 JEB

modelCovOld<-function(model, data, na.rm=T, ...){
  model<-lavaanify(model)
  vars<-unique(c(model$lhs, model$rhs))
  data<-data[,which(names(data) %in% vars)]
  cov(data, na.rm=na.rm)
}


modelCov<-function(model, data, ...){
  ov<-vnames(flatten.model.syntax(model), "ov")
  samp<-Sample(data=data, ov.names=ov, data.type="full", ...)
  cmat<-samp@cov[[1]]
  rownames(cmat)<-colnames(cmat)<-samp@ov.names
  cmat
}