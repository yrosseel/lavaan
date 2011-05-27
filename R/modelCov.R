#takes a model in lavaan syntax and the user's data and returns the covariance matrix
#of observed variables.  Useful so that the user can do things like diagnose errors
#in the cov matrix, use cov2cor to look at the correlation matrix, try and invert the
#sample covariance matrix, etc.

#last updated 5/27/2011 JEB
# changelog: using sem and inspect to get output.  This way, all arguments such as groups, etc, can be used

#modelCovOld<-function(model, data, na.rm=T, ...){
#  model<-lavaanify(model)
#  vars<-unique(c(model$lhs, model$rhs))
#  data<-data[,which(names(data) %in% vars)]
#  cov(data, na.rm=na.rm)
#}


#modelCovOld<-function(model, data, ...){
#  ov<-vnames(flatten.model.syntax(model), "ov")
#  samp<-Sample(data=data, ov.names=ov, data.type="full", ...)
#  cmat<-samp@cov[[1]]
#  rownames(cmat)<-colnames(cmat)<-samp@ov.names
#  cmat
#}

modelCov<-function(model, data, ...) {
    fit <- sem(model, data=data, ..., se="none", do.fit=FALSE)
    inspect(fit, "sampstat")
} 