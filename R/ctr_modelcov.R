# takes a model in lavaan syntax and the user's data and returns the covariance
# matrix of observed variables.  Useful so that the user can do things like
# diagnose errors in the cov matrix, use cov2cor to look at the correlation
# matrix, try and invert the sample covariance matrix, etc.

# update 5/27/2011 JEB
# changelog: using sem and inspect to get output.
#            This way, all arguments such as groups, etc, can be used
# update 3 june 2011 YR: removed se="none" (since now implied by do.fit=FALSE)
# update 13 dec 2011 YR: changed name (to avoid confusion with the
#                        model-implied cov)
inspectSampleCov <- function(model, data, ...) {
  fit <- sem(model, data = data, ..., do.fit = FALSE)
  inspect(fit, "sampstat")
}
