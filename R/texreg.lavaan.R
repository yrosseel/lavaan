# takes a fitted lavaan SEM  and (using texreg)
# returns latex syntax for a nicely-formatted regression 
# ``mod`` a fitted Lavaan SEM object with regressions
# ``ci`` Logical. should the confidence intervals from the SEM be reported instead of the standard errors?
# ``custom.model.names`` character vector. Should be the same length as the number of regressions in the SEM. 
#    The model names in the header of the regression table. Defaults to the names of the endogenous varibles. 
# ``...`` Other input to the texreg function


texreg.lavaan <- function(mod,ci=FALSE,custom.model.names,...){
    require(texreg)

    ctab <- parameterEstimates(mod)
    ctab <- droplevels(ctab[ctab$op=='~',])

    models <- split(ctab[,-c(1:2)],ctab$lhs)
    cmn <- if(missing(custom.model.names)) names(models) else custom.model.names

    tr <- lapply(models,function(m){
                     args <-  list(coef.names=m$rhs,
                                  coef=m$est,
                                  se=m$se,
                                   pvalues=m$pvalue)
                     if(ci){
                         args$ci.low=m$ci.lower
                         args$ci.up=m$ci.upper
                     }
                     do.call("createTexreg",args)})

    texreg(tr,custom.model.names=cmn,...)
}
