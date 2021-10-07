# collect information about the model that we can use
# (eg. is theta diagonal or not, is the structurual model recursive or not,
#  is the model just a regression model, etc)
#
# initial version: YR 15 March 2021
# - YR 05 Oct 2021: use det(I - B) to check if B is acyclic
 

lav_model_properties <- function(GLIST, lavpartable = NULL, lavpta = NULL,
                                 m.free.idx = NULL) {

    if(is.null(lavpta)) {
        lavpta <- lav_partable_attributes(lavpartable)
    }

    ngroups <- lavpta$ngroups


    # is the model a univariate/multivariate linear multiple regression
    # model (per group)?
    uvreg <- logical(ngroups)
    uvord <- logical(ngroups)
    mvreg <- logical(ngroups)
    acyclic <- logical(ngroups)
    nexo <- integer(ngroups)

    for(g in seq_len(ngroups)) {
       
       # at least 1 regression
       if(length(lavpta$vnames$eqs.y[[g]]) == 0L) {
           next
       }

       # acyclic?
       B <- GLIST$beta
       beta.idx <- which(names(GLIST) == "beta")
       # keep fixed values (if any); fill in 1 in all 'free' positions
       B[ m.free.idx[[beta.idx]] ] <- 1
       IB <- diag(nrow(B)) - B
       # if B is acyclic, we should be able to permute the rows/cols of B
       # so that B is upper/lower triangular, and so det(IB) = 1
       if(det(IB) == 1) {
           acyclic[g] <- TRUE
       }  


       # no latent variables, at least 1 dependent variable
       if(lavpta$nfac[[g]] > 0L) {
           next
       }

       # no mediators
       if(length(lavpta$vnames$eqs.y[[g]]) != length(lavpta$vnames$ov.y[[g]])) {
           next
       }

       # categorical y?
       if(length(lavpta$vnames$ov.ord[[g]]) > 0L) {
           # we only flag the univariate version
           if(length(lavpta$vnames$ov.ord[[g]]) == 1L &&
              length(lavpta$vnames$ov.y[[g]])   == 1L &&
              lavpta$vnames$ov.ord[[g]][1] == lavpta$vnames$ov.y[[g]][1]) {
               uvord[g] <- TRUE
           }
       } else {
           if(length(lavpta$vnames$ov.y[[g]]) > 1L) {
               mvreg[g] <- TRUE
           } else {
               uvreg[g] <- TRUE
           }
       }

       nexo[g] <- length(lavpta$vnames$eqs.x[[g]])

    } # g

    modprop <- list( uvreg = uvreg, uvord = uvord, mvreg = mvreg,
                     nexo = nexo, acyclic = acyclic )

    modprop
}
