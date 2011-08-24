# in time, we will add here various utilities to extract/compute
# information about the latent factors (in a CFA context)
#
# eg. (model-based) internal consistency reliability
#
# YR. 24 aug 2011

reliability <- function(object) {

    # only LISREL representation for now
    if(object@Model@representation != "LISREL") {
        stop("This function needs revision; it only works with the LISREL represenation!")
    }

    # check if we have latent variables
    nfac <- length(lavaanNames(object, type="lv"))
    if(nfac < 1) {
        stop("lavaanError: no latent factors in this model")
    }
    
    G <- object@Sample@ngroups
    theta.idx <- which( names(object@Model@GLIST) == "theta" )

    OUT <- vector("list", length=G)
    for(g in 1:G) {
               S <- object@Sample@cov[[g]]
            nvar <- object@Sample@nvar
        SigmaHat <- object@Fit@Sigma.hat[[g]]
           THETA <- object@Model@GLIST[[ theta.idx[g] ]]

        # 1. cronbach's alpha
        alpha <- nvar/(nvar-1) * (1.0 - sum(diag(S))/sum(S))

        # 2. model-based reliability (Bentler 2010)
        rho <- 1.0 - sum(THETA)/sum(SigmaHat)

        # 3. TODO... add more for the one-factor case!

        OUT[[g]] <- list(alpha = alpha,
                         rho   = rho)
    }

    if(G == 1) {
        OUT <- OUT[[1]]
    } else {
        names(OUT) <- object@Sample@group.label
    }

    OUT
}
