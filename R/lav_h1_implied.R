# compute sample statistics for the unrestricted (h1) model
lav_h1_implied <- function(lavdata        = NULL,
                           lavsamplestats = NULL,
                           lavoptions     = NULL) {

    if(lavdata@nlevels == 1L) {
        if(lavsamplestats@missing.flag) {
            if(lavoptions$conditional.x) {
                implied <- list() # not available yet
            } else {
                implied <- list(cov  = lapply(lavsamplestats@missing.h1,
                                              "[[", "sigma"),
                                mean = lapply(lavsamplestats@missing.h1,
                                              "[[", "mu"))
            }
        } else {
            if(lavoptions$conditional.x) {
                implied <- list(res.cov    = lavsamplestats@res.cov,
                                res.int    = lavsamplestats@res.int,
                                res.slopes = lavsamplestats@res.slopes,
                                res.th     = lavsamplestats@res.th)
            } else {
                implied <- list(cov    = lavsamplestats@cov,
                                mean   = lavsamplestats@mean,
                                th     = lavsamplestats@th)
            }
        } # complete data
    } else {
        # estimate Mu.B, Mu.W, Sigma.B and Sigma.W for unrestricted model
        # TODO
        implied <- list()
    }

    implied
}

lav_h1_implied_cluster <- function(lavdata        = NULL,
                                   lavsamplestats = NULL,
                                   lavoptions     = NULL) {

}
