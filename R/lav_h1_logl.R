# compute logl for the unrestricted (h1) model -- per group
lav_h1_logl <- function(lavdata        = NULL,
                        lavsamplestats = NULL,
                        lavoptions     = NULL) {

    # number of groups
    ngroups <- lavdata@ngroups
    logl <- rep( as.numeric(NA, ngroups) )

    # only for ML (for now)
    if(lavoptions$estimator == "ML") {
        for(g in seq_len(ngroups) ) {
            if(lavdata@nlevels > 1L) {
                # not ready yet
                #logl[g] <- lav_mvnorm_cluster_loglik_samplestats_2l()
            } else if(lavsamplestats@missing.flag) {
                logl[g] <-
                    lav_mvnorm_missing_loglik_samplestats(
                        Yp    = lavsamplestats@missing[[g]],
                        Mu    = lavsamplestats@missing.h1[[g]]$mu,
                        Sigma = lavsamplestats@missing.h1[[g]]$sigma)
            } else { # single-level, complete data
                if(lavoptions$conditional.x) {
                    # TODO
                    #logl[g] <- lav_mvreg_loglik_samplestats(
                    #    sample.res.beta =
                    #    sample.res.cov=
                    #    sample.XX=
                    #    sample.nobs=
                    #    Beta=
                    #    Sigma=)
                } else {
                    logl[g] <-
                        lav_mvnorm_h1_loglik_samplestats(
                            sample.mean = lavsamplestats@mean[[g]],
                            sample.cov.inv = lavsamplestats@icov[[g]],
                            sample.nobs = lavsamplestats@nobs[[g]])
                }
            } # complete
        } # g
    } # ML

    logl
}
