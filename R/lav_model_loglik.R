# compute the loglikelihood of the data, given the model
lav_model_loglik <- function(lavdata        = NULL,
                             lavsamplestats = NULL,
                             lavimplied        = NULL,
                             lavoptions     = NULL) {

    ngroups <- lavdata@ngroups

    logl.group <- rep(as.numeric(NA), ngroups)

    # should compute logl, or return NA?
    logl.ok <- FALSE
    if(lavoptions$estimator %in% c("ML", "MML")) {
        # check if everything is numeric, OR if we have exogenous
        # factor with 2 levels only
        if(all(lavdata@ov$type == "numeric")) {
            logl.ok <- TRUE
        } else {
            not.idx <- which(lavdata@ov$type != "numeric")
            for(i in not.idx) {
                if(lavdata@ov$type[i] == "factor" &&
                   lavdata@ov$exo[i] == 1L &&
                   lavdata@ov$nlev[i] == 2L) {
                    logl.ok <- TRUE
                } else {
                    logl.ok <- FALSE
                    break
                }
            }
        }
    }
  
    # lavsamplestats filled in? (not if no data...)
    if(length(lavsamplestats@ntotal) == 0L) {
        logl.ok <- FALSE
    }

    if(logl.ok) {
        for(g in seq_len(ngroups) ) {
            if(lavdata@nlevels > 1L) {
                # not ready yet
                #logl.group[g] <- lav_mvnorm_cluster_loglik_samplestats_2l()
            } else if(lavsamplestats@missing.flag) {
                logl.group[g] <-
                    lav_mvnorm_missing_loglik_samplestats(
                        Yp    = lavsamplestats@missing[[g]],
                        Mu    = lavimplied$mean[[g]],
                        Sigma = lavimplied$cov[[g]])
            } else { # single-level, complete data
                if(lavoptions$conditional.x) {
                    sample.beta <- t( cbind(lavsamplestats@res.int[[g]],
                                            lavsamplestats@res.slopes[[g]]) )
                    sample.XX   <- crossprod(cbind(1,lavdata@eXo[[g]]))
                    Beta        <- t( cbind(lavimplied$res.int[[g]],
                                            lavimplied$res.slopes[[g]]) )
                    logl.group[g] <- lav_mvreg_loglik_samplestats(
                        sample.res.beta = sample.beta,
                        sample.res.cov  = lavsamplestats@res.cov[[g]],
                        sample.XX       = sample.XX,
                        sample.nobs     = lavsamplestats@nobs[[g]],
                        Beta            = Beta,
                        Sigma           = lavimplied$res.cov[[g]],
                        Sinv.method     = "eigen",
                        Sigma.inv       = NULL)
                } else {
                    if(lavoptions$meanstructure) {
                        Mu <- lavimplied$mean[[g]]
                    } else {
                        Mu <- lavsamplestats@mean[[g]]
                    }
                    logl.group[g] <- lav_mvnorm_loglik_samplestats(
                        sample.mean = lavsamplestats@mean[[g]],
                        sample.cov  = lavsamplestats@cov[[g]],
                        sample.nobs = lavsamplestats@nobs[[g]],
                        Mu          = Mu,
                        Sigma       = lavimplied$cov[[g]],
                        Sinv.method = "eigen",
                        Sigma.inv   = NULL)
                }
            } # complete
        } # g
    } # logl.ok is TRUE



    out <- list(loglik        = sum(logl.group),
                loglik.group  = logl.group,
                npar          = 0L,
                AIC           = as.numeric(NA),
                BIC           = as.numeric(NA),
                estimator     = lavoptions$estimator,
                conditional.x = lavoptions$conditional.x,
                fixed.x       = lavoptions$fixed.x)
    out
}
