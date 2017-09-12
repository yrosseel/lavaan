# compute the loglikelihood of the data, given the model
lav_model_loglik <- function(lavdata        = NULL,
                             lavsamplestats = NULL,
                             lavimplied     = NULL,
                             lavmodel       = NULL,
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
        # nlevels + fiml
        if(lavdata@nlevels > 1L && lavsamplestats@missing.flag) {
            logl.ok <- FALSE
        }
    }
  
    # lavsamplestats filled in? (not if no data...)
    if(length(lavsamplestats@ntotal) == 0L) {
        logl.ok <- FALSE
    }

    if(logl.ok) {
        for(g in seq_len(ngroups) ) {
            if(lavdata@nlevels > 1L) {
                # here, we assume only 2 levels, at [[1]] and [[2]]
                Sigma.W <- lavimplied$cov[[  (g-1)*2 + 1]]
                Mu.W    <- lavimplied$mean[[ (g-1)*2 + 1]]
                Sigma.B <- lavimplied$cov[[  (g-1)*2 + 2]]
                Mu.B    <- lavimplied$mean[[ (g-1)*2 + 2]]

                logl.group[g] <- lav_mvnorm_cluster_loglik_samplestats_2l(
                    YLp          = lavsamplestats@YLp[[g]],
                    Lp           = lavdata@Lp[[g]],
                    Mu.W         = Mu.W,
                    Sigma.W      = Sigma.W,
                    Mu.B         = Mu.B,
                    Sigma.B      = Sigma.B,
                    Sinv.method  = "eigen",
                    log2pi       = TRUE,
                    minus.two    = FALSE)

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

    # logl
    logl <- sum(logl.group)

    # number of parameters, taking into account any equality constraints
    npar <- lavmodel@nx.free
    if(nrow(lavmodel@con.jac) > 0L) {
        ceq.idx <- attr(lavmodel@con.jac, "ceq.idx")
        if(length(ceq.idx) > 0L) {
            neq <- qr(lavmodel@con.jac[ceq.idx,,drop=FALSE])$rank
            npar <- npar - neq
        }
    }

    # logl
    logl <- sum(logl.group)

    if(logl.ok) {
        # AIC
        AIC <- (-2 * logl) + (2 * npar)

        # BIC
        BIC <- (-2 * logl) + (npar * log(lavsamplestats@ntotal))

        # BIC2
        N.star <- (lavsamplestats@ntotal + 2) / 24
        BIC2 <- (-2 * logl) + (npar * log(N.star))   
    } else {
        AIC <- BIC <- BIC2 <- as.numeric(NA)
    }

    out <- list(loglik        = logl,
                loglik.group  = logl.group,
                npar          = npar,
                ntotal        = lavsamplestats@ntotal,
                AIC           = AIC,
                BIC           = BIC,
                BIC2          = BIC2,
                estimator     = lavoptions$estimator,
                conditional.x = lavoptions$conditional.x,
                fixed.x       = lavoptions$fixed.x)
    out
}
