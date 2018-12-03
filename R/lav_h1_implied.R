# compute sample statistics for the unrestricted (h1) model
# and also the logl (if available)
lav_h1_implied_logl <- function(lavdata        = NULL,
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

        logl <- lav_h1_logl(lavdata        = lavdata,
                            lavsamplestats = lavsamplestats,
                            lavoptions     = lavoptions)
    } else {
        # estimate Mu.B, Mu.W, Sigma.B and Sigma.W for unrestricted model
        ngroups <- lavdata@ngroups
        nlevels <- lavdata@nlevels
        implied <- list(cov  = vector("list", length = ngroups * nlevels),
                        mean = vector("list", length = ngroups * nlevels))
        loglik.group <- numeric(lavdata@ngroups)

        for(g in 1:lavdata@ngroups) {
            if(lavoptions$verbose) {
                cat("\nFitting unrestricted (H1) model in group ", g, "\n")
            }
            OUT <- lav_mvnorm_cluster_em_sat(YLp      = lavsamplestats@YLp[[g]],
                                             Lp       = lavdata@Lp[[g]],
                                             verbose  = lavoptions$verbose,
                                             tol      = 1e-04,     # option?
                                             min.variance = 1e-05, # option?
                                             max.iter = 5000L)     # option?
            if(lavoptions$verbose) {
                cat("\n")
            }

            # if any near-zero within variance(s), produce warning here
            zero.var <- which(diag(OUT$Sigma.W) <= 1e-05)
            if(length(zero.var)) {
                gtxt <- if(ngroups > 1L) {
                                paste(" in group ", g, ".", sep = "")
                            } else { " " }
                txt <- c("H1 estimation resulted in a within covariance matrix",
                          gtxt, "with (near) zero variances for some of the
                          level-1 variables: ",
                         lavdata@ov.names.l[[g]][[1]][zero.var])
                warning(lav_txt2message(txt))
            }

            # store in implied
            implied$cov[[(g-1)*nlevels + 1L]] <- OUT$Sigma.W
            implied$cov[[(g-1)*nlevels + 2L]] <- OUT$Sigma.B
            implied$mean[[(g-1)*nlevels + 1L]] <- OUT$Mu.W
            implied$mean[[(g-1)*nlevels + 2L]] <- OUT$Mu.B

            # store logl per group
            loglik.group[g] <- OUT$logl
        }

        logl <- list(loglik = sum(loglik.group), loglik.group = loglik.group)
    }

    list(implied = implied, logl = logl)
}

