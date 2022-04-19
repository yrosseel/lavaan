# Browne's residual test statistic
# see Browne (1984) eq 2.20a

# T.B = (N-1) * t(RES) %*% Delta.c %*%
#                          solve(t(Delta.c) %*% Gamma %*% Delta.c) %*%
#                          t(Delta.c) %*% RES
#
#     = (N-1) * t(RES) %*% (Gamma.inv -
#                           Gamma.inv %*% Delta %*%
#                                  solve(t(Delta) %*% Gamma.inv %*% Delta) %*%
#                           t(Delta) %*% Gamma.inv) %*% RES


# TODo: - allow for 'structured' (model-based) version of Gamma
#       - allow for unbiased version of Gamma
#       - allow for non-linear equality constraints
#         (see Browne, 1982, eq 1.7.19)

lav_test_browne <- function(lavobject = NULL, n.minus.one = "default",
                            ADF = TRUE) {

    # slots
    lavdata <- lavobject@Data
    lavsamplestats <- lavobject@SampleStats
    lavmodel <- lavobject@Model

    # check input
    if(!inherits(lavobject, "lavaan")) {
        stop("lavaan ERROR: object is not a lavaan object.")
    }
    if(!ADF && lavmodel@categorical) {
        stop("lavaan ERROR: normal theory version not available in the categorical setting.")
    }
    if(lavdata@nlevels > 1L) {
        stop("lavaan ERROR: Browne's test is not available when data is multilevel.")
    }
    if(length(lavmodel@ceq.nonlinear.idx) > 0L) {
        stop("lavaan ERROR: Browne's test is not available (yet) when nonlinear equality constraints are involved.")
    }

    if(!is.logical(n.minus.one)) {
        if(lavobject@Options$estimator == "ML") {
            n.minus.one <- FALSE
        } else {
            n.minus.one <- TRUE
        }
    }

    # linear equality constraints?
    lineq.flag <- FALSE
    if(lavmodel@eq.constraints) {
        lineq.flag <- TRUE
    } else if(.hasSlot(lavmodel, "ceq.simple.only") &&
              lavmodel@ceq.simple.only) {
        lineq.flag <- TRUE
    }

    # ingredients
    Delta <- computeDelta(lavmodel)
    if(ADF) {
        # ADF version
        if(!is.null(lavsamplestats@NACOV[[1]])) {
            Gamma <- lavsamplestats@NACOV
        } else {
            Gamma <- lavGamma(lavobject)
        }
    } else {
        # NT version
        Gamma <- lavGamma(lavobject, ADF = FALSE,
                          NT.rescale = lavobject@Options$estimator == "ML")
    }
    WLS.obs <- lavsamplestats@WLS.obs
    WLS.est <- lav_model_wls_est(lavmodel)
    nobs <-  lavsamplestats@nobs
    ntotal <- lavsamplestats@ntotal

    # compute T.B per group
    ngroups <- lavdata@ngroups
    stat.group <- numeric(ngroups)


    # 1. standard setting: no equality constraints
    if(!lineq.flag) {
        for(g in seq_len(ngroups)) {
            RES <- WLS.obs[[g]] - WLS.est[[g]]
            Delta.g <- Delta[[g]]
            Delta.c <- lav_matrix_orthogonal_complement(Delta.g)
            tDGD <- solve(crossprod(Delta.c, Gamma[[g]]) %*% Delta.c)
            if(n.minus.one) {
                Ng <- nobs[[g]] - 1L
            } else {
                Ng <- nobs[[g]]
            }
            tResDelta.c <- crossprod(RES, Delta.c)
            stat.group[g] <- Ng * drop(tResDelta.c %*% tDGD %*% t(tResDelta.c))
        }
        STAT <- sum(stat.group)

    # 2. linear equality constraint
    } else if(lineq.flag) {
        RES.all <- do.call("c", WLS.obs) - do.call("c", WLS.est)
        Delta.all <- do.call("rbind", Delta)
        if(lavmodel@eq.constraints) {
            Delta.g <- Delta.all %*% lavmodel@eq.constraints.K
        } else if(.hasSlot(lavmodel, "ceq.simple.only") &&
                  lavmodel@ceq.simple.only) {
            Delta.g <- Delta.all %*% lavmodel@ceq.simple.K
        }
        Gamma.inv.weighted <- vector("list", ngroups)
        for(g in seq_len(ngroups)) {
            if(n.minus.one) {
                Ng <- nobs[[g]] - 1L
            } else {
                Ng <- nobs[[g]]
            }
            Gamma.inv.weighted[[g]] <- solve(Gamma[[g]]) * Ng/ntotal
        }
        GI <- lav_matrix_bdiag(Gamma.inv.weighted)
        q1 <- drop( t(RES.all) %*% GI %*% RES.all)
        q2 <- drop( t(RES.all) %*%
                    GI %*% Delta.g %*%
                    solve( t(Delta.g) %*% GI %*% Delta.g) %*%
                    t(Delta.g) %*% GI %*%
                    RES.all )
        STAT <- ntotal * (q1 - q2)
        stat.group <- STAT * unlist(nobs) / ntotal # proxy only

    # 3. nonlinear equality constraints
    } else {
        # TODO
    }


    DF <- lavobject@test[[1]]$df
    if(ADF) {
        NAME <- "browne.residual.adf"
    } else {
        NAME <- "browne.residual.nt"
    }
    out <- list(test       = NAME,
                stat       = STAT,
                stat.group = stat.group,
                df         = DF,
                pvalue     = 1 - pchisq(STAT, DF))
    out
}
