# utility functions needed to compute various (robust) fit measures:
#
# - lav_fit_catml_dwls (for 'robust' RMSEA/CFI if data is cateogrical)
# - lav_fit_fiml_corrected (correct RMSEA/CFI if data is incomplete)

# compute scaling-factor (c.hat3) for fit.dwls, using fit.catml ingredients
# see:
#     Savalei, V. (2021) Improving Fit Indices In SEM with categorical data.
#     Multivariate Behavioral Research, 56(3), 390-407.
#
lav_fit_catml_dwls <- function(lavobject) {

    # empty list
    empty.list <- list(XX3 = as.numeric(NA), df3 = as.numeric(NA),
                       c.hat3 = as.numeric(NA), XX3.scaled = as.numeric(NA),
                       XX3.null = as.numeric(NA), df3.null = as.numeric(NA),
                       c.hat3.null = as.numeric(NA))

    # limitations
    if( !lavobject@Model@categorical ||
         lavobject@Options$conditional.x ||
         length(unlist(lavobject@pta$vnames$ov.num)) > 0L ) {
        return(empty.list)
    } else {
        lavdata        <- lavobject@Data
        lavsamplestats <- lavobject@SampleStats
    }

    # 'refit' using estimator = "catML"
    fit.catml <- try(lav_object_catml(lavobject), silent = TRUE)
    if(inherits(fit.catml, "try-error")) {
        return(empty.list)
    }

    XX3 <- fit.catml@test[[1]]$stat
    df3 <- fit.catml@test[[1]]$df


    # compute 'k'
    V      <- lavTech(fit.catml, "wls.v") # NT-ML weight matrix

    W.dwls <- lavTech(lavobject, "wls.v") # DWLS weight matrix
    Gamma  <- lavTech(lavobject, "gamma") # acov of polychorics
    Delta  <- lavTech(lavobject, "delta")
    E.inv  <- lavTech(lavobject, "inverted.information")

    fg <- unlist(lavsamplestats@nobs)/lavsamplestats@ntotal

    # Fixme: as we only need the trace, perhaps we could do this
    # group-specific? (see lav_test_satorra_bentler_trace_original)
    V.g <- V; W.dwls.g <- W.dwls; Gamma.f <- Gamma; Delta.g <- Delta
    for(g in seq_len(lavdata@ngroups)) {
        ntotal <- nrow(Gamma[[g]])
        nvar <- lavobject@Model@nvar[[g]]
        pstar <- nvar * (nvar - 1) / 2
        rm.idx <- seq_len(ntotal - pstar)

        # reduce
        Delta.g[[g]]  <- Delta[[g]][-rm.idx,,drop=FALSE]
        # reduce and weight
        W.dwls.g[[g]] <- fg[g] * W.dwls[[g]][-rm.idx, -rm.idx]
        V.g[[g]]      <- fg[g] * V[[g]] # should already have the right dims
        Gamma.f[[g]]  <- 1/fg[g] * Gamma[[g]][-rm.idx, -rm.idx]
    }
    # create 'big' matrices
    W.dwls.all <- lav_matrix_bdiag(W.dwls.g)
    V.all      <- lav_matrix_bdiag(V.g)
    Gamma.all  <- lav_matrix_bdiag(Gamma.f)
    Delta.all  <- do.call("rbind", Delta.g)

    # compute trace
    WiU.all <- diag(nrow(W.dwls.all)) - Delta.all %*% E.inv %*% t(Delta.all) %*% W.dwls.all
    ks <- sum(diag(t(WiU.all) %*% V.all %*% WiU.all %*% Gamma.all))

    # convert to lavaan 'scaling.factor'
    c.hat3 <- ks/df3
    XX3.scaled <- XX3/c.hat3

    # baseline model
    XX3.null <- fit.catml@baseline$test[[1]]$stat
    df3.null <- fit.catml@baseline$test[[1]]$df

    kbs <- sum(diag(Gamma.all))
    c.hat3.null <- kbs/df3.null

    # return values
    list(XX3 = XX3, df3 = df3, c.hat3 = c.hat3, XX3.scaled = XX3.scaled,
         XX3.null = XX3.null, df3.null = df3.null, c.hat3.null = c.hat3.null)
}
