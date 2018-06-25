# numeric approximation of the Hessian
# using an analytic gradient
lav_model_hessian <- function(lavmodel       = NULL,
                              lavsamplestats = NULL,
                              lavdata        = NULL,
                              lavoptions     = NULL,
                              lavcache       = NULL,
                              group.weight   = TRUE) {

    estimator <- lavmodel@estimator

    # computing the Richardson extrapolation
    Hessian <- matrix(0, lavmodel@nx.free, lavmodel@nx.free)
    x <- lav_model_get_parameters(lavmodel = lavmodel)
    for(j in 1:lavmodel@nx.free) {
        h.j <- 1e-06
        x.left <- x.left2 <- x.right <- x.right2 <- x
        x.left[j]  <- x[j] - h.j; x.left2[j]  <- x[j] - 2*h.j
        x.right[j] <- x[j] + h.j; x.right2[j] <- x[j] + 2*h.j

        g.left <-
            lav_model_gradient(lavmodel       = lavmodel,
                               GLIST          = lav_model_x2GLIST(lavmodel =
                                                            lavmodel, x.left),
                               lavsamplestats = lavsamplestats,
                               lavdata        = lavdata,
                               lavcache       = lavcache,
                               type           = "free",
                               group.weight   = group.weight)
        g.left2 <-
            lav_model_gradient(lavmodel       = lavmodel,
                               GLIST          = lav_model_x2GLIST(lavmodel =
                                                            lavmodel, x.left2),
                               lavsamplestats = lavsamplestats,
                               lavdata        = lavdata,
                               lavcache       = lavcache,
                               type           = "free",
                               group.weight   = group.weight)

        g.right <-
            lav_model_gradient(lavmodel       = lavmodel,
                               GLIST          = lav_model_x2GLIST(lavmodel =
                                                            lavmodel, x.right),
                               lavsamplestats = lavsamplestats,
                               lavdata        = lavdata,
                               lavcache       = lavcache,
                               type           = "free",
                               group.weight   = group.weight)

        g.right2 <-
            lav_model_gradient(lavmodel       = lavmodel,
                               GLIST          = lav_model_x2GLIST(lavmodel =
                                                            lavmodel, x.right2),
                               lavsamplestats = lavsamplestats,
                               lavdata        = lavdata,
                               lavcache       = lavcache,
                               type           = "free",
                               group.weight   = group.weight)

        Hessian[,j] <- (g.left2 - 8*g.left + 8*g.right - g.right2)/(12*h.j)
    }

    # make symmetric (NEEDED? probably not)
    Hessian <- ( Hessian + t(Hessian) )/2.0

    Hessian
}

# if only chol would accept a complex matrix...
lav_model_hessian_complex <- function(lavmodel       = NULL,
                                      lavsamplestats = NULL,
                                      lavdata        = NULL,
                                      lavcache       = NULL,
                                      group.weight   = TRUE) {

    gradf <- function(x) {
        GLIST <-  lav_model_x2GLIST(lavmodel = lavmodel, x = x)
        dx <- lav_model_gradient(lavmodel       = lavmodel,
                                 GLIST          = GLIST,
                                 lavsamplestats = lavsamplestats,
                                 lavdata        = lavdata,
                                 lavcache       = lavcache,
                                 type           = "free",
                                 group.weight   = group.weight)
        dx
    }

    x <- lav_model_get_parameters(lavmodel = lavmodel)
    Hessian <- lav_func_jacobian_complex(func = gradf, x = x)

    Hessian
}
