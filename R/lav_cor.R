# user-visible routine to 
# compute polychoric/polyserial/... correlations
#
# YR 17 Sept 2013
# 
# - YR 26 Nov 2013: big change - make it wrapper around lavaan()
#                   estimator = "none" means two.step (starting values)

lavCor <- function(object, 
                   # lav.data options
                   ordered    = NULL, 
                   group      = NULL, 
                   missing    = "listwise",
                   ov.names.x = NULL, 
                   # partable options
                   meanstructure = FALSE,
                   fixed.x       = TRUE,
                   # lavaan options
                   se         = "none", 
                   estimator  = "two.step", 
                   # other options (for lavaan)
                   ...,
                   return.type = "cor") {

    # check estimator
    estimator <- tolower(estimator)
    if(estimator %in% c("two.step", "two.stage")) {
        estimator <- "none"
    }

    # check object class
    if(inherits(object, "lavaan")) {
        lav.data <- object@Data
    } else if(inherits(object, "lavData")) {
        lav.data <- object
    } else if(inherits(object, "data.frame")) {
        NAMES <- names(object)
        if(!is.null(group)) {
            NAMES <- NAMES[- match(group, NAMES)]
        }
        lav.data <- lavData(data = object, group = group, 
                            ov.names = NAMES, ordered = ordered,
                            ov.names.x = ov.names.x,
                            missing = missing)
    } else {
        stop("lavaan ERROR: lavCor can not handle objects of class ", 
             class(object))
    }

    # generate partable for unrestricted model
    PT.un <- 
        lav_partable_unrestricted(ov.names      = lav.data@ov.names,
                                  ov            = lav.data@ov,
                                  ov.names.x    = lav.data@ov.names.x,
                                  sample.cov    = NULL,
                                  meanstructure = meanstructure,
                                  sample.mean   = NULL,
                                  sample.th     = NULL,
                                  fixed.x       = fixed.x)

    fit <- lavaan(slotParTable = PT.un, slotData = lav.data,
                  model.type = "unrestricted",
                  se = se, estimator = estimator, ...)

    # check return.type
    if(return.type %in% c("cor","cov")) {
        out <- inspect(fit, "sampstat")
        if(fit@Data@ngroups == 1L) {
            out <- out$cov
        } else {
            out <- lapply(out, "[[", "cov")
        }
    } else if(return.type %in% c("sampstat")) {
        out <- inspect(fit, "sampstat")
    } else if(return.type %in% c("parameterEstimates", "pe", 
              "parameterestimates", "est")) {
        out <- parameterEstimates(fit)
    } else {
        out <- fit
    }

    out
}

