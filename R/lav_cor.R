# user-visible routine to
# compute polychoric/polyserial/... correlations
#
# YR 17 Sept 2013
#
# - YR 26 Nov 2013: big change - make it a wrapper around lavaan()
#                   estimator = "none" means two.step (starting values)

lavCor <- function(object,
                   # lav.data options
                   ordered    = NULL,
                   group      = NULL,
                   missing    = "listwise",
                   ov.names.x = NULL,
                   # lavaan options
                   se         = "none",
                   estimator  = "two.step",
                   # other options (for lavaan)
                   ...,
                   output = "cor") {

    # shortcut if object = lavaan object
    if(inherits(object, "lavaan")) {
        out <- lav_cor_output(object, output = output)
        return(out)
    }

    # check estimator
    estimator <- tolower(estimator)
    if(estimator %in% c("two.step", "two.stage")) {
        estimator <- "none"
    }

    # se?
    se <- tolower(se); output <- tolower(output)
    if(se != "none") {
        if(output %in% c("cor","cov","sampstat","th","thresholds")) {
            warning("lavaan WARNING: argument `se' is ignored since standard erros are not needed for the requested `output'")
            se <- "none"
        }
    }

    # check object class
    if(inherits(object, "lavData")) {
        lav.data <- object
    } else if(inherits(object, "data.frame")) {
        NAMES <- names(object)
        if(!is.null(group)) {
            NAMES <- NAMES[- match(group, NAMES)]
        }
        lav.data <- lavData(data = object, group = group,
                            ov.names = NAMES, ordered = ordered,
                            ov.names.x = ov.names.x,
                            lavoptions = list(missing = missing))
    } else {
        stop("lavaan ERROR: lavCor can not handle objects of class ",
             paste(class(object), collapse= " "))
    }

    # set default estimator if se != "none"
    categorical <- any(lav.data@ov$type == "ordered")
    if(se != "none" && estimator == "none") {
        if(categorical) {
            estimator <- "WLSMV"
        } else {
            estimator <- "ML"
        }
    }

    # extract partable options from dots
    dots <- list(...)
    meanstructure <- FALSE
    fixed.x       <- FALSE
    mimic         <- "lavaan"
    conditional.x <- FALSE
    if(!is.null(dots$meanstructure)) {
        meanstructure <- dots$meanstructure
    }
    if(categorical) {
        meanstructure <- TRUE
    }
    if(!is.null(dots$fixed.x)) {
        fixed.x <- dots$fixed.x
    }
    if(!is.null(dots$mimic)) {
        mimic <- dots$mimic
    }
    if(!is.null(dots$conditional.x)) {
        conditional.x <- dots$conditional.x
    }

    # override, only for backwards compatibility (eg moments() in JWileymisc)
    if(missing %in% c("ml", "fiml")) {
        meanstructure = TRUE
    }

    # generate partable for unrestricted model
    PT.un <-
        lav_partable_unrestricted(lavobject     = NULL,
                                  lavdata       = lav.data,
                                  lavoptions    = list(meanstructure = meanstructure,
                                                       fixed.x = fixed.x,
                                                       conditional.x = conditional.x,
                                                       mimic = mimic),
                                  sample.cov    = NULL,
                                  sample.mean   = NULL,
                                  sample.th     = NULL)


    FIT <- lavaan(slotParTable = PT.un, slotData = lav.data,
                  model.type = "unrestricted",
                  missing = missing,
                  se = se, estimator = estimator, ...)

    out <- lav_cor_output(FIT, output = output)

    out
}

lav_cor_output <- function(object, output = "cor") {

    # check output
    if(output %in% c("cor","cov")) {
        out <- inspect(object, "sampstat")
        if(object@Data@ngroups == 1L) {
            if(object@Model@conditional.x) {
                out <- out$res.cov
            } else {
                out <- out$cov
            }
            if(output == "cor") {
                out <- cov2cor(out)
            }
        } else {
            if(object@Model@conditional.x) {
                out <- lapply(out, "[[", "res.cov")
            } else {
                out <- lapply(out, "[[", "cov")
            }
            if(output == "cor") {
                out <- lapply(out, cov2cor)
            }
        }
    } else if(output %in% c("th","thresholds")) {
        out <- inspect(object, "sampstat")
        if(object@Data@ngroups == 1L) {
            if(object@Model@conditional.x) {
                out <- out$res.th
            } else {
                out <- out$th
            }
        } else {
            if(object@Model@conditional.x) {
                out <- lapply(out, "[[", "res.th")
            } else {
                out <- lapply(out, "[[", "th")
            }
        }
    } else if(output %in% c("sampstat")) {
        out <- inspect(object, "sampstat")
    } else if(output %in% c("parameterEstimates", "pe",
              "parameterestimates", "est")) {
        #out <- parameterEstimates(object)
        out <- standardizedSolution(object)
    } else {
        out <- object
    }

    out
}

