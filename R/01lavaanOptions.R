# initial version YR 02/08/2010

# this function collects and checks the user-provided options/arguments, 
# and fills in the "default" values, or changes them in an attempt to
# produce a consistent set of values...
#
# returns a list with the named options 

setLavaanOptions <- function(opt = formals(lavaan))
{
    # everything lowercase
    opt.old <- opt
    opt <- lapply(opt, function(x) { if(is.character(x)) tolower(x) else x})
    # except group.partial, which may contain capital letters
    opt$group.partial <- opt.old$group.partial

    # do.fit implies se="none and test="none" (unless not default)
    if(!opt$do.fit) {
        if(opt$se == "default") {
            opt$se <- "none"
        }
        if(opt$test == "default") {
            opt$test <- "none"
        }
    }

    # group.equal and group.partial
    if(is.null(opt$group.equal) || nchar(opt$group.equal) == 0L) {
        opt$group.equal <- character(0)
    } else if(length(opt$group.equal) == 0) {
        # nothing to do
    } else if(all(opt$group.equal %in% c("loadings", "intercepts", "means",
                                         "regressions", "residuals",
                                         "residual.covariances",
                                         "lv.variances", "lv.covariances"))) {
        # nothing to do 
    } else {
        stop("unknown value for `group.equal' argument: ",
             opt$group.equal, "\n")
    }
    if(is.null(opt$group.partial) || nchar(opt$group.partial) == 0L) {
        opt$group.partial <- character(0)
    } else if(length(opt$group.partial) == 0) {
        # nothing to do
    }

    # mimic
    if(opt$mimic == "default") {
        # WARNING: this will likely change soon
        # for now, we use mimic=Mplus as the default, but since
        # there are an increasing number of Mplus oddities, we will
        # make are own decisions in the future.
        opt$mimic <- "Mplus"
    } else if(opt$mimic == "mplus") {
        opt$mimic <- "Mplus"
    } else if(opt$mimic == "eqs") {
        opt$mimic <- "EQS"
    } else if(opt$mimic == "lisrel") {
        cat("Warning: mimic=\"LISREL\" is not ready yet. Using EQS instead.\n")
        opt$mimic <- "EQS"
    } else {
        stop("mimic must be \"Mplus\" or \"EQS\" \n")
    }

    # representation
    if(opt$representation == "default") {
        opt$representation <- "LISREL"
    } else if(opt$representation == "lisrel") {
        opt$representation <- "LISREL"
    } else if(opt$representation == "eqs" || 
              opt$representation == "bentler-weeks") {
        opt$representation <- "EQS"
    } else {
        stop("representation must be \"LISREL\" or \"EQS\" \n")
    }


    # missing
    if(opt$missing == "default") {
        #if(opt$mimic == "Mplus") { # since version 5?
        #    opt$missing <- "ml"
        #} else {
            opt$missing <- "listwise"
        #}
    } else if(opt$missing %in% c("ml", "direct", "fiml")) {
        opt$missing <- "ml"
    } else if(opt$missing %in% c("two.stage", "listwise")) {
        # nothing to do
    } else {
        stop("unknown value for `missing' argument: ", opt$missing, "\n")
    }

    if(opt$missing == "ml" && opt$data.type == "moment") {
        stop("value for `missing' is \"ml\" but only sample statistics are provided")
    }

    # default test statistic
    if(opt$test == "default") {
        opt$test <- "standard"
    } else if(opt$test %in% c("none", "standard", 
                              "satorra.bentler", "yuan.bentler",
                              "bollen.stine", "bootstrap")) {
        # nothing to do
    } else if(opt$test == "satorra" || opt$test == "satorra-bentler") {
        opt$test <- "satorra.bentler"
    } else if(opt$test == "yuan" || opt$test == "yuan-bentler") {
        opt$test <- "yuan.bentler"
    } else {
        stop("`test' argument must one of \"none\", \"standard\"",
             "\"satorra.bentler\", \"yuan.bentler\", \"bollen.stine\"",
             " or \"bootstrap\"")
    }

    # meanstructure
    if(is.logical(opt$meanstructure)) {
        if(opt$meanstructure == FALSE) {
            # user explicitly wants meanstructure == FALSE
            # check for conflicting arguments
            if(opt$estimator %in% c("mlm", "mlr", "mlf"))
                warning("lavaan WARNING: estimator forces meanstructure = TRUE")
            if(opt$missing == "ml")
                warning("lavaan WARNING: missing argument forces meanstructure = TRUE")
        }
    } else if(opt$meanstructure == "default") {
        # by default: no meanstructure!
        opt$meanstructure <- FALSE
    } else {
        stop("meanstructure must be TRUE, FALSE or \"default\"\n")
    }

    # default estimator and se
    if(opt$se == "bootstrap") {
        stop("bootstrapped standard errors are not implemented yet!\n")
    }

    if(opt$estimator == "default" || opt$estimator == "ml") {
        opt$estimator <- "ML"
        if(opt$se == "default") {
            opt$se <- "standard"
        } else if(opt$se == "first.order" || 
                  opt$se == "none"        || 
                  opt$se == "standard"    ||
                  opt$se == "robust.mlr"  || 
                  opt$se == "robust.mlm") {
            # nothing to do
        } else if(opt$se == "robust") {
            if(opt$missing == "ml") {
                opt$se <- "robust.mlr"
            } else {
                opt$se <- "robust.mlm"
            }
        } else {
            stop("unknown value for `se' argument when estimator is ML: ", 
                 opt$se, "\n")
        }
    } else if(opt$estimator == "mlm") {
        opt$estimator <- "ML"
        opt$information <- "expected"
        opt$meanstructure <- TRUE
        opt$se <- "robust.mlm"
        opt$test <- "satorra.bentler"
        if(opt$missing == "ml") {
            stop("the MLM estimator can not be used when data are incomplete\n")
        }
    } else if(opt$estimator == "mlf") {
        opt$estimator <- "ML"
        opt$meanstructure <- TRUE
        opt$se <- "first.order"
    } else if(opt$estimator == "mlr") {
        opt$estimator <- "ML"
        opt$se <- "robust.mlr"
        opt$test <- "yuan.bentler"
        opt$meanstructure <- TRUE
    } else if(opt$estimator == "gls") {
        opt$estimator <- "GLS"
        if(opt$se == "default" || opt$se == "standard") {
            opt$se <- "standard"
        } else if(opt$se == "none") {
            # nothing to do
        } else {
            stop("invalid value for `se' argument when estimator is GLS: ", 
                 opt$se, "\n")
        }
        if(opt$test != "standard") {
            stop("invalid value for `test' argument when estimator is GLS: ", 
                 opt$test, "\n")
        }
    } else if(opt$estimator == "wls") {
        opt$estimator <- "WLS"
        if(opt$se == "default" || opt$se == "standard") {
            opt$se <- "standard"
        } else if(opt$se == "none") {
            # nothing to do
        } else {
            stop("invalid value for `se' argument when estimator is WLS: ", 
                 opt$se, "\n")
        }
        if(opt$test != "standard") {
            stop("invalid value for `test' argument when estimator is GLS: ", 
                 opt$test, "\n")
        }
    } else if(opt$estimator == "uls") {
        opt$estimator <- "none"
        if(opt$se == "default" || opt$se == "standard") {
            opt$se <- "standard"
        } else if(opt$se == "none") {
            # nothing to do
        } else {
            stop("invalid value for `se' argument when estimator is ULS: ", 
                 opt$se, "\n")
        }
        if(opt$test != "standard") {
            stop("invalid value for `test' argument when estimator is ULS: ",
                 opt$test, "\n")
        }
    } else {
        stop("unknown value for `estimator' argument: ", opt$estimator, "\n")
    }

    # likelihood approach (wishart or normal)
    if(opt$estimator != "ML") {
        if(opt$likelihood != "default") {
            stop("likelihood argument is only relevant if estimator = ML")
        } 
    } else { # ml
        if(opt$likelihood == "default") {
           opt$likelihood <- "normal"
            if(opt$mimic == "EQS"    || 
               opt$mimic == "LISREL" || 
               opt$mimic == "AMOS") {
                opt$likelihood <- "wishart"
            }
        } else if(opt$likelihood == "wishart" || opt$likelihood == "normal") {
            # nothing to do
        } else {
            stop("invalid value for `likelihood' argument: ", 
                 opt$likelihood, "\n")
        }
    }

    # information
    if(opt$information == "default") {
        if(opt$missing == "ml"     || 
           opt$se == "robust.mlr"  || 
           opt$se == "first.order" ||
           nchar(opt$constraints) > 0L) {
            opt$information <- "observed"
        } else {
            opt$information <- "expected"
        }
    } else if(opt$information %in% c("observed", "expected")) {
        # nothing to do
    } else {
        stop("information must be either \"expected\" or \"observed\"\n")
    }

    # fixed.x
    if(is.logical(opt$fixed.x)) {
        # nothing to do
    } else if(opt$fixed.x == "default") {
        if(opt$estimator == "ML" && opt$mimic == "Mplus") {
            opt$fixed.x <- TRUE
        } else {
            opt$fixed.x <- FALSE
        }
    } else {
        stop("fixed.x must be TRUE, FALSE or \"default\"\n")
    }


    # meanstructure
    if(opt$missing == "ml" || opt$model.type == "growth") {
        opt$meanstructure <- TRUE
    }
    if("intercepts" %in% opt$group.equal ||
       "means" %in% opt$group.equal) {
        opt$meanstructure <- TRUE
    }
    stopifnot(is.logical(opt$meanstructure))
    stopifnot(is.logical(opt$verbose))
    stopifnot(is.logical(opt$warn))

    if(opt$debug) {
        opt$verbose <- opt$warn <- TRUE
    }

    if(opt$debug) { cat("lavaan DEBUG: lavaanOptions OUT\n"); str(opt) }

    opt
}
