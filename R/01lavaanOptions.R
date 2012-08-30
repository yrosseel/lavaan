# initial version YR 02/08/2010

# this function collects and checks the user-provided options/arguments, 
# and fills in the "default" values, or changes them in an attempt to
# produce a consistent set of values...
#
# returns a list with the named options 

setLavaanOptions <- function(opt = formals(lavaan))
{

    if(opt$debug) { cat("lavaan DEBUG: lavaanOptions IN\n"); str(opt) }

    # everything lowercase
    opt.old <- opt
    opt <- lapply(opt, function(x) { if(is.character(x)) tolower(x) else x})
    # except group,group.partial, which may contain capital letters
    opt$group <- opt.old$group
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

    # mimic
    if(opt$mimic == "default" || opt$mimic == "lavaan") {
        opt$mimic <- "lavaan"
    } else if(opt$mimic == "mplus") {
        opt$mimic <- "Mplus"
    } else if(opt$mimic == "eqs") {
        opt$mimic <- "EQS"
    } else if(opt$mimic == "lisrel") {
        cat("Warning: mimic=\"LISREL\" is not ready yet. Using EQS instead.\n")
        opt$mimic <- "EQS"
    } else {
        stop("mimic must be \"lavaan\", \"Mplus\" or \"EQS\" \n")
    }

    # group.equal and group.partial
    if(opt$group.equal[1] == "none") {
        opt$group.equal <- character(0)
    } else if(is.null(opt$group.equal) || nchar(opt$group.equal) == 0L) {
        if(opt$mimic == "Mplus" && !is.null(opt$group)) {
            if(opt$categorical) {
                opt$group.equal <- c("loadings", "thresholds")
            } else {
                opt$group.equal <- c("loadings", "intercepts")
            }
        } else {
            opt$group.equal <- character(0)
        }
    } else if(length(opt$group.equal) == 0) {
        # nothing to do
    } else if(all(opt$group.equal %in% c("loadings", "intercepts", "means",
                                         "regressions", "residuals",
                                         "residual.covariances", "thresholds",
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
    } else {
        # strip white space
        opt$group.partial <- gsub("[[:space:]]+", "", opt$group.partial)
    }

    # if categorical, and group.equal contains "intercepts", also add
    # thresholds (and vice versa)
    if(opt$categorical && "intercepts" %in% opt$group.equal) {
        opt$group.equal <- unique(c(opt$group.equal, "thresholds"))
    }
    if(opt$categorical && "thresholds" %in% opt$group.equal) {
        opt$group.equal <- unique(c(opt$group.equal, "intercepts"))
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
        if(opt$mimic == "Mplus" && 
           opt$estimator %in% c("default", "ml", "mlr")) { 
            # since version 5?
            opt$missing <- "ml" 
            # check later if this is ok
        } else {
            opt$missing <- "listwise"
        }
    } else if(opt$missing %in% c("ml", "direct", "fiml")) {
        opt$missing <- "ml"
        if(opt$estimator %in% c("mlm", "mlmv", "gls", "wls", "wlsm", "wlsmv",
                                "uls", "ulsm", "ulsmv")) {
            stop("lavaan ERROR: missing=\"ml\" is not allowed for estimator MLM, MLMV, GLS, ULS, ULSM, ULSMV, DWLS, WLS, WLSM, WLSMV")
        }
    } else if(opt$missing %in% c("two.stage", "listwise")) {
        # nothing to do
    } else {
        stop("unknown value for `missing' argument: ", opt$missing, "\n")
    }

    # default test statistic
    if(opt$test == "default") {
        opt$test <- "standard"
    } else if(opt$test %in% c("none", "standard")) {
        # nothing to do
    } else if(opt$test == "satorra" || 
              opt$test == "sb" || 
              opt$test == "SB" ||
              opt$test == "satorra.bentler" ||
              opt$test == "satorra-bentler") {
        opt$test <- "satorra.bentler"
    } else if(opt$test == "yuan" || 
              opt$test == "yb" ||
              opt$test == "YB" ||
              opt$test == "yuan.bentler" ||
              opt$test == "yuan-bentler") {
        opt$test <- "yuan.bentler"
    } else if(opt$test == "m.adjusted" ||
              opt$test == "m" ||
              opt$test == "mean.adjusted" ||
              opt$test == "mean-adjusted") {
        opt$test <- "satorra.bentler"
    } else if(opt$test == "mean.var.adjusted" ||
              opt$test == "mean-var-adjusted" ||
              opt$test == "mv" ||
              opt$test == "second.order" ||
              opt$test == "satterthwaite" ||
              opt$test == "Satterthwaite" ||
              opt$test == "mv.adjusted") {
        opt$test <- "mean.var.adjusted"
    } else if(opt$test == "mplus6" ||
              opt$test == "scale.shift" ||
              opt$test == "scaled.shifted") {
        opt$test <- "scaled.shifted"
    } else if(opt$test == "bootstrap" || 
              opt$test == "boot" ||
              opt$test == "bollen.stine" || 
              opt$test == "bollen-stine") {
        opt$test <- "bollen.stine"
    } else {
        stop("`test' argument must one of \"none\", \"standard\",
            \"satorra.bentler\", \"yuan.bentler\",
            \"mean.var.adjusted\", \"scaled.shifted\",
            \"bollen.stine\", or \"bootstrap\"")
    }
 
    # check missing
    if(opt$missing == "ml" && 
       opt$test %in% c("satorra.bentler", 
                       "mean.var.adjusted", "scaled.shifted")) {
        opt$missing <- "listwise"
    }

    # meanstructure
    if(is.logical(opt$meanstructure)) {
        if(opt$meanstructure == FALSE) {
            # user explicitly wants meanstructure == FALSE
            # check for conflicting arguments
            if(opt$estimator %in% c("mlm", "mlmv", "mlr", "mlf", "ulsm", "ulsmv", "wlsm", "wlsmv"))
                warning("lavaan WARNING: estimator forces meanstructure = TRUE")
            if(opt$missing == "ml")
                warning("lavaan WARNING: missing argument forces meanstructure = TRUE")
        }
    } else if(opt$meanstructure == "default") {
        # by default: no meanstructure!
        opt$meanstructure <- FALSE
        # unless there is a group argument? (added since 0.4-10)
        if(!is.null(opt$group)) opt$meanstructure <- TRUE
    } else {
        stop("meanstructure must be TRUE, FALSE or \"default\"\n")
    }

    # estimator and se
    if(opt$se == "boot" || opt$se == "bootstrap") {
        opt$se <- "bootstrap"
        opt$information <- "observed"
        opt$bootstrap <- as.integer(opt$bootstrap)
        stopifnot(opt$bootstrap > 0L)
    }

    # default estimator
    if(opt$estimator == "default") {
        if(opt$categorical) 
            opt$estimator <- "wlsmv"
        else
            opt$estimator <- "ml"
    }

    if(opt$estimator == "ml") {
        opt$estimator <- "ML"
        if(opt$se == "default") {
            opt$se <- "standard"
        } else if(opt$se == "first.order" || 
                  opt$se == "bootstrap"   ||
                  opt$se == "none"        || 
                  opt$se == "standard"    ||
                  opt$se == "robust.huber.white"  || 
                  opt$se == "robust.sem") {
            # nothing to do
        } else if(opt$se == "robust") {
            if(opt$missing == "ml") {
                opt$se <- "robust.huber.white"
            } else {
                opt$se <- "robust.sem"
            }
        } else {
            stop("unknown value for `se' argument when estimator is ML: ", 
                 opt$se, "\n")
        }

    } else if(opt$estimator == "mlm" || opt$estimator == "mlmv") {
        if(opt$test != "none") {
            if(opt$estimator == "mlm") {
                opt$test <- "satorra.bentler"
            } else if(opt$estimator == "mlmv") {          
                opt$test <- "scaled.shifted"
            }
        }
        opt$estimator <- "ML"
        opt$information <- "expected"
        opt$meanstructure <- TRUE
        if(opt$se == "bootstrap") stop("use ML estimator for bootstrap")
        if(opt$se != "none") opt$se <- "robust.sem"
        opt$missing <- "listwise"
    } else if(opt$estimator == "mlf") {
        opt$estimator <- "ML"
        opt$meanstructure <- TRUE
        if(opt$se == "bootstrap") stop("use ML estimator for bootstrap")
        if(opt$se != "none") opt$se <- "first.order"
    } else if(opt$estimator == "mlr") {
        opt$estimator <- "ML"
        if(opt$se == "bootstrap") stop("use ML estimator for bootstrap")
        if(opt$se != "none") opt$se <- "robust.huber.white"
        if(opt$test != "none") opt$test <- "yuan.bentler"
        opt$meanstructure <- TRUE
    } else if(opt$estimator == "gls") {
        opt$estimator <- "GLS"
        if(opt$se == "default" || opt$se == "standard") {
            opt$se <- "standard"
        } else if(opt$se == "none" || opt$se == "bootstrap") {
            # nothing to do
        } else {
            stop("invalid value for `se' argument when estimator is GLS: ", 
                 opt$se, "\n")
        }
        if(!opt$test %in% c("standard","none")) {
            stop("invalid value for `test' argument when estimator is GLS: ", 
                 opt$test, "\n")
        }
        opt$missing <- "listwise"       
    } else if(opt$estimator == "wls") {
        opt$estimator <- "WLS"
        if(opt$se == "default" || opt$se == "standard") {
            opt$se <- "standard"
        } else if(opt$se == "none" || opt$se == "bootstrap") {
            # nothing to do
        } else if(opt$se == "robust.sem") {
            # nothing to do
        } else if(opt$se == "robust") {
            opt$se <- "robust.sem"
        } else {
            stop("invalid value for `se' argument when estimator is WLS: ", 
                 opt$se, "\n")
        }
        if(!opt$test %in% c("standard","none")) {
            stop("invalid value for `test' argument when estimator is WLS: ", 
                 opt$test, "\n")
        }
        opt$missing <- "listwise"
    } else if(opt$estimator == "dwls") {
        opt$estimator <- "DWLS"
        if(opt$se == "default" || opt$se == "standard") {
            opt$se <- "standard"
        } else if(opt$se == "none" || opt$se == "bootstrap") {
            # nothing to do
        } else if(opt$se == "robust.sem") {
            # nothing to do
        } else if(opt$se == "robust") {
            opt$se <- "robust.sem"
        } else {
            stop("invalid value for `se' argument when estimator is DWLS: ",
                 opt$se, "\n")
        }
        if(!opt$test %in% c("standard","none","satorra.bentler", 
                            "mean.adjusted",
                            "mean.var.adjusted","scaled.shifted")) {
            stop("invalid value for `test' argument when estimator is DWLS: ",
                 opt$test, "\n")
        }
        opt$missing <- "listwise"
    } else if(opt$estimator == "wlsm") {
        opt$estimator <- "DWLS"
        if(opt$se == "bootstrap") stop("use (D)WLS estimator for bootstrap")
        if(opt$se != "none") opt$se <- "robust.sem"
        if(opt$test != "none") opt$test <- "satorra.bentler"
        opt$missing <- "listwise"
     } else if(opt$estimator == "wlsmv") {
        opt$estimator <- "DWLS"
        if(opt$se == "bootstrap") stop("use (D)WLS estimator for bootstrap")
        if(opt$se != "none") opt$se <- "robust.sem"
        if(opt$test != "none") opt$test <- "scaled.shifted"
        opt$missing <- "listwise"
    } else if(opt$estimator == "uls") {
        opt$estimator <- "ULS"
        if(opt$se == "default" || opt$se == "standard") {
            opt$se <- "standard"
        } else if(opt$se == "none" || opt$se == "bootstrap") {
            # nothing to do
        } else if(opt$se == "robust.sem") {
            # nothing to do
        } else if(opt$se == "robust") {
            opt$se <- "robust.sem"
        } else {
            stop("invalid value for `se' argument when estimator is ULS: ", 
                 opt$se, "\n")
        }
        if(!opt$test %in% c("standard","none", "satorra.bentler",
                            "mean.adjusted",
                            "mean.var.adjusted","scaled.shifted")) {
            stop("invalid value for `test' argument when estimator is ULS: ",
                 opt$test, "\n")
        }
        opt$missing <- "listwise"
    } else if(opt$estimator == "ulsm") {
        opt$estimator <- "ULS"
        if(opt$se == "bootstrap") stop("use ULS estimator for bootstrap")
        if(opt$se != "none") opt$se <- "robust.sem"
        if(opt$test != "none") opt$test <- "satorra.bentler"
        opt$missing <- "listwise"
    } else if(opt$estimator == "ulsmv") {
        opt$estimator <- "ULS"
        if(opt$se == "bootstrap") stop("use ULS estimator for bootstrap")
        if(opt$se != "none") opt$se <- "robust.sem"
        if(opt$test != "none") opt$test <- "scaled.shifted"
        opt$missing <- "listwise"
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
           opt$se == "robust.huber.white"  || 
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
        if(opt$estimator == "ML" && (opt$mimic == "Mplus" || 
                                     opt$mimic == "lavaan")) {
            opt$fixed.x <- TRUE
        } else if(opt$categorical) {
            opt$fixed.x <- TRUE
        } else {
            opt$fixed.x <- FALSE
        }
    } else {
        stop("fixed.x must be TRUE, FALSE or \"default\"\n")
    }


    # meanstructure again
    if(opt$missing == "ml" || opt$model.type == "growth") {
        opt$meanstructure <- TRUE
    }
    if("intercepts" %in% opt$group.equal ||
       "means" %in% opt$group.equal) {
        opt$meanstructure <- TRUE
    }
    if(opt$se == "robust.huber.white" || 
       opt$se == "robust.sem" ||
       opt$test == "satorra.bentler" ||
       opt$test == "mean.var.adjusted" ||
       opt$test == "scaled.shifted" ||
       opt$test == "yuan.bentler") {
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
