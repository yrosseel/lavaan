# EFA: exploratory factor analysis
#
# EFA is implemented as a special version of ESEM
# - it is therefore a wrapper around the lavaan() function to simplify
#   the input
# - a lavaan model is generated with a single 'block' that can be rotated
# - the 'default' output produces output that is more in line with traditional
#   EFA software (in R) like factanal() and fa() from the psych package

# YR 20 Sept 2022 - first version

efa <- function(data           = NULL,
                nfactors       = 1L,
                sample.cov     = NULL,
                sample.nobs    = NULL,
                rotation       = "geomin",
                rotation.args  = list(),
                ov.names       = names(data),
                se             = "none",
                ...,
                output         = "lavaan.efa") {

    # handle ov.names
    if(!is.null(data) && inherits(data, "data.frame")) {
        if(length(ov.names) > 0L) {
            data <- data[, ov.names, drop = FALSE]
        } else {
            ov.names <- names(data)
        }
    }

    # unless it is specified in rotation.args, set order.lv.by to "sumofsquares"
    if(is.null(rotation.args$order.lv.by)) {
        rotation.args$order.lv.by = "sumofsquares"
    }

    # output
    output <- tolower(output)
    if(!output %in% c("lavaan", "efa", "lavaan.efa")) {
        stop("lavaan ERROR: output= must be either \"lavaan\" or \"lavaan.efa\"")
    }
    if(output == "efa") {
        output <- "lavaan.efa"
    }

    # handle dotdotdot
    dotdotdot <- list(...)

    # check for unallowed arguments
    if(!is.null(dotdotdot$group)) {
        stop("lavaan ERROR: efa has no support for multiple groups (for now)")
    }
    if(!is.null(dotdotdot$cluster)) {
        stop("lavaan ERROR: efa has no support for clustered data (for now)")
    }
    if(!is.null(dotdotdot$ordered)) {
        stop("lavaan ERROR: efa has no support for ordered data (for now)")
    }

    # switch off unneeded actions
    if(is.null(dotdotdot$baseline)) {
        dotdotdot$baseline <- FALSE
    }
    if(is.null(dotdotdot$loglik)) {
        dotdotdot$loglik <- FALSE
    }
    # implied is needed for lavResiduals -> fitMeasures
    # h1? needed for missing?
    # check.post?
    # ...


    # generate parameter table for an EFA model
    # note: we generate the full parameter table; perhaps it is more efficient
    #       to just creat a FLAT table (=output of  lavParseModelString)
    #       and add model=FLAT to the lavaan call?
    PT.efa <- lav_partable_generate_efa(ov.names = ov.names,
                                        nfactors = nfactors,
                                        ordered  = NULL)

    # call lavaan
    FIT <- do.call("lavaan",
                   args = c(list(slotParTable  = PT.efa,
                                 data          = data,
                                 sample.cov    = sample.cov,
                                 sample.nobs   = sample.nobs,
                                 rotation      = rotation,
                                 rotation.args = rotation.args),
                            dotdotdot))

    if(output == "lavaan") {
        out <- FIT
    } else if(output == "lavaan.efa") {
        # create a list
        out <- list(version = FIT@version,
                    lavoptions = FIT@Options,
                    lavpartable = FIT@ParTable,
                    standardized = standardizedSolution(FIT),
                    #lavsamplestats = FIT@SampleStats,
                    fit = fitMeasures(FIT),
                    glist.std = lavInspect(FIT, "std"),
                    converged = lavInspect(FIT, "converged"))
        class(out) <- c("lavaan.efa", "list")
    }

    out
}

print.lavaan.efa <- function(object, nd = 3L, cutoff = 0.3,
                             dot.cutoff = 0.1, ...) {
    if(!object$converged) {
        cat("** WARNING ** Optimizer did not end normally\n")
        cat("** WARNING ** Estimates below are most likely unreliable\n")
    }
    cat("\n")
    LAMBDA <- unclass(object$glist.std$lambda)
    lav_print_loadings(LAMBDA, nd = nd, cutoff = cutoff,
                       dot.cutoff = dot.cutoff)
    cat("\n")

    invisible(LAMBDA)
}

# more verbose output
summary.lavaan.efa <- function(object, nd = 3L, cutoff = 0.3,
                               dot.cutoff = 0.1, ...) {

    cat("This is ",
                sprintf("lavaan %s", object$version),
                " -- exploratory factor analysis\n", sep = "")
    if(!object$converged) {
        cat("** WARNING ** Optimizer did not end normally\n")
        cat("** WARNING ** Estimates below are most likely unreliable\n")
    }
    cat("\n")

    # fit measures
    print(round(object$fit[c("chisq", "df", "pvalue", "cfi",
                             "rmsea", "srmr")], nd))
    cat("\n")

    cat("Standardized factor loadings:\n\n")
    LAMBDA <- unclass(object$glist.std$lambda)
    THETA <- unname(unclass(object$glist.std$theta))
    lav_print_loadings(LAMBDA, nd = nd, cutoff = cutoff,
                       dot.cutoff = dot.cutoff,
                       resvar = diag(THETA))
    cat("\n")

    sumsq <- colSums(LAMBDA * LAMBDA)
    propvar <- sumsq/nrow(LAMBDA)
    cumvar <- cumsum(propvar)
    var.table <- rbind(sumsq, propvar, cumvar)
    #colnames(var.table) <- colnames(LAMBDA)
    rownames(var.table) <- c("Sum of squared loadings",
                             "Proportion Var",
                             "Cumulative Var")
    class(var.table) <- c("lavaan.matrix", "matrix", "array")
    print(var.table, nd = nd)
    cat("\n")

    # factor correlations:
    if(ncol(LAMBDA) > 1L && !object$lavoptions$rotation.args$orthogonal) {
        cat("Factor correlations:\n\n")
        print(object$glist.std$psi)
        cat("\n")
    }

    invisible(object)
}
