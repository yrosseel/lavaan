# initial version: YR 03/05/2017

# header
lav_object_print_header <- function(object) {

    cat(sprintf("lavaan %s ",
                packageDescription("lavaan", fields="Version")))

    # catch FAKE run
    FAKE <- FALSE
    if(object@Options$optim.method == "none") {
        FAKE <- TRUE
    }

    # Convergence or not?
    if(FAKE) {
        cat("-- DRY RUN with 0 iterations --\n")
    } else if(object@optim$iterations > 0) {
        if(object@optim$converged) {
        cat(sprintf("ended normally after %i iterations\n",
                    object@optim$iterations))
        } else {
            cat(sprintf("did NOT end normally after %i iterations\n",
                object@optim$iterations))
            cat("** WARNING ** Estimates below are most likely unreliable\n")
        }
    } else {
        cat("did not run (perhaps do.fit = FALSE)?\n")
        cat("** WARNING ** Estimates below are simply the starting values\n")
    }
    cat("\n")

    # number of free parameters
    #t0.txt <- sprintf("  %-40s", "Number of free parameters")
    #t1.txt <- sprintf("  %10i", object@optim$npar)
    #t2.txt <- ""
    #cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
    #cat("\n")
}

# optim
lav_object_print_optim <- function(object) {

    #cat("Optimization information:\n\n")

    t0.txt <- sprintf("  %-40s", "Optimization method")
    t1.txt <- sprintf("  %10s", toupper(object@Options$optim.method))
    t2.txt <- ""
    cat(t0.txt, t1.txt, t2.txt, "\n", sep="")

    t0.txt <- sprintf("  %-40s", "Number of free parameters")
    t1.txt <- sprintf("  %10i",   object@optim$npar)
    t2.txt <- ""
    cat(t0.txt, t1.txt, t2.txt, "\n", sep="")

    if(object@Model@eq.constraints) {
        t0.txt <- sprintf("  %-40s", "Number of equality constraints")
        t1.txt <- sprintf("  %10i", nrow(object@Model@ceq.JAC))
        t2.txt <- ""
        cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
    }
    if(nrow(object@Model@cin.JAC) > 0L) {
        t0.txt <- sprintf("  %-40s", "Number of inequality constraints")
        t1.txt <- sprintf("  %10i", nrow(object@Model@cin.JAC))
        t2.txt <- ""
        cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
    }

    cat("\n")
}

# rotation
lav_object_print_rotation <- function(object) {

    #cat("Rotation information:\n\n")

    t0.txt <- sprintf("  %-20s", "Rotation method")
    if(object@Options$rotation.args$orthogonal) {
        MM <- paste(toupper(object@Options$rotation), " ", "ORTHOGONAL",
                    sep = "")
    } else {
        MM <- paste(toupper(object@Options$rotation), " ", "OBLIQUE",
                    sep = "")
    }
    t1.txt <- sprintf("  %30s", MM)
    t2.txt <- ""
    cat(t0.txt, t1.txt, t2.txt, "\n", sep="")


    if(object@Options$rotation == "geomin") {
        t0.txt <- sprintf("  %-40s", "Geomin epsilon")
        t1.txt <- sprintf("  %10.6g",
            object@Options$rotation.args$geomin.epsilon)
        cat(t0.txt, t1.txt, "\n", sep="")
    } else if(object@Options$rotation == "orthomax") {
        t0.txt <- sprintf("  %-40s", "Orthomax gamma")
        t1.txt <- sprintf("  %10.6g",
            object@Options$rotation.args$orthomax.gamma)
        cat(t0.txt, t1.txt, "\n", sep="")
    } else if(object@Options$rotation == "cf") {
        t0.txt <- sprintf("  %-40s", "Crawford-Ferguson gamma")
        t1.txt <- sprintf("  %10.6g",
            object@Options$rotation.args$cf.gamma)
        cat(t0.txt, t1.txt, "\n", sep="")
    } else if(object@Options$rotation == "oblimin") {
        t0.txt <- sprintf("  %-40s", "Oblimin gamma")
        t1.txt <- sprintf("  %10.6g",
            object@Options$rotation.args$oblimin.gamma)
        cat(t0.txt, t1.txt, "\n", sep="")
    }

    t0.txt <- sprintf("  %-30s", "Rotation algorithm (rstarts)")
    tmp <- paste(toupper(object@Options$rotation.args$algorithm),
                 " (", object@Options$rotation.args$rstarts, ")", sep = "")
    t1.txt <- sprintf("  %20s", tmp)
    cat(t0.txt, t1.txt, t2.txt, "\n", sep="")

    t0.txt <- sprintf("  %-40s", "Standardized metric")
    if(object@Options$rotation.args$std.ov) {
        t1.txt <- sprintf("  %10s", "TRUE")
    } else {
        t1.txt <- sprintf("  %10s", "FALSE")
    }
    t2.txt <- ""
    cat(t0.txt, t1.txt, t2.txt, "\n", sep="")

    t0.txt <- sprintf("  %-40s", "Row weights")
    tmp.txt <- object@Options$rotation.args$row.weights
    t1.txt <- sprintf("  %10s", paste(toupper(substring(tmp.txt, 1, 1)),
                                      substring(tmp.txt, 2), sep = ""))
    t2.txt <- ""
    cat(t0.txt, t1.txt, t2.txt, "\n", sep="")

    cat("\n")
}


