# initial version: YR 03/05/2017

# header
lav_object_print_header <- function(object) {
    
    # catch FAKE run
    FAKE <- FALSE
    if(object@Options$optim.method == "none") {
        FAKE <- TRUE
    }

    # Convergence or not?
    if(FAKE) {
        cat(sprintf("lavaan (%s) -- DRY RUN with 0 iterations\n",
                    packageDescription("lavaan", fields="Version")))
    } else if(object@optim$iterations > 0) {
        if(object@optim$converged) {
        cat(sprintf("lavaan (%s) converged normally after %3i iterations\n",
                    packageDescription("lavaan", fields="Version"),
                    object@optim$iterations))
        } else {
            cat(sprintf("** WARNING ** lavaan (%s) did NOT converge after %i iterations\n",
                packageDescription("lavaan", fields="Version"),
                object@optim$iterations))
            cat("** WARNING ** Estimates below are most likely unreliable\n")
        }
    } else {
        cat(sprintf("** WARNING ** lavaan (%s) model has NOT been fitted\n",
                    packageDescription("lavaan", fields="Version")))
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



