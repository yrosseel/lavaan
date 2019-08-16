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

}

# optim
lav_object_print_optim <- function(object, nd = 3L) {

    #cat("Optimization information:\n\n")

    # first column
    c1 <- c("Estimator", "Optimization method", "Number of free parameters")

    # second column
    c2 <- c(toupper(object@Options$estimator),
            toupper(object@Options$optim.method), object@optim$npar)

    # optional output
    if(object@Model@eq.constraints) {
        c1 <- c(c1, "Number of equality constraints")
        c2 <- c(c2, nrow(object@Model@ceq.JAC))
    }
    if(nrow(object@Model@cin.JAC) > 0L) {
        c1 <- c(c1, "Number of inequality constraints")
        c2 <- c(c2, nrow(object@Model@cin.JAC))
    }
    if(nrow(object@Model@con.jac) > 0L) {
        c1 <- c(c1, "Row rank of the constraints matrix")
        c2 <- c(c2, qr(object@Model@con.jac)$rank)
    }

    # empty last row
    c1 <- c(c1, ""); c2 <- c(c2, "")

    # format
    c1 <- format(c1, width = 40L)
    c2 <- format(c2, width = 11L + max(0, (nd - 3L)) * 4L, justify = "right")

    # character matrix
    M <- cbind(c1, c2, deparse.level = 0)
    colnames(M) <- rep("",  ncol(M))
    rownames(M) <- rep(" ", nrow(M))

    # print
    write.table(M, row.names = TRUE, col.names = FALSE, quote = FALSE)

    invisible(M)
}

# print rotation information
lav_object_print_rotation <- function(object, header = FALSE, nd = 3L) {

    # header
    if(header) {
        cat("Rotation information:\n\n")
    }

    # container
    c1 <- c2 <- character(0L)

    # rotation method
    c1 <- c(c1, "Rotation method")
    if(object@Options$rotation.args$orthogonal) {
        MM <- paste(toupper(object@Options$rotation), " ", "ORTHOGONAL",
                    sep = "")
    } else {
        MM <- paste(toupper(object@Options$rotation), " ", "OBLIQUE",
                    sep = "")
    }
    c2 <- c(c2, MM)

    # method options
    if(object@Options$rotation == "geomin") {
        c1 <- c(c1, "Geomin epsilon")
        c2 <- c(c2, object@Options$rotation.args$geomin.epsilon)
    } else if(object@Options$rotation == "orthomax") {
        c1 <- c(c1, "Orthomax gamma")
        c2 <- c(c2, object@Options$rotation.args$orthomax.gamma)
    } else if(object@Options$rotation == "cf") {
        c1 <- c(c1, "Crawford-Ferguson gamma")
        c2 <- c(c2, object@Options$rotation.args$cf.gamma)
    } else if(object@Options$rotation == "oblimin") {
        c1 <- c(c1, "Oblimin gamma")
        c2 <- c(c2, object@Options$rotation.args$oblimin.gamma)
    }

    # rotation algorithm
    c1 <- c(c1, "Rotation algorithm (rstarts)")
    tmp <- paste(toupper(object@Options$rotation.args$algorithm),
                 " (", object@Options$rotation.args$rstarts, ")", sep = "")
    c2 <- c(c2, tmp)

    # Standardized metric (or not)
    c1 <- c(c1, "Standardized metric")
    if(object@Options$rotation.args$std.ov) {
        c2 <- c(c2, "TRUE")
    } else {
        c2 <- c(c2, "FALSE")
    }

    # Row weights
    c1 <- c(c1, "Row weights")
    tmp.txt <- object@Options$rotation.args$row.weights
    c2 <- c(c2, paste(toupper(substring(tmp.txt, 1, 1)),
                      substring(tmp.txt, 2), sep = ""))

    # empty row
    c1 <- c(c1, ""); c2 <- c(c2, "")

    # format c1/c2
    c1 <- format(c1, width = 33L)
    c2 <- format(c2, width = 18L + max(0, (nd - 3L)) * 4L, justify = "right")

    # create character matrix
    M <- cbind(c1, c2, deparse.level = 0)
    colnames(M) <- rep("",  ncol(M))
    rownames(M) <- rep(" ", nrow(M))

    # print
    write.table(M, row.names = TRUE, col.names = FALSE, quote = FALSE)

    invisible(M)
}

# test statistics
#
# print 'blocks' of test statistics
# blocks with 'scaling.factors' come first
#
lav_object_print_test_statistics <- function(object, nd = 3L) {

    num.format  <- paste("%", max(8L, nd + 5L), ".", nd, "f", sep = "")

    # check if test == "none", @test is empty, or estimator = "MML"
    if(object@Options$test[1] == "none" || length(object@test) == 0L ||
       object@Options$estimator == "MML") {
        return( character(0L) )
    }

    # header
    cat("Model Test User Model:\n")

    # extract tests
    TEST <- object@test

    # locate 'robust' tests (here: having a scaling factor)
    has.no.scaling <- unname(sapply((lapply(TEST, "[[", "scaling.factor")),
                             is.null))
    robust.idx     <- which(!has.no.scaling)
    non.robust.idx <- which(has.no.scaling)
    if(length(robust.idx) > 0L) {
        non.robust.idx <- non.robust.idx[-1] # remove 'standard', because it
                                             # is shown together with robust
    }
    BLOCKS <- c(robust.idx, non.robust.idx)
    nBlocks <- length(BLOCKS)

    # print out blocks
    for(block in BLOCKS) {

        # one or two-columns for this block?
        if(length(robust.idx) > 0L && block %in% robust.idx) {
            twocolumn <- TRUE
        } else {
            twocolumn <- FALSE
        }

        if(!twocolumn) {
            if(is.na(TEST[[block]]$df) || TEST[[block]]$df == 0L) {
                c1 <- c("Test statistic", "Degrees of freedom")
                c2 <- c(sprintf(num.format, TEST[[block]]$stat),
                        ifelse(TEST[[block]]$df %% 1 == 0, # integer
                               TEST[[block]]$df,
                               sprintf(num.format, TEST[[block]]$df)))
                c3 <- c("", "")
            } else {
                PLABEL <- "P-value"
                if(!is.null(TEST[[block]]$refdistr)) {
                    if(TEST[[block]]$refdistr == "chisq") {
                        PLABEL <- "P-value (Chi-square)"
                    } else if(TEST[[block]]$refdistr == "unknown") {
                        PLABEL <- "P-value (Unknown)"
                    } else if(TEST[[block]]$refdistr == "bootstrap") {
                        PLABEL <- "P-value (Bollen-Stine bootstrap)"
                    }
                }
                c1 <- c("Test statistic", "Degrees of freedom", PLABEL)
                c2 <- c(sprintf(num.format, TEST[[block]]$stat),
                        ifelse(TEST[[block]]$df %% 1 == 0, # integer          
                               TEST[[block]]$df,                    
                               sprintf(num.format, TEST[[block]]$df)),
                        sprintf(num.format, TEST[[block]]$pvalue))
                c3 <- c("", "", "")
            }
        } else {
            if(is.na(TEST[[block]]$df) || TEST[[block]]$df == 0L) {
                c1 <- c("Test Statistic", "Degrees of freedom")
                c2 <- c(sprintf(num.format, TEST[[1L]]$stat),
                        ifelse(TEST[[1L]]$df %% 1 == 0, # integer          
                               TEST[[1L]]$df,                    
                               sprintf(num.format, TEST[[1L]]$df)))
                c3 <- c(sprintf(num.format, TEST[[block]]$stat),
                        ifelse(TEST[[block]]$df %% 1 == 0, # integer          
                               TEST[[block]]$df,                    
                               sprintf(num.format, TEST[[block]]$df)))
            } else {
                if(!is.null(TEST[[1L]]$refdistr)) {
                    if(TEST[[1L]]$refdistr == "chisq") {
                        PLABEL <- "P-value (Chi-square)"
                    } else if(TEST[[1L]]$refdistr == "unknown") {
                        PLABEL <- "P-value (Unknown)"
                    } else {
                        PLABEL <- "P-value"
                    }
                }
                c1 <- c("Test Statistic", "Degrees of freedom", PLABEL,
                        "Scaling correction factor")
                c2 <- c(sprintf(num.format, TEST[[1L]]$stat),
                        ifelse(TEST[[1L]]$df %% 1 == 0, # integer          
                               TEST[[1L]]$df,                    
                               sprintf(num.format, TEST[[1L]]$df)),
                        sprintf(num.format, TEST[[1L]]$pvalue), "")
                c3 <- c(sprintf(num.format, TEST[[block]]$stat),
                        ifelse(TEST[[block]]$df %% 1 == 0, # integer          
                               TEST[[block]]$df,                    
                               sprintf(num.format, TEST[[block]]$df)),
                        sprintf(num.format, TEST[[block]]$pvalue),
                        sprintf(num.format, TEST[[block]]$scaling.factor))
                if(!is.null(TEST[[block]]$shift.parameter)) {
                    if(object@Data@ngroups == 1L) {
                        c1 <- c(c1, "Shift parameter")
                        c2 <- c(c2, "")
                        c3 <- c(c3, 
                            sprintf(num.format, TEST[[block]]$shift.parameter))
                    } else {
                        c1 <- c(c1, "Shift parameter for each group:")
                        c2 <- c(c2, "")
                        c3 <- c(c3, "")
                        for(g in 1:object@Data@ngroups) {
                            c1 <- c(c1, sprintf("    %-38s", 
                                            object@Data@group.label[[g]]))
                            c2 <- c(c2, "")
                            c3 <- c(c3, sprintf(num.format, 
                                          TEST[[block]]$shift.parameter[g]))
                        } 
                    } 
                } # shift
            }
        }

        # if twocolumn, label first row
        if(twocolumn && block == BLOCKS[1]) {
            c1 <- c("", c1); c2 <- c("Standard", c2); c3 <- c("Robust", c3)
        } else {
            # empty row
            c1 <- c("", c1); c2 <- c("", c2); c3 <- c("", c3)
        }

        # format c1/c2
        c1 <- format(c1, width = 43L)
        c2 <- format(c2, width = 8L + max(0, (nd - 3L)) * 4L, justify = "right")
        c3 <- format(c3, width = 8L + nd, justify = "right")

        # create character matrix
        if(twocolumn) {
            M <- cbind(c1, c2, c3, deparse.level = 0)
        } else {
            M <- cbind(c1, c2, deparse.level = 0)
        }
        colnames(M) <- rep("",  ncol(M))
        rownames(M) <- rep(" ", nrow(M))

        # print
        write.table(M, row.names = TRUE, col.names = FALSE, quote = FALSE)

        # additional comment if twocolumn
        if(twocolumn && TEST[[1]]$df > 0L) {
            cat("    for the", TEST[[block]]$label, "\n")
        }

        # multiple groups?
        ngroups <- object@Data@ngroups
        if(ngroups > 1L) {
            c1 <- c2 <- c3 <- character(ngroups)
            for(g in 1:object@Data@ngroups) {
                tmp <- sprintf("  %-40s", object@Data@group.label[[g]])
                c1[g] <- format(tmp, width = 43L)
                if(!twocolumn) {
                    tmp <- sprintf(num.format, TEST[[block]]$stat.group[g])
                    c2[g] <- format(tmp, width = 8L + max(0, (nd - 3L)) * 4L, 
                                    justify = "right")
                } else {
                    tmp <- sprintf(num.format, TEST[[1]]$stat.group[g])
                    c2[g] <- format(tmp, width = 8L + max(0, (nd - 3L)) * 4L,
                                    justify = "right")
                    tmp <- sprintf(num.format, TEST[[block]]$stat.group[g])
                    c3[g] <- format(tmp, width = 8L + nd, justify = "right")
                }
            }
            if(twocolumn) {
                M <- cbind(c1, c2, c3, deparse.level = 0)
            } else {
                M <- cbind(c1, c2, deparse.level = 0)
            }
            colnames(M) <- rep("",  ncol(M))
            rownames(M) <- rep(" ", nrow(M))
            cat("  Test statistic for each group:\n")
            write.table(M, row.names = TRUE, col.names = FALSE, quote = FALSE)
        }
    
    } # blocks

    #invisible(M)
}

lav_object_print_short_summary <- function(object, nd = 3L) {

    # 1. print header
    lav_object_print_header(object)

    # 2. print optim info (including estimator)
    lav_object_print_optim(object, nd = nd)

    # 3. if EFA/ESEM, print rotation info
    if(.hasSlot(object@Model, "nefa") && object@Model@nefa > 0L) {
        lav_object_print_rotation(object, nd = nd)
    }

    # 4. print lavdata
    lav_data_print_short(object@Data, nd = nd)

    # 5. print test statistics
    lav_object_print_test_statistics(object, nd = nd)

    # 5b. only if MML was used?
    if(object@Options$estimator == "MML") {
        fm <- fitMeasures(object, c("logl", "npar", "aic", "bic", "bic2"),
                          output = "text")
        print.lavaan.fitMeasures(fm, nd = nd)
    }

}

