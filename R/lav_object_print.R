# initial version: YR 03/05/2017

# header
lav_object_print_header <- function(object) {

    # sam or sem?
    if(.hasSlot(object, "internal") && !is.null(object@internal$sam.method)) {
        cat("This is ",
            sprintf("lavaan %s",
                        packageDescription("lavaan", fields="Version")),
            " -- using the SAM approach to SEM", sep = "")
    } else {
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
    }
    cat("\n")
}

# sam
lav_object_print_sam_header <- function(object,  nd = 3L) {

    # grab 'SAM' information from @internal slot
    SAM <- object@internal

    # sam method
    c1 <- c("SAM method")
    c2 <- toupper(SAM$sam.method)

    # options
    if(SAM$sam.method == "local") {
        c1 <- c(c1, "Mapping matrix M method")
        c2 <- c(c2, SAM$sam.local.options$M.method)
        # TODo: more!
    }

    # number of measurement blocks
    c1 <- c(c1, "Number of measurement blocks")
    c2 <- c(c2, length(SAM$sam.mm.list))

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

# optim
lav_object_print_optim <- function(object, nd = 3L) {

    #cat("Optimization information:\n\n")

    c1 <- c("Estimator")

    # second column
    tmp.est <- toupper(object@Options$estimator)
    if(tmp.est == "DLS") {
        dls.first.letter <- substr(object@Options$estimator.args$dls.GammaNT,
                                   1L, 1L)
        tmp.est <- paste("DLS-", toupper(dls.first.letter), sep = "")
    }
    c2 <- tmp.est

    # additional estimator args
    if(!is.null(object@Options$estimator.args) &&
       length(object@Options$estimator.args) > 0L) {
        if(object@Options$estimator == "DLS") {
            c1 <- c(c1, "Estimator DLS value for a")
            c2 <- c(c2, object@Options$estimator.args$dls.a)
        }
    }

    # optimization method + npar
    c1 <- c(c1, "Optimization method", "Number of model parameters")
    c2 <- c(c2, toupper(object@Options$optim.method),
                lavInspect(object, "npar"))

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
        con.jac.rank <- qr(object@Model@con.jac)$rank
        if(con.jac.rank == (nrow(object@Model@ceq.JAC) +
                            nrow(object@Model@cin.JAC)) ) {
            # nothing to do (don't print, as this is redundant information)
        } else {
            c1 <- c(c1, "Row rank of the constraints matrix")
            c2 <- c(c2, qr(object@Model@con.jac)$rank)
        }
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
    if(object@Options$rotation == "none") {
        MM <- toupper(object@Options$rotation)
    } else if(object@Options$rotation.args$orthogonal) {
        MM <- paste(toupper(object@Options$rotation), " ", "ORTHOGONAL",
                    sep = "")
    } else {
        MM <- paste(toupper(object@Options$rotation), " ", "OBLIQUE",
                    sep = "")
    }
    c2 <- c(c2, MM)


    if(object@Options$rotation != "none") {

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

    }

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
    lavoptions <- object@Options

    # check if test == "none", @test is empty, or estimator = "MML"
    if(lavoptions$test[1] == "none" || length(object@test) == 0L ||
       lavoptions$estimator == "MML") {
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
                if(TEST[[block]]$test == "scaled.shifted") {
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
                #c1 <- c(c1, paste("    for the", TEST[[block]]$label))
                c1 <- c(c1, paste("    ", TEST[[block]]$label))
                c2 <- c(c2, "")
                c3 <- c(c3, "")
            }
        }

        # if twocolumn, label first row
        if(twocolumn && block == BLOCKS[1]) {
            c1 <- c("", c1); c2 <- c("Standard", c2); c3 <- c("Robust", c3)
        } else {
            # empty row
            c1 <- c("", c1); c2 <- c("", c2); c3 <- c("", c3)
        }

        # if information type is different from 'se', print it
        if(length(lavoptions$information) > 1L &&
           lavoptions$information[1] != lavoptions$information[2]) {
            c1 <- c(c1, "Information")
            tmp.txt <- lavoptions$information[2]
            c2 <- c(c2, paste(toupper(substring(tmp.txt,1,1)),
                                      substring(tmp.txt, 2), sep = ""))
            c3 <- c(c3, "")
        }
        # if h1.information type is different from 'se', print it
        if(length(lavoptions$h1.information) > 1L &&
           lavoptions$h1.information[1] != lavoptions$h1.information[2]) {
            c1 <- c(c1, "Information saturated (h1) model")
            tmp.txt <- lavoptions$h1.information[2]
            c2 <- c(c2, paste(toupper(substring(tmp.txt,1,1)),
                                      substring(tmp.txt, 2), sep = ""))
            c3 <- c(c3, "")
        }
        # if observed.information type is different from 'se', print it
        if(length(lavoptions$observed.information) > 1L &&
           lavoptions$information[2] == "observed" &&
           (lavoptions$observed.information[1] !=
            lavoptions$observed.information[2]) ) {
            c1 <- c(c1, "Observed information based on")
            tmp.txt <- lavoptions$observed.information[2]
            c2 <- c(c2, paste(toupper(substring(tmp.txt,1,1)),
                                      substring(tmp.txt, 2), sep = ""))
            c3 <- c(c3, "")
        }


        # format c1/c2/c3 (note: fitMeasures uses 35/16/8)
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

lav_object_print_sam_test_statistics <- function(object, nd = 3L) {

    # format for numeric values
    num.format  <- paste("%", max(8L, nd + 5L), ".", nd, "f", sep = "")
    int.format  <- paste("%", max(8L, nd + 5L), "d", sep = "")
    char.format <- paste("%", max(8L, nd + 5L), "s", sep = "")

    # measurement
    tmp <- object@internal$sam.mm.table
    COLNAMES <- colnames(tmp); colnames(tmp) <- NULL
    M <- lapply(tmp, function(x) {
                if(is.integer(x)) {
                    sprintf(int.format, x)
                } else if(is.character(x)) {
                    sprintf(char.format, x)
                } else if(is.numeric(x)) {
                    sprintf(num.format, x)
                } else {
                    sprintf(char.format, as.character(x))
                }  })
    # make matrix
    M <- do.call("cbind", M)
    # add column names as first row
    M <- rbind(sprintf(char.format, COLNAMES), M, deparse.level = 0L)
    # add margin
    #M <- cbind(rep("aaa", nrow(M)), M, deparse.level = 0L)
    rownames(M) <- rep(" ", nrow(M))
    write.table(M, row.names = TRUE, col.names = FALSE, quote = FALSE)
}

lav_object_print_short_summary <- function(object, nd = 3L) {

    # 1. print header
    lav_object_print_header(object)

    # 1.b sam or sem?
    if(.hasSlot(object, "internal") && !is.null(object@internal$sam.method)) {
        cat("\n")
        lav_object_print_sam_header(object)

        # print lavdata
        lav_data_print_short(object@Data, nd = nd)

        # local test statistics
        lav_object_print_sam_test_statistics(object, nd = nd)

        # global test statistics (for global only)
        if(object@internal$sam.method == "global") {
            lav_object_print_test_statistics(object, nd = nd)
        }

    } else {
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
            fm <- fitMeasures(object, c("logl", "aic", "bic", "bic2"),
                              output = "text")
            print.lavaan.fitMeasures(fm, nd = nd, add.h0 = FALSE)
        }
    } # regular sem
}

