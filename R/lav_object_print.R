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

    # empty third column
    c3 <- rep("", length(c1))

    # empty last row
    c1 <- c(c1, ""); c2 <- c(c2, ""); c3 <- c(c3, "")

    # format
    c1 <- format(c1, width = 40L)
    c2 <- format(c2, width = 11L + max(0, (nd - 3L)) * 4L, justify = "right")
    c3 <- format(c3, width = 11L + nd, justify = "right")

    # character matrix
    M <- cbind(c1, c2, c3, deparse.level = 0)
    colnames(M) <- rep("",  ncol(M))
    rownames(M) <- rep(" ", nrow(M))

    # print
    write.table(M, row.names = TRUE, col.names = FALSE, quote = FALSE)

    invisible(M)
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

# short summary
lav_object_print_test_statistics <- function(object) {
    
    # robust/scaled statistics?
    if(object@Options$test %in% c("satorra.bentler",
                                  "yuan.bentler",
                                  "yuan.bentler.mplus",
                                  "mean.var.adjusted",
                                  "scaled.shifted") &&
       length(object@test) > 1L) {
        scaled <- TRUE
        if(object@Options$test == "scaled.shifted")
            shifted <- TRUE
        else
            shifted <- FALSE
    } else {
        scaled <- FALSE
        shifted <- FALSE
    }

    # 0. heading
    #h.txt <- sprintf("\nChi-square test user model (h0)",
    #                 object@Options$estimator)
    #t0.txt <- sprintf("  %-40s", "Estimator")
    #t1.txt <- sprintf("  %10s", object@Options$estimator)
    #t2.txt <- ifelse(scaled,
    #          sprintf("  %10s", "Robust"), "")
    #cat(t0.txt, t1.txt, t2.txt, "\n", sep="")

    # check if test == "none"
    if(object@Options$test != "none" && object@Options$estimator != "MML") {

        # 1. chi-square values
        t0.txt <- sprintf("  %-40s", "Model Fit Test Statistic")
        t1.txt <- sprintf("  %10.3f", object@test[[1]]$stat)
        t2.txt <- ifelse(scaled,
                  sprintf("  %10.3f", object@test[[2]]$stat), "")
        cat(t0.txt, t1.txt, t2.txt, "\n", sep="")

        # 2. degrees of freedom
        t0.txt <- sprintf("  %-40s", "Degrees of freedom")
        t1.txt <- sprintf("  %10i",   object@test[[1]]$df)
        t2.txt <- ifelse(scaled,
                         ifelse(round(object@test[[2]]$df) ==
                                object@test[[2]]$df,
                                sprintf("  %10i",   object@test[[2]]$df),
                                sprintf("  %10.3f", object@test[[2]]$df)),
                         "")
        cat(t0.txt, t1.txt, t2.txt, "\n", sep="")

        # 3. P-value
        if(is.na(object@test[[1]]$df)) {
            t0.txt <- sprintf("  %-40s", "P-value")
            t1.txt <- sprintf("  %10.3f", object@test[[1]]$pvalue)
            t2.txt <- ifelse(scaled,
                      sprintf("  %10.3f", object@test[[2]]$pvalue), "")
            cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
        } else if(object@test[[1]]$df > 0) {
            if(object@test[[1]]$refdistr == "chisq") {
                t0.txt <- sprintf("  %-40s", "P-value (Chi-square)")
            } else if(length(object@test) == 1L &&
                      object@test[[1]]$refdistr == "unknown") {
                t0.txt <- sprintf("  %-40s", "P-value (Unknown)")
            } else {
                t0.txt <- sprintf("  %-40s", "P-value")
            }
            t1.txt <- sprintf("  %10.3f", object@test[[1]]$pvalue)
            t2.txt <- ifelse(scaled,
                      sprintf("  %10.3f", object@test[[2]]$pvalue), "")
            cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
        } else {
            # FIXME: should we do this? To warn that exact 0.0 was not obtained?
            if(object@optim$fx > 0) {
                t0.txt <- sprintf("  %-35s", "Minimum Function Value")
                t1.txt <- sprintf("  %15.13f", object@optim$fx)
                t2.txt <- ""
                cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
            }
        }

        # 3b. Do we have a Bollen-Stine p-value?
        if(object@Options$test == "bollen.stine") {
            t0.txt <- sprintf("  %-40s", "P-value (Bollen-Stine Bootstrap)")
            t1.txt <- sprintf("  %10.3f", object@test[[2]]$pvalue)
            cat(t0.txt, t1.txt, "\n", sep="")
        }

        # 4. Scaling correction factor
        if(scaled) {
            t0.txt <- sprintf("  %-40s", "Scaling correction factor")
            t1.txt <- sprintf("  %10s", "")
            t2.txt <- sprintf("  %10.3f", object@test[[2]]$scaling.factor)
            cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
            if(object@Options$test == "yuan.bentler") {
                cat("    for the Yuan-Bentler correction\n")
            } else if(object@Options$test == "yuan.bentler.mplus") {
                cat("    for the Yuan-Bentler correction (Mplus variant)\n")
            } else if(object@Options$test == "satorra.bentler") {
                if(object@Options$mimic == "Mplus" &&
                   object@Options$estimator == "ML") {
                    cat("    for the Satorra-Bentler correction (Mplus variant)\n")
                } else if(object@Options$mimic == "Mplus" &&
                          object@Options$estimator == "DWLS") {
                    cat("    for the Satorra-Bentler correction (WLSM)\n")
                } else if(object@Options$mimic == "Mplus" &&
                          object@Options$estimator == "ULS") {
                    cat("    for the Satorra-Bentler correction (ULSM)\n")
                } else {
                    cat("    for the Satorra-Bentler correction\n")
                }
            } else if(object@Options$test == "mean.var.adjusted") {
                if(object@Options$mimic == "Mplus" &&
                   object@Options$estimator == "ML") {
                    cat("    for the mean and variance adjusted correction (MLMV)\n")
                } else if(object@Options$mimic == "Mplus" &&
                          object@Options$estimator == "DWLS") {
                    cat("    for the mean and variance adjusted correction (WLSMV)\n")
                } else if(object@Options$mimic == "Mplus" &&
                          object@Options$estimator == "ULS") {
                    cat("    for the mean and variance adjusted correction (ULSMV)\n")
                } else {
                    cat("    for the mean and variance adjusted correction\n")
                }
            }
        }

        # 4b. Shift parameter?
        if(shifted) {
            if(object@Data@ngroups == 1L) {
                t0.txt <- sprintf("  %-40s", "Shift parameter")
                t1.txt <- sprintf("  %10s", "")
                t2.txt <- sprintf("  %10.3f",
                                  object@test[[2]]$shift.parameter)
                cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
            } else { # multiple groups, multiple shift values!
                cat("  Shift parameter for each group:\n")
                for(g in 1:object@Data@ngroups) {
                    t0.txt <- sprintf("    %-38s", object@Data@group.label[[g]])
                    t1.txt <- sprintf("  %10s", "")
                    t2.txt <- sprintf("  %10.3f",
                                     object@test[[2]]$shift.parameter[g])
                    cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
                }
            }
            if(object@Options$mimic == "Mplus" &&
               object@Options$estimator == "DWLS") {
                cat("    for simple second-order correction (WLSMV)\n")
            } else {
                cat("    for simple second-order correction (Mplus variant)\n")
            }
        }
        if(object@Data@ngroups > 1L) {
            cat("\n")
            cat("Chi-square for each group:\n\n")
            for(g in 1:object@Data@ngroups) {
                t0.txt <- sprintf("  %-40s", object@Data@group.label[[g]])
                t1.txt <- sprintf("  %10.3f", object@test[[1]]$stat.group[g])
                t2.txt <- ifelse(scaled, sprintf("  %10.3f",
                                 object@test[[2]]$stat.group[g]), "")
                cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
            }
        }
    } # test != none

}

lav_object_print_short_summary <- function(object, nd = 3L) {

    # 1. print header
    lav_object_print_header(object)

    # 2. print optim info (including estimator)
    lav_object_print_optim(object, nd = nd)

    # 3. if EFA/ESEM, print rotation info
    if(.hasSlot(object@Model, "nefa") && object@Model@nefa > 0L) {
        lav_object_print_rotation(object)
    }

    # 4. print lavdata
    lav_data_print_short(object@Data)

    # 5. print test statistics
    lav_object_print_test_statistics(object)

    # 5b. only if MML was used?
    if(object@Options$estimator == "MML") {
        fm <- fitMeasures(object, c("logl", "npar", "aic", "bic", "bic2"))
        print.lavaan.fitMeasures(fm)
    }
}

