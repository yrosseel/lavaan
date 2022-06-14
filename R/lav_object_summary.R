# initial version: YR 03/05/2017

# major change: YR 14/06/2022
# - summary() is now silent if not printed
# - here, we only collect the necessary ingredients, and store them in a
#   a list
# - the result is a S3 class lavaan.summary
# - the actual printing is done by print.lavaan.summary (see lav_print.R)

# create summary of a lavaan object
lav_object_summary <- function(object, header       = TRUE,
                                       fit.measures = FALSE,
                                       estimates    = TRUE,
                                       ci           = FALSE,
                                       fmi          = FALSE,
                                       std          = FALSE,
                                       standardized = FALSE,
                                       remove.step1 = TRUE,
                                       cov.std      = TRUE,
                                       rsquare      = FALSE,
                                       std.nox      = FALSE,
                                       modindices   = FALSE) {

    # return a list with the main ingredients
    res <- list()

    # this is to avoid partial matching of 'std' with std.nox
    standardized <- std || standardized

    if(std.nox) {
        standardized <- TRUE
    }

    # create the 'short' summary
    if(header) {

        # 1. create header/greeting including version information
        res$header <- lav_object_summary_header(object)

        # sam or sem?
        if(.hasSlot(object, "internal") &&
        !is.null(object@internal$sam.method)) {

            # SAM version

            # 2. sam header
            res$sam.header <- lav_object_summary_sam_header(object)

            # 3. no EFA (for now)?

            # 4. summarize lavdata
            res$data <- lav_data_summary_short(object@Data)

            # 5a. sam local test statistics
            res$sam <- list(sam.method    = object@internal$sam.method,
                            sam.mm.table  = object@internal$sam.mm.table,
                            sam.mm.rel    = object@internal$sam.mm.rel,
                            sam.struc.fit = object@internal$sam.struc.fit,
                            ngroups       = object@Data@ngroups,
                            group.label   = object@Data@group.label)

            # 5b. global test statistics (for global only)
            if(object@internal$sam.method == "global") {
                res$test <- object@test
            }

        } else {

            # SEM version

            # 2. summarize optim info (including estimator)
            res$optim <- lav_object_summary_optim(object)

            # 3. if EFA/ESEM, summarize rotation info
            if(.hasSlot(object@Model, "nefa") && object@Model@nefa > 0L) {
                res$efa <- lav_object_summary_rotation(object)
            }

            # 4. summarize lavdata
            res$data <- lav_data_summary_short(object@Data)

            # 5. test statistics
            TEST <- object@test
            # double check if we have attr(TEST, "info") (perhaps old object?)
            if(is.null(attr(TEST, "info"))) {
                lavdata <- object@Data
                lavoptions <- object@Options
                attr(TEST, "info") <-
                    list(ngroups = lavdata@ngroups,
                         group.label = lavdata@group.label,
                         information = lavoptions$information,
                         h1.information = lavoptions$h1.information,
                         observed.information = lavoptions$observed.information)
            }
            res$test <- TEST



        } # regular sem
    } # header

    # only if requested, add the additional fit measures
    if(fit.measures) {
        if(length(object@Options$test) == 1L && object@Options$test == "none") {
            warning("lavaan WARNING: fit measures not available if test = \"none\"\n\n")
        } else if(object@optim$npar > 0L && !object@optim$converged) {
            warning("lavaan WARNING: fit measures not available if model did not converge\n\n")
        } else {
            FIT <- fitMeasures(object, fit.measures="default")
            res$fit = FIT
        }
    }

    # main ingredient: the parameter table
    if(estimates) {
        PE <- parameterEstimates(object, ci = ci, standardized = standardized,
                                 rsquare = rsquare, fmi = fmi,
                                 cov.std = cov.std,
                                 remove.eq = FALSE, remove.system.eq = TRUE,
                                 remove.ineq = FALSE, remove.def = FALSE,
                                 remove.nonfree = FALSE,
                                 remove.step1 = remove.step1,
                                 #remove.nonfree.scales = TRUE,
                                 output = "text",
                                 header = TRUE)
        if(standardized && std.nox) {
            PE$std.all <- NULL
        }
        res$pe <- as.data.frame(PE)
    }

    # modification indices?
    if(modindices) {
        MI <- modificationIndices(object, standardized=TRUE, cov.std = cov.std)
        res$mi <- MI
    }

    # create lavaan.summary S3 class
    class(res) <- c("lavaan.summary", "list")

    res
}


# header
lav_object_summary_header <- function(object) {

    # sam or sem?
    if(.hasSlot(object, "internal") && !is.null(object@internal$sam.method)) {
        txt1 <- paste("This is ", sprintf("lavaan %s",
                      packageDescription("lavaan", fields="Version")),
                      " -- using the SAM approach to SEM", sep = "")
        txt2 <- ""
    } else {
        txt1 <- paste(sprintf("lavaan %s ",
                      packageDescription("lavaan", fields = "Version")),
                      sep = "")

        # catch FAKE run
        FAKE <- FALSE
        if(object@Options$optim.method == "none") {
            FAKE <- TRUE
        }

        # Convergence or not?
        txt2 <- ""
        if(FAKE) {
            txt1 <- paste(txt1, "-- DRY RUN with 0 iterations --", sep = "")
        } else if(object@optim$iterations > 0) {
            if(object@optim$converged) {
            txt1 <- paste(txt1, sprintf("ended normally after %i iterations",
                                object@optim$iterations), sep = "")
            } else {
                txt1 <-
                    paste(txt1,
                          sprintf("did NOT end normally after %i iterations",
                          object@optim$iterations), sep = "")
                txt2 <-
                    "** WARNING ** Estimates below are most likely unreliable"
            }
        } else {
            txt1 <- paste(txt1, "did not run (perhaps do.fit = FALSE)?",
                          sep = "")
            txt2 <-
                "** WARNING ** Estimates below are simply the starting values"
        }
    }

    txt <- c(txt1, txt2)
    txt
}

# optim
lav_object_summary_optim <- function(object) {

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

    list(c1 = c1, c2 = c2)
}

# print rotation information
lav_object_summary_rotation <- function(object) {

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

    list(c1 = c1, c2 = c2)
}

# sam header info
lav_object_summary_sam_header <- function(object) {

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

    # estimator measurement blocks
    c1 <- c(c1, "Estimator measurement part")
    c2 <- c(c2, SAM$sam.mm.estimator)

    # estimator structural part
    c1 <- c(c1, "Estimator  structural part")
    c2 <- c(c2, SAM$sam.struc.estimator)

    list(c1 = c1, c2 = c2)
}

