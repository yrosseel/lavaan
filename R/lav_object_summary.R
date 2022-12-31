# initial version: YR 03/05/2017

# major change: YR 14/06/2022 for 0.6-12
# - summary() is now silent if not printed
# - here, we only collect the necessary ingredients, and store them in a
#   a list
# - the result is a S3 class lavaan.summary
# - the actual printing is done by print.lavaan.summary (see lav_print.R)

# YR 26 July 2022: add fm.args= argument to change the way (some) fit measures
#                  are computed
# YR 24 Sept 2022: add efa= argument

# create summary of a lavaan object
lav_object_summary <- function(object, header       = TRUE,
                                       fit.measures = FALSE,
                                       fm.args      =
                                           list(standard.test     = "default",
                                                scaled.test       = "default",
                                                rmsea.ci.level       = 0.90,
                                                rmsea.h0.closefit    = 0.05,
                                                rmsea.h0.notclosefit = 0.08),
                                       estimates    = TRUE,
                                       ci           = FALSE,
                                       fmi          = FALSE,
                                       std          = FALSE,
                                       standardized = FALSE,
                                       remove.step1 = TRUE,
                                       cov.std      = TRUE,
                                       rsquare      = FALSE,
                                       std.nox      = FALSE,
                                       efa          = FALSE,
                                       efa.args     =
                                           list(lambda           = TRUE,
                                                theta            = TRUE,
                                                psi              = TRUE,
                                                eigenvalues      = TRUE,
                                                sumsq.table      = TRUE,
                                                lambda.structure = FALSE,
                                                se               = FALSE,
                                                zstat            = FALSE,
                                                pvalue           = FALSE),
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

        # 1. collect header information
        if(.hasSlot(object, "version")) {
            VERSION <- object@version
        } else {
            VERSION <- "pre 0.6"
        }
        res$header <- list(lavaan.version   = VERSION,
                           sam.approach     = (.hasSlot(object, "internal") &&
                                       !is.null(object@internal$sam.method)),
                           optim.method     = object@Options$optim.method,
                           optim.iterations = object@optim$iterations,
                           optim.converged  = object@optim$converged)

        # sam or sem?
        if(.hasSlot(object, "internal") &&
        !is.null(object@internal$sam.method)) {

            # SAM version

            # 2. sam header
            res$sam.header <-
                list(sam.method          = object@internal$sam.method,
                     sam.local.options   = object@internal$sam.local.options,
                     sam.mm.list         = object@internal$sam.mm.list,
                     sam.mm.estimator    = object@internal$sam.mm.estimator,
                     sam.struc.estimator = object@internal$sam.struc.estimator)

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
            res$optim <- list(estimator      = object@Options$estimator,
                              estimator.args = object@Options$estimator.args,
                              optim.method   = object@Options$optim.method,
                              npar           = object@Model@nx.free,
                              eq.constraints = object@Model@eq.constraints,
                              nrow.ceq.jac   = nrow(object@Model@ceq.JAC),
                              nrow.cin.jac   = nrow(object@Model@cin.JAC),
                              nrow.con.jac   = nrow(object@Model@con.jac),
                              con.jac.rank   = qr(object@Model@con.jac)$rank)


            # 3. if EFA/ESEM, summarize rotation info
            if(.hasSlot(object@Model, "nefa") && object@Model@nefa > 0L) {
                res$rotation <-
                    list(rotation      = object@Options$rotation,
                         rotation.args = object@Options$rotation.args)
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

    # efa-related info
    if(efa) {
        res$efa <- lav_summary_efa(object, efa.args = efa.args)
    } # efa

    # only if requested, add the additional fit measures
    if(fit.measures) {
        # some early warnings (to avoid a hard stop)
        if(object@Data@data.type == "none") {
            warning("lavaan WARNING: fit measures not available if there is no data\n\n")
        } else if(length(object@Options$test) == 1L &&
                         object@Options$test == "none") {
            warning("lavaan WARNING: fit measures not available if test = \"none\"\n\n")
        } else if(object@optim$npar > 0L && !object@optim$converged) {
            warning("lavaan WARNING: fit measures not available if model did not converge\n\n")
        } else {
            FIT <- lav_fit_measures(object, fit.measures = "default",
                                    fm.args = fm.args)
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
