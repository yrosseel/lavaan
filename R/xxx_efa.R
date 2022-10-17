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
                bounds         = "pos.var",
                baseline       = FALSE,
                ...,
                output         = "efa") {

    # handle dotdotdot
    dotdotdot <- list(...)

    # twolevel?
    twolevel.flag <- !is.null(dotdotdot$cluster)

    # check for unallowed arguments
    if(!is.null(dotdotdot$group)) {
        stop("lavaan ERROR: efa has no support for multiple groups (for now)")
    }

    # handle ov.names
    if(!is.null(data) && inherits(data, "data.frame")) {
        if(length(ov.names) > 0L) {
            if(twolevel.flag) {
                data <- data[, c(ov.names, dotdotdot$cluster)]
            } else {
                data <- data[, ov.names, drop = FALSE]
            }
        } else {
            ov.names <- names(data)
        }
    }

    # check nfactors
    if(any(nfactors < 1L)) {
        stop("lavaan ERROR: nfactors must be greater than zero.")
    } else {
        # check for maximum number of factors
        # Fixme: can we do this more efficiently? also holds for categorical?
        nvar <- length(ov.names)
        p.star <- nvar * (nvar + 1)/2
        nfac.max <- 0L
        for(nfac in seq_len(nvar)) {
            # compute number of free parameters
            npar <- nfac*nvar + nfac*(nfac+1L)/2 + nvar - nfac^2
            if(npar > p.star) {
                nfac.max <- nfac - 1L
                break
            }
        }
        if(any(nfactors > nfac.max)) {
            stop("lavaan ERROR: when nvar = ", nvar,
                 " the maximum number of factors is ", nfac.max, sep = "")
        }
    }

    # output
    output <- tolower(output)
    if(!output %in% c("lavaan", "efa")) {
        stop("lavaan ERROR: output= must be either \"lavaan\" or \"efa\"")
    }
    if(output == "lavaan" && length(nfactors) > 1L) {
        stop("lavaan ERROR: when output = \"lavaan\", nfactors must be a single (integer) number.")
    }

    # fit models
    nfits <- length(nfactors)
    out <- vector("list", length = nfits)
    for(f in seq_len(nfits)) {
        # generate model syntax
        model.syntax <- lav_syntax_efa(ov.names = ov.names,
                                       nfactors = nfactors[f],
                                       twolevel = twolevel.flag)
        # call lavaan (using sem())
        FIT <- do.call("sem",
                       args = c(list(model         = model.syntax,
                                     data          = data,
                                     sample.cov    = sample.cov,
                                     sample.nobs   = sample.nobs,
                                     rotation      = rotation,
                                     rotation.args = rotation.args,
                                     bounds        = bounds,
                                     baseline      = baseline),
                                dotdotdot))

        if(output == "efa") {
            FIT@Options$model.type <- "efa"
        }

        out[[f]] <- FIT
    }

    # class
    if(nfits == 1L && output == "lavaan") {
        out <- out[[1]]
    } else {
        class(out) <- c("efaList", "list")
    }

    out
}

# generate 'efa' syntax for a single block of factors
lav_syntax_efa <- function(ov.names = NULL, nfactors = 1L, twolevel = FALSE) {

    if(twolevel) {
        tmp <- lav_syntax_efa(ov.names = ov.names, nfactors = nfactors)
        model <- c("level: 1", tmp, "level: 2", tmp)
    } else {
        model <- character(nfactors)
        for(f in seq_len(nfactors)) {
            txt <- paste('efa("efa")*f', f, " =~ ",
                         paste(ov.names, collapse = " + "), sep = "")
            model[f] <- txt
        }
    }

    model
}

lav_summary_efa <- function(object) {

    nblocks <- object@Model@nblocks


    STD <- lavTech(object, "std", add.class = TRUE, add.labels = TRUE,
                   list.by.group = FALSE)
    lambda.idx <- which(names(STD) == "lambda")
    theta.idx  <- which(names(STD) == "theta")
    psi.idx    <- which(names(STD) == "psi")
    LAMBDA <- STD[lambda.idx]
    THETA  <- STD[theta.idx]
    PSI    <- STD[psi.idx]

    # eigenvalues correlation matrix
    std.ov <- object@Options$rotation.args$std.ov
    COV <- object@h1$implied$cov # h1
    if(std.ov) {
        COV <- lapply(COV, cov2cor)
    }
    eigvals <- lapply(seq_len(nblocks), function(b) {
                        tmp <- eigen(COV[[b]], only.values = TRUE)$values
                        names(tmp) <- paste("ev", 1:nrow(LAMBDA[[b]]), sep = "")
                        class(tmp) <- c("lavaan.vector", "numeric")
                        tmp
                      })

    # sum-of-squares table
    sumsq.table <- lapply(seq_len(nblocks), function(b) {
                              sumsq <- diag(PSI[[b]] %*% crossprod(LAMBDA[[b]]))
                              propvar <- sumsq/nrow(LAMBDA[[b]])
                              cumvar <- cumsum(propvar)
                              tmp <- rbind(sumsq, propvar, cumvar)
                              rownames(tmp) <- c("Sum of squared loadings",
                                                 "Proportion Var",
                                                 "Cumulative Var")
                              class(tmp) <- c("lavaan.matrix", "matrix")
                              tmp
                          })

    # (factor) structure coefficients
    lambda.structure <- lapply(seq_len(nblocks), function(b) {
                                   tmp <- LAMBDA[[b]] %*% PSI[[b]]
                                   class(tmp) <- c("lavaan.matrix", "matrix")
                                   tmp
                               })

    # standard errors (if any)
    lambda.se <- psi.se <- theta.se <- NULL
    if(object@Options$se != "none") {
        SE <- lavTech(object, "std.se", add.class = TRUE, add.labels = TRUE,
                   list.by.group = FALSE)
        lambda.se <- SE[lambda.idx]
        theta.se  <- SE[theta.idx]
        psi.se    <- SE[psi.idx]
    }

    block.names <- character(0L)
    if(length(object@Data@group.label) > 0L) {
        block.names <- object@Data@group.label
    } else if(length(object@Data@level.label) > 0L) {
        block.names <- object@Data@level.label
    } else if(length(object@Data@group.label) > 0L &&
              length(object@Data@level.label) > 0L) {
        block.names <- paste(object@Data@group.label,
                             object@Data@level.label, sep = " -- ")
    }

    res <- list(nblocks     = nblocks,
                block.names = block.names,
                lambda      = LAMBDA,
                theta       = THETA,
                psi         = PSI,
                std.ov      = std.ov,
                eigvals     = eigvals,
                sumsq.table = sumsq.table,
                orthogonal  = object@Options$rotation.args$orthogonal,
                lambda.structure = lambda.structure,
                lambda.se   = lambda.se,
                psi.se      = psi.se,
                theta.se    = theta.se)

    res
}

# print only (standardized) loadings
print.lavaan.efa <- function(x, nd = 3L, cutoff = 0.3,
                          dot.cutoff = 0.1, ...) {
    # unclass
    y <- unclass(x)

    if(!y$header$optim.converged) {
        cat("** WARNING ** Optimizer did not end normally\n")
        cat("** WARNING ** Estimates below are most likely unreliable\n")
    }

    # loadings per block
    for(b in seq_len(y$efa$nblocks)) {
        cat("\n")
        if(length(y$efa$block.names) > 0L) {
            cat(y$efa$block.names[[b]], ":\n\n", sep = "")
        }
        LAMBDA <- unclass(y$efa$lambda[[b]])
        lav_print_loadings(LAMBDA, nd = nd, cutoff = cutoff,
                           dot.cutoff = dot.cutoff,
                           x.se = y$efa$lambda.se[[b]])
        cat("\n")
    }

    invisible(LAMBDA)
}

# print efaList
print.efaList <- function(x, nd = 3L, cutoff = 0.3,
                          dot.cutoff = 0.1, ...) {
    # unclass
    y <- unclass(x)

    nfits <- length(y)
    for(f in seq_len(nfits)) {
        res <-  lav_object_summary(y[[f]], fit.measures = FALSE,
                                           estimates    = FALSE,
                                           modindices   = FALSE,
                                           efa          = TRUE)
        print.lavaan.efa(res, nd = nd, cutoff = cutoff,
                         dot.cutoff = dot.cutoff, ...)
    }

    invisible(y)
}

# summary efaList
summary.efaList <- function(object, ...) {
    # unclass
    y <- unclass(object)

    # dotdotdot
    dotdotdot <- list(...)
    nd <- 3L; cutoff <- 0.3; dot.cutoff <- 0.1
    if(!is.null(dotdotdot$nd)) {
        nd <- dotdotdot$nd
    }
    if(!is.null(dotdotdot$cutoff)) {
        cutoff <- dotdotdot$cutoff
    }
    if(!is.null(dotdotdot$dot.cutoff)) {
        dot.cutoff <- dotdotdot$dot.cutoff
    }

    # extract useful info from first model
    out <- lav_object_summary(y[[1]], header = TRUE, estimates = FALSE,
                              efa = FALSE)
    lavaan.version <- out$header$lavaan.version
    estimator      <- out$optim$estimator
    cat("This is ",
        sprintf("lavaan %s", lavaan.version),
        " -- running exploratory factor analysis\n", sep = "")

    # everything converged?
    if(!all(sapply(y, lavInspect, "converged"))) {
        cat("lavaan WARNING: not all models did converge!\n")
    }
    cat("\n")

    # robust test statistics?
    robust.flag <- FALSE
    if(any(c("satorra.bentler", "yuan.bentler", "yuan.bentler.mplus",
             "mean.var.adjusted", "scaled.shifted") %in%
       unlist(sapply(slot(y[[1]], "test"), "[", "test")))) {
        robust.flag <- TRUE
    }

    # categorical?
    categorical.flag <- FALSE
    if(y[[1]]@Model@categorical) {
        categorical.flag <- TRUE
    }

    # fit.measures
    if(robust.flag) {
        if(categorical.flag) {
            fit.measures <- c("chisq.scaled", "df.scaled", "pvalue.scaled",
                              "rmsea.scaled")
        } else {
            fit.measures <- c("chisq.scaled", "df.scaled", "pvalue.scaled",
                              "rmsea.robust")
        }
    } else {
        fit.measures <- c("chisq", "df", "pvalue", "rmsea")
    }
    fit.names <- c("chisq", "df", "pvalue", "rmsea")

    # estimator
    estimator      <- out$optim$estimator
    estimator.args <- out$optim$estimator.args
    rotation       <- out$rotation
    rotation.args  <- out$rotation.args

        c1 <- c("Estimator")
        # second column
        tmp.est <- toupper(estimator)
        if(tmp.est == "DLS") {
            dls.first.letter <- substr(estimator.args$dls.GammaNT,
                                       1L, 1L)
            tmp.est <- paste("DLS-", toupper(dls.first.letter), sep = "")
        }
        c2 <- tmp.est

        # additional estimator args
        if(!is.null(estimator.args) &&
           length(estimator.args) > 0L) {
            if(estimator == "DLS") {
                c1 <- c(c1, "Estimator DLS value for a")
                c2 <- c(c2, estimator.args$dls.a)
            }
        }

# rotation method
        c1 <- c(c1, "Rotation method")
        if(rotation$rotation == "none") {
            MM <- toupper(rotation$rotation)
        } else if(rotation$rotation.args$orthogonal) {
            MM <- paste(toupper(rotation$rotation), " ", "ORTHOGONAL",
                    sep = "")
        } else {
            MM <- paste(toupper(rotation$rotation), " ", "OBLIQUE",
                        sep = "")
        }
        c2 <- c(c2, MM)

        if(rotation$rotation != "none") {

            # method options
            if(rotation$rotation == "geomin") {
                c1 <- c(c1, "Geomin epsilon")
                c2 <- c(c2, rotation$rotation.args$geomin.epsilon)
            } else if(rotation$rotation == "orthomax") {
                c1 <- c(c1, "Orthomax gamma")
                c2 <- c(c2, rotation$rotation.args$orthomax.gamma)
            } else if(rotation$rotation == "cf") {
                c1 <- c(c1, "Crawford-Ferguson gamma")
                c2 <- c(c2, rotation$rotation.args$cf.gamma)
            } else if(rotation$rotation == "oblimin") {
                c1 <- c(c1, "Oblimin gamma")
                c2 <- c(c2, rotation$rotation.args$oblimin.gamma)
            }

            # rotation algorithm
            c1 <- c(c1, "Rotation algorithm (rstarts)")
            tmp <- paste(toupper(rotation$rotation.args$algorithm),
                         " (", rotation$rotation.args$rstarts, ")", sep = "")
            c2 <- c(c2, tmp)

            # Standardized metric (or not)
           c1 <- c(c1, "Standardized metric")
            if(rotation$rotation.args$std.ov) {
                c2 <- c(c2, "TRUE")
            } else {
                c2 <- c(c2, "FALSE")
            }

            # Row weights
            c1 <- c(c1, "Row weights")
            tmp.txt <- rotation$rotation.args$row.weights
            c2 <- c(c2, paste(toupper(substring(tmp.txt, 1, 1)),
                              substring(tmp.txt, 2), sep = ""))
        }

        # format c1/c2
        c1 <- format(c1, width = 33L)
        c2 <- format(c2, width = 18L + max(0, (nd - 3L)) * 4L,
                     justify = "right")

        # create character matrix
        M <- cbind(c1, c2, deparse.level = 0)
        colnames(M) <- rep("",  ncol(M))
        rownames(M) <- rep(" ", nrow(M))

        # print
        write.table(M, row.names = TRUE, col.names = FALSE, quote = FALSE)

    # data
    if(!is.null(out$data)) {
        cat("\n")
        lav_data_print_short(out$data, nd = nd)
    }
    cat("\n")

    nfits <- length(y)
    RES <- lapply(y, lav_object_summary, header = FALSE,
               fit.measures = FALSE, estimates = TRUE, efa = TRUE)

    # number of factors
    nfactors <- sapply(RES, function(x) ncol(x$efa$lambda[[1]]))

    # fit measures
    if(nfits > 1L) {
        cat("Overview models:\n")
    } else {
        cat("Fit measures:\n")
    }
    Table <- do.call("rbind", (lapply(y, fitMeasures, fit.measures)))
    colnames(Table) <- fit.names
    rownames(Table) <- paste("nfactors = ", nfactors, sep = "")
    class(Table) <- c("lavaan.matrix", "matrix")
    print(Table, nd = nd, shift = 2L)

    # eigenvalues
    cat("\n")
    if(RES[[1]]$efa$std.ov) {
        cat("Eigenvalues correlation matrix:\n")
    } else {
        cat("Eigenvalues covariance matrix:\n")
    }
    for(b in seq_len(RES[[1]]$efa$nblocks)) {
        cat("\n")
        if(length(RES[[1]]$efa$block.names) > 0L) {
            cat(RES[[1]]$efa$block.names[[b]], ":\n\n", sep = "")
        }
        print(RES[[1]]$efa$eigvals[[b]], nd = nd, shift = 2L)
    } # blocks
    cat("\n")

    # print summary for each model
    for(f in seq_len(nfits)) {
        res <- RES[[f]]
        attr(res, "nd") <- nd
        attr(res, "cutoff") <- cutoff
        attr(res, "dot.cutoff") <- dot.cutoff

        cat("Number of factors: ", nfactors[f], "\n")
        print(res)
    }

    invisible(y)
}
