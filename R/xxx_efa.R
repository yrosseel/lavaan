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

# summary information for a single (lavaan) model
lav_summary_efa <- function(object) {

    stopifnot(inherits(object, "lavaan"))

    nblocks <- object@Model@nblocks
    orthogonal.flag <- object@Options$rotation.args$orthogonal

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
        nvar <- nrow(LAMBDA[[b]])
        nfactor <- ncol(LAMBDA[[b]])

        # sum of squares:
        # - if orthogonal, this is really the sum of the squared factor loadings
        # - if oblique, we need to take the correlation into account
        sumsq <- diag(PSI[[b]] %*% crossprod(LAMBDA[[b]]))

        # reorder
        if(nfactor > 1L) {
            # determine order
            order.idx <- sort.int(sumsq, decreasing = TRUE,                                                       index.return = TRUE)$ix
            # re-order from large to small
            sumsq <- sumsq[order.idx]
        }

        # Proportion var (= sumsq/nvar)
        propvar <- sumsq/nrow(LAMBDA[[b]])

        # Cumulative var
        cumvar <- cumsum(propvar)

        # Proportion 'explained' (= proportion of total sumsq)
        #   note: sum(sumsq) == sum(communalities)
        propexpl <- sumsq/sum(sumsq)

        # construct table
        tmp <- rbind(sumsq, propvar, cumvar, propexpl)

        # total + colnames
        if(nfactor > 1L) {
            # add total column
            tmp <- cbind(tmp, rowSums(tmp))
            tmp[3, ncol(tmp)] <- tmp[2, ncol(tmp)]
            colnames(tmp) <-
                c(colnames(LAMBDA[[b]])[order.idx],
                  "total")
        } else {
            colnames(tmp) <-
                colnames(LAMBDA[[b]])[1]
        }

        # rownames
        if(nfactor == 1L) {
            ssq.label <- "Sum of squared loadings"
        } else if(orthogonal.flag) {
          ssq.label <- "Sum of sq (ortho) loadings"
        } else {
          ssq.label <- "Sum of sq (obliq) loadings"
        }
        rownames(tmp) <- c(ssq.label,
                           "Proportion var",
                           "Cumulative var",
                           "Proportion of total")

        # class
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
    RES <- vector("list", nfits)
    for(ff in seq_len(nfits)) {
        res <-  lav_object_summary(y[[ff]], fit.measures = FALSE,
                                            estimates    = FALSE,
                                            modindices   = FALSE,
                                            efa          = TRUE)
        RES[[ff]] <- print.lavaan.efa(res, nd = nd, cutoff = cutoff,
                                      dot.cutoff = dot.cutoff, ...)
    }

    invisible(RES)
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

    # header information
    lavaan.version <- out$header$lavaan.version
    converged.flag <- all(sapply(y, lavInspect, "converged"))

    # estimator
    estimator      <- out$optim$estimator
    estimator.args <- out$optim$estimator.args

    # rotation
    rotation       <- out$rotation$rotation
    rotation.args  <- out$rotation$rotation.args

    # data
    lavdata <- out$data

    # main part: lav_object_summary information per model
    RES <- lapply(y, lav_object_summary, header = FALSE,
               fit.measures = FALSE, estimates = TRUE, efa = TRUE)

    # number of factors
    nfactors <- sapply(RES, function(x) ncol(x$efa$lambda[[1]]))

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

    # table with fitmeasures
    Table <- do.call("rbind", (lapply(y, fitMeasures, fit.measures)))
    colnames(Table) <- fit.names
    rownames(Table) <- paste("nfactors = ", nfactors, sep = "")
    class(Table) <- c("lavaan.matrix", "matrix")

    # create return object
    out <- list(lavaan.version = lavaan.version,
                converged.flag = converged.flag,
                estimator      = estimator,
                estimator.args = estimator.args,
                rotation       = rotation,
                rotation.args  = rotation.args,
                lavdata        = lavdata,
                fit.table      = Table,
                model.list     = RES)

    # add nd, cutoff and dot.cutoff as attributes
    attr(out, "nd")         <- nd
    attr(out, "cutoff")     <- cutoff
    attr(out, "dot.cutoff") <- dot.cutoff

    # create class
    class(out) <- c("efaList.summary", "list")

    out
}


# print summary efaList
print.efaList.summary <- function(x, ..., nd = 3L, cutoff = 0.3,
                                  dot.cutoff = 0.1) {

    # unclass
    y <- unclass(x)

    # get nd, if it is stored as an attribute
    ND <- attr(y, "nd")
    if(!is.null(ND) && is.numeric(ND)) {
        nd <- as.integer(ND)
    }
    # get cutoff, if it is stored as an attribute
    CT <- attr(y, "cutoff")
    if(!is.null(CT) && is.numeric(CT)) {
        cutoff <- CT
    } else {
        cutoff <- 0.3
    }
    # get dot.cutoff, if it is stored as an attribute
    DC <- attr(y, "dot.cutoff")
    if(!is.null(DC) && is.numeric(DC)) {
        dot.cutoff <- DC
    } else {
        dot.cutoff <- 0.1
    }

    cat("This is ",
        sprintf("lavaan %s", x$lavaan.version),
        " -- running exploratory factor analysis\n", sep = "")

    # everything converged?
    if(!x$converged.flag) {
        cat("lavaan WARNING: not all models did converge!\n")
    }
    cat("\n")


    # estimator
    c1 <- c("Estimator")
    # second column
    tmp.est <- toupper(x$estimator)
    if(tmp.est == "DLS") {
        dls.first.letter <- substr(x$estimator.args$dls.GammaNT,
                                   1L, 1L)
        tmp.est <- paste("DLS-", toupper(dls.first.letter), sep = "")
    }
    c2 <- tmp.est

    # additional estimator args
    if(!is.null(x$estimator.args) &&
       length(x$estimator.args) > 0L) {
        if(x$estimator == "DLS") {
            c1 <- c(c1, "Estimator DLS value for a")
            c2 <- c(c2, x$estimator.args$dls.a)
        }
    }

    # rotation method
    c1 <- c(c1, "Rotation method")
    if(x$rotation == "none") {
        MM <- toupper(x$rotation)
    } else if(x$rotation.args$orthogonal) {
        MM <- paste(toupper(x$rotation), " ", "ORTHOGONAL",
                sep = "")
    } else {
       MM <- paste(toupper(x$rotation), " ", "OBLIQUE",
                    sep = "")
    }
    c2 <- c(c2, MM)

    if(x$rotation != "none") {

        # method options
        if(x$rotation == "geomin") {
            c1 <- c(c1, "Geomin epsilon")
            c2 <- c(c2, x$rotation.args$geomin.epsilon)
        } else if(x$rotation == "orthomax") {
            c1 <- c(c1, "Orthomax gamma")
            c2 <- c(c2, x$rotation.args$orthomax.gamma)
        } else if(x$rotation == "cf") {
            c1 <- c(c1, "Crawford-Ferguson gamma")
            c2 <- c(c2, x$rotation.args$cf.gamma)
        } else if(x$rotation == "oblimin") {
            c1 <- c(c1, "Oblimin gamma")
            c2 <- c(c2, x$rotation.args$oblimin.gamma)
        }

        # rotation algorithm
        c1 <- c(c1, "Rotation algorithm (rstarts)")
        tmp <- paste(toupper(x$rotation.args$algorithm),
                     " (", x$rotation.args$rstarts, ")", sep = "")
        c2 <- c(c2, tmp)

        # Standardized metric (or not)
       c1 <- c(c1, "Standardized metric")
        if(x$rotation.args$std.ov) {
            c2 <- c(c2, "TRUE")
        } else {
            c2 <- c(c2, "FALSE")
        }

        # Row weights
        c1 <- c(c1, "Row weights")
        tmp.txt <- x$rotation.args$row.weights
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
    if(!is.null(x$lavdata)) {
        cat("\n")
        lav_data_print_short(x$lavdata, nd = nd)
    }
    cat("\n")

    # number of models
    nfits <- length(x$model.list)

    # number of factors
    nfactors <- sapply(x$model.list, function(xx) ncol(xx$efa$lambda[[1]]))

    # fit measures
    if(nfits > 1L) {
        cat("Overview models:\n")
    } else {
        cat("Fit measures:\n")
    }
    print(x$fit.table, nd = nd, shift = 2L)

    # eigenvalues
    cat("\n")
    if(x$model.list[[1]]$efa$std.ov) {
        cat("Eigenvalues correlation matrix:\n")
    } else {
        cat("Eigenvalues covariance matrix:\n")
    }
    for(b in seq_len(x$model.list[[1]]$efa$nblocks)) {
        cat("\n")
        if(length(x$model.list[[1]]$efa$block.names) > 0L) {
            cat(x$model.list[[1]]$efa$block.names[[b]], ":\n\n", sep = "")
        }
        print(x$model.list[[1]]$efa$eigvals[[b]], nd = nd, shift = 2L)
    } # blocks
    cat("\n")

    # print summary for each model
    for(f in seq_len(nfits)) {
        res <- x$model.list[[f]]
        attr(res, "nd") <- nd
        attr(res, "cutoff") <- cutoff
        attr(res, "dot.cutoff") <- dot.cutoff

        cat("Number of factors: ", nfactors[f], "\n")
        print(res)
    }

    invisible(y)
}
