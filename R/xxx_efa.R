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
#
# workflow:
# - summary() first calls summary.efaList()
# - for each model, summary.efaList() calls lav_object_summary() with
#   efa = TRUE and efa.args
# - for each model, lav_object_summary() calls
#   lav_summary_efa(object, efa.args = efa.args) to populate the $efa slot
#
lav_summary_efa <- function(object,
                            efa.args = list(lambda           = TRUE,
                                            theta            = TRUE,
                                            psi              = TRUE,
                                            eigenvalues      = TRUE,
                                            sumsq.table      = TRUE,
                                            lambda.structure = FALSE,
                                            se               = FALSE,
                                            zstat            = FALSE,
                                            pvalue           = FALSE)) {

    stopifnot(inherits(object, "lavaan"))

    nblocks <- object@Model@nblocks
    orthogonal.flag <- object@Options$rotation.args$orthogonal

    # get standardized solution
    LAMBDA <- THETA <- PSI <- NULL
    STD <- lavTech(object, "std", add.class = TRUE, add.labels = TRUE,
                   list.by.group = FALSE)
    lambda.idx <- which(names(STD) == "lambda")
    theta.idx  <- which(names(STD) == "theta")
    psi.idx    <- which(names(STD) == "psi")

    # LAMBDA
    LAMBDA <- STD[lambda.idx]
    names(LAMBDA) <- NULL

    # THETA
    THETA  <- STD[theta.idx]
    # make THETA diagonal
    THETA <- lapply(seq_len(nblocks), function(b) {
                 tmp <- diag(THETA[[b]])
                 class(tmp) <- c("lavaan.vector", "numeric")
                 tmp
             })

    # PSI
    PSI    <- STD[psi.idx]
    names(PSI)    <- NULL

    # eigenvalues correlation matrix
    std.ov <- object@Options$rotation.args$std.ov
    COV <- object@h1$implied$cov # h1
    if(std.ov) {
        COV <- lapply(COV, cov2cor)
    }
    eigvals <- NULL
    if(efa.args$eigenvalues) {
        eigvals <- lapply(seq_len(nblocks), function(b) {
                        tmp <- eigen(COV[[b]], only.values = TRUE)$values
                        names(tmp) <- paste("ev", 1:nrow(LAMBDA[[b]]), sep = "")
                        class(tmp) <- c("lavaan.vector", "numeric")
                        tmp
                      })
    }

    # sum-of-squares table
    sumsq.table <- NULL
    if(efa.args$sumsq.table) {
        sumsq.table <- lapply(seq_len(nblocks), function(b) {
            nvar <- nrow(LAMBDA[[b]])
            nfactor <- ncol(LAMBDA[[b]])

            # sum of squares:
            # - if orthogonal, this is really the sum of the squared factor
            #   loadings
            # - if oblique, we need to take the correlation into account
            sumsq <- diag(PSI[[b]] %*% crossprod(LAMBDA[[b]]))

            # reorder
            if(nfactor > 1L) {
                # determine order
                order.idx <- sort.int(sumsq, decreasing = TRUE,                                                       index.return = TRUE)$ix
                # re-order from large to small
                sumsq <- sumsq[order.idx]
            }

            # Proportion 'explained' (= proportion of total sumsq)
            #   note: sum(sumsq) == sum(communalities)
            propexpl <- sumsq/sum(sumsq)

            # Proportion var (= sumsq/nvar)
            propvar <- sumsq/nrow(LAMBDA[[b]])

            # Cumulative var
            cumvar <- cumsum(propvar)

            # construct table
            tmp <- rbind(sumsq, propexpl, propvar, cumvar)

            # total + colnames
            if(nfactor > 1L) {
                # add total column
                tmp <- cbind(tmp, rowSums(tmp))
                tmp[4, ncol(tmp)] <- tmp[3, ncol(tmp)]
                 colnames(tmp) <-
                    c(colnames(LAMBDA[[b]])[order.idx],
                      "total")
            } else {
                colnames(tmp) <- colnames(LAMBDA[[b]])[1]
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
                               "Proportion of total",
                               "Proportion var",
                               "Cumulative var")

            # class
            class(tmp) <- c("lavaan.matrix", "matrix")

            tmp
        })
    } # sumsq.table

    # (factor) structure coefficients
    if(efa.args$lambda.structure) {
        lambda.structure <- lapply(seq_len(nblocks), function(b) {
                                   tmp <- LAMBDA[[b]] %*% PSI[[b]]
                                   class(tmp) <- c("lavaan.matrix", "matrix")
                                   tmp
                                   })
    } else {
        lambda.structure <- NULL
    }

    # standard errors (if any)
    lambda.se    <- theta.se    <- psi.se    <- NULL
    lambda.zstat <- theta.zstat <- psi.zstat <- NULL
    lambda.pval  <- theta.pval  <- psi.pval  <- NULL
    if(object@Options$se != "none") {
        SE <- lavTech(object, "std.se", add.class = TRUE, add.labels = TRUE,
                      list.by.group = FALSE)

        se.flag <- ( efa.args$se || efa.args$zstat || efa.args$pvalue )

        # ALWAYS use lambda.se
        if(efa.args$lambda) {
            lambda.se <- SE[lambda.idx]
            names(lambda.se) <- NULL
        }

        # theta.se
        if(se.flag && efa.args$theta) {
            theta.se  <- SE[theta.idx]
            # make theta.se diagonal
            theta.se <- lapply(seq_len(nblocks), function(b) {
                            tmp <- diag(theta.se[[b]])
                            class(tmp) <- c("lavaan.vector", "numeric")
                            tmp
                        })
        }

        # psi.se
        if(se.flag && efa.args$psi) {
            psi.se    <- SE[psi.idx]
            names(psi.se)   <- NULL
        }

        # compute zstat
        if(efa.args$zstat || efa.args$pvalue) {
            if(efa.args$lambda) {
                lambda.zstat <- lapply(seq_len(nblocks), function(b) {
                                tmp.se <- lambda.se[[b]]
                                tmp.se[ tmp.se < sqrt(.Machine$double.eps) ] <-
                                    as.numeric(NA)
                                tmp <- LAMBDA[[b]] / tmp.se
                                class(tmp) <- c("lavaan.matrix", "matrix")
                                tmp
                                })
            }
            if(efa.args$theta) {
                theta.zstat <- lapply(seq_len(nblocks), function(b) {
                                tmp.se <- theta.se[[b]]
                                tmp.se[ tmp.se < sqrt(.Machine$double.eps) ] <-
                                    as.numeric(NA)
                                tmp <- THETA[[b]] / tmp.se
                                class(tmp) <- c("lavaan.vector", "numeric")
                                tmp
                               })
            }
            if(efa.args$psi) {
                psi.zstat <- lapply(seq_len(nblocks), function(b) {
                                tmp.se <- psi.se[[b]]
                                tmp.se[ tmp.se < sqrt(.Machine$double.eps) ] <-
                                   as.numeric(NA)
                                tmp <- PSI[[b]] / tmp.se
                                class(tmp) <- c("lavaan.matrix.symmetric",
                                                "matrix")
                                tmp
                              })
            }
        }

        # compute pval
        if(efa.args$pvalue) {
            if(efa.args$lambda) {
                lambda.pval <- lapply(seq_len(nblocks), function(b) {
                                tmp <- 2 * (1 - pnorm( abs(lambda.zstat[[b]]) ))
                                class(tmp) <- c("lavaan.matrix", "matrix")
                                tmp
                               })
            }
            if(efa.args$theta) {
                theta.pval <- lapply(seq_len(nblocks), function(b) {
                                tmp <- 2 * (1 - pnorm( abs(theta.zstat[[b]]) ))
                                class(tmp) <- c("lavaan.vector", "numeric")
                                tmp
                               })
            }
            if(efa.args$psi) {
                psi.pval <- lapply(seq_len(nblocks), function(b) {
                                tmp <- 2 * (1 - pnorm( abs(psi.zstat[[b]]) ))
                                class(tmp) <- c("lavaan.matrix.symmetric",
                                                "matrix")
                                tmp
                            })
            }
        }
    } # se/zstat/pvalue

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

    # we remove them here; we may have needed them for other parts
    if(!efa.args$lambda) {
        LAMBDA <- NULL
    }
    if(!efa.args$theta) {
        THETA <- NULL
    }
    if(!efa.args$psi) {
        PSI <- NULL
    }
    if(!efa.args$se) {
        # always keep lambda.se (for the signif stars)
        theta.se <- psi.se <- NULL
    }
    if(!efa.args$zstat) {
        lambda.zstat <- theta.zstat <- psi.zstat <- NULL
    }

    res <- list(nblocks          = nblocks,
                block.names      = block.names,
                std.ov           = std.ov,
                eigvals          = eigvals,
                sumsq.table      = sumsq.table,
                orthogonal       = object@Options$rotation.args$orthogonal,
                lambda.structure = lambda.structure,
                lambda           = LAMBDA,
                theta            = THETA,
                psi              = PSI,
                lambda.se        = lambda.se,
                lambda.zstat     = lambda.zstat,
                lambda.pvalue    = lambda.pval,
                psi.se           = psi.se,
                psi.zstat        = psi.zstat,
                psi.pvalue       = psi.pval,
                theta.se         = theta.se,
                theta.zstat      = theta.zstat,
                theta.pvalue     = theta.pval)

    res
}

# print only (standardized) loadings
print.lavaan.efa <- function(x, nd = 3L, cutoff = 0.3,
                          dot.cutoff = 0.1, alpha.level = 0.01, ...) {
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
                           alpha.level = alpha.level,
                           x.se = y$efa$lambda.se[[b]])
        cat("\n")
    }

    invisible(LAMBDA)
}

# print efaList
print.efaList <- function(x, nd = 3L, cutoff = 0.3,
                          dot.cutoff = 0.1, alpha.level = 0.01, ...) {
    # unclass
    y <- unclass(x)

    nfits <- length(y)
    RES <- vector("list", nfits)
    for(ff in seq_len(nfits)) {
        res <-  lav_object_summary(y[[ff]], fit.measures = FALSE,
                                            estimates    = FALSE,
                                            modindices   = FALSE,
                                            efa          = TRUE,
                                            efa.args = list(
                                                lambda           = TRUE,
                                                theta            = FALSE,
                                                psi              = FALSE,
                                                eigenvalues      = FALSE,
                                                sumsq.table      = FALSE,
                                                lambda.structure = FALSE,
                                                se               = FALSE,
                                                zstat            = FALSE,
                                                pvalue           = FALSE))
        RES[[ff]] <- print.lavaan.efa(res, nd = nd, cutoff = cutoff,
                                      dot.cutoff = dot.cutoff,
                                      alpha.level = alpha.level, ...)
    }

    invisible(RES)
}

# summary efaList
summary.efaList <- function(object, nd = 3L, cutoff = 0.3, dot.cutoff = 0.1,
                            alpha.level = 0.01,
                            lambda = TRUE, theta = TRUE, psi = TRUE,
                            fit.table = TRUE,
                            eigenvalues = TRUE, sumsq.table = TRUE,
                            lambda.structure = FALSE, se = FALSE,
                            zstat = FALSE, pvalue = FALSE, ...) {
    # unclass
    y <- unclass(object)

    # construct efa.args
    efa.args <- list(lambda = lambda, theta = theta, psi = psi,
                     eigenvalues = eigenvalues, sumsq.table = sumsq.table,
                     lambda.structure = lambda.structure,
                     se = se, zstat = zstat, pvalue = pvalue)

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
               fit.measures = FALSE, estimates = TRUE, efa = TRUE,
               efa.args = efa.args)

    # number of factors (for ALL blocks)
    nfactors <- sapply(y, function(x) x@pta$nfac[[1]])

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
    Table <- NULL
    if(fit.table) {
        # categorical?
        categorical.flag <- FALSE
        if(y[[1]]@Model@categorical) {
            categorical.flag <- TRUE
        }
        if(robust.flag) {
            if(categorical.flag) {
                fit.measures <- c("chisq.scaled", "df.scaled", "pvalue.scaled")
            } else {
                fit.measures <- c("chisq.scaled", "df.scaled", "pvalue.scaled")
            }
        } else {
            fit.measures <- c("chisq", "df", "pvalue")
        }
        fit.names <- c("chisq", "df", "pvalue")

        # if ML, add AIC and BIC
        if(y[[1]]@Model@estimator == "ML") {
            fit.measures <- c("aic", "bic", fit.measures)
            fit.names <- c("aic", "bic", fit.names)
        }

        # table with fitmeasures
        Table <- do.call("rbind", (lapply(y, fitMeasures, fit.measures)))
        colnames(Table) <- fit.names
        rownames(Table) <- paste("nfactors = ", nfactors, sep = "")
        class(Table) <- c("lavaan.matrix", "matrix")
    }

    # create return object
    out <- list(lavaan.version = lavaan.version,
                converged.flag = converged.flag,
                estimator      = estimator,
                estimator.args = estimator.args,
                rotation       = rotation,
                rotation.args  = rotation.args,
                lavdata        = lavdata,
                fit.table      = Table,
                nfactors       = nfactors,
                model.list     = RES)

    # add nd, cutoff, dot.cutoff, ... as attributes (for printing)
    attr(out, "nd")          <- nd
    attr(out, "cutoff")      <- cutoff
    attr(out, "dot.cutoff")  <- dot.cutoff
    attr(out, "alpha.level") <- alpha.level

    # create class
    class(out) <- c("efaList.summary", "list")

    out
}


# print summary efaList
print.efaList.summary <- function(x, nd = 3L, cutoff = 0.3,
                                  dot.cutoff = 0.1, alpha.level = 0.01,
                                  ...) {

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
    }
    # get dot.cutoff, if it is stored as an attribute
    DC <- attr(y, "dot.cutoff")
    if(!is.null(DC) && is.numeric(DC)) {
        dot.cutoff <- DC
    }
    # get alpha.level, if it is stored as an attribute
    AL <- attr(y, "alpha.level")
    if(!is.null(AL) && is.numeric(AL)) {
        alpha.level <- AL
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

    # number of models
    nfits <- length(x$model.list)

    # number of factors
    nfactors <- x$nfactors

    # fit measures
    if(!is.null(x$fit.table)) {
        cat("\n")
        if(nfits > 1L) {
            cat("Overview models:\n")
        } else {
            cat("Fit measures:\n")
        }
        print(x$fit.table, nd = nd, shift = 2L)
    }

    # eigenvalues
    if(!is.null(x$model.list[[1]]$efa$eigvals[[1]])) {
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
    }

    # print summary for each model
    for(f in seq_len(nfits)) {
        res <- x$model.list[[f]]
        attr(res, "nd") <- nd
        attr(res, "cutoff") <- cutoff
        attr(res, "dot.cutoff") <- dot.cutoff
        attr(res, "alpha.level") <- alpha.level

        if(nfits > 1L) {
            if(f == 1L) {
                cat("\n")
            }
            cat("Number of factors: ", nfactors[f], "\n")
        }
        print(res)
    }

    invisible(y)
}
