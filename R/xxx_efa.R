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

    # output
    output <- tolower(output)
    if(!output %in% c("lavaan", "efa")) {
        stop("lavaan ERROR: output= must be either \"lavaan\" or \"efa\"")
    }

    # set order.lv.by= to "sumofsquares" by default
    # except for bi-factor, where we set order.lv.by= to "none"
    if(tolower(rotation) %in% c("bi-geomin", "bigeomin",
                                "bi-quartimin", "biquartimin")) {
        if(is.null(rotation.args$order.lv.by)) {
            rotation.args$order.lv.by <- "none"
        }
    }

    # generate model syntax
    model.syntax <- lav_syntax_efa(ov.names = ov.names, nfactors = nfactors,
                                   twolevel = twolevel.flag)
    # call lavaan (using sem())
    FIT <- do.call("sem",
                   args = c(list(model         = model.syntax,
                                 data          = data,
                                 sample.cov    = sample.cov,
                                 sample.nobs   = sample.nobs,
                                 rotation      = rotation,
                                 rotation.args = rotation.args),
                            dotdotdot))

    if(output == "efa") {
        FIT@Options$model.type <- "efa"
    }

    FIT
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
    COV <- object@SampleStats@cov
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
                              sumsq <- colSums(LAMBDA[[b]] * LAMBDA[[b]])
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
print.lavaan.efa <- function(y, nd = 3L, cutoff = 0.3,
                          dot.cutoff = 0.1, ...) {
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

