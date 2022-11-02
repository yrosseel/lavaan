# This function will (eventually) contain functions that can be used to
# 'predict' the values of outcome variables (y), given the values of
# input variables (x).

# first version YR 2 Nov 2022

# method = "predictive.distribution" is based on the following article:

# Mark de Rooij, Julian D. Karch, Marjolein Fokkema, Zsuzsa Bakk, Bunga Citra
# Pratiwi & Henk Kelderman (2022) SEM-Based Out-of-Sample Predictions,
# StructuralEquation Modeling: A Multidisciplinary Journal
# DOI:10.1080/10705511.2022.2061494

# main function
lavPredictY <- function(object,
                        newdata  = NULL,
                        ynames   = lavNames(object, "ov.y"),
                        xnames   = lavNames(object, "ov.x"),
                        method   = "predictive.distribution",
                        label    = TRUE,
                        assemble = TRUE) {

    stopifnot(inherits(object, "lavaan"))
    lavmodel       <- object@Model
    lavdata        <- object@Data
    #lavsamplestats <- object@SampleStats
    lavimplied     <- object@implied

    # need full data set
    if(is.null(newdata)) {
        # use internal copy:
        if(lavdata@data.type != "full") {
            stop("lavaan ERROR: sample statistics were used for fitting and newdata is empty")
        } else if(is.null(lavdata@X[[1]])) {
            stop("lavaan ERROR: no local copy of data; FIXME!")
        } else {
            data.obs <- lavdata@X
            ov.names <- lavdata@ov.names
        }
        # eXo <- lavdata@eXo
    } else {
        OV <- lavdata@ov
        newData <- lavData(data        = newdata,
                           group       = lavdata@group,
                           ov.names    = lavdata@ov.names,
                           ov.names.x  = lavdata@ov.names.x,
                           ordered     = OV$name[ OV$type == "ordered" ],
                           lavoptions  = list(std.ov = lavdata@std.ov,
                                              group.label = lavdata@group.label,
                                              missing = lavdata@missing,
                                              warn = TRUE),
                           allow.single.case = TRUE)
        # if ordered, check if number of levels is till the same (new in 0.6-7)
        if(lavmodel@categorical) {
            orig.ordered.idx <- which(lavdata@ov$type == "ordered")
            orig.ordered.lev <- lavdata@ov$nlev[orig.ordered.idx]
            match.new.idx <- match(lavdata@ov$name[orig.ordered.idx],
                                    newData@ov$name)
            new.ordered.lev <- newData@ov$nlev[match.new.idx]
            if(any(orig.ordered.lev - new.ordered.lev != 0)) {
                stop("lavaan ERROR: ",
                     "mismatch number of categories for some ordered variables",
                     "\n\t\tin newdata compared to original data.")
            }
        }

        data.obs <- newData@X
        # eXo <- newData@eXo
        ov.names <- newData@ov.names
    }

    # check ynames
    if(length(ynames) == 0L) {
        stop("lavaan ERROR: please specifify the y-variables in the ynames= argument")
    } else if(!is.list(ynames)) {
        ynames <- rep(list(ynames), lavdata@ngroups)
    }

    # check xnames
    if(length(xnames) == 0L) {
        stop("lavaan ERROR: please specifify the x-variables in the xnames= argument")
    } else if(!is.list(xnames)) {
        xnames <- rep(list(xnames), lavdata@ngroups)
    }

    # create y.idx and x.idx
    y.idx <- x.idx <- vector("list", lavdata@ngroups)
    for(g in seq_len(lavdata@ngroups)) {
        # ynames in ov.names for this group?
        missing.idx <- which(!ynames[[g]] %in% ov.names[[g]])
        if(length(missing.idx) > 0L) {
            stop("lavaan ERROR: some variable names in ynames do not appear in the dataset:\n\t\t", paste(ynames[[g]][missing.idx], collapse = " "))
        } else {
            y.idx[[g]] <- which(ov.names[[g]] %in% ynames[[g]])
        }

        # xnames in ov.names for this group?
        missing.idx <- which(!xnames[[g]] %in% ov.names[[g]])
        if(length(missing.idx) > 0L) {
            stop("lavaan ERROR: some variable names in xnames do not appear in the dataset:\n", paste(xnames[[g]][missing.idx], collapse = " "))
        } else {
            x.idx[[g]] <- which(ov.names[[g]] %in% xnames[[g]])
        }
    }

    # prediction method
    if(method == "predictive.distribution") {
        out <- lav_predict_y_conditional(lavobject = NULL, lavmodel = lavmodel,
                   lavdata = lavdata, lavimplied = lavimplied,
                   data.obs = data.obs, y.idx = y.idx, x.idx = x.idx)
    } else {
        stop("lavaan ERROR: method must be \"predictive.distribution\" (for now).")
    }

    # label?
    if(label) {
        # column names
        for(g in seq_len(lavdata@ngroups)) {
            colnames(out[[g]]) <- ynames[[g]]
        }

        # group.labels
        if(lavdata@ngroups > 1L) {
            names(out) <- lavdata@group.label
        }
    }

    # lavaan.matrix
    out <- lapply(out, "class<-", c("lavaan.matrix", "matrix"))

    if(lavdata@ngroups == 1L) {
        res <- out[[1L]]
    } else {
        res <- out
    }

    # assemble multiple groups into a single data.frame?
    if(lavdata@ngroups > 1L && assemble) {
        if(!is.null(newdata)) {
            lavdata <- newData
        }
        DATA <- matrix(as.numeric(NA), nrow = sum(unlist(lavdata@norig)),
                                       ncol = ncol(out[[1L]])) # assume == per g
        colnames(DATA) <- colnames(out[[1L]])
        for(g in seq_len(lavdata@ngroups)) {
            DATA[ lavdata@case.idx[[g]], ] <- out[[g]]
        }
        DATA <- as.data.frame(DATA, stringsAsFactors = FALSE)

        if(!is.null(newdata)) {
            DATA[, lavdata@group] <- newdata[, lavdata@group ]
        } else {
            # add group
            DATA[, lavdata@group ] <- rep(as.character(NA), nrow(DATA))
            if(lavdata@missing == "listwise") {
                # we will loose the group label of omitted variables!
                DATA[unlist( lavdata@case.idx ), lavdata@group ] <-
                    rep( lavdata@group.label, unlist( lavdata@nobs ) )
            } else {
                DATA[unlist( lavdata@case.idx ), lavdata@group ] <-
                    rep( lavdata@group.label, unlist( lavdata@norig ) )
            }
        }

        res <- DATA
    }

    res
}


# method = "predictive.distribution"
lav_predict_y_conditional <- function(lavobject      = NULL,  # for convenience
                                      # object ingredients
                                      lavmodel       = NULL,
                                      lavdata        = NULL,
                                      lavimplied     = NULL,
                                      # new data
                                      data.obs       = NULL,
                                      # y and x
                                      y.idx          = NULL,
                                      x.idx          = NULL,
                                      # options
                                      level          = 1L) { # not used for now

    # full object?
    if(inherits(lavobject, "lavaan")) {
        lavmodel       <- lavobject@Model
        lavdata        <- lavobject@Data
        #lavsamplestats <- lavobject@SampleStats
        lavimplied     <- lavobject@implied
    } else {
        stopifnot(!is.null(lavmodel), !is.null(lavdata),
                  # !is.null(lavsamplestats),
                  !is.null(lavimplied))
    }

    # data.obs?
    if(is.null(data.obs)) {
        data.obs <- lavdata@X
    }

    # checks
    if(lavmodel@conditional.x) {
        stop("lavaan ERROR: no support for contidional.x = TRUE (yet).")
    }
    if(lavmodel@categorical) {
        stop("lavaan ERROR: no support for categorical data (yet).")
    }
    if(lavdata@nlevels > 1L) {
        stop("lavaan ERROR: no support for multilevel data (yet).")
    }
    #if(lavdata@missing %in% c("ml", "ml.x")) {
    #    stop("lavaan ERROR: no support for missing = \"ml\" (yet).")
    #}

    # output container
    YPRED <- vector("list", length = lavdata@ngroups)

    # run over all groups
    for(g in 1:lavdata@ngroups) {

        # multiple levels?
        if(lavdata@nlevels > 1L) {
            # TODO!
            stop("lavaan ERROR: no support for multilevel data (yet)!")
        } else {
            data.obs.g <- data.obs[[g]]

            # we assume conditional.x = FALSE (for now)
            cov.g   <- lavimplied$cov[[g]]

            if(lavmodel@meanstructure) {
                mean.g <- lavimplied$mean[[g]]
            } else {
                # use unrestricted mean vector
                mean.g <- colMeans(data.obs.g, na.rm = TRUE)
                # or zero vector?
                # mean.g <- rep(0, ncol(data.obs.g))
            }

            y.idx.g <- y.idx[[g]]
            x.idx.g <- x.idx[[g]]

            # partition y/x
            Sxx <- cov.g[x.idx.g, x.idx.g, drop = FALSE]
            Sxy <- cov.g[x.idx.g, y.idx.g, drop = FALSE]

            # x-data only
            Xtest <- data.obs.g[, x.idx.g, drop = FALSE]

            # mx/my
            mx <- mean.g[x.idx.g]
            my <- mean.g[y.idx.g]

            # center using mx
            Xtest <- t( t(Xtest) - mx )

            # prediction rule
            tmp <- Xtest %*% solve(Sxx, Sxy)
            YPRED[[g]] <- t( t(tmp) + my )

        } # single level

    } # g

    YPRED
}

