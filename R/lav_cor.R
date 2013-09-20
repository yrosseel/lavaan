# user-visible routine to 
# compute polychoric/polyserial/... correlations
#
# YR 17 Sept 2013

lavCor <- function(object, ordered = NULL, ov.names.x = NULL, 
                   # method="two.step", 
                   group = NULL, missing = "listwise",
                   WLS.W = FALSE, details = FALSE,
                   labels = TRUE, verbose = FALSE) {

    # check object class
    if(inherits(object, "lavaan")) {
        lav.data <- object@lavData
    } else if(inherits(object, "lavData")) {
        lav.data <- object
    } else if(inherits(object, "data.frame")) {
        NAMES <- names(object)
        if(!is.null(group)) {
            NAMES <- NAMES[- match(group, NAMES)]
        }
        lav.data <- lavData(data = object, group = group, 
                            ov.names = NAMES, ordered = ordered,
                            ov.names.x = ov.names.x,
                            missing = missing)
    } else {
        stop("lavaan ERROR: lavCor can not handle objects of class ", 
             class(object))
    }

    out <- lav_cor(lav.data = lav.data, WLS.W = WLS.W, details = details, 
                   labels = labels, verbose = verbose)

    out
}

# internal version
lav_cor <- function(lav.data = NULL, WLS.W = FALSE, details = FALSE, 
                    labels = FALSE, verbose = FALSE) {

    # shortcuts
    vartable   <- lav.data@ov
    ov.names   <- lav.data@ov.names
    ngroups    <- lav.data@ngroups

    COR <- vector("list", length=ngroups)
    for(g in 1:ngroups) {
        ov.types  <- vartable$type[ match(ov.names[[g]], vartable$name) ]
        ov.levels <- vartable$nlev[ match(ov.names[[g]], vartable$name) ]
        CAT <- muthen1984(Data       = lav.data@X[[g]],
                          ov.names   = ov.names[[g]],
                          ov.types   = ov.types,
                          ov.levels  = ov.levels,
                          ov.names.x = lav.data@ov.names.x[[g]],
                          eXo        = lav.data@eXo[[g]],
                          group      = g, # for error messages only
                          missing    = lav.data@missing, # listwise or pairwise?
                          WLS.W      = WLS.W,
                          verbose    = verbose)
        COR[[g]] <- unname(CAT$COV)
        if(details) {
            attr(COR[[g]], "ov.names") <- ov.names[[g]]
            attr(COR[[g]], "ov.types") <- ov.types
            attr(COR[[g]], "TH")       <- CAT$TH
            attr(COR[[g]], "TH.IDX")   <- CAT$TH.IDX
            attr(COR[[g]], "TH.NAMES") <- CAT$TH.NAMES
            attr(COR[[g]], "SLOPES")   <- CAT$SLOPES
            attr(COR[[g]], "VAR")      <- CAT$VAR
            attr(COR[[g]], "COV")      <- CAT$COV
        }
        if(WLS.W) {
            th <- unlist(CAT$TH)
            th[ov.types == "numeric"] <- -1*th[ov.types == "numeric"]
            WLS.obs <- c(th,
                         vec(CAT$SLOPES), # FIXME
                         unlist(CAT$VAR[ov.types == "numeric"]),
                         vech(CAT$COV, diagonal=FALSE))
            attr(COR[[g]], "WLS.obs") <- WLS.obs
            attr(COR[[g]], "WLS.W") <- CAT$WLS.W
        }

        if(labels) {
            rownames(COR[[g]]) <- colnames(COR[[g]]) <- ov.names[[g]]
            class(COR[[g]]) <- c("lavaan.matrix.symmetric", "matrix")
        }
    }

    if(ngroups == 1L) {
        out <- COR[[1]]
    } else {
        out <- COR
    }

    out
}

# test for bivariate normality
lavCorCellFit <- function(object, ordered = NULL, group = NULL, 
                          ov.names.x = NULL, missing = "listwise") {

    # check object class
    if(inherits(object, "lavaan")) {
        lav.data <- object@lavData
    } else if(inherits(object, "lavData")) {
        lav.data <- object
    } else if(inherits(object, "data.frame")) {
        NAMES <- names(object)
        if(!is.null(group)) {
            NAMES <- NAMES[- match(group, NAMES)]
        }
        lav.data <- lavData(data = object, group = group,
                            ov.names = NAMES, ordered = ordered,
                            ov.names.x = ov.names.x,
                            missing = missing)
    } else {
        stop("lavaan ERROR: lavCor can not handle objects of class ",
             class(object))
    }

    PI <- unlist(lav_pairwise_tables_sample_pi(object = lav.data))

    out <- lav_pairwise_tables(object    = object,
                               vartable  = vartable,
                               X         = lav.data@X, 
                               ov.names  = lav.data@ov.names, 
                               Pi        = PI,
                               average   = TRUE,
                               std.resid = TRUE,
                               method    = "LR",
                               collapse  = FALSE)

    #out$est.freq <- unlist(PI) * out$nobs

    out
}

# test for bivariate normality
lavCorTableFit <- function(object, ordered = NULL, group = NULL, 
                           ov.names.x = NULL, missing = "listwise") {

    # check object class
    if(inherits(object, "lavaan")) {
        lav.data <- object@lavData
    } else if(inherits(object, "lavData")) {
        lav.data <- object
    } else if(inherits(object, "data.frame")) {
        NAMES <- names(object)
        if(!is.null(group)) {
            NAMES <- NAMES[- match(group, NAMES)]
        }
        lav.data <- lavData(data = object, group = group,
                            ov.names = NAMES, ordered = ordered,
                            ov.names.x = ov.names.x,
                            missing = missing)
    } else {
        stop("lavaan ERROR: lavCor can not handle objects of class ",
             class(object))
    }

    COR <- lav_cor(lav.data = lav.data, WLS.W = FALSE, details = TRUE,
                   labels = FALSE, verbose = FALSE)

    # relist
    if(!is.list(COR)) {
        COR <- list(COR)
    }
    TH <- lapply(lapply(COR, attr, "TH"), unlist)
    TH.IDX <- lapply(lapply(COR, attr, "TH.IDX"), unlist)

    out <- lav_pairwise_tables_freq(vartable = lav.data@ov,
                                    X        = lav.data@X,
                                    ov.names = lav.data@ov.names,
                                    as.data.frame. = TRUE)

    # LR
    # not defined if out$obs.prop is (close to) zero
    # we use freq=0.5 for these empty cells
    zero.idx <- which(obs.prop < .Machine$double.eps)
    obs.prop[zero.idx] <- 0.5/out$nobs[zero.idx]
    out$std.resid <- 2*out$nobs*(obs.prop*log(obs.prop/est.prop))
    

    # only 1 row per table
    row.idx <- which(!duplicated(out$id))
    out <- out[row.idx,,drop=FALSE]

    # remove some cell-specific columns
    out$row <- NULL; out$col <- NULL
    out$obs.freq <- NULL; out$est.freq <- NULL

    # add table-wise info
    # FIXME: we need to filter out 'numeric' variables
    #out$cors <- unlist( lapply(COR, vech, diag=FALSE) )

    out
}
