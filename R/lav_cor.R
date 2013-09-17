# user-visible routine to 
# compute polychoric/polyserial/... correlations
#
# YR 17 Sept 2013

lavCor <- function(object, ordered = NULL, # method="two.step", 
                   missing = "listwise",
                   WLS.W = FALSE, details = FALSE,
                   labels = TRUE, verbose = FALSE) {

    lav_data <- lav_data_extract(object = object, ordered = ordered)
    vartable   <- lav_data$vartable
    X          <- lav_data$X
    ov.names   <- lav_data$ov.names
    eXo        <- lav_data$eXo
    ov.names.x <- lav_data$ov.names.x

    # number of groups
    ngroups <- length(X)

    COR <- vector("list", length=ngroups)
    for(g in 1:ngroups) {
        ov.types  <- vartable$type[ match(ov.names[[g]], vartable$name) ]
        ov.levels <- vartable$nlev[ match(ov.names[[g]], vartable$name) ]
        CAT <- muthen1984(Data       = X[[g]],
                          ov.names   = ov.names[[g]],
                          ov.types   = ov.types,
                          ov.levels  = ov.levels,
                          ov.names.x = ov.names.x[[g]],
                          eXo        = eXo[[g]],
                          group      = g, # for error messages only
                          missing    = missing, # listwise or pairwise?
                          WLS.W      = WLS.W,
                          verbose    = verbose)
        COR[[g]] <- unname(CAT$COV)
        if(details) {
            attr(COR[[g]], "ov.names") <- ov.names[[g]]
            attr(COR[[g]], "ov.types") <- ov.types
            attr(COR[[g]], "TH")       <- CAT$TH
            attr(COR[[g]], "TH.names") <- CAT$TH.NAMES
            attr(COR[[g]], "slopes")   <- CAT$SLOPES
            attr(COR[[g]], "var")      <- CAT$VAR
            attr(COR[[g]], "cov")      <- CAT$COV
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
