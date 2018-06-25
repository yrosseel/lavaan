# return 'attributes' of a lavaan partable -- generate a new set if necessary
lav_partable_attributes <- function(partable, pta = NULL) {

    if(is.null(pta)) {
        # attached to partable?
        pta <- attributes(partable)
        if(!is.null(pta$vnames) && !is.null(pta$nvar)) {
            # looks like a pta
            return(pta)
        } else {
            pta <- list()
        }
    }

    # vnames
    pta$vnames <- lav_partable_vnames(partable, type="all")

    # vidx
    OV <- pta$vnames$ov
    LV <- pta$vnames$lv
    nblocks <- length(pta$vnames$ov)
    pta$vidx <- lapply(names(pta$vnames), function(v) {
                    lapply(seq_len(nblocks), function(g) {
                        if(grepl("lv", v)) {
                            match(pta$vnames[[v]][[g]], LV[[g]])
                        } else if(grepl("th", v)) {
                            # thresholds have '|t' pattern
                            TH <-  sapply(strsplit(pta$vnames[[v]][[g]],
                                          "|t", fixed = TRUE), "[[", 1L)
                            match(TH, OV[[g]])
                        } else if(grepl("eqs", v)){
                            # mixture of OV/LV
                            integer(0L)
                        } else {
                            match(pta$vnames[[v]][[g]], OV[[g]])
                        }
                    })
                })
    names(pta$vidx) <- names(pta$vnames)

    # meanstructure
    pta$meanstructure <- any(partable$op == "~1")

    # nblocks
    pta$nblocks <- nblocks

    # ngroups
    pta$ngroups <- lav_partable_ngroups(partable)

    # nlevels
    pta$nlevels <- lav_partable_nlevels(partable)

    # nvar
    pta$nvar <- lapply(pta$vnames$ov, length)

    # nfac
    pta$nfac <- lapply(pta$vnames$lv, length)

    # nfac.nonnormal - for numerical integration
    pta$nfac.nonnormal <- lapply(pta$vnames$lv.nonnormal, length)

    # th.idx (new in 0.6-1)
    pta$th.idx <- lapply(seq_len(pta$nblocks), function(g) {
                            out <- numeric( length(pta$vnames$th.mean[[g]]) )
                            idx <- ( pta$vnames$th.mean[[g]] %in%
                                     pta$vnames$th[[g]] )
                            out[idx] <- pta$vidx$th[[g]]
                            out
                        })

    pta
}

