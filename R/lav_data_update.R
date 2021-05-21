# update lavdata object with new dataset
# - assuming everything else stays the same
# - optionally, also provide boot.idx (per group) to adapt internal slots

# YR - 18 Jan 2021 (so we don't need to export lav_data_*_patterns functions)

lav_data_update <- function(lavdata = NULL, newX = NULL, BOOT.idx = NULL,
                            lavoptions = NULL) {

    stopifnot(length(newX) == lavdata@ngroups)
    stopifnot(!is.null(lavoptions))
    newdata <- lavdata

    # replace data 'X' slot for each group
    for(g in 1:lavdata@ngroups) {

        # replace raw data
        newdata@X[[g]] <- newX[[g]]

        # Mp
        if(lavoptions$missing != "listwise") {
            newdata@Mp[[g]] <- lav_data_missing_patterns(newX[[g]],
                                  sort.freq = FALSE, coverage = FALSE)
        }

        # Rp
        if(length(lavdata@ov.names.x[[g]]) == 0L &&
           all(lavdata@ov.names[[g]] %in%
               lavdata@ov$name[lavdata@ov$type == "ordered"])) {
            newdata@Rp[[g]] <- lav_data_resp_patterns(newX[[g]])
        }

        # Lp
        if(lavdata@nlevels > 1L) {
            # extract cluster variable(s), for this group
            clus <- matrix(0, nrow(newX[[g]]), lavdata@nlevels - 1L)
            for(l in 2:lavdata@nlevels) {
                clus[,(l-1L)] <- lavdata@Lp[[g]]$cluster.idx[[l]]
            }
            newdata@Lp[[g]] <- lav_data_cluster_patterns(Y = newX[[g]],
                                clus = clus,
                                cluster = lavdata@cluster,
                                ov.names = lavdata@ov.names[[g]],
                                ov.names.l = lavdata@ov.names.l[[g]])
        }
    }

    # if boot.idx if provided, also adapt eXo and WT
    if(!is.null(BOOT.idx)) {

        boot.idx <- BOOT.idx[[g]]

        # eXo
        if(!is.null(lavdata@eXo[[g]])) {
            newdata@eXo[[g]] <- lavdata@eXo[[g]][boot.idx,,drop=FALSE]
        }

        # sampling weights
        if(!is.null(lavdata@weights[[g]])) {
            newdata@weights[[g]] <- lavdata@weights[[g]][boot.idx]
        }
    } # g

    # return update data object
    newdata
}

