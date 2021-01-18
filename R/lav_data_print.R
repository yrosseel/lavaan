# print object from lavData class
#

setMethod("show", "lavData",
function(object) {
    # print 'lavData' object
    lav_data_print_short(object)
})

lav_data_print_short <- function(object, nd = 3L) {

    lavdata <- object
    num.format  <- paste("%", max(8L, nd + 5L), ".", nd, "f", sep = "")

    # listwise deletion?
    listwise <- twocolumn <- FALSE
    for(g in 1:lavdata@ngroups) {
       if(lavdata@nobs[[1L]] != lavdata@norig[[1L]]) {
           listwise <- twocolumn <- TRUE
           break
       }
    }

    # header? no, for historical reasons only
    #cat("Data information:\n\n")

    c1 <- c2 <- c3 <- character(0L)

    # number of observations
    if(lavdata@ngroups == 1L) {
        if(listwise) {
            c1 <- c(c1, ""); c2 <- c(c2, "Used"); c3 <- c(c3, "Total")
        }
        c1 <- c(c1, "Number of observations")
        c2 <- c(c2, lavdata@nobs[[1L]])
        c3 <- c(c3, ifelse(listwise, lavdata@norig[[1L]], ""))
    } else {
        c1 <- c(c1, "Number of observations per group:");
        if(listwise) {
            c2 <- c(c2, "Used"); c3 <- c(c3, "Total")
        } else {
            c2 <- c(c2, ""); c3 <- c(c3, "")
        }
        for(g in 1:lavdata@ngroups) {
            c1 <- c(c1, sprintf("  %-40s", lavdata@group.label[[g]]))
            c2 <- c(c2, lavdata@nobs[[g]])
            c3 <- c(c3, ifelse(listwise, lavdata@norig[[g]], ""))
        } # g
    }

    # number of clusters
    if(lavdata@ngroups == 1L) {
        if( (.hasSlot(lavdata, "nlevels")) && # in case we have an old obj
            (lavdata@nlevels > 1L) ) {
            for(l in 2:lavdata@nlevels) {
                c1 <- c(c1,
                        paste("Number of clusters [", lavdata@cluster[l-1], "]",
                        sep = ""))
                c2 <- c(c2, lavdata@Lp[[1]]$nclusters[[l]])
                c3 <- c(c3, "")
            }
        } else if( (.hasSlot(lavdata, "cluster")) &&
                   (length(lavdata@cluster) > 0L) ) {
            c1 <- c(c1, paste("Number of clusters [", lavdata@cluster, "]",
                              sep = ""))
            c2 <- c(c2, lavdata@Lp[[1]]$nclusters[[2]])
            c3 <- c(c3, "")
        }
    } else {
        if( (.hasSlot(lavdata, "nlevels")) && (lavdata@nlevels > 1L) ) {
            for(l in 2:lavdata@nlevels) {
                c1 <- c(c1,
                  paste("Number of clusters [", lavdata@cluster[l-1], "]:",
                        sep = ""))
                c2 <- c(c2, ""); c3 <- c(c3, "")
                for(g in 1:lavdata@ngroups) {
                    c1 <- c(c1, sprintf("  %-40s", lavdata@group.label[[g]]))
                    c2 <- c(c2, lavdata@Lp[[g]]$nclusters[[l]])
                    c3 <- c(c3, "")
                }
            }
        } else if( (.hasSlot(lavdata, "cluster")) &&
               (length(lavdata@cluster) > 0L) ) {
            c1 <- c(c1,
             paste("Number of clusters [", lavdata@cluster, "]:", sep = ""))
            c2 <- c(c2, ""); c3 <- c(c3, "")
            for(g in 1:lavdata@ngroups) {
                c1 <- c(c1, sprintf("  %-40s", lavdata@group.label[[g]]))
                c2 <- c(c2, lavdata@Lp[[g]]$nclusters[[2]])
                c3 <- c(c3, "")
            }
        }
    }

    # missing patterns?
    if(!is.null(lavdata@Mp[[1L]])) {
        if(lavdata@ngroups == 1L) {
            c1 <- c(c1, "Number of missing patterns")
            c2 <- c(c2, lavdata@Mp[[1L]]$npatterns)
            c3 <- c(c3, "")
        } else {
            c1 <- c(c1, "Number of missing patterns per group:")
            c2 <- c(c2, ""); c3 <- c(c3, "")
            for(g in 1:lavdata@ngroups) {
                c1 <- c(c1, sprintf("  %-40s", lavdata@group.label[[g]]))
                c2 <- c(c2, lavdata@Mp[[g]]$npatterns)
                c3 <- c(c3, "")
            }
        }
    }

    # sampling weights?
    if( (.hasSlot(lavdata, "weights")) && # in case we have an old object
        (!is.null(lavdata@weights[[1L]])) ) {
        c1 <- c(c1, "Sampling weights variable")
        c2 <- c(c2, lavdata@sampling.weights)
        c3 <- c(c3, "")
    }

    # empty row
    c1 <- c(c1, ""); c2 <- c(c2, ""); c3 <- c(c3, "")

    # format c1/c2
    c1 <- format(c1, width = 43L)
    c2 <- format(c2, width = 8L + max(0, (nd - 3L)) * 4L, justify = "right")
    c3 <- format(c3, width = 8L + nd, justify = "right")

    # create character matrix
    if(twocolumn) {
        M <- cbind(c1, c2, c3, deparse.level = 0)
    } else {
        M <- cbind(c1, c2, deparse.level = 0)
    }
    colnames(M) <- rep("",  ncol(M))
    rownames(M) <- rep(" ", nrow(M))

    # print
    write.table(M, row.names = TRUE, col.names = FALSE, quote = FALSE)

    invisible(M)
}

