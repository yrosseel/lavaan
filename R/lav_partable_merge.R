# merge two parameter tables
# - but allow different number of columns
lav_partable_merge <- function(pt1 = NULL, pt2 = NULL, 
                               remove.duplicated = FALSE,
                               fromLast=FALSE,
                               warn = TRUE) {

    pt1 <- as.data.frame(pt1, stringsAsFactors = FALSE)
    pt2 <- as.data.frame(pt2, stringsAsFactors = FALSE)

    # check minimum requirements: lhs, op, rhs
    stopifnot( !is.null(pt1$lhs), !is.null(pt1$op), !is.null(pt1$rhs),
               !is.null(pt2$lhs), !is.null(pt2$op), !is.null(pt2$rhs) )

    # both should have group (or not)
    if(is.null(pt1$group) && is.null(pt2$group)) {
        TMP <- rbind(pt1[, c("lhs","op","rhs","group")],
                     pt2[, c("lhs","op","rhs","group")])
    } else {
        if(is.null(pt1$group) && !is.null(pt2$group)) {
            pt1$group <- rep(1L, length(pt1$lhs))
        } else if(is.null(pt2$group) && !is.null(pt1$group)) {
            pt2$group <- rep(1L, length(pt2$lhs))
        }
        TMP <- rbind(pt1[, c("lhs","op","rhs","group")],
                     pt2[, c("lhs","op","rhs","group")])
    }

    # check for duplicated elements
    if(remove.duplicated) {
        # if fromLast = TRUE, idx is in pt1
        # if fromLast = FALSE, idx is in pt2
        idx <- which(duplicated(TMP, fromLast=fromLast)) 
    
        if(length(idx)) {
            if(warn) {
                warning("lavaan WARNING: duplicated parameters are ignored:\n",
                paste(apply(pt1[idx, c("lhs","op","rhs")], 1,
                      paste, collapse=" "), collapse="\n"))
            }
            if(fromLast) {
                pt1 <- pt1[-idx,]
            } else {
                idx <- idx - nrow(pt1)
                pt2 <- pt2[-idx,]
            }
        }
    } else if(!is.null(pt1$start) && !is.null(pt2$start)) {
        # copy start values from pt1 to pt2
        for(i in 1:length(pt1$lhs)) {
            idx <- which(pt2$lhs == pt1$lhs[i] &
                         pt2$op  == pt1$op[i] &
                         pt2$rhs == pt1$rhs[i] &
                         pt2$group == pt1$group[i])

            pt2$start[idx] <- pt1$start[i]
        }
    }

    # nicely merge, using 'id' column (if it comes first)
    if(is.null(pt1$id) && !is.null(pt2$id)) {
        nid <- max(pt2$id)
        pt1$id <- (nid+1L):(nid+nrow(pt1))
    } else if(is.null(pt2$id) && !is.null(pt1$id)) {
        nid <- max(pt1$id)
        pt2$id <- (nid+1L):(nid+nrow(pt2))
    }

    NEW <- base::merge(pt1, pt2, all = TRUE, sort = FALSE)

    NEW
}
