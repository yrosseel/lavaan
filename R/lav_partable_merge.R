# merge two parameter tables
# - but allow different number of columns
lav_partable_merge <- function(pt1 = NULL, pt2 = NULL, 
                               remove.duplicated = FALSE,
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
    idx <- which(duplicated(TMP, fromLast=TRUE)) # idx should be in pt1
    if(length(idx)) {
        if(warn) {
            warning("lavaan WARNING: duplicated parameters are ignored:\n",
            paste(apply(pt1[idx, c("lhs","op","rhs")], 1,
                  paste, collapse=" "), collapse="\n"))
        l}
        pt1 <- pt1[-idx,]
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
