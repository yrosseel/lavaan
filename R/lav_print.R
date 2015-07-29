## NOTE:
## round(1.2355, 3) = 1.236
## but 
## round(1.2345, 3) = 1.234
##
## perhaps we should add 0.0005 or something to avoid this?

print.lavaan.data.frame <- function(x, ..., nd=3) {

    ROW.NAMES <- rownames(x)
    y <- as.data.frame(lapply(x, function(x) {
                              if(is.numeric(x)) round(x, nd) else x}))
    rownames(y) <- ROW.NAMES

    if(!is.null(attr(x, "header"))) {
        cat("\n", attr(x, "header"), "\n\n", sep = "")
    }

    print(y, ...)
    invisible(x)
}


# prints only lower triangle of a symmetric matrix
print.lavaan.matrix.symmetric <- function(x, ..., nd=3) {
    # print only lower triangle of a symmetric matrix
    # this function was inspired by the `print.correlation' function
    # in package nlme 
    y <- x; y <- unclass(y)
    ll <- lower.tri(x, diag=TRUE)
    y[ll] <- format(round(x[ll], digits=nd)); y[!ll] <- ""
    if (!is.null(colnames(x))) {
      colnames(y) <- abbreviate(colnames(x), minlength = nd + 3)
    }
    print(y, ..., quote = FALSE)
    invisible(x)
}


print.lavaan.matrix <- function(x, ..., nd=3) {
    y <- unclass(x)
    if (!is.null(colnames(x))) {
      colnames(y) <- abbreviate(colnames(x), minlength = nd + 3)
    }
    print( round(y, nd), ... )
    invisible(x)
}

print.lavaan.vector <- function(x, ..., nd=3) {
    y <- unclass(x)
    #if(!is.null(names(x))) {
    #    names(y) <- abbreviate(names(x), minlength = nd + 3)
    #}
    print( round(y, nd), ... )
    invisible(x)
}

print.lavaan.character <- function(x) {
    cat(x)
    invisible(x)
}

print.lavaan.parameterEstimates <- function(x, ..., nd = 3L) {

    # format for numeric values
    num.format  <- paste("%", max(8, nd + 5), ".", nd, "f", sep = "")
    char.format <- paste("%", max(8, nd + 5), "s", sep="") 

    # output sections
    GSECTIONS <- c("Latent Variables", 
                   "Composites", 
                   "Regressions", 
                   "Covariances",
                   "Intercepts", 
                   "Thresholds", 
                   "Variances", 
                   "Scales y*",
                   "Group Weight",
                   "R-Square")
    ASECTIONS <- c("Defined Parameters", 
                   "Constraints")

    cat("\nParameter Estimates:\n\n")

    # header
    # 1. information
    # 2. se
    # 3. bootstrap requested/successful draws

    # 1.
    t0.txt <- sprintf("  %-40s", "Information")
    tmp.txt <- attr(x, "information")
    t1.txt <- sprintf("  %10s", paste(toupper(substring(tmp.txt,1,1)),
                     substring(tmp.txt,2), sep=""))
    cat(t0.txt, t1.txt, "\n", sep="")

    # 2.
    t0.txt <- sprintf("  %-31s", "Standard Errors")
    tmp.txt <- attr(x, "se")
    t1.txt <- sprintf("  %19s", paste(toupper(substring(tmp.txt,1,1)),
                                      substring(tmp.txt,2), sep=""))
    cat(t0.txt, t1.txt, "\n", sep="")

    # 3.
    if(attr(x, "se") == "bootstrap" && !is.null(attr(x, "bootstrap"))) {
        t0.txt <- sprintf("  %-40s", "Number of requested bootstrap draws")
        t1.txt <- sprintf("  %10i", attr(x, "bootstrap"))
        cat(t0.txt, t1.txt, "\n", sep="")
        t0.txt <- sprintf("  %-40s", "Number of successful bootstrap draws")
        t1.txt <- sprintf("  %10i", attr(x, "bootstrap.successful"))
        cat(t0.txt, t1.txt, "\n", sep="")
    }
    
    # number of groups
    if(is.null(x$group)) {
        ngroups <- 1L
        x$group <- rep(1L, length(x$lhs))
    } else {
        ngroups <- max(x$group)
    }

    # round to 3 digits after the decimal point
    y <- as.data.frame(
           lapply(x, function(x) {
               if(is.numeric(x)) { 
                   sprintf(num.format, x)   
               } else {
                   x
               }
           }),
           stringsAsFactors = FALSE)

    # always remove group/op/rhs/label/exo columns
    y$op <- y$group <- y$rhs <- y$label <- y$exo <- NULL

    # if standardized, remove std.nox column (space reasons only)
    y$std.nox <- NULL

    # convert to character matrix
    m <- as.matrix(format.data.frame(y, na.encode = FALSE, 
                                         justify = "right"))
   
    # use empty row names
    rownames(m) <- rep("", nrow(m))

    # handle se == 0.0
    if(!is.null(x$se)) {
        se.idx <- which(x$se == 0)
        if(length(se.idx) > 0L) {
            m[se.idx, "se"] <- ""
            if(!is.null(x$z)) {
                m[se.idx, "z"] <- ""
            }
            if(!is.null(x$pvalue)) {
                m[se.idx, "pvalue"] <- ""
            }
        }

        # handle se == NA
        se.idx <- which(is.na(x$se))
        if(length(se.idx) > 0L) {
            if(!is.null(x$z)) {
                m[se.idx, "z"] <- ""
            }
            if(!is.null(x$pvalue)) {
                m[se.idx, "pvalue"] <- ""
            }
        }
    }

    # rename some column names
    colnames(m)[ colnames(m) ==    "lhs" ] <- ""
    colnames(m)[ colnames(m) ==     "op" ] <- ""
    colnames(m)[ colnames(m) ==    "rhs" ] <- ""
    colnames(m)[ colnames(m) ==    "est" ] <- "Estimate"
    colnames(m)[ colnames(m) ==     "se" ] <- "Std.Err"
    colnames(m)[ colnames(m) ==      "z" ] <- "Z-value"
    colnames(m)[ colnames(m) == "pvalue" ] <- "P(>|z|)"
    colnames(m)[ colnames(m) == "std.lv" ] <- "Std.lv"
    colnames(m)[ colnames(m) == "std.all"] <- "Std.all"
    colnames(m)[ colnames(m) == "std.nox"] <- "Std.nox"
    colnames(m)[ colnames(m) == "prior"  ] <- "Prior"
 
    # format column names
    colnames(m) <- sprintf(char.format, colnames(m))

    # first the group-specific sections
    for(g in 1:ngroups) {

        # ov/lv names
        ov.names <- lavNames(x, "ov", group = g)
        lv.names <- lavNames(x, "lv", group = g)

        # group header
        if(ngroups > 1L) {
            group.label <- attr(x, "group.label")
            #if(g > 1L) { 
                cat("\n\n")
            #} else {
            #    cat("\n")
            #}
            cat("Group ", g, " [", group.label[g], "]:\n", sep="")
        }

        for(s in GSECTIONS) {
            if(s == "Latent Variables") {
                row.idx <- which( x$op == "=~" & !x$lhs %in% ov.names & 
                                  x$group == g)
                if(length(row.idx) == 0L) next
                m[row.idx,1] <- .makeNames(x$rhs[row.idx], x$label[row.idx])
            } else if(s == "Composites") {
                row.idx <- which( x$op == "<~" & x$group == g)
                if(length(row.idx) == 0L) next
                m[row.idx,1] <- .makeNames(x$rhs[row.idx], x$label[row.idx])
            } else if(s == "Regressions") {
                row.idx <- which( x$op == "~" & x$group == g)
                if(length(row.idx) == 0L) next
                m[row.idx,1] <- .makeNames(x$rhs[row.idx], x$label[row.idx])
            } else if(s == "Covariances") {
                row.idx <- which(x$op == "~~" & x$lhs != x$rhs & !x$exo &
                                 x$group == g)
                if(length(row.idx) == 0L) next
                m[row.idx,1] <- .makeNames(x$rhs[row.idx], x$label[row.idx])
            } else if(s == "Intercepts") {
                row.idx <- which(x$op == "~1" & !x$exo & x$group == g)
                if(length(row.idx) == 0L) next
                m[row.idx,1] <- .makeNames(x$lhs[row.idx], x$label[row.idx])
            } else if(s == "Thresholds") {
                row.idx <- which(x$op == "|" & x$group == g)
                if(length(row.idx) == 0L) next
                m[row.idx,1] <- .makeNames(paste(x$lhs[row.idx], "|", 
                    x$rhs[row.idx], sep=""), x$label[row.idx])
            } else if(s == "Variances") {
                row.idx <- which(x$op == "~~" & x$lhs == x$rhs & !x$exo &
                                 x$group == g)
                if(length(row.idx) == 0L) next
                m[row.idx,1] <- .makeNames(x$rhs[row.idx], x$label[row.idx])
            } else if(s == "Scales y*") {
                row.idx <- which(x$op == "~*~" & x$group == g)
                if(length(row.idx) == 0L) next
                m[row.idx,1] <- .makeNames(x$rhs[row.idx], x$label[row.idx])
            } else if(s == "Group Weight") {
                row.idx <- which(x$lhs == "group" & x$op == "%" & x$group == g)
                if(length(row.idx) == 0L) next
                m[row.idx,1] <- .makeNames(x$rhs[row.idx], x$label[row.idx])
            } else if(s == "R-Square") {
                row.idx <- which(x$op == "r2" & x$group == g)
                if(length(row.idx) == 0L) next
                m[row.idx,1] <- .makeNames(x$rhs[row.idx], x$label[row.idx])
            } else {
                row.idx <- integer(0L)
            }

            if(s %in% c("Latent Variables", "Composites", 
                        "Regressions", "Covariances")) {
                nel <- length(row.idx)
                M <- matrix("", nrow = nel*2, ncol = ncol(m))
                colnames(M) <- colnames(m)
                rownames(M) <- rep("", NROW(M))
                #colnames(M)[1] <- sprintf("%-17s", paste(s, ":", sep = ""))
                LHS <- paste(x$lhs[row.idx], x$op[row.idx])
                lhs.idx <- seq(1, nel*2L, 2L)
                rhs.idx <- seq(1, nel*2L, 2L) + 1L
                M[lhs.idx, 1] <- sprintf(" %-15s", LHS)
                M[rhs.idx,  ] <- m[row.idx,]
                # avoid duplicated LHS labels
                if(nel > 1L) {
                    del.idx <- integer(0)
                    old.lhs <- ""
                    for(i in 1:nel) {
                        if(LHS[i] == old.lhs) {
                            del.idx <- c(del.idx, lhs.idx[i])
                        }
                        old.lhs <- LHS[i]
                    }
                    M <- M[-del.idx,,drop=FALSE]
                }
                cat("\n", s, ":\n", sep = "")
                #cat("\n")
                print(M, quote = FALSE)
            } else if(s == "R-Square") {
                M <- m[row.idx,1:2,drop=FALSE]
                colnames(M) <- colnames(m)[1:2]
                rownames(M) <- rep("", NROW(M))
                #colnames(M)[1] <- sprintf("%-17s", paste(s, ":", sep = ""))
                cat("\n", s, ":\n", sep = "")
                #cat("\n")
                print(M, quote = FALSE)
            } else {
                #M <- rbind(matrix("", nrow = 1L, ncol = ncol(m)),
                #           m[row.idx,])
                M <- m[row.idx,,drop=FALSE]
                colnames(M) <- colnames(m)
                rownames(M) <- rep("", NROW(M))
                #colnames(M)[1] <- sprintf("%-17s", paste(s, ":", sep = ""))
                cat("\n", s, ":\n", sep = "")
                #cat("\n")
                print(M, quote = FALSE)
            }
        }

    }    

    # asections
    for(s in ASECTIONS) {
        if(s == "Defined Parameters") {
            row.idx <- which(x$op == ":=")
            m[row.idx,1] <- .makeNames(x$lhs[row.idx], "")
            M <- m[row.idx,,drop=FALSE]
            colnames(M) <- colnames(m)
        } else if(s == "Constraints") {
            row.idx <- which(x$op %in% c("==", "<", ">"))
            if(length(row.idx) == 0) next
            m[row.idx,1] <- .makeConNames(x$lhs[row.idx], x$op[row.idx], 
                                          x$rhs[row.idx], nd = nd)
            m[row.idx,2] <- sprintf(num.format, abs(x$est[row.idx]))
            M <- m[row.idx,1:2,drop=FALSE]
            colnames(M) <- c("", sprintf(char.format, "|Slack|"))
        } else {
            row.idx <- integer(0L)
        }

        if(length(row.idx) == 0L) {
            next
        }

        rownames(M) <- rep("", NROW(M))
        #colnames(M)[1] <- sprintf("%-17s", paste(s, ":", sep = ""))
        #cat("\n")
        cat("\n", s, ":\n", sep = "")
        print(M, quote = FALSE)
    }
    cat("\n")

    invisible(m)
}

.makeNames <- function(NAMES, LABELS) {

    W <- 14

    multiB <- FALSE
    if(any(nchar(NAMES) != nchar(NAMES, "bytes")))
        multiB <- TRUE
    if(any(nchar(LABELS) != nchar(LABELS, "bytes")))
        multiB <- TRUE
    # labels?
    l.idx <- which(nchar(LABELS) > 0L)
    if(length(l.idx) > 0L) {
        if(!multiB) {
            LABELS <- abbreviate(LABELS, 4)
            LABELS[l.idx] <- paste(" (", LABELS[l.idx], ")", sep="")
            MAX.L <- max(nchar(LABELS))
            NAMES <- abbreviate(NAMES, minlength = (W - MAX.L),
                                strict = TRUE)
        } else {
            # do not abbreviate anything (eg in multi-byte locales)
            MAX.L <- 4L
        }
        NAMES <- sprintf(paste("%-", (W - MAX.L), "s%", MAX.L, "s",
                               sep=""), NAMES, LABELS)
    } else {
        if(!multiB) {
            NAMES <- abbreviate(NAMES, minlength = W, strict = TRUE)
        } else {
            NAMES <- sprintf(paste("%-", W, "s", sep=""), NAMES)
        }
    }

    char.format <- paste("   %-", W, "s", sep = "")
    sprintf(char.format, NAMES)
}

.makeConNames <- function(lhs, op, rhs, nd) {

    nd <- max(nd, 3)
    W <- 40 + (nd - 3)*3

    nel <- length(lhs)
    if(length(nel) == 0L) return(character(0))
    NAMES <- character(nel)
    for(i in 1:nel) {
        if(rhs == "0" && op == ">") {
            con.string <- paste(lhs, " - 0", sep="")
        } else if(rhs == "0" && op == "<") {
            con.string <- paste(rhs, " - (", lhs, ")", sep="")
        } else if(rhs != "0" && op == ">") {
            con.string <- paste(lhs, " - (", rhs, ")", sep="")
        } else if(rhs != "0" && op == "<") {
            con.string <- paste(rhs, " - (", lhs, ")", sep="")
        } else if(rhs == "0" && op == "==") {
            con.string <- paste(lhs, " - 0", sep="")
        } else if(rhs != "0" && op == "==") {
            con.string <- paste(lhs, " - (", rhs, ")", sep="")
        }
        con.string <- abbreviate(con.string, W, strict = TRUE)
        char.format <- paste("    %-", W, "s", sep = "")
        NAMES[i] <- sprintf(char.format, con.string)
    }

    NAMES
}
