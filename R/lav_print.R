## NOTE:
## round(1.2355, 3) = 1.236
## but
## round(1.2345, 3) = 1.234
##
## perhaps we should add 0.0005 or something to avoid this?

print.lavaan.data.frame <- function(x, ..., nd = 3L) {

    ROW.NAMES <- rownames(x)
    y <- as.data.frame(lapply(x, function(x) {
                              if(is.numeric(x)) round(x, nd) else x}))
    rownames(y) <- ROW.NAMES

    if(!is.null(attr(x, "header"))) {
        cat("\n", attr(x, "header"), "\n\n", sep = "")
    }

    print(y, ...)

    if(!is.null(attr(x, "footer"))) {
        cat("\n", attr(x, "footer"), "\n\n", sep = "")
    }

    invisible(x)
}

print.lavaan.list <- function(x, ...) {

    y <- unclass(x)
    attr(y, "header") <- NULL

    header <- attr(x, "header")
    if(!is.null(header)) {
        if(is.character(header)) {
            cat("\n", header, "\n\n", sep = "")
        } else {
            print(header); cat("\n")
        }
    }

    print(y, ...)
    invisible(x)
}


# prints only lower triangle of a symmetric matrix
print.lavaan.matrix.symmetric <- function(x, ..., nd = 3L) {
    # print only lower triangle of a symmetric matrix
    # this function was inspired by the `print.correlation' function
    # in package nlme
    y <- x; y <- unclass(y)
    ll <- lower.tri(x, diag=TRUE)
    y[ll] <- format(round(x[ll], digits=nd)); y[!ll] <- ""
    if (!is.null(colnames(x))) {
      colnames(y) <- abbreviate(colnames(x), minlength = nd + 3L)
    }
    print(y, ..., quote = FALSE)
    invisible(x)
}


print.lavaan.matrix <- function(x, ..., nd = 3L) {
    y <- unclass(x)
    if (!is.null(colnames(x))) {
      colnames(y) <- abbreviate(colnames(x), minlength = nd + 3L)
    }
    print( round(y, nd), ... )
    invisible(x)
}

print.lavaan.vector <- function(x, ..., nd = 3L) {
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
    num.format  <- paste("%", max(8L, nd + 5L), ".", nd, "f", sep = "")
    int.format  <- paste("%", max(8L, nd + 5L), "d", sep = "")
    char.format <- paste("%", max(8L, nd + 5L), "s", sep = "")

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

    # header?
    header <- attr(x, "header")
    if(is.null(header)) {
        header <- FALSE
    }

    if(header) {
        cat("\nParameter Estimates:\n\n")

        # info about standard errors (if we have x$se only)
        # 1. information
        # 2. se
        # 3. bootstrap requested/successful draws
        if(!is.null(x$se)) {

            # container
            c1 <- c2 <- character(0L)

            # which type of standard errors?
            c1 <- c(c1, "Standard errors")
            if(attr(x, "se") == "robust.huber.white") {
                tmp.txt <- "sandwich" # since 0.6-6
            } else {
                tmp.txt <- attr(x, "se")
            }
            c2 <- c(c2, paste(toupper(substring(tmp.txt, 1, 1)),
                                      substring(tmp.txt, 2), sep = ""))
            # information
            if(attr(x, "se") != "bootstrap") {

                # type for information
                if(attr(x, "se") == "robust.huber.white") {
                    c1 <- c(c1, "Information bread")
                } else {
                    c1 <- c(c1, "Information")
                }
                tmp.txt <- attr(x, "information")
                c2 <- c(c2, paste(toupper(substring(tmp.txt, 1, 1)),
                                  substring(tmp.txt, 2), sep = ""))

                # if observed, which type? (hessian of h1)
                if(attr(x, "information") == "observed") {
                    c1 <- c(c1, "Observed information based on")
                    tmp.txt <- attr(x, "observed.information")
                    c2 <- c(c2, paste(toupper(substring(tmp.txt, 1, 1)),
                                      substring(tmp.txt, 2), sep = ""))
                }

                # if h1 is involved, structured or unstructured?
                if(attr(x, "information") %in% c("expected", "first.order") ||
                   attr(x, "observed.information") == "h1") {
                    if(attr(x, "se") == "robust.huber.white" &&
                       attr(x, "h1.information") !=
                       attr(x, "h1.information.meat")) {
                        c1 <- c(c1, "Information bread saturated (h1) model")
                    } else {
                        c1 <- c(c1, "Information saturated (h1) model")
                    }
                    tmp.txt <- attr(x, "h1.information")
                    c2 <- c(c2, paste(toupper(substring(tmp.txt,1,1)),
                                      substring(tmp.txt, 2), sep = ""))
                }

                # if sandwich, which information for the meat? (first.order)
                # only print if it is NOT first.order
                if(attr(x, "se") == "robust.huber.white" &&
                   attr(x, "information.meat") != "first.order") {
                    c1 <- c(c1, "Information meat")
                    tmp.txt <- attr(x, "information.meat")
                    c2 <- c(c2, paste(toupper(substring(tmp.txt,1,1)),
                                      substring(tmp.txt, 2), sep = ""))
                }

                # if sandwich, structured or unstructured for the meat?
                # only print if it is not the same as h1.information
                if(attr(x, "se") == "robust.huber.white" &&
                   attr(x, "h1.information.meat") !=
                   attr(x, "h1.information")) {
                    c1 <- c(c1, "Information meat saturated (h1) model")
                    tmp.txt <- attr(x, "h1.information.meat")
                    c2 <- c(c2, paste(toupper(substring(tmp.txt,1,1)),
                                      substring(tmp.txt, 2), sep = ""))
                }
            } # no bootstrap

            # 4.
            if(attr(x, "se") == "bootstrap" && !is.null(attr(x, "bootstrap"))) {
                c1 <- c(c1, "Number of requested bootstrap draws")
                c2 <- c(c2, attr(x, "bootstrap"))
                c1 <- c(c1, "Number of successful bootstrap draws")
                c2 <- c(c2, attr(x, "bootstrap.successful"))
            }

            # format c1/c2
            c1 <- format(c1, width = 38L)
            c2 <- format(c2,
                   width = 13L + max(0, (nd - 3L)) * 4L, justify = "right")

            # create character matrix
            M <- cbind(c1, c2, deparse.level = 0)
            colnames(M) <- rep("",  ncol(M))
            rownames(M) <- rep(" ", nrow(M))

            # print
            write.table(M, row.names = TRUE, col.names = FALSE, quote = FALSE)
        }
    }

    # number of groups
    if(is.null(x$group)) {
        ngroups <- 1L
        x$group <- rep(1L, length(x$lhs))
    } else {
        ngroups <- lav_partable_ngroups(x)
    }

    # number of levels
    if(is.null(x$level)) {
        nlevels <- 1L
        x$level <- rep(1L, length(x$lhs))
    } else {
        nlevels <- lav_partable_nlevels(x)
    }

    # block column
    if(is.null(x$block)) {
        x$block <- rep(1L, length(x$lhs))
    }

    # step column (SAM)
    #if(!is.null(x$step)) {
    #    tmp.LABEL <- rep("", length(x$lhs))
    #    p1.idx <- which(x$step == 1L)
    #    p2.idx <- which(x$step == 2L)
    #    tmp.LABEL[p1.idx] <- "1"
    #    tmp.LABEL[p2.idx] <- "2"
    #
    #    if(is.null(x$label)) {
    #        x$label <- tmp.LABEL
    #    } else {
    #        x$label <- paste(x$label, tmp.LABEL, sep = "")
    #    }
    #
    #    x$step <- NULL
    #}

    # round to 3 digits after the decimal point
    y <- as.data.frame(
           lapply(x, function(x) {
               if(is.integer(x)) {
                   sprintf(int.format, x)
               } else if(is.numeric(x)) {
                   sprintf(num.format, x)
               } else {
                   x
               }
           }),
           stringsAsFactors = FALSE)

    # always remove /block/level/group/op/rhs/label/exo columns
    y$op <- y$group <- y$rhs <- y$label <- y$exo <- NULL
    y$block <- y$level <- NULL
    y$efa <- NULL

    # if standardized, remove std.nox column (space reasons only)
    # unless, std.all is already removed
    if(!is.null(y$std.all)) {
        y$std.nox <- NULL
    }

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
            ## for lavaan.mi-class objects (semTools)
            if(!is.null(x$t)) {
                m[se.idx, "t"] <- ""
            }
            if(!is.null(x$df)) {
                m[se.idx, "df"] <- ""
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
            ## for lavaan.mi-class objects (semTools)
            if(!is.null(x$t)) {
                m[se.idx, "t"] <- ""
            }
            if(!is.null(x$df)) {
                m[se.idx, "df"] <- ""
            }
        }
    }

    # handle lower/upper boundary points
    if(!is.null(x$lower)) {
        b.idx <- which( abs(x$lower - x$est) < sqrt(.Machine$double.eps) &
                        (is.na(x$se) | (is.finite(x$se) & x$se != 0.0)) )
        if(length(b.idx) > 0L && !is.null(x$pvalue)) {
            m[b.idx, "pvalue"] <- ""
            if(is.null(x$label)) {
                x$label <- rep("", length(x$lhs))
            }
            x$label[b.idx] <- ifelse(nchar(x$label[b.idx]) > 0L,
                                     paste(x$label[b.idx], "+lb", sep = ""),
                                     "lb")
        }
        # remove lower column
        m <- m[, colnames(m) != "lower"]
    }
    if(!is.null(x$upper)) {
        b.idx <- which( abs(x$upper - x$est) < sqrt(.Machine$double.eps) &
                        is.finite(x$se) & x$se != 0.0)
        if(length(b.idx) > 0L && !is.null(x$pvalue)) {
            m[b.idx, "pvalue"] <- ""
            if(is.null(x$label)) {
                x$label <- rep("", length(x$lhs))
            }
            x$label[b.idx] <- ifelse(nchar(x$label[b.idx]) > 0L,
                                     paste(x$label[b.idx], "+ub", sep = ""),
                                     "ub")
        }
        # remove upper column
        m <- m[, colnames(m) != "upper"]
    }


    # handle fmi
    if(!is.null(x$fmi)) {
        se.idx <- which(x$se == 0)
        if(length(se.idx) > 0L) {
            m[se.idx, "fmi"] <- ""
            ## for lavaan.mi-class objects (semTools)
            if (!is.null(x$riv)) m[se.idx, "riv"] <- ""
        }

        not.idx <- which(x$op %in% c(":=", "<", ">", "=="))
        if(length(not.idx) > 0L) {
            if(!is.null(x$fmi)) {
                m[not.idx, "fmi"] <- ""
                ## for lavaan.mi-class objects (semTools)
                if (!is.null(x$riv)) m[not.idx, "riv"] <- ""
            }
        }
    }

    # for blavaan, handle Post.SD and PSRF
    if(!is.null(x$Post.SD)) {
        se.idx <- which(x$Post.SD == 0)
        if(length(se.idx) > 0L) {
            m[se.idx, "Post.SD"] <- ""
            if(!is.null(x$psrf)) {
                m[se.idx, "psrf"] <- ""
            }
            if(!is.null(x$PSRF)) {
                m[se.idx, "PSRF"] <- ""
            }
        }

        # handle psrf for defined parameters
        not.idx <- which(x$op %in% c(":=", "<", ">", "=="))
        if(length(not.idx) > 0L) {
            if(!is.null(x$psrf)) {
                m[not.idx, "psrf"] <- ""
            }
            if(!is.null(x$PSRF)) {
                m[not.idx, "PSRF"] <- ""
            }
        }
    }

    # rename some column names
    colnames(m)[ colnames(m) ==    "lhs" ] <- ""
    colnames(m)[ colnames(m) ==     "op" ] <- ""
    colnames(m)[ colnames(m) ==    "rhs" ] <- ""
    colnames(m)[ colnames(m) ==   "step" ] <- "Step"
    colnames(m)[ colnames(m) ==    "est" ] <- "Estimate"
    colnames(m)[ colnames(m) ==     "se" ] <- "Std.Err"
    colnames(m)[ colnames(m) ==      "z" ] <- "z-value"
    colnames(m)[ colnames(m) == "pvalue" ] <- "P(>|z|)"
    colnames(m)[ colnames(m) == "std.lv" ] <- "Std.lv"
    colnames(m)[ colnames(m) == "std.all"] <- "Std.all"
    colnames(m)[ colnames(m) == "std.nox"] <- "Std.nox"
    colnames(m)[ colnames(m) == "prior"  ] <- "Prior"
    colnames(m)[ colnames(m) == "fmi"    ] <- "FMI"
    ## for lavaan.mi-class objects (semTools)
    if ("t" %in% colnames(m)) {
      colnames(m)[ colnames(m) == "t"      ] <- "t-value"
      colnames(m)[ colnames(m) == "P(>|z|)"] <- "P(>|t|)"
      colnames(m)[ colnames(m) == "riv"    ] <- "RIV"
    }

    # format column names
    colnames(m) <- sprintf(char.format, colnames(m))

    # exceptions for blavaan: Post.Mean (width = 9), Prior (width = 14)
    if(!is.null(x$Post.Mean)) {
        tmp <- gsub("[ \t]+", "", colnames(m), perl=TRUE)

        # reformat "Post.Mean" column
        col.idx <- which(tmp == "Post.Mean")
        if(length(col.idx) > 0L) {
            tmp.format <- paste("%", max(9, nd + 5), "s", sep="")
            colnames(m)[col.idx] <- sprintf(tmp.format, colnames(m)[col.idx])
            m[,col.idx] <- sprintf(tmp.format, m[,col.idx])
        }

        # reformat "Prior" column
        col.idx <- which(tmp == "Prior")
        if(length(col.idx) > 0L) {
            MAX <- max( nchar( m[,col.idx] ) ) + 1L
            tmp.format <- paste("%", max(MAX, nd + 5), "s", sep="")
            colnames(m)[col.idx] <- sprintf(tmp.format, colnames(m)[col.idx])
            m[,col.idx] <- sprintf(tmp.format, m[,col.idx])
        }
    }

    b <- 0L
    # group-specific sections
    for(g in 1:ngroups) {

        # group header
        if(ngroups > 1L) {
            group.label <- attr(x, "group.label")
            cat("\n\n")
            cat("Group ", g, " [", group.label[g], "]:\n", sep="")
        }

        for(l in 1:nlevels) {

            # block number
            b <- b + 1L

            # ov/lv names
            ov.names <- lavNames(x, "ov", block = b)
            lv.names <- lavNames(x, "lv", block = b)

            # level header
            if(nlevels > 1L) {
                level.label <- attr(x, "level.label")
                cat("\n\n")
                cat("Level ", l, " [", level.label[l], "]:\n", sep="")
            }

            # group-specific sections
            for(s in GSECTIONS) {
                if(s == "Latent Variables") {
                        row.idx <- which( x$op == "=~" & !x$lhs %in% ov.names &
                                      x$block == b)
                    if(length(row.idx) == 0L) next
                    m[row.idx,1] <- .makeNames(x$rhs[row.idx], x$label[row.idx])
                } else if(s == "Composites") {
                    row.idx <- which( x$op == "<~" & x$block == b)
                    if(length(row.idx) == 0L) next
                    m[row.idx,1] <- .makeNames(x$rhs[row.idx], x$label[row.idx])
                } else if(s == "Regressions") {
                    row.idx <- which( x$op == "~" & x$block == b)
                    if(length(row.idx) == 0L) next
                    m[row.idx,1] <- .makeNames(x$rhs[row.idx], x$label[row.idx])
                } else if(s == "Covariances") {
                    row.idx <- which(x$op == "~~" & x$lhs != x$rhs & !x$exo &
                                     x$block == b)
                    if(length(row.idx) == 0L) next
                    # make distinction between residual and plain
                    y.names <- unique( c(lavNames(x, "eqs.y"),
                                         lavNames(x, "ov.ind"),
                                         lavNames(x, "lv.ind")) )
                    PREFIX <- rep("", length(row.idx))
                    PREFIX[ x$rhs[row.idx] %in% y.names ] <- "  ."
                    m[row.idx,1] <- .makeNames(x$rhs[row.idx], x$label[row.idx],
                                               PREFIX = PREFIX)
                    #m[row.idx,1] <- .makeNames(x$rhs[row.idx], x$label[row.idx])
                } else if(s == "Intercepts") {
                    row.idx <- which(x$op == "~1" & !x$exo & x$block == b)
                    if(length(row.idx) == 0L) next
                    # make distinction between intercepts and means
                    y.names <- unique( c(lavNames(x, "eqs.y"),
                                         lavNames(x, "ov.ind"),
                                         lavNames(x, "lv.ind")) )
                    PREFIX <- rep("", length(row.idx))
                    PREFIX[ x$lhs[row.idx] %in% y.names ] <- "  ."
                        m[row.idx,1] <- .makeNames(x$lhs[row.idx], x$label[row.idx],
                                               PREFIX = PREFIX)
                    #m[row.idx,1] <- .makeNames(x$lhs[row.idx], x$label[row.idx])
                } else if(s == "Thresholds") {
                    row.idx <- which(x$op == "|" & x$block == b)
                    if(length(row.idx) == 0L) next
                    m[row.idx,1] <- .makeNames(paste(x$lhs[row.idx], "|",
                        x$rhs[row.idx], sep=""), x$label[row.idx])
                } else if(s == "Variances") {
                    row.idx <- which(x$op == "~~" & x$lhs == x$rhs & !x$exo &
                                     x$block == b)
                    if(length(row.idx) == 0L) next
                    # make distinction between residual and plain
                    y.names <- unique( c(lavNames(x, "eqs.y"),
                                         lavNames(x, "ov.ind"),
                                         lavNames(x, "lv.ind")) )
                    PREFIX <- rep("", length(row.idx))
                    PREFIX[ x$rhs[row.idx] %in% y.names ] <- "  ."
                    m[row.idx,1] <- .makeNames(x$rhs[row.idx], x$label[row.idx],
                                               PREFIX = PREFIX)
                } else if(s == "Scales y*") {
                    row.idx <- which(x$op == "~*~" & x$block == b)
                    if(length(row.idx) == 0L) next
                    m[row.idx,1] <- .makeNames(x$rhs[row.idx], x$label[row.idx])
                } else if(s == "Group Weight") {
                        row.idx <- which(x$lhs == "group" & x$op == "%" & x$block == b)
                    if(length(row.idx) == 0L) next
                    m[row.idx,1] <- .makeNames(x$rhs[row.idx], x$label[row.idx])
                } else if(s == "R-Square") {
                        row.idx <- which(x$op == "r2" & x$block == b)
                    if(length(row.idx) == 0L) next
                    m[row.idx,1] <- .makeNames(x$rhs[row.idx], x$label[row.idx])
                } else {
                        row.idx <- integer(0L)
                }

                # do we need special formatting for this section?
                # three types:
                #  - regular (nothing to do, except row/colnames)
                #  - R-square
                #  - Latent Variables (and Composites), Regressions and Covariances
                #    'bundle' the output per lhs element

                # bundling
                if(s %in% c("Latent Variables", "Composites",
                            "Regressions", "Covariances")) {
                    nel <- length(row.idx)
                    M <- matrix("", nrow = nel*2, ncol = ncol(m))
                    colnames(M) <- colnames(m)
                    rownames(M) <- rep("", NROW(M))
                    #colnames(M)[1] <- sprintf("%-17s", paste(s, ":", sep = ""))
                    if(is.null(x$efa)) {
                        LHS <- paste(x$lhs[row.idx], x$op[row.idx])
                    } else {
                        LHS <- paste(x$lhs[row.idx], x$op[row.idx],
                                     x$efa[row.idx])
                    }
                    lhs.idx <- seq(1, nel*2L, 2L)
                    rhs.idx <- seq(1, nel*2L, 2L) + 1L
                    if(s == "Covariances") {
                        # make distinction between residual and plain
                        y.names <- unique( c(lavNames(x, "eqs.y"),
                                             lavNames(x, "ov.ind"),
                                             lavNames(x, "lv.ind")) )
                            PREFIX <- rep("", length(row.idx))
                        PREFIX[ x$lhs[row.idx] %in% y.names ] <- "."
                    } else {
                        PREFIX <- rep("", length(LHS))
                    }
                    M[lhs.idx, 1] <- sprintf("%1s%-15s", PREFIX, LHS)
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
                            if(length(del.idx) > 0L) {
                            M <- M[-del.idx,,drop=FALSE]
                        }
                    }
                    cat("\n", s, ":\n", sep = "")
                        #cat("\n")
                    print(M, quote = FALSE)

                # R-square
                } else if(s == "R-Square") {
                    M <- m[row.idx,1:2,drop=FALSE]
                    colnames(M) <- colnames(m)[1:2]
                    rownames(M) <- rep("", NROW(M))
                    #colnames(M)[1] <- sprintf("%-17s", paste(s, ":", sep = ""))
                    cat("\n", s, ":\n", sep = "")
                    #cat("\n")
                    print(M, quote = FALSE)

                # Regular
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
            } # GSECTIONS

        } # levels

    } # groups

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

.makeNames <- function(NAMES, LABELS, PREFIX = NULL) {

    W <- 14
    if(is.null(PREFIX)) {
        PREFIX <- rep("", length(NAMES))
    }

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

    char.format <- paste("%3s%-", W, "s", sep = "")
    sprintf(char.format, PREFIX, NAMES)
}

.makeConNames <- function(lhs, op, rhs, nd) {

    nd <- max(nd, 3)
    W <- 41 + (nd - 3)*3

    nel <- length(lhs)
    if(length(nel) == 0L) return(character(0))
    NAMES <- character(nel)
    for(i in 1:nel) {
        if(rhs[i] == "0" && op[i] == ">") {
            con.string <- paste(lhs[i], " - 0", sep="")
        } else if(rhs[i] == "0" && op[i] == "<") {
            con.string <- paste(rhs[i], " - (", lhs[i], ")", sep="")
        } else if(rhs[i] != "0" && op[i] == ">") {
            con.string <- paste(lhs[i], " - (", rhs[i], ")", sep="")
        } else if(rhs[i] != "0" && op[i] == "<") {
            con.string <- paste(rhs[i], " - (", lhs[i], ")", sep="")
        } else if(rhs[i] == "0" && op[i] == "==") {
            con.string <- paste(lhs[i], " - 0", sep="")
        } else if(rhs[i] != "0" && op[i] == "==") {
            con.string <- paste(lhs[i], " - (", rhs[i], ")", sep="")
        }
        con.string <- abbreviate(con.string, W, strict = TRUE)
        char.format <- paste("   %-", W, "s", sep = "")
        NAMES[i] <- sprintf(char.format, con.string)
    }

    NAMES
}

summary.lavaan.fsr <- function(object, ...) {

    dotdotdot <- list(...)
    if(!is.null(dotdotdot$nd)) {
        nd <- dotdotdot$nd
    } else {
        nd <- 3L
    }

    print.lavaan.fsr(x = object, nd = nd, mm = TRUE, struc = TRUE)
}

print.lavaan.fsr <- function(x, ..., nd = 3L, mm = FALSE, struc = FALSE) {

    y <- unclass(x)

    # print header
    if(!is.null(y$header)) {
        cat(y$header)
        cat("\n")
    }

    if(mm && !is.null(y$MM.FIT)) {
        cat("\n")
        nblocks <- length(y$MM.FIT)
        for(b in seq_len(nblocks)) {
            cat("Measurement block for latent variable(s):",
                paste(lavNames(y$MM.FIT[[b]], "lv")), "\n")

            # fit measures?
            b.options <- lavInspect(y$MM.FIT[[b]], "options")
            if(!(length(b.options$test) == 1L && b.options$test == "none")) {
                cat("\n")
                print(fitMeasures(y$MM.FIT[[b]], c("chisq", "df", "pvalue", "cfi", "rmsea", "srmr")))
            }

            # parameter estimates
            PE <- parameterEstimates(y$MM.FIT[[b]], ci = FALSE,
                                     output = "text", header = TRUE)
            print.lavaan.parameterEstimates(PE, ..., nd = nd)
            cat("\n")
        }
    }

    # print PE
    if(struc) {
        cat("Structural Part\n")
        cat("\n")
        #print.lavaan.parameterEstimates(y$PE, ..., nd = nd)

        print(summary(y$STRUC.FIT, fit.measures = FALSE, estimates = FALSE,
                                    modindices = FALSE))
        FIT <- fitMeasures(y$STRUC.FIT, fit.measures="default")
        if(FIT["df"] > 0) {
            print.lavaan.fitMeasures( FIT, add.h0 = FALSE )
        }
    }
    PE <- parameterEstimates(y$STRUC.FIT, ci = FALSE,
                             remove.eq = FALSE, remove.system.eq = TRUE,
                             remove.ineq = FALSE, remove.def = FALSE,
                             remove.nonfree = FALSE,
                             output = "text", header = TRUE)
    print.lavaan.parameterEstimates(PE, ..., nd = nd)

    invisible(y)
}


# print warnings/errors in a consistent way
# YR 12 July 2018

lav_txt2message <- function(txt, header = "lavaan WARNING:",
                            footer = "", txt.width = 70L, shift = 3L) {
    # make sure we only have a single string
    txt <- paste(txt, collapse = "")

    # split the txt in little chunks
    chunks <- strsplit(txt, "\\s+", fixed = FALSE)[[1]]

    # chunk size (number of characters)
    chunk.size <- nchar(chunks)

    # remove empty chunks (needed?)
    empty.idx <- which(chunk.size == 0)
    if(length(empty.idx) > 0L) {
        chunks <- chunks[-empty.idx]
        chunk.size <- chunk.size[-empty.idx]
    }

    # insert "\n" so the txt never gets wider than txt.width
    num.lines <- floor((sum(chunk.size) + length(chunk.size))/txt.width + 0.5)
    target <- character(num.lines)

    line.size <- shift
    line.num <- 1L
    start.chunk <- 1L
    end.chunck <- 1L
    for(ch in seq_len( length(chunks) )) {
        line.size <- line.size + chunk.size[ch] + 1L
        if(line.size > txt.width) {
            end.chunk <- ch - 1L
            target[line.num] <- paste(c(rep(" ", (shift-1)),
                                        chunks[ start.chunk:end.chunk ]),
                                      collapse = " ")
            line.num <- line.num + 1L
            start.chunk <- ch
            line.size <- shift + chunk.size[ch] + 1L
        }
    }
    # last line
    target[line.num] <- paste(c(rep(" ", (shift-1)),
                              chunks[ start.chunk:ch ]), collapse = " ")

    body <- paste(target, collapse = "\n")
    if(nchar(footer) == 0L) {
       out <- paste(c(header, body), collapse = "\n")
    } else {
       out <- paste(c(header, body, footer), collapse = "\n")
    }

    out
}

# new in 0.6-12
print.lavaan.summary <- function(x, ..., nd = 3L) {

    y <- unclass(x) # change to ordinary list

    # get nd, if it is stored as an attribute
    ND <- attr(y, "nd")
    if(!is.null(ND) && is.numeric(ND)) {
        nd <- as.integer(ND)
    }

    # header
    if(!is.null(y$header)) {
        lavaan.version   <- y$header$lavaan.version
        sam.approach     <- y$header$sam.approach
        optim.method     <- y$header$optim.method
        optim.iterations <- y$header$optim.iterations
        optim.converged  <- y$header$optim.converged

        # sam or sem?
        if(sam.approach) {
            cat("This is ",
                sprintf("lavaan %s", lavaan.version),
                " -- using the SAM approach to SEM\n", sep = "")
        } else {
            cat(sprintf("lavaan %s ", lavaan.version))

            # Convergence or not?
            if(optim.method == "none") {
                cat("-- DRY RUN with 0 iterations --\n")
            } else if(optim.iterations > 0) {
                if(optim.converged) {
                    if(optim.iterations == 1L) {
                        cat("ended normally after 1 iteration\n")
                    } else {
                        cat(sprintf("ended normally after %i iterations\n",
                            optim.iterations))
                    }
                } else {
                    if(optim.iterations == 1L) {
                        cat("did NOT end normally after 1 iteration\n")
                    } else {
                      cat(sprintf("did NOT end normally after %i iterations\n",
                            optim.iterations))
                    }
                    cat("** WARNING ** Estimates below are most likely unreliable\n")
                }
            } else {
                cat("did not run (perhaps do.fit = FALSE)?\n")
                cat("** WARNING ** Estimates below are simply the starting values\n")
            }
        }
    }

    # optim
    if(!is.null(y$optim)) {
        estimator      <- y$optim$estimator
        estimator.args <- y$optim$estimator.args
        optim.method   <- y$optim$optim.method
        npar           <- y$optim$npar
        eq.constraints <- y$optim$eq.constraints
        nrow.ceq.jac   <- y$optim$nrow.ceq.jac
        nrow.cin.jac   <- y$optim$nrow.cin.jac
        nrow.con.jac   <- y$optim$nrow.con.jac
        con.jac.rank   <- y$optim$con.jac.rank

        cat("\n")
        # cat("Optimization information:\n\n")

        c1 <- c("Estimator")
        # second column
        tmp.est <- toupper(estimator)
        if(tmp.est == "DLS") {
            dls.first.letter <- substr(estimator.args$dls.GammaNT,
                                       1L, 1L)
            tmp.est <- paste("DLS-", toupper(dls.first.letter), sep = "")
        }
        c2 <- tmp.est

        # additional estimator args
        if(!is.null(estimator.args) &&
           length(estimator.args) > 0L) {
            if(estimator == "DLS") {
                c1 <- c(c1, "Estimator DLS value for a")
                c2 <- c(c2, estimator.args$dls.a)
            }
        }

        # optimization method + npar
        c1 <- c(c1, "Optimization method", "Number of model parameters")
        c2 <- c(c2, toupper(optim.method), npar)

        # optional output
        if(eq.constraints) {
            c1 <- c(c1, "Number of equality constraints")
            c2 <- c(c2, nrow.ceq.jac)
        }
        if(nrow.cin.jac > 0L) {
            c1 <- c(c1, "Number of inequality constraints")
            c2 <- c(c2, nrow.cin.jac)
        }
        if(nrow.con.jac > 0L) {
            if(con.jac.rank == (nrow.ceq.jac + nrow.cin.jac)) {
            # nothing to do (don't print, as this is redundant information)
            } else {
                c1 <- c(c1, "Row rank of the constraints matrix")
                c2 <- c(c2, con.jac.rank)
            }
        }

        # format
        c1 <- format(c1, width = 40L)
        c2 <- format(c2, width = 11L + max(0, (nd - 3L)) * 4L,
                     justify = "right")

        # character matrix
        M <- cbind(c1, c2, deparse.level = 0)
        colnames(M) <- rep("",  ncol(M))
        rownames(M) <- rep(" ", nrow(M))

        # print
        write.table(M, row.names = TRUE, col.names = FALSE, quote = FALSE)
    }

    # sam header
    if(!is.null(y$sam.header)) {
        cat("\n")
        sam.method          <- y$sam.header$sam.method
        sam.local.options   <- y$sam.header$sam.local.options
        sam.mm.list         <- y$sam.header$sam.mm.list
        sam.mm.estimator    <- y$sam.header$sam.mm.estimator
        sam.struc.estimator <- y$sam.header$sam.struc.estimator

        # sam method
        c1 <- c("SAM method")
        c2 <- toupper(sam.method)

        # options
        if(sam.method == "local") {
            c1 <- c(c1, "Mapping matrix M method")
            c2 <- c(c2, sam.local.options$M.method)
            # TODo: more!
        }

        # number of measurement blocks
        c1 <- c(c1, "Number of measurement blocks")
        c2 <- c(c2, length(sam.mm.list))

        # estimator measurement blocks
        c1 <- c(c1, "Estimator measurement part")
        c2 <- c(c2, sam.mm.estimator)

        # estimator structural part
        c1 <- c(c1, "Estimator  structural part")
        c2 <- c(c2, sam.struc.estimator)

        # format
        c1 <- format(c1, width = 40L)
        c2 <- format(c2, width = 11L + max(0, (nd - 3L)) * 4L,
                     justify = "right")

        # character matrix
        M <- cbind(c1, c2, deparse.level = 0)
        colnames(M) <- rep("",  ncol(M))
        rownames(M) <- rep(" ", nrow(M))

        # print
        write.table(M, row.names = TRUE, col.names = FALSE, quote = FALSE)
    }

    # efa/rotation
    if(!is.null(y$rotation)) {
        cat("\n")
        rotation      <- y$rotation
        rotation.args <- y$rotation.args

        #cat("Rotation information:\n\n")
        # container
        c1 <- c2 <- character(0L)

        # rotation method
        c1 <- c(c1, "Rotation method")
        if(rotation$rotation == "none") {
            MM <- toupper(rotation$rotation)
        } else if(rotation$rotation.args$orthogonal) {
            MM <- paste(toupper(rotation$rotation), " ", "ORTHOGONAL",
                    sep = "")
        } else {
            MM <- paste(toupper(rotation$rotation), " ", "OBLIQUE",
                        sep = "")
        }
        c2 <- c(c2, MM)

        if(rotation$rotation != "none") {

            # method options
            if(rotation$rotation == "geomin") {
                c1 <- c(c1, "Geomin epsilon")
                c2 <- c(c2, rotation$rotation.args$geomin.epsilon)
            } else if(rotation$rotation == "orthomax") {
                c1 <- c(c1, "Orthomax gamma")
                c2 <- c(c2, rotation$rotation.args$orthomax.gamma)
            } else if(rotation$rotation == "cf") {
                c1 <- c(c1, "Crawford-Ferguson gamma")
                c2 <- c(c2, rotation$rotation.args$cf.gamma)
            } else if(rotation$rotation == "oblimin") {
                c1 <- c(c1, "Oblimin gamma")
                c2 <- c(c2, rotation$rotation.args$oblimin.gamma)
            }

            # rotation algorithm
            c1 <- c(c1, "Rotation algorithm (rstarts)")
            tmp <- paste(toupper(rotation$rotation.args$algorithm),
                         " (", rotation$rotation.args$rstarts, ")", sep = "")
            c2 <- c(c2, tmp)

            # Standardized metric (or not)
            c1 <- c(c1, "Standardized metric")
            if(rotation$rotation.args$std.ov) {
                c2 <- c(c2, "TRUE")
            } else {
                c2 <- c(c2, "FALSE")
            }

            # Row weights
            c1 <- c(c1, "Row weights")
            tmp.txt <- rotation$rotation.args$row.weights
            c2 <- c(c2, paste(toupper(substring(tmp.txt, 1, 1)),
                              substring(tmp.txt, 2), sep = ""))
        }

        # format c1/c2
        c1 <- format(c1, width = 33L)
        c2 <- format(c2, width = 18L + max(0, (nd - 3L)) * 4L,
                     justify = "right")

        # create character matrix
        M <- cbind(c1, c2, deparse.level = 0)
        colnames(M) <- rep("",  ncol(M))
        rownames(M) <- rep(" ", nrow(M))

        # print
        write.table(M, row.names = TRUE, col.names = FALSE, quote = FALSE)
    }

    # data object
    if(!is.null(y$data)) {
        cat("\n")
        lav_data_print_short(y$data, nd = nd)
    }

    # sam local stats: measurement blocks + structural part
    if(!is.null(y$sam)) {
        cat("\n")
        sam.method    <- y$sam$sam.method
        sam.mm.table  <- y$sam$sam.mm.table
        sam.mm.rel    <- y$sam$sam.mm.rel
        sam.struc.fit <- y$sam$sam.struc.fit
        ngroups       <- y$sam$ngroups
        group.label   <- y$sam$group.label

        # measurement
        tmp <- sam.mm.table
        if(sam.method == "global") {
            cat("Summary Information Measurement Part:\n\n")
        } else {
            cat("Summary Information Measurement + Structural:\n\n")
        }
        print(tmp, row.names = rep(" ", nrow(tmp)), nd = nd)

        if(sam.method == "local") {
            # reliability information
            c1 <- c2 <- character(0L)
            if(ngroups == 1L) {
                cat("\n")
                cat("  Model-based reliability latent variables:\n\n")
                tmp <- data.frame(as.list(sam.mm.rel[[1]]))
                class(tmp) <- c("lavaan.data.frame", "data.frame")
                print(tmp, row.names = rep(" ", nrow(tmp)), nd = nd)
            } else {
                cat("\n")
                cat("  Model-based reliability latent variables (per group):\n")
                for(g in 1:ngroups) {
                    cat("\n")
                    cat("  Group ", g, " [", group.label[g], "]:\n\n",
                        sep = "")
                    tmp <- data.frame(as.list(sam.mm.rel[[g]]))
                    class(tmp) <- c("lavaan.data.frame", "data.frame")
                    print(tmp, row.names = rep(" ", nrow(tmp)), nd = nd)
                }
            }

            cat("\n")
            cat("  Summary Information Structural part:\n\n")
            tmp <- data.frame(as.list(sam.struc.fit))
            class(tmp) <- c("lavaan.data.frame", "data.frame")
            print(tmp, row.names = rep(" ", nrow(tmp)), nd = nd)
        }
    }

    # test statistics
    if(!is.null(y$test)) {
        cat("\n")
        lav_test_print(y$test, nd = nd)
    }

    # extra fit measures (if present)
    if(!is.null(y$fit)) {
        print.lavaan.fitMeasures(y$fit, nd = nd, add.h0 = FALSE )
    }

    # parameter table
    if(!is.null(y$pe)) {
        PE <- y$pe
        class(PE) <- c("lavaan.parameterEstimates", "lavaan.data.frame",
                       "data.frame")
        print(PE, nd = nd)
    }

    # modification indices
    if(!is.null(y$mi)) {
        cat("Modification Indices:\n\n")
        MI <- y$mi
        rownames(MI) <- NULL
        print(MI, nd = nd)
    }

    invisible(y)
}
