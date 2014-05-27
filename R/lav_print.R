## NOTE:
## round(1.2355, 3) = 1.236
## but 
## round(1.2345, 3) = 1.234
##
## perhaps we should add 0.0005 or something to avoid this?

print.lavaan.data.frame <- function(x, ..., nd=3) {

    y <- as.data.frame(lapply(x, function(x) {
                              if(is.numeric(x)) round(x, nd) else x}))

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

