# utility functions
#
# initial version: YR 25/03/2009


# invert positive definite symmetric matrix (eg cov matrix)
# using choleski decomposition
# return log determinant as an attribute
inv.chol <- function(S, logdet=FALSE) {
    cS <- chol(S)
    #if( inherits(cS, "try-error") ) {
    #    print(S)
    #    warning("lavaan WARNING: symmetric matrix is not positive symmetric!")
    #}
    S.inv <- chol2inv( cS )
    if(logdet) {
        attr(S.inv, "logdet") <- sum(log(diag(cS)^2))
    }
    S.inv
}

# convert correlation matrix + standard deviations to covariance matrix
# based on cov2cor in package:stats
cor2cov <- function(R, sds, names=NULL) {

    p <- (d <- dim(R))[1L]
    if(!is.numeric(R) || length(d) != 2L || p != d[2L]) 
        stop("'V' is not a square numeric matrix")

    if(any(!is.finite(sds))) 
        warning("sds had 0 or NA entries; non-finite result is doubtful")

    #if(sum(diag(R)) != p) 
    #    stop("The diagonal of a correlation matrix should be all ones.")

    if(p != length(sds)) 
        stop("The standard deviation vector and correlation matrix have a different number of variables")

    S <- R
    S[] <- sds * R * rep(sds, each=p)

    # optionally, add names
    if(!is.null(names)) {
        stopifnot(length(names) == p)
        rownames(S) <- colnames(S) <- names
    }

    S
} 

# convert characters within single quotes to numeric vector
# eg. s <- '3 4.3 8e-3 2.0'
#     x <- char2num(s)
char2num <- function(s = '') {
    # first, strip all ',' or ';'
    s. <- gsub(","," ", s); s. <- gsub(";"," ", s.)
    tc <- textConnection(s.)
    x <- scan(tc, quiet=TRUE)
    close(tc)
    x 
}

# create full matrix based on lower.tri or upper.tri elements; add names
# always ROW-WISE!!
getCov <- function(x, lower = TRUE, diagonal = TRUE, sds = NULL,
                   names = paste("V", 1:nvar, sep="")) {

    # check x and sds
    if(is.character(x)) x <- char2num(x)
    if(is.character(sds)) sds <- char2num(sds)
   
    nels <- length(x)
    if(lower) {
        COV <- lav_matrix_lower2full(x, diagonal = diagonal)
    } else {
        COV <- lav_matrix_upper2full(x, diagonal = diagonal)
    }
    nvar <- ncol(COV)

    # if diagonal is false, assume unit diagonal
    if(!diagonal) diag(COV) <- 1

    # check if we have a sds argument
    if(!is.null(sds)) {
        stopifnot(length(sds) == nvar)
        COV <- cor2cov(COV, sds)
    }

    # names
    stopifnot(length(names) == nvar)
    rownames(COV) <- colnames(COV) <- names

    COV
}


# translate row+col matrix indices to vec idx
rowcol2vec <- function(row.idx, col.idx, nrow, symmetric=FALSE) {

    idx <- row.idx + (col.idx-1)*nrow
    if(symmetric) {
        idx2 <- col.idx + (row.idx-1)*nrow
        idx <- unique(sort(c(idx, idx2)))
    }
    idx
}

# dummy function to 'pretty' print a vector with fixed width
pprint.vector <- function(x, 
                          digits.after.period=3,
                          ncols=NULL, max.col.width=11,
                          newline=TRUE) {

    n <- length(x)
    var.names <- names(x)

    total.width = getOption("width")

    max.width <- max(nchar(var.names))
    if( max.width < max.col.width) { # shrink
        max.col.width <- max( max.width, digits.after.period+2)
    }

    # automatic number of columns
    if(is.null(ncols)) {
        ncols <- floor( (total.width-2) / (max.col.width+2) )
    }
    nrows <- ceiling(n / ncols)

    if(digits.after.period >= (max.col.width-3)) {
        max.col.width <- digits.after.period + 3
    }
    string.format <- paste(" %", max.col.width, "s", sep="")
    number.format <- paste(" %", max.col.width, ".", digits.after.period, "f", sep="")

    for(nr in 1:nrows) {
        rest <- min(ncols, n)
        if(newline) cat("\n")
        # labels
        for(nc in 1:rest) {
            vname <- substr(var.names[(nr-1)*ncols + nc], 1, max.col.width)
            cat(sprintf(string.format, vname))
        }        
        cat("\n")
        for(nc in 1:rest) {
            cat(sprintf(number.format, x[(nr-1)*ncols + nc]))
        }
        cat("\n")
        n <- n - ncols
    }
    if(newline) cat("\n")
}

# print only lower half of symmetric matrix
pprint.matrix.symm <- function(x, 
                               digits.after.period=3,
                               ncols=NULL, max.col.width=11,
                               newline=TRUE) {

    n <- ncol <- ncol(x); nrow <- nrow(x)
    stopifnot(ncol == nrow)
    var.names <- rownames(x)

    total.width = getOption("width")

    max.width <- max(nchar(var.names)) 
    if( max.width < max.col.width) { # shrink
        max.col.width <- max( max.width, digits.after.period+2)
    }

    # automatic number of columns
    if(is.null(ncols)) {
        ncols <- floor( (total.width-2) / (max.col.width+2) )
    }

    nblocks <- ceiling(n / ncols)

    if(digits.after.period >= (max.col.width-3)) {
        max.col.width <- digits.after.period + 3
    }
    fc.format     <- paste(" %", min(max.width, max.col.width), "s", sep="")
    string.format <- paste(" %", max.col.width, "s", sep="")
    number.format <- paste(" %", max.col.width, ".", digits.after.period, "f", sep="")

    for(nb in 1:nblocks) {
        rest <- min(ncols, n)
        if(newline) cat("\n")
        # empty column
        cat(sprintf(fc.format, ""))
        # labels
        for(nc in 1:rest) {
            vname <- substr(var.names[(nb-1)*ncols + nc], 1, max.col.width)
            cat(sprintf(string.format, vname))
        }        
        cat("\n")
        row.start <- (nb-1)*ncols + 1
        for(nr in row.start:nrow) {
            # label
            vname <- substr(var.names[nr], 1, max.col.width)
            cat(sprintf(fc.format, vname))
            col.rest <- min(rest, (nr - row.start + 1))
            for(nc in 1:col.rest) {
                value <- x[nr, (nb-1)*ncols + nc]
                cat(sprintf(number.format, value))
            }
            cat("\n")
        }
        n <- n - ncols
    }
    if(newline) cat("\n")
}

# elimination of rows/cols symmetric matrix
eliminate.rowcols <- function(x, el.idx=integer(0)) {

    if(length(el.idx) == 0) { 
        return( x )
    }
    stopifnot(ncol(x) == nrow(x))
    stopifnot(min(el.idx) > 0 && max(el.idx) <= ncol(x))

    x[-el.idx, -el.idx]
}

# elimination of rows/cols pstar symmetric matrix
eliminate.pstar.idx <- function(nvar=1, el.idx=integer(0), 
                                meanstructure=FALSE, type="all") {

    if(length(el.idx) > 0) {
        stopifnot(min(el.idx) > 0 && max(el.idx) <= nvar)
    }
 
    XX <- utils::combn(1:(nvar+1),2)
    XX[2,] <- XX[2,] - 1
    
    if(type == "all") {
        idx <- !(apply(apply(XX, 2, function(x) {x %in% el.idx}), 2, all))
    } else {
        idx <- !(apply(apply(XX, 2, function(x) {x %in% el.idx}), 2, any))
    }

    if(meanstructure) {
        idx <- c(!(1:nvar %in% el.idx), idx)
        #idx <- c(rep(TRUE, nvar), idx)
        
    }

    idx
}


# construct 'augmented' covariance matrix
# based on the covariance matrix and the mean vector
augmented.covariance <- function(S., mean) {

    S <- as.matrix(S.)
    m <- as.matrix(mean)
    p <- ncol(S)

    if(nrow(m) != p) {
        stop("incompatible dimension of mean vector")
    }

    out <- matrix(0, ncol=(p+1), nrow=(p+1))

    out[1:p,1:p] <- S + m %*% t(m)
    out[p+1,1:p] <- t(m)
    out[1:p,p+1] <- m
    out[p+1,p+1] <- 1
    
    out
}


# R code for forming block diagonal matrix posted by Berten Gunter on R-help:
# Berton Gunter gunter.berton at gene.com Fri Sep 2 01:08:58 CEST 2005
# slight changes by YR
bdiag <- function(...) {
    mlist <- list(...)
    # handle case in which list of matrices is given
    if(length(mlist) == 1) mlist <- unlist(mlist, recursive=FALSE)
    # remove/ignore NULL blocks
    null.idx <- sapply(mlist, is.null); mlist <- mlist[!null.idx]
    csdim <- rbind(c(0,0), apply(sapply(mlist,dim), 1, cumsum ))
    ret <- array(0, dim=csdim[length(mlist)+1,])
    add1 <- matrix(rep(1:0,2), ncol=2)
    for(i in seq(along=mlist)) {
        indx <- apply(csdim[i:(i+1),]+add1,2, function(x) x[1]:x[2])
          ## non-square matrix
        if(is.null(dim(indx))) ret[indx[[1]], indx[[2]]] <- mlist[[i]]
            ## square matrix
        else ret[indx[,1],indx[,2]] <- mlist[[i]]
    }
    ret
}

# rbind list of matrices (same number of columns)
rbindlist <- function(...) {
     mlist <- list(...)
     if(length(mlist) == 1) mlist <- unlist(mlist, recursive=FALSE)
     nrows <- sapply(mlist, nrow)
     end.row <- cumsum(nrows)
     start.row <- c(1, end.row + 1)
     ret <- matrix(0, nrow=sum(sapply(mlist, nrow)), ncol=ncol(mlist[[1]]))
     for(i in seq(along=mlist)) {
         ret[start.row[i]:end.row[i],] <- mlist[[i]]    
     }   
     ret
}



# linesearch using 'armijo' backtracking
# to find a suitable `stepsize' (alpha)
linesearch.backtracking.armijo <- function(f.alpha, s.alpha, alpha=10) {

    tau <- 0.5
    ftol <- 0.001

    f.old <- f.alpha(0)
    s.old <- s.alpha(0)

    armijo.condition <- function(alpha) {
        f.new <- f.alpha(alpha)

        # condition
        f.new > f.old + ftol * alpha * s.old
    }

    i <- 1
    while(armijo.condition(alpha)) {
        alpha <- alpha * tau
        f.new <- f.alpha(alpha)
        cat("... backtracking: ", i, "alpha = ", alpha, "f.new = ", f.new, "\n")
        i <- i + 1
    }

    alpha
}

steepest.descent <- function(start, objective, gradient, iter.max, verbose) {

    x <- start

    if(verbose) {
        cat("Steepest descent iterations\n")
        cat("iter        function  abs.change  rel.change     step.size       norm.gx\n")
        gx <- gradient(x)
        norm.gx <- sqrt( gx %*% gx )
        fx <- objective(x)
        cat(sprintf("%4d   %11.7E                         %11.5E %11.5E",
                    0, fx, 0, norm.gx), "\n")
    }

    for(iter in 1:iter.max) {
        fx.old <- objective(x)

        # normalized gradient
        gx <- gradient(x)
        old.gx <- gx
        norm.gx <- sqrt( gx %*% gx )
        gradient.old <- gx / norm.gx
        direction.vector <- (-1) * gradient.old

        f.alpha <- function(alpha) {
            new.x <- x + alpha * direction.vector
            fx <- objective(new.x)
            #cat("  [stepsize]  iter ", iter, " step size = ", alpha,
            #    " fx = ", fx, "\n", sep="")
            # for optimize only
            if(is.infinite(fx)) {
                fx <- .Machine$double.xmax
            }
            fx
        }

        #s.alpha <- function(alpha) {
        #    new.x <- x + alpha * direction.vector
        #    gradient.new <- gradient(new.x)
        #    norm.gx <- sqrt( gradient.new %*% gradient.new)
        #    gradient.new <- gradient.new/norm.gx
        #    as.numeric(gradient.new %*% direction.vector)
        #}

        # find step size
        #alpha <- linesearch.backtracking.armijo(f.alpha, s.alpha, alpha=1)

        if(iter == 1) {
            alpha <- 0.1
        } else {
            alpha <- optimize(f.alpha, lower=0.0, upper=1)$minimum
            if( f.alpha(alpha) > fx.old ) {
                alpha <- optimize(f.alpha, lower=-1, upper=0.0)$minimum
            }
        }
        

        # steepest descent step
        old.x <- x
        x <- x + alpha * direction.vector
        gx.old <- gx
        gx <- gradient(x)
        dx.max <- max(abs( gx ))

        # verbose
        if(verbose) {
            fx <- fx.old
            fx.new <- objective(x)
            abs.change <- fx.new - fx.old
            rel.change <- abs.change / fx.old
            norm.gx <- sqrt(gx %*% gx)
            if(verbose) {
                cat(sprintf("%4d   %11.7E %10.7f %10.7f %11.5E %11.5E",
                            iter, fx.new, abs.change, rel.change, alpha, norm.gx), 
                    "\n")
            }
        }

        # convergence check
        if( dx.max < 1e-05 )
            break
    }

    x
}

