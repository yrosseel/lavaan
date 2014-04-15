# tools for the multivariate Bernoulli distribution
#
# see:
#
# Maydeu-Olivares & Joe (2005). Limited- and Full-Information Estimation and
# Goodness-of-Fit Testing in 2^n Contingency Tables: A Unified Framework.
# Journal of the American Statistical Association, 100, 1009--1020.

# YR. 15 April 2014 -- first version

# compute higher-order joint moments (Teugels 1991)
# PROP must be an array, with dim = rep(2L, nitems)
lav_tables_mvb_getPiDot <- function(PROP, order. = nitems) {

    # number of items/dimensions
    nitems <- length(dim(PROP))

    # compute 'pi dot' up to order = order.
    pidot <- unlist(
        lapply(1:order., function(Order) {
            IDX <- utils::combn(1:nitems, Order)
            tmp <- apply(IDX, 2L, function(idx)
                as.numeric(apply(PROP, idx, sum))[1L])
            tmp 
        })
    )

    pidot
}

# compute 'T' matrix, so that pidot = T %*% prop
lav_tables_mvb_getT <- function(nitems = 3L, order. = nitems, rbind. = FALSE) {

    # index matrix
    INDEX <- array(1:(2^nitems), dim = rep(2L, nitems))

    T.r <- lapply(1:order., function(Order) {
        IDX <- utils::combn(1:nitems, Order)
        TT <- matrix(0L, ncol(IDX), 2^nitems)
        TT <- do.call("rbind", 
            lapply(1:ncol(IDX), function(i) {
                TRue <- as.list(rep(TRUE, nitems)); TRue[ IDX[,i] ] <- 1L
                ARGS <- c(list(INDEX), TRue)
                T1 <- integer( 2^nitems )
                T1[ as.vector(do.call("[", ARGS)) ] <- 1L
                T1
            }))
        TT
    })

    if(rbind.) {
        T.r <- do.call("rbind", T.r)
    }

    T.r
}

# simple test function to check that  pidot = T %*% prop
lav_tables_mvb_test <- function(nitems = 3L, verbose = FALSE) {
    
    freq <- sample( 5:50, 2^nitems, replace=TRUE)
    prop <- freq/sum(freq)
    TABLE <- array(freq, dim=rep(2, nitems))
    PROP  <- array(prop, dim=rep(2, nitems))
    # note: freq is always as.numeric(TABLE)
    #       prop is always as.numeric(PROP)

    pidot <- lav_tables_mvb_getPiDot(PROP)
    T.r <- lav_tables_mvb_getT(nitems = nitems, order. = nitems, rbind. = TRUE)

    if(verbose) {
        out <- cbind(as.numeric(T.r %*% prop), pidot)
        colnames(out) <- c("T * prop", "pidot")
        print(out)
    }

    all.equal(pidot, as.numeric(T.r %*% prop))
}

