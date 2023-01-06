\name{summary.efaList}
\alias{summary.efaList}
\alias{efaList}
\alias{print.efaList.summary}
\title{Summarizing EFA Fits}
\description{
S3 summary and print methods for class \code{efaList}.}
\usage{
\method{summary}{efaList}(object,
        nd = 3L, cutoff = 0.3, dot.cutoff = 0.1, alpha.level = 0.01,
        lambda = TRUE, theta = TRUE, psi = TRUE, fit.table = TRUE,
        fs.determinacy = FALSE, eigenvalues = TRUE, sumsq.table = TRUE, 
        lambda.structure = FALSE, se = FALSE, zstat = FALSE,
        pvalue = FALSE, \dots)

\method{print}{efaList.summary}(x, nd = 3L, cutoff = 0.3, dot.cutoff = 0.1,
               alpha.level = 0.01, \dots)
}
\arguments{
\item{object}{An object of class \code{efaList}, usually, a result of a call
    to \code{\link{efa}} with (the default) \code{output = "efa"}.}
\item{x}{An object of class \code{summary.efaList}, usually, a result of a 
    call to \code{summary.efaList}.}
\item{nd}{Integer. The number of digits that are printed after the decimal 
    point in the output.}
\item{cutoff}{Numeric. Factor loadings smaller that this value (in absolute
    value) are not printed
    (even if they are significantly different from zero). The idea is that
    only medium to large factor loadings are printed, to better see the
    overall structure.}
\item{dot.cutoff}{Numeric. Factor loadings larger (in absolute value) 
    than this value, but smaller (in absolute value) than the cutoff value 
    are shown as a dot. They represent small loadings that may still need your
    attention.}
\item{alpha.level}{Numeric. If the the p-value of a factor loading is smaller
    than this value, a significance star is printed to the right of the 
    factor loading. To switch this off, use \code{alpha.level = 0}.}
\item{lambda}{Logical. If \code{TRUE}, include the (standardized) factor 
    loadings in the summary.}
\item{theta}{Logical. If \code{TRUE}, include the unique variances and
    the communalities in the table of factor loadings.}
\item{psi}{Logical. If \code{TRUE}, include the factor correlations in the 
    summary. Ignored if only a single factor is used.}
\item{fit.table}{Logical. If \code{TRUE}, show fit information for each
    model.}
\item{fs.determinacy}{Logical. If \code{TRUE}, show the factor score
    determinacy values per factor (assuming regression factor scores are used)
    and their squared values.}
\item{eigenvalues}{Logical. If \code{TRUE}, include the eigenvalues of the
    sample variance-covariance matrix in the summary.}
\item{sumsq.table}{Logical. If \code{TRUE}, include a table including sums
    of squares of factor loadings (and related measures) in the summary.
    The sums of squares are computed as the
    diagonal elements of Lambda times Psi (where Psi is the matrix of
    factor correlations.).
    If orthogonal rotation was used, Psi is diagonal and the sums of 
    squares are identical to the sums of the squared column elements of the 
    Lambda matrix (i.e., the factor loadings). This is no longer the case
    when obique rotation has been used. But in both cases (orthgonal or
    oblique), the (total) sum of the sums of squares equals the sum of the
    communalities. In the second row of the table (Proportion of total), 
    the sums of squares are divided by the total. In the third row of the 
    table (Proportion var), the sums of squares are divided by the 
    number of items.}
\item{lambda.structure}{Logical. If \code{TRUE}, show the structure matrix
    (i.e., the factor loadings multiplied by the factor correlations).}
\item{se}{Logical. If \code{TRUE}, include the standard errors of the 
    standardized lambda, theta and psi elements in the summary.}
\item{zstat}{Logical. If \code{TRUE}, include the Z-statistics of the 
    standardized lambda, theta and psi elements in the summary.}    
\item{pvalue}{Logical. If \code{TRUE}, include the P-values of the 
    standardized lambda, theta and psi elements in the summary.}
\item{...}{Further arguments passed to or from other methods.}
}
\value{
    The function \code{summary.efaList} computes and returns a list of 
    summary statistics for the list of EFA models in \code{object}.
}
\examples{
## The famous Holzinger and Swineford (1939) example
fit <- efa(data = HolzingerSwineford1939, 
           ov.names = paste("x", 1:9, sep = ""),
           nfactors = 1:3,
           rotation = "geomin",
           rotation.args = list(geomin.epsilon = 0.01, rstarts = 1))
summary(fit, nd = 3L, cutoff = 0.2, dot.cutoff = 0.05,
        lambda.structure = TRUE, pvalue = TRUE)
}
