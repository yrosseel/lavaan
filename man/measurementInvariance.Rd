\name{measurementInvariance}
\alias{measurementInvariance}
\alias{measurementinvariance}
\title{
Measurement Invariance Tests
}
\description{
Testing measurement invariance across groups using a typical sequence of
model comparison tests. 
}
\usage{
measurementInvariance(..., strict = FALSE, quiet = FALSE)
}
\arguments{
  \item{...}{The same arguments as for any lavaan model. 
    See \code{\link{cfa}} for more information.}
  \item{strict}{If \code{TRUE}, the sequence requires `strict' invariance.
    See details for more information.}
  \item{quiet}{If \code{FALSE}, a summary is printed out containing an
    overview of the different models that are fitted, together with some
    model comparison tests.}
}
\details{
If \code{strict = FALSE}, the following four models are tested in order:
\enumerate{
  \item{Model 1: configural invariance. The same factor structure is imposed
        on all groups.}
  \item{Model 2: weak invariance. The factor loadings are constrained to be
        equal across groups.}
  \item{Model 3: strong invariance. The factor loadings and intercepts are
        constrained to be equal across groups.}
  \item{Model 4: The factor loadings, intercepts and means are constrained to
        be equal across groups.}
}
Each time a more restricted model is fitted, a chi-square difference test
is reported, comparing the current model with the previous one, and comparing
the current model to the baseline model (Model 1). In addition, the difference
in cfi is also reported (delta.cfi).

If \code{strict = TRUE}, the following five models are tested in order:
\enumerate{
  \item{Model 1: configural invariance. The same factor structure is imposed
        on all groups.}
  \item{Model 2: weak invariance. The factor loadings are constrained to be
        equal across groups.}
  \item{Model 3: strong invariance. The factor loadings and intercepts are
        constrained to be equal across groups.}
  \item{Model 4: strict invariance. The factor loadings, intercepts and 
        residual variances are constrained to be equal across groups.}
  \item{Model 5: The factor loadings, intercepts, residual variances and means 
        are constrained to be equal across groups.}
}

Note that if the chi-square test statistic is scaled (eg. a
Satorra-Bentler or Yuan-Bentler test statistic), a special version of the 
chi-square difference test is used as described in 
\url{http://www.statmodel.com/chidiff.shtml}
}

\value{
Invisibly, all model fits in the sequence are returned as a list.
}
\references{
Vandenberg, R. J., and Lance, C. E. (2000). A review and synthesis of the 
measurement invariance literature: Suggestions, practices, and recommendations 
for organizational research. Organizational Research Methods, 3, 4-70.
}
\seealso{
\code{\link{cfa}}
}
\examples{
HW.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

measurementInvariance(HW.model, data=HolzingerSwineford1939, group="school")
}
