\name{modificationIndices}
\alias{modificationIndices}
\alias{modificationindices}
\alias{modindices}
\title{Modification Indices}
\description{
Modification indices of a latent variable model.}
\usage{
modificationIndices(object, standardized = TRUE, power = FALSE, 
                    delta = 0.1, alpha = 0.05, high.power = 0.75)
modindices(object, standardized = TRUE, power = FALSE,
                    delta = 0.1, alpha = 0.05, high.power = 0.75)
}
\arguments{
\item{object}{An object of class \code{\linkS4class{lavaan}}.}
\item{standardized}{If \code{TRUE}, two extra columns (sepc.lv and sepc.all) 
will contain standardized values for the epc's. In the first column (sepc.lv),
standardizization is based on the variances of the (continuous) latent
variables. In the second column (sepc.all), standardization is based
on both the variances of both (continuous) observed and latent variables.
(Residual) covariances are standardized using (residual) variances.}
\item{power}{If \code{TRUE}, the (post-hoc) power is computed for each
modification index, using the values of \code{delta} and \code{alpha}.}
\item{delta}{The value of the effect size, as used in the post-hoc power 
computation, currently using the unstandardized metric of the epc column.}
\item{alpha}{The significance level used for deciding if the modification
index is statistically significant or not.}
\item{high.power}{If the computed power is higher than this cutoff value,
the power is considered `high'. If not, the power is considered `low'. 
This affects the values in the 'decision' column in the output.}
}
\value{
  A data.frame containing modification indices and EPC's.
}
\examples{
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- cfa(HS.model, data=HolzingerSwineford1939)
modindices(fit)
}
