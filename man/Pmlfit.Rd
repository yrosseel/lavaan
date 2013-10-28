\name{lav_tables_fit}
\alias{lav_tables_fit_Cp}
\alias{lav_tables_fit_CF}
\alias{lav_tables_fit_CM}
\title{Pairwise maximum likelihood fit statistics}

\description{
Three measures of fit for the pairwise maximum likelihood estimation method that are based on likelihood ratios (LR) are defined: 
\eqn{C_F}, \eqn{C_M}, and \eqn{C_P}. Subscripts \eqn{F}, \eqn{M}, and \eqn{P} signify comparisons of model-implied proportions of response 
patterns with observed proportions of full response patterns (\eqn{C_F}), expected proportions under just the assumption of multivariate 
normality (\eqn{C_M}), and observed proportions of pairs of item responses (\eqn{C_P}).
}

\usage{
lav_tables_fit_CF(object, est = "h0")
lav_tables_fit_Cp(object, alpha = 0.05)
lav_tables_fit_CM(object)
}

\arguments{
\item{object}{An object of class \code{\linkS4class{lavaan}}.}
\item{alpha}{The nominal signifiance.}
\item{est}{If \code{"h0"}, the user-specified model is used; if \code{"h1"},
the unrestricted model is used.}
}

\references{
Barendse, M. T., Ligtvoet, R., Timmerman, M. E., & Oort, F. J. (under review). Structural Equation Modeling of Discrete data: 
Model Fit after Pairwise Maximum Likelihood.

Joreskog, K. G., & Moustaki, I. (2001). Factor analysis of ordinal variables: A comparison of three approaches. 
\emph{Mulitvariate Behavioral Research, 36}, 347-387. 
}

\details{
\subsection{\eqn{C_f}}{
\eqn{C_f} compares the log-likelihood of the model-implied proportions of multivariate response patterns(\eqn{/pi_r})with respect to the 
correlations and thresholds given in the PML estimation method with the observed sample proportions of response patterns (\eqn{p_r}). 
Multiplied by two times the sample size, we obtain
\deqn{
C_F = 2N\sum_{r}p_{r}\ln[p_{r}/\hat{\pi}_{r}],
}
which asymptotically has a chi-square distribution with
\deqn{
df_F = m^k - n - 1,
}
where \eqn{k} denote items with discrete response scales, \eqn{m} denote the response options, and \eqn{n} denote the number of parameters 
to be estimated. 
}

\subsection{\eqn{C_m}}{
The \eqn{C_M} is based on the \eqn{C_f} it, and compares the log-likelihood of the model-implied proportions of the model-of-interest with 
the model-implied proportions of the model that only assumes an underlying multivariate normal distribution, without any further 
restrictions:
\deqn{
C_M = C_{F1} - C_{F0},
}
where \eqn{C_{F1}} is \eqn{C_F} for Model 1, the model of interest, and \eqn{C_{F0}} is \eqn{C_F} for Model 0, the model that assumes 
underlying multivariate normality and that has all polychoric correlations and all thresholds as its parameters. Statistic \eqn{C_M} has a 
chi-square distribution with 
\deqn{
df_M = k^{2}/2 - k/2 + k(m-1) - n_{1},
}
where \eqn{k} denote items with discrete response scales, \eqn{m} denote the response options, \eqn{k^{2}/2 - k/2} denote the number of 
polychoric correlations, \eqn{k(m-1)} denote the number of thresholds, and \eqn{n_1} is the number of parameters of the model of interest. 
}

\subsection{\eqn{C_p}}{
The \eqn{C_P} only consider pairs of responses, and compare observed (\eqn{p}) and model-implied proportions of pairs of responses(\eqn{/pi}). 
For items \eqn{i} and \eqn{j} 
\deqn{
C_P=2N\sum_{c_i=1}^m \sum_{c_j=1}^m
p_{c_i,c_j}\ln[p_{c_i,c_j}/\hat{\pi}_{c_i,c_j}],
}
which has an asymptotic chi-square distribution with degrees of freedom equal to the information \eqn{(m^2 -1)} minus the number of 
parameters (2(m-1) thresholds and 1 correlation),
\deqn{
df_P = m^{2} - 2(m - 1) - 2.
}
As there are \eqn{k(k-1)/2} possible pairs of items, this \eqn{C_P} should be applied with a Bonferroni adjusted level of significance 
\eqn{\alpha^*}, with
\deqn{
\alpha^*= \alpha / k(k-1)/2),
}
to keep the family-wise error rate at \eqn{\alpha}. The hypothesis of overall goodness-of-fit is tested at \eqn{\alpha} and rejected as 
soon as \eqn{C_P} is significant at \eqn{\alpha^*} for at least one pair of items. Notice that with dichotomous items, \eqn{m = 2}, 
and \eqn{df_P = 0}, 
so that hypothesis can not be tested. 
}
}

\seealso{
\code{\link{lavTables}, \link{lavaan}}
}

\examples{
HS9 <- HolzingerSwineford1939[,c("x1","x2","x3","x4","x5",
                                 "x6","x7","x8","x9")]
HSbinary <- as.data.frame( lapply(HS9, cut, 2, labels=FALSE) )
HSbinary$school <- HolzingerSwineford1939$school
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- cfa(HS.model, data=HSbinary, ordered=names(HSbinary),group = "school", estimator="PML")
#Fit statistics for a multigroup model
lav_tables_fit_CM(fit)
lav_tables_fit_Cp(fit)
lav_tables_fit_CF(fit)

#HS.model <- ' visual  =~ x1 + x2 + x3
#              textual =~ x4 + x5 + x6
#              speed   =~ x7 + x8 + x9 
#              new     =~ x7 + x8'
#fit2 <- cfa(HS.model, data=HSbinary, ordered=names(HSbinary),group = "school", estimator="PML")
#Fit statistics
#lav_tables_fit_Cm(fit)
#lav_tables_fit_Cp(fit)
#lav_tables_fit_Cf(fit)

#model comparison
#anova(lav_tables_fit_Cf(fit),lav_tables_fit_Cf(fit2))
}

\keyword{
pairwise maximum likelihood analysis, discrete data, goodness of fit, lavaan
}

