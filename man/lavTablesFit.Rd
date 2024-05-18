\name{lavTablesFitCp}
\alias{lavTablesFit}
\alias{lavTablesFitCp}
\alias{lavTablesFitCf}
\alias{lavTablesFitCm}
\title{Pairwise maximum likelihood fit statistics}

\description{
Three measures of fit for the pairwise maximum likelihood estimation method that are based on likelihood ratios (LR) are defined: 
\eqn{C_F}, \eqn{C_M}, and \eqn{C_P}. Subscript \eqn{F} signifies a comparison of model-implied proportions of full response 
patterns with observed sample proportions, subscript \eqn{M} signifies a comparison of model-implied proportions of full response 
patterns with the proportions implied by the assumption of multivariate normality, and subscript \eqn{P} signifies
a comparison of model-implied proportions of pairs of item responses with the observed proportions of pairs of item responses. 
}

\usage{
lavTablesFitCf(object)
lavTablesFitCp(object, alpha = 0.05)
lavTablesFitCm(object)
}

\arguments{
\item{object}{An object of class \code{\linkS4class{lavaan}}.}
\item{alpha}{The nominal level of signifiance of global fit.}
}

\references{
Barendse, M. T., Ligtvoet, R., Timmerman, M. E., & Oort, F. J. (2016). Structural Equation Modeling of Discrete data: 
Model Fit after Pairwise Maximum Likelihood. \emph{Frontiers in psychology, 7}, 1-8. 

Joreskog, K. G., & Moustaki, I. (2001). Factor analysis of ordinal variables: A comparison of three approaches. 
\emph{Multivariate Behavioral Research, 36}, 347-387. 
}

\details{
\subsection{\eqn{C_F}}{
The \eqn{C_F} statistic compares the log-likelihood of the model-implied proportions (\eqn{\pi_r}) with the observed proportions (\eqn{p_r})
of the full multivariate responses patterns:
\deqn{
C_F = 2N\sum_{r}p_{r}\ln[p_{r}/\hat{\pi}_{r}],
}
which asymptotically has a chi-square distribution with
\deqn{
df_F = m^k - n - 1,
}
where \eqn{k} denotes the number of items with discrete response scales, \eqn{m} denotes the number of response options, and \eqn{n} denotes 
the number of parameters to be estimated. Notice that \eqn{C_F} results may be biased because of large numbers of empty cells in the multivariate 
contingency table.
}

\subsection{\eqn{C_M}}{
The \eqn{C_M} statistic is based on the \eqn{C_F} statistic, and compares the proportions implied by the model of interest (Model 1) 
with proportions implied by the assumption of an underlying multivariate normal distribution (Model 0):
\deqn{
C_M = C_{F1} - C_{F0},
}
where \eqn{C_{F0}} is \eqn{C_F} for Model 0 and \eqn{C_{F1}} is \eqn{C_F} for Model 1. Statistic \eqn{C_M} has a chi-square distribution with 
degrees of freedom
\deqn{
df_M = k(k-1)/2 + k(m-1) - n_{1},
}
where \eqn{k} denotes the number of items with discrete response scales, \eqn{m} denotes the number of response options, and \eqn{k(k-1)/2} 
denotes the number of polychoric correlations, \eqn{k(m-1)} denotes the number of thresholds, and \eqn{n_1} is the number of parameters of the 
model of interest. Notice that \eqn{C_M}  results may be biased because of large numbers of empty cells in the multivariate contingency table. However, 
bias may cancels out as both Model 1 and Model 0 contain the same pattern of empty responses.
}

\subsection{\eqn{C_P}}{
With the \eqn{C_P} statistic we only consider pairs of responses, and compare observed sample proportions (\eqn{p}) with model-implied proportions 
of pairs of responses(\eqn{\pi}). For items \eqn{i} and \eqn{j} we obtain a pairwise likelihood ratio test statistic \eqn{C_{P_{ij}}}
\deqn{
C_{P_{ij}}=2N\sum_{c_i=1}^m \sum_{c_j=1}^m
p_{c_i,c_j}\ln[p_{c_i,c_j}/\hat{\pi}_{c_i,c_j}],
}
where \eqn{m} denotes the number of response options and \eqn{N} denotes sample size. The \eqn{C_P} statistic has an asymptotic chi-square distribution 
with degrees of freedom equal to the information \eqn{(m^2 -1)} minus the number of parameters (2(m-1) thresholds and 1 correlation),
\deqn{
df_P = m^{2} - 2(m - 1) - 2.
}
As \eqn{k} denotes the number of items, there are \eqn{k(k-1)/2} possible pairs of items. The \eqn{C_P} statistic should therefore be applied with 
a Bonferroni adjusted level of significance \eqn{\alpha^*}, with
\deqn{
\alpha^*= \alpha /(k(k-1)/2)),
}
to keep the family-wise error rate at \eqn{\alpha}. The hypothesis of overall goodness-of-fit is tested at \eqn{\alpha} and rejected as 
soon as \eqn{C_P} is significant at \eqn{\alpha^*} for at least one pair of items. Notice that with dichotomous items, \eqn{m = 2}, 
and \eqn{df_P = 0}, so that hypothesis can not be tested. 
}
}

\seealso{
\code{\link{lavTables}, \link{lavaan}}
}

\examples{
# Data
HS9 <- HolzingerSwineford1939[,c("x1","x2","x3","x4","x5",
                                 "x6","x7","x8","x9")]
HSbinary <- as.data.frame( lapply(HS9, cut, 2, labels=FALSE) )

# Single group example with one latent factor
HS.model <- ' trait =~ x1 + x2 + x3 + x4 '
fit <- cfa(HS.model, data=HSbinary[,1:4], ordered=names(HSbinary[,1:4]),
           estimator="PML")
lavTablesFitCm(fit)
lavTablesFitCp(fit)
lavTablesFitCf(fit)
}


