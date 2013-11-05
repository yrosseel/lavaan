\name{lav_tables_fit}
\alias{lav_tables_fit_Cp}
\alias{lav_tables_fit_Cf}
\alias{lav_tables_fit_Cm}
\title{Pairwise maximum likelihood fit statistics}

\description{
Three measures of fit for the pairwise maximum likelihood estimation method that are based on likelihood ratios (LR) are defined: 
\eqn{C_F}, \eqn{C_M}, and \eqn{C_P}. Subscript \eqn{F} signifies a comparison of model-implied proportions of full response 
patterns with observed sample proportions, subscript \eqn{M} signifies a comparison of model-implied proportions of full response 
patterns with the model-implied proportions under just the assumption of multivariate normality, and subscript \eqn{P} signifies
a comparison of model-implied proportions of pairs of item responses with the observed proportions of pairs of item responses. 
}

\usage{
lav_tables_fit_Cf(object, est = "h0")
lav_tables_fit_Cp(object, alpha = 0.05)
lav_tables_fit_Cm(object)
}

\arguments{
\item{object}{An object of class \code{\linkS4class{lavaan}}.}
\item{alpha}{The nominal level of signifiance of global fit.}
\item{est}{If \code{"h0"}, a user-specified model is used; if \code{"h1"},
the unrestricted model of multivariate normality is used.}
}

\references{
Barendse, M. T., Ligtvoet, R., Timmerman, M. E., & Oort, F. J. (under review). Structural Equation Modeling of Discrete data: 
Model Fit after Pairwise Maximum Likelihood.

Joreskog, K. G., & Moustaki, I. (2001). Factor analysis of ordinal variables: A comparison of three approaches. 
\emph{Mulitvariate Behavioral Research, 36}, 347-387. 
}

\details{
\subsection{\eqn{C_F}}{
The \eqn{C_F} compares the log-likelihood of the model-implied proportions of multivariate response patterns(\eqn{\pi_r}) calculated with respect 
to the estimated PML correlations and thresholds with the log-likelihood of the observed sample proportions of full response patterns (\eqn{p_r}). 
Multiplied by two times the sample size (\eqn{N}), we obtain a likelihood ratio test statistic
\deqn{
C_F = 2N\sum_{r}p_{r}\ln[p_{r}/\hat{\pi}_{r}],
}
which asymptotically has a chi-square distribution with
\deqn{
df_F = m^k - n - 1,
}
where \eqn{k} denotes the number of items with discrete response scales, \eqn{m} denotes the number of response options, and \eqn{n} denotes 
the number of parameters to be estimated. Notice that CF results may be biased because of large numbers of empty cells in the multivariate 
contingency table.
}

\subsection{\eqn{C_M}}{
The \eqn{C_M} is based on the \eqn{C_F} it, and compares the fit of the model-implied proportions of the model-of-interest (Model 1) with 
the model-implied proportions of a model that only assumes an underlying multivariate normal distribution, without any further 
restrictions (Model 0):
\deqn{
C_M = C_{F0} - C_{F1},
}
where \eqn{C_{F0}} is \eqn{C_F} for Model 0, the model of interest, and \eqn{C_{F1}} is \eqn{C_F} for Model 1, the model that assumes 
underlying multivariate normality. Statistic \eqn{C_M} has a chi-square distribution with degrees of freedom
\deqn{
df_M = k(k-1)/2 + k(m-1) - n_{1},
}
where \eqn{k} denotes the number of items with discrete response scales, \eqn{m} denotes the number of response options, and \eqn{k(k-1)/2} 
denotes the number of polychoric correlations, \eqn{k(m-1)} denote the number of thresholds, and \eqn{n_1} is the number of parameters of the 
model of interest. Notice that CM results may be biased because of large numbers of empty cells in the multivariate contingency table. However, 
bias may cancels out as both Model 1 and Model 0 contain the same pattern of empty responses.
}

\subsection{\eqn{C_P}}{
With the \eqn{C_P} we only consider pairs of responses, and compare observed sample proportions (\eqn{p}) with model-implied proportions 
of pairs of responses(\eqn{\pi}). For items \eqn{i} and \eqn{j} we obtain a pairwise likelihood ratio test statistic \eqn{C_{P_{ij}}}
\deqn{
C_{P_{ij}}=2N\sum_{c_i=1}^m \sum_{c_j=1}^m
p_{c_i,c_j}\ln[p_{c_i,c_j}/\hat{\pi}_{c_i,c_j}],
}
where \eqn{m} denotes the number of response options and \eqn{N} denotes sample size. The (\eqn{C_P}) has an asymptotic chi-square distribution 
with degrees of freedom equal to the information \eqn{(m^2 -1)} minus the number of parameters (2(m-1) thresholds and 1 correlation),
\deqn{
df_P = m^{2} - 2(m - 1) - 2.
}
As \eqn{k} denotes the number of items, there are \eqn{k(k-1)/2} possible pairs of items. The \eqn{C_P} should therefore be applied with a Bonferroni 
adjusted level of significance \eqn{\alpha^*}, with
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
#Data
HS9 <- HolzingerSwineford1939[,c("x1","x2","x3","x4","x5",
                                 "x6","x7","x8","x9")]
HSbinary <- as.data.frame( lapply(HS9, cut, 2, labels=FALSE) )
HSbinary$school <- HolzingerSwineford1939$school

#Single group example with one latent factor
HS.model1 <- ' trait =~ x1 + x2 + x3 + x4 + x5 + x6 '
fit <- cfa(HS.model1, data=HSbinary[,1:6], ordered=names(HSbinary))
lavaan:::lav_tables_fit_CM(fit)
lavaan:::lav_tables_fit_Cp(fit)
lavaan:::lav_tables_fit_CF(fit, est="h0")

#Single group example with a two latent factors
#HS.model2 <- ' visual  =~ x1 + x2 + x3
#              textual =~ x4 + x5 + x6'
#fit2 <- cfa(HS.model2, data=HSbinary[,1:6], ordered=names(HSbinary), estimator="PML")
#Fit statistics
#lav_tables_fit_Cm(fit2)
#lav_tables_fit_Cp(fit2)
#lav_tables_fit_Cf(fit2)

#model comparison between the two different models
#anova(lav_tables_fit_Cf(fit1),lav_tables_fit_Cf(fit2))

#Multigroup example
#HS.model <- ' visual  =~ x1 + x2 + x3
#              textual =~ x4 + x5 + x6
#              speed   =~ x7 + x8 + x9 '
#fit <- cfa(HS.model, data=HSbinary, ordered=names(HSbinary),group = "school", estimator="PML")
#lav_tables_fit_CM(fit)
#lav_tables_fit_Cp(fit)
#lav_tables_fit_CF(fit)
}

\keyword{
pairwise maximum likelihood analysis, discrete data, goodness of fit, lavaan
}

