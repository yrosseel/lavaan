# the multivariate normal distribution, unrestricted (h1), missing values

# 1) loglikelihood --> same as h0 but where Mu and Sigma are unrestricted
# 2) 3) 4) 5)      --> (idem)

# YR 26 March 2016: first version

# here, we estimate Mu and Sigma from Y with missing values, assuming normality
# this is a rewrite of the 'estimate.moments.EM' function in <= 0.5-22
lav_mvnorm_missing_h1_estimate_moments <- function(Y        = NULL,
                                                   Mp       = NULL,
                                                   verbose  = FALSE,
                                                   max.iter = 500L,
                                                   tol      = 1e-05) {

}
