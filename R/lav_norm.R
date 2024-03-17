# simple derivatives of the normal distribution

# dnorm
dnorm_dummy <- function(y, mu = 0, sigma2 = 1) {
  sigma <- sqrt(sigma2)
  1 / (sigma * sqrt(2 * pi)) * exp(-0.5 * ((y - mu) / sigma * (y - mu) / sigma))
}

# dnorm_dmu_x <- function(x, y, sigma2 = 1) {
#     dnorm_dummy(y = y, mu = x, sigma2 = sigma2)
# }
# numDeriv:::grad(func=dnorm_dmu_x, x=0.3, y=2.3, sigma2=16)

# partial derivative - mu
dnorm_dmu <- function(y, mu = 0, sigma2 = 1) {
  dy <- dnorm(x = y, mean = mu, sd = sqrt(sigma2))
  (y - mu) / sigma2 * dy
}

# dnorm_dsigma2_x <- function(x, y, mu = 0) {
#    dnorm_dummy(y = y, mu = mu, sigma2 = x)
# }
# numDeriv:::grad(func=dnorm_dsigma2_x, x=16, y=2.3, mu=0.3)

# partial derivative - sigma2
dnorm_dsigma2 <- function(y, mu = 0, sigma2 = 1) {
  dy <- dnorm(x = y, mean = mu, sd = sqrt(sigma2))
  (1 / (2 * sigma2 * sigma2) * (y - mu) * (y - mu) - 1 / (2 * sigma2)) * dy
}

# dnorm_dy_x <- function(x, mu = 0, sigma2 = 1) {
#    dnorm_dummy(y = x, mu = mu, sigma2 = sigma2)
# }
# numDeriv:::grad(func=dnorm_dy_x, x=2.3, mu=0.3, sigma2=16)


# partial derivative - y
dnorm_dy <- function(y, mu = 0, sigma2 = 1) {
  dy <- dnorm(x = y, mean = mu, sd = sqrt(sigma2))
  -(y - mu) / sigma2 * dy
}


#### d log dnorm ####
#
# d log dnorm() / d theta   = 1/dy d dnorm() / d theta
dlogdnorm <- function(y, mu = 0, sigma2 = 1) {
  sigma <- sqrt(sigma2)
  -log(sigma * sqrt(2 * pi)) + (-0.5 * ((y - mu) / sigma * (y - mu) / sigma))
}

# dlogdnorm_dmu_x <- function(x, y, sigma2 = 1) {
#    dlogdnorm(y = y, mu = x, sigma2 = sigma2)
# }
# numDeriv:::grad(func=dlogdnorm_dmu_x, x=0.3, y=2.3, sigma2=16)

# partial derivative - mu
dlogdnorm_dmu <- function(y, mu = 0, sigma2 = 1) {
  (y - mu) / sigma2
}
# dlogdnorm_dmu(y = 2.3, mu = 0.3, sigma2 = 16)

# dlogdnorm_dsigma2_x <- function(x, y, mu = 0) {
#    dlogdnorm(y = y, mu = mu, sigma2 = x)
# }
# numDeriv:::grad(func=dlogdnorm_dsigma2_x, x=16, y=2.3, mu=0.3)

# partial derivative - sigma2
dlogdnorm_dsigma2 <- function(y, mu = 0, sigma2 = 1) {
  1 / (2 * sigma2 * sigma2) * (y - mu) * (y - mu) - 1 / (2 * sigma2)
}
# dlogdnorm_dsigma2(y = 2.3, mu = 0.3, sigma2 = 16)

# dlogdnorm_dy_x <- function(x, mu = 0, sigma2 = 1) {
#    dlogdnorm(y = x, mu = mu, sigma2 = sigma2)
# }
# numDeriv:::grad(func=dlogdnorm_dy_x, x=2.3, mu=0.3, sigma2=16)


# partial derivative - y
dlogdnorm_dy <- function(y, mu = 0, sigma2 = 1) {
  -(y - mu) / sigma2
}
# dlogdnorm_dy(y = 2.3, mu = 0.3, sigma2 = 16)
