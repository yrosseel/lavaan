### Terrence D. Jorgensen
### Last updated: 25 January 2017
### test lav_mvnorm_* functions

context("lav_mvnorm_*")


## complete data
varnames <- paste("x", 1:9, sep = "")
H9 <- HolzingerSwineford1939[ , varnames]
## impose missingness
H9miss <- H9
H9miss$x5 <- ifelse(H9miss$x1 <= quantile(H9miss$x1, .3), NA, H9miss$x5)
H9miss$x9 <- ifelse(H9miss$x4 <= quantile(H9miss$x4, .3), NA, H9miss$x9)
## fit model to complete and incomplete data
HS.model <- '
    visual  =~ x1 + x2 + x3
    textual =~ x4 + x5 + x6
    speed   =~ x7 + x8 + x9
'


#########################
## Test Complete cases ##
#########################

cfit <- cfa(HS.model, data = H9, meanstructure = TRUE)

## save summary statistics
cM <- lavInspect(cfit, "sampstat", add.class = FALSE)$mean # matches round(colMeans(H9), 3)
cS <- lavInspect(cfit, "sampstat", add.class = FALSE)$cov  # matches round(cov(H9)*300/301, 3)
## model-implied moments
cMu <- lavInspect(cfit, "mean.ov", add.class = FALSE)
cSigma <- lavInspect(cfit, "cov.ov", add.class = FALSE)


## sum casewise log-likelihoods under saturated model
cLL1 <- fitMeasures(cfit)[["unrestricted.logl"]]
cLL2 <- sum(mnormt::dmnorm(H9, mean = cM, varcov = cS, log = TRUE))
#cLL3 <- sum(mvtnorm::dmvnorm(H9, mean = cM, sigma = cS, log = TRUE))
## functions of actual interest
cLL4 <- sum(lav_mvnorm_dmvnorm(Y = as.matrix(H9), Mu = cM, Sigma = cS))
cLL5 <- lav_mvnorm_h1_loglik_data(as.matrix(H9), casewise = FALSE)
cLL6 <- sum(lav_mvnorm_h1_loglik_data(as.matrix(H9), casewise = TRUE))

test_that("6 saturated log-likelihoods match for complete data", {
  expect_equal(cLL1, cLL2)
#  expect_equal(cLL1, cLL3)
  expect_equal(cLL1, cLL4)
  expect_equal(cLL1, cLL5)
  expect_equal(cLL1, cLL6)
})
rm(cLL1, cLL2,
   #cLL3, 
   cLL4, cLL5, cLL6)


## sum casewise log-likelihoods under target model
cLL1 <- fitMeasures(cfit)[["logl"]]
cLL2 <- sum(mnormt::dmnorm(H9, mean = cMu, varcov = cSigma, log = TRUE))
#cLL3 <- sum(mvtnorm::dmvnorm(H9, mean = cMu, sigma = cSigma, log = TRUE))
cLL4 <- sum(lav_mvnorm_dmvnorm(Y = as.matrix(H9), Mu = cMu, Sigma = cSigma))
cLL5 <- lav_mvnorm_loglik_samplestats(sample.mean = cM, sample.cov = cS,
                                      sample.nobs = nobs(cfit),
                                      Mu = cMu, Sigma = cSigma)

test_that("5 target-model log-likelihoods match for complete data", {
  expect_equal(cLL1, cLL2)
#  expect_equal(cLL1, cLL3)
  expect_equal(cLL1, cLL4)
  expect_equal(cLL1, cLL5)
})
rm(cLL1, cLL2,
   #cLL3, 
   cLL4, cLL5)


##################
## Missing Data ##
##################

mfit <- cfa(HS.model, data = H9miss, meanstructure = TRUE, missing = "fiml")

## list per missind-data pattern
lavInspect(mfit, "coverage")
pattern <- lavInspect(mfit, "pattern")
H9logic <- !is.na(H9miss)

## indicators for which pattern each row belongs to   # lav_data_missing_patterns(H9miss)$case.idx
indPatterns <- sapply(1:4, function(pp) {
  apply(H9logic, 1, function(x) all(x == pattern[pp, ]))
})
all(rowSums(indPatterns) == 1) # check exactly 1 pattern per person

## lists of sample stats per pattern
# (mN <- colSums(indPatterns)) # N per pattern
# mM <- lapply(1:4, function(pp) {
#   colMeans(H9miss[indPatterns[,pp], varnames[ pattern[pp,] ] ])
# })
# mS <- lapply(1:4, function(pp) {
#   cov(H9miss[indPatterns[,pp], varnames[pattern[pp,]]]) * (mN[pp] - 1) / mN[pp]
# })
## lists of model-implied moments
mMu <- lavInspect(mfit, "mean.ov", add.class = FALSE)
mSigma <- lavInspect(mfit, "cov.ov", add.class = FALSE)




## sum casewise log-likelihoods under saturated model for each pattern
mLL1 <- fitMeasures(mfit)[["logl"]]
mLL2 <- sum(sapply(1:4, function(pp) {
  sum(apply(H9miss[indPatterns[,pp], varnames[pattern[pp,]]], 1,
            mnormt::dmnorm, mean = mMu[varnames[pattern[pp,]]],
            varcov = mSigma[varnames[pattern[pp,]], varnames[pattern[pp,]]], log = TRUE))
}))
#mLL3 <- sum(sapply(1:4, function(pp) {
#  sum(apply(H9miss[indPatterns[,pp], varnames[pattern[pp,]]], 1,
#            mvtnorm::dmvnorm, mean = mMu[varnames[pattern[pp,]]],
#            sigma = mSigma[varnames[pattern[pp,]], varnames[pattern[pp,]]], log = TRUE))
#}))
## functions of actual interest
mLL4 <- lav_mvnorm_missing_loglik_data(H9miss, Mu = mMu, Sigma = mSigma, pattern = FALSE)
mLL5 <- lav_mvnorm_missing_loglik_data(H9miss, Mu = mMu, Sigma = mSigma, pattern = TRUE)
## from sample stats
mLL6 <- lav_mvnorm_missing_loglik_samplestats(mfit@SampleStats@missing[[1]], Mu = mMu, Sigma = mSigma)

test_that("6 target-model log-likelihoods match for missing data", {
  expect_equal(mLL1, mLL2)
#  expect_equal(mLL1, mLL3)
  expect_equal(mLL1, mLL4)
  expect_equal(mLL1, mLL5)
  expect_equal(mLL1, mLL6)
})
rm(mLL1, mLL2,
   #mLL3, 
   mLL4, mLL5, mLL6)


#########################
## run tests in this file
# test_file("tests/testthat/test_lav_mvnorm.R")

