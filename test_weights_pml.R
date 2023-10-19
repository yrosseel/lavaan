library(tidyverse)
devtools::load_all(".")

# Generate binary data
set.seed(15123)
pop.model <- '
f1 =~ 0.8 * y1 + 0.7 * y2 + 0.47 * y3 + 0.38 * y4 + 0.34 * y5

y1 | -1.43 * t1
y2 | -0.55 * t1
y3 | -0.13 * t1
y4 | -0.72 * t1
y5 | -1.13 * t1
'
Data <-
  simulateData(pop.model, sample.nobs = 1000) %>%
  as_tibble() %>%
  mutate(across(everything(), ordered)) %>%
  mutate(score = apply(., 1, function(x) sum(as.numeric(x))),
         wt = score + rchisq(1000, 1),
         wt = wt / sum(wt) * 1000)

# Fit model
mod <- '
f1 =~ y1 + y2 + y3 + y4 + y5
'
fit1 <- sem(model = mod, data = Data, estimator = "PML", std.lv = TRUE)
fit2 <- sem(model = mod, data = Data, estimator = "PML", std.lv = TRUE,
            sampling.weights = "wt")
# > coef(fit2)
# f1=~y1 f1=~y2 f1=~y3 f1=~y4 f1=~y5  y1|t1  y2|t1  y3|t1  y4|t1  y5|t1
# 0.713  0.740  0.499  0.327  0.415 -1.515 -0.680 -0.214 -0.810 -1.222
# > parTable(fit2)$se[1:10] |> round(3)
# [1] 0.075 0.070 0.062 0.065 0.070 0.062 0.043 0.040 0.045 0.053

# Check derivatives. Delta below gives the derivatives of standardised
# thresholds and polychoric correlations (rows of Delta in this order) with
# respect to model parameter vector theta. In lavaan in a factor analysis model
# the order of the individual parameters is: loadings, unstandardized
# thresholds, factor correlations
Delta1 <- lavaan:::computeDelta(lavmodel = fit1@Model)[[1]]
Delta2 <- lavaan:::computeDelta(lavmodel = fit2@Model)[[1]]
# > Delta2
# [,1]      [,2]      [,3] [,4] [,5] [,6]
# [1,] 0.0000000 0.0000000 0.0000000    1    0    0
# [2,] 0.0000000 0.0000000 0.0000000    0    1    0
# [3,] 0.0000000 0.0000000 0.0000000    0    0    1
# [4,] 0.7724621 0.6878226 0.0000000    0    0    0
# [5,] 0.4310165 0.0000000 0.6878226    0    0    0
# [6,] 0.0000000 0.4310165 0.7724621    0    0    0
