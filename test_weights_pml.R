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
  mutate(wt = 1)  # abs(rnorm(1000, mean = 1, sd = 2))
                  # wt = wt / sum(wt) * 1000\

# Fit model
mod <- '
f1 =~ y1 + y2 + y3 + y4 + y5
'
fit1 <- sem(model = mod, data = Data, estimator = "PML", std.lv = TRUE)
fit2 <- sem(model = mod, data = Data, estimator = "PML", std.lv = TRUE,
            sampling.weights = "wt")








