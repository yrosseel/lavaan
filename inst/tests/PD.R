library(lavaan)

model <- '
  # measurement model
    dem60 =~ y1 + a*y2 + b*y3 + c*y4
    dem65 =~ y5 + a*y6 + b*y7 + c*y8
    ind60 =~ x1 + x2 + x3
  # regressions
    dem60 ~ b1*ind60 
    dem65 ~ b2*ind60 + b3*dem60
  # residual correlations
    y1 ~~ y5 
    y2 ~~ y4 
    y2 ~~ y6
    y3 ~~ y7
    y4 ~~ y8
    y6 ~~ y8 
'
fit <- sem(model, data=PoliticalDemocracy)

# default extract functions
source("common.srcR", echo = TRUE)
