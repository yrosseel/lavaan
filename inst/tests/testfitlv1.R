 library(lavaan)
options(warn = 1L)

# one-factor + eXo
set.seed(1234)
pop.model <- ' f =~ 0.7*y1 + 0.7*y2 + 0.7*y3 + 0.7*y4 + 0.7*y5
               f ~ (-2.3)*x1 + 0.8*x2
               y1 ~ 0.2*x2 
               y3 ~ 0.7*x1 '
Data <- simulateData(pop.model, sample.nobs=100)

model <- ' f =~ y1 + y2 + y3 + y4 + y5
           f ~ x1 + x2
           y1 ~ x2
           y3 ~ x1 '
fit <- sem(model, data=Data, fixed.x=TRUE)

# default extract functions
source("common.srcR", echo = TRUE)


# create missing values
missing.per.var <- floor(nrow(Data) / 10)
Data.missing <- as.data.frame(lapply(Data, function(x) {
    idx <- sample(1:length(x), missing.per.var); x[idx] <- NA; x}))

# listwise deletion
fit1 <- sem(model, data=Data.missing, fixed.x=FALSE, missing="listwise")
# FIML
fit2 <- sem(model, data=Data.missing, fixed.x=FALSE, missing="ml")
fit <- fit2

# default extract functions
source("common.srcR", echo = TRUE)


# create binary version 
Data.binary <- Data
Data.binary$y1 <- cut(Data$y1, 2L, labels=FALSE)
Data.binary$y2 <- cut(Data$y2, 2L, labels=FALSE)
Data.binary$y3 <- cut(Data$y3, 2L, labels=FALSE)
Data.binary$y4 <- cut(Data$y4, 2L, labels=FALSE)
Data.binary$y5 <- cut(Data$y5, 2L, labels=FALSE)


fit <- sem(model, data=Data.binary, estimator="WLSMV", 
           ordered=c("y1","y2","y3","y4","y5"))

# default extract functions
source("common.srcR", echo = TRUE)


# create missing values
missing.per.var <- floor(nrow(Data) / 10)
Data.missing <- as.data.frame(lapply(Data.binary, function(x) {
    idx <- sample(1:length(x), missing.per.var); x[idx] <- NA; x}))

# listwise deletion
fit1 <- sem(model, data=Data.missing, missing="listwise", estimator="WLSMV",
            ordered=c("y1","y2","y3","y4","y5"))
# pairwise
fit2 <- sem(model, data=Data.missing, missing="pairwise", estimator="WLSMV",
            ordered=c("y1","y2","y3","y4","y5"))
fit <- fit2

# default extract functions
source("common.srcR", echo = TRUE)





