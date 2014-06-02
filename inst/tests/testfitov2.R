 library(lavaan)
options(warn = 1L)

# 2-variable model: simple regression
set.seed(1234)
x1 = rnorm(100)
y1 = 0.5 + 2*x1 + rnorm(100)
Data <- data.frame(y1 = y1, x1 = x1)

model <- ' y1 ~ x1 '
fit <- sem(model, data=Data, fixed.x=TRUE)

# default extract functions
source("common.srcR", echo = TRUE)

fit <- sem(model, data=Data, fixed.x=FALSE)

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

# default extract functions
source("common.srcR", echo = TRUE)




# create binary version
Data.binary <- Data
Data.binary$y1 <- cut(Data$y1, 2L, labels=FALSE)

model <- 'y1 ~ x1'
fit <- sem(model, data=Data.binary, estimator="WLSMV", ordered=c("y1"))

# default extract functions
source("common.srcR", echo = TRUE)

