 library(lavaan)
options(warn = 1L)

# 1-variable model
set.seed(1234)
Data <- data.frame(y1 = rnorm(100))

# variance only
model <- 'y1 ~~ y1'
fit <- lavaan(model, data=Data)

# default extract functions
source("common.srcR", echo = TRUE)


# intercept + variance only
model <- 'y1 ~ 1; y1 ~~ y1' 
fit <- lavaan(model, data=Data)

# default extract functions
source("common.srcR", echo = TRUE)

# create missing values
missing.per.var <- floor(nrow(Data) / 10)
Data.missing <- as.data.frame(lapply(Data, function(x) {
    idx <- sample(1:length(x), missing.per.var); x[idx] <- NA; x}))

model <- 'y1 ~ 1; y1 ~~ y1'

# listwise deletion
fit1 <- lavaan(model, data=Data.missing, fixed.x=FALSE, missing="listwise")
# FIML
fit2 <- lavaan(model, data=Data.missing, fixed.x=FALSE, missing="ml")
# should be the same here, since there is only 1 column, and all
# empty cases are removed
all.equal(coef(fit1), coef(fit2))



# create binary version
Data.binary <- as.data.frame( lapply(Data, cut, 2, labels=FALSE) )
Data.binary <- as.data.frame( lapply(Data.binary, ordered) )

model <- 'y1 | t1'
fit <- sem(model, data=Data.binary, estimator="WLSMV")

# default extract functions
source("common.srcR", echo = TRUE)




