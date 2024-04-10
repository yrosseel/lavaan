# for blavaan
# TDJ: add "..." to make the generic actually generic, for lavaan.mi objects

# S3 generic for S3 dispatch
fitMeasures <- function(object, fit.measures = "all",
                        baseline.model = NULL, h1.model = NULL,
                        fm.args = list(
                          standard.test = "default",
                          scaled.test = "default",
                          rmsea.ci.level = 0.90,
                          rmsea.close.h0 = 0.05,
                          rmsea.notclose.h0 = 0.08,
                          robust = TRUE,
                          cat.check.pd = TRUE
                        ),
                        output = "vector", ...) {
  UseMethod("fitMeasures", object)
}
fitmeasures <- function(object, fit.measures = "all",
                        baseline.model = NULL, h1.model = NULL,
                        fm.args = list(
                          standard.test = "default",
                          scaled.test = "default",
                          rmsea.ci.level = 0.90,
                          rmsea.close.h0 = 0.05,
                          rmsea.notclose.h0 = 0.08,
                          robust = TRUE,
                          cat.check.pd = TRUE
                        ),
                        output = "vector", ...) {
  UseMethod("fitmeasures", object)
}


# S4 generic for S4 dispatch
setGeneric(
  "fitMeasures",
  function(object, fit.measures = "all",
           baseline.model = NULL, h1.model = NULL,
           fm.args = list(
             standard.test = "default",
             scaled.test = "default",
             rmsea.ci.level = 0.90,
             rmsea.close.h0 = 0.05,
             rmsea.notclose.h0 = 0.08,
             robust = TRUE,
             cat.check.pd = TRUE
           ),
           output = "vector", ...) {
    standardGeneric("fitMeasures")
  }
)
setGeneric(
  "fitmeasures",
  function(object, fit.measures = "all",
           baseline.model = NULL, h1.model = NULL,
           fm.args = list(
             standard.test = "default",
             scaled.test = "default",
             rmsea.ci.level = 0.90,
             rmsea.close.h0 = 0.05,
             rmsea.notclose.h0 = 0.08,
             robust = TRUE,
             cat.check.pd = TRUE
           ),
           output = "vector", ...) {
    standardGeneric("fitmeasures")
  }
)


# S3 generics
inspect <- function(object, what = "free", ...) {
  UseMethod("inspect", object)
}

lavInspect <- function(object, what = "free",
                       add.labels = TRUE,
                       add.class = TRUE,
                       list.by.group = TRUE,
                       drop.list.single.group = TRUE) {
  UseMethod("lavInspect", object)
}

lavTech <- function(object, what = "free",
                    add.labels = FALSE,
                    add.class = FALSE,
                    list.by.group = FALSE,
                    drop.list.single.group = FALSE) {
  UseMethod("lavTech", object)
}
