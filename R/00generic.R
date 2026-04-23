# for blavaan
# TDJ: add "..." to make the generic actually generic, for lavaan.mi objects

# S3 generic for S3 dispatch
fitMeasures <- function(object, fit.measures = "all",           # nolint start
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
                        output = "vector", ...) {               # nolint end
  UseMethod("fitMeasures", object)
}
fitmeasures <- function(object, fit.measures = "all",           # nolint start
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
                        output = "vector", ...) {               # nolint end
  UseMethod("fitmeasures", object)
}


# S4 generic for S4 dispatch
setGeneric(
  "fitMeasures",                                              # nolint start
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
           output = "vector", ...) {                          # nolint end
    standardGeneric("fitMeasures")
  }
)
setGeneric(
  "fitmeasures",
  function(object, fit.measures = "all",                      # nolint start
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
           output = "vector", ...) {                     # nolint end
    standardGeneric("fitmeasures")
  }
)


# S3 generics
lavInspect <- function(object, what = "free",           # nolint start
                       add.labels = TRUE,
                       add.class = TRUE,
                       list.by.group = TRUE,
                       drop.list.single.group = TRUE) { # nolint end
  UseMethod("lavInspect", object)
}

inspect <- function(object, what = "free", ...) {
    UseMethod("inspect", object)
}

lavTech <- function(object, what = "free",              # nolint start
                    add.labels = FALSE,
                    add.class = FALSE,
                    list.by.group = FALSE,
                    drop.list.single.group = FALSE) {   # nolint end
  UseMethod("lavTech", object)
}
