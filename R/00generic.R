# for blavaan
# TDJ: add "..." to make the generic actually generic, for lavaan.mi objects

# S3 generic for S3 dispatch
fitMeasures <- function(object, fit_measures = "all",           # nolint
                        baseline_model = NULL, h1_model = NULL,
                        fm_args = list(
                          standard.test = "default",
                          scaled.test = "default",
                          rmsea.ci.level = 0.90,
                          rmsea.close.h0 = 0.05,
                          rmsea.notclose.h0 = 0.08,
                          robust = TRUE,
                          cat.nonpd = "na"
                        ),
                        output = "vector", level = NULL, ...) {
  UseMethod("fitMeasures", object)
}
fitmeasures <- function(object, fit_measures = "all",
                        baseline_model = NULL, h1_model = NULL,
                        fm_args = list(
                          standard.test = "default",
                          scaled.test = "default",
                          rmsea.ci.level = 0.90,
                          rmsea.close.h0 = 0.05,
                          rmsea.notclose.h0 = 0.08,
                          robust = TRUE,
                          cat.nonpd = "na"
                        ),
                        output = "vector", level = NULL, ...) {
  UseMethod("fitmeasures", object)
}


# S4 generic for S4 dispatch
setGeneric(
  "fitMeasures",
  function(object, fit_measures = "all",
           baseline_model = NULL, h1_model = NULL,
           fm_args = list(
             standard.test = "default",
             scaled.test = "default",
             rmsea.ci.level = 0.90,
             rmsea.close.h0 = 0.05,
             rmsea.notclose.h0 = 0.08,
             robust = TRUE,
             cat.nonpd = "na"
           ),
           output = "vector", level = NULL, ...) {
    standardGeneric("fitMeasures")
  }
)
setGeneric(
  "fitmeasures",
  function(object, fit_measures = "all",
           baseline_model = NULL, h1_model = NULL,
           fm_args = list(
             standard.test = "default",
             scaled.test = "default",
             rmsea.ci.level = 0.90,
             rmsea.close.h0 = 0.05,
             rmsea.notclose.h0 = 0.08,
             robust = TRUE,
             cat.nonpd = "na"
           ),
           output = "vector", level = NULL, ...) {
    standardGeneric("fitmeasures")
  }
)


# S3 generics
lavInspect <- function(object, what = "free",
                       add_labels = TRUE,
                       add_class = TRUE,
                       list_by_group = TRUE,
                       drop_list_single_group = TRUE,
                       ...) {
  UseMethod("lavInspect", object)
}

inspect <- function(object, what = "free", ...) {
    UseMethod("inspect", object)
}

lavTech <- function(object, what = "free",
                    add_labels = FALSE,
                    add_class = FALSE,
                    list_by_group = FALSE,
                    drop_list_single_group = FALSE,
                    ...) {
  UseMethod("lavTech", object)
}
