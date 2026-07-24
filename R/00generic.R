# for blavaan
# TDJ: add "..." to make the generic actually generic, for lavaan.mi objects

# S3 generic for S3 dispatch
#
# NOTE: the generic only lists 'object' and '...'; the full argument list
# (fit_measures, baseline_model, h1_model, fm_args, output, level) lives on
# the methods, where the defaults are defined once. This keeps the generic
# minimal and avoids having to keep the generic and method signatures in sync.
fitMeasures <- function(object, ...) {                          # nolint
  UseMethod("fitMeasures", object)
}
fitmeasures <- function(object, ...) {
  UseMethod("fitmeasures", object)
}


# S4 generic for S4 dispatch
setGeneric("fitMeasures", function(object, ...) {
  standardGeneric("fitMeasures")
})
setGeneric("fitmeasures", function(object, ...) {
  standardGeneric("fitmeasures")
})


# S3 generics
#
# As for fitMeasures() above, the generics only list 'object' and '...'; the
# full argument list (what, add_labels, add_class, list_by_group,
# drop_list_single_group) with its defaults lives on the methods.
lavInspect <- function(object, ...) {
  UseMethod("lavInspect", object)
}

inspect <- function(object, ...) {
    UseMethod("inspect", object)
}

lavTech <- function(object, ...) {
  UseMethod("lavTech", object)
}
