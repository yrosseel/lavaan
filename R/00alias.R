# create a thin forwarding alias for the function `target` (a name):
# the alias has the same formals as the target -- so documented usage,
# argument matching and missing() semantics are unaffected -- but its
# body merely re-dispatches the call to the target.
#
# rationale: a plain copy (`alias <- target`) stores the complete
# function a *second* time in the lazy-loading DB of the installed
# package; for the larger user-facing functions and their synonyms
# (parameterEstimates/parameterestimates, modificationIndices/..., ...)
# these copies added up to more than 100 Kb of installed size.
#
# note: the target must already exist when the alias is created, so the
# alias must come later in the collation order (00alias.R itself sorts
# before all lav_*.R files).
# the modindices() default `maximum_number = nrow(list_1)` refers to a
# variable that is only created inside the body; the thin forwarders
# below (which share the formals but not the body) would otherwise
# trigger a false-positive codetools NOTE
utils::globalVariables("list_1")

lav_alias <- function(target) {
  pkg_env <- parent.frame() # the package environment under construction
  fun <- eval(
    bquote(function() {
      sc <- sys.call()
      sc[[1L]] <- .(as.name(target)) # resolved in the lavaan namespace
      eval(sc, envir = parent.frame())
    }),
    envir = pkg_env
  )
  formals(fun) <- formals(get(target, envir = pkg_env, mode = "function"))
  fun
}
