# create (if not already created) an environment to put cached objects in
# this is executed when the package is 'compiled' !
if (!exists("lavaan_cache_env")) lavaan_cache_env <- new.env(parent = emptyenv())

# tracing possibility in functions defined below, an example of use :
#
# in the function where you want to trace add a line
#     ldw_trace(x)
# where x is a characterstring you want to show in the trace
#
# thereafter execute a script like this:
# library(lavaan)
# lavaan:::set_trace(TRUE)
# model <- '
# # latent variable definitions
#    ind60 =~ x1 + x2 + x3
#    dem60 =~ y1 + a*y2 + b*y3 + c*y4
#    dem65 =~ y5 + a*y6 + b*y7 + c*y8
#
# # regressions
#   dem60 ~ ind60
#   dem65 ~ ind60 + dem60
#
# # residual correlations
#   y1 ~~ y5
#   y2 ~~ y4 + y6
#   y3 ~~ y7
#   y4 ~~ y8
#   y6 ~~ y8
# '
# fit <- sem(model, data = PoliticalDemocracy)
# summary(fit)
# lavaan:::set_trace(FALSE)
# lavaan:::print_trace("PolDem_trace.txt")
#

ldw_trace <- function(content = "") {
  ignore.in.stack <- c(
    "eval", "try", "tryCatch", "tryCatchList", "tryCatchOne", "doTryCatch",
    "which", "unique", "as.list", "as.character", "unlist", "ldw_trace",
    "source", "withVisible", "tryCatch.W.E", "withCallingHandlers", "do.call"
  )
  if (!exists("TRACE", lavaan_cache_env)) {
    return(invisible(NULL))
  }
  if (!exists("TRACENR", lavaan_cache_env)) assign("TRACENR", 1L, lavaan_cache_env)
  tracenr <- get("TRACENR", lavaan_cache_env)
  x <- sub("[() ].*$", "", as.character(sys.calls()))
  if (length(x) == 0) {
    return(invisible(NULL))
  }
  a <- paste0("trc", formatC(tracenr, format = "d", width = 5, flag = "0"))
  x <- x[!(x %in% ignore.in.stack)]
  if (length(x) > 0) {
    assign(a, list(stack = x, content = content, time = Sys.time()), lavaan_cache_env)
    assign("TRACENR", tracenr + 1L, lavaan_cache_env)
  }

  invisible(NULL)
}

set_trace <- function(state = NULL, silent = FALSE) {
  traceon <- exists("TRACE", lavaan_cache_env)
  msg <- ""
  if (is.null(state)) {
    rm(list = ls(lavaan_cache_env, pattern = "^trc"), envir = lavaan_cache_env)
    if (exists("TRACENR", lavaan_cache_env)) rm("TRACENR", envir = lavaan_cache_env)
    msg <- "Traces removed."
  } else if (state) {
    if (traceon) {
      msg <- "Trace already active!"
    } else {
      assign("TRACE", TRUE, lavaan_cache_env)
      msg <- "Trace on."
    }
  } else {
    if (traceon) {
      rm("TRACE", envir = lavaan_cache_env)
      msg <- "Trace off."
    } else {
      msg <- "Trace not active!"
    }
  }
  if (!silent) cat(msg, "\n", sep = "")

  invisible(NULL)
}

get_trace <- function() {
  traceobjects <- ls(lavaan_cache_env, pattern = "^trc")
  if (length(traceobjects) == 0) {
    return(list())
  }
  x <- mget(traceobjects, envir = lavaan_cache_env)
  x <- x[order(names(x))]

  x
}

print_trace <- function(file = "", clean_after = (file != "")) {
  cat("Trace print on ", format(Sys.time(), format = "%F"), "\n\n", file = file)
  x <- get_trace()
  for (x1 in x) {
    cat(format(x1$time, format = "%T"),
      paste(x1$stack, collapse = ">"), ":", x1$content,
      "\n",
      sep = " ", file = file, append = TRUE
    )
  }
  if (clean_after) set_trace(NULL, TRUE)
}
