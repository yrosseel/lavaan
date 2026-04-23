# create (if not already created) an environment to put cached objects in
# this is executed when the package is 'compiled' !
if (!exists("lavaan_cache_env"))
    lavaan_cache_env <- new.env(parent = emptyenv())

# tracing possibility in functions defined below, an example of use :
#
# in the function where you want to trace add a line
#     lav_trace(x)
# where x is a characterstring you want to show in the trace
#
# thereafter execute a script like this:
# library(lavaan)
# lavaan:::lav_trace_set(TRUE)
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
# lavaan:::lav_trace_set(FALSE)
# lavaan:::lav_trace_print("PolDem_trace.txt")
#

lav_trace <- function(content = "") {
  ignore_in_stack <- c(
    "eval", "try", "tryCatch", "tryCatchList", "tryCatchOne", "doTryCatch",
    "which", "unique", "as.list", "as.character", "unlist", "lav_trace",
    "source", "withVisible", "tryCatch.W.E", "withCallingHandlers", "do.call"
  )
  if (!exists("trace", lavaan_cache_env)) {
    return(invisible(NULL))
  }
  if (!exists("tracenr", lavaan_cache_env))
              assign("tracenr", 1L, lavaan_cache_env)
  tracenr <- get("tracenr", lavaan_cache_env)
  x <- sub("[() ].*$", "", as.character(sys.calls()))
  if (length(x) == 0) {
    return(invisible(NULL))
  }
  a <- paste0("trc", formatC(tracenr, format = "d", width = 5, flag = "0"))
  x <- x[!(x %in% ignore_in_stack)]
  if (length(x) > 0) {
    assign(a, list(stack = x, content = content, time = Sys.time()),
                        lavaan_cache_env)
    assign("tracenr", tracenr + 1L, lavaan_cache_env)
  }

  invisible(NULL)
}

lav_trace_set <- function(state = NULL, silent = FALSE) {
  traceon <- exists("trace", lavaan_cache_env)
  msg <- ""
  if (is.null(state)) {
    rm(list = ls(lavaan_cache_env, pattern = "^trc"), envir = lavaan_cache_env)
    if (exists("tracenr", lavaan_cache_env))
        rm("tracenr", envir = lavaan_cache_env)
    msg <- "Traces removed."
  } else if (state) {
    if (traceon) {
      msg <- "Trace already active!"
    } else {
      assign("trace", TRUE, lavaan_cache_env)
      msg <- "Trace on."
    }
  } else {
    if (traceon) {
      rm("trace", envir = lavaan_cache_env)
      msg <- "Trace off."
    } else {
      msg <- "Trace not active!"
    }
  }
  if (!silent) cat(msg, "\n", sep = "")

  invisible(NULL)
}

lav_trace_get <- function() {
  traceobjects <- ls(lavaan_cache_env, pattern = "^trc")
  if (length(traceobjects) == 0) {
    return(list())
  }
  x <- mget(traceobjects, envir = lavaan_cache_env)
  x <- x[order(names(x))]

  x
}

lav_trace_print <- function(file = "", clean_after = (file != "")) {
  cat("Trace print on ", format(Sys.time(), format = "%F"), "\n\n", file = file)
  x <- lav_trace_get()
  for (x1 in x) {
    cat(format(x1$time, format = "%T"),
      paste(x1$stack, collapse = ">"), ":", x1$content,
      "\n",
      sep = " ", file = file, append = TRUE
    )
  }
  if (clean_after) lav_trace_set(NULL, TRUE)
}

lav_trace_summary <- function(file = "", clean_after = FALSE) {
  x <- lav_trace_get()
  temp <- new.env(parent = emptyenv())
  for (x1 in x) {
    nn <- length(x1$stack)
    mm <- paste(x1$stack[nn], paste(x1$stack[seq_len(nn - 1L)],
          collapse = ">"), sep = "\t")
    assign(mm, 1L + get0(mm, temp, ifnotfound = 0L), temp)
  }
  objects <- sort(ls(temp))
  for (i in seq_along(objects)) {
    cat(objects[i], get(objects[i], temp), "\n", file = file, append = TRUE)
  }
  if (clean_after) lav_trace_set(NULL, TRUE)
}
