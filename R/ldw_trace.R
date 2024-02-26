# create an environment to put cached object in
# this is executed when the package is 'compiled' !
trace_env <- new.env(parent = emptyenv())

ldw_trace <- function(content = "") {
  if (!exists("TRACE", trace_env)) return(invisible(NULL))
  if (!exists("TRACENR", trace_env)) assign("TRACENR", 1L, trace_env)
  tracenr <- get("TRACENR", trace_env)
  x <- sub("[() ].*$", "", as.character(sys.calls()))
  if (length(x) == 0) return(invisible(NULL))
  a <- paste0("trc", formatC(tracenr, format="d", width = 5, flag = "0"))
  x <- x[!(x %in% c("eval", "try", "tryCatch", "tryCatchList", "tryCatchOne", "doTryCatch",
                    "which", "unique", "as.list", "as.character", "unlist", "ldw_trace",
                    "source", "withVisible"))]
  if (length(x) > 0) {
    assign(a, list(stack = x, content = content, time = Sys.time()), trace_env)
    assign("TRACENR", tracenr + 1L, trace_env)
  }
  return(invisible(NULL))
}
set_trace <- function(state = NULL, silent = FALSE) {
  traceon <- exists("TRACE", trace_env)
  msg <- ""
  if (is.null(state)) {
    rm(list = ls(trace_env, pattern = "^trc"), envir = trace_env)
    if (exists("TRACENR", trace_env)) rm("TRACENR", envir = trace_env)
    msg <- "Traces removed." 
  } else if (state) {
    if (traceon) {
      msg <- "Trace already active!"
    } else {
      assign("TRACE", TRUE, trace_env)
      msg <- "Trace on."
    }
  } else {
    if (traceon) {
      rm("TRACE", envir = trace_env)
      msg <- "Trace off." 
    } else {
      msg <- "Trace not active!"
    }
  }
  if (!silent) cat(msg, "\n", sep = "") 
  return(invisible(NULL))
}
get_trace <- function() {
  traceobjects <- ls(trace_env, pattern = "^trc")
  if (length(traceobjects) == 0) return(list())
  x <- mget(traceobjects, envir = trace_env)
  x <- x[order(names(x))]
  return(x)
}
print_trace <- function(file = "", clean_after = (file != "")) {
  cat("Trace print on ", format(Sys.time(), format = "%F"), "\n\n", file = file)
  x <- get_trace()
  for (x1 in x) {
    cat( format(x1$time, format = "%T"), 
      paste(x1$stack, collapse = ">"), ":", x1$content, 
        "\n", sep = " ", file = file, append = TRUE)
  }
  if (clean_after) set_trace(NULL, TRUE)
}
