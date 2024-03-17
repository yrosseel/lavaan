# data.frame utilities
# Y.R. 11 April 2013

# - 10 nov 2019: * removed lav_dataframe_check_vartype(), as we can simply use
#                  sapply(lapply(frame, class), "[", 1L) (unused anyway)
#                * removed lav_dataframe_check_ordered() as we can simply use
#                  any(sapply(frame[, ov.names], inherits, "ordered"))

# construct vartable, but allow 'ordered/factor' argument to intervene
# we do NOT change the data.frame
lav_dataframe_vartable <- function(frame = NULL, ov.names = NULL,
                                   ov.names.x = NULL,
                                   ordered = NULL,
                                   factor = NULL,
                                   as.data.frame. = FALSE) {
  if (missing(ov.names)) {
    var.names <- names(frame)
  } else {
    ov.names <- unlist(ov.names, use.names = FALSE)
    ov.names.x <- unlist(ov.names.x, use.names = FALSE)
    var.names <- unique(c(ov.names, ov.names.x))
  }
  nvar <- length(var.names)
  var.idx <- match(var.names, names(frame))


  nobs <- integer(nvar)
  type <- character(nvar)
  user <- integer(nvar)
  exo <- ifelse(var.names %in% ov.names.x, 1L, 0L)
  mean <- numeric(nvar)
  var <- numeric(nvar)
  nlev <- integer(nvar)
  lnam <- character(nvar)
  for (i in seq_len(nvar)) {
    x <- frame[[var.idx[i]]]

    type.x <- class(x)[1L]

    # correct for matrix with 1 column
    if (inherits(x, "matrix") && (is.null(dim(x)) ||
      (!is.null(dim) && ncol(x) == 1L))) {
      type.x <- "numeric"
    }

    # correct for integers
    if (inherits(x, "integer")) {
      type.x <- "numeric"
    }

    # handle the 'labelled' type from the haven package
    # - if the variable name is not in 'ordered', we assume
    #   it is numeric (for now) 11 March 2018
    if (inherits(x, "labelled") && !(var.names[i] %in% ordered)) {
      type.x <- "numeric"
    }

    # handle ordered/factor
    if (!is.null(ordered) && var.names[i] %in% ordered) {
      type.x <- "ordered"
      lev <- sort(unique(x)) # we assume integers!
      nlev[i] <- length(lev)
      lnam[i] <- paste(lev, collapse = "|")
      user[i] <- 1L
    } else if (!is.null(factor) && var.names[i] %in% factor) {
      type.x <- "factor"
      lev <- sort(unique(x)) # we assume integers!
      nlev[i] <- length(lev)
      lnam[i] <- paste(lev, collapse = "|")
      user[i] <- 1L
    } else {
      nlev[i] <- nlevels(x)
      lnam[i] <- paste(levels(x), collapse = "|")
    }

    type[i] <- type.x
    nobs[i] <- sum(!is.na(x))
    mean[i] <- ifelse(type.x == "numeric", mean(x, na.rm = TRUE),
      as.numeric(NA)
    )
    var[i] <- ifelse(type.x == "numeric", var(x, na.rm = TRUE),
      as.numeric(NA)
    )
  }

  VAR <- list(
    name = var.names, idx = var.idx, nobs = nobs, type = type, exo = exo,
    user = user, mean = mean, var = var, nlev = nlev, lnam = lnam
  )

  if (as.data.frame.) {
    VAR <- as.data.frame(VAR,
      stringsAsFactors = FALSE,
      row.names = 1:length(VAR$name)
    )
    class(VAR) <- c("lavaan.data.frame", "data.frame")
  }

  VAR
}
