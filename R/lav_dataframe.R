# data.frame utilities
# Y.R. 11 April 2013

# - 10 nov 2019: * removed lav_dataframe_check_vartype(), as we can simply use
#                  sapply(lapply(frame, class), "[", 1L) (unused anyway)
#                * removed lav_dataframe_check_ordered() as we can simply use
#                  any(sapply(frame[, ov.names], inherits, "ordered"))

# construct vartable, but allow 'ordered/factor' argument to intervene
# we do NOT change the data.frame
lav_dataframe_vartable <- function(frame = NULL, ov_names = NULL,
                                   ov_names_x = NULL,
                                   ordered = NULL,
                                   factor = NULL,
                                   as_data_frame = FALSE,
                                   allow_empty_cell = FALSE) {
  if (missing(ov_names)) {
    var_names <- names(frame)
  } else {
    ov_names <- unlist(ov_names, use.names = FALSE)
    ov_names_x <- unlist(ov_names_x, use.names = FALSE)
    var_names <- unique(c(ov_names, ov_names_x))
  }
  nvar <- length(var_names)
  var_idx <- match(var_names, names(frame))


  nobs <- integer(nvar)
  type <- character(nvar)
  user <- integer(nvar)
  exo <- ifelse(var_names %in% ov_names_x, 1L, 0L)
  mean <- numeric(nvar)
  var <- numeric(nvar)
  nlev <- integer(nvar)
  lnam <- character(nvar)
  for (i in seq_len(nvar)) {
    x <- frame[[var_idx[i]]]

    type_x <- class(x)[1L]

    # correct for matrix with 1 column
    if (inherits(x, "matrix") && (is.null(dim(x)) ||
      (!is.null(dim) && ncol(x) == 1L))) {
      type_x <- "numeric"
    }

    # correct for integers
    if (inherits(x, "integer")) {
      type_x <- "numeric"
    }

    # handle the 'labelled' type from the haven package
    # - if the variable name is not in 'ordered', we assume
    #   it is numeric (for now) 11 March 2018
    # - haven >= 2.0 renamed the class to 'haven_labelled' (and
    #   'haven_labelled_spss' for SPSS user-defined missings); the old
    #   'labelled' class only matched pre-2.0 haven, so labelled numeric
    #   variables were left with type = "haven_labelled" (which, e.g.,
    #   sent lavPredict() down the numerical EBM path). Detect a labelled
    #   *numeric* variable directly instead (July 2026).
    if (is.numeric(x) &&
        inherits(x, c("labelled", "haven_labelled")) &&
        !(var_names[i] %in% ordered)) {
      type_x <- "numeric"
    }

    # handle ordered/factor
    if (!is.null(ordered) && var_names[i] %in% ordered) {
      type_x <- "ordered"
      if (allow_empty_cell) {
        if (inherits(x, 'factor')) {
          nlev[i] <- nlevels(x)
          lnam[i] <- paste(levels(x), collapse = "|")
        } else {
          nlev[i] <- max(as.numeric(x), na.rm = TRUE)
          lnam[i] <- paste(1:nlev[i], collapse = "|")
        }
      } else {
        lev <- sort(unique(x)) # we assume integers!
        nlev[i] <- length(lev)
        lnam[i] <- paste(lev, collapse = "|")
      }
      user[i] <- 1L
    } else if (!is.null(factor) && var_names[i] %in% factor) {
      type_x <- "factor"
      lev <- sort(unique(x)) # we assume integers!
      nlev[i] <- length(lev)
      lnam[i] <- paste(lev, collapse = "|")
      user[i] <- 1L
    } else {
      nlev[i] <- nlevels(x)
      lnam[i] <- paste(levels(x), collapse = "|")
    }

    type[i] <- type_x
    nobs[i] <- sum(!is.na(x))
    mean[i] <- ifelse(type_x == "numeric", mean(x, na.rm = TRUE),
      as.numeric(NA)
    )
    var[i] <- ifelse(type_x == "numeric", var(x, na.rm = TRUE),
      as.numeric(NA)
    )
  }

  var_1 <- list(
    name = var_names, idx = var_idx, nobs = nobs, type = type, exo = exo,
    user = user, mean = mean, var = var, nlev = nlev, lnam = lnam
  )

  if (as_data_frame) {
    var_1 <- as.data.frame(var_1,
      stringsAsFactors = FALSE,
      row.names = seq_along(var_1$name)
    )
    class(var_1) <- c("lavaan.data.frame", "data.frame")
  }

  var_1
}
