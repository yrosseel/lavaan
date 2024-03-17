# generate syntax for an independence model
lav_syntax_independence <- function(ov.names = character(0),
                                    ov.names.x = character(0),
                                    sample.cov = NULL) {
  ov.names.nox <- ov.names[!ov.names %in% ov.names.x]
  nvar <- length(ov.names.nox)
  lv.names <- paste("f", 1:nvar, sep = "")

  # check sample.cov
  if (!is.null(sample.cov)) {
    if (is.list(sample.cov)) {
      ngroups <- length(sample.cov)
    } else {
      ngroups <- 1L
      sample.cov <- list(sample.cov)
    }
    stopifnot(is.matrix(sample.cov[[1]]))
    # stopifnot(length(ov.names) == nrow(sample.cov[[1]]))
    # FIXME: check rownames and reorder...
  }

  # construct lavaan syntax for an independence model
  txt <- "# independence model\n"

  # =~ lines (each observed variables has its own latent variable)
  # excepct for ov's that are in ov.names.x
  txt <- paste(txt, paste(lv.names, " =~ 1*", ov.names.nox,
    "\n",
    sep = "", collapse = ""
  ), sep = "")

  # residual ov variances fixed to zero
  txt <- paste(txt, paste(ov.names.nox, " ~~ 0*", ov.names.nox,
    "\n",
    sep = "", collapse = ""
  ), sep = "")

  # latent variances
  if (is.null(sample.cov)) {
    txt <- paste(txt, paste(lv.names, " ~~ ", lv.names,
      "\n",
      sep = "", collapse = ""
    ), sep = "")
  } else {
    # fill in sample values
    ov.idx <- match(ov.names.nox, ov.names)

    start.txt <- paste("start(c(",
      apply(matrix(
        unlist(lapply(sample.cov, function(x) {
          diag(x)[ov.idx]
        })),
        ncol = ngroups
      ), 1, paste, collapse = ","), "))",
      sep = ""
    )
    txt <- paste(txt, paste(lv.names, " ~~ ", start.txt, " * ",
      lv.names,
      "\n",
      sep = "", collapse = ""
    ), sep = "")
  }

  # latent *covariances* fixed to zero (= independence!)
  if (length(lv.names) > 1L) {
    tmp <- utils::combn(lv.names, 2)
    txt <- paste(txt, paste(tmp[1, ], " ~~ 0*", tmp[2, ], "\n",
      sep = "",
      collapse = ""
    ), sep = "")
  }

  # if 'independent x' variables, add an 'empty' regression
  if ((nx <- length(ov.names.x)) > 0) {
    # dummy regression line
    txt <- paste(txt, paste("f1 ~ 0*", ov.names.x,
      "\n",
      sep = "", collapse = ""
    ), sep = "")
  }

  # Note: no need to pass starting values here, lavaanStart will
  # use the sample statistics anyway....

  txt
}
