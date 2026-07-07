## NOTE:
## round(1.2355, 3) = 1.236
## but
## round(1.2345, 3) = 1.234
##
## perhaps we should add 0.0005 or something to avoid this?

lav_dataframe_print <- function(x, ..., nd = 3L) {
  row_names <- rownames(x)
  y <- as.data.frame(lapply(x, function(x) {
    if (is.numeric(x)) round(x, nd) else x
  }))
  rownames(y) <- row_names

  if (!is.null(attr(x, "header"))) {
    cat("\n", attr(x, "header"), "\n\n", sep = "")
  }

  print(y, ...)

  if (!is.null(attr(x, "footer"))) {
    cat("\n", attr(x, "footer"), "\n\n", sep = "")
  }

  invisible(x)
}

lav_lavaanlist_print <- function(x, ...) {
  y <- unclass(x)
  attr(y, "header") <- NULL

  header <- attr(x, "header")
  if (!is.null(header)) {
    if (is.character(header)) {
      cat("\n", header, "\n\n", sep = "")
    } else {
      print(header)
      cat("\n")
    }
  }

  print(y, ...)
  invisible(x)
}


# prints only lower triangle of a symmetric matrix
lav_mat_sym_print <- function(x, ..., nd = 3L, shift = 0L,
                                          diag_na_dot = TRUE) {
  # print only lower triangle of a symmetric matrix
  # this function was inspired by the `print.correlation' function
  # in package nlme
  x <- as.matrix(x) # just in case
  y <- x
  y <- unclass(y)
  attributes(y)[c("header", "footer")] <- NULL
  ll <- lower.tri(x, diag = TRUE)
  y[ll] <- format(round(x[ll], digits = nd))
  y[!ll] <- ""
  if (diag_na_dot) {
    # print a "." instead of NA on the main diagonal (eg lav_efalist_summary)
    diag_idx <- lav_mat_diag_idx(ncol(x))
    tmp <- x[diag_idx]
    if (all(is.na(tmp))) {
      y[diag_idx] <- paste(strrep(" ", nd + 2L), ".", sep = "")
    }
  }
  if (!is.null(colnames(x))) {
    colnames(y) <- abbreviate(colnames(x), minlength = nd + 3L)
  }
  if (shift > 0L) {
    empty_string <- rep(strrep(x = " ", times = shift), times = nrow(x))
    if (!is.null(rownames(x))) {
      rownames(y) <- paste(empty_string, rownames(x), sep = "")
    } else {
      rownames(y) <- empty_string
    }
  }

  if (!is.null(attr(x, "header"))) {
    cat("\n", attr(x, "header"), "\n\n", sep = "")
  }

  print(y, ..., quote = FALSE, right = TRUE)

  if (!is.null(attr(x, "footer"))) {
    cat("\n", attr(x, "footer"), "\n\n", sep = "")
  }

  invisible(x)
}


lav_mat_print <- function(x, ..., nd = 3L, shift = 0L) {
  x <- as.matrix(x) # just in case
  y <- unclass(x)
  attributes(y)[c("header", "footer")] <- NULL
  if (!is.null(colnames(x))) {
    colnames(y) <- abbreviate(colnames(x), minlength = nd + 3L)
  }
  if (shift > 0L) {
    empty_string <- rep(strrep(x = " ", times = shift), times = nrow(x))
    if (!is.null(rownames(x))) {
      rownames(y) <- paste(empty_string, rownames(x), sep = "")
    } else {
      rownames(y) <- empty_string
    }
  }
  if (!is.null(attr(x, "header"))) {
    cat("\n", attr(x, "header"), "\n\n", sep = "")
  }

  print(round(y, nd), right = TRUE, ...)

  if (!is.null(attr(x, "footer"))) {
    cat("\n", attr(x, "footer"), "\n\n", sep = "")
  }

  invisible(x)
}

lav_vector_print <- function(x, ..., nd = 3L, shift = 0L) {
  y <- unclass(x)
  attributes(y)[c("header", "footer")] <- NULL
  # if(!is.null(names(x))) {
  #    names(y) <- abbreviate(names(x), minlength = nd + 3)
  # }
  if (!is.null(attr(x, "header"))) {
    cat("\n", attr(x, "header"), "\n\n", sep = "")
  }

  if (shift > 0L) {
    empty_string <- strrep(x = " ", times = shift)
    tmp <- format(y, digits = nd, width = 2L + nd)
    tmp[1] <- paste(empty_string, tmp[1], sep = "")
    print(tmp, quote = FALSE, ...)
  } else {
    print(round(y, nd), right = TRUE, ...)
  }

  if (!is.null(attr(x, "footer"))) {
    cat("\n", attr(x, "footer"), "\n\n", sep = "")
  }

  invisible(x)
}

lav_parameterestimates_print <- function(x, ..., nd = 3L) {
  # format for numeric values
  num_format <- paste("%", max(8L, nd + 5L), ".", nd, "f", sep = "")
  int_format <- paste("%", max(8L, nd + 5L), "d", sep = "")
  char_format <- paste("%", max(8L, nd + 5L), "s", sep = "")

  # output sections
  gsections <- c(
    "Latent Variables",
    "Composites",
    "Regressions",
    "Covariances",
    "Intercepts",
    "Thresholds",
    "Variances",
    "Scales y*",
    "Group Weight",
    "R-Square"
  )
  asections <- c(
    "Defined Parameters",
    "Constraints"
  )

  # header?
  header <- attr(x, "header")
  if (is.null(header)) {
    header <- FALSE
  }

  if (header) {
    cat("\nParameter Estimates:\n\n")

    # info about parameterization (if categorical only)
    categorical_flag <- attr(x, "categorical")
    if (categorical_flag) {
      # container
      c1 <- c2 <- character(0L)

      # which parameterization?
      c1 <- c(c1, "Parameterization")
      tmp_txt <- attr(x, "parameterization")
      c2 <- c(c2, paste(toupper(substring(tmp_txt, 1, 1)),
        substring(tmp_txt, 2),
        sep = ""
      ))

      # format c1/c2
      c1 <- format(c1, width = 37L)
      c2 <- format(c2,
        width = 14L + max(0, (nd - 3L)) * 4L, justify = "right"
      )

      # create character matrix
      m_1 <- cbind(c1, c2, deparse.level = 0)
      colnames(m_1) <- rep("", ncol(m_1))
      rownames(m_1) <- rep(" ", nrow(m_1))

      # print
      write.table(m_1, row.names = TRUE, col.names = FALSE, quote = FALSE)
    }

    # info about standard errors (if we have x$se only)
    # 1. information
    # 2. se
    # 3. bootstrap requested/successful draws
    if (!is.null(x$se)) {
      # container
      c1 <- c2 <- character(0L)

      # which type of standard errors?
      c1 <- c(c1, "Standard errors")
      if (attr(x, "se") == "robust.huber.white") {
        tmp_txt <- "sandwich" # since 0.6-6
      } else if (identical(attr(x, "estimator"), "IV") &&
                 attr(x, "se") != "bootstrap") {
        tmp_txt <- "two-stage" # estimator = "IV"
      } else {
        tmp_txt <- attr(x, "se")
      }
      c2 <- c(c2, paste(toupper(substring(tmp_txt, 1, 1)),
        substring(tmp_txt, 2),
        sep = ""
      ))
      # information
      if (attr(x, "se") != "bootstrap" &&
          identical(attr(x, "estimator"), "IV")) {
        # estimator = "IV": report the two-stage vcov methods and the type of
        # Gamma matrix instead of the information matrix
        c1 <- c(c1, "Standard errors stage 1")
        c2 <- c(c2, attr(x, "iv.vcov.stage1"))
        c1 <- c(c1, "Standard errors stage 2")
        c2 <- c(c2, attr(x, "iv.vcov.stage2"))
        if (!is.null(attr(x, "iv.gamma"))) {
          c1 <- c(c1, "Gamma matrix")
          c2 <- c(c2, attr(x, "iv.gamma"))
        }
      } else if (attr(x, "se") != "bootstrap") {
        # type for information
        if (attr(x, "se") == "robust.huber.white") {
          c1 <- c(c1, "Information bread")
        } else {
          c1 <- c(c1, "Information")
        }
        tmp_txt <- attr(x, "information")
        c2 <- c(c2, paste(toupper(substring(tmp_txt, 1, 1)),
          substring(tmp_txt, 2),
          sep = ""
        ))

        # if observed, which type? (hessian of h1)
        if (attr(x, "information") == "observed") {
          c1 <- c(c1, "Observed information based on")
          tmp_txt <- attr(x, "observed.information")
          c2 <- c(c2, paste(toupper(substring(tmp_txt, 1, 1)),
            substring(tmp_txt, 2),
            sep = ""
          ))
        }

        # if h1 is involved, structured or unstructured?
        if (attr(x, "information") %in% c("expected", "first.order") ||
          attr(x, "observed.information") == "h1") {
          if (attr(x, "se") == "robust.huber.white" &&
            attr(x, "h1.information") !=
              attr(x, "h1.information.meat")) {
            c1 <- c(c1, "Information bread saturated (h1) model")
          } else {
            c1 <- c(c1, "Information saturated (h1) model")
          }
          tmp_txt <- attr(x, "h1.information")
          c2 <- c(c2, paste(toupper(substring(tmp_txt, 1, 1)),
            substring(tmp_txt, 2),
            sep = ""
          ))
        }

        # if sandwich, which information for the meat? (first.order)
        # only print if it is NOT first.order
        if (attr(x, "se") == "robust.huber.white" &&
          attr(x, "information.meat") != "first.order") {
          c1 <- c(c1, "Information meat")
          tmp_txt <- attr(x, "information.meat")
          c2 <- c(c2, paste(toupper(substring(tmp_txt, 1, 1)),
            substring(tmp_txt, 2),
            sep = ""
          ))
        }

        # if sandwich, structured or unstructured for the meat?
        # only print if it is not the same as h1.information
        if (attr(x, "se") == "robust.huber.white" &&
          attr(x, "h1.information.meat") !=
            attr(x, "h1.information")) {
          c1 <- c(c1, "Information meat saturated (h1) model")
          tmp_txt <- attr(x, "h1.information.meat")
          c2 <- c(c2, paste(toupper(substring(tmp_txt, 1, 1)),
            substring(tmp_txt, 2),
            sep = ""
          ))
        }
      } # no bootstrap

      #TDJ: Pooling options for lavaan.mi-class objects (which NEVER bootstrap)
      if (isTRUE(attr(x, "pooled"))) {
        ## add an empty element for a space before pooling section
        c1 <- c(c1, "", "Pooled across imputations")
        c2 <- c(c2, "", "Rubin's (1987) rules")
      }
      if (!is.null(attr(x, "scale.W"))) {
        c1 <- c(c1, "Augment within-imputation variance")
        if (attr(x, "scale.W")) {
          c2 <- c(c2, "Scale by average RIV")
        } else {
          c2 <- c(c2, "Add between component")
        }
      }
      if (!is.null(attr(x, "asymptotic"))) {
        c1 <- c(c1, "Wald test for pooled parameters")
        if (attr(x, "asymptotic")) {
          c2 <- c(c2, "Normal (z) distribution")
        } else {
          c2 <- c(c2, "t(df) distribution")
        }
      }

      # 4.
      if (attr(x, "se") == "bootstrap" && !is.null(attr(x, "bootstrap"))) {
        c1 <- c(c1, "Number of requested bootstrap draws")
        c2 <- c(c2, attr(x, "bootstrap"))
        c1 <- c(c1, "Number of successful bootstrap draws")
        c2 <- c(c2, attr(x, "bootstrap.successful"))
      }

      # format c1/c2
      c1 <- format(c1, width = 37L)
      c2 <- format(c2,
        width = 14L + max(0, (nd - 3L)) * 4L, justify = "right"
      )

      # create character matrix
      m_1 <- cbind(c1, c2, deparse.level = 0)
      colnames(m_1) <- rep("", ncol(m_1))
      rownames(m_1) <- rep(" ", nrow(m_1))

      # print
      write.table(m_1, row.names = TRUE, col.names = FALSE, quote = FALSE)

      #TDJ: Message for lavaan.mi-class objects when df > 1000 for t test
      if (isTRUE(attr(x, "infDF"))) {
        cat(c("\n  Pooled t statistics with df >= 1000 are displayed with",
              "\n  df = Inf(inity) to save space. Although the t distribution",
              "\n  with large df closely approximates a standard normal",
              "\n  distribution, exact df for reporting these t tests can be",
              "\n  obtained from parameterEstimates.mi() \n\n"), sep = "")
      }

    }
  }

  # number of groups
  if (is.null(x$group)) {
    ngroups <- 1L
    x$group <- rep(1L, length(x$lhs))
  } else {
    ngroups <- lav_pt_ngroups(x)
  }

  # number of levels
  if (is.null(x$level)) {
    nlevels <- 1L
    x$level <- rep(1L, length(x$lhs))
  } else {
    nlevels <- lav_pt_nlevels(x)
  }

  # block column
  if (is.null(x$block)) {
    x$block <- rep(1L, length(x$lhs))
  }

  # step column (SAM)
  # if(!is.null(x$step)) {
  #    tmp.LABEL <- rep("", length(x$lhs))
  #    p1.idx <- which(x$step == 1L)
  #    p2.idx <- which(x$step == 2L)
  #    tmp.LABEL[p1.idx] <- "1"
  #    tmp.LABEL[p2.idx] <- "2"
  #
  #    if(is.null(x$label)) {
  #        x$label <- tmp.LABEL
  #    } else {
  #        x$label <- paste(x$label, tmp.LABEL, sep = "")
  #    }
  #
  #    x$step <- NULL
  # }

  # round to 3 digits after the decimal point
  y <- as.data.frame(
    lapply(x, function(x) {
      if (is.integer(x)) {
        sprintf(int_format, x)
      } else if (is.character(x)) { # perhaps plabel
        sprintf(char_format, x)
      } else if (is.numeric(x)) {
        sprintf(num_format, x)
      } else {
        x
      }
    }),
    stringsAsFactors = FALSE
  )

  # always remove /block/level/group/op/rhs/label/exo columns
  y$op <- y$group <- y$rhs <- y$label <- y$exo <- NULL
  y$block <- y$level <- NULL
  y$efa <- NULL

  # if standardized, remove std.nox column (space reasons only)
  # unless, std.all is already removed
  if (!is.null(y$std.all)) {
    y$std.nox <- NULL
  }

  # convert to character matrix
  m <- as.matrix(format.data.frame(y,
    na.encode = FALSE,
    justify = "right"
  ))

  # use empty row names
  rownames(m) <- rep("", nrow(m))

  # handle se == 0.0
  if (!is.null(x$se)) {
    se_idx <- which(x$se == 0)
    if (length(se_idx) > 0L) {
      m[se_idx, "se"] <- ""
      if (!is.null(x$z)) {
        m[se_idx, "z"] <- ""
      }
      if (!is.null(x$pvalue)) {
        m[se_idx, "pvalue"] <- ""
      }
      ## for lavaan.mi-class objects (semTools)
      if (!is.null(x$t)) {
        m[se_idx, "t"] <- ""
      }
      if (!is.null(x$df)) {
        m[se_idx, "df"] <- ""
      }
    }

    # handle se == NA
    se_idx <- which(is.na(x$se))
    if (length(se_idx) > 0L) {
      if (!is.null(x$z)) {
        m[se_idx, "z"] <- ""
      }
      if (!is.null(x$pvalue)) {
        m[se_idx, "pvalue"] <- ""
      }
      ## for lavaan.mi-class objects (semTools)
      if (!is.null(x$t)) {
        m[se_idx, "t"] <- ""
      }
      if (!is.null(x$df)) {
        m[se_idx, "df"] <- ""
      }
    }
  }

  # handle lower/upper boundary points
  # attention: to distinguish fixed values from free but bounded values,
  #            we have added a little bit of fuzz to the bounds of
  #            free parameters only
  if (!is.null(x$lower)) {
    diff_lower <- abs(x$lower - x$est)
    b_idx <- which(diff_lower > 0 & diff_lower < sqrt(.Machine$double.eps))
    if (length(b_idx) > 0L) {
      if (!is.null(x$pvalue)) {
        m[b_idx, "pvalue"] <- ""
      }
      if (is.null(x$label)) {
        x$label <- rep("", length(x$lhs))
      }
      x$label[b_idx] <- ifelse(nchar(x$label[b_idx]) > 0L,
        paste(x$label[b_idx], "+lb", sep = ""),
        "lb"
      )
    }
    # remove lower column
    m <- m[, colnames(m) != "lower"]
  }
  if (!is.null(x$upper)) {
    diff_upper <- abs(x$upper - x$est)
    b_idx <- which(diff_upper > 0 & diff_upper < sqrt(.Machine$double.eps))
    if (length(b_idx) > 0L) {
      if (!is.null(x$pvalue)) {
        m[b_idx, "pvalue"] <- ""
      }
      if (is.null(x$label)) {
        x$label <- rep("", length(x$lhs))
      }
      x$label[b_idx] <- ifelse(nchar(x$label[b_idx]) > 0L,
        paste(x$label[b_idx], "+ub", sep = ""),
        "ub"
      )
    }
    # remove upper column
    m <- m[, colnames(m) != "upper"]
  }


  # handle fmi
  if (!is.null(x$fmi)) {
    se_idx <- which(x$se == 0)
    if (length(se_idx) > 0L) {
      m[se_idx, "fmi"] <- ""
      ## for lavaan.mi-class objects (semTools)
      if (!is.null(x$riv)) m[se_idx, "riv"] <- ""
    }

    not_idx <- which(x$op %in% c(":=", "<", ">", "=="))
    if (length(not_idx) > 0L) {
      if (!is.null(x$fmi)) {
        m[not_idx, "fmi"] <- ""
        ## for lavaan.mi-class objects (semTools)
        if (!is.null(x$riv)) m[not_idx, "riv"] <- ""
      }
    }
  }

  # for blavaan, handle Post.SD and PSRF
  if (!is.null(x$Post.SD)) {
    se_idx <- which(x$Post.SD == 0)
    if (length(se_idx) > 0L) {
      m[se_idx, "Post.SD"] <- ""
      if (!is.null(x$psrf)) {
        m[se_idx, "psrf"] <- ""
      }
      if (!is.null(x$PSRF)) {
        m[se_idx, "PSRF"] <- ""
      }
    }

    # handle psrf for defined parameters
    not_idx <- which(x$op %in% c(":=", "<", ">", "=="))
    if (length(not_idx) > 0L) {
      if (!is.null(x$psrf)) {
        m[not_idx, "psrf"] <- ""
      }
      if (!is.null(x$PSRF)) {
        m[not_idx, "PSRF"] <- ""
      }
    }
  }

  # rename some column names
  colnames(m)[colnames(m) == "lhs"] <- ""
  colnames(m)[colnames(m) == "op"] <- ""
  colnames(m)[colnames(m) == "rhs"] <- ""
  colnames(m)[colnames(m) == "step"] <- "Step"
  colnames(m)[colnames(m) == "est"] <- "Estimate"
  colnames(m)[colnames(m) == "se"] <- "Std.Err"
  colnames(m)[colnames(m) == "z"] <- "z-value"
  colnames(m)[colnames(m) == "pvalue"] <- "P(>|z|)"
  colnames(m)[colnames(m) == "std.lv"] <- "Std.lv"
  colnames(m)[colnames(m) == "std.all"] <- "Std.all"
  colnames(m)[colnames(m) == "std.nox"] <- "Std.nox"
  colnames(m)[colnames(m) == "std.user"] <- "Std.usr"
  colnames(m)[colnames(m) == "prior"] <- "Prior"
  colnames(m)[colnames(m) == "fmi"] <- "FMI"
  ## for lavaan.mi-class objects (semTools)
  if ("t" %in% colnames(m)) {
    colnames(m)[colnames(m) == "t"] <- "t-value"
    colnames(m)[colnames(m) == "P(>|z|)"] <- "P(>|t|)"
    colnames(m)[colnames(m) == "riv"] <- "RIV"
  }

  # format column names
  colnames(m) <- sprintf(char_format, colnames(m))

  # exceptions for blavaan: Post.Mean (width = 9), Prior (width = 14)
  if (!is.null(x$Post.Mean)) {
    tmp <- gsub("[ \t]+", "", colnames(m), perl = TRUE)

    # reformat "Post.Mean" column
    col_idx <- which(tmp == "Post.Mean")
    if (length(col_idx) > 0L) {
      tmp_format <- paste("%", max(9, nd + 5), "s", sep = "")
      colnames(m)[col_idx] <- sprintf(tmp_format, colnames(m)[col_idx])
      m[, col_idx] <- sprintf(tmp_format, m[, col_idx])
    }

    # reformat "Prior" column
    col_idx <- which(tmp == "Prior")
    if (length(col_idx) > 0L) {
      max_1 <- max(nchar(m[, col_idx])) + 1L
      tmp_format <- paste("%", max(max_1, nd + 5), "s", sep = "")
      colnames(m)[col_idx] <- sprintf(tmp_format, colnames(m)[col_idx])
      m[, col_idx] <- sprintf(tmp_format, m[, col_idx])
    }
  }

  b <- 0L
  # group-specific sections
  for (g in 1:ngroups) {
    # group header
    if (ngroups > 1L) {
      group_label <- attr(x, "group.label")
      cat("\n\n")
      cat("Group ", g, " [", group_label[g], "]:\n", sep = "")
    }

    for (l in 1:nlevels) {
      # block number
      b <- b + 1L

      # ov/lv names
      ov_names <- lav_object_vnames(x, "ov", block = b)
      # lv_names <- lav_object_vnames(x, "lv", block = b)

      # level header
      if (nlevels > 1L) {
        level_label <- attr(x, "level.label")
        cat("\n\n")
        cat("Level ", l, " [", level_label[l], "]:\n", sep = "")
      }

      # group-specific sections
      for (s in gsections) {
        if (s == "Latent Variables") {
          row_idx <- which(x$op == "=~" & !x$lhs %in% ov_names &
            x$block == b)
          if (length(row_idx) == 0L) next
          m[row_idx, 1] <- lav_print_format_names(x$rhs[row_idx],
                                                 x$label[row_idx])
        } else if (s == "Composites") {
          row_idx <- which(x$op == "<~" & x$block == b)
          if (length(row_idx) == 0L) next
          m[row_idx, 1] <- lav_print_format_names(x$rhs[row_idx],
                                                 x$label[row_idx])
        } else if (s == "Regressions") {
          row_idx <- which(x$op == "~" & x$block == b)
          if (length(row_idx) == 0L) next
          m[row_idx, 1] <- lav_print_format_names(x$rhs[row_idx],
                                                  x$label[row_idx])
        } else if (s == "Covariances") {
          row_idx <- which(x$op == "~~" & x$lhs != x$rhs & !x$exo &
            x$block == b)
          if (length(row_idx) == 0L) next
          # make distinction between residual and plain
          y_names <- unique(c(
            lav_object_vnames(x, "eqs.y"),
            lav_object_vnames(x, "ov.ind"),
            lav_object_vnames(x, "lv.ind")
          ))
          prefix <- rep("", length(row_idx))
          prefix[x$rhs[row_idx] %in% y_names] <- "  ."
          m[row_idx, 1] <- lav_print_format_names(x$rhs[row_idx],
                                      x$label[row_idx], prefix = prefix)
          # m[row.idx,1] <- lav_print_format_names(x$rhs[row.idx],
          #                                            x$label[row.idx])
        } else if (s == "Intercepts") {
          row_idx <- which(x$op == "~1" & !x$exo & x$block == b)
          if (length(row_idx) == 0L) next
          # make distinction between intercepts and means
          y_names <- unique(c(
            lav_object_vnames(x, "eqs.y"),
            lav_object_vnames(x, "ov.ind"),
            lav_object_vnames(x, "lv.ind")
          ))
          prefix <- rep("", length(row_idx))
          prefix[x$lhs[row_idx] %in% y_names] <- "  ."
          m[row_idx, 1] <- lav_print_format_names(x$lhs[row_idx],
                           x$label[row_idx], prefix = prefix)
          # m[row.idx,1] <- lav_print_format_names(x$lhs[row.idx],
          #                                             x$label[row.idx])
        } else if (s == "Thresholds") {
          row_idx <- which(x$op == "|" & x$block == b)
          if (length(row_idx) == 0L) next
          m[row_idx, 1] <- lav_print_format_names(paste(x$lhs[row_idx], "|",
            x$rhs[row_idx],
            sep = ""
          ), x$label[row_idx])
        } else if (s == "Variances") {
          row_idx <- which(x$op == "~~" & x$lhs == x$rhs & !x$exo &
            x$block == b)
          if (length(row_idx) == 0L) next
          # make distinction between residual and plain
          y_names <- unique(c(
            lav_object_vnames(x, "eqs.y"),
            lav_object_vnames(x, "ov.ind"),
            lav_object_vnames(x, "lv.ind")
          ))
          prefix <- rep("", length(row_idx))
          prefix[x$rhs[row_idx] %in% y_names] <- "  ."
          m[row_idx, 1] <- lav_print_format_names(x$rhs[row_idx],
                                 x$label[row_idx], prefix = prefix)
        } else if (s == "Scales y*") {
          row_idx <- which(x$op == "~*~" & x$block == b)
          if (length(row_idx) == 0L) next
          m[row_idx, 1] <- lav_print_format_names(x$rhs[row_idx],
                                                     x$label[row_idx])
        } else if (s == "Group Weight") {
          row_idx <- which(x$lhs == "group" & x$op == "%" & x$block == b)
          if (length(row_idx) == 0L) next
          m[row_idx, 1] <- lav_print_format_names(x$rhs[row_idx],
                                                     x$label[row_idx])
        } else if (s == "R-Square") {
          row_idx <- which(x$op == "r2" & x$block == b)
          if (length(row_idx) == 0L) next
          m[row_idx, 1] <- lav_print_format_names(x$rhs[row_idx],
                                                     x$label[row_idx])
        } else {
          row_idx <- integer(0L)
        }

        # do we need special formatting for this section?
        # three types:
        #  - regular (nothing to do, except row/colnames)
        #  - R-square
        #  - Latent Variables (and Composites), Regressions and Covariances
        #    'bundle' the output per lhs element

        # bundling
        if (s %in% c(
          "Latent Variables", "Composites",
          "Regressions", "Covariances"
        )) {
          nel <- length(row_idx)
          m_1 <- matrix("", nrow = nel * 2, ncol = ncol(m))
          colnames(m_1) <- colnames(m)
          rownames(m_1) <- rep("", NROW(m_1))
          # colnames(M)[1] <- sprintf("%-17s", paste(s, ":", sep = ""))
          if (is.null(x$efa)) {
            lhs <- paste(x$lhs[row_idx], x$op[row_idx])
          } else {
            lhs <- paste(
              x$lhs[row_idx], x$op[row_idx],
              x$efa[row_idx]
            )
          }
          lhs_idx <- seq(1, nel * 2L, 2L)
          rhs_idx <- seq(1, nel * 2L, 2L) + 1L
          if (s == "Covariances") {
            # make distinction between residual and plain
            y_names <- unique(c(
              lav_object_vnames(x, "eqs.y"),
              lav_object_vnames(x, "ov.ind"),
              lav_object_vnames(x, "lv.ind")
            ))
            prefix <- rep("", length(row_idx))
            prefix[x$lhs[row_idx] %in% y_names] <- "."
          } else {
            prefix <- rep("", length(lhs))
          }
          m_1[lhs_idx, 1] <- sprintf("%1s%-15s", prefix, lhs)
          m_1[rhs_idx, ] <- m[row_idx, ]
          # avoid duplicated LHS labels
          if (nel > 1L) {
            del_idx <- integer(0)
            old_lhs <- ""
            for (i in 1:nel) {
              if (lhs[i] == old_lhs) {
                del_idx <- c(del_idx, lhs_idx[i])
              }
              old_lhs <- lhs[i]
            }
            if (length(del_idx) > 0L) {
              m_1 <- m_1[-del_idx, , drop = FALSE]
            }
          }
          cat("\n", s, ":\n", sep = "")
          # cat("\n")
          print(m_1, quote = FALSE)

          # R-square
        } else if (s == "R-Square") {
          m_1 <- m[row_idx, 1:2, drop = FALSE]
          colnames(m_1) <- colnames(m)[1:2]
          rownames(m_1) <- rep("", NROW(m_1))
          # colnames(M)[1] <- sprintf("%-17s", paste(s, ":", sep = ""))
          cat("\n", s, ":\n", sep = "")
          # cat("\n")
          print(m_1, quote = FALSE)

          # Regular
        } else {
          # M <- rbind(matrix("", nrow = 1L, ncol = ncol(m)),
          #           m[row.idx,])
          m_1 <- m[row_idx, , drop = FALSE]
          colnames(m_1) <- colnames(m)
          rownames(m_1) <- rep("", NROW(m_1))
          # colnames(M)[1] <- sprintf("%-17s", paste(s, ":", sep = ""))
          cat("\n", s, ":\n", sep = "")
          # cat("\n")
          print(m_1, quote = FALSE)
        }
      } # GSECTIONS
    } # levels
  } # groups

  # asections
  for (s in asections) {
    if (s == "Defined Parameters") {
      row_idx <- which(x$op == ":=")
      m[row_idx, 1] <- lav_print_format_names(x$lhs[row_idx], "")
      m_1 <- m[row_idx, , drop = FALSE]
      colnames(m_1) <- colnames(m)
    } else if (s == "Constraints") {
      row_idx <- which(x$op %in% c("==", "<", ">"))
      if (length(row_idx) == 0) next
      m[row_idx, 1] <- lav_print_format_con(x$lhs[row_idx],
                          x$op[row_idx], x$rhs[row_idx], nd = nd)
      m[row_idx, 2] <- sprintf(num_format, abs(x$est[row_idx]))
      m_1 <- m[row_idx, 1:2, drop = FALSE]
      colnames(m_1) <- c("", sprintf(char_format, "|Slack|"))
    } else {
      row_idx <- integer(0L)
    }

    if (length(row_idx) == 0L) {
      next
    }

    rownames(m_1) <- rep("", NROW(m_1))
    # colnames(M)[1] <- sprintf("%-17s", paste(s, ":", sep = ""))
    # cat("\n")
    cat("\n", s, ":\n", sep = "")
    print(m_1, quote = FALSE)
  }
  cat("\n")

  invisible(m)
}

lav_print_format_names <- function(names_1, labels_1, prefix = NULL) {
  m_w <- 14
  if (is.null(prefix)) {
    prefix <- rep("", length(names_1))
  }

  multi_b <- FALSE
  if (any(nchar(names_1) != nchar(names_1, "bytes"))) {
    multi_b <- TRUE
  }
  if (any(nchar(labels_1) != nchar(labels_1, "bytes"))) {
    multi_b <- TRUE
  }
  # labels?
  l_idx <- which(nchar(labels_1) > 0L)
  if (length(l_idx) > 0L) {
    if (!multi_b) {
      labels_1 <- abbreviate(labels_1, 4)
      labels_1[l_idx] <- paste(" (", labels_1[l_idx], ")", sep = "")
      max_l <- max(nchar(labels_1))
      names_1 <- abbreviate(names_1,
        minlength = (m_w - max_l),
        strict = TRUE
      )
    } else {
      # do not abbreviate anything (eg in multi-byte locales)
      max_l <- 4L
    }
    names_1 <- sprintf(paste("%-", (m_w - max_l), "s%", max_l, "s",
      sep = ""
    ), names_1, labels_1)
  } else {
    if (!multi_b) {
      names_1 <- abbreviate(names_1, minlength = m_w, strict = TRUE)
    } else {
      names_1 <- sprintf(paste("%-", m_w, "s", sep = ""), names_1)
    }
  }

  char_format <- paste("%3s%-", m_w, "s", sep = "")
  sprintf(char_format, prefix, names_1)
}

lav_print_format_con <- function(lhs, op, rhs, nd) {
  nd <- max(nd, 3)
  m_w <- 41 + (nd - 3) * 3

  nel <- length(lhs)
  if (length(nel) == 0L) {
    return(character(0))
  }
  names_1 <- character(nel)
  for (i in 1:nel) {
    if (rhs[i] == "0" && op[i] == ">") {
      con_string <- paste(lhs[i], " - 0", sep = "")
    } else if (rhs[i] == "0" && op[i] == "<") {
      con_string <- paste(rhs[i], " - (", lhs[i], ")", sep = "")
    } else if (rhs[i] != "0" && op[i] == ">") {
      con_string <- paste(lhs[i], " - (", rhs[i], ")", sep = "")
    } else if (rhs[i] != "0" && op[i] == "<") {
      con_string <- paste(rhs[i], " - (", lhs[i], ")", sep = "")
    } else if (rhs[i] == "0" && op[i] == "==") {
      con_string <- paste(lhs[i], " - 0", sep = "")
    } else if (rhs[i] != "0" && op[i] == "==") {
      con_string <- paste(lhs[i], " - (", rhs[i], ")", sep = "")
    }
    con_string <- abbreviate(con_string, m_w, strict = TRUE)
    char_format <- paste("   %-", m_w, "s", sep = "")
    names_1[i] <- sprintf(char_format, con_string)
  }

  names_1
}

# new in 0.6-12
lav_summary_print <- function(x, ..., nd = 3L) {
  y <- unclass(x) # change to ordinary list

  # get nd, if it is stored as an attribute
  nd_1 <- attr(y, "nd")
  if (!is.null(nd_1) && is.numeric(nd_1)) {
    nd <- as.integer(nd_1)
  }

  # header
  if (!is.null(y$header)) {
    lavaan_version <- y$header$lavaan.version
    sam_approach <- y$header$sam.approach
    optim_method <- y$header$optim.method
    optim_iterations <- y$header$optim.iterations
    optim_converged <- y$header$optim.converged

    # sam or sem?
    if (sam_approach) {
      cat("This is ",
        sprintf("lavaan %s", lavaan_version),
        " -- using the SAM approach to SEM\n",
        sep = ""
      )
    } else {
      cat(sprintf("lavaan %s ", lavaan_version))

      # Convergence or not?
      if (optim_method == "none") {
        cat("-- DRY RUN with 0 iterations --\n")
      } else if (optim_iterations > 0) {
        if (optim_converged) {
          if (optim_iterations == 1L) {
            cat("ended normally after 1 iteration\n")
          } else {
            cat(sprintf(
              "ended normally after %i iterations\n",
              optim_iterations
            ))
          }
        } else {
          if (optim_iterations == 1L) {
            cat("did NOT end normally after 1 iteration\n")
          } else {
            cat(sprintf(
              "did NOT end normally after %i iterations\n",
              optim_iterations
            ))
          }
          cat("** WARNING ** Estimates below are most likely unreliable\n")
        }
      } else {
        cat("did not run (perhaps do.fit = FALSE)?\n")
        cat("** WARNING ** Estimates below are simply the starting values\n")
      }
    }
  }

  #TDJ: print header for lavaan.mi object (or nothing when NULL)
  cat(y$top_of_lavaanmi)

  # optim
  if (!is.null(y$optim)) {
    estimator <- y$optim$estimator
    estimator_args <- y$optim$estimator.args
    optim_method <- y$optim$optim.method
    npar <- y$optim$npar
    eq_constraints <- y$optim$eq.constraints
    nrow_ceq_jac <- y$optim$nrow.ceq.jac
    nrow_cin_jac <- y$optim$nrow.cin.jac
    nrow_con_jac <- y$optim$nrow.con.jac
    con_jac_rank <- y$optim$con.jac.rank

    cat("\n")
    # cat("Optimization information:\n\n")

    c1 <- c("Estimator")
    # second column
    tmp_est <- toupper(estimator)
    if (tmp_est == "DLS") {
      dls_first_letter <- substr(
        estimator_args$dls.GammaNT,
        1L, 1L
      )
      tmp_est <- paste("DLS-", toupper(dls_first_letter), sep = "")
    }
    # reduced-bias M-estimation: show IRBM / ERBM (estimator is ML internally)
    if (!is.null(estimator_args$rbm.method)) {
      tmp_est <- switch(estimator_args$rbm.method,
        implicit = "IRB-ML",
        explicit = "ERB-ML",
        none = "ML",
        "RBM"
      )
    }
    c2 <- tmp_est

    # additional estimator args
    if (!is.null(estimator_args) &&
      length(estimator_args) > 0L) {
      if (estimator == "DLS") {
        c1 <- c(c1, "Estimator DLS value for a")
        c2 <- c(c2, estimator_args$dls.a)
      }
    }

    # optimization method + npar
    c1 <- c(c1, "Optimization method", "Number of model parameters")
    c2 <- c(c2, toupper(optim_method), npar)

    # optional output
    if (eq_constraints) {
      c1 <- c(c1, "Number of equality constraints")
      c2 <- c(c2, nrow_ceq_jac)
    }
    if (nrow_cin_jac > 0L) {
      c1 <- c(c1, "Number of inequality constraints")
      c2 <- c(c2, nrow_cin_jac)
    }
    if (nrow_con_jac > 0L) {
      if (con_jac_rank == (nrow_ceq_jac + nrow_cin_jac)) {
        # nothing to do (don't print, as this is redundant information)
      } else {
        c1 <- c(c1, "Row rank of the constraints matrix")
        c2 <- c(c2, con_jac_rank)
      }
    }

    # format
    c1 <- format(c1, width = 40L)
    c2 <- format(c2,
      width = 11L + max(0, (nd - 3L)) * 4L,
      justify = "right"
    )

    # character matrix
    m <- cbind(c1, c2, deparse.level = 0)
    colnames(m) <- rep("", ncol(m))
    rownames(m) <- rep(" ", nrow(m))

    # print
    write.table(m, row.names = TRUE, col.names = FALSE, quote = FALSE)
  }

  # sam header
  if (!is.null(y$sam.header)) {
    cat("\n")
    sam_method <- y$sam.header$sam.method
    sam_local_options <- y$sam.header$sam.local.options
    sam_mm_list <- y$sam.header$sam.mm.list
    sam_mm_estimator <- y$sam.header$sam.mm.estimator
    sam_struc_estimator <- y$sam.header$sam.struc.estimator

    # sam method
    c1 <- c("SAM method")
    c2 <- toupper(sam_method)

    # options
    if (sam_method == "local") {
      c1 <- c(c1, "Mapping matrix M method")
      c2 <- c(c2, sam_local_options$M.method)
      # TODo: more!
    }

    # number of measurement blocks
    c1 <- c(c1, "Number of measurement blocks")
    c2 <- c(c2, length(sam_mm_list))

    # estimator measurement blocks
    c1 <- c(c1, "Estimator measurement part")
    c2 <- c(c2, sam_mm_estimator)

    # estimator structural part
    c1 <- c(c1, "Estimator  structural part")
    c2 <- c(c2, sam_struc_estimator)

    # format
    c1 <- format(c1, width = 40L)
    c2 <- format(c2,
      width = 11L + max(0, (nd - 3L)) * 4L,
      justify = "right"
    )

    # character matrix
    m <- cbind(c1, c2, deparse.level = 0)
    colnames(m) <- rep("", ncol(m))
    rownames(m) <- rep(" ", nrow(m))

    # print
    write.table(m, row.names = TRUE, col.names = FALSE, quote = FALSE)
  }

  # efa/rotation
  if (!is.null(y$rotation)) {
    cat("\n")
    rotation <- y$rotation
    # rotation_args <- y$rotation.args

    # cat("Rotation information:\n\n")
    # container
    c1 <- c2 <- character(0L)

    # rotation method
    c1 <- c(c1, "Rotation method")
    if (rotation$rotation == "none") {
      mm <- toupper(rotation$rotation)
    } else if (rotation$rotation.args$orthogonal) {
      mm <- paste(toupper(rotation$rotation), " ", "ORTHOGONAL",
        sep = ""
      )
    } else {
      mm <- paste(toupper(rotation$rotation), " ", "OBLIQUE",
        sep = ""
      )
    }
    c2 <- c(c2, mm)

    if (rotation$rotation != "none") {
      # method options
      if (rotation$rotation == "geomin") {
        c1 <- c(c1, "Geomin epsilon")
        c2 <- c(c2, rotation$rotation.args$geomin_epsilon)
      } else if (rotation$rotation == "orthomax") {
        c1 <- c(c1, "Orthomax gamma")
        c2 <- c(c2, rotation$rotation.args$orthomax_gamma)
      } else if (rotation$rotation == "cf") {
        c1 <- c(c1, "Crawford-Ferguson gamma")
        c2 <- c(c2, rotation$rotation.args$cf_gamma)
      } else if (rotation$rotation == "oblimin") {
        c1 <- c(c1, "Oblimin gamma")
        c2 <- c(c2, rotation$rotation.args$oblimin_gamma)
      } else if (rotation$rotation == "promax") {
        c1 <- c(c1, "Promax kappa")
        c2 <- c(c2, rotation$rotation.args$promax_kappa)
      }

      # rotation algorithm
      c1 <- c(c1, "Rotation algorithm (rstarts)")
      tmp <- paste(toupper(rotation$rotation.args$algorithm),
        " (", rotation$rotation.args$rstarts, ")",
        sep = ""
      )
      c2 <- c(c2, tmp)

      # Standardized metric (or not)
      c1 <- c(c1, "Standardized metric")
      if (rotation$rotation.args$std_ov) {
        c2 <- c(c2, "TRUE")
      } else {
        c2 <- c(c2, "FALSE")
      }

      # Row weights
      c1 <- c(c1, "Row weights")
      tmp_txt <- rotation$rotation.args$row_weights
      c2 <- c(c2, paste(toupper(substring(tmp_txt, 1, 1)),
        substring(tmp_txt, 2),
        sep = ""
      ))
    }

    # format c1/c2
    c1 <- format(c1, width = 33L)
    c2 <- format(c2,
      width = 18L + max(0, (nd - 3L)) * 4L,
      justify = "right"
    )

    # create character matrix
    m <- cbind(c1, c2, deparse.level = 0)
    colnames(m) <- rep("", ncol(m))
    rownames(m) <- rep(" ", nrow(m))

    # print
    write.table(m, row.names = TRUE, col.names = FALSE, quote = FALSE)
  }

  # data object
  if (!is.null(y$data)) {
    cat("\n")
    lav_data_print_short(y$data, nd = nd)
  }

  # sam local stats: measurement blocks + structural part
  if (!is.null(y$sam)) {
    cat("\n")
    sam_method <- y$sam$sam.method
    sam_mm_table <- y$sam$sam.mm.table
    sam_mm_rel <- y$sam$sam.mm.rel
    sam_struc_fit <- y$sam$sam.struc.fit
    ngroups <- y$sam$ngroups
    nlevels <- y$sam$nlevels
    group_label <- y$sam$group.label
    level_label <- y$sam$level.label
    block_label <- y$sam$block.label

    # measurement
    tmp <- sam_mm_table
    if (sam_method == "global") {
      cat("Summary Information Measurement Part:\n\n")
    } else {
      cat("Summary Information Measurement + Structural:\n\n")
    }
    print(tmp, row.names = rep(" ", nrow(tmp)), nd = nd)

    if (sam_method == "local") {
      # reliability information
      c1 <- c2 <- character(0L)
      if (ngroups == 1L && nlevels == 1L) {
        cat("\n")
        cat("  Model-based reliability latent variables:\n\n")
        tmp <- data.frame(as.list(sam_mm_rel[[1]]))
        class(tmp) <- c("lavaan.data.frame", "data.frame")
        print(tmp, row.names = rep(" ", nrow(tmp)), nd = nd)
      } else if (ngroups > 1L && nlevels == 1L) {
        cat("\n")
        cat("  Model-based reliability latent variables (per group):\n")
        for (g in 1:ngroups) {
          cat("\n")
          cat("  Group ", g, " [", group_label[g], "]:\n\n",
            sep = ""
          )
          tmp <- data.frame(as.list(sam_mm_rel[[g]]))
          class(tmp) <- c("lavaan.data.frame", "data.frame")
          print(tmp, row.names = rep(" ", nrow(tmp)), nd = nd)
        }
      } else if (ngroups == 1L && nlevels > 1L) {
        cat("\n")
        cat("  Model-based reliability latent variables (per level):\n")
        for (g in 1:nlevels) {
          cat("\n")
          cat("  Level ", g, " [", level_label[g], "]:\n\n",
            sep = ""
          )
          tmp <- data.frame(as.list(sam_mm_rel[[g]]))
          class(tmp) <- c("lavaan.data.frame", "data.frame")
          print(tmp, row.names = rep(" ", nrow(tmp)), nd = nd)
        }
      } else if (ngroups > 1L && nlevels > 1L) {
        cat("\n")
        cat("  Model-based reliability latent variables (per group/level):\n")
        for (g in seq_along(block_label)) {
          cat("\n")
          cat("  Group/Level ", g, " [", block_label[g], "]:\n\n",
            sep = ""
          )
          tmp <- data.frame(as.list(sam_mm_rel[[g]]))
          class(tmp) <- c("lavaan.data.frame", "data.frame")
          print(tmp, row.names = rep(" ", nrow(tmp)), nd = nd)
        }
      }

      cat("\n")
      cat("  Summary Information Structural part:\n\n")
      tmp <- data.frame(as.list(sam_struc_fit))
      class(tmp) <- c("lavaan.data.frame", "data.frame")
      print(tmp, row.names = rep(" ", nrow(tmp)), nd = nd)
    }
  }

  # test statistics
  if (!is.null(y$test)) {
    cat("\n")
    lav_test_print(y$test, nd = nd)
  }

  # extra fit measures (if present)
  if (!is.null(y$fit)) {
    add_h0 <- FALSE
    if (!is.null(attr(y$fit, "add.h0"))) {
      add_h0 <- isTRUE(attr(y$fit, "add.h0"))
    }
    lav_fitmeasures_print(y$fit, nd = nd, add.h0 = add_h0)
  }

  # efa output
  if (!is.null(y$efa)) {
    # get cutoff, if it is stored as an attribute
    ct <- attr(y, "cutoff")
    if (!is.null(ct) && is.numeric(ct)) {
      cutoff <- ct
    } else {
      cutoff <- 0.3
    }
    # get dot.cutoff, if it is stored as an attribute
    dc <- attr(y, "dot.cutoff")
    if (!is.null(dc) && is.numeric(dc)) {
      dot_cutoff <- dc
    } else {
      dot_cutoff <- 0.1
    }
    # get alpha.level, if it is stored as an attribute
    al <- attr(y, "alpha.level")
    if (!is.null(al) && is.numeric(al)) {
      alpha_level <- al
    } else {
      alpha_level <- 0.01
    }

    for (b in seq_len(y$efa$nblocks)) {
      if (length(y$efa$block.label) > 0L) {
        cat("\n\n")
        cat(y$efa$block.label[[b]], ":\n\n", sep = "")
      }
      if (!is.null(y$efa$lambda[[b]])) {
        cat("\n")
        if (!is.null(y$efa$lambda.se[[b]]) && alpha_level > 0) {
          cat("Standardized loadings: (* = significant at ",
            round(alpha_level * 100),
            "% level)\n\n",
            sep = ""
          )
        } else {
          cat("Standardized loadings:\n\n")
        }
        mm_lambda <- unclass(y$efa$lambda[[b]])
        mm_theta <- unname(unclass(y$efa$theta[[b]]))
        lav_print_loadings(mm_lambda,
          nd = nd, cutoff = cutoff,
          dot_cutoff = dot_cutoff,
          alpha_level = alpha_level,
          resvar = mm_theta, # diag elements only
          x_se = y$efa$lambda.se[[b]]
        )
      }

      if (!is.null(y$efa$sumsq_table[[b]])) {
        cat("\n")
        print(y$efa$sumsq_table[[b]], nd = nd)
      }

      # factor correlations:
      if (!y$efa$orthogonal && !is.null(y$efa$psi[[b]]) &&
        ncol(y$efa$psi[[b]]) > 1L) {
        cat("\n")
        if (!is.null(y$efa$psi.se[[b]]) && alpha_level > 0) {
          cat("Factor correlations: (* = significant at ",
            round(alpha_level * 100),
            "% level)\n\n",
            sep = ""
          )
        } else {
          cat("Factor correlations:\n\n")
        }
        lav_print_psi(y$efa$psi[[b]],
          nd = nd,
          alpha_level = alpha_level,
          x_se = y$efa$psi.se[[b]]
        )
      }

      # factor score determinacy (for regression scores only!)
      if (!is.null(y$efa$fs_determinacy[[b]])) {
        cat("
Correlation regression factor scores and factors (determinacy):\n\n")
        print(y$efa$fs_determinacy[[b]], nd = nd)
        cat("
R2 regression factor scores (= squared correlations):\n\n")
        tmp <- y$efa$fs_determinacy[[b]]
        tmp2 <- tmp * tmp
        class(tmp2) <- c("lavaan.vector", "numeric")
        print(tmp2, nd = nd)
      }

      # lambda_structure
      if (!is.null(y$efa$lambda_structure[[b]])) {
        cat("\n")
        cat("Standardized structure (= LAMBDA %*% PSI):\n\n")
        print(y$efa$lambda_structure[[b]], nd = nd)
      }

      # standard errors lambda
      if (!is.null(y$efa$theta.se[[b]])) { # we check for theta.se
        # as lambda.se is needed for '*'
        cat("\n")
        cat("Standard errors standardized loadings:\n\n")
        print(y$efa$lambda.se[[b]], nd = nd)
      }

      # z-statistics lambda
      if (!is.null(y$efa$lambda.zstat[[b]])) {
        cat("\n")
        cat("Z-statistics standardized loadings:\n\n")
        print(y$efa$lambda.zstat[[b]], nd = nd)
      }

      # pvalues lambda
      if (!is.null(y$efa$lambda.pvalue[[b]])) {
        cat("\n")
        cat("P-values standardized loadings:\n\n")
        print(y$efa$lambda.pvalue[[b]], nd = nd)
      }

      # standard errors theta
      if (!is.null(y$efa$theta.se[[b]])) {
        cat("\n")
        cat("Standard errors unique variances:\n\n")
        print(y$efa$theta.se[[b]], nd = nd)
      }

      # z-statistics theta
      if (!is.null(y$efa$theta.zstat[[b]])) {
        cat("\n")
        cat("Z-statistics unique variances:\n\n")
        print(y$efa$theta.zstat[[b]], nd = nd)
      }

      # pvalues theta
      if (!is.null(y$efa$theta.pvalue[[b]])) {
        cat("\n")
        cat("P-values unique variances:\n\n")
        print(y$efa$theta.pvalue[[b]], nd = nd)
      }

      # standard errors psi
      if (!is.null(y$efa$theta.se[[b]])) { # we check for theta.se
        # as psi.se is needed for '*'
        cat("\n")
        cat("Standard errors factor correlations:\n\n")
        print(y$efa$psi.se[[b]], nd = nd)
      }

      # z-statistics psi
      if (!is.null(y$efa$psi.zstat[[b]])) {
        cat("\n")
        cat("Z-statistics factor correlations:\n\n")
        print(y$efa$psi.zstat[[b]], nd = nd)
      }

      # pvalues psi
      if (!is.null(y$efa$psi.pvalue[[b]])) {
        cat("\n")
        cat("P-values factor correlations:\n\n")
        print(y$efa$psi.pvalue[[b]], nd = nd)
      }
    } # blocks
    cat("\n")
  } # efa

  # residuals (printed in the lavResiduals(output = "text") style), after the
  # fit measures and just before the parameter table
  if (!is.null(y$residuals)) {
    res_text <- attr(y$residuals, "text")
    if (!is.null(res_text)) {
      lav_residuals_summary_print(res_text, nd = nd)
    }
  }

  # parameter table
  if (!is.null(y$pe) && is.null(y$efa)) {
    pe <- y$pe
    class(pe) <- c(
      "lavaan.parameterEstimates", "lavaan.data.frame",
      "data.frame"
    )
    print(pe, nd = nd)

    # note for nonlinear defined (:=) parameters under the delta method
    if (isTRUE(attr(pe, "delta.note"))) {
      lav_msg_note(gettext(
        "Standard errors and confidence intervals of the (nonlinear) defined
         (:=) parameters are based on the first-order delta method; for
         strongly nonlinear definitions, se.def = \"mc\" (Monte Carlo) or
         se = \"bootstrap\" may be more accurate."))
    }
  }

  # modification indices
  if (!is.null(y$mi)) {
    cat("Modification Indices:\n\n")
    mi <- y$mi
    rownames(mi) <- NULL
    print(mi, nd = nd)
  }

  invisible(y)
}

# helper function to print the loading matrix, masking small loadings
lav_print_loadings <- function(x, nd = 3L, cutoff = 0.3, dot_cutoff = 0.1,
                               alpha_level = 0.01, resvar = NULL, x_se = NULL) {
  # unclass
  y <- unclass(x)

  # round, and create a character matrix
  y <- format(round(y, nd), width = 3L + nd, justify = "right")

  # right-align column names
  colnames(y) <- format(colnames(y), width = 3L + nd, justify = "right")

  # create dot/empty string
  dot_string <- format(".", width = 3L + nd, justify = "right")
  empty_string <- format(" ", width = 3L + nd)

  # print a 'dot' if dot.cutoff < |loading| < cutoff
  if (dot_cutoff < cutoff) {
    y[abs(x) < cutoff & abs(x) > dot_cutoff] <- dot_string
  }

  # print nothing if |loading| < dot.cutoff
  y[abs(x) < min(dot_cutoff, cutoff)] <- empty_string

  # add 'star' for significant loadings (if provided) using alpha = 0.01
  if (!is.null(x_se) && !any(is.na(x_se))) {
    col_names <- colnames(y)
    row_names <- rownames(y)
    x_se[x_se < sqrt(.Machine$double.eps)] <- 1 # to avoid NA
    zstat <- x / x_se
    z_cutoff <- qnorm(1 - (alpha_level / 2))
    zstat_string <- ifelse(abs(zstat) > z_cutoff, "*", " ")
    y <- matrix(paste(y, zstat_string, sep = ""), nrow(y), ncol(y))
    colnames(y) <- col_names
    rownames(y) <- row_names
  }

  # add resvar
  if (!is.null(resvar)) {
    names_1 <- colnames(y)
    y <- cbind(y, format(round(cbind(resvar, 1 - resvar), nd),
      width = 12L + nd, justify = "right"
    ))
    resvar_names <- format(c("unique.var", "communalities"),
      width = 12L + nd, justify = "right"
    )
    colnames(y) <- c(names_1, resvar_names)
  }

  # print
  print(y, quote = FALSE)
}

# helper function to print the psi matrix, showing signif stars
lav_print_psi <- function(x, nd = 3L, alpha_level = 0.01, x_se = NULL) {
  # unclass
  y <- unclass(x)

  # round, and create a character matrix
  y <- format(round(y, nd), width = 3L + nd, justify = "right")

  # right-align column names
  colnames(y) <- format(colnames(y), width = 3L + nd, justify = "right")

  # add 'star' for significant loadings (if provided) using alpha = 0.01
  if (!is.null(x_se) && !any(is.na(x_se))) {
    col_names <- colnames(y)
    row_names <- rownames(y)
    x_se[x_se < sqrt(.Machine$double.eps)] <- 1 # to avoid NA
    zstat <- x / x_se
    z_cutoff <- qnorm(1 - (alpha_level / 2))
    zstat_string <- ifelse(abs(zstat) > z_cutoff, "*", " ")
    y <- matrix(paste(y, zstat_string, sep = ""), nrow(y), ncol(y))
    colnames(y) <- col_names
    rownames(y) <- row_names
  }

  # remove upper part
  ll <- upper.tri(x, diag = FALSE)
  y[ll] <- ""

  # print
  print(y, quote = FALSE)
}
