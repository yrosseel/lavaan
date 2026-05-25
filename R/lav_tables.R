# construct 1D, 2D or pattern-based frequency tables
# YR. 10 April 2013
# Notes:
# - we do NOT make a distinction here between unordered and ordered categorical
#   variables
# - object can be a matrix (most likely with integers), a full data frame,
#   a fitted lavaan object, or a lavData object
# - 11 May 2013: added collapse=TRUE, min.std.resid options (suggested
#   by Myrsini Katsikatsou
# - 11 June 2013: added dimension, to get one-way and two-way (three-way?)
#   tables
# - 20 Sept 2013: - allow for sample-based or model-based cell probabilities
#                   re-organize/re-name to provide a more consistent interface
#                   rows in the output can be either: cells, tables or patterns
#                 - dimension=0 equals type="pattern
#                 - collapse=TRUE is replaced by type="table"
#                 - changed names of statistics: std.resid is now GR.average
#                 - added many more statistics; some based on the model, some
#                   on the unrestricted model
# - 8 Nov 2013:   - skip empty cells for G2, instead of adding 0.5 to obs
# - 7 Feb 2016:   - take care of conditional.x = TRUE

lavTables <- function(object,                         # nolint start
                      # what type of table?
                      dimension = 2L,
                      type = "cells",
                      # if raw data, additional attributes
                      categorical = NULL,
                      group = NULL,
                      # which statistics / fit indices?
                      statistic = "default",
                      G2.min = 3.0, # needed for G2.{p/n}large
                      X2.min = 3.0, # needed for X2.{p/n}large
                      # pvalues for statistics?
                      p.value = FALSE,
                      # Bonferroni
                      # alpha.adj    = FALSE,
                      # output format
                      output = "data.frame",
                      patternAsString = TRUE) {      # nolint end
  # check object
  object <- lav_object_check_version(object)

  # check input
  if (!(dimension == 0L || dimension == 1L || dimension == 2L)) {
    lav_msg_stop(gettext(
      "dimension must be 0, 1 or 2 for pattern, one-way or two-way tables"))
  }
  stopifnot(type %in% c("cells", "table", "pattern"))
  if (type == "pattern") {
    dimension <- 0L
  }

  # extract or create lavdata
  lavdata <- lav_lavdata(object, ordered = categorical, group = group)

  # is 'object' a lavaan object?
  lavobject <- NULL
  if (inherits(object, "lavaan")) {
    lavobject <- object
  }

  # case 1: response patterns
  if (dimension == 0L) {
    out <- lav_tables_pattern(
      lavobject = lavobject, lavdata = lavdata,
      statistic = statistic,
      pattern_as_string = patternAsString
    )
    # output format
    if (output == "data.frame") {
      class(out) <- c("lavaan.data.frame", "data.frame")
    } else {
      lav_msg_warn(gettextf("output option `%s' is not available; ignored.",
                            output))
    }

    # case 2: one-way/univariate
  } else if (dimension == 1L) {
    out <- lav_tables_oneway(
      lavobject = lavobject, lavdata = lavdata,
      statistic = statistic
    )

    # output format
    if (output == "data.frame") {
      class(out) <- c("lavaan.data.frame", "data.frame")
    } else {
      lav_msg_warn(gettextf("output option `%s' is not available; ignored.",
                            output))
    }

    # case 3a: two-way/pairwise/bivariate + cells
  } else if (dimension == 2L && type == "cells") {
    out <- lav_tables_pairwise_cells(
      lavobject = lavobject,
      lavdata = lavdata,
      statistic = statistic
    )
    # output format
    if (output == "data.frame") {
      class(out) <- c("lavaan.data.frame", "data.frame")
    } else if (output == "table") {
      out <- lav_tables_cells_format(out, lavdata = lavdata)
    } else {
      lav_msg_warn(gettextf("output option `%s' is not available; ignored.",
                            output))
    }

    #  case 3b: two-way/pairwise/bivariate + collapsed table
  } else if (dimension == 2L && (type == "table" || type == "tables")) {
    out <- lav_tables_pairwise_table(
      lavobject = lavobject,
      lavdata = lavdata,
      statistic = statistic,
      g2_min = G2.min,
      x2_min = X2.min,
      p_value = p.value
    )
    # output format
    if (output == "data.frame") {
      class(out) <- c("lavaan.data.frame", "data.frame")
    } else if (output == "table") {
      out <- lav_tables_table_format(out,
        lavdata = lavdata,
        lavobject = lavobject
      )
    } else {
      lav_msg_warn(gettext("output option `%s' is not available; ignored.",
                           output))
    }
  }

  if ((is.data.frame(out) && nrow(out) == 0L) ||
    (is.list(out) && length(out) == 0L)) {
    # empty table (perhaps, no categorical variables)
    return(invisible(out))
  }

  out
}

# shortcut, always dim=2, type="cells"
# lavTablesFit <- function(object,
#                         # if raw data, additional attributes
#                         categorical  = NULL,
#                         group        = NULL,
#                         # which statistics / fit indices?
#                         statistic    = "default",
#                         G2.min       = 3.0,
#                         X2.min       = 3.0,
#                         # pvalues for statistics?
#                         p.value      = FALSE,
#                         # output format
#                         output       = "data.frame") {
#
#    lavTables(object = object, dimension = 2L, type = "table",
#              categorical = categorical, group = group,
#              statistic = statistic,
#              G2.min = G2.min, X2.min = X2.min, p.value = p.value,
#              output = output, patternAsString = FALSE)
# }

# lavTables1D <- function(object,
#                        # if raw data, additional attributes
#                        categorical  = NULL,
#                       group        = NULL,
#                       # which statistics / fit indices?
#                        statistic    = "default",
#                        # output format
#                        output = "data.frame") {
#
#   lavTables(object = object, dimension = 1L,
#              categorical = categorical, group = group,
#              statistic = statistic, p.value = FALSE,
#              output = output, patternAsString = FALSE)
# }


lav_tables_pattern <- function(lavobject = NULL, lavdata = NULL,
                               statistic = NULL, pattern_as_string = TRUE) {
  # this only works if we have 'categorical' variables
  cat_idx <- which(lavdata@ov$type %in% c("ordered", "factor"))
  if (length(cat_idx) == 0L) {
    lav_msg_warn(gettext("no categorical variables are found"))
    return(data.frame(
      pattern = character(0L), nobs = integer(0L),
      obs.freq = integer(0L), obs.prop = numeric(0L)
    ))
  }
  # no support yet for mixture of endogenous ordered + numeric variables
  if (!is.null(lavobject) &&
    length(lav_object_vnames(lavobject, "ov.nox")) > length(cat_idx)) {
    lav_msg_warn(gettext("some endogenous variables are not categorical"))
    return(data.frame(
      pattern = character(0L), nobs = integer(0L),
      obs.freq = integer(0L), obs.prop = numeric(0L)
    ))
  }

  # default statistics
  if (!is.null(lavobject)) {
    if (length(statistic) == 1L && statistic == "default") {
      statistic <- c("G2", "X2")
    } else {
      stopifnot(statistic %in% c("G2.un", "X2.un", "G2", "X2"))
    }
  } else {
    # only data
    if (length(statistic) == 1L && statistic == "default") {
      # if data, none by default
      statistic <- character(0L)
    } else {
      stopifnot(statistic %in% c("G2.un", "X2.un"))
    }
  }

  # first, create basic table with response patterns
  for (g in 1:lavdata@ngroups) {
    pat <- lav_data_resp_patterns(lavdata@X[[g]])$pat
    obs_freq <- as.integer(rownames(pat))
    if (pattern_as_string) {
      pat <- data.frame(
        pattern = apply(pat, 1, paste, collapse = ""),
        stringsAsFactors = FALSE
      )
    } else {
      pat <- as.data.frame(pat, stringsAsFactors = FALSE)
      names(pat) <- lavdata@ov.names[[g]]
    }
    # pat$id <- 1:nrow(pat)
    if (lavdata@ngroups > 1L) {
      pat$group <- rep(g, nrow(pat))
    }
    nobs_1 <- sum(obs_freq)
    pat$nobs <- rep(nobs_1, nrow(pat))
    pat$obs.freq <- obs_freq
    rownames(pat) <- NULL
    if (g == 1L) {
      out <- pat
    } else {
      out <- rbind(out, pat)
    }
  }

  out$obs.prop <- out$obs.freq / out$nobs

  if (any(c("X2.un", "G2.un") %in% statistic)) {
    # not a good statistic... we only have uni+bivariate information
    lav_msg_warn(gettext(
      "limited information used for thresholds and correlations;
      but X2/G2 assumes full information"))
    pi0 <- lav_tables_resp_pi(
      lavobject = lavobject, lavdata = lavdata,
      est = "h1"
    )

    out$est.prop.un <- unlist(pi0)
    if ("G2.un" %in% statistic) {
      out$G2.un <- lav_tables_stat_g2(
        out$obs.prop, out$est.prop.un,
        out$nobs
      )
    }
    if ("X2.un" %in% statistic) {
      out$X2.un <- lav_tables_stat_x2(
        out$obs.prop, out$est.prop.un,
        out$nobs
      )
    }
  }
  if (any(c("X2", "G2") %in% statistic)) {
    if (lavobject@Options$estimator %in% c("FML")) {
      # ok, nothing to say
    } else if (lavobject@Options$estimator %in%
      c("WLS", "DWLS", "PML", "ULS")) {
      lav_msg_warn(gettextf(
        "estimator %s is not using full information while est.prop is
        using full information", lavobject@Options$estimator))
    } else {
      lav_msg_stop(gettextf(
        "estimator %s is not supported.", lavobject@Options$estimator))
    }

    pi0 <- lav_tables_resp_pi(
      lavobject = lavobject, lavdata = lavdata,
      est = "h0"
    )

    out$est.prop <- unlist(pi0)
    if ("G2" %in% statistic) {
      out$G2 <- lav_tables_stat_g2(
        out$obs.prop, out$est.prop,
        out$nobs
      )
    }
    if ("X2" %in% statistic) {
      out$X2 <- lav_tables_stat_x2(
        out$obs.prop, out$est.prop,
        out$nobs
      )
    }
  }

  # remove nobs?
  # out$nobs <- NULL

  out
}

# pairwise tables, rows = table cells
lav_tables_pairwise_cells <- function(lavobject = NULL, lavdata = NULL,
                                      statistic = character(0L)) {
  # this only works if we have at least two 'categorical' variables
  cat_idx <- which(lavdata@ov$type %in% c("ordered", "factor"))
  if (length(cat_idx) == 0L) {
    lav_msg_warn(gettext("no categorical variables are found"))
    return(data.frame(
      id = integer(0L), lhs = character(0L), rhs = character(0L),
      nobs = integer(0L), row = integer(0L), col = integer(0L),
      obs.freq = integer(0L), obs.prop = numeric(0L)
    ))
  }
  if (length(cat_idx) == 1L) {
    lav_msg_warn(gettext("at least two categorical variables are needed"))
    return(data.frame(
      id = integer(0L), lhs = character(0L), rhs = character(0L),
      nobs = integer(0L), row = integer(0L), col = integer(0L),
      obs.freq = integer(0L), obs.prop = numeric(0L)
    ))
  }

  # default statistics
  if (!is.null(lavobject)) {
    if (length(statistic) == 1L && statistic == "default") {
      statistic <- c("X2")
    } else {
      stopifnot(statistic %in% c(
        "cor", "th", "X2", "G2",
        "cor.un", "th.un", "X2.un", "G2.un"
      ))
    }
  } else {
    if (length(statistic) == 1L && statistic == "default") {
      # if data, none by default
      statistic <- character(0L)
    } else {
      stopifnot(statistic %in% c("cor.un", "th.un", "X2.un", "G2.un"))
    }
  }

  # initial table, observed cell frequencies
  out <- lav_tables_pairwise_freq_cell(
    lavdata = lavdata,
    as_data_frame = TRUE
  )
  out$obs.prop <- out$obs.freq / out$nobs

  if (any(c("cor.un", "th.un", "X2.un", "G2.un") %in% statistic)) {
    pi0 <- lav_tables_pairwise_sample_pi(
      lavobject = lavobject,
      lavdata = lavdata
    )
    out$est.prop.un <- unlist(pi0)
    if ("G2.un" %in% statistic) {
      out$G2.un <- lav_tables_stat_g2(
        out$obs.prop, out$est.prop.un,
        out$nobs
      )
    }
    if ("X2.un" %in% statistic) {
      out$X2.un <- lav_tables_stat_x2(
        out$obs.prop, out$est.prop.un,
        out$nobs
      )
    }

    if ("cor.un" %in% statistic) {
      cor_1 <- attr(pi0, "COR")
      cor_all <- unlist(lapply(cor_1, function(x) {
        x[lower.tri(x, diag = FALSE)]
      }))
      out$cor.un <- cor_all[out$id]
    }
  }
  if (any(c("cor", "th", "X2", "G2") %in% statistic)) {
    pi0 <- lav_tables_pairwise_model_pi(lavobject = lavobject)
    out$est.prop <- unlist(pi0)
    if ("G2" %in% statistic) {
      out$G2 <- lav_tables_stat_g2(
        out$obs.prop, out$est.prop,
        out$nobs
      )
    }
    if ("X2" %in% statistic) {
      out$X2 <- lav_tables_stat_x2(
        out$obs.prop, out$est.prop,
        out$nobs
      )
    }
    if ("cor" %in% statistic) {
      cor_1 <- attr(pi0, "COR")
      cor_all <- unlist(lapply(cor_1, function(x) {
        x[lower.tri(x, diag = FALSE)]
      }))
      out$cor <- cor_all[out$id]
    }
  }

  out
}

# G2 statistic
lav_tables_stat_g2 <- function(obs_prop = NULL, est_prop = NULL, nobs = NULL) {
  # not defined if out$obs_prop is (close to) zero
  zero_idx <- which(obs_prop < .Machine$double.eps)
  if (length(zero_idx)) {
    obs_prop[zero_idx] <- as.numeric(NA)
  }
  # the usual G2 formula
  g2 <- 2 * nobs * (obs_prop * log(obs_prop / est_prop))
  g2
}

# X2 (aka X2) statistic
lav_tables_stat_x2 <- function(obs_prop = NULL, est_prop = NULL, nobs = NULL) {
  res_prop <- obs_prop - est_prop
  x2 <- nobs * (res_prop * res_prop) / est_prop
  x2
}

# pairwise tables, rows = tables
lav_tables_pairwise_table <- function(lavobject = NULL, lavdata = NULL,
                                      statistic = character(0L),
                                      g2_min = 3.0,
                                      x2_min = 3.0,
                                      p_value = FALSE) {
  # default statistics
  if (!is.null(lavobject)) {
    if (length(statistic) == 1L && statistic == "default") {
      statistic <- c("X2", "X2.average")
    } else {
      stopifnot(statistic %in% c(
        "X2", "G2", "X2.un", "G2.un",
        "cor", "cor.un",
        "RMSEA.un", "RMSEA",
        "G2.average",
        "G2.nlarge",
        "G2.plarge",
        "X2.average",
        "X2.nlarge",
        "X2.plarge"
      ))
    }
  } else {
    if (length(statistic) == 1L && statistic == "default") {
      # if data, none by default
      statistic <- character(0L)
    } else {
      stopifnot(statistic %in% c(
        "cor.un", "X2.un", "G2.un",
        "RMSEA.un"
      ))
    }
  }

  # identify 'categorical' variables
  # cat.idx <- which(lavdata@ov$type %in% c("ordered","factor"))

  # pairwise tables
  # pairwise.tables <- utils::combn(vartable$name[cat.idx], m=2L)
  # pairwise.tables <- rbind(seq_len(ncol(pairwise.tables)),
  #                         pairwise.tables)
  # ntables <- ncol(pairwise.tables)

  # initial table, observed cell frequencies
  # out <- as.data.frame(t(pairwise.tables))
  # names(out) <- c("id", "lhs", "rhs")

  # collapse approach
  stat_cell <- character(0)
  if (any(c("G2", "G2.average", "G2.plarge", "G2.nlarge") %in% statistic)) {
    stat_cell <- c(stat_cell, "G2")
  }
  if (any(c("X2", "X2.average", "X2.plarge", "X2.nlarge") %in% statistic)) {
    stat_cell <- c(stat_cell, "X2")
  }
  if ("G2" %in% statistic || "RMSEA" %in% statistic) {
    stat_cell <- c(stat_cell, "G2")
  }
  if ("X2.un" %in% statistic) {
    stat_cell <- c(stat_cell, "X2.un")
  }
  if ("G2.un" %in% statistic || "RMSEA.un" %in% statistic) {
    stat_cell <- c(stat_cell, "G2.un")
  }
  if ("cor.un" %in% statistic) {
    stat_cell <- c(stat_cell, "cor.un")
  }
  if ("cor" %in% statistic) {
    stat_cell <- c(stat_cell, "cor")
  }

  # get table with table cells
  out_cell <- lav_tables_pairwise_cells(
    lavobject = lavobject,
    lavdata = lavdata,
    statistic = stat_cell
  )
  # only 1 row per table
  row_idx <- which(!duplicated(out_cell$id))
  if (is.null(out_cell$group)) {
    out <- out_cell[row_idx, c("lhs", "rhs", "nobs"), drop = FALSE]
  } else {
    out <- out_cell[row_idx, c("lhs", "rhs", "group", "nobs"), drop = FALSE]
  }

  # df
  if (length(statistic) > 0L) {
    nrow <- tapply(out_cell$row, INDEX = out_cell$id, FUN = max)
    ncol <- tapply(out_cell$col, INDEX = out_cell$id, FUN = max)
    out$df <- nrow * ncol - nrow - ncol
  }

  # cor
  if ("cor" %in% statistic) {
    out$cor <- out_cell[row_idx, "cor"]
  }

  # cor.un
  if ("cor.un" %in% statistic) {
    out$cor.un <- out_cell[row_idx, "cor.un"]
  }

  # X2
  if ("X2" %in% statistic) {
    out$X2 <- tapply(out_cell$X2,
      INDEX = out_cell$id, FUN = sum,
      na.rm = TRUE
    )
    if (p_value) {
      out$X2.pval <- pchisq(out$X2, df = out$df, lower.tail = FALSE)
    }
  }
  if ("X2.un" %in% statistic) {
    out$X2.un <- tapply(out_cell$X2.un,
      INDEX = out_cell$id, FUN = sum,
      na.rm = TRUE
    )
    if (p_value) {
      out$X2.un.pval <- pchisq(out$X2.un, df = out$df, lower.tail = FALSE)
    }
  }

  # G2
  if ("G2" %in% statistic) {
    out$G2 <- tapply(out_cell$G2,
      INDEX = out_cell$id, FUN = sum,
      na.rm = TRUE
    )
    if (p_value) {
      out$G2.pval <- pchisq(out$G2, df = out$df, lower.tail = FALSE)
    }
  }
  if ("G2.un" %in% statistic) {
    out$G2.un <- tapply(out_cell$G2.un,
      INDEX = out_cell$id, FUN = sum,
      na.rm = TRUE
    )
    if (p_value) {
      out$G2.un.pval <- pchisq(out$G2.un, df = out$df, lower.tail = FALSE)
    }
  }

  if ("RMSEA" %in% statistic) {
    g2 <- tapply(out_cell$G2, INDEX = out_cell$id, FUN = sum, na.rm = TRUE)
    # note: there seems to be a mistake in Appendix 1 eqs 43/44 of Joreskog
    # SSI paper (2005) 'SEM with ordinal variables using LISREL'
    # 2*N*d should N*d
    out$RMSEA <- sqrt(pmax(0, (g2 - out$df) / (out$nobs * out$df)))
    if (p_value) {
      # note: MUST use 1 - pchisq (instead of lower.tail = FALSE)
      # because for ncp > 80, routine only computes lower tail
      out$RMSEA.pval <- 1.0 - pchisq(g2,
        ncp = 0.1 * 0.1 * out$nobs * out$df,
        df = out$df, lower.tail = TRUE
      )
    }
  }
  if ("RMSEA.un" %in% statistic) {
    g2 <- tapply(out_cell$G2.un,
      INDEX = out_cell$id, FUN = sum,
      na.rm = TRUE
    )
    # note: there seems to be a mistake in Appendix 1 eqs 43/44 of Joreskog
    # SSI paper (2005) 'SEM with ordinal variables using LISREL'
    # 2*N*d should N*d
    out$RMSEA.un <- sqrt(pmax(0, (g2 - out$df) / (out$nobs * out$df)))
    if (p_value) {
      # note: MUST use 1 - pchisq (instead of lower.tail = FALSE)
      # because for ncp > 80, routine only computes lower tail
      out$RMSEA.un.pval <- 1.0 - pchisq(g2,
        ncp = 0.1 * 0.1 * out$nobs * out$df,
        df = out$df, lower.tail = TRUE
      )
    }
  }

  if ("G2.average" %in% statistic) {
    out$G2.average <- tapply(out_cell$G2,
      INDEX = out_cell$id, FUN = mean,
      na.rm = TRUE
    )
  }

  if ("G2.nlarge" %in% statistic) {
    out$G2.min <- rep(g2_min, length(out$lhs))
    out$G2.nlarge <- tapply(out_cell$G2,
      INDEX = out_cell$id,
      FUN = function(x) sum(x > g2_min, na.rm = TRUE)
    )
  }

  if ("G2.plarge" %in% statistic) {
    out$G2.min <- rep(g2_min, length(out$lhs))
    out$G2.plarge <- tapply(out_cell$G2,
      INDEX = out_cell$id,
      FUN = function(x) sum(x > g2_min, na.rm = TRUE) / length(x)
    )
  }

  if ("X2.average" %in% statistic) {
    out$X2.average <- tapply(out_cell$X2,
      INDEX = out_cell$id, FUN = mean,
      na.rm = TRUE
    )
  }

  if ("X2.nlarge" %in% statistic) {
    out$X2.min <- rep(x2_min, length(out$lhs))
    out$X2.nlarge <- tapply(out_cell$X2,
      INDEX = out_cell$id,
      FUN = function(x) sum(x > x2_min, na.rm = TRUE)
    )
  }

  if ("X2.plarge" %in% statistic) {
    out$X2.min <- rep(x2_min, length(out$lhs))
    out$X2.plarge <- tapply(out_cell$X2,
      INDEX = out_cell$id,
      FUN = function(x) sum(x > x2_min, na.rm = TRUE) / length(x)
    )
  }

  out
}


lav_tables_oneway <- function(lavobject = NULL, lavdata = NULL,
                              statistic = NULL) {
  # shortcuts
  vartable <- lavdata@ov
  xx <- lavdata@X

  # identify 'categorical' variables
  cat_idx <- which(vartable$type %in% c("ordered", "factor"))
  ncat <- length(cat_idx)

  # do we have any categorical variables?
  if (length(cat_idx) == 0L) {
    lav_msg_warn(gettext("no categorical variables are found"))
    return(data.frame(
      id = integer(0L), lhs = character(0L), rhs = character(0L),
      nobs = integer(0L),
      obs.freq = integer(0L), obs.prop = numeric(0L),
      est.prop = numeric(0L), X2 = numeric(0L)
    ))
  } else {
    labels <- strsplit(vartable$lnam[cat_idx], "\\|")
  }

  # ok, we have an overview of all categorical variables in the data
  ngroups <- length(xx)

  # for each group, for each categorical variable, collect information
  tables <- vector("list", length = ngroups)
  for (g in 1:ngroups) {
    tables[[g]] <- lapply(seq_len(ncat),
      FUN = function(x) {
        idx <- cat_idx[x]
        nrow <- vartable$nlev[idx]
        ncell <- nrow
        nvar <- length(lavdata@ov.names[[g]])
        id <- (g - 1) * nvar + x

        # compute observed frequencies
        freq <- tabulate(xx[[g]][, idx], nbins = ncell)

        list(
          id = rep.int(id, ncell),
          lhs = rep.int(vartable$name[idx], ncell),
          # op = rep.int("freq", ncell),
          rhs = labels[[x]],
          group = rep.int(g, ncell),
          nobs = rep.int(sum(freq), ncell),
          obs.freq = freq,
          obs.prop = freq / sum(freq)
        )
      }
    )
  }

  for (g in 1:ngroups) {
    table_1 <- tables[[g]]
    table_1 <- lapply(table_1, as.data.frame, stringsAsFactors = FALSE)
    if (g == 1L) {
      out <- do.call(rbind, table_1)
    } else {
      out <- rbind(out, do.call(rbind, table_1))
    }
  }
  if (g == 1) {
    # remove group column
    out$group <- NULL
  }

  # default statistics
  if (!is.null(lavobject)) {
    if (length(statistic) == 1L && statistic == "default") {
      statistic <- c("X2")
    } else {
      stopifnot(statistic %in% c(
        "th.un",
        "th", "G2", "X2"
      ))
    }

    # sample based
    # note, there is no G2.un or X2.un: always saturated!
    if ("th.un" %in% statistic) {
      # sample based
      th_un <- unlist(lapply(1:lavdata@ngroups, function(x) {
        if (lavobject@Model@conditional.x) {
          th <- lavobject@SampleStats@res.th[[x]][
            lavobject@SampleStats@th.idx[[x]] > 0
          ]
        } else {
          th <- lavobject@SampleStats@th[[x]][
            lavobject@SampleStats@th.idx[[x]] > 0
          ]
        }
        th_idx <- lavobject@SampleStats@th.idx[[x]][
          lavobject@SampleStats@th.idx[[x]] > 0
        ]
        unname(unlist(tapply(th,
          INDEX = th_idx,
          function(y) c(y, Inf)
        )))
      }))
      # overwrite obs.prop
      # NOTE: if we have exogenous variables, obs.prop will NOT
      #       correspond with qnorm(th_un)
      out$obs.prop <- unname(unlist(tapply(th_un,
        INDEX = out$id,
        FUN = function(x) {
          (pnorm(c(x, Inf)) -
            pnorm(c(-Inf, x)))[-(length(x) + 1)]
        }
      )))

      out$th.un <- th_un
    }

    # model based
    if (any(c("th", "G2", "X2") %in% statistic)) {
      # model based
      th_h0 <- unlist(lapply(1:lavdata@ngroups, function(x) {
        if (lavobject@Model@conditional.x) {
          th <- lavobject@implied$res.th[[x]][
            lavobject@SampleStats@th.idx[[x]] > 0
          ]
        } else {
          th <- lavobject@implied$th[[x]][
            lavobject@SampleStats@th.idx[[x]] > 0
          ]
        }
        th_idx <- lavobject@SampleStats@th.idx[[x]][
          lavobject@SampleStats@th.idx[[x]] > 0
        ]
        unname(unlist(tapply(th,
          INDEX = th_idx,
          function(x) c(x, Inf)
        )))
      }))

      est_prop <- unname(unlist(tapply(th_h0,
        INDEX = out$id,
        FUN = function(x) {
          (pnorm(c(x, Inf)) -
            pnorm(c(-Inf, x)))[-(length(x) + 1)]
        }
      )))
      out$est.prop <- est_prop

      if ("th" %in% statistic) {
        out$th <- th_h0
      }
      if ("G2" %in% statistic) {
        out$G2 <- lav_tables_stat_g2(
          out$obs.prop, out$est.prop,
          out$nobs
        )
      }
      if ("X2" %in% statistic) {
        out$X2 <- lav_tables_stat_x2(
          out$obs.prop, out$est.prop,
          out$nobs
        )
      }
    }
  } else {
    if (length(statistic) == 1L && statistic == "default") {
      # if data, none by default
      statistic <- character(0L)
    } else {
      stopifnot(statistic %in% c("th.un"))
    }

    if ("th.un" %in% statistic) {
      out$th.un <- unlist(tapply(out$obs.prop,
        INDEX = out$id,
        FUN = function(x) qnorm(cumsum(x))
      ))
    }
  }

  out
}

# HJ 15/1/2023 MODIFIED to add sampling weights
# compute pairwise (two-way) frequency tables
lav_tables_pairwise_freq_cell <- function(lavdata = NULL,
                                          as_data_frame = TRUE) {
  # shortcuts
  vartable <- as.data.frame(lavdata@ov, stringsAsFactors = FALSE)
  xx <- lavdata@X
  # ov_names <- lavdata@ov.names
  ngroups <- lavdata@ngroups
  # wt <- lavdata@weights

  # identify 'categorical' variables
  cat_idx <- which(vartable$type %in% c("ordered", "factor"))

  # do we have any categorical variables?
  if (length(cat_idx) == 0L) {
    lav_msg_stop(gettext("no categorical variables are found"))
  } else if (length(cat_idx) == 1L) {
    lav_msg_stop(gettext("at least two categorical variables are needed"))
  }

  # pairwise tables
  pairwise_tables <- utils::combn(vartable$name[cat_idx], m = 2L)
  pairwise_tables <- rbind(pairwise_tables, seq_len(ncol(pairwise_tables)))
  ntables <- ncol(pairwise_tables)

  # for each group, for each pairwise table, collect information
  tables <- vector("list", length = ngroups)
  for (g in 1:ngroups) {
    tables[[g]] <- apply(pairwise_tables,
      MARGIN = 2,
      FUN = function(x) {
        idx1 <- which(vartable$name == x[1])
        idx2 <- which(vartable$name == x[2])
        id <- (g - 1) * ntables + as.numeric(x[3])
        nrow <- vartable$nlev[idx1]
        ncol <- vartable$nlev[idx2]
        ncell <- nrow * ncol

        # compute two-way observed frequencies
        y1 <- xx[[g]][, idx1]
        y2 <- xx[[g]][, idx2]
        # freq <- table(y1, y2) # we lose missings; useNA is ugly
        freq <- lav_bvord_freq(y1, y2)

        # >>>>>>>> HJ/MK PML CODE >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        # If we want to use weighted frequencies we can use the code
        # below. However, it will probably make sense only when the
        # weights are normalised. If they're not, we may get quite ugly
        # and nonsensical numbers here. So for now, just keep the
        # lavtables as is (using non-weighted frequencies).
        #
        # freq <- lav_bvord_freq(y1, Y2, wt[[g]])

        # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


        list(
          id = rep.int(id, ncell),
          lhs = rep.int(x[1], ncell),
          # op = rep.int("table", ncell),
          rhs = rep.int(x[2], ncell),
          group = rep.int(g, ncell),
          nobs = rep.int(sum(freq), ncell),
          row = rep.int(seq_len(nrow), times = ncol),
          col = rep(seq_len(ncol), each = nrow),
          obs.freq = lav_mat_vec(freq) # col by col!
        )
      }
    )
  }

  if (as_data_frame) {
    for (g in 1:ngroups) {
      table_1 <- tables[[g]]
      table_1 <- lapply(table_1, as.data.frame, stringsAsFactors = FALSE)
      if (g == 1) {
        out <- do.call(rbind, table_1)
      } else {
        out <- rbind(out, do.call(rbind, table_1))
      }
    }
    if (g == 1) {
      # remove group column
      out$group <- NULL
    }
  } else {
    if (ngroups == 1L) {
      out <- tables[[1]]
    } else {
      out <- tables
    }
  }

  out
}


# low-level function to compute expected proportions per cell
# object
lav_tables_pairwise_model_pi <- function(lavobject = NULL) {
  stopifnot(lavobject@Model@categorical)

  # shortcuts
  lavmodel <- lavobject@Model
  implied <- lavobject@implied
  ngroups <- lavobject@Data@ngroups
  ov_types <- lavobject@Data@ov$type
  th_idx <- lavobject@Model@th.idx
  sigma_hat <- if (lavmodel@conditional.x) implied$res.cov else implied$cov
  th <- if (lavmodel@conditional.x) implied$res.th else implied$th

  pi0 <- vector("list", length = ngroups)
  for (g in 1:ngroups) {
    sigmahat <- sigma_hat[[g]]
    cors <- sigmahat[lower.tri(sigmahat)]
    if (any(abs(cors) > 1)) {
      lav_msg_warn(gettext(
        "some model-implied correlations are larger than 1.0"))
    }
    nvar <- nrow(sigmahat)

    # shortcut for all ordered - tablewise
    if (all(ov_types == "ordered") && !is.null(lavobject@Cache[[g]]$long)) {
      # FREQ.OBS <- c(FREQ.OBS, lavobject@Cache[[g]]$bifreq)
      long2 <- lav_pml_longvec_th_rho(
        no_x = nvar,
        all_thres = th[[g]],
        index_var_of_thres = th_idx[[g]],
        rho_xixj = cors
      )
      # get expected probability per table, per pair
      pi0[[g]] <- lav_pml_expprob_vec(
        ind_vec = lavobject@Cache[[g]]$long,
        th_rho_vec = long2
      )
    } else {
      pi_group <- integer(0)
      # order! first i, then j, lav_mat_vec(table)!
      for (i in seq_len(nvar - 1L)) {
        for (j in (i + 1L):nvar) {
          if (ov_types[i] == "ordered" && ov_types[j] == "ordered") {
            pi_table <- lav_bvord_noexo_pi(
              rho = sigmahat[i, j],
              th_y1 = th[[g]][th_idx[[g]] == i],
              th_y2 = th[[g]][th_idx[[g]] == j]
            )
            pi_group <- c(pi_group, lav_mat_vec(pi_table))
          }
        }
      }
      pi0[[g]] <- pi_group
    }
  } # g

  # add  COR/th/TH.IDX
  attr(pi0, "COR") <- sigma_hat
  attr(pi0, "TH") <- th
  attr(pi0, "TH.IDX") <- th_idx

  pi0
}

# low-level function to compute expected proportions per cell
# using sample-based correlations + thresholds
#
# object can be either lavData or lavaan class
lav_tables_pairwise_sample_pi <- function(lavobject = NULL, lavdata = NULL) {
  # get COR, TH and th.idx
  if (!is.null(lavobject)) {
    if (lavobject@Model@conditional.x) {
      cor_1 <- lavobject@SampleStats@res.cov
      th <- lavobject@SampleStats@res.th
    } else {
      cor_1 <- lavobject@SampleStats@cov
      th <- lavobject@SampleStats@th
    }
    th_idx <- lavobject@SampleStats@th.idx
  } else if (!is.null(lavdata)) {
    fit_un <- lav_object_cor(object = lavdata, se = "none", output = "fit")
    if (fit_un@Model@conditional.x) {
      cor_1 <- fit_un@SampleStats@res.cov
      th <- fit_un@SampleStats@res.th
    } else {
      cor_1 <- fit_un@SampleStats@cov
      th <- fit_un@SampleStats@th
    }
    th_idx <- fit_un@SampleStats@th.idx
  } else {
    lav_msg_stop(gettext("both lavobject and lavdata are NULL"))
  }

  lav_tables_pairwise_sample_pi_cor(
    cor_1 = cor_1, th = th,
    th_idx = th_idx
  )
}

# low-level function to compute expected proportions per cell
lav_tables_pairwise_sample_pi_cor <- function(cor_1 = NULL, th = NULL, # nolint
                                              th_idx = NULL) {
  ngroups <- length(cor_1)

  pi0 <- vector("list", length = ngroups)
  for (g in 1:ngroups) {
    sigmahat <- cor_1[[g]]
    cors <- sigmahat[lower.tri(sigmahat)]
    if (any(abs(cors) > 1)) {
      lav_msg_warn(gettext(
        "some model-implied correlations are larger than 1.0"))
    }
    nvar <- nrow(sigmahat)
    th_idx_1 <- th_idx[[g]]

    # reconstruct ov.types
    ov_types <- rep("numeric", nvar)
    ord_idx <- unique(th_idx_1[th_idx_1 > 0])
    ov_types[ord_idx] <- "ordered"

    pi_group <- integer(0)
    # order! first i, then j, lav_mat_vec(table)!
    for (i in seq_len(nvar - 1L)) {
      for (j in (i + 1L):nvar) {
        if (ov_types[i] == "ordered" && ov_types[j] == "ordered") {
          pi_table <- lav_bvord_noexo_pi(
            rho = sigmahat[i, j],
            th_y1 = th[[g]][th_idx_1 == i],
            th_y2 = th[[g]][th_idx_1 == j]
          )
          pi_group <- c(pi_group, lav_mat_vec(pi_table))
        }
      }
    }
    pi0[[g]] <- pi_group
  } # g

  # add  COR/TH/TH.IDX
  attr(pi0, "COR") <- cor_1
  attr(pi0, "TH") <- th
  attr(pi0, "TH.IDX") <- th_idx

  pi0
}

# low-level function to compute expected proportions per PATTERN
# using sample-based correlations + thresholds
#
# object can be either lavData or lavaan class
#
# only valid if estimator = FML, POM or NOR
#
lav_tables_resp_pi <- function(lavobject = NULL, lavdata = NULL,
                               est = "h0") {
  # shortcuts
  if (!is.null(lavobject)) {
    lavmodel <- lavobject@Model
    implied <- lavobject@implied
  }
  ngroups <- lavdata@ngroups

  # h0 or unrestricted?
  if (est == "h0") {
    sigma_hat <- if (lavmodel@conditional.x) implied$res.cov else implied$cov
    th <- if (lavmodel@conditional.x) implied$res.th else implied$th
    th_idx_1 <- lavobject@SampleStats@th.idx
  } else {
    if (is.null(lavobject)) {
      fit_un <- lav_object_cor(object = lavdata, se = "none", output = "fit")
      sigma_hat <-
        if (fit_un@Model@conditional.x) fit_un@implied$res.cov else fit_un@implied$cov
      th <-
        if (fit_un@Model@conditional.x) fit_un@implied$res.th else fit_un@implied$th
      th_idx_1 <- fit_un@SampleStats@th.idx
    } else {
      if (lavobject@Model@conditional.x) {
        sigma_hat <- lavobject@SampleStats@res.cov
        th <- lavobject@SampleStats@res.th
      } else {
        sigma_hat <- lavobject@SampleStats@cov
        th <- lavobject@SampleStats@th
      }
      th_idx_1 <- lavobject@SampleStats@th.idx
    }
  }

  pi0 <- vector("list", length = ngroups)
  for (g in 1:ngroups) {
    sigmahat <- sigma_hat[[g]]
    cors <- sigmahat[lower.tri(sigmahat)]
    if (any(abs(cors) > 1)) {
      lav_msg_warn(gettext(
        "some model-implied correlations are larger than 1.0"))
    }
    nvar <- nrow(sigmahat)
    th_idx <- th_idx_1[[g]]
    mean_1 <- rep(0, nvar)

    # reconstruct ov.types
    ov_types <- rep("numeric", nvar)
    ord_idx <- unique(th_idx[th_idx > 0])
    ov_types[ord_idx] <- "ordered"

    if (all(ov_types == "ordered")) {
      # get patterns ## FIXME GET it
      if (!is.null(lavdata@Rp[[g]]$pat)) {
        pat_1 <- lavdata@Rp[[g]]$pat
      } else {
        pat_1 <- lav_data_resp_patterns(lavdata@X[[g]])$pat
      }
      npatterns <- nrow(pat_1)
      freq <- as.numeric(rownames(pat_1))
      pi_group <- numeric(npatterns)
      th_var <- lapply(
        1:nvar,
        function(x) c(-Inf, th[[g]][th_idx == x], +Inf)
      )
      # FIXME!!! ok to set diagonal to 1.0?
      diag(sigmahat) <- 1.0
      for (r in 1:npatterns) {
        # compute probability for each pattern
        lower <- sapply(1:nvar, function(x) th_var[[x]][pat_1[r, x]])
        upper <- sapply(1:nvar, function(x) th_var[[x]][pat_1[r, x] + 1L])
        # handle missing values
        na_idx <- which(is.na(pat_1[r, ]))
        if (length(na_idx) > 0L) {
          lower <- lower[-na_idx]
          upper <- upper[-na_idx]
          mean_r <- mean_1[-na_idx]
          sigmahat_r <- sigmahat[-na_idx, -na_idx, drop = FALSE]
        } else {
          mean_r <- mean_1
          sigmahat_r <- sigmahat
        }
        pi_group[r] <- sadmvn(lower, upper,
          mean = mean_r,
          varcov = sigmahat_r
        )
      }
    } else { # case-wise
      pi_group <- rep(as.numeric(NA), lavdata@nobs[[g]])
      lav_msg_warn(gettext("casewise PI not implemented"))
    }

    pi0[[g]] <- pi_group
  } # g

  pi0
}

lav_tables_table_format <- function(out, lavdata = lavdata,
                                    lavobject = lavobject) {
  # determine column we need
  names_1 <- names(out)
  stat_idx <- which(names_1 %in% c(
    "cor", "cor.un",
    "G2", "G2.un",
    "X2", "X2.un",
    "RMSEA", "RMSEA.un",
    "G2.average", "G2.plarge", "G2.nlarge",
    "X2.average", "X2.plarge", "X2.nlarge"
  ))
  if (length(stat_idx) == 0) {
    if (!is.null(out$obs.freq)) {
      stat_idx <- which(names_1 == "obs.freq")
    } else if (!is.null(out$nobs)) {
      stat_idx <- which(names_1 == "nobs")
    }
    uni <- NULL
  } else if (length(stat_idx) > 1) {
    lav_msg_stop(gettext(
      "more than one statistic for table output:"),
      paste(names_1[stat_idx], collapse = " ")
    )
  } else {
    # univariate version of same statistic
    if (names_1[stat_idx] == "G2.average") {
      uni <- lavTables(lavobject, dimension = 1L, statistic = "G2")
    } else if (names_1[stat_idx] == "X2.average") {
      uni <- lavTables(lavobject, dimension = 1L, statistic = "X2")
    } else {
      uni <- NULL
    }
  }

  out_1 <- vector("list", length = lavdata@ngroups)
  for (g in 1:lavdata@ngroups) {
    if (lavdata@ngroups == 1L) { # no group column
      stat <- out[[stat_idx]]
    } else {
      stat <- out[[stat_idx]][out$group == g]
    }
    rn <- lavdata@ov.names[[g]]
    out_1[[g]] <- lav_getcov(stat, diagonal = FALSE, lower = FALSE, names = rn)
    # change diagonal elements: replace by univariate stat
    # if possible
    diag(out_1[[g]]) <- as.numeric(NA)
    if (!is.null(uni)) {
      if (!is.null(uni$group)) {
        idx <- which(uni$group == g)
      } else {
        idx <- seq_along(uni$lhs)
      }
      if (names_1[stat_idx] == "G2.average") {
        diag(out_1[[g]]) <- tapply(uni$G2[idx],
          INDEX = uni$id[idx],
          FUN = mean
        )
      } else if (names_1[stat_idx] == "X2.average") {
        diag(out_1[[g]]) <- tapply(uni$X2[idx],
          INDEX = uni$id[idx],
          FUN = mean
        )
      }
    } else if (names_1[stat_idx] %in% c("cor", "cor.un")) {
      diag(out_1[[g]]) <- 1
    }
    class(out_1[[g]]) <- c("lavaan.matrix.symmetric", "matrix")
  }
  if (lavdata@ngroups > 1L) {
    names(out_1) <- lavdata@group.label
    out <- out_1
  } else {
    out <- out_1[[1]]
  }

  out
}

lav_tables_cells_format <- function(out, lavdata = lavdata,
                                    drop_list_single_group = FALSE) {
  out_1 <- vector("list", length = lavdata@ngroups)
  if (is.null(out$group)) {
    out$group <- rep(1L, length(out$lhs))
  }
  # do we have a statistic?
  # determine column we need
  names_1 <- names(out)
  stat_idx <- which(names_1 %in% c(
    "cor", "cor.un",
    "G2", "G2.un",
    "X2", "X2.un",
    "RMSEA", "RMSEA.un",
    "G2.average", "G2.plarge", "G2.nlarge",
    "X2.average", "X2.plarge", "X2.nlarge"
  ))
  if (length(stat_idx) == 0) {
    statistic <- "obs.freq"
  } else if (length(stat_idx) > 1) {
    lav_msg_stop(gettext(
      "more than one statistic for table output:"),
      paste(names_1[stat_idx], collapse = " ")
    )
  } else {
    statistic <- names_1[stat_idx]
  }

  for (g in 1:lavdata@ngroups) {
    case_idx <- which(out$group == g)
    id_group <- unique(out$id[out$group == g])
    tmp <- lapply(id_group, function(x) {
      tx <- out[out$id == x, ]
      m <- matrix(
        tx[, statistic],
        max(tx$row), max(tx$col)
      )
      rownames(m) <- unique(tx$row)
      colnames(m) <- unique(tx$col)
      class(m) <- c("lavaan.matrix", "matrix")
      m
    })
    names(tmp) <- unique(paste(out$lhs[case_idx], out$rhs[case_idx],
      sep = "_"
    ))
    out_1[[g]] <- tmp
  }

  if (lavdata@ngroups == 1L && drop_list_single_group) {
    out_1 <- out_1[[1]]
  } else {
    if (length(lavdata@group.label) > 0L) {
      names(out_1) <- unlist(lavdata@group.label)
    }
  }

  out_1
}



# The function lav_tables_uni_freq_cell computes the univariate (one-way)
# frequency tables.
# The function closely follows the "logic" of the lavaan function
# lav_tables_pairwise_freq_cell.
# The output is either a list or a data.frame depending on the value the logical
# input argument as.data.frame. Either way, the same information is contained
# which is:
# a) the observed (univariate) frequencies f_ia, i=1,...,p (variables),
#    a=1,...,ci (response categories), with a index running faster than i index.
# b) an index vector with the name varb which indicates which variable each
#    frequency refers to.
# c) an index vector with the name group which indicates which group each
#    frequency refers to when multi-group analysis.
# d) an index vector with the name level which indicates which level within
#    each ordinal variable each frequency refers to.
# e) a vector nobs which gives how many cases were considered to compute the
#    corresponding frequency. Since we use the available data for each variable
#    when missing=="available_cases" we expect these numbers to differ when
#    missing values are present.
# f) an index vector with the name id indexing each univariate table,
#    1 goes to first variable in the first group, 2 to 2nd variable in the
#    second group and so on. The last table has the index equal to
#                                   (no of groups) * (no of variables).

lav_tables_uni_freq_cell <- function(lavdata = NULL,
                                            as_data_frame = TRUE) {
  # shortcuts
  vartable <- as.data.frame(lavdata@ov, stringsAsFactors = FALSE)
  xx <- lavdata@X
  # ov_names <- lavdata@ov.names
  ngroups <- lavdata@ngroups

  # identify 'categorical' variables
  cat_idx <- which(vartable$type %in% c("ordered", "factor"))

  # do we have any categorical variables?
  if (length(cat_idx) == 0L) {
    lav_msg_stop(gettext("no categorical variables are found"))
  }

  # univariate tables
  univariate_tables <- vartable$name[cat_idx]
  univariate_tables <- rbind(univariate_tables,
    seq_along(univariate_tables),
    deparse.level = 0
  )
  ntables <- ncol(univariate_tables)

  # for each group, for each pairwise table, collect information
  uni_tables <- vector("list", length = ngroups)
  for (g in 1:ngroups) {
    uni_tables[[g]] <- apply(univariate_tables,
      MARGIN = 2,
      FUN = function(x) {
        idx1 <- which(vartable$name == x[1])
        id <- (g - 1) * ntables + as.numeric(x[2])
        ncell <- vartable$nlev[idx1]

        # compute one-way observed frequencies
        y1 <- xx[[g]][, idx1]
        uni_freq <- tabulate(y1, nbins = max(y1, na.rm = TRUE))

        list(
          id = rep.int(id, ncell),
          varb = rep.int(x[1], ncell),
          group = rep.int(g, ncell),
          nobs = rep.int(sum(uni_freq), ncell),
          level = seq_len(ncell),
          obs.freq = uni_freq
        )
      }
    )
  }

  if (as_data_frame) {
    for (g in 1:ngroups) {
      uni_table <- uni_tables[[g]]
      uni_table <- lapply(uni_table, as.data.frame,
        stringsAsFactors = FALSE
      )
      if (g == 1) {
        out <- do.call(rbind, uni_table)
      } else {
        out <- rbind(out, do.call(rbind, uni_table))
      }
    }
    if (g == 1) {
      # remove group column
      out$group <- NULL
    }
  } else {
    if (ngroups == 1L) {
      out <- uni_tables[[1]]
    } else {
      out <- uni_tables
    }
  }

  out
}
