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

lavTables <- function(object,
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
                      # Bonferonni
                      # alpha.adj    = FALSE,
                      # output format
                      output = "data.frame",
                      patternAsString = TRUE) {
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
  lavdata <- lavData(object, ordered = categorical, group = group)

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
      patternAsString = patternAsString
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
      G2.min = G2.min,
      X2.min = X2.min,
      p.value = p.value
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
                               statistic = NULL, patternAsString = TRUE) {
  # this only works if we have 'categorical' variables
  cat.idx <- which(lavdata@ov$type %in% c("ordered", "factor"))
  if (length(cat.idx) == 0L) {
    lav_msg_warn(gettext("no categorical variables are found"))
    return(data.frame(
      pattern = character(0L), nobs = integer(0L),
      obs.freq = integer(0L), obs.prop = numeric(0L)
    ))
  }
  # no support yet for mixture of endogenous ordered + numeric variables
  if (!is.null(lavobject) &&
    length(lavNames(lavobject, "ov.nox")) > length(cat.idx)) {
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
    obs.freq <- as.integer(rownames(pat))
    if (patternAsString) {
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
    NOBS <- sum(obs.freq)
    pat$nobs <- rep(NOBS, nrow(pat))
    pat$obs.freq <- obs.freq
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
    PI <- lav_tables_resp_pi(
      lavobject = lavobject, lavdata = lavdata,
      est = "h1"
    )

    out$est.prop.un <- unlist(PI)
    if ("G2.un" %in% statistic) {
      out$G2.un <- lav_tables_stat_G2(
        out$obs.prop, out$est.prop.un,
        out$nobs
      )
    }
    if ("X2.un" %in% statistic) {
      out$X2.un <- lav_tables_stat_X2(
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

    PI <- lav_tables_resp_pi(
      lavobject = lavobject, lavdata = lavdata,
      est = "h0"
    )

    out$est.prop <- unlist(PI)
    if ("G2" %in% statistic) {
      out$G2 <- lav_tables_stat_G2(
        out$obs.prop, out$est.prop,
        out$nobs
      )
    }
    if ("X2" %in% statistic) {
      out$X2 <- lav_tables_stat_X2(
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
  cat.idx <- which(lavdata@ov$type %in% c("ordered", "factor"))
  if (length(cat.idx) == 0L) {
    lav_msg_warn(gettext("no categorical variables are found"))
    return(data.frame(
      id = integer(0L), lhs = character(0L), rhs = character(0L),
      nobs = integer(0L), row = integer(0L), col = integer(0L),
      obs.freq = integer(0L), obs.prop = numeric(0L)
    ))
  }
  if (length(cat.idx) == 1L) {
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
    as.data.frame. = TRUE
  )
  out$obs.prop <- out$obs.freq / out$nobs

  if (any(c("cor.un", "th.un", "X2.un", "G2.un") %in% statistic)) {
    PI <- lav_tables_pairwise_sample_pi(
      lavobject = lavobject,
      lavdata = lavdata
    )
    out$est.prop.un <- unlist(PI)
    if ("G2.un" %in% statistic) {
      out$G2.un <- lav_tables_stat_G2(
        out$obs.prop, out$est.prop.un,
        out$nobs
      )
    }
    if ("X2.un" %in% statistic) {
      out$X2.un <- lav_tables_stat_X2(
        out$obs.prop, out$est.prop.un,
        out$nobs
      )
    }

    if ("cor.un" %in% statistic) {
      COR <- attr(PI, "COR")
      cor.all <- unlist(lapply(COR, function(x) {
        x[lower.tri(x, diag = FALSE)]
      }))
      out$cor.un <- cor.all[out$id]
    }
  }
  if (any(c("cor", "th", "X2", "G2") %in% statistic)) {
    PI <- lav_tables_pairwise_model_pi(lavobject = lavobject)
    out$est.prop <- unlist(PI)
    if ("G2" %in% statistic) {
      out$G2 <- lav_tables_stat_G2(
        out$obs.prop, out$est.prop,
        out$nobs
      )
    }
    if ("X2" %in% statistic) {
      out$X2 <- lav_tables_stat_X2(
        out$obs.prop, out$est.prop,
        out$nobs
      )
    }
    if ("cor" %in% statistic) {
      COR <- attr(PI, "COR")
      cor.all <- unlist(lapply(COR, function(x) {
        x[lower.tri(x, diag = FALSE)]
      }))
      out$cor <- cor.all[out$id]
    }
  }

  out
}

# G2 statistic
lav_tables_stat_G2 <- function(obs.prop = NULL, est.prop = NULL, nobs = NULL) {
  # not defined if out$obs.prop is (close to) zero
  zero.idx <- which(obs.prop < .Machine$double.eps)
  if (length(zero.idx)) {
    obs.prop[zero.idx] <- as.numeric(NA)
  }
  # the usual G2 formula
  G2 <- 2 * nobs * (obs.prop * log(obs.prop / est.prop))
  G2
}

# X2 (aka X2) statistic
lav_tables_stat_X2 <- function(obs.prop = NULL, est.prop = NULL, nobs = NULL) {
  res.prop <- obs.prop - est.prop
  X2 <- nobs * (res.prop * res.prop) / est.prop
  X2
}

# pairwise tables, rows = tables
lav_tables_pairwise_table <- function(lavobject = NULL, lavdata = NULL,
                                      statistic = character(0L),
                                      G2.min = 3.0,
                                      X2.min = 3.0,
                                      p.value = FALSE) {
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
  stat.cell <- character(0)
  if (any(c("G2", "G2.average", "G2.plarge", "G2.nlarge") %in% statistic)) {
    stat.cell <- c(stat.cell, "G2")
  }
  if (any(c("X2", "X2.average", "X2.plarge", "X2.nlarge") %in% statistic)) {
    stat.cell <- c(stat.cell, "X2")
  }
  if ("G2" %in% statistic || "RMSEA" %in% statistic) {
    stat.cell <- c(stat.cell, "G2")
  }
  if ("X2.un" %in% statistic) {
    stat.cell <- c(stat.cell, "X2.un")
  }
  if ("G2.un" %in% statistic || "RMSEA.un" %in% statistic) {
    stat.cell <- c(stat.cell, "G2.un")
  }
  if ("cor.un" %in% statistic) {
    stat.cell <- c(stat.cell, "cor.un")
  }
  if ("cor" %in% statistic) {
    stat.cell <- c(stat.cell, "cor")
  }

  # get table with table cells
  out.cell <- lav_tables_pairwise_cells(
    lavobject = lavobject,
    lavdata = lavdata,
    statistic = stat.cell
  )
  # only 1 row per table
  row.idx <- which(!duplicated(out.cell$id))
  if (is.null(out.cell$group)) {
    out <- out.cell[row.idx, c("lhs", "rhs", "nobs"), drop = FALSE]
  } else {
    out <- out.cell[row.idx, c("lhs", "rhs", "group", "nobs"), drop = FALSE]
  }

  # df
  if (length(statistic) > 0L) {
    nrow <- tapply(out.cell$row, INDEX = out.cell$id, FUN = max)
    ncol <- tapply(out.cell$col, INDEX = out.cell$id, FUN = max)
    out$df <- nrow * ncol - nrow - ncol
  }

  # cor
  if ("cor" %in% statistic) {
    out$cor <- out.cell[row.idx, "cor"]
  }

  # cor.un
  if ("cor.un" %in% statistic) {
    out$cor.un <- out.cell[row.idx, "cor.un"]
  }

  # X2
  if ("X2" %in% statistic) {
    out$X2 <- tapply(out.cell$X2,
      INDEX = out.cell$id, FUN = sum,
      na.rm = TRUE
    )
    if (p.value) {
      out$X2.pval <- pchisq(out$X2, df = out$df, lower.tail = FALSE)
    }
  }
  if ("X2.un" %in% statistic) {
    out$X2.un <- tapply(out.cell$X2.un,
      INDEX = out.cell$id, FUN = sum,
      na.rm = TRUE
    )
    if (p.value) {
      out$X2.un.pval <- pchisq(out$X2.un, df = out$df, lower.tail = FALSE)
    }
  }

  # G2
  if ("G2" %in% statistic) {
    out$G2 <- tapply(out.cell$G2,
      INDEX = out.cell$id, FUN = sum,
      na.rm = TRUE
    )
    if (p.value) {
      out$G2.pval <- pchisq(out$G2, df = out$df, lower.tail = FALSE)
    }
  }
  if ("G2.un" %in% statistic) {
    out$G2.un <- tapply(out.cell$G2.un,
      INDEX = out.cell$id, FUN = sum,
      na.rm = TRUE
    )
    if (p.value) {
      out$G2.un.pval <- pchisq(out$G2.un, df = out$df, lower.tail = FALSE)
    }
  }

  if ("RMSEA" %in% statistic) {
    G2 <- tapply(out.cell$G2, INDEX = out.cell$id, FUN = sum, na.rm = TRUE)
    # note: there seems to be a mistake in Appendix 1 eqs 43/44 of Joreskog
    # SSI paper (2005) 'SEM with ordinal variables using LISREL'
    # 2*N*d should N*d
    out$RMSEA <- sqrt(pmax(0, (G2 - out$df) / (out$nobs * out$df)))
    if (p.value) {
      # note: MUST use 1 - pchisq (instead of lower.tail = FALSE)
      # because for ncp > 80, routine only computes lower tail
      out$RMSEA.pval <- 1.0 - pchisq(G2,
        ncp = 0.1 * 0.1 * out$nobs * out$df,
        df = out$df, lower.tail = TRUE
      )
    }
  }
  if ("RMSEA.un" %in% statistic) {
    G2 <- tapply(out.cell$G2.un,
      INDEX = out.cell$id, FUN = sum,
      na.rm = TRUE
    )
    # note: there seems to be a mistake in Appendix 1 eqs 43/44 of Joreskog
    # SSI paper (2005) 'SEM with ordinal variables using LISREL'
    # 2*N*d should N*d
    out$RMSEA.un <- sqrt(pmax(0, (G2 - out$df) / (out$nobs * out$df)))
    if (p.value) {
      # note: MUST use 1 - pchisq (instead of lower.tail = FALSE)
      # because for ncp > 80, routine only computes lower tail
      out$RMSEA.un.pval <- 1.0 - pchisq(G2,
        ncp = 0.1 * 0.1 * out$nobs * out$df,
        df = out$df, lower.tail = TRUE
      )
    }
  }

  if ("G2.average" %in% statistic) {
    out$G2.average <- tapply(out.cell$G2,
      INDEX = out.cell$id, FUN = mean,
      na.rm = TRUE
    )
  }

  if ("G2.nlarge" %in% statistic) {
    out$G2.min <- rep(G2.min, length(out$lhs))
    out$G2.nlarge <- tapply(out.cell$G2,
      INDEX = out.cell$id,
      FUN = function(x) sum(x > G2.min, na.rm = TRUE)
    )
  }

  if ("G2.plarge" %in% statistic) {
    out$G2.min <- rep(G2.min, length(out$lhs))
    out$G2.plarge <- tapply(out.cell$G2,
      INDEX = out.cell$id,
      FUN = function(x) sum(x > G2.min, na.rm = TRUE) / length(x)
    )
  }

  if ("X2.average" %in% statistic) {
    out$X2.average <- tapply(out.cell$X2,
      INDEX = out.cell$id, FUN = mean,
      na.rm = TRUE
    )
  }

  if ("X2.nlarge" %in% statistic) {
    out$X2.min <- rep(X2.min, length(out$lhs))
    out$X2.nlarge <- tapply(out.cell$X2,
      INDEX = out.cell$id,
      FUN = function(x) sum(x > X2.min, na.rm = TRUE)
    )
  }

  if ("X2.plarge" %in% statistic) {
    out$X2.min <- rep(X2.min, length(out$lhs))
    out$X2.plarge <- tapply(out.cell$X2,
      INDEX = out.cell$id,
      FUN = function(x) sum(x > X2.min, na.rm = TRUE) / length(x)
    )
  }

  out
}


lav_tables_oneway <- function(lavobject = NULL, lavdata = NULL,
                              statistic = NULL) {
  # shortcuts
  vartable <- lavdata@ov
  X <- lavdata@X

  # identify 'categorical' variables
  cat.idx <- which(vartable$type %in% c("ordered", "factor"))
  ncat <- length(cat.idx)

  # do we have any categorical variables?
  if (length(cat.idx) == 0L) {
    lav_msg_warn(gettext("no categorical variables are found"))
    return(data.frame(
      id = integer(0L), lhs = character(0L), rhs = character(0L),
      nobs = integer(0L),
      obs.freq = integer(0L), obs.prop = numeric(0L),
      est.prop = numeric(0L), X2 = numeric(0L)
    ))
  } else {
    labels <- strsplit(vartable$lnam[cat.idx], "\\|")
  }

  # ok, we have an overview of all categorical variables in the data
  ngroups <- length(X)

  # for each group, for each categorical variable, collect information
  TABLES <- vector("list", length = ngroups)
  for (g in 1:ngroups) {
    TABLES[[g]] <- lapply(seq_len(ncat),
      FUN = function(x) {
        idx <- cat.idx[x]
        nrow <- vartable$nlev[idx]
        ncell <- nrow
        nvar <- length(lavdata@ov.names[[g]])
        id <- (g - 1) * nvar + x

        # compute observed frequencies
        FREQ <- tabulate(X[[g]][, idx], nbins = ncell)

        list(
          id = rep.int(id, ncell),
          lhs = rep.int(vartable$name[idx], ncell),
          # op = rep.int("freq", ncell),
          rhs = labels[[x]],
          group = rep.int(g, ncell),
          nobs = rep.int(sum(FREQ), ncell),
          obs.freq = FREQ,
          obs.prop = FREQ / sum(FREQ)
        )
      }
    )
  }

  for (g in 1:ngroups) {
    TABLE <- TABLES[[g]]
    TABLE <- lapply(TABLE, as.data.frame, stringsAsFactors = FALSE)
    if (g == 1L) {
      out <- do.call(rbind, TABLE)
    } else {
      out <- rbind(out, do.call(rbind, TABLE))
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
      th <- unlist(lapply(1:lavdata@ngroups, function(x) {
        if (lavobject@Model@conditional.x) {
          TH <- lavobject@SampleStats@res.th[[x]][
            lavobject@SampleStats@th.idx[[x]] > 0
          ]
        } else {
          TH <- lavobject@SampleStats@th[[x]][
            lavobject@SampleStats@th.idx[[x]] > 0
          ]
        }
        TH.IDX <- lavobject@SampleStats@th.idx[[x]][
          lavobject@SampleStats@th.idx[[x]] > 0
        ]
        unname(unlist(tapply(TH,
          INDEX = TH.IDX,
          function(y) c(y, Inf)
        )))
      }))
      # overwrite obs.prop
      # NOTE: if we have exogenous variables, obs.prop will NOT
      #       correspond with qnorm(th)
      out$obs.prop <- unname(unlist(tapply(th,
        INDEX = out$id,
        FUN = function(x) {
          (pnorm(c(x, Inf)) -
            pnorm(c(-Inf, x)))[-(length(x) + 1)]
        }
      )))

      out$th.un <- th
    }

    # model based
    if (any(c("th", "G2", "X2") %in% statistic)) {
      # model based
      th.h0 <- unlist(lapply(1:lavdata@ngroups, function(x) {
        if (lavobject@Model@conditional.x) {
          TH <- lavobject@implied$res.th[[x]][
            lavobject@SampleStats@th.idx[[x]] > 0
          ]
        } else {
          TH <- lavobject@implied$th[[x]][
            lavobject@SampleStats@th.idx[[x]] > 0
          ]
        }
        TH.IDX <- lavobject@SampleStats@th.idx[[x]][
          lavobject@SampleStats@th.idx[[x]] > 0
        ]
        unname(unlist(tapply(TH,
          INDEX = TH.IDX,
          function(x) c(x, Inf)
        )))
      }))

      est.prop <- unname(unlist(tapply(th.h0,
        INDEX = out$id,
        FUN = function(x) {
          (pnorm(c(x, Inf)) -
            pnorm(c(-Inf, x)))[-(length(x) + 1)]
        }
      )))
      out$est.prop <- est.prop

      if ("th" %in% statistic) {
        out$th <- th.h0
      }
      if ("G2" %in% statistic) {
        out$G2 <- lav_tables_stat_G2(
          out$obs.prop, out$est.prop,
          out$nobs
        )
      }
      if ("X2" %in% statistic) {
        out$X2 <- lav_tables_stat_X2(
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
                                          as.data.frame. = TRUE) {
  # shortcuts
  vartable <- as.data.frame(lavdata@ov, stringsAsFactors = FALSE)
  X <- lavdata@X
  ov.names <- lavdata@ov.names
  ngroups <- lavdata@ngroups
  wt <- lavdata@weights

  # identify 'categorical' variables
  cat.idx <- which(vartable$type %in% c("ordered", "factor"))

  # do we have any categorical variables?
  if (length(cat.idx) == 0L) {
    lav_msg_stop(gettext("no categorical variables are found"))
  } else if (length(cat.idx) == 1L) {
    lav_msg_stop(gettext("at least two categorical variables are needed"))
  }

  # pairwise tables
  pairwise.tables <- utils::combn(vartable$name[cat.idx], m = 2L)
  pairwise.tables <- rbind(pairwise.tables, seq_len(ncol(pairwise.tables)))
  ntables <- ncol(pairwise.tables)

  # for each group, for each pairwise table, collect information
  TABLES <- vector("list", length = ngroups)
  for (g in 1:ngroups) {
    TABLES[[g]] <- apply(pairwise.tables,
      MARGIN = 2,
      FUN = function(x) {
        idx1 <- which(vartable$name == x[1])
        idx2 <- which(vartable$name == x[2])
        id <- (g - 1) * ntables + as.numeric(x[3])
        nrow <- vartable$nlev[idx1]
        ncol <- vartable$nlev[idx2]
        ncell <- nrow * ncol

        # compute two-way observed frequencies
        Y1 <- X[[g]][, idx1]
        Y2 <- X[[g]][, idx2]
        # FREQ <- table(Y1, Y2) # we loose missings; useNA is ugly
        FREQ <- lav_bvord_freq(Y1, Y2)

        # >>>>>>>> HJ/MK PML CODE >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        # If we want to use weighted frequencies we can use the code
        # below. However, it will probably make sense only when the
        # weights are normalised. If they're not, we may get quite ugly
        # and nonsensical numbers here. So for now, just keep the
        # lavtables as is (using non-weighted frequencies).
        #
        # FREQ <- lav_bvord_freq(Y1, Y2, wt[[g]])

        # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


        list(
          id = rep.int(id, ncell),
          lhs = rep.int(x[1], ncell),
          # op = rep.int("table", ncell),
          rhs = rep.int(x[2], ncell),
          group = rep.int(g, ncell),
          nobs = rep.int(sum(FREQ), ncell),
          row = rep.int(seq_len(nrow), times = ncol),
          col = rep(seq_len(ncol), each = nrow),
          obs.freq = lav_matrix_vec(FREQ) # col by col!
        )
      }
    )
  }

  if (as.data.frame.) {
    for (g in 1:ngroups) {
      TABLE <- TABLES[[g]]
      TABLE <- lapply(TABLE, as.data.frame, stringsAsFactors = FALSE)
      if (g == 1) {
        out <- do.call(rbind, TABLE)
      } else {
        out <- rbind(out, do.call(rbind, TABLE))
      }
    }
    if (g == 1) {
      # remove group column
      out$group <- NULL
    }
  } else {
    if (ngroups == 1L) {
      out <- TABLES[[1]]
    } else {
      out <- TABLES
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
  ov.types <- lavobject@Data@ov$type
  th.idx <- lavobject@Model@th.idx
  Sigma.hat <- if (lavmodel@conditional.x) implied$res.cov else implied$cov
  TH <- if (lavmodel@conditional.x) implied$res.th else implied$th

  PI <- vector("list", length = ngroups)
  for (g in 1:ngroups) {
    Sigmahat <- Sigma.hat[[g]]
    cors <- Sigmahat[lower.tri(Sigmahat)]
    if (any(abs(cors) > 1)) {
      lav_msg_warn(gettext(
        "some model-implied correlations are larger than 1.0"))
    }
    nvar <- nrow(Sigmahat)

    # shortcut for all ordered - tablewise
    if (all(ov.types == "ordered") && !is.null(lavobject@Cache[[g]]$long)) {
      # FREQ.OBS <- c(FREQ.OBS, lavobject@Cache[[g]]$bifreq)
      long2 <- LongVecTH.Rho(
        no.x = nvar,
        all.thres = TH[[g]],
        index.var.of.thres = th.idx[[g]],
        rho.xixj = cors
      )
      # get expected probability per table, per pair
      PI[[g]] <- pairwiseExpProbVec(
        ind.vec = lavobject@Cache[[g]]$long,
        th.rho.vec = long2
      )
    } else {
      PI.group <- integer(0)
      # order! first i, then j, lav_matrix_vec(table)!
      for (i in seq_len(nvar - 1L)) {
        for (j in (i + 1L):nvar) {
          if (ov.types[i] == "ordered" && ov.types[j] == "ordered") {
            PI.table <- lav_bvord_noexo_pi(
              rho = Sigmahat[i, j],
              th.y1 = TH[[g]][th.idx[[g]] == i],
              th.y2 = TH[[g]][th.idx[[g]] == j]
            )
            PI.group <- c(PI.group, lav_matrix_vec(PI.table))
          }
        }
      }
      PI[[g]] <- PI.group
    }
  } # g

  # add  COR/TH/TH.IDX
  attr(PI, "COR") <- Sigma.hat
  attr(PI, "TH") <- TH
  attr(PI, "TH.IDX") <- th.idx

  PI
}

# low-level function to compute expected proportions per cell
# using sample-based correlations + thresholds
#
# object can be either lavData or lavaan class
lav_tables_pairwise_sample_pi <- function(lavobject = NULL, lavdata = NULL) {
  # get COR, TH and th.idx
  if (!is.null(lavobject)) {
    if (lavobject@Model@conditional.x) {
      COR <- lavobject@SampleStats@res.cov
      TH <- lavobject@SampleStats@res.th
    } else {
      COR <- lavobject@SampleStats@cov
      TH <- lavobject@SampleStats@th
    }
    TH.IDX <- lavobject@SampleStats@th.idx
  } else if (!is.null(lavdata)) {
    fit.un <- lavCor(object = lavdata, se = "none", output = "fit")
    if (fit.un@Model@conditional.x) {
      COR <- fit.un@SampleStats@res.cov
      TH <- fit.un@SampleStats@res.th
    } else {
      COR <- fit.un@SampleStats@cov
      TH <- fit.un@SampleStats@th
    }
    TH.IDX <- fit.un@SampleStats@th.idx
  } else {
    lav_msg_stop(gettext("both lavobject and lavdata are NULL"))
  }

  lav_tables_pairwise_sample_pi_cor(
    COR = COR, TH = TH,
    TH.IDX = TH.IDX
  )
}

# low-level function to compute expected proportions per cell
lav_tables_pairwise_sample_pi_cor <- function(COR = NULL, TH = NULL,
                                              TH.IDX = NULL) {
  ngroups <- length(COR)

  PI <- vector("list", length = ngroups)
  for (g in 1:ngroups) {
    Sigmahat <- COR[[g]]
    cors <- Sigmahat[lower.tri(Sigmahat)]
    if (any(abs(cors) > 1)) {
      lav_msg_warn(gettext(
        "some model-implied correlations are larger than 1.0"))
    }
    nvar <- nrow(Sigmahat)
    th.idx <- TH.IDX[[g]]

    # reconstruct ov.types
    ov.types <- rep("numeric", nvar)
    ord.idx <- unique(th.idx[th.idx > 0])
    ov.types[ord.idx] <- "ordered"

    PI.group <- integer(0)
    # order! first i, then j, lav_matrix_vec(table)!
    for (i in seq_len(nvar - 1L)) {
      for (j in (i + 1L):nvar) {
        if (ov.types[i] == "ordered" && ov.types[j] == "ordered") {
          PI.table <- lav_bvord_noexo_pi(
            rho = Sigmahat[i, j],
            th.y1 = TH[[g]][th.idx == i],
            th.y2 = TH[[g]][th.idx == j]
          )
          PI.group <- c(PI.group, lav_matrix_vec(PI.table))
        }
      }
    }
    PI[[g]] <- PI.group
  } # g

  # add  COR/TH/TH.IDX
  attr(PI, "COR") <- COR
  attr(PI, "TH") <- TH
  attr(PI, "TH.IDX") <- TH.IDX

  PI
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
    Sigma.hat <- if (lavmodel@conditional.x) implied$res.cov else implied$cov
    TH <- if (lavmodel@conditional.x) implied$res.th else implied$th
    TH.IDX <- lavobject@SampleStats@th.idx
  } else {
    if (is.null(lavobject)) {
      fit.un <- lavCor(object = lavdata, se = "none", output = "fit")
      Sigma.hat <- if (fit.un@Model@conditional.x) fit.un@implied$res.cov else fit.un@implied$cov
      TH <- if (fit.un@Model@conditional.x) fit.un@implied$res.th else fit.un@implied$th
      TH.IDX <- fit.un@SampleStats@th.idx
    } else {
      if (lavobject@Model@conditional.x) {
        Sigma.hat <- lavobject@SampleStats@res.cov
        TH <- lavobject@SampleStats@res.th
      } else {
        Sigma.hat <- lavobject@SampleStats@cov
        TH <- lavobject@SampleStats@th
      }
      TH.IDX <- lavobject@SampleStats@th.idx
    }
  }

  PI <- vector("list", length = ngroups)
  for (g in 1:ngroups) {
    Sigmahat <- Sigma.hat[[g]]
    cors <- Sigmahat[lower.tri(Sigmahat)]
    if (any(abs(cors) > 1)) {
      lav_msg_warn(gettext(
        "some model-implied correlations are larger than 1.0"))
    }
    nvar <- nrow(Sigmahat)
    th.idx <- TH.IDX[[g]]
    MEAN <- rep(0, nvar)

    # reconstruct ov.types
    ov.types <- rep("numeric", nvar)
    ord.idx <- unique(th.idx[th.idx > 0])
    ov.types[ord.idx] <- "ordered"

    if (all(ov.types == "ordered")) {
      # get patterns ## FIXME GET it
      if (!is.null(lavdata@Rp[[g]]$pat)) {
        PAT <- lavdata@Rp[[g]]$pat
      } else {
        PAT <- lav_data_resp_patterns(lavdata@X[[g]])$pat
      }
      npatterns <- nrow(PAT)
      freq <- as.numeric(rownames(PAT))
      PI.group <- numeric(npatterns)
      TH.VAR <- lapply(
        1:nvar,
        function(x) c(-Inf, TH[[g]][th.idx == x], +Inf)
      )
      # FIXME!!! ok to set diagonal to 1.0?
      diag(Sigmahat) <- 1.0
      for (r in 1:npatterns) {
        # compute probability for each pattern
        lower <- sapply(1:nvar, function(x) TH.VAR[[x]][PAT[r, x]])
        upper <- sapply(1:nvar, function(x) TH.VAR[[x]][PAT[r, x] + 1L])
        # handle missing values
        na.idx <- which(is.na(PAT[r, ]))
        if (length(na.idx) > 0L) {
          lower <- lower[-na.idx]
          upper <- upper[-na.idx]
          MEAN.r <- MEAN[-na.idx]
          Sigmahat.r <- Sigmahat[-na.idx, -na.idx, drop = FALSE]
        } else {
          MEAN.r <- MEAN
          Sigmahat.r <- Sigmahat
        }
        PI.group[r] <- sadmvn(lower, upper,
          mean = MEAN.r,
          varcov = Sigmahat.r
        )
      }
    } else { # case-wise
      PI.group <- rep(as.numeric(NA), lavdata@nobs[[g]])
      lav_msg_warn(gettext("casewise PI not implemented"))
    }

    PI[[g]] <- PI.group
  } # g

  PI
}

lav_tables_table_format <- function(out, lavdata = lavdata,
                                    lavobject = lavobject) {
  # determine column we need
  NAMES <- names(out)
  stat.idx <- which(NAMES %in% c(
    "cor", "cor.un",
    "G2", "G2.un",
    "X2", "X2.un",
    "RMSEA", "RMSEA.un",
    "G2.average", "G2.plarge", "G2.nlarge",
    "X2.average", "X2.plarge", "X2.nlarge"
  ))
  if (length(stat.idx) == 0) {
    if (!is.null(out$obs.freq)) {
      stat.idx <- which(NAMES == "obs.freq")
    } else if (!is.null(out$nobs)) {
      stat.idx <- which(NAMES == "nobs")
    }
    UNI <- NULL
  } else if (length(stat.idx) > 1) {
    lav_msg_stop(gettext(
      "more than one statistic for table output:"),
      paste(NAMES[stat.idx], collapse = " ")
    )
  } else {
    # univariate version of same statistic
    if (NAMES[stat.idx] == "G2.average") {
      UNI <- lavTables(lavobject, dimension = 1L, statistic = "G2")
    } else if (NAMES[stat.idx] == "X2.average") {
      UNI <- lavTables(lavobject, dimension = 1L, statistic = "X2")
    } else {
      UNI <- NULL
    }
  }

  OUT <- vector("list", length = lavdata@ngroups)
  for (g in 1:lavdata@ngroups) {
    if (lavdata@ngroups == 1L) { # no group column
      STAT <- out[[stat.idx]]
    } else {
      STAT <- out[[stat.idx]][out$group == g]
    }
    RN <- lavdata@ov.names[[g]]
    OUT[[g]] <- getCov(STAT, diagonal = FALSE, lower = FALSE, names = RN)
    # change diagonal elements: replace by univariate stat
    # if possible
    diag(OUT[[g]]) <- as.numeric(NA)
    if (!is.null(UNI)) {
      if (!is.null(UNI$group)) {
        idx <- which(UNI$group == g)
      } else {
        idx <- 1:length(UNI$lhs)
      }
      if (NAMES[stat.idx] == "G2.average") {
        diag(OUT[[g]]) <- tapply(UNI$G2[idx],
          INDEX = UNI$id[idx],
          FUN = mean
        )
      } else if (NAMES[stat.idx] == "X2.average") {
        diag(OUT[[g]]) <- tapply(UNI$X2[idx],
          INDEX = UNI$id[idx],
          FUN = mean
        )
      }
    } else if (NAMES[stat.idx] %in% c("cor", "cor.un")) {
      diag(OUT[[g]]) <- 1
    }
    class(OUT[[g]]) <- c("lavaan.matrix.symmetric", "matrix")
  }
  if (lavdata@ngroups > 1L) {
    names(OUT) <- lavdata@group.label
    out <- OUT
  } else {
    out <- OUT[[1]]
  }

  out
}

lav_tables_cells_format <- function(out, lavdata = lavdata,
                                    drop.list.single.group = FALSE) {
  OUT <- vector("list", length = lavdata@ngroups)
  if (is.null(out$group)) {
    out$group <- rep(1L, length(out$lhs))
  }
  # do we have a statistic?
  # determine column we need
  NAMES <- names(out)
  stat.idx <- which(NAMES %in% c(
    "cor", "cor.un",
    "G2", "G2.un",
    "X2", "X2.un",
    "RMSEA", "RMSEA.un",
    "G2.average", "G2.plarge", "G2.nlarge",
    "X2.average", "X2.plarge", "X2.nlarge"
  ))
  if (length(stat.idx) == 0) {
    statistic <- "obs.freq"
  } else if (length(stat.idx) > 1) {
    lav_msg_stop(gettext(
      "more than one statistic for table output:"),
      paste(NAMES[stat.idx], collapse = " ")
    )
  } else {
    statistic <- NAMES[stat.idx]
  }

  for (g in 1:lavdata@ngroups) {
    case.idx <- which(out$group == g)
    ID.group <- unique(out$id[out$group == g])
    TMP <- lapply(ID.group, function(x) {
      Tx <- out[out$id == x, ]
      M <- matrix(
        Tx[, statistic],
        max(Tx$row), max(Tx$col)
      )
      rownames(M) <- unique(Tx$row)
      colnames(M) <- unique(Tx$col)
      class(M) <- c("lavaan.matrix", "matrix")
      M
    })
    names(TMP) <- unique(paste(out$lhs[case.idx], out$rhs[case.idx],
      sep = "_"
    ))
    OUT[[g]] <- TMP
  }

  if (lavdata@ngroups == 1L && drop.list.single.group) {
    OUT <- OUT[[1]]
  } else {
    if (length(lavdata@group.label) > 0L) {
      names(OUT) <- unlist(lavdata@group.label)
    }
  }

  OUT
}
