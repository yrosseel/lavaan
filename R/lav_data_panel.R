# ------------------------ long panel data ----------------------------------- #
# YR 17 Oct 2026: panel = c(id = "pid", time = "wave")
#
# Long data (one row per person and occasion) are reshaped to wide data (one
# row per person), so that a model for repeated measures can be written for
# the wide variables, by hand (x1, x2, x3) or with indexed names (x[t in 1:3]):
# the variable x at time 3 becomes the column x3 (or x_3 with sep = "_").
#
# - the time variable must contain non-negative integers; missing occasions
#   become missing values in the wide data (handled by missing = "ml")
# - only the wide variables that the model uses are created; a model name
#   matches x<sep><time>, with or without zero padding (x3, x03)
# - columns that are constant within each id (time-invariant covariates,
#   group, cluster, sampling weights) are kept as a single column
# - base names in ordered= and aux= stand for all their wide variables
# ---------------------------------------------------------------------------- #

lav_data_panel <- function(data = NULL, panel = NULL, flat_model = NULL,
                           ordered = NULL, aux = NULL, group = NULL,
                           cluster = NULL, sampling_weights = NULL) {
  spec <- lav_data_panel_spec(panel)
  temp <- lav_data_panel_check(data, spec)
  data <- temp$data
  id <- temp$id
  time_int <- temp$time
  id_name <- spec$id
  time_name <- spec$time
  sep <- spec$sep

  # one row per id: positions of the rows per (id, time)
  ids <- unique(id)
  times <- sort(unique(time_int))
  id_idx <- match(id, ids)
  time_idx <- match(time_int, times)
  key <- (id_idx - 1L) * length(times) + time_idx
  if (anyDuplicated(key)) {
    dup <- which(duplicated(key))[1L]
    lav_msg_stop(gettextf(
      "panel=: the data contain more than one row for id %1$s at time %2$s;
      every id can have at most one row per occasion.",
      dQuote(as.character(id[dup]), FALSE), time_int[dup]))
  }
  row_of <- matrix(NA_integer_, nrow = length(ids), ncol = length(times))
  row_of[cbind(id_idx, time_idx)] <- seq_along(id)

  # which columns are constant within each id? (missing values ignored)
  other <- setdiff(names(data), c(id_name, time_name))
  first_row <- vapply(seq_along(ids), function(i) {
    r <- row_of[i, ]
    r[!is.na(r)][1L]
  }, integer(1L))
  constant <- vapply(other, function(v) {
    x <- data[[v]]
    lav_data_panel_constant(x, id_idx, length(ids))
  }, logical(1L))

  # the names used by the model
  model_names <- lav_data_panel_model_names(flat_model)

  # candidate wide names: column x, time t -> x<sep>t (also zero padded)
  width <- max(nchar(as.character(times)))
  cand_name <- character(0L)
  cand_var <- character(0L)
  cand_time <- integer(0L)
  for (v in other) {
    plain <- paste0(v, sep, times)
    padded <- paste0(v, sep, formatC(times, width = width, flag = "0"))
    cand_name <- c(cand_name, plain, padded)
    cand_var <- c(cand_var, rep(v, 2L * length(times)))
    cand_time <- c(cand_time, seq_along(times), seq_along(times))
  }
  keep <- !duplicated(cand_name)
  cand_name <- cand_name[keep]
  cand_var <- cand_var[keep]
  cand_time <- cand_time[keep]

  # expand base names in ordered= and aux= (y -> y1, y2, ...)
  expand_base <- function(x) {
    if (is.null(x) || is.logical(x)) return(x)
    out <- character(0L)
    for (v in x) {
      if (v %in% other && !constant[[v]]) {
        out <- c(out, paste0(v, sep, times))
      } else {
        out <- c(out, v)
      }
    }
    unique(out)
  }
  ordered <- expand_base(ordered)
  aux <- expand_base(aux)

  wanted <- unique(c(model_names, if (is.character(ordered)) ordered,
                     if (is.character(aux)) aux))

  # a model name that is a time-varying column of the long data
  varying_used <- intersect(wanted, other[!constant])
  if (length(varying_used) > 0L) {
    v <- varying_used[1L]
    lav_msg_stop(gettextf(
      "panel=: the variable %1$s changes over time within an id; in the
      model, refer to it at a given time, as %2$s (or with an index, as
      %3$s[t in %4$s:%5$s]).", dQuote(v, FALSE),
      paste0(v, sep, times[1L]), paste0(v, sep), times[1L],
      times[length(times)]))
  }
  for (arg in list(list("group", group), list("cluster", cluster),
                   list("sampling_weights", sampling_weights))) {
    v <- arg[[2L]]
    if (is.character(v) && length(v) == 1L && v %in% other && !constant[[v]]) {
      lav_msg_stop(gettextf(
        "panel=: the %1$s variable %2$s must be constant within each id.",
        arg[[1L]], dQuote(v, FALSE)))
    }
  }

  # build the wide data
  wide <- data.frame(row.names = seq_along(ids))
  wide[[id_name]] <- ids
  for (v in other[constant]) {
    wide[[v]] <- data[[v]][first_row]
  }
  wide_names <- intersect(wanted, cand_name)
  ambiguous <- intersect(wide_names, other[constant])
  if (length(ambiguous) > 0L) {
    v <- ambiguous[1L]
    k <- match(v, cand_name)
    lav_msg_stop(gettextf(
      "panel=: the name %1$s is both a column of the data and the variable
      %2$s at time %3$s; use another separator, for example
      panel = c(id = \"%4$s\", time = \"%5$s\", sep = \"_\").",
      dQuote(v, FALSE), dQuote(cand_var[k], FALSE), times[cand_time[k]],
      id_name, time_name))
  }
  for (nm in wide_names) {
    k <- match(nm, cand_name)
    rows <- row_of[, cand_time[k]]
    wide[[nm]] <- data[[cand_var[k]]][rows]
  }

  list(data = wide, ordered = ordered, aux = aux)
}

# the panel= argument: c(id = , time = ) or list(id = , time = , sep = )
lav_data_panel_spec <- function(panel) {
  if (is.list(panel)) panel <- unlist(panel)
  if (!is.character(panel) || length(panel) < 2L) {
    lav_msg_stop(gettext(
      "panel= should be a character vector with (at least) the names of the
      id and the time variable, as in panel = c(id = \"pid\", time = \"wave\")
      (optionally with sep = \"_\")."))
  }
  nms <- names(panel)
  if (is.null(nms) || all(nms == "")) {
    nms <- c("id", "time", "sep")[seq_along(panel)]
  }
  nms <- tolower(nms)
  bad <- setdiff(nms, c("id", "time", "sep"))
  if (length(bad) > 0L || anyDuplicated(nms) ||
      !all(c("id", "time") %in% nms)) {
    lav_msg_stop(gettext(
      "panel= should contain the elements id and time (and optionally sep),
      as in panel = c(id = \"pid\", time = \"wave\")."))
  }
  names(panel) <- nms
  sep <- if ("sep" %in% nms) panel[["sep"]] else ""
  list(id = panel[["id"]], time = panel[["time"]], sep = sep)
}

# is x constant within each id (ignoring missing values)?
lav_data_panel_constant <- function(x, id_idx, nid) {
  ok <- !is.na(x)
  if (!any(ok)) return(TRUE)
  xx <- if (is.factor(x)) as.integer(x) else x
  first <- rep(NA, nid)
  pos <- which(ok)
  first_pos <- pos[!duplicated(id_idx[pos])]
  first[id_idx[first_pos]] <- xx[first_pos]
  all(xx[pos] == first[id_idx[pos]])
}

# the (observed) names a flat model refers to, interaction parts included
lav_data_panel_model_names <- function(flat_model) {
  if (is.null(flat_model) || is.null(flat_model$lhs)) {
    return(character(0L))
  }
  rel <- !flat_model$op %in% c(":", "==", "<", ">", ":=", ":~")
  nms <- c(flat_model$lhs[rel], flat_model$rhs[rel])
  nms <- unlist(strsplit(nms, ":", fixed = TRUE))
  unique(nms[nzchar(nms)])
}

# check the id and time variables of panel=; remove rows with a missing id
# or time; returns the data and the id and (integer) time vectors
lav_data_panel_check <- function(data = NULL, spec = NULL) {
  if (is.null(data) || !is.data.frame(data)) {
    lav_msg_stop(gettext(
      "the panel= argument needs (long) raw data in the data= argument."))
  }
  id_name <- spec$id
  time_name <- spec$time
  for (v in c(id_name, time_name)) {
    if (!v %in% names(data)) {
      lav_msg_stop(gettextf(
        "panel= variable %s not found in the data.", dQuote(v, FALSE)))
    }
  }

  # rows without id or time are useless
  id <- data[[id_name]]
  time <- data[[time_name]]
  if (is.factor(time)) time <- as.character(time)
  time_num <- suppressWarnings(as.numeric(time))
  bad_time <- !is.na(time) & (is.na(time_num) | time_num < 0 |
                                time_num != round(time_num))
  if (any(bad_time)) {
    lav_msg_stop(gettextf(
      "the time variable %1$s of panel= must contain non-negative integers
      (found %3$s); recode it first, for example with
      as.integer(factor(%2$s)).", dQuote(time_name, FALSE), time_name,
      dQuote(as.character(time[which(bad_time)[1L]]), FALSE)))
  }
  drop_idx <- which(is.na(id) | is.na(time_num))
  if (length(drop_idx) > 0L) {
    lav_msg_warn(gettextf(
      "%1$s row(s) of the data have a missing value for the id or the time
      variable of panel=, and are removed.", length(drop_idx)))
    data <- data[-drop_idx, , drop = FALSE]
    id <- id[-drop_idx]
    time_num <- time_num[-drop_idx]
  }
  if (nrow(data) == 0L) {
    lav_msg_stop(gettext("no rows left in the data for panel=."))
  }
  list(data = data, id = id, time = as.integer(time_num))
}

# ------------------------ long panel data, two-level ------------------------ #
# panel= in a two-level model: the data stay long (one row per occasion),
# the id variable is the cluster variable, and the lagged variables of the
# model (y[t-1] in the syntax, y.lag1 in the parameter table; see
# lav_expand_index_long) are created: y.lag1 is the value of y at time - 1
# for the same id, and missing if there is no such occasion (no row, or a
# missing value). The lags follow the values of the time variable, so gaps
# in the measurement occasions are respected.
lav_data_panel_long <- function(data = NULL, panel = NULL, lagged = NULL,
                                cluster = NULL, flat_model = NULL) {
  spec <- lav_data_panel_spec(panel)
  temp <- lav_data_panel_check(data, spec)
  data <- temp$data
  id <- temp$id
  time_int <- temp$time
  id_name <- spec$id

  # the cluster variable is the id variable
  if (!is.null(cluster) && !identical(cluster, id_name)) {
    lav_msg_stop(gettextf(
      "panel=: in a two-level model, the cluster variable is the id
      variable %1$s (found cluster = %2$s).", dQuote(id_name, FALSE),
      dQuote(paste(cluster, collapse = " "), FALSE)))
  }

  # one row per (id, time)
  key <- paste(match(id, unique(id)), time_int, sep = ":")
  if (anyDuplicated(key)) {
    dup <- which(duplicated(key))[1L]
    lav_msg_stop(gettextf(
      "panel=: the data contain more than one row for id %1$s at time %2$s;
      every id can have at most one row per occasion.",
      dQuote(as.character(id[dup]), FALSE), time_int[dup]))
  }

  # the lagged variables
  for (i in seq_len(NROW(lagged))) {
    nm <- lagged$name[i]
    v <- lagged$var[i]
    if (nm %in% names(data)) {
      lav_msg_stop(gettextf(
        "panel=: the data already contain a variable %1$s; this name is
        used for %2$s at an earlier (or later) occasion.",
        dQuote(nm, FALSE), dQuote(v, FALSE)))
    }
    if (!v %in% names(data)) {
      lav_msg_stop(gettextf(
        "panel=: the lagged variable %s is not found in the data.",
        dQuote(v, FALSE)))
    }
    other_key <- paste(match(id, unique(id)), time_int - lagged$lag[i],
                       sep = ":")
    data[[nm]] <- data[[v]][match(other_key, key)]
  }

  # occasions without the value of a lagged predictor (the first occasion
  # of each id, an occasion after a gap) can not be used: the likelihood is
  # conditional on the predictors (Mplus drops them silently as well);
  # remove them here, rather than with a warning about missing values in
  # exogenous variables
  if (NROW(lagged) > 0L) {
    endo <- character(0L)
    if (!is.null(flat_model)) {
      endo <- c(flat_model$lhs[flat_model$op == "~"],
                flat_model$rhs[flat_model$op == "=~"])
    }
    pred <- setdiff(unique(lagged$name), endo)
    if (length(pred) > 0L) {
      keep <- stats::complete.cases(data[, pred, drop = FALSE])
      data <- data[keep, , drop = FALSE]
    }
  }

  list(data = data, cluster = id_name)
}
