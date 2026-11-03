# compare the measurement parameters of a multiple-group model across the
# groups: one table per factor (loadings: indicators by groups), or one table
# of intercepts/thresholds (indicators by groups), with the range of every
# parameter over the groups, so that (possibly) non-invariant items stand
# out. Meant for rotated (ESEM, with or without the agreement rotation) and
# aligned models, but works for any multiple-group model.
#
# YR 03 Nov 2026: first version

lavCompareGroups <- function(object, what = "loadings", # nolint
                             standardized = FALSE, add_range = TRUE,
                             sort = FALSE) {
  if (!inherits(object, "lavaan")) {
    lav_msg_stop(gettext("object is not a lavaan object."))
  }
  ngroups <- object@Data@ngroups
  if (ngroups < 2L) {
    lav_msg_stop(gettext("lavCompareGroups() requires a multiple-group model."))
  }
  if (object@Data@nlevels > 1L) {
    lav_msg_stop(gettext("lavCompareGroups() is not available for multilevel
                          models."))
  }
  what <- tolower(what)
  if (length(what) != 1L || !what %in% c("loadings", "intercepts")) {
    lav_msg_stop(gettext("what must be \"loadings\" or \"intercepts\"."))
  }
  group_label <- object@Data@group.label
  if (length(group_label) == 0L) {
    group_label <- as.character(seq_len(ngroups))
  }

  # the estimates: the model metric, or completely standardized
  pe <- parameterEstimates(object, se = FALSE, zstat = FALSE, pvalue = FALSE,
                           ci = FALSE, standardized = standardized,
                           remove.system.eq = FALSE, remove.eq = FALSE,
                           remove.ineq = FALSE, remove.def = TRUE)
  est_col <- if (standardized) "std.all" else "est"
  if (!est_col %in% names(pe)) {
    lav_msg_stop(gettextf("column %s is not available.", dQuote(est_col)))
  }

  # one table: parameters selected by op (and lhs for the loadings), rows
  # in the order of the first group, columns the groups
  make_table <- function(op, lhs = NULL) {
    idx1 <- which(pe$op == op & pe$group == 1L &
                  (if (is.null(lhs)) TRUE else pe$lhs == lhs))
    if (length(idx1) == 0L) {
      return(NULL)
    }
    row_label <- if (op == "=~") {
      pe$rhs[idx1]
    } else if (op == "|") {
      paste0(pe$lhs[idx1], "|", pe$rhs[idx1])
    } else {
      pe$lhs[idx1]
    }
    out <- matrix(NA_real_, length(idx1), ngroups,
                  dimnames = list(row_label, group_label))
    for (g in seq_len(ngroups)) {
      for (i in seq_along(idx1)) {
        idx_g <- which(pe$lhs == pe$lhs[idx1[i]] & pe$op == op &
                       pe$rhs == pe$rhs[idx1[i]] & pe$group == g)
        if (length(idx_g) == 1L) {
          out[i, g] <- pe[[est_col]][idx_g]
        }
      }
    }
    if (add_range) {
      rng <- apply(out, 1L, function(x) {
        if (all(is.na(x))) NA_real_ else diff(range(x, na.rm = TRUE))
      })
      out <- cbind(out, range = rng)
    }
    if (sort) {
      rng <- apply(out[, seq_len(ngroups), drop = FALSE], 1L, function(x) {
        if (all(is.na(x))) -Inf else diff(range(x, na.rm = TRUE))
      })
      out <- out[order(rng, decreasing = TRUE), , drop = FALSE]
    }
    # a data frame (the group labels are printed in full)
    col_names <- colnames(out)
    out <- as.data.frame(out)
    names(out) <- col_names
    class(out) <- c("lavaan.data.frame", "data.frame")
    out
  }

  if (what == "loadings") {
    lv_names <- unique(pe$lhs[pe$op == "=~"])
    out <- lapply(lv_names, function(lv) make_table("=~", lhs = lv))
    names(out) <- lv_names
    out <- out[!vapply(out, is.null, logical(1L))]
  } else {
    # the intercepts of the (non-ordered) indicators and the thresholds of
    # the ordered ones; the intercepts of the factors are left out
    ov_names <- unique(pe$rhs[pe$op == "=~"])
    nu <- make_table("~1")
    if (!is.null(nu)) {
      nu <- nu[rownames(nu) %in% ov_names, , drop = FALSE]
      ordered_ov <- unique(pe$lhs[pe$op == "|"])
      nu <- nu[!rownames(nu) %in% ordered_ov, , drop = FALSE]
      if (nrow(nu) == 0L) {
        nu <- NULL
      }
    }
    th <- make_table("|")
    out <- rbind(nu, th)
    if (is.null(out)) {
      lav_msg_stop(gettext("the model has no intercepts or thresholds."))
    }
    class(out) <- c("lavaan.data.frame", "data.frame")
  }
  out
}
