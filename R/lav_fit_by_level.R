# level-specific fit measures for multilevel models, using the
# 'partially saturated model' approach (Ryu & West, 2009):
#
# - to assess the fit of the WITHIN part of the model, we fit an auxiliary
#   model where the between part is saturated (and vice versa); any misfit
#   of this partially saturated model can then be attributed to the
#   within part only
# - for the incremental fit indices (CFI/TLI), we also need level-specific
#   baseline models: the independence model at the target level, again
#   combined with a saturated model at the other level
#
# the (fitted) partially saturated models are stored in the @internal slot
# (internal$fit.by.level), and are used by fitMeasures(..., level = ) and
# by summary(..., fit.measures = TRUE)
#
# YR 11 July 2026 (initial version)

# generate 'saturated' rows for all blocks at the other level(s)
# (level != level.keep)
#
# important: any exogenous (fixed.x) covariates at the saturated level
# must REMAIN exogenous, because we re-use the sample statistics (and h1)
# slots of the original model: the saturated part is therefore saturated
# *conditionally* on the covariates (free y-(co)variances, free y~x
# slopes), just as if the user had written out a saturated level in the
# model syntax
lav_fit_by_level_sat_rows <- function(lavobject, level.keep = 1L) {
  lavdata <- lavobject@Data
  lavpta <- lavobject@pta
  lavh1 <- lavobject@h1
  nlevels <- lavdata@nlevels
  ngroups <- lavdata@ngroups
  fixed_x <- lavobject@Options$fixed.x

  lhs <- op <- rhs <- character(0L)
  block <- group <- level <- free <- exo <- integer(0L)
  ustart <- numeric(0L)

  add_rows <- function(lhs1, op1, rhs1, b, g, l, free1, exo1, ustart1) {
    n <- length(lhs1)
    lhs <<- c(lhs, lhs1)
    op <<- c(op, rep(op1, n))
    rhs <<- c(rhs, rhs1)
    block <<- c(block, rep(b, n))
    group <<- c(group, rep(g, n))
    level <<- c(level, rep(l, n))
    free <<- c(free, rep(free1, length.out = n))
    exo <<- c(exo, rep(exo1, length.out = n))
    ustart <<- c(ustart, rep(ustart1, length.out = n))
    invisible(NULL)
  }

  for (g in seq_len(ngroups)) {
    for (l in seq_len(nlevels)) {
      if (l == level.keep) {
        next
      }
      b <- (g - 1L) * nlevels + l

      ov_names <- lavpta$vnames$ov[[b]]
      if (fixed_x) {
        ov_x <- lavpta$vnames$ov.x[[b]]
      } else {
        ov_x <- character(0L)
      }
      ov_nox <- ov_names[!ov_names %in% ov_x]

      # h1 (unrestricted) moments of this block, for starting values
      sample_cov <- NULL
      sample_mean <- NULL
      if (length(lavh1) > 0L && !is.null(lavh1$implied$cov[[b]])) {
        sample_cov <- lavh1$implied$cov[[b]]
        sample_mean <- lavh1$implied$mean[[b]]
      }
      s_val <- function(v1, v2) {
        if (is.null(sample_cov)) {
          return(as.numeric(NA))
        }
        sample_cov[cbind(match(v1, ov_names), match(v2, ov_names))]
      }
      m_val <- function(v1) {
        if (is.null(sample_mean)) {
          return(as.numeric(NA))
        }
        sample_mean[match(v1, ov_names)]
      }

      # 1. variances of ALL observed variables first, in the same order
      # as lavpta$vnames$ov[[b]]: the first-appearance order of the ov's
      # in the parameter table must match the original model, because the
      # sample statistics (which are re-used) are stored in that order
      is_x <- ov_names %in% ov_x
      add_rows(ov_names, "~~", ov_names, b, g, l,
               free1 = ifelse(is_x, 0L, 1L),
               exo1 = ifelse(is_x, 1L, 0L),
               ustart1 = s_val(ov_names, ov_names))

      # 2. free covariances among the non-exogenous variables
      if (length(ov_nox) > 1L) {
        tmp <- utils::combn(ov_nox, 2L)
        add_rows(tmp[1L, ], "~~", tmp[2L, ], b, g, l,
                 free1 = 1L, exo1 = 0L,
                 ustart1 = s_val(tmp[1L, ], tmp[2L, ]))
      }

      # 3. free slopes of all y's on all (fixed) x's, plus the fixed
      # covariances among the exogenous variables
      if (length(ov_x) > 0L) {
        beta_start <- try({
          s_xx <- sample_cov[match(ov_x, ov_names),
                             match(ov_x, ov_names), drop = FALSE]
          s_xy <- sample_cov[match(ov_x, ov_names),
                             match(ov_nox, ov_names), drop = FALSE]
          solve(s_xx, s_xy) # nx x ny
        }, silent = TRUE)
        for (x in ov_x) {
          st <- rep(as.numeric(NA), length(ov_nox))
          if (!inherits(beta_start, "try-error") &&
              !is.null(sample_cov)) {
            st <- beta_start[match(x, ov_x), ]
          }
          add_rows(ov_nox, "~", rep(x, length(ov_nox)), b, g, l,
                   free1 = 1L, exo1 = 0L, ustart1 = st)
        }

        if (length(ov_x) > 1L) {
          tmp <- utils::combn(ov_x, 2L)
          add_rows(tmp[1L, ], "~~", tmp[2L, ], b, g, l,
                   free1 = 0L, exo1 = 1L,
                   ustart1 = s_val(tmp[1L, ], tmp[2L, ]))
        }
      }

      # 4. intercepts
      if (lavobject@Options$meanstructure) {
        if (l == 1L) {
          # within level: intercepts fixed to zero, except for
          # within-only variables (the same convention as lavaanify)
          b2 <- (g - 1L) * nlevels + 2L
          ov_l2 <- lavpta$vnames$ov[[b2]]
          within_only <- !ov_nox %in% ov_l2
          add_rows(ov_nox, "~1", rep("", length(ov_nox)), b, g, l,
                   free1 = ifelse(within_only, 1L, 0L), exo1 = 0L,
                   ustart1 = ifelse(within_only, m_val(ov_nox), 0))
        } else {
          add_rows(ov_nox, "~1", rep("", length(ov_nox)), b, g, l,
                   free1 = 1L, exo1 = 0L, ustart1 = m_val(ov_nox))
        }
        if (length(ov_x) > 0L) {
          add_rows(ov_x, "~1", rep("", length(ov_x)), b, g, l,
                   free1 = 0L, exo1 = 1L, ustart1 = m_val(ov_x))
        }
      }
    } # levels
  } # groups

  data.frame(
    lhs = lhs, op = op, rhs = rhs,
    user = rep(1L, length(lhs)),
    block = block, group = group, level = level,
    free = free, ustart = ustart, exo = exo,
    stringsAsFactors = FALSE
  )
}

# construct the parameter table for a partially saturated model:
# - the rows of the target level (level.keep) are taken from the user model
#   (or from the independence model, if independent = TRUE), while
# - the other level(s) are saturated
#
# returns NULL if the model contains constraints involving parameters from
# different levels (in that case, the misfit cannot be attributed to a
# single level)
lav_pt_partial_saturation <- function(lavobject, level.keep = 1L,
                                      independent = FALSE) {
  stopifnot(inherits(lavobject, "lavaan"))
  pt_user <- as.data.frame(lavobject@ParTable, stringsAsFactors = FALSE)

  # rows for the target level
  if (independent) {
    pt_keep <- as.data.frame(
      lav_pt_independence(lavobject = lavobject),
      stringsAsFactors = FALSE
    )
  } else {
    pt_keep <- pt_user
    # use the fitted values as starting values
    if (!is.null(pt_keep$est)) {
      pt_keep$ustart <- pt_keep$est
    }
  }

  # saturated rows for the other level(s); note: exogenous (fixed.x)
  # covariates remain exogenous (see lav_fit_by_level_sat_rows)
  pt_sat <- lav_fit_by_level_sat_rows(lavobject, level.keep = level.keep)

  # the columns we retain in the new parameter table
  core <- c(
    "lhs", "op", "rhs", "user", "block", "group", "level",
    "free", "ustart", "exo", "label", "plabel"
  )
  chr_core <- c("lhs", "op", "rhs", "label", "plabel")
  for (cc in core) {
    if (is.null(pt_keep[[cc]])) {
      if (cc == "group") {
        # single group: all (block > 0) rows belong to group 1
        pt_keep$group <- ifelse(pt_keep$block > 0L, 1L, 0L)
      } else {
        pt_keep[[cc]] <- if (cc %in% chr_core) "" else 0L
      }
    }
    if (is.null(pt_sat[[cc]])) {
      if (cc == "group") {
        pt_sat$group <- ifelse(pt_sat$block > 0L, 1L, 0L)
      } else {
        pt_sat[[cc]] <- if (cc %in% chr_core) "" else 0L
      }
    }
  }

  # the user parameter table may use labels in the group and level
  # columns (e.g., "g1", "within"), while the generated saturated and
  # independence partables use integer indices; map the indices to the
  # user values, so the combined partable is consistent
  group_values <- lav_pt_group_values(pt_user)
  level_values <- lav_pt_level_values(pt_user)
  level_keep_val <- level_values[level.keep]
  pt_sat$group <- group_values[pt_sat$group]
  pt_sat$level <- level_values[pt_sat$level]
  if (independent) {
    pt_keep$group <- group_values[pt_keep$group]
    pt_keep$level <- level_values[pt_keep$level]
  }

  con_ops <- c("==", "<", ">", ":=")
  keep_rows <- pt_keep[pt_keep$level %in% level_keep_val &
                       !pt_keep$op %in% con_ops, core, drop = FALSE]
  sat_rows <- pt_sat[!pt_sat$level %in% level_keep_val, core, drop = FALSE]

  # handle constraint/definition rows of the user model (they carry
  # level = 0): keep them if they only involve parameters of the target
  # level; silently drop them if they only involve parameters of the
  # other level(s); if they mix both, level-specific fit is not defined
  con_rows <- NULL
  if (!independent) {
    con_idx <- which(pt_user$op %in% con_ops)
    if (length(con_idx) > 0L) {
      mod_flag <- !pt_user$op %in% con_ops
      lab_all <- character(0L)
      if (!is.null(pt_user$plabel)) {
        lab_all <- pt_user$plabel
      } else {
        lab_all <- rep("", nrow(pt_user))
      }
      lab_user <- pt_user$label
      if (is.null(lab_user)) {
        lab_user <- rep("", nrow(pt_user))
      }
      keep_flag <- mod_flag & pt_user$level %in% level_keep_val
      drop_flag <- mod_flag & !pt_user$level %in% level_keep_val
      keep_syms <- c(lab_all[keep_flag], lab_user[keep_flag])
      drop_syms <- c(lab_all[drop_flag], lab_user[drop_flag])
      keep_syms <- keep_syms[nchar(keep_syms) > 0L]
      drop_syms <- drop_syms[nchar(drop_syms) > 0L]
      # a label may be attached to parameters at both levels (implying a
      # cross-level equality constraint)
      if (any(keep_syms %in% drop_syms)) {
        return(NULL)
      }

      keep_con <- logical(length(con_idx))
      for (ii in seq_along(con_idx)) {
        i <- con_idx[ii]
        refs <- character(0L)
        for (txt in c(pt_user$lhs[i], pt_user$rhs[i])) {
          if (nchar(txt) > 0L && is.na(suppressWarnings(as.numeric(txt)))) {
            refs <- c(refs, all.vars(parse(
              text = txt,
              keep.source = FALSE
            )))
          }
        }
        in_keep <- refs %in% keep_syms
        in_drop <- refs %in% drop_syms
        if (pt_user$op[i] == ":=") {
          if (length(refs) > 0L && all(in_keep)) {
            keep_con[ii] <- TRUE
            # the defined parameter becomes available for later rows
            keep_syms <- c(keep_syms, pt_user$lhs[i])
          } else {
            drop_syms <- c(drop_syms, pt_user$lhs[i])
          }
        } else {
          if (any(in_keep) && any(in_drop)) {
            # cross-level constraint
            return(NULL)
          }
          if (any(in_keep) && !any(in_drop)) {
            keep_con[ii] <- TRUE
          }
        }
      }
      if (any(keep_con)) {
        # pt_keep == pt_user here (independent = FALSE), but with all
        # core columns present
        con_rows <- pt_keep[con_idx[keep_con], core, drop = FALSE]
      }
    }
  }

  # combine, per block; constraint rows at the end
  pt_new <- rbind(keep_rows, sat_rows)
  pt_new <- pt_new[order(pt_new$block), , drop = FALSE]
  if (!is.null(con_rows)) {
    pt_new <- rbind(pt_new, con_rows)
  }

  # renumber the free parameters
  idx_free <- which(pt_new$free > 0L)
  pt_new$free <- 0L
  pt_new$free[idx_free] <- seq_along(idx_free)
  pt_new$id <- seq_len(nrow(pt_new))
  rownames(pt_new) <- NULL

  as.list(pt_new)
}

# fit the partially saturated models (one per level), together with their
# level-specific baseline models; returns a (possibly empty) named list
# with one fitted lavaan object per level; each object's @baseline slot
# contains the level-specific baseline, so that (for example) CFI/TLI are
# automatically computed with respect to the correct baseline model
lav_object_fit_by_level <- function(object) {
  out <- list()
  lavdata <- object@Data
  if (lavdata@nlevels < 2L) {
    return(out)
  }
  lavoptions <- object@Options

  # we need a converged fit and a test statistic
  if (!lavoptions$do.fit || !object@optim$converged ||
      object@optim$npar == 0L) {
    return(out)
  }
  if (length(lavoptions$test) == 1L && lavoptions$test == "none") {
    return(out)
  }
  if (length(object@h1) == 0L) {
    return(out)
  }
  # not (yet) for random slopes or efa blocks
  if (!is.null(object@ParTable$rv) &&
      any(nchar(object@ParTable$rv) > 0L)) {
    return(out)
  }
  if (!is.null(object@ParTable$efa) &&
      any(nchar(object@ParTable$efa) > 0L)) {
    return(out)
  }
  # not (yet) for categorical, conditional.x, correlation structures,
  # or free group weights
  if (object@Model@categorical || lavoptions$conditional.x ||
      object@Model@correlation || lavoptions$group.w.free) {
    return(out)
  }

  # options for the auxiliary fits
  lavoptions$se <- "none"
  lavoptions$h1 <- FALSE # already provided by sloth1
  lavoptions$baseline <- FALSE # we construct our own baselines
  lavoptions$fit.by.level <- FALSE # no recursion
  lavoptions$loglik <- TRUE
  lavoptions$implied <- TRUE
  lavoptions$check.start <- FALSE
  lavoptions$check.gradient <- FALSE
  lavoptions$check.post <- FALSE
  lavoptions$check.vcov <- FALSE
  lavoptions$bounds <- "none"
  lavoptions$optim.bounds <- list()
  lavoptions$start <- "default" # the ustart values are already good
  lavoptions$rstarts <- 0L
  lavoptions$do.fit <- TRUE
  if (lavoptions$optim.method != "noniter") {
    lavoptions$optim.method <- "nlminb"
  }

  level_names <- c("within", "between")

  for (l in seq_len(lavdata@nlevels)) {
    # partially saturated model for this level
    pt_ps <- lav_pt_partial_saturation(object,
      level.keep = l,
      independent = FALSE
    )
    if (is.null(pt_ps)) {
      lav_msg_warn(gettext(
        "level-specific fit measures are not available: the model contains
         (equality) constraints involving parameters from different
         levels."))
      return(list())
    }
    fit_ps <- try(
      lavaan(pt_ps,
        slot_options      = lavoptions,
        slot_sample_stats = object@SampleStats,
        slot_data         = lavdata,
        slot_cache        = object@Cache,
        sloth1            = object@h1
      ),
      silent = TRUE
    )
    if (inherits(fit_ps, "try-error") || !fit_ps@optim$converged) {
      lav_msg_warn(gettextf(
        "estimation of the partially saturated model for the %s level
         failed.", level_names[l]))
      next
    }

    # level-specific baseline model
    pt_bl <- lav_pt_partial_saturation(object,
      level.keep = l,
      independent = TRUE
    )
    fit_bl <- try(
      lavaan(pt_bl,
        slot_options      = lavoptions,
        slot_sample_stats = object@SampleStats,
        slot_data         = lavdata,
        slot_cache        = object@Cache,
        sloth1            = object@h1
      ),
      silent = TRUE
    )
    if (inherits(fit_bl, "try-error") || !fit_bl@optim$converged) {
      lav_msg_warn(gettextf(
        "estimation of the level-specific baseline model for the %s level
         failed.", level_names[l]))
    } else {
      fit_ps@baseline <- list(
        partable = fit_bl@ParTable,
        test = fit_bl@test
      )
    }

    out[[level_names[l]]] <- fit_ps
  }

  out
}

# resolve a user-specified 'level' argument (1/2, "within"/"between", or
# the level label used in the model syntax) to "within" or "between"
lav_fit_by_level_target <- function(object, level = NULL) {
  if (object@Data@nlevels < 2L) {
    lav_msg_stop(gettext(
      "level-specific fit measures are only available for multilevel
       models."))
  }

  level_names <- c("within", "between")
  target <- NULL
  if (is.numeric(level) && length(level) == 1L &&
      level %in% seq_len(object@Data@nlevels)) {
    target <- level_names[as.integer(level)]
  } else if (is.character(level) && length(level) == 1L) {
    idx <- pmatch(tolower(level), level_names)
    if (is.na(idx)) {
      # perhaps the label of a level, as used in the model syntax
      idx <- pmatch(tolower(level), tolower(object@Data@level.label))
    }
    if (!is.na(idx)) {
      target <- level_names[idx]
    }
  }
  if (is.null(target)) {
    lav_msg_stop(gettextf(
      "level must be a single level number (1 or 2) or level name
       (%s).", lav_msg_view(c(level_names, object@Data@level.label),
                            log_sep = "or")))
  }

  target
}

# retrieve the partially saturated fit for a given level, computing it
# on the fly if it was not stored (e.g., fit.by.level = FALSE)
lav_fit_by_level_get <- function(object, level = NULL) {
  target <- lav_fit_by_level_target(object, level = level)

  fbl <- object@internal$fit.by.level
  if (is.null(fbl) || length(fbl) == 0L) {
    # not stored: compute on the fly
    fbl <- lav_object_fit_by_level(object)
  }
  fit_l <- NULL
  if (length(fbl) > 0L) {
    fit_l <- fbl[[target]]
  }
  if (is.null(fit_l)) {
    lav_msg_stop(gettextf(
      "no partially saturated model is available for the %s level.",
      target))
  }

  fit_l
}

# build a small (measures x levels) matrix with level-specific fit
# measures, based on the STORED partially saturated fits (no on-the-fly
# computation); returns NULL if they are not available
# (used by summary(fit, fit.measures = TRUE))
lav_fit_by_level_fm <- function(object) {
  fbl <- object@internal$fit.by.level
  if (is.null(fbl) || length(fbl) == 0L) {
    return(NULL)
  }

  # scaled test statistics available?
  scaled <- any(!sapply(
    lapply(object@test, "[[", "scaling.factor"),
    is.null
  ))

  measures <- c("chisq", "df", "pvalue")
  if (scaled) {
    measures <- c(measures, "chisq.scaled", "df.scaled", "pvalue.scaled")
  }
  measures <- c(measures, "cfi", "tli", "rmsea")
  if (scaled) {
    measures <- c(measures, "cfi.robust", "tli.robust", "rmsea.robust")
  }

  level_names <- c("within", "between")
  srmr_names <- c("srmr_within", "srmr_between")
  res <- matrix(NA_real_,
    nrow = length(measures) + 1L,
    ncol = length(level_names)
  )
  rownames(res) <- c(measures, "srmr")
  colnames(res) <- level_names

  for (i in seq_along(level_names)) {
    fit_l <- fbl[[level_names[i]]]
    if (is.null(fit_l)) {
      next
    }
    fm <- try(
      lav_fit(
        object = fit_l,
        fit_measures = list(fit.measures = c(measures, srmr_names[i])),
        output = "vector"
      ),
      silent = TRUE
    )
    if (inherits(fm, "try-error")) {
      next
    }
    idx <- match(c(measures, srmr_names[i]), names(fm))
    res[, i] <- unname(fm[idx])
  }

  # drop all-NA robust rows (e.g., when a baseline model failed)
  class(res) <- c("lavaan.matrix", "matrix", "array")
  res
}

# print the (measures x levels) matrix produced by lav_fit_by_level_fm()
# the layout mirrors lav_fitmeasures_print(): the Within column lines up
# with the regular (single/standard) column, and the Between column lines
# up with the second column that is used when a scaled test statistic is
# requested
lav_fit_by_level_print <- function(x, nd = 3L) {
  if (is.null(x) || !is.matrix(x) || nrow(x) == 0L) {
    return(invisible(x))
  }

  labels <- c(
    chisq          = gettext("Test statistic"),
    df             = gettext("Degrees of freedom"),
    pvalue         = gettext("P-value"),
    chisq.scaled   = gettext("Test statistic (scaled)"),
    df.scaled      = gettext("Degrees of freedom (scaled)"),
    pvalue.scaled  = gettext("P-value (scaled)"),
    cfi            = gettext("Comparative Fit Index (CFI)"),
    tli            = gettext("Tucker-Lewis Index (TLI)"),
    cfi.robust     = gettext("Robust Comparative Fit Index (CFI)"),
    tli.robust     = gettext("Robust Tucker-Lewis Index (TLI)"),
    rmsea          = gettext("RMSEA"),
    rmsea.robust   = gettext("Robust RMSEA"),
    srmr           = gettext("SRMR")
  )

  num_format <- paste("%", max(8L, nd + 5L), ".", nd, "f", sep = "")

  cat("\n")
  cat(gettext("Level-specific Fit Measures (other level saturated):"),
      "\n\n", sep = "")

  # container three columns (same layout as lav_fitmeasures_print)
  c1 <- ""
  c2 <- gettext("Within")
  c3 <- gettext("Between")

  fmt1 <- function(val, m) {
    if (is.na(val)) {
      return("NA")
    }
    if (m %in% c("df", "df.scaled") && val %% 1 == 0) {
      return(as.character(as.integer(val)))
    }
    sprintf(num_format, val)
  }

  for (r in seq_len(nrow(x))) {
    m <- rownames(x)[r]
    label <- labels[m]
    if (is.na(label)) {
      label <- m
    }
    c1 <- c(c1, label)
    c2 <- c(c2, fmt1(x[r, 1L], m))
    c3 <- c(c3, fmt1(x[r, 2L], m))
  }

  # format c1/c2/c3 (identical widths as lav_fitmeasures_print)
  c1 <- format(c1, width = 35L)
  c2 <- format(c2, width = 16L + max(0, (nd - 3L)) * 4L, justify = "right")
  c3 <- format(c3, width = 8L + nd, justify = "right")

  # create character matrix
  m <- cbind(c1, c2, c3, deparse.level = 0)
  colnames(m) <- rep("", ncol(m))
  rownames(m) <- rep(" ", nrow(m))

  # print
  write.table(m, row.names = TRUE, col.names = FALSE, quote = FALSE)

  invisible(x)
}
