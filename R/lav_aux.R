# auxiliary (aux=) variable handling
#
# auxiliary variables are observed variables that are NOT part of the
# user-specified model. They are used to make the missing-at-random (MAR)
# assumption more plausible under missing data. Two mechanisms are supported:
#
#  - missing = "two.stage" / "robust.two.stage": the auxiliary variables are
#    included in the first-stage EM run that estimates the saturated moments
#    the model is fitted to (see lav_samplestats.R); the model estimates are
#    therefore affected. (Handled via the lavData@aux channel + the EM code.)
#
#  - missing = "ml" / "fiml" (full-information ML): the auxiliary variables
#    are added to the model as a block of "saturated correlates" (Graham,
#    2003): each auxiliary variable is freely correlated with every other
#    auxiliary variable, with every model observed variable, and has a free
#    mean. The (enlarged) model is then estimated with casewise FIML, so the
#    auxiliary variables affect both the point estimates and the standard
#    errors, while leaving the model degrees of freedom unchanged. The
#    auxiliary parameters are nuisance parameters and are hidden from the
#    default output. (Handled by the functions below.)

# clean up the user-supplied aux= names: drop names that are not found in the
# data, that are model (observed) variables, or that are binary or
# ordered/categorical (not supported). Returns a character vector of the
# (continuous) auxiliary variable names that can be used.
lav_aux_clean_names <- function(aux      = NULL,
                                data     = NULL,
                                ov_names = NULL,
                                ordered  = NULL) {
  if (is.null(aux)) {
    return(character(0L))
  }
  if (!is.character(aux)) {
    lav_msg_stop(gettext("aux= argument must be a character vector."))
  }
  aux <- unique(aux[nchar(aux) > 0L])
  if (length(aux) == 0L) {
    return(character(0L))
  }

  # we need raw data
  if (is.null(data) || !(is.data.frame(data) || is.matrix(data))) {
    return(character(0L))
  }
  data_names <- if (is.data.frame(data)) names(data) else colnames(data)

  # 1) drop aux names that are not found in the data
  bad_idx <- which(!aux %in% data_names)
  if (length(bad_idx) > 0L) {
    lav_msg_warn(gettextf(
      "auxiliary (aux) variable(s) not found in the data and ignored: %s",
      paste(aux[bad_idx], collapse = " ")))
    aux <- aux[-bad_idx]
  }
  if (length(aux) == 0L) {
    return(character(0L))
  }

  # 2) drop aux names that are model (observed) variables
  model_ov <- unique(unlist(ov_names))
  bad_idx <- which(aux %in% model_ov)
  if (length(bad_idx) > 0L) {
    lav_msg_note(gettextf(
      "auxiliary (aux) variable(s) are also model variables and ignored
       as auxiliary: %s", paste(aux[bad_idx], collapse = " ")))
    aux <- aux[-bad_idx]
  }
  if (length(aux) == 0L) {
    return(character(0L))
  }

  # 3) drop binary or ordered/categorical aux variables
  cat_flag <- logical(length(aux))
  for (i in seq_along(aux)) {
    if (!is.null(ordered) && aux[i] %in% ordered) {
      cat_flag[i] <- TRUE
      next
    }
    vv <- if (is.data.frame(data)) data[[aux[i]]] else data[, aux[i]]
    if (is.factor(vv) || is.ordered(vv) || is.character(vv) ||
        is.logical(vv)) {
      cat_flag[i] <- TRUE
    } else {
      # numeric: binary?
      if (length(unique(vv[!is.na(vv)])) <= 2L) {
        cat_flag[i] <- TRUE
      }
    }
  }
  if (any(cat_flag)) {
    lav_msg_warn(gettextf(
      "auxiliary (aux) variable(s) are binary or ordered/categorical;
       this is not supported and they are removed: %s",
      paste(aux[cat_flag], collapse = " ")))
    aux <- aux[!cat_flag]
  }

  aux
}

# append the 'saturated correlates' block for the auxiliary variables to a
# flat model (the parsed-syntax object produced by lavParseModelString). The
# new entries are all free parameters (mod.idx = 0, no modifier):
#   aux_i ~~ aux_j   (auxiliary variances and covariances)
#   aux_i ~ 1        (auxiliary means)
#   aux_i ~~ ov_k    (auxiliary correlated with every model observed variable)
# lavParTable() expands these per group/block, just like ordinary syntax.
lav_aux_satcor_flat <- function(flat_model = NULL,
                                aux        = character(0L),
                                ov_names   = character(0L)) {
  if (length(aux) == 0L) {
    return(flat_model)
  }
  ov <- setdiff(unique(unlist(ov_names)), aux)
  na <- length(aux)

  lhs <- character(0L)
  op  <- character(0L)
  rhs <- character(0L)

  # aux ~~ aux (variances + covariances)
  for (i in seq_len(na)) {
    for (j in i:na) {
      lhs <- c(lhs, aux[i]); op <- c(op, "~~"); rhs <- c(rhs, aux[j])
    }
  }
  # aux ~ 1 (means)
  lhs <- c(lhs, aux); op <- c(op, rep("~1", na)); rhs <- c(rhs, rep("", na))
  # aux ~~ model observed variables
  for (a in aux) {
    for (v in ov) {
      lhs <- c(lhs, a); op <- c(op, "~~"); rhs <- c(rhs, v)
    }
  }
  nn <- length(lhs)

  # append to each column of flat_model (handle the known columns explicitly,
  # pad any other columns by type)
  for (cn in names(flat_model)) {
    val <- switch(cn,
      lhs     = lhs,
      op      = op,
      rhs     = rhs,
      block   = rep(1L, nn),
      mod.idx = rep(0L, nn),
      NULL
    )
    if (is.null(val)) {
      old <- flat_model[[cn]]
      if (is.character(old)) {
        val <- rep("", nn)
      } else if (is.integer(old)) {
        val <- rep(0L, nn)
      } else {
        val <- rep(as.numeric(NA), nn)
      }
    }
    flat_model[[cn]] <- c(flat_model[[cn]], val)
  }
  # modifiers/constraints attributes are unchanged (all new rows are free)
  attr(flat_model, "aux") <- aux

  flat_model
}

# augment the baseline (independence) parameter table when auxiliary variables
# are present: the model observed variables stay mutually uncorrelated (as in
# the independence model), but the auxiliary variables become a saturated block
# (freely correlated with each other and with every model observed variable).
# We start from the ordinary independence partable (which already contains the
# variances and means in the correct variable order -- crucial, so that the
# saturated reference stays aligned with the sample statistics) and merge in the
# free auxiliary covariances. This matches the baseline used by semTools'
# auxiliary() so that CFI/TLI/NFI are computed correctly.
lav_aux_satcor_baseline <- function(indep_partable = NULL,
                                    aux            = character(0L),
                                    ov_names       = character(0L),
                                    ngroups        = 1L) {
  if (length(aux) == 0L) {
    return(indep_partable)
  }
  ov <- setdiff(unique(unlist(ov_names)), aux)
  na <- length(aux)
  blk <- character(0L)
  # auxiliary covariances (off-diagonal; variances already in independence pt)
  if (na > 1L) {
    for (i in 1:(na - 1L)) {
      for (j in (i + 1L):na) {
        blk <- c(blk, paste(aux[i], "~~", aux[j]))
      }
    }
  }
  # auxiliary ~~ model observed variables
  for (a in aux) {
    for (v in ov) {
      blk <- c(blk, paste(a, "~~", v))
    }
  }
  if (length(blk) == 0L) {
    return(indep_partable)
  }

  sat_pt <- lavaanify(paste(blk, collapse = "\n"), ngroups = ngroups)
  sat_pt <- sat_pt[, c("lhs", "op", "rhs", "user", "block", "group")]

  pt_df <- as.data.frame(indep_partable, stringsAsFactors = FALSE)
  merged <- lav_partable_merge(pt_df, sat_pt,
    remove.duplicated = TRUE, warn = FALSE
  )
  n_new <- nrow(merged) - nrow(pt_df)
  if (n_new > 0L) {
    new_rows <- seq_len(n_new) + nrow(pt_df)
    merged$free[new_rows]   <- seq_len(n_new) + max(pt_df$free)
    merged$plabel[new_rows] <- paste0(".p", new_rows, ".")
  }
  as.list(merged)
}

# two-stage estimation with auxiliary variables (Savalei & Bentler, 2009):
# the stage-1 saturated moments are estimated over the augmented set
# [model variables, auxiliary variables]; the model is then fitted to the
# model-variable subset. For the (robust) standard errors, the stage-1
# asymptotic covariance Omega must therefore be computed over the augmented
# set and then reduced to the model-variable moments. This function returns
# the indices, within the augmented moment vector [mu (p_aug),
# vech(Sigma) (p_aug*(p_aug+1)/2)], that correspond to the model-variable
# moments [mu_model, vech(Sigma_model)] (the model variables come first, the
# auxiliary variables last). Selecting this sub-block of Omega gives the
# correct covariance of the model-variable stage-1 estimates (the sub-matrix
# of the inverse information, NOT the inverse of the model-variable
# information).
lav_aux_moment_idx <- function(p_model = 0L, p_aug = 0L) {
  mean_idx <- seq_len(p_model)
  ri <- lav_mat_vech_row_idx(p_aug) # incl. diagonal
  ci <- lav_mat_vech_col_idx(p_aug)
  vech_keep <- which(ri <= p_model & ci <= p_model)
  c(mean_idx, p_aug + vech_keep)
}

# given a (fitted) lavaan parameter table and the auxiliary variable names,
# return a logical vector marking the rows that belong to the auxiliary
# (saturated-correlates) block (i.e. any row that involves an auxiliary
# variable). Used to hide these nuisance parameters from the default output.
lav_aux_idx <- function(partable = NULL, aux = character(0L)) {
  if (length(aux) == 0L) {
    return(rep(FALSE, length(partable$lhs)))
  }
  (partable$lhs %in% aux) | (partable$rhs %in% aux)
}
