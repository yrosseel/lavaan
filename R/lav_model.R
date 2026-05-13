# constructor of the matrix lavoptions$representation
#
# initial version: YR 22/11/2010
# - YR 14 Jan 2014: moved to lav_model.R
# - YR 18 Nov 2014: more efficient handling of linear equality constraints
# - YR 02 Dec 2014: allow for bare-minimum parameter tables
# - YR 25 Jan 2017: collect options in lavoptions
# - YR 12 Mar 2021: add lavpta as argument; create model attributes (ma)
# - YR 21 Jan 2025: add composites

# construct MATRIX lavoptions$representation of the model
lav_model <- function(lavpartable = NULL,                          # nolint
                      lavoptions = NULL,
                      th_idx = list()) {
  # handle bare-minimum partables
  lavpartable <- lav_partable_complete(lavpartable)
  lavpta <- lav_partable_attributes(lavpartable)
  lavpartable <- lav_partable_set_cache(lavpartable, lavpta)

  # global info from user model
  nblocks <- lav_partable_nblocks(lavpartable)
  ngroups <- lav_partable_ngroups(lavpartable)
  meanstructure <- any(lavpartable$op == "~1")
  correlation <- lavoptions$correlation
  if (is.null(correlation)) {
    correlation <- FALSE
  }
  composites_option <- lavoptions$composites
  if (is.null(composites_option)) {
    composites_option <- TRUE
  }
  composites <- any(lavpartable$op == "<~") && composites_option
  categorical <- any(lavpartable$op == "|")
  if (categorical) {
    meanstructure <- TRUE

    # handle th.idx if length(th.idx) != nblocks
    if (nblocks != length(th_idx)) {
      th_idx <- rep(th_idx, each = nblocks)
    }
  }
  group_w_free <- any(lavpartable$lhs == "group" & lavpartable$op == "%")
  multilevel <- FALSE
  if (!is.null(lavpartable$level)) {
    nlevels <- lav_partable_nlevels(lavpartable)
    if (nlevels > 1L) {
      multilevel <- TRUE
    }
  } else {
    nlevels <- 1L
  }

  nefa <- lav_partable_nefa(lavpartable)
  if (nefa > 0L) {
    efa_values <- lav_partable_efa_values(lavpartable)
  }

  # check for simple equality constraints
  eq_simple <- any(lavpartable$free > 0L & duplicated(lavpartable$free))
  if (eq_simple) {
    # just like in <0.5-18, add (temporary) 'unco' column
    # so we can fill in x.unco.idx
    lavpartable$unco <- integer(length(lavpartable$id))
    idx_free <- which(lavpartable$free > 0L)
    lavpartable$unco[idx_free] <- seq_along(idx_free)
  }

  # handle variable definitions and (in)equality constraints
  tmp_con <- lav_constraints_parse(
    partable = lavpartable,
    constraints = NULL
  )

  # handle *linear* equality constraints special
  if (tmp_con$ceq.linear.only.flag) {
    con_jac <- tmp_con$ceq.JAC
    con_lambda <- numeric(nrow(tmp_con$ceq.JAC))
    attr(con_jac, "inactive.idx") <- integer(0L)
    attr(con_jac, "ceq.idx") <- seq_len(nrow(tmp_con$ceq.JAC))
  } else {
    con_jac <- matrix(0, 0, 0)
    con_lambda <- numeric(0)
  }

  # select model matrices
  if (lavoptions$representation == "LISREL") {
    tmp_rep <- lav_lisrel(lavpartable, target = NULL, extra = TRUE,
                          allow_composites = composites)
  } else if (lavoptions$representation == "RAM") {
    tmp_rep <- lav_ram(lavpartable, target = NULL, extra = TRUE)
  } else {
    lav_msg_stop(gettextf(
      "%1$s argument must be either %2$s or %3$s",
      "representation", "LISREL", "RAM"))
  }
  if (lav_debug()) print(tmp_rep)

  # FIXME: check for non-existing parameters
  bad_idx <-
    which((tmp_rep$mat == "" | is.na(tmp_rep$row) | is.na(tmp_rep$col)) &
          !lavpartable$op %in% c("==", "<", ">", ":=", "|~"))

  if (length(bad_idx) > 0L) {
    this_formula <- paste(lavpartable$lhs[bad_idx[1]],
      lavpartable$op[bad_idx[1]],
      lavpartable$rhs[bad_idx[1]],
      sep = " "
    )
    if (lavoptions$representation == "LISREL") {
      lav_msg_stop(gettextf(
        "a model parameter is not defined in the LISREL representation %s.
        Upgrade to latent variables or consider using representation = 'RAM'.",
        this_formula))
    } else {
      lav_msg_stop(
        gettextf("parameter is not defined: %s", this_formula)
      )
    }
  }

  # prepare nG-sized slots
  tmp_ng <- sum(unlist(attr(tmp_rep, "mmNumber")))
  tmp_glist <- vector(mode = "list", tmp_ng)
  names(tmp_glist) <- unlist(attr(tmp_rep, "mmNames"))
  dim_names <- vector(mode = "list", length = tmp_ng)
  is_symmetric <- logical(tmp_ng)
  mm_size <- integer(tmp_ng)

  m_free_idx <- m_user_idx <- vector(mode = "list", length = tmp_ng)
  x_free_idx <- x_unco_idx <- x_user_idx <- vector(
    mode = "list",
    length = tmp_ng
  )

  # prepare nblocks-sized slots
  nvar <- integer(nblocks)
  nmat <- unlist(attr(tmp_rep, "mmNumber"))
  mm_idx <- lav_model_group_mm_indices(nmat)
  num_idx <- vector("list", length = nblocks)
  nexo <- integer(nblocks)
  ov_x_dummy_ov_idx <- vector(mode = "list", length = nblocks)
  ov_x_dummy_lv_idx <- vector(mode = "list", length = nblocks)
  ov_y_dummy_ov_idx <- vector(mode = "list", length = nblocks)
  ov_y_dummy_lv_idx <- vector(mode = "list", length = nblocks)
  ov_efa_idx <- vector(mode = "list", length = nblocks)
  lv_efa_idx <- vector(mode = "list", length = nblocks)

  offset <- 0L
  # keep track of ov.names across blocks
  for (g in 1:nblocks) {
    # observed and latent variables for this block
    ov_names <- lav_partable_vnames(lavpartable, "ov", block = g)
    ov_names_nox <- lav_partable_vnames(lavpartable, "ov.nox", block = g)
    ov_names_x <- lav_partable_vnames(lavpartable, "ov.x", block = g)
    ov_num <- lav_partable_vnames(lavpartable, "ov.num", block = g)
    if (lavoptions$conditional.x) {
      if (nlevels > 1L) {
        if (ngroups == 1L) {
          other_block_names <- lav_partable_vnames(lavpartable, "ov",
            block = seq_len(nblocks)[-g]
          )
        } else {
          # TEST ME!
          # which group is this?
          this_group <- ceiling(g / nlevels)
          blocks_within_group <- (this_group - 1L) * nlevels + seq_len(nlevels)
          other_block_names <- lav_partable_vnames(lavpartable, "ov",
            block = blocks_within_group[-g]
          )
        }


        if (length(ov_names_x) > 0L) {
          idx <- which(ov_names_x %in% other_block_names)
          if (length(idx) > 0L) {
            ov_names_nox <- unique(c(ov_names_nox, ov_names_x[idx]))
            ov_names_x <- ov_names_x[-idx]
            ov_names <- ov_names_nox
          }
        }
      }
      nvar[g] <- length(ov_names_nox)
      if (correlation) {
        num_idx[[g]] <- integer(0L)
      } else {
        num_idx[[g]] <- which(ov_names_nox %in% ov_num)
      }
    } else {
      nvar[g] <- length(ov_names)
      if (correlation) {
        num_idx[[g]] <- integer(0L)
      } else {
        num_idx[[g]] <- which(ov_names %in% ov_num)
      }
    }
    nexo[g] <- length(ov_names_x)

    if (nefa > 0L) {
      lv_names <- lav_partable_vnames(lavpartable, "lv", block = g)
    }

    # model matrices for this block
    mm_number <- attr(tmp_rep, "mmNumber")[[g]]
    mm_names <- attr(tmp_rep, "mmNames")[[g]]
    mm_symmetric <- attr(tmp_rep, "mmSymmetric")[[g]]
    mm_dim_names <- attr(tmp_rep, "mmDimNames")[[g]]
    mm_rows <- attr(tmp_rep, "mmRows")[[g]]
    mm_cols <- attr(tmp_rep, "mmCols")[[g]]

    for (mm in 1:mm_number) {
      # offset in tmp.glist
      offset <- offset + 1L

      # matrix size, symmetric, dim.names
      if (mm_symmetric[mm]) {
        tmp_n <- mm_rows[mm]
        mm_size <- as.integer(tmp_n * (tmp_n + 1) / 2)
      } else {
        mm_size <- as.integer(mm_rows[mm] * mm_cols[mm])
      }
      mm_size[offset] <- mm_size
      is_symmetric[offset] <- mm_symmetric[mm]
      dim_names[[offset]] <- mm_dim_names[[mm]]

      # select elements for this matrix
      idx <- which(lavpartable$block == g & tmp_rep$mat == mm_names[mm])

      # 1. first assign free values only, to get vector index
      #    -> to be used in lav_model_objective
      free_idx <- lav_matrix_rowcol_idx(
        row = tmp_rep$row[idx],
        col = tmp_rep$col[idx],
        value = lavpartable$free[idx],
        nrow = mm_rows[mm],
        ncol = mm_cols[mm],
        symmetric = mm_symmetric[mm]
      )
      m_free_idx[[offset]] <- free_idx$m.idx
      x_free_idx[[offset]] <- free_idx$x.idx

      # 2. if simple equality constraints, unconstrained free parameters
      #    -> to be used in lav_model_gradient
      if (eq_simple) {
        unco_idx <- lav_matrix_rowcol_idx(
          row = tmp_rep$row[idx],
          col = tmp_rep$col[idx],
          value = lavpartable$unco[idx],
          nrow = mm_rows[mm],
          ncol = mm_cols[mm],
          symmetric = mm_symmetric[mm]
        )
        # m.unco.idx[[offset]] <-     which(tmp > 0)
        x_unco_idx[[offset]] <- unco_idx$x.idx
      } else {
        # m.unco.idx[[offset]] <- m.free.idx[[offset]]
        x_unco_idx[[offset]] <- x_free_idx[[offset]]
      }

      # 3. general mapping between user and tmp.glist
      user_idx <- lav_matrix_rowcol_idx(
        row = tmp_rep$row[idx],
        col = tmp_rep$col[idx],
        value = lavpartable$id[idx],
        nrow = mm_rows[mm],
        ncol = mm_cols[mm],
        symmetric = mm_symmetric[mm]
      )
      m_user_idx[[offset]] <- user_idx$m.idx
      x_user_idx[[offset]] <- user_idx$x.idx

      # 4. now assign starting/fixed values
      # create empty matrix
      # FIXME: again, we may want to use sparse matrices here...
      tmp <- matrix(0.0,
        nrow = mm_rows[mm],
        ncol = mm_cols[mm]
      )
      tmp[cbind(tmp_rep$row[idx], tmp_rep$col[idx])] <- lavpartable$start[idx]
      if (mm_symmetric[mm]) {
        tmp_tt <- t(tmp)
        tmp[lower.tri(tmp)] <- tmp_tt[lower.tri(tmp_tt)]
      }

      # 4b. override with cov.x (if conditional.x = TRUE)
      # new in 0.6-1
      # shouldn't be needed, if lavpartable$start contains cov.x values
      # if(mm.names[mm] == "cov.x") {
      #    tmp <- cov.x[[g]]
      # }
      # 4c. override with mean.x (if conditional.x = TRUE)
      # new in 0.6-1
      # shouldn't be needed, if lavpartable$start contains mean.x values
      # if(mm.names[mm] == "mean.x") {
      #    tmp <- as.matrix(mean.x[[g]])
      # }

      # representation specific stuff
      if (lavoptions$representation == "LISREL" &&
        mm_names[mm] == "lambda") {
        ov_dummy_names_nox <- attr(tmp_rep, "ov.dummy.names.nox")[[g]]
        ov_dummy_names_x <- attr(tmp_rep, "ov.dummy.names.x")[[g]]
        ov_dummy_names <- c(ov_dummy_names_nox, ov_dummy_names_x)
        # define dummy latent variables
        if (length(ov_dummy_names)) {
          # in this case, lv.names will be extended with the dummys
          tmp_lv_names <- mm_dim_names$psi[[1]]
          row_tmp_idx <- match(ov_dummy_names, ov_names)
          col_tmp_idx <- match(ov_dummy_names, tmp_lv_names)
          # Fix lambda values to 1.0
          tmp[cbind(row_tmp_idx, col_tmp_idx)] <- 1.0

          ov_x_dummy_ov_idx[[g]] <- match(ov_dummy_names_x, ov_names)
          ov_x_dummy_lv_idx[[g]] <- match(ov_dummy_names_x, tmp_lv_names)
          ov_y_dummy_ov_idx[[g]] <- match(ov_dummy_names_nox, ov_names)
          ov_y_dummy_lv_idx[[g]] <- match(ov_dummy_names_nox, tmp_lv_names)
        }
      }

      # representation specific
      if (lavoptions$representation == "LISREL" && mm_names[mm] == "delta") {
        # only categorical values are listed in the lavpartable
        # but all remaining values should be 1.0
        idx <- which(tmp[, 1L] == 0.0)
        tmp[idx, 1L] <- 1.0
      }

      # representation specific
      if (lavoptions$representation == "RAM" && mm_names[mm] == "ov.idx") {
        tmp[1, ] <- attr(tmp_rep, "ov.idx")[[g]]
      }

      # assign matrix to tmp.glist
      tmp_glist[[offset]] <- tmp
    } # mm

    # efa related info
    if (nefa > 0L) {
      ov_efa_idx[[g]] <- vector("list", length = nefa)
      lv_efa_idx[[g]] <- vector("list", length = nefa)
      for (set in seq_len(nefa)) {
        # determine ov idx for this set
        ov_efa <-
          unique(lavpartable$rhs[lavpartable$op == "=~" &
            lavpartable$block == g &
            lavpartable$efa == efa_values[set]])
        ov_efa_idx[[g]][[set]] <- match(ov_efa, ov_names)

        lv_efa <-
          unique(lavpartable$lhs[lavpartable$op == "=~" &
            lavpartable$block == g &
            lavpartable$efa == efa_values[set]])
        lv_efa_idx[[g]][[set]] <- match(lv_efa, lv_names)
      }
      names(ov_efa_idx[[g]]) <- efa_values
      names(lv_efa_idx[[g]]) <- efa_values
    } # efa

    # set variances composites (new in 0.6-20)
    if (composites) {
      mm_in_group <- mm_idx[[g]]
      tmp_glist[mm_in_group] <-
        lav_lisrel_comp_set_intresvar(tmp_glist[mm_in_group])
    }
  } # g

  # fixed.x parameters?
  # fixed.x <- any(lavpartable$exo > 0L & lavpartable$free == 0L)
  # if(categorical) {
  #    fixed.x <- TRUE
  # }

  # dirty hack to mimic MUML
  if (!is.null(lavoptions$tech.muml.scale)) {
    lav_msg_warn(gettext("using muml scale in group 2"))

    # find matrix
    lambda_idx <- which(names(tmp_glist) == "lambda")[2L]

    # find rows/cols
    b_names <- paste0("b", ov_names) ## ad-hoc assumption!!!
    tmp_cols <- match(b_names, tmp_lv_names)
    tmp_rows <- seq_len(nvar[2])
    stopifnot(length(tmp_cols) == length(tmp_rows))
    tmp_glist[[lambda_idx]][cbind(tmp_rows, tmp_cols)] <-
      lavoptions$tech.muml.scale
  }

  # which free parameters are observed variances?
  ov_names <- lav_partable_vnames(lavpartable, "ov")
  x_free_var_idx <- lavpartable$free[lavpartable$free &
    # !duplicated(lavpartable$free) &
    lavpartable$lhs %in% ov_names &
    lavpartable$op == "~~" &
    lavpartable$lhs == lavpartable$rhs]

  rv_lv <- rv_ov <- list()
  if (multilevel) {
    # store information about random slopes (if any)
    lv_names <- lav_partable_vnames(lavpartable, "lv")
    # we should also add split-y names (x) to lv.names
    # FIXME: make this work for multiple work multilevel
    level_values <- lav_partable_level_values(lavpartable)
    ovx1 <- lav_object_vnames(lavpartable, "ov.x", level = level_values[1])
    ovx2 <- lav_object_vnames(lavpartable, "ov.x", level = level_values[2])
    ovx12 <- ovx2[ovx2 %in% ovx1]
    lv_names <- c(lv_names, ovx12)

    # RV LV
    rv_idx <- which(nchar(lavpartable$rv) > 0L &
      lavpartable$level == level_values[1] &
      lavpartable$rhs %in% lv_names)
    if (length(rv_idx)) {
      rv_lv <- lapply(rv_idx, function(x) {
        c(lavpartable$lhs[x], lavpartable$rhs[x])
      })
      names(rv_lv) <- lavpartable$rv[rv_idx]
    }

    # RV OV
    rv_idx <- which(nchar(lavpartable$rv) > 0L &
      lavpartable$level == level_values[1] &
      !lavpartable$rhs %in% lv_names)
    if (length(rv_idx)) {
      rv_ov <- lapply(rv_idx, function(x) {
        c(lavpartable$lhs[x], lavpartable$rhs[x])
      })
      names(rv_ov) <- lavpartable$rv[rv_idx]
    }
  } # multilevel

  # new in 0.6-9: model properties
  modprop <- lav_model_properties(
    glist = tmp_glist,
    lavpartable = lavpartable,
    nmat = nmat,
    m_free_idx = m_free_idx
  )

  tmp_model <- new("lavModel",
    GLIST = tmp_glist,
    dimNames = dim_names,
    isSymmetric = is_symmetric,
    mmSize = mm_size,
    representation = lavoptions$representation,
    modprop = modprop,
    meanstructure = meanstructure,
    correlation = correlation,
    composites = composites,
    categorical = categorical,
    multilevel = multilevel,
    link = lavoptions$link,
    nblocks = nblocks,
    ngroups = ngroups, # breaks rsem????
    nefa = nefa,
    group.w.free = group_w_free,
    nmat = nmat,
    mm.idx = mm_idx,
    nvar = nvar,
    num.idx = num_idx,
    th.idx = th_idx,
    nx.free = max(lavpartable$free),
    nx.unco = if (is.null(lavpartable$unco)) {
      max(lavpartable$free)
    } else {
      max(lavpartable$unco)
    },
    nx.user = max(lavpartable$id),
    m.free.idx = m_free_idx,
    x.free.idx = x_free_idx,
    x.free.var.idx = x_free_var_idx,
    # m.unco.idx=m.unco.idx,
    x.unco.idx = x_unco_idx,
    m.user.idx = m_user_idx,
    x.user.idx = x_user_idx,
    x.def.idx = which(lavpartable$op == ":="),
    x.ceq.idx = which(lavpartable$op == "=="),
    x.cin.idx = which(lavpartable$op == ">" | lavpartable$op == "<"),
    ceq.simple.only = tmp_con$ceq.simple.only,
    ceq.simple.K = tmp_con$ceq.simple.K,
    eq.constraints = tmp_con$ceq.linear.only.flag,
    eq.constraints.K = tmp_con$ceq.JAC.NULL,
    eq.constraints.k0 = tmp_con$ceq.rhs.NULL,
    def.function = tmp_con$def.function,
    ceq.function = tmp_con$ceq.function,
    ceq.JAC = tmp_con$ceq.JAC,
    ceq.rhs = tmp_con$ceq.rhs,
    ceq.jacobian = tmp_con$ceq.jacobian,
    ceq.linear.idx = tmp_con$ceq.linear.idx,
    ceq.nonlinear.idx = tmp_con$ceq.nonlinear.idx,
    cin.simple.only = tmp_con$cin.simple.only,
    cin.function = tmp_con$cin.function,
    cin.JAC = tmp_con$cin.JAC,
    cin.rhs = tmp_con$cin.rhs,
    cin.jacobian = tmp_con$cin.jacobian,
    cin.linear.idx = tmp_con$cin.linear.idx,
    cin.nonlinear.idx = tmp_con$cin.nonlinear.idx,
    con.jac = con_jac,
    con.lambda = con_lambda,
    nexo = nexo,
    fixed.x = lavoptions$fixed.x,
    conditional.x = lavoptions$conditional.x,
    parameterization = lavoptions$parameterization,
    ov.x.dummy.ov.idx = ov_x_dummy_ov_idx,
    ov.x.dummy.lv.idx = ov_x_dummy_lv_idx,
    ov.y.dummy.ov.idx = ov_y_dummy_ov_idx,
    ov.y.dummy.lv.idx = ov_y_dummy_lv_idx,
    ov.efa.idx = ov_efa_idx,
    lv.efa.idx = lv_efa_idx,
    rv.lv = rv_lv,
    rv.ov = rv_ov,
    estimator = lavoptions$estimator,
    estimator.args = lavoptions$estimator.args
  )

  if (lav_debug()) {
    cat("lavaan debug: lavaanModel\n")
    print(str(tmp_model))
    print(tmp_model@GLIST)
  }

  tmp_model
}

# for backwards compatibility
# tmp.model <- lav_model
