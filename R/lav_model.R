# constructor of the matrix lavoptions$representation
#
# initial version: YR 22/11/2010
# - YR 14 Jan 2014: moved to lav_model.R
# - YR 18 Nov 2014: more efficient handling of linear equality constraints
# - YR 02 Dec 2014: allow for bare-minimum parameter tables
# - YR 25 Jan 2017: collect options in lavoptions
# - YR 12 Mar 2021: add lavpta as argument; create model attributes (ma)

# construct MATRIX lavoptions$representation of the model
lav_model <- function(lavpartable = NULL,                          # nolint
                      lavoptions = NULL,
                      th.idx = list()) {
  # handle bare-minimum partables
  lavpartable <- lav_partable_complete(lavpartable)
  lavpta = lav_partable_attributes(lavpartable)
  lavpartable <- lav_partable_set_cache(lavpartable, lavpta)

  # global info from user model
  nblocks <- lav_partable_nblocks(lavpartable)
  ngroups <- lav_partable_ngroups(lavpartable)
  meanstructure <- any(lavpartable$op == "~1")
  correlation <- lavoptions$correlation
  if (is.null(correlation)) {
    correlation <- FALSE
  }
  categorical <- any(lavpartable$op == "|")
  if (categorical) {
    meanstructure <- TRUE

    # handle th.idx if length(th.idx) != nblocks
    if (nblocks != length(th.idx)) {
      th.idx <- rep(th.idx, each = nblocks)
    }
  }
  group.w.free <- any(lavpartable$lhs == "group" & lavpartable$op == "%")
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
    efa.values <- lav_partable_efa_values(lavpartable)
  }

  # check for simple equality constraints
  eq.simple <- any(lavpartable$free > 0L & duplicated(lavpartable$free))
  if (eq.simple) {
    # just like in <0.5-18, add (temporary) 'unco' column
    # so we can fill in x.unco.idx
    lavpartable$unco <- integer(length(lavpartable$id))
    idx.free <- which(lavpartable$free > 0L)
    lavpartable$unco[idx.free] <- seq_along(idx.free)
  }

  # handle variable definitions and (in)equality constraints
  tmp.con <- lav_constraints_parse(
    partable = lavpartable,
    constraints = NULL,
    debug = lavoptions$debug
  )

  # handle *linear* equality constraints special
  if (tmp.con$ceq.linear.only.flag) {
    con.jac <- tmp.con$ceq.JAC
    con.lambda <- numeric(nrow(tmp.con$ceq.JAC))
    attr(con.jac, "inactive.idx") <- integer(0L)
    attr(con.jac, "ceq.idx") <- seq_len(nrow(tmp.con$ceq.JAC))
  } else {
    con.jac <- matrix(0, 0, 0)
    con.lambda <- numeric(0)
  }

  # select model matrices
  if (lavoptions$representation == "LISREL") {
    tmp.rep <- lav_lisrel(lavpartable, target = NULL, extra = TRUE)
  } else if (lavoptions$representation == "RAM") {
    tmp.rep <- lav_ram(lavpartable, target = NULL, extra = TRUE)
  } else {
    lav_msg_stop(gettextf(
      "%1$s argument must be either %2$s or %3$s",
      "representation", "LISREL", "RAM"))
  }
  if (lavoptions$debug) print(tmp.rep)

  # FIXME: check for non-existing parameters
  bad.idx <- which((tmp.rep$mat == "" | is.na(tmp.rep$row) | is.na(tmp.rep$col)) &
    !lavpartable$op %in% c("==", "<", ">", ":="))

  if (length(bad.idx) > 0L) {
    this.formula <- paste(lavpartable$lhs[bad.idx[1]],
      lavpartable$op[bad.idx[1]],
      lavpartable$rhs[bad.idx[1]],
      sep = " "
    )
    if (lavoptions$representation == "LISREL") {
      lav_msg_stop(gettextf(
        "a model parameter is not defined in the LISREL representation %s.
        Upgrade to latent variables or consider using representation = 'RAM'.",
        this.formula)      )
    } else {
      lav_msg_stop(
        gettextf("parameter is not defined: %s", this.formula)
      )
    }
  }

  # prepare nG-sized slots
  tmp.ng <- sum(unlist(attr(tmp.rep, "mmNumber")))
  tmp.glist <- vector(mode = "list", tmp.ng)
  names(tmp.glist) <- unlist(attr(tmp.rep, "mmNames"))
  dim.names <- vector(mode = "list", length = tmp.ng)
  is.symmetric <- logical(tmp.ng)
  mm.size <- integer(tmp.ng)

  m.free.idx <- m.user.idx <- vector(mode = "list", length = tmp.ng)
  x.free.idx <- x.unco.idx <- x.user.idx <- vector(
    mode = "list",
    length = tmp.ng
  )

  # prepare nblocks-sized slots
  nvar <- integer(nblocks)
  nmat <- unlist(attr(tmp.rep, "mmNumber"))
  num.idx <- vector("list", length = nblocks)
  nexo <- integer(nblocks)
  ov.x.dummy.ov.idx <- vector(mode = "list", length = nblocks)
  ov.x.dummy.lv.idx <- vector(mode = "list", length = nblocks)
  ov.y.dummy.ov.idx <- vector(mode = "list", length = nblocks)
  ov.y.dummy.lv.idx <- vector(mode = "list", length = nblocks)
  ov.efa.idx <- vector(mode = "list", length = nblocks)
  lv.efa.idx <- vector(mode = "list", length = nblocks)

  offset <- 0L
  # keep track of ov.names across blocks
  for (g in 1:nblocks) {
    # observed and latent variables for this block
    ov.names <- lav_partable_vnames(lavpartable, "ov", block = g)
    ov.names.nox <- lav_partable_vnames(lavpartable, "ov.nox", block = g)
    ov.names.x <- lav_partable_vnames(lavpartable, "ov.x", block = g)
    ov.num <- lav_partable_vnames(lavpartable, "ov.num", block = g)
    if (lavoptions$conditional.x) {
      if (nlevels > 1L) {
        if (ngroups == 1L) {
          other.block.names <- lav_partable_vnames(lavpartable, "ov",
            block = seq_len(nblocks)[-g]
          )
        } else {
          # TEST ME!
          # which group is this?
          this.group <- ceiling(g / nlevels)
          blocks.within.group <- (this.group - 1L) * nlevels + seq_len(nlevels)
          other.block.names <- lav_partable_vnames(lavpartable, "ov",
            block = blocks.within.group[-g]
          )
        }


        if (length(ov.names.x) > 0L) {
          idx <- which(ov.names.x %in% other.block.names)
          if (length(idx) > 0L) {
            ov.names.nox <- unique(c(ov.names.nox, ov.names.x[idx]))
            ov.names.x <- ov.names.x[-idx]
            ov.names <- ov.names.nox
          }
        }
      }
      nvar[g] <- length(ov.names.nox)
      if (correlation) {
        num.idx[[g]] <- integer(0L)
      } else {
        num.idx[[g]] <- which(ov.names.nox %in% ov.num)
      }
    } else {
      nvar[g] <- length(ov.names)
      if (correlation) {
        num.idx[[g]] <- integer(0L)
      } else {
        num.idx[[g]] <- which(ov.names %in% ov.num)
      }
    }
    nexo[g] <- length(ov.names.x)

    if (nefa > 0L) {
      lv.names <- lav_partable_vnames(lavpartable, "lv", block = g)
    }

    # model matrices for this block
    mm.number <- attr(tmp.rep, "mmNumber")[[g]]
    mm.names <- attr(tmp.rep, "mmNames")[[g]]
    mm.symmetric <- attr(tmp.rep, "mmSymmetric")[[g]]
    mm.dim.names <- attr(tmp.rep, "mmDimNames")[[g]]
    mm.rows <- attr(tmp.rep, "mmRows")[[g]]
    mm.cols <- attr(tmp.rep, "mmCols")[[g]]

    for (mm in 1:mm.number) {
      # offset in tmp.glist
      offset <- offset + 1L

      # matrix size, symmetric, dim.names
      if (mm.symmetric[mm]) {
        tmp.n <- mm.rows[mm]
        mm.size <- as.integer(tmp.n * (tmp.n + 1) / 2)
      } else {
        mm.size <- as.integer(mm.rows[mm] * mm.cols[mm])
      }
      mm.size[offset] <- mm.size
      is.symmetric[offset] <- mm.symmetric[mm]
      dim.names[[offset]] <- mm.dim.names[[mm]]

      # select elements for this matrix
      idx <- which(lavpartable$block == g & tmp.rep$mat == mm.names[mm])

      # create empty `pattern' matrix
      # FIXME: one day, we may want to use sparse matrices...
      #        but they should not slow things down!
      tmp <- matrix(0L,
        nrow = mm.rows[mm],
        ncol = mm.cols[mm]
      )

      # 1. first assign free values only, to get vector index
      #    -> to be used in lav_model_objective
      tmp[cbind(tmp.rep$row[idx], tmp.rep$col[idx])] <- lavpartable$free[idx]
      if (mm.symmetric[mm]) {
        # NOTE: we assume everything is in the UPPER tri!
        tmp.tt <- t(tmp)
        tmp[lower.tri(tmp)] <- tmp.tt[lower.tri(tmp.tt)]
      }
      m.free.idx[[offset]] <- which(tmp > 0)
      x.free.idx[[offset]] <- tmp[which(tmp > 0)]

      # 2. if simple equality constraints, unconstrained free parameters
      #    -> to be used in lav_model_gradient
      if (eq.simple) {
        tmp[cbind(
          tmp.rep$row[idx],
          tmp.rep$col[idx]
        )] <- lavpartable$unco[idx]
        if (mm.symmetric[mm]) {
          # NOTE: we assume everything is in the UPPER tri!
          tmp.tt <- t(tmp)
          tmp[lower.tri(tmp)] <- tmp.tt[lower.tri(tmp.tt)]
        }
        # m.unco.idx[[offset]] <-     which(tmp > 0)
        x.unco.idx[[offset]] <- tmp[which(tmp > 0)]
      } else {
        # m.unco.idx[[offset]] <- m.free.idx[[offset]]
        x.unco.idx[[offset]] <- x.free.idx[[offset]]
      }

      # 3. general mapping between user and tmp.glist
      tmp[cbind(tmp.rep$row[idx], tmp.rep$col[idx])] <- lavpartable$id[idx]
      if (mm.symmetric[mm]) {
        tmp.tt <- t(tmp)
        tmp[lower.tri(tmp)] <- tmp.tt[lower.tri(tmp.tt)]
      }
      m.user.idx[[offset]] <- which(tmp > 0)
      x.user.idx[[offset]] <- tmp[which(tmp > 0)]

      # 4. now assign starting/fixed values
      # create empty matrix
      # FIXME: again, we may want to use sparse matrices here...
      tmp <- matrix(0.0,
        nrow = mm.rows[mm],
        ncol = mm.cols[mm]
      )
      tmp[cbind(tmp.rep$row[idx], tmp.rep$col[idx])] <- lavpartable$start[idx]
      if (mm.symmetric[mm]) {
        tmp.tt <- t(tmp)
        tmp[lower.tri(tmp)] <- tmp.tt[lower.tri(tmp.tt)]
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
        mm.names[mm] == "lambda") {
        ov.dummy.names.nox <- attr(tmp.rep, "ov.dummy.names.nox")[[g]]
        ov.dummy.names.x <- attr(tmp.rep, "ov.dummy.names.x")[[g]]
        ov.dummy.names <- c(ov.dummy.names.nox, ov.dummy.names.x)
        # define dummy latent variables
        if (length(ov.dummy.names)) {
          # in this case, lv.names will be extended with the dummys
          tmp.lv.names <- mm.dim.names$psi[[1]]
          row.tmp.idx <- match(ov.dummy.names, ov.names)
          col.tmp.idx <- match(ov.dummy.names, tmp.lv.names)
          # Fix lambda values to 1.0
          tmp[cbind(row.tmp.idx, col.tmp.idx)] <- 1.0

          ov.x.dummy.ov.idx[[g]] <- match(ov.dummy.names.x, ov.names)
          ov.x.dummy.lv.idx[[g]] <- match(ov.dummy.names.x, tmp.lv.names)
          ov.y.dummy.ov.idx[[g]] <- match(ov.dummy.names.nox, ov.names)
          ov.y.dummy.lv.idx[[g]] <- match(ov.dummy.names.nox, tmp.lv.names)
        }
      }

      # representation specific
      if (lavoptions$representation == "LISREL" && mm.names[mm] == "delta") {
        # only categorical values are listed in the lavpartable
        # but all remaining values should be 1.0
        idx <- which(tmp[, 1L] == 0.0)
        tmp[idx, 1L] <- 1.0
      }

      # representation specific
      if (lavoptions$representation == "RAM" && mm.names[mm] == "ov.idx") {
        tmp[1, ] <- attr(tmp.rep, "ov.idx")[[g]]
      }

      # assign matrix to tmp.glist
      tmp.glist[[offset]] <- tmp
    } # mm

    # efa related info
    if (nefa > 0L) {
      ov.efa.idx[[g]] <- vector("list", length = nefa)
      lv.efa.idx[[g]] <- vector("list", length = nefa)
      for (set in seq_len(nefa)) {
        # determine ov idx for this set
        ov.efa <-
          unique(lavpartable$rhs[lavpartable$op == "=~" &
            lavpartable$block == g &
            lavpartable$efa == efa.values[set]])
        ov.efa.idx[[g]][[set]] <- match(ov.efa, ov.names)

        lv.efa <-
          unique(lavpartable$lhs[lavpartable$op == "=~" &
            lavpartable$block == g &
            lavpartable$efa == efa.values[set]])
        lv.efa.idx[[g]][[set]] <- match(lv.efa, lv.names)
      }
      names(ov.efa.idx[[g]]) <- efa.values
      names(lv.efa.idx[[g]]) <- efa.values
    } # efa
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
    lambda.idx <- which(names(tmp.glist) == "lambda")[2L]

    # find rows/cols
    b.names <- paste0("b", ov.names) ## ad-hoc assumption!!!
    tmp.cols <- match(b.names, tmp.lv.names)
    tmp.rows <- seq_len(nvar[2])
    stopifnot(length(tmp.cols) == length(tmp.rows))
    tmp.glist[[lambda.idx]][cbind(tmp.rows, tmp.cols)] <-
      lavoptions$tech.muml.scale
  }

  # which free parameters are observed variances?
  ov.names <- vnames(lavpartable, "ov")
  x.free.var.idx <- lavpartable$free[lavpartable$free &
    # !duplicated(lavpartable$free) &
    lavpartable$lhs %in% ov.names &
    lavpartable$op == "~~" &
    lavpartable$lhs == lavpartable$rhs]

  rv.lv <- rv.ov <- list()
  if (multilevel) {
    # store information about random slopes (if any)
    lv.names <- lav_partable_vnames(lavpartable, "lv")
    # we should also add splitted-y names (x) to lv.names
    # FIXME: make this work for multiple work multilevel
    level.values <- lav_partable_level_values(lavpartable)
    ovx1 <- lavNames(lavpartable, "ov.x", level = level.values[1])
    ovx2 <- lavNames(lavpartable, "ov.x", level = level.values[2])
    ovx12 <- ovx2[ovx2 %in% ovx1]
    lv.names <- c(lv.names, ovx12)

    # RV LV
    rv.idx <- which(nchar(lavpartable$rv) > 0L &
      lavpartable$level == level.values[1] &
      lavpartable$rhs %in% lv.names)
    if (length(rv.idx)) {
      rv.lv <- lapply(rv.idx, function(x) {
        c(lavpartable$lhs[x], lavpartable$rhs[x])
      })
      names(rv.lv) <- lavpartable$rv[rv.idx]
    }

    # RV OV
    rv.idx <- which(nchar(lavpartable$rv) > 0L &
      lavpartable$level == level.values[1] &
      !lavpartable$rhs %in% lv.names)
    if (length(rv.idx)) {
      rv.ov <- lapply(rv.idx, function(x) {
        c(lavpartable$lhs[x], lavpartable$rhs[x])
      })
      names(rv.ov) <- lavpartable$rv[rv.idx]
    }
  } # multilevel

  # new in 0.6-9: model properties
  modprop <- lav_model_properties(
    GLIST = tmp.glist,
    lavpartable = lavpartable,
    nmat = nmat,
    m.free.idx = m.free.idx
  )

  tmp.model <- new("lavModel",
    GLIST = tmp.glist,
    dimNames = dim.names,
    isSymmetric = is.symmetric,
    mmSize = mm.size,
    representation = lavoptions$representation,
    modprop = modprop,
    meanstructure = meanstructure,
    correlation = correlation,
    categorical = categorical,
    multilevel = multilevel,
    link = lavoptions$link,
    nblocks = nblocks,
    ngroups = ngroups, # breaks rsem????
    nefa = nefa,
    group.w.free = group.w.free,
    nmat = nmat,
    nvar = nvar,
    num.idx = num.idx,
    th.idx = th.idx,
    nx.free = max(lavpartable$free),
    nx.unco = if (is.null(lavpartable$unco)) {
      max(lavpartable$free)
    } else {
      max(lavpartable$unco)
    },
    nx.user = max(lavpartable$id),
    m.free.idx = m.free.idx,
    x.free.idx = x.free.idx,
    x.free.var.idx = x.free.var.idx,
    # m.unco.idx=m.unco.idx,
    x.unco.idx = x.unco.idx,
    m.user.idx = m.user.idx,
    x.user.idx = x.user.idx,
    x.def.idx = which(lavpartable$op == ":="),
    x.ceq.idx = which(lavpartable$op == "=="),
    x.cin.idx = which(lavpartable$op == ">" | lavpartable$op == "<"),
    ceq.simple.only = tmp.con$ceq.simple.only,
    ceq.simple.K = tmp.con$ceq.simple.K,
    eq.constraints = tmp.con$ceq.linear.only.flag,
    eq.constraints.K = tmp.con$ceq.JAC.NULL,
    eq.constraints.k0 = tmp.con$ceq.rhs.NULL,
    def.function = tmp.con$def.function,
    ceq.function = tmp.con$ceq.function,
    ceq.JAC = tmp.con$ceq.JAC,
    ceq.rhs = tmp.con$ceq.rhs,
    ceq.jacobian = tmp.con$ceq.jacobian,
    ceq.linear.idx = tmp.con$ceq.linear.idx,
    ceq.nonlinear.idx = tmp.con$ceq.nonlinear.idx,
    cin.function = tmp.con$cin.function,
    cin.JAC = tmp.con$cin.JAC,
    cin.rhs = tmp.con$cin.rhs,
    cin.jacobian = tmp.con$cin.jacobian,
    cin.linear.idx = tmp.con$cin.linear.idx,
    cin.nonlinear.idx = tmp.con$cin.nonlinear.idx,
    con.jac = con.jac,
    con.lambda = con.lambda,
    nexo = nexo,
    fixed.x = lavoptions$fixed.x,
    conditional.x = lavoptions$conditional.x,
    parameterization = lavoptions$parameterization,
    ov.x.dummy.ov.idx = ov.x.dummy.ov.idx,
    ov.x.dummy.lv.idx = ov.x.dummy.lv.idx,
    ov.y.dummy.ov.idx = ov.y.dummy.ov.idx,
    ov.y.dummy.lv.idx = ov.y.dummy.lv.idx,
    ov.efa.idx = ov.efa.idx,
    lv.efa.idx = lv.efa.idx,
    rv.lv = rv.lv,
    rv.ov = rv.ov,
    estimator = lavoptions$estimator,
    estimator.args = lavoptions$estimator.args
  )

  if (lavoptions$debug) {
    cat("lavaan lavoptions$debug: lavaanModel\n")
    print(str(tmp.model))
    print(tmp.model@GLIST)
  }

  tmp.model
}

# for backwards compatibility
# tmp.model <- lav_model
