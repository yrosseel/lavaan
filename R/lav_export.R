# export `lavaan' lav model description to third-party software
#

lav_export <- function(object, target = "lavaan", prefix = "sem",
                      dir_name = "lav_export", export = TRUE, ...) {
  dotdotdot <- list(...)
  lav_adapt_func(environment(), dotdotdot, NULL)

  stopifnot(inherits(object, "lavaan"))
  # check object
  object <- lav_object_check_version(object)
  target <- tolower(target)

  # check for conditional.x = TRUE
  # if(object@Model@conditional.x) {
  #    stop("lavaan ERROR: this function is not (yet)
  #           available if conditional.x = TRUE")
  # }

  ngroups <- object@Data@ngroups

  # 2. create syntax file
  if (target == "lavaan") {
    out <- lav_export_lavaan(object)
  } else if (target == "mplus") {
    # a single combined data file (with GRP/CLUS columns when needed)
    data_file <- paste(prefix, ".", target, ".raw", sep = "")
    out <- lav_export_mplus(object, data_file = data_file)
  } else if (target == "lisrel") {
    out <- lav_export_lisrel(object)
  } else if (target == "eqs") {
    out <- lav_export_eqs(object)
  } else if (target == "openmx") {
    out <- lav_export_openmx(object)
  } else {
    lav_msg_stop(gettextf("target %s has not been implemented yet", target))
  }

  # export to file?
  if (export) {
    if (!dir.exists(dir_name)) {
      dir.create(path = dir_name, recursive = TRUE)
    }
    input_file <- paste(dir_name, "/", prefix, ".", target, ".in", sep = "")
    cat(out, file = input_file, sep = "")

    # write data (if available)
    if (identical(object@Data@data.type, "full")) {
      data_combined <- lav_export_data(object)
      write.table(data_combined,
        file = paste(dir_name, "/", data_file, sep = ""),
        na = "-999999",
        col.names = FALSE, row.names = FALSE, quote = FALSE
      )
    } else if (identical(object@Data@data.type, "moment")) {
      data_1 <- object@SampleStats@cov[[1L]]
      write.table(data_1,
        file = paste(dir_name, "/", data_file, sep = ""),
        na = "-999999",
        col.names = FALSE, row.names = FALSE, quote = FALSE
      )
    } else {
      lav_msg_warn(gettext("no data available"))
    }
    return(invisible(out))
  }

  out
}

# build a single combined data frame from a fitted object, with columns in the
# order expected by the generated input: observed variables, [weight], [GRP],
# [CLUS]. Variable names are sanitized (dots -> underscores).
lav_export_data <- function(object) {
  ngroups <- object@Data@ngroups
  nlevels <- object@Data@nlevels
  # single-level clustered data (cluster= but no level: structure): the cluster
  # column is needed for Mplus TYPE=COMPLEX (see lav_export_mplus_header)
  clustered <- nlevels == 1L && length(object@Data@cluster) > 0L
  ov.names <- lav_pt_vnames(object@ParTable, "ov")
  weight.name <- object@Data@sampling.weights

  parts <- vector("list", ngroups)
  clus.offset <- 0L
  for (g in seq_len(ngroups)) {
    X <- object@Data@X[[g]]
    colnames(X) <- object@Data@ov.names[[g]]
    df <- as.data.frame(X, stringsAsFactors = FALSE)
    eXo <- object@Data@eXo[[g]]
    if (!is.null(eXo) && ncol(eXo) > 0L) {
      colnames(eXo) <- object@Data@ov.names.x[[g]]
      newc <- setdiff(colnames(eXo), colnames(df))
      if (length(newc)) {
        df <- cbind(df, as.data.frame(eXo[, newc, drop = FALSE]))
      }
    }
    # reorder observed columns to match the partable ov order
    df <- df[, ov.names, drop = FALSE]
    if (length(weight.name) && !is.null(object@Data@weights[[g]])) {
      df[[weight.name]] <- object@Data@weights[[g]]
    }
    if (ngroups > 1L) df[["GRP"]] <- g
    if (nlevels > 1L || clustered) {
      # cluster.idx restarts at 1 within each group; offset so that cluster ids
      # are unique across groups (Mplus requires this)
      cidx <- object@Data@Lp[[g]]$cluster.idx[[2L]]
      df[["CLUS"]] <- cidx + clus.offset
      clus.offset <- clus.offset + max(cidx)
    }
    parts[[g]] <- df
  }

  do.call(rbind, parts)
}
lavExport <- lav_alias("lav_export") # synonym #nolint


lav_export_check <- function(lav) {
  if (inherits(lav, "lavaan")) {
    lav <- lav@ParTable
  } else if (is.list(lav)) {
    # nothing to do
  } else {
    lav_msg_stop(gettext("lav must be of class `lavaan' or a parTable"))
  }

  # check syntax
  if (is.null(lav$ustart)) lav$ustart <- lav$est

  # check if free is missing
  if (is.null(lav$free)) lav$free <- rep(0L, length(lav$ustart))

  # check if label is missing
  if (is.null(lav$label)) lav$label <- rep("", length(lav$ustart))

  # check if group is missing
  if (is.null(lav$group)) lav$group <- rep(1L, length(lav$ustart))

  # if eq.id not all zero, create labels instead
  # if(!is.null(lav$eq.id) && !all(lav$eq.id == 0L)) {
  #    lav$label <- paste("p",as.character(lav$eq.id), sep="")
  #    lav$label[lav$label == "p0"] <- ""
  # }

  lav
}

# Regenerate (compact) lavaan model syntax from a fitted lavaan object
# (or a parameter table).
#
# Design goals:
#  - the regenerated syntax, when refitted with the same data and the same
#    fitting options, reproduces an equivalent model (same free/fixed pattern,
#    same equality constraints, same fixed values)
#  - the syntax is as compact as possible: rhs terms that share the same lhs and
#    operator are combined on a single line (eg `f =~ x1 + x2 + x3`)
#  - multiple groups are written using `group:` blocks; multilevel models using
#    `level:` blocks (nested within `group:` blocks if both are present)
lav_export_lavaan <- function(object, ...) {
  # accept a fitted lavaan object or a parameter table
  if (inherits(object, "lavaan")) {
    pt        <- as.data.frame(object@ParTable, stringsAsFactors = FALSE)
    ngroups   <- object@Data@ngroups
    nlevels   <- object@Data@nlevels
  } else if (is.list(object)) {
    pt <- as.data.frame(object, stringsAsFactors = FALSE)
    ngroups   <- if (is.null(pt$group)) 1L else max(pt$group)
    nlevels   <- if (is.null(pt$level)) 1L else max(pt$level)
  } else {
    lav_msg_stop(gettext("object must be of class \"lavaan\" or a parameter table"))
  }
  nn <- length(pt$lhs)

  # ensure the columns we rely on exist
  if (is.null(pt$free))   pt$free   <- rep(0L, nn)
  if (is.null(pt$ustart)) {
    pt$ustart <- if (!is.null(pt$est)) pt$est else rep(as.numeric(NA), nn)
  }
  if (is.null(pt$label))  pt$label  <- rep("", nn)
  if (is.null(pt$plabel)) pt$plabel <- rep("", nn)
  if (is.null(pt$block))  pt$block  <- rep(1L, nn)
  if (is.null(pt$group))  pt$group  <- ifelse(pt$block == 0L, 0L, 1L)
  if (is.null(pt$level))  pt$level  <- ifelse(pt$block == 0L, 0L, 1L)
  pt$label[is.na(pt$label)] <- ""

  # operator categories
  con.op <- c("==", "<", ">", ":=")
  is.con <- pt$op %in% con.op

  # ---- 1. resolve equality classes -> clean display labels -----------------
  display.label <- lav_export_lavaan_labels(pt, is.con)

  # ---- 2. emit the model, block by block -----------------------------------
  body <- character(0L)
  for (g in seq_len(ngroups)) {
    if (ngroups > 1L) {
      # use the numeric block index: group labels may contain characters that
      # are not valid in the block header syntax (eg "Grant-White")
      body <- c(body, "", paste0("group: ", g))
    }
    for (l in seq_len(nlevels)) {
      if (nlevels > 1L) {
        body <- c(body, paste0("level: ", l))
      }
      block.idx <- which(!is.con & pt$group == g & pt$level == l)
      body <- c(body,
                lav_export_lavaan_block(pt, block.idx, display.label))
    }
  }

  # ---- 3. constraints and defined parameters (block-independent) ------------
  con.lines <- lav_export_lavaan_constraints(pt, is.con, display.label)
  if (length(con.lines)) {
    body <- c(body, "", con.lines)
  }

  header <- "# model syntax autogenerated by lavExport()\n"
  out <- paste0(header, paste(body, collapse = "\n"), "\n")
  class(out) <- c("lavaan.character", "character")
  out
}

# format a fixed numeric value compactly but without loss of precision
lav_export_num <- function(x) {
  if (is.na(x)) {
    return("NA")
  }
  out <- formatC(x, digits = 17L, format = "g")
  # round-trip check: use the shortest representation that recovers x
  for (d in 1:16) {
    short <- formatC(x, digits = d, format = "g")
    if (isTRUE(all.equal(as.numeric(short), x))) {
      out <- short
      break
    }
  }
  trimws(out)
}

# determine, for each parameter row, the clean label to display (or "")
# parameters that are constrained equal (via a shared label or via a pairwise
# `==` constraint) get a single shared display label
lav_export_lavaan_labels <- function(pt, is.con) {
  nn <- length(pt$lhs)
  par.idx <- which(!is.con)

  # union-find over parameter rows (indexed by position in pt)
  parent <- seq_len(nn)
  find <- function(i) {
    while (parent[i] != i) {
      parent[i] <<- parent[parent[i]]
      i <- parent[i]
    }
    i
  }
  union <- function(i, j) {
    ri <- find(i)
    rj <- find(j)
    if (ri != rj) parent[rj] <<- ri
  }

  # lookup: plabel -> row, (non-empty) label -> rows
  plab2row <- integer(0L)
  if (any(nzchar(pt$plabel))) {
    idx <- which(nzchar(pt$plabel))
    plab2row <- idx
    names(plab2row) <- pt$plabel[idx]
  }
  resolve <- function(token) {
    # return parameter row(s) referenced by a token (plabel or label)
    if (token %in% names(plab2row)) {
      return(unname(plab2row[token]))
    }
    which(!is.con & pt$label == token)
  }

  # (a) union rows that share the same non-empty label
  lab.tab <- pt$label[par.idx]
  for (lab in unique(lab.tab[nzchar(lab.tab)])) {
    rows <- par.idx[lab.tab == lab]
    if (length(rows) > 1L) {
      for (k in 2:length(rows)) union(rows[1L], rows[k])
    }
  }

  # (b) union rows linked by a pairwise equality constraint between two
  #     parameter references (plabel or label); these are folded into a shared
  #     label and the `==` line itself is dropped later
  token.re <- "^[A-Za-z.][A-Za-z0-9._]*$"
  eq.idx <- which(pt$op == "==")
  for (i in eq.idx) {
    lhs <- pt$lhs[i]
    rhs <- pt$rhs[i]
    if (grepl(token.re, lhs) && grepl(token.re, rhs)) {
      lrow <- resolve(lhs)
      rrow <- resolve(rhs)
      if (length(lrow) && length(rrow)) {
        union(lrow[1L], rrow[1L])
      }
    }
  }

  # build classes
  root <- vapply(seq_len(nn), find, integer(1L))
  display.label <- rep("", nn)

  # collect names already in use (variable names + user labels) to avoid
  # collisions when generating fresh labels
  is.plabel.style <- function(x) grepl("^\\.", x)
  user.lab <- unique(pt$label[!is.con & nzchar(pt$label) &
                              !is.plabel.style(pt$label)])
  reserved <- unique(c(pt$lhs, pt$rhs, user.lab,
                       pt$lhs[pt$op == ":="]))
  gen.counter <- 0L
  gen.name <- function() {
    repeat {
      gen.counter <<- gen.counter + 1L
      cand <- paste0("p", gen.counter)
      if (!cand %in% reserved) {
        reserved <<- c(reserved, cand)
        return(cand)
      }
    }
  }

  for (r in unique(root[par.idx])) {
    members <- which(root == r & !is.con)
    if (length(members) == 0L) next
    # prefer an existing user label
    member.user.lab <- pt$label[members]
    member.user.lab <- member.user.lab[nzchar(member.user.lab) &
                                       !is.plabel.style(member.user.lab)]
    if (length(member.user.lab)) {
      lab <- member.user.lab[1L]
    } else if (length(members) > 1L) {
      lab <- gen.name()
    } else {
      lab <- ""
    }
    display.label[members] <- lab
  }

  display.label
}

# emit the (compact) model lines for one block
lav_export_lavaan_block <- function(pt, block.idx, display.label) {
  if (length(block.idx) == 0L) {
    return(character(0L))
  }

  out <- character(0L)

  # modifier-decorated rhs term for a single row; returns NA to drop the term
  term <- function(i) {
    rhs <- pt$rhs[i]
    lab <- display.label[i]
    if (pt$free[i] == 0L) {
      if (is.na(pt$ustart[i])) {
        return(NA_character_) # eg fixed.x (co)variance/mean: refit re-adds it
      }
      val <- lav_export_num(pt$ustart[i])
      if (nzchar(lab)) {
        # fixed value *and* a referenced label: stack both modifiers
        return(paste0(val, "*", lab, "*", rhs))
      }
      return(paste0(val, "*", rhs))
    }
    if (nzchar(lab)) {
      return(paste0(lab, "*", rhs))
    }
    rhs
  }

  # ordered operator groups for stable, readable output
  op.order <- c("=~", "<~", "~", "~1", "~~", "|", "~*~")
  ops <- intersect(op.order, unique(pt$op[block.idx]))
  ops <- c(ops, setdiff(unique(pt$op[block.idx]), op.order))

  for (op in ops) {
    rows <- block.idx[pt$op[block.idx] == op]
    if (op == "~1") {
      # one line per lhs: `x1 ~ <mod>1`
      for (lhs in unique(pt$lhs[rows])) {
        i <- rows[pt$lhs[rows] == lhs][1L]
        lab <- display.label[i]
        if (pt$free[i] == 0L) {
          if (is.na(pt$ustart[i])) next
          mod <- paste0(lav_export_num(pt$ustart[i]), "*")
          if (nzchar(lab)) mod <- paste0(mod, lab, "*")
        } else {
          mod <- if (nzchar(lab)) paste0(lab, "*") else ""
        }
        out <- c(out, paste0(lhs, " ~ ", mod, "1"))
      }
      next
    }
    # group by lhs and combine rhs terms with `+`
    for (lhs in unique(pt$lhs[rows])) {
      lrows <- rows[pt$lhs[rows] == lhs]
      if (op == "~~") {
        # list the variance (diagonal) before any covariances
        diag.first <- order(pt$rhs[lrows] != lhs)
        lrows <- lrows[diag.first]
      }
      terms <- vapply(lrows, term, character(1L))
      terms <- terms[!is.na(terms)]
      if (length(terms) == 0L) next
      out <- c(out, paste0(lhs, " ", op, " ", paste(terms, collapse = " + ")))
    }
  }

  out
}

# determine the defined parameters (:=) and constraints (==, <, >) to emit, with
# all references to plabels / old labels rewritten to the clean display labels.
# pairwise equality constraints between two parameters are dropped (they are
# folded into shared labels). returns a data.frame with columns lhs, op, rhs.
lav_export_constraints_df <- function(pt, is.con, display.label) {
  con.idx <- which(is.con)
  empty <- data.frame(lhs = character(0L), op = character(0L),
                      rhs = character(0L), stringsAsFactors = FALSE)
  if (length(con.idx) == 0L) {
    return(empty)
  }

  # token translation map: plabel/old-label -> display label
  trans <- character(0L)
  has.disp <- which(nzchar(display.label))
  if (length(has.disp)) {
    if (any(nzchar(pt$plabel[has.disp]))) {
      pidx <- has.disp[nzchar(pt$plabel[has.disp])]
      trans[pt$plabel[pidx]] <- display.label[pidx]
    }
    if (any(nzchar(pt$label[has.disp]))) {
      lidx <- has.disp[nzchar(pt$label[has.disp])]
      trans[pt$label[lidx]] <- display.label[lidx]
    }
  }
  rewrite <- function(expr) {
    if (length(trans) == 0L || is.na(expr) || !nzchar(expr)) {
      return(expr)
    }
    for (from in names(trans)) {
      to <- trans[[from]]
      if (identical(from, to)) next
      pat <- paste0("(?<![A-Za-z0-9._])", lav_export_escape(from),
                    "(?![A-Za-z0-9._])")
      expr <- gsub(pat, to, expr, perl = TRUE)
    }
    expr
  }

  # which `==` rows are pairwise equalities already folded into labels?
  token.re <- "^[A-Za-z.][A-Za-z0-9._]*$"
  is.param.ref <- function(token) {
    token %in% pt$plabel[!is.con] ||
      token %in% pt$label[!is.con & nzchar(pt$label)]
  }

  lhs <- rhs <- op <- character(0L)
  for (i in con.idx) {
    o <- pt$op[i]
    if (o == "==" &&
        grepl(token.re, pt$lhs[i]) && grepl(token.re, pt$rhs[i]) &&
        is.param.ref(pt$lhs[i]) && is.param.ref(pt$rhs[i])) {
      next # folded into a shared label
    }
    op  <- c(op, o)
    lhs <- c(lhs, if (o == ":=") pt$lhs[i] else rewrite(pt$lhs[i]))
    rhs <- c(rhs, rewrite(pt$rhs[i]))
  }
  data.frame(lhs = lhs, op = op, rhs = rhs, stringsAsFactors = FALSE)
}

# emit lavaan-syntax lines for defined parameters and constraints
lav_export_lavaan_constraints <- function(pt, is.con, display.label) {
  con.df <- lav_export_constraints_df(pt, is.con, display.label)
  if (nrow(con.df) == 0L) {
    return(character(0L))
  }
  paste0(con.df$lhs, " ", con.df$op, " ", con.df$rhs)
}

# escape regex metacharacters in a literal token (eg the dots in plabels)
lav_export_escape <- function(x) {
  gsub("([.\\\\+*?\\[\\]^$(){}=!<>|:-])", "\\\\\\1", x, perl = TRUE)
}

lav_export_lisrel <- function(lav) {
  lav <- lav_export_check(lav)
  lav_msg_stop(gettext("this function needs revision"))
}

lav_export_eqs <- function(lav) {
  lav <- lav_export_check(lav)
  lav_msg_stop(gettext("this function needs revision"))
}

lav_export_openmx <- function(lav) {
  lav <- lav_export_check(lav)
  lav_msg_stop(gettext("this function needs revision"))
}
