# inspect a lavaanList object

lav_lavaanlist_inspect <- function(object, what = "free", ...) {
  dotdotdot <- list(...)
  if (length(dotdotdot) > 0L) {
    for (j in seq_along(dotdotdot)) {
      lav_msg_warn(gettextf(
        "Unknown argument %s for %s", sQuote(names(dotdotdot)[j]),
        sQuote("inspect"))
      )
    }
  }
  lavListInspect(
    object = object,
    what = what,
    add.labels = TRUE,
    add.class = TRUE,
    drop.list.single.group = TRUE
  )
}

# the `tech' version: no labels, full matrices, ... for further processing
lav_lavaanlist_lavtech <- function(object,
                               what = "free",                            # nolint start
                               add.labels = FALSE,
                               add.class = FALSE,
                               list.by.group = FALSE,
                               drop.list.single.group = FALSE) {         # nolint end
  lavListInspect(
    object = object, what = what,
    add.labels = add.labels, add.class = add.class,
    list.by.group = list.by.group,
    drop.list.single.group = drop.list.single.group
  )
}

lavListTech <- function(object,                                          # nolint start
                        what = "free",
                        add.labels = FALSE,
                        add.class = FALSE,
                        list.by.group = FALSE,
                        drop.list.single.group = FALSE) {                # nolint end
  lavListInspect(
    object = object, what = what,
    add.labels = add.labels, add.class = add.class,
    list.by.group = list.by.group,
    drop.list.single.group = drop.list.single.group
  )
}

# just in case some uses lavInspect on a lavaanList object
lav_lavaanlist_lavinspect <- function(object,
                                  what = "free",
                                  add.labels = TRUE,                     # nolint start
                                  add.class = TRUE,
                                  list.by.group = TRUE,
                                  drop.list.single.group = TRUE) {       # nolint end
  lavListInspect(
    object = object, what = what,
    add.labels = add.labels, add.class = add.class,
    list.by.group = list.by.group,
    drop.list.single.group = drop.list.single.group
  )
}

lavListInspect <- function(object,                                       # nolint start
                           what = "free",
                           add.labels = TRUE,
                           add.class = TRUE,
                           list.by.group = TRUE,
                           drop.list.single.group = TRUE) {              # nolint end
  # object must inherit from class lavaanList
  stopifnot(inherits(object, "lavaanList"))

  # only a single argument
  if (length(what) > 1) {
    lav_msg_stop(gettext(
      "`what' arguments contains multiple arguments; only one is allowed"))
  }

  # be case insensitive
  what <- tolower(what)


  #### model matrices, with different contents ####
  if (what == "free") {
    lav_lavaanlist_inspect_mms(object,
      what = "free",
      type = "free", add_labels = add.labels, add_class = add.class,
      list_by_group = list.by.group,
      drop_list_single_group = drop.list.single.group
    )
  } else if (what == "partable" || what == "user") {
    lav_lavaanlist_inspect_mms(object,
      what = "free",
      type = "partable", add_labels = add.labels, add_class = add.class,
      list_by_group = list.by.group,
      drop_list_single_group = drop.list.single.group
    )
  } else if (what == "start" || what == "starting.values") {
    lav_lavaanlist_inspect_mms(object,
      what = "start",
      add_labels = add.labels, add_class = add.class,
      list_by_group = list.by.group,
      drop_list_single_group = drop.list.single.group
    )


    #### parameter table ####
  } else if (what == "list") {
    parTable(object)

    #### data + missingness ####
  } else if (what == "ngroups") {
    object@Data@ngroups
  } else if (what == "group") {
    object@Data@group
  } else if (what == "cluster") {
    object@Data@cluster
  } else if (what == "nlevels") {
    object@Data@nlevels
  } else if (what == "nclusters") {
    lav_object_inspect_cluster_info(object,
      level = 2L,
      what = "nclusters",
      drop.list.single.group = drop.list.single.group
    )
  } else if (what == "ncluster.size") {
    lav_object_inspect_cluster_info(object,
      level = 2L,
      what = "ncluster.size",
      drop.list.single.group = drop.list.single.group
    )
  } else if (what == "cluster.size") {
    lav_object_inspect_cluster_info(object,
      level = 2L,
      what = "cluster.size",
      drop.list.single.group = drop.list.single.group
    )
  } else if (what == "cluster.id") {
    lav_object_inspect_cluster_info(object,
      level = 2L,
      what = "cluster.id",
      drop.list.single.group = drop.list.single.group
    )
  } else if (what == "cluster.idx") {
    lav_object_inspect_cluster_info(object,
      level = 2L,
      what = "cluster.idx",
      drop.list.single.group = drop.list.single.group
    )
  } else if (what == "cluster.label") {
    lav_object_inspect_cluster_info(object,
      level = 2L,
      what = "cluster.label",
      drop.list.single.group = drop.list.single.group
    )
  } else if (what == "cluster.sizes") {
    lav_object_inspect_cluster_info(object,
      level = 2L,
      what = "cluster.sizes",
      drop.list.single.group = drop.list.single.group
    )
  } else if (what == "average.cluster.size") {
    lav_object_inspect_cluster_info(object,
      level = 2L,
      what = "average.cluster.size",
      drop.list.single.group = drop.list.single.group
    )
  } else if (what == "ordered") {
    object@Data@ordered
  } else if (what == "group.label") {
    object@Data@group.label
  } else if (what == "level.label") {
    object@Data@level.label
  } else if (what == "nobs") { # only for original!
    unlist(object@Data@nobs)
  } else if (what == "norig") { # only for original!
    unlist(object@Data@norig)
  } else if (what == "ntotal") { # only for original!
    sum(unlist(object@Data@nobs))

    #### from the model object (but stable) over datasets? ####
  } else if (what == "th.idx") {
    lav_lavaanlist_inspect_th_idx(object,
      add_labels = add.labels, add_class = add.class,
      drop_list_single_group = drop.list.single.group
    )


    #### meanstructure, categorical ####
  } else if (what == "meanstructure") {
    object@Model@meanstructure
  } else if (what == "categorical") {
    object@Model@categorical
  } else if (what == "fixed.x") {
    object@Model@fixed.x
  } else if (what == "parameterization") {
    object@Model@parameterization

    # options
  } else if (what == "options" || what == "lavoptions") {
    object@Options

    # call
  } else if (what == "call") {
    as.list(object@call)

    #### not found ####
  } else {
    lav_msg_stop(gettextf(
      "unknown `what' argument in inspect function: `%s'", what))
  }
}


lav_lavaanlist_inspect_start <- function(object) {
  # from 0.5-19, they are in the partable
  if (!is.null(object@ParTable$start)) {
    out <- object@ParTable$start
  } else {
    # in < 0.5-19, we should look in @Fit@start
    out <- object@Fit@start
  }

  out
}

lav_lavaanlist_inspect_mms <- function(
    object, what = "free",
    type = "free", add_labels = FALSE, add_class = FALSE,
    list_by_group = FALSE,
    drop_list_single_group = FALSE) {
  glist <- object@Model@GLIST

  for (mm in seq_along(glist)) {
    if (add_labels) {
      dimnames(glist[[mm]]) <- object@Model@dimNames[[mm]]
    }

    if (what == "free") {
      # fill in free parameter counts
      if (type == "free") {
        m_el_idx <- object@Model@m.free.idx[[mm]]
        x_el_idx <- object@Model@x.free.idx[[mm]]
        # } else if(type == "unco") {
        #    m.el.idx <- object@Model@m.unco.idx[[mm]]
        #    x.el.idx <- object@Model@x.unco.idx[[mm]]
      } else if (type == "partable") {
        m_el_idx <- object@Model@m.user.idx[[mm]]
        x_el_idx <- object@Model@x.user.idx[[mm]]
      } else {
        lav_msg_stop(gettextf("unknown type argument: %s", type))
      }
      # erase everything
      glist[[mm]][, ] <- 0.0
      glist[[mm]][m_el_idx] <- x_el_idx
    } else if (what == "start") {
      # fill in starting values
      m_user_idx <- object@Model@m.user.idx[[mm]]
      x_user_idx <- object@Model@x.user.idx[[mm]]
      start_1 <- lav_lavaanlist_inspect_start(object)
      glist[[mm]][m_user_idx] <- start_1[x_user_idx]
    }

    # class
    if (add_class) {
      if (object@Model@isSymmetric[mm]) {
        class(glist[[mm]]) <- c("lavaan.matrix.symmetric", "matrix")
      } else {
        class(glist[[mm]]) <- c("lavaan.matrix", "matrix")
      }
    }
  }

  # try to reflect `equality constraints'
  con_flag <- FALSE
  if (what == "free" && object@Model@eq.constraints) {
    # extract constraints from parameter table
    pt_1 <- parTable(object)
    con_1 <- pt_1[pt_1$op %in% c("==", "<", ">"), c("lhs", "op", "rhs")]
    rownames(con_1) <- NULL

    # replace 'labels' by parameter numbers
    id <- lav_partable_constraints_label_id(pt_1)
    label <- names(id)
    for (con in seq_len(nrow(con_1))) {
      # lhs
      lhs_labels <- all.vars(as.formula(paste("~", con_1[con, "lhs"])))

      if (length(lhs_labels) > 0L) {
        # par id
        lhs_freeid <- id[match(lhs_labels, label)]

        # substitute
        tmp <- con_1[con, "lhs"]
        for (pat in seq_along(lhs_labels)) {
          tmp <- sub(lhs_labels[pat], lhs_freeid[pat], tmp)
        }
        con_1[con, "lhs"] <- tmp
      }

      # rhs
      rhs_labels <- all.vars(as.formula(paste("~", con_1[con, "rhs"])))

      if (length(rhs_labels) > 0L) {
        # par id
        rhs_freeid <- id[match(rhs_labels, label)]
        # substitute
        tmp <- con_1[con, "rhs"]
        for (pat in seq_along(rhs_labels)) {
          tmp <- sub(rhs_labels[pat], rhs_freeid[pat], tmp)
        }
        con_1[con, "rhs"] <- tmp
      }
    } # con

    # add this info at the top
    # GLIST <- c(constraints = list(CON), GLIST)
    # no, not a good idea, it does not work with list.by.group

    # add it as a 'header' attribute?
    attr(con_1, "header") <- "Note: model contains equality constraints:"
    con_flag <- TRUE
  }

  # should we group them per group?
  if (list_by_group) {
    lavmodel <- object@Model
    nmat <- lavmodel@nmat

    out <- vector("list", length = object@Data@ngroups)
    for (g in 1:object@Data@ngroups) {
      # which mm belong to group g?
      mm_in_group <- 1:nmat[g] + cumsum(c(0, nmat))[g]

      out[[g]] <- glist[mm_in_group]
    }

    if (object@Data@ngroups == 1L && drop_list_single_group) {
      out <- out[[1]]
    } else {
      if (length(object@Data@group.label) > 0L) {
        names(out) <- unlist(object@Data@group.label)
      }
    }
  } else {
    out <- glist
  }

  # header
  if (con_flag) {
    attr(out, "header") <- con_1
  }

  # lavaan.list
  if (add_class) {
    class(out) <- c("lavaan.list", "list")
  }

  out
}

lav_lavaanlist_inspect_th_idx <- function(
    object,
    add_labels = FALSE, add_class = FALSE, drop_list_single_group = FALSE) {
  # thresholds idx -- usually, we get it from SampleStats
  # but fortunately, there is a copy in Model, but no names...
  out <- object@Model@th.idx

  # nblocks
  nblocks <- length(out)

  # labels + class
  for (b in seq_len(nblocks)) {
    # if(add.labels && length(OUT[[b]]) > 0L) {
    #    names(OUT[[b]]) <- object@SampleStats@th.names[[b]]
    # }
    if (add_class && !is.null(out[[b]])) {
      class(out[[b]]) <- c("lavaan.vector", "numeric")
    }
  }

  if (nblocks == 1L && drop_list_single_group) {
    out <- out[[1]]
  } else {
    if (object@Data@nlevels == 1L &&
      length(object@Data@group.label) > 0L) {
      names(out) <- unlist(object@Data@group.label)
    } else if (object@Data@nlevels > 1L &&
      length(object@Data@group.label) == 0L) {
      names(out) <- object@Data@level.label
    }
  }

  out
}
