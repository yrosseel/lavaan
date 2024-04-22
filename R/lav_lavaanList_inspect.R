# inspect a lavaanList object

inspect.lavaanList <- function(object, what = "free", ...) {
  lavListInspect(
    object = object,
    what = what,
    add.labels = TRUE,
    add.class = TRUE,
    drop.list.single.group = TRUE
  )
}

# the `tech' version: no labels, full matrices, ... for further processing
lavTech.lavaanList <- function(object,
                               what = "free",
                               add.labels = FALSE,
                               add.class = FALSE,
                               list.by.group = FALSE,
                               drop.list.single.group = FALSE) {
  lavListInspect(
    object = object, what = what,
    add.labels = add.labels, add.class = add.class,
    list.by.group = list.by.group,
    drop.list.single.group = drop.list.single.group
  )
}

lavListTech <- function(object,
                        what = "free",
                        add.labels = FALSE,
                        add.class = FALSE,
                        list.by.group = FALSE,
                        drop.list.single.group = FALSE) {
  lavListInspect(
    object = object, what = what,
    add.labels = add.labels, add.class = add.class,
    list.by.group = list.by.group,
    drop.list.single.group = drop.list.single.group
  )
}

# just in case some uses lavInspect on a lavaanList object
lavInspect.lavaanList <- function(object,
                                  what = "free",
                                  add.labels = TRUE,
                                  add.class = TRUE,
                                  list.by.group = TRUE,
                                  drop.list.single.group = TRUE) {
  lavListInspect(
    object = object, what = what,
    add.labels = add.labels, add.class = add.class,
    list.by.group = list.by.group,
    drop.list.single.group = drop.list.single.group
  )
}

lavListInspect <- function(object,
                           what = "free",
                           add.labels = TRUE,
                           add.class = TRUE,
                           list.by.group = TRUE,
                           drop.list.single.group = TRUE) {
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
    lav_lavaanList_inspect_modelmatrices(object,
      what = "free",
      type = "free", add.labels = add.labels, add.class = add.class,
      list.by.group = list.by.group,
      drop.list.single.group = drop.list.single.group
    )
  } else if (what == "partable" || what == "user") {
    lav_lavaanList_inspect_modelmatrices(object,
      what = "free",
      type = "partable", add.labels = add.labels, add.class = add.class,
      list.by.group = list.by.group,
      drop.list.single.group = drop.list.single.group
    )
  } else if (what == "start" || what == "starting.values") {
    lav_lavaanList_inspect_modelmatrices(object,
      what = "start",
      add.labels = add.labels, add.class = add.class,
      list.by.group = list.by.group,
      drop.list.single.group = drop.list.single.group
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
    lav_lavaanList_inspect_th_idx(object,
      add.labels = add.labels, add.class = add.class,
      drop.list.single.group = drop.list.single.group
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


lav_lavaanList_inspect_start <- function(object) {
  # from 0.5-19, they are in the partable
  if (!is.null(object@ParTable$start)) {
    OUT <- object@ParTable$start
  } else {
    # in < 0.5-19, we should look in @Fit@start
    OUT <- object@Fit@start
  }

  OUT
}

lav_lavaanList_inspect_modelmatrices <- function(
    object, what = "free",
    type = "free", add.labels = FALSE, add.class = FALSE,
    list.by.group = FALSE,
    drop.list.single.group = FALSE) {
  GLIST <- object@Model@GLIST

  for (mm in 1:length(GLIST)) {
    if (add.labels) {
      dimnames(GLIST[[mm]]) <- object@Model@dimNames[[mm]]
    }

    if (what == "free") {
      # fill in free parameter counts
      if (type == "free") {
        m.el.idx <- object@Model@m.free.idx[[mm]]
        x.el.idx <- object@Model@x.free.idx[[mm]]
        # } else if(type == "unco") {
        #    m.el.idx <- object@Model@m.unco.idx[[mm]]
        #    x.el.idx <- object@Model@x.unco.idx[[mm]]
      } else if (type == "partable") {
        m.el.idx <- object@Model@m.user.idx[[mm]]
        x.el.idx <- object@Model@x.user.idx[[mm]]
      } else {
        lav_msg_stop(gettextf("unknown type argument: %s", type))
      }
      # erase everything
      GLIST[[mm]][, ] <- 0.0
      GLIST[[mm]][m.el.idx] <- x.el.idx
    } else if (what == "start") {
      # fill in starting values
      m.user.idx <- object@Model@m.user.idx[[mm]]
      x.user.idx <- object@Model@x.user.idx[[mm]]
      START <- lav_lavaanList_inspect_start(object)
      GLIST[[mm]][m.user.idx] <- START[x.user.idx]
    }

    # class
    if (add.class) {
      if (object@Model@isSymmetric[mm]) {
        class(GLIST[[mm]]) <- c("lavaan.matrix.symmetric", "matrix")
      } else {
        class(GLIST[[mm]]) <- c("lavaan.matrix", "matrix")
      }
    }
  }

  # try to reflect `equality constraints'
  con.flag <- FALSE
  if (what == "free" && object@Model@eq.constraints) {
    # extract constraints from parameter table
    PT <- parTable(object)
    CON <- PT[PT$op %in% c("==", "<", ">"), c("lhs", "op", "rhs")]
    rownames(CON) <- NULL

    # replace 'labels' by parameter numbers
    ID <- lav_partable_constraints_label_id(PT)
    LABEL <- names(ID)
    for (con in 1:nrow(CON)) {
      # lhs
      LHS.labels <- all.vars(as.formula(paste("~", CON[con, "lhs"])))

      if (length(LHS.labels) > 0L) {
        # par id
        LHS.freeid <- ID[match(LHS.labels, LABEL)]

        # substitute
        tmp <- CON[con, "lhs"]
        for (pat in 1:length(LHS.labels)) {
          tmp <- sub(LHS.labels[pat], LHS.freeid[pat], tmp)
        }
        CON[con, "lhs"] <- tmp
      }

      # rhs
      RHS.labels <- all.vars(as.formula(paste("~", CON[con, "rhs"])))

      if (length(RHS.labels) > 0L) {
        # par id
        RHS.freeid <- ID[match(RHS.labels, LABEL)]
        # substitute
        tmp <- CON[con, "rhs"]
        for (pat in 1:length(RHS.labels)) {
          tmp <- sub(RHS.labels[pat], RHS.freeid[pat], tmp)
        }
        CON[con, "rhs"] <- tmp
      }
    } # con

    # add this info at the top
    # GLIST <- c(constraints = list(CON), GLIST)
    # no, not a good idea, it does not work with list.by.group

    # add it as a 'header' attribute?
    attr(CON, "header") <- "Note: model contains equality constraints:"
    con.flag <- TRUE
  }

  # should we group them per group?
  if (list.by.group) {
    lavmodel <- object@Model
    nmat <- lavmodel@nmat

    OUT <- vector("list", length = object@Data@ngroups)
    for (g in 1:object@Data@ngroups) {
      # which mm belong to group g?
      mm.in.group <- 1:nmat[g] + cumsum(c(0, nmat))[g]
      mm.names <- names(GLIST[mm.in.group])

      OUT[[g]] <- GLIST[mm.in.group]
    }

    if (object@Data@ngroups == 1L && drop.list.single.group) {
      OUT <- OUT[[1]]
    } else {
      if (length(object@Data@group.label) > 0L) {
        names(OUT) <- unlist(object@Data@group.label)
      }
    }
  } else {
    OUT <- GLIST
  }

  # header
  if (con.flag) {
    attr(OUT, "header") <- CON
  }

  # lavaan.list
  if (add.class) {
    class(OUT) <- c("lavaan.list", "list")
  }

  OUT
}

lav_lavaanList_inspect_th_idx <- function(
    object,
    add.labels = FALSE, add.class = FALSE, drop.list.single.group = FALSE) {
  # thresholds idx -- usually, we get it from SampleStats
  # but fortunately, there is a copy in Model, but no names...
  OUT <- object@Model@th.idx

  # nblocks
  nblocks <- length(OUT)

  # labels + class
  for (b in seq_len(nblocks)) {
    # if(add.labels && length(OUT[[b]]) > 0L) {
    #    names(OUT[[b]]) <- object@SampleStats@th.names[[b]]
    # }
    if (add.class && !is.null(OUT[[b]])) {
      class(OUT[[b]]) <- c("lavaan.vector", "numeric")
    }
  }

  if (nblocks == 1L && drop.list.single.group) {
    OUT <- OUT[[1]]
  } else {
    if (object@Data@nlevels == 1L &&
      length(object@Data@group.label) > 0L) {
      names(OUT) <- unlist(object@Data@group.label)
    } else if (object@Data@nlevels > 1L &&
      length(object@Data@group.label) == 0L) {
      names(OUT) <- object@Data@level.label
    }
  }

  OUT
}
