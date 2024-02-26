lav_lavaan_step01_ovnames_1 <- function(slotParTable, model, dotdotdot.parser) {
  # als slotPartable niet null -> copieer in flat.model
  # anders
  #   als model character 
  #     parseer model naar flat.model
  #   anders
  #     als het een "formula" is (geef warning)
  #       zet "~ x1 + x2 + x3" om naar character "f =~ x1 + x2 + x3" en parseer naar flatmodel 
  #       zet "y "~ x1 + x2 + x3" om naar character en parseer naar flatmodel 
  #       iets anders --> ***error***
  #     anders 
  #       als het een lavaan object is, flat.model uithalen via functie parTable
  #       anders 
  #         als model een list is
  #           als minimum aanwezig (kolommen lhs, op, rhs, free)
  #             flat.model = model
  #             kolom blmock eventueel overschrijven met kolom group
  #           anders --> ***error***
  #         anders
  #           als model NULL is --> ***error***
  #           TODO: wat als model bvb een function is ???  Geeft nu een fout object 'flat.model' not found
  #           
  # 1a. get ov.names and ov.names.x (per group) -- needed for lavData()
  if (!is.null(slotParTable)) {
    flat.model <- slotParTable
  } else if (is.character(model)) {
    if (is.null(dotdotdot.parser)) {
      flat.model <- lavParseModelString(model, parser = "old") # for now
    } else {
      flat.model <- lavParseModelString(model, parser = dotdotdot.parser)
    }
  } else if (inherits(model, "formula")) {
    # two typical cases:
    # 1. regression type formula
    # 2. no quotes, eg f =~ x1 + x2 + x3     TODO: this isn't a valid formula !!!
    tmp <- as.character(model)
    if (tmp[1] == "~" && length(tmp) == 2L) {
      # looks like an unquoted single factor model f =~ something
      warning("lavaan WARNING: model seems to be a formula; please enclose the model syntax between quotes")
      # create model and hope for the best
      model.bis <- paste("f =", paste(tmp, collapse = " "), sep = "")
      flat.model <- lavParseModelString(model.bis)
    } else if (tmp[1] == "~" && length(tmp) == 3L) {
      # looks like a (unquoted) regression formula
      warning("lavaan WARNING: model seems to be a formula; please enclose the model syntax between quotes")
      # create model and hope for the best
      model.bis <- paste(tmp[2], tmp[1], tmp[3])
      flat.model <- lavParseModelString(model.bis)
    } else {
      stop("lavaan ERROR: model seems to be a formula; please enclose the model syntax between quotes")
    }
  } else if (inherits(model, "lavaan")) {
    # hm, a lavaan model; let's try to extract the parameter table
    # and see what happens
    flat.model <- parTable(model)
  } else if (is.list(model)) {
    # two possibilities: either model is already lavaanified
    # or it is something else...
    
    # look for the bare minimum columns: lhs - op - rhs
    if (!is.null(model$lhs) && !is.null(model$op)  &&
        !is.null(model$rhs) && !is.null(model$free)) {
      
      # ok, we have something that looks like a parameter table
      # FIXME: we need to check for redundant arguments
      # (but if cfa/sem was used, we can not trust the call)
      # redundant <- c("meanstructure", "int.ov.free", "int.lv.free",
      #        "fixed.x", "orthogonal", "std.lv", "parameterization",
      #        "auto.fix.first", "auto.fix.single", "auto.var",
      #        "auto.cov.lv.x", "auto.cov.y", "auto.th", "auto.delta")
      flat.model <- model
      
      # fix semTools issue here? for auxiliary() which does not use
      # block column yet
      if (!is.null(flat.model$block)) {
        nn <- length(flat.model$lhs)
        if (length(flat.model$block) != nn) {
          flat.model$block <- flat.model$group
        }
        if (any(is.na(flat.model$block))) {
          flat.model$block <- flat.model$group
        }
      } else if (!is.null(flat.model$group)) {
        flat.model$block <- flat.model$group
      }
      
    } else {
      bare.minimum <- c("lhs", "op", "rhs", "free")
      missing.idx <- is.na(match(bare.minimum, names(model)))
      missing.txt <- paste(bare.minimum[missing.idx], collapse = ", ")
      stop("lavaan ERROR: model is a list, but not a parameterTable?",
           "\n  lavaan  NOTE: ",
           "missing column(s) in parameter table: [", missing.txt, "]")
    }
  } else if (is.null(model)) {
    stop("lavaan ERROR: model is NULL!")
  }
  return(flat.model)
}

lav_lavaan_step01_ovnames_2 <- function(flat.model, ov.order, data, sample.cov, slotData) {
  # ov.order naar lowercase, check "data" of "model"
  # indien "data" flat.model proberen aanpassen met lav_partable_ov_from_data (warning indien mislukt)
  #           (datanames halen uit data of sample.cov of slotData)
  # TODO: naar lowercase en check kan veel vroeger gebeuren
  # TODO: in lav_partable_ov_from_data en waar gebruikt: stockeren in attributes ipv kunstgreep?
  
  # Ok, we got a flattened model; usually this a flat.model object, but it could
  # also be an already lavaanified parTable, or a bare-minimum list with
  # lhs/op/rhs/free elements
  
  # new in 0.6-14
  # if ov.order = "data", it would seem we need to intervene here;
  # we do this by 'injecting' dummy lhs da rhs statement in flat.model, to
  # 'trick' lav_partable_vnames() (which only sees the model!)
  ov.order <- tolower(ov.order)
  if (ov.order == "data") {
    flat.model.orig <- flat.model
    try(flat.model <- lav_partable_ov_from_data(flat.model, data = data,
                                                sample.cov = sample.cov,
                                                slotData   = slotData),
        silent = TRUE)
    if (inherits(flat.model, "try-error")) {
      warning("lavaan WARNING: ov.order = \"data\" setting failed; switching back to ov.order = \"model\"")
      flat.model <- flat.model.orig
    }
  } else if (ov.order != "model") {
    stop("lavaan ERROR: ov.order= argument should be \"model\" (default) or \"data\"")
  }
  return(flat.model)
}

lav_lavaan_step01_ovnames_3 <- function(flat.model, ov.names, ngroups) {
  # indien "group :" voorkomt in flat.model
  #   tmp.group.values : rechterzijden
  #   kopieer model zonder attributen en lavaanify dat -> tmp.lav
  #   ov.names, ov.names.y, ov.names.x, lv.names uit tmp.lav via lav_partable_vnames
  # anders
  #   als flat.model$group niet null
  #     group.values er uit halen via lav_partable_group_values (indien méér dan 1)
  #     ov.names, ov.names.y, ov.names.x, lv.names uit flat.model halen via lav_partable_vnames
  #   anders
  #     ov.names, ov.names.y, ov.names.x, lv.names uit flat.model halen via lav_partable_vnames
  #   
  # Check length(ov.names) > 1  TODO: als ov.names een list is, check voor alle elementen ???  
  #   
  # group blocks?
  # TODO: call lav_partable_vnames only ones and not for each type, when vnames is corrected
  #
  flat.model.2 <- NULL
  tmp.lav <- NULL
  group.values <- NULL
  
  if (any(flat.model$op == ":" & tolower(flat.model$lhs) == "group")) {
    # here, we only need to figure out:
    # - ngroups
    # - ov's per group
    #
    # - FIXME: we need a more efficient way, avoiding lavaanify/lav_partable_vnames
    #
    group.idx <- which(flat.model$op == ":" & tolower(flat.model$lhs) == "group")
    # replace by 'group' (in case we got 'Group'):
    flat.model$lhs[group.idx] <- "group"
    tmp.group.values <- unique(flat.model$rhs[group.idx])
    tmp.ngroups <- length(tmp.group.values)
    
    flat.model.2 <- flat.model
    attr(flat.model.2, "modifiers") <- NULL
    attr(flat.model.2, "constraints") <- NULL
    tmp.lav <- lavaanify(flat.model.2, ngroups = tmp.ngroups, warn = FALSE)
    ov.names <- ov.names.y <- ov.names.x <- lv.names <- vector("list",
                                                               length = tmp.ngroups)
    for (g in seq_len(tmp.ngroups)) {
      ov.names[[g]]   <- unique(unlist(lav_partable_vnames(tmp.lav,
                                                           type = "ov", group = tmp.group.values[g])))
      ov.names.y[[g]] <- unique(unlist(lav_partable_vnames(tmp.lav,
                                                           type = "ov.nox", group = tmp.group.values[g])))
      ov.names.x[[g]] <- unique(unlist(lav_partable_vnames(tmp.lav,
                                                           type = "ov.x", group = tmp.group.values[g])))
      lv.names[[g]] <- unique(unlist(lav_partable_vnames(tmp.lav,
                                                         type = "lv", group = tmp.group.values[g])))
    }
  } else if (!is.null(flat.model$group)) {
    # user-provided full partable with group column!
    ngroups <- lav_partable_ngroups(flat.model)
    if (ngroups > 1L) {
      group.values <- lav_partable_group_values(flat.model)
      ov.names <- ov.names.y <- ov.names.x <- lv.names <- vector("list",
                                                                 length = ngroups)
      for (g in seq_len(ngroups)) {
        # collapsed over levels (if any)
        ov.names[[g]]   <- unique(unlist(lav_partable_vnames(flat.model,
                                                             type = "ov", group = group.values[g])))
        ov.names.y[[g]] <- unique(unlist(lav_partable_vnames(flat.model,
                                                             type = "ov.nox", group = group.values[g])))
        ov.names.x[[g]] <- unique(unlist(lav_partable_vnames(flat.model,
                                                             type = "ov.x", group = group.values[g])))
        lv.names[[g]]   <- unique(unlist(lav_partable_vnames(flat.model,
                                                             type = "lv", group = group.values[g])))
      }
    } else {
      ov.names   <- lav_partable_vnames(flat.model, type = "ov")
      ov.names.y <- lav_partable_vnames(flat.model, type = "ov.nox")
      ov.names.x <- lav_partable_vnames(flat.model, type = "ov.x")
      lv.names   <- lav_partable_vnames(flat.model, type = "lv")
    }
  } else {
    # collapse over levels (if any)
    ov.names   <- unique(unlist(lav_partable_vnames(flat.model, type = "ov")))
    ov.names.y <- unique(unlist(lav_partable_vnames(flat.model, type = "ov.nox")))
    ov.names.x <- unique(unlist(lav_partable_vnames(flat.model, type = "ov.x")))
    lv.names   <- unique(unlist(lav_partable_vnames(flat.model, type = "lv")))
  }
  
  # sanity check (new in 0.6-8): do we have any ov.names?
  # detect early
  if (length(ov.names) == 0L) {
    stop("lavaan ERROR: ov.names is empty: model does not refer to any observed variables; check your syntax.")
  }
  return(list(flat.model = flat.model,
              ov.names = ov.names,
              ov.names.x = ov.names.x,
              ov.names.y = ov.names.y,
              lv.names = lv.names,
              group.values = group.values,
              ngroups = ngroups))
}

lav_lavaan_step01_ovnames_4 <- function(lv.names, data, sample.cov, dotdotdot, slotOptions) {
  # latente variabelen kunnen niet in data zitten --> *** error *** (behalve indien expliciet gevraagd)
  # latente interacties niet toegelaten ---> *** error ***
  
  # sanity check: ov.names.x should NOT appear in ov.names.y
  # this may happen if 'x' is exogenous in one block, but not in another...
  #endo.idx <- which(ov.names.x %in% ov.names.y)
  #if (length(endo.idx) > 0L) {
  #    # remove from x! (new in 0.6-8)
  #    ov.names.x <- ov.names.x[-endo.idx]
  #}
  
  
  # handle for lv.names that are also observed variables (new in 0.6-6)
  lv.lv.names <- unique(unlist(lv.names))
  if (length(lv.lv.names) > 0L) {
    
    # check for lv.names in data/cov
    if (!is.null(data)) {
      bad.idx <- which(lv.lv.names %in% names(data))
    } else if (!is.null(sample.cov)) {
      bad.idx <- which(lv.lv.names %in% rownames(sample.cov))
    } else {
      bad.idx <- integer(0L)
    }
    
    # if found, hard stop
    if (length(bad.idx) > 0L) {
      if (!is.null(dotdotdot$check.lv.names) &&
          !dotdotdot$check.lv.names) {
        # ignore it, user switched this check off -- new in 0.6-7
      } else {
        stop("lavaan ERROR: some latent variable names collide ",
             "with observed\n\t\tvariable names: ",
             paste(lv.lv.names[bad.idx], collapse = " "))
      }
      
      # rename latent variables (by adding 'lat')
      #flat.model.idx <- which(flat.model$op == "=~" &
      #                  flat.model$lhs %in% lv.names[bad.idx])
      #flat.model$lhs[flat.model.idx] <- paste(flat.model$lhs[flat.model.idx], "lat", sep = "")
      
      # add names to ov.names
      #ov.names <- c(ov.names, lv.names[bad.idx])
      # what about ov.names.y and ov.names.x?
    }
  }
  
  # sanity check: we do not support latent interaction yet (using the :)
  lv.int.idx <- which(grepl(":", lv.lv.names))
  if (length(lv.int.idx) > 0L) {
    if (!is.null(dotdotdot$check.lv.interaction) &&
        !dotdotdot$check.lv.interaction) {
      # ignore, user (or sam) switched this check off - new in 0.6-16
    } else if (!is.null(slotOptions) && !slotOptions$check.lv.interaction) {
      # ignore
    } else {
      txt <- c("Interaction terms involving latent variables (",
               lv.lv.names[lv.int.idx[1]], ") are not supported.",
               " You may consider creating product indicators to define ",
               "the latent interaction term. See the indProd() function ",
               "in the semTools package.")
      stop(lav_txt2message(txt, header = "lavaan ERROR:"))
    }
  }
}

lav_lavaan_step01_ovnames_5 <- function(data, cluster, flat.model, group.values, ngroups) {
  # indien "level :" voorkomt in flat.model 
  #   cluster moet vermeld zijn indien data opgegeven, zoniet *** error ***
  #   bepaal tmp.group.values en tmp.level.values
  #   er moeten minstens 2 levels gedefinieerd zijn, zoniet *** error ***
  #   kopieer flat.model zonder attributen en lavaanify dat -> tmp.lav
  #   check minstens 2 levels voor tmp.lav
  #   bepaal ov.names.l per group en per level (via lav_partable_vnames op tmp.lav)
  # anders
  #   als lav_partable_nlevels(flat.model) > 0 
  #     cluster met vermeld zijn indien data opgegeven, zoniet *** error ***
  #     bepaal ov.names.l per group en per level via lavNames op flat.model
  #   anders
  #     geen levels (ov.names.l = list())
  #   
  # handle ov.names.l
  if (any(flat.model$op == ":" & tolower(flat.model$lhs) == "level")) {
    
    # check for cluster argument
    if (!is.null(data) && is.null(cluster)) {
      stop("lavaan ERROR: cluster argument is missing.")
    }
    
    # here, we only need to figure out:
    # - nlevels
    # - ov's per level
    # - FIXME: we need a more efficient way, avoiding lavaanify/vnames
    
    group.idx <- which(flat.model$op == ":" & flat.model$lhs == "group")
    tmp.group.values <- unique(flat.model$rhs[group.idx])
    tmp.ngroups <- max(c(length(tmp.group.values), 1))
    
    level.idx <- which(flat.model$op == ":" & tolower(flat.model$lhs) == "level")
    # replace by "level" (in case we got 'Level')
    flat.model$lhs[level.idx] <- "level"
    tmp.level.values <- unique(flat.model$rhs[level.idx])
    tmp.nlevels <- length(tmp.level.values)
    
    # we need at least 2 levels (for now)
    if (tmp.nlevels < 2L) {
      stop("lavaan ERROR: when data is clustered, you must specify a model\n",
           "  for each level in the model syntax (for now); see example(Demo.twolevel)")
    }
    
    flat.model.2 <- flat.model
    attr(flat.model.2, "modifiers") <- NULL
    attr(flat.model.2, "constraints") <- NULL
    tmp.lav <- lavaanify(flat.model.2, ngroups = tmp.ngroups, warn = FALSE)
    # check for empty levels
    if (max(tmp.lav$level) < 2L) {
      stop("lavaan ERROR: at least one level has no model syntax;",
           "you must specify a model for each level in the model syntax (for now); see example(Demo.twolevel)")
    }
    ov.names.l <- vector("list", length = tmp.ngroups) # per group
    
    for (g in seq_len(tmp.ngroups)) {
      ov.names.l[[g]] <- vector("list", length = tmp.nlevels)
      for (l in seq_len(tmp.nlevels)) {
        if (tmp.ngroups > 1L) {
          ov.names.l[[g]][[l]] <-
            unique(unlist(lav_partable_vnames(tmp.lav,
                                              type = "ov",
                                              group = tmp.group.values[g],
                                              level = tmp.level.values[l])))
        } else {
          ov.names.l[[g]][[l]] <-
            unique(unlist(lav_partable_vnames(tmp.lav,
                                              type = "ov",
                                              level = tmp.level.values[l])))
        }
      } # levels
    } # groups
  } else {
    # perhaps model is already a parameter table
    nlevels <- lav_partable_nlevels(flat.model)
    if (nlevels > 1L) {
      
      # check for cluster argument (only if we have data)
      if (!is.null(data) && is.null(cluster)) {
        stop("lavaan ERROR: cluster argument is missing.")
      }
      
      ngroups <- lav_partable_ngroups(flat.model)
      group.values <- lav_partable_group_values(flat.model)
      ov.names.l <- vector("list", length = ngroups)
      for (g in 1:ngroups) {
        # note: lavNames() will return a list if any level:
        ov.names.l[[g]] <- lavNames(flat.model, "ov", group = group.values[g])
      }
    } else {
      # no level: in model syntax
      ov.names.l <- list()
    }
  }
  return(list(
    flat.model = flat.model,
    ov.names.l = ov.names.l
  ))
}

lav_lavaan_step01_ovnames_6 <- function(ordered, flat.model, data) {
  # interpretatie en check ordered parameter, eventueel aanpassen
  
  # sanity check ordered argument (just in case, add lhs variables names)
  if (!is.null(ordered)) { # new in 0.6-4
    if (is.logical(ordered) && ordered) {  # ordered = TRUE
      # assume the user means: ordered = names(Data)
      ordered <- lavNames(flat.model, "ov.nox") # new in 0.6-6: changed from ov
    } else if (is.logical(ordered) && !ordered) {
      ordered <- character(0L)
    } else if (!is.character(ordered)) {
      stop("lavaan ERROR: ordered argument must be a character vector")
    } else if (length(ordered) == 1L && nchar(ordered) == 0L) {
      ordered <- character(0L)
    } else {
      # check if all names in "ordered" occur in the dataset?
      if (!is.null(data)) {
        if (inherits(data, "data.frame")) {
          NAMES <- names(data)
        } else if (inherits(data, "matrix")) {
          NAMES <- colnames(data)
        }
        
        missing.idx <- which(!ordered %in% NAMES)
        if (length(missing.idx) > 0L) { # FIXme: warn = FALSE has no eff
          warning("lavaan WARNING: ordered variable(s): ",
                  paste(ordered[missing.idx], collapse = " "),
                  "\n  could not be found in the data and will be ignored")
        }
      }
    }
  }
  # add the variable names that were treated as ordinal
  # in the model syntax
  ordered <- unique(c(ordered, lavNames(flat.model, "ov.ord")))
  return(ordered)
}
