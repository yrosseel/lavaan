# this function is written by Myrsini Katsikatsou

############################## pairwiseTables FUNCTION ########################
# This function can be public. It gets as an input a raw data set of ordinal
# variables and it returns a list of all pairwise frequency tables.
#
# The input arguments of the function:
# data : matrix or data frame containing the data. The rows correspond to
#       different observations and the columns to different observed categorical
#       (ordinal or nominal) variables. No continuous variables or covariates
#       should be contained in data. If the variables contained in the data are
#       distinguished into indicators of exogenous latent variables (lv) and
#       indicators of endogenous latent variables, those for exogenous lv should
#       be presented first (in the first columns of data) followed by the
#       indicators for endogenous lv.
# var.levels: NULL or vector or list, specifies the levels (response categories)
#            for each categorical variable contained in data.
#            If NULL, the levels encoutered in data are used. If a response
#            category is not observed in the data, then var.levels should be
#            defined.
#            If vector, that implies that all variables have the same levels as
#            given in the vector.
#            If list, the components of the list are vectors, as many as the
#            number of variables in data. Each vector gives the levels of
#            the corresponding categorical variable in data.
# no.x : NULL or integer, gives the number of indicators for exogenous lv.
#        The default value is NULL indicating that data contains only
#        indicators of exogenous latent variables.
# perc : TRUE/FALSE. If FALSE the observed frequencies are reported, otherwise
#        the observed percentages are given.
# na.exclude : TRUE/FALSE. If TRUE, listwise deletion is applied to data.
#              Otherwise, cases with missing values are preserved and and an
#              extra level with label NA is included in the tables.

# The output of the function:
# It is a list of three components: $pairTables, $VarLevels and $Ncases_del.
# pairTables : a list of so many tables as the number of variable pairs formed
#              by data. If there are indicators of both exogenous and endogenous
#              variables, then first all the matrices referring to pairs of
#              indicators of exogenous lv are reported, followed by all the
#              matrices referring to pairs of indicators of endogenous lv, which
#              in turn folowed by all the matrices of pairs: one indicator of an
#              exogenous - one indicator of an endogenous lv.
# VarLevels : a list of as many vectors as the number of variables in the data.
#             Each vector gives the levels/ response categories of each variable
# Ncases_del : An integer reporting the number of cases deleted by data because
#              of missing values (listwise deletion) when na.exclude=TRUE.



pairwiseTables <- function(data, var.levels = NULL, no.x = NULL,
                           perc = FALSE, na.exclude = TRUE) {
  # data in right format?
  if ((!is.matrix(data)) & (!is.data.frame(data))) {
    lav_msg_stop(gettext("data is neither a matrix nor a data.frame"))
  }

  # at least two variables
  no.var <- dim(data)[2]
  if (no.var < 2) {
    lav_msg_stop(gettext("there are less than 2 variables"))
  }

  # no.x < no.var  ?
  if (no.x > no.var) {
    lav_msg_stop(gettext(
      "number of indicators for exogenous latent variables is larger than
      the total number of variables in data"))
  }


  # if data as matrix, transforma as data.frame
  if (is.matrix(data)) {
    data <- as.data.frame(data)
  }

  # listwise deletion
  if (na.exclude) {
    old.data <- data
    data <- na.omit(data)
  }

  # all columns of data.frame should be of class factor so that function levels
  # can be applied
  if (!all(sapply(data, class) == "factor")) {
    if (nrow(data) > 1) {
      data <- data.frame(sapply(data, factor))
    } else {
      data <- apply(data, 2, factor)
      data <- as.data.frame(matrix(data, nrow = 1))
    }
  }

  # the levels observed for each variable, obs.levels is a list
  obs.levels <- lapply(data, levels)

  # number of variables in data same as number of vectors in var.levels
  if (is.list(var.levels) && no.var != length(var.levels)) {
    lav_msg_stop(gettext(
      "the length of var.levels does not match the number of variables of
      the given data set"))
  }

  # create var.levels if a list is not given
  old.var.levels <- var.levels
  if (!is.list(old.var.levels)) {
    if (is.null(old.var.levels)) {
      var.levels <- obs.levels
    } else {
      var.levels <- vector("list", no.var)
      var.levels <- lapply(var.levels, function(x) {
        x <- old.var.levels
      })
    }
  }
  names(var.levels) <- names(data)

  # also check that obs.levels exist in the object var.levels given by the user, i.e. old.var.levels
  if (is.list(old.var.levels)) {
    for (i in 1:no.var) {
      if (!all(obs.levels[[i]] %in% old.var.levels[[i]])) {
        lav_msg_stop(gettext(
          "levels observed in data are not mentioned in var.levels"))
      }
    }
  } else if (is.vector(old.var.levels)) {
    if (!all(apply(na.omit(data), 2, function(x) {
      x %in% old.var.levels
    }))) {
      lav_msg_stop(gettext("levels observed in data are not mentioned
                           in var.levels"))
    }
  }

  no.given.levels <- sapply(var.levels, length)

  # assign the right levels for each variable as given in object var.levels if it is not the case
  # it is not the case when the observed levels are a subgroup of the var.levels given
  if (!is.null(old.var.levels)) {
    no.obs.levels <- sapply(obs.levels, length)
    if (!all(no.obs.levels == no.given.levels)) {
      index <- c(1:no.var)[no.obs.levels != no.given.levels]
      for (i in index) {
        data[, i] <- factor(data[, i], levels = var.levels[[i]])
      }
    }
  }

  # compute the bivariate frequency tables
  # Split first into two cases: a) only indicators of exogenous latent variables
  # b) otherwise
  if (is.null(no.x) || no.x == no.var) {
    pairs.index <- utils::combn(no.var, 2)
    no.pairs <- dim(pairs.index)[2]
    res <- vector("list", no.pairs)
    for (i in 1:no.pairs) {
      res[[i]] <- table(data[, pairs.index[, i]], useNA = "ifany")
    }
  } else {
    no.y <- no.var - no.x
    pairs.xixj.index <- utils::combn(no.x, 2) # row 1 gives i index, row 2 j index, j runs faster than i
    pairs.yiyj.index <- utils::combn(no.y, 2)
    pairs.xiyj.index <- expand.grid(1:no.y, 1:no.x)
    pairs.xiyj.index <- rbind(pairs.xiyj.index[, 2], pairs.xiyj.index[, 1]) # row 1 gives i index, row 2 j index, j runs faster than i

    no.pairs.xixj <- dim(pairs.xixj.index)[2]
    no.pairs.yiyj <- dim(pairs.yiyj.index)[2]
    no.pairs.xiyj <- dim(pairs.xiyj.index)[2]
    no.all.pairs <- no.pairs.xixj + no.pairs.yiyj + no.pairs.xiyj

    data.x <- data[, 1:no.x]
    data.y <- data[, (no.x + 1):no.var]

    res <- vector("list", no.all.pairs)
    for (i in 1:no.pairs.xixj) {
      res[[i]] <- table(data.x[, pairs.xixj.index[, i]], useNA = "ifany")
    }

    j <- 0
    for (i in (no.pairs.xixj + 1):(no.pairs.xixj + no.pairs.yiyj)) {
      j <- j + 1
      res[[i]] <- table(data.y[, pairs.yiyj.index[, j]], useNA = "ifany")
    }

    j <- 0
    for (i in (no.pairs.xixj + no.pairs.yiyj + 1):no.all.pairs) {
      j <- j + 1
      res[[i]] <- table(
        cbind(
          data.x[, pairs.xiyj.index[1, j], drop = FALSE],
          data.y[, pairs.xiyj.index[2, j], drop = FALSE]
        ),
        useNA = "ifany"
      )
    }
  }

  # if percentages are asked
  if (perc) {
    Nobs <- dim(data)[1]
    res <- lapply(res, function(x) {
      x / Nobs
    })
  }

  # Ncases_del = the number of cases deleted because they had missing values
  if (na.exclude) {
    Ncases_deleted <- dim(old.data)[1] - dim(data)[1]
  } else {
    Ncases_deleted <- 0
  }

  list(pairTables = res, VarLevels = var.levels, Ncases_del = Ncases_deleted)
}
