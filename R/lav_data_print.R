# print object from lavData class
#

setMethod(
  "show", "lavData",
  function(object) {
    # print 'lavData' object
    res <- lav_data_summary_short(object)
    lav_data_print_short(res, nd = 3L)
  }
)

# create summary information for @lavdata slot
lav_data_summary_short <- function(object) {
  # which object?
  if (inherits(object, "lavaan")) {
    lavdata <- object@Data
  } else if (inherits(object, "lavData")) {
    lavdata <- object
  } else {
    lav_msg_stop(gettext("object must be lavaan or lavData object"))
  }

  # two or three columns (depends on nobs/norig)
  threecolumn <- FALSE
  for (g in 1:lavdata@ngroups) {
    if (lavdata@nobs[[g]] != lavdata@norig[[g]]) {
      threecolumn <- TRUE
      break
    }
  }

  # clustered data?
  clustered <- FALSE
  if (.hasSlot(lavdata, "cluster") && # in case we have an old obj
    length(lavdata@cluster) > 0L) {
    clustered <- TRUE
  }

  # multilevel data?
  multilevel <- FALSE
  if (.hasSlot(lavdata, "nlevels") && # in case we have an old obj
    lavdata@nlevels > 1L) {
    multilevel <- TRUE
  }

  # extract summary information
  datasummary <- list(
    ngroups = lavdata@ngroups,
    nobs = unlist(lavdata@nobs)
  )

  # norig?
  if (threecolumn) {
    datasummary$norig <- unlist(lavdata@norig)
  }

  # multiple groups?
  if (lavdata@ngroups > 1L) {
    datasummary$group.label <- lavdata@group.label
  }

  # sampling weights?
  if ((.hasSlot(lavdata, "weights")) && # in case we have an old object
    (!is.null(lavdata@weights[[1L]]))) {
    datasummary$sampling.weights <- lavdata@sampling.weights
  }

  # clustered/multilevel data?
  if (clustered) {
    if (multilevel) {
      datasummary$nlevels <- lavdata@nlevels
    }
    datasummary$cluster <- lavdata@cluster

    if (lavdata@ngroups == 1L) {
      datasummary$nclusters <- unlist(lavdata@Lp[[1]]$nclusters)
    } else {
      tmp <- vector("list", length = lavdata@ngroups)
      for (g in seq_len(lavdata@ngroups)) {
        tmp[[g]] <- unlist(lavdata@Lp[[g]]$nclusters)
      }
      datasummary$nclusters <- tmp
    }
  }

  # missing data?
  if (!is.null(lavdata@Mp[[1L]])) {
    datasummary$npatterns <- sapply(lavdata@Mp, "[[", "npatterns")
    if (multilevel && !is.null(lavdata@Mp[[1L]]$Zp)) {
      datasummary$npatterns2 <- sapply(lapply(
        lavdata@Mp,
        "[[", "Zp"
      ), "[[", "npatterns")
    }
  }

  datasummary
}

lav_data_print_short <- function(object, nd = 3L) {
  # object should data summary
  if (inherits(object, "lavaan")) {
    object <- lav_data_summary_short(object)
  }
  datasummary <- object

  num.format <- paste("%", max(8L, nd + 5L), ".", nd, "f", sep = "")

  # threecolumn
  threecolumn <- !is.null(datasummary$norig)

  # multilevel?
  multilevel <- !is.null(datasummary$nlevels)

  # clustered?
  clustered <- !is.null(datasummary$cluster) && is.null(datasummary$nlevels)

  # header? no, for historical reasons only
  # cat("Data information:\n\n")

  c1 <- c2 <- c3 <- character(0L)

  # number of observations
  if (datasummary$ngroups == 1L) {
    if (threecolumn) {
      c1 <- c(c1, "")
      c2 <- c(c2, "Used")
      c3 <- c(c3, "Total")
    }
    c1 <- c(c1, "Number of observations")
    c2 <- c(c2, datasummary$nobs)
    c3 <- c(c3, ifelse(threecolumn, datasummary$norig, ""))
  } else {
    c1 <- c(c1, "Number of observations per group:")
    if (threecolumn) {
      c2 <- c(c2, "Used")
      c3 <- c(c3, "Total")
    } else {
      c2 <- c(c2, "")
      c3 <- c(c3, "")
    }
    for (g in 1:datasummary$ngroups) {
      c1 <- c(c1, sprintf("  %-40s", datasummary$group.label[g]))
      c2 <- c(c2, datasummary$nobs[g])
      c3 <- c(c3, ifelse(threecolumn, datasummary$norig[g], ""))
    } # g
  }

  # number of clusters
  if (datasummary$ngroups == 1L) {
    if (multilevel) {
      for (l in 2:datasummary$nlevels) {
        c1 <- c(
          c1,
          paste("Number of clusters [",
            datasummary$cluster[l - 1], "]",
            sep = ""
          )
        )
        c2 <- c(c2, datasummary$nclusters[l])
        c3 <- c(c3, "")
      }
    } else if (clustered) {
      c1 <- c(c1, paste("Number of clusters [", datasummary$cluster, "]",
        sep = ""
      ))
      c2 <- c(c2, datasummary$nclusters[2])
      c3 <- c(c3, "")
    }
  } else {
    if (multilevel) {
      for (l in 2:datasummary$nlevels) {
        c1 <- c(
          c1,
          paste("Number of clusters [", datasummary$cluster[l - 1], "]:",
            sep = ""
          )
        )
        c2 <- c(c2, "")
        c3 <- c(c3, "")
        for (g in 1:datasummary$ngroups) {
          c1 <- c(c1, sprintf("  %-40s", datasummary$group.label[g]))
          c2 <- c(c2, datasummary$nclusters[[g]][l])
          c3 <- c(c3, "")
        }
      }
    } else if (clustered) {
      c1 <- c(
        c1,
        paste("Number of clusters [", datasummary$cluster, "]:", sep = "")
      )
      c2 <- c(c2, "")
      c3 <- c(c3, "")
      for (g in 1:datasummary$ngroups) {
        c1 <- c(c1, sprintf("  %-40s", datasummary$group.label[g]))
        c2 <- c(c2, datasummary$nclusters[[g]][2])
        c3 <- c(c3, "")
      }
    }
  }

  # missing patterns?
  if (!is.null(datasummary$npatterns)) {
    if (datasummary$ngroups == 1L) {
      if (multilevel) {
        c1 <- c(c1, "Number of missing patterns -- level 1")
        c2 <- c(c2, datasummary$npatterns)
        c3 <- c(c3, "")
        if (!is.null(datasummary$npatterns2)) {
          c1 <- c(c1, "Number of missing patterns -- level 2")
          c2 <- c(c2, datasummary$npatterns2)
          c3 <- c(c3, "")
        }
      } else {
        c1 <- c(c1, "Number of missing patterns")
        c2 <- c(c2, datasummary$npatterns)
        c3 <- c(c3, "")
      }
    } else {
      if (multilevel) {
        c1 <- c(c1, "Number of missing patterns per group:")
        c2 <- c(c2, "")
        c3 <- c(c3, "")
        for (g in 1:datasummary$ngroups) {
          c1 <- c(
            c1,
            paste(sprintf(
              "  %-40s",
              datasummary$group.label[g]
            ), "-- level 1")
          )
          c2 <- c(c2, datasummary$npatterns[g])
          c3 <- c(c3, "")
          if (!is.null(datasummary$npatterns2)) {
            c1 <- c(
              c1,
              paste(sprintf(
                "  %-40s",
                datasummary$group.label[g]
              ), "-- level 2")
            )
            c2 <- c(c2, datasummary$npatterns2[g])
            c3 <- c(c3, "")
          }
        }
      } else {
        c1 <- c(c1, "Number of missing patterns per group:")
        c2 <- c(c2, "")
        c3 <- c(c3, "")
        for (g in 1:datasummary$ngroups) {
          c1 <- c(c1, sprintf("  %-40s", datasummary$group.label[g]))
          c2 <- c(c2, datasummary$npatterns[g])
          c3 <- c(c3, "")
        }
      }
    }
  }

  # sampling weights?
  if (!is.null(datasummary$sampling.weights)) {
    c1 <- c(c1, "Sampling weights variable")
    c2 <- c(c2, datasummary$sampling.weights)
    c3 <- c(c3, "")
  }

  # format c1/c2
  c1 <- format(c1, width = 43L)
  c2 <- format(c2, width = 8L + max(0, (nd - 3L)) * 4L, justify = "right")
  c3 <- format(c3, width = 8L + nd, justify = "right")

  # create character matrix
  if (threecolumn) {
    M <- cbind(c1, c2, c3, deparse.level = 0)
  } else {
    M <- cbind(c1, c2, deparse.level = 0)
  }
  colnames(M) <- rep("", ncol(M))
  rownames(M) <- rep(" ", nrow(M))

  # print
  write.table(M, row.names = TRUE, col.names = FALSE, quote = FALSE)

  invisible(M)
}
