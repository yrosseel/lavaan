# user visible function to add 'matrix' entries in the parameter table
lavMatrixRepresentation <- function(partable, representation = "LISREL", # nolint start
                                    allow.composites = TRUE, # new in 0.6-20
                                    add.attributes = FALSE,
                                    as.data.frame. = TRUE) {             # nolint end
  # check parameter table
  partable <- lav_partable_complete(partable)

  # get model matrices
  if (representation == "LISREL") {
    rep_1 <- lav_lisrel(partable, target = NULL, extra = add.attributes,
                      allow.composites = allow.composites)
  } else if (representation == "RAM") {
    rep_1 <- lav_ram(partable, target = NULL, extra = add.attributes)
  } else {
    lav_msg_stop(gettext("representation must either \"LISREL\" or \"RAM\"."))
  }

  partable$mat <- rep_1$mat
  partable$row <- rep_1$row
  partable$col <- rep_1$col

  if (as.data.frame.) {
    partable <- as.data.frame(partable, stringsAsFactors = FALSE)
    class(partable) <- c("lavaan.data.frame", "data.frame")
  }

  if (add.attributes) {
    if (representation == "LISREL") {
      attr(partable, "ov.dummy.names.nox") <- attr(rep_1, "ov.dummy.names.nox")
      attr(partable, "ov.dummy.names.x") <- attr(rep_1, "ov.dummy.names.x")
    } else if (representation == "RAM") {
      attr(partable, "ov.idx") <- attr(rep_1, "ov.idx")
    }
    attr(partable, "mmNames") <- attr(rep_1, "mmNames")
    attr(partable, "mmNumber") <- attr(rep_1, "mmNumber")
    attr(partable, "mmRows") <- attr(rep_1, "mmRows")
    attr(partable, "mmCols") <- attr(rep_1, "mmCols")
    attr(partable, "mmDimNames") <- attr(rep_1, "mmDimNames")
    attr(partable, "mmSymmetric") <- attr(rep_1, "mmSymmetric")
  }

  partable
}
