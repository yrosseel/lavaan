# handle ov.order = "data" by adding attribute "ovda" to FLAT
lav_partable_ov_from_data <- function(FLAT = NULL, # nolint
                                      data = NULL,
                                      sample.cov = NULL,
                                      slotData = NULL) { # nolint
  # current model-based ov.names
  ov.names <- lav_partable_vnames(FLAT, type = "ov")

  # get data-based ov.names
  data.names <- NULL
  if (!is.null(data)) {
    data.names <- names(data)
  } else if (!is.null(sample.cov)) {
    # multiple group/blocks?
    if (is.list(sample.cov)) {
      data.names <- unique(unlist(lapply(sample.cov, colnames)))
      if (is.null(data.names)) {
        # try again with rows
        data.names <- unique(unlist(lapply(sample.cov, rownames)))
      }
    } else {
      data.names <- colnames(sample.cov)
      if (is.null(data.names)) {
        # try again with rows
        data.names <- rownames(sample.cov)
      }
    }
  } else if (!is.null(slotData)) {
    data.names <- unique(unlist(slotData@ov.names))
  }

  if (is.null(data.names) || length(data.names) == 0L) {
    lav_msg_stop(gettext("could not find variable names in data/sample.cov"))
  }

  # extract needed ov.names in the same order as the data
  ov.names.data <- data.names[data.names %in% ov.names]

  # check if we have all of them
  if (length(ov.names.data) != length(ov.names)) {
    idx.missing <- which(!(ov.names %in% ov.names.data))
    lav_msg_stop(gettextf(
      "some (observed) variables specified in the model are not found
      in the data: %s",
      lav_msg_view(ov.names[idx.missing], "none")))
  }

  # check if the order is the same
  if (!identical(ov.names, ov.names.data)) {
    attr(FLAT, "ovda") <- ov.names.data # nolint
  }
  return(FLAT)
}
