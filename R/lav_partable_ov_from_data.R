# handle ov.order = "data" by adding attribute "ovda" to FLAT
lav_partable_ov_from_data <- function(flat = NULL, # nolint
                                      data = NULL,
                                      sample_cov = NULL,
                                      slot_data = NULL) { # nolint
  # current model-based ov.names
  ov_names <- lav_partable_vnames(flat, type = "ov")

  # get data-based ov.names
  data_names <- NULL
  if (!is.null(data)) {
    data_names <- names(data)
  } else if (!is.null(sample_cov)) {
    # multiple group/blocks?
    if (is.list(sample_cov)) {
      data_names <- unique(unlist(lapply(sample_cov, colnames)))
      if (is.null(data_names)) {
        # try again with rows
        data_names <- unique(unlist(lapply(sample_cov, rownames)))
      }
    } else {
      data_names <- colnames(sample_cov)
      if (is.null(data_names)) {
        # try again with rows
        data_names <- rownames(sample_cov)
      }
    }
  } else if (!is.null(slot_data)) {
    data_names <- unique(unlist(slot_data@ov.names))
  }

  if (is.null(data_names) || length(data_names) == 0L) {
    lav_msg_stop(gettext("could not find variable names in data/sample.cov"))
  }

  # extract needed ov.names in the same order as the data
  ov_names_data <- data_names[data_names %in% ov_names]

  # check if we have all of them
  if (length(ov_names_data) != length(ov_names)) {
    idx_missing <- which(!(ov_names %in% ov_names_data))
    lav_msg_stop(gettextf(
      "some (observed) variables specified in the model are not found
      in the data: %s",
      lav_msg_view(ov_names[idx_missing], "none")))
  }

  # check if the order is the same
  if (!identical(ov_names, ov_names_data)) {
    attr(flat, "ovda") <- ov_names_data # nolint
  }
  flat
}
