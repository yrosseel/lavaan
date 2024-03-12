# store pta in attributes of partable
lav_partable_set_cache <- function(partable, pta = NULL, force = FALSE) {
  if (is.null(pta)) {
    if (force) attr(partable, "vnames") <- NULL
    pta <- lav_partable_attributes(partable)
  }

  for (n in names(pta)) {
    attr(partable, n) <- pta[[n]]
  }

  partable
}

lav_partable_remove_cache <- function(partable) {
  attributelist <- names(attributes(partable))

  for (n in attributelist) {
    if (n != "ovda" && n != "names") attr(partable, n) <- NULL
  }

  partable
}
