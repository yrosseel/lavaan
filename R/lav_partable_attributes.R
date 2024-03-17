# return 'attributes' of a lavaan partable -- generate a new set if necessary
lav_partable_attributes <- function(partable, pta = NULL) {
  if (is.null(pta)) {
    # attached to partable?
    pta <- attributes(partable)
    if (!is.null(pta$vnames) && !is.null(pta$nvar)) {
      # looks like a pta
      pta$ovda <- NULL
      return(pta)
    } else {
      pta <- list()
    }
  }

  # vnames
  pta$vnames <- lav_partable_vnames(partable, type = "*")

  # vidx
  tmp.ov <- pta$vnames$ov
  tmp.lv <- pta$vnames$lv
  nblocks <- length(pta$vnames$ov)
  pta$vidx <- lapply(names(pta$vnames), function(v) {
    lapply(seq_len(nblocks), function(b) {
      if (v == "lv.marker") {
        match(pta$vnames[[v]][[b]], tmp.ov[[b]])
      } else if (grepl("lv", v)) {
        match(pta$vnames[[v]][[b]], tmp.lv[[b]])
      } else if (grepl("th", v)) {
        # thresholds have '|t' pattern
        tmp.th <- sapply(strsplit(pta$vnames[[v]][[b]],
          "|t",
          fixed = TRUE
        ), "[[", 1L)
        match(tmp.th, tmp.ov[[b]])
      } else if (grepl("eqs", v)) {
        # mixture of tmp.ov/tmp.lv
        integer(0L)
      } else {
        match(pta$vnames[[v]][[b]], tmp.ov[[b]])
      }
    })
  })
  names(pta$vidx) <- names(pta$vnames)

  # meanstructure
  pta$meanstructure <- any(partable$op == "~1")

  # nblocks
  pta$nblocks <- nblocks

  # ngroups
  pta$ngroups <- lav_partable_ngroups(partable)

  # nlevels
  pta$nlevels <- lav_partable_nlevels(partable)

  # nvar
  pta$nvar <- lapply(pta$vnames$ov, length)

  # nfac
  pta$nfac <- lapply(pta$vnames$lv, length)

  # nfac.nonnormal - for numerical integration
  pta$nfac.nonnormal <- lapply(pta$vnames$lv.nonnormal, length)

  # th.idx (new in 0.6-1)
  pta$th.idx <- lapply(seq_len(pta$nblocks), function(b) {
    out <- numeric(length(pta$vnames$th.mean[[b]]))
    idx <- (pta$vnames$th.mean[[b]] %in%
      pta$vnames$th[[b]])
    out[idx] <- pta$vidx$th[[b]]
    out
  })

  pta
}
