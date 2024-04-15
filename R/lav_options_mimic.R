# LDW 11/4/24 : overwrite defaults depending on mimic in separate function
#
lav_options_mimic <- function(opt) {
  mlr.test <- "yuan.bentler.mplus" # for now
  if (opt$mimic == "lavaan") {
    if (is.character(opt$conditional.x)) { # = "default"
      if (lav_options_estimatorgroup(opt$estimator) == "ML") {
        opt$conditional.x <- FALSE
      }
    }
    if (opt$fixed.x == "default") {
      if (any(lav_options_estimatorgroup(opt$estimator) == c("MML", "ML")) &&
          is.character(opt$start) && opt$start != "simple") { # new in 0.6-12
        opt$fixed.x <- TRUE
      }
    }
    if (is.character(opt$zero.keep.margins)) { # = "default"
      opt$zero.keep.margins <- TRUE
    }
  } else if (opt$mimic == "Mplus") {
    if (length(opt$group.equal) == 0L || all(nchar(opt$group.equal) == 0L)) {
      if (opt$.categorical) {
        opt$group.equal <- c("loadings", "thresholds")
      } else {
        if (is.logical(opt$meanstructure) && !opt$meanstructure) {
          opt$group.equal <- "loadings"
        } else {
          opt$group.equal <- c("loadings", "intercepts")
        }
      }
    }
    if (opt$missing == "default") {
      if (!opt$.categorical && any(opt$estimator == c("ml", "mlr")))
      {
        # since version 5?
        opt$missing <- "ml"
        # check later if this is ok
      }
    }
    if (opt$estimator != "pml") {
      if (opt$meanstructure == "default") opt$meanstructure <- TRUE
    }
    if (opt$estimator == "mlr") mlr.test <- "yuan.bentler.mplus"
    if (any(lav_options_estimatorgroup(opt$estimator) ==
            c("ML", "REML", "NTRLS", "catML"))) {
      }
    if (is.character(opt$conditional.x)) { # = "default"
      if (lav_options_estimatorgroup(opt$estimator) == "ML") {
        opt$conditional.x <- FALSE
      }
    }
    if (opt$fixed.x == "default") {
      if (any(lav_options_estimatorgroup(opt$estimator) == c("MML", "ML")) &&
          is.character(opt$start) && opt$start != "simple") { # new in 0.6-12
        opt$fixed.x <- TRUE
      }
    }
    if (is.character(opt$zero.keep.margins)) { # = "default"
      opt$zero.keep.margins <- TRUE
    }
    opt$baseline.conditional.x.free.slopes <- FALSE
  } else if (opt$mimic == "EQS") {
    if (opt$estimator == "mlr") mlr.test <- "yuan.bentler"
    if (any(lav_options_estimatorgroup(opt$estimator) ==
            c("ML", "REML", "NTRLS", "catML"))) {
      if (opt$likelihood == "default") opt$likelihood <- "wishart"
    }
  } else if (opt$mimic == "LISREL") {
    if (any(lav_options_estimatorgroup(opt$estimator) ==
            c("ML", "REML", "NTRLS", "catML"))) {
      if (opt$likelihood == "default") opt$likelihood <- "wishart"
    }
  }
  if (opt$estimator == "mlr") {
    if (opt$test[1] == "default") {
      opt$test <- mlr.test
    } else {
      opt$test <- union(mlr.test, opt$test)
    }
  }
  opt
}