load_lavaan_data <- function(name) {
  env <- new.env(parent = emptyenv())
  utils::data(list = name, package = "lavaan", envir = env)
  env[[name]]
}

holzinger_cfa_model <- function() {
  "
    visual  =~ x1 + x2 + x3
    textual =~ x4 + x5 + x6
    speed   =~ x7 + x8 + x9
  "
}

political_democracy_model <- function() {
  "
    ind60 =~ x1 + x2 + x3
    dem60 =~ y1 + y2 + y3 + y4
    dem65 =~ y5 + y6 + y7 + y8

    dem60 ~ ind60
    dem65 ~ ind60 + dem60

    y1 ~~ y5
    y2 ~~ y4 + y6
    y3 ~~ y7
    y4 ~~ y8
    y6 ~~ y8
  "
}

growth_model <- function() {
  "
    i =~ 1*t1 + 1*t2 + 1*t3 + 1*t4
    s =~ 0*t1 + 1*t2 + 2*t3 + 3*t4

    i ~ x1 + x2
    s ~ x1 + x2

    t1 ~ c1
    t2 ~ c2
    t3 ~ c3
    t4 ~ c4
  "
}

twolevel_model <- function() {
  "
    level: 1
      fw =~ y1 + y2 + y3
      fw ~ x1 + x2 + x3

    level: 2
      fb =~ y1 + y2 + y3
      fb ~ w1 + w2
  "
}

make_ordinal_holzinger <- function() {
  data <- load_lavaan_data("HolzingerSwineford1939")
  for (name in paste0("x", seq_len(9L))) {
    data[[name]] <- ordered(cut(data[[name]], breaks = 3L, labels = FALSE))
  }
  data
}

make_missing_holzinger <- function(seed = 20260511L, missing_rate = 0.08) {
  data <- load_lavaan_data("HolzingerSwineford1939")
  withr::with_seed(seed, {
    n_missing <- max(1L, round(nrow(data) * missing_rate))
    for (name in paste0("x", seq_len(9L))) {
      data[sample.int(nrow(data), size = n_missing), name] <- NA
    }
  })
  data
}

make_large_cfa_data <- function(seed = 20260511L, n = 600L,
                                factors = 6L, items_per_factor = 5L) {
  withr::with_seed(seed, {
    eta <- matrix(stats::rnorm(n * factors), nrow = n, ncol = factors)
    out <- data.frame(matrix(NA_real_, nrow = n,
                             ncol = factors * items_per_factor))
    names(out) <- paste0("x", seq_len(ncol(out)))

    loading <- 0.75
    resid_sd <- sqrt(1 - loading^2)
    for (factor_id in seq_len(factors)) {
      for (item_id in seq_len(items_per_factor)) {
        col_id <- (factor_id - 1L) * items_per_factor + item_id
        out[[col_id]] <- loading * eta[, factor_id] +
          stats::rnorm(n, sd = resid_sd)
      }
    }
    out
  })
}

large_cfa_model <- function(factors = 6L, items_per_factor = 5L) {
  paste(vapply(seq_len(factors), function(factor_id) {
    item_ids <- ((factor_id - 1L) * items_per_factor + 1L):
      (factor_id * items_per_factor)
    paste0("f", factor_id, " =~ ",
           paste0("x", item_ids, collapse = " + "))
  }, character(1L)), collapse = "\n")
}

lavaan_profiling_scenarios <- function() {
  list(
    list(
      label = "cfa_holzinger_ml",
      title = "Holzinger-Swineford CFA",
      family = "cfa",
      estimator = "ML",
      data = "HolzingerSwineford1939",
      notes = "Three-factor CFA with public example data and default SE/test.",
      fit = function() {
        lavaan::cfa(
          holzinger_cfa_model(),
          data = load_lavaan_data("HolzingerSwineford1939"),
          std.lv = TRUE
        )
      }
    ),
    list(
      label = "sem_political_democracy_ml",
      title = "Political Democracy SEM",
      family = "sem",
      estimator = "ML",
      data = "PoliticalDemocracy",
      notes = "Bollen political democracy SEM from the package README.",
      fit = function() {
        lavaan::sem(
          political_democracy_model(),
          data = load_lavaan_data("PoliticalDemocracy")
        )
      }
    ),
    list(
      label = "growth_demo_ml",
      title = "Demo Growth Model",
      family = "growth",
      estimator = "ML",
      data = "Demo.growth",
      notes = "Linear growth model with predictors and time-varying covariates.",
      fit = function() {
        lavaan::growth(
          growth_model(),
          data = load_lavaan_data("Demo.growth")
        )
      }
    ),
    list(
      label = "twolevel_demo_ml",
      title = "Demo Two-Level SEM",
      family = "twolevel",
      estimator = "ML",
      data = "Demo.twolevel",
      notes = "Multilevel model using the public clustered demo dataset.",
      fit = function() {
        lavaan::sem(
          twolevel_model(),
          data = load_lavaan_data("Demo.twolevel"),
          cluster = "cluster"
        )
      }
    ),
    list(
      label = "ordinal_holzinger_wlsmv",
      title = "Ordinal Holzinger-Swineford CFA",
      family = "ordinal",
      estimator = "WLSMV",
      data = "HolzingerSwineford1939 transformed to ordered indicators",
      notes = "Ordered three-category indicators stress WLS sample stats.",
      fit = function() {
        lavaan::cfa(
          holzinger_cfa_model(),
          data = make_ordinal_holzinger(),
          ordered = paste0("x", seq_len(9L)),
          estimator = "WLSMV",
          std.lv = TRUE
        )
      }
    ),
    list(
      label = "fiml_holzinger_ml",
      title = "FIML Missing-Data CFA",
      family = "missing",
      estimator = "ML",
      data = "HolzingerSwineford1939 with deterministic missingness",
      notes = "CFA with 8 percent deterministic missingness and missing='ml'.",
      fit = function() {
        lavaan::cfa(
          holzinger_cfa_model(),
          data = make_missing_holzinger(),
          missing = "ml",
          std.lv = TRUE
        )
      }
    ),
    list(
      label = "multigroup_holzinger_loadings",
      title = "Multigroup Loading-Invariance CFA",
      family = "multigroup",
      estimator = "ML",
      data = "HolzingerSwineford1939 grouped by school",
      notes = "Two-group CFA with equality-constrained loadings.",
      fit = function() {
        lavaan::cfa(
          holzinger_cfa_model(),
          data = load_lavaan_data("HolzingerSwineford1939"),
          group = "school",
          group.equal = "loadings",
          std.lv = TRUE
        )
      }
    ),
    list(
      label = "sample_cov_holzinger_ml",
      title = "Sample-Covariance CFA",
      family = "sample-cov",
      estimator = "ML",
      data = "HolzingerSwineford1939 sample moments",
      notes = "CFA using sample.cov/sample.mean instead of a raw data frame.",
      fit = function() {
        data <- load_lavaan_data("HolzingerSwineford1939")
        ov <- paste0("x", seq_len(9L))
        lavaan::cfa(
          holzinger_cfa_model(),
          sample.cov = stats::cov(data[ov]),
          sample.mean = colMeans(data[ov]),
          sample.nobs = nrow(data),
          meanstructure = TRUE,
          std.lv = TRUE
        )
      }
    ),
    list(
      label = "synthetic_large_cfa_ml",
      title = "Synthetic Larger CFA",
      family = "synthetic",
      estimator = "ML",
      data = "Deterministic 30-indicator synthetic CFA",
      notes = "Six-factor CFA with 600 rows and 30 indicators.",
      fit = function() {
        lavaan::cfa(
          large_cfa_model(),
          data = make_large_cfa_data(),
          std.lv = TRUE
        )
      }
    )
  )
}

profiling_scenario_metadata <- function(scenarios) {
  rows <- lapply(scenarios, function(scenario) {
    data.frame(
      label = scenario$label,
      title = scenario$title,
      family = scenario$family,
      estimator = scenario$estimator,
      data = scenario$data,
      notes = scenario$notes,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}
