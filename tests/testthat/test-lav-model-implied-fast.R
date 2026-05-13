test_that("lav_model_implied_fast matches individual LISREL moment helpers", {
  fit <- cfa(
    "visual =~ x1 + x2 + x3",
    data = HolzingerSwineford1939,
    meanstructure = TRUE,
    se = "none",
    test = "none"
  )

  lavmodel <- fit@Model
  glist <- lavmodel@GLIST
  fast <- lavaan:::lav_model_implied_fast(
    lavmodel = lavmodel,
    glist = glist,
    need_sigma = TRUE,
    need_mu = TRUE,
    extra = TRUE
  )

  mm_idx <- lavaan:::lav_model_get_mm_idx(lavmodel)
  ref_sigma <- vector("list", lavmodel@nblocks)
  ref_mu <- vector("list", lavmodel@nblocks)
  for (g in seq_len(lavmodel@nblocks)) {
    mlist <- glist[mm_idx[[g]]]
    ref_sigma[[g]] <- lavaan:::lav_model_sigma_extra(
      lavaan:::lav_lisrel_sigma(mlist = mlist),
      check_sigma_pd = get0(
        "opt_check_sigma_pd",
        lavaan:::lavaan_cache_env,
        ifnotfound = "chol"
      )
    )
    ref_mu[[g]] <- lavaan:::lav_lisrel_mu(mlist = mlist)
  }

  expect_equal(fast$sigma, ref_sigma, tolerance = 1e-12)
  expect_equal(fast$mu, ref_mu, tolerance = 1e-12)
})

test_that("lav_model_implied_fast matches threshold and slope helpers", {
  dat <- HolzingerSwineford1939
  dat$x1 <- ordered(cut(dat$x1, breaks = 3L, labels = FALSE))
  dat$x2 <- ordered(cut(dat$x2, breaks = 3L, labels = FALSE))
  dat$x3 <- ordered(cut(dat$x3, breaks = 3L, labels = FALSE))

  categorical_fit <- cfa(
    "visual =~ x1 + x2 + x3",
    data = dat,
    ordered = c("x1", "x2", "x3"),
    parameterization = "delta",
    se = "none",
    test = "none"
  )
  categorical_model <- categorical_fit@Model
  categorical_glist <- categorical_model@GLIST
  categorical_fast <- lavaan:::lav_model_implied_fast(
    lavmodel = categorical_model,
    glist = categorical_glist,
    need_th = TRUE
  )

  mm_idx <- lavaan:::lav_model_get_mm_idx(categorical_model)
  ref_th <- vector("list", categorical_model@nblocks)
  for (g in seq_len(categorical_model@nblocks)) {
    ref_th[[g]] <- lavaan:::lav_lisrel_th(
      mlist = categorical_glist[mm_idx[[g]]],
      th_idx = categorical_model@th.idx[[g]]
    )
  }
  expect_equal(categorical_fast$th, ref_th, tolerance = 1e-12)

  regression_fit <- sem(
    "x1 ~ x2 + x3",
    data = HolzingerSwineford1939,
    meanstructure = TRUE,
    fixed.x = TRUE,
    conditional.x = TRUE,
    se = "none",
    test = "none"
  )
  regression_model <- regression_fit@Model
  regression_glist <- regression_model@GLIST
  regression_fast <- lavaan:::lav_model_implied_fast(
    lavmodel = regression_model,
    glist = regression_glist,
    need_pi = TRUE
  )

  mm_idx <- lavaan:::lav_model_get_mm_idx(regression_model)
  ref_pi <- vector("list", regression_model@nblocks)
  for (g in seq_len(regression_model@nblocks)) {
    ref_pi[[g]] <- lavaan:::lav_lisrel_pi(
      mlist = regression_glist[mm_idx[[g]]]
    )
  }
  expect_equal(regression_fast$pi, ref_pi, tolerance = 1e-12)
})
