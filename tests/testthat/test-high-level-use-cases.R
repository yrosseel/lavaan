test_that("cfa fits the Holzinger-Swineford three-factor model", {
  data("HolzingerSwineford1939")

  model <- paste(
    "visual =~ x1 + x2 + x3",
    "textual =~ x4 + x5 + x6",
    "speed =~ x7 + x8 + x9",
    sep = "\n"
  )

  fit <- cfa(model, data = HolzingerSwineford1939)

  expect_true(methods::is(fit, "lavaan"))
  expect_true(lavInspect(fit, "converged"))
  expect_equal(as.numeric(fitMeasures(fit, "df")), 24)
  expect_equal(round(as.numeric(fitMeasures(fit, "cfi")), 3), 0.931)
})

test_that("sem fits the Political Democracy model", {
  data("PoliticalDemocracy")

  model <- paste(
    "ind60 =~ x1 + x2 + x3",
    "dem60 =~ y1 + a*y2 + b*y3 + c*y4",
    "dem65 =~ y5 + a*y6 + b*y7 + c*y8",
    "dem60 ~ ind60",
    "dem65 ~ ind60 + dem60",
    "y1 ~~ y5",
    "y2 ~~ y4 + y6",
    "y3 ~~ y7",
    "y4 ~~ y8",
    "y6 ~~ y8",
    sep = "\n"
  )

  fit <- sem(model, data = PoliticalDemocracy)
  paths <- parameterEstimates(fit)
  dem65_from_dem60 <- paths$est[
    paths$op == "~" & paths$lhs == "dem65" & paths$rhs == "dem60"
  ]

  expect_true(lavInspect(fit, "converged"))
  expect_equal(as.numeric(fitMeasures(fit, "df")), 38)
  expect_equal(round(as.numeric(fitMeasures(fit, "cfi")), 3), 0.997)
  expect_equal(round(dem65_from_dem60, 3), 0.865)
})

test_that("growth fits a linear growth curve model", {
  data("Demo.growth")

  model <- paste(
    "i =~ 1*t1 + 1*t2 + 1*t3 + 1*t4",
    "s =~ 0*t1 + 1*t2 + 2*t3 + 3*t4",
    sep = "\n"
  )

  fit <- growth(model, data = Demo.growth)
  means <- parameterEstimates(fit)
  means <- means[
    means$op == "~1" & means$lhs %in% c("i", "s"),
    c("lhs", "est")
  ]

  expect_true(lavInspect(fit, "converged"))
  expect_equal(as.numeric(fitMeasures(fit, "df")), 5)
  expect_equal(round(as.numeric(fitMeasures(fit, "rmsea")), 3), 0.039)
  expect_equal(
    stats::setNames(round(means$est, 3), means$lhs),
    c(i = 0.615, s = 1.006)
  )
})
