#basic working example
testthat::test_that("Returns TRUE when no errors present", {
  data = HolzingerSwineford1939

  #should return a list if it works
  fit <- efa(data = data,
             ov.names = paste("x", 1:9, sep = ""),
             nfactors = 1:3,
             rotation = "geomin",
             rotation.args = list(geomin.epsilon = 0.01, rstarts = 1))

  expect_true(is.list(fit))
})


#no rotation
testthat::test_that("Returns warning when errors present - rotation", {
  data = HolzingerSwineford1939

  expect_error(efa(data = data,
                   ov.names = paste("x", 1:9, sep = ""),
                   nfactors = 1:3,
                   rotation = "",
                   rotation.args = list(geomin.epsilon = 0.01, rstarts = 1)))
})

#no valid factors
testthat::test_that("Returns error message when errors present - nfactors - 1", {
  data = HolzingerSwineford1939

  #Logical nfactors input
  expect_error(efa(data = data,
                   ov.names = paste("x", 1:9, sep = ""),
                   nfactors = TRUE,
                   rotation = "",
                   rotation.args = list(geomin.epsilon = 0.01, rstarts = 1)),
               label = "invalid 'length' argument")
})

testthat::test_that("Returns error message when errors present - nfactors - 2", {
  data = HolzingerSwineford1939

  #NA nfactors input
  expect_error(efa(data = data,
                   ov.names = paste("x", 1:9, sep = ""),
                   nfactors = NA,
                   rotation = "",
                   rotation.args = list(geomin.epsilon = 0.01, rstarts = 1)))
})

testthat::test_that("Returns error message when errors present - nfactors - 3", {
  data = HolzingerSwineford1939

  #0 nfactors input
  expect_error(efa(data = data,
      ov.names = paste("x", 1:9, sep = ""),
      nfactors = 0,
      rotation = "geomin",
      rotation.args = list(geomin.epsilon = 0.01, rstarts = 1)),
    label = "nfactors must be greater than zero")
})

testthat::test_that("Returns error message when errors present - nfactors - 4", {
  data = HolzingerSwineford1939

  ov.names <- paste("x", 1:9, sep = "")
  nvar <- length(ov.names)
  p.star <- nvar * (nvar + 1)/2
  nfac.max <- 0L
  for(nfac in seq_len(nvar)) {
    # compute number of free parameters
    npar <- nfac*nvar + nfac*(nfac+1L)/2 + nvar - nfac^2
    if(npar > p.star) {
      nfac.max <- nfac - 1L
      break
    }
  }

  #nmax nfactors  < input
  expect_error(efa(data = data,
                   ov.names = ov.names,
                   nfactors = 6,
                   rotation = "geomin",
                   rotation.args = list(geomin.epsilon = 0.01, rstarts = 1)),
               label = paste("lavaan ERROR: when nvar = ", nvar,
                             " the maximum number of factors is ", nfac.max, sep = ""))
})

#no ov.names
testthat::test_that("Returns error message when errors present - ov.names - 1", {
  data = HolzingerSwineford1939

  #ov.name undefined
  expect_error(efa(data = data,
                    ov.names = "",
                    nfactors = 1:3,
                    rotation = "geomin",
                    rotation.args = list(geomin.epsilon = 0.01, rstarts = 1)),
               label = "undefined columns selected")
})

testthat::test_that("Returns error message when errors present - ov.names - 2", {
  data = NULL

  #ov.name length == 0L
  expect_error(efa(data = data,
                   ov.names = NULL,
                   nfactors = 1:3,
                   rotation = "geomin",
                   rotation.args = list(geomin.epsilon= 0.01, rstarts = 1)),
               label = "could not extract variable names from data or sample.cov")
})

#invalid groups input
testthat::test_that("Returns error message when errors present - dotdotdot$groups", {
  data = HolzingerSwineford1939

  #invalid dotdotdot$groups
  expect_error(efa(data = data,
                   ov.names = paste("x", 1:9, sep = ""),
                   nfactors = 1:3,
                   rotation = "geomin",
                   rotation.args = list(geomin.epsilon = 0.01, rstarts = 1),
                   group = ""),
               label = "efa has no support for multiple groups (for now)")
})

#invalid output
testthat::test_that("Returns error message when errors present - output - 1", {
  data = HolzingerSwineford1939

  #invalid output
  expect_error(efa(data = data,
                   ov.names = paste("x", 1:9, sep = ""),
                   nfactors = 1:3,
                   rotation = "geomin",
                   rotation.args = list(geomin.epsilon = 0.01, rstarts = 1),
                   output = ""),
               label = "output= must be either \"lavaan\" or \"efa\""
  )
})

#output = 'lavaan' error
testthat::test_that("Returns error message when errors present - output - 2", {
  data = HolzingerSwineford1939

  #invalid nfactors with output lavaan
  expect_error(efa(data = data,
                   ov.names = paste("x", 1:9, sep = ""),
                   nfactors = 1:3,
                   rotation = "geomin",
                   rotation.args = list(geomin.epsilon = 0.01, rstarts = 1),
                   output = "lavaan"),
               label = "when output = \"lavaan\", nfactors must be a single (integer) number.")
})

#invalid data input
testthat::test_that("Returns error message when errors present - data", {
  data = numeric()

  #invalid input data
  expect_error(efa(data = data,
                   ov.names = paste("x", 1:9, sep = ""),
                   nfactors = 1:3,
                   rotation = "geomin",
                   rotation.args = list(geomin.epsilon = 0.01, rstarts = 1))
  )
})


