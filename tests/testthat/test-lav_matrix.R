context("lav_matrix")

A     <- matrix(1:16, nrow=2)
A_sqr <- matrix(1:16, nrow=4)
A_sym <- matrix(1:9, nrow=3); A_sym[upper.tri(A_sym)] <- t(A_sym)[upper.tri(A_sym)]

test_that("lav_matrix_vech matches using lower.tri", {
  expect_identical(
    lav_matrix_vech(A),
    A[lower.tri(A, diag = TRUE)]
  )
})

test_that("lav_matrix_vech without diagonal matches using lower.tri", {
  expect_identical(
    lav_matrix_vech(A, diagonal = FALSE),
    A[lower.tri(A)])
})

test_that("lav_matrix_vech and lav_matrix_vechru are identical on a symmetric matrix", {
  for (diagonal in c(TRUE, FALSE))
    expect_identical(lav_matrix_vech(A_sym, diagonal), lav_matrix_vechru(A_sym, diagonal))
})

test_that("lav_matrix_vechr and lav_matrix_vechu are identical on a symmetric matrix", {
  for (diagonal in c(TRUE, FALSE))
    expect_identical(lav_matrix_vechr(A_sym, diagonal), lav_matrix_vechu(A_sym, diagonal))
})