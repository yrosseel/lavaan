# small utility functions to deal with PRELIS
# Y.R.: 11 dec 2012
prelis.read.cor <- function(file = "") {
  # read in numbers as characters
  txt <- scan(file, what = "character", quiet = TRUE)

  # convert to numbers
  txt <- gsub("D", "e", txt)
  x <- as.numeric(txt)

  # create COR/COR matrix
  COR <- lav_matrix_lower2full(x, diagonal = TRUE)
  COR
}


prelis.read.acm <- function(file = "", rescale = 1e-3) {
  # read in raw data -- ignore first three elements
  # first element: 123.456789 (check?)
  # second element: 2.72 version number of prelis
  # third element: almost zero??
  zz <- file(file, "rb")
  raw <- readBin(zz, what = "double", n = 1e+05)[-c(1, 2, 3)]
  close(zz)

  # scale numbers
  raw <- raw * rescale

  ACM <- lav_matrix_lower2full(raw, diagonal = TRUE)

  # elements are divided by 2??
  ACM <- ACM * 2
  ACM
}

prelis.write.data <- function(data, file = "prelis", na.rm = TRUE,
                              labels = FALSE, std.ov = FALSE) {
  dfile <- paste(file, ".raw", sep = "")
  write.table(data,
    file = dfile, na = "-999999", col.names = FALSE,
    row.names = FALSE, quote = FALSE
  )
  if (labels) {
    lfile <- paste(file, ".lab", sep = "")
    write.table(unique(names(data)),
      file = lfile, row.names = F,
      col.names = F, quote = F
    )
  }
}

prelis.run <- function(X, type = "OR", keep.files = FALSE) {
  label <- names(X)
  nvar <- ncol(X)

  # write raw data
  prelis.write.data(X, file = "prelistmp")

  # write syntax
  txt <- paste("DA NI=", nvar, " NO=0 MI=-999999\n", sep = "")
  txt <- paste(txt, "LA", sep = "")
  tmp <- 0
  for (i in 1:nvar) {
    if (tmp %% 6 == 0) txt <- paste(txt, "\n", sep = "")
    txt <- paste(txt, label[i], " ", sep = "")
    tmp <- tmp + 1
  }
  txt <- paste(txt, "\n")
  txt <- paste(txt, "RA FI=prelistmp.raw\n", sep = "")
  txt <- paste(txt, type, " ALL\n", sep = "")
  txt <- paste(txt, "OU MA=PM SA=prelistmp.acm SM=prelistmp.cor\n", sep = "")
  writeLines(txt, con = "prelistmp.in")

  # run prelis
  system("prelis prelistmp.in prelistmp.out")

  # read in acm and cor
  ACM <- prelis.read.acm(file = "prelistmp.acm")
  COR <- prelis.read.cor(file = "prelistmp.cor")

  # clean up
  if (!keep.files) {
    unlink(c(
      "prelistmp.in", "prelistmp.out", "prelistmp.acm",
      "prelistmp.cor", "prelistmp.FREQ", "prelistmp.raw"
    ))
  }

  list(COR = COR, ACM = ACM)
}
