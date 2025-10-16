.onAttach <- function(libname, pkgname) {
  version <- read.dcf(
    file = system.file("DESCRIPTION", package = pkgname),
    fields = "Version"
  )[1]
  packageStartupMessage(
    "This is ", paste(pkgname, version), "\n",
    pkgname, " is FREE software! Please report any bugs."
  )
}
