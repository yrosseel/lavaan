\name{lav_plotinfo_tikzcode}
\alias{lav_plotinfo_tikzcode}
\title{Creates code to show a diagram in tikz format.}
\usage{
lav_plotinfo_tikzcode(plotinfo,
                      outfile = "",
                      cex = 1.3,
                      sloped.labels = TRUE,
                      standalone = FALSE,
                      mlovcolors = c("lightgreen", "lightblue"),
                      lightness = 1,
                      italic = TRUE,
                      auto.subscript = TRUE
                     )
}
\arguments{
\item{plotinfo}{A plotinfo structure as returned from \code{lav_plotinfo_positions}.}
\item{outfile}{A connection or a character string, the file to store the code.}
\item{cex}{Minimum distance between nodes in cm.}
\item{sloped.labels}{logical, sloped labels for the edges.}
\item{standalone}{logical, add code to make teX file standalone?}
\item{mlovcolors}{Array of two colors for distinguishing ov nodes with same name
between the levels in a  multilevel model.}
\item{lightness}{A scalar factor to modify the distances between nodes.}
\item{italic}{Are labels to be in italic font? If FALSE the labels will be
shown in mathrm font.}
\item{auto.subscript}{Logical, see \code{\link{lav_label_code}}.}
}
\value{
NULL, invisible
}
\description{
Creates the code to make a diagram in tikz.
}
\seealso{
\code{\link{lav_plotinfo_positions}}
}
\examples{
model <- 'alpha  =~ x1 + x2 + x3       # latent variable
          beta <~ x4 + x5 + x6         # composite
          gamma =~ x7 + x8 + x9        # latent variable
          Xi =~ x10 + x11 + x12 + x13 # latent variable
          # regressions
          Xi ~ v * alpha + t * beta + 1
          alpha ~ yy * Theta1 + tt * beta + ss * gamma
          '
test <- lav_model_plotinfo(model)
test <- lav_plotinfo_positions(test)
lav_plotinfo_tikzcode(test) # no file given, so output to R console
modelml <- '
        level: 1
        fw =~ 1*y_1 + y_2 + y_3 + y_5
        level: 2
        fb =~ 1*y_1 + y_2 + y_3 + y_4
        y_2 ~~ cv24 * y_4
  '
testml <- lav_model_plotinfo(modelml)
testml <- lav_plotinfo_positions(testml)
# in the line hereunder no file is given, so output to R console
lav_plotinfo_tikzcode(testml, cex = 1.4, sloped.labels = FALSE,
        standalone = TRUE, auto.subscript = FALSE)
\dontrun{
# example creating tex file with the above diagrams
zz <- file("testtikz.tex", open="w")
writeLines(c(
  '\\documentclass{article}',
  '\\usepackage{amsmath, amssymb}',
  '\\usepackage{amsfonts}',
  '\\usepackage[utf8]{inputenc}',
  '\\usepackage[english]{babel}',
  '\\usepackage{xcolor}',
  '\\usepackage{color}',
  '\\usepackage{tikz}',
  '\\usetikzlibrary {shapes.geometric}',
  '\\begin{document}'),
  zz)
lav_plotinfo_tikzcode(test, outfile = "temp.tex")
tmp <- readLines("tmp.tex")
writeLines(tmp, zz)
writeLines(" ", zz)
lav_plotinfo_tikzcode(testml, outfile = "temp.tex", cex = 1.4,
              sloped.labels = FALSE, auto.subscript = FALSE)
tmp <- readLines("tmp.tex")
writeLines(tmp, zz)
writeLines("\\end{document}", zz)
close(zz)
openPDF <- function(f) {
  os <- .Platform$OS.type
  if (os=="windows")
    shell.exec(normalizePath(f))
  else {
    pdf <- getOption("pdfviewer", default='')
    if (nchar(pdf)==0)
      stop("The 'pdfviewer' option is not set. Use options(pdfviewer=...)")
    system2(pdf, args=c(f))
  }
}
tools::texi2dvi("testtikz.tex", pdf = TRUE, clean = TRUE)
openPDF("testtikz.pdf")
}
}
