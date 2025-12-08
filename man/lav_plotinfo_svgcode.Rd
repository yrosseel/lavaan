\name{lav_plotinfo_svgcode}
\alias{lav_plotinfo_svgcode}
\title{Creates code to show a diagram in svg format.}
\usage{
lav_plotinfo_svgcode(plotinfo,
                     outfile = "",
                     sloped.labels = TRUE,
                     standalone = FALSE,
                     stroke.width = 2L,
                     font.size = 20L,
                     idx.font.size = 15L,
                     dy = 5L,
                     mlovcolors = c("lightgreen", "lightblue"),
                     lightness= 1,
                     font.family = "Latin Modern Math, arial, Aerial, sans",
                     italic = TRUE,
                     auto.subscript = TRUE
                     )
}
\arguments{
\item{plotinfo}{A plotinfo structure as returned from \code{lav_plotinfo_positions}.}
\item{outfile}{A connection or a character string, the file to store the code.}
\item{sloped.labels}{logical, sloped labels for the edges.}
\item{standalone}{logical, add code to make a html file with svg in it? If
FALSE (default) the outfile - if specified - must have a file extension
'svg', if TRUE the file extension must be 'htm' or 'html'.}
\item{stroke.width}{Value for stroke-width parameter in svg.}
\item{font.size}{An integer specifying normal font size to use.}
\item{idx.font.size}{An integer specifying font size to use for a subscript.}
\item{dy}{An integer specifying the distance to move the baseline of the subscript.}
\item{mlovcolors}{Array of two colors for distinguishing ov nodes with same name
between the levels in a  multilevel model.}
\item{lightness}{A scalar factor to modify the distances between nodes.}
\item{font.family}{Fonts to be tried for rendering the labels. The first one,
Latin Modern Math, is close to the default font used in tikz for mathematical
expressions and can be downloaded from
\href{https://ctan.org/tex-archive/fonts/lm-math}{CTAN}.}
\item{italic}{Are labels to be in italic font? If FALSE the labels will be
shown in mathrm font.}
\item{auto.subscript}{Logical, see \code{\link{lav_label_code}}.}
}
\value{
NULL (invisible)}
\description{
Create svg code to show a diagram of the model.
}
\seealso{
\code{\link{lav_plotinfo_positions}}
}
\examples{
model <- 'alpha  =~ x1 + x2 + x3       # latent variable
          beta <~ x4 + x5 + x6         # composite
          gamma =~ x7 + x8 + x9        # latent variable
          Xi =~ x10 + x11 + x12 + x13  # latent variable
          # regressions
          Xi ~ v * alpha + t * beta + 1
          alpha ~ yy * Theta1 + tt * beta + ss * gamma
          '
test <- lav_model_plotinfo(model)
test <- lav_plotinfo_positions(test)
lav_plotinfo_svgcode(test) # no file given, so output to R console
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
lav_plotinfo_svgcode(testml, sloped.labels = FALSE, standalone = TRUE,
            auto.subscript = FALSE)
\dontrun{
# example creating html file with the above diagrams
zz <- file("demosvg.html", open="w")
writeLines(c(
  '<!DOCTYPE html>',
  '<html>',
  '<body>',
  '<h2>SVG diagrams created by lav_plot R package</h2>'),
  zz)
lav_plotinfo_svgcode(test, outfile = "temp.svg")
tmp <- readLines("tmp.svg")
writeLines(tmp, zz)
writeLines("<br />", zz)
lav_plotinfo_svgcode(testml, outfile = "temp.svg", sloped.labels = FALSE,
             standalone = TRUE)
tmp <- readLines("tmp.svg")
writeLines(tmp, zz)
writeLines(c("</body>", "</html>"), zz)
close(zz)
browseURL("demosvg.html")
}
}
