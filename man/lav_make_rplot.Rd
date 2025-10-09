\name{lav_make_rplot}
\alias{lav_make_rplot}
\title{Plots a diagram.}
\usage{
lav_make_rplot(nodes.edges,
               sloped.labels = TRUE,
               outfile = "",
               addgrid = TRUE,
               mlovcolors = c("lightgreen", "lightblue"),
               lightness = 1,
               italic = TRUE,
               auto.subscript = TRUE
               )
}
\arguments{
\item{nodes.edges}{A list with data.frames \code{nodes} and \code{edges} and
integer \code{mlrij} as returned from \code{lav_position_nodes}.}
\item{sloped.labels}{Logical, sloped labels above (or under) the edges.}
\item{outfile}{Character string naming the file to store the diagram in PNG or
  NA to show the plot in R.}
\item{addgrid}{Logical, add a grid with indicated 'positions' to the graph?}
\item{mlovcolors}{Array of two colors for distinguishing ov nodes with same name
between the levels in a  multilevel model.}
\item{lightness}{A scalar factor to modify the distances between nodes.}
\item{italic}{Are labels to be in italic font? \strong{Attention:} The switch 
to an italic font will only be possible if the label value isn't an expression!}
\item{auto.subscript}{Logical, see \code{\link{lav_format_label}}.}
}
\value{
NULL (invisible)}
\description{
Creates a plot in R showing a simple diagram of the model.
}
\seealso{
\code{\link{lav_position_nodes}}
}
\examples{
model <- 'alpha  =~ x1 + x2 + x_3     # latent variable
          beta <~ x4 + x5 + x6        # composite
          gamma =~ x7 + x8 + x9       # latent variable
          Xi =~ x10 + x11 + x12 + x13 # latent variable
          # regressions
          Xi ~ v * alpha + t * beta + c * 1
          alpha ~ yy * Theta1 + tt * beta + ss * gamma
          '
test <- lav_get_model_info(model)
test1 <- lav_position_nodes(test)
lav_make_rplot(test1, lightness = 1.1)
# better position for constant in regressen Xi, no sloped labels
test2 <- lav_position_nodes(test, placenodes = list(`1vanXi` = c(2, 5)))
lav_make_rplot(test2, FALSE, lightness = 1.1)
modelml <- '
        level: 1
        fw =~ 1*y_1 + y_2 + y_3 + y_5
        level: 2
        fb =~ 1*y_1 + y_2 + y_3 + y_4
        y_2 ~~ cv24 * y_4
  '
test <- lav_get_model_info(modelml)
test <- lav_position_nodes(test)
lav_make_rplot(test, sloped.labels = FALSE, addgrid = FALSE, 
                auto.subscript = FALSE)
\dontrun{
# example where plot is stored in a PNG file
lav_make_rplot(test, sloped.labels = FALSE, addgrid = FALSE,
                auto.subscript = FALSE, outfile="demo_rplot.png")
}
}
