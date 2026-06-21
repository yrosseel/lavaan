\name{lav_plotinfo_positions}
\alias{lav_plotinfo_positions}
\title{Position the nodes in the diagram}
\usage{
lav_plotinfo_positions(plotinfo,
                       placenodes = NULL,
                       edgelabelsbelow = NULL,
                       group_covar_indicators = FALSE,
                       debug = FALSE)
}
\arguments{
\item{plotinfo}{The plotinfo structure as returned from \code{lav_model_plotinfo}.}
\item{placenodes}{optional list with members \code{\var{nodename} = c(row, column)},
\code{row} and \code{column} don't have to be integers.}
\item{edgelabelsbelow}{optional list with members \code{c(\var{nodename1},
\var{nodename2})}.}
\item{group_covar_indicators}{logical, should items with indicators which
have an explicit covariance link be placed in the same group, i.e. forced
to be on the same side of the diagram?}
\item{debug}{logical, print debug information?}
}
\value{
A plotinfo structure with modified \code{nodes} and \code{edges} 
data.frames, and an integer \code{mlrij} giving the position at which a line
should be drawn for multilevel models.

The data.frames have new columns defined as follows:
\enumerate{
\item \code{nodes}
  \describe{
    \item{rij}{index of the row where the node will be placed, initialized NA.}
    \item{kolom}{index of the column where the node will be placed,
    initialized NA.}
  }
\item \code{edges}
  \describe{
    \item{vananker}{character, anchor point for starting node, initialized NA.}
    \item{naaranker}{character, anchor point for destination node,
    initialized NA.}
    \item{controlpt.kol}{real, column position of control point if the
    edge has to be drawn as a quadratic Bezier curve.}
    \item{controlpt.row}{real, row position of control point if the
    edge has to be drawn as a quadratic Bezier curve.}
    \item{labelbelow}{logical, TRUE if label has to be positioned under 
    the line, initialized FALSE.}
  }
}
}
\description{
Computes the positions for the nodes and anchors and control points for the 
edges in the diagram.
}
\details{
This function \strong{tries} to arrange the nodes and the anchor points for the
edges in a way that keeps the diagram reasonably compact, while placing the
dependent variable of a regression to the right of the variables it depends on.
If the result is not what you want, you can inspect the plot with the
\code{lav_plotinfo_rgraph} or \code{lav_plot} functions and then place some
nodes at a different location (\code{placenodes}) and/or move some edge labels
to the other side of the edge (\code{labelsbelow}).

For multilevel models -- maximum two levels -- the nodes in block 2 are grouped
at the top of the diagram and those in block 1 at the bottom. The item 
\code{mlrij} gives the position of the separation between the blocks.
}

\seealso{
\code{\link{lav_model_plotinfo}}, \code{\link{lav_plotinfo_rgraph}},
\code{\link{lav_plot}}
}

\examples{
model <- 'alfa  =~ 1 * x1 + x2 + x3         # latent variable
          beta <~ x4 + x5 + x6              # composite
          gamma =~ 1 * x7 + x8 + x9         # latent variable
          Xi =~ 1 * x10 + x11 + x12 + x13   # latent variable
          # regressions
          Xi ~ v * alfa + t * beta + cc * 1
          alfa ~ tt * beta + ss * gamma + yy * Theta1
          # variances and covariances
          x2 ~~ cc25 * x5
          x3 ~~ cc36 * x6
          x3 ~~ cc34 * x4
          gamma ~~ 0.55 * gamma
          '
test <- lav_model_plotinfo(model)
(test_positioned <- lav_plotinfo_positions(test))
}
