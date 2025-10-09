\name{lav_get_model_info}
\alias{lav_get_model_info}
\title{Get model information}
\usage{
lav_get_model_info(model = NULL, infile = NULL, varlv = FALSE)
}
\arguments{
\item{model}{A character vector specifying the model in lavaan syntax or a list
(or data.frame) with at least members lhs, op, rhs, label and fixed or a fitted
lavaan object (in which case the \code{ParTable} object is extracted and column
\code{est} is used as value to show).
Should be \code{NULL} if infile is given.}
\item{infile}{A character string specifying the file which contains the model
syntax.}
\item{varlv}{A logical indicating that the (residual) variance of a variable
should be plotted as a seperate latent variable (with a smaller circle then
ordinary latent variables). In this case covariances between two such variables
will be plotted as covariance between the latent variables for the variance.}
}
\value{
A list with members nodes and edges, containing the necessary data to create
a diagram. The return value consists of a list of two data.frames:
\enumerate{
\item nodes
  \describe{
    \item{id}{integer, autoincrement identification of the node.}
    \item{naam}{character, name of the node as specified in the model.
    For intercepts the name is "1vanXXXX", with XXXX the name
    of the regressed variable.}
    \item{tiepe}{character, type of node: ov (observed variable),
    lv (latent variable), varlv (variance as latent variable),
    cv (composite variable), wov (within level variable in multilevel model),
    bov (between level variable in multilevel model), 
    const (intercept of regression).}
    \item{blok}{integer, level (0 if not a multilevel model).}
    \item{voorkeur}{character, preferred place in the plot, l(left), r(right)
    or m = (middle).}
  }
\item edges
  \describe{
    \item{id}{integer, autoincrement identification of the edge.}
    \item{label}{character, label for the edge, made from the label specified in
    the model and the fixed (or estimated) value if present.}
    \item{van}{integer, id of the starting node.}
    \item{naar}{integer, id of the destination node.}
    \item{tiepe}{character, lavaan operator, except (residual)
      variance, coded here as '~~~'.}
  }
}
}

\description{
Gets information from model necessary to produce a plot.
}

\examples{
model <- 'alpha  =~ 1 * x1 + x2 + x3        # latent variable
          beta <~ x4 + x5 + x6              # composite
          gamma =~ 1 * x7 + x8 + x9         # latent variable
          Xi =~ 1 * x10 + x11 + x12 + x13   # latent variable
          # regressions
          Xi ~ v * alpha + t * beta + cc * 1
          alpha ~ tt * beta + ss * gamma + yy * Theta1
          # variances and covariances
          x2 ~~ cc25 * x5
          x3 ~~ cc36 * x6
          x3 ~~ cc34 * x4
          gamma ~~ 0.55 * gamma
          '
(test <- lav_get_model_info(model))
}
