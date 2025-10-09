\name{lav_format_label}
\alias{lav_format_label}
\title{Provides formatted label for use in tikz, svg or rplot}
\usage{
lav_format_label(label = "", value = "", show =FALSE,
                 idx.font.size = 20L, dy = 7L,
                 italic = TRUE, auto.subscript = TRUE)
}
\arguments{
\item{label}{A character string in one of the formats
\itemize{
\item \var{name}
\item \var{name_index}
\item \var{name=value}
\item \var{name_index=value}
}
If \var{value} is specified in parameter \code{label} and parameter
\code{value} is empty, the value in \code{label} will be used.}
\item{value}{A character string specifying a value or empty.}
\item{show}{A logical indicating that the result should be plotted in an Rplot.}
\item{idx.font.size}{An integer specifying font size to use for the
subscript in svg.}
\item{dy}{An integer specifying the distance to move the baseline of
the subscript in svg.}
\item{italic}{A logical to indicate if the labels font should be cursive.
Only used for tikz output and in shown Rplot if result isn't an expression.}
\item{auto.subscript}{Logical, if TRUE and \code{label} starts with one or more
alhabetic chars followed by one or more digits, an underscore will be
inserted between these parts.}
}
\value{
a list with members \code{svg}, \code{tikz} and \code{r}, giving the result.
}
\description{
Creates the code to use in tikz, svg or R to show the label. The label is
printed in an xy-plot if show = TRUE.
}
\details{
If both \code{label} and \code{value} are empty, the resulting codes are
also empty and nothing will be shown.

If \code{label} is empty and \code{value} is not, processing is done as if
\code{label} was specified and \code{value} was empty.

If \code{label} contains the string "1van", the label value is set to "1", this
is to allow distinct names for regression intercepts and still label them as "1".

If \var{name} in \code{label} is a Greek character or \code{varepsilon},
it is attempted to generate code which shows the Greek symbol. If there is an
\var{index} part in \code{label}, it is attempted to generate code which
displays this value as a subscript. In the svg code the values
\code{idx.font.size} and \code{dy} are used for the subscript in the
generated code. If a \var{value} is present, it is attempted
to show this value after the label, with an equal sign in between both.
}
\examples{
lav_format_label("x3")
lav_format_label("beta10", 0.65, show = TRUE)
lav_format_label("A_i,j=0.45", show = TRUE)
lav_format_label("Gamma")
lav_format_label(value="1.2345")
}
