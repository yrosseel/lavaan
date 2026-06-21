\name{lav_label_code}
\alias{lav_label_code}
\title{Provides a formatted label for use in tikz, svg or Rplot}
\usage{
lav_label_code(label = "", value = "", show =FALSE,
               idx_font_size = 20L, dy = 7L,
               italic = TRUE, auto_subscript = TRUE)
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
\item{idx_font_size}{An integer specifying font size to use for the
subscript in svg.}
\item{dy}{An integer specifying the distance to move the baseline of
the subscript in svg.}
\item{italic}{A logical indicating whether the label's font should be italic.
Only used for tikz output and in the displayed Rplot when the result is not an expression.}
\item{auto_subscript}{Logical, if TRUE and \code{label} starts with one or more
alphabetic chars followed by one or more digits, an underscore will be
inserted between these parts.}
}
\value{
a list with members \code{svg}, \code{tikz} and \code{r}, giving the result.
}
\description{
Creates the code used to show the label in tikz, svg or R. The label is
plotted in an xy-plot if show = TRUE.
}
\details{
If both \code{label} and \code{value} are empty, the resulting codes are
also empty and nothing will be shown.

If \code{label} is empty and \code{value} is not, processing is done as if
\code{label} was specified and \code{value} was empty.

If \code{label} contains the string "1van", the label value is set to "1".
This allows distinct names for regression intercepts while still labeling
them as "1".

If \var{name} in \code{label} is a Greek character or \code{varepsilon},
the function attempts to generate code that shows the Greek symbol. If
\code{label} contains an \var{index} part, the function attempts to generate
code that displays this value as a subscript. In the svg code, the values
\code{idx_font_size} and \code{dy} are used for the subscript. If a
\var{value} is present, the function attempts to show this value after the
label, with an equal sign between the two.
}
\examples{
lav_label_code("x3")
lav_label_code("beta10", 0.65, show = TRUE)
lav_label_code("A_i,j=0.45", show = TRUE)
lav_label_code("Gamma")
lav_label_code(value="1.2345")
}
