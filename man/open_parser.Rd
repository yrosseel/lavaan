\name{open_parser}
\alias{open_parser}
\alias{lav_parse_options}
\title{The open parser configuration}
\description{
The lavaan model syntax parser can be extended with new operators, modifiers and groupings.
This extensions are only used in the 'open' parser.} 
\usage{
lav_parse_options(
    pkgname   = NULL,
    operators = NULL,
    modifiers = NULL,
    groupings  = NULL
  )
}
\arguments{
\item{pkgname}{The name of the package that extends the elements.}
\item{operators}{A character vector with new operators. These operators must start and end with
a percentage sign and be at least 3 characters long.}
\item{modifiers}{A data.frame with columns \code{mod}, \code{side} and \code{expect}. The elements in \code{mod}
must be character strings in snake_case. The elements in \code{side} must be \code{l} or \code{r}, indicating
whether the modifier can be expected on the left hand of the operator or on the right hand side. The elements in
\code{expect} must be \code{char}, \code{num} or \code{raw}, specifying the value for the modifier to expect, where
\code{raw} means the value is not evaluated in the parser and returned as an uninterpreted string.}
\item{groupings}{A character vector with new groupings, they must be snake_case.}
}
\value{
  A \code{list} with \code{data.frames} containing the currently active parser configuration.
}
\examples{
curconfig <- lav_parse_options("testPKG", 
    operators = c('\%w\%', '\%q\%'),
    modifiers = data.frame(mod = c("linksc", "linksn", "rechtsn"),
                           side = c("l", "l", "r"),
                           expect = c("char", "num", "raw")),
    groupings = c("color", "colour")
)
for (j in seq_along(curconfig)) {
    cat(names(curconfig)[j], ":\\n")
    print(curconfig[[j]])
}
model <- '
     ind60 =~ x1 + x2 + x3
     dem60 =~ y1 + a*y2 + start(0.3)*b*y3 + upper(0.99)*c*y4
     dem65 =~ y5 + a*y6 + b*y7 + rechtsn(2 * pi / 3)*c*y8
     dem60 \%w\% aa *0.2*y4

    dem60 ~ ind60
    dem65 ~ ind60 + dem60
    linksc("klopt") * dem60 =~ aaa * y8
    linksn(sqrt(23.5)) * dem60 =~ bbb * y9
    efa("efatje") * dem65 =~ xyz * x3

    y1 ~~ y5
    y2 ~~ y4 + y6
    y3 ~~ y7
    y4 ~~ y8
    y6 ~~ y8
'
print(result <- lavParseModelString(model, TRUE, parser = "open"))
print(attributes(result))

model2 <- '
  color: green  
    level: 1
        fw =~ y1 + y2 + y3
        fw ~ x1 + x2 + x3
    level: 2
        fb =~ y1 + y2 + y3
        fb ~ w1 + w2
  color: blue  
    level: 1
        fw =~ y1 + y2 + y3
        fw ~ x1 + x2 + x3
    level: 2
        fb =~ y1 + y2 + y3
        fb ~ w1 + w2
'
print(result2 <- lavParseModelString(model2, TRUE, parser = "open"))
print(attributes(result2))
}
