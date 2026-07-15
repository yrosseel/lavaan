# lavaan

lavaan is a free, open source R package for latent variable analysis. You can
use lavaan to estimate a large variety of multivariate statistical models,
including path analysis, confirmatory factor analysis, structural equation
modeling and growth curve models.

The lavaan package is developed to provide useRs, researchers and teachers a
free open-source, but commercial-quality package for latent variable modeling.
The long-term goal of lavaan is to implement all the state-of-the-art
capabilities that are currently available in commercial packages. However,
lavaan is still under development, and much work still needs to be done.

## Job ad

As part of an NWO-funded project for Dutch open-science infrastructural
development, we are looking for a (second) 
**researcher/junior postdoc for the lavaan project**. 
The job ad can be found on the University of [Twente](https://utwentecareers.nl/en/vacancies/2525/researcher-in-methods-and-statistics-specializing-in-programming-and-structural-equation-modelling/).

The deadline is *31 July 2026*.

The proposal for the full project
can be found on [Zenodo](https://doi.org/10.5281/zenodo.19552779). The post-doc
position above is part of Work Package 2 (WP2).


## Installation

Install the stable version from CRAN:

```R
install.packages("lavaan")
```

Or the development version from GitHub:

```R
# install.packages("remotes")
remotes::install_github("yrosseel/lavaan")
```

## Usage

To get a first impression of how lavaan works in practice, consider the
following example of a SEM model (the Political Democracy Example from 
Bollen's 1989 book):

```R
library(lavaan)

model <- '
   # latent variable definitions
     ind60 =~ x1 + x2 + x3
     dem60 =~ y1 + y2 + y3 + y4
     dem65 =~ y5 + y6 + y7 + y8
   # regressions
     dem60 ~ ind60
     dem65 ~ ind60 + dem60
   # residual covariances
     y1 ~~ y5
     y2 ~~ y4 + y6
     y3 ~~ y7
     y4 ~~ y8
     y6 ~~ y8
'
fit <- sem(model, data = PoliticalDemocracy)
summary(fit)
```

More information can be found on the website: https://lavaan.org

