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

## News:

As part of an NWO-funded project for Dutch open-science infrastructural development, we are looking for a **post-doc/programmer for the lavaan project**. The position is initially for 1 year, but conditional on positive first-year review, it is designed to be extended to as many as 4 years.
The job ad can be found on the University of [Amsterdam (UVA) website](https://werkenbij.uva.nl/en/vacancies/researcher-in-methods-and-statistics-specializing-in-programming-and-structural-equation-modeling-netherlands-15033),
but it is also published on [AcademicTransfer.com](https://www.academictransfer.com/en/jobs/360501/researcher-in-methods-and-statistics-specializing-in-programming-and-structural-equation-modeling).
And the same keywords can be searched to find the job posting on LinkedIn.
The deadline is already *17 May 2026*, but the deadline might be extended if sufficient candidates are not already found.
The proposal for the full project can be found on [Zenodo](https://doi.org/10.5281/zenodo.19552779). The post-doc position above is part of Work Package 1 (WP1).


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

