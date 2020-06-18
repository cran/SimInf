

![](https://raw.githubusercontent.com/stewid/SimInf/master/logo/logo.png)

> A flexible and efficient framework for data-driven stochastic disease spread simulations

[![Build Status](https://dev.azure.com/stefanwidgren/SimInf/_apis/build/status/stewid.SimInf?branchName=master)](https://dev.azure.com/stefanwidgren/SimInf/_build)
[![CRAN status](https://www.r-pkg.org/badges/version/SimInf)](https://CRAN.R-project.org/package=SimInf)
[![CRAN RStudio mirror downloads](https://cranlogs.r-pkg.org/badges/last-month/SimInf)](https://CRAN.R-project.org/package=SimInf)
[![Code coverage](https://img.shields.io/azure-devops/coverage/stefanwidgren/SimInf/1)](https://dev.azure.com/stefanwidgren/SimInf/_build/latest?definitionId=1&branchName=master)

# SimInf

The package provides an efficient and very flexible framework to
conduct data-driven epidemiological modeling in realistic large scale
disease spread simulations. The framework integrates infection
dynamics in subpopulations as continuous-time Markov chains using the
Gillespie stochastic simulation algorithm and incorporates available
data such as births, deaths and movements as scheduled events at
predefined time-points. Using C code for the numerical solvers and
'OpenMP' (if available) to divide work over multiple processors
ensures high performance when simulating a sample outcome. One of our
design goals was to make the package extendable and enable usage of
the numerical solvers from other R extension packages in order to
facilitate complex epidemiological research. The package contains
template models and can be extended with user-defined models.

## Getting started

You can use one of the predefined compartment models in SimInf, for
example, SEIR. But you can also define a custom model 'on the fly'
using the model parser method `mparse`. The method takes a character
vector of transitions in the form of `X -> propensity -> Y` and
automatically generates the C and R code for the model. The left hand
side of the first arrow (`->`) is the initial state, the right hand
side of the last arrow (`->`) is the final state, and the propensity
is written between the two arrows. The flexibility of the `mparse`
approach allows for quick prototyping of new models or features. To
illustrate the `mparse` functionality, let us consider an SIR model in
a closed population i.e., no births or deaths. Let `beta` denote the
transmission rate of spread between a susceptible individual and an
infectious individual and `gamma` the recovery rate from infection
(`gamma` = 1 / average duration of infection). The model can be
described as:


```r
library(SimInf)

transitions <- c("S -> beta*S*I/(S+I+R) -> I",
                 "I -> gamma*I -> R")
compartments <- c("S", "I", "R")
```

The `transitions` and `compartments` variables together with the
constants `beta` and `gamma` can now be used to generate a model with
`mparse`. The model also needs to be initialised with the initial
condition `u0` and `tspan`, a vector of time points where the state of
the system is to be returned. Let us create a model that consists of
1000 replicates of a population, denoted a *node* in SimInf, that each
starts with 99 susceptibles, 5 infected and 0 recovered individuals.


```r
n <- 1000
u0 <- data.frame(S = rep(99, n), I = rep(5, n), R = rep(0, n))

model <- mparse(transitions = transitions,
                compartments = compartments,
                gdata = c(beta = 0.16, gamma = 0.077),
                u0 = u0,
                tspan = 1:180)
```

To generate data from the model and then print some basic information
about the outcome, run the following commands:




```r
result <- run(model)
result
```

```
#> Model: SimInf_model
#> Number of nodes: 1000
#> Number of transitions: 2
#> Number of scheduled events: 0
#> 
#> Global data
#> -----------
#>  Parameter Value
#>  beta      0.160
#>  gamma     0.077
#> 
#> Compartments
#> ------------
#>      Min. 1st Qu. Median   Mean 3rd Qu.   Max.
#>  S   1.00   17.00  27.00  37.91   53.00  99.00
#>  I   0.00    0.00   2.00   5.74    9.00  47.00
#>  R   0.00   36.00  73.00  60.34   85.00 103.00
```

There are several functions in SimInf to facilitate analysis and
post-processing of simulated data, for example, `trajectory`,
`prevalence` and `plot`. The default `plot` will display the median
count in each compartment across nodes as a colored line together with
the inter-quartile range using the same color, but with transparency.


```r
plot(result)
```

![](https://siminf.org/img/mparse-SIR.png)

Most modeling and simulation studies require custom data analysis once
the simulation data has been generated.  To support this, SimInf
provides the `trajectory` method to obtain a `data.frame` with the
number of individuals in each compartment at the time points specified
in `tspan`. Below is the first 10 lines of the `data.frame` with
simulated data.


```r
trajectory(result)
```

```
#>    node time  S I R
#> 1     1    1 98 6 0
#> 2     2    1 98 6 0
#> 3     3    1 98 6 0
#> 4     4    1 99 5 0
#> 5     5    1 97 7 0
#> 6     6    1 98 5 1
#> 7     7    1 99 5 0
#> 8     8    1 99 5 0
#> 9     9    1 97 7 0
#> 10   10    1 97 6 1
...
```

Finally, let us use the `prevalence` method to explore the proportion
of infected individuals across all nodes. It takes a model object and
a formula specification, where the left hand side of the formula
specifies the compartments representing cases i.e., have an attribute
or a disease and the right hand side of the formula specifies the
compartments at risk. Below is the first 10 lines of the `data.frame`.


```r
prevalence(result, I ~ S + I + R)
```

```
#>    time prevalence
#> 1     1 0.05196154
#> 2     2 0.05605769
#> 3     3 0.06059615
#> 4     4 0.06516346
#> 5     5 0.06977885
#> 6     6 0.07390385
#> 7     7 0.07856731
#> 8     8 0.08311538
#> 9     9 0.08794231
#> 10   10 0.09321154
...
```

## Learn more

See the
[vignette](https://CRAN.R-project.org/package=SimInf/vignettes/SimInf.pdf)
to learn more about special features that the SimInf R package
provides, for example, how to:

- use continuous state variables

- use the SimInf framework from another R package

- incorporate available data such as births, deaths and movements as
  scheduled events at predefined time-points.

## Installation

You can install the released version of `SimInf` from
[CRAN](https://CRAN.R-project.org/package=SimInf)


```r
install.packages("SimInf")
```

or use the `remotes` package to install the development version from
[GitHub](https://github.com/stewid/SimInf)


```r
library(remotes)
install_github("stewid/SimInf")
```

We refer to section 3.1 in the
[vignette](https://CRAN.R-project.org/package=SimInf/vignettes/SimInf.pdf)
for detailed installation instructions.

## Authors

In alphabetical order: Pavol Bauer [![ORCID
iD](https://orcid.org/sites/default/files/images/orcid_16x16.gif)](https://orcid.org/0000-0003-4328-7171),
Robin Eriksson [![ORCID
iD](https://orcid.org/sites/default/files/images/orcid_16x16.gif)](https://orcid.org/0000-0002-4291-712X),
Stefan Engblom [![ORCID
iD](https://orcid.org/sites/default/files/images/orcid_16x16.gif)](https://orcid.org/0000-0002-3614-1732),
and Stefan Widgren [![ORCID
iD](https://orcid.org/sites/default/files/images/orcid_16x16.gif)](https://orcid.org/0000-0001-5745-2284)
**(Maintainer)**

Any suggestions, bug reports, forks and pull requests are
appreciated. Get in touch.

## Citation

SimInf is research software. To cite SimInf in publications, please
use:

- Widgren S, Bauer P, Eriksson R, Engblom S (2019) SimInf: An R
  Package for Data-Driven Stochastic Disease Spread Simulations.
  Journal of Statistical Software, 91(12), 1--42. doi:
  [10.18637/jss.v091.i12](https://doi.org/10.18637/jss.v091.i12)

- Bauer P, Engblom S, Widgren S (2016) Fast event-based
  epidemiological simulations on national scales. International
  Journal of High Performance Computing Applications, 30(4),
  438--453. doi:
  [10.1177/1094342016635723](https://doi.org/10.1177/1094342016635723)

## Acknowledgments

This work was financially supported by the Swedish Research Council
within the UPMARC Linnaeus centre of Excellence (Pavol Bauer, Robin
Eriksson and Stefan Engblom), the Swedish Research Council Formas
(Stefan Engblom and Stefan Widgren), the Swedish Board of Agriculture
(Stefan Widgren), and by the Swedish strategic research program
eSSENCE (Stefan Widgren).

## Versioning

The `SimInf` package uses [semantic versioning](https://semver.org/).

## License

The `SimInf` package is licensed under the
[GPLv3](https://github.com/stewid/SimInf/blob/master/LICENSE).
