% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SISe3.R
\docType{data}
\name{u0_SISe3}
\alias{u0_SISe3}
\title{Example data to initialize the \sQuote{SISe3} model}
\format{
A \code{data.frame}
}
\usage{
data(u0_SISe3)
}
\description{
Example data to initialize a population of 1600 nodes and
demonstrate the \code{\linkS4class{SISe3}} model.
}
\details{
A \code{data.frame} with the number of individuals in the
\sQuote{S_1}, \sQuote{S_2}, \sQuote{S_3}, \sQuote{I_1},
\sQuote{I_2} and \sQuote{I_3} compartments in 1600 nodes. Note
that the \sQuote{I_1}, \sQuote{I_2} and \sQuote{I_3} compartments
are zero.
}
\examples{
\dontrun{
## For reproducibility, call the set.seed() function and specify
## the number of threads to use. To use all available threads,
## remove the set_num_threads() call.
set.seed(123)
set_num_threads(1)

## Create an 'SISe3' model with 1600 nodes and initialize it to
## run over 4*365 days and record data at weekly time-points.

## Load the initial population and add ten infected individuals to
## I_1 in the first node.
u0 <- u0_SISe3
u0$I_1[1] <- 10

## Define 'tspan' to run the simulation over 4*365 and record the
## state of the system at weekly time-points.
tspan <- seq(from = 1, to = 4*365, by = 7)

## Load scheduled events for the population of nodes with births,
## deaths and between-node movements of individuals.
events <- events_SISe3

## Create a 'SISe3' model
model <- SISe3(u0 = u0, tspan = tspan, events = events,
               phi = rep(0, nrow(u0)), upsilon_1 = 1.8e-2,
               upsilon_2 = 1.8e-2, upsilon_3 = 1.8e-2,
               gamma_1 = 0.1, gamma_2 = 0.1, gamma_3 = 0.1,
               alpha = 1, beta_t1 = 1.0e-1, beta_t2 = 1.0e-1,
               beta_t3 = 1.25e-1, beta_t4 = 1.25e-1, end_t1 = 91,
               end_t2 = 182, end_t3 = 273, end_t4 = 365, epsilon = 0)

## Run the model to generate a single stochastic trajectory.
result <- run(model)

## Summarize trajectory
summary(result)

## Plot the proportion of nodes with at least one infected
## individual.
plot(result, I_1 + I_2 + I_3 ~ ., level = 2, type = "l")
}
}
\keyword{dataset}
