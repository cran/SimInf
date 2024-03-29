% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SISe.R
\name{u0_SISe}
\alias{u0_SISe}
\title{Example data to initialize the \sQuote{SISe} model}
\usage{
u0_SISe()
}
\value{
A \code{data.frame}
}
\description{
Example data to initialize a population of 1600 nodes and
demonstrate the \code{\linkS4class{SISe}} model.
}
\details{
A \code{data.frame} with the number of individuals in the
\sQuote{S} and \sQuote{I} compartments in 1600 nodes. Note that
the \sQuote{I} compartment is zero.
}
\examples{
\dontrun{
## For reproducibility, call the set.seed() function and specify
## the number of threads to use. To use all available threads,
## remove the set_num_threads() call.
set.seed(123)
set_num_threads(1)

## Create an 'SISe' model with 1600 nodes and initialize it to
## run over 4*365 days and record data at weekly time-points.

## Load the initial population and add ten infected individuals to
## the first node.
u0 <- u0_SISe()
u0$I[1] <- 10

## Define 'tspan' to run the simulation over 4*365 and record the
## state of the system at weekly time-points.
tspan <- seq(from = 1, to = 4*365, by = 7)

## Load scheduled events for the population of nodes with births,
## deaths and between-node movements of individuals.
events <- events_SISe()

## Create an 'SISe' model
model <- SISe(u0 = u0, tspan = tspan, events = events_SISe(),
              phi = 0, upsilon = 1.8e-2, gamma = 0.1, alpha = 1,
              beta_t1 = 1.0e-1, beta_t2 = 1.0e-1, beta_t3 = 1.25e-1,
              beta_t4 = 1.25e-1, end_t1 = 91, end_t2 = 182,
              end_t3 = 273, end_t4 = 365, epsilon = 0)

## Run the model to generate a single stochastic trajectory.
result <- run(model)

## Summarize trajectory
summary(result)

## Plot the proportion of nodes with at least one infected
## individual.
plot(result, I~S+I, level = 2, type = "l")
}
}
