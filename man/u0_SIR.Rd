% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SIR.R
\name{u0_SIR}
\alias{u0_SIR}
\title{Example data to initialize the \sQuote{SIR} model}
\usage{
u0_SIR()
}
\value{
A \code{data.frame}
}
\description{
Example data to initialize a population of 1600 nodes and
demonstrate the \code{\linkS4class{SIR}} model.
}
\details{
A \code{data.frame} with the number of individuals in the
\sQuote{S}, \sQuote{I} and \sQuote{R} compartments in 1600
nodes. Note that the \sQuote{I} and \sQuote{R} compartments are
zero.
}
\examples{
\dontrun{
## For reproducibility, call the set.seed() function and specify
## the number of threads to use. To use all available threads,
## remove the set_num_threads() call.
set.seed(123)
set_num_threads(1)

## Create an 'SIR' model with 1600 nodes and initialize
## it to run over 4*365 days. Add one infected individual
## to the first node.
u0 <- u0_SIR()
u0$I[1] <- 1
tspan <- seq(from = 1, to = 4*365, by = 1)
model <- SIR(u0     = u0,
             tspan  = tspan,
             events = events_SIR(),
             beta   = 0.16,
             gamma  = 0.01)

## Run the model to generate a single stochastic trajectory.
result <- run(model)
plot(result)

## Summarize trajectory
summary(result)
}
}
