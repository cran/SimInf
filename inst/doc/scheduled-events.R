## ----pressure, echo=FALSE, fig.align="left", fig.cap="**Figure 1.** Illustration of movements between nodes. Each time step depicts movements during one time unit, for example, a day. The network has *N=4* nodes where node *1* is infected and nodes *2*--*4* are non-infected. Arrows indicate movements of individuals from a source node to a destination node and labels denote the size of the shipment. Here, infection may spread from node *1* to node *3* at *t=2* and then from node *3* to node *2* at *t=3*.", out.width = '100%'----
knitr::include_graphics("img/temporal-network.svg")

## ---- eval = TRUE, echo = TRUE, message = FALSE---------------------
events <- data.frame(
    event      = rep("extTrans", 6),  ## Event "extTrans" is a movement between nodes
    time       = c(1, 1, 2, 2, 3, 3), ## The time that the event happens
    node       = c(3, 3, 1, 4, 3, 4), ## In which node does the event occur
    dest       = c(4, 2, 3, 3, 2, 2), ## Which node is the destination node
    n          = c(9, 2, 8, 3, 5, 4), ## How many individuals are moved
    proportion = c(0, 0, 0, 0, 0, 0), ## This is not used when n > 0
    select     = c(4, 4, 4, 4, 4, 4), ## Use the 4th column in the model select matrix
    shift      = c(0, 0, 0, 0, 0, 0)) ## Not used in this example

## -------------------------------------------------------------------
events

## -------------------------------------------------------------------
library(SimInf)

model <- SIR(u0 = data.frame(S = c(10, 15, 20, 25),
                             I = c( 5,  0,  0,  0),
                             R = c( 0,  0,  0,  0)),
             tspan = 0:3,
             beta = 0,
             gamma = 0,
             events = events)

## -------------------------------------------------------------------
model@events@E

## -------------------------------------------------------------------
set.seed(1)
result <- run(model)

## ---- fig.width=7, fig.height=4, fig.align="left", fig.cap="**Figure 2.** Number of susceptible, infected and recovered individuals in each node."----
plot(result, range = FALSE)

## -------------------------------------------------------------------
trajectory(result)

