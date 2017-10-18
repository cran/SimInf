### R code from vignette source 'SimInf.Rnw'

###################################################
### code chunk number 1: SimInf.Rnw:100-101
###################################################
options(prompt = "R> ", continue = "+    ", width = 74, useFancyQuotes = FALSE)


###################################################
### code chunk number 2: load-SimInf
###################################################
library("SimInf")


###################################################
### code chunk number 3: SIR-u0
###################################################
u0 <- data.frame(S = 99, I = 1, R = 0)


###################################################
### code chunk number 4: SIR-tspan
###################################################
tspan <- seq(1, 6 * 30, by = 7)


###################################################
### code chunk number 5: SIR-model
###################################################
model <- SIR(u0 = u0, tspan = tspan, beta = 0.16, gamma = 0.077)


###################################################
### code chunk number 6: SIR-run
###################################################
result <- run(model, threads = 1, seed = 22)


###################################################
### code chunk number 7: SIR-U
###################################################
U(result)[, 1:10]


###################################################
### code chunk number 8: SimInf.Rnw:1033-1034 (eval = FALSE)
###################################################
## plot(result)


###################################################
### code chunk number 9: SIR-plot-proportion
###################################################
plot(result)


###################################################
### code chunk number 10: SIR-multiple-nodes
###################################################
u0 <- data.frame(S = c(99, 90), I = c( 1,  2), R = c( 0,  0))
model <- SIR(u0 = u0, tspan = tspan, beta = 0.16, gamma = 0.077)
result <- run(model, threads = 1, seed = 22)


###################################################
### code chunk number 11: SIR-U-multiple-nodes
###################################################
U(result)[, 1:10]


###################################################
### code chunk number 12: SIR-run-outer
###################################################
u0 <- data.frame(S = rep(99, 500), I = rep(1, 500), R = rep(0, 500))
model <- SIR(u0, 1:75, beta = 0.16, gamma = 0.077)
x <- seq(from = 0.2, to = 1.8, by = 0.05)
y <- seq(from = 0.2, to = 1.1, by = 0.05)
pop <- run_outer(x, y, model, gamma ~ beta,
function(model) {prevalence(run(model), "pop")[75]})
bnp <- run_outer(x, y, model, gamma ~ beta,
function(model) {prevalence(run(model), "bnp")[75]})


###################################################
### code chunk number 13: SIR-run-outer-plot
###################################################
opar <- par(mfrow = c(1, 2), oma = c(1, 1, 3, 0), mar = c(4, 3, 1, 1))
contour(x * model@gdata["beta"], y * model@gdata["gamma"],
        pop, method = "edge", bty = "l")
title("Population prevalence")
mtext(expression(beta), side = 1, line = 3)
mtext(expression(gamma), side = 2, line = 2)
contour(x * model@gdata["beta"], y * model@gdata["gamma"],
bnp, method = "edge", bty = "l")
title("Between-node prevalence")
mtext(expression(beta), side = 1, line = 3)
mtext(expression(gamma), side = 2, line = 2)
par(opar)


###################################################
### code chunk number 14: SIR-u0-multiple-nodes
###################################################
u0 <- data.frame(S = c(10, 0, 0, 0, 0), I = c(10, 0, 0, 0, 0),
R = c(10, 0, 0, 0, 0))


###################################################
### code chunk number 15: SIR-events
###################################################
events <- structure(list(
    event      = c(1, 3, 0, 3, 3, 3),
    time       = c(2, 3, 4, 4, 4, 5),
    node       = c(2, 1, 2, 1, 1, 1),
    dest       = c(0, 3, 0, 3, 4, 5),
    n          = c(2, 1, 1, 1, 1, 2),
    proportion = c(0, 0, 0, 0, 0, 0),
    select     = c(1, 2, 2, 2, 2, 2),
    shift      = c(0, 0, 0, 0, 0, 0)),
    .Names = c("event", "time", "node", "dest", "n",
               "proportion", "select", "shift"),
    row.names = c(NA, -6L), class = "data.frame")


###################################################
### code chunk number 16: SIR-events-show
###################################################
events


###################################################
### code chunk number 17: U-model-events
###################################################
model <- SIR(u0, 1:5, events = events, beta = 0.16, gamma = 0.077)
U(run(model, threads = 1, seed = 22))


###################################################
### code chunk number 18: SIR-u0-multiple-nodes-II
###################################################
u0 <- data.frame(S = c(90, 100, 100, 100, 100), I = c(10, 0, 0, 0, 0),
R = c(0, 0, 0, 0, 0))


###################################################
### code chunk number 19: SIR-events-II
###################################################
events <- structure(list(
    event      = c(3, 3, 3, 3),
    time       = c(3, 25, 5, 10),
    node       = c(1, 1, 1, 1),
    dest       = c(3, 2, 3, 3),
    n          = c(1, 9, 1, 5),
    proportion = c(0, 0, 0, 0),
    select     = c(2, 2, 2, 2),
    shift      = c(0, 0, 0, 0)),
                    .Names = c("event", "time", "node", "dest",
                        "n", "proportion", "select", "shift"),
    row.names = c(NA, -4L), class = "data.frame")


###################################################
### code chunk number 20: SIR-events-II-show
###################################################
events


###################################################
### code chunk number 21: SIR-events-II-run
###################################################
model <- SIR(u0, 1:75, events = events, beta = 0.16, gamma = 0.077)
result <- run(model, threads = 1, seed = 5)


###################################################
### code chunk number 22: SIR-spaghetti (eval = FALSE)
###################################################
## plot(result, N = TRUE, spaghetti = TRUE, compartments = "I")


###################################################
### code chunk number 23: SIR-spaghetti-plot
###################################################
plot(result, N = TRUE, spaghetti = TRUE, compartments = "I")


###################################################
### code chunk number 24: load-SISe3sp-events
###################################################
data("u0_SISe3", package = "SimInf")
data("events_SISe3", package = "SimInf")


###################################################
### code chunk number 25: load-nodes
###################################################
data("nodes", package = "SimInf")
d <- distance_matrix(x = nodes$x, y = nodes$y, cutoff = 2500)


###################################################
### code chunk number 26: create-SISe3sp-tspan
###################################################
tspan <- seq(1, 4 * 365, by = 7)


###################################################
### code chunk number 27: create-SISe3sp-u0
###################################################
set.seed(123)
i <- sample(1:1600, 160)
u0_SISe3$I_1[i] <- as.integer(u0_SISe3$S_1[i] * 0.05)
u0_SISe3$I_2[i] <- as.integer(u0_SISe3$S_2[i] * 0.05)
u0_SISe3$I_3[i] <- as.integer(u0_SISe3$S_3[i] * 0.05)
u0_SISe3$S_1[i] <- u0_SISe3$S_1[i] - u0_SISe3$I_1[i]
u0_SISe3$S_2[i] <- u0_SISe3$S_2[i] - u0_SISe3$I_2[i]
u0_SISe3$S_3[i] <- u0_SISe3$S_3[i] - u0_SISe3$I_3[i]


###################################################
### code chunk number 28: create-SISe3sp
###################################################
model <- SISe3_sp(u0 = u0_SISe3, tspan = tspan, phi = rep(0, 1600),
events = events_SISe3, upsilon_1 = 0.013, upsilon_2 = 0.013, upsilon_3
= 0.013, gamma_1 = 0.1, gamma_2 = 0.1, gamma_3 = 0.1, alpha = 1,
beta_t1 = 0.095, beta_t2 = 0.12, beta_t3 = 0.10, beta_t4 = 0.15,
end_t1 = 91, end_t2 = 182, end_t3 = 273, end_t4 = 365, distance = d,
coupling = 0.1)


###################################################
### code chunk number 29: SISe3sp-model-summary
###################################################
summary(model)


###################################################
### code chunk number 30: display-SISe3sp-model
###################################################
plot(events(model))


###################################################
### code chunk number 31: SimInf.Rnw:1563-1573 (eval = FALSE)
###################################################
## result <- run(model, threads = 1, seed = 1)
## N1 <- colSums(susceptible(result, age = 1) + infected(result, age = 1))
## N2 <- colSums(susceptible(result, age = 2) + infected(result, age = 2))
## N3 <- colSums(susceptible(result, age = 3) + infected(result, age = 3))
## opar <- par(mfrow = c(3, 1), mar = c(2, 5, 1, 1))
## plot(N3 ~ tspan, ylab = "N", type = "l", bty = "l", lty = 1)
## plot(N2 ~ tspan, ylab = "N", type = "l", bty = "l", lty = 2)
## plot(N1 ~ tspan, ylab = "N", type = "l", bty = "l", lty = 3)
## mtext("Day", side = 1, line = 1, at = 750)
## par(opar)


###################################################
### code chunk number 32: population-size
###################################################
result <- run(model, threads = 1, seed = 1)
N1 <- colSums(susceptible(result, age = 1) + infected(result, age = 1))
N2 <- colSums(susceptible(result, age = 2) + infected(result, age = 2))
N3 <- colSums(susceptible(result, age = 3) + infected(result, age = 3))
opar <- par(mfrow = c(3, 1), mar = c(2, 5, 1, 1))
plot(N3 ~ tspan, ylab = "N", type = "l", bty = "l", lty = 1)
plot(N2 ~ tspan, ylab = "N", type = "l", bty = "l", lty = 2)
plot(N1 ~ tspan, ylab = "N", type = "l", bty = "l", lty = 3)
mtext("Day", side = 1, line = 1, at = 750)
par(opar)


###################################################
### code chunk number 33: SimInf.Rnw:1609-1610
###################################################
options(continue = "+  ")


###################################################
### code chunk number 34: plot-infected-nodes
###################################################
plot_nodes <- function(x) {
  wnp <- prevalence(x, type = "wnp") > 0
  opar <- par(mfrow = c(3, 6), mar = c(2, 0.3, 0.3, 0.3))
  on.exit(par(opar))
  for(i in seq(1, dim(wnp)[2], length.out = 18)) {
    Si <- which(!wnp[, i])
    Ii <- which(wnp[, i])
    plot(nodes$x[Si], nodes$y[Si], col = "yellow", pch = 20, xlab = "",
      ylab = "", xaxt = "n", yaxt = "n", ann = FALSE)
    day <- 7 * as.integer(i - 1) + 1
    mtext(sprintf("Day: %i", day), side = 1, line = 0.5)
    points(nodes$x[Ii], nodes$y[Ii], col = "blue", pch = 20)
  }
}


###################################################
### code chunk number 35: SimInf.Rnw:1630-1631
###################################################
options(continue = "+    ")


###################################################
### code chunk number 36: SimInf.Rnw:1638-1640 (eval = FALSE)
###################################################
## model@gdata["coupling"] <- 0
## plot_nodes(run(model, threads = 1, seed = 1))


###################################################
### code chunk number 37: plot-SISe3sp-I
###################################################
png("SimInf-plot-SISe3sp-I.png", height = 480, width = 960)
plot_nodes(run(model, threads = 1, seed = 1))
dev.off()


###################################################
### code chunk number 38: plot-SISe3sp-II
###################################################
model@gdata["coupling"] <- 0
png("SimInf-plot-SISe3sp-II.png", height = 480, width = 960)
plot_nodes(run(model, threads = 1, seed = 1))
dev.off()


###################################################
### code chunk number 39: plot-SISe3sp-III
###################################################
model@gdata["coupling"] <- 0.1
model@events@N <- cbind(model@events@N,
"3" = c(0L, -1L, 0L, -1L, 0L, -1L))
i <- which(model@events@event == 3)
model@events@shift[i] <- 3L
png("SimInf-plot-SISe3sp-III.png", height = 480, width = 960)
plot_nodes(run(model, threads = 1, seed = 1))
dev.off()


###################################################
### code chunk number 40: SimInf.Rnw:1713-1719 (eval = FALSE)
###################################################
## model@gdata["coupling"] <- 0.1
## model@events@N <- cbind(model@events@N,
## "3" = c(0L, -1L, 0L, -1L, 0L, -1L))
## i <- which(model@events@event == 3)
## model@events@shift[i] <- 3L
## plot_nodes(run(model, threads = 1, seed = 1))


###################################################
### code chunk number 41: SIR-mparse-I
###################################################
transitions <- c("S -> b*S*I/(S+I+R) -> I", "I -> g*I -> R")
compartments <- c("S", "I", "R")


###################################################
### code chunk number 42: SIR-mparse-II
###################################################
m <- mparse(transitions, compartments, b = 0.16, g = 0.077)
model <- init(m, u0 = cbind(S = 99, I = 1, R = 0), tspan = 1:180)


###################################################
### code chunk number 43: SIR-mparse-III (eval = FALSE)
###################################################
## result <- run(model, threads = 1, seed = 22)
## plot(result)


###################################################
### code chunk number 44: SIR-mparse-IV
###################################################
result <- run(model, threads = 1, seed = 22)
plot(result)


###################################################
### code chunk number 45: SIR-mparse-incidence (eval = FALSE)
###################################################
## m <- mparse(c("S -> b*S*I/(S+I+R) -> I + Icum", "I -> g*I -> R"),
## c("S", "I", "Icum", "R"), b = 0.16, g = 0.077)
## model <- init(m, cbind(S = 99, I = 1, Icum = 0, R = 0), 1:180)
## result <- run(model, threads = 1, seed = 22)
## plot(stepfun(result@tspan[-1], diff(c(0, U(result)["Icum",]))),
## main = "", xlab = "Time", ylab = "Number of cases",
## do.points = FALSE)


###################################################
### code chunk number 46: SIR-mparse-incidence-plot
###################################################
m <- mparse(c("S -> b*S*I/(S+I+R) -> I + Icum", "I -> g*I -> R"),
c("S", "I", "Icum", "R"), b = 0.16, g = 0.077)
model <- init(m, cbind(S = 99, I = 1, Icum = 0, R = 0), 1:180)
result <- run(model, threads = 1, seed = 22)
plot(stepfun(result@tspan[-1], diff(c(0, U(result)["Icum",]))),
main = "", xlab = "Time", ylab = "Number of cases",
do.points = FALSE)


###################################################
### code chunk number 47: mparse-predator-prey (eval = FALSE)
###################################################
## m <- mparse(transitions = c("@ -> bR*R -> R",
## "R -> (dR+(bR-dR)*R/K)*R -> @", "R -> alpha/(1+w*R)*R*F -> @",
## "@ -> bF*alpha/(1+w*R)*R*F -> F", "F -> dF*F -> @"),
## compartments = c("R", "F"), bR = 2, bF = 2, dR = 1, K = 1000,
## alpha = 0.007, w = 0.0035, dF = 2)
## model <- init(m, cbind(R = 1000, F = 100), 1:100)
## result <- run(model, threads = 1, seed = 1)
## plot(result, N = TRUE)


###################################################
### code chunk number 48: mparse-predator-prey-plot
###################################################
m <- mparse(transitions = c("@ -> bR*R -> R",
"R -> (dR+(bR-dR)*R/K)*R -> @", "R -> alpha/(1+w*R)*R*F -> @",
"@ -> bF*alpha/(1+w*R)*R*F -> F", "F -> dF*F -> @"),
compartments = c("R", "F"), bR = 2, bF = 2, dR = 1, K = 1000,
alpha = 0.007, w = 0.0035, dF = 2)
model <- init(m, cbind(R = 1000, F = 100), 1:100)
result <- run(model, threads = 1, seed = 1)
plot(result, N = TRUE)


###################################################
### code chunk number 49: create-package-skeleton (eval = FALSE)
###################################################
## path = tempdir()
## package_skeleton(m, "PredatorPrey", path = path)


###################################################
### code chunk number 50: install-package-skeleton (eval = FALSE)
###################################################
## install.packages(file.path(path, "PredatorPrey"), repos = NULL)


###################################################
### code chunk number 51: load-package-skeleton (eval = FALSE)
###################################################
## library("PredatorPrey")


