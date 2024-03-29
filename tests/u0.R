## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
## Copyright (C) 2015 -- 2022 Stefan Widgren
##
## SimInf is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## SimInf is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <https://www.gnu.org/licenses/>.

library(SimInf)
library(tools)
source("util/check.R")

model <- SIR(u0     = data.frame(S = 1:3, I = 4:6, R = 7:9),
             tspan  = 1:10,
             events = NULL,
             beta   = 0,
             gamma  = 0)

u0(model) <- data.frame(S = 10:12, I = 13:15, R = 16:18)
stopifnot(identical(model@u0,
                    matrix(c(10L, 11L, 12L,
                             13L, 14L, 15L,
                             16L, 17L, 18L),
                           nrow = 3,
                           ncol = 3,
                           byrow = TRUE,
                           dimnames = list(c("S", "I", "R"), NULL))))

stopifnot(identical(u0(model), data.frame(S = 10:12, I = 13:15, R = 16:18)))

res <- assertError(
    u0(model) <- data.frame(S = 10:13, I = 14:17, R = 18:21))
check_error(res, "The number of rows in 'u0' must match nodes in 'model'.")
