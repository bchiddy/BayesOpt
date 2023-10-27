# Analysis of benchmark algorithms across chosen test functions
# Author: Ben Chiddy
# Date: 24/10/2023
# Reference: (R Core Team, 2023)

## Scaled Goldstein-Price

## Analysis variables
ninit       <- 12
end         <- 50
pmin        <- 0
pmax        <- 1
num_par     <- 2
n           <- 50
runs        <- 1000
true_min    <- -3.12
noise       <- TRUE
f1          <- goldprsc
f2          <- goldprsc_opt

## Nelder-Mead
params_data <- data.frame()
for(r in 1:runs) {
  y <- c()
  xx1 <- c()
  xx2 <-c()
  os <- optim(runif(num_par, pmin, pmax), 
              f=f2, 
              method="Nelder-Mead") 
  params_data <- rbind(params_data, data.frame(
    x1 = xx1[1:50],
    x2 = xx2[1:50],
    y = y[1:50],
    bov = bov(y, n),
    run = r
  ))
}
write.table(params_data, file="data/goldstein_price/params_nlm.csv", 
            row.names = FALSE, col.names = FALSE, sep=",")

## L-BFGS-B
params_data <- data.frame()
for(r in 1:runs) {
  y <- c()
  xx1 <- c()
  xx2 <-c()
  os <- optim(runif(num_par, pmin, pmax), 
              f=f2, 
              lower=pmin, 
              upper=pmax, 
              method="L-BFGS-B") 
  params_data <- rbind(params_data, data.frame(
    x1 = xx1[1:50],
    x2 = xx2[1:50],
    y = y[1:50],
    bov = bov(y, n),
    run = r
  ))
}
write.table(params_data, file="data/goldstein_price/params_lfbgs.csv", 
            row.names = FALSE, col.names = FALSE, sep=",")

## SANN
params_data <- data.frame()
for(r in 1:runs) {
  y <- c()
  xx1 <- c()
  xx2 <-c()
  os <- optim(runif(num_par, pmin, pmax), 
              f=f2, 
              method="SANN") 
  params_data <- rbind(params_data, data.frame(
    x1 = xx1[1:50],
    x2 = xx2[1:50],
    y = y[1:50],
    bov = bov(y, n),
    run = r
  ))
}
write.table(params_data, file="data/goldstein_price/params_SANN.csv", 
            row.names = FALSE, col.names = FALSE, sep=",")

## RS
params_data <- data.frame()
for(r in 1:runs) {
  vals <- random_search(f=f1, n=50, names=c("x1","x2"), ninit = 12,
                        start_vals = initialize_values(2,pmin,pmax,12),
                        pmin=pmin, pmax=pmax)
  params_data <- rbind(params_data, data.frame(
    x1 = vals[,1],
    x2 = vals[,2],
    y = vals[,3],
    bov = bov(y, n),
    run = r
  ))
}
write.table(params_data, file="data/goldstein_price/params_RS.csv", 
            row.names = FALSE, col.names = FALSE, sep=",")

## Ackley

## Analysis variables
ninit       <- 12
end         <- 50
pmin        <- -5
pmax        <- 5
num_par     <- 2
n           <- 50
runs        <- 1000
true_min    <- 0
noise       <- FALSE
f1          <- ackley
f2          <- ackley_opt

## Nelder-Mead
params_data <- data.frame()
for(r in 1:runs) {
  y <- c()
  xx1 <- c()
  xx2 <-c()
  os <- optim(runif(num_par, pmin, pmax), 
              f=f2, 
              method="Nelder-Mead") 
  params_data <- rbind(params_data, data.frame(
    x1 = xx1[1:50],
    x2 = xx2[1:50],
    y = y[1:50],
    bov = bov(y, n),
    run = r
  ))
}
write.table(params_data, file="data/ackley/params_nlm.csv", 
            row.names = FALSE, col.names = FALSE, sep=",")

## L-BFGS-B
params_data <- data.frame()
for(r in 1:runs) {
  y <- c()
  xx1 <- c()
  xx2 <-c()
  os <- optim(runif(num_par, pmin, pmax), 
              f=f2, 
              lower=pmin, 
              upper=pmax, 
              method="L-BFGS-B") 
  params_data <- rbind(params_data, data.frame(
    x1 = xx1[1:50],
    x2 = xx2[1:50],
    y = y[1:50],
    bov = bov(y, n),
    run = r
  ))
}
write.table(params_data, file="data/ackley/params_lfbgs.csv", 
            row.names = FALSE, col.names = FALSE, sep=",")

## SANN
params_data <- data.frame()
for(r in 1:runs) {
  y <- c()
  xx1 <- c()
  xx2 <-c()
  os <- optim(runif(num_par, pmin, pmax), 
              f=f2, 
              method="SANN") 
  params_data <- rbind(params_data, data.frame(
    x1 = xx1[1:50],
    x2 = xx2[1:50],
    y = y[1:50],
    bov = bov(y, n),
    run = r
  ))
}
write.table(params_data, file="data/ackley/params_SANN.csv", 
            row.names = FALSE, col.names = FALSE, sep=",")

# RS
params_data <- data.frame()
for(r in 1:runs) {
  vals <- random_search(f=f1, n=50, names=c("x1","x2"), ninit = 12,
                        start_vals = initialize_values(2,pmin,pmax,12),
                        pmin=pmin, pmax=pmax)
  params_data <- rbind(params_data, data.frame(
    x1 = vals[,1],
    x2 = vals[,2],
    y = vals[,3],
    bov = bov(y, n),
    run = r
  ))
}
write.table(params_data, file="data/ackley/params_RS.csv", 
            row.names = FALSE, col.names = FALSE, sep=",")

## Bukin number 6

## Analysis variables
ninit       <- 12
end         <- 50
pmin        <- -5
pmax        <- 5
num_par     <- 2
n           <- 50
runs        <- 1000
true_min    <- 0
noise       <- FALSE
f1          <- bukin6
f2          <- bukin6_opt

## Nelder-Mead
params_data <- data.frame()
for(r in 1:runs) {
  y <- c()
  xx1 <- c()
  xx2 <-c()
  os <- optim(runif(num_par, pmin, pmax), 
              f=f2, 
              method="Nelder-Mead") 
  params_data <- rbind(params_data, data.frame(
    x1 = xx1[1:50],
    x2 = xx2[1:50],
    y = y[1:50],
    bov = bov(y, n),
    run = r
  ))
}
write.table(params_data, file="data/bukin6/params_nlm.csv", 
            row.names = FALSE, col.names = FALSE, sep=",")

## L-BFGS-B
params_data <- data.frame()
for(r in 1:runs) {
  y <- c()
  xx1 <- c()
  xx2 <-c()
  os <- optim(runif(num_par, pmin, pmax), 
              f=f2, 
              lower=pmin, 
              upper=pmax, 
              method="L-BFGS-B") 
  params_data <- rbind(params_data, data.frame(
    x1 = xx1[1:50],
    x2 = xx2[1:50],
    y = y[1:50],
    bov = bov(y, n),
    run = r
  ))
}
write.table(params_data, file="data/bukin6/params_lfbgs.csv", 
            row.names = FALSE, col.names = FALSE, sep=",")

## SANN
params_data <- data.frame()
for(r in 1:runs) {
  y <- c()
  xx1 <- c()
  xx2 <-c()
  os <- optim(runif(num_par, pmin, pmax), 
              f=f2, 
              method="SANN") 
  params_data <- rbind(params_data, data.frame(
    x1 = xx1[1:50],
    x2 = xx2[1:50],
    y = y[1:50],
    bov = bov(y, n),
    run = r
  ))
}
write.table(params_data, file="data/bukin6/params_SANN.csv", 
            row.names = FALSE, col.names = FALSE, sep=",")

## RS
params_data <- data.frame()
for(r in 1:runs) {
  vals <- random_search(f=f1, n=50, names=c("x1","x2"), ninit = 12,
                        start_vals = initialize_values(2,pmin,pmax,12),
                        pmin=pmin, pmax=pmax)
  params_data <- rbind(params_data, data.frame(
    x1 = vals[,1],
    x2 = vals[,2],
    y = vals[,3],
    bov = bov(y, n),
    run = r
  ))
}
write.table(params_data, file="data/bukin6/params_RS.csv", 
            row.names = FALSE, col.names = FALSE, sep=",")

