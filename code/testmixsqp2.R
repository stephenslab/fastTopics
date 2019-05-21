# This is a short script used to verify convergence of the SQP updates
# to the correct solution for a more challenging example.

# SET UP ENVIRONMENT
# ------------------
library(Rcpp)
source("misc.R")
source("mixsqp.R")
sourceCpp("mixsqp.cpp")

# LOAD DATA
# ---------
cat("Loading data set.\n")
load("../data/tacks.RData")
L <- tacks$L
w <- tacks$w

# Initialize the solution estimate.
m  <- ncol(L)
x0 <- rep(1/m,m)

# FIT MODEL USING EM & SQP
# ------------------------
cat("Fitting model by iterating SQP updates.\n")
fit1 <- mixem(L,w,x0,numiter = 4)
fit2 <- mixsqp(L,w,fit1$x,numiter = 24,zero.threshold = 1e-8,verbose = TRUE)

# SUMMARIZE RESULTS
# -----------------
cat(sprintf("Objective at true solution: %0.12f\n",
            mixobjective(L,w,tacks$x,0)))
cat(sprintf("Objective at SQP solution:  %0.12f\n",fit2$value))
cat(sprintf("Largest difference between solutions: %0.1e\n",
            max(abs(tacks$x - fit2$x))))
