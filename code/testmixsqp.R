# This is a short script to verify the implementation of the SQP and
# EM updates for computing maximum-likelihood estimates of the mixture
# proportions in a mixture model.

# SET UP ENVIRONMENT
# ------------------
library(Rcpp)
# sourceCpp("mixsqp.cpp",verbose = TRUE)
source("misc.R")
source("mixsqp.R")

# LOAD DATA
# ---------
cat("Loading data set.\n")
mixdata <- readRDS("../data/mixdata.rds")
L       <- mixdata$L
w       <- mixdata$w

# Initialize the solution estimate.
m  <- ncol(L)
x0 <- rep(1/m,m)

# FIT MODEL USING SQP
# -------------------
cat("Fitting model by iterating SQP updates.\n")
fit1 <- mixsqp(L,w,x0,numiter = 14,verbose = TRUE)

# FIT MODEL USING EM
# ------------------
cat("Fitting model by iterating EM updates.\n")
fit2 <- mixem(L,w,x0,numiter = 1000)

# SUMMARIZE RESULTS
# -----------------
cat(sprintf("Objective at true solution: %0.12f\n",mixdata$value))
cat(sprintf("Objective at SQP solution:  %0.12f\n",fit1$value))
cat(sprintf("Objective at EM solution:   %0.12f\n",fit2$value))
cat(sprintf("Largest difference between true and SQP solutions: %0.1e\n",
            max(abs(mixdata$x - fit1$x))))
cat(sprintf("Largest difference between true and EM solutions: %0.1e\n",
            max(abs(mixdata$x - fit2$x))))
cat(sprintf("EM updates yield improvement between %0.1e and %0.1e.\n",
            min(-diff(fit2$progress$obj)),max(-diff(fit2$progress$obj))))
