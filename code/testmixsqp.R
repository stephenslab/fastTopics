# TO DO: Explain here what this script does, and how to use it.

# SET UP ENVIRONMENT
# ------------------
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

# FIT MODEL USING EM
# ------------------
cat("Fitting model by iterating EM updates.\n")
fit1 <- mixem(L,w,x0,numiter = 1000)

# FIT MODEL USING SQP
# -------------------
cat("Fitting model by iterating SQP updates.\n")
# TO DO.

# SUMMARIZE RESULTS
# -----------------
cat(sprintf("Objective at true solution: %0.12f\n",mixdata$value))
cat(sprintf("Objective at EM solution:   %0.12f\n",fit1$value))
cat(sprintf("Objective at SQP solution:  %0.12f\n",0))
cat(sprintf("Largest difference between true and EM solutions:  %0.1e\n",
            max(abs(mixdata$x - fit1$x))))
cat(sprintf("Largest difference between true and SQP solutions: %0.1e\n",0))
cat(sprintf("Each EM update yields an improvement between %0.1e and %0.1e.\n",
            min(-diff(fit1$progress$obj)),max(-diff(fit1$progress$obj))))
cat(sprintf("Each SQP update yields an improvement between %0.1e and %0.1e.\n",
            0,0))
