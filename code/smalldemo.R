# Small script illustrating application of the multiplicative updates
# and the alternating SQP method for fitting a Poisson topic model.

# SCRIPT PARAMETERS
# -----------------
# Number of factors (topics).
k <- 13  

# SET UP ENVIRONMENT
# ------------------
library(parallel)
library(daarem)
library(Rcpp)
library(readr)
library(ggplot2)
library(cowplot)
source("misc.R")
source("betanmf.R")
source("mixsqp.R")
source("altsqp.R")
sourceCpp("mixsqp.cpp")

# Initialize the sequence of pseudorandom numbers.
set.seed(1)

# LOAD DATA
# ---------
# I remove columns in which all the counts are zero.
cat("Loading count data.\n")
counts        <- suppressMessages(read_csv("../data/droplet_small.csv.gz"))
class(counts) <- "data.frame"
counts        <- as.matrix(counts)
n             <- nrow(counts)
m             <- ncol(counts)
cat(sprintf("Loaded %d x %d counts matrix.\n",n,m))

# GENERATE INITIAL ESTIMATES
# --------------------------
# Generate initial estimates of the factors (stored as an m x k
# matrix) and loadings (stored as an n x k matrix).
F <- matrix(runif(m*k),m,k)
L <- matrix(runif(n*k),n,k)

# RUN MULTIPLICATIVE UPDATES
# --------------------------
cat("Fitting Poisson topic model by iterating multiplicative updates.\n")
fit1 <- betanmf(counts,L,t(F),numiter = 100)

# RUN ALTERNATING SQP METHOD
# --------------------------
cat("Fitting Poisson topic model by iterating SQP updates.\n")
fit2 <- altsqp(counts,F,L,numiter = 100,method = "accelerated",
               control = list(nc = 2,order = 10))

# PLOT IMPROVEMENT IN SOLUTIONS OVER TIME
# ---------------------------------------
bestf <- -251269.912745
pdat <- rbind(cbind(data.frame(iter      = 1:100,objective = fit1$value,
                               method    = "betanmf")),
              cbind(data.frame(iter      = 1:100,
                               objective = fit2$value,
                               method    = "altsqp")))
p1 <- ggplot(pdat,aes(x = iter,y = objective - bestf,color = method)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("darkorange","darkblue")) +
  scale_y_continuous(breaks = 10^seq(0,6,0.5),trans = "log10") +
  labs(x = "iteration",y = "distance from minimum")
print(p1)
