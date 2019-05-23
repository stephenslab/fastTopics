# Small script illustrating application of the multiplicative updates
# and the alternating SQP method for fitting a Poisson topic model.

# SCRIPT PARAMETERS
# -----------------
# Number of factors (topics).
K <- 13  

# SET UP ENVIRONMENT
# ------------------
library(Matrix)
library(readr)
library(ggplot2)
library(cowplot)
source("../code/misc.R")
source("../code/betanmf.R")

# Initialize the sequence of pseudorandom numbers.
set.seed(1)

# LOAD DATA
# ---------
cat("Loading count data.\n")
counts        <- suppressMessages(read_csv("../data/droplet_small.csv.gz"))
class(counts) <- "data.frame"
counts        <- as.matrix(counts)
counts        <- Matrix(counts,sparse = TRUE)
n             <- nrow(counts)
p             <- ncol(counts)
cat(sprintf("Loaded %d x %d counts matrix.\n",n,p))

# GENERATE INITIAL ESTIMATES
# --------------------------
# Generate initial estimates of the factors (stored as an p x K
# matrix) and loadings (stored as an n x K matrix).
F <- matrix(runif(p*K),p,K)
L <- matrix(runif(n*K),n,K)

# RUN MULTIPLICATIVE UPDATES
# --------------------------
cat("Fitting Poisson topic model by iterating multiplicative updates.\n")
fit.betanmf <- betanmf(counts,L,t(F),numiter = 100)

# RUN ALTERNATING SQP METHOD
# --------------------------
cat("Fitting Poisson topic model by iterating SQP updates.\n")
# TO DO: Revise this line of code.
fit.altsqp <- altsqp(counts,F,L,numiter = 50)
    
# PLOT IMPROVEMENT IN SOLUTIONS OVER TIME
# ---------------------------------------
bestf <- -250300.7483213344530668
p1 <- ggplot(fit.betanmf$progress,
             aes(x = iter,y = objective - bestf + 1e-8)) +
  geom_line(color = "darkblue",size = 1) +
  geom_line(color = "darkorange",data = fit.altsqp$progress,size = 1) +
  scale_y_continuous(trans = "log10") +
  labs(x = "iteration",y = "distance from minimum")

