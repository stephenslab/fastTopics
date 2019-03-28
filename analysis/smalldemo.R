# Small script illustrating application of the betanmf and altsqp
# methods for fitting a Poisson topic model.

# SCRIPT PARAMETERS
# -----------------
K  <- 13  # Number of factors (topics).
nc <- 2   # Number of threads to use.

# SET UP ENVIRONMENT
# ------------------
library(parallel)
library(Matrix)
library(quadprog)
library(readr)
library(ggplot2)
library(cowplot)
source("../code/misc.R")
source("../code/altsqp.R")
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
cat("Fitting Poisson topic model using multiplicative updates.\n")
fit.betanmf <- betanmf(counts,L,t(F),numiter = 100)

# RUN ALTERNATING SQP METHOD
# --------------------------
cat("Fitting Poisson topic model using alternating SQP method.\n")
fit.altsqp <- altsqp(counts,F,L,nc = 2,numiter = 10)
    
# PLOT IMPROVEMENT IN SOLUTIONS OVER TIME
# ---------------------------------------
bestf <- -250300.7483213344530668
p1 <- ggplot(fit.betanmf$progress,
             aes(x = iter,y = objective - bestf + 1e-8)) +
  geom_line(color = "darkblue",size = 1) +
  geom_line(color = "darkorange",data = fit.altsqp$progress,size = 1) +
  scale_y_continuous(trans = "log10") +
  labs(x = "iteration",y = "distance from minimum")

