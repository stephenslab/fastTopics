# TO DO: Explain here what this function does, and how to use it.

# SCRIPT PARAMETERS
# -----------------

# Number of factors (topics).
K <- 13

# SET UP ENVIRONMENT
# ------------------
library(readr)
source("../code/misc.R")
source("../code/betanmf.R")

# Initialize the sequence of pseudorandom numbers.
set.seed(1)

# LOAD DATA
# ---------
cat("Loading count data.\n")
counts        <- suppressMessages(read_csv("../data/droplet.csv.gz"))
class(counts) <- "data.frame"
counts        <- as.matrix(counts)
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
fit.betanmf <- betanmf(counts,L,t(F),numiter = 200)
    
# RUN ALTERNATING SQP METHOD
# --------------------------
cat("Fitting Poisson topic model using alternating SQP method.\n")
# TO DO.
