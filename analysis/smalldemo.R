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
cat(sprintf("Loaded %d x %d counts matrix.\n",nrow(counts),ncol(counts)))

# GENERATE INITIAL ESTIMATES
# --------------------------

# RUN MULTIPLICATIVE UPDATES
# --------------------------
cat("Fitting Poisson topic model using multiplicative updates.\n")
# TO DO.

# RUN ALTERNATING SQP METHOD
# --------------------------
cat("Fitting Poisson topic model using alternating SQP method.\n")
