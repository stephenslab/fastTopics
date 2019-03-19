# TO DO: Explain here what this function does, and how to use it.

# SCRIPT PARAMETERS
# -----------------

# Number of factors (topics).
K <- 13

# SET UP ENVIRONMENT
# ------------------
library(readr)

# Initialize the sequence of pseudorandom numbers.
set.seed(1)

# LOAD DATA
# ---------
cat("Loading count data.\n")
counts        <- suppressMessages(read_csv("../data/droplet.csv.gz"))
class(counts) <- "data.frame"
counts        <- as.matrix(counts)
cat(sprintf("Loaded %d x %d counts matrix.\n",nrow(counts),ncol(counts)))
