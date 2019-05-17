# TO DO: Explain here what this script does, and how to use it.

# LOAD DATA
# ---------
cat("Loading data set.\n")
mixdata <- readRDS("../data/mixdata.rds")
L       <- mixdata$L
w       <- mixdata$w
