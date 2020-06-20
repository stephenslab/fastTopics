# TO DO: Explain here what this script is for, and how to use it.

# Simulate data.
set.seed(1)
dat <- simulate_poisson_gene_data(800,1000,4)
X   <- dat$X
