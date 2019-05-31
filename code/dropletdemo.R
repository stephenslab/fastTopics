# TO DO: Explain here what this script does, and how to use it.

# SET UP ENVIRONMENT
# ------------------
library(parallel)
library(Rcpp)
library(readr)
library(ggplot2)
library(cowplot)
source("misc.R")
source("betanmf.R")
source("mixsqp.R")
source("altsqp.R")
sourceCpp("mixsqp.cpp")

# LOAD DATA
# ---------
cat("Loading count data.\n")
counts <- read_csv("../data/droplet.csv.gz",col_names = FALSE,progress = FALSE,
                   col_types = cols(.default = col_double()))
class(counts) <- "data.frame"
counts        <- as.matrix(counts)
n             <- nrow(counts)
p             <- ncol(counts)
cat(sprintf("Loaded %d x %d counts matrix.\n",n,p))

# LOAD INITIAL ESTIMATES
# ----------------------
cat("Loading initial estimates of factors and loadings.\n")
F0 <- read_csv("../data/droplet_factors_rough.csv.gz",col_names = FALSE,
               progress = FALSE,col_types = cols(.default = col_double()))
L0 <- read_csv("../data/droplet_loadings_rough.csv.gz",col_names = FALSE,
               progress = FALSE,col_types = cols(.default = col_double()))
class(F0) <- "data.frame"
class(L0) <- "data.frame"
F0        <- as.matrix(F0)
L0        <- as.matrix(L0)

stop()

# RUN ALTERNATING SQP METHOD
# --------------------------
cat("Fitting Poisson topic model using alternating SQP method.\n")
fit <- altsqp(counts,F0,L0,nc = 2,numiter = 10)
