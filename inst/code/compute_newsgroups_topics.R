# Analyze the "20 Newsgroups" data using fastTopics.
#
# sinteractive --mem=24G -c 8 --time=24:00:00
# module load R/4.2.0
# .libPaths()[1]
# /home/pcarbo/R_libs_4_20
library(tools)
library(Matrix)
library(fastTopics)
load("../datafiles/newsgroups.RData")
set.seed(1)

# Remove words that appear in fewer than 10 documents.
x <- colSums(counts > 0)
j <- which(x > 9)
counts <- counts[,j]

# Fit a Poisson NMF using fastTopics, with k = 10 factors/topics.
pnmf <- fit_poisson_nmf(counts,k = 10,numiter = 100,method = "em",
                        control = list(numiter = 4,nc = 8,extrapolate = FALSE),
                        init.method = "random",verbose = "detailed")
pnmf <- fit_poisson_nmf(counts,fit0 = pnmf,numiter = 100,method = "scd",
                        control = list(numiter = 4,nc = 8,extrapolate = TRUE),
                        verbose = "detailed")

# Perform the "grade of membership" differential expression analysis
# using the fitted Poisson NMF model.
de_le <- de_analysis(pnmf,counts,shrink.method = "ash",
                     lfc.stat = "le",pseudocount = 0.1,
                     control = list(ns = 1e4,nc = 8,nsplit = 1000))
de_vsnull <- de_analysis(pnmf,counts,shrink.method = "ash",
                         lfc.stat = "vsnull",pseudocount = 0.1,
                         control = list(ns = 1e4,nc = 8,nsplit = 1000))

# Save the outputs to an .Rdata file.
session_info <- sessionInfo()
save(list = c("pnmf","de_le","de_vsnull","session_info"),
     file = "newsgroups_topics.RData")
resaveRdaFiles("newsgroups_topics.RData")
