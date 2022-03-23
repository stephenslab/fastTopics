# TO DO: Explain here briefly what this script is for.
library(Matrix)
set.seed(1)
X <- simulate_count_data(80,100,k = 4,sparse = TRUE)$X
fit0 <- fit_topic_model(X,k = 4)
fit <- fit_topic_model_map(X,fit0,numiter = 20)
