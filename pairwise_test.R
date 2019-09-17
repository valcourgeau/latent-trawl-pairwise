library('zeallot')
library(compiler)
library(parallel)
library(assertthat)

pollution_data <- read.csv('data/clean_pollution_data.csv')

test_params <- noven_example_params
test_params[1] <- 1/noven_example_params[1]
test_params[2] <- (noven_example_params[2] + noven_example_params[3])/abs(noven_example_params[1])

std_param_pl <- TrawlPL(data = pollution_data$NO2[1:1000], depth = 5, parametrisation='standard', parallel = F)
std_param_pl(test_params)

noven_param_pl <- TrawlPL(data = pollution_data$NO[1:100], depth = 3, parametrisation='noven')
noven_param_pl(noven_example_params)


# library('optimParallel')
# cl <- makeCluster(detectCores() - 1)
# optimParallel(std_param_pl, parallel = list(cl=cl), par = test_params, method = 'L-BFGS-B', lower=c(0.05, 1, 1, 1e-2), upper=c(0.35, 6, 20, 1), control = list(trace=3))
# parallel::stopCluster(cl)

optim(std_param_pl, par = test_params, method = 'L-BFGS-B', lower=c(0.1, 0.5, 1, 1e-2), upper=c(2, 6, 20, 1), control = list(trace=3))

