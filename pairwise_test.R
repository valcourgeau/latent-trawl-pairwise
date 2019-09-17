library('zeallot')
library(compiler)
library(parallel)
library(assertthat)

pollution_data <- read.csv('data/clean_pollution_data.csv')

test_params <- noven_example_params
test_params[1] <- 1/noven_example_params[1]
test_params[2] <- (noven_example_params[2] + noven_example_params[3])/abs(noven_example_params[1])

std_param_pl <- TrawlPL(data = pollution_data$NO2[1:100000], depth = 6, parametrisation='standard', parallel = T)
std_param_pl(test_params)


data_tmp <- pollution_data$PM10[1:50000]
init_params <- DupuisSimplified(data_tmp)
init_params[3:4] <- init_params[4:3]
init_params <- init_params[1:4]
names(init_params) <- c('xi', 'sigma', 'kappa', 'rho')

std_param_pl_blocked <- function(trawl_params){
return(std_param_pl(c(init_params[1:3], trawl_params)))  
}

noven_param_pl <- TrawlPL(data = pollution_data$NO[1:100], depth = 3, parametrisation='noven')
noven_param_pl(noven_example_params)


# library('optimParallel')
# cl <- makeCluster(detectCores() - 1)
# optimParallel(std_param_pl, parallel = list(cl=cl), par = test_params, method = 'L-BFGS-B', lower=c(0.05, 1, 1, 1e-2), upper=c(0.35, 6, 20, 1), control = list(trace=3))
# parallel::stopCluster(cl)

optim_pm10_full <- optim(std_param_pl, par = init_params, method = 'L-BFGS-B', lower=c(0.1, 0.1, 1, 1e-2), upper=c(1.0, 6, 5, 1), control = list(trace=3))
change_rate <- 0.3
optim_pm10_restricted <- optim(std_param_pl, par = init_params, method = 'L-BFGS-B',
                               lower=init_params*(1-change_rate), upper=init_params*(1+change_rate),
                               control = list(trace=3))



optim_pm10 <- optim(std_param_pl_blocked, par = init_params[4], method = 'L-BFGS-B', lower=c(0.01), upper=c(2), control = list(trace=3))
optim_pm10$par

h_compute <- 0:40
acf_tmp <- acf_trawl_num_approx(h_compute, alpha=1/init_params[1], beta = init_params[2]/init_params[1]-init_params[3], kappa = init_params[3], rho = optim_pm10_full$par)
acf(data_tmp)
lines(h_compute, acf_tmp)
acf_tmp <- acf_trawl_num_approx(h_compute, alpha=1/optim_pm10_restricted$par[1], beta = optim_pm10_restricted$par[2]/optim_pm10_restricted$par[1]-optim_pm10_restricted$par[3], kappa = optim_pm10_restricted$par[3], rho = optim_pm10_restricted$par[4], delta = 0.01)
lines(h_compute, acf_tmp)

eva::dgpd(data)

