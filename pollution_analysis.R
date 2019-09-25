pollution_data <- read.csv('data/clean_pollution_data.csv')
# std_param_pl <- TrawlPL(data = pollution_data$NO2[1:100000],
#                         depth = 6,
#                         parametrisation='standard',
#                         parallel = T,
#                         type = 'exp')


pd_2456 <- apply(pollution_data[,-1][,c(2,4,5,6)],
        function(x){
          EVTrawlFit(x[1:100000],
                 depth = 5,
                 parametrisation = 'standard',
                 type = 'exp', parallel = T)
          }, MARGIN = 2)
pd_2456 <- t(pd_2456)
pd_2456
# write.csv(pd_2456, file = 'analysis/pollution/pollution_data_2456_100k.csv')
ev_fit
acf(pollution_data$O3)

# NO2 [1] 0.1274843 0.5843454 1.4550261 0.2659386

# check
acf_indices <- 0:50

par <- ParametrisationTranslator(apply(res[2:1001,], MARGIN = 2, mean), 'standard')
alpha <- par[1]
beta <- par[2]
kappa <- par[3]
trawl_param <- par[4:length(par)]
acf(pollution_data$NO2[1:100000])
lines(acf_indices, acf_trawl_num_approx(acf_indices, alpha = alpha, beta = beta, kappa = kappa, rho = trawl_param))



res <- SubSampleFit(pollution_data$NO2[1:100000], depth = 5, sub_length = 2500, trials = 1000, parallel = T, file_csv = 'no2_5_2500_100.csv', type = 'exp')
for(i in 3:nrow(res)){
  print(apply(res[2:i,], MARGIN = 2, sd)/sqrt(i))
  print(apply(res[2:i,], MARGIN = 2, mean))
}

