pollution_data <- read.csv('data/clean_pollution_data.csv')
# std_param_pl <- TrawlPL(data = pollution_data$NO2[1:100000],
#                         depth = 6,
#                         parametrisation='standard',
#                         parallel = T,
#                         type = 'exp')

small <- (length(pollution_data$O3)-25000+1):length(pollution_data$O3)

pd_small <- apply(pollution_data[,-1],
        function(x){
          EVTrawlFit(x[small],
                 depth = 5,
                 parametrisation = 'standard',
                 type = 'exp', parallel = T)
          }, MARGIN = 2)
pd_small <- t(pd_small)
pd_small

pd_std_trf <- apply(pollution_data[,-1],
                  function(x){
                    EVTrawlFit(x[small],
                               depth = 5,
                               parametrisation = 'std_trf',
                               type = 'exp', parallel = T, bounds = 'bd', n_trials=20, acf_depth=20)
                  }, MARGIN = 2)
pd_std_trf <- t(pd_std_trf)
pd_std_trf

# write.csv(pd_std_trf, file = 'analysis/pollution/pollution_std_trf_all_depth_5_25k.csv')

pd_sum_exp <- apply(pollution_data[,-1],
            function(x){
              EVTrawlFit(x[small],
                         depth = 5,
                         parametrisation = 'standard',
                         type = 'sum_exp', parallel = T)
            }, MARGIN = 2)
pd_sum_exp <- t(pd_sum_exp)
pd_sum_exp

pd_gamma <- apply(pollution_data[,-1],
                    function(x){
                      EVTrawlFit(x[small],
                                 depth = 5,
                                 parametrisation = 'standard',
                                 type = 'gamma', parallel = T)
                    }, MARGIN = 2)
pd_gamma <- t(pd_gamma)
pd_gamma

# write.csv(pd, file = 'analysis/pollution/pollution_data_all_depth_4_25k.csv')


par <- ParametrisationTranslator(pd_std_trf[3,], 'standard', target = 'noven')
par
alpha <- par[1]
beta <- par[2]
kappa <- par[3]
trawl_param <- par[4:length(par)]
acf_indices <- 0:50
acf(pollution_data$NO[small])
lines(acf_indices, acf_trawl_num_approx(acf_indices, alpha = alpha, beta = beta, kappa = kappa, rho = trawl_param, delta = 0.1, type = 'exp'))


DupuisSimplified(data_u = pollution_data$O3[small], n_trials = 20, acf_depth = 20)


par <- ParametrisationTranslator(pd[5,], 'standard', target = 'noven')
alpha <- par[1]
beta <- par[2]
kappa <- par[3]
trawl_param <- par[4:length(par)]
acf_indices <- 0:50
acf(pollution_data$PM10)
lines(acf_indices, acf_trawl_num_approx(acf_indices, alpha = alpha, beta = beta, kappa = kappa, rho = trawl_param))



pd_std_trf <- apply(pollution_data[,-1],
                function(x){
                  EVTrawlFit(x[1:10000],
                             depth = 5,
                             parametrisation = 'std_trf',
                             type = 'exp', parallel = T)
                }, MARGIN = 2)
pd_std_trf <- t(pd_std_trf)
pd_std_trf
# write.csv(pd_std_trf, file = 'analysis/pollution/pollution_data_all_trf_100k.csv')


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

