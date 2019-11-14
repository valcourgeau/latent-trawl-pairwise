source('dupuis_regression.R')
source('sim_study_utils.R')

pollution_data <- read.csv("data/clean_pollution_data.csv")
pollution_data <- pollution_data[,-1]
head(pollution_data)

max_length <- 100000
marginal_fit <- lapply(1:ncol(pollution_data), 
       function(i){
         SubSampleFit(
           data=pollution_data[1:max_length,i],
           depth = 5,
           sub_length = max_length-1,
           method = 'PL',
           trials = 1,
           parallel = F,
           file_csv = paste(colnames(pollution_data)[i], '_full_pl.csv', sep=''),
           subfolder = 'analysis/pollution/results/',
           n_trials=50,
           seed=41,
           acf_depth=20,
           bounds='ow'
         )
       }
)

results_numeric <- marginal_fit[[1]][1,]
for(res in marginal_fit){
  results_numeric <- rbind(results_numeric, res[2,])
}
colnames(results_numeric) <- c('xi', 'sigma', 'kappa', 'rho')
results <- list(params = t(as.matrix(marginal_fit[[1]][1,])), numerics = results_numeric[-1,])
colnames(results$params) <- c('N', 'delta', 'sub_length', 'trials')
rownames(results$numerics) <- colnames(pollution_data)
results
rlist::list.save(results, 'analysis/pollution/results/results.RData')

max_length <- 50000
marginal_fit <- lapply(1:ncol(pollution_data), 
                       function(i){
                         SubSampleFit(
                           data =pollution_data[1:max_length,i],
                           depth = 5,
                           sub_length = 5000,
                           method = 'GMM',
                           trials = 100,
                           file_csv = paste(colnames(pollution_data)[i], '_sub_sample_1.csv', sep=''),
                           subfolder = 'analysis/pollution/results/',
                           parallel = T,
                           n_trials=30,
                           seed=41
                         )
                       }
)



SubSampleFit(
  data = pollution_data[,4],
  depth = 5,
  sub_length = 10000,
  method = 'GMM',
  trials = 3,
  file_csv = paste('o3_test', '.csv', sep=''),
  parallel = T,
  n_trials=20
)
