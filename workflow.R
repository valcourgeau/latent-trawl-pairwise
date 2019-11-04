source('dupuis_regression.R')
source('sim_study_utils.R')

pollution_data <- read.csv("data/clean_pollution_data.csv")
pollution_data <- pollution_data[,-1]
head(pollution_data)

max_length <- 50000
marginal_fit <- lapply(1:ncol(pollution_data), 
       function(i){
         SubSampleFit(
           data =pollution_data[1:50000,i],
           depth = 5,
           sub_length = 5000,
           method = 'GMM',
           trials = 100,
           file_csv = paste(colnames(pollution_data)[i], '_sub_sample_100.csv', sep=''),
           subfolder = 'analysis/pollution/results/',
           parallel = T,
           n_trials=20,
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
