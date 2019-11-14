
# air pollution
pollution_data <- read.csv("data/clean_pollution_data.csv")
pollution_data <- pollution_data[,-1]
head(pollution_data)


max_length <- 100000

for(i in 1:length(colnames(pollution_data))){
  plot_data <- pollution_data[1:max_length,i]
  pl_results <- SubSampleFit(
    data=pollution_data[1:max_length,i],
    depth = 5,
    sub_length = max_length-1,
    method = 'PL',
    trials = 1,
    parallel = F,
    file_csv = paste(colnames(pollution_data)[i], '_full_1_PL.csv', sep=''),
    subfolder = 'analysis/pollution/results/',
    bounds='ow',
    n_trials=20,
    seed=41,
    acf_depth=20
  )
}

for(i in 1:length(colnames(pollution_data))){
  plot_data <- pollution_data[1:max_length,i]
  pl_results <- SubSampleFit(
    data=pollution_data[1:max_length,i],
    depth = 5,
    sub_length = 5000,
    method = 'PL',
    trials = 200,
    parallel = T,
    file_csv = paste(colnames(pollution_data)[i], '_full_200_PL.csv', sep=''),
    subfolder = 'analysis/pollution/results/',
    bounds='ow',
    n_trials=30,
    seed=41,
    acf_depth=20
  )
}

