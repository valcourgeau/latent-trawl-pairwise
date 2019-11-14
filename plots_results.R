library(vioplot)
pollution_data <- read.csv("data/clean_pollution_data.csv")
pollution_data <- pollution_data[,-1]
head(pollution_data)

for(i in 1:4){
  for(j in 1:ncol(pollution_data)){
    data_main <- read.csv(paste('analysis/pollution/results/', colnames(pollution_data)[j],'_full_gmm.csv', sep=''))
    data_main <- data_main[,-1]
    data_gmm <- read.csv(paste('analysis/pollution/results/', colnames(pollution_data)[j],'_sub_sample_gmm_200.csv', sep=''))
    data_gmm <- data_gmm[,-1]
    data_gmm <- data_gmm[-1,]
    paste(xi_tmp, ' ', sigma_tmp, ' ', kappa_tmp, ' ', rho_tmp)
    print(apply(data_gmm, 2, mean))
    print(apply(data_gmm, 2, sd)/sqrt(length(data_gmm[,1])))
    vioplot(data_gmm[,i])
    xi_tmp <- data_main[2,1]
    sigma_tmp <- data_main[2,2]
    kappa_tmp <- data_main[2,3]
    rho_tmp <- data_main[2,4]
    print(xi_tmp)
  }
}

  
results <- rep(0, 4)
for(i in 1:ncol(pollution_data)){
    data_main <- read.csv(paste('analysis/pollution/results/', colnames(pollution_data)[i],'_full_gmm.csv', sep=''))
    data_main <- data_main[,-1]
    data_gmm <- read.csv(paste('analysis/pollution/results/', colnames(pollution_data)[i],'_sub_sample_gmm_200.csv', sep=''))
    data_gmm <- data_gmm[,-1]
    data_gmm <- data_gmm[-1,]
    results <- rbind(results, (data_main[2,] %>% as.vector))
    results <- rbind(results, (apply(data_gmm, 2, mean)) %>% as.vector)
    results <- rbind(results, (apply(data_gmm, 2, sd)/sqrt(length(data_gmm[,1]))) %>% as.vector)
}
results <- results[-1,]
rownames(results) <- NULL
rownames(results) <- c(
  rbind(
    colnames(pollution_data),
    vapply(colnames(pollution_data), function(x){paste(x,'_mean',sep='')}, '1'),
    vapply(colnames(pollution_data), function(x){paste(x,'_sd',sep='')}, '1')
  )
)
results  
write.csv(results, 'analysis/pollution/results/results_100k_gmm.csv')
hist(data_gmm[,4], breaks = 20)

read.csv('analysis/pollution/results/PM10_.csv')

results_pl <- rep(0, 4)
for(i in 1:ncol(pollution_data)){
  data_main <- read.csv(paste('analysis/pollution/results/', colnames(pollution_data)[i],'_full_pl.csv', sep=''))
  data_main <- data_main[,-1]
  data_gmm <- read.csv(paste('analysis/pollution/results/', colnames(pollution_data)[i],'_full_200_PL.csv', sep=''))
  data_gmm <- data_gmm[,-1]
  data_gmm <- data_gmm[-1,]
  results_pl <- rbind(results_pl, (data_main[2,] %>% as.vector))
  results_pl <- rbind(results_pl, (apply(data_gmm, 2, mean)) %>% as.vector)
  results_pl <- rbind(results_pl, (apply(data_gmm, 2, sd)/sqrt(length(data_gmm[,1]))) %>% as.vector)
}
results_pl <- results_pl[-1,]
rownames(results_pl) <- NULL
rownames(results_pl) <- c(
  rbind(
    colnames(pollution_data),
    vapply(colnames(pollution_data), function(x){paste(x,'_mean',sep='')}, '1'),
    vapply(colnames(pollution_data), function(x){paste(x,'_sd',sep='')}, '1')
  )
)
results - results_pl

# Add histograms / violin plots comparing the two!
