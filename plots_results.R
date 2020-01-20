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
    
    print(apply(data_gmm, 2, mean))
    print(apply(data_gmm, 2, sd)/sqrt(length(data_gmm[,1])))
    vioplot(data_gmm[,i])
    xi_tmp <- data_main[2,1]
    sigma_tmp <- data_main[2,2]
    kappa_tmp <- data_main[2,3]
    rho_tmp <- data_main[2,4]
    print(xi_tmp)
    paste(xi_tmp, ' ', sigma_tmp, ' ', kappa_tmp, ' ', rho_tmp)
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
par(mfrow=c(3,2), mai=c(0.6,0.7,0.6,0.1))

ylims <- c(12,14,14,10,17,20)
for(i in 1:ncol(pollution_data)){
  data_main <- read.csv(paste('analysis/pollution/results/', colnames(pollution_data)[i],'_full_pl.csv', sep=''))
  data_main <- data_main[,-1]
  data_pl <- read.csv(paste('analysis/pollution/results/', colnames(pollution_data)[i],'_full_200_PL.csv', sep=''))
  data_pl <- data_pl[,-1]
  data_pl <- data_pl[-1,]
  data_main <- read.csv(paste('analysis/pollution/results/', colnames(pollution_data)[i],'_full_gmm.csv', sep=''))
  data_main <- data_main[,-1]
  data_gmm <- read.csv(paste('analysis/pollution/results/', colnames(pollution_data)[i],'_sub_sample_gmm_200.csv', sep=''))
  data_gmm <- data_gmm[,-1]
  data_gmm <- data_gmm[-1,]
  
  col_name <- colnames(pollution_data)[i]
  
  dark_blue <- rgb(0.019, 0.06019608, 0.87058824,0.65)
  light_blue <- rgb(0.01960784, 0.61960784, 0.90196078, 0.5)
  hist(data_pl[,4][data_pl[,4]>1e-2], breaks=20, col=dark_blue,
       probability = T, ylim=c(0,ylims[i]), xlim=c(0,1),
       main=bquote("Bootstrap" ~ .(col_name) ~ rho),
       xlab = 'Value',
       cex.axis=1.9, cex.lab=2.2, cex.main=2.5)
  hist(pmax(data_gmm[,4]*(1+rnorm(length(data_gmm[,4]), sd = 0.1)), 0.0),
       add=T, breaks=10, col=light_blue, probability = T)
  legend(0.6,ylims[i]-2, legend=c('PL', 'GMM'), cex=2.,
         fill=c(dark_blue, light_blue))
}


# extreme propagation

EmpiricalPropagation(pollution_data[,4], 0, 3)

par(mfrow=c(3,2), mai=c(0.6,0.7,0.6,0.1))
ylims <- c(12,14,14,10,17,20)
for(i in 1:ncol(pollution_data)){
  results_pl
  col_name <- colnames(pollution_data)[i]
  
  dark_blue <- rgb(0.019, 0.06019608, 0.87058824,0.65)
  light_blue <- rgb(0.01960784, 0.61960784, 0.90196078, 0.5)
  plot(0:10,
    PropagationProbability(x_now = 0,
                           x_future = 0,
                           h = 0:10,
                           xi = results[3*(i-1)+1, 1],
                           sigma = results[3*(i-1)+1, 2],
                           kappa = results[3*(i-1)+1, 3],
                           rho = results[3*(i-1)+1, 4]),
    ylim=c(0,1.2), ylab='Conditional Probability',
    cex.axis=1.9, cex.lab=2.2, cex.main=2.5
  )
  legend(0.6,ylims[i]-2, legend=c('PL', 'GMM'), cex=2.,
         fill=c(dark_blue, light_blue))
}



