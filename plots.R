
# air pollution
par(mai=c(0.6,0.55,0.1,0.2), mfrow=c(6,4))
pollution_data <- read.csv("data/clean_pollution_data.csv")
pollution_data <- pollution_data[,-1]
head(pollution_data)

air_pollution_results <- rlist::list.load('analysis/pollution/results/results.RData')
air_pollution_results$numerics


max_length <- 100000
i <- 1
for(i in 1:length(colnames(pollution_data))){
  plot_data <- pollution_data[1:max_length,i]
  plot_nums <- air_pollution_results$numerics
  # pl_results <- SubSampleFit(
  #   data=pollution_data[1:max_length,i],
  #   depth = 5,
  #   sub_length = max_length-1,
  #   method = 'PL',
  #   trials = 1,
  #   parallel = F,
  #   file_csv = paste(colnames(pollution_data)[i], '_full_1_PL.csv', sep=''),
  #   subfolder = 'analysis/pollution/results/',
  #   n_trials=20,
  #   seed=41,
  #   acf_depth=20
  # )
  
  pl_results <- read.csv(paste('analysis/pollution/results/',
                               colnames(pollution_data)[i],
                         '_full_pl.csv', sep=''))
  DupuisSimplified(
    data_u = plot_data,
    n_trials = 30,
    acf_depth = 10,
    plot.it=T,
    mult_fac = c(pl_results[2,5], 0)
  )
  
  # histograms
  alpha_tmp <- if(plot_nums[i,1] < 0){2}else{1/plot_nums[i,1]}
  beta_tmp <- plot_nums[i,2] / abs(plot_nums[i,1]) - plot_nums[i,3]
  acf_depth <- 30
  acf_vals <- vapply(c(0.01, 1:(acf_depth-1)), function(h){
    acf_trawl(h, alpha = alpha_tmp, beta = beta_tmp, kappa = plot_nums[i,3], 
              rho = plot_nums[i,4], delta = 0.5, end_seq = 50)}, 1)
  acf_vals_pl <- vapply(c(0.01, 1:(acf_depth-1)), function(h){
    acf_trawl(h, alpha = alpha_tmp, beta = beta_tmp, kappa = plot_nums[i,3], 
              rho = pl_results[2,5], delta = 0.5, end_seq = 50)}, 1)
  
  acf(plot_data, lag.max = 30, bty='n',
      cex.axis=1.5,
      cex.lab=1.5,
      main='')
  lines(0:(acf_depth-1), acf_vals,
        lty=2, lwd=3, col='dodgerblue3',
        cex=1.3, cex.axis=1.5, cex.lab=1.3)
  lines(0:(acf_depth-1), acf_vals_pl,
        lty=4, lwd=3, col='darkblue',
        cex=1.3, cex.axis=1.5, cex.lab=1.3)
  legend(15, 0.9,
         legend=c(expression(paste(rho, ' PL')), expression(paste(rho, ' GMM'))),
         col=c('darkblue', "dodgerblue3"),
         pch = c(NA, NA), lty=c(2, 4), cex=1.1, lwd=2,
         bty='n')
  
  hist(plot_data[plot_data>0], breaks = 50, probability = T,
       xlab = paste(colnames(pollution_data)[i],' Value'),
       ylab='Density',
       cex.axis=1.5,
       cex.lab=1.5,
       main='')
  plotting_x <- seq(0,max(plot_data[plot_data>0]),length.out=100)
  lines(plotting_x, 
        evir::dgpd(plotting_x, xi = plot_nums[i,1], beta=plot_nums[i,2]),
        lwd=3, lty=2, col='dodgerblue3', cex=1.3)
  fExtremes::qqparetoPlot(plot_data[plot_data>0], xi = plot_nums[i,1],
                          labels=F, bty='n', cex.axis=1.3)
  abline(v=quantile(plot_data[plot_data>0], 0.95), col='black', lwd=2, lty=2)
  title(main="", sub="",
        xlab=paste(colnames(pollution_data)[i], '(ordered)'), ylab="GPD Quantiles",
        cex.axis=1.5,
        cex.lab=1.5)
}

