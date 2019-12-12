# devtools::install_github("vinecopulib/rvinecopulib")
library(rvinecopulib)
library(ggraph)
pollution_data <- read.csv("data/clean_pollution_data.csv")
pollution_data <- pollution_data[,-1]
origin_pollution_data <- read.csv("data/non_gpd_clean_data.csv")
origin_pollution_data <- origin_pollution_data[,-1]

pollution_gmm <- read.csv('analysis/pollution/results/results_100k_gmm.csv')
pollution_gmm <- pollution_gmm[(0:(nrow(pollution_gmm)/4+1)) * 3 + 1,]
row.names(pollution_gmm) <- pollution_gmm[,1]
pollution_gmm <- pollution_gmm[,-1]

unif_pollution_data <- UniformFromGPDForMatrix(
  dataset_origin = origin_pollution_data,
  dataset = pollution_data,
  params = pollution_gmm)

vine_config <- list()
vine_config[['family_set']] =  c("clayton", "gumbel", "indep")
vine_config[['trunc_lvl']] = NA
vine_config[['selcrit']] = 'aic'
vine_config[['core']] = parallel::detectCores()-1
vine_config[['show_trace']] = FALSE
# full_vine_mbic <- rvinecopulib::vinecop(data = unif_pollution_data,
#                                         family_set = c("clayton", "gumbel", "indep"),
#                                         trunc_lvl = NA,
#                                         selcrit='mbicv',
#                                         core=parallel::detectCores()-1)
# print(full_vine_mbic)
# full_vine_mbic$structure$struct_array
# plot(full_vine_mbic)
# 
# full_vine_bic <- rvinecopulib::vinecop(data = unif_pollution_data,
#                                         family_set = c("clayton", "gumbel", "indep"),
#                                         trunc_lvl = NA,
#                                         selcrit='bic',
#                                         core=parallel::detectCores()-1)
# print(full_vine_bic)

evc <- ExtremeVineCollection(pollution_data, unif_pollution_data, horizons = c(1), vine_config = vine_config, rescaling = T)
evc$O3[[1]]$quantile_values


col_number_tmp <- 1
cond_value <- 0.96
vine_tmp <- evc[[col_number_tmp]][[3]]$vine_fit
vine_struc <- evc[[col_number_tmp]][[3]]$vine_fit$structure

samples <- t(vapply(950:1000/1000, function(x){ apply(ExtremeVineConditionalSimulation(vine_tmp, col_number = col_number_tmp, value = x, n = 1), 2, mean)},
                    rep(0, ncol(pollution_data))))

par(mai=c(0.9,0.6,0.4,0.1), mfrow=c(6,1))
plot(950:1000/1000, samples[,col_number_tmp],
     ylab=paste('Probability of', colnames(pollution_data)[col_number_tmp]),
     xlab = paste(colnames(pollution_data)[col_number_tmp], ' quantiles'), main = 'Probability of extremes',
     cex.lab=1.5, cex.axis=1.4)
for(i in c(1:ncol(pollution_data))[-col_number_tmp]){
  plot(950:1000/1000, samples[,i],
       ylab=paste('Probability of', colnames(pollution_data)[i]),
       xlab=paste(colnames(pollution_data)[col_number_tmp], ' quantiles'),
       cex.lab=1.5, cex.axis=1.4)
  abline(h = evc[[col_number_tmp]][[1]]$quantile_values[i], lty=2, lwd=2)
}


ExtremeVineExtractConditional(evc$O3[[1]], 4)
