# devtools::install_github("vinecopulib/rvinecopulib")
library(rvinecopulib)
library(ggraph)

#############################
#############################

# DATA
max_depth <- 100000
pollution_data <- read.csv("data/clean_pollution_data.csv")
pollution_data <- pollution_data[,-1]

# extreme data
test_pollution_data <- pollution_data[(max_depth+1):nrow(pollution_data),]
pollution_data <- pollution_data[1:max_depth,]

# empirical quantiles
empirical_quantiles <- 1.0-apply(pollution_data, 2, function(x){mean(x>0)})
test_empirical_quantiles <- 1.0-apply(test_pollution_data, 2, function(x){mean(x>0)})

# Booleans
bool_pollution_data <- t(apply(pollution_data, 1, function(x){x>0}))
test_bool_pollution_data <- t(apply(test_pollution_data, 1, function(x){x>0}))


origin_pollution_data <- read.csv("data/non_gpd_clean_data.csv")
origin_pollution_data <- origin_pollution_data[,-1]
test_origin_pollution_data <- origin_pollution_data[(max_depth+1):nrow(origin_pollution_data),]
origin_pollution_data <- origin_pollution_data[1:max_depth,]

pollution_gmm <- read.csv('analysis/pollution/results/results_100k_gmm.csv')
pollution_gmm <- pollution_gmm[(0:(nrow(pollution_gmm)/4+1)) * 3 + 1,]
row.names(pollution_gmm) <- pollution_gmm[,1]
pollution_gmm <- pollution_gmm[,-1]

unif_pollution_data <- UniformFromGPDForMatrix(
  dataset_origin = origin_pollution_data,
  dataset = pollution_data,
  params = pollution_gmm)

test_unif_pollution_data <- UniformFromGPDForMatrixTestData(
  dataset_origin = origin_pollution_data,
  dataset = pollution_data,
  test_dataset_origin = test_origin_pollution_data,
  test_dataset = test_pollution_data,
  params = pollution_gmm)


#############################
#############################

# VINE CFG
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


#############################
#############################

# Fitting Vine
evc <- ExtremeVineCollection(dataset = pollution_data, uniform_dataset = unif_pollution_data,
                             horizons = c(1, 72), vine_config = vine_config, rescaling = T)
evc$O3[[1]]$quantiles

#############################
#############################

# Shortcuts
col_cond_on <- 3
final_col <- ncol(pollution_data) + 1
vine_tmp <- evc[[col_cond_on]][[1]]$vine_fit
vine_struc <- evc[[col_cond_on]][[1]]$vine_fit$structure
vine_quantiles <- evc[[col_cond_on]][[1]]$quantiles

# col_number is the column on which we condition on.
# Usually, it is implicitly set to be the last column
samples <- t(vapply(1:100/100,
                    function(x){
                      apply(ExtremeVineConditionalSimulation(
                        vine_tmp, col_number = final_col, value = x, n = 1),
                        MARGIN = 2,
                        FUN = mean)},
                    rep(0, ncol(pollution_data)+1)))

cor_mat <- cor(unif_pollution_data)

par(mai=c(0.55,0.6,0.4,0.1), mfrow=c(4,2))
plot(1:100/100, samples[,final_col],
     ylab=paste('(Control)'),
     xlab = paste(colnames(pollution_data)[col_cond_on], ' quantiles'), main = 'Quantiles as function of time t extreme',
     cex.lab=1.5, cex.axis=1.7, cex.main = 1.8)
for(i in c(1:ncol(pollution_data))){
  plot(1:100/100, samples[,i],
       ylab=paste('Quantiles', colnames(pollution_data)[i], 'at t+1'),
       xlab=paste(colnames(pollution_data)[col_cond_on], ' quantiles at time t'),
       main=paste(colnames(pollution_data)[i], 'with corr',  round(cor_mat[col_cond_on, i], 2)),
       cex.lab=1.5, cex.axis=1.7)
  abline(h = evc[[col_cond_on]][[1]]$quantiles[i], lty=2, lwd=2)
  abline(v = evc[[col_cond_on]][[1]]$quantiles[col_cond_on], lty=3, lwd=2)
}

cor(unif_pollution_data)

ExtremeVineExtractConditional(evc$O3[[1]], 4)

#############################
#############################

# Predicting tests
ExtremeVineConditionalPredict(
  vine = vine_tmp,
  quantile_values = vine_quantiles,
  col_number = final_col,
  values = seq(from=vine_quantiles[col_cond_on], to=0.99, length.out = 5),
  n = 1000
)

ExtremeVineConditionalIndicatorPredict(
  vine = vine_tmp,
  quantile_values = vine_quantiles,
  col_number = final_col,
  values = seq(from=vine_quantiles[col_cond_on], to=0.99, length.out = 5),
  n = 1000
)


#############################
#############################

# Preducting (actual)

ExtremeVineConditionalPredict(
  vine = vine_tmp,
  quantile_values = vine_quantiles,
  col_number = final_col,
  values = test_unif_pollution_data[1:1000,col_cond_on],
  n = 10
)





