library(rvinecopulib)
library(ggraph)

#############################
#############################
# CONFIGURATION

# VINE CFG
vine_config <- list()
vine_config[['family_set']] =  c("clayton", "gumbel", "indep")
vine_config[['trunc_lvl']] = NA
vine_config[['selcrit']] = 'aic'
vine_config[['core']] = parallel::detectCores()-1
vine_config[['show_trace']] = FALSE

vine_horizons_set <- c(1, 2) # c(1, 2, 3, 4, 6, 12, 24, 48)

n_tron <- 100

n_ctron <- 10
predict_train_index <- 1:10
predict_test_index <- 1:10


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

subdataset_collection <- list()
subdataset_collection[['pollution_data']] <- pollution_data
subdataset_collection[['test_pollution_data']] <- test_pollution_data
subdataset_collection[['empirical_quantiles']] <- empirical_quantiles
subdataset_collection[['test_empirical_quantiles']] <- test_empirical_quantiles
subdataset_collection[['bool_pollution_data']] <- bool_pollution_data
subdataset_collection[['test_bool_pollution_data']] <- test_bool_pollution_data
subdataset_collection[['origin_pollution_data']] <- origin_pollution_data
subdataset_collection[['test_origin_pollution_data']] <- test_origin_pollution_data
subdataset_collection[['pollution_gmm']] <- pollution_gmm
subdataset_collection[['unif_pollution_data']] <- unif_pollution_data
subdataset_collection[['test_unif_pollution_data']] <- test_unif_pollution_data

rlist::list.save(subdataset_collection,
     paste('data/end_to_end_pollution/subdatasets/', Sys.Date(), '_subdataset_collection.RData', sep = ''))

#############################
#############################

evc <- ExtremeVineCollection(dataset = pollution_data, uniform_dataset = unif_pollution_data,
                             horizons = vine_horizons_set, vine_config = vine_config, rescaling = T)
rlist::list.save(evc,
    paste('analysis/pollution/end_to_end_pollution/vine_checkpoints/', Sys.Date(), '_vines.RData', sep = ''))

#############################
#############################

# TRON

tron_compute <- list()
for(horizon_number in 1:length(vine_horizons_set)){
  tron_for_one_horizon <- list()
  for(col_cond_on in 1:ncol(pollution_data)){
    final_col <- ncol(pollution_data) + 1
    vine_tmp <- evc[[col_cond_on]][[1]]$vine_fit
    vine_struc <- evc[[col_cond_on]][[1]]$vine_fit$structure
    vine_quantiles <- evc[[col_cond_on]][[1]]$quantiles
    
    input_data_tron <- ExtremeVineData(
      dataset = pollution_data,
      uniform_dataset = unif_pollution_data,
      horizon = 1,
      col_number = col_cond_on,
      rescaling = T
    )
    
    # TRONs
    ev_trons <- ExtremeVineTRON(vine = vine_tmp, extreme_quantile = empirical_quantiles[col_cond_on],
                                quantile_values = vine_quantiles,
                                col_number = final_col, cond_threshold = 0.0,
                                ecdf_rescaling = input_data_tron$ecdf_rescaling[[final_col]],
                                xi = pollution_gmm[col_cond_on, 1], sigma = pollution_gmm[col_cond_on, 2],
                                n=n_tron, seed=42)
    ev_trons[['horizon']] <- vine_horizons_set[horizon_number]
    tron_for_one_horizon[[colnames(pollution_data)[col_cond_on]]] <- ev_trons
  }
  tron_compute[[horizon_number]] <- tron_for_one_horizon  
}

rlist::list.save(tron_compute,
                 paste('analysis/pollution/end_to_end_pollution/tron/', Sys.Date(), '_TRON.RData', sep = ''))

#############################
#############################

# C-TRON

# Prediction results
train_prediction_output <- ExtremeVinePredictData(pollution_data, col_cond_on, 1)
test_prediction_output <- ExtremeVinePredictData(test_pollution_data, col_cond_on, 1)

conditional_tron_compute <- list()
for(horizon_number in 1:length(vine_horizons_set)){
  ctron_for_one_horizon <- list()
  for(col_cond_on in 1:ncol(pollution_data)){
    print(col_cond_on)
    ev_ctrons <- list()
    final_col <- ncol(pollution_data) + 1
    vine_tmp <- evc[[col_cond_on]][[1]]$vine_fit
    vine_struc <- evc[[col_cond_on]][[1]]$vine_fit$structure
    vine_quantiles <- evc[[col_cond_on]][[1]]$quantiles
    
    # true values
    train_prediction_output <- ExtremeVinePredictData(test_pollution_data, col_cond_on, 1)
    test_prediction_output <- ExtremeVinePredictData(test_pollution_data, col_cond_on, 1)
    
    # input data
    input_data <- ExtremeVineTestData(
      dataset = pollution_data,
      uniform_dataset = unif_pollution_data,
      test_dataset = test_pollution_data,
      test_uniform_dataset = test_unif_pollution_data,
      horizon = 1,
      col_number = col_cond_on,
      rescaling = T
    )
    
    xvine_train_test_data <- list(
      'train'=input_data$xvine_data[predict_train_index, final_col],
      'test'=input_data$xvine_test_data[predict_test_index, final_col]
    )
    
    two_phases_predict <- list()
    for(tag in c('train', 'test')){
      print(tag)
      prediction_results <- list()
      pred_1 <- ExtremeVineConditionalPredict(
        vine = vine_tmp,
        quantile_values = vine_quantiles,
        col_number = final_col,
        values = xvine_train_test_data[[tag]],
        n = n_ctron
      )
      prediction_results[['avg_pred']] <- pred_1
      
      pred_2 <- ExtremeVineConditionalIndicatorPredict(
        vine = vine_tmp,
        quantile_values = vine_quantiles,
        col_number = final_col,
        values = input_data$xvine_data[1:100],
        n = n_ctron
      )
      prediction_results[['indicator_pred']] <- pred_2
      
      two_phases_predict[[tag]] <- prediction_results
    }
    
    ev_ctrons[['predictions']] <- two_phases_predict
    ev_ctrons[['horizon']] <- vine_horizons_set[horizon_number]
    
    ctron_for_one_horizon[[colnames(pollution_data)[col_cond_on]]] <- ev_ctrons
  }
  conditional_tron_compute[[horizon_number]] <- ctron_for_one_horizon  
}

rlist::list.save(conditional_tron_compute,
                 paste('analysis/pollution/end_to_end_pollution/pollution/end_to_end_pollution/conditional_tron/', Sys.Date(), '_TRON.RData', sep = ''))


#############################
#############################

# Quantiles vs conditional values

# Shortcuts

for(horizon_number in 1:length(vine_horizons_set)){
  for(col_cond_on in 1:ncol(pollution_data)){
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
    
    multiplier <- 1.4
    pic_name <- paste("analysis/pollution/end_to_end_pollution/images/response_quantiles/",  vine_horizons_set[horizon_number], "/response_",
                      colnames(pollution_data)[col_cond_on], "_horizon", vine_horizons_set[horizon_number], ".png", sep = "")
    png(pic_name, width = 800, height = 800)
    layout(matrix(c(1,1,2,3,4,5,6,7), 4, 2, byrow = TRUE),
           widths=c(1,1), heights=c(1, multiplier, multiplier, multiplier))
    plot(1:100/100, samples[,final_col],
         ylab=paste('(Control)'),
         xlab = paste(colnames(pollution_data)[col_cond_on], ' quantiles'),
         type='l', lwd=3,
         main = paste('Response Vine Quantiles as function', colnames(pollution_data)[col_cond_on], 'quantiles'),
         cex.lab=1.5, cex.axis=1.7, cex.main = 1.5)
    for(i in c(1:ncol(pollution_data))){
      plot(1:100/100, samples[,i],
           ylab=paste('Quantiles ', colnames(pollution_data)[i], ' at t+', vine_horizons_set[horizon_number], sep = ''),
           xlab=paste(colnames(pollution_data)[col_cond_on], ' quantiles at time t'), 
           type='l', lwd=3,
           main=paste('(Horizon ', vine_horizons_set[horizon_number], ') ',
                      colnames(pollution_data)[i], ' with corr ',  round(cor_mat[col_cond_on, i], 2), sep=''),
           cex.lab=1.5, cex.axis=1.7, cex.main=1.7)
      abline(h = evc[[col_cond_on]][[horizon_number]]$quantiles[i], lty=2, lwd=2)
      abline(v = evc[[col_cond_on]][[horizon_number]]$quantiles[col_cond_on], lty=3, lwd=2)
    }
    dev.off()
  }
}
