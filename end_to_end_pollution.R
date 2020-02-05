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

n_ctron <- 100
predict_train_index <- 1:100
predict_test_index <- 1:100


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
    paste('data/end_to_end_pollution/vine_checkpoints/', Sys.Date(), '_vines.RData', sep = ''))

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
                 paste('data/end_to_end_pollution/tron/', Sys.Date(), '_TRON.RData', sep = ''))

#############################
#############################

# C-TRON

# Prediction results
train_prediction_output <- ExtremeVinePredictData(pollution_data, col_cond_on, 1)
test_prediction_output <- ExtremeVinePredictData(test_pollution_data, col_cond_on, 1)



pred_2$pred

conditiona_tron_compute <- list()
for(horizon_number in 1:length(vine_horizons_set)){
  ctron_for_one_horizon <- list()
  for(col_cond_on in 1:ncol(pollution_data)){
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
    
    xvine_train_test_data <- c(
      input_data$xvine_data[predict_train_index, final_col],
      input_data$xvine_test_data[predict_test_index, final_col]
    )
    
    two_phases_predict <- list()
    tag_index <- 1
    for(tag in c('train', 'test')){
      prediction_results <- list()
      pred_1 <- ExtremeVineConditionalPredict(
        vine = vine_tmp,
        quantile_values = vine_quantiles,
        col_number = final_col,
        values = xvine_train_test_data[tag_index],
        n = n_ctron
      )
      prediction_results[['conditional']] <- pred_1
      
      pred_2 <- ExtremeVineConditionalIndicatorPredict(
        vine = vine_tmp,
        quantile_values = vine_quantiles,
        col_number = final_col,
        values = input_data$xvine_data[1:100],
        n = n_ctron
      )
      prediction_results[['indicator']] <- pred_2
      
      two_phases_predict[[tag]] <- prediction_results
    }
    
    ev_ctrons[['predictions']] <- two_phases_predict
    ev_ctrons[['horizon']] <- vine_horizons_set[horizon_number]
    
    ctron_for_one_horizon[[colnames(pollution_data)[col_cond_on]]] <- ev_trons
  }
  conditiona_tron_compute[[horizon_number]] <- ctron_for_one_horizon  
}

