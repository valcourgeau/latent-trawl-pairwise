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

vine_horizons_set <- c(1,2,3,4,6,12,24,48)

evc <- ExtremeVineCollection(dataset = pollution_data, uniform_dataset = unif_pollution_data,
                             horizons = vine_horizons_set, vine_config = vine_config, rescaling = T)

