
UniformsFromGPD <- function(data_origin, data, xi, sigma){
  # take oirinal data and zero/GPD-distributed data
  # and transform into uniform(0,1)
  positive_index <- which(data > 0.0)  
  p_zero <- 1 - length(positive_index) / length(data)
  data[positive_index] <- p_zero +  (1-p_zero) * 
    evir::pgpd(q = data[positive_index],
               xi = xi,
               beta = sigma)
  ecdf_before_threshold <- ecdf(data_origin[-positive_index])
  data[-positive_index] <- p_zero * ecdf_before_threshold(data_origin[-positive_index])
  
  return(data)
}

UniformDataFromGPD <- function(dataset_origin, dataset, params){
  # datasets are organised in cols
  # params in rows (xi, sigma)
  return(vapply(
    1:ncol(dataset_origin),
    FUN = function(i){
      UniformsFromGPD(data_origin = dataset_origin[,i],
                      data = dataset[,i],
                      xi = params[i,1],
                      sigma = params[i,2])
    },
    rep(0, nrow(dataset_origin))
  ))
}

pollution_data <- read.csv("data/clean_pollution_data.csv")
pollution_data <- pollution_data[,-1]
origin_pollution_data <- read.csv("data/non_gpd_clean_data.csv")
origin_pollution_data <- origin_pollution_data[,-1]

pollution_gmm <- read.csv('analysis/pollution/results/results_100k_gmm.csv')
pollution_gmm <- pollution_gmm[(0:(nrow(pollution_gmm)/4+1)) * 3 + 1,]
row.names(pollution_gmm) <- pollution_gmm[,1]
pollution_gmm <- pollution_gmm[,-1]

UniformsFromGPD(data_origin = origin_pollution_data[,1],
                data = pollution_data[,1],
                xi = pollution_gmm[1,1],
                sigma = pollution_gmm[1,2])

unif_pollution_data <- UniformDataFromGPD(dataset_origin = origin_pollution_data,
                   dataset = pollution_data,
                   params = pollution_gmm)

