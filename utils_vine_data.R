
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
