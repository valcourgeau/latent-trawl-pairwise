ExtremeCompleteCDF <- function(data_origin, data, xi, sigma){
  positive_index <- which(data > 0.0)  
  p_zero <- 1 - length(positive_index) / length(data)
  origin_quantile <- quantile(data_origin, p_zero)
  function(x){
    return(evir::pgpd(q = x,
                      xi = xi,
                      beta = sigma))
  }
  
  print(summary(pmax(data_origin-origin_quantile, 0)))
  print(summary(data))
}

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

UniformFromGPDForMatrix <- function(dataset_origin, dataset, params){
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

UniformsFromGPDTestData <- function(data_origin, data, test_data_origin, test_data, xi, sigma){
  # take oirinal data and zero/GPD-distributed data
  # and transform into uniform(0,1)
  positive_index <- which(data > 0.0)  
  p_zero <- 1 - length(positive_index) / length(data)
  ecdf_before_threshold <- ecdf(data_origin[-positive_index])
  
  test_positive_index <- which(test_data > 0)
  test_data[test_positive_index] <- p_zero +  (1-p_zero) * 
    evir::pgpd(q = test_data[test_positive_index],
               xi = xi,
               beta = sigma)
  test_data[-test_positive_index] <- p_zero * ecdf_before_threshold(test_data_origin[-test_positive_index])
  
  return(test_data)
}

UniformFromGPDForMatrixTestData <- function(dataset_origin, dataset, test_dataset_origin, test_dataset, params){
  # datasets are organised in cols
  # params in rows (xi, sigma)
  return(vapply(
    1:ncol(dataset_origin),
    FUN = function(i){
      UniformsFromGPDTestData(data_origin = dataset_origin[,i],
                      data = dataset[,i],
                      test_data_origin = test_dataset_origin[,i],
                      test_data = test_dataset[,i],
                      xi = params[i,1],
                      sigma = params[i,2])
    },
    rep(0, nrow(test_dataset_origin))
  ))
}

FilterExtremeIndex <- function(dataset, col_number, horizon){
  n_elems <- nrow(dataset)
  index_pick <- which(dataset[,col_number] > 0.0)
  index_pick <- index_pick[which(index_pick <= n_elems-horizon)]
  index_origin <- index_pick
  index_pick <- index_pick+horizon
  index_pick <- matrix(
    rep(index_pick,each=ncol(dataset)),
    ncol=ncol(dataset),
    byrow=TRUE)
  
 
  index_pick <- cbind(index_pick, index_origin)
  return(index_pick)
}

ExtremeVineData <- function(dataset, uniform_dataset, col_number, horizon, rescaling=F){
  # concatenate the variable at time t as the last column
  assertthat::are_equal(nrow(dataset), nrow(uniform_dataset))
  assertthat::are_equal(ncol(dataset), ncol(uniform_dataset))
  
  index_pick <- FilterExtremeIndex(dataset, col_number, horizon)
  # t + horizon data
  xvine_data <- vapply(
    1:ncol(dataset),
    function(i){uniform_dataset[index_pick[,i],i]},
    index_pick[,1])
  
  # adding origin (at time t) data
  xvine_data <- cbind(
    xvine_data,
    uniform_dataset[index_pick[,ncol(dataset)+1],col_number]
  )
  
  ecdf_rescaling <- apply(xvine_data, MARGIN = 2,
                          function(x){ecdf(x)})
  
  if(rescaling){
    xvine_data <- apply(xvine_data, MARGIN = 2,
                        function(x){ecdf_tmp <- ecdf(x); return(ecdf_tmp(x))})
  }
  
  extreme_data <- vapply(
    1:ncol(dataset),
    function(i){dataset[index_pick[,i],i]},
    index_pick[,1])
  xvne_proba_zero <- 1-apply(extreme_data > 0, 2, mean)
  
  quantiles <- xvne_proba_zero
  quantile_values <- vapply(1:ncol(uniform_dataset),
                            function(i){quantile(xvine_data[,i], quantiles[i])}, 1.0)
  return(list(xvine_data=xvine_data, quantiles=quantiles, 
              quantile_values=quantile_values, ecdf_rescaling=ecdf_rescaling))
}

# ExtremeVineTestData <- function(dataset, uniform_dataset, test_dataset, test_uniform_dataset, horizon, col_number, rescaling=F){
#   # concatenate the variable at time t as the last column
#   unif_extract <- uniform_dataset[which(dataset[,col_number]>0), col_number]
#   base_ecdf <- ecdf(unif_extract)
#   test_unif_extract <- test_uniform_dataset[which(test_dataset[,col_number]>0), col_number]
#   
#   if(rescaling){
#     return(base_ecdf(test_unif_extract))  
#   }else{
#     return(test_unif_extract)
#   }
# }


ExtremeVinePredictData <- function(dataset, col_number, horizon){
  # concatenate the variable at time t as the last column
  
  index_pick <- FilterExtremeIndex(dataset, col_number, horizon)
  # t + horizon data
  xvine_data <- vapply(
    1:ncol(dataset),
    function(i){dataset[index_pick[,i],i]>0},
    index_pick[,1])
  
  # adding origin (at time t) data
  xvine_data <- cbind(
    xvine_data,
    dataset[index_pick[,ncol(dataset)+1],col_number] > 0
  )
  
  return(xvine_data)
}

ExtremeVineTestData <- function(dataset, test_dataset, uniform_dataset, test_uniform_dataset, col_number, horizon, rescaling=F){
  # concatenate the variable at time t as the last column
  assertthat::are_equal(nrow(dataset), nrow(uniform_dataset))
  assertthat::are_equal(ncol(dataset), ncol(uniform_dataset))
  
  assertthat::are_equal(nrow(test_dataset), nrow(test_uniform_dataset))
  assertthat::are_equal(ncol(test_dataset), ncol(test_uniform_dataset))
  
  index_pick <- FilterExtremeIndex(dataset, col_number, horizon)
  index_pick_test <- FilterExtremeIndex(test_dataset, col_number, horizon)
  
  # t + horizon data
  xvine_data <- vapply(
    1:ncol(dataset),
    function(i){uniform_dataset[index_pick[,i],i]},
    index_pick[,1])
  
  # adding origin (at time t) data
  xvine_data <- cbind(
    xvine_data,
    uniform_dataset[index_pick[,ncol(dataset)+1],col_number]
  )
  
  pred_data <- vapply(
    1:ncol(dataset),
    function(i){dataset[index_pick[,i],i] > 0},
    index_pick[,1])
  
  # adding origin (at time t) data
  pred_data <- cbind(
    pred_data,
    dataset[index_pick[,ncol(dataset)+1],col_number] > 0
  )
  
  # same for test
  xvine_test_data <- vapply(
    1:ncol(dataset),
    function(i){test_uniform_dataset[index_pick_test[,i],i]},
    index_pick_test[,1])
  
  xvine_test_data <- cbind(
    xvine_test_data,
    test_uniform_dataset[index_pick_test[,ncol(dataset)+1],col_number]
  )
  
  pred_test_data <- vapply(
    1:ncol(dataset),
    function(i){test_dataset[index_pick_test[,i],i]>0},
    index_pick_test[,1])
  
  pred_test_data <- cbind(
    pred_test_data,
    test_dataset[index_pick_test[,ncol(dataset)+1],col_number]>0
  )
  
  ecdf_rescaling <- apply(xvine_data, MARGIN = 2,
                      function(x){ecdf(x)})
  if(rescaling){
    xvine_data <- apply(xvine_data, MARGIN = 2,
                        function(x){ecdf_tmp <- ecdf(x); return(ecdf_tmp(x))})
    xvine_test_data <- vapply(1:ncol(xvine_test_data), function(i){ecdf_rescaling[[i]](xvine_test_data[,i])}, rep(0, nrow(xvine_test_data)))
  }
  
  extreme_data <- vapply(
    1:ncol(dataset),
    function(i){dataset[index_pick[,i],i]},
    index_pick[,1])
  xvne_proba_zero <- 1-apply(extreme_data > 0, 2, mean)
  
  quantiles_zero <- xvne_proba_zero
  quantile_values <- vapply(1:ncol(uniform_dataset),
                            function(i){quantile(xvine_data[,i], quantiles_zero[i])}, 1.0)
  return(list(xvine_data=xvine_data, quantiles=quantiles_zero, quantile_values=quantile_values,
              xvine_test_data=xvine_test_data, ecdf_rescaling=ecdf_rescaling,
              pred_data=pred_data, pred_test_data=pred_test_data))
}


ExtremeVineFit <- function(dataset, uniform_dataset, col_number, horizon, vine_config, rescaling=F){
  xvine_data <- ExtremeVineData(dataset, uniform_dataset, col_number, horizon, rescaling)
  time_before <- Sys.time()
  cat('Starting fit on column', col_number, 'with horizon', horizon, '\n')
  target_name <- colnames(dataset)[col_number]
  colnames(xvine_data$xvine_data) <- c(colnames(dataset), paste('ex-ante', target_name))
  vine_fit <- rvinecopulib::vinecop(
    data = xvine_data$xvine_data,
    family_set = vine_config[['family_set']],
    trunc_lvl = vine_config[['trunc_lvl']],
    selcrit = vine_config[['selcrit']],
    core = vine_config[['core']],
    show_trace = vine_config[['show_trace']])
  print(Sys.time() - time_before)
  return(list(vine_fit=vine_fit, quantiles=xvine_data$quantiles, quantile_values=xvine_data$quantile_values))
}

ExtremeVineCollection <- function(dataset, uniform_dataset, horizons, vine_config, rescaling=F){
  # Structure: colnames -> horizons
  n_col <- ncol(dataset)
  result_vines <- lapply(
    1:n_col,
    function(col_number){
        return(lapply(
          horizons,
          function(horizon){
            return(
              ExtremeVineFit(dataset, uniform_dataset, col_number, horizon, vine_config, rescaling)
            )
          }
        )
      )
    }
  )
  names(result_vines) <- colnames(dataset)
  return(result_vines)
}

ExtremeVineSimulation <- function(vine_collection, n){
  #vine_collection: vars -> horizon
  return(
    lapply(vine_collection,
         FUN = function(vine_per_var){
           lapply(vine_per_var,
                  function(vine){
                    rvinecopulib::rvinecop(n = n, vine=vine)
                  })
         }))
}

ExtremeVineAIC <- function(vine_collection){
  return(
    lapply(vine_collection,
           FUN = function(vine_per_var){
             lapply(vine_per_var,
                    function(vine){
                      return(2.0 * (vine$npars - vine$loglik))
                    })
           }))
}

ExtremeVineBIC <- function(vine_collection, n){
  return(
    lapply(vine_collection,
           FUN = function(vine_per_var){
             lapply(vine_per_var,
                    function(vine){
                      return(2.0 * (log(n)*vine$npars - vine$loglik))
                    })
           }))
}

ExtremeVineMBIC <- function(vine_collection){
  return(
    lapply(vine_collection,
           FUN = function(vine_per_var){
             lapply(vine_per_var,
                    function(vine){
                      return(rvinecopulib::mBICV(vine))
                    })
           }))
}

ExtremeVineExtractLink <- function(vine, level){
  # in rows
  vine_structure <- vine$structure
  n_vars <- vine_structure$d
  origin <- vine_structure$order[1:(n_vars-level)]
  target <- rvinecopulib::as_rvine_matrix(vine_structure)[level,1:(n_vars-level)]
  return(cbind(origin, target))
}

ExtremeVineExtractBicop <- function(vine, link){
  #returns list with level, link number and bicop
  assertthat::assert_that(!assertthat::are_equal(link[1], link[2]), msg=paste('Link nodes should be different. Received', link[1], link[2]))
  vine_structure <- vine$structure
  n_vars <- vine_structure$d
  list_links <- lapply(1:n_vars, function(lvl){ExtremeVineExtractLink(vine, lvl)})
  is_this_in_level <- lapply(list_links,
                             function(list_lks){
                               apply(list_lks, MARGIN = 1, function(lk){all(link %in% lk)})
                              })
  any_at_each_level <- vapply(is_this_in_level, any, T)
  if(any(any_at_each_level)){
    link_level <- which(any_at_each_level)
    link_number <- which(is_this_in_level[[link_level]])
    link_bicop <- vine$pair_copulas[[link_level]][[link_number]]
  }else{
    link_level <- NA
    link_number <- NA
    link_bicop <- NA
  }
  
  return(list(level=link_level, number=link_number, bicop=link_bicop))
}


ExtremeVineExtractConditional <- function(vine, level){
  # in rows
  vine_structure <- vine$structure
  n_vars <- vine_structure$d
  origin <- vine_structure$order[1:(n_vars-level)]
  if(level <= 1){return(numeric())}
  else{
    cond <- rvinecopulib::as_rvine_matrix(vine_structure)[1:(level-1),1:(n_vars-level)]
    if(level == 2){
      return(as.matrix(cond))
    }
    return(t(cond))
  }
}




