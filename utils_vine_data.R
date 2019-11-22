
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

FilterExtremeIndex <- function(dataset, col_number, horizon){
  n_elems <- nrow(dataset)
  index_pick <- which(dataset[,col_number] > 0.0)
  index_pick <- index_pick[which(index_pick <= n_elems-horizon)]
  index_pick <- matrix(
    rep(index_pick,each=ncol(dataset)),
    ncol=ncol(dataset),
    byrow=TRUE)
  index_pick[,-col_number] <- index_pick[,-col_number]+horizon
  
  return(index_pick)
}

ExtremeVineData <- function(dataset, uniform_dataset, col_number, horizon){
  index_pick <- FilterExtremeIndex(dataset, col_number, horizon)
  xvine_data <- vapply(
    1:ncol(dataset),
    function(i){uniform_dataset[index_pick[,i],i]},
    index_pick[,1])
  return(xvine_data)
}

ExtremeVineFit <- function(dataset, uniform_dataset, col_number, horizon, vine_config){
  xvine_data <- ExtremeVineData(dataset, uniform_dataset, col_number, horizon)
  time_before <- Sys.time()
  vine_fit <- rvinecopulib::vinecop(
    data = xvine_data,
    family_set = vine_config[['family_set']],
    trunc_lvl = vine_config[['trunc_lvl']],
    selcrit=vine_config[['selcrit']],
    core=vine_config[['core']])
  print(Sys.time() - time_before)
  return(vine_fit)
}

ExtremeVineCollection <- function(dataset, uniform_dataset, horizons, vine_config){
  # Structure: colnames -> horizons
  n_col <- ncol(dataset)
  result_vines <- lapply(
    1:n_col,
    function(col_number){
        return(lapply(
          horizons,
          function(horizon){return(ExtremeVineFit(dataset, uniform_dataset, col_number, horizon, vine_config))}
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


ExtremeVineExtractConditional <- function(vine, level){
  # in rows
  vine_structure <- vine$structure
  n_vars <- vine_structure$d
  origin <- vine_structure$order[1:(n_vars-level)]
  if(level <= 1){return(numeric())}
  else{
    cond <- rvinecopulib::as_rvine_matrix(vine_structure)[1:(level-1),1:(n_vars-level)]
    return(t(cond))
  }
}

ExtremeVineConditionalSimulation <- function(vine, value, col_number){
  n_vars <- evc$O3[[1]]$structure$d
  
  sim_vals <- rep(0, n_vars)
  sim_vals[i] <- value
  # while(prod(sim_vals) == 0.0){
  #     
  # }  
}



