
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

ExtremeVineData <- function(dataset, uniform_dataset, col_number, horizon, rescaling=F){
  assertthat::are_equal(nrow(dataset), nrow(uniform_dataset))
  assertthat::are_equal(ncol(dataset), ncol(uniform_dataset))
  
  index_pick <- FilterExtremeIndex(dataset, col_number, horizon)
  xvine_data <- vapply(
    1:ncol(dataset),
    function(i){uniform_dataset[index_pick[,i],i]},
    index_pick[,1])
  
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
  return(list(xvine_data=xvine_data, quantiles=quantiles, quantile_values=quantile_values))
}

ExtremeVineFit <- function(dataset, uniform_dataset, col_number, horizon, vine_config, rescaling=F){
  xvine_data <- ExtremeVineData(dataset, uniform_dataset, col_number, horizon, rescaling)
  
  time_before <- Sys.time()
  cat('Starting fit on column', col_number, '\n')
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
          function(horizon){return(ExtremeVineFit(dataset, uniform_dataset, col_number, horizon, vine_config, rescaling))}
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

# handle when one computation has to be handle before the other
ExtremeVineCondSimSingle <- function(vine, col_number, value){
  sim_vals <- rep(0, vine$structure$d)
  sim_vals[col_number] <- value
  random_start <- runif(vine$structure$d) # uncorrelated unif(0,1)
  
  for(level in 1:(vine$structure$trunc_lvl-1)){
    links <- ExtremeVineExtractLink(vine, level)
    already_simulated_vars <- which(sim_vals > 0.0)
    links_with_only_one <- apply(
      links,
      MARGIN = 1,
      function(link){sum(as.numeric(already_simulated_vars %in% link))==1})
    
    while(any(links_with_only_one)){
      which_to_choose <- which(links_with_only_one)
      
      links <- links[links_with_only_one,]
      cond_vars <- NA
      if(level > 1){
        cond_vars <- ExtremeVineExtractConditional(vine, level = level)
        is_link_usable_cond <- apply(cond_vars,
                                     MARGIN = 1,
                                     FUN = function(lk_cond_var){
                                       all(lk_cond_var %in% already_simulated_vars)
                                     })
        cond_vars <- cond_vars[is_link_usable_cond,]
        
        if(is.null(dim(links))){
          links <- links[is_link_usable_cond[links_with_only_one]]
        }else{
          links <- links[is_link_usable_cond[links_with_only_one],]
        }
        which_to_choose <- which_to_choose[is_link_usable_cond]
      }
      target_vars <- apply(
        as.matrix(links),
        MARGIN = 1,
        function(link){!(link %in% already_simulated_vars)})
      target_vars <- t(target_vars)
      already_simed_vars <- links[!target_vars]
      target_vars <- links[target_vars]
      
      if(is.null(dim(links))){
        links <- matrix(links, ncol=2)
      }
      
      if(all(vapply(which_to_choose, is.na, FUN.VALUE = NA))){break}
      tmp <- t(vapply(1:nrow(links),
                      function(i){
                        bicop_tmp <- vine$pair_copulas[[level]][[which_to_choose[i]]]
                        bicop_data_with_cond <- ExtremeVinePlacingConditional(
                          rdm_data = sim_vals[already_simed_vars[i]],
                          cond_data = random_start[target_vars[i]],
                          cond_on =  which(links[i,]==target_vars[i])) # first level of interaction, we use uniforms as per Bevacqua et al. (2017)
                        val <- rvinecopulib::hbicop(u=bicop_data_with_cond,
                                                    cond_var = which(links[i,]==target_vars[i]),
                                                    bicop_tmp,
                                                    inverse = T)
                        
                        # we know that cond_vars variables have been simulated
                        if(level > 1){
                          conditional_level <- level
                          current_target <- target_vars[i]
                          current_cond <- NA
                          if(is.null(dim(cond_vars))){
                            conditional_vals <- cond_vars
                          }else{
                            conditional_vals <- cond_vars[i,]
                          }
                          
                          while(conditional_level > 1){
                            more_than_one <- F
                            for(cond_v in conditional_vals){
                              if(!more_than_one){
                                extracted_bicop <- ExtremeVineExtractBicop(
                                  vine, c(current_target, cond_v))
                                current_cond <- cond_v
                                # checking if the link does exist
                                if(!is.na(extracted_bicop$level)){
                                  if(extracted_bicop$level+1 == conditional_level){
                                    more_than_one <- TRUE
                                    conditional_level <- extracted_bicop$level
                                    
                                    assertthat::assert_that(sim_vals[cond_v] > 0.0)
                                    val <- rvinecopulib::hbicop(u=c(val, sim_vals[cond_v]),
                                                                cond_var = 2,
                                                                extracted_bicop$bicop,
                                                                inverse = T)
                                  }else{
                                    # warning('extracted_bicop$level+1 != conditional_level')
                                  }
                                }
                              }
                            }
                            current_links <- ExtremeVineExtractLink(vine, conditional_level)
                            conditional_vals <- ExtremeVineExtractConditional(vine, conditional_level)
                            
                            is_condi_link <- apply(current_links, 1,
                                                   function(lk){all(c(current_target, current_cond) %in% lk)})
                            if(is.null(dim(conditional_vals))){
                              conditional_vals <- conditional_vals[is_condi_link]
                            }else{
                              conditional_vals <- conditional_vals[is_condi_link,]
                            }
                            
                          }
                        }
                        return(c(target_vars[i], val))},
                      c(1,2)))
      sim_vals[tmp[,1]] <- tmp[,2]
      
      # updating links_with_only_one which might have unlock the access to a variable
      already_simulated_vars <- which(sim_vals > 0.0)
      links <- ExtremeVineExtractLink(vine, level)
      links_with_only_one <- apply(
        links,
        MARGIN = 1,
        function(link){sum(as.numeric(already_simulated_vars %in% link))==1})
    }
  }
  
  return(sim_vals)
}

ExtremeVineConditionalSimulation <- function(vine, col_number, value, n, seed=42){
  set.seed(42)
  return(
    t(
      vapply(1:n, FUN = function(i){ExtremeVineCondSimSingle(vine, col_number, value)}, rep(0, vine_tmp$structure$d))
      )
    )
}

ExtremeVinePlacingConditional <- function(rdm_data, cond_data, cond_on){
  if(cond_on == 1){
    return(c(cond_data, rdm_data))
  }else if(cond_on == 2){
    return(c(rdm_data, cond_data))
  }else{
    warning('cond_on should in c(1,2).')
  }
}



