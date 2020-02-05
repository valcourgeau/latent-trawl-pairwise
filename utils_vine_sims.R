ExtremeVinePlacingConditional <- function(rdm_data, cond_data, cond_on){
  if(cond_on == 1){
    return(c(cond_data, rdm_data))
  }else if(cond_on == 2){
    return(c(rdm_data, cond_data))
  }else{
    warning('cond_on should in c(1,2).')
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
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  return(
    t(
      vapply(1:n, FUN = function(i){ExtremeVineCondSimSingle(vine, col_number, value)}, rep(0, vine_tmp$structure$d))
    )
  )
}

ExtremeVineConditionalPredict <- function(vine, quantile_values, col_number, values, n, seed=42){
  n_vars <- vine$structure$d
  predict_vals <-  vapply(
    1:length(values), 
    function(i){ExtremeVineConditionalSimulation(vine, col_number, value = values[i], n = n, seed = seed) %>% as.matrix},
    matrix(0, ncol = n_vars, nrow=n))
  results <- list(
    mean_pred =  t(apply(predict_vals, c(2,3), mean)),
    sd_pred =  t(apply(predict_vals, c(2,3), sd)),
    sd_pred_sample =  t(apply(predict_vals, c(2,3), sd))/sqrt(n),
    quantile_values = quantile_values
  )
  
  results[['pred']] <- t(apply(results[['mean_pred']][,-col_number], 1, function(x){x>quantile_values}))
  
  return(results)
}

ExtremeVineConditionalIndicatorPredict <- function(vine, quantile_values, col_number, values, n, seed=42){
  n_vars <- vine$structure$d
  predict_vals <-  vapply(
    1:length(values), 
    function(i){ExtremeVineConditionalSimulation(vine, col_number, value = values[i], n = n, seed = seed) %>% as.matrix},
    matrix(0, ncol = n_vars, nrow=n))
  
  # cut the unnecessary values
  tmp_vals <- apply(
    predict_vals[,-col_number,],
    c(1,3),
    function(x){x>quantile_values})
  tmp_vals <- aperm(tmp_vals, c(2,1,3))
  
  results <- list(
    mean_pred =  t(apply(tmp_vals, c(2,3), mean)),
    sd_pred =  t(apply(tmp_vals, c(2,3), sd)),
    sd_pred_sample =  t(apply(tmp_vals, c(2,3), sd))/sqrt(n),
    quantile_values = quantile_values
  )
  results[['pred']] <- t(apply(results[['mean_pred']][,-col_number], 1, function(x){x>quantile_values}))
  return(results)
}

ExtremeVineMarginalCondSim <- function(xi, sigma, u, n){
  # sim from X | X > u (~ GPD(u; xi, sigma+u*xi)) where X | X > 0 ~ GPD(xi, sigma)
  # no reset of seed
  return(u + evir::rgpd(n = n, xi = xi, mu = 0, beta = sigma + u*xi))
}

ExtremeVineTRON <- function(vine, quantile_values, extreme_quantile, col_number, cond_threshold, xi, sigma, n, ecdf_rescaling, seed=42){
  # cond_theshold is in the original scale: threshold for extreme data (default: 0.0)
  # Computes the (non-conditional) TRONs
  if(!is.null(seed)){set.seed(seed)}
  
  cond_values <- ExtremeVineMarginalCondSim(xi = xi, sigma = sigma, u = cond_threshold, n = n)
  cond_values <- extreme_quantile + (1-extreme_quantile) * evir::pgpd(cond_values, xi = xi, mu = 0.0, beta = sigma)
  cond_values <- ecdf_rescaling(cond_values)
  # return one sims for each value incond_values
  
  values_from_vine <- t(
      vapply(cond_values,
             FUN = function(cond_val){
               ExtremeVineConditionalSimulation(vine = vine, col_number = col_number, value = cond_val, n = 1, seed = NULL)
             },
             rep(0, vine_tmp$structure$d))
    )
  above_threshold <- t(apply(values_from_vine[,-col_number], 1, function(x){x>quantile_values}))
  
  result <- list()
  result[['TRON']] <- apply(above_threshold, 2, mean)
  result[['sd']] <- apply(above_threshold, 2, sd)
  result[['sd_sample']] <- result[['sd']] / sqrt(n)
  result[['quantile_values']] <- quantile_values
  
  return(result)
}
