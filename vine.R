# devtools::install_github("vinecopulib/rvinecopulib")
library(rvinecopulib)
library(ggraph)
pollution_data <- read.csv("data/clean_pollution_data.csv")
pollution_data <- pollution_data[,-1]
origin_pollution_data <- read.csv("data/non_gpd_clean_data.csv")
origin_pollution_data <- origin_pollution_data[,-1]

pollution_gmm <- read.csv('analysis/pollution/results/results_100k_gmm.csv')
pollution_gmm <- pollution_gmm[(0:(nrow(pollution_gmm)/4+1)) * 3 + 1,]
row.names(pollution_gmm) <- pollution_gmm[,1]
pollution_gmm <- pollution_gmm[,-1]

unif_pollution_data <- UniformFromGPDForMatrix(
  dataset_origin = origin_pollution_data,
  dataset = pollution_data,
  params = pollution_gmm)

vine_config <- list()
vine_config[['family_set']] =  c("clayton", "gumbel", "indep")
vine_config[['trunc_lvl']] = NA
vine_config[['selcrit']] = 'mbicv'
vine_config[['core']] = parallel::detectCores()-1

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

evc <- ExtremeVineCollection(pollution_data, unif_pollution_data, horizons = c(1,2), vine_config = vine_config)
evc$O3[[1]]$structure

vine_tmp <- evc$O3[[1]]
vine_struc <- evc$O3[[1]]$structure
col_number_tmp <- 3
sim_vals <- rep(0, vine_tmp$structure$d)
sim_vals[col_number_tmp] <- 0.98
random_start <- runif(vine_tmp$structure$d)

for(level in 1:(vine_tmp$structure$trunc_lvl-1)){
  already_simulated_vars <- which(sim_vals > 0.0)
  
  links <- ExtremeVineExtractLink(vine_tmp, level)
  links_with_only_one <- apply(
    links,
    MARGIN = 1,
    function(link){sum(as.numeric(already_simulated_vars %in% link))==1})
  
  which_to_choose <- which(links_with_only_one)
  
  links <- links[links_with_only_one,]
  cond_vars <- NA
  if(level > 1){
    cond_vars <- ExtremeVineExtractConditional(vine_tmp, level = level)
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
    bicop_tmp <- vine_tmp$pair_copulas[[level]][[which_to_choose[i]]]
    bicop_data_with_cond <- ExtremeVinePlacingConditional(
      rdm_data = sim_vals[already_simed_vars[i]],
      cond_data = random_start[target_vars[i]],
      cond_on =  which(links[i,]==target_vars[i]))
    val <- rvinecopulib::hbicop(u=bicop_data_with_cond,
                                cond_var = which(links[i,]==target_vars[i]),
                                bicop_tmp,
                                inverse = F)

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
              vine_tmp, c(current_target, cond_v))
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
                                            inverse = F)
              }else{
                # warning('extracted_bicop$level+1 != conditional_level')
              }
            }
          }
        }
        current_links <- ExtremeVineExtractLink(vine_tmp, conditional_level)
        conditional_vals <- ExtremeVineExtractConditional(vine_tmp, conditional_level)

        is_condi_link <- apply(current_links, 1,
                            function(lk){all(c(current_target, current_cond) %in% lk)})
        if(is.null(dim(conditional_vals))){
          conditional_vals <- conditional_vals[is_condi_link]
        }else{
          conditional_vals <- conditional_vals[is_condi_link,]
        }
        
      }
    }
    # val <- val[which(!simed_var_bicop[i,])]
    return(c(target_vars[i], val))},
    c(1,2)))
  sim_vals[tmp[,1]] <- tmp[,2]
  print('sim_vals')
  print(sim_vals)
}
  
sim_vals

# simulating

ExtremeVineExtractConditional(evc$O3[[1]], 4)
