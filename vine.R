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
sim_vals[col_number_tmp] <- 0.35

for(level in 1:vine_tmp$structure$trunc_lvl){
  already_simulated <- which(sim_vals > 0.0)
  links <- as.matrix(ExtremeVineExtractLink(vine_tmp, level))
  print(links)
  links_with_only_one <- apply(
    links,
    MARGIN = 1,
    function(link){sum(as.numeric(already_simulated %in% link))==1})
  which_to_choose <- which(links_with_only_one)
  links <- links[links_with_only_one,]
  simed_var_bicop <- t(apply(
    links,
    MARGIN = 1,
    function(link){
      return(link %in% already_simulated)
    }))
  
  already_simed <- links[simed_var_bicop]
  print(already_simed)
  to_sim <- links[!simed_var_bicop]
  print(to_sim)
  print(nrow(links))
  print(which_to_choose)
  tmp <- t(vapply(1:nrow(links),
                function(i){
    bicop_tmp <- vine_tmp$pair_copulas[[level]][[which_to_choose[i]]]
    val <- rvinecopulib::hbicop(u=sim_vals[already_simed],
                                cond_var = which(simed_var_bicop[i,]),
                                bicop_tmp)
    print(val)
    val <- val[which(!simed_var_bicop[i,])]
    return(c(to_sim[i], val))},
    c(1,2)))
  print(tmp)
}
  
  # simulating

ExtremeVineExtractConditional(evc$O3[[1]], 4)
