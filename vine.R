# devtools::install_github("vinecopulib/rvinecopulib")
library(rvinecopulib)
library(ggraph)

max_data_length <- 100000
pollution_data <- read.csv("~/GitHub/latent-trawl-pairwise/data/clean_pollution_data.csv")
pollution_data <- pollution_data[1:max_data_length,-1]
origin_pollution_data <- read.csv("~/GitHub/latent-trawl-pairwise/data/non_gpd_clean_data.csv")
origin_pollution_data <- origin_pollution_data[1:max_data_length,-1]

pollution_gmm <- read.csv('~/GitHub/latent-trawl-pairwise/analysis/pollution/results/results_100k_gmm.csv')
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

ExtremeVine(dataset = pollution_data,
            uniform_dataset = unif_pollution_data,
            col_number = 2,
            horizon = 1,
            vine_config = vine_config)

evc <- ExtremeVineCollection(
  dataset = pollution_data,
  uniform_dataset = unif_pollution_data,
  horizons=c(1,2),
  vine_config = vine_config
)


full_vine_mbic <- rvinecopulib::vinecop(data = unif_pollution_data,
                                        family_set = c("clayton", "gumbel", "indep"),
                                        trunc_lvl = NA,
                                        selcrit='mbicv',
                                        core=parallel::detectCores()-1)
print(full_vine_mbic)
full_vine_mbic$structure$struct_array
plot(full_vine_mbic)

full_vine_bic <- rvinecopulib::vinecop(data = unif_pollution_data,
                                        family_set = c("clayton", "gumbel", "indep"),
                                        trunc_lvl = NA,
                                        selcrit='bic',
                                        core=parallel::detectCores()-1)
print(full_vine_bic)