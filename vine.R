# devtools::install_github("vinecopulib/rvinecopulib")
library(rvinecopulib)

pollution_data <- read.csv("data/clean_pollution_data.csv")
pollution_data <- pollution_data[,-1]
origin_pollution_data <- read.csv("data/non_gpd_clean_data.csv")
origin_pollution_data <- origin_pollution_data[,-1]

pollution_gmm <- read.csv('analysis/pollution/results/results_100k_gmm.csv')
pollution_gmm <- pollution_gmm[(0:(nrow(pollution_gmm)/4+1)) * 3 + 1,]
row.names(pollution_gmm) <- pollution_gmm[,1]
pollution_gmm <- pollution_gmm[,-1]

unif_pollution_data <- UniformDataFromGPD(
  dataset_origin = origin_pollution_data,
  dataset = pollution_data,
  params = pollution_gmm)

full_vine_mbic <- rvinecopulib::vinecop(data = unif_pollution_data,
                                        family_set = c("clayton", "gumbel", "indep"),
                                        trunc_lvl = NA,
                                        selcrit='mbicv',
                                        core=parallel::detectCores()-1)
print(full_vine_mbic)


full_vine_bic <- rvinecopulib::vinecop(data = unif_pollution_data,
                                        family_set = c("clayton", "gumbel", "indep"),
                                        trunc_lvl = NA,
                                        selcrit='bic',
                                        core=parallel::detectCores()-1)
print(full_vine_bic)