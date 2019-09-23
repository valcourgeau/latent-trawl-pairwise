
EVTrawlFit <- function(data, parametrisation='standard'){
  custom_mle_results <- CustomMarginalMLE(data)
  kappa <- GetKappa(data = data, params = custom_mle_results, parametrisation = parametrisation)
  return(c(custom_mle_results, kappa))
}

EVTrawlFit(o)
