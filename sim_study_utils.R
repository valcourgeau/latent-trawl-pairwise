
EVTrawlFit <- function(data, depth, parametrisation='standard', type='exp', parallel=F){
  custom_mle_results <- CustomMarginalMLE(data)
  kappa <- GetKappa(data = data, params = custom_mle_results, parametrisation = 'standard')
  params <- c(custom_mle_results, kappa)

  if(parametrisation == 'noven'){
    params <- ParametrisationTranslator(params, parametrisation = 'standard')
  }
  
  trawl_pl_basic <- TrawlPL(data = data, depth = depth, parametrisation = parametrisation, type = type, parallel = parallel)
  trawl_pl_restricted <- function(trawl_params){
    return(trawl_pl_basic(c(params, trawl_params)))
  }
  
  trawl_cfg <- GetTrawlParamsConfig(type)
  
  # rdm start
  #TODO add parallel L-BFGS-B
  trawl_res <- stats::optim(fn = trawl_pl_restricted,
                            par = vapply(1:trawl_cfg$n_params, function(i){
                              runif(n = 1,
                                    min = trawl_cfg$lower[i],
                                    max = trawl_cfg$upper[i])
                              }, 1.),
                            lower = trawl_cfg$lower,
                            upper = trawl_cfg$upper,
                            method = 'L-BFGS-B',
                            control = list(trace=3))
  
  return(c(params, trawl_res$par))
}


