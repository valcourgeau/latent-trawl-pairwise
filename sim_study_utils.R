
EVTrawlFit <- function(data, depth, method, parametrisation='standard', type='exp',
                       parallel=F, bounds='config', ...){
  # method 'PL' or 'GMM'
 
  trawl_cfg <- GetTrawlParamsConfig(type)
  # rdm start
  #TODO add parallel L-BFGS-B
  if(method == 'PL'){
    custom_mle_results <- CustomMarginalMLE(data, parametrisation)
    kappa <- GetKappa(data = data, params = custom_mle_results, parametrisation = 'standard')
    params <- c(custom_mle_results, kappa)
    
    if(parametrisation == 'noven'){
      params <- ParametrisationTranslator(params, parametrisation = 'standard')
    }
    

    trawl_pl_basic <- TrawlPL(data = data, depth = depth, parametrisation = parametrisation,
                              type = type, parallel = parallel)
    trawl_pl_restricted <- function(trawl_params){
      return(trawl_pl_basic(c(params, trawl_params)))
    }
    if(bounds == 'config'){
      for(i in 1:100){
        print(parallel)
      }
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
    }else{
      trawl_p <- DupuisSimplified(data_u = data, ...)[1:(4+trawl_cfg$n_params-1)]
      trawl_res <- stats::optim(fn = trawl_pl_restricted,
                                par = vapply(1:trawl_cfg$n_params, function(i){
                                  runif(n = 1,
                                        min = trawl_cfg$lower[i],
                                        max = trawl_cfg$upper[i])
                                }, 1.),
                                lower = trawl_p*0.8,
                                upper = trawl_p*1.2,
                                method = 'L-BFGS-B',
                                control = list(trace=3))
    }
    return(c(params, trawl_res$par))
  }else{
    if(method == 'GMM'){
      trawl_p <- DupuisSimplified(data_u = data, ...)[1:(4+trawl_cfg$n_params-1)]
      return(trawl_p)
    }
  }
}


SubSampleFit <- function(data, depth, sub_length, method, trials, file_csv,
                         parametrisation='standard', type='exp', parallel=F, seed=42, subfolder='simulation/results/', ...){
  # method 'PL' or 'GMM'
  n <- length(data)
  set.seed(seed)
  start_points <- sample(1:(n-sub_length), size = trials, replace = F)
  results <- matrix(0, nrow=trials, ncol=3+GetTrawlParamsConfig(type)$n_params)
  
  if(parallel){
    cores <- parallel::detectCores(logical = TRUE)
    cl <- parallel::makeCluster(cores)
    parallel::clusterExport(cl, c('CaseSeparator',
                        'CppCaseSeparator',
                        # 'data',
                        'DupuisSimplified',
                        'CustomLikelihood',
                        'CheckAllNonpositive',
                        'TransformationMap',
                        'acf_trawl',
                        'ParametrisationTranslator',
                        'CaseOneOne',
                        'CaseOneZero',
                        'CaseZeroZero',
                        'CppCaseOneOne',
                        'CppCaseOneZero',
                        'CppCaseZeroZero',
                        'CheckAllPositive',
                        'StandTrawlTerms',
                        'EVTrawlFit',
                        'CustomMarginalMLE',
                        'GetKappa',
                        'TrawlPL',
                        'GetTrawlParamsConfig',
                        'GetTrawlEnv',
                        'TrawlPLFunctional',
                        'TrawlPLStandard',
                        'TrawlPLNoven',
                        'GetTrawlFunctions',
                        'PairPDFConstructor',
                        'PLConstructor',
                        GetTrawlEnvsList()))
    # parallel::clusterEvalQ(cl, library(zeallot))
    # parallel::clusterEvalQ(cl, library(compiler))
    
    sub_sample_time <- Sys.time()
    results <- parallel::parLapply(X = start_points,
                        cl = cl,
                        fun = function(start_pt){
                          EVTrawlFit(data = data[start_pt:(start_pt+sub_length)],
                                     depth = depth,
                                     parametrisation = parametrisation,
                                     method = method,
                                     type = type,
                                     parallel = F,
                                     ...)
                        }, chunk.size = 5)
    print(Sys.time() - sub_sample_time)
    parallel::stopCluster(cl)
    results <- matrix(unlist(results), ncol=length(results[[1]]), byrow = T)
  }else{
    print('No parallel trials but parallel PL.')
    sub_sample_time <- Sys.time()
    results<- t(vapply(start_points,
                       FUN = function(start_pt){
                         tmp_pr <- EVTrawlFit(data = data[start_pt:(start_pt+sub_length)],
                                    depth = depth,
                                    parametrisation = parametrisation,
                                    type = type,
                                    method = method,
                                    parallel = T,
                                    ...)
                         print(tmp_pr)
                       },
                       FUN.VALUE = rep(0, ncol(results))))
    print(Sys.time() - sub_sample_time)
  }
  
  results <- rbind(c(n, depth, sub_length, trials, rep(0, GetTrawlParamsConfig(type)$n_params-1)), results)
  write.csv(results, file = paste(subfolder, file_csv, sep = ''))
  return(results)
}

